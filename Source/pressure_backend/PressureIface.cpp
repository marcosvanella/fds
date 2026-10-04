// Selector and the full solve_pressure() wrapper (common layer around a backend).
#include "PressureIface.H"
#include "PressureBackend.H"
#include "CommonLayer.H"
#include "Composite.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <sstream>

namespace pb {

using namespace amrex;

namespace {
bool has_nonzero (iMultiFab const* m)
{
    if (!m) { return false; }
    return m->max(0) != 0 || m->min(0) != 0;
}
bool has_zero (iMultiFab const* m)
{
    if (!m) { return false; }
    return m->min(0) == 0;
}

std::string validate (PressureProblem const& p)
{
    if (!p.rhs || !p.phi) { return "rhs and phi are required"; }
    if (p.ba.empty()) { return "empty BoxArray"; }
    if (p.rhs->boxArray() != p.ba || p.phi->boxArray() != p.ba) { return "rhs/phi BoxArray differ from problem BoxArray"; }
    if (!(p.rhs->DistributionMap() == p.dm) || !(p.phi->DistributionMap() == p.dm)) { return "rhs/phi DistributionMapping differ"; }
    if (p.rhs->nComp() < 1 || p.phi->nComp() < 1) { return "rhs/phi need at least one component"; }
    if (!p.geom.Domain().contains(p.ba.minimalBox())) { return "BoxArray outside the geometry domain"; }
    if (p.ba.minimalBox() != p.geom.Domain() || p.ba.numPts() != p.geom.Domain().numPts()) {
        return "BoxArray does not cover the domain exactly (box domain required)";
    }
    if (p.gauge_weight && (p.gauge_weight->boxArray() != p.ba || !(p.gauge_weight->DistributionMap() == p.dm) || p.gauge_weight->nComp() < 1)) {
        return "gauge_weight BoxArray/DistributionMapping differ from the problem or it has no component";
    }
    if (p.gauge_offset && (p.gauge_offset->boxArray() != p.ba || !(p.gauge_offset->DistributionMap() == p.dm) || p.gauge_offset->nComp() < 1)) {
        return "gauge_offset BoxArray/DistributionMapping differ from the problem or it has no component";
    }
    for (int d = 0; d < 3; ++d) {
        if (!p.cell_width[d].empty()) {
            if (static_cast<int>(p.cell_width[d].size()) != p.geom.Domain().length(d)) { return "cell_width size differs from the domain length"; }
            for (auto w : p.cell_width[d]) {
                if (std::abs(w - p.geom.CellSize(d)) > 1.0e-12*p.geom.CellSize(d)) { return "uniform cell_width disagrees with the geometry cell size"; }
            }
        }
    }
    for (int d = 0; d < 3; ++d) {
        const bool lo = p.bc[face_index(d,0)] == BC::Periodic, hi = p.bc[face_index(d,1)] == BC::Periodic;
        if (lo != hi) { return "periodic BC must be set on both faces of a direction"; }
        if (lo != (p.geom.isPeriodic(d) != 0)) { return "Geometry periodicity does not match the BC"; }
    }
    return "";
}
}

Selection select_backend (PressureProblem const& p, BackendKind requested)
{
    Selection s;
    if (!p.levels.empty()) { return select_composite(p, requested); }
    if (p.nlevels != 1) {
        s.message = "composite (multi-level) pressure solve is not built"; return s;
    }
    if (p.cylindrical || !p.geom.IsCartesian()) {
        s.message = "cylindrical (or other non-Cartesian) geometry is not built (FDS CYLINDRICAL scales rows and RHS by the radius factor)"; return s;
    }
    for (int d = 0; d < 3; ++d) {
        auto const& w = p.cell_width[d];
        if (w.empty()) { continue; }
        const auto mm = std::minmax_element(w.begin(), w.end());
        if (*mm.second - *mm.first > 1.0e-12 * std::abs(*mm.second)) {
            s.message = "non-uniform cell widths (stretched mesh) are not built: FFT::Poisson and MLPoisson assume uniform spacing";
            return s;
        }
    }
    if (p.cell_coef_a || p.face_coef_b[0] || p.face_coef_b[1] || p.face_coef_b[2]) {
        s.message = "variable coefficients (masked or non-unit operator) are not built"; return s;
    }
    if (p.component_id) { s.message = "driver-supplied component ids are not built (masked branch)"; return s; }
    if (has_nonzero(p.cell_class)) { s.message = "masked cells (cell_class != 0: obstructed/solid, known-value or pinned cells) are not built (masked branch)"; return s; }
    if (has_zero(p.uncovered)) { s.message = "covered cells (uncovered == 0) are not built"; return s; }
    int nopen = 0;
    for (int f = 0; f < 6; ++f) { nopen += (p.bc[f] == BC::Dirichlet) ? 1 : 0; }
    if (nopen != 0 && nopen != 6) {
        s.message = "single level with mixed open/closed domain faces is not built"; return s;
    }
    s.ok = true;
    s.kind = (requested == BackendKind::Auto) ? BackendKind::FFT : requested;
    return s;
}

PressureResult solve_pressure (PressureProblem const& p, PressureOptions const& o, PressureWorkspace* ws)
{
    PressureResult R;
    auto fail = [&] (Status st, std::string const& msg) {
        R.status = st; R.message = msg;
        if (o.verbose > 0) { Print() << "PRESSURE ERROR (" << to_string(st) << "): " << msg << "\n"; }
        return R;
    };
    // Selector first: "not built" must be reported for unsupported requests whatever else is unset.
    Selection sel = select_backend(p, o.backend);
    if (!sel.ok) { return fail(Status::NotBuilt, sel.message); }
    if (!p.levels.empty()) {
        std::string badc = validate_composite(p);
        if (!badc.empty()) { return fail(Status::InvalidInput, badc); }
        R = solve_composite(p, o, ws);
        if (R.status == Status::NotBuilt || R.status == Status::InvalidInput) {
            if (o.verbose > 0) { Print() << "PRESSURE ERROR (" << to_string(R.status) << "): " << R.message << "\n"; }
        }
        return R;
    }
    std::string bad = validate(p);
    if (!bad.empty()) { return fail(Status::InvalidInput, bad); }

    std::unique_ptr<PressureBackend> be = (sel.kind == BackendKind::FFT) ? make_fft_backend() : make_mlmg_backend();
    R.backend = be->name();

    const Real vol = p.geom.CellSize(0) * p.geom.CellSize(1) * p.geom.CellSize(2);
    ComponentMap cm = label_components(p);

    MultiFab b(p.ba, p.dm, 1, 0);
    MultiFab::Copy(b, *p.rhs, 0, 0, 1, 0);
    if (o.remove_mean) {
        remove_mean(b, cm, p.uncovered, vol, nullptr, p.mean_kind);
    }
    MultiFab work(p.ba, p.dm, 1, 1);
    work.setVal(0.0);
    if (o.use_initial_guess) { MultiFab::Copy(work, *p.phi, 0, 0, 1, 0); }

    R.backend_status = be->solve(p, o, work, b);

    // Gauge, then the true residual of the problem the backend was given (mean-removed b).
    apply_gauge(work, cm, p.uncovered, vol, nullptr, p.gauge_weight, p.gauge_offset);
    if (o.check_residual) {
        ResidualNorms rn = true_residual(p, work, b);
        R.residual_checked = true;
        R.residual_rel2 = rn.rel2; R.residual_relmax = rn.relmax;
        R.residual_ok = (rn.rel2 <= o.residual_tol);
        if (!R.residual_ok) {
            std::ostringstream m;
            m << "true residual ||b-L*H||_2/||b||_2 = " << rn.rel2 << " exceeds " << o.residual_tol;
            R.warnings.push_back(m.str());
        }
    }
    for (auto const& c : cm.comps) {
        if (c.singular && c.removed_rel > o.removed_mean_warn) {
            std::ostringstream m;
            m << "removed mean of component " << c.id << " is " << c.removed_rel << " of rms(b), above " << o.removed_mean_warn;
            R.warnings.push_back(m.str());
        }
    }
    if (o.verbose > 0) { for (auto const& w : R.warnings) { Print() << "PRESSURE WARNING: " << w << "\n"; } }
    R.components = cm.comps;
    if (!R.backend_status.converged) {
        R.status = Status::NotConverged; R.message = "backend did not converge";
        if (o.verbose > 0) { Print() << "PRESSURE WARNING: backend did not converge\n"; }
    }
    MultiFab::Copy(*p.phi, work, 0, 0, 1, 0);
    return R;
}

} // namespace pb
