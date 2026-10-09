// Composite (multi-level) pressure solve: amrex::MLMG over MLPoisson on the whole hierarchy (milestone M2, first
// version). The common layer (mean removal, gauge, true residual, diagnostics) runs from outside on the
// uncovered cells of all levels with one fixed-point scale (ExactSum.H). See frozen/composite-notes.md.
#include "CommonLayer.H"
#include "Composite.H"
#include "ExactSum.H"
#include "PbWorkspaceImpl.H"

#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace pb {

using namespace amrex;

// ---------------------------------------------------------------------------------------------------------------
// Workspace
// ---------------------------------------------------------------------------------------------------------------
PressureWorkspace::PressureWorkspace () = default;
PressureWorkspace::~PressureWorkspace () = default;
PressureWorkspace::PressureWorkspace (PressureWorkspace&&) noexcept = default;
PressureWorkspace& PressureWorkspace::operator= (PressureWorkspace&&) noexcept = default;
#ifdef PB_WITH_HYPRE
bool PressureWorkspace::built () const { return m_impl && (m_impl->nlev > 0 || m_impl->fft || m_impl->hypre); }
#else
bool PressureWorkspace::built () const { return m_impl && (m_impl->nlev > 0 || m_impl->fft); }
#endif
int PressureWorkspace::num_levels () const { return m_impl ? m_impl->nlev : 0; }
bool PressureWorkspace::fft_plan_built () const { return m_impl && m_impl->fft; }
#ifdef PB_WITH_HYPRE
bool PressureWorkspace::hypre_built () const { return m_impl && (m_impl->hypre || m_impl->hsys); }
#else
bool PressureWorkspace::hypre_built () const { return false; }
#endif
void PressureWorkspace::discard () { m_impl.reset(); }
PressureWorkspace::Impl& PressureWorkspace::ensure_impl ()
{
    if (!m_impl) { m_impl = std::make_unique<Impl>(); }
    return *m_impl;
}

namespace {

LinOpBCType lin_bc (BC b)
{
    switch (b) {
    case BC::Neumann: return LinOpBCType::Neumann;
    case BC::Periodic: return LinOpBCType::Periodic;
    case BC::Dirichlet: return LinOpBCType::Dirichlet;
    }
    return LinOpBCType::Neumann;
}

bool has_nonzero_mask (iMultiFab const* m) { return m && (m->max(0) != 0 || m->min(0) != 0); }
bool has_zero_mask (iMultiFab const* m) { return m && m->min(0) == 0; }

int hidden_direction (PressureProblem const& p)
{
    Box const& dom = p.levels[0].geom.Domain();
    for (int d = 0; d < 3; ++d) { if (dom.length(d) == 1) { return d; } }
    return -1;
}

Real cell_volume (Geometry const& g) { return g.CellSize(0)*g.CellSize(1)*g.CellSize(2); }

// Covered mask of level l: the cells under the next finer level coarsened by its ratio are 0.
std::unique_ptr<iMultiFab> make_uncovered (BoxArray const& ba, DistributionMapping const& dm,
                                           BoxArray const* fine_ba, IntVect const& ratio)
{
    auto m = std::make_unique<iMultiFab>(ba, dm, 1, 0);
    m->setVal(1);
    if (fine_ba) {
        BoxArray cba = *fine_ba;
        cba.coarsen(ratio);
        for (MFIter mfi(*m); mfi.isValid(); ++mfi) {
            auto const& a = m->array(mfi);
            for (auto const& is : cba.intersections(mfi.validbox())) {
                amrex::LoopOnCpu(is.second, [&] (int i, int j, int k) { a(i,j,k) = 0; });
            }
        }
    }
    return m;
}

void build_impl (PressureWorkspace::Impl& W, PressureProblem const& p)
{
    const int n = static_cast<int>(p.levels.size());
    W.nlev = n;
    W.bc = p.bc;
    W.hidden = hidden_direction(p);
    W.ba.clear(); W.dm.clear(); W.geom.clear(); W.ratio.clear(); W.unc.clear(); W.vol.clear(); W.nunc.clear();
    for (int l = 0; l < n; ++l) {
        PressureLevel const& L = p.levels[l];
        W.ba.push_back(L.ba); W.dm.push_back(L.dm); W.geom.push_back(L.geom);
        W.ratio.push_back(l == 0 ? IntVect(1) : L.ref_ratio);
    }
    for (int l = 0; l < n; ++l) {
        const bool has_fine = (l + 1 < n);
        W.unc.push_back(make_uncovered(W.ba[l], W.dm[l], has_fine ? &W.ba[l+1] : nullptr, has_fine ? W.ratio[l+1] : IntVect(1)));
        W.vol.push_back(cell_volume(W.geom[l]));
        Long ncov = 0;
        if (has_fine) { BoxArray c = W.ba[l+1]; c.coarsen(W.ratio[l+1]); ncov = c.numPts(); }
        W.nunc.push_back(W.ba[l].numPts() - ncov);
    }
    W.eba = W.ba; W.egeom = W.geom;
    std::array<BC,6> ebc = p.bc;
    if (W.hidden >= 0) {
        const int d = W.hidden;
        const int lo = W.geom[0].Domain().smallEnd(d);
        int od = (d == 0) ? 1 : 0;                       // a non-hidden direction sets the isotropic cell size
        const Real h0 = W.geom[0].CellSize(od);
        int R = 1;
        for (int l = 0; l < n; ++l) {
            if (l > 0) { for (int t = 0; t < 3; ++t) { if (t != d) { R *= W.ratio[l][t]; break; } } }
            const int nyl = W.ext_n * R;
            BoxList bl;
            for (int i = 0; i < static_cast<int>(W.ba[l].size()); ++i) {
                Box b = W.ba[l][i];
                b.setSmall(d, lo); b.setBig(d, lo + nyl - 1);
                bl.push_back(b);
            }
            W.eba[l] = BoxArray(bl);
            Box dom = W.geom[l].Domain();
            dom.setSmall(d, lo); dom.setBig(d, lo + nyl - 1);
            RealBox rb = W.geom[l].ProbDomain();
            Real plo[3] = {rb.lo(0), rb.lo(1), rb.lo(2)}, phi_[3] = {rb.hi(0), rb.hi(1), rb.hi(2)};
            phi_[d] = plo[d] + W.ext_n * h0;
            Array<int,3> per{W.geom[l].isPeriodic(0), W.geom[l].isPeriodic(1), W.geom[l].isPeriodic(2)};
            per[d] = 1;
            W.egeom[l] = Geometry(dom, RealBox({plo[0],plo[1],plo[2]}, {phi_[0],phi_[1],phi_[2]}), CoordSys::cartesian, per);
        }
        ebc[face_index(d,0)] = BC::Periodic; ebc[face_index(d,1)] = BC::Periodic;     // the hidden direction's term is 0 whatever its BC (D-057)
    }
    // Q1: masked composite is assembled HYPRE only — no MLMG hierarchy (C/F interp cannot exclude Solid/Known).
    if (composite_is_masked(p)) {
        W.mlmg.reset();
        W.mlp.reset();
        return;
    }
    LPInfo info;
    Vector<Geometry> geoms(W.egeom.begin(), W.egeom.end());
    Vector<BoxArray> bas(W.eba.begin(), W.eba.end());
    Vector<DistributionMapping> dms(W.dm.begin(), W.dm.end());
    W.mlmg.reset();
    W.mlp = std::make_unique<MLPoisson>(geoms, bas, dms, info);
    W.mlp->setMaxOrder(2);
    Array<LinOpBCType,AMREX_SPACEDIM> lo{lin_bc(ebc[face_index(0,0)]), lin_bc(ebc[face_index(1,0)]), lin_bc(ebc[face_index(2,0)])};
    Array<LinOpBCType,AMREX_SPACEDIM> hi{lin_bc(ebc[face_index(0,1)]), lin_bc(ebc[face_index(1,1)]), lin_bc(ebc[face_index(2,1)])};
    W.mlp->setDomainBC(lo, hi);
    // Homogeneous Dirichlet values (nullptr) on every level. The coarse-fine boundary of level l > 0 is
    // filled by MLMG from level l-1 at every solve.
    for (int l = 0; l < n; ++l) { W.mlp->setLevelBC(l, nullptr); }
    W.mlmg = std::make_unique<MLMG>(*W.mlp);
}

bool impl_matches (PressureWorkspace::Impl const& W, PressureProblem const& p)
{
    const int n = static_cast<int>(p.levels.size());
    if (W.nlev != n || W.bc != p.bc) { return false; }
    for (int l = 0; l < n; ++l) {
        PressureLevel const& L = p.levels[l];
        if (!(W.ba[l] == L.ba) || !(W.dm[l] == L.dm) || !same_geometry(W.geom[l], L.geom)) { return false; }
        if (l > 0 && W.ratio[l] != L.ref_ratio) { return false; }
    }
    return true;
}

} // anonymous

bool PressureWorkspace::matches (PressureProblem const& p) const
{
#ifdef PB_WITH_HYPRE
    if (p.levels.empty()) { return m_impl && ((m_impl->fft && m_impl->fft->plan_matches(p)) || (m_impl->hypre && m_impl->hypre->plan_matches(p))); }
#else
    if (p.levels.empty()) { return m_impl && m_impl->fft && m_impl->fft->plan_matches(p); }
#endif
    return m_impl && m_impl->nlev > 0 && impl_matches(*m_impl, p);
}

Status PressureWorkspace::rebuild (PressureProblem const& p, std::string* message)
{
    if (p.levels.empty()) { return rebuild_single_level(*this, p, message); }
    m_impl.reset();
    mark_regrid();
    auto fail = [&] (Status st, std::string const& m) { if (message) { *message = m; } return st; };
    Selection sel = select_composite(p, BackendKind::Auto);
    if (!sel.ok) { return fail(Status::NotBuilt, sel.message); }
    std::string bad = validate_composite(p);
    if (!bad.empty()) { return fail(Status::InvalidInput, bad); }
    auto w = std::make_unique<Impl>();
    build_impl(*w, p);
    m_impl = std::move(w);
    if (message) { message->clear(); }
    return Status::Ok;
}

// ---------------------------------------------------------------------------------------------------------------
// Selector and validation
// ---------------------------------------------------------------------------------------------------------------
Selection select_composite (PressureProblem const& p, BackendKind requested)
{
    Selection s;
    if (requested == BackendKind::FFT) {
        s.message = "composite (multi-level) solve on the FFT backend is not built: only MLMG and HYPRE do composite"; return s;
    }
    for (auto const& L : p.levels) {
        if (!L.geom.IsCartesian()) {
            s.message = "cylindrical (or other non-Cartesian) geometry is not built (FDS CYLINDRICAL scales rows and RHS by the radius factor)"; return s;
        }
    }
    if (p.cylindrical) {
        s.message = "cylindrical (or other non-Cartesian) geometry is not built (FDS CYLINDRICAL scales rows and RHS by the radius factor)"; return s;
    }
    for (int d = 0; d < 3; ++d) {
        auto const& w = p.cell_width[d];
        if (w.empty()) { continue; }
        const auto mm = std::minmax_element(w.begin(), w.end());
        if (*mm.second - *mm.first > 1.0e-12 * std::abs(*mm.second)) {
            s.message = "non-uniform cell widths (stretched mesh) are not built: MLPoisson assumes uniform spacing"; return s;
        }
    }
    if (p.cell_coef_a || p.face_coef_b[0] || p.face_coef_b[1] || p.face_coef_b[2]) {
        s.message = "variable coefficients (masked or non-unit operator) are not built"; return s;
    }
    if (p.component_id) { s.message = "composite with driver-supplied component ids is not built (masked branch)"; return s; }
    if (p.cell_class || p.known_value) {
        s.message = "PressureProblem::cell_class/known_value are single-level fields; with `levels` set use PressureLevel::cell_class/known_value"; return s;
    }
    if (has_zero_mask(p.uncovered)) { s.message = "composite with a caller-supplied covered-cell mask is not built: covered cells are derived from the level BoxArrays"; return s; }
    const bool masked = composite_is_masked(p);
    if (masked && requested == BackendKind::MLMG) {
        s.message = "composite with masked cells is not built on MLMG (C/F interpolation cannot exclude Solid/Known per stencil; Q1: use assembled HYPRE)"; return s;
    }
    if (masked && requested == BackendKind::FFT) {
        s.message = "composite with masked cells cannot go through FFT"; return s;
    }
#ifndef PB_WITH_HYPRE
    if (masked) { s.message = kHypreNotBuilt; return s; }
#endif
    if (p.gauge_weight || p.gauge_offset) {
        s.message = "PressureProblem::gauge_weight/gauge_offset are single-level fields; with `levels` set use PressureLevel::gauge_weight/gauge_offset"; return s;
    }
    // Ratios and the one-cell direction.
    Box const& dom0 = p.levels[0].geom.Domain();
    int nhidden = 0, hd = -1;
    for (int d = 0; d < 3; ++d) { if (dom0.length(d) == 1) { ++nhidden; hd = d; } }
    if (nhidden > 1) { s.message = "composite with more than one one-cell direction is not built"; return s; }
    if (hd >= 0) {
        // D-057: only a one-cell y is the FDS TWO_D case (its term is dropped, Dirichlet faces act as Neumann: effective_bc).
        if (hd != 1 && (p.bc[face_index(hd,0)] == BC::Dirichlet || p.bc[face_index(hd,1)] == BC::Dirichlet)) {
            s.message = "composite with a Dirichlet face in a one-cell x or z direction is not built (the direction would carry a non-zero operator term)"; return s;
        }
    }
    for (std::size_t l = 1; l < p.levels.size(); ++l) {
        IntVect const& r = p.levels[l].ref_ratio;
        int rr = 0;
        for (int d = 0; d < 3; ++d) {
            if (d == hd) {
                if (r[d] != 1) { s.message = "composite refinement in a one-cell direction is not built (ratio must be 1 there)"; return s; }
                continue;
            }
            if (rr == 0) { rr = r[d]; }
            if (r[d] != rr) { s.message = "composite with anisotropic refinement ratios is not built (one ratio for all refined directions)"; return s; }
        }
        if (rr != 2 && rr != 4) { s.message = "composite refinement ratios other than 2 and 4 are not built"; return s; }
    }
    s.ok = true;
    if (requested == BackendKind::HYPRE) { s.kind = BackendKind::HYPRE; }
    else if (requested == BackendKind::Auto && composite_is_masked(p)) { s.kind = BackendKind::HYPRE; }  // D-092 product: Auto → HYPRE
    else { s.kind = BackendKind::MLMG; }
    return s;
}

std::string validate_composite (PressureProblem const& p, bool structural)
{
    const int n = static_cast<int>(p.levels.size());
    for (int l = 0; l < n; ++l) {
        PressureLevel const& L = p.levels[l];
        const std::string tag = "level " + std::to_string(l) + ": ";
        if (!L.rhs || !L.phi) { return tag + "rhs and phi are required"; }
        if (L.ba.empty()) { return tag + "empty BoxArray"; }
        if (L.dm.size() != L.ba.size()) { return tag + "DistributionMapping does not match the BoxArray"; }
        if (!(L.rhs->boxArray() == L.ba) || !(L.phi->boxArray() == L.ba)) { return tag + "rhs/phi BoxArray differ from the level BoxArray"; }
        if (!(L.rhs->DistributionMap() == L.dm) || !(L.phi->DistributionMap() == L.dm)) { return tag + "rhs/phi DistributionMapping differ from the level"; }
        if (L.rhs->nComp() < 1 || L.phi->nComp() < 1) { return tag + "rhs/phi need at least one component"; }
        if ((L.gauge_weight != nullptr) != (p.levels[0].gauge_weight != nullptr)) { return tag + "gauge_weight must be given on every level or on none"; }
        if ((L.gauge_offset != nullptr) != (p.levels[0].gauge_offset != nullptr)) { return tag + "gauge_offset must be given on every level or on none"; }
        for (MultiFab const* g : {L.gauge_weight, L.gauge_offset, L.known_value}) {
            if (g && (!(g->boxArray() == L.ba) || !(g->DistributionMap() == L.dm) || g->nComp() < 1)) {
                return tag + "gauge_weight/gauge_offset/known_value BoxArray/DistributionMapping differ from the level or they have no component";
            }
        }
        if (L.cell_class) {
            if (L.cell_class->boxArray() != L.ba || !(L.cell_class->DistributionMap() == L.dm) || L.cell_class->nComp() < 1) {
                return tag + "cell_class BoxArray/DistributionMapping differ from the level or it has no component";
            }
            if (L.cell_class->min(0) < CellGas || L.cell_class->max(0) > CellKnown) {
                return tag + "cell_class has a code other than 0 (gas), 1 (solid) or 2 (known value)";
            }
        }
        if (structural) {
            if (!L.geom.Domain().contains(L.ba.minimalBox())) { return tag + "BoxArray outside the level domain"; }
            if (!L.ba.isDisjoint()) { return tag + "BoxArray is not disjoint"; }
        }
        for (int d = 0; d < 3; ++d) {
            const bool lo = p.bc[face_index(d,0)] == BC::Periodic, hi = p.bc[face_index(d,1)] == BC::Periodic;
            if (lo != hi) { return "periodic BC must be set on both faces of a direction"; }
            if (lo != (L.geom.isPeriodic(d) != 0)) { return tag + "Geometry periodicity does not match the BC"; }
        }
        if (!structural) { continue; }     // cheap path with a matching workspace: the layout was validated when the workspace was built
        if (l == 0) {
            if (L.ba.minimalBox() != L.geom.Domain() || L.ba.numPts() != L.geom.Domain().numPts()) {
                return tag + "level 0 BoxArray does not cover the domain exactly";
            }
            for (int d = 0; d < 3; ++d) {
                if (!p.cell_width[d].empty()) {
                    if (static_cast<int>(p.cell_width[d].size()) != L.geom.Domain().length(d)) { return "cell_width size differs from the level 0 domain length"; }
                    for (auto w : p.cell_width[d]) {
                        if (std::abs(w - L.geom.CellSize(d)) > 1.0e-12*L.geom.CellSize(d)) { return "uniform cell_width disagrees with the level 0 geometry cell size"; }
                    }
                }
            }
            continue;
        }
        PressureLevel const& C = p.levels[l-1];
        IntVect const& rr = L.ref_ratio;
        if (rr.min() < 1) { return tag + "ref_ratio must be positive"; }
        if (L.geom.Domain() != amrex::refine(C.geom.Domain(), rr)) { return tag + "domain is not the coarser domain refined by ref_ratio"; }
        for (int d = 0; d < 3; ++d) {
            if (std::abs(L.geom.CellSize(d)*rr[d] - C.geom.CellSize(d)) > 1.0e-10*C.geom.CellSize(d)) { return tag + "cell size is not the coarser cell size divided by ref_ratio"; }
            if (std::abs(L.geom.ProbLo(d) - C.geom.ProbLo(d)) > 1.0e-10*C.geom.ProbLength(d)) { return tag + "problem domain differs from the coarser level's"; }
            if (std::abs(L.geom.ProbHi(d) - C.geom.ProbHi(d)) > 1.0e-10*C.geom.ProbLength(d)) { return tag + "problem domain differs from the coarser level's"; }
        }
        if (!L.ba.coarsenable(rr)) { return tag + "BoxArray is not aligned to ref_ratio (boxes not coarsenable)"; }
        BoxArray cba = L.ba;
        cba.coarsen(rr);
        if (!C.ba.contains(cba)) { return tag + "not nested in the coarser level"; }
        // Proper nesting with one coarse cell of buffer (except at non-periodic domain faces): the coarse cells
        // next to the coarsened level, taken through periodic images, must belong to the coarser BoxArray.
        const Box cdom = C.geom.Domain();
        std::vector<IntVect> shifts{IntVect(0)};
        for (auto const& sv : C.geom.periodicity().shiftIntVect()) { shifts.push_back(sv); }
        for (int i = 0; i < static_cast<int>(cba.size()); ++i) {
            const Box bg = amrex::grow(cba[i], 1);
            for (auto const& sv : shifts) {
                const Box piece = amrex::shift(bg, sv) & cdom;
                if (piece.ok() && !C.ba.contains(piece)) { return tag + "not properly nested: needs one cell of the coarser level around it"; }
            }
        }
    }
    return "";
}

// ---------------------------------------------------------------------------------------------------------------
// Helpers on the hierarchy
// ---------------------------------------------------------------------------------------------------------------
namespace {

using LevelMFs = std::vector<std::unique_ptr<MultiFab>>;

std::vector<SumTerm> terms_of (PressureWorkspace::Impl const& W, LevelMFs const& f, bool volume_weighted)
{
    std::vector<SumTerm> t;
    for (int l = 0; l < W.nlev; ++l) {
        SumTerm s;
        s.mf = f[l].get(); s.comp = 0; s.weight = volume_weighted ? double(W.vol[l]) : 1.0; s.uncovered = W.unc[l].get();
        t.push_back(s);
    }
    return t;
}

double max_abs_uncovered (PressureWorkspace::Impl const& W, LevelMFs const& f)
{
    double m = 0.0;
    for (int l = 0; l < W.nlev; ++l) {
        for (MFIter mfi(*f[l]); mfi.isValid(); ++mfi) {
            auto const& a = f[l]->const_array(mfi);
            auto const& u = W.unc[l]->const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { if (u(i,j,k) != 0) { m = std::max(m, std::abs(double(a(i,j,k)))); } });
        }
    }
    ParallelDescriptor::ReduceRealMax(m);
    return m;
}

LevelMFs squares_of (PressureWorkspace::Impl const& W, LevelMFs const& f)
{
    LevelMFs q;
    for (int l = 0; l < W.nlev; ++l) {
        q.push_back(std::make_unique<MultiFab>(W.ba[l], W.dm[l], 1, 0));
        for (MFIter mfi(*q[l]); mfi.isValid(); ++mfi) {
            auto const& o = q[l]->array(mfi);
            auto const& a = f[l]->const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { o(i,j,k) = a(i,j,k)*a(i,j,k); });
        }
    }
    return q;
}

double total_volume (PressureWorkspace::Impl const& W)
{
    double v = 0.0;
    for (int l = 0; l < W.nlev; ++l) { v += double(W.vol[l]) * double(W.nunc[l]); }
    return v;
}

// Volume-weighted exact mean (and, if wanted, rms) over the uncovered cells of the hierarchy.
void hierarchy_mean (PressureWorkspace::Impl const& W, LevelMFs const& f, double& mean, double& rms, bool want_rms = true)
{
    const double V = total_volume(W);
    ExactSumResult s = exact_sum_multi(terms_of(W, f, true), 1);
    mean = (V > 0.0) ? s.sum[0] / V : 0.0;
    rms = 0.0;
    if (!want_rms) { return; }
    LevelMFs q = squares_of(W, f);
    ExactSumResult s2 = exact_sum_multi(terms_of(W, q, true), 1);
    rms = (V > 0.0) ? std::sqrt(s2.sum[0] / V) : 0.0;
}

void average_down_all (PressureWorkspace::Impl const& W, LevelMFs& f)
{
    for (int l = W.nlev - 1; l >= 1; --l) {
        amrex::average_down(*f[l], *f[l-1], 0, 1, W.ratio[l]);
    }
}

// Copy between the original (one cell in the hidden direction d) and the extruded layout (same box indices, same
// distribution mapping, so local boxes correspond one to one).
void replicate_to_ext (PressureWorkspace::Impl const& W, MultiFab const& orig, MultiFab& ext)
{
    const int d = W.hidden;
    const int lo = W.geom[0].Domain().smallEnd(d);
    for (MFIter mfi(ext); mfi.isValid(); ++mfi) {
        auto const& e = ext.array(mfi);
        auto const& o = orig.const_array(mfi.index());
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            int iv[3] = {i, j, k}; iv[d] = lo;
            e(i,j,k) = o(iv[0], iv[1], iv[2]);
        });
    }
}

void plane_from_ext (PressureWorkspace::Impl const& W, MultiFab const& ext, MultiFab& orig)
{
    (void)W;
    for (MFIter mfi(orig); mfi.isValid(); ++mfi) {
        auto const& o = orig.array(mfi);
        auto const& e = ext.const_array(mfi.index());
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { o(i,j,k) = e(i,j,k); });
    }
}

template <class V> Vector<MultiFab*> ptrs (V const& v)
{
    Vector<MultiFab*> r;
    for (auto const& m : v) { r.push_back(m.get()); }
    return r;
}

} // anonymous

// ---------------------------------------------------------------------------------------------------------------
// Masked composite solve (Q1: assembled HYPRE; D3/D4/D6)
// ---------------------------------------------------------------------------------------------------------------
namespace {

void average_down_gas (PressureWorkspace::Impl const& W, PressureProblem const& p, LevelMFs& f)
{
    for (int l = W.nlev - 1; l >= 1; --l) {
        // Gas-only average-down (D3e): covered coarse = mean of fine gas cells under it; else 0.
        MultiFab fine_gas(W.ba[l], W.dm[l], 1, 0);
        MultiFab fine_cnt(W.ba[l], W.dm[l], 1, 0);
        fine_gas.setVal(0.0); fine_cnt.setVal(0.0);
        iMultiFab const* cls = p.levels[l].cell_class;
        if (!cls) {
            MultiFab::Copy(fine_gas, *f[l], 0, 0, 1, 0);
            fine_cnt.setVal(1.0);
        } else {
            for (MFIter mfi(fine_gas); mfi.isValid(); ++mfi) {
                auto const& o = fine_gas.array(mfi);
                auto const& c = fine_cnt.array(mfi);
                auto const& a = f[l]->const_array(mfi);
                auto const& cl = cls->const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    if (cl(i,j,k) == CellGas) { o(i,j,k) = a(i,j,k); c(i,j,k) = Real(1); }
                });
            }
        }
        MultiFab crse_sum(W.ba[l-1], W.dm[l-1], 1, 0);
        MultiFab crse_cnt(W.ba[l-1], W.dm[l-1], 1, 0);
        crse_sum.setVal(0.0); crse_cnt.setVal(0.0);
        amrex::average_down(fine_gas, crse_sum, 0, 1, W.ratio[l]);
        amrex::average_down(fine_cnt, crse_cnt, 0, 1, W.ratio[l]);
        for (MFIter mfi(*f[l-1]); mfi.isValid(); ++mfi) {
            auto const& dest = f[l-1]->array(mfi);
            auto const& s = crse_sum.const_array(mfi);
            auto const& c = crse_cnt.const_array(mfi);
            auto const& u = W.unc[l-1]->const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                if (u(i,j,k) != 0) { return; }
                if (c(i,j,k) > Real(0)) { dest(i,j,k) = s(i,j,k) / c(i,j,k); }
                else { dest(i,j,k) = Real(0); }
            });
        }
    }
}

PressureResult solve_composite_masked (PressureProblem const& p, PressureOptions const& o, PressureWorkspace* ws, bool full)
{
    MFIter::allowMultipleMFIters(true);
    PressureResult R;
    R.backend = "HYPRE";
#ifndef PB_WITH_HYPRE
    R.status = Status::NotBuilt; R.message = kHypreNotBuilt; R.backend = ""; return R;
#else
    PressureWorkspace local;
    if (!ws) { ws = &local; }
    if (!ws->matches(p)) {
        R.workspace_rebuilt = ws->num_levels() > 0;
        std::string msg;
        Status st = ws->rebuild(p, &msg);
        if (st != Status::Ok) { R.status = st; R.message = msg; return R; }
    }
    PressureWorkspace::Impl& W = *ws->impl();
    const int nlev = W.nlev;
    R.nlevels = nlev;

    // Label (D4); Q2 Known-on-C/F → InvalidInput.
    CompositeComponentMap cm;
    std::string lab_err = label_components_composite_masked(p, W.unc, W.geom, W.ratio, cm);
    if (!lab_err.empty()) {
        R.status = Status::InvalidInput; R.message = lab_err; R.backend = ""; return R;
    }
    R.mask_label_rounds = cm.rounds;
    R.mask_label_seconds = cm.seconds;
    if (cm.comps.empty()) {
        R.status = Status::InvalidInput; R.message = "masked composite without an uncovered gas cell"; R.backend = ""; return R;
    }

    Long nunc_gas = 0;
    for (auto const& c : cm.comps) { nunc_gas += c.ncells; }
    R.ncells_uncovered = nunc_gas;

    LevelMFs b, phi;
    const int ngp = 1;
    for (int l = 0; l < nlev; ++l) {
        b.push_back(std::make_unique<MultiFab>(W.ba[l], W.dm[l], 1, 0));
        MultiFab::Copy(*b[l], *p.levels[l].rhs, 0, 0, 1, 0);
        // Zero covered / Solid / Known RHS (not part of the gas problem).
        for (MFIter mfi(*b[l]); mfi.isValid(); ++mfi) {
            auto const& a = b[l]->array(mfi);
            auto const& lab = cm.label[l]->const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { if (lab(i,j,k) < 0) { a(i,j,k) = Real(0); } });
        }
        phi.push_back(std::make_unique<MultiFab>(W.ba[l], W.dm[l], 1, ngp));
        phi[l]->setVal(0.0);
        if (o.use_initial_guess) { MultiFab::Copy(*phi[l], *p.levels[l].phi, 0, 0, 1, 0); }
    }

    // Per-component mean removal (D4): level-volume weights, one mean per singular multi-level component.
    if (o.remove_mean) {
        const int nc = int(cm.comps.size());
        std::vector<double> sum_vb(nc, 0.0), sum_v(nc, 0.0), sum_vb2(nc, 0.0);
        for (int l = 0; l < nlev; ++l) {
            const double v = double(W.vol[l]);
            for (MFIter mfi(*b[l]); mfi.isValid(); ++mfi) {
                auto const& a = b[l]->const_array(mfi);
                auto const& lab = cm.label[l]->const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    const int id = lab(i,j,k);
                    if (id < 0) { return; }
                    const double x = double(a(i,j,k));
                    if (p.mean_kind == MeanKind::Volume) {
                        sum_vb[id] += v * x; sum_v[id] += v; sum_vb2[id] += v * x * x;
                    } else {
                        sum_vb[id] += v * x; sum_v[id] += 1.0; sum_vb2[id] += (v*x)*(v*x);
                    }
                });
            }
        }
        ParallelDescriptor::ReduceRealSum(sum_vb.data(), nc);
        ParallelDescriptor::ReduceRealSum(sum_v.data(), nc);
        if (full) { ParallelDescriptor::ReduceRealSum(sum_vb2.data(), nc); }
        std::vector<double> shift(nc, 0.0);
        for (int q = 0; q < nc; ++q) {
            if (!cm.comps[q].singular || sum_v[q] == 0.0) { continue; }
            double mean = 0.0, removed = 0.0, rms = 0.0;
            if (p.mean_kind == MeanKind::Volume) {
                mean = sum_vb[q] / sum_v[q];
                rms = (full && sum_v[q] > 0) ? std::sqrt(sum_vb2[q] / sum_v[q]) : 0.0;
                const double floor_ = std::ldexp(std::abs(mean) + rms, -52);
                if (std::abs(mean) > floor_) { removed = mean; shift[q] = mean; }
            } else {
                mean = sum_vb[q] / sum_v[q]; // arithmetic mean of F=v*b over cells (sum_v = ncells)
                rms = (full && sum_v[q] > 0) ? std::sqrt(sum_vb2[q] / sum_v[q]) : 0.0;
                const double floor_ = std::ldexp(std::abs(mean) + rms, -52);
                if (std::abs(mean) > floor_) { removed = mean; /* b -= mean/v per level below */ }
            }
            cm.comps[q].removed_mean = removed;
            cm.comps[q].removed_rel = (full && rms > 0.0) ? std::abs(mean) / rms : 0.0;
            if (p.mean_kind != MeanKind::Volume && std::abs(removed) > 0.0) {
                // store mean of F; shift applied as mean/v_l
                shift[q] = removed; // flag: ScaledArithmetic uses shift as mean(F)
            }
        }
        for (int l = 0; l < nlev; ++l) {
            const double v = double(W.vol[l]);
            for (MFIter mfi(*b[l]); mfi.isValid(); ++mfi) {
                auto const& a = b[l]->array(mfi);
                auto const& lab = cm.label[l]->const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    const int id = lab(i,j,k);
                    if (id < 0 || shift[id] == 0.0) { return; }
                    if (p.mean_kind == MeanKind::Volume) { a(i,j,k) -= shift[id]; }
                    else { a(i,j,k) -= shift[id] / v; }
                });
            }
        }
    }

    // HYPRE system with sealed pins.
    std::vector<HypreLayoutLevel> hl;
    for (int l = 0; l < nlev; ++l) {
        HypreLayoutLevel L;
        L.ba = W.ba[l]; L.dm = W.dm[l]; L.geom = W.geom[l]; L.ratio = W.ratio[l];
        L.cell_class = p.levels[l].cell_class;
        L.known_value = p.levels[l].known_value;
        hl.push_back(L);
    }
    std::array<BC,6> const ebc = effective_bc(p.bc, p.levels[0].geom.Domain());
    std::vector<std::pair<int, IntVect>> pins;
    for (auto const& cp : cm.pins) { pins.push_back({cp.level, cp.cell}); }
    const bool reused = W.hsys && W.hsys->ok() && W.hsys->masked() && W.hsys->matches_masked(hl, ebc, o.hypre, pins);
    if (!reused) { W.hsys = std::make_unique<HypreSystem>(hl, ebc, o.hypre, pins); }
    if (!W.hsys->ok()) {
        R.status = Status::NotBuilt; R.message = W.hsys->message(); R.backend = "";
        W.hsys.reset();
        return R;
    }
    std::vector<MultiFab*> ph; std::vector<MultiFab const*> bb;
    for (int l = 0; l < nlev; ++l) { ph.push_back(phi[l].get()); bb.push_back(b[l].get()); }
    BackendStatus bs = W.hsys->solve(ph, bb, o.tol_rel, o.max_iter, o.use_initial_guess);
    bs.plan_reused = reused;
    bs.setup_seconds = reused ? 0.0 : W.hsys->setup_seconds();
    R.backend_status = bs;

    // Fill Solid=0, Known=g on every level.
    for (int l = 0; l < nlev; ++l) {
        if (!p.levels[l].cell_class) { continue; }
        for (MFIter mfi(*phi[l]); mfi.isValid(); ++mfi) {
            auto const& f = phi[l]->array(mfi);
            auto const& cl = p.levels[l].cell_class->const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                if (cl(i,j,k) == CellSolid) { f(i,j,k) = Real(0); }
                else if (cl(i,j,k) == CellKnown) {
                    f(i,j,k) = p.levels[l].known_value ? p.levels[l].known_value->const_array(mfi)(i,j,k) : Real(0);
                }
            });
        }
    }

    // Per-component gauge (D4): volume-weighted mean of phi over the component's uncovered gas cells.
    bool any_sing = false;
    for (auto const& c : cm.comps) { any_sing = any_sing || c.singular; }
    LevelMFs phr;
    if (any_sing && o.check_residual && full) {
        for (int l = 0; l < nlev; ++l) {
            phr.push_back(std::make_unique<MultiFab>(W.ba[l], W.dm[l], 1, phi[l]->nGrow()));
            MultiFab::Copy(*phr[l], *phi[l], 0, 0, 1, 0);
        }
        // mean-remove per singular component on the copy
        const int nc = int(cm.comps.size());
        std::vector<double> sum_vp(nc, 0.0), sum_v(nc, 0.0);
        for (int l = 0; l < nlev; ++l) {
            const double v = double(W.vol[l]);
            for (MFIter mfi(*phr[l]); mfi.isValid(); ++mfi) {
                auto const& a = phr[l]->const_array(mfi);
                auto const& lab = cm.label[l]->const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    const int id = lab(i,j,k);
                    if (id < 0 || !cm.comps[id].singular) { return; }
                    sum_vp[id] += v * double(a(i,j,k)); sum_v[id] += v;
                });
            }
        }
        ParallelDescriptor::ReduceRealSum(sum_vp.data(), nc);
        ParallelDescriptor::ReduceRealSum(sum_v.data(), nc);
        for (int l = 0; l < nlev; ++l) {
            for (MFIter mfi(*phr[l]); mfi.isValid(); ++mfi) {
                auto const& a = phr[l]->array(mfi);
                auto const& lab = cm.label[l]->const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    const int id = lab(i,j,k);
                    if (id < 0 || !cm.comps[id].singular || sum_v[id] == 0.0) { return; }
                    a(i,j,k) -= sum_vp[id] / sum_v[id];
                });
            }
        }
        average_down_gas(W, p, phr);
    }
    // Gauge on the returned solution.
    {
        const int nc = int(cm.comps.size());
        std::vector<double> sum_num(nc, 0.0), sum_den(nc, 0.0);
        const bool has_w = (p.levels[0].gauge_weight != nullptr);
        const bool has_g = (p.levels[0].gauge_offset != nullptr);
        for (int l = 0; l < nlev; ++l) {
            const double v = double(W.vol[l]);
            for (MFIter mfi(*phi[l]); mfi.isValid(); ++mfi) {
                auto const& a = phi[l]->const_array(mfi);
                auto const& lab = cm.label[l]->const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    const int id = lab(i,j,k);
                    if (id < 0 || !cm.comps[id].singular) { return; }
                    double rho = 1.0, g0 = 0.0;
                    if (has_w) { rho = double(p.levels[l].gauge_weight->const_array(mfi)(i,j,k)); }
                    if (has_g) { g0 = double(p.levels[l].gauge_offset->const_array(mfi)(i,j,k)); }
                    sum_num[id] += v * rho * (double(a(i,j,k)) - g0);
                    sum_den[id] += v * rho;
                });
            }
        }
        ParallelDescriptor::ReduceRealSum(sum_num.data(), nc);
        ParallelDescriptor::ReduceRealSum(sum_den.data(), nc);
        for (int q = 0; q < nc; ++q) {
            if (!cm.comps[q].singular || sum_den[q] == 0.0) { continue; }
            cm.comps[q].gauge_shift = sum_num[q] / sum_den[q];
        }
        for (int l = 0; l < nlev; ++l) {
            for (MFIter mfi(*phi[l]); mfi.isValid(); ++mfi) {
                auto const& a = phi[l]->array(mfi);
                auto const& lab = cm.label[l]->const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    const int id = lab(i,j,k);
                    if (id < 0 || !cm.comps[id].singular) { return; }
                    a(i,j,k) -= cm.comps[id].gauge_shift;
                });
            }
        }
    }
    average_down_gas(W, p, phi);

    // Residual via HYPRE apply: A = -V L ⇒ L phi = -A phi / V; res = b - L phi = b + Aphi/V (Known already in A/RHS).
    if (o.check_residual && full) {
        LevelMFs& phv = phr.empty() ? phi : phr;
        LevelMFs Ax;
        for (int l = 0; l < nlev; ++l) {
            Ax.push_back(std::make_unique<MultiFab>(W.ba[l], W.dm[l], 1, 0));
            Ax[l]->setVal(0.0);
        }
        std::vector<MultiFab const*> px; std::vector<MultiFab*> py;
        for (int l = 0; l < nlev; ++l) { px.push_back(phv[l].get()); py.push_back(Ax[l].get()); }
        W.hsys->apply(px, py);
        LevelMFs res;
        for (int l = 0; l < nlev; ++l) {
            res.push_back(std::make_unique<MultiFab>(W.ba[l], W.dm[l], 1, 0));
            const double v = double(W.vol[l]);
            for (MFIter mfi(*res[l]); mfi.isValid(); ++mfi) {
                auto const& r = res[l]->array(mfi);
                auto const& bb = b[l]->const_array(mfi);
                auto const& ax = Ax[l]->const_array(mfi);
                auto const& lab = cm.label[l]->const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    if (lab(i,j,k) < 0) { r(i,j,k) = Real(0); return; }
                    r(i,j,k) = Real(double(bb(i,j,k)) + double(ax(i,j,k)) / v);
                });
            }
        }
        ResidualSums S;
        S.singular = any_sing; S.what = "masked composite ";
        {
            LevelMFs rq = squares_of(W, res), bq = squares_of(W, b), pq = squares_of(W, phv);
            // Restrict sums to gas labels: zero non-gas in copies used for sums? terms_of uses unc; also zero Solid via label.
            // Zero residual on pin cells for nopin path below; for r2 use gas uncovered via unc (Solid covered by label zeroing already in res).
            S.r2 = exact_sum_multi(terms_of(W, rq, true), 1).sum[0];
            S.b2 = exact_sum_multi(terms_of(W, bq, true), 1).sum[0];
            S.phi2 = exact_sum_multi(terms_of(W, pq, true), 1).sum[0];
            S.sumr = exact_sum_multi(terms_of(W, res, true), 1).sum[0];
            S.sumw = 0.0;
            for (auto const& c : cm.comps) { S.sumw += double(c.ncells); } // cell count; volumes differ — use volume sum:
            S.sumw = 0.0;
            for (int l = 0; l < nlev; ++l) {
                Long ng = 0;
                for (MFIter mfi(*cm.label[l]); mfi.isValid(); ++mfi) {
                    auto const& lab = cm.label[l]->const_array(mfi);
                    amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { if (lab(i,j,k) >= 0) { ++ng; } });
                }
                ParallelDescriptor::ReduceLongSum(ng);
                S.sumw += double(W.vol[l]) * double(ng);
            }
            S.rmax = max_abs_uncovered(W, res); S.bmax = max_abs_uncovered(W, b);
            std::vector<SumTerm> ta; ta.reserve(nlev);
            double A_finest = 0.0;
            for (int l = 0; l < nlev; ++l) {
                const double al = operator_anorm(W.geom[l].CellSizeArray(), W.geom[l].Domain());
                A_finest = al;
                SumTerm t; t.mf = pq[l].get(); t.comp = 0; t.weight = double(W.vol[l])*al*al; t.uncovered = W.unc[l].get();
                ta.push_back(t);
            }
            const double a2phi2 = exact_sum_multi(ta, 1).sum[0];
            const double anorm_eff = (S.phi2 > 0.0) ? std::sqrt(a2phi2 / S.phi2) : 0.0;
            S.anorm = std::min(2.0 * anorm_eff, A_finest);
        }
        S.r2_nopin = S.r2; S.rmax_nopin = S.rmax;
        if (bs.pin_applied && !bs.pins.empty()) {
            double pin2 = 0.0, pin_abs = 0.0;
            for (auto const& pr : bs.pins) {
                const int pl = pr.first; const IntVect pc = pr.second;
                double rp = 0.0;
                for (MFIter mfi(*res[pl]); mfi.isValid(); ++mfi) {
                    if (mfi.validbox().contains(pc)) {
                        auto const& a = res[pl]->array(mfi);
                        rp += double(a(pc[0], pc[1], pc[2]));
                        a(pc[0], pc[1], pc[2]) = Real(0);
                    }
                }
                ParallelDescriptor::ReduceRealSum(rp);
                pin2 += double(W.vol[pl]) * rp * rp;
                pin_abs = std::max(pin_abs, std::abs(rp));
            }
            S.pin2 = pin2; S.pin_abs = pin_abs;
            LevelMFs rq2 = squares_of(W, res);
            S.r2_nopin = exact_sum_multi(terms_of(W, rq2, true), 1).sum[0];
            S.rmax_nopin = max_abs_uncovered(W, res);
            S.pin_excluded = true;
        }
        evaluate_residual(R, o, S);
    }

    if (full) {
        for (auto const& c : cm.comps) {
            if (c.singular && c.removed_rel > o.removed_mean_warn) {
                std::ostringstream m;
                m << "removed mean of component " << c.id << " is " << c.removed_rel << " of rms(b), above " << o.removed_mean_warn;
                R.warnings.push_back(m.str());
            }
        }
    }
    if (o.verbose > 0) { for (auto const& w : R.warnings) { Print() << "PRESSURE WARNING: " << w << "\n"; } }
    R.components = cm.comps;
    if (!bs.converged) {
        R.status = Status::NotConverged; R.message = "backend did not converge";
        if (o.verbose > 0) { Print() << "PRESSURE WARNING: backend did not converge\n"; }
    }
    for (int l = 0; l < nlev; ++l) { MultiFab::Copy(*p.levels[l].phi, *phi[l], 0, 0, 1, 0); }
    return R;
#endif
}

} // namespace

// ---------------------------------------------------------------------------------------------------------------
// Solve
// ---------------------------------------------------------------------------------------------------------------
PressureResult solve_composite (PressureProblem const& p, PressureOptions const& o, PressureWorkspace* ws, bool full)
{
    PressureResult R;
    const bool masked = composite_is_masked(p);
    Selection const sel0 = select_composite(p, o.backend);
    const bool use_hypre = (sel0.kind == BackendKind::HYPRE) || (o.backend == BackendKind::HYPRE);
    R.backend = use_hypre ? "HYPRE" : "MLMG";
    if (masked && !use_hypre) {
        R.status = Status::NotBuilt;
        R.message = "composite with masked cells is not built on MLMG (C/F interpolation cannot exclude Solid/Known per stencil; Q1: use assembled HYPRE)";
        R.backend = "";
        return R;
    }
    PressureWorkspace local;
    if (!ws) { ws = &local; }
    double setup_s = 0.0;   // workspace (MLMG operator) build; 0 when it was reused
    if (!ws->matches(p)) {
        const double t0 = amrex::second();
        R.workspace_rebuilt = ws->num_levels() > 0;
        std::string msg;
        Status st = ws->rebuild(p, &msg);
        if (st != Status::Ok) { R.status = st; R.message = msg; return R; }
        setup_s = amrex::second() - t0;
    }
    PressureWorkspace::Impl& W = *ws->impl();
    const int nlev = W.nlev;
    R.nlevels = nlev;

    if (masked) {
        return solve_composite_masked(p, o, ws, full);
    }

    bool open = false;
    {
        std::array<BC,6> const eb = effective_bc(p.bc, p.levels[0].geom.Domain());
        for (int f = 0; f < 6; ++f) { open = open || (eb[f] == BC::Dirichlet); }
    }
    ComponentInfo ci;
    ci.id = 0; ci.singular = !open;
    for (int l = 0; l < nlev; ++l) { ci.ncells += W.nunc[l]; }
    R.ncells_uncovered = ci.ncells;
    if (full) {   // pin: lowest index (x fastest) uncovered cell of the coarsest level that has one; recorded, not applied
        for (int l = 0; l < nlev; ++l) {
            if (W.nunc[l] == 0) { continue; }
            Box const& dom = W.geom[l].Domain();
            const Long nx = dom.length(0), ny = dom.length(1);
            Long lowest = std::numeric_limits<Long>::max();
            for (MFIter mfi(*W.unc[l]); mfi.isValid(); ++mfi) {
                auto const& u = W.unc[l]->const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    if (u(i,j,k) != 0) {
                        lowest = std::min(lowest, Long(i - dom.smallEnd(0)) + nx*(Long(j - dom.smallEnd(1)) + ny*Long(k - dom.smallEnd(2))));
                    }
                });
            }
            ParallelDescriptor::ReduceLongMin(lowest);
            ci.pin = IntVect(int(lowest % nx) + dom.smallEnd(0), int((lowest / nx) % ny) + dom.smallEnd(1), int(lowest / (nx*ny)) + dom.smallEnd(2));
            break;
        }
    }

    // Right-hand side copy (covered cells set to zero: they are replaced by the fine residual inside MLMG), phi work fields.
    LevelMFs b, phi;
    const int ngp = 1;
    for (int l = 0; l < nlev; ++l) {
        b.push_back(std::make_unique<MultiFab>(W.ba[l], W.dm[l], 1, 0));
        MultiFab::Copy(*b[l], *p.levels[l].rhs, 0, 0, 1, 0);
        for (MFIter mfi(*b[l]); mfi.isValid(); ++mfi) {
            auto const& a = b[l]->array(mfi);
            auto const& u = W.unc[l]->const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { if (u(i,j,k) == 0) { a(i,j,k) = Real(0.0); } });
        }
        phi.push_back(std::make_unique<MultiFab>(W.ba[l], W.dm[l], 1, ngp));
        phi[l]->setVal(0.0);
        if (o.use_initial_guess) { MultiFab::Copy(*phi[l], *p.levels[l].phi, 0, 0, 1, 0); }
    }

    // Mean removal (D-067): composite compatibility is sum over uncovered cells of v*b = 0 (v = cell volume of the level).
    //  Volume (default): the same constant sum(v*b)/sum(v) is subtracted on every level.
    //  ScaledArithmetic (FDS parity): the arithmetic mean of F = v*b is removed from F, b_k -= mean(F)/v_l.
    if (o.remove_mean) {
        std::vector<double> shift_lev(nlev, 0.0);     // constant subtracted from b on each level
        double mean = 0.0, rms = 0.0, removed = 0.0, floor_ = 0.0;
        if (p.mean_kind == MeanKind::Volume) {
            hierarchy_mean(W, b, mean, rms, full);
            floor_ = std::ldexp(max_abs_uncovered(W, b), -52);   // idempotence: below round-off of b itself, leave alone
            removed = mean;
            if (std::abs(mean) > floor_ && ci.singular) { for (int l = 0; l < nlev; ++l) { shift_lev[l] = mean; } }
            else { removed = 0.0; }
        } else {
            LevelMFs F;
            for (int l = 0; l < nlev; ++l) {
                F.push_back(std::make_unique<MultiFab>(W.ba[l], W.dm[l], 1, 0));
                MultiFab::Copy(*F[l], *b[l], 0, 0, 1, 0);
                F[l]->mult(W.vol[l], 0, 1, 0);
            }
            ExactSumResult sF = exact_sum_multi(terms_of(W, F, false), 1);
            const double n = double(sF.count[0]);
            mean = (n > 0) ? sF.sum[0] / n : 0.0;                 // arithmetic mean of F
            if (full) {
                LevelMFs Fq = squares_of(W, F);
                ExactSumResult sF2 = exact_sum_multi(terms_of(W, Fq, false), 1);
                rms = (n > 0) ? std::sqrt(sF2.sum[0] / n) : 0.0;  // rms of F
            }
            floor_ = std::ldexp(max_abs_uncovered(W, F), -52);
            if (std::abs(mean) > floor_ && ci.singular) {
                removed = mean;
                for (int l = 0; l < nlev; ++l) { shift_lev[l] = mean / W.vol[l]; }
            }
        }
        if (ci.singular) {
            ci.removed_mean = removed;      // Volume: constant subtracted from b; ScaledArithmetic: mean of F = v*b
            ci.removed_rel = (full && rms > 0.0) ? std::abs(mean) / rms : 0.0;
        }
        for (int l = 0; l < nlev; ++l) {
            if (shift_lev[l] == 0.0) { continue; }
            const double shift = shift_lev[l];
            for (MFIter mfi(*b[l]); mfi.isValid(); ++mfi) {
                auto const& a = b[l]->array(mfi);
                auto const& u = W.unc[l]->const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { if (u(i,j,k) != 0) { a(i,j,k) -= shift; } });
            }
        }
    }

    // Backend.
    MLMG& mlmg = *W.mlmg;
    mlmg.setVerbose(o.verbose >= 3 ? o.verbose - 2 : 0);   // verbose >= 3 prints the MLMG iteration history
    mlmg.setMaxIter(o.max_iter);
    mlmg.setConvergenceNormType(MLMGNormType::bnorm);
    mlmg.setThrowException(true);
    // MLMG works on the original layout, or on the extruded copy when a direction has one cell.
    const bool ext = (W.hidden >= 0);
    LevelMFs phix, bx;
    if (ext) {
        for (int l = 0; l < nlev; ++l) {
            phix.push_back(std::make_unique<MultiFab>(W.eba[l], W.dm[l], 1, ngp));
            bx.push_back(std::make_unique<MultiFab>(W.eba[l], W.dm[l], 1, 0));
            phix[l]->setVal(0.0);
            replicate_to_ext(W, *b[l], *bx[l]);
            replicate_to_ext(W, *phi[l], *phix[l]);
        }
    }
    Vector<MultiFab*> pphi = ptrs(ext ? phix : phi);
    Vector<MultiFab const*> pb_;
    for (auto const& m : (ext ? bx : b)) { pb_.push_back(m.get()); }
    BackendStatus bs;
    if (!use_hypre) {
        bs = mlmg_solve(mlmg, pphi, pb_, o.mlmg_krylov, Real(o.tol_rel), o.max_iter);
        bs.setup_seconds = setup_s; bs.plan_reused = (setup_s == 0.0);
        R.backend_status = bs;
        if (ext) { for (int l = 0; l < nlev; ++l) { plane_from_ext(W, *phix[l], *phi[l]); } }
    }
#ifdef PB_WITH_HYPRE
    else {
        // HYPRE: assembled matrix on the original layout (a one-cell direction simply has no term), cached in the workspace.
        std::vector<HypreLayoutLevel> hl;
        for (int l = 0; l < nlev; ++l) { HypreLayoutLevel L; L.ba = W.ba[l]; L.dm = W.dm[l]; L.geom = W.geom[l]; L.ratio = W.ratio[l]; hl.push_back(L); }
        std::array<BC,6> const ebc = effective_bc(p.bc, p.levels[0].geom.Domain());
        const bool reused = W.hsys && W.hsys->ok() && W.hsys->matches(hl, ebc, o.hypre, ci.singular);
        if (!reused) { W.hsys = std::make_unique<HypreSystem>(hl, ebc, o.hypre, ci.singular); }
        if (!W.hsys->ok()) {
            R.status = Status::NotBuilt; R.message = W.hsys->message(); R.backend = "";
            W.hsys.reset();
            return R;
        }
        std::vector<MultiFab*> ph; std::vector<MultiFab const*> bb;
        for (int l = 0; l < nlev; ++l) { ph.push_back(phi[l].get()); bb.push_back(b[l].get()); }
        bs = W.hsys->solve(ph, bb, o.tol_rel, o.max_iter, o.use_initial_guess);
        bs.plan_reused = reused;
        bs.setup_seconds = reused ? 0.0 : W.hsys->setup_seconds();
        R.backend_status = bs;
    }
#endif

    // Gauge (D-067): sum(V*rho*(phi - KRES)) / sum(V*rho) over the uncovered cells of the hierarchy (exact sums, one
    // fixed-point scale per sum) is removed from phi on all levels; with rho = 1 and KRES = 0 this is the plain exact
    // volume-weighted mean. Then the covered coarse cells take the average-down of the fine solution. A component
    // that is not singular is not shifted.
    // E2: the residual verdict is evaluated on the solver's own solution with the composite volume-weighted mean removed (before this gauge),
    // so it and the round-off floor cannot depend on the gauge constant (see mean_removed_copy in CommonLayer.H).
    LevelMFs phr;
    if (ci.singular && o.check_residual && full) {
        for (int l = 0; l < nlev; ++l) {
            phr.push_back(std::make_unique<MultiFab>(W.ba[l], W.dm[l], 1, phi[l]->nGrow()));
            phr[l]->setVal(0.0);
            MultiFab::Copy(*phr[l], *phi[l], 0, 0, 1, 0);
        }
        double sh = 0.0, rm = 0.0;
        hierarchy_mean(W, phr, sh, rm, false);
        for (int l = 0; l < nlev; ++l) { phr[l]->plus(Real(-sh), 0, 1, 0); }
        average_down_all(W, phr);
    }
    LevelMFs& phv = phr.empty() ? phi : phr;     // the field the residual check sees
    if (ci.singular) {
        const bool has_w = (p.levels[0].gauge_weight != nullptr), has_g = (p.levels[0].gauge_offset != nullptr);
        double shift = 0.0;
        if (!has_w && !has_g) {
            double rms = 0.0;
            hierarchy_mean(W, phi, shift, rms, false);
        } else {
            LevelMFs X;                                               // phi - KRES
            std::vector<SumTerm> num, den;
            for (int l = 0; l < nlev; ++l) {
                X.push_back(std::make_unique<MultiFab>(W.ba[l], W.dm[l], 1, 0));
                MultiFab::Copy(*X[l], *phi[l], 0, 0, 1, 0);
                if (has_g) { MultiFab::Subtract(*X[l], *p.levels[l].gauge_offset, 0, 0, 1, 0); }
                SumTerm tn;
                tn.mf = X[l].get(); tn.weight = double(W.vol[l]); tn.uncovered = W.unc[l].get();
                tn.wfield = has_w ? p.levels[l].gauge_weight : nullptr;
                num.push_back(tn);
                if (has_w) {
                    SumTerm td;
                    td.mf = p.levels[l].gauge_weight; td.weight = double(W.vol[l]); td.uncovered = W.unc[l].get();
                    den.push_back(td);
                }
            }
            const double sn = exact_sum_multi(num, 1).sum[0];
            const double sd = has_w ? exact_sum_multi(den, 1).sum[0] : total_volume(W);
            shift = (sd > 0.0) ? sn / sd : 0.0;
        }
        ci.gauge_shift = shift;
        for (int l = 0; l < nlev; ++l) { phi[l]->plus(Real(-shift), 0, 1, 0); }
    }
    average_down_all(W, phi);

    // True residual: composite b - L phi on the uncovered cells, computed fresh from the final phi
    // (MLMG::compResidual: level residuals, reflux of the fine fluxes to the coarse cells, C/F ghost values from the
    // coarse solution), independent of the iteration's own residual.
    if (o.check_residual && full) {
        LevelMFs res, resx;
        for (int l = 0; l < nlev; ++l) {
            res.push_back(std::make_unique<MultiFab>(W.ba[l], W.dm[l], 1, 0)); res[l]->setVal(0.0);
            if (ext) {
                resx.push_back(std::make_unique<MultiFab>(W.eba[l], W.dm[l], 1, 0)); resx[l]->setVal(0.0);
                replicate_to_ext(W, *phv[l], *phix[l]);          // mean-removed phi (E2), covered cells = average-down
            }
        }
        Vector<MultiFab*> pres = ptrs(ext ? resx : res);
        Vector<MultiFab*> pphi_r = ptrs(ext ? phix : phv);
        mlmg.compResidual(pres, pphi_r, pb_);
        if (ext) { for (int l = 0; l < nlev; ++l) { plane_from_ext(W, *resx[l], *res[l]); } }
        ResidualSums S;
        S.singular = ci.singular; S.what = "composite ";
        {
            LevelMFs rq = squares_of(W, res), bq = squares_of(W, b), pq = squares_of(W, phv);
            S.r2 = exact_sum_multi(terms_of(W, rq, true), 1).sum[0];
            S.b2 = exact_sum_multi(terms_of(W, bq, true), 1).sum[0];
            S.phi2 = exact_sum_multi(terms_of(W, pq, true), 1).sum[0];
            S.sumr = exact_sum_multi(terms_of(W, res, true), 1).sum[0];
            S.sumw = total_volume(W);
            S.rmax = max_abs_uncovered(W, res); S.bmax = max_abs_uncovered(W, b);
            // E3 / D-087 (docs/pressure/08 §8.1): anorm_eff with level-wise A_l weighted like phi2; restore margin with
            // anorm = min(2*anorm_eff, A_finest). Same reduction as before (weight vol[l]*A_l^2); no new reduction.
            std::vector<SumTerm> ta; ta.reserve(nlev);
            double A_finest = 0.0;
            for (int l = 0; l < nlev; ++l) {
                const double al = operator_anorm(W.geom[l].CellSizeArray(), W.geom[l].Domain());
                A_finest = al;   // levels are coarse->fine; last is finest
                SumTerm t; t.mf = pq[l].get(); t.comp = 0; t.weight = double(W.vol[l])*al*al; t.uncovered = W.unc[l].get();
                ta.push_back(t);
            }
            const double a2phi2 = exact_sum_multi(ta, 1).sum[0];
            const double anorm_eff = (S.phi2 > 0.0) ? std::sqrt(a2phi2 / S.phi2) : 0.0;
            S.anorm = std::min(2.0 * anorm_eff, A_finest);
        }
        S.r2_nopin = S.r2; S.rmax_nopin = S.rmax;
        if (bs.pin_applied) {
            // Pin rows are identity, not equations of L: exclude every listed pin (D6); fall back to the single pin_cell.
            auto const& plist = bs.pins;
            const int npin = plist.empty() ? 1 : int(plist.size());
            double pin2 = 0.0, pin_abs = 0.0;
            for (int q = 0; q < npin; ++q) {
                const int pl = plist.empty() ? bs.pin_level : plist[q].first;
                const IntVect pc = plist.empty() ? bs.pin_cell : plist[q].second;
                double rp = 0.0;
                for (MFIter mfi(*res[pl]); mfi.isValid(); ++mfi) {
                    if (mfi.validbox().contains(pc)) {
                        auto const& a = res[pl]->array(mfi);
                        rp += double(a(pc[0], pc[1], pc[2]));
                        a(pc[0], pc[1], pc[2]) = Real(0);
                    }
                }
                ParallelDescriptor::ReduceRealSum(rp);
                pin2 += double(W.vol[pl]) * rp * rp;
                pin_abs = std::max(pin_abs, std::abs(rp));
            }
            S.pin2 = pin2; S.pin_abs = pin_abs;
            LevelMFs rq2 = squares_of(W, res);
            S.r2_nopin = exact_sum_multi(terms_of(W, rq2, true), 1).sum[0];
            S.rmax_nopin = max_abs_uncovered(W, res);
            S.pin_excluded = true;
        }
        evaluate_residual(R, o, S);
        if (o.verbose >= 2 && bs.pin_applied) {
            Print() << "PRESSURE INFO: composite residual full " << R.residual_rel2 << " (max " << R.residual_relmax << "), pin cell excluded " << R.residual_rel2_nopin
                    << " (max " << R.residual_relmax_nopin << "), pin row " << R.residual_pin_rel << "\n";
        }
    }
    if (full && ci.singular && ci.removed_rel > o.removed_mean_warn) {
        std::ostringstream m;
        m << "removed mean of component " << ci.id << " is " << ci.removed_rel << " of rms(b), above " << o.removed_mean_warn;
        R.warnings.push_back(m.str());
    }
    if (o.verbose > 0) { for (auto const& w : R.warnings) { Print() << "PRESSURE WARNING: " << w << "\n"; } }
    R.components.push_back(ci);
    if (!bs.converged) {
        R.status = Status::NotConverged; R.message = "backend did not converge";
        if (o.verbose > 0) { Print() << "PRESSURE WARNING: backend did not converge\n"; }
    }
    for (int l = 0; l < nlev; ++l) { MultiFab::Copy(*p.levels[l].phi, *phi[l], 0, 0, 1, 0); }
    return R;
}

// ---------------------------------------------------------------------------------------------------------------
// Face-centred gradient
// ---------------------------------------------------------------------------------------------------------------
PressureResult face_gradient_composite (PressureProblem const& p, std::vector<std::array<MultiFab*,3>>& grad_out)
{
    PressureResult R;
    R.backend = "MLMG";
    auto fail = [&] (Status st, std::string const& msg) { R.status = st; R.message = msg; return R; };
    if (p.levels.empty()) { return fail(Status::InvalidInput, "face_gradient_composite needs a problem with levels"); }
    Selection sel = select_composite(p, BackendKind::Auto);
    if (!sel.ok) { return fail(Status::NotBuilt, sel.message); }
    std::string bad = validate_composite(p);
    if (!bad.empty()) { return fail(Status::InvalidInput, bad); }
    const int nlev = static_cast<int>(p.levels.size());
    if (static_cast<int>(grad_out.size()) != nlev) { return fail(Status::InvalidInput, "grad_out needs one entry per level"); }
    for (int l = 0; l < nlev; ++l) {
        for (int d = 0; d < 3; ++d) {
            MultiFab* g = grad_out[l][d];
            if (!g) { return fail(Status::InvalidInput, "grad_out has a null MultiFab"); }
            if (g->nComp() < 1) { return fail(Status::InvalidInput, "grad_out needs one component"); }
            if (!(g->boxArray() == amrex::convert(p.levels[l].ba, IntVect::TheDimensionVector(d))) || !(g->DistributionMap() == p.levels[l].dm)) {
                return fail(Status::InvalidInput, "grad_out BoxArray/DistributionMapping differ from the face-centred level layout");
            }
        }
    }
    if (composite_is_masked(p)) {
        return fail(Status::NotBuilt, "face_gradient_composite with masked cells is not built");
    }
    PressureWorkspace::Impl W;
    build_impl(W, p);
    const bool ext = (W.hidden >= 0);
    LevelMFs phi;
    Vector<MultiFab*> pphi;
    std::vector<std::array<std::unique_ptr<MultiFab>,3>> flux(nlev);
    Vector<Array<MultiFab*,AMREX_SPACEDIM>> pflux(nlev);
    for (int l = 0; l < nlev; ++l) {
        phi.push_back(std::make_unique<MultiFab>(W.eba[l], W.dm[l], 1, 1));
        phi[l]->setVal(0.0);
        if (ext) {
            MultiFab tmp(W.ba[l], W.dm[l], 1, 0);
            MultiFab::Copy(tmp, *p.levels[l].phi, 0, 0, 1, 0);
            replicate_to_ext(W, tmp, *phi[l]);
        } else {
            MultiFab::Copy(*phi[l], *p.levels[l].phi, 0, 0, 1, 0);
        }
        pphi.push_back(phi[l].get());
        for (int d = 0; d < 3; ++d) {
            flux[l][d] = std::make_unique<MultiFab>(amrex::convert(W.eba[l], IntVect::TheDimensionVector(d)), W.dm[l], 1, 0);
            flux[l][d]->setVal(0.0);
            pflux[l][d] = flux[l][d].get();
        }
    }
    W.mlp->prepareForSolve();
    for (int l = nlev - 1; l >= 0; --l) { W.mlp->prepareForFluxes(l, l > 0 ? pphi[l-1] : nullptr); }
    W.mlp->getFluxes(pflux, pphi, MLLinOp::Location::FaceCenter);
    // flux = -dphi/dn; copy with the sign of the gradient. The one-cell direction has no gradient.
    for (int l = 0; l < nlev; ++l) {
        for (int d = 0; d < 3; ++d) {
            MultiFab& g = *grad_out[l][d];
            if (d == W.hidden) { g.setVal(0.0, 0, 1, 0); continue; }
            if (ext) {
                for (MFIter mfi(g); mfi.isValid(); ++mfi) {
                    auto const& o = g.array(mfi);
                    auto const& e = flux[l][d]->const_array(mfi.index());
                    amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { o(i,j,k) = -e(i,j,k); });
                }
            } else {
                MultiFab::Copy(g, *flux[l][d], 0, 0, 1, 0);
                g.mult(Real(-1.0), 0, 1, 0);
            }
        }
    }
    // D-032 (3): coarse faces covered by fine faces take the average of the fine gradients.
    for (int l = nlev - 1; l >= 1; --l) {
        Array<MultiFab const*,3> fine{grad_out[l][0], grad_out[l][1], grad_out[l][2]};
        Array<MultiFab*,3> crse{grad_out[l-1][0], grad_out[l-1][1], grad_out[l-1][2]};
        amrex::average_down_faces(fine, crse, W.ratio[l], 0);
        // A face on a periodic domain boundary exists twice (index lo and index hi+1). If one copy is under the
        // fine level and the other is not, the fine average must be carried over to the other copy.
        for (int d = 0; d < 3; ++d) {
            if (d == W.hidden || !W.geom[l-1].isPeriodic(d)) { continue; }
            BoxArray fc = W.ba[l];
            fc.coarsen(W.ratio[l]);                                       // covered coarse cells
            MultiFab& g = *grad_out[l-1][d];
            MultiFab cov(g.boxArray(), g.DistributionMap(), 1, 0);
            cov.setVal(0.0);
            for (MFIter mfi(cov); mfi.isValid(); ++mfi) {
                auto const& c = cov.array(mfi);
                for (int i = 0; i < static_cast<int>(fc.size()); ++i) {
                    const Box fb = amrex::surroundingNodes(fc[i], d) & mfi.validbox();
                    if (fb.ok()) { amrex::LoopOnCpu(fb, [&] (int ii, int jj, int kk) { c(ii,jj,kk) = 1.0; }); }
                }
            }
            MultiFab tv(g.boxArray(), g.DistributionMap(), 1, 0), rv(g.boxArray(), g.DistributionMap(), 1, 0), rc(g.boxArray(), g.DistributionMap(), 1, 0);
            MultiFab::Copy(tv, g, 0, 0, 1, 0);
            MultiFab::Multiply(tv, cov, 0, 0, 1, 0);
            rv.setVal(0.0); rc.setVal(0.0);
            rv.ParallelAdd(tv, 0, 0, 1, 0, 0, W.geom[l-1].periodicity());
            rc.ParallelAdd(cov, 0, 0, 1, 0, 0, W.geom[l-1].periodicity());
            for (MFIter mfi(g); mfi.isValid(); ++mfi) {
                auto const& o = g.array(mfi); auto const& v = rv.const_array(mfi); auto const& n = rc.const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { if (n(i,j,k) > 0.0) { o(i,j,k) = v(i,j,k)/n(i,j,k); } });
            }
        }
    }
    return R;
}

} // namespace pb
