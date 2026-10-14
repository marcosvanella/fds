// Composite (multi-level) pressure solve: amrex::MLMG over MLPoisson on the whole hierarchy (milestone M2, first
// version). The common layer (mean removal, gauge, true residual, diagnostics) runs from outside on the
// uncovered cells of all levels with one fixed-point scale (ExactSum.H). See frozen/composite-notes.md.
#include "CommonLayer.H"
#include "Composite.H"
#include "ExactSum.H"
#include "PbWorkspaceImpl.H"

#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_FFT_Poisson.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#if defined(PB_WITH_HYPRE) && defined(AMREX_USE_HYPRE)
#include "HYPRE_utilities.h"
#endif

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

// MLPoisson with an FFT bottom (PressureOptions::mlmg_fft_bottom): BottomSolver::custom solves the coarsest MG level of AMR level 0 exactly with
// amrex::FFT::Poisson (same 7-point operator and homogeneous BCs). The plan is built on the first bottom solve.
class FftBottomPoisson final : public MLPoisson {
public:
    using MLPoisson::MLPoisson;
    [[nodiscard]] bool supportCustomBottomSolver () const override { return true; }
    // Bottom level: at most one box per rank, and periodic / Neumann / Dirichlet faces only.
    bool fft_bottom_ok () const
    {
        const int mg = NMGLevels(0) - 1;
        if (static_cast<int>(m_grids[0][mg].size()) > ParallelDescriptor::NProcs()) { return false; }
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            for (auto t : {m_lobc[0][d], m_hibc[0][d]}) {
                if (t != LinOpBCType::Periodic && t != LinOpBCType::Neumann && t != LinOpBCType::Dirichlet) { return false; }
            }
        }
        return true;
    }
    void customBottomSolve (MLMG* mlmg, MultiFab& x, MultiFab const& b, Real, Real, int) override
    {
        if (!m_fft) {
            Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM> bc;
            auto fb = [] (LinOpBCType t) { return t == LinOpBCType::Periodic ? FFT::Boundary::periodic : t == LinOpBCType::Neumann ? FFT::Boundary::even : FFT::Boundary::odd; };
            for (int d = 0; d < AMREX_SPACEDIM; ++d) { bc[d] = {fb(m_lobc[0][d]), fb(m_hibc[0][d])}; }
            m_fft = std::make_unique<FFT::Poisson<MultiFab>>(Geom(0, NMGLevels(0) - 1), bc);
        }
        m_fft->solve(x, b);
        mlmg->postCG(0, 1);
    }
private:
    std::unique_ptr<FFT::Poisson<MultiFab>> m_fft;
};

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
    // Masked composite: no MLPoisson hierarchy (HYPRE assembles its own matrix; masked MLMG builds MaskedMlmg on its first solve).
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
    W.face_neumann = p.face_neumann;
    const bool robin = has_face_neumann(p);
    if (robin) {
        // Mixed domain faces: MLABecLaplacian with A = 0, B = -1 (the MLPoisson operator) and Robin a phi + b dphi/dn = 0 on the flagged
        // faces: a = 1, b = 0 (Dirichlet) except b = 1, a = 0 (Neumann) over the face_neumann cells of levels >= 1.
        auto op = std::make_unique<MLABecLaplacian>(geoms, bas, dms, info);
        op->setScalars(0.0, -1.0);
        for (int l = 0; l < n; ++l) { op->setACoeffs(l, 0.0); op->setBCoeffs(l, 1.0); }
        W.mlp = std::move(op);
    } else {
        W.mlp = std::make_unique<FftBottomPoisson>(geoms, bas, dms, info);
    }
    W.mlp->setMaxOrder(2);
    Array<LinOpBCType,AMREX_SPACEDIM> lo{lin_bc(ebc[face_index(0,0)]), lin_bc(ebc[face_index(1,0)]), lin_bc(ebc[face_index(2,0)])};
    Array<LinOpBCType,AMREX_SPACEDIM> hi{lin_bc(ebc[face_index(0,1)]), lin_bc(ebc[face_index(1,1)]), lin_bc(ebc[face_index(2,1)])};
    for (int d = 0; d < 3 && robin; ++d) {
        if (!p.face_neumann[face_index(d,0)].empty()) { lo[d] = LinOpBCType::Robin; }
        if (!p.face_neumann[face_index(d,1)].empty()) { hi[d] = LinOpBCType::Robin; }
    }
    W.mlp->setDomainBC(lo, hi);
    // Homogeneous Dirichlet values (nullptr) on every level. The coarse-fine boundary of level l > 0 is
    // filled by MLMG from level l-1 at every solve.
    for (int l = 0; l < n; ++l) {
        if (!robin) { W.mlp->setLevelBC(l, nullptr); continue; }
        MultiFab ra(W.eba[l], W.dm[l], 1, 1), rb(W.eba[l], W.dm[l], 1, 1), rf(W.eba[l], W.dm[l], 1, 1);
        ra.setVal(1.0); rb.setVal(0.0); rf.setVal(0.0);
        const Box dom = W.egeom[l].Domain();
        for (int f = 0; f < 6 && l > 0; ++f) {
            if (p.face_neumann[f].empty()) { continue; }
            const int d = f % 3, side = f / 3, e = side ? dom.bigEnd(d) + 1 : dom.smallEnd(d) - 1;
            for (MFIter mfi(ra); mfi.isValid(); ++mfi) {
                Box g = mfi.growntilebox(1);
                if (e < g.smallEnd(d) || e > g.bigEnd(d)) { continue; }
                g.setRange(d, e);
                auto const& A = ra.array(mfi); auto const& B = rb.array(mfi);
                amrex::LoopOnCpu(g, [&] (int i, int j, int k) {
                    IntVect iv(i,j,k); iv[d] = side ? dom.bigEnd(d) : dom.smallEnd(d);
                    if (face_neumann_at(p, l, f, iv)) { A(i,j,k) = 0.0; B(i,j,k) = 1.0; }
                });
            }
        }
        W.mlp->setLevelBC(l, nullptr, &ra, &rb, &rf);
    }
    W.mlmg = std::make_unique<MLMG>(*W.mlp);
}

// MLABecLaplacian folds the Robin terms into its A coefficients at every prepareForSolve and needs setScalars before each one.
void refresh_robin (PressureWorkspace::Impl& W)
{
    if (auto* op = dynamic_cast<MLABecLaplacian*>(W.mlp.get())) { op->setScalars(0.0, -1.0); }
}

bool impl_matches (PressureWorkspace::Impl const& W, PressureProblem const& p)
{
    const int n = static_cast<int>(p.levels.size());
    if (W.nlev != n || W.bc != p.bc || W.face_neumann != p.face_neumann) { return false; }
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
    if (has_face_neumann(p) && (masked || requested == BackendKind::HYPRE)) {
        s.message = "composite with a mixed domain face (face_neumann) is built on unmasked MLMG only (Robin BC); HYPRE and masked composite are not built"; return s;
    }
#if !defined(PB_WITH_HYPRE) || !defined(AMREX_USE_HYPRE)
    if (masked && requested == BackendKind::MLMG) {
        s.message = "composite with masked cells on MLMG needs the HYPRE bottom solver (pressure/13 ruling 6): built without HYPRE"; return s;
    }
#endif
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
    if (hd >= 0 && has_face_neumann(p)) { s.message = "composite with a mixed domain face (face_neumann) and a one-cell direction is not built"; return s; }
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
    for (int f = 0; f < 6 && n > 0; ++f) {
        if (p.face_neumann[f].empty()) { continue; }
        const Box dom = p.levels[0].geom.Domain();
        const int d = f % 3, t1 = (d == 0) ? 1 : 0, t2 = (d == 2) ? 1 : 2;
        if (p.bc[f] != BC::Dirichlet) { return "face_neumann is set on a face that is not Dirichlet"; }
        if (static_cast<Long>(p.face_neumann[f].size()) != static_cast<Long>(dom.length(t1)) * dom.length(t2)) { return "face_neumann size differs from the level-0 face cell count"; }
    }
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
// Masked composite MLMG (P3, docs/pressure/13): MLABecLaplacian, A = 1, B = 1 (a phi - div(beta grad phi), the negative of L),
// beta = 0 on every face touching a Solid or Known cell, 1 elsewhere; a = 1 on Solid/Known cells (identity rows), a = sum 2/dx^2 over
// the Known faces of a gas cell (as MLMGBackend.cpp, solve_masked). Rulings: (1) every coarse face beta is the arithmetic average
// of the fine faces (AMReX average_down_faces, both between AMR levels and between MG levels); (2) level-0 coefficients under level 1
// are the averaged level-1 ones (AMReX averageDownCoeffsToCoarseAmrLevel overwrites the covered cells and faces, the level-0 mask
// is not used there); (3) mechanism: a = 1 only on cells whose whole volume is excluded at that level; after AMReX's averaging, a is
// re-set from the averaged excluded fraction f (f == 1: a = 1, else a = the averaged Known term, 0 for a Solid-only problem), so a
// partly solid cell keeps a pure Laplacian row (the AMReX overset mask was not used: its coarsening stops at the first partly masked
// cell); (4) max_coarsening_level from the thinnest fluid-separating Solid run; (5) a Solid run shorter than the ratio between gas
// cells on a level >= 1 is NotBuilt; (6) bottom solver HYPRE BoomerAMG (IJ interface, options namespace pb_mlmg_bottom, HYPRE row 0
// pinned when the bottom is singular), gauge: the mean is removed per singular fluid region of each level-0 MG level (MaskedABecLap).
// Smoother: AMReX GSRB, V(2,2). C/F: MLMG's own interpolation reads the coarse cell under an open C/F face and its tangential
// neighbours, so a Solid/Known coarse cell there (HYPRE drops it from the stencil) is NotBuilt: the obstacle may not cross the C/F
// interface. With that, MLMG and the assembled HYPRE operator are the same discretisation.
// ---------------------------------------------------------------------------------------------------------------
#if defined(PB_WITH_HYPRE) && defined(AMREX_USE_HYPRE)
namespace {
constexpr int kRunCap = 4;                // ruling 4: Solid runs are measured up to this many level-0 cells (ghost width)

class MaskedABecLap final : public MLABecLaplacian {
public:
    using MLABecLaplacian::MLABecLaplacian;
    Vector<Vector<MultiFab>> frac, akn;   // per AMR and MG level: excluded volume fraction, averaged Known a-term
    void set_masks (int amrlev, MultiFab const& f, MultiFab const& k)
    {
        if (frac.empty()) { frac.resize(NAMRLevels()); akn.resize(NAMRLevels()); }
        for (auto* v : {&frac, &akn}) {
            (*v)[amrlev].clear();
            for (int mg = 0; mg < NMGLevels(amrlev); ++mg) {
                MultiFab const& a = m_a_coeffs[amrlev][mg];
                (*v)[amrlev].emplace_back(a.boxArray(), a.DistributionMap(), 1, 0);
                (*v)[amrlev].back().setVal(0.0);
            }
        }
        MultiFab::Copy(frac[amrlev][0], f, 0, 0, 1, 0);
        MultiFab::Copy(akn[amrlev][0], k, 0, 0, 1, 0);
        m_chained = false; m_regions = false;
    }
    // Ruling 6 gauge (D-067 per connected fluid region): a region of an MG level of AMR level 0 (cells not fully excluded, joined by
    // faces with beta > 0) without a Dirichlet face or Known term carries a constant null vector. MLMG then calls the two functions
    // below on level 0 (its residual each cycle, the bottom right-hand side, and the bottom correction, which BoomerAMG returns with an
    // arbitrary constant per region): they remove the mean per singular region over its cells, and leave excluded cells alone. On a
    // finer AMR level (only the start of the solve, where the caller has already made the rhs compatible) nothing is done.
    bool isSingular (int amrlev) const override { return amrlev == 0 && m_any_sing; }
    bool isBottomSingular () const override { return m_any_sing; }
    Vector<Real> getSolvabilityOffset (int amrlev, int mglev, MultiFab const& rhs) const override
    {
        if (amrlev != 0) { return Vector<Real>(1, Real(0)); }
        Reg const& R = m_reg[mglev];
        std::vector<Real> s(2*R.n, Real(0));
        for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
            auto const& r = rhs.const_array(mfi);
            auto const& g = R.id.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const int q = g(i,j,k);
                if (q >= 0 && R.sing[q]) { s[2*q] += r(i,j,k); s[2*q+1] += Real(1); }
            });
        }
        // coarse MG levels may live on a sub-communicator (AMReX consolidation): reduce over it, as MLCellLinOp does
        if (R.n > 0) { ParallelAllReduce::Sum(s.data(), 2*R.n, ParallelContext::CommunicatorSub()); }
        R.off.assign(R.n, Real(0));
        for (int q = 0; q < R.n; ++q) { if (s[2*q+1] > 0) { R.off[q] = s[2*q] / s[2*q+1]; } }
        return Vector<Real>(1, Real(0));
    }
    void fixSolvabilityByOffset (int amrlev, int mglev, MultiFab& rhs, Vector<Real> const&) const override
    {
        if (amrlev != 0) { return; }
        Reg const& R = m_reg[mglev];
        for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
            auto const& r = rhs.array(mfi);
            auto const& g = R.id.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const int q = g(i,j,k);
                if (q >= 0 && R.sing[q]) { r(i,j,k) -= R.off[q]; }
            });
        }
    }
    // BoomerAMG's coarsest-grid elimination needs a nonsingular bottom matrix. AMReX (pb_mlmg_bottom.adjust_singular_matrix = 1) pins
    // HYPRE row 0, the low corner of the first box of rank 0; that covers one singular region if the corner lies in it. Empty = Ok.
    std::string bottom_pin_check () const
    {
        Reg const& R = m_reg.back();
        int ns = 0;
        for (char v : R.sing) { ns += v ? 1 : 0; }
        if (ns == 0) { return std::string(); }
        if (ns > 1) { return std::to_string(ns) + " separate singular fluid regions on the bottom MG level"; }
        BoxArray const& ba = R.id.boxArray();
        DistributionMapping const& dm = R.id.DistributionMap();
        int ib = 0;
        while (ib < int(ba.size()) && dm[ib] != 0) { ++ib; }
        int q = -1;
        if (ParallelDescriptor::MyProc() == 0 && ib < int(ba.size())) { { Box const bx = ba[ib]; q = R.id[ib](bx.smallEnd()); } }
        ParallelDescriptor::ReduceIntMax(q);
        if (q < 0 || !R.sing[q]) { return "the pinned bottom row (HYPRE row 0) is not in the singular fluid region"; }
        return std::string();
    }
    void prepareForSolve () override
    {
        MLABecLaplacian::prepareForSolve();
        if (!m_chained) { chain(); m_chained = true; }
        for (int al = 0; al < NAMRLevels(); ++al) {
            for (int mg = 0; mg < NMGLevels(al); ++mg) {
                for (MFIter mfi(m_a_coeffs[al][mg]); mfi.isValid(); ++mfi) {
                    auto const& a = m_a_coeffs[al][mg].array(mfi);
                    auto const& f = frac[al][mg].const_array(mfi);
                    auto const& k = akn[al][mg].const_array(mfi);
                    amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int kk) { a(i,j,kk) = (f(i,j,kk) >= Real(1) - Real(1e-12)) ? Real(1) : k(i,j,kk); });
                }
            }
        }
        if (!m_regions) { regions_build_all(); m_regions = true; }     // MLMG::solve calls this again; the regions do not change
    }
private:
    bool m_chained = false, m_regions = false;
    struct Reg {
        int n = 0;
        iMultiFab id;                     // region of each cell (-1: fully excluded)
        std::vector<char> sing;           // per region: no Dirichlet face and no Known term
        mutable std::vector<Real> off;
    };
    Vector<Reg> m_reg;                    // per MG level of AMR level 0
    bool m_any_sing = false;
    // Connected fluid regions of MG level `mg` of AMR level 0: minimum-label propagation, then compact numbering by repeated global minimum.
    void regions_build_all ()
    {
        m_reg.resize(NMGLevels(0));
        m_any_sing = false;
        for (int mg = 0; mg < NMGLevels(0); ++mg) {
            regions_build(mg, m_reg[mg]);
            for (char v : m_reg[mg].sing) { m_any_sing = m_any_sing || v; }
        }
    }
    void regions_build (int mg, Reg& Rg)
    {
        iMultiFab& lab = Rg.id;
        MultiFab const& a = m_a_coeffs[0][mg];
        MultiFab const& f = frac[0][mg];
        auto const& bc = m_b_coeffs[0][mg];
        Geometry const& geom = m_geom[0][mg];
        Box const dom = geom.Domain();
        const int nx = dom.length(0), ny = dom.length(1);
        lab.define(a.boxArray(), a.DistributionMap(), 1, 1);
        lab.setVal(-1);
        for (MFIter mfi(lab); mfi.isValid(); ++mfi) {
            auto const& g = lab.array(mfi);
            auto const& ff = f.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                if (ff(i,j,k) < Real(1) - Real(1e-12)) { g(i,j,k) = (i - dom.smallEnd(0)) + nx*((j - dom.smallEnd(1)) + ny*(k - dom.smallEnd(2))); }
            });
        }
        for (bool changed = true; changed; ) {
            lab.FillBoundary(geom.periodicity());
            changed = false;
            for (MFIter mfi(lab); mfi.isValid(); ++mfi) {
                auto const& g = lab.array(mfi);
                Array4<Real const> b[3] = {bc[0].const_array(mfi), bc[1].const_array(mfi), bc[2].const_array(mfi)};
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    if (g(i,j,k) < 0) { return; }
                    for (int d = 0; d < 3; ++d) {
                        for (int side = 0; side < 2; ++side) {
                            IntVect q(i,j,k), fc(i,j,k); q[d] += side ? 1 : -1; fc[d] += side;
                            if (!geom.isPeriodic(d) && !dom.contains(q)) { continue; }
                            const int v = g(q[0],q[1],q[2]);
                            if (v >= 0 && v < g(i,j,k) && b[d](fc[0],fc[1],fc[2]) > Real(0)) { g(i,j,k) = v; changed = true; }
                        }
                    }
                });
            }
            ParallelDescriptor::ReduceBoolOr(changed);
        }
        // compact numbering: region q = the q-th smallest root label
        Rg.n = 0; Rg.sing.clear();
        iMultiFab id(lab.boxArray(), lab.DistributionMap(), 1, 0);
        id.setVal(-1);
        for (;;) {
            int mn = std::numeric_limits<int>::max();
            for (MFIter mfi(lab); mfi.isValid(); ++mfi) {
                auto const& g = lab.const_array(mfi); auto const& o = id.const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { if (g(i,j,k) >= 0 && o(i,j,k) < 0) { mn = std::min(mn, g(i,j,k)); } });
            }
            ParallelDescriptor::ReduceIntMin(mn);
            if (mn == std::numeric_limits<int>::max()) { break; }
            bool fixed = false;                     // a Known term or an open Dirichlet domain face
            for (MFIter mfi(lab); mfi.isValid(); ++mfi) {
                auto const& g = lab.const_array(mfi); auto const& o = id.array(mfi); auto const& aa = a.const_array(mfi);
                Array4<Real const> b[3] = {bc[0].const_array(mfi), bc[1].const_array(mfi), bc[2].const_array(mfi)};
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    if (g(i,j,k) != mn) { return; }
                    o(i,j,k) = Rg.n;
                    if (aa(i,j,k) > Real(0)) { fixed = true; }
                    for (int d = 0; d < 3; ++d) {
                        if (geom.isPeriodic(d)) { continue; }
                        if (i*(d==0) + j*(d==1) + k*(d==2) == dom.smallEnd(d) && m_lobc[0][d] == LinOpBCType::Dirichlet && b[d](i,j,k) > Real(0)) { fixed = true; }
                        IntVect fc(i,j,k); fc[d] += 1;
                        if (fc[d] == dom.bigEnd(d) + 1 && m_hibc[0][d] == LinOpBCType::Dirichlet && b[d](fc[0],fc[1],fc[2]) > Real(0)) { fixed = true; }
                    }
                });
            }
            ParallelDescriptor::ReduceBoolOr(fixed);
            Rg.sing.push_back(fixed ? 0 : 1);
            ++Rg.n;
        }
        for (MFIter mfi(lab); mfi.isValid(); ++mfi) {
            auto const& g = lab.array(mfi); auto const& o = id.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { g(i,j,k) = o(i,j,k); });
        }
    }
    // The same averaging chain as MLABecLaplacian::averageDownCoeffs (finest AMR level first: its MG levels, then into the covered
    // cells of the next coarser level's MG level 0; level 0's MG levels last).
    void chain ()
    {
        for (int al = NAMRLevels() - 1; al >= 0; --al) {
            for (int mg = 1; mg < NMGLevels(al); ++mg) {
                const IntVect r = (al > 0) ? IntVect(this->mg_coarsen_ratio) : this->mg_coarsen_ratio_vec[mg-1];
                amrex::average_down(frac[al][mg-1], frac[al][mg], 0, 1, r);
                amrex::average_down(akn[al][mg-1], akn[al][mg], 0, 1, r);
            }
            if (al > 0) {
                amrex::average_down(frac[al].back(), frac[al-1].front(), 0, 1, IntVect(this->mg_coarsen_ratio));
                amrex::average_down(akn[al].back(), akn[al-1].front(), 0, 1, IntVect(this->mg_coarsen_ratio));
            }
        }
    }
};

std::uint64_t composite_class_hash (PressureProblem const& p)
{
    std::uint64_t g = 1469598103934665603ULL;
    for (auto const& L : p.levels) {
        PressureProblem q; q.cell_class = L.cell_class;
        g ^= class_content_hash(q); g *= 1099511628211ULL; g ^= g >> 29;
    }
    return g;
}
} // namespace

struct MaskedMlmg {
    std::uint64_t hash = 0;
    bool has_known = false;
    std::unique_ptr<MaskedABecLap> op;
    std::unique_ptr<MLMG> mlmg;           // refers to *op: declared after it, destroyed first
    std::vector<std::unique_ptr<iMultiFab>> cls1;   // cell codes with one ghost cell (Known terms of the right-hand side)
    BackendStatus diag;                   // the mg_* fields
    double setup_s = 0.0;
};

namespace {
// Cell codes of every level with `ng` ghost cells: Gas 0, Solid 1, Known 2, covered by the next finer level -2, outside the domain -1.
// Ghost cells of a level > 0 outside that level take the code of the coarse cell under them (C/F), then the level's own boxes.
std::vector<std::unique_ptr<iMultiFab>> class_codes (PressureWorkspace::Impl const& W, PressureProblem const& p, std::vector<int> const& ng)
{
    std::vector<std::unique_ptr<iMultiFab>> cg(W.nlev);
    for (int l = 0; l < W.nlev; ++l) {
        cg[l] = std::make_unique<iMultiFab>(W.ba[l], W.dm[l], 1, ng[l]);
        cg[l]->setVal(-1);
        iMultiFab const* cls = p.levels[l].cell_class;
        for (MFIter mfi(*cg[l]); mfi.isValid(); ++mfi) {
            auto const& o = cg[l]->array(mfi);
            auto const& u = W.unc[l]->const_array(mfi);
            Array4<int const> c;
            if (cls) { c = cls->const_array(mfi); }
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { o(i,j,k) = (u(i,j,k) == 0) ? -2 : (cls ? c(i,j,k) : int(CellGas)); });
        }
        if (l > 0) {
            const IntVect r = W.ratio[l];
            BoxArray cba(W.ba[l]); cba.grow(ng[l]); cba.coarsen(r);
            iMultiFab ct(cba, W.dm[l], 1, 0);
            ct.setVal(-1);
            ct.ParallelCopy(*cg[l-1], 0, 0, 1, 0, 0, W.geom[l-1].periodicity());
            for (MFIter mfi(*cg[l]); mfi.isValid(); ++mfi) {
                auto const& o = cg[l]->array(mfi);
                auto const& c = ct.const_array(mfi);
                Box const vb = mfi.validbox();
                amrex::LoopOnCpu(amrex::grow(vb, ng[l]), [&] (int i, int j, int k) {
                    if (vb.contains(IntVect(i,j,k))) { return; }
                    o(i,j,k) = c(amrex::coarsen(i, r[0]), amrex::coarsen(j, r[1]), amrex::coarsen(k, r[2]));
                });
            }
        }
        cg[l]->FillBoundary(W.geom[l].periodicity());
        // FillBoundary leaves the ghost cells of a level > 0 that lie outside the periodic domain image; coarse-filled values there are kept.
    }
    return cg;
}

// Build the masked operator and its MLMG. Empty string = built; otherwise the NotBuilt reason.
std::string build_masked_mlmg (PressureWorkspace::Impl const& W, PressureProblem const& p, std::array<BC,6> const& ebc, MaskedMlmg& M)
{
    const int n = W.nlev;
    if (W.hidden >= 0) { return "masked composite MLMG with a one-cell direction is not built (use HYPRE)"; }
    std::vector<int> R(n, 1), ng(n, 1);
    for (int l = 1; l < n; ++l) { R[l] = R[l-1] * W.ratio[l][0]; }
    for (int l = 0; l < n; ++l) { ng[l] = kRunCap * R[l] + 1; }
    auto cg = class_codes(W, p, ng);

    // Rulings 4 and 5: Solid runs along each direction that start and end at a gas cell (fluid-separating), on uncovered cells.
    double tmin = std::numeric_limits<double>::max();
    int nlong = 0;
    Long thin = 0;
    std::string thin_at;
    for (int l = 0; l < n; ++l) {
        for (MFIter mfi(*cg[l]); mfi.isValid(); ++mfi) {
            auto const& c = cg[l]->const_array(mfi);
            Box const gb = amrex::grow(mfi.validbox(), ng[l]);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                if (c(i,j,k) != CellSolid) { return; }
                for (int d = 0; d < 3; ++d) {
                    IntVect lo(i,j,k); lo[d] -= 1;
                    if (c(lo[0],lo[1],lo[2]) != CellGas) { continue; }        // not the start of a run bounded by gas
                    int len = 0; bool sep = false, lng = true;
                    for (IntVect q(i,j,k); gb.contains(q); q[d] += 1) {
                        const int v = c(q[0],q[1],q[2]);
                        if (v == CellSolid) { ++len; continue; }
                        lng = false; sep = (v == CellGas); break;
                    }
                    if (lng) { ++nlong; continue; }
                    if (!sep) { continue; }
                    if (l > 0 && len < W.ratio[l][d]) {
                        if (thin++ == 0) { std::ostringstream m; m << "level " << l << " cell (" << i << "," << j << "," << k << ") direction " << d << ", " << len << " cell(s)"; thin_at = m.str(); }
                    }
                    tmin = std::min(tmin, double(len) / double(R[l]));
                }
            });
        }
    }
    ParallelDescriptor::ReduceLongSum(thin);
    ParallelDescriptor::ReduceIntSum(nlong);
    ParallelDescriptor::ReduceRealMin(tmin);
    if (thin > 0) {
        // the first offending cell is only known on its rank; the count is global
        return "masked composite MLMG is not built: " + std::to_string(thin) + " Solid wall(s) thinner than the refinement ratio between gas cells"
               + (thin_at.empty() ? std::string() : " (e.g. " + thin_at + ")") + "; averaging would erase them (pressure/13 ruling 5): use HYPRE";
    }
    // C/F: MLMG's coarse/fine interpolation (order 2) reads the coarse cell under an open C/F face and its tangential neighbours. A
    // Solid or Known coarse cell there enters with its MLMG value 0 (HYPRE drops it from the stencil), which is wrong for a
    // pressure: the obstacle may not cross the C/F interface of a level >= 1.
    Long cfbad = 0;
    std::string cf_at;
    for (int l = 1; l < n; ++l) {
        const IntVect r = W.ratio[l];
        iMultiFab inf(W.ba[l], W.dm[l], 1, 1);           // 1 on this level's cells (valid and periodic images), 0 elsewhere
        inf.setVal(0); inf.setVal(1, 0, 1, 0);
        inf.FillBoundary(W.geom[l].periodicity());
        BoxArray cba(W.ba[l]); cba.coarsen(r); cba.grow(2);
        iMultiFab cc(cba, W.dm[l], 1, 0);
        cc.setVal(-1);
        cc.ParallelCopy(*cg[l-1], 0, 0, 1, 0, 0, W.geom[l-1].periodicity());
        Box const fdom = W.geom[l].Domain();
        for (MFIter mfi(inf); mfi.isValid(); ++mfi) {
            auto const& m = inf.const_array(mfi);
            auto const& c = cg[l]->const_array(mfi);
            auto const& k0 = cc.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                if (c(i,j,k) != CellGas) { return; }
                for (int d = 0; d < 3; ++d) {
                    for (int side = 0; side < 2; ++side) {
                        IntVect q(i,j,k); q[d] += side ? 1 : -1;
                        if (m(q[0],q[1],q[2]) == 1) { continue; }                                 // fine neighbour
                        if (!W.geom[l].isPeriodic(d) && !fdom.contains(q)) { continue; }         // domain face
                        IntVect Q = amrex::coarsen(q, r);
                        if (k0(Q[0],Q[1],Q[2]) != CellGas) { continue; }                        // closed face
                        for (int t = 0; t < 3; ++t) {
                            if (t == d) { continue; }
                            for (int st = -1; st <= 1; st += 2) {
                                IntVect T = Q; T[t] += st;
                                const int v = k0(T[0],T[1],T[2]);
                                if (v == CellSolid || v == CellKnown) {
                                    if (cfbad++ == 0) { std::ostringstream m2; m2 << "level " << l << " cell (" << i << "," << j << "," << k << ") face " << d << (side ? "+" : "-"); cf_at = m2.str(); }
                                }
                            }
                        }
                    }
                }
            });
        }
    }
    ParallelDescriptor::ReduceLongSum(cfbad);
    if (cfbad > 0) {
        return "masked composite MLMG is not built: a Solid/Known coarse cell lies in the coarse/fine interpolation stencil of "
               + std::to_string(cfbad) + " open C/F face(s)" + (cf_at.empty() ? std::string() : " (e.g. " + cf_at + ")")
               + "; the obstacle crosses the C/F interface (cover it with the fine level): use HYPRE";
    }

    int maxc = -1;
    if (nlong > 0) { tmin = std::min(tmin, double(kRunCap)); }
    if (tmin < std::numeric_limits<double>::max()) { maxc = std::max(0, int(std::floor(std::log2(tmin) + 1e-12))); }

    // Level coefficients (MG level 0 of each AMR level).
    std::vector<MultiFab> a0(n), f0(n), k0(n);
    std::vector<std::array<MultiFab,3>> b0(n);
    M.has_known = false;
    for (int l = 0; l < n; ++l) {
        Box const dom = W.geom[l].Domain();
        const auto dx = W.geom[l].CellSizeArray();
        a0[l].define(W.ba[l], W.dm[l], 1, 0); f0[l].define(W.ba[l], W.dm[l], 1, 0); k0[l].define(W.ba[l], W.dm[l], 1, 0);
        for (MFIter mfi(a0[l]); mfi.isValid(); ++mfi) {
            auto const& a = a0[l].array(mfi); auto const& f = f0[l].array(mfi); auto const& kn = k0[l].array(mfi);
            auto const& c = cg[l]->const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const int v = c(i,j,k);
                if (v == CellSolid || v == CellKnown) { a(i,j,k) = Real(1); f(i,j,k) = Real(1); kn(i,j,k) = Real(0); M.has_known = M.has_known || v == CellKnown; return; }
                Real s = 0;                                     // gas (or covered: overwritten by the fine average)
                for (int d = 0; d < 3; ++d) {
                    for (int side = 0; side < 2; ++side) {
                        IntVect q(i,j,k); q[d] += side ? 1 : -1;
                        if (c(q[0],q[1],q[2]) == CellKnown) { s += Real(2) / (dx[d]*dx[d]); }
                    }
                }
                a(i,j,k) = s; f(i,j,k) = Real(0); kn(i,j,k) = s;
            });
        }
        for (int d = 0; d < 3; ++d) {
            b0[l][d].define(amrex::convert(W.ba[l], IntVect::TheDimensionVector(d)), W.dm[l], 1, 0);
            for (MFIter mfi(b0[l][d]); mfi.isValid(); ++mfi) {
                auto const& b = b0[l][d].array(mfi);
                auto const& c = cg[l]->const_array(mfi);
                const bool per = W.geom[l].isPeriodic(d);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    IntVect hi(i,j,k), lo(i,j,k); lo[d] -= 1;
                    auto open = [] (int v) { return v == CellGas || v == -2; };     // covered faces are overwritten by the fine average
                    bool g;
                    if (!per && hi[d] == dom.smallEnd(d)) { g = open(c(hi[0],hi[1],hi[2])); }
                    else if (!per && lo[d] == dom.bigEnd(d)) { g = open(c(lo[0],lo[1],lo[2])); }
                    else { g = open(c(hi[0],hi[1],hi[2])) && open(c(lo[0],lo[1],lo[2])); }
                    b(i,j,k) = g ? Real(1) : Real(0);
                });
            }
        }
    }
    ParallelDescriptor::ReduceBoolOr(M.has_known);

    LPInfo info;
    if (maxc >= 0) { info.setMaxCoarseningLevel(maxc); }
    Vector<Geometry> geoms(W.geom.begin(), W.geom.end());
    Vector<BoxArray> bas(W.ba.begin(), W.ba.end());
    Vector<DistributionMapping> dms(W.dm.begin(), W.dm.end());
    M.mlmg.reset();
    M.op = std::make_unique<MaskedABecLap>(geoms, bas, dms, info);
    M.op->setMaxOrder(2);
    Array<LinOpBCType,AMREX_SPACEDIM> lo{lin_bc(ebc[face_index(0,0)]), lin_bc(ebc[face_index(1,0)]), lin_bc(ebc[face_index(2,0)])};
    Array<LinOpBCType,AMREX_SPACEDIM> hi{lin_bc(ebc[face_index(0,1)]), lin_bc(ebc[face_index(1,1)]), lin_bc(ebc[face_index(2,1)])};
    M.op->setDomainBC(lo, hi);
    M.op->setScalars(Real(1), Real(1));
    for (int l = 0; l < n; ++l) {
        M.op->setLevelBC(l, nullptr);
        M.op->setACoeffs(l, a0[l]);
        M.op->setBCoeffs(l, Array<MultiFab const*,AMREX_SPACEDIM>{&b0[l][0], &b0[l][1], &b0[l][2]});
        M.op->set_masks(l, f0[l], k0[l]);
    }
    M.op->prepareForSolve();
    {
        const std::string e = M.op->bottom_pin_check();
        if (!e.empty()) { return "masked composite MLMG is not built: " + e + " (pressure/13 ruling 6 gauge for the BoomerAMG bottom): use HYPRE"; }
    }

    // Diagnostics of rulings 1-4 on the built hierarchy.
    BackendStatus& D = M.diag;
    D = BackendStatus();
    D.mg_max_coarsening = maxc;
    D.mg_levels0 = M.op->NMGLevels(0);
    for (int al = 0; al < n; ++al) {
        for (int mg = 0; mg < M.op->NMGLevels(al); ++mg) {
            for (MFIter mfi(M.op->m_a_coeffs[al][mg]); mfi.isValid(); ++mfi) {
                auto const& a = M.op->m_a_coeffs[al][mg].const_array(mfi);
                auto const& f = M.op->frac[al][mg].const_array(mfi);
                auto const& k = M.op->akn[al][mg].const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int kk) {
                    if (f(i,j,kk) >= Real(1) - Real(1e-12)) { ++D.mg_solid_cells; if (a(i,j,kk) != Real(1)) { ++D.mg_solid_bad; } }
                    else if (f(i,j,kk) > Real(0)) { ++D.mg_partial_cells; if (a(i,j,kk) != k(i,j,kk)) { ++D.mg_partial_shifted; } }
                });
            }
        }
    }
    for (int d = 0; d < 3; ++d) {
        for (MFIter mfi(M.op->m_b_coeffs[0][0][d]); mfi.isValid(); ++mfi) {
            auto const& b = M.op->m_b_coeffs[0][0][d].const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { if (b(i,j,k) > Real(0) && b(i,j,k) < Real(1)) { ++D.mg_half_faces; } });
        }
    }
    Long sums[5] = {Long(D.mg_solid_cells), Long(D.mg_solid_bad), Long(D.mg_partial_cells), Long(D.mg_partial_shifted), Long(D.mg_half_faces)};
    ParallelDescriptor::ReduceLongSum(sums, 5);
    D.mg_solid_cells = sums[0]; D.mg_solid_bad = sums[1]; D.mg_partial_cells = sums[2]; D.mg_partial_shifted = sums[3]; D.mg_half_faces = sums[4];

    // Ruling 6: BoomerAMG bottom (pressure/13 family; ParmParse pb_mlmg_bottom.bamg_* overrides).
    {
        ParmParse pp("pb_mlmg_bottom");
        int ct = 8, it = 6, rd = 18, ru = 18, rc = 9, ml = 25; Real st = 0.25;
        pp.queryAdd("bamg_coarsen_type", ct); pp.queryAdd("bamg_interp_type", it); pp.queryAdd("bamg_strong_threshold", st);
        pp.queryAdd("bamg_down_relax_type", rd); pp.queryAdd("bamg_up_relax_type", ru); pp.queryAdd("bamg_coarse_relax_type", rc);
        pp.queryAdd("bamg_max_levels", ml);
        int adj = 1; pp.queryAdd("adjust_singular_matrix", adj);
    }
    M.mlmg = std::make_unique<MLMG>(*M.op);
    {
        std::string bot = "hypre";              // diagnostics only: pb_mlmg_bottom.solver = smoother | bicgstab
        ParmParse("pb_mlmg_bottom").query("solver", bot);
        M.mlmg->setBottomSolver(bot == "smoother" ? BottomSolver::smoother : bot == "bicgstab" ? BottomSolver::bicgstab : BottomSolver::hypre);
    }
    M.mlmg->setHypreInterface(Hypre::Interface::ij);
    M.mlmg->setHypreOptionsNamespace("pb_mlmg_bottom");
    M.cls1.clear();
    for (int l = 0; l < n; ++l) {
        M.cls1.push_back(std::make_unique<iMultiFab>(W.ba[l], W.dm[l], 1, 1));
        iMultiFab::Copy(*M.cls1[l], *cg[l], 0, 0, 1, 1);
    }
    return std::string();
}
} // namespace
#endif

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

PressureResult solve_composite_masked (PressureProblem const& p, PressureOptions const& o, PressureWorkspace* ws, bool full, bool use_mlmg)
{
    MFIter::allowMultipleMFIters(true);
    PressureResult R;
    R.backend = use_mlmg ? "MLMG" : "HYPRE";
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

    std::array<BC,6> const ebc = effective_bc(p.bc, p.levels[0].geom.Domain());
    BackendStatus bs;
    std::vector<std::unique_ptr<MultiFab>> mrhs;          // masked MLMG: its right-hand side s (residual check)
    if (use_mlmg) {
#ifdef AMREX_USE_HYPRE
        const std::uint64_t h = composite_class_hash(p);
        const bool reused = W.mmg && W.mmg->hash == h;
        if (!reused) {
            const double t0 = amrex::second();
            auto M = std::make_shared<MaskedMlmg>();
            std::string err = build_masked_mlmg(W, p, ebc, *M);
            if (!err.empty()) { W.mmg.reset(); R.status = Status::NotBuilt; R.message = err; R.backend = ""; return R; }
            M->hash = h; M->setup_s = amrex::second() - t0;
            W.mmg = M;
        }
        MaskedMlmg& M = *W.mmg;
        // s = -b + sum over Known faces 2 g/dx^2 on uncovered gas cells, 0 elsewhere (A = -L); excluded cells start at 0.
        for (int l = 0; l < nlev; ++l) {
            mrhs.push_back(std::make_unique<MultiFab>(W.ba[l], W.dm[l], 1, 0));
            std::unique_ptr<MultiFab> kv;
            if (M.has_known) {
                kv = std::make_unique<MultiFab>(W.ba[l], W.dm[l], 1, 1);
                kv->setVal(0.0);
                if (p.levels[l].known_value) { MultiFab::Copy(*kv, *p.levels[l].known_value, 0, 0, 1, 0); }
                kv->FillBoundary(W.geom[l].periodicity());
            }
            const auto dx = W.geom[l].CellSizeArray();
            for (MFIter mfi(*mrhs[l]); mfi.isValid(); ++mfi) {
                auto const& sr = mrhs[l]->array(mfi);
                auto const& bb = b[l]->const_array(mfi);
                auto const& lab = cm.label[l]->const_array(mfi);
                auto const& c = M.cls1[l]->const_array(mfi);
                auto const& f = phi[l]->array(mfi);
                Array4<Real const> g;
                if (kv) { g = kv->const_array(mfi); }
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    if (c(i,j,k) == CellSolid || c(i,j,k) == CellKnown) { f(i,j,k) = Real(0); }
                    if (lab(i,j,k) < 0) { sr(i,j,k) = Real(0); return; }
                    Real v = -bb(i,j,k);
                    if (kv) {
                        for (int d = 0; d < 3; ++d) {
                            for (int side = 0; side < 2; ++side) {
                                IntVect q(i,j,k); q[d] += side ? 1 : -1;
                                if (c(q[0],q[1],q[2]) == CellKnown) { v += Real(2) * g(q[0],q[1],q[2]) / (dx[d]*dx[d]); }
                            }
                        }
                    }
                    sr(i,j,k) = v;
                });
            }
        }
        MLMG& mlmg = *M.mlmg;
        mlmg.setVerbose(o.verbose >= 3 ? o.verbose - 2 : 0);
        mlmg.setMaxIter(o.max_iter);
        mlmg.setConvergenceNormType(MLMGNormType::bnorm);
        mlmg.setThrowException(true);
        Vector<MultiFab const*> ps;
        for (auto const& m : mrhs) { ps.push_back(m.get()); }
        // Auto = plain V-cycles (pressure/13 ruling 7), run like None to tol/10 of the own residual (mlmg_solve).
        const Krylov kk = (o.mlmg_krylov == Krylov::Auto) ? Krylov::None : o.mlmg_krylov;
        bs = mlmg_solve(mlmg, ptrs(phi), ps, kk, Real(o.tol_rel), o.max_iter, true);
        HYPRE_ClearAllErrors();     // the BoomerAMG bottom stops at its own tolerance or iteration cap; its error flag is sticky
        const BackendStatus& D = M.diag;
        bs.mg_max_coarsening = D.mg_max_coarsening; bs.mg_levels0 = D.mg_levels0;
        bs.mg_solid_cells = D.mg_solid_cells; bs.mg_solid_bad = D.mg_solid_bad; bs.mg_partial_cells = D.mg_partial_cells;
        bs.mg_partial_shifted = D.mg_partial_shifted; bs.mg_half_faces = D.mg_half_faces;
        bs.plan_reused = reused;
        bs.setup_seconds = reused ? 0.0 : M.setup_s;
        R.backend_status = bs;
#endif
    } else {
        // HYPRE system with sealed pins.
        std::vector<HypreLayoutLevel> hl;
        for (int l = 0; l < nlev; ++l) {
            HypreLayoutLevel L;
            L.ba = W.ba[l]; L.dm = W.dm[l]; L.geom = W.geom[l]; L.ratio = W.ratio[l];
            L.cell_class = p.levels[l].cell_class;
            L.known_value = p.levels[l].known_value;
            hl.push_back(L);
        }
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
        bs = W.hsys->solve(ph, bb, o.tol_rel, o.max_iter, o.use_initial_guess);
        bs.plan_reused = reused;
        bs.setup_seconds = reused ? 0.0 : W.hsys->setup_seconds();
        R.backend_status = bs;
    }

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
        LevelMFs res;
        if (use_mlmg) {
#ifdef AMREX_USE_HYPRE
            // compResidual on a copy (it fills ghost cells): s - A phi = -(b - L phi) on gas cells, 0 on excluded cells.
            LevelMFs tmp, rr;
            for (int l = 0; l < nlev; ++l) {
                tmp.push_back(std::make_unique<MultiFab>(W.ba[l], W.dm[l], 1, 1)); tmp[l]->setVal(0.0);
                MultiFab::Copy(*tmp[l], *phv[l], 0, 0, 1, 0);
                rr.push_back(std::make_unique<MultiFab>(W.ba[l], W.dm[l], 1, 0)); rr[l]->setVal(0.0);
            }
            Vector<MultiFab const*> ps;
            for (auto const& m : mrhs) { ps.push_back(m.get()); }
            W.mmg->mlmg->compResidual(ptrs(rr), ptrs(tmp), ps);
            for (int l = 0; l < nlev; ++l) {
                res.push_back(std::make_unique<MultiFab>(W.ba[l], W.dm[l], 1, 0));
                for (MFIter mfi(*res[l]); mfi.isValid(); ++mfi) {
                    auto const& r = res[l]->array(mfi);
                    auto const& q = rr[l]->const_array(mfi);
                    auto const& lab = cm.label[l]->const_array(mfi);
                    amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { r(i,j,k) = (lab(i,j,k) < 0) ? Real(0) : -q(i,j,k); });
                }
            }
#endif
        } else {
            LevelMFs Ax;
            for (int l = 0; l < nlev; ++l) {
                Ax.push_back(std::make_unique<MultiFab>(W.ba[l], W.dm[l], 1, 0));
                Ax[l]->setVal(0.0);
            }
            // The assembled matrix pins one cell per singular component (identity row, column eliminated). Applied to a field whose
            // pin value is not 0 it misses that coupling in the pin's neighbours, so apply it to phv shifted per component to phi_pin = 0
            // (L of the true operator is unchanged by a constant on a singular component).
            LevelMFs pz;
            {
                const int nc = int(cm.comps.size());
                std::vector<double> pv(nc, 0.0);
                for (auto const& cp : cm.pins) {
                    for (MFIter mfi(*phv[cp.level]); mfi.isValid(); ++mfi) {
                        if (!mfi.validbox().contains(cp.cell)) { continue; }
                        const int id = cm.label[cp.level]->const_array(mfi)(cp.cell);
                        if (id >= 0) { pv[id] = double(phv[cp.level]->const_array(mfi)(cp.cell)); }
                    }
                }
                if (nc > 0) { ParallelDescriptor::ReduceRealSum(pv.data(), nc); }
                for (int l = 0; l < nlev; ++l) {
                    pz.push_back(std::make_unique<MultiFab>(W.ba[l], W.dm[l], 1, phv[l]->nGrow()));
                    MultiFab::Copy(*pz[l], *phv[l], 0, 0, 1, phv[l]->nGrow());
                    for (MFIter mfi(*pz[l]); mfi.isValid(); ++mfi) {
                        auto const& a = pz[l]->array(mfi);
                        auto const& lab = cm.label[l]->const_array(mfi);
                        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                            const int id = lab(i,j,k);
                            if (id >= 0 && cm.comps[id].singular) { a(i,j,k) -= Real(pv[id]); }
                        });
                    }
                }
            }
            std::vector<MultiFab const*> px; std::vector<MultiFab*> py;
            for (int l = 0; l < nlev; ++l) { px.push_back(pz[l].get()); py.push_back(Ax[l].get()); }
            W.hsys->apply(px, py);
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
        return solve_composite_masked(p, o, ws, full, !use_hypre);   // explicit MLMG: masked composite MLMG (P3); Auto stays HYPRE (D-092)
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
        refresh_robin(W);
        auto* fbp = dynamic_cast<FftBottomPoisson*>(W.mlp.get());
        // Krylov only: under plain V-cycles the FFT bottom stagnated at true residual ~1e-7 on vortex_ring_48_4m (200 iterations; harness cases converge).
        const bool krylov = o.mlmg_krylov == Krylov::BiCGSTAB || o.mlmg_krylov == Krylov::GMRES;
        const bool fft_bottom = o.mlmg_fft_bottom && krylov && fbp && !ext && fbp->fft_bottom_ok();
        mlmg.setBottomSolver(fft_bottom ? BottomSolver::custom : BottomSolver::Default);
        bs = mlmg_solve(mlmg, pphi, pb_, o.mlmg_krylov, Real(o.tol_rel), o.max_iter, fft_bottom);   // the FFT bottom is a fixed linear map: kept under Krylov
        if (fft_bottom) { bs.method += "+FFT bottom"; }
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
        refresh_robin(W);
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
    // Solid-only mask: the gradient of the unmasked operator, zeroed on every face touching a Solid cell (beta = 0 there). C/F faces
    // never touch a Solid cell (NotBuilt in the solve), so they are unchanged. Known cells stay NotBuilt.
    for (auto const& L : p.levels) { if (L.known_value) { return fail(Status::NotBuilt, "face_gradient_composite with Known cells is not built"); } }
    PressureProblem pu = p;
    for (auto& L : pu.levels) { L.cell_class = nullptr; }
    PressureWorkspace::Impl W;
    build_impl(W, pu);
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
    for (int l = 0; l < nlev; ++l) {
        if (!p.levels[l].cell_class) { continue; }
        iMultiFab c(p.levels[l].ba, p.levels[l].dm, 1, 1);
        c.setVal(0);
        iMultiFab::Copy(c, *p.levels[l].cell_class, 0, 0, 1, 0);
        c.FillBoundary(p.levels[l].geom.periodicity());
        for (int d = 0; d < 3; ++d) {
            const IntVect e = IntVect::TheDimensionVector(d);
            for (MFIter mfi(*grad_out[l][d]); mfi.isValid(); ++mfi) {
                auto const& g = grad_out[l][d]->array(mfi);
                auto const& s = c.const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { if (s(i,j,k) == 1 || s(i-e[0],j-e[1],k-e[2]) == 1) { g(i,j,k) = 0.0; } });
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
