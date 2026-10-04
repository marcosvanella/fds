// LevelRegistry.H: the per-level objects of the driver (S9; Role 3 plan section 4 items 1, 2, 5; RegridInterface.H).
//
// One LevelRegistry per run. It holds, for every AMR level, the Level (layout), the Fields (registered field MultiFabs on that level's BoxArray/DistributionMapping,
// D-031 ghost widths), the SideData (solid flag, face mask, source flag) and the covered-by-finer mask. It implements fdsrt::LevelListener, so Role 3's AmrCore callbacks
// (make_level / remake_level / clear_level) construct and destroy the per-level objects at regrid. Level 0 is the FDS &MESH set: TimeLoop adopts the objects it owns
// (adopt_level0, non-owning); level 0 layout is fixed by the meshes and cannot be remade.
//
// What a level > 0 has: layout, Fields, SideData (rebuilt from the layout: no obstructions, every face between boxes of the level or towards the coarser level is
// open, the domain edge is a wall), the covered mask. The FDS binding is made by TimeLoop::bind_level (fds_fine_level.f90: a fine-level MESH_TYPE with its metrics, cell table and
// aliased fields; patch 0007/0008 route), after which Level::fds_bound() is true; a level that was not bound still aborts a stage request with a clear message. Not built yet for a fine
// level: wall cells at the domain edge (a fine box on the edge has no wall cells; its edge ghost layer is filled by the mirror rule of BcStep, periodic and Neumann/mirror faces only), wall cells of obstructions,
// the INTERPOLATED_MESH mask, zone sums (notes/level-binding.md, notes/two-level-run.md).
//
// New fine level, derived stage fields (architect ruling): RegistryTransfer transfers the conserved state (rho*Z, faces); the stage arrays D, DS, RSUM, MU, KRES, H, HS of a new level are
// NOT transferred. FDS skips the density/EOS update in the first cycle, so a new level would divide by RSUM = 0 in DIVERGENCE_PART_1. Rule: RegistryTransfer::derive injects the PARENT value
// (piecewise constant, the coarse cell under each fine cell) into those seven arrays of the new level, before the first stage runs; the stage then recomputes them from the fine state.
// Not derived yet: TMP/RSUM from the fine ZZ by the EOS (the injected parent values stand in until the first density update). TwoLevelRun.cpp sets the hook; notes/two-level-run.md.
//
// Kernel-facing rules (M2a): (a) passive scalars are the Fields' business (ncomp = N_TOTAL_SCALARS); (b) only uniform Cartesian metrics (Level::dx).
#ifndef FDS_LEVELREGISTRY_H_
#define FDS_LEVELREGISTRY_H_

#include <AMReX_iMultiFab.H>

#include <functional>
#include <map>
#include <memory>
#include <vector>

#include "Fields.H"
#include "RegridInterface.H"
#include "SideData.H"

namespace fdsamr {

// SideData provider for a layout-only level: SOLID = 0 everywhere; a face of a valid cell is kWall at the domain edge (periodic edges included, as in FDS),
// kInterface (open: SideData turns it into a no-wall face) wherever the neighbouring cell is inside the domain but in another box: a same-level box-to-box face and
// a coarse-fine face are both open faces. Obstructions at a fine level are not representable (no OBST, Role 3 plan item 5).
CellWallProvider layout_cell_walls(const Level& lev);

// Mask of the valid cells of `ba` that are covered by the finer level: 1 where a cell of ba, refined by `ratio`, is inside a box of `fine_ba` (the fine level's BoxArray).
// iMultiFab on (ba, dm), 1 component, no ghost cells.
amrex::iMultiFab build_covered_mask(const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, const amrex::BoxArray& fine_ba, const amrex::IntVect& ratio);

class LevelRegistry : public fdsrt::LevelListener {
public:
    // dom: the domain data of level 0 (periodicity, scalar counts); nscalars: N_TOTAL_SCALARS; names: fields to allocate per fine level (default_field_names() if empty).
    LevelRegistry(const DomainInfo& dom, int nscalars, const std::vector<std::string>& names = {});

    // Level 0: non-owning view of the objects TimeLoop owns. Must be called once, before any make_level.
    void adopt_level0(const Level& l0, Fields& F, SideData& sd);

    // fdsrt::LevelListener (called between time steps only). make_level(0) / remake_level(0) are accepted only if the layout equals the adopted level-0 layout.
    void make_level(const fdsrt::LevelLayout& l) override;
    void remake_level(const fdsrt::LevelLayout& l) override;   // the previous level objects stay readable through retired(level) until release_retired(level)
    void clear_level(int level) override;
    // D-058 (R4): bracket of one dynamic regrid. Between begin_regrid() and end_regrid() nothing that make_level / remake_level / clear_level replaces or removes is freed: remake_level
    // keeps the previous objects of the level (retired_*), clear_level keeps the removed levels (retired_level / retired_fields work for them as well), and end_regrid() frees
    // everything retired, all levels, so that Role 3 can copy from every old field before any of them goes away. Outside the bracket the behaviour is the S9 one (remake_level
    // retires until release_retired(l); clear_level frees at once). Nesting is refused.
    void begin_regrid() override;
    void end_regrid() override;
    bool in_regrid() const { return m_in_regrid; }

    int num_levels() const { return static_cast<int>(m_lv.size()); }   // levels 0..finest
    bool has_level(int l) const { return l >= 0 && l < num_levels() && m_lv[l].F != nullptr; }
    const Level& level(int l) const;
    Fields& fields(int l);
    const Fields& fields(int l) const;
    SideData& side_data(int l);

    // Per-level SideData rebuild, callable at regrid (make_level/remake_level call it; Role 3 may call it again after it changed the boxes of the finer level).
    // Uses the provider set with set_side_provider (a fine-level FDS provider once such objects exist) or layout_cell_walls. Level 0: the provider of the FDS set-up
    // is not available here; a level-0 rebuild is refused (TimeLoop owns it).
    void rebuild_side_data(int l);
    void set_side_provider(int l, CellWallProvider p) { ensure(l); m_lv[l].provider = std::move(p); }

    // Fine mask of level l: nonzero = covered by level l+1 (nullptr for the finest level, or when l+1 does not exist). Recomputed whenever a level is made, remade
    // or cleared. Use with exact_sum_uncovered.
    const amrex::iMultiFab* covered_mask(int l) const;

    // True when level l has FDS MESH_TYPE objects behind its boxes (always for level 0). false for every level > 0 today.
    bool fds_bound(int l) const { return has_level(l) && m_lv[l].lev->fds_bound(); }

    // After remake_level: the previous Fields/Level of that level, for Role 3's transfer. release_retired frees them. nullptr when nothing is retired.
    const Level* retired_level(int l) const;
    Fields* retired_fields(int l);
    void release_retired(int l);   // also frees what clear_level retired inside a regrid bracket for that level

    // D-058 (R4): "initial fields on level L". t = 0 hierarchy: the initial conditions are evaluated DIRECTLY at the cell centres and face centres of the new level (not interpolated);
    // the caller then averages down once (the Role 3 DriverAdapter average_down_registry) and refills ghosts. `cell` gives rho, T and the N_TOTAL_SCALARS mass fractions at a cell
    // centre, `velocity(dir, x, y, z)` the normal velocity at a face centre (optional: faces stay 0). Fills the valid region of RHO, RHOS (= RHO), TMP, ZZ, ZZS (= ZZ) and, with
    // `velocity`, U/V/W and US/VS/WS (= U/V/W) of level l >= 1 (level 0 is the FDS initialisation); fields that the registry does not carry are skipped. Ghost cells are not touched.
    // x, y, z are physical coordinates: ProbLo + (index + 1/2) dx for a cell, ProbLo + index dx in the normal direction of a face (AMReX face index = low face of the cell of that index).
    struct InitialFill {
        std::function<void(double x, double y, double z, double* rho, double* tmp, double* zz /* [N_TOTAL_SCALARS] */)> cell;
        std::function<double(int dir, double x, double y, double z)> velocity;
    };
    void fill_initial_level(int l, const InitialFill& f);

    // Counters for tests and logs.
    int n_make = 0, n_remake = 0, n_clear = 0, n_side_rebuilds = 0, n_begin_regrid = 0, n_end_regrid = 0;

private:
    struct Slot {
        std::unique_ptr<Level> lev_own;
        std::unique_ptr<Fields> F_own;
        std::unique_ptr<SideData> sd_own;
        const Level* lev = nullptr;
        Fields* F = nullptr;
        SideData* sd = nullptr;
        std::unique_ptr<amrex::iMultiFab> covered;
        CellWallProvider provider;
        // retired by remake_level
        std::unique_ptr<Level> old_lev;
        std::unique_ptr<Fields> old_F;
        std::unique_ptr<SideData> old_sd;
    };
    void ensure(int l) { if (l >= num_levels()) m_lv.resize(l + 1); }
    void build_slot(const fdsrt::LevelLayout& l, Slot& s);
    void refresh_covered();
    DomainInfo m_dom;
    int m_ns;
    std::vector<std::string> m_names;
    std::vector<Slot> m_lv;
    bool m_have0 = false;
    bool m_in_regrid = false;
    struct Cleared { std::unique_ptr<Level> lev; std::unique_ptr<Fields> F; std::unique_ptr<SideData> sd; };
    std::map<int, Cleared> m_cleared;   //!< levels removed by clear_level inside a regrid bracket, kept readable until end_regrid
};

}  // namespace fdsamr

#endif
