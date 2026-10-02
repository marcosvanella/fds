// LevelRegistry.H: the per-level objects of the driver (S9; Role 3 plan section 4 items 1, 2, 5; RegridInterface.H).
//
// One LevelRegistry per run. It holds, for every AMR level, the Level (layout), the Fields (registered field MultiFabs on that level's BoxArray/DistributionMapping,
// D-031 ghost widths), the SideData (solid flag, face mask, source flag) and the covered-by-finer mask. It implements fdsrt::LevelListener, so Role 3's AmrCore callbacks
// (make_level / remake_level / clear_level) construct and destroy the per-level objects at regrid. Level 0 is the FDS &MESH set: TimeLoop adopts the objects it owns
// (adopt_level0, non-owning); level 0 layout is fixed by the meshes and cannot be remade.
//
// What a level > 0 has today: layout, Fields, SideData (rebuilt from the layout: no obstructions, every face between boxes of the level or towards the coarser level is
// open, the domain edge is a wall), the covered mask. What it does NOT have is an FDS binding: the kernels take an FDS mesh number, and a fine-level FDS MESH_TYPE
// object (metrics, wall tables, CELL_INDEX, POINT_TO_BOX route of patch 0005 DRAFT) does not exist (notes/level-interface.md section 3). Level::fds_bound() is false for it,
// and TimeLoop aborts with that message if a stage is requested on it. No fine-level physics is faked here.
//
// Kernel-facing rules (M2a): (a) passive scalars are the Fields' business (ncomp = N_TOTAL_SCALARS); (b) only uniform Cartesian metrics (Level::dx).
#ifndef FDS_LEVELREGISTRY_H_
#define FDS_LEVELREGISTRY_H_

#include <AMReX_iMultiFab.H>

#include <functional>
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
    void release_retired(int l);

    // Counters for tests and logs.
    int n_make = 0, n_remake = 0, n_clear = 0, n_side_rebuilds = 0;

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
};

}  // namespace fdsamr

#endif
