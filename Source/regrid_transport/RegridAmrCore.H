// RegridAmrCore.H: the AmrCore subclass of Role 3 (R2a). AmrCore is the AMReX class that owns the grids of all levels (BoxArray and
// DistributionMapping per level) and calls back to make, remake and clear a level and to tag cells for refinement.
//
// R2a builds the STATIC hierarchy: the grids come from fdsrt::Hierarchy (R1: finer &MESH boxes plus nesting boxes), not from tagging.
// Every created, changed or removed level is reported to a fdsrt::LevelListener (RegridInterface.H; the driver's LevelRegistry implements it).
// R4 adds the dynamic part: ErrorEst runs the caller's tag function, force-tags the static finer &MESH footprints and clips the tags to the refinable region
// (IR-008); init_from_tags builds the t = 0 hierarchy level by level (ruling: initial conditions evaluated directly on each new level); regrid_dynamic
// rebuilds the finer levels between full steps. The data of new and remade levels are filled through a LevelDataTransfer.
#ifndef FDSRT_REGRIDAMRCORE_H_
#define FDSRT_REGRIDAMRCORE_H_

#include <AMReX_AmrCore.H>
#include <AMReX_iMultiFab.H>

#include <functional>
#include <vector>

#include "Hierarchy.H"
#include "RegridInterface.H"

namespace fdsrt {

// amrex::AmrInfo (ratios, blocking factors, max grid size, nesting width, per direction) from the hierarchy and the &AMR parameters.
amrex::AmrInfo make_amr_info(const Hierarchy& h, const AmrParams& p);
// Geometry of level 0 (physical box = bounding box of the level-0 meshes, periodicity from the hierarchy).
amrex::Geometry make_level0_geometry(const Hierarchy& h);
amrex::Box to_amrex(const IBox& b);

// Counters of the dynamic part (per level index of the level that was tagged; the tags create/keep level+1).
struct RegridStats {
    int n_regrids = 0;                                  // calls of regrid_dynamic
    int n_changed = 0;                                  // of those, regrids that changed some level
    std::vector<long> tags_set;                         // tags set after the criteria, the forced boxes and the clip (cells), summed over regrids
    std::vector<long> discarded_outside_region;         // criterion tags removed because they lie outside the refinable region (cells; IR-008 R1-R5)
    std::vector<long> buffer_blocks_discarded;          // blocks added by the tag buffer outside the region and removed (blocks of bf/ratio cells)
};

class RegridAmrCore : public amrex::AmrCore {
public:
    // Fills `tags` of level `lev` from the fields of that level (field ghost layer filled by the caller). Static boxes and the region clip are applied afterwards.
    using TagFn = std::function<void(int lev, amrex::TagBoxArray& tags, amrex::Real time)>;

    // Distribution of the grids of `level` over ranks. Default (empty function): AMReX space-filling-curve mapping, deterministic.
    using DmFn = std::function<amrex::DistributionMapping(int level, const amrex::BoxArray&)>;

    RegridAmrCore(const Hierarchy& h, const AmrParams& p);

    // Installs the static hierarchy: SetBoxArray / SetDistributionMap / SetFinestLevel and calls listener.make_level for levels 1..top, and
    // for level 0 when notify_level0 is true (the driver owns level 0 and adopts it, so the driver passes false and level0_dm = its mapping).
    void init_static(LevelListener& listener, bool notify_level0 = false, const DmFn& dm_fn = DmFn(),
                     const amrex::DistributionMapping* level0_dm = nullptr);

    // ---- R4 dynamic hierarchy -------------------------------------------------------------------------------------------------------------------
    void set_tag_function(TagFn f) { tag_fn_ = std::move(f); }
    void set_data_transfer(LevelDataTransfer* t) { transfer_ = t; }

    // t = 0: level 0 is given (the driver's meshes; `level0_dm` its mapping, default AMReX mapping); finer levels are built one at a time: tag the finest
    // level, make the next one (make_level + fill_initial_level), repeat until MAX_LEVEL or no tags; then transfer.hierarchy_done(true).
    void init_from_tags(LevelListener& listener, amrex::Real time, bool notify_level0 = false, const DmFn& dm_fn = DmFn(),
                        const amrex::DistributionMapping* level0_dm = nullptr);

    // Between full steps: re-tag all levels, rebuild the grids of levels 1.. and fill the data (begin_regrid, make/remake/clear, end_regrid, hierarchy_done).
    // Returns true if the grids of some level changed. Does nothing (false) when MAX_LEVEL = 0.
    bool regrid_dynamic(amrex::Real time);

    // 1 where the cell of level `lev` lies under level lev+1 (for the TAG_KEEP hysteresis); all zero if lev is the finest level.
    amrex::iMultiFab covered_mask(int lev) const;

    const RegridStats& stats() const { return stats_; }

    LevelLayout layout(int level) const;     // what the registry gets for `level`
    const Hierarchy& hierarchy() const { return h_; }

protected:
    // AmrCore callbacks. R2a forwards to the listener; the data transfer of a real regrid is R4.
    void MakeNewLevelFromScratch(int lev, amrex::Real time, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) override;
    void MakeNewLevelFromCoarse(int lev, amrex::Real time, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) override;
    void RemakeLevel(int lev, amrex::Real time, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm) override;
    void ClearLevel(int lev) override;
    void ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;   // R4: tag function + static boxes + region clip
public:
    // AmrMesh calls this after the tag buffer and the coarsening by blocking factor / ratio: the buffer must not reach outside the refinable region.
    void ManualTagsPlacement(int lev, amrex::TagBoxArray& tags, const amrex::Vector<amrex::IntVect>& bf_lev) override;
protected:
    amrex::DistributionMapping MakeDistributionMap(int lev, amrex::BoxArray const& ba) override;

private:
    Hierarchy h_;
    LevelListener* listener_ = nullptr;
    LevelDataTransfer* transfer_ = nullptr;
    TagFn tag_fn_;
    DmFn dm_fn_;
    amrex::DistributionMapping level0_dm_;
    bool have_level0_dm_ = false;
    RegridStats stats_;
    std::vector<bool> reported_discard_;
};

}  // namespace fdsrt

#endif
