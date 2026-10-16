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
#include "TagOps.H"
#include "../pressure_backend/ObstMask.H"

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

// What a regrid did, for the post-regrid hook (ruling D-063).
struct RegridOutcome {
    bool changed = false;               // the grids of some level changed (or the number of levels)
    bool retained_and_new = false;      // some level >= 1 kept part of its old grids AND gained new cells: new fine faces next to retained ones
    long retained_cells = 0;            // cells of levels >= 1 present in both the old and the new grids
    long new_cells = 0;                 // cells of levels >= 1 in the new grids that were not in the old ones
    int finest = 0;
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
    // Preferred path: AmrParams TAG_* keys + set_tag_fields (ErrorEst calls apply_tag_criteria). set_tag_function is the fallback for unit tests
    // that build a TagCriterion by hand (used only when no TAG_* criterion is set).
    void set_tag_function(TagFn f) { tag_fn_ = std::move(f); }
    // Also apply set_tag_function on top of the TAG_* criteria (driver: step-0 forced initial patch, --init-patch; off again after the initial regrid)
    void set_force_tag_function(bool on) { force_tag_fn_ = on; }
    // --tag-obst NBUF: force-tag every OBST of pb::obst_list() (ObstMask.H) on each level below the top, rounded outward and grown by NBUF level-0 cells (<0: off)
    void set_tag_obst(int nbuf) { tag_obst_nbuf_ = nbuf; }
    // --tag-obst-snap (default on for two or more levels): tag, on each level below the top, every OBST whose FDS snap at the next level is not the refinement of its snap at this level
    // (including OBSTs that are thin at this level), rounded outward and grown by 1 cell, so a coarse-fine interface never cuts an OBST whose representation differs between the levels
    void set_tag_obst_snap(bool on) { tag_obst_snap_ = on; }
    // D-010: no tags within NBUF level-0 cells of a stateful wall (pb::nofine_list; -1 off). Clustering then takes only tagged cells
    // (grid efficiency 1), else a box around a tagged region would cover the cleared hole.
    void set_nofine(int nbuf) { nofine_nbuf_ = nbuf; if (nbuf >= 0 && !pb::nofine_list().empty()) SetGridEff(1.0); }
    using TagFieldsFn = std::function<TagFields(int lev)>;
    void set_tag_fields(TagFieldsFn f) { tag_fields_fn_ = std::move(f); }
    const AmrParams& amr_params() const { return params_; }
    void set_data_transfer(LevelDataTransfer* t) { transfer_ = t; }

    // D-080: by default we buffer tags ourselves (ghost exchange + dilation) and pass n_error_buf = 0 to AMReX so it does not
    // double-grow. set_use_amrex_tag_buffer(true) restores AMReX TagBoxArray::buffer (incomplete across box edges; negative-control only).
    void set_use_amrex_tag_buffer(bool use);

    // t = 0: level 0 is given (the driver's meshes; `level0_dm` its mapping, default AMReX mapping); finer levels are built one at a time: tag the finest
    // level, make the next one (make_level + fill_initial_level), repeat until MAX_LEVEL or no tags; then transfer.hierarchy_done(true).
    void init_from_tags(LevelListener& listener, amrex::Real time, bool notify_level0 = false, const DmFn& dm_fn = DmFn(),
                        const amrex::DistributionMapping* level0_dm = nullptr);

    // Between full steps: re-tag all levels, rebuild the grids of levels 1.. and fill the data (begin_regrid, make/remake/clear, end_regrid, hierarchy_done).
    // Returns true if the grids of some level changed. Does nothing (false) when MAX_LEVEL = 0.
    bool regrid_dynamic(amrex::Real time);

    // Post-regrid hook (D-063): called at the end of regrid_dynamic, after hierarchy_done(false) and before control returns to the time loop (so before the next dt),
    // when POST_REGRID_PROJECTION asks for it: AUTO (default) and ON whenever the grids changed, OFF never. The driver installs
    // the composite projection here (PostRegridProjection.H, DriverAdapter.H: install_post_regrid_projection).
    using PostRegridHook = std::function<void(const RegridOutcome&)>;
    void set_post_regrid_hook(PostRegridHook f) { post_hook_ = std::move(f); }
    const RegridOutcome& last_outcome() const { return last_outcome_; }
    int hook_calls() const { return hook_calls_; }
    // The decision alone (also used by the tests): does the setting ask for the projection after this outcome?
    static bool projection_wanted(PostRegridProjection mode, const RegridOutcome& o);

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
    AmrParams params_;
    LevelListener* listener_ = nullptr;
    LevelDataTransfer* transfer_ = nullptr;
    TagFn tag_fn_;
    bool force_tag_fn_ = false;
    int tag_obst_nbuf_ = -1;
    bool tag_obst_snap_ = false;
    int nofine_nbuf_ = -1;
    TagFieldsFn tag_fields_fn_;
    DmFn dm_fn_;
    amrex::DistributionMapping level0_dm_;
    bool have_level0_dm_ = false;
    RegridStats stats_;
    int n_error_buf_ = 1;                 // &AMR N_ERROR_BUF (semantics unchanged; applied in ErrorEst via buffer_tags)
    bool own_tag_buffer_ = true;          // false = AMReX TagBoxArray::buffer path (negative control)
    PostRegridProjection post_mode_ = PostRegridProjection::Auto;
    PostRegridHook post_hook_;
    RegridOutcome last_outcome_;
    int hook_calls_ = 0;
    std::vector<bool> reported_discard_;
};

}  // namespace fdsrt

#endif
