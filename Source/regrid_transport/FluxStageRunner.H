// FluxStageRunner.H: the interface flux overwrite (D-050, D-061) as a step of one time stage, Role 3 side. It sits between the driver's stage pieces: after the
// stage fluxes of a kind have been read out on every level it builds the override lists of every level from the level above it (finest level first) and hands them
// to the driver; with `overwrite` off it clears them instead (negative control, and the FDS behaviour of the FR-016 baseline). It never runs a kernel itself.
//
// One stage, in the order FDS has (the caller runs the pieces; see TimeLoopStages.H for the version on the driver's TimeLoop):
//   1. viscosity + mass finite differences, all levels;
//   2. ADV read-out, all levels                                    -> runner.set_overrides(fa, FluxKind::Adv)
//   3. density (cell update with the ADV overrides), all levels; exchange / boundary / wall steps (the coarse-fine ghost hook runs here)
//   4. DIVERGENCE_PART_1 read-out of the diffusive fluxes, all levels  -> runner.set_overrides(fa, FluxKind::Dif)
//   5. DIVERGENCE_PART_1 again where an override was set (DEL_RHO_D_DEL_Z from the overwritten fluxes).
// Plain words: a "face flux" is the amount of a quantity that crosses one cell face per unit area and time. A coarse face on the outline of a fine patch gets the
// area-weighted mean of the fine face fluxes that tile it; that makes the amount leaving the coarse cell equal to the amount entering the fine cells, so nothing is
// created or lost at the interface and no flux register (refluxing) is needed when every level takes the same time step.
#ifndef FDSRT_FLUXSTAGERUNNER_H_
#define FDSRT_FLUXSTAGERUNNER_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>

#include <vector>

#include "FluxOverrideOps.H"
#include "RegridInterface.H"

namespace fdsrt {

// Grids of one level, as the runner needs them (the driver's Level / the AmrCore data).
struct StageLevel {
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    amrex::IntVect ratio_from_parent{1, 1, 1};   // 1 in a single-cell (hidden) direction; (1,1,1) for level 0
};

struct StageStats {
    long adv_entries = 0, dif_entries = 0;       // face entries set on this rank, summed over levels and boxes (a face shared by two coarse boxes counts twice)
    long calls = 0;                              // calls of set_overrides
};

class FluxStageRunner {
public:
    // levels[0] is the coarsest. overwrite = false: every call clears the lists (the empty set; the kernels are then bitwise the plain FDS ones).
    FluxStageRunner(std::vector<StageLevel> levels, bool overwrite = true) : lv_(std::move(levels)), overwrite_(overwrite) {}

    int num_levels() const { return static_cast<int>(lv_.size()); }
    bool overwrite() const { return overwrite_; }
    void set_overwrite(bool on) { overwrite_ = on; }
    const StageStats& stats() const { return stats_; }

    // Reads fa.stage_flux(l+1, kind, d) for l = finest-1 .. 0 and calls fa.set_flux_override(l, kind, lists) for EVERY level 0..finest (the finest level gets an
    // empty set, which also clears the lists of the previous stage). Valid after the read-out of `kind` on all levels (see the order above).
    void set_overrides(FluxAccess& fa, FluxKind kind);

    // The lists a level would get from level l+1's fluxes, without handing them over (tests).
    std::vector<std::vector<FluxOverride>> lists_for_level(const FluxAccess& fa, int level, FluxKind kind, OverrideStats* st = nullptr) const;

private:
    std::vector<StageLevel> lv_;
    bool overwrite_;
    StageStats stats_;
};

}  // namespace fdsrt

#endif
