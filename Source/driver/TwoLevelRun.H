// TwoLevelRun.H: `fds_amr <case.fds> --two-level-run`: multi-level run (S14 / D-053). Level 0 is the FDS case; levels 1..max_level are static nested refined patches made by Role 3's
// RegridAmrCore (Hierarchy / AmrInput, listener = LevelRegistry, data = RegistryTransfer), bound through TimeLoop::bind_level; the time loop runs the composite step (Role 2
// PressureIface over all bound levels, interface flux overwrite, face average-down) and the post-regrid projection (D-063) makes the initial hierarchy velocity satisfy div u = D.
// Default max_level = 1 (two levels, as before). With --levels 3 (max_level = 2) a nested inner patch (patch2) is tagged on level 1.
#ifndef FDS_TWOLEVELRUN_H_
#define FDS_TWOLEVELRUN_H_

#include <cmath>
#include <limits>
#include <string>

#include "FdsAmr.H"
#include "TimeLoop.H"

namespace fdsamr {

struct TwoLevelOptions {
    int steps = 20;
    int max_level = 1;               // finest level index (1 = two levels; 2 = three levels). Levels are bound 1..max_level.
    int ratio = 2;                   // refinement ratio in the directions with more than one cell (same at every level)
    int patch[4] = {4, 11, 4, 11};   // level-0 cells tagged to create level 1: x lo, x hi, z lo, z hi (y: the whole domain)
    int patch2[4] = {6, 9, 6, 9};    // level-0 cells of the nested region tagged on levels >= 1 (must lie inside patch); used when max_level >= 2
    int maxsize = 16;                // MAX_GRID_SIZE of fine levels (fine cells)
    int blocking = 4;                // BLOCKING_FACTOR
    std::string projection = "AUTO"; // POST_REGRID_PROJECTION of the initial projection: ON / AUTO / OFF
    bool overwrite = true;           // interface flux overwrite
    std::string outdir = ".";
    std::string chid = "two_level";
    int log_every = 1;
    int boundary_test = 0;           // 1: call stage_boundary(1,3) and (1,6) on the bound fine level and check level 0 is untouched; 2: same with the pre-S14.1 level-0-only guard (must stop the run)
    bool bind_only = false;          // true: build hierarchy, bind all levels, report layout, skip advance (registry/bind gate when a full step is not wanted)
    int regrid_every = 0;            // CLI --regrid-every N -> AmrParams.regrid_interval (0 = static after setup; >0 mid-run gate)
    std::string plotfile;            // non-empty: WriteMultiLevelPlotfile here after the run (RHO + cell-centred U,V,W); CLI --plotfile [DIR]
    int plot_every = 0;              // CLI --plot-every N: when N>0 and plotfile set, after step s with s%N==0 write plt_<chid>_s<s> (keeps mid-regrid + final dumps)
    // TAG_VORT path (vortex-ring owner-show). NaN = keep blob TAG_SPEC_* emit when regrid_every>0.
    // tag_vort_auto: thr = 0.45 * max|VORT| on level 0 after the first undivided-|ω| fill (before initial regrid); TAG_KEEP=0.85 (tight core).
    double tag_vort = std::numeric_limits<double>::quiet_NaN();
    bool tag_vort_auto = false;
};

// Returns the number of failed checks (0 = completed). Prints one line per step: composite mass, drift, max |div u - D| by level.
// With o.max_level >= 2 the same path builds and binds three or more levels (D-053).
int two_level_run(TimeLoop& loop, const Level0& l0, double dt_setup, const TwoLevelOptions& o);

}  // namespace fdsamr

#endif
