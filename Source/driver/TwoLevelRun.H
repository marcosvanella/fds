// TwoLevelRun.H: `fds_amr <case.fds> --two-level-run`: the real two-level run (S14). Level 0 is the FDS case, level 1 a static refined patch made by Role 3's RegridAmrCore
// (the hierarchy of the Hierarchy / AmrInput code, listener = LevelRegistry, data = RegistryTransfer), bound through TimeLoop::bind_level; the time loop runs the composite step
// (composite pressure solve, interface flux overwrite, face average-down) and the post-regrid projection (D-063, Role 3's PostRegridProjection with the pressure_backend composite solver)
// makes the initial two-level velocity satisfy div u = D.
#ifndef FDS_TWOLEVELRUN_H_
#define FDS_TWOLEVELRUN_H_

#include <string>

#include "FdsAmr.H"
#include "TimeLoop.H"

namespace fdsamr {

struct TwoLevelOptions {
    int steps = 20;
    int ratio = 2;                   // refinement ratio in the directions with more than one cell
    int patch[4] = {4, 11, 4, 11};   // coarse cells of the fine patch: x lo, x hi, z lo, z hi (y: the whole domain)
    int maxsize = 16;                // MAX_GRID_SIZE of level 1 (fine cells)
    int blocking = 4;                // BLOCKING_FACTOR
    std::string projection = "AUTO"; // POST_REGRID_PROJECTION of the initial projection: ON / AUTO / OFF
    bool overwrite = true;           // interface flux overwrite
    std::string outdir = ".";
    std::string chid = "two_level";
    int log_every = 1;
};

// Returns the number of failed checks (0 = completed). Prints one line per step: composite mass, drift, max |div u - D| by level.
int two_level_run(TimeLoop& loop, const Level0& l0, double dt_setup, const TwoLevelOptions& o);

}  // namespace fdsamr

#endif
