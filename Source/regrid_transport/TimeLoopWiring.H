// TimeLoopWiring.H: installs Role 3's hooks into Role 1's TimeLoop (Source/driver/TimeLoop.H, read only). Kept apart from DriverAdapter because it names TimeLoop,
// whose implementation needs the FDS Fortran objects: the unit tests use DriverAdapter and the real BcStep only; this file is compiled by the driver build.
#ifndef FDSRT_TIMELOOPWIRING_H_
#define FDSRT_TIMELOOPWIRING_H_

#include "DriverAdapter.H"
#include "TimeLoop.H"

namespace fdsrt {

// Hierarchy installed (RegridAmrCore::init_static done): set the FDS ghost-rule hook on levels 1..finest. The hook runs inside TimeLoop::stage_exchange for every exchange
// code, after the same-level fill and before the boundary routines. Lifetime: `th` and `stats` are copied/stored by pointer; `stats` must outlive the loop.
void install_cf_ghost_hooks(fdsamr::TimeLoop& loop, int finest_level, const ThermoProvider& th = ThermoProvider(), CfHookStats* stats = nullptr);

// After the last stage of a step (or before outputs): volume average-down of every transfer scalar, finest level first (the covered cells then hold the fine data).
void average_down_hierarchy(fdsamr::TimeLoop& loop);

}  // namespace fdsrt

#endif
