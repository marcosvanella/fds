// TimeLoopWiring.cpp: see TimeLoopWiring.H.
#include "TimeLoopWiring.H"

namespace fdsrt {

void install_cf_ghost_hooks(fdsamr::TimeLoop& loop, int finest_level, const ThermoProvider& th, CfHookStats* stats)
{
    for (int l = 1; l <= finest_level; ++l) loop.set_cf_ghost_hook(l, make_cf_ghost_hook(loop.registry(), th, stats));
}

void average_down_hierarchy(fdsamr::TimeLoop& loop) { average_down_registry(loop.registry()); }

}  // namespace fdsrt
