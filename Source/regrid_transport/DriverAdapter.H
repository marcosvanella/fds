// DriverAdapter.H: glue between Role 3's level operators (LevelOps.H) and Role 1's per-level objects (Source/driver: LevelRegistry, CfGhostHook).
// Includes driver headers (read only); nothing in Source/driver is edited. Cell-centred scalar fields only (R2a): face velocities and H/HS belong to later steps
// (face-velocity transfer, pressure boundary values), they are skipped here on purpose.
#ifndef FDSRT_DRIVERADAPTER_H_
#define FDSRT_DRIVERADAPTER_H_

#include <string>
#include <vector>

#include "GhostExchange.H"
#include "LevelRegistry.H"

namespace fdsrt {

// True for a field the R2a transfers handle: cell-centred and not H/HS (pressure ghosts are the pressure backend's).
bool is_transfer_scalar(const fdsamr::Fields& F, const std::string& name);

// The hook to install with TimeLoop::set_cf_ghost_hook(level, hook) for every level > 0. For each field of the request: fill nlayers (default 2) ghost layers
// of that field on the level from level-1, piecewise-constant (fill_cf_ghosts_pc). Level 0 and a missing coarser level: no-op.
fdsamr::CfGhostHook make_cf_ghost_hook(fdsamr::LevelRegistry& reg, int nlayers = 2);

// Fine -> coarse volume average of the given fields (all transfer scalars if names is empty) for levels top-1 .. 0, finest first, so that data propagates
// through the whole hierarchy. Covered coarse cells only.
void average_down_registry(fdsamr::LevelRegistry& reg, const std::vector<std::string>& names = {});

}  // namespace fdsrt

#endif
