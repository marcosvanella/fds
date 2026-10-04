// DriverAdapter.H: glue between Role 3's level operators (LevelOps.H) and Role 1's per-level objects (Source/driver: LevelRegistry, CfGhostHook).
// Includes driver headers (read only); nothing in Source/driver is edited. Cell-centred scalar fields only (R2a): face velocities and H/HS belong to later steps
// (face-velocity transfer, pressure boundary values), they are skipped here on purpose.
#ifndef FDSRT_DRIVERADAPTER_H_
#define FDSRT_DRIVERADAPTER_H_

#include <functional>
#include <string>
#include <vector>

#include "GhostExchange.H"
#include "LevelRegistry.H"

namespace fdsrt {

// True for a field the R2a transfers handle: cell-centred and not H/HS (pressure ghosts are the pressure backend's).
bool is_transfer_scalar(const fdsamr::Fields& F, const std::string& name);

// Equation of state for the ghost cells of TMP and RSUM (FDS: GET_SPECIFIC_GAS_CONSTANT, PBAR_P). Optional: without it TMP and RSUM ghosts are injected/averaged
// (documented fallback in LevelOps.H). pbar gets the level and the (ghost) cell.
struct ThermoProvider {
    int n_tracked = 0;
    std::function<double(const double* zz)> rsum;
    std::function<double(int level, const amrex::IntVect& cell)> pbar;
};

// Totals of what the hook did (for logs and tests).
struct CfHookStats {
    long fine_ghost_cells = 0, covered_cells = 0, conflicts = 0, calls = 0;
};

// The hook to install with TimeLoop::set_cf_ghost_hook(level, hook) for every level > 0. FDS ghost rules of LevelOps.H, in this order, for the fields of the request:
//  1. volume average-down of the present fields into the covered cells of level-1 (so layers >= 3 are consistent),
//  2. covered coarse cells next to a face, layers 1 and 2, by the FDS rule (fine -> coarse), then the same-level ghost refresh of those coarse fields,
//  3. fine ghost cells, layers 1 and 2, by the FDS rule (coarse -> fine).
// Cell-centred scalars only: RHO/RHOS+ZZ/ZZS, TMP, RSUM, MU, KRES, D, DS; H/HS and the face fields are skipped. Level 0 and a missing coarser level: no-op.
fdsamr::CfGhostHook make_cf_ghost_hook(fdsamr::LevelRegistry& reg, const ThermoProvider& th = ThermoProvider(), CfHookStats* stats = nullptr);

// Fine -> coarse restriction of the given fields (all transfer scalars if names is empty) for levels top-1 .. 0, finest first, so that data propagates
// through the whole hierarchy. Covered coarse cells only. Species: RHO with ZZ (and RHOS with ZZS) are restricted by average_down_species (rho and rho*Z averaged,
// Z rebuilt as rho*Z/rho); the mass fractions are NEVER averaged linearly (that breaks the species mass, ~1e-4 relative at a flame front). Other fields: volume mean.
void average_down_registry(fdsamr::LevelRegistry& reg, const std::vector<std::string>& names = {});

}  // namespace fdsrt

#endif
