// CompositeOps.H: the AMReX-only pieces of the two-level (composite) step that are not FDS kernels (S14). Free functions on a LevelRegistry, unit-testable without the FDS objects.
//
// Plain words. A fine level sits on top of a coarse one. Where they meet (the coarse-fine interface) the step needs five things the single-level code never had:
//   * the coarse faces under a fine face carry the area mean of the fine faces (average_down), for the velocities and for the momentum-equation forcing FVX, FVY, FVZ;
//   * the fine ghost faces of the velocity are interpolated from the coarse level (fill_fine_face_ghosts);
//   * the fine level's Poisson right-hand side is built in C++ from the same terms FDS uses (fine_pressure_rhs);
//   * the pressure H in the fine ghost cells is set so that the velocity-corrector kernel, which differences H across the face, reproduces the composite face gradient (set_fine_h_ghosts).
// Face fields are nodal MultiFabs of the Fields; AMReX face a is the low face of cell a (FDS face I = a).
#ifndef FDS_COMPOSITEOPS_H_
#define FDS_COMPOSITEOPS_H_

#include <AMReX_MultiFab.H>

#include <array>

#include "LevelRegistry.H"

namespace fdsamr {

// Fine -> coarse area mean of the faces named in `names` (U,V,W or US,VS,WS or FVX,FVY,FVZ), finest level first, for levels 0..nl-1 of the registry.
void composite_average_down_faces(LevelRegistry& reg, int nl, const std::array<const char*, 3>& names);

// Ghost faces of level l >= 1 for the three face fields `names`: faces outside the fine boxes that no fine box owns get the coarse-interpolated value (normal-linear, constant
// across the coarse face tangentially: fdsrt::prolong_faces_level), then the same-level fill overwrites those that a neighbouring fine box owns. The level's valid faces are untouched.
// The same for FVX, FVY, FVZ, which are cell-shaped arrays holding the flux at the upper face of the cell (FDS index): the nodal AMReX face n is FV(n-1).
void composite_average_down_flux(LevelRegistry& reg, int nl);

void composite_fill_fine_face_ghosts(LevelRegistry& reg, int l, const std::array<const char*, 3>& names);

// PRHS of a fine level, the interior formula of PRESSURE_SOLVER_COMPUTE_RHS (a fine box has no wall cell): (FVX(I-1)-FVX(I))/dx + (FVY(J-1)-FVY(J))/dy + (FVZ(K-1)-FVZ(K))/dz - DDDT.
void composite_fine_pressure_rhs(Fields& F, const Level& lev, amrex::MultiFab& rhs);

// H (or HS) ghost cells of level l >= 1 next to the coarse level: H(ghost) = H(inside) -/+ g * dx so that the kernel difference (H(I+1)-H(I))/dx on the interface face equals the
// composite face gradient g = grad[d] (nodal MultiFab on the level's BoxArray). Ghost cells that another fine box owns (same level, periodic images included) are not touched.
// Returns the number of ghost cells set on this rank.
long composite_set_fine_h_ghosts(const Level& lev, amrex::MultiFab& H, const std::array<amrex::MultiFab*, 3>& grad);

}  // namespace fdsamr

#endif
