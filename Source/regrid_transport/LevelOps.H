// LevelOps.H: level-to-level operators on cell-centred MultiFabs (R2a). AMReX only, no FDS and no driver types.
//
//  average_down_cells : fine -> coarse volume average on the coarse cells covered by the fine level (conservative: the integral over the
//                       covered region is unchanged up to round-off).
//  fill_cf_ghosts_pc  : ghost cells of the fine level that belong to no fine grid (the cells across a coarse-fine face) get the value of the
//                       coarse cell that contains them (piecewise-constant, no slope). With ratio >= 2 and two layers, layer 2 is in the same
//                       coarse cell as layer 1, so layer 2 = layer 1, as the FDS scheme of the numeric-scheme note prescribes (Phase 3).
//                       Ghost cells filled by a neighbouring fine grid (same-level exchange) and ghost cells outside a non-periodic domain
//                       edge are left untouched.
#ifndef FDSRT_LEVELOPS_H_
#define FDSRT_LEVELOPS_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace fdsrt {

// Components scomp .. scomp+ncomp-1 of `fine` are averaged into the same components of `crse` on every coarse cell that lies under a fine grid.
// Other coarse cells are not touched. Valid cells only.
void average_down_cells(const amrex::MultiFab& fine, amrex::MultiFab& crse, const amrex::Geometry& fine_geom, const amrex::Geometry& crse_geom,
                        const amrex::IntVect& ratio, int scomp, int ncomp);

// nlayers (<= fine.nGrow(), <= ratio) ghost layers of `fine`, components scomp .. scomp+ncomp-1, from `crse`. Run the same-level ghost fill first.
// Returns the number of ghost cells written on this rank.
long fill_cf_ghosts_pc(amrex::MultiFab& fine, const amrex::MultiFab& crse, const amrex::Geometry& fine_geom, const amrex::Geometry& crse_geom,
                       const amrex::IntVect& ratio, int nlayers, int scomp, int ncomp);

}  // namespace fdsrt

#endif
