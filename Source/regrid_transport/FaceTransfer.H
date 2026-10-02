// FaceTransfer.H: prolongation of face-centred velocity for a new fine patch (R4). AMReX only.
//
// prolong_faces_normal_linear: the fallback for the cases the AMReX interpolator FaceDivFree cannot do (ratio 1 in a single-cell direction,
// see notes/r4-facedivfree-findings.md). A fine face takes the coarse face value of the coarse cell that contains it in the two transverse
// directions, and is interpolated linearly between the two coarse faces of that cell in the normal direction. Therefore
//   * a fine face that lies on a coarse face has exactly the coarse value (the interface faces of the design note);
//   * the mean of the fine faces on a coarse face is the coarse value (average-down of faces gives the coarse field back);
//   * every fine cell has exactly the divergence of its coarse parent cell (each direction contributes (U_hi - U_lo)/dx_coarse),
//     so a divergence-free coarse field gives a divergence-free fine field.
// The price against FaceDivFree is a velocity that is constant across a coarse face in the transverse directions (less smooth).
// Written with amrex::ParallelFor (glue code, not a physics kernel; D-043).
#ifndef FDSRT_FACETRANSFER_H_
#define FDSRT_FACETRANSFER_H_

#include <AMReX_FArrayBox.H>

#include <array>

namespace fdsrt {

// crse[d]: coarse faces of direction d covering (for the normal direction) the coarse faces of the cells that contain `fine_cells`;
// fine[d]: fine faces of direction d over `fine_cells` (cell-centred box of the fine patch), which are written (all of them, d = 0..2).
// `fine_cells` must be a multiple of `ratio` on both ends (a union of whole coarse cells).
void prolong_faces_normal_linear(const std::array<const amrex::FArrayBox*, 3>& crse, const std::array<amrex::FArrayBox*, 3>& fine,
                                 const amrex::Box& fine_cells, const amrex::IntVect& ratio);

}  // namespace fdsrt

#endif
