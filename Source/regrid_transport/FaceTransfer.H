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
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <array>

namespace fdsrt {

// crse[d]: coarse faces of direction d covering (for the normal direction) the coarse faces of the cells that contain `fine_cells`;
// fine[d]: fine faces of direction d over `fine_cells` (cell-centred box of the fine patch), which are written (all of them, d = 0..2).
// `fine_cells` must be a multiple of `ratio` on both ends (a union of whole coarse cells).
void prolong_faces_normal_linear(const std::array<const amrex::FArrayBox*, 3>& crse, const std::array<amrex::FArrayBox*, 3>& fine,
                                 const amrex::Box& fine_cells, const amrex::IntVect& ratio);

// ---- D-060 level path (the single face-velocity prolongation of a regrid) ---------------------------------------------------------------------------
// Face-centred MultiFabs of a level: fine[d] / crse[d] / old_fine[d] are nodal in direction d, one component, any ghost width (only valid faces are read and written).
// prolong_faces_level fills the VALID faces of fine[0..2] (all of them) from the coarse level by prolong_faces_normal_linear, in three steps:
//   1. the coarse faces under each fine box are gathered into temporary coarse-aligned arrays (ParallelCopy with the periodicity of the coarse level);
//   2. prolong_faces_normal_linear; fine faces on a coarse face carry exactly the coarse value (the interface faces of a patch included);
//   3. if old_fine is given, the valid faces of the old fine level are copied over the overlapping faces (bitwise), so only faces without fine data are new.
// Every fine box must be a union of whole coarse cells (blocking factor >= ratio). The fine MultiFab arrays must already exist on the new BoxArray.
void prolong_faces_level(const std::array<amrex::MultiFab*, 3>& fine, const std::array<const amrex::MultiFab*, 3>& crse, const amrex::Geometry& crse_geom,
                         const amrex::IntVect& ratio, const std::array<const amrex::MultiFab*, 3>& old_fine = {nullptr, nullptr, nullptr});

// Fine -> coarse mean of the faces: every coarse face that carries fine faces gets the mean of the (r_t1 * r_t2) fine faces that lie on it (area average; the
// normal-velocity flux through the coarse face is then the sum of the fine ones). Only coarse faces that coincide with fine faces are written (interface included).
void average_down_faces(const std::array<const amrex::MultiFab*, 3>& fine, const std::array<amrex::MultiFab*, 3>& crse, const amrex::IntVect& ratio);

// max over the valid cells of `vel` (three face MultiFabs, same cell BoxArray) of |div u - D|, div u = sum_d (U_d(hi) - U_d(lo)) / dx_d; D = nullptr: target 0.
// Local to this rank, the caller reduces. Directions with one cell are included like any other (their faces are real faces).
double max_divergence_error_local(const std::array<const amrex::MultiFab*, 3>& vel, const amrex::Geometry& geom, const amrex::MultiFab* D = nullptr);

}  // namespace fdsrt

#endif
