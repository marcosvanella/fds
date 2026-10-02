// CellTransfer.H: prolongation of conserved cell data to a new fine patch, with a positivity-preserving limiter and a parent-sum-preserving clip (R4).
// AMReX only, host loops (regrid-time work may stay on the host until Phase 11, NFR-043), C++ glue (D-043).
//
// Terms. A *parent* is the coarse cell and its *children* are the r_x*r_y*r_z fine cells inside it. The data are conserved densities, e.g. rho*Z_n of every
// tracked species n (mass per volume); rho itself is rebuilt by the caller as the sum over species, so the mass fractions Y_n = rho*Z_n / rho lie in [0,1]
// and sum to 1 as soon as every rho*Z_n >= 0.
//
// prolong_conserved (per component, per parent):
//   1. slope in each direction with r_d > 1: monotonised-central limited difference of the coarse neighbours (zero at a local extremum);
//   2. child = q_parent + sum_d slope_d * x_d, with x_d the child centre offset in coarse-cell units, symmetric about the parent centre, so the children
//      average to the parent value (conservation to round-off);
//   3. all slopes of the parent are scaled by one factor alpha <= 1 so that every child stays inside [max(min of the 3x3x3 neighbourhood, floor),
//      max of the neighbourhood] (no new extremes, no negative value from a non-negative parent). With data >= floor no clip ever fires;
//   4. clip (only if a child is still below `floor`, e.g. a parent that is negative itself or round-off): children are raised to the floor and the others
//      scaled so that the children's sum equals r^3 * q_parent again (parent sum kept); a parent below the floor gets the piecewise-constant children
//      q_parent (sum kept, value still below the floor, counted as `unfixable`). Every clip is counted (FR-025).
#ifndef FDSRT_CELLTRANSFER_H_
#define FDSRT_CELLTRANSFER_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <vector>

namespace fdsrt {

struct ProlongOpts {
    bool use_floor = true;     // enforce value >= floor (species densities). false: no floor (no clip), only the extremum bounds of step 3
    double floor = 0.0;
};

struct ProlongStats {
    long parents = 0;          // coarse parents processed (all ranks)
    long limited = 0;          // parents whose slopes were scaled by alpha < 1 or set to zero by the limiter in at least one component
    long clips = 0;            // parents where step 4 changed at least one child (component-parent pairs)
    long unfixable = 0;        // of those, parents that are below the floor themselves
};

// Step 4 alone, on one parent and one component: `children` (count n) are modified in place; returns 0 = nothing to do, 1 = clipped and rescaled,
// 2 = parent below the floor (children set to the parent value). Exposed for the unit test.
int clip_children(double* children, int n, double parent, double floor);

// Fills `fine` (ncomp components, valid cells) from `crse` (same ncomp, one valid ghost layer: same-level exchange and physical boundary already done).
// Every fine box must be a union of whole parents (its edges multiples of `ratio`). If `old_fine` is given, its valid cells are then copied over the
// overlapping part of `fine` (bitwise), so only the cells that had no fine data are interpolated.
ProlongStats prolong_conserved(amrex::MultiFab& fine, const amrex::MultiFab& crse, const amrex::Geometry& crse_geom, const amrex::IntVect& ratio,
                               const ProlongOpts& opt, const amrex::MultiFab* old_fine = nullptr);

}  // namespace fdsrt

#endif
