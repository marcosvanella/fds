// CellTransfer.H: prolongation of conserved cell data to a new fine patch by injection (owner directive Oct 14).
// AMReX only, host loops (regrid-time work may stay on the host until Phase 11, NFR-043), C++ glue (D-043).
//
// A *parent* is the coarse cell and its *children* are the r_x*r_y*r_z fine cells inside it. prolong_conserved gives every child the parent's value in
// every component (piecewise constant: no slope, no limiter, no floor, no renormalisation). Each component is conserved per parent exactly, and since
// rho, rho*Z_n and T are all injected the children keep the parent's equation-of-state state and mass fractions. (Until Oct 14: MC-limited linear
// slopes per component with a positivity clip; a slope that differs per species moved rho*Z between species against the injected T.)
#ifndef FDSRT_CELLTRANSFER_H_
#define FDSRT_CELLTRANSFER_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <vector>

namespace fdsrt {

struct ProlongOpts {           // kept for the callers; injection uses neither field
    bool use_floor = true;
    double floor = 0.0;
};

struct ProlongStats {
    long parents = 0;          // coarse parents processed (all ranks)
    long limited = 0;          // always 0 (injection; kept for the callers' statistics)
    long clips = 0;            // always 0
    long unfixable = 0;        // always 0
};

// Fills `fine` (ncomp components, valid cells) from the valid cells of `crse` (same ncomp).
// Every fine box must be a union of whole parents (its edges multiples of `ratio`). If `old_fine` is given, its valid cells are then copied over the
// overlapping part of `fine` (bitwise), so only the cells that had no fine data are injected.
ProlongStats prolong_conserved(amrex::MultiFab& fine, const amrex::MultiFab& crse, const amrex::Geometry& crse_geom, const amrex::IntVect& ratio,
                               const ProlongOpts& opt, const amrex::MultiFab* old_fine = nullptr);

}  // namespace fdsrt

#endif
