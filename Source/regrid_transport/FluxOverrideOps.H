// FluxOverrideOps.H: the interface flux overwrite (D-050), Role 3 side: from the stage face fluxes of a fine level to the override lists of the
// coarse level below it (FluxOverride of RegridInterface.H, shapes in notes/flux-override-interface.md). AMReX only, no driver types, no hooks needed:
// the fine fluxes come in as MultiFabs (what FluxAccess::stage_flux returns), the result is what FluxAccess::set_flux_override takes.
//
// A coarse face is an interface face when exactly one of the two coarse cells it separates is covered by the fine level (periodic wrap respected; a face
// whose other side is outside a non-periodic domain edge is not an interface). Its value is the area-weighted sum of the fine faces that tile it:
//   F_coarse(n) = sum over fine faces f of (A_f / A_c) F_f(n)       (uniform cells: the plain mean of the r1*r2 fine face values; ARO of wall.f90)
// summed in the order of FDS (k outer, j, i inner). A face shared by several coarse boxes is listed in every box that holds it (Role 1 requirement).
#ifndef FDSRT_FLUXOVERRIDEOPS_H_
#define FDSRT_FLUXOVERRIDEOPS_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <vector>

#include "RegridInterface.H"

namespace fdsrt {

struct OverrideStats {
    long entries = 0;   // face entries on this rank, summed over local boxes (a shared face counts once per box)
};

// fine_flux[d]: the stage flux of the fine level, nodal in direction d, nscal components, no ghost layer, on (fine_ba, fine_dm) (FluxAccess::stage_flux).
// Returns one vector per LOCAL coarse box (MFIter order over crse_ba/crse_dm), each with one FluxOverride per direction that has at least one face; boxes without
// an interface face give an empty vector. An empty fine_ba gives all-empty vectors (the empty set: the kernels must then be bitwise unchanged).
std::vector<std::vector<FluxOverride>> build_flux_overrides(const amrex::BoxArray& crse_ba, const amrex::DistributionMapping& crse_dm, const amrex::Geometry& crse_geom,
                                                           const amrex::BoxArray& fine_ba, const amrex::DistributionMapping& fine_dm, const amrex::IntVect& ratio,
                                                           const amrex::MultiFab* const fine_flux[3], OverrideStats* stats = nullptr);

}  // namespace fdsrt

#endif
