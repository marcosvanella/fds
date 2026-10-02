// RegridInterface.H: PROPOSAL (R0) for the interface between regrid_transport (Role 3) and the driver registry (Role 1).
// Declarations only; for review by Role 1 and the Integration Lead. Nothing here is implemented or linked yet.
//
// Division of work (ruling 1): Role 3 owns the hierarchy (an AmrCore subclass: the AMReX class that holds the grids of every level
// and the callbacks "make / remake / clear level") and the level operators. Role 1 owns the field registry and the kernel shim and
// generalises Fields/TimeLoop to a level index. The two sides meet in three places:
//   1. layouts: Role 3 tells the registry when a level is created, changed or removed (LevelListener);
//   2. fluxes: the registry exposes the stage face fluxes of each level and accepts an override list (see notes/flux-override-interface.md);
//   3. the time step: one global dt, the minimum over all levels and ranks, the same in both stages (ADR-002 v1.2).
#ifndef FDSRT_REGRID_INTERFACE_H_
#define FDSRT_REGRID_INTERFACE_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

#include <array>
#include <vector>

namespace fdsrt {

// Everything the registry needs to build the fields of one level. Boxes are in the index space of that level; the level-0 boxes are the
// level-0 &MESH boxes in input order followed by the boxes under finer meshes, so box i < n_level0_meshes is FDS mesh i+1.
struct LevelLayout {
    int level = 0;
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    amrex::IntVect ref_ratio_from_parent{1, 1, 1};  // 1 in a single-cell (hidden) direction
};

// Role 1 implements this (one object for the whole registry). Called from Role 3's AmrCore callbacks, between time steps only.
class LevelListener {
public:
    virtual ~LevelListener() = default;
    virtual void make_level(const LevelLayout& l) = 0;    // new level: allocate Fields/SideData; data is filled by Role 3 afterwards
    virtual void remake_level(const LevelLayout& l) = 0;  // same level, new grids: old fields stay readable until Role 3 has copied them
    virtual void clear_level(int level) = 0;
};

enum class FluxKind { Adv, Dif, Heat };

// One override list per box and direction (shapes and units: notes/flux-override-interface.md). Empty list = no change, bitwise.
struct FluxOverride {
    int dir = 0;                              // normal direction 0,1,2
    int nscal = 0;                            // values per face
    std::vector<std::array<int, 3>> face;     // AMReX face index (low face of cell a), sorted by (k,j,i)
    std::vector<double> value;                // face.size()*nscal, face-major
};

// Role 1 implements this for the kernel loop.
class FluxAccess {
public:
    virtual ~FluxAccess() = default;
    // Stage face fluxes of a level: nodal in `dir`, component n-1 = scalar n, no ghost layer, same BoxArray/DistributionMapping as the level.
    virtual const amrex::MultiFab& stage_flux(int level, FluxKind kind, int dir) const = 0;
    // Replaces the override lists for the next stage on `level` (all boxes owned by this rank); an empty vector clears them.
    virtual void set_flux_override(int level, FluxKind kind, const std::vector<std::vector<FluxOverride>>& per_local_box) = 0;
    // Smallest stable dt on this rank for `level`; Role 3 / TimeLoop reduce it over levels and ranks.
    virtual double local_dt(int level) const = 0;
};

}  // namespace fdsrt

#endif
