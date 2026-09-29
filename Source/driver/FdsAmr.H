// FdsAmr.H: level-0 AMReX layout built from the FDS mesh set-up data (M2a, S1).
//
// Kernel-facing rules (M2a): (a) passive scalars (N_TOTAL_SCALARS beyond the tracked species) are handled by Fields.cpp (S2),
// not here; this file carries only the two counts. (b) Only uniform Cartesian metrics are used; R(I)/RRN(I) are dropped,
// and CYLINDRICAL or TRN* meshes are rejected (IR-002).
#ifndef FDS_AMR_H_
#define FDS_AMR_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>

namespace fdsamr {

// One entry per FDS mesh NM = 1..nmeshes (index nm-1 in the vectors).
struct MeshInfo {
    int ijk[3];
    double xb[6];
    int rank;
};

struct DomainInfo {
    int periodic[3];
    int cylindrical;
    int n_tracked;
    int n_total;
    int nranks;
};

// Level 0: one box per FDS &MESH, box i = FDS mesh i+1, owner rank = FDS PROCESS(NM). No re-boxing.
struct Level0 {
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    amrex::Vector<MeshInfo> mesh;
    DomainInfo dom;
    double dx[3];
};

// Queries the Fortran set-up data (fds_mesh_query.f90) and builds Level0. Aborts with a message for meshes outside the M2a scope
// (nonuniform grid, cylindrical, unequal cell sizes, overlapping boxes).
Level0 build_level0();

// Prints one line per box (rank 0 only).
void print_level0(const Level0& l0);

}  // namespace fdsamr

#endif
