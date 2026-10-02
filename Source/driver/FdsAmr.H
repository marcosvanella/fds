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

// One AMR level (S9, per-level interface). Level 0 is one box per FDS &MESH, box i = FDS mesh i+1, owner rank = FDS PROCESS(NM), no re-boxing. For a level > 0
// the boxes come from the regrid (fdsrt::LevelLayout) and `mesh` stays empty until fine-level FDS MESH_TYPE objects exist (notes/level-interface.md section 3):
// box i of level l is FDS mesh `fds_mesh_offset + i + 1` once it is bound. Level0 is an alias, so existing code and the single-level path are unchanged.
struct Level {
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    amrex::Vector<MeshInfo> mesh;
    DomainInfo dom;
    double dx[3];
    int level = 0;                                  // AMR level index (0 = coarsest)
    amrex::IntVect ref_ratio_from_parent{1, 1, 1};  // refinement ratio to level-1 (1 for level 0)
    int fds_mesh_offset = 0;                        // FDS mesh number of box i is fds_mesh_offset + i + 1 (level 0: 0). Meaningful only when mesh is bound
    bool fds_bound() const { return !mesh.empty(); }   // true when every box has an FDS MESH_TYPE behind it (always at level 0)
};
using Level0 = Level;

// Builds Level0 from the mesh list (index i = FDS mesh i+1). Aborts with a message for meshes outside the M2a scope (cylindrical,
// unequal cell sizes, corners off the cell lattice, overlapping boxes, rank out of range). Pure C++: FdsSetup.cpp fills the mesh list
// from the FDS set-up, the unit tests fill it synthetically.
Level0 assemble_level0(const amrex::Vector<MeshInfo>& meshes, const DomainInfo& dom);

// Layout-only level (level > 0 or a test hierarchy): geometry, boxes, ranks and cell sizes, no FDS mesh binding (mesh stays empty). `ba`/`dm`/`geom` come from the
// regrid (fdsrt::LevelLayout). dx is taken from geom.CellSize(). Pure C++, no FDS symbols.
Level make_layout_level(int level, const amrex::Geometry& geom, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                        const amrex::IntVect& ref_ratio_from_parent, const DomainInfo& dom);

// Prints one line per box (rank 0 only).
void print_level0(const Level0& l0);

}  // namespace fdsamr

#endif
