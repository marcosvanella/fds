// GhostShare.H: D-059 diagnostic. A covered coarse cell that is a ghost cell of two or more faces of the fine patches (a patch corner or edge, or a patch thinner than
// four coarse cells) holds ONE value in AMReX where FDS has one array per mesh; the ruling accepts this as a documented limitation on condition that the driver counts
// those cells and prints the count at each regrid. This is the count, from the grids alone (no field data, no MPI communication; every rank gets the same numbers).
// "Serves a face" follows LevelOps fill_covered_ghosts_fds exactly: layer 1 = covered cell with an uncovered coarse neighbour across the face; layer 2 = covered cell
// whose neighbour towards the face is a layer-1 cell (covered) and whose next neighbour is uncovered. Periodic wrap respected; a neighbour outside a non-periodic edge
// is not an interface.
#ifndef FDSRT_GHOSTSHARE_H_
#define FDSRT_GHOSTSHARE_H_

#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>

#include <string>

namespace fdsrt {

struct GhostShareCount {
    long serving_one = 0;         // covered coarse cells that are the ghost cell of exactly one face (and one layer)
    long serving_several = 0;     // covered coarse cells that serve two or more (face, layer) candidates: the limitation of D-059
    long serving_two_faces = 0;   // of those, cells whose candidates lie on two or more different faces (corner/edge cells); the rest serve one face in layers 1 and 2
    std::string to_string() const;
};

// fine_ba: grids of the fine level; ratio: refinement ratio to the coarse level whose cells are counted; crse_geom: geometry of the coarse level.
GhostShareCount count_shared_ghost_cells(const amrex::BoxArray& fine_ba, const amrex::IntVect& ratio, const amrex::Geometry& crse_geom);

}  // namespace fdsrt

#endif
