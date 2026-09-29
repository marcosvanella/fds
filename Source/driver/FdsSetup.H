// FdsSetup.H: level-0 layout and cell wall data from the FDS set-up (needs the FDS objects; the unit tests do not use this file).
// Kernel-facing rules (M2a): (a) passive scalars are handled by Fields.cpp; (b) only uniform Cartesian metrics are used.
#ifndef FDS_SETUP_H_
#define FDS_SETUP_H_

#include <AMReX_Array4.H>
#include <AMReX_Box.H>

#include "FdsAmr.H"

namespace fdsamr {

// Queries the Fortran set-up data (fds_mesh_query.f90) and builds Level0 (assemble_level0).
Level0 build_level0();

// Cell wall provider for SideData (SideData.H) reading MESHES(nm)%CELL: comp 0 = SOLID, comps 1..6 = WALL_INDEX(-1,1,-2,2,-3,3).
void fds_cell_walls(int nm, const amrex::Box& vbox, amrex::Array4<int> const& a);

}  // namespace fdsamr

#endif
