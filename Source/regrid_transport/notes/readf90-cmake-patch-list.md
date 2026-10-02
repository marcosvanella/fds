# Patch list for the Chief Architect: `&AMR` input and CMake (WITH_AMREX only)

Status: R0 proposal. Nothing here is applied; Role 3 does not edit `read.f90` or the top-level CMake.

## Current state
`AmrInput.H/.cpp` parse `&AMR` and `&AMR_REGION` from the input text on their own (`parse_amr_params`). They need no change to
FDS to be tested. Names are working names (IR-003). The parser rejects unknown names, REF_RATIO other than 2 or 4, bad list
lengths, non-power-of-2 blocking factors, blocking factors that grow faster than the ratio, and out-of-range values.

## Needed in existing files
1. `Source/read.f90`, behind `#ifdef WITH_AMREX`: FDS stops on an unknown namelist group, so `&AMR` and `&AMR_REGION` must be
   skipped by the FDS namelist reader (as other groups it does not read), otherwise every AMR input is rejected before the C++
   parser runs. Exact place: where the reader dispatches on the group name. The C++ side then reads the same input file text.
2. No Fortran change (ruling of 2026-10-02): the C++ driver reads the input file itself and calls `parse_amr_params`, then
   `build_hierarchy_from_meshes` with the mesh list from the FDS set-up (`MeshInput` = ijk, xb, rank per mesh).
3. Top-level `CMakeLists.txt`, inside the `USE_AMREX` block: `add_subdirectory(Source/regrid_transport)` and link
   `fds_regrid_transport` into the `fds` target (include path and static library are already provided by that directory).
4. Decided (ruling of 2026-10-02): finer `&MESH` entries without an `&AMR` line are an error naming the missing line and a mesh
   pair; AMR mode is never inferred. Implemented in `group_meshes`. Consequence: the unmodified `ns2d_16_int_1to2_refinement`
   input needs an `&AMR MAX_LEVEL=1 /` line (the test case in `tests/cases/` has it).

## Tests
`ctest` in this directory: `regrid_transport_amr_input`, `regrid_transport_hierarchy` (pure C++, no AMReX).
