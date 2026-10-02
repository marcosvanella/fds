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
2. `Source/read.f90` (or the set-up entry): hand the input file name or text to the C++ driver so that `parse_amr_params` can run
   at set-up. Alternative with no Fortran change: the C++ driver reads the file itself (path is already known to FDS_SETUP).
3. Top-level `CMakeLists.txt`, inside the `USE_AMREX` block: `add_subdirectory(Source/regrid_transport)` and link
   `fds_regrid_transport` into the `fds` target (include path and static library are already provided by that directory).
4. Decision needed: when an input has finer `&MESH` entries but no `&AMR` line, is AMR mode inferred (MAX_LEVEL from the
   mesh cell sizes) or is that an error naming the missing line? The unmodified `ns2d_16_int_1to2_refinement` input has no `&AMR`.

## Tests
`ctest` in this directory: `regrid_transport_amr_input` (pure C++, no AMReX).
