# Patch 0002: `USE_AMREX` option in the top-level CMakeLists.txt

Apply after 0001. Touches only `CMakeLists.txt` (+8 lines): the option `USE_AMREX` (default OFF) and, before `install(TARGETS fds)`,
`if(USE_AMREX) include(Source/driver/CMakeLists.txt) endif()`. With the option OFF nothing else changes: the added option is unused, no
language is enabled, no target is touched (CMake configure output differs only by the extra cache entry).

With `USE_AMREX=ON`, `Source/driver/CMakeLists.txt` turns the `fds` target into `fds_amr`: Fortran sources compiled with `-DWITH_AMREX`,
C++ `main.cpp` and `FdsAmr.cpp`, `fds_mesh_query.f90` added, linked with C++ against AMReX (`AMREX_INSTALL_DIR`, default
`/workspace/amrex-install`) and OpenMP. The Fortran-only compile options and macros of the `fds` target are restricted to the Fortran
compiler with generator expressions.

Configure example (out of tree): `cmake -S <src> -B <build> $FDS_CMAKE_COMMON -DUSE_AMREX=ON` (see `../README.md`).

## Evidence
- `git apply --check` passes on the current tree and on a scratch copy with 0001 applied.
- USE_AMREX=ON build of the scratch copy with 0001+0002 and the driver files: builds, no warnings from the driver files
  (`-Wall -Wextra`), `tests/check_setup_amr.sh`: `shunn3_32` 1 rank -> 1 box, `shunn3_4mesh_32` 4 ranks -> 4 boxes (see `../README.md`).
- USE_AMREX=OFF with 0001+0002 both applied (`tests/check_off_bitwise.sh`): `shunn3_32` 1 rank PASS (16 files), `shunn3_4mesh_32`
  4 ranks PASS (47 files), all bitwise identical to the baseline. Further evidence for the source change itself is in `0001-setup-callable.md`.
