# Patch 0005 (DRAFT): `mesh.f90`, standard-conforming `POINT_TO_BOX` hook

**Status: draft for the Architect, due end of S5. Not validated with the oneAPI compilers (`ifx`/`ifort` are not available on this box). Validated only with
`git apply --check` and a gfortran 14.2 build of the patched tree (`-std=f2018`, the FDS flags), plus the unchanged-output checks below.**

Apply after 0001 and 0002. Touches only `Source/mesh.f90` (+37 lines, all inside `#ifdef WITH_AMREX`; preprocessor output with the macro undefined is unchanged).

## What it does
- `MESH_POINTERS` gets `USE FDS_AMREX_HOOKS, ONLY: BOX_VIEW` (module `Source/driver/fds_amrex_hooks.f90`, compiled only with `USE_AMREX=ON`).
- `SUBROUTINE POINT_TO_BOX(NM)`: calls `POINT_TO_MESH(NM)` first, then re-points the module pointers U,V,W,US,VS,WS,D,DS,H,HS,KRES,FVX,FVY,FVZ,RHO,RHOS,MU,TMP,Q,RSUM,ZZ,ZZS
  at `BOX_VIEW(NM)%...` when that view is associated. A routine that needs the box data calls `POINT_TO_BOX(NM)` instead of `POINT_TO_MESH(NM)`; everything else is
  pointed exactly as before.
- The views are built by `FDS_HOOK_SET_VIEW` from the C address of a FAB with `C_F_POINTER` (rank-1 flat array) and pointer bounds remapping
  `P(lb1:ub1,lb2:ub2,lb3:ub3) => FLAT` (Fortran 2008): no `MESHES(NM)%X` allocatable is involved, no compiler extension, no `LOC`/`TRANSFER`, no Cray pointer.
  The S2 alias (`fds_alias.c`, `fds_box_shim.f90`) that makes the `MESHES(NM)%X` allocatable descriptors describe FAB memory is a non-standard descriptor
  trick (W2); `POINT_TO_BOX` is its standard-conforming replacement: once the kernels call `POINT_TO_BOX` the alias can be dropped.

## Why it is a draft
- The S3/S4 kernels still run through the alias (bitwise checks pass); `POINT_TO_BOX` and `FDS_HOOK_SET_VIEW` are compiled and linked in the `USE_AMREX=ON`
  build, but no kernel calls them yet. The hook needs S5 (time loop) to be exercised with real routines; `FDS_HOOK_SET_VIEW` is not yet unit-tested.
- Not compiled with `ifx`/`ifort`; pointer bounds remapping of a rank-1 target onto rank 3/4 is standard but some oneAPI versions have had defects with it: the
  Architect should build it with the project's oneAPI toolchain before accepting it.
- The routines FDS calls through `POINT_TO_MESH` per module (e.g. `WALL_ROUTINES` local pointers `RHOP,UU,VV,WW,ZZP`) are assigned in those routines from the module
  pointers, so they see the box data automatically; module-level copies kept across calls would not, there are none in the M2a kernels.

## Evidence
- `git apply --check` passes on the tree at HEAD and in order with 0001..0004 (scratch tree `role1-s4-work/s4tree`).
- `USE_AMREX=OFF` (all of 0003, 0004, 0005 applied): `shunn3_32` 1 rank PASS (16 files), `shunn3_4mesh_32` 4 ranks PASS (47 files), bitwise to the baseline.
- `USE_AMREX=ON`, gfortran 14.2: patched tree builds, kernel checks and driver tests pass as in 0003.

## Not done here
oneAPI validation, a unit test of `FDS_HOOK_SET_VIEW`/`POINT_TO_BOX` with data from a FAB, the switch of the S3 kernels from the alias to `POINT_TO_BOX` (S5/S6).
