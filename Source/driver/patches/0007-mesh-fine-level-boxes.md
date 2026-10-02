# Patch 0007 (DRAFT): `mesh.f90`, fine-level mesh objects for `POINT_TO_BOX` (D-056, option B)

**Status: DRAFT (applied to HEAD as a draft, commit 10c2a3b7b0), not validated with the oneAPI compilers (`ifx`/`ifort` are not available on this box), NOT used for physics.** Condition of D-056: patch 0005 must pass oneAPI validation
first (the Architect is arranging it); this patch grows 0005 and has the same validation requirement plus the points under "What oneAPI has to confirm". If 0005 or 0007 fails, fix it or
carry fine boxes on the alias route (`fds_shim_bind`); option A (a longer `MESHES`) is not an alternative without asking. 0005 itself is not edited: 0007 applies on top of it (the file
index line `172704f` of the tree with 0005 applied is the base).

Apply after 0005 (and 0003/0004/0006, independent). Touches only `Source/mesh.f90`, every added or changed line inside `#ifdef WITH_AMREX`; with the macro undefined the
preprocessed source is the unchanged FDS source (`USE_AMREX=OFF` bitwise check: below).

## What it does
- `MESH_POINTERS`: `USE FDS_AMREX_HOOKS, ONLY: BOX_VIEW,BOX_VIEW_TYPE,FDS_HOOK_FINE_ABORT`; one new module array
  `TYPE(FINE_LEVEL_TYPE), ALLOCATABLE, TARGET, SAVE :: FINE_LEVEL(:)`, **indexed by refinement level** (`FINE_LEVEL(1)` = level 1), with `NM0`, `N_BOXES`, `BOX(:)` (the `MESH_TYPE` objects of
  the boxes of the level) and `VIEW(:)` (`BOX_VIEW_TYPE`, the pointer views of the AMReX data of the box, as `BOX_VIEW(NM)` for level 0). The driver allocates and fills it (draft driver
  file `Source/driver/draft/fds_fine_mesh_b.f90`); no FDS routine does. A box of level L has the FDS mesh number `NM = FINE_LEVEL(L)%NM0 + IB`, always `> SIZE(MESHES)`
  (`Level::fds_mesh_offset` of the driver).
- `MESHES`, `PROCESS`, `OMESH` allocation, MPI maps and the level-0 mesh-count loops are not touched.
- `POINT_TO_MESH(NM)` is split in two: the unchanged statements `U=>M%U ... ` become `POINT_TO_MESH_OBJECT(M)` (a `TYPE(MESH_TYPE), POINTER, INTENT(IN)` dummy), `POINT_TO_MESH(NM)` is
  `M=>MESHES(NM); CALL POINT_TO_MESH_OBJECT(M)`. This is a cut of the routine in two places and a rename of its `END` (under the macro), no statement of the body changes.
- `POINT_TO_MESH(NM)` with `NM > SIZE(MESHES)` (or `NM < 1`) **stops the run with a clear message** through `FDS_HOOK_FINE_ABORT` (hook module: message on unit 0 naming the routine and the
  number, `ERROR STOP 1`): a routine that has not been made fine-ready cannot silently read `MESHES(NM)`.
- `POINT_TO_BOX(NM)`: for `NM <= SIZE(MESHES)` as in 0005; for a larger number it finds the level whose range contains `NM`, calls `POINT_TO_MESH_OBJECT(FINE_LEVEL(L)%BOX(IB))` and then re-points
  the listed arrays at `FINE_LEVEL(L)%VIEW(IB)` exactly as 0005 does at `BOX_VIEW(NM)`; a number that belongs to no level stops the run with the same message.

## Kernel wrappers and direct `MESHES(NM)` reads (driver side, in this commit)
`fds_k_*` wrappers (`fds_kernels.f90`, 16 entries with a mesh number) call `FDS_HOOK_FINE_GUARD(name, NM, NMESHES)` first: a fine number aborts with the message unless
`fds_hook_set_fine_ready(1)` was called (not called anywhere: **fine boxes are not usable for physics**). This covers every route the C++ driver has into the kernels, including the direct
`MESHES(NM)%...` reads that exist in kernel-side code, which a bounds-free build would read out of range without notice: `mass.f90` (2, `M_DOT_PPP`, same in the generated
`fds_density_split.f90`), `velo.f90` (2), `wall.f90` (2), `turb.f90` (9), the `OMESH` accesses (`velo.f90` 24 and `wall.f90` 9: `EXTERNAL_GHOSTS_FILLED` skips them at level interfaces, D-055),
`CHANGE_TIME_STEP_INDEX(NM)`/`DT_NEW_K(NM)` (sized `NMESHES`; the driver keeps per-box dt vectors by level offset). A build with `-fcheck=bounds` also stops at each of them. Before a wrapper is
made fine-ready, each of these reads is converted to the module pointer `M`/`POINT_TO_BOX` data (a small guarded kernel patch, numbered then).

## Contents of a fine box (not in this patch: the driver draft `draft/fds_fine_box_b.f90` builds them; list of `level-interface.md` section 3)
Metrics (`X,Y,Z,XC..,DX..,RDX..,DXN..,RDXN..,R=RRN=1`, `IBAR..`), `CELL`/`CELL_INDEX` (all gas), `WALL`/`EXTERNAL_WALL`/`WALL_INDEX` with the domain-edge walls and the
interface walls without `OMESH` neighbour (`INTERPOLATED_BOUNDARY`, `NIC=1`), box state arrays through `VIEW`. No obstructions at fine levels for now. Only the metrics and the pointer routing exist in the
draft driver file so far.

## Evidence (gfortran 14.2, Open MPI 5.0.7; scratch tree = HEAD tree + 0007, `-DFDS_AMR_FINE_B_DRAFT=ON`)
See "Results" below.

## What oneAPI has to confirm (for the Architect)
1. `POINT_TO_MESH_OBJECT(M)` with a `POINTER, INTENT(IN)` dummy and the assignments `U=>M%U` of allocatable components (the same statements as before, now through a pointer dummy rather than a local pointer).
2. `ALLOCATABLE` components of a derived type that are `TARGET` through the `TARGET` module variable `FINE_LEVEL` (`FINE_LEVEL(L)%BOX(IB)%U` as a pointer target; `ASSOCIATED(U, FINE_LEVEL(1)%BOX(IB)%U)`).
3. `MOVE_ALLOC` of components when the level array grows (draft driver file), and the rank-1 to rank-3/4 pointer bounds remapping already required by 0005.
4. `ERROR STOP` inside a routine called from a C++ main under `mpirun`: the run must end with a non-zero exit code and the message on stderr (checked here with gfortran/Open MPI).

## Results
- `git apply --check` passes on the tree at HEAD (patches 0003-0006 already applied there; `mesh.f90` index `172704f`, the result of 0005).
- `USE_AMREX=OFF`: the preprocessed `mesh.f90` (`gfortran -cpp -E -P`, macro undefined, blank lines dropped) is byte-identical with and without the patch (851 lines each); full OFF run of the patched tree: `check_off_bitwise.sh` `shunn3_32` 1 rank PASS (16 output files bitwise identical to the baseline), `shunn3_4mesh_32` 4 ranks PASS (47 files).
- `USE_AMREX=ON` with 0007 applied (scratch tree, gfortran 14.2, `-DFDS_AMR_FINE_B_DRAFT=ON`): builds; `tests/run_fine_guard_check.sh <build> <work> 1`: wrapper guard PASS, `--fine-b-selftest` PASS (metrics and `U` of two fine boxes reached through `POINT_TO_BOX`, a view replaces the FDS array with its own bounds, level-0 `POINT_TO_BOX(1)` and `POINT_TO_MESH(1)` unchanged), `--fine-b-abort 1` (POINT_TO_MESH on a fine number) and `2` (POINT_TO_BOX on an unknown number) stop with the message and a non-zero exit, `--fine-b-abort 3` (POINT_TO_BOX on a fine box) returns.
- Same patched ON build: final fields and step logs of `dec2`, `dec4` (4 ranks) and `dec1` byte-identical to the pre-S9 binary with `--exact-zone-sums` (`bitcmp.sh`), as for the unpatched ON build.
- Default ON build (without 0007): all driver checks pass (`run_driver_tests.sh`, decomposition, kernel checks, csmag, outputs, M2a baselines, periodic regression, zone sums, solid dump, wrapper guard).

## Fine-level kernels (item (1) of the D-056 option B work; needs 0008 too)
- Fine-ready wrappers: with 0007 + 0008 applied the driver's `fds_p_*` (TimeLoop path) and `fds_k_*` (C++ path) wrappers `visc`, `dens`, `dens_pre`, `dens_post`, `vflux`, `div1`, `div2`, `vpred`,
  `vcorr`, `p_mfd`, `p_dens`, `p_dens_pre`, `baroclinic`, `noflux`, flag set/clear call `POINT_TO_BOX` through the patched files and use `BOX_OBJ(NM)` (driver module `fds_box_obj.f90`:
  `MESHES(NM)` or `FINE_LEVEL(L)%BOX(IB)`) instead of `MESHES(NM)`; every other NM-taking `fds_p_*` wrapper calls `FDS_HOOK_L0_ONLY` and aborts for a fine number (this closes the gap that
  the S9.6 guard covered only the `fds_k_*` entry names). The CMake auto-detection defines `FDS_FINE_B` for the driver Fortran when `Source/mesh.f90` declares `FINE_LEVEL_TYPE`.
  `fds_hook_set_fine_ready(1)` is still not called by any production code path: fine boxes are not used for physics before 0007 and 0008 are validated.
- `FDS_FINE_B_SET_VIEW(L,IB,WHICH,LB,EXT,P)` (C name `fds_fine_b_set_view`) binds fields of `FINE_LEVEL(L)%VIEW(IB)` to AMReX data, as `fds_hook_set_view` does for level 0.
- Contents of a built box (`BUILD_FINE_BOX` in the draft): metrics with the FDS arithmetic, `CELL_INDEX`/`CELL` with a solid-exterior ghost ring, `EDGE(0:0)`, zero-size wall tables,
  `INTERPOLATED_MESH`, all state, field and work arrays, LES filter width and `CSD2`, background pressure and atmosphere. Not built: obstructions (none at fine levels), `OMESH`, interface walls
  without `OMESH` neighbour and domain-edge walls (see the finding below), particles and radiation (abort), slices and other output arrays.
- Test `tests/run_fine_b_shadow_check.sh <build> [steps] [work]` (needs the draft build, `-DFDS_AMR_FINE_B_DRAFT=ON`): `--run --fine-b-shadow` repeats every kernel call of a single-mesh run on a
  clone of the level-0 box at mesh number `NMESHES+1`; `--fine-b-build` does it with a box made by `BUILD_FINE_BOX` (only state arrays copied; the box's own views are bound through
  `FDS_FINE_B_SET_VIEW`). Comparison is bit by bit, including ghost layers and the global sums (`DSUM`, `USUM`, `PSUM`, `P_ZONE`). The fine box has the cells and spacing of the level-0 box,
  i.e. it is a level-1 box covering the whole domain compared with the single-level run at that resolution. Clone mode: all stages, every array, `dec1` and a 16^3 periodic LES case (`tg16`):
  bitwise. Built mode (`tg16`): the core (cells at least two from every box face) is bitwise for every stage; the shell is reported with its largest relative difference. The level-0 final
  fields of the shadowed runs are byte-identical to a plain run.
- **Finding for the Architect (design):** a fine box has no walls, so at a periodic or physical domain-edge face the cells next to the face differ from FDS's treatment (the level-0 box has
  `INTERPOLATED` periodic walls there: wall terms in `DIVERGENCE_PART_1`, `NO_FLUX`, `BAROCLINIC`, ghost values of `VISCOSITY_BC`). The driver fills those ghost layers for interior interfaces from
  the coarse level (D-055); for the domain-edge faces of a fine level the walls (and the `INTERPOLATED_MESH` coverage mask of the coarse cells) still have to be built.
- Pressure stages on a fine box (RHS, residual check, velocity error, `H_GHOST`) are not fine-ready (pass for a pressure pass).
