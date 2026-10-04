# Patch 0008 (VALIDATED, with 0007): `velo.f90`, `mass.f90`, `divg.f90`, `wall.f90`, `turb.f90`, kernels find their box through `POINT_TO_BOX`

**Status: VALIDATED, goes with 0007 (oneAPI and GNU Debug).** Evidence: oneAPI (ifx 2026.1.1, Intel MPI 2021.18.1): 0 warnings (also with -check all -traceback -fpe0 -init=snan), OFF output bitwise equal to the unpatched tree, driver tests pass (including the checked-bounds build). GNU Debug (gfortran 14.2, -fcheck=all, FP traps): fine guard, fine-b shadow and flux hook checks pass. The fine-b shadow check runs on a single level-0 mesh, 1 rank only. The amended 0007 (DT_NEW(*) dummies in velo.f90) is commit a72d491fb9. It needs only `POINT_TO_BOX` of patch 0005; it is independent of
0007 but is only useful with it (without 0007 `POINT_TO_BOX` is the 0005 routine and a level-0 mesh number gives the same result as `POINT_TO_MESH`).

## What it does
Five lines per file (a macro, with its comment), placed before the `MODULE` statement of each of the five files:

    #ifdef WITH_AMREX
    #define POINT_TO_MESH POINT_TO_BOX
    #endif

Every `CALL POINT_TO_MESH(NM)` of the file (velo 12, mass 3, divg 4, wall 3, turb 16: 38 call sites, no statement is edited) is preprocessed to `CALL POINT_TO_BOX(NM)`. The kernels of these files
therefore reach the arrays of a fine-level box (`NM > SIZE(MESHES)`, `FINE_LEVEL(L)%BOX(IB)` with the views `FINE_LEVEL(L)%VIEW(IB)`, 0007) the same way they reach a level-0 box. For a level-0 number
without `BOX_VIEW` data `POINT_TO_BOX` is `POINT_TO_MESH`; with `BOX_VIEW` data the viewed arrays are re-pointed at the AMReX data (the 0005 behaviour).

Why a macro and not a rename of 38 lines: the patch is five hunks of one block each, it is easy to review and to re-apply when FDS changes the files, and with `WITH_AMREX` undefined it vanishes.

## Not routed (level 0 only; documented, each aborts or is skipped at a fine number)
- the `OMESH`/`MESHES(NOM)` branches (`velo.f90`, `wall.f90`): skipped by `EXTERNAL_GHOSTS_FILLED` (D-055) and by the zero-wall fine boxes; a fine box never reads a `MESHES(NOM)`;
- `SETTLING_VELOCITY`, cut-cell, radiation and particle paths (not reachable on a fine box: the builder aborts for particles and radiation);
- `pres.f90` (`PRESSURE_SOLVER_COMPUTE_RHS`, `PRESSURE_SOLVER_CHECK_RESIDUALS`, `COMPUTE_VELOCITY_ERROR`): they need a fine `PRHS`/`IPS`/boundary-condition set-up and are for a later pressure pass; the
  driver wrappers of these stages call `FDS_HOOK_L0_ONLY` and abort for a fine number.

## `USE_AMREX=OFF` bitwise check
`WITH_AMREX` is undefined for OFF builds, so the preprocessed source is the unchanged FDS source.
- `gfortran -cpp -E -P` of each of the five files with the macro undefined, patched vs unpatched, blank lines dropped: identical for all five files.
- Full OFF run of the patched tree (0003-0008 applied): see "Results".

## Results (gfortran 14.2, Open MPI 5.0.7)
- `git apply --check` passes on the HEAD files.
- Preprocessed sources identical with the macro undefined (above).
- OFF run `check_off_bitwise.sh`: with 0003-0008 applied (copy of the tree outside the reference, rebuilt after the 0008 files were added; HEAD already carries 0007): `shunn3_32` 1 rank PASS (16 output files bitwise identical to the baseline), `shunn3_4mesh_32` 4 ranks PASS (47 files).
- ON with 0007 + 0008 (scratch tree, `-DFDS_AMR_FINE_B_DRAFT=ON`): `tests/run_fine_b_shadow_check.sh` (see `0007-mesh-fine-level-boxes.md`, "Fine-level kernels"): every kernel stage that the
  `FDS_P_*`/`FDS_K_*` wrappers expose runs on `FINE_LEVEL(1)%BOX(1)` through `POINT_TO_BOX` and gives bitwise the result of the level-0 run, for a clone and for a box built from scratch.
- Default ON build (without 0007/0008): every driver check passes and the dec2/dec4/dec1 outputs are byte-identical to the previous binary with `--exact-zone-sums`.
