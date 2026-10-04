# Patch 0009 (VALIDATED): `divg.f90`, diffusive face-flux read-out and override hook

**Status: VALIDATED (oneAPI and GNU Debug).** Evidence: oneAPI (ifx 2026.1.1, Intel MPI 2021.18.1): 0 warnings (also with -check all -traceback -fpe0 -init=snan), OFF output bitwise equal to the unpatched tree, driver tests pass (including the checked-bounds build). GNU Debug (gfortran 14.2, -fcheck=all, FP traps): fine guard, fine-b shadow and flux hook checks pass. The fine-b shadow check runs on a single level-0 mesh, 1 rank only. The amended 0007 (DT_NEW(*) dummies in velo.f90) is commit a72d491fb9. Applies to the HEAD `Source/divg.f90` and on top of 0008. Independent of 0007, but the hook is only useful with the AMR driver.
Design: `notes/flux-hooks-design.md`. This patch was first announced as 0008; 0008 is the kernel-routing patch, so the flux hook is 0009.

## What it does
Two additions in `DIVERGENCE_PART_1`, both inside `#ifdef WITH_AMREX`:
- `USE FDS_FLUX_HOOKS, ONLY: FDS_HOOK_DIF_FLUX` in the declaration part;
- `CALL FDS_HOOK_DIF_FLUX(NM,LBOUND(RHO_D_DZDX,4),RHO_D_DZDX,RHO_D_DZDY,RHO_D_DZDZ)` right after the species-sum fix and `SET_EXIMDIFFLX_3D`, before "Store diffusive flux for output".

`RHO_D_DZDX/Y/Z` are the species diffusive face fluxes exactly as `DEL_RHO_D_DEL_Z` reads them (wall corrections and the species-sum fix applied). The hook (driver module `fds_flux_hooks.f90`)
returns at once when the box has no registered flux array and no override. With a registered array (mode 1) it copies the faces; with an override list (mode 2) it replaces the listed faces in the
kernel's own arrays before the divergence is formed. The patch edits no existing statement.

## `USE_AMREX=OFF` bitwise check
`WITH_AMREX` is undefined for OFF builds, so the preprocessed source is the unchanged FDS source.
- `gfortran -cpp -E -P` of `divg.f90` with the macro undefined, HEAD vs HEAD + 0008 + 0009, blank lines dropped: identical (1431 lines).
- Full OFF run of the patched tree (0003-0009 applied): see "Results".

## Results (gfortran 14.2, Open MPI 5.0.7)
- `patch -p1` applies to HEAD and on top of 0008.
- ON, `tests/run_flux_hook_check.sh` (cases `dec1` 1 rank, `dec4_np4` 4 ranks, two species): hooks registered with an empty list, and with the read-out list fed back, give final fields and step log bitwise equal to the
  plain run for ADV (generated density copy) and DIF (this patch, re-run of `DIVERGENCE_PART_1`); scaled lists change the result, so the hooks act.
- OFF run `check_off_bitwise.sh` with 0003-0009 applied (rebuilt): `shunn3_32` 1 rank PASS (16 output files bitwise identical to the baseline), `shunn3_4mesh_32` 4 ranks PASS (47 files).
