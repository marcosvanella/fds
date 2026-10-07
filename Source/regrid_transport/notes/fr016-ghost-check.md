# FR-016 ghost check (ns2d_16_int_1to2_refinement)

Compares the coarse-fine ghost values that the C++ hook (`make_cf_ghost_hook`, `LevelOps`) writes with the values FDS itself leaves in its mesh arrays.

## Baseline data
The frozen baseline directory holds plain outputs only (slice files, device files, restart); it has no ghost-cell dumps. The ghost values come from an instrumented FDS build
(reference-dump hook of `Source/driver/tests/refdump/`, applied to the baseline commit), run on 13 ranks (one file per mesh: `<prefix>`, `<prefix>.2`, ...):
- steps 2 and 3: `FDSREF_FILE=<prefix> FDSREF_STEPS=2,3` (the same dumps serve kernelcheck);
- step 1: `FDSREF_FILE=<prefix1> FDSREF_STEPS=1` with a SEPARATE prefix (never give it to kernelcheck, which needs steps 2,3). The record count is written when MAIN_LOOP ends,
  so let the run finish (a copy of the case with T_END=0.1 is enough); a killed run leaves a truncated dump. Step 1 holds the initial-pass records.

## Run
- Standalone (no D/DS): `-DRT_FR016_DUMP=<prefix>` adds ctest `regrid_transport_fr016_ghosts` (steps 2,3), `-DRT_FR016A_DUMP=<prefix1>` adds `regrid_transport_fr016a_ghosts_step1`;
  by hand `fr016_ghost_check <cases dir> <prefix> [step]`. The tool rebuilds the two-level state from the dump (valid cells of all meshes), applies the C++ fill and compares layers
  1 and 2 on the coarse meshes (inside the hole) and the fine mesh.
- Through a real TimeLoop stage, incl. D and DS (E4 of `tests/run_e2e_driver.sh`, needs a driver build): `RT_FR016_DUMP=<prefix>` runs steps 2,3, `RT_FR016A_DUMP=<prefix1>` step 1;
  by hand `fds_amr ns2d_16_l0.fds --rt-e2e ghost --dumps <prefix> --step <n>`. Records: DENS_P (RHO ZZ TMP RSUM), DENS_C (RHOS ZZS), VISC_P (MU KRES), DIV1_P (D), DIV1_C (DS).

## Result (steps 2 and 3)
- Fine side (RHO, ZZ, TMP, RSUM, RHOS, ZZS, MU, KRES): bitwise equal.
- Coarse side: RHO and TMP agree to 3e-15 relative, ZZ, RSUM, MU, KRES equal. RHO, ZZ, MU are uniform in this case, KRES varies and is the strong check.
- Exception: 4 of 16 corner-zone covered cells of KRES differ (see Limitation).

## Result (step 1, separate FDSREF_STEPS=1 dump; D and DS compared through the stage)
- Fine side, 64 cells per field and layer: RHO ZZ TMP RSUM RHOS ZZS MU D DS bitwise; KRES to 1.3e-14 relative.
- Coarse side, 16 cells per field and layer: ZZ RSUM ZZS MU D DS equal; RHO RHOS TMP to 3.3e-15 (4.7e-15 in the hole edge/corner cells); KRES to 4.6e-15.
- D-059 exceptions (coarse layer 1, hole edge/corner cells, 16 per field): KRES 2 cells (max rel 0.42), D 2 cells (0.5), DS 2 cells (0.52); every other field 0.
  The standalone tool (no D/DS) gives the same KRES exception.

## FDS rules implemented (routine, line in `Source/`)
- `wall.f90` ASSIGN_GHOST_VALUE, 282-388: the coarse side reads one layer of the fine mesh next to the face (probe points `init.f90` 3161-3163). The value is the area-weighted mean
  (ARO = min(1, area ratio), 330-335) of the r_t*r_t fine cells of that layer, not the r^3 volume mean of average-down. ZZ: mass-weighted, clipped to [0,1] (349+); RSUM rebuilt from ZZ;
  TMP = PBAR_P/(RSUM*RHO). Layer 2 = layer 1 for RHO, ZZ, TMP (the second-order branch, 357+, needs equal-size cells). The fine side takes the coarse value.
- `velo.f90` 545 (VISCOSITY_BC): MU, KRES, D, DS get layer 1 only, the mean over the face layer. Same N_INT_CELLS pattern at 1407 (VELOCITY_BC) and 1899.
- Early returns under EXTERNAL_GHOSTS_FILLED: `wall.f90` 306, `velo.f90` 522.
- Average-down keeps its volume value only in covered cells no ghost rule touches.

## Limitation (Architect decision)
A covered coarse cell that is the ghost cell of two faces (edge or corner of the fine patch, or a patch thinner than 4 coarse cells) holds one value in AMReX; FDS has one array per mesh
with the same collision but a different winner per mesh. The code takes the first candidate (x, y, z; layer 1 before 2) and counts `conflicts` in `CfStats`.
Options: accept, or per-face ghost storage in the fine-level kernels. Edge/corner ghost cells on the fine side (written by the piecewise-constant fill, not defined by FDS) are unverified.
