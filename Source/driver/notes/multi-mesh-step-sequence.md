# Step sequence on a multi-mesh input (shunn3_4mesh_32): driver against plain FDS

Report: on `shunn3_4mesh_32` (four 16x1x16 meshes, 2x2) the driver's T after step 1 is 6.189e-4 s, plain FDS 6.195e-4 s (step 2: 1.2998e-3 against 1.3009e-3), the same
at 1 and 4 ranks; the one-mesh `shunn3_32` agrees.

Findings (plain FDS, USE_AMREX=OFF build of the committed tree, the same input with only the `&MESH` lines changed):

| meshes | T after step 1 | T after step 2 |
|---|---|---|
| 1 (32x1x32) | 6.188e-4 | 1.2994e-3 |
| 2 (split in x) | 6.188e-4 | 1.2994e-3 |
| 2 (split in z), 2 meshes in the reverse order | 6.188e-4 | 1.2994e-3 |
| 4 (2x2, the input) | 6.195e-4 | 1.3009e-3 |
| 4 (2x2), `VELOCITY_TOLERANCE=1e-10` (33 iterations instead of 1) | 6.194e-4 | 1.3007e-3 |
| driver, 4 meshes (np 1 and np 4) and the 1-mesh variant | 6.189e-4 | 1.2998e-3 |

- The time-step rules of the driver are FDS's (`TimeStep.H`; the step-start rule, the repeat-pass rule, the DT increase of 1.1 per step): every later DT is the
  FDS value times the same factor, and the whole offset (-8.5e-4 relative in T, constant over the steps compared) comes from the DT of step 1, which FDS finds by
  the repeated first-pass reduction from the initial DT (0.312 s).
- Plain FDS itself changes the step-1 DT only for the 2x2 arrangement, where four meshes share a corner; the two-mesh splits give the one-mesh value, and so does the
  driver, whose level-0 data has no mesh interfaces (shared faces, ghost cells and corner ghost cells filled from the neighbours by AMReX). Tightening the pressure
  iteration of FDS does not remove the 2x2 difference (not a pressure-iteration effect). The cause is therefore in FDS's own treatment of the cells at an
  interior mesh corner, not in the driver's dt candidates, mesh set for the minimum or the order of the pressure pass; where exactly in FDS has not been traced.
- The remaining difference of the driver from the 1-mesh FDS (1.6e-4 in the step-1 T) is present for the 1-mesh variant too, i.e. it is not a multi-mesh effect.

Ruling: single-mesh FDS is the ground truth and a multi-mesh run must equal it; the driver is unchanged. The multi-mesh FDS first-dt offset above (6.195e-4 against
6.188e-4 s) is a known reference quirk; the candidate cause (unconfirmed) is interpolated-boundary UVW_SAVE data in the first trial CFL pass.

`tests/run_step_sequence_check.sh` compares the driver (4-mesh input at 1 and 4 ranks, and the 1-mesh input) with the SINGLE-MESH FDS baseline
`shunn3_4mesh_32__1mesh` only: the three driver runs must be identical, DT within the .csv rounding (6e-3) and T within 5e-4 over the first 40 steps. The multi-mesh FDS
run is not a reference. It is in the regular list (`run_driver_tests.sh`).
