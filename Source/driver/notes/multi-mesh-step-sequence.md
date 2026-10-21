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

- With EXTERNAL_GHOSTS_FILLED, MATCH_VELOCITY returns at once; `fds_p_save_uvw` (after exchange 3 and 6) stands in for it and also sets U_GHOST/V_GHOST/W_GHOST, the ghost-cell midpoint velocities that COMPUTE_VISCOSITY (Deardorff UP/VP/WP) and VELOCITY_FLUX put in the interface ghost cells. Without them the ghosts keep the initial wind profile from INITIALIZE_MESH_VARIABLES_1.
- Known differences from stock multi-mesh (accepted, Architect ruling: the driver keeps the current-value ghosts). Stock VELOCITY_BC fills the tangential interface ghost velocities in its own edge loop (the z-edges of 'ghost cell values at edge of computational domain'), so (1) the ground-edge (vegetation, BOUNDARY_FUEL_MODEL_BC) branch, whose VEL_T reads U,V of the interface ghost cell, gets the value left by the previous call whenever the writing edge comes later in the edge list; the sub-ground ghost face (k=0) it sets is then per-mesh. The driver fills the ghosts before the routine (current values). On the 2.5 m 2x2 UGLMAT case at step 1: 60 interface cells at k=1 over all meshes, up to 4.85e-2 rel in the corrector MU (V&V, physics-matrix/RESULTS.md), spreading to V/W within the step. (2) The diagonal ghosts across a mesh-corner line (outside the box in two or three directions) are not refreshed by stock (cell fields, US/VS/WS) or a step late (U/V/W, corrector branch); the driver fills them; at step 1 FVX/FVY 1.6e-6 rel in the corner columns.
  (3) The OPEN-boundary branch of the same loop decides wind vs free-slip from U_NORM, which reads interface ghosts (W at i=0 / j=0, etc.) that a later edge or a later orientation of the same edge refreshes: at the top OPEN boundary next to an x-/y-interface (II=0 or JJ=0 side only) stock reads the previous call's value and can pick free-slip where the driver picks the wind (US(0,j,KBP1) 1.6e-3 rel at step 2, FVZ 8.4e-3).
- `FDSTL_VERIFY_STOCK_EDGE_LAG=1` (verification only, off by default, AMR build only; fds_amrex_hooks.f90, fds_step.f90 FDS_P_LAG_SAVE, GhostExchange.cpp) reproduces (1)-(3): U,V,W are saved before the fill, every ghost written by stock's 'ghost cell values at edge' block gets the key IE*8+P of its first write (P = position of the orientation in SIGN_LOOP x ORIENTATION_LOOP), and the ground-edge and OPEN-branch reads take the saved value when that key is larger than the reader's; the diagonal ghosts keep their values through the fill. With it (and the --ls-fine AREA_ADJUST fix below), steps 1-10 of the 2.5 m 2x2 UGLMAT case agree with stock 36975d7 (stage dumps, interior) to <=1e-10 rel to the field maximum in every state field except D, DS (<=2.5e-8) and DDDT (<=7e-7), which are cancellation: DIV = sum of velocity differences of O(|U|/dx ~ 3/s) that cancel to O(1e-5), so a 1e-11 m/s velocity round-off gives 1.5e-11 abs in DDDT (= 1e-12 of the term size); DT agrees to 7.6e-14 rel. Waiver (Architect-accepted, multi-mesh item 1): FVX reaches 1.1e-10 rel against the 1e-10 target, accepted as round-off growth: it appears only from step 8, its inputs (U,V,W,US,VS,WS,MU,RHO) are <=6.7e-11, and D/DS/DDDT are shown to be cancellation. The ghosts that still differ at the end of a step (376 points, corner lines: H, HS, MU, KRES, U, W) are unread: NaN in all of them after step 1 leaves step 2 bitwise unchanged.
- --ls-fine (fds_ls_grid.f90 LS_APPLY_L0) used to set the wall AREA_ADJUST to max(1,total flux fraction); DIVERGENCE_PART_1 also scales Q_CON_F by AREA_ADJUST, so the ignition cells got a different wall heat source (DS 3.4e-4 rel at step 2). It now keeps the mean LS AREA_ADJUST (the vege.f90 value at ratio 1) and encodes the fraction as total/AREA_ADJUST.
- Not reproduced, no effect up to step 10: stock exchanges PHI_LS (code 14) only after the predictor (PHI1_LS), so OMESH%PHI_LS stays at its initial value and the predictor level set of a mesh sees PHI_LS_MIN in the interface ghosts; the --ls-fine global grid sees the real neighbour. It matters once the front crosses a mesh interface.
