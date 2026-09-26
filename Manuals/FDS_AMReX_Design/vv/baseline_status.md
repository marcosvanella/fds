# Baseline status: FireX 36975d7, GNU/Open MPI reference build

- **Build:** FDS `FDS-6.11.1-1244-g36975d765f-AMReX` (FireX worktree HEAD 36975d765fcead401e14b094a04f910ac42eab8a), GNU Fortran + Open MPI 5.0.7 (Debian package). The binaries are Release `build/firex-36975d7/ompi_gnu_rel/fds` and Debug (the `_db` build; sha256 below).
- **HYPRE:** commit 63331f19c7 = **v2.32.0-24-g63331f19c**. The FDS run header prints "Hypre library version: 3.0.0" because CMakeLists hard-codes HYPRE_GIT_VERSION. That label is wrong; every manifest records the real version.
- **Written:** 2026-09-25. Author: V&V lead. The generator is `vv-runs/scripts/make_report.py`, which reads `vv-runs/analysis/gnu_ompi_firex-36975d7/analysis.json` and every `manifest.json`.
- **Baseline directory:** `vv-runs/baseline/gnu_ompi_firex-36975d7/`. It is read-only (`chmod -R a-w`) and has `SHA256SUMS` over every file, using relative paths. A-09 Debug runs and the Tier 2 `emb_1to1` run are in `vv-runs/a09/gnu_ompi_firex-36975d7/`, G0 is in `vv-runs/g0/…` and the restart-repeat calibration is in `vv-runs/calib/…`. Those three directories are not part of the frozen baseline.
- **Run protocol:**
  - OMP_NUM_THREADS=1.
  - Launch with `timeout -k 60 T mpirun --bind-to none -np N`, where T = max(600 s, 3× estimate).
  - Never more than 8 concurrent ranks. Runs are strictly sequential, one at a time.
  - Pre-launch gate on load and RAM (see Deviations).
- **Class note (requirements v0.3.2):**
  - Whole-run comparisons against FDS are **T2**. Only kernels are bitwise.
  - The T0/T1 results below come from repeat runs of the same binary (determinism, rank count, restart). Their only purpose is to calibrate tolerances.
  - Multi-mesh cases are **T2 references**, not T1 references.
  - Single-mesh runs are also T2-class for whole-run comparison against FDS. They are labelled "single mesh" so the A-24 single-mesh references can be told apart.

## G0: PASS

| check | result |
|---|---|
| Release binary sha256 | 9d3b5991447a70b8cd1d067e7db6fec8506479b32fd62961339f8d387e8050ec (matches expected) |
| Debug binary sha256 | f21f8e4d4cf85eabbc9a88191f64a875c77c52d8dda1145c9bb270b284f15dcf |
| `fds -V` | FDS-6.11.1-1244-g36975d765f-AMReX; Open MPI v5.0.7 |
| ns2d_16, 1 rank | exit 0, "STOP: FDS completed successfully", wall 0.91 s (FDS elapsed 0.681 s) |
| obst_activation_default, 4 ranks | exit 0, normal stop, wall 0.70 s (FDS elapsed 0.47 s) |
| Worktree | HEAD 36975d765f… before and after; `git status --short` empty before and after |
| Setup-only sweep (T_END=0) over all 941 inputs | **NOT RUN**: pending, per instruction |

## A-09: AMR prototype cases (Release, plus a Debug no-trap check)

UGLMAT derivation for `ns2d_16_int_1to2_refinement_uglmat`:
- The source is `Adaptive_Mesh_Refinement/ns2d_16_int_1to2_refinement.fds`.
- The whole record `&PRES VELOCITY_TOLERANCE=1.E-6, MAX_PRESSURE_ITERATIONS=100/` is replaced by `&PRES SOLVER='UGLMAT HYPRE' /`.
- CHID becomes `ns2d_16_int_1to2_refinement_uglmat`.
- Nothing else changes. The diff is in `derivation.diff` in the run folder.

| run | binary | ranks | wall s | est s | exit/normal stop | trap lines | pressure it. mean | max | steps |
|---|---|---|---|---|---|---|---|---|---|
| ns2d_16_emb_1to2_refinement | rel | 2 | 14.8 | 17 | 0/yes | 0 | 6.35 | 16 | 8.72e+03 |
| ns2d_16_emb_1to2_refinement__sf17 | rel | 2 | 14.7 | 17 | 0/yes | 0 | — | — | — |
| ns2d_16_int_1to2_refinement | rel | 1 | 90.9 | 31 | 0/yes | 0 | 35.75 | 46 | 4.55e+03 |
| ns2d_16_int_1to2_refinement__sf17 | rel | 1 | 93.7 | 31 | 0/yes | 0 | — | — | — |
| ns2d_16_int_1to2_refinement_uglmat | rel | 1 | 29.7 | 9 | 0/yes | 0 | 1.00 | 1 | 4.44e+03 |
| ns2d_16_int_1to2_refinement_uglmat__sf17 | rel | 1 | 30.4 | 9 | 0/yes | 0 | — | — | — |
| ns2d_16_int_1to2_refinement_uglmat__sf17_rep2 | rel | 1 | 30 | 9 | 0/yes | 0 | — | — | — |
| ns2d_16_emb_1to1_refinement | rel | 2 | 3.92 | 8 | 0/yes | 0 | 13.39 | 17 | 1.41e+03 |
| ns2d_16_emb_1to1_refinement__debug | db | 2 | 33.4 | 31.4 | 0/yes | 0 | — | — | — |
| ns2d_16_emb_1to1_refinement__sf17 | rel | 2 | 3.92 | 8 | 0/yes | 0 | — | — | — |
| ns2d_16_emb_1to2_refinement__debug | db | 2 | 102 | 119 | 0/yes | 0 | — | — | — |
| ns2d_16_int_1to2_refinement__debug | db | 1 | 254 | 727 | 0/yes | 0 | — | — | — |
| ns2d_16_int_1to2_refinement_uglmat__debug | db | 1 | 130 | 237 | 0/yes | 0 | — | — | — |

Results:
- **Status:** every A-09 run (Release and Debug) exited 0 with "FDS completed successfully". There were **0** NaN, array-bounds or FPE trap lines in the Debug runs.
- **int_1to2 wall time:** 90.9 s against a 31 s estimate.
- **Pressure iterations:**
  - The default solver on the interior 1:2 patch needs a mean of 35.75 and a max of 46 pressure iterations per step, over 4552 steps.
  - UGLMAT needs 1 iteration per step, and the run is about 3× faster.
  - emb_1to2 has a mean of 6.3 and a max of 16. emb_1to1 has a mean of 13.4 and a max of 17.
- **RMS u error against the analytic solution (ns2d.py method):**

| run | 0–30 s | 0–2π |
|---|---|---|
| int_1to2 | 1.023 | 0.2981 |
| int_1to2 UGLMAT | 1.131 | 0.2981 |
| emb_1to2 | 0.554 | 0.1307 |
| emb_1to1 (Tier 2) | 0.586 | 0.1360 |
| uniform ns2d_16 (reference) | — | 0.1307 |
| uniform ns2d_32 (reference) | — | 0.0319 |

  - The interior 1:2 refined patch is therefore **worse than the uniform coarse grid** (0.298 vs 0.131 up to 2π). This matches the known interior-refinement defect.
  - An analytic-pressure metric was tried and dropped. It did not converge with N on the uniform ns2d series (0.597, 0.684, 0.702 for N = 16, 32, 64), so the formula or normalization is not validated.
- **int_1to2, default solver vs UGLMAT (SIG_FIGS=17):**
  - The time bases differ, so the formal T1 check fails on structure.
  - Interpolated over 0–30 s:
    - UVEL: max |Δ| 3.53 (0.861 of ‖ref‖∞) at t = 18.6 s.
    - PRES: max |Δ| 6.04 (1.72 of ‖ref‖∞).
    - VISC: identical.
  - Over 0–2π:
    - UVEL: max |Δ| 2.3e-5 (1.1e-5 of ‖ref‖∞).
    - PRES: 0.66 (0.25 of ‖ref‖∞), which is consistent with a pressure-level offset.
  - UVEL divergence times: |Δ| > 1e-6 at 0.75 s, > 1e-3 at 11.38 s, > 0.1 at 15.54 s and > 1 at 16.71 s. This is the classic exponential growth of a round-off-level difference in a marginally resolved flow.
- **UGLMAT repeat determinism** (sf17 vs sf17_rep2): T0 PASS and T1 PASS (bitwise identical).

Release vs Debug (T2-style, interpolated on release times; SIG_FIGS as committed = 8 digits; max |Δ| and max |Δ|/‖ref‖∞):

| case | quantity | max abs | max abs/‖ref‖∞ | t at max (s) |
|---|---|---|---|---|
| ns2d_16_int_1to2_refinement | UVEL | 0.144 | 0.0351 | 29.85 |
| ns2d_16_int_1to2_refinement | PRES | 0.172 | 0.0488 | 29.85 |
| ns2d_16_int_1to2_refinement | VISC | 0 | 0 | 0.00 |
| ns2d_16_int_1to2_refinement_uglmat | UVEL | 0.143 | 0.0289 | 29.79 |
| ns2d_16_int_1to2_refinement_uglmat | PRES | 0.592 | 0.151 | 29.95 |
| ns2d_16_int_1to2_refinement_uglmat | VISC | 0 | 0 | 0.00 |
| ns2d_16_emb_1to2_refinement | UVEL | 0.00351 | 0.0016 | 19.53 |
| ns2d_16_emb_1to2_refinement | PRES | 0.00379 | 0.00185 | 21.30 |
| ns2d_16_emb_1to2_refinement | VISC | 0 | 0 | 0.00 |
| ns2d_16_emb_1to1_refinement | UVEL | 0 | 0 | 0.00 |
| ns2d_16_emb_1to1_refinement | PRES | 1e-08 | 4.9e-09 | 20.79 |
| ns2d_16_emb_1to1_refinement | VISC | 0 | 0 | 0.00 |

Release and Debug differ at the 0.14 level in UVEL for int_1to2 (see the table above). Debug uses different optimization and floating-point contraction, so the time step sequence differs. The difference then grows along the same chaotic path as the default-vs-UGLMAT pair. emb_1to2 differs by 1.6e-3 of ‖ref‖∞, and emb_1to1 matches to 5e-9. None of this is a Debug-detected fault: there are no traps.

## Tier 1 capture

Scope: the 40 Tier 1 inventory rows, which are 39 cases plus the DETERMINISM row. They include the 3 A-09 Tier 1 cases, the restart set and the two derived `VV/` inputs.
- Each case was run as committed and as a SIG_FIGS=17 copy.
- The DETERMINISM row became the repeat and 1-rank runs.
- On top of that come 4 alternate-rank runs and 7 GLMAT (FR-002) copies.
- That makes 91 runs in the baseline directory: 84 from the Tier 1 driver plus the 7 A-09 Release runs. **Every run exited 0 with normal completion. None timed out and none hit a trap.**

"Est s" is the inventory estimate. The FDS criterion is evaluated on the as-committed row with a numpy port of `fdsplotlib` dataplot (metric, error, tolerance), or with the case's own FDS script criterion.

| run | variant | ranks | meshes | wall s | est s | exit/normal | ref class (v0.3.2) | FDS criterion (as-committed rows) | T0/T1 repeat checks |
|---|---|---|---|---|---|---|---|---|---|
| 1_step_2_step_compare | as committed | 4 | 4 | 51.1 | 41 | 0/yes | T2 ref (multi-mesh) | PASS (Q2S: relative err 1.51e-04 vs tol 0.005, metric end); PASS (Q4S: relative err 1.01e-05 vs tol 0.005, metric end) | — |
| 1_step_2_step_compare__sf17 | SIG_FIGS=17 copy | 4 | 4 | 50.5 | 41 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | 1_step_2_step_compare 4 vs 1 rank: T0 FAIL / T1 FAIL |
| 1_step_2_step_compare__sf17_np1 | SIG_FIGS=17 copy at 1 rank(s): chaotic case: alternate rank  | 1 | 4 | 187 | 41 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | 1_step_2_step_compare 4 vs 1 rank: T0 FAIL / T1 FAIL |
| dancing_eddies_1mesh | as committed | 1 | 1 | 285 | 46 | 0/yes | T2 ref (single mesh) | no FDS criterion for this case | — |
| dancing_eddies_1mesh__sf17 | SIG_FIGS=17 copy | 1 | 1 | 301 | 46 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |
| dancing_eddies_default | as committed | 4 | 4 | 84.9 | 18 | 0/yes | T2 ref (multi-mesh) | PASS (pres: absolute err 1.74e-04 vs tol 0.01, metric end) | — |
| dancing_eddies_default__sf17 | SIG_FIGS=17 copy | 4 | 4 | 82.7 | 18 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | dancing_eddies_default repeat (4 ranks): T0 PASS / T1 PASS |
| dancing_eddies_default__sf17_glmat | GLMAT copy (FR-002): SIG_FIGS=17 + &PRES SOLVER='GLMAT' + CH | 4 | 4 | 118 | 18 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | — |
| dancing_eddies_default__sf17_rep2 | SIG_FIGS=17 copy, determinism repeat #2 (inventory row 24) | 4 | 4 | 82.2 | 18 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | dancing_eddies_default repeat (4 ranks): T0 PASS / T1 PASS |
| dancing_eddies_uglmat_refine | as committed | 4 | 4 | 99.3 | 18 | 0/yes | T2 ref (multi-mesh) | plot-only dataplot rows (no pass/fail); they reference dancing_eddies_tight_devc.csv, which is not in Tier 1 | — |
| dancing_eddies_uglmat_refine__sf17 | SIG_FIGS=17 copy | 4 | 4 | 98.5 | 18 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | — |
| divergence_test_2 | as committed | 1 | 1 | 0.98 | 3 | 0/yes | T2 ref (single mesh) | PASS (Vdot\|Vdot: absolute err 4.83e-16, 5.17e-16 vs tol 1e-10, metric end) | — |
| divergence_test_2__sf17 | SIG_FIGS=17 copy | 1 | 1 | 0.94 | 3 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |
| divergence_test_3 | as committed | 4 | 4 | 84.5 | 13 | 0/yes | T2 ref (multi-mesh) | PASS (div_min\|div_max: absolute err 1.43e-12, 6.36e-13 vs tol 1e-10, metric end) | — |
| divergence_test_3__sf17 | SIG_FIGS=17 copy | 4 | 4 | 79.7 | 13 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | — |
| divergence_test_3__sf17_glmat | GLMAT copy (FR-002): SIG_FIGS=17 + &PRES SOLVER='GLMAT' + CH | 4 | 4 | 135 | 13 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | — |
| energy_budget_tmix | as committed | 1 | 1 | 1.97 | 3 | 0/yes | T2 ref (single mesh) | PASS (Temp: relative err 2.37e-04 vs tol 0.01, metric end) | — |
| energy_budget_tmix__sf17 | SIG_FIGS=17 copy | 1 | 1 | 2.02 | 3 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |
| lapse_rate | as committed | 4 | 4 | 0.95 | 3 | 0/yes | T2 ref (multi-mesh) | PASS (T: relative err 1.23e-05 vs tol 0.001, metric end); PASS (P: relative err 5.82e-07 vs tol 0.001, metric end) | — |
| lapse_rate__sf17 | SIG_FIGS=17 copy | 4 | 4 | 0.95 | 3 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | lapse_rate 4 vs 1 rank: T0 FAIL / T1 PASS |
| lapse_rate__sf17_glmat | GLMAT copy (FR-002): SIG_FIGS=17 + &PRES SOLVER='GLMAT' + CH | 4 | 4 | 0.95 | 3 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | — |
| lapse_rate__sf17_np1 | SIG_FIGS=17 copy at 1 rank(s): firebot rank count (inventory | 1 | 4 | 1.03 | 3 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | lapse_rate 4 vs 1 rank: T0 FAIL / T1 PASS |
| ns2d_16 | as committed | 1 | 1 | 0.85 | 3 | 0/yes | T2 ref (single mesh) | no FDS pass/fail (ns2d_error plot N/A); ns2d.py RMS u error 1.3075e-01 | — |
| ns2d_16__sf17 | SIG_FIGS=17 copy | 1 | 1 | 0.85 | 3 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |
| ns2d_16_emb_1to2_refinement | as committed | 2 | 2 | 14.8 | 17 | 0/yes | T2 ref (multi-mesh) | no FDS criterion (commented out, A-09); V&V: RMS u err vs analytic 0.5545 (0–30 s), 0.1307 (0–2π) | — |
| ns2d_16_emb_1to2_refinement__sf17 | SIG_FIGS=17 copy | 2 | 2 | 14.7 | 17 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | — |
| ns2d_16_int_1to2_refinement | as committed | 1 | 13 | 90.9 | 31 | 0/yes | T2 ref (multi-mesh) | no FDS criterion (commented out, A-09); V&V: RMS u err vs analytic 1.0229 (0–30 s), 0.2981 (0–2π) | — |
| ns2d_16_int_1to2_refinement__sf17 | SIG_FIGS=17 copy | 1 | 13 | 93.7 | 31 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | — |
| ns2d_16_int_1to2_refinement_uglmat | derived VV/ns2d_16_int_1to2_refinement_uglmat | 1 | 13 | 29.7 | 9 | 0/yes | T2 ref (multi-mesh) | no FDS criterion (commented out, A-09); V&V: RMS u err vs analytic 1.1307 (0–30 s), 0.2981 (0–2π) | — |
| ns2d_16_int_1to2_refinement_uglmat__sf17 | derived uglmat + SIG_FIGS=17 | 1 | 13 | 30.4 | 9 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | A09 uglmat repeat: T0 PASS / T1 PASS |
| ns2d_16_int_1to2_refinement_uglmat__sf17_rep2 | derived uglmat + SIG_FIGS=17, determinism repeat #2 | 1 | 13 | 30 | 9 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | A09 uglmat repeat: T0 PASS / T1 PASS |
| ns2d_16_nupt1 | as committed | 1 | 1 | 0.72 | 3 | 0/yes | T2 ref (single mesh) | plot-only rows (N/A / Convergent Series): no FDS pass/fail; ns2d.py RMS u error 5.2227e-02 | — |
| ns2d_16_nupt1__sf17 | SIG_FIGS=17 copy | 1 | 1 | 0.71 | 3 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |
| ns2d_32 | as committed | 1 | 1 | 4.02 | 4 | 0/yes | T2 ref (single mesh) | no FDS pass/fail (ns2d_error plot N/A); ns2d.py RMS u error 3.1949e-02 | — |
| ns2d_32__sf17 | SIG_FIGS=17 copy | 1 | 1 | 4.12 | 4 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |
| ns2d_32_nupt1 | as committed | 1 | 1 | 3.28 | 4 | 0/yes | T2 ref (single mesh) | plot-only rows (N/A / Convergent Series): no FDS pass/fail; ns2d.py RMS u error 1.3062e-02 | — |
| ns2d_32_nupt1__sf17 | SIG_FIGS=17 copy | 1 | 1 | 3.24 | 4 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |
| ns2d_64 | as committed | 1 | 1 | 28.8 | 9 | 0/yes | T2 ref (single mesh) | no FDS pass/fail (ns2d_error plot N/A); ns2d.py RMS u error 8.0479e-03; observed order 32→64 1.99 (V&V gate ≥1.8 PASS) | — |
| ns2d_64__sf17 | SIG_FIGS=17 copy | 1 | 1 | 28.5 | 9 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |
| ns2d_64_nupt1 | as committed | 1 | 1 | 23.6 | 9 | 0/yes | T2 ref (single mesh) | PASS (UVEL: relative err 1.70e-04 vs tol 0.01, metric mean); ns2d.py RMS u error 3.3741e-03; observed order 32→64 1.95 (V&V gate ≥1.8 PASS) | — |
| ns2d_64_nupt1__sf17 | SIG_FIGS=17 copy | 1 | 1 | 23.6 | 9 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |
| ns2d_8 | as committed | 1 | 1 | 0.38 | 3 | 0/yes | T2 ref (single mesh) | no FDS pass/fail (ns2d_error plot N/A); ns2d.py RMS u error 4.8220e-01 | — |
| ns2d_8__sf17 | SIG_FIGS=17 copy | 1 | 1 | 0.38 | 3 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |
| ns2d_8_nupt1 | as committed | 1 | 1 | 0.34 | 3 | 0/yes | T2 ref (single mesh) | plot-only rows (N/A / Convergent Series): no FDS pass/fail; ns2d.py RMS u error 2.0111e-01 | — |
| ns2d_8_nupt1__sf17 | SIG_FIGS=17 copy | 1 | 1 | 0.35 | 3 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |
| obst_activation_default | as committed | 4 | 4 | 0.71 | 3 | 0/yes | T2 ref (multi-mesh) | PASS (D_max: absolute err 3.71e-14 vs tol 1e-13, metric max) | — |
| obst_activation_default__sf17 | SIG_FIGS=17 copy | 4 | 4 | 0.71 | 3 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | — |
| obst_activation_default__sf17_glmat | GLMAT copy (FR-002): SIG_FIGS=17 + &PRES SOLVER='GLMAT' + CH | 4 | 4 | 0.95 | 3 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | — |
| obst_activation_ulmat | as committed | 4 | 4 | 0.83 | 3 | 0/yes | T2 ref (multi-mesh) | FAIL (D_max: absolute err 1.77e-12 vs tol 1e-13, metric max) | — |
| obst_activation_ulmat__sf17 | SIG_FIGS=17 copy | 4 | 4 | 0.81 | 3 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | — |
| obst_activation_ulmat__sf17_glmat | GLMAT copy (FR-002): SIG_FIGS=17 + &PRES SOLVER='GLMAT' + CH | 4 | 4 | 0.9 | 3 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | — |
| obst_coarse_fine_interface | as committed | 2 | 2 | 1.03 | 3 | 0/yes | T2 ref (multi-mesh) | PASS (DP: absolute err 1.47e+00 vs tol 2, metric end) | — |
| obst_coarse_fine_interface__sf17 | SIG_FIGS=17 copy | 2 | 2 | 1.01 | 3 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | obst_coarse_fine_interface 2 vs 1 rank: T0 FAIL / T1 PASS |
| obst_coarse_fine_interface__sf17_np1 | SIG_FIGS=17 copy at 1 rank(s): firebot rank count (inventory | 1 | 2 | 0.99 | 3 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | obst_coarse_fine_interface 2 vs 1 rank: T0 FAIL / T1 PASS |
| random_meshes | as committed | 4 | 5 | 59.2 | 7 | 0/yes | T2 ref (multi-mesh) | PASS (Vdot: relative err 2.46e-03 vs tol 0.02, metric end) | — |
| random_meshes__sf17 | SIG_FIGS=17 copy | 4 | 5 | 59.9 | 7 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | — |
| restart_test1_continuous | derived VV/restart_test1_continuous (copy of 1a, T_END=10, n | 1 | 1 | 56.8 | 17 | 0/yes | T2 ref (single mesh) | no FDS criterion (firebot restart set, no dataplot row); see restart section | — |
| restart_test1_continuous__sf17 | derived VV/restart_test1_continuous (copy of 1a, T_END=10, n | 1 | 1 | 56.9 | 17 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |
| restart_test1a | as committed | 1 | 1 | 36.5 | 10 | 0/yes | T2 ref (single mesh) | no FDS criterion (firebot restart set, no dataplot row); see restart section | — |
| restart_test1a__sf17 | SIG_FIGS=17 copy | 1 | 1 | 35.6 | 10 | 0/yes | T2 ref (single mesh) | (see as-committed row) | RESTART repeat: T0 PASS / T1 PASS |
| restart_test1b | as committed; restart from restart_test1a outputs (CHID rest | 1 | 1 | 43.7 | 17 | 0/yes | T2 ref (single mesh) | no FDS criterion (firebot restart set, no dataplot row); see restart section | — |
| restart_test1b__sf17 | SIG_FIGS=17 copy; restart from restart_test1a__sf17 outputs  | 1 | 1 | 44.2 | 17 | 0/yes | T2 ref (single mesh) | (see as-committed row) | RESTART repeat: T0 PASS / T1 PASS |
| saad_512_cfl_1 | as committed | 1 | 1 | 2.85 | 3 | 0/yes | T2 ref (single mesh) | not used by saad_mms_temporal_error.py (only p25/p125/p0625 enter the order estimate); no FDS criterion | — |
| saad_512_cfl_1__sf17 | SIG_FIGS=17 copy | 1 | 1 | 2.99 | 3 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |
| saad_512_cfl_p0625 | as committed | 1 | 1 | 30.3 | 10 | 0/yes | T2 ref (single mesh) | saad_mms_temporal_error.py (p25/p125/p0625 triple): PASS (L2 order rho=2.0053, Z=2.0053, need ≥ 1.99) | — |
| saad_512_cfl_p0625__sf17 | SIG_FIGS=17 copy | 1 | 1 | 29 | 10 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |
| saad_512_cfl_p125 | as committed | 1 | 1 | 16.1 | 6 | 0/yes | T2 ref (single mesh) | saad_mms_temporal_error.py (p25/p125/p0625 triple): PASS (L2 order rho=2.0053, Z=2.0053, need ≥ 1.99) | — |
| saad_512_cfl_p125__sf17 | SIG_FIGS=17 copy | 1 | 1 | 16.6 | 6 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |
| saad_512_cfl_p25 | as committed | 1 | 1 | 8.72 | 5 | 0/yes | T2 ref (single mesh) | saad_mms_temporal_error.py (p25/p125/p0625 triple): PASS (L2 order rho=2.0053, Z=2.0053, need ≥ 1.99) | — |
| saad_512_cfl_p25__sf17 | SIG_FIGS=17 copy | 1 | 1 | 8.65 | 5 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |
| saad_512_cfl_p5 | as committed | 1 | 1 | 4.78 | 4 | 0/yes | T2 ref (single mesh) | not used by saad_mms_temporal_error.py (only p25/p125/p0625 enter the order estimate); no FDS criterion | — |
| saad_512_cfl_p5__sf17 | SIG_FIGS=17 copy | 1 | 1 | 5 | 4 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |
| shunn3_128 | as committed | 1 | 1 | 40.2 | 20 | 0/yes | T2 ref (single mesh) | shunn_mms.py limits apply only at N=512 (Tier 2): not applicable. V&V: L2 e_rho=2.201e-03, e_Z=5.987e-04, e_u=3.014e-04, e_H=5.058e-03; observed order 32→128: rho 1.95, Z 1.69, u 2.00, H 0.92 (V&V gate ≥1.8: FAIL ['e_Z', 'e_H']) | — |
| shunn3_128__sf17 | SIG_FIGS=17 copy | 1 | 1 | 40.7 | 20 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |
| shunn3_32 | as committed | 1 | 1 | 1.39 | 3 | 0/yes | T2 ref (single mesh) | shunn_mms.py limits apply only at N=512 (Tier 2): not applicable. V&V: L2 e_rho=3.300e-02, e_Z=6.245e-03, e_u=4.798e-03, e_H=1.815e-02 | — |
| shunn3_32__sf17 | SIG_FIGS=17 copy | 1 | 1 | 1.39 | 3 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |
| shunn3_4mesh_128 | as committed | 4 | 4 | 8.51 | 9 | 0/yes | T2 ref (multi-mesh) | no FDS script uses it. V&V: L2 e_rho=2.915e-03, e_u=5.405e-04 vs shunn3_128 e_rho=2.201e-03 | — |
| shunn3_4mesh_128__sf17 | SIG_FIGS=17 copy | 4 | 4 | 8.18 | 9 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | — |
| shunn3_4mesh_128__sf17_glmat | GLMAT copy (FR-002): SIG_FIGS=17 + &PRES SOLVER='GLMAT' + CH | 4 | 4 | 11.1 | 9 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | — |
| shunn3_64 | as committed | 1 | 1 | 5.58 | 5 | 0/yes | T2 ref (single mesh) | shunn_mms.py limits apply only at N=512 (Tier 2): not applicable. V&V: L2 e_rho=8.852e-03, e_Z=2.091e-03, e_u=1.362e-03, e_H=9.745e-03 | — |
| shunn3_64__sf17 | SIG_FIGS=17 copy | 1 | 1 | 5.55 | 5 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |
| soborot_superbee_square_wave_128 | as committed | 4 | 4 | 10.7 | 5 | 0/yes | T2 ref (multi-mesh) | soborot_mass_transport.py: PASS (L1(4 mesh)=8.568478e-03 ≤ 1e-2; \|L1(4)−L1(1)\|=0.00e+00 ≤ 1e-10) | — |
| soborot_superbee_square_wave_128_1mesh | as committed | 1 | 1 | 38.1 | 10 | 0/yes | T2 ref (single mesh) | soborot_mass_transport.py: L1(1 mesh)=8.568478e-03 (pair criterion above) | — |
| soborot_superbee_square_wave_128_1mesh__sf17 | SIG_FIGS=17 copy | 1 | 1 | 37.3 | 10 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |
| soborot_superbee_square_wave_128__sf17 | SIG_FIGS=17 copy | 4 | 4 | 11 | 5 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | soborot_superbee_128 4 vs 1 rank: T0 FAIL / T1 FAIL |
| soborot_superbee_square_wave_128__sf17_glmat | GLMAT copy (FR-002): SIG_FIGS=17 + &PRES SOLVER='GLMAT' + CH | 4 | 4 | 11.1 | 5 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | — |
| soborot_superbee_square_wave_128__sf17_np1 | SIG_FIGS=17 copy at 1 rank(s): determinism: 1-vs-4 ranks (in | 1 | 4 | 40.6 | 5 | 0/yes | T2 ref (multi-mesh) | (see as-committed row) | soborot_superbee_128 4 vs 1 rank: T0 FAIL / T1 FAIL |
| species_conservation_1 | as committed | 1 | 1 | 18 | 6 | 0/yes | T2 ref (single mesh) | PASS (PROPANE: relative err 9.12e-03 vs tol 0.01, metric end) | — |
| species_conservation_1__sf17 | SIG_FIGS=17 copy | 1 | 1 | 18 | 6 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |
| species_conservation_2 | as committed | 1 | 1 | 2.97 | 4 | 0/yes | T2 ref (single mesh) | PASS (M1\|M2: absolute err 0.00e+00, 3.78e-03 vs tol 0.01, metric end) | — |
| species_conservation_2__sf17 | SIG_FIGS=17 copy | 1 | 1 | 2.88 | 4 | 0/yes | T2 ref (single mesh) | (see as-committed row) | — |

**Totals: 91 runs in the Tier 1 table (84 Tier 1 capture runs plus 7 A-09 Release runs), 0 failed.**

- Sum of wall time for all of them: 3017 s (50.3 min), against a summed estimate of 913 s.
- Excluding A-09: 2712 s (45.2 min).
- As-committed runs only (39 runs, one per case): 1193 s (19.9 min), against the summed inventory estimate of 390 s and the test-plan figure of ~8 min (480 s). That is a factor of 2.5.

**Caveat on wall times:**
- Other jobs ran on the development machine throughout (AMReX prototype runs and Intel oneAPI/cmake builds), with load average between 3 and 18.
- 33 of 113 launches went ahead with load1 ≥ 2 under the busy-core fallback (see Deviations).
- The measured times are therefore upper bounds.
- The large overruns against the estimates do not come from contention, though; they are real and should be re-tiered (see Re-tier suggestions).
  - One example: the 1-rank dancing_eddies runs took 261–301 s, while the 4-rank run took 83 s. That is consistent with 4 × 83 = 332 CPU-s, so the estimate is what is low.
  - Repeat pairs agree within a few % (for example dancing_eddies_default: 84.9, 82.7 and 82.2 s), so the noise is small next to the 3–8× overruns.

### FDS criteria summary (as-committed runs)

- **PASS:**
  - 1_step_2_step_compare: M rel. error 1.51e-4, H 1.01e-5 (tolerance 5e-3).
  - dancing_eddies_default pres: 1.74e-4 (tolerance 1e-2).
  - divergence_test_2: 4.8e-16 / 5.2e-16 (tolerance 1e-10).
  - divergence_test_3: 1.43e-12 / 6.4e-13 (tolerance 1e-10).
  - energy_budget_tmix: 2.37e-4 (tolerance 1e-2).
  - lapse_rate: T 1.23e-5, P 5.8e-7 (tolerance 1e-3).
  - ns2d_64_nupt1: 1.70e-4 (tolerance 1e-2).
  - obst_activation_default: 3.7e-14 (tolerance 1e-13).
  - obst_coarse_fine_interface DP: abs 1.47 (tolerance 2).
  - random_meshes: 2.46e-3 (tolerance 2e-2).
  - species_conservation_1: **9.12e-3 (tolerance 1e-2, only 9 % margin)**.
  - species_conservation_2: 0 / 3.78e-3 (tolerance 1e-2).
- **FAIL: obst_activation_ulmat.** D_max abs error 1.77e-12 against tolerance 1e-13, with a max divergence of about 1.9e-12. The ULMAT solver on this GNU/HYPRE build does not reach the 1e-13 divergence the FDS criterion expects. The default (FFT) case passes at 3.7e-14. This is a finding against the reference build itself, and FDS-AMR must not be held to 1e-13 on ULMAT until it is explained.
- **Script criteria:**
  - soborot_mass_transport: L1(4 mesh) = L1(1 mesh) = 8.568478e-3 ≤ 1e-2, |Δ| = 0 ≤ 1e-10. PASS.
  - saad temporal L2 order: ρ 2.0053, Z 2.0053 (≥ 1.99). PASS.
  - ns2d spatial order (32→64), V&V gate ≥ 1.8: ν=0 1.99, ν=0.1 1.95. PASS. The RMS u error for ν=0 is 0.482 / 0.1307 / 0.0319 / 0.00805 for N = 8 / 16 / 32 / 64.
- **shunn3 MMS spatial L2 (V&V order gate ≥ 1.8, 32→128):** ρ 1.95, u 2.00, **Z 1.69, H 0.92 (FAIL)**.
  - The FDS `shunn_mms.py` limits only apply at N=512, which is Tier 2, so they are not applicable here.
  - The gate is too strict for Z and H at these resolutions: H is still pre-asymptotic, with order 0.95 for 64→128, while Z reaches 1.80 for 64→128. It needs re-calibration, not a code fix.
  - Errors:

| N | e_rho | e_Z | e_u | e_H |
|---|---|---|---|---|
| 32 | 3.30e-2 | 6.25e-3 | 4.80e-3 | 1.82e-2 |
| 64 | 8.85e-3 | 2.09e-3 | 1.36e-3 | 9.75e-3 |
| 128 | 2.20e-3 | 5.99e-4 | 3.01e-4 | 5.06e-3 |

- **Plot-only dataplot rows** (no pass/fail): ns2d_8/16/32(_nupt1) "Convergent Series" and the dancing_eddies iteration and CPU rows. The latter reference `dancing_eddies_tight_devc.csv`, which is not in Tier 1.

## Determinism, rank-count and repeat checks (T0/T1; these calibrate the tolerances)

- T0 means bitwise identical.
- T1 is `compare_csv.py --class T1`.
- When the time base differs (the dt sequence is not identical), the formal T1 check fails on structure. In that case the last column gives the interpolated max |Δ|.

| comparison | file | T0 | T1 | T1 max abs | T1 max rel | T0 failing cols | T1 failing cols | note | interp max abs (if time base differs) |
|---|---|---|---|---|---|---|---|---|---|
| A09 uglmat repeat | _devc.csv | PASS | PASS | 0 | 0 | — | — | — | — |
| A09 uglmat repeat | _hrr.csv | PASS | PASS | 0 | 0 | — | — | — | — |
| A09 uglmat repeat | _steps.csv | PASS | PASS | 0 | 0 | — | — | — | — |
| DET dancing_eddies_default repeat (4 ranks) | _devc.csv | PASS | PASS | 0 | 0 | — | — | — | — |
| DET dancing_eddies_default repeat (4 ranks) | _hrr.csv | PASS | PASS | 0 | 0 | — | — | — | — |
| DET dancing_eddies_default repeat (4 ranks) | _steps.csv | PASS | PASS | 0 | 0 | — | — | — | — |
| DET soborot_superbee_128 4 vs 1 rank | _devc.csv | PASS | PASS | 0 | 0 | — | — | — | — |
| DET soborot_superbee_128 4 vs 1 rank | _hrr.csv | FAIL | FAIL | 6.93e-15 | 1.87 | Q_CONV, Q_COND, Q_PRES, Q_TOTAL | Q_PRES | — | — |
| DET soborot_superbee_128 4 vs 1 rank | _steps.csv | PASS | PASS | 0 | 0 | — | — | — | — |
| RANK 1_step_2_step_compare 4 vs 1 rank | _devc.csv | FAIL | FAIL | — | — | — | — | FAIL time base differs: ref 21 rows [0,10], test 21 rows [0,10] | 2.38e-10 |
| RANK 1_step_2_step_compare 4 vs 1 rank | _hrr.csv | FAIL | PASS | 1.19e-07 | 0.826 | HRR, Q_RADI, Q_CONV, Q_COND, Q_DIFF | — | — | — |
| RANK 1_step_2_step_compare 4 vs 1 rank | _steps.csv | PASS | PASS | 0 | 0 | — | — | — | — |
| RANK obst_coarse_fine_interface 2 vs 1 rank | _devc.csv | FAIL | PASS | 2.91e-11 | 2.43e-16 | P1, P2 | — | — | — |
| RANK obst_coarse_fine_interface 2 vs 1 rank | _hrr.csv | FAIL | PASS | 1.46e-11 | 3.65e-12 | Q_RADI, Q_CONV, Q_COND, Q_PRES, Q_ENTH | — | — | — |
| RANK obst_coarse_fine_interface 2 vs 1 rank | _steps.csv | PASS | PASS | 0 | 0 | — | — | — | — |
| RANK obst_coarse_fine_interface 2 vs 1 rank | _ctrl.csv | PASS | PASS | 0 | 0 | — | — | — | — |
| RANK lapse_rate 4 vs 1 rank | _hrr.csv | FAIL | PASS | 2.44e-15 | 4.26e-14 | Q_RADI, Q_CONV, Q_COND, Q_ENTH, Q_TOTAL | — | — | — |
| RANK lapse_rate 4 vs 1 rank | _line.csv | PASS | PASS | 0 | 0 | — | — | — | — |
| RANK lapse_rate 4 vs 1 rank | _steps.csv | PASS | PASS | 0 | 0 | — | — | — | — |
| RESTART repeat: restart_test1a+1b (sf17) run twice | _devc.csv | PASS | PASS | 0 | 0 | — | — | — | — |
| RESTART repeat: restart_test1a+1b (sf17) run twice | _hrr.csv | PASS | PASS | 0 | 0 | — | — | — | — |
| RESTART repeat: restart_test1a+1b (sf17) run twice | _steps.csv | PASS | PASS | 0 | 0 | — | — | — | — |
| RESTART repeat: restart_test1a (sf17) run twice | _devc.csv | PASS | PASS | 0 | 0 | — | — | — | — |
| RESTART repeat: restart_test1a (sf17) run twice | _hrr.csv | PASS | PASS | 0 | 0 | — | — | — | — |
| RESTART repeat: restart_test1a (sf17) run twice | _steps.csv | PASS | PASS | 0 | 0 | — | — | — | — |

Reading:
- **Same binary, same rank count, repeated:** bitwise identical. That covers dancing_eddies_default at 4 ranks (all CSVs, T0 PASS) and the UGLMAT int_1to2 run (T0 PASS). Run-to-run determinism holds at fixed decomposition, including HYPRE.
- **Different rank count:** not bitwise identical, because MPI reduction order changes the round-off.
  - lapse_rate (4 vs 1): T1 PASS; hrr differences are ≤ 2.4e-15 abs.
  - obst_coarse_fine_interface (2 vs 1): T1 PASS; devc max |Δ| 2.9e-11 on P ~ 1e5, relative 2.4e-16.
  - soborot_128 (4 vs 1): devc tracers bitwise identical. hrr differs only in near-zero columns, e.g. Q_PRES max |Δ| 6.9e-15. The rel. 1.87 on a near-zero value is what fails T1 there. **T1 needs an absolute floor**, around 1e-13 × the column scale.
  - 1_step_2_step_compare (4 vs 1): the time stamps differ at up to about 3e-12 relative, because dt depends on global reductions. Interpolated hrr differences are ≤ 1.2e-7 abs, and T1 passes on hrr.

## Restart (restart_test1a → 1b vs continuous)

- restart_test1a runs to T_END=5 s with DT_RESTART=5.
- restart_test1b restarts it (same CHID `restart_test1a`, RESTART=T) to T_END=10 s.
- The continuous reference is derived from 1a: T_END=10, no DT_RESTART, CHID `restart_test1_continuous`.
- Structural finding:
  - FDS clips 1a's last step to land exactly on T_END: step 344 has dt 8.03e-3, against about 1.3e-2 in the continuous run.
  - The restarted trajectory is therefore a different discrete trajectory from the continuous one, by construction.
  - The step logs are identical through step 300 (t = 4.43734), and the continuous run reaches t=10 in 727 steps against 729.
  - With this case pair, T0/T1 row matching is not applicable to "restart vs continuous".
- A second confounder: FDS DEVC and HRR outputs are **time-averaged over the output interval** (DEVC TIME_AVERAGED defaults to .TRUE.). The interval is T_END/NFRAMES: 0.0125 s in 1a, 0.025 s in the continuous run. 1b inherits 1a's 0.0125 s, because DT_DEVC is read back from the restart file (dump.f90). Rows at the same time stamp therefore average different windows, even while the trajectories are still identical.
- What is measured:
  - (a) 1a vs continuous, at time stamps that occur exactly in both files.
  - (b) 1b over 5–10 s vs continuous, T2-style interpolated.
  - (c) Restart repeatability (the T0 calibration): 1a+1b run twice.

| run set | file | 1a vs continuous at exactly common times (t ≤ 5 s) | 1b (5–10 s) vs continuous, worst 3 columns, interpolated |
|---|---|---|---|
| restart_test1__sf17 | _devc.csv | 60/200 rows identical; last common t=4.9785; first diff t=0.422; max\|Δ\| 0.28 | null#7: 0.00765 (2.34 of ‖ref‖∞); null#9: 0.0197 (1.72 of ‖ref‖∞); null#17: 2.19 (1.71 of ‖ref‖∞) |
| restart_test1__sf17 | _hrr.csv | 60/200 rows identical; last common t=4.9785; first diff t=0.422; max\|Δ\| 45.5 | Q_ENTH: 51.1 (1.23 of ‖ref‖∞); Q_TOTAL: 50.3 (0.945 of ‖ref‖∞); Q_COND: 4.91 (0.656 of ‖ref‖∞) |
| restart_test1 (committed) | _devc.csv | 60/200 rows identical; last common t=4.9785; first diff t=0.422; max\|Δ\| 0.28 | null#6: 0.00321 (2.46 of ‖ref‖∞); null#7: 0.00765 (2.34 of ‖ref‖∞); null#5: 0.00144 (2.01 of ‖ref‖∞) |
| restart_test1 (committed) | _hrr.csv | 60/200 rows identical; last common t=4.9785; first diff t=0.422; max\|Δ\| 45.5 | Q_ENTH: 51.1 (1.23 of ‖ref‖∞); Q_TOTAL: 50.3 (0.945 of ‖ref‖∞); Q_COND: 4.91 (0.656 of ‖ref‖∞) |

Reading:
- **Restart repeatability** (c): T0 PASS. Running the whole 1a → restart → 1b chain twice gives bitwise-identical devc, hrr and steps CSVs, and 1a alone is also bitwise identical. Reading and writing the restart file is deterministic.
- **1a vs continuous** (a):
  - 60 of the 200 common time stamps are bitwise identical in devc and hrr, and the first difference is at t = 0.42 s.
  - The step logs are identical through step 300, so these differences come from the averaging windows, not from a trajectory difference.
  - Up to t ≈ 4.98 s the common time stamps match exactly, which shows the dt sequence is identical until the T_END clip.
- **1b (5–10 s) vs continuous** (b), T2-style:
  - HRR max |Δ| 65.9 kW (0.22 of ‖ref‖∞ = 299 kW).
  - Q_RADI 0.22, Q_CONV 0.12, Q_TOTAL 0.95 and Q_ENTH 1.23 of ‖ref‖∞. The last two are net terms, small next to HRR.
  - Point velocities are O(1) of ‖ref‖∞: max |Δ| up to 2.2 m/s.
  - The plume-whirl case is turbulent and chaotic, so a one-step perturbation (the clip) grows to O(1) within a few seconds.
- **Conclusion:**
  - This case pair can only calibrate restart reproducibility, which is bitwise.
  - It cannot calibrate restart-vs-continuous equivalence.
  - For a T0/T1 restart-equivalence test, FDS-AMR needs a restart point that the continuous run hits exactly without clipping (for example DT fixed with LOCK_TIME_STEP, and T_END a multiple of DT), and matching NFRAMES/DT_DEVC in both halves. Proposed as a VV-owned case.

## GLMAT copies (FR-002, informative)

Each copy is the SIG_FIGS=17 input with `&PRES SOLVER='GLMAT'` added and CHID suffixed with `_glmat`. It is compared against the FFT sf17 run, T2-style.

| case | file | quantity (worst 4 by abs/‖ref‖∞) | max abs | ‖ref‖∞ | max abs/‖ref‖∞ |
|---|---|---|---|---|---|
| obst_activation_default | _devc.csv | D_min | 1.79e-12 | 1.04e-13 | 17.2 |
| obst_activation_default | _devc.csv | D_max | 1.77e-12 | 1.37e-13 | 12.9 |
| obst_activation_default | _hrr.csv | Q_COND | 3.01e-14 | 1.78e-15 | 16.9 |
| obst_activation_default | _hrr.csv | Q_CONV | 0.388 | 0.388 | 1 |
| obst_activation_default | _hrr.csv | Q_TOTAL | 0.388 | 0.388 | 1 |
| obst_activation_default | _hrr.csv | Q_ENTH | 2.6 | 42.1 | 0.0618 |
| obst_activation_ulmat | _devc.csv | D_max | 1.08e-12 | 1.87e-12 | 0.576 |
| obst_activation_ulmat | _devc.csv | D_min | 1.2e-12 | 2.09e-12 | 0.574 |
| obst_activation_ulmat | _hrr.csv | Q_COND | 8.51e-14 | 8.45e-14 | 1.01 |
| obst_activation_ulmat | _hrr.csv | Q_CONV | 0.534 | 0.534 | 1 |
| obst_activation_ulmat | _hrr.csv | Q_TOTAL | 0.534 | 0.534 | 1 |
| obst_activation_ulmat | _hrr.csv | Q_ENTH | 22 | 45.2 | 0.486 |
| divergence_test_3 | _devc.csv | div_max | 4.33e-11 | 1.1e-12 | 39.4 |
| divergence_test_3 | _devc.csv | div_min | 3.98e-11 | 1.13e-12 | 35.2 |
| divergence_test_3 | _hrr.csv | Q_CONV | 6.96e-14 | 2.97e-16 | 234 |
| divergence_test_3 | _hrr.csv | Q_TOTAL | 8.3e-14 | 1.31e-15 | 63.4 |
| divergence_test_3 | _hrr.csv | Q_COND | 3.16e-14 | 1.05e-15 | 30 |
| divergence_test_3 | _hrr.csv | Q_ENTH | 0.00437 | 0.0726 | 0.0602 |
| dancing_eddies_default | _devc.csv | error | 0.00046 | 0.000492 | 0.937 |
| dancing_eddies_default | _devc.csv | iter | 5 | 6 | 0.833 |
| dancing_eddies_default | _devc.csv | pres | 0.00387 | 0.401 | 0.00966 |
| dancing_eddies_default | _hrr.csv | Q_CONV | 7.67e-14 | 2.55e-16 | 301 |
| dancing_eddies_default | _hrr.csv | Q_TOTAL | 8.75e-14 | 1.39e-15 | 62.9 |
| dancing_eddies_default | _hrr.csv | Q_COND | 1.64e-14 | 1.15e-15 | 14.2 |
| dancing_eddies_default | _hrr.csv | Q_ENTH | 2.16e-13 | 3.12e-14 | 6.94 |
| lapse_rate | _hrr.csv | Q_CONV | 0.00591 | 0.0179 | 0.33 |
| lapse_rate | _hrr.csv | Q_COND | 6.13e-06 | 0.000127 | 0.0483 |
| lapse_rate | _hrr.csv | Q_TOTAL | 0.00591 | 0.998 | 0.00592 |
| lapse_rate | _hrr.csv | Q_ENTH | 3.58e-06 | 0.999 | 3.58e-06 |
| lapse_rate | _line.csv | T | 1.26e-08 | 20 | 6.31e-10 |
| lapse_rate | _line.csv | P | 0 | 1.01e+05 | 0 |
| soborot_superbee_square_wave_128 | _devc.csv | Y_TRACER-1 | 0 | 0 | 0 |
| soborot_superbee_square_wave_128 | _devc.csv | Y_TRACER-2 | 0 | 0 | 0 |
| soborot_superbee_square_wave_128 | _devc.csv | Y_TRACER-3 | 0 | 0 | 0 |
| soborot_superbee_square_wave_128 | _devc.csv | Y_TRACER-4 | 0 | 0 | 0 |
| soborot_superbee_square_wave_128 | _hrr.csv | HRR | 0 | 0 | 0 |
| soborot_superbee_square_wave_128 | _hrr.csv | HRR_OX | 0 | 0 | 0 |
| soborot_superbee_square_wave_128 | _hrr.csv | Q_RADI | 0 | 0 | 0 |
| soborot_superbee_square_wave_128 | _hrr.csv | Q_CONV | 0 | 0.004 | 0 |
| shunn3_4mesh_128 | _hrr.csv | Q_ENTH | 0.0318 | 0.0319 | 0.996 |
| shunn3_4mesh_128 | _hrr.csv | ZONE_1 | 1.15e-05 | 0.00397 | 0.00289 |
| shunn3_4mesh_128 | _hrr.csv | Q_PRES | 4.3e-07 | 0.000158 | 0.00272 |
| shunn3_4mesh_128 | _hrr.csv | Q_TOTAL | 4.3e-07 | 0.000158 | 0.00272 |
| shunn3_4mesh_128 | _mass.csv | SCALAR | 5.09e-06 | 0.2 | 2.54e-05 |
| shunn3_4mesh_128 | _mass.csv | LUMPED SCALAR | 5.09e-06 | 0.2 | 2.54e-05 |
| shunn3_4mesh_128 | _mass.csv | BACKGROUND | 2.54e-05 | 1 | 2.54e-05 |
| shunn3_4mesh_128 | _mass.csv | LUMPED BACKGROUND | 2.54e-05 | 1 | 2.54e-05 |

Reading:
- GLMAT and FFT agree where the physics is insensitive:
  - soborot: identical.
  - lapse_rate: T profile within 1.3e-8 abs.
  - shunn3_4mesh_128: mass within 2.5e-5 rel.
- Divergence-type outputs (D_max, div_max) differ at the 1e-12 to 1e-11 level. These are round-off quantities whose ‖ref‖∞ is itself only 1e-13 to 1e-12, which is why the relative numbers are large. The same holds for the Q_CONV, Q_COND and Q_TOTAL entries whose ‖ref‖∞ is about 1e-15.
- dancing_eddies: pres differs by 0.97 % of ‖ref‖∞ (3.9e-3 abs). The iteration count and error columns differ by construction, since they are solver diagnostics. The 4-mesh GLMAT result equals the A-24 single-mesh result to 1.3e-10, so this difference is entirely the multi-mesh FFT iteration error at VELOCITY_TOLERANCE.
- GLMAT wall time vs FFT:

| case | FFT | GLMAT |
|---|---|---|
| divergence_test_3 | 80 s | 135 s |
| dancing_eddies_default | 83 s | 118 s |
| shunn3_4mesh_128 | 8.2 s | 11.1 s |

## A-24: derived single-mesh references (requirements v0.3.2)

- Each derived copy covers the same domain and cell size, with all MESH lines (including MULT expansion) merged into one `&MESH` of the same total IJK.
- Nothing else is changed, including the CHID. The copy runs on 1 rank.
- Each copy is run as committed and as a SIG_FIGS=17 copy.
- The multi-mesh originals are **T2 references**. dancing_eddies_default was already captured in Tier 1. layer_1mesh is used as is.

| run | variant | ranks | meshes | wall s | exit/normal | reference class |
|---|---|---|---|---|---|---|
| dancing_eddies_default | as committed | 4 | 4 | 84.9 | 0/yes | T2 reference (multi-mesh) |
| dancing_eddies_default__1mesh | derived single-mesh copy (A-24): all MESH lines merged into one mesh,  | 1 | 1 | 262 | 0/yes | single-mesh reference (T2 class per v0.3.2) |
| dancing_eddies_default__1mesh_sf17 | derived single-mesh copy (A-24) + SIG_FIGS=17 | 1 | 1 | 272 | 0/yes | single-mesh reference (T2 class per v0.3.2) |
| dancing_eddies_default__sf17 | SIG_FIGS=17 copy | 4 | 4 | 82.7 | 0/yes | T2 reference (multi-mesh) |
| layer_1mesh | as committed (A-24 single-mesh reference used as is) | 1 | 1 | 94 | 0/yes | single-mesh reference (T2 class per v0.3.2) |
| layer_1mesh__sf17 | SIG_FIGS=17 copy | 1 | 1 | 89.8 | 0/yes | single-mesh reference (T2 class per v0.3.2) |
| shunn3_4mesh_32 | as committed (multi-mesh original, A-24) | 4 | 4 | 0.6 | 0/yes | T2 reference (multi-mesh) |
| shunn3_4mesh_32__1mesh | derived single-mesh copy (A-24): all MESH lines merged into one mesh,  | 1 | 1 | 1.2 | 0/yes | single-mesh reference (T2 class per v0.3.2) |
| shunn3_4mesh_32__1mesh_sf17 | derived single-mesh copy (A-24) + SIG_FIGS=17 | 1 | 1 | 1.19 | 0/yes | single-mesh reference (T2 class per v0.3.2) |
| shunn3_4mesh_32__sf17 | SIG_FIGS=17 copy (multi-mesh original, A-24) | 4 | 4 | 0.59 | 0/yes | T2 reference (multi-mesh) |
| symmetry_test_mpi | as committed (multi-mesh original, A-24) | 8 | 8 | 36.1 | 0/yes | T2 reference (multi-mesh) |
| symmetry_test_mpi__1mesh | derived single-mesh copy (A-24): all MESH lines merged into one mesh,  | 1 | 1 | 15.4 | 0/yes | single-mesh reference (T2 class per v0.3.2) |
| symmetry_test_mpi__1mesh_sf17 | derived single-mesh copy (A-24) + SIG_FIGS=17 | 1 | 1 | 15 | 0/yes | single-mesh reference (T2 class per v0.3.2) |
| symmetry_test_mpi__sf17 | SIG_FIGS=17 copy (multi-mesh original, A-24) | 8 | 8 | 42.2 | 0/yes | T2 reference (multi-mesh) |

Exact edits (also in each run folder as `derivation.diff` and `manifest.json` → derivation_diff):

`shunn3_4mesh_32__1mesh`:
```diff
--- Verification/Scalar_Analytical_Solution/shunn3_4mesh_32.fds
+++ shunn3_4mesh_32.fds
@@ -2,10 +2,7 @@

 MESH IJK=32,1,32, XB=-1,1,-0.05,0.05,-1,1/

-&MESH IJK=16,1,16, XB=-1,0,-0.05,0.05,-1,0/
-&MESH IJK=16,1,16, XB= 0,1,-0.05,0.05,-1,0/
-&MESH IJK=16,1,16, XB= 0,1,-0.05,0.05, 0,1/
-&MESH IJK=16,1,16, XB=-1,0,-0.05,0.05, 0,1/
+&MESH IJK=32,1,32, XB=-1,1,-0.05,0.05,-1,1 /

 &TIME T_END=1/
```

`dancing_eddies_default__1mesh`:
```diff
--- Verification/Pressure_Solver/dancing_eddies_default.fds
+++ dancing_eddies_default.fds
@@ -1,7 +1,6 @@
 &HEAD CHID='dancing_eddies_default', TITLE='Simple 2D Tunnel'/

-&MULT ID='mesh',DX=0.075,I_UPPER=3 /
-&MESH IJK=75,1,40, XB=0.0,0.075,-0.001,0.001,-0.02,0.020, MULT_ID='mesh' /
+&MESH IJK=300,1,40, XB=0.0,0.3,-0.001,0.001,-0.02,0.020 /

 &TIME T_END=2.0 /
```

`symmetry_test_mpi__1mesh`:
```diff
--- Verification/Flowfields/symmetry_test_mpi.fds
+++ symmetry_test_mpi.fds
@@ -1,7 +1,6 @@
 &HEAD CHID='symmetry_test_mpi', TITLE='Test symmetry of flow solver, MPI version' /

-&MESH IJK=10,10,10, XB=0.0,0.5,0.0,0.5,0.0,0.5, MULT_ID='mult' /
-&MULT ID='mult', DX=0.5, DY=0.5, DZ=0.5, I_UPPER=1, J_UPPER=1, K_UPPER=1 /
+&MESH IJK=20,20,20, XB=0.0,1.0,0.0,1.0,0.0,1.0 /

 &TIME T_END=500. /
```

Multi-mesh (reference) vs single-mesh copy, SIG_FIGS=17 pair, interpolated onto the multi-mesh output times. Every non-timing column; "rel∞" = max|Δ|/‖multi-mesh‖∞, "pt rel" = pointwise max |Δ|/|ref| over points with |ref| ≥ 1e-3·‖ref‖∞:

**A24 shunn3_4mesh_32: multi-mesh vs single-mesh (sf17)**, columns identically zero in both runs omitted (formal T1 on raw rows: _hrr.csv FAIL, _mass.csv FAIL)

| file | quantity | max abs | ‖ref‖∞ | rel∞ | pt rel | t at max | \|Δ\| at end |
|---|---|---|---|---|---|---|---|
| _hrr.csv | Q_PRES | 4.78e-05 | 0.0011 | 0.0435 | 0.781 | 0.9742 | 1.15e-06 |
| _hrr.csv | Q_ENTH | 0.0373 | 0.0857 | 0.435 | 9.24 | 0.529 | 0.000232 |
| _hrr.csv | Q_TOTAL | 4.78e-05 | 0.0011 | 0.0435 | 0.781 | 0.9742 | 1.15e-06 |
| _hrr.csv | ZONE_1 | 0.000218 | 0.0098 | 0.0222 | 0.347 | 0.9541 | 0.000149 |
| _mass.csv | Total | 2.93e-05 | 1.2 | 2.44e-05 | 2.44e-05 | 0.4592 | 9.43e-06 |
| _mass.csv | BACKGROUND | 3.67e-05 | 1 | 3.67e-05 | 3.67e-05 | 0.4592 | 1.18e-05 |
| _mass.csv | SCALAR | 7.33e-06 | 0.2 | 3.67e-05 | 3.67e-05 | 0.4592 | 2.36e-06 |
| _mass.csv | LUMPED BACKGROUND | 3.67e-05 | 1 | 3.67e-05 | 3.67e-05 | 0.4592 | 1.18e-05 |
| _mass.csv | LUMPED SCALAR | 7.33e-06 | 0.2 | 3.67e-05 | 3.67e-05 | 0.4592 | 2.36e-06 |

**A24 shunn3_4mesh_32: multi-mesh vs single-mesh (committed)**, columns identically zero in both runs omitted (formal T1 on raw rows: _hrr.csv FAIL, _mass.csv FAIL)

| file | quantity | max abs | ‖ref‖∞ | rel∞ | pt rel | t at max | \|Δ\| at end |
|---|---|---|---|---|---|---|---|
| _hrr.csv | Q_PRES | 4.78e-05 | 0.0011 | 0.0435 | 0.781 | 0.9742 | 1.15e-06 |
| _hrr.csv | Q_ENTH | 0.0373 | 0.0857 | 0.435 | 9.24 | 0.529 | 0.000232 |
| _hrr.csv | Q_TOTAL | 4.78e-05 | 0.0011 | 0.0435 | 0.781 | 0.9742 | 1.15e-06 |
| _hrr.csv | ZONE_1 | 0.000218 | 0.0098 | 0.0222 | 0.347 | 0.9541 | 0.000149 |
| _mass.csv | Total | 2.93e-05 | 1.2 | 2.44e-05 | 2.44e-05 | 0.4592 | 9.4e-06 |
| _mass.csv | BACKGROUND | 3.67e-05 | 1 | 3.67e-05 | 3.67e-05 | 0.4592 | 1.18e-05 |
| _mass.csv | SCALAR | 7.34e-06 | 0.2 | 3.67e-05 | 3.67e-05 | 0.4592 | 2.35e-06 |
| _mass.csv | LUMPED BACKGROUND | 3.67e-05 | 1 | 3.67e-05 | 3.67e-05 | 0.4592 | 1.18e-05 |
| _mass.csv | LUMPED SCALAR | 7.34e-06 | 0.2 | 3.67e-05 | 3.67e-05 | 0.4592 | 2.35e-06 |

**A24 dancing_eddies_default: multi-mesh vs single-mesh (sf17)**, columns identically zero in both runs omitted (formal T1 on raw rows: _devc.csv FAIL, _hrr.csv FAIL)

| file | quantity | max abs | ‖ref‖∞ | rel∞ | pt rel | t at max | \|Δ\| at end |
|---|---|---|---|---|---|---|---|
| _devc.csv | pres | 0.00387 | 0.401 | 0.00966 | 0.594 | 1.82 | 0.000174 |
| _devc.csv | error | 0.00046 | 0.000492 | 0.937 | 0.976 | 0.8304 | 0.000138 |
| _devc.csv | iter | 5 | 6 | 0.833 | 0.833 | 1.316 | 0.853 |
| _hrr.csv | Q_CONV | 2.85e-16 | 2.55e-16 | 1.12 | 7.95 | 1.598 | 4.22e-17 |
| _hrr.csv | Q_COND | 2.21e-16 | 1.15e-15 | 0.192 | 12.2 | 0.01197 | 7.77e-17 |
| _hrr.csv | Q_ENTH | 2.94e-14 | 3.12e-14 | 0.944 | 9 | 0.09865 | 0 |
| _hrr.csv | Q_TOTAL | 2.31e-16 | 1.39e-15 | 0.166 | 9.86 | 1.598 | 3.55e-17 |
| _hrr.csv | MLR_LJ AIR | 2.91e-06 | 0.0231 | 0.000126 | 0.00017 | 0.8945 | 1.62e-10 |

**A24 dancing_eddies_default: multi-mesh vs single-mesh (committed)**, columns identically zero in both runs omitted (formal T1 on raw rows: _devc.csv FAIL, _hrr.csv FAIL)

| file | quantity | max abs | ‖ref‖∞ | rel∞ | pt rel | t at max | \|Δ\| at end |
|---|---|---|---|---|---|---|---|
| _devc.csv | pres | 0.00387 | 0.401 | 0.00966 | 0.595 | 1.82 | 0.000174 |
| _devc.csv | error | 0.00046 | 0.000492 | 0.937 | 0.976 | 0.8304 | 0.000138 |
| _devc.csv | iter | 5 | 6 | 0.833 | 0.833 | 1.316 | 0.853 |
| _hrr.csv | Q_CONV | 2.85e-16 | 2.55e-16 | 1.12 | 7.95 | 1.598 | 4.22e-17 |
| _hrr.csv | Q_COND | 2.21e-16 | 1.15e-15 | 0.192 | 12.2 | 0.01197 | 7.77e-17 |
| _hrr.csv | Q_ENTH | 2.94e-14 | 3.12e-14 | 0.944 | 9 | 0.09865 | 0 |
| _hrr.csv | Q_TOTAL | 2.31e-16 | 1.39e-15 | 0.166 | 9.86 | 1.598 | 3.55e-17 |
| _hrr.csv | MLR_LJ AIR | 2.91e-06 | 0.0231 | 0.000126 | 0.00017 | 0.8945 | 0 |

**A24 symmetry_test_mpi: multi-mesh vs single-mesh (sf17)**, columns identically zero in both runs omitted (formal T1 on raw rows: _devc.csv FAIL, _hrr.csv FAIL)

| file | quantity | max abs | ‖ref‖∞ | rel∞ | pt rel | t at max | \|Δ\| at end |
|---|---|---|---|---|---|---|---|
| _devc.csv | u_1 | 8.62e-06 | 0.000771 | 0.0112 | 0.156 | 9.5 | 7.29e-07 |
| _devc.csv | u_2 | 1.3e-05 | 0.000826 | 0.0157 | 0.197 | 1 | 5.02e-07 |
| _devc.csv | v_1 | 8.62e-06 | 0.000771 | 0.0112 | 0.156 | 9.5 | 7.29e-07 |
| _devc.csv | v_2 | 1.3e-05 | 0.000826 | 0.0157 | 0.197 | 1 | 5.02e-07 |
| _devc.csv | w_1 | 8.62e-06 | 0.000771 | 0.0112 | 0.156 | 9.5 | 7.29e-07 |
| _devc.csv | w_2 | 1.3e-05 | 0.000826 | 0.0157 | 0.197 | 1 | 5.02e-07 |
| _hrr.csv | Q_CONV | 5.51e-05 | 0.0171 | 0.00323 | 0.00757 | 4.5 | 5.24e-09 |
| _hrr.csv | Q_PRES | 1.4e-07 | 0.397 | 3.52e-07 | 3.58e-07 | 490 | 5.05e-08 |
| _hrr.csv | Q_ENTH | 1.1e-05 | 0.38 | 2.88e-05 | 6.67e-05 | 6.5 | 1.08e-06 |
| _hrr.csv | Q_TOTAL | 5.51e-05 | 0.38 | 0.000145 | 0.000337 | 4.5 | 5.57e-08 |
| _hrr.csv | MLR_AIR | 5.27e-11 | 0.00336 | 1.57e-08 | 1.57e-08 | 500 | 5.27e-11 |
| _hrr.csv | ZONE_1 | 0.00373 | 1.36e+05 | 2.75e-08 | 2.75e-08 | 500 | 0.00373 |

**A24 symmetry_test_mpi: multi-mesh vs single-mesh (committed)**, columns identically zero in both runs omitted (formal T1 on raw rows: _devc.csv FAIL, _hrr.csv FAIL)

| file | quantity | max abs | ‖ref‖∞ | rel∞ | pt rel | t at max | \|Δ\| at end |
|---|---|---|---|---|---|---|---|
| _devc.csv | u_1 | 8.62e-06 | 0.000771 | 0.0112 | 0.156 | 9.5 | 7.29e-07 |
| _devc.csv | u_2 | 1.3e-05 | 0.000826 | 0.0157 | 0.197 | 1 | 5.02e-07 |
| _devc.csv | v_1 | 8.62e-06 | 0.000771 | 0.0112 | 0.156 | 9.5 | 7.29e-07 |
| _devc.csv | v_2 | 1.3e-05 | 0.000826 | 0.0157 | 0.197 | 1 | 5.02e-07 |
| _devc.csv | w_1 | 8.62e-06 | 0.000771 | 0.0112 | 0.156 | 9.5 | 7.29e-07 |
| _devc.csv | w_2 | 1.3e-05 | 0.000826 | 0.0157 | 0.197 | 1 | 5.02e-07 |
| _hrr.csv | Q_CONV | 5.51e-05 | 0.0171 | 0.00323 | 0.00757 | 4.5 | 5e-09 |
| _hrr.csv | Q_PRES | 1.4e-07 | 0.397 | 3.52e-07 | 3.58e-07 | 490 | 5e-08 |
| _hrr.csv | Q_ENTH | 1.1e-05 | 0.38 | 2.88e-05 | 6.68e-05 | 6.5 | 1.09e-06 |
| _hrr.csv | Q_TOTAL | 5.51e-05 | 0.38 | 0.000145 | 0.000337 | 4.5 | 6e-08 |
| _hrr.csv | MLR_AIR | 1e-10 | 0.00336 | 2.98e-08 | 6.92e-08 | 5 | 1e-10 |
| _hrr.csv | ZONE_1 | 0.01 | 1.36e+05 | 7.36e-08 | 9.89e-08 | 423 | 0 |

MMS error on the lower-left quadrant (the only mesh the multi-mesh `_mms.csv` holds):

| run | region | MMS time | e_rho | e_Z | e_u | e_H |
|---|---|---|---|---|---|---|
| shunn3_4mesh_128 | multi-mesh run: _mms.csv holds mesh 1 only = lower-left quadrant (x,z in [-1,0]) | 1.000 | 0.00291 | 0.00144 | 0.00054 | 0.00162 |
| shunn3_4mesh_32 | multi-mesh run: _mms.csv holds mesh 1 only = lower-left quadrant (x,z in [-1,0]) | 1.000 | 0.0322 | 0.0137 | 0.00508 | 0.00905 |
| shunn3_128 (quadrant) | lower-left quadrant (x,z in [-1,0]) of the single-mesh field | 0.904 | 0.0022 | 0.000599 | 0.000301 | 0.00506 |
| shunn3_4mesh_32__1mesh (quadrant) | lower-left quadrant (x,z in [-1,0]) of the single-mesh field | 1.000 | 0.0323 | 0.0137 | 0.00505 | 0.0101 |

shunn3_128 writes its MMS slice at t = 0.904 (MMS_TIMER=0.9), while the 4-mesh cases write at t = 1.0. Its row is therefore not directly comparable with shunn3_4mesh_128. The A-24 comparison is shunn3_4mesh_32 against shunn3_4mesh_32__1mesh (quadrant), both at t = 1.0.

**A-24 summary (multi-mesh T2 reference vs single-mesh copy, per quantity, SIG_FIGS=17):**
- Committed and sf17 pairs give the same differences to 3 digits.
- For shunn3_4mesh_32 and dancing_eddies_default, the multi-mesh and single-mesh time stamps differ (dt comes from different reductions), so the formal row-matched T1 fails on structure. For symmetry_test_mpi the time base matches and T1 fails on the values. All numbers are interpolated onto the multi-mesh times.
- The single-mesh copies are also 1-rank runs, so they carry no MPI effect.

**shunn3_4mesh_32** (4 meshes / 4 ranks vs 1 mesh; 32×32, T_END = 1):
- Mass: Total max |Δ| 2.9e-5 (2.4e-5 of ‖ref‖∞). BACKGROUND 3.7e-5 abs (3.7e-5 rel). SCALAR 7.3e-6 abs (3.7e-5 rel).
- HRR file:
  - Q_PRES and Q_TOTAL: 4.8e-5 abs (4.35 % of ‖ref‖∞ = 1.1e-3).
  - Q_ENTH: 3.7e-2 abs (43.5 % of 8.6e-2; a small net term).
  - ZONE_1 (background pressure): 2.2e-4 Pa (2.2 % of 9.8e-3).
- MMS L2 error on the lower-left quadrant, 4-mesh vs 1-mesh:

| quantity | 4-mesh | 1-mesh | difference |
|---|---|---|---|
| e_rho | 3.2249e-2 | 3.2255e-2 | 0.02 % |
| e_Z | 1.37092e-2 | 1.37094e-2 | 1e-5 relative |
| e_u | 5.08e-3 | 5.05e-3 | 0.6 % |
| e_H | 9.05e-3 | 1.008e-2 | 10 % |

  The multi-mesh `_mms.csv` only holds mesh 1, so a full-domain 4-mesh error is not available from FDS output.

**dancing_eddies_default** (4 meshes / 4 ranks vs 1 mesh; 300×40, T_END = 2 s):
- devc `pres`: max |Δ| 3.87e-3 Pa (0.97 % of ‖ref‖∞ = 0.401), at t = 1.82 s.
- `error` (the Poisson velocity-error diagnostic) and `iter` differ by construction: the multi-mesh run iterates to VELOCITY_TOLERANCE, and the single mesh gets a direct FFT solution.
- HRR file: MLR 2.9e-6 (1.3e-4 of ‖ref‖∞). The energy terms are round-off, ≤ 3e-14 abs.
- Cross-checks:
  - The single-mesh copy agrees with the 4-mesh GLMAT copy (a global direct solve) to 1.3e-10 in `pres`. Merging the meshes and removing the inter-mesh iteration error are equivalent.
  - The derived copy agrees with FDS's committed `dancing_eddies_1mesh` to 1.5e-11. The committed file also sets TUNNEL_PRECONDITIONER=F, which makes no material difference on 1 mesh. The committed 1mesh case is therefore not a verbatim A-24 derivation, but it is numerically equivalent.

**symmetry_test_mpi** (8 meshes / 8 ranks vs 1 mesh; 20³, T_END = 500 s):
- devc point velocities u/v/w at 2 points: max |Δ| 8.6e-6 m/s (1.1 % of ‖ref‖∞ = 7.7e-4) and 1.3e-5 m/s (1.6 % of 8.3e-4).
- In both runs, u = v = w to ≤ 1.3e-18, so the symmetry is preserved by both decompositions.
- HRR file:
  - Q_CONV: 5.5e-5 (0.32 %).
  - Q_TOTAL: 1.45e-4 of ‖ref‖∞.
  - Q_ENTH: 2.9e-5.
  - Q_PRES: 3.5e-7.
  - ZONE_1 pressure: 2.8e-8 relative.

**layer_1mesh**: captured as is, committed and sf17 (94 s and 90 s). Its multi-mesh partner layer_4mesh is not in A-24 scope, so no difference is reported.

**Implication for T2 tolerances:**
- Mesh decomposition alone moves "smooth" integrated quantities by 1e-5 to 1e-4 relative (mass, pressure level).
- It moves local point values in weak flows by about 1 % of their scale (quiescent symmetry case, dancing eddies `pres`).
- It moves near-zero balance terms (Q_ENTH, Q_PRES in shunn) by up to 40 % of their small scale.
- T2 tolerances for FDS-AMR-vs-FDS on multi-mesh references should not be tighter than these numbers. They should use an absolute floor tied to each column's scale, not pointwise relative error.

## As-committed vs SIG_FIGS=17 sanity

Every run with a sf17 copy was checked. SIG_FIGS should only affect output formatting, so the step sequence must be identical, with the same row count, and every value must agree within 8-significant-figure rounding (pointwise relative ≤ 5e-8).

88 file pairs, 0 not PASS.

| case | file | status | rows | max pointwise rel | max abs/col ‖‖∞ | worst column |
|---|---|---|---|---|---|---|
| 1_step_2_step_compare | _devc.csv | PASS | 21 | 4.91e-08 | 3.05e-08 | Q4S |
| 1_step_2_step_compare | _hrr.csv | PASS | 21 | 3.75e-08 | 2.52e-08 | Q_PRES |
| dancing_eddies_1mesh | _devc.csv | PASS | 1001 | 4.85e-08 | 2.5e-08 | Time |
| dancing_eddies_1mesh | _hrr.csv | PASS | 1001 | 4.76e-08 | 3.67e-08 | Q_CONV |
| dancing_eddies_default | _devc.csv | PASS | 1001 | 4.77e-08 | 2.49e-08 | Time |
| dancing_eddies_default | _hrr.csv | PASS | 1001 | 4.92e-08 | 4.34e-08 | Q_COND |
| dancing_eddies_uglmat_refine | _devc.csv | PASS | 1001 | 4.66e-08 | 2.5e-08 | Time |
| dancing_eddies_uglmat_refine | _hrr.csv | PASS | 1001 | 4.87e-08 | 3.55e-08 | Q_COND |
| divergence_test_2 | _devc.csv | PASS | 58 | 3.64e-08 | 1.58e-08 | Vdot |
| divergence_test_2 | _hrr.csv | PASS | 58 | 4.71e-08 | 2.74e-08 | Q_RADI |
| divergence_test_3 | _devc.csv | PASS | 1001 | 4.78e-08 | 3.28e-08 | div_min |
| divergence_test_3 | _hrr.csv | PASS | 1001 | 4.98e-08 | 4.75e-08 | Q_COND |
| energy_budget_tmix | _devc.csv | PASS | 1001 | 4.91e-08 | 1.3e-08 | Temp |
| energy_budget_tmix | _hrr.csv | PASS | 1001 | 4.71e-08 | 2.33e-08 | Time |
| lapse_rate | _hrr.csv | PASS | 5 | 3.65e-08 | 3.65e-08 | Q_COND |
| layer_1mesh | _devc.csv | PASS | 31 | 2.83e-08 | 2.19e-08 | z_int |
| layer_1mesh | _hrr.csv | PASS | 1001 | 4.88e-08 | 4.68e-08 | Q_RADI |
| ns2d_16 | _devc.csv | PASS | 296 | 4.98e-08 | 2.42e-08 | PRES |
| ns2d_16 | _hrr.csv | PASS | 296 | 4.42e-08 | 7.95e-09 | Time |
| ns2d_16_emb_1to2_refinement | _devc.csv | PASS | 1001 | 4.7e-08 | 2.43e-08 | PRES |
| ns2d_16_emb_1to2_refinement | _hrr.csv | PASS | 1001 | 4.58e-08 | 1.66e-08 | Time |
| ns2d_16_int_1to2_refinement | _devc.csv | PASS | 1001 | 4.79e-08 | 1.66e-08 | Time |
| ns2d_16_int_1to2_refinement | _hrr.csv | PASS | 1001 | 4.63e-08 | 3.63e-08 | Q_ENTH |
| ns2d_16_int_1to2_refinement_uglmat | _devc.csv | PASS | 1001 | 4.9e-08 | 1.67e-08 | Time |
| ns2d_16_int_1to2_refinement_uglmat | _hrr.csv | PASS | 1001 | 4.77e-08 | 2.6e-08 | Q_PRES |
| ns2d_16_nupt1 | _devc.csv | PASS | 228 | 4.71e-08 | 2.55e-08 | PRES |
| ns2d_16_nupt1 | _hrr.csv | PASS | 228 | 3.88e-08 | 7.9e-09 | Time |
| ns2d_32 | _devc.csv | PASS | 569 | 4.75e-08 | 2.46e-08 | PRES |
| ns2d_32 | _hrr.csv | PASS | 569 | 4.4e-08 | 2.91e-08 | Q_PRES |
| ns2d_32_nupt1 | _devc.csv | PASS | 441 | 4.98e-08 | 2.56e-08 | UVEL |
| ns2d_32_nupt1 | _hrr.csv | PASS | 441 | 4.85e-08 | 2.91e-08 | Q_PRES |
| ns2d_64 | _devc.csv | PASS | 1001 | 4.89e-08 | 2.48e-08 | PRES |
| ns2d_64 | _hrr.csv | PASS | 1001 | 4.67e-08 | 7.93e-09 | Time |
| ns2d_64_nupt1 | _devc.csv | PASS | 937 | 4.83e-08 | 2.62e-08 | UVEL |
| ns2d_64_nupt1 | _hrr.csv | PASS | 937 | 4.29e-08 | 7.96e-09 | Time |
| ns2d_8 | _devc.csv | PASS | 157 | 4.88e-08 | 2.43e-08 | PRES |
| ns2d_8 | _hrr.csv | PASS | 157 | 4.05e-08 | 7.92e-09 | Time |
| ns2d_8_nupt1 | _devc.csv | PASS | 121 | 4.79e-08 | 3.26e-08 | PRES |
| ns2d_8_nupt1 | _hrr.csv | PASS | 121 | 4.74e-08 | 7.87e-09 | Time |
| obst_activation_default | _devc.csv | PASS | 145 | 4.7e-08 | 4.28e-08 | D_min |
| obst_activation_default | _hrr.csv | PASS | 145 | 4.75e-08 | 2.72e-08 | Q_COND |
| obst_activation_ulmat | _devc.csv | PASS | 145 | 4.83e-08 | 2.67e-08 | D_max |
| obst_activation_ulmat | _hrr.csv | PASS | 145 | 4.64e-08 | 2.5e-08 | Time |
| obst_coarse_fine_interface | _devc.csv | PASS | 44 | 4.75e-08 | 4.16e-08 | P2 |
| obst_coarse_fine_interface | _hrr.csv | PASS | 44 | 4.55e-08 | 4.37e-08 | Q_COND |
| random_meshes | _devc.csv | PASS | 1001 | 4.88e-08 | 4.12e-08 | Vdot |
| random_meshes | _hrr.csv | PASS | 1001 | 4.94e-08 | 2.11e-08 | Q_ENTH |
| restart_test1_continuous | _devc.csv | PASS | 401 | 4.91e-08 | 3.72e-08 | null#17 |
| restart_test1_continuous | _hrr.csv | PASS | 401 | 4.97e-08 | 4.2e-08 | Q_DIFF |
| restart_test1a | _devc.csv | PASS | 345 | 4.94e-08 | 4.19e-08 | null#14 |
| restart_test1a | _hrr.csv | PASS | 345 | 4.79e-08 | 3.68e-08 | Q_COND |
| restart_test1b | _devc.csv | PASS | 727 | 4.94e-08 | 4.25e-08 | null#17 |
| restart_test1b | _hrr.csv | PASS | 727 | 4.79e-08 | 3.88e-08 | Q_DIFF |
| saad_512_cfl_1 | _devc.csv | PASS | 101 | 4.99e-08 | 4.98e-08 | u_max |
| saad_512_cfl_1 | _hrr.csv | PASS | 101 | 4.31e-08 | 1.28e-08 | Q_PRES |
| saad_512_cfl_1 | _mass.csv | PASS | 101 | 8.75e-09 | 8.75e-09 | SCALAR |
| saad_512_cfl_p0625 | _devc.csv | PASS | 1001 | 4.83e-08 | 1.36e-08 | p0_min |
| saad_512_cfl_p0625 | _hrr.csv | PASS | 1001 | 4.87e-08 | 3.46e-08 | Q_ENTH |
| saad_512_cfl_p0625 | _mass.csv | PASS | 1001 | 4.83e-08 | 1.28e-08 | Time |
| saad_512_cfl_p125 | _devc.csv | PASS | 801 | 4.83e-08 | 1.36e-08 | p0_min |
| saad_512_cfl_p125 | _hrr.csv | PASS | 801 | 4.83e-08 | 1.28e-08 | Time |
| saad_512_cfl_p125 | _mass.csv | PASS | 801 | 4.83e-08 | 1.28e-08 | Time |
| saad_512_cfl_p25 | _devc.csv | PASS | 401 | 4.83e-08 | 1.36e-08 | p0_min |
| saad_512_cfl_p25 | _hrr.csv | PASS | 401 | 4.83e-08 | 2.01e-08 | Q_ENTH |
| saad_512_cfl_p25 | _mass.csv | PASS | 401 | 4.83e-08 | 1.28e-08 | Time |
| saad_512_cfl_p5 | _devc.csv | PASS | 201 | 4.83e-08 | 1.36e-08 | p0_min |
| saad_512_cfl_p5 | _hrr.csv | PASS | 201 | 4.83e-08 | 2.89e-08 | Q_ENTH |
| saad_512_cfl_p5 | _mass.csv | PASS | 201 | 4.83e-08 | 1.28e-08 | Time |
| shunn3_128 | _hrr.csv | PASS | 187 | 4.58e-08 | 2.83e-08 | Q_PRES |
| shunn3_128 | _mass.csv | PASS | 187 | 3.73e-08 | 4.99e-09 | Time |
| shunn3_32 | _hrr.csv | PASS | 80 | 4.42e-08 | 3.71e-08 | Q_PRES |
| shunn3_32 | _mass.csv | PASS | 80 | 3.59e-08 | 2.49e-08 | SCALAR |
| shunn3_4mesh_128 | _hrr.csv | PASS | 186 | 4.55e-08 | 2.83e-08 | Q_PRES |
| shunn3_4mesh_128 | _mass.csv | PASS | 186 | 4.99e-08 | 4.99e-08 | BACKGROUND |
| shunn3_4mesh_32 | _hrr.csv | PASS | 81 | 4.21e-08 | 3.71e-08 | Q_PRES |
| shunn3_4mesh_32 | _mass.csv | PASS | 81 | 4.93e-08 | 4.93e-08 | BACKGROUND |
| shunn3_64 | _hrr.csv | PASS | 100 | 4.68e-08 | 2.28e-08 | ZONE_1 |
| shunn3_64 | _mass.csv | PASS | 100 | 3.74e-08 | 5e-09 | Time |
| soborot_superbee_square_wave_128 | _devc.csv | PASS | 11 | 4.61e-08 | 3.4e-08 | Y_TRACER-20 |
| soborot_superbee_square_wave_128 | _hrr.csv | PASS | 11 | 3.77e-08 | 2.15e-08 | Q_COND |
| soborot_superbee_square_wave_128_1mesh | _devc.csv | PASS | 11 | 4.61e-08 | 3.4e-08 | Y_TRACER-20 |
| soborot_superbee_square_wave_128_1mesh | _hrr.csv | PASS | 11 | 4.69e-08 | 2.15e-08 | Q_COND |
| species_conservation_1 | _hrr.csv | PASS | 547 | 4.73e-08 | 2.97e-08 | Q_TOTAL |
| species_conservation_1 | _mass.csv | PASS | 61 | 3.42e-08 | 1.38e-08 | PROPANE |
| species_conservation_2 | _devc.csv | PASS | 188 | 4.34e-08 | 2e-08 | Time |
| species_conservation_2 | _hrr.csv | PASS | 188 | 4.92e-08 | 4.29e-08 | ZONE_2 |
| symmetry_test_mpi | _devc.csv | PASS | 1001 | 4.7e-08 | 6.44e-09 | u_1 |
| symmetry_test_mpi | _hrr.csv | PASS | 1001 | 4.72e-08 | 3.68e-08 | ZONE_1 |

## Findings and failures

1. **obst_activation_ulmat: FDS criterion FAIL.**
   - D_max abs error 1.77e-12 against tolerance 1e-13, with max divergence about 1.9e-12.
   - The default-solver sibling passes at 3.7e-14, and the GLMAT copies of both land at about 1e-12.
   - The reference build itself does not meet this FDS criterion. Carry it as a known reference deviation until it is checked against a firebot or Intel build. Do not hold FDS-AMR to 1e-13 here.
2. **shunn3 V&V order gate (≥ 1.8, 32→128) fails for Z (1.69) and H (0.92).** ρ (1.95) and u (2.00) pass. The FDS script limits only apply at N=512 (Tier 2). This is a gate calibration issue: re-calibrate to ≥ 1.6 for Z and drop H, or move to 64→256.
3. **species_conservation_1** passes with only a 9 % margin: 9.12e-3 against 1e-2.
4. **A-09:**
   - The interior 1:2 refinement (int_1to2) gives about 2.3× the u error of the uniform coarse grid (0.298 vs 0.131 to 2π). It is also 3× slower than UGLMAT because of 36 pressure iterations per step.
   - The default and UGLMAT solutions diverge chaotically after about 11 s.
   - These are expected prototype defects, now quantified.
5. **Restart vs continuous is not T0/T1-comparable with this case pair**, because of the T_END clip and the time-averaged DEVC windows. Restart reproducibility is bitwise.
6. **Rank-count changes break bitwise identity** at the 1e-16 to 1e-13 relative level, as expected. T1 as currently specified needs an absolute floor for near-zero columns (the soborot Q_PRES case).
7. No run failed, timed out or trapped: 113 manifests, all exit 0 with normal stop, and 0 trap lines, including all 4 Debug runs.

## Re-tier suggestions (measured, as-committed, single run)

| case | ranks | measured s | estimate s | ratio |
|---|---|---|---|---|
| dancing_eddies_1mesh | 1 | 285 | 46 | 6.2× |
| dancing_eddies_uglmat_refine | 4 | 99 | 18 | 5.5× |
| ns2d_16_int_1to2_refinement | 1 | 91 | 31 | 2.9× |
| dancing_eddies_default | 4 | 85 | 18 | 4.7× |
| divergence_test_3 | 4 | 85 | 13 | 6.5× |
| random_meshes | 4 | 59 | 7 | 8.5× |
| restart_test1_continuous / 1b / 1a | 1 | 57 / 44 / 37 | 17 / 17 / 10 | 3–4× |
| 1_step_2_step_compare | 4 | 51 | 41 | 1.2× |
| shunn3_128 | 1 | 40 | 20 | 2.0× |
| soborot_superbee_square_wave_128_1mesh | 1 | 38 | 10 | 3.8× |
| saad_512_cfl_p0625 | 1 | 30 | 10 | 3.0× |
| ns2d_64 / ns2d_64_nupt1 | 1 | 29 / 24 | 9 / 9 | 3× |

Proposal:
- The as-committed set takes 1193 s (19.9 min) serial, or 2411 rank-seconds. That does not fit the ~8 min Tier 1 budget.
- The 26 cases with ≤ 30 s each take 193 s together. The 34 cases with ≤ 60 s take 549 s.
- Move dancing_eddies_1mesh (285 s) to Tier 2. It is also numerically equivalent to the A-24 derived single-mesh copy, to 1.5e-11.
- Consider moving dancing_eddies_uglmat_refine, divergence_test_3 and random_meshes to Tier 2 as well, or run Tier 1 with 2–3 cases in parallel within the 8-rank cap. Either way Tier 1 comes back to about 8–10 min.
- The estimates in `case_inventory.*` should be regenerated from `docs/vv/measured_runtimes.csv`. The inventory is not edited here, per instruction.

## Deviations from the run protocol (all recorded per run in `manifest.json` → `preflight`)

- **Load gate.** The rule "load1 < 2 before launch" was implemented as: pass if load1 < 2, **or** if the instantaneous busy-core count, sampled from /proc/stat over 3 s, is < 1.0. There was a re-check every 30 s and a 30-min cap; the cap was never reached.
  - 80 of 113 gated launches (G0, A-09, Tier 1, A-24 and calibration) passed on load1 < 2. **33 launched with load1 ≥ 2 on the busy-core fallback.**
  - The 1-min load average stayed ≥ 2 for long stretches because of other short build bursts, including an unrelated `vv-runs/A-38` FDS job that ran for more than 40 min on the same binary. Without the fallback the capture would have stalled.
  - Total gate waiting time was about 60 min.
- **RAM gate.** It uses MemAvailable, not MemFree. MemFree was below ranks×0.5 + 1 GB for 99 of 113 launches (the page cache was full), while MemAvailable was always ≥ 4.5 GB. Both values are recorded in every manifest.
- **Restart files are written in every run.** In this FireX build DT_RESTART defaults to T_END/20, so every run wrote `.restart` files.
- **The "drop binaries > 5 MB" rule is applied per file.** Many runs keep dozens of small `.sf`/`.restart`/`.s3d` files; int_1to2, for example, keeps 117 `.sf` and 13 `.restart` files, about 28 MB. The frozen baseline directory is 1.9 GB. The dropped files are listed below with size and sha256.
- **restart_test1a large binaries** were kept until restart_test1b had copied them. They were then dropped (`prune_later`) and the 1a manifest was updated.
- **restart_test1_continuous** (inventory "VV/") is a derived copy of 1a (T_END=10, DT_RESTART removed, CHID changed). The diff is in its folder.
- **Serial execution.** All runs were strictly sequential, which is below the 8-concurrent-rank cap. Timeouts were max(600 s, 3× estimate); none triggered.
- **Wall times** include `mpirun` start-up. "fds_elapsed_s" in the CSV is FDS's own elapsed time.

## Pending / not done

- G0 setup-only sweep (T_END=0) over all 941 inputs: not run, pending per instruction.
- Intel oneAPI toolchain baseline: not started, since this capture covers the GNU/Open MPI reference only.
- A full-domain MMS error for the multi-mesh shunn3 runs is not possible from FDS output, because `_mms.csv` holds mesh 1 only.
- Suggested follow-ups, not done here:
  - A VV-owned restart-equivalence case with a fixed DT.
  - Repeating obst_activation_ulmat on another build.

## Integrity

- `vv-runs/baseline/gnu_ompi_firex-36975d7/SHA256SUMS` covers every file in the baseline directory, with paths relative to it. It lists 5622 files and was written before the directory was made read-only (`chmod -R a-w`, 2026-09-25). The sha256 of SHA256SUMS itself is 892069296d03215c1212781340d5b89e48c645a138fa1bffdf14173eea3134fe. A verification run after the chmod passed.
- Verify with: `cd vv-runs/baseline/gnu_ompi_firex-36975d7 && sha256sum -c --quiet SHA256SUMS`.
- this repository and `(local FDS master checkout)` were never written. `test-plan.md`, `case_inventory.*` and `requirements.md` were not edited.

## Dropped large binaries (> 5 MB, sha256 recorded in each manifest)

420 files, 6.63 GB total not kept.

| run dir | file | MB | sha256 (prefix) |
|---|---|---|---|
| baseline/gnu_ompi_firex-36975d7/1_step_2_step_compare | 1_step_2_step_compare_1.restart | 10.7 | d3c6fd71afedc68e… |
| baseline/gnu_ompi_firex-36975d7/1_step_2_step_compare | 1_step_2_step_compare_2.restart | 10.7 | 47a214eba28559a5… |
| baseline/gnu_ompi_firex-36975d7/1_step_2_step_compare | 1_step_2_step_compare_3.restart | 10.7 | cdde695d4cecf002… |
| baseline/gnu_ompi_firex-36975d7/1_step_2_step_compare | 1_step_2_step_compare_4.restart | 10.7 | a44c213b5eda0657… |
| baseline/gnu_ompi_firex-36975d7/1_step_2_step_compare__sf17 | 1_step_2_step_compare_1.restart | 10.7 | d3c6fd71afedc68e… |
| baseline/gnu_ompi_firex-36975d7/1_step_2_step_compare__sf17 | 1_step_2_step_compare_2.restart | 10.7 | 47a214eba28559a5… |
| baseline/gnu_ompi_firex-36975d7/1_step_2_step_compare__sf17 | 1_step_2_step_compare_3.restart | 10.7 | cdde695d4cecf002… |
| baseline/gnu_ompi_firex-36975d7/1_step_2_step_compare__sf17 | 1_step_2_step_compare_4.restart | 10.7 | a44c213b5eda0657… |
| baseline/gnu_ompi_firex-36975d7/1_step_2_step_compare__sf17_np1 | 1_step_2_step_compare_1.restart | 10.7 | 6401b09cadf6e19e… |
| baseline/gnu_ompi_firex-36975d7/1_step_2_step_compare__sf17_np1 | 1_step_2_step_compare_2.restart | 10.7 | dd15e86efb8df39b… |
| baseline/gnu_ompi_firex-36975d7/1_step_2_step_compare__sf17_np1 | 1_step_2_step_compare_3.restart | 10.7 | 037bc3d8280debf3… |
| baseline/gnu_ompi_firex-36975d7/1_step_2_step_compare__sf17_np1 | 1_step_2_step_compare_4.restart | 10.7 | 72bc757a32baa7c3… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_1mesh | dancing_eddies_1mesh_1.restart | 60.4 | e632e0002906a9e4… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_1mesh | dancing_eddies_1mesh_1_1.sf | 48.1 | 46f4b63cbf81b897… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_1mesh | dancing_eddies_1mesh_1_2.sf | 48.1 | 1e7e9389da34ed72… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_1mesh | dancing_eddies_1mesh_1_3.sf | 48.1 | 49fd4fcec5eada27… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_1mesh | dancing_eddies_1mesh_1_4.sf | 48.1 | dd364b4b09556ddb… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_1mesh | dancing_eddies_1mesh_1_5.sf | 48.1 | 02051c7b77eafb60… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_1mesh__sf17 | dancing_eddies_1mesh_1.restart | 60.4 | e632e0002906a9e4… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_1mesh__sf17 | dancing_eddies_1mesh_1_1.sf | 48.1 | 46f4b63cbf81b897… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_1mesh__sf17 | dancing_eddies_1mesh_1_2.sf | 48.1 | 1e7e9389da34ed72… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_1mesh__sf17 | dancing_eddies_1mesh_1_3.sf | 48.1 | 49fd4fcec5eada27… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_1mesh__sf17 | dancing_eddies_1mesh_1_4.sf | 48.1 | dd364b4b09556ddb… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_1mesh__sf17 | dancing_eddies_1mesh_1_5.sf | 48.1 | 02051c7b77eafb60… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_1.restart | 15.4 | 50ad51b4e3d54bc4… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_1_1.sf | 12.2 | cefb3ecc637f4a34… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_1_2.sf | 12.2 | 4349ef09ef680d4d… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_1_3.sf | 12.2 | 00eee9e89fa12f14… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_1_4.sf | 12.2 | d90d74aada942e7b… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_1_5.sf | 12.2 | ae0379bda74c3cb7… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_2.restart | 15.4 | dd817f2d7c7f74ee… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_2_1.sf | 12.2 | fec0a392a08fd0cc… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_2_2.sf | 12.2 | 801a6148499bd710… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_2_3.sf | 12.2 | 0cd987d36a2b8cb5… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_2_4.sf | 12.2 | 13ad150dd7b8ef3f… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_2_5.sf | 12.2 | bc68dc315a2bf956… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_3.restart | 15.4 | 41e97a4874ce0dd8… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_3_1.sf | 12.2 | 3f420d91bf0bf778… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_3_2.sf | 12.2 | 82e21a8f6013b470… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_3_3.sf | 12.2 | edf8e8c47ea428d3… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_3_4.sf | 12.2 | 02efa676ef3d6ed0… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_3_5.sf | 12.2 | 3d3d427c504e15e9… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_4.restart | 15.4 | 2395a06438297d44… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_4_1.sf | 12.2 | f568b74c3b50b460… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_4_2.sf | 12.2 | a529cb1ad542e60b… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_4_3.sf | 12.2 | 8c9da95ffc00e6e5… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_4_4.sf | 12.2 | 7470395bca8d23ee… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default | dancing_eddies_default_4_5.sf | 12.2 | 452a0383fdf0494c… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__1mesh | dancing_eddies_default_1.restart | 60.4 | e3169a9c095b5167… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__1mesh | dancing_eddies_default_1_1.sf | 48.1 | cee6f33df463b99e… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__1mesh | dancing_eddies_default_1_2.sf | 48.1 | efcade469d3b00b8… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__1mesh | dancing_eddies_default_1_3.sf | 48.1 | 4830913df7509fab… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__1mesh | dancing_eddies_default_1_4.sf | 48.1 | d85abd8cf8c6d5c1… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__1mesh | dancing_eddies_default_1_5.sf | 48.1 | 57beb09eaa91683b… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__1mesh_sf17 | dancing_eddies_default_1.restart | 60.4 | e3169a9c095b5167… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__1mesh_sf17 | dancing_eddies_default_1_1.sf | 48.1 | cee6f33df463b99e… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__1mesh_sf17 | dancing_eddies_default_1_2.sf | 48.1 | efcade469d3b00b8… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__1mesh_sf17 | dancing_eddies_default_1_3.sf | 48.1 | 4830913df7509fab… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__1mesh_sf17 | dancing_eddies_default_1_4.sf | 48.1 | d85abd8cf8c6d5c1… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__1mesh_sf17 | dancing_eddies_default_1_5.sf | 48.1 | 57beb09eaa91683b… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_1.restart | 15.4 | 50ad51b4e3d54bc4… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_1_1.sf | 12.2 | cefb3ecc637f4a34… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_1_2.sf | 12.2 | 4349ef09ef680d4d… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_1_3.sf | 12.2 | 00eee9e89fa12f14… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_1_4.sf | 12.2 | d90d74aada942e7b… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_1_5.sf | 12.2 | ae0379bda74c3cb7… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_2.restart | 15.4 | dd817f2d7c7f74ee… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_2_1.sf | 12.2 | fec0a392a08fd0cc… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_2_2.sf | 12.2 | 801a6148499bd710… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_2_3.sf | 12.2 | 0cd987d36a2b8cb5… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_2_4.sf | 12.2 | 13ad150dd7b8ef3f… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_2_5.sf | 12.2 | bc68dc315a2bf956… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_3.restart | 15.4 | 41e97a4874ce0dd8… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_3_1.sf | 12.2 | 3f420d91bf0bf778… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_3_2.sf | 12.2 | 82e21a8f6013b470… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_3_3.sf | 12.2 | edf8e8c47ea428d3… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_3_4.sf | 12.2 | 02efa676ef3d6ed0… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_3_5.sf | 12.2 | 3d3d427c504e15e9… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_4.restart | 15.4 | 2395a06438297d44… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_4_1.sf | 12.2 | f568b74c3b50b460… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_4_2.sf | 12.2 | a529cb1ad542e60b… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_4_3.sf | 12.2 | 8c9da95ffc00e6e5… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_4_4.sf | 12.2 | 7470395bca8d23ee… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17 | dancing_eddies_default_4_5.sf | 12.2 | 452a0383fdf0494c… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_1.restart | 15.4 | a03beda6db357648… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_1_1.sf | 12.2 | c048d80c36df188d… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_1_2.sf | 12.2 | 1a05b0d9d094067b… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_1_3.sf | 12.2 | 517ba804e387cc65… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_1_4.sf | 12.2 | 691383258b1fca11… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_1_5.sf | 12.2 | c87edb8bf6f90e94… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_2.restart | 15.4 | 037660b7afddcb68… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_2_1.sf | 12.2 | a3a612f39d0c7822… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_2_2.sf | 12.2 | 055074349e7b4694… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_2_3.sf | 12.2 | 555b4eb65c6f9255… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_2_4.sf | 12.2 | d0391e46fd7aff4c… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_2_5.sf | 12.2 | e0b9e18395eb2847… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_3.restart | 15.4 | 2f35813566dda27c… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_3_1.sf | 12.2 | 67ce838084185738… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_3_2.sf | 12.2 | 6602f72c818f9c9a… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_3_3.sf | 12.2 | 43e81e9664ba894a… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_3_4.sf | 12.2 | 78738fdbfe7a9fdc… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_3_5.sf | 12.2 | c76e53a8cfb70028… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_4.restart | 15.4 | ae48e3dc74409708… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_4_1.sf | 12.2 | 0d07c88640e1f942… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_4_2.sf | 12.2 | 4d45154f61092682… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_4_3.sf | 12.2 | 4fc99a3bc449a617… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_4_4.sf | 12.2 | 3c2a169cd92e1dd2… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_glmat | dancing_eddies_default_glmat_4_5.sf | 12.2 | 51ceeeaf04943b19… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_1.restart | 15.4 | 50ad51b4e3d54bc4… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_1_1.sf | 12.2 | cefb3ecc637f4a34… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_1_2.sf | 12.2 | 4349ef09ef680d4d… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_1_3.sf | 12.2 | 00eee9e89fa12f14… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_1_4.sf | 12.2 | d90d74aada942e7b… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_1_5.sf | 12.2 | ae0379bda74c3cb7… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_2.restart | 15.4 | dd817f2d7c7f74ee… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_2_1.sf | 12.2 | fec0a392a08fd0cc… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_2_2.sf | 12.2 | 801a6148499bd710… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_2_3.sf | 12.2 | 0cd987d36a2b8cb5… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_2_4.sf | 12.2 | 13ad150dd7b8ef3f… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_2_5.sf | 12.2 | bc68dc315a2bf956… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_3.restart | 15.4 | 41e97a4874ce0dd8… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_3_1.sf | 12.2 | 3f420d91bf0bf778… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_3_2.sf | 12.2 | 82e21a8f6013b470… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_3_3.sf | 12.2 | edf8e8c47ea428d3… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_3_4.sf | 12.2 | 02efa676ef3d6ed0… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_3_5.sf | 12.2 | 3d3d427c504e15e9… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_4.restart | 15.4 | 2395a06438297d44… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_4_1.sf | 12.2 | f568b74c3b50b460… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_4_2.sf | 12.2 | a529cb1ad542e60b… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_4_3.sf | 12.2 | 8c9da95ffc00e6e5… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_4_4.sf | 12.2 | 7470395bca8d23ee… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_default__sf17_rep2 | dancing_eddies_default_4_5.sf | 12.2 | 452a0383fdf0494c… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine | dancing_eddies_uglmat_refine_1.restart | 15.6 | 0512ca6042d7bb65… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine | dancing_eddies_uglmat_refine_1_1.sf | 12.3 | aabfa5d1ccb06b49… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine | dancing_eddies_uglmat_refine_1_2.sf | 12.3 | 4b0c3aab2426f4c1… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine | dancing_eddies_uglmat_refine_1_3.sf | 12.3 | 6c571e54a1f6d8bf… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine | dancing_eddies_uglmat_refine_1_4.sf | 12.3 | 16169571f2562f0b… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine | dancing_eddies_uglmat_refine_1_5.sf | 12.3 | ee4158065b32517c… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine | dancing_eddies_uglmat_refine_2.restart | 15.6 | 4c975060d038ff9e… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine | dancing_eddies_uglmat_refine_2_1.sf | 12.3 | 185342118359a433… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine | dancing_eddies_uglmat_refine_2_2.sf | 12.3 | f5d3fcc091c039df… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine | dancing_eddies_uglmat_refine_2_3.sf | 12.3 | 4a73874190a8f4c0… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine | dancing_eddies_uglmat_refine_2_4.sf | 12.3 | f2cdb0689bdea320… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine | dancing_eddies_uglmat_refine_2_5.sf | 12.3 | 8dc1c79f9be3bb7a… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine__sf17 | dancing_eddies_uglmat_refine_1.restart | 15.6 | b94aabafa62edf68… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine__sf17 | dancing_eddies_uglmat_refine_1_1.sf | 12.3 | aabfa5d1ccb06b49… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine__sf17 | dancing_eddies_uglmat_refine_1_2.sf | 12.3 | 4b0c3aab2426f4c1… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine__sf17 | dancing_eddies_uglmat_refine_1_3.sf | 12.3 | 6c571e54a1f6d8bf… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine__sf17 | dancing_eddies_uglmat_refine_1_4.sf | 12.3 | 16169571f2562f0b… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine__sf17 | dancing_eddies_uglmat_refine_1_5.sf | 12.3 | ee4158065b32517c… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine__sf17 | dancing_eddies_uglmat_refine_2.restart | 15.6 | f5740b322d569b56… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine__sf17 | dancing_eddies_uglmat_refine_2_1.sf | 12.3 | 185342118359a433… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine__sf17 | dancing_eddies_uglmat_refine_2_2.sf | 12.3 | f5d3fcc091c039df… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine__sf17 | dancing_eddies_uglmat_refine_2_3.sf | 12.3 | 4a73874190a8f4c0… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine__sf17 | dancing_eddies_uglmat_refine_2_4.sf | 12.3 | f2cdb0689bdea320… |
| baseline/gnu_ompi_firex-36975d7/dancing_eddies_uglmat_refine__sf17 | dancing_eddies_uglmat_refine_2_5.sf | 12.3 | 8dc1c79f9be3bb7a… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3 | divergence_test_3_1.restart | 15.4 | 42bad8191f1bb10f… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3 | divergence_test_3_1_1.sf | 12.2 | 256d99eeb7b68b7e… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3 | divergence_test_3_1_2.sf | 12.2 | b3d3425d383cc1e0… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3 | divergence_test_3_1_3.sf | 12.2 | fef8844d61c9969f… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3 | divergence_test_3_1_4.sf | 12.2 | c636e7b8a1189c6d… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3 | divergence_test_3_2.restart | 15.4 | 925d2d527a70d170… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3 | divergence_test_3_2_1.sf | 12.2 | 81a4986d87e8ea49… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3 | divergence_test_3_2_2.sf | 12.2 | 1bb84028a70d4b14… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3 | divergence_test_3_2_3.sf | 12.2 | b00a33888a43ea9d… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3 | divergence_test_3_2_4.sf | 12.2 | 259affeb71ae6c31… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3 | divergence_test_3_3.restart | 15.4 | 65d32e2470d9b21a… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3 | divergence_test_3_3_1.sf | 12.2 | 013eb9fe166d0788… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3 | divergence_test_3_3_2.sf | 12.2 | ec6eb22974882a24… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3 | divergence_test_3_3_3.sf | 12.2 | 5d2973ba6f630f13… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3 | divergence_test_3_3_4.sf | 12.2 | 3a92c6d2c480d52c… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3 | divergence_test_3_4.restart | 15.4 | 94ca2fd5a87521ba… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3 | divergence_test_3_4_1.sf | 12.2 | 3f7cdc3ce5f1bad8… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3 | divergence_test_3_4_2.sf | 12.2 | 4844cb0b4070c94a… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3 | divergence_test_3_4_3.sf | 12.2 | 5ffd1c15cf1c54e2… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3 | divergence_test_3_4_4.sf | 12.2 | 32cd3fbaccd3e800… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17 | divergence_test_3_1.restart | 15.4 | 42bad8191f1bb10f… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17 | divergence_test_3_1_1.sf | 12.2 | 256d99eeb7b68b7e… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17 | divergence_test_3_1_2.sf | 12.2 | b3d3425d383cc1e0… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17 | divergence_test_3_1_3.sf | 12.2 | fef8844d61c9969f… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17 | divergence_test_3_1_4.sf | 12.2 | c636e7b8a1189c6d… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17 | divergence_test_3_2.restart | 15.4 | 925d2d527a70d170… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17 | divergence_test_3_2_1.sf | 12.2 | 81a4986d87e8ea49… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17 | divergence_test_3_2_2.sf | 12.2 | 1bb84028a70d4b14… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17 | divergence_test_3_2_3.sf | 12.2 | b00a33888a43ea9d… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17 | divergence_test_3_2_4.sf | 12.2 | 259affeb71ae6c31… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17 | divergence_test_3_3.restart | 15.4 | 65d32e2470d9b21a… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17 | divergence_test_3_3_1.sf | 12.2 | 013eb9fe166d0788… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17 | divergence_test_3_3_2.sf | 12.2 | ec6eb22974882a24… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17 | divergence_test_3_3_3.sf | 12.2 | 5d2973ba6f630f13… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17 | divergence_test_3_3_4.sf | 12.2 | 3a92c6d2c480d52c… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17 | divergence_test_3_4.restart | 15.4 | 94ca2fd5a87521ba… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17 | divergence_test_3_4_1.sf | 12.2 | 3f7cdc3ce5f1bad8… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17 | divergence_test_3_4_2.sf | 12.2 | 4844cb0b4070c94a… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17 | divergence_test_3_4_3.sf | 12.2 | 5ffd1c15cf1c54e2… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17 | divergence_test_3_4_4.sf | 12.2 | 32cd3fbaccd3e800… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17_glmat | divergence_test_3_glmat_1.restart | 15.4 | 2c53c23c3e79725b… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17_glmat | divergence_test_3_glmat_1_1.sf | 12.2 | 1b569630218b0d66… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17_glmat | divergence_test_3_glmat_1_2.sf | 12.2 | d2a9d421f681c7ee… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17_glmat | divergence_test_3_glmat_1_3.sf | 12.2 | 4363d44b6b419cac… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17_glmat | divergence_test_3_glmat_1_4.sf | 12.2 | c7050a7d30701c49… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17_glmat | divergence_test_3_glmat_2.restart | 15.4 | 91d0caab3001e1d5… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17_glmat | divergence_test_3_glmat_2_1.sf | 12.2 | 3694291200c8377b… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17_glmat | divergence_test_3_glmat_2_2.sf | 12.2 | c046165f4c67cd6e… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17_glmat | divergence_test_3_glmat_2_3.sf | 12.2 | 48ab54a15e50e07a… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17_glmat | divergence_test_3_glmat_2_4.sf | 12.2 | f36e1efbbe48e18e… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17_glmat | divergence_test_3_glmat_3.restart | 15.4 | 7edd47737f102fcc… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17_glmat | divergence_test_3_glmat_3_1.sf | 12.2 | 0d869851ce0f62f6… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17_glmat | divergence_test_3_glmat_3_2.sf | 12.2 | c0b0d79f9dfea6ea… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17_glmat | divergence_test_3_glmat_3_3.sf | 12.2 | 890a78450d0d0d19… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17_glmat | divergence_test_3_glmat_3_4.sf | 12.2 | e2d9d44077010008… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17_glmat | divergence_test_3_glmat_4.restart | 15.4 | 15fef9d164ffc7ab… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17_glmat | divergence_test_3_glmat_4_1.sf | 12.2 | 3234d9651fdb614f… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17_glmat | divergence_test_3_glmat_4_2.sf | 12.2 | 9257e565a12c6583… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17_glmat | divergence_test_3_glmat_4_3.sf | 12.2 | 98063ec01c6af648… |
| baseline/gnu_ompi_firex-36975d7/divergence_test_3__sf17_glmat | divergence_test_3_glmat_4_4.sf | 12.2 | 80755f6d165051c1… |
| baseline/gnu_ompi_firex-36975d7/layer_1mesh | layer_1mesh_1.restart | 9.5 | fd93b770ede86b75… |
| baseline/gnu_ompi_firex-36975d7/layer_1mesh | layer_1mesh_1_1.s3d | 8.4 | 0801fb04ca5bc798… |
| baseline/gnu_ompi_firex-36975d7/layer_1mesh | layer_1mesh_1_3.s3d | 5.6 | 6e632e0e13b2938d… |
| baseline/gnu_ompi_firex-36975d7/layer_1mesh__sf17 | layer_1mesh_1.restart | 9.5 | fd93b770ede86b75… |
| baseline/gnu_ompi_firex-36975d7/layer_1mesh__sf17 | layer_1mesh_1_1.s3d | 8.4 | 0801fb04ca5bc798… |
| baseline/gnu_ompi_firex-36975d7/layer_1mesh__sf17 | layer_1mesh_1_3.s3d | 5.6 | 6e632e0e13b2938d… |
| baseline/gnu_ompi_firex-36975d7/ns2d_32 | ns2d_32_1.restart | 5.7 | f1d87f6ff49e6577… |
| baseline/gnu_ompi_firex-36975d7/ns2d_32__sf17 | ns2d_32_1.restart | 5.7 | f1d87f6ff49e6577… |
| baseline/gnu_ompi_firex-36975d7/ns2d_32_nupt1 | ns2d_32_nupt1_1.restart | 5.7 | ca5b8b66fb049ab2… |
| baseline/gnu_ompi_firex-36975d7/ns2d_32_nupt1__sf17 | ns2d_32_nupt1_1.restart | 5.7 | ca5b8b66fb049ab2… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64 | ns2d_64_1.restart | 21.5 | cee0899ff8afbbda… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64 | ns2d_64_1_1.sf | 16.9 | 2932f23ede94bd8d… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64 | ns2d_64_1_2.sf | 16.9 | 52f801eed74141c8… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64 | ns2d_64_1_3.sf | 16.9 | 9d365ac7c00a6927… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64 | ns2d_64_1_4.sf | 16.9 | ee2f30b5ec4c546a… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64 | ns2d_64_1_5.sf | 16.9 | af45393f03a92a4e… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64 | ns2d_64_1_6.sf | 16.9 | ae48c634b25ff592… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64 | ns2d_64_1_7.sf | 16.9 | 4b9da157c29718f6… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64 | ns2d_64_1_8.sf | 16.9 | 975b62ada58d8cb9… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64__sf17 | ns2d_64_1.restart | 21.5 | cee0899ff8afbbda… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64__sf17 | ns2d_64_1_1.sf | 16.9 | 2932f23ede94bd8d… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64__sf17 | ns2d_64_1_2.sf | 16.9 | 52f801eed74141c8… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64__sf17 | ns2d_64_1_3.sf | 16.9 | 9d365ac7c00a6927… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64__sf17 | ns2d_64_1_4.sf | 16.9 | ee2f30b5ec4c546a… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64__sf17 | ns2d_64_1_5.sf | 16.9 | af45393f03a92a4e… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64__sf17 | ns2d_64_1_6.sf | 16.9 | ae48c634b25ff592… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64__sf17 | ns2d_64_1_7.sf | 16.9 | 4b9da157c29718f6… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64__sf17 | ns2d_64_1_8.sf | 16.9 | 975b62ada58d8cb9… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64_nupt1 | ns2d_64_nupt1_1.restart | 21.5 | 77f81cf8346573ad… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64_nupt1 | ns2d_64_nupt1_1_1.sf | 15.7 | 994d2b51a98fa5b2… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64_nupt1 | ns2d_64_nupt1_1_2.sf | 15.7 | bb0a8f7b6c072dbe… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64_nupt1 | ns2d_64_nupt1_1_3.sf | 15.7 | ade0cfa1d8269c7b… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64_nupt1 | ns2d_64_nupt1_1_4.sf | 15.7 | 5abf8b295fd22397… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64_nupt1 | ns2d_64_nupt1_1_5.sf | 15.7 | 5e4c66820ab8e7b7… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64_nupt1 | ns2d_64_nupt1_1_6.sf | 15.7 | 64e52b86b3e19fd1… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64_nupt1 | ns2d_64_nupt1_1_7.sf | 15.7 | 9bd58a6adaf9706d… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64_nupt1 | ns2d_64_nupt1_1_8.sf | 15.7 | 673f50a305dd3d84… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64_nupt1__sf17 | ns2d_64_nupt1_1.restart | 21.5 | 77f81cf8346573ad… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64_nupt1__sf17 | ns2d_64_nupt1_1_1.sf | 15.7 | 994d2b51a98fa5b2… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64_nupt1__sf17 | ns2d_64_nupt1_1_2.sf | 15.7 | bb0a8f7b6c072dbe… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64_nupt1__sf17 | ns2d_64_nupt1_1_3.sf | 15.7 | ade0cfa1d8269c7b… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64_nupt1__sf17 | ns2d_64_nupt1_1_4.sf | 15.7 | 5abf8b295fd22397… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64_nupt1__sf17 | ns2d_64_nupt1_1_5.sf | 15.7 | 5e4c66820ab8e7b7… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64_nupt1__sf17 | ns2d_64_nupt1_1_6.sf | 15.7 | 64e52b86b3e19fd1… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64_nupt1__sf17 | ns2d_64_nupt1_1_7.sf | 15.7 | 9bd58a6adaf9706d… |
| baseline/gnu_ompi_firex-36975d7/ns2d_64_nupt1__sf17 | ns2d_64_nupt1_1_8.sf | 15.7 | 673f50a305dd3d84… |
| baseline/gnu_ompi_firex-36975d7/restart_test1_continuous | restart_test1_continuous_1.prt5 | 6.1 | 8ab220d995ec9e4b… |
| baseline/gnu_ompi_firex-36975d7/restart_test1_continuous | restart_test1_continuous_1.restart | 18.4 | 776eb88638c41f5c… |
| baseline/gnu_ompi_firex-36975d7/restart_test1_continuous | restart_test1_continuous_1_1.iso | 24.4 | 7e5fb0bc220445b4… |
| baseline/gnu_ompi_firex-36975d7/restart_test1_continuous | restart_test1_continuous_1_2.iso | 8.0 | 882237cf5cb9b4b9… |
| baseline/gnu_ompi_firex-36975d7/restart_test1_continuous | restart_test1_continuous_1_4.iso | 12.2 | fe87e18cc53e6e4d… |
| baseline/gnu_ompi_firex-36975d7/restart_test1_continuous__sf17 | restart_test1_continuous_1.prt5 | 6.1 | 8ab220d995ec9e4b… |
| baseline/gnu_ompi_firex-36975d7/restart_test1_continuous__sf17 | restart_test1_continuous_1.restart | 18.4 | 13c99a96e9295d44… |
| baseline/gnu_ompi_firex-36975d7/restart_test1_continuous__sf17 | restart_test1_continuous_1_1.iso | 24.4 | 7e5fb0bc220445b4… |
| baseline/gnu_ompi_firex-36975d7/restart_test1_continuous__sf17 | restart_test1_continuous_1_2.iso | 8.0 | 882237cf5cb9b4b9… |
| baseline/gnu_ompi_firex-36975d7/restart_test1_continuous__sf17 | restart_test1_continuous_1_4.iso | 12.2 | fe87e18cc53e6e4d… |
| baseline/gnu_ompi_firex-36975d7/restart_test1a | restart_test1a_1.prt5 | 5.3 | 00295bbafab89881… |
| baseline/gnu_ompi_firex-36975d7/restart_test1a | restart_test1a_1.restart | 18.9 | 62ce30105e00e00a… |
| baseline/gnu_ompi_firex-36975d7/restart_test1a | restart_test1a_1_1.iso | 18.9 | 15616269fd0e9a54… |
| baseline/gnu_ompi_firex-36975d7/restart_test1a | restart_test1a_1_2.iso | 6.7 | b81affeb73249854… |
| baseline/gnu_ompi_firex-36975d7/restart_test1a | restart_test1a_1_4.iso | 10.0 | 8c8dd59454badfad… |
| baseline/gnu_ompi_firex-36975d7/restart_test1a__sf17 | restart_test1a_1.prt5 | 5.3 | 00295bbafab89881… |
| baseline/gnu_ompi_firex-36975d7/restart_test1a__sf17 | restart_test1a_1.restart | 18.9 | 5327500b56de31af… |
| baseline/gnu_ompi_firex-36975d7/restart_test1a__sf17 | restart_test1a_1_1.iso | 18.9 | 15616269fd0e9a54… |
| baseline/gnu_ompi_firex-36975d7/restart_test1a__sf17 | restart_test1a_1_2.iso | 6.7 | b81affeb73249854… |
| baseline/gnu_ompi_firex-36975d7/restart_test1a__sf17 | restart_test1a_1_4.iso | 10.0 | 8c8dd59454badfad… |
| baseline/gnu_ompi_firex-36975d7/restart_test1b | restart_test1a_1.prt5 | 11.2 | 309260df5854d0b4… |
| baseline/gnu_ompi_firex-36975d7/restart_test1b | restart_test1a_1.restart | 18.4 | 90896dd219d685de… |
| baseline/gnu_ompi_firex-36975d7/restart_test1b | restart_test1a_1_1.iso | 44.2 | de8f2468e54650c1… |
| baseline/gnu_ompi_firex-36975d7/restart_test1b | restart_test1a_1_1.s3d | 5.3 | 33b25019050f9b19… |
| baseline/gnu_ompi_firex-36975d7/restart_test1b | restart_test1a_1_2.iso | 13.9 | 1cbebb2be1a9c8da… |
| baseline/gnu_ompi_firex-36975d7/restart_test1b | restart_test1a_1_3.s3d | 5.6 | a99801bc5c5ac4bc… |
| baseline/gnu_ompi_firex-36975d7/restart_test1b | restart_test1a_1_4.iso | 21.1 | fa4bce1c8a5e9678… |
| baseline/gnu_ompi_firex-36975d7/restart_test1b__sf17 | restart_test1a_1.prt5 | 11.2 | 309260df5854d0b4… |
| baseline/gnu_ompi_firex-36975d7/restart_test1b__sf17 | restart_test1a_1.restart | 18.4 | 5879df7fb7e9efa0… |
| baseline/gnu_ompi_firex-36975d7/restart_test1b__sf17 | restart_test1a_1_1.iso | 44.2 | de8f2468e54650c1… |
| baseline/gnu_ompi_firex-36975d7/restart_test1b__sf17 | restart_test1a_1_1.s3d | 5.3 | 33b25019050f9b19… |
| baseline/gnu_ompi_firex-36975d7/restart_test1b__sf17 | restart_test1a_1_2.iso | 13.9 | 1cbebb2be1a9c8da… |
| baseline/gnu_ompi_firex-36975d7/restart_test1b__sf17 | restart_test1a_1_3.s3d | 5.6 | a99801bc5c5ac4bc… |
| baseline/gnu_ompi_firex-36975d7/restart_test1b__sf17 | restart_test1a_1_4.iso | 21.1 | fa4bce1c8a5e9678… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_1 | saad_512_cfl_1_1.restart | 12.6 | 5a15f46da3ea7a34… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_1__sf17 | saad_512_cfl_1_1.restart | 12.6 | 5a15f46da3ea7a34… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p0625 | saad_512_cfl_p0625_1.restart | 12.6 | 03be9f05abaf0562… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p0625 | saad_512_cfl_p0625_1_1.sf | 10.3 | a1c0103a5bb1b030… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p0625 | saad_512_cfl_p0625_1_10.sf | 10.3 | e6177cd0c6acab95… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p0625 | saad_512_cfl_p0625_1_2.sf | 10.3 | 1154ec431eac0891… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p0625 | saad_512_cfl_p0625_1_3.sf | 10.3 | b8fde9c48e23e5b2… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p0625 | saad_512_cfl_p0625_1_4.sf | 10.3 | 1ee4901036ebc940… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p0625 | saad_512_cfl_p0625_1_5.sf | 10.3 | ffe629471087f084… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p0625 | saad_512_cfl_p0625_1_6.sf | 10.3 | dc79029128e2ee12… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p0625 | saad_512_cfl_p0625_1_7.sf | 10.3 | 056f63941f2f0c87… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p0625 | saad_512_cfl_p0625_1_8.sf | 10.3 | 126e0d3e3e1e3df4… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p0625 | saad_512_cfl_p0625_1_9.sf | 10.3 | c4ee5a7cabee1c82… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p0625__sf17 | saad_512_cfl_p0625_1.restart | 12.6 | 03be9f05abaf0562… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p0625__sf17 | saad_512_cfl_p0625_1_1.sf | 10.3 | a1c0103a5bb1b030… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p0625__sf17 | saad_512_cfl_p0625_1_10.sf | 10.3 | e6177cd0c6acab95… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p0625__sf17 | saad_512_cfl_p0625_1_2.sf | 10.3 | 1154ec431eac0891… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p0625__sf17 | saad_512_cfl_p0625_1_3.sf | 10.3 | b8fde9c48e23e5b2… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p0625__sf17 | saad_512_cfl_p0625_1_4.sf | 10.3 | 1ee4901036ebc940… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p0625__sf17 | saad_512_cfl_p0625_1_5.sf | 10.3 | ffe629471087f084… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p0625__sf17 | saad_512_cfl_p0625_1_6.sf | 10.3 | dc79029128e2ee12… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p0625__sf17 | saad_512_cfl_p0625_1_7.sf | 10.3 | 056f63941f2f0c87… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p0625__sf17 | saad_512_cfl_p0625_1_8.sf | 10.3 | 126e0d3e3e1e3df4… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p0625__sf17 | saad_512_cfl_p0625_1_9.sf | 10.3 | c4ee5a7cabee1c82… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p125 | saad_512_cfl_p125_1.restart | 12.6 | 1a038412cf549336… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p125 | saad_512_cfl_p125_1_1.sf | 8.2 | d05df2fc791f3270… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p125 | saad_512_cfl_p125_1_10.sf | 8.2 | dc0f4714266a4d9a… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p125 | saad_512_cfl_p125_1_2.sf | 8.2 | 9ea0c6ee8e9f522c… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p125 | saad_512_cfl_p125_1_3.sf | 8.2 | 6bee4e2644cd8b17… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p125 | saad_512_cfl_p125_1_4.sf | 8.2 | 859adced21531940… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p125 | saad_512_cfl_p125_1_5.sf | 8.2 | 87b25b492c94edb6… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p125 | saad_512_cfl_p125_1_6.sf | 8.2 | 01bb1ae9e5ce0e58… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p125 | saad_512_cfl_p125_1_7.sf | 8.2 | 25be9996271762ea… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p125 | saad_512_cfl_p125_1_8.sf | 8.2 | ea601ff536dc20dc… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p125 | saad_512_cfl_p125_1_9.sf | 8.2 | 053de75ea2d7de98… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p125__sf17 | saad_512_cfl_p125_1.restart | 12.6 | 1a038412cf549336… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p125__sf17 | saad_512_cfl_p125_1_1.sf | 8.2 | d05df2fc791f3270… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p125__sf17 | saad_512_cfl_p125_1_10.sf | 8.2 | dc0f4714266a4d9a… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p125__sf17 | saad_512_cfl_p125_1_2.sf | 8.2 | 9ea0c6ee8e9f522c… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p125__sf17 | saad_512_cfl_p125_1_3.sf | 8.2 | 6bee4e2644cd8b17… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p125__sf17 | saad_512_cfl_p125_1_4.sf | 8.2 | 859adced21531940… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p125__sf17 | saad_512_cfl_p125_1_5.sf | 8.2 | 87b25b492c94edb6… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p125__sf17 | saad_512_cfl_p125_1_6.sf | 8.2 | 01bb1ae9e5ce0e58… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p125__sf17 | saad_512_cfl_p125_1_7.sf | 8.2 | 25be9996271762ea… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p125__sf17 | saad_512_cfl_p125_1_8.sf | 8.2 | ea601ff536dc20dc… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p125__sf17 | saad_512_cfl_p125_1_9.sf | 8.2 | 053de75ea2d7de98… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p25 | saad_512_cfl_p25_1.restart | 12.6 | 97cebe10966ec572… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p25__sf17 | saad_512_cfl_p25_1.restart | 12.6 | 97cebe10966ec572… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p5 | saad_512_cfl_p5_1.restart | 12.6 | 3e0607f03d3b49dc… |
| baseline/gnu_ompi_firex-36975d7/saad_512_cfl_p5__sf17 | saad_512_cfl_p5_1.restart | 12.6 | 3e0607f03d3b49dc… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128 | shunn3_128_1.restart | 70.1 | bba082fca317d73e… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128 | shunn3_128_1_1.sf | 12.5 | 072c7e035c093444… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128 | shunn3_128_1_10.sf | 12.5 | 46e597284fd2951d… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128 | shunn3_128_1_11.sf | 12.5 | b579e61aeacd9fba… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128 | shunn3_128_1_12.sf | 12.5 | ae857a40b03d18e9… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128 | shunn3_128_1_2.sf | 12.5 | dc42d2f476669619… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128 | shunn3_128_1_3.sf | 12.5 | 4b196340d8156f36… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128 | shunn3_128_1_4.sf | 12.5 | 770f17c99ed9a037… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128 | shunn3_128_1_5.sf | 12.5 | c161ff807cb94dff… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128 | shunn3_128_1_6.sf | 12.5 | 2fa4b1db855605f8… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128 | shunn3_128_1_7.sf | 12.5 | b899cf10f6664b5f… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128 | shunn3_128_1_8.sf | 12.5 | 055938515455de1f… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128 | shunn3_128_1_9.sf | 12.5 | a01facf230452e73… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128__sf17 | shunn3_128_1.restart | 70.1 | bba082fca317d73e… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128__sf17 | shunn3_128_1_1.sf | 12.5 | 072c7e035c093444… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128__sf17 | shunn3_128_1_10.sf | 12.5 | 46e597284fd2951d… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128__sf17 | shunn3_128_1_11.sf | 12.5 | b579e61aeacd9fba… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128__sf17 | shunn3_128_1_12.sf | 12.5 | ae857a40b03d18e9… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128__sf17 | shunn3_128_1_2.sf | 12.5 | dc42d2f476669619… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128__sf17 | shunn3_128_1_3.sf | 12.5 | 4b196340d8156f36… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128__sf17 | shunn3_128_1_4.sf | 12.5 | 770f17c99ed9a037… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128__sf17 | shunn3_128_1_5.sf | 12.5 | c161ff807cb94dff… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128__sf17 | shunn3_128_1_6.sf | 12.5 | 2fa4b1db855605f8… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128__sf17 | shunn3_128_1_7.sf | 12.5 | b899cf10f6664b5f… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128__sf17 | shunn3_128_1_8.sf | 12.5 | 055938515455de1f… |
| baseline/gnu_ompi_firex-36975d7/shunn3_128__sf17 | shunn3_128_1_9.sf | 12.5 | a01facf230452e73… |
| baseline/gnu_ompi_firex-36975d7/shunn3_4mesh_128 | shunn3_4mesh_128_1.restart | 19.1 | 7e340550bcc63da7… |
| baseline/gnu_ompi_firex-36975d7/shunn3_4mesh_128 | shunn3_4mesh_128_2.restart | 19.1 | 759b0a16d6fc58e9… |
| baseline/gnu_ompi_firex-36975d7/shunn3_4mesh_128 | shunn3_4mesh_128_3.restart | 19.1 | 5ceeb8b221845360… |
| baseline/gnu_ompi_firex-36975d7/shunn3_4mesh_128 | shunn3_4mesh_128_4.restart | 19.1 | b8143a6693eba1ac… |
| baseline/gnu_ompi_firex-36975d7/shunn3_4mesh_128__sf17 | shunn3_4mesh_128_1.restart | 19.1 | 7e340550bcc63da7… |
| baseline/gnu_ompi_firex-36975d7/shunn3_4mesh_128__sf17 | shunn3_4mesh_128_2.restart | 19.1 | 759b0a16d6fc58e9… |
| baseline/gnu_ompi_firex-36975d7/shunn3_4mesh_128__sf17 | shunn3_4mesh_128_3.restart | 19.1 | 5ceeb8b221845360… |
| baseline/gnu_ompi_firex-36975d7/shunn3_4mesh_128__sf17 | shunn3_4mesh_128_4.restart | 19.1 | b8143a6693eba1ac… |
| baseline/gnu_ompi_firex-36975d7/shunn3_4mesh_128__sf17_glmat | shunn3_4mesh_128_glmat_1.restart | 19.1 | 49d9188e3172cc0d… |
| baseline/gnu_ompi_firex-36975d7/shunn3_4mesh_128__sf17_glmat | shunn3_4mesh_128_glmat_2.restart | 19.1 | 27c7f82bfc973018… |
| baseline/gnu_ompi_firex-36975d7/shunn3_4mesh_128__sf17_glmat | shunn3_4mesh_128_glmat_3.restart | 19.1 | 6aa402272be9e7a5… |
| baseline/gnu_ompi_firex-36975d7/shunn3_4mesh_128__sf17_glmat | shunn3_4mesh_128_glmat_4.restart | 19.1 | 920e6f08de20d6db… |
| baseline/gnu_ompi_firex-36975d7/shunn3_64 | shunn3_64_1.restart | 18.0 | f100d442b93ce0aa… |
| baseline/gnu_ompi_firex-36975d7/shunn3_64__sf17 | shunn3_64_1.restart | 18.0 | f100d442b93ce0aa… |
| baseline/gnu_ompi_firex-36975d7/soborot_superbee_square_wave_128 | soborot_superbee_square_wave_128_1.restart | 17.1 | bb7f22a4c595bbe2… |
| baseline/gnu_ompi_firex-36975d7/soborot_superbee_square_wave_128 | soborot_superbee_square_wave_128_2.restart | 17.1 | bc490f278ed2ce66… |
| baseline/gnu_ompi_firex-36975d7/soborot_superbee_square_wave_128 | soborot_superbee_square_wave_128_3.restart | 17.1 | f00058ad2ffef71c… |
| baseline/gnu_ompi_firex-36975d7/soborot_superbee_square_wave_128 | soborot_superbee_square_wave_128_4.restart | 17.1 | 0d99ea27fbe45839… |
| baseline/gnu_ompi_firex-36975d7/soborot_superbee_square_wave_128_1mesh | soborot_superbee_square_wave_128_1mesh_1.restart | 66.0 | 4786607bce14f158… |
| baseline/gnu_ompi_firex-36975d7/soborot_superbee_square_wave_128_1mesh__sf17 | soborot_superbee_square_wave_128_1mesh_1.restart | 66.0 | 4786607bce14f158… |
| baseline/gnu_ompi_firex-36975d7/soborot_superbee_square_wave_128__sf17 | soborot_superbee_square_wave_128_1.restart | 17.1 | ad8434116e001d74… |
| baseline/gnu_ompi_firex-36975d7/soborot_superbee_square_wave_128__sf17 | soborot_superbee_square_wave_128_2.restart | 17.1 | f34c07502409a85a… |
| baseline/gnu_ompi_firex-36975d7/soborot_superbee_square_wave_128__sf17 | soborot_superbee_square_wave_128_3.restart | 17.1 | a711214c8aaf1583… |
| baseline/gnu_ompi_firex-36975d7/soborot_superbee_square_wave_128__sf17 | soborot_superbee_square_wave_128_4.restart | 17.1 | 3b667ec240ed6e23… |
| baseline/gnu_ompi_firex-36975d7/soborot_superbee_square_wave_128__sf17_glmat | soborot_superbee_square_wave_128_glmat_1.restart | 17.1 | 4f1dc2cc28ff05f6… |
| baseline/gnu_ompi_firex-36975d7/soborot_superbee_square_wave_128__sf17_glmat | soborot_superbee_square_wave_128_glmat_2.restart | 17.1 | 2f8ca5323cb48ddf… |
| baseline/gnu_ompi_firex-36975d7/soborot_superbee_square_wave_128__sf17_glmat | soborot_superbee_square_wave_128_glmat_3.restart | 17.1 | 22b92aaeb039fde0… |
| baseline/gnu_ompi_firex-36975d7/soborot_superbee_square_wave_128__sf17_glmat | soborot_superbee_square_wave_128_glmat_4.restart | 17.1 | d5c4d9a843f5e0f3… |
| baseline/gnu_ompi_firex-36975d7/soborot_superbee_square_wave_128__sf17_np1 | soborot_superbee_square_wave_128_1.restart | 17.1 | 6468b09eca4f51e0… |
| baseline/gnu_ompi_firex-36975d7/soborot_superbee_square_wave_128__sf17_np1 | soborot_superbee_square_wave_128_2.restart | 17.1 | 5de811842848274b… |
| baseline/gnu_ompi_firex-36975d7/soborot_superbee_square_wave_128__sf17_np1 | soborot_superbee_square_wave_128_3.restart | 17.1 | 6e09409f6c498125… |
| baseline/gnu_ompi_firex-36975d7/soborot_superbee_square_wave_128__sf17_np1 | soborot_superbee_square_wave_128_4.restart | 17.1 | e0c992619c3905e7… |
| baseline/gnu_ompi_firex-36975d7/species_conservation_1 | species_conservation_1_1.restart | 8.8 | 36a09d844003bc1a… |
| baseline/gnu_ompi_firex-36975d7/species_conservation_1__sf17 | species_conservation_1_1.restart | 8.8 | 36a09d844003bc1a… |
| baseline/gnu_ompi_firex-36975d7/symmetry_test_mpi__1mesh | symmetry_test_mpi_1.restart | 6.7 | bb66379e5b7eb419… |
| baseline/gnu_ompi_firex-36975d7/symmetry_test_mpi__1mesh_sf17 | symmetry_test_mpi_1.restart | 6.7 | 36b8ff6bfb4e744f… |
| calib/gnu_ompi_firex-36975d7/restart_test1a__sf17_rep2 | restart_test1a_1.prt5 | 5.3 | 00295bbafab89881… |
| calib/gnu_ompi_firex-36975d7/restart_test1a__sf17_rep2 | restart_test1a_1.restart | 18.9 | 63c7243c1b310227… |
| calib/gnu_ompi_firex-36975d7/restart_test1a__sf17_rep2 | restart_test1a_1_1.iso | 18.9 | 15616269fd0e9a54… |
| calib/gnu_ompi_firex-36975d7/restart_test1a__sf17_rep2 | restart_test1a_1_2.iso | 6.7 | b81affeb73249854… |
| calib/gnu_ompi_firex-36975d7/restart_test1a__sf17_rep2 | restart_test1a_1_4.iso | 10.0 | 8c8dd59454badfad… |
| calib/gnu_ompi_firex-36975d7/restart_test1b__sf17_rep2 | restart_test1a_1.prt5 | 11.2 | 309260df5854d0b4… |
| calib/gnu_ompi_firex-36975d7/restart_test1b__sf17_rep2 | restart_test1a_1.restart | 18.4 | 1a53d7fc10f2fa23… |
| calib/gnu_ompi_firex-36975d7/restart_test1b__sf17_rep2 | restart_test1a_1_1.iso | 44.2 | de8f2468e54650c1… |
| calib/gnu_ompi_firex-36975d7/restart_test1b__sf17_rep2 | restart_test1a_1_1.s3d | 5.3 | 33b25019050f9b19… |
| calib/gnu_ompi_firex-36975d7/restart_test1b__sf17_rep2 | restart_test1a_1_2.iso | 13.9 | 1cbebb2be1a9c8da… |
| calib/gnu_ompi_firex-36975d7/restart_test1b__sf17_rep2 | restart_test1a_1_3.s3d | 5.6 | a99801bc5c5ac4bc… |
| calib/gnu_ompi_firex-36975d7/restart_test1b__sf17_rep2 | restart_test1a_1_4.iso | 21.1 | fa4bce1c8a5e9678… |
