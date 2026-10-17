# 15 — Effect of a looser predictor pressure tolerance on the solution

Status: **MEASURED**, owner-provided NVIDIA test machine (CPU run, 1 rank, 1 OpenMP thread pinned: `OMP_PROC_BIND=close`, `OMP_PLACES=cores`; GCC 13, Open MPI 5, HYPRE 2.32, AMReX 26.09). Source: tip `d09c02b93c` plus a **private scratch patch, not committed and not in any product file** (two environment switches, below). No default was changed; the product default stays `tol_rel` 1e-12 for both solves.

## 1. Question
Section 8.4 of `14-composite-solve-profile.md` measured a predictor-only tolerance (corrector unchanged) as cheaper, with the divergence gate unchanged because the corrector decides the final divergence. The objection to using it: the predicted velocity advances density, temperature and species into the corrector, so a looser predictor solve changes the predicted state, not only the pressure. This document measures how much the **state** changes.

## 2. Method
- Cases (all 3 levels except CSIRO, 2 levels): **B** near-burner 3-level (`tests/cases/burner_obst_3level.fds`, 49 steps, no reaction: HRR is identically 0, so the HRR rows are not exercised there); **S** `simple_test` 3-level (`--tag-vort auto`, 49 steps, burning propane); **C** CSIRO `Case_C064_LS` (24x24x8 level 0, `--levels 2 --init-patch --patch 10 19 0 23 --blocking 2 --maxsize 8 --regrid-every 2 --tag-phi --nofine 1`, 60 steps, about 15% of the 120 s case, lean).
- Predictor tolerance 1e-12 (baseline), 1e-10, 1e-8; the corrector solve is 1e-12 in every run. Baseline run **twice** (a, b) for the repeat-run noise floor. Same binary, 1 rank, 1 thread, same input, same command in all runs; runs sequential.
- Override mechanism (scratch only): environment variable `FDSTL_PRED_TOL` sets `tol_rel` (and `residual_tol`, so the MLMG polish pass does not tighten it back) for the predictor solve only; the existing `--pressure-tol` flag is not usable for this because it applies to both solves and turns on the study path. The earlier throwaway switch of section 7.7 did not change the cycle count (void); this one does (the solve time falls with the tolerance, table 3).
- Diagnostic (scratch only): after every step the driver wrote every **leaf gas cell** (uncovered by a finer level, not a fine OBST cell) with TMP, RHO, Q (heat release per volume) and ZZ(1..ns); a script compared run vs baseline cell by cell, same cell index, every step. Layouts (box sets) were identical in all runs at every step. Quantities per step: max abs and max relative difference of TMP, RHO, each species mass fraction Z, Q; total HRR (sum Q dV); max-cell Q; total mass (sum RHO dV). Pointwise relative = |a-b| / max(|b|, 1e-6 max|b| of that field), the floor keeps cells where a species or Q is zero from dividing by about 0. Max is over cells and over all steps.

## 3. Results

Repeat-run noise floor (baseline a vs baseline b): **exactly 0** in every quantity, every step, all three cases (bitwise identical; divergence gate values identical to all digits). So every difference below is caused by the predictor tolerance alone, not by run-to-run noise.

### 3.1 Maximum differences against the 1e-12 baseline (over leaf cells and all steps)
| Case | pred tol | TMP abs (K) | TMP rel | RHO abs (kg/m3) | RHO rel | Z abs (worst of species) | Z rel (worst) | HRR total abs (W) / rel | max-cell Q abs (W/m3) / rel | total mass abs (kg) / rel |
|---|---|---|---|---|---|---|---|---|---|---|
| B | 1e-10 | 9.5e-8 | 3.2e-10 | 3.8e-10 | 3.1e-10 | 1.6e-10 | 7.1e-9 | 0 (no reaction) | 0 | 2.1e-14 / 8.6e-16 |
| B | 1e-8 | 8.5e-6 | 2.9e-8 | 3.4e-8 | 2.8e-8 | 1.4e-8 | 6.7e-7 | 0 | 0 | 1.9e-12 / 7.5e-14 |
| S | 1e-10 | 3.9e-6 | 1.4e-9 | 2.0e-10 | 1.4e-9 | 6.1e-10 | 2.6e-8 | 2.6e-6 / 1.7e-11 | 8.0e-4 / 5.7e-11 | 1.3e-12 / 5.1e-14 |
| S | 1e-8 | 2.1e-4 | 7.3e-8 | 1.8e-8 | 7.4e-8 | 4.8e-8 | 3.1e-6 | 6.8e-5 / 4.5e-10 | 5.1e-2 / 3.2e-9 | 1.8e-11 / 7.2e-13 |
| C | 1e-10 | 7.3e-6 | 2.3e-8 | 2.5e-8 | 2.3e-8 | 4.3e-10 | 2.0e-6 | 0.21 / 1.4e-9 | 1.2e-4 / 1.2e-9 | 3.2e-6 / 4.9e-12 |
| C | 1e-8 | 1.2e-3 | 3.9e-6 | 4.3e-6 | 3.9e-6 | 8.3e-8 | 5.8e-4 | 36.5 / 2.4e-7 | 2.1e-2 / 2.2e-7 | 6.6e-4 / 1.0e-9 |

Pointwise-relative Z is dominated by cells where the species is nearly absent (floor 1e-6 of the field maximum); normalized by the field maximum the species differences are about as large as the RHO relative differences (column "Z abs" over Z of order 1). Per-species values (B: Z1, Z2; S and C: Z1, Z2, Z3 in FDS ZZ index order) are in the scratch series files; the worst species is shown. HRR "total" and "max cell" for S and C are the only heat-release rows with data.

### 3.2 Time series of the total leaf mass difference (run minus baseline, kg; relative in parentheses)
| Case | pred tol | step 1 | step 5 | step 10 | step 20 | step 30 | last step |
|---|---|---|---|---|---|---|---|
| B | 1e-10 | 0 | 0 | 0 | 0 | 0 | -2.1e-14 (-8.6e-16) |
| B | 1e-8 | 0 | 0 | -3.6e-15 | -3.6e-14 | -1.7e-13 | -1.9e-12 (-7.5e-14) |
| S | 1e-10 | 0 | 1.1e-14 | 1.1e-14 | -5.0e-14 | -3.9e-14 | 1.3e-12 (5.1e-14) |
| S | 1e-8 | 0 | 5.5e-12 | 5.5e-12 | -3.1e-12 | -1.9e-12 | -1.8e-11 (-7.2e-13) |
| C | 1e-10 | 0 | 5.4e-7 | 2.8e-7 | -3.1e-7 | -7.2e-7 | -3.2e-6 (-4.9e-12) |
| C | 1e-8 | 0 | -3.4e-5 | -7.5e-5 | -1.8e-4 | -2.7e-4 | -6.6e-4 (-1.0e-9) |

Total mass is the quantity least affected (it is set by the boundary flux and the corrected velocity). In C at 1e-8 the difference grows steadily (about -2e-11 relative per step) rather than staying at round-off.

Maximum temperature difference grows with time in every case (B 1e-10: 5.8e-11 K at step 5, 3.6e-8 K at step 49; S 1e-8: 1.6e-8 K at step 5, 9.5e-5 K at step 49; C 1e-8: 2.0e-6 K at step 5, 1.0e-3 K at step 60).

### 3.3 Divergence gate and cost
| Case | pred tol | worst fine div u - D (gate 1e-9) | gate | t_loop (s) | pressure (s) | pressure share | t_loop vs baseline |
|---|---|---|---|---|---|---|---|
| B | 1e-12 a / b | 3.59e-13 / 3.59e-13 | PASS | 18.56 / 18.54 | 8.83 / 8.85 | 47.6% / 47.7% | - |
| B | 1e-10 | 3.62e-13 | PASS | 18.07 | 8.36 | 46.2% | -2.6% |
| B | 1e-8 | 3.63e-13 | PASS | 17.58 | 7.87 | 44.8% | -5.2% |
| S | 1e-12 a / b | 7.26e-11 / 7.26e-11 | PASS | 88.24 / 88.13 | 54.49 / 54.44 | 61.8% / 61.8% | - |
| S | 1e-10 | 7.26e-11 | PASS | 84.75 | 51.15 | 60.4% | -3.9% |
| S | 1e-8 | 7.26e-11 | PASS | 81.70 | 47.74 | 58.4% | -7.4% |
| C | 1e-12 a / b | 3.98e-13 / 3.98e-13 | PASS | 3.77 / 3.75 | 1.37 / 1.36 | 36.3% / 36.4% | - |
| C | 1e-10 | 3.97e-13 | PASS | 3.66 | 1.26 | 34.6% | -2.9% |
| C | 1e-8 | 3.97e-13 | PASS | 3.57 | 1.17 | 32.8% | -5.4% |
The share is pressure time over t_loop (driver perf file). The mass and species balance gates of the runs also pass in all 12 runs. Baseline a/b times agree within 1.3% (these are not the timing-study runs of section 8; S is the 49-step version with the established flame, hence the larger pressure share than the 30-step runs of section 8).

## 4. Verdict
- **Not round-off.** Round-off is about 1e-13 relative. Even at the tighter 1e-10 predictor tolerance the maximum relative differences are 3e-10 (B), 1.4e-9 (S), 2.3e-8 (C) in TMP and RHO, i.e. 3 to 5 orders above round-off; at 1e-8 they are 3e-8 (B), 7e-8 (S), 4e-6 (C). Only the total mass difference in B and S reaches round-off-like levels (1e-15 to 1e-13 relative); in C it does not (5e-12 and 1e-9). The baseline repeat is bitwise identical, so the differences are caused by the tolerance and are deterministic.
- The divergence gate passes in every run with the same margin as the baseline (the corrector decides the final divergence, as in section 8.4). The gate therefore does **not** detect these state differences; it is the wrong instrument for this question.
- The differences are small (TMP at most 1.2e-3 K, 3.9e-6 relative, at 1e-8 on C after 60 steps; heat release 2.4e-7 relative) but they grow with time and with flow complexity, and they are largest in the case with the level-set fire (C). The effect is a state difference of a size that scales with the tolerance (about 50x to 170x from 1e-10 to 1e-8 in B, S, C): consistent with the predicted state carrying the pressure error into the density, temperature and species update.
- Cost: -2.6 to -3.9% t_loop at 1e-10 and -5.2 to -7.4% at 1e-8 (1 thread, 1 rank). This is a modest gain against a state difference that is not at round-off, so the objection to loosening the predictor tolerance by default is supported by the data. Neither value is recommended as a default here; the default is unchanged at 1e-12. If a looser predictor tolerance is ever proposed, it needs a justified-tolerance argument on the predicted state (differences in this table are the evidence) and a longer run on the level-set case than the 60 steps used here.

## 5. Limits
60 steps (C) and 49 steps (B, S) only; 1 rank, 1 thread; CPU build; B has no reaction (HRR rows empty); S and C only cover the species that react or are emitted in those inputs; the leaf comparison needs identical layouts, which held in all runs here (a looser tolerance could change a tag-driven layout in a longer run). The scratch patch (the two switches and the leaf dump) and the raw dumps stay in a private scratch directory on the test machine; nothing was pushed.
