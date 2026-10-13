# M3: pyramid solver cost table (P4 + docs/pressure/13 section 3 rules)

Solver-only: P1 velocity is not in yet, so the right-hand side is a fixed analytic field, not physics. Harness: `pb_harness mode=comp_pyramid n=48 mgs=16 cost=1 nrep=5`
(`run_pyr_cost`, harness/composite_modes.cpp). Case: the pyramid_48_4m geometry. L0 48^3, L1 ratio 2 over the pyramid, three-step OBST at odd fine
indices, periodic x/y, Neumann z, max_grid_size 16. Masked composite arms share one rhs, mask, gauge (mean removed) and layout. FFT is not applicable to a
masked or composite problem, so FFT and HYPRE PCG are timed on the unmasked single-level L0 baseline (same L0 rhs), "where applicable" per the rules.

Rules followed: tol_rel 1e-10 for every arm. Fresh workspace per arm per rep: first solve with zero guess, then re-solve with zero guess, then warm
solve (previous solution as initial guess for an rhs perturbed by 1%, standing in for the next step). Arms are interleaved within each rep. Rep 0 is an
untimed warm-up, then median of 5. Times are the max over ranks of the backend's own clocks. Setup: HYPRE = IJ assembly plus BoomerAMG setup, FFT = plan.
MLMG = first solve minus re-solve (rules item 4); MLMG rebuilds its operator on every solve (no workspace cache), so its re-solve includes the build.
Amortized per step = setup/N + 2 x re-solve. True residual is the gauge-free check with pin rows excluded.

OBST-face divergence max: 0 by construction for every masked arm. The masked operator has no solid-gas faces (Neumann closure), and there is no
velocity field yet (P1) to take a divergence of. To be revisited with P1.

GPU rows: skipped. This build has no CUDA-enabled HYPRE or AMReX GPU pressure path.

P4: MLMG-preconditioned BiCGStab and CG (single-level symmetric) at tol 1e-10 are covered by ctest `pb_krylov` (single-level MLMG CG and BiCGStab,
composite BiCGStab, composite CG = InvalidInput). The masked composite BiCGStab check is in `pb_comp_masked_mlmg`.

### CPU, 1 rank

| Arm | Method | Iters (zero / warm) | Setup ms | First solve ms | Re-solve ms | Warm solve ms | Amortized N=1 / 10 / 100 ms | True rel residual |
|---|---|---|---|---|---|---|---|---|
| FFT (single-level unmasked L0) | FFT | 0 / 0 | 0.4 | 3.4 | 3.0 | 2.8 | 6.5 / 6.1 / 6.1 | 5.2e-14 |
| HYPRE PCG+BoomerAMG (single-level unmasked L0) | PCG | 18 / 15 | 1037.8 | 1296.0 | 185.1 | 197.4 | 1408.1 / 474.0 / 380.6 | 6.8e-12 |
| HYPRE Auto Krylov+BoomerAMG (masked composite) | GMRES | 18 / 17 | 2506.5 | 3100.2 | 553.9 | 396.1 | 3614.2 / 1358.4 / 1132.8 | 1.2e-11 |
| HYPRE GMRES(30)+BoomerAMG (masked composite) | GMRES | 18 / 17 | 2209.0 | 2571.0 | 478.0 | 580.7 | 3165.0 / 1176.9 / 978.1 | 1.2e-11 |
| BoomerAMG alone (masked composite) | BoomerAMG | 44 / 37 | 2155.7 | 3112.0 | 1349.1 | 780.5 | 4853.9 / 2913.8 / 2719.8 | 1.0e-11 |
| MLMG V-cycle (masked composite) | V-cycle | 15 / 12 | 95.5 | 501.3 | 471.7 | 348.9 | 1038.8 / 952.9 / 944.3 | 1.8e-13 |
| MLMG-preconditioned BiCGStab (masked composite) | BiCGSTAB | 5 / 6 | 177.6 | 550.7 | 531.6 | 796.8 | 1240.8 / 1081.0 / 1065.0 | 3.2e-12 |

Peak RSS (process high-water mark over all arms): max per rank 221 MB, sum 221 MB.

### CPU, 4 ranks (box has 1 core: oversubscribed, times NOT valid; iterations and residuals are)

| Arm | Method | Iters (zero / warm) | Setup ms | First solve ms | Re-solve ms | Warm solve ms | Amortized N=1 / 10 / 100 ms | True rel residual |
|---|---|---|---|---|---|---|---|---|
| FFT (single-level unmasked L0) | FFT | 0 / 0 | 0.5 | 2.5 | 6.1 | 5.1 | 12.7 / 12.3 / 12.2 | 5.2e-14 |
| HYPRE PCG+BoomerAMG (single-level unmasked L0) | PCG | 18 / 15 | 1253.6 | 2563.0 | 841.7 | 281.2 | 2937.0 / 1808.7 / 1695.9 | 6.4e-12 |
| HYPRE Auto Krylov+BoomerAMG (masked composite) | GMRES | 19 / 17 | 2706.0 | 3888.6 | 999.3 | 2741.7 | 4704.6 / 2269.3 / 2025.7 | 4.5e-12 |
| HYPRE GMRES(30)+BoomerAMG (masked composite) | GMRES | 19 / 17 | 2765.7 | 4297.3 | 2030.7 | 1897.9 | 6827.1 / 4338.1 / 4089.1 | 4.5e-12 |
| BoomerAMG alone (masked composite) | BoomerAMG | 44 / 37 | 2878.7 | 6234.1 | 2641.1 | 1671.7 | 8160.9 / 5570.0 / 5311.0 | 1.3e-11 |
| MLMG V-cycle (masked composite) | V-cycle | 15 / 12 | 0.0 | 2536.1 | 4492.2 | 3335.5 | 8984.5 / 8984.5 / 8984.5 | 1.3e-13 |
| MLMG-preconditioned BiCGStab (masked composite) | BiCGSTAB | 5 / 5 | 45.6 | 1676.4 | 2242.4 | 3559.5 | 4530.4 / 4489.3 / 4485.2 | 3.2e-12 |

Peak RSS (process high-water mark over all arms): max per rank 98 MB, sum 367 MB.

## Reading

- On 1 rank, MLMG V-cycle is the cheapest masked-composite arm at every regrid interval: about 0.94 s per step versus 0.98 s for HYPRE GMRES at
  N=100, and versus 3.2 to 3.6 s for HYPRE at N=1. HYPRE's 2.2 to 2.5 s setup dominates unless the grids live for 10 or more steps.
- MLMG-preconditioned BiCGStab cuts iterations from 15 to 5 but not time (2 V-cycle applications per iteration).
- BoomerAMG alone needs 44 iterations, about 2.5x the Krylov-wrapped arms.
- Warm start saves 2 to 7 iterations; it does not change the ranking.
- The unmasked single-level FFT baseline costs about 6 ms per step, roughly 150x below any masked composite arm.
- The 4-rank rows were run on a 1-core box and are oversubscribed. Valid 4-rank timings need a multi-core machine (pop-os is approved).
