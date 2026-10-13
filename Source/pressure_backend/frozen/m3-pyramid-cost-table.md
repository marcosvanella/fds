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

### Box (1 core, idle: load 1.7), 1 rank

| Arm | Method | Iters (zero / warm) | Setup ms | First solve ms | Re-solve ms | Warm solve ms | Amortized N=1 / 10 / 100 ms | True rel residual |
|---|---|---|---|---|---|---|---|---|
| FFT (single-level unmasked L0) | FFT | 0 / 0 | 0.4 | 3.3 | 2.8 | 2.7 | 5.9 / 5.6 / 5.5 | 5.2e-14 |
| HYPRE PCG+BoomerAMG (single-level unmasked L0) | PCG | 18 / 15 | 661.3 | 828.5 | 164.9 | 142.7 | 991.1 / 395.9 / 336.4 | 6.8e-12 |
| HYPRE Auto Krylov+BoomerAMG (masked composite) | GMRES | 18 / 17 | 1619.8 | 1951.5 | 333.9 | 319.6 | 2287.6 / 829.8 / 684.0 | 1.2e-11 |
| HYPRE GMRES(30)+BoomerAMG (masked composite) | GMRES | 18 / 17 | 1607.4 | 1936.5 | 318.6 | 302.3 | 2244.6 / 798.0 / 653.3 | 1.2e-11 |
| BoomerAMG alone (masked composite) | BoomerAMG | 44 / 37 | 1606.0 | 2234.6 | 663.6 | 540.1 | 2933.1 / 1487.7 / 1343.2 | 1.0e-11 |
| HYPRE PCG+BoomerAMG, symmetric C/F (masked composite, hypre.symmetric_cf) | PCG | 18 / 17 | 1610.6 | 1877.2 | 276.3 | 255.6 | 2163.2 / 713.6 / 568.7 | 1.1e-11 |
| MLMG V-cycle (masked composite) | V-cycle | 15 / 12 | 11.7 | 371.6 | 363.3 | 278.5 | 738.3 / 727.8 / 726.7 | 1.8e-13 |
| MLMG-preconditioned BiCGStab (masked composite) | BiCGSTAB | 5 / 6 | 29.6 | 359.0 | 329.9 | 368.7 | 689.4 / 662.8 / 660.1 | 3.2e-12 |

Peak RSS (process high-water mark over all arms): max per rank 223 MB, sum 223 MB.
### pop-os (i9-13900HX, 24 cores / 32 threads, P+E hybrid; shared, load average 19-21 during the run), 1 rank

| Arm | Method | Iters (zero / warm) | Setup ms | First solve ms | Re-solve ms | Warm solve ms | Amortized N=1 / 10 / 100 ms | True rel residual |
|---|---|---|---|---|---|---|---|---|
| FFT (single-level unmasked L0) | FFT | 0 / 0 | 0.4 | 4.3 | 3.6 | 3.6 | 7.6 / 7.2 / 7.1 | 5.2e-14 |
| HYPRE PCG+BoomerAMG (single-level unmasked L0) | PCG | 18 / 15 | 805.8 | 1255.8 | 442.1 | 375.7 | 1690.0 / 964.9 / 892.3 | 6.8e-12 |
| HYPRE Auto Krylov+BoomerAMG (masked composite) | GMRES | 18 / 17 | 2006.7 | 2823.7 | 837.9 | 787.1 | 3682.5 / 1876.4 / 1695.8 | 1.2e-11 |
| HYPRE GMRES(30)+BoomerAMG (masked composite) | GMRES | 18 / 17 | 2001.8 | 2844.8 | 837.3 | 737.5 | 3676.4 / 1874.8 / 1694.6 | 1.2e-11 |
| BoomerAMG alone (masked composite) | BoomerAMG | 44 / 37 | 1993.4 | 3776.0 | 1718.9 | 1455.6 | 5431.1 / 3637.1 / 3457.7 | 1.0e-11 |
| HYPRE PCG+BoomerAMG, symmetric C/F (masked composite, hypre.symmetric_cf) | PCG | 18 / 17 | 1957.1 | 2662.6 | 716.4 | 681.2 | 3389.9 / 1628.5 / 1452.3 | 1.1e-11 |
| MLMG V-cycle (masked composite) | V-cycle | 15 / 12 | 64.8 | 683.5 | 618.7 | 502.2 | 1302.2 / 1243.8 / 1238.0 | 1.8e-13 |
| MLMG-preconditioned BiCGStab (masked composite) | BiCGSTAB | 5 / 6 | 44.4 | 679.8 | 638.5 | 713.0 | 1321.5 / 1281.5 / 1277.5 | 3.2e-12 |

Peak RSS (process high-water mark over all arms): max per rank 227 MB, sum 227 MB.
### pop-os, 4 ranks (--bind-to core)

| Arm | Method | Iters (zero / warm) | Setup ms | First solve ms | Re-solve ms | Warm solve ms | Amortized N=1 / 10 / 100 ms | True rel residual |
|---|---|---|---|---|---|---|---|---|
| FFT (single-level unmasked L0) | FFT | 0 / 0 | 0.4 | 2.5 | 1.8 | 2.0 | 4.0 / 3.7 / 3.7 | 5.2e-14 |
| HYPRE PCG+BoomerAMG (single-level unmasked L0) | PCG | 18 / 15 | 425.0 | 587.3 | 166.0 | 145.0 | 757.0 / 374.5 / 336.2 | 6.4e-12 |
| HYPRE Auto Krylov+BoomerAMG (masked composite) | GMRES | 19 / 17 | 922.4 | 1261.1 | 338.7 | 327.5 | 1599.9 / 769.7 / 686.7 | 4.5e-12 |
| HYPRE GMRES(30)+BoomerAMG (masked composite) | GMRES | 19 / 17 | 911.3 | 1270.2 | 345.9 | 316.0 | 1603.0 / 782.9 / 700.8 | 4.5e-12 |
| BoomerAMG alone (masked composite) | BoomerAMG | 44 / 37 | 901.2 | 1606.4 | 728.4 | 627.2 | 2358.1 / 1547.0 / 1465.9 | 1.3e-11 |
| HYPRE PCG+BoomerAMG, symmetric C/F (masked composite, hypre.symmetric_cf) | PCG | 19 / 17 | 895.9 | 1199.2 | 313.8 | 282.2 | 1523.6 / 717.3 / 636.6 | 3.7e-12 |
| MLMG V-cycle (masked composite) | V-cycle | 15 / 12 | 27.5 | 424.3 | 400.8 | 320.8 | 829.2 / 804.4 / 802.0 | 1.3e-13 |
| MLMG-preconditioned BiCGStab (masked composite) | BiCGSTAB | 5 / 5 | 25.2 | 360.7 | 339.3 | 349.2 | 703.9 / 681.2 / 678.9 | 3.2e-12 |

Peak RSS (process high-water mark over all arms): max per rank 96 MB, sum 360 MB.

## Symmetric composite operator (hypre.symmetric_cf, selectable; default unchanged)

FDS multi-mesh style C/F coupling: every fine C/F face carries the flux (phi_f - phi_C)/h_eff, h_eff = (r+1)/2 h_f. The fine row and the coarse
(reflux) row use the same weight, with no interpolation slopes. The volume-scaled matrix passes the symmetry check (||A - A^T||_inf < 1e-13 ||A||_inf),
so Auto Krylov picks PCG. Accuracy cost on the pyramid (48^3, 1 rank; `comp_pyramid n=48 mgs=16 symcmp=1`). Differences are mean-removed and
taken over uncovered gas cells, relative to max|ref| and ||ref||_2 globally. "Near C/F" = cells within 2 cells of the interface on their level.

| Pair | max (global) | L2 (global) | max (near C/F) | L2 (near C/F) |
|---|---|---|---|---|
| symmetric vs current HYPRE (non-symmetric) | 7.9e-3 | 2.9e-3 | 7.9e-3 | 1.8e-3 |
| symmetric vs MLMG | 7.9e-3 | 2.9e-3 | 7.9e-3 | 1.8e-3 |
| current HYPRE vs MLMG (reference) | 4.9e-11 | 1.8e-12 | 2.6e-12 | 5.1e-13 |

The 0.8 % maximum difference sits at the C/F interface; this is first-order coupling versus the second-order (slope) interpolation of MLMG and
default HYPRE. Cost: PCG 18 iterations, re-solve 10-17 % cheaper than GMRES, the same setup. Test: `pb_comp_masked_mlmg` (np 1 and 2:
PCG chosen, max diff vs MLMG < 2 %).

## Reading

- Box, 1 rank (idle machine): MLMG-preconditioned BiCGStab is the cheapest masked-composite arm when grids change often: 0.69 / 0.66 s per step
  at N = 1 / 10, versus 0.74 / 0.73 s for the MLMG V-cycle and 2.2 / 0.80 s for HYPRE GMRES. At N = 100, HYPRE wins: PCG with the symmetric
  C/F at 0.57 s and default GMRES at 0.65 s, versus BiCGStab at 0.66 s. HYPRE's 1.6 s setup only pays off when grids live for about 50 steps or more.
- MLMG-preconditioned BiCGStab cuts iterations from 15 to 5. Its time is close to the V-cycle (2 V-cycle applications per iteration).
- BoomerAMG alone needs 44 iterations, about 2.5x the Krylov-wrapped arms.
- Warm start saves 2 to 7 iterations; it does not change the ranking.
- The unmasked single-level FFT baseline costs about 6 ms per step, roughly 150x below any masked composite arm.
- 4 ranks (pop-os): every HYPRE arm speeds up by about 2.2-2.5x versus pop-os 1 rank. MLMG V-cycle speeds up only 1.5x and BiCGStab 1.9x.
  At 4 ranks MLMG-preconditioned BiCGStab is the cheapest masked arm at N = 1 and 10: 0.70 / 0.68 s per step, versus 1.52 / 0.72 s for the
  best HYPRE arm. At N = 100 they are about equal: HYPRE PCG with the symmetric C/F at 0.64 s (with the 0.8 % accuracy cost above), BiCGStab
  at 0.68 s, default HYPRE GMRES at 0.69 s.
- pop-os was shared (load 19-21 from other jobs) on hybrid P/E cores, so its absolute times are pessimistic. Its 1-rank times are slower than
  the idle box. The ratios between arms in one run are the reliable part (interleaved, median of 5).
