# P2 findings: pressure (Poisson for H) with AMReX FFT::Poisson, composite MLMG and HYPRE

Scope: prototype P2 of the FDS-to-AMReX refactor. The checks follow `docs/pressure/01-amr-mapping-spec.md` §C.1, §C.4–C.5, §I and REC-I3, `docs/requirements.md` §2.1a (eps_H), FR-005 and FR-037/FR-039, and the refined acceptance list from the Pressure Solver Lead (2026-09-25).
- Code: `prototypes/p2_pressure_mlmg/`. The README there has the exact build and run commands.
- `scripts/run_all.sh` reruns every check and exits nonzero on failure.
- Logs: `prototypes/p2_pressure_mlmg/logs/`. `check_results.md` is the machine-generated version of the tables below.

Environment:
- AMReX 99ddfda (banner "26.09") from `(local AMReX install built against HYPRE 2.32)`: 3D, MPI+OMP, FFT via FFTW 3.3.10, Fortran interfaces.
- Linked **statically** against the FireX-pinned HYPRE `(local GNU third-party library tree)/libs/hypre/63331f19/lib/libHYPRE.a`. The CMake cache has `HYPRE_LIBRARIES` set to that `.a`, `ldd` shows no HYPRE shared object, and `HYPRE_BoomerAMGSolve` is resolved inside the executable.
- GNU 14.2 and OpenMPI 5.0.7, Release build.
- Runs used `mpirun --bind-to none`, at most 4 ranks, `OMP_NUM_THREADS=1` (one repeat set used 2 threads) and `OMP_DYNAMIC=false`.
- All runs were on 2026-09-25. Load average was 4–9 on 8 cores because other FDS jobs were running.
- The final full rerun of `run_all.sh` took 109 s and exited 0.

## 0. Summary

| # | Requirement | Result | Key numbers |
|---|---|---|---|
| 1 | FFT::Poisson vs MLMG (`setMaxOrder(2)`), single level, gauge = mean removal, rel. L2 ≤ eps_H (max-norm also reported) | **PASS** | Worst rel. L2 is 7.0e-14 (eps_H 1e-8 at N=64, 3.9e-8 at N=128). Worst max\|ΔH\| is 6.2e-16 (relative 7.5e-14). Covers Neumann, periodic and Dirichlet at 64³, Neumann at 128³ and a 34×18×32 domain. MLMG ‖r‖₂/‖b‖₂ ≤ 4.2e-13, which meets 01 §I.2 (≤ 1e-12) |
| 2 | C/F order ≥ 1.8 at maxorder 2 and 3 (ratio 2 and 4, ≥ 3 resolutions); C/F flux matching; compare with UGLMAT | **PASS** (order, flux). **UGLMAT comparison not done** | Minimum observed order (max-norm, C/F-adjacent cells) at maxorder 2: 1.94 (ratio 2), 1.93 (ratio 4), 1.92 (Dirichlet, ratio 2). At maxorder 3: 2.04 / 2.07 / 1.97. Global orders are ≥ 1.88. Flux identity (coarse C/F flux = fine-face sum inside the composite operator): ≤ 1.8e-14 relative. UGLMAT reference cannot be reproduced in P2 (§2.4) |
| 3 | Max H difference on 1/2/4 ranks × max_grid_size 32/16, per path | **PASS** | Worst rel. L2 across decompositions is 2.4e-16 and worst max\|ΔH\| is 4.1e-17, for every path (FFT, MLMG, MLMG+HYPRE bottom, HYPRE full solve). Results are not bitwise identical across layouts |
| 4 | 3 repeats per path at fixed ranks/threads/layout, OMP_DYNAMIC=false: bitwise? | **PASS (bitwise identical)** | All 4 paths, at 4 ranks × 1 thread and at 2 ranks × 2 threads, max_grid_size 16 |
| 5 | HYPRE version and per-solve times | recorded | `HYPRE (v2.32.0-24-g63331f19c - HEAD branch) initialized`. MLMG+HYPRE bottom agrees with the default bottom to ≤ 1.3e-15 rel. L2. Times: §5, indicative only |
| 6 | Fortran interface (`amrex_poisson` + `amrex_multigrid`) vs C++ | **PASS** | rel. L2 1.6e-16 vs C++ MLMG and 4.9e-14 vs FFT (not bitwise) |
| 7 | MLMG coarsest level for 34×18×32 | **confirmed** | 34×18×32 → 17×9×16 (2 MG levels) for max_grid_size 32, 16 and 8 |

Nothing is still running on the development machine. All P2 processes ended.

## 1. Requirement 1: FFT vs MLMG on a single level

**Setup.**
- RHS: `exp(-r²/0.02) + 0.5 sin(3πx) cos(2πy+0.4)(z+0.25)`, smooth and non-separable (`p2_main.cpp:rhs_func`).
- Singular BCs: the RHS mean is removed in application code (REC-C2).
- Solvers:
  - `FFT::Poisson`: BC `even` for Neumann, `odd` for Dirichlet, `periodic`.
  - `MLPoisson` + MLMG with `setMaxOrder(2)` (the default in the code), `tol_rel = 1e-12` and B-norm convergence (`setConvergenceNormType(MLMGNormType::bnorm)`).
- Gauge: volume-weighted mean removal after each singular solve.
- Metric: `‖ΔH‖₂/‖H_ref‖₂` after mean removal, with the FFT as reference. `max|ΔH|` and `max|ΔH|/max|H_ref|` are also reported.
- All runs: 4 ranks, max_grid_size 32 (34×18×32: 2 ranks).

| case | comparison | rel. L2 | max\|ΔH\| | rel. max | eps_H | result |
|---|---|---|---|---|---|---|
| Neumann 64³ | MLMG vs FFT | 4.93e-14 | 5.48e-16 | 4.78e-14 | 1e-8 | PASS |
| Neumann 64³ | MLMG+HYPRE bottom vs FFT | 4.93e-14 | 5.48e-16 | 4.78e-14 | 1e-8 | PASS |
| Neumann 64³ | HYPRE full solve (MLMG, max_coarsening_level=0) vs FFT | 9.5e-15 | 1.3e-16 | 1.1e-14 | 1e-8 | PASS |
| periodic 64³ | MLMG vs FFT | 4.00e-14 | 2.41e-16 | 3.66e-14 | 1e-8 | PASS |
| periodic 64³ | MLMG+HYPRE bottom vs FFT | 4.00e-14 | 2.41e-16 | 3.67e-14 | 1e-8 | PASS |
| Dirichlet (homog.) 64³ | MLMG vs FFT | 8.6e-15 | 5.5e-17 | 7.1e-15 | 1e-8 | PASS |
| Dirichlet 64³ | **MLMG maxorder 3** vs FFT (info) | **1.8e-3** | 1.6e-5 | 2.1e-3 | 1e-8 | differs, as expected |
| Neumann 128³ | MLMG vs FFT | 6.96e-14 | 6.21e-16 | 5.42e-14 | 3.9e-8 | PASS |
| Neumann 34×18×32 (mgs 32/16/8) | MLMG vs FFT | 1.29e-14 | 2.06e-16 | 1.99e-14 | 1e-8 | PASS |
| Neumann 34×18×32 | MLMG+HYPRE bottom vs FFT | 1.29e-14 | 2.1e-16 | 2.0e-14 | 1e-8 | PASS |

- **maxorder 2 is required for FFT/MLMG equivalence at Dirichlet (OPEN) faces.** The measurement confirms FR-039 and 01 §C.1: at maxorder 3 the Dirichlet ghost stencil differs, and H differs from the FFT by 1.8e-3. That is five orders of magnitude above eps_H.
- For homogeneous Neumann and periodic single-level problems, maxorder has no effect: the results are bitwise equal to maxorder 2.
- Final MLMG residuals ‖r‖₂/‖b‖₂ were 3e-14 to 4.2e-13. The 01 §I.2 condition (≤ 1e-12) is met, so the 1e-12 tolerance is reachable at N ≤ 128.

## 2. Requirement 2: coarse-fine order, flux matching, UGLMAT

**Setup** (`mode=mms`).
- Unit cube with a static two-level hierarchy, ratio 2 or 4.
- The fine patch covers half the domain in each direction. It is shifted off-centre by N/8 in y and z and sits away from the domain boundary.
- Manufactured solution: `sin(2πx+0.3) sin(2πy+0.7) sin(2πz+1.1) + 0.5 cos(4πx+0.2) cos(2π(y+z))`.
- The RHS is the point value of the analytic Laplacian.
- BCs:
  - Periodic (no boundary error, so only the C/F error is seen).
  - Inhomogeneous Dirichlet, with exact face values in the ghost cells via `setLevelBC`.
- Periodic runs: the composite volume-weighted RHS mean over uncovered cells is removed (REC-C2), then the RHS is averaged down onto covered coarse cells (01 §C.2).
- One composite `MLMG::solve` with `tol_rel = 1e-11`, 4 ranks, max_grid_size 32.
- Error metric: max-norm of (φ − φ_exact), with the composite mean removed for periodic runs. It is computed globally (uncovered coarse cells plus fine cells) and restricted to **C/F-adjacent cells**: fine cells touching the C/F boundary and uncovered coarse cells sharing a face with the patch.

### 2.1 Observed order (p = log₂(e_N/e_2N))

Max-norm, C/F-adjacent cells (both sides):

| bc | ratio | maxorder | N=16 | N=32 | N=64 | N=128 | p (16→32, 32→64, 64→128) | min p global | result |
|---|---|---|---|---|---|---|---|---|---|
| periodic | 2 | 2 | 3.22e-2 | 8.41e-3 | 2.14e-3 | 5.38e-4 | 1.94, 1.98, 1.99 | 1.94 | **PASS** |
| periodic | 2 | 3 | 2.37e-2 | 5.47e-3 | 1.31e-3 | 3.18e-4 | 2.11, 2.07, 2.04 | 1.99 | ≥ 1.8 |
| periodic | 4 | 2 | 3.40e-2 | 8.94e-3 | 2.28e-3 | – | 1.93, 1.97 | 1.93 | **PASS** |
| periodic | 4 | 3 | 2.36e-2 | 5.50e-3 | 1.31e-3 | – | 2.10, 2.07 | 1.99 | ≥ 1.8 |
| Dirichlet | 2 | 2 | 2.93e-2 | 7.77e-3 | 1.98e-3 | – | 1.92, 1.97 | 1.88 | **PASS** |
| Dirichlet | 2 | 3 | 2.19e-2 | 5.41e-3 | 1.38e-3 | – | 2.02, 1.97 | 1.94 | ≥ 1.8 |

Findings:
- **maxorder 2 meets the ≥ 1.8 target at the C/F interface for both ratio 2 and ratio 4**, and the order approaches 2 as the grid is refined. REC-I3 therefore accepts maxorder 2 for multi-level runs, and the FR-039 fallback to maxorder 3 is not needed.
- **Where the largest error sits.** At maxorder 2 the largest error is in the uncovered **coarse** cells next to the patch (err_max_all = err_max_CF).
- **Effect of maxorder 3.**
  - It lowers the C/F error by about 40% (2.14e-3 → 1.31e-3 at N=64, ratio 2) and the global max error by about 15%.
  - It raises the error in the fine cells next to the interface (7.8e-4 → 1.19e-3).
  - The L2 error is almost unchanged (6.43e-4 vs 6.33e-4).
- **Fine-side-only order.** Taken alone, the fine-side C/F error has order 1.77–1.84 on the coarsest pair (16→32) and 1.92–1.98 after that. The gate uses the maximum over both sides, as specified. The full table is in `logs/check_results.md`.
- **Tangential interpolation.** It stays at order 3 regardless of maxorder (hard-coded `IBD_max_order_DEF = 3`, per 01 §C.1). It was not varied here.
- **HYPRE bottom on two levels.** MLMG with the HYPRE bottom solver on the two-level case (N=64, ratio 2, maxorder 2) gives the same error as the default bottom: relative difference 0.0 at 7 significant digits.

### 2.2 Flux matching across the C/F interface

The fluxes are taken from `MLMG::getFluxes` after the composite solve (F = −∇φ per level). The fine face fluxes are area-averaged onto the coarse faces with `average_down_faces`, which gives (sum of fine fluxes)/rr². Three quantities were measured in every uncovered coarse cell next to the patch:

| quantity | meaning | result (max over all 21 two-level runs) |
|---|---|---|
| **flux identity** | MLMG's own composite residual `r = b − Lφ`, minus `(b + div F)` built with the **fine-face sum** on C/F faces, relative to max\|b\| | **≤ 1.8e-14 (machine precision) → PASS**. The composite operator uses exactly the fine-flux sum as the coarse C/F face flux |
| coarse-cell balance | `div F` (fine-sum fluxes on C/F faces) − b, relative | 2.8e-13 … 1.5e-10. This equals the composite residual, so it is set by the solver tolerance (tol_rel 1e-11) and is not a mismatch |
| raw per-level coarse flux vs fine average | `getFluxes` on the coarse level at C/F faces (computed from coarse data incl. averaged-down covered cells) vs fine average, relative to max\|F\| | 1.2e-1 (N=16) → 1.6e-2 (N=128), i.e. **first order, not zero** |

Consequence for the port:
- The coarse-level flux or face velocity that `MLMG::getFluxes` returns on C/F faces must **not** be used as it is.
- The projection must overwrite coarse C/F faces (and covered faces) with `average_down_faces` of the fine fluxes, as 01 §C.2 already prescribes (and `MacProjector` does). After that step, the interface mismatch is zero by construction and the coarse-cell divergence closes to the solver residual (above).

### 2.3 Singular composite RHS

In the periodic two-level runs, the composite RHS has a nonzero discrete mean before the application-level fix (−5.1e-3 at N=16, −7.9e-5 at N=128, O(h²)). Without the fix, MLMG still reported convergence in 9 iterations. The true composite residual, however, stayed at 2.6e-5·max\|b\| (N=16) and 6.3e-6·max\|b\| (N=32): MLMG's internal singular fix converges to a shifted problem. Evidence runs: `logs/info_nofix_periodic_rr2_mo2_n{16,32}.out` (`compat_fix=0`). This is evidence for REC-C2: remove the composite, volume-weighted RHS mean over uncovered cells in application code before every singular solve.

### 2.4 UGLMAT reference: not reproduced

REC-I2 criterion 2 asks for `‖H_AMR − H_exact‖ ≤ 1.1 × ‖H_UGLMAT − H_exact‖` on `Adaptive_Mesh_Refinement/ns2d_16_int_1to2_refinement` / `ns2d_16_emb_1to2_refinement`. That requires FDS baseline runs of derived `SOLVER='UGLMAT HYPRE'` inputs, which is the FR-016 V&V work (see R-34 on whether baseline even passes there). The spec defines no UGLMAT reference that can be reproduced on P2's manufactured 3-D problem. I also did not re-implement the UGLMAT two-point C/F flux, because that would not be the FDS reference. **As instructed, this item stops here: not done in P2.**

## 3. Requirement 3: decomposition

Setup: 64³ grid, Neumann and periodic BCs, 1/2/4 ranks × max_grid_size 32/16, which gives 6 runs per BC. Each run solves all four paths. Every run is compared with the 1-rank, max_grid_size 32 run of the same path. H is gauge-fixed and written as raw float64.

| bc | path | worst rel. L2 | worst max\|ΔH\| | worst rel. max | bitwise vs ref | result (eps_H 1e-8) |
|---|---|---|---|---|---|---|
| Neumann | FFT | 2.1e-16 | 6.9e-18 | 6.1e-16 | no (all 5) | PASS |
| Neumann | MLMG (default bottom) | 2.3e-16 | 5.2e-18 | 4.5e-16 | no | PASS |
| Neumann | MLMG + HYPRE bottom | 1.7e-16 | 3.5e-18 | 3.0e-16 | no | PASS |
| Neumann | HYPRE full solve | 1.8e-16 | 4.1e-17 | 3.6e-15 | no | PASS |
| periodic | FFT | 2.4e-16 | 3.5e-18 | 5.3e-16 | no | PASS |
| periodic | MLMG | 1.3e-16 | 1.7e-18 | 2.6e-16 | no | PASS |
| periodic | MLMG + HYPRE bottom | 1.2e-16 | 1.7e-18 | 2.6e-16 | no | PASS |
| periodic | HYPRE full solve | 1.2e-16 | 2.0e-17 | 3.0e-15 | no | PASS |

The differences are at round-off level. As requirements.md FR-005(iii) anticipates, none of the paths is bitwise invariant to rank count or box layout. That includes the FFT path, with differences from 1 rank/mgs 16 vs 1 rank/mgs 32 onward, since the gauge and RHS-mean reductions change order.

## 4. Requirement 4: bitwise repeat runs

Setup: 64³ Neumann, max_grid_size 16, `OMP_DYNAMIC=false`, 3 repeats per set. The table gives the md5 of the gauge-fixed H.

| set | FFT | MLMG | MLMG+HYPRE bottom | HYPRE full |
|---|---|---|---|---|
| 4 ranks × 1 thread | 767205d7f7dd ×3 | 23846706f158 ×3 | 47d226e07ae9 ×3 | 26481654cf8d ×3 |
| 2 ranks × 2 threads | 9a5c2061929b ×3 | b05dc2a0885c ×3 | 159aef8e8fae ×3 | 3c5bd4986855 ×3 |

**Bitwise identical run to run for every path.** The same md5 values also came out of three separate invocations of `run_all.sh`. This supports FR-005(iv) for the pressure solve at fixed ranks, threads and layout. As a negative control, flipping one value by 1e-15 in one repeat file makes `check.py` report FAIL and exit 1.

## 5. HYPRE version and per-solve times

- **HYPRE version actually used:** `HYPRE (v2.32.0-24-g63331f19c - HEAD branch) initialized`. This is printed by AMReX at init (`AMReX.cpp:771-781`) in every run log, and `check.py` checks it. Together with the static link above and the agreement checks, this closes **R-30**.
- **HYPRE was really called.** With `bottom_verbose=2` (`logs/r1_hypre_verbose.out`), HYPRE reported "solver = BoomerAMG; preconditioner = none" and 16 `HYPRE BoomerAMG: Num. iterations = …` lines, one per bottom solve or full solve.
- **HYPRE vs default bottom.** MLMG+HYPRE bottom vs MLMG default bottom: rel. L2 ≤ 1.9e-16 on power-of-two grids and ≤ 1.3e-15 on 34×18×32. On power-of-two grids with agglomeration, MLMG coarsens to 2×2×2, so HYPRE only solves 8 unknowns. On 34×18×32, HYPRE solves the 17×9×16 = 4896-cell bottom.
  - To give HYPRE a real test, the "HYPRE full solve" path (`max_coarsening_level=0` + HYPRE bottom) was added. It agrees with FFT to 9.5e-15 (Neumann) and 9.5e-14 (Dirichlet).

**Per-solve wall time.** These are **indicative only**: the development machine was shared with other FDS jobs (load 4–9 on 8 cores). Runs used 4 ranks × 1 thread and 64³ unless noted. FFT times include `FFT::Poisson` construction (planning). MLMG times are `MLMG::solve` only; operator setup is excluded.

| case | FFT | MLMG (default bottom) | MLMG + HYPRE bottom | HYPRE full solve |
|---|---|---|---|---|
| Neumann 64³ | 0.0051 s | 0.032 s (12 it) | 0.030 s (12 it) | 1.60 s (4 MLMG it) |
| periodic 64³ | 0.0033 s | 0.030 s (11 it) | 0.029 s (11 it) | 2.95 s (5 it) |
| Dirichlet 64³ | 0.0040 s | 0.031 s (12 it) | 0.030 s (12 it) | 0.93 s (3 it) |
| Neumann 128³ | 0.041 s | 0.26 s (13 it) | 0.26 s (13 it) | – |
| Neumann 34×18×32 (2 ranks) | 0.0017 s | 0.013 s (9 it) | 0.054 s (9 it) | – |

The values are from the final rerun. Across the three full reruns, the same solve varied by up to about 3×; for example, 128³ MLMG+HYPRE took 0.73 s in one rerun and 0.26 s in another.

MLMG costs about 6–9× the FFT per solve here. That is consistent with the 02 §7.6 estimate (FFT 2–8× cheaper), but it must be re-measured on an idle machine (NFR-030 rules).

## 6. Fortran interface (feasibility)

`src/p2_fortran.F90` runs the same 64³ all-Neumann solve through `amrex_poisson` (`set_maxorder(2)`, `set_domain_bc`, `set_level_bc`) and `amrex_multigrid` (`set_always_use_bnorm`, `solve`, `comp_residual`). It uses the same RHS formula, compatibility fix and gauge as the C++ code, on 1 rank with mgs 32 and on 2 ranks with mgs 16.

| Fortran run | vs C++ MLMG (1 rank, mgs 32) | vs C++ FFT | vs C++ MLMG (2 ranks, mgs 16) | ‖r‖₂/‖b‖₂ |
|---|---|---|---|---|
| 1 rank, mgs 32 | 1.6e-16 | 4.9e-14 | 2.5e-16 | 2.95e-13 |
| 2 ranks, mgs 16 | 1.6e-16 | 4.9e-14 | 2.2e-16 | 2.95e-13 |

All PASS. The results are not bitwise equal to C++ at the same layout. The likely cause, not proven, is that the RHS is evaluated with gfortran `exp`/`sin` rather than the C++ libm calls. Two pitfalls:
- `comp_residual`/`solve` need array arguments with INTENT(INOUT), so array constructors like `[phi]` don't compile.
- `amrex_multifab%norm0/norm2/sum` are collective, so calling them inside an `if (amrex_parallel_ioprocessor())` block deadlocks. My first version hung this way.

## 7. MLMG coarsening for 34×18×32

The MG hierarchy printed from `MLPoisson::NMGLevels`/`Geom(0,m)` is **34×18×32 → 17×9×16**, i.e. 2 MG levels, for max_grid_size 32, 16 and 8. MLMG stops at 17×9×16 because 17 and 9 are odd and domain-based coarsening needs the whole domain to coarsen by 2. The bottom solve is therefore a 4896-cell problem: bicgstab by default, or HYPRE. For comparison, 64³ and 128³ coarsen to 2×2×2 (6 and 7 levels). Iteration counts were not affected here (9 iterations for 34×18×32 vs 12 for 64³). Poorly factorable domain sizes still raise the bottom-solver share, which matters for 01 REC-D3 (tunnels). HYPRE as the bottom solver is the planned mitigation.

## 8. Surprises and pitfalls

1. **BoomerAMG used as the whole solver diverges on singular (all-Neumann or periodic) matrices.** MLMG with `max_coarsening_level=0` + HYPRE bottom failed after 1 iteration with resid/bnorm = 2e+117 (64³) and 3e+110 (32³), then aborted with "MLMG failing so lets stop here". Evidence: `logs/info_hypre_full_nosingfix.out`, which `run_all.sh` reruns and expects to fail. It works only with `hypre.adjust_singular_matrix=1`, which pins one row (`AMReX_HypreABecLap3.cpp:486`). PCG+BoomerAMG (`hypre.hypre_solver=PCG hypre.hypre_preconditioner=BoomerAMG`, close to the FDS UGLMAT setup) also failed without it, stalling at resid/bnorm 0.31 after 200 MLMG iterations, and works with it (4 iterations, agreeing with FFT to 9.5e-15). That was an ad-hoc run; it is not in `run_all.sh` and no log was kept. For HYPRE used as the bottom of a normal MLMG hierarchy on the 2×2×2 or 17×9×16 bottom, the flag was not needed, but `run_all.sh` sets it for every run. **Any HYPRE-full reference configuration (01 REC-G2.2) must set `hypre.adjust_singular_matrix=1` for closed domains.**
2. **Raw coarse fluxes at C/F faces are first-order wrong.** `MLMG::getFluxes` computes each level's fluxes independently. At C/F faces the coarse flux is not the fine sum (§2.2), so `average_down_faces` is mandatory in the projection.
3. **MLMG converges happily on an incompatible singular RHS.** Its internal fix shifts the problem instead (§2.3), so the application must remove the RHS mean.
4. **`ParmParse::queryarr` never shrinks a pre-filled `std::vector`.** It only grows it (`AMReX_ParmParse.cpp:1313-1316`), so a default list longer than the input keeps its tail. This bit my `solvers=` option: a 2-element input ran 5 solvers. The driver must start from an empty vector.
5. `MLMG::setAlwaysUseBNorm` is deprecated in 99ddfda. Use `setConvergenceNormType(MLMGNormType::bnorm)`. The Fortran interface still exposes only `set_always_use_bnorm`.
6. The FFT agrees with MLMG at round-off level (≤ 7e-14) even though the FFT is a direct solve. So the FFT→MLMG switch in FR-039 has no measurable H jump on single-level input.

## 9. What remains

- **UGLMAT comparison (REC-I2 criterion 2):** needs FDS UGLMAT-HYPRE runs of the `ns2d_16_*` cases (FR-016 V&V). Not done in P2.
- **Timings on an idle machine** (NFR-030 rules; TinyProfiler regions `MLMG::solve()` and `FFT::Poisson::solve`). The numbers above are indicative only.
- **Not covered by P2:**
  - frozen-state FDS input (01 §H.2), i.e. an RHS taken from an FDS step;
  - inhomogeneous Dirichlet/OPEN values lifted into the FFT RHS (02 §7.2);
  - mixed faces (Robin/per-cell BCs);
  - obstructions (E-1/E-2) and per-zone gauge;
  - `MacProjector` (AMReX-Hydro) and the `∇·U = S` sign check (01 §A.1 [VERIFY]);
  - 3-level hierarchies;
  - tangential C/F interpolation order other than 3.
- **Decomposition runs at 8 ranks** and V&V's byte-identical `SIG_FIGS=17` checks (FR-005(i)) are out of P2 scope; the development machine limit tonight was 4 ranks.
- `runs/` holds 149 MB of raw fields. They can be deleted, and `run_all.sh` regenerates them.
