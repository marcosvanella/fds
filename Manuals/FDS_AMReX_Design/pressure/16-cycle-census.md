# 16 — Cycle census of the composite pressure solve, and what drives the cycle count

Status: **MEASURED** (iteration counts only; no timings), tip `bea3ec428a` plus a **private scratch patch that is not committed and changes no default** (environment switches, section 2). 1 rank, 1 OpenMP thread, GCC 13, Open MPI 5, HYPRE 2.32, AMReX 26.09, CPU. Nothing pushed. Counts do not depend on the machine.

> **Erratum 2026-10-21 (`pb_mlmg_bottom.*` keys).** Root cause confirmed in source: `pressure_flag` (`Source/driver/main.cpp` 57) forwarded only `hypre.*` and `pressure.*` tokens into ParmParse (since `25640298df`, 2026-10-08), and `amrex::Initialize(argc, argv, false, ...)` (main.cpp 159, since `004b593534`) does not parse argv, so `pb_mlmg_bottom.{solver,bamg_*,adjust_singular_matrix,recompute_preconditioner}`, read only by the masked composite MLMG (`CompositeSolve.cpp` `build_masked_mlmg` 1078-1095; keys present since `0b0c36d823`/`d2c134e36d`), were never set from the command line. Consequence for doc 14: the bottom smoother/BiCGStab rows of §7.4 ran the default BoomerAMG bottom (outputs byte-identical to the baseline); they are marked invalid there. Fix (src, 2026-10-21): `pb_mlmg_bottom.` added to the forwarded prefixes, and a rank-0 line `masked composite MLMG: bottom solver <name> (pb_mlmg_bottom.solver)` printed when the key is set; checked on the box (burner 3-level case, 2 steps): the line appears and the multi-level CSV now differs from the default run. Default behaviour unchanged.

**Where it ran.** The owner-provided NVIDIA test machine was listed and connected, but it holds no built tree of this tip with the matching HYPRE/AMReX configuration, so every run here was made on the project's Linux machine (8 vCPU, CPU only). Threads were not pinned because up to six runs shared the box; this only affects time, which is not reported here. Re-run on the test machine if wanted: the counts must be identical (single rank, deterministic).

## 1. Summary
- Default composite path for every masked case here (OBST plus OPEN vent): **masked composite MLMG**, plain V(2,2) cycles to `tol_rel` 1e-12 (stops at the own relative residual about 1e-13), **bottom solver HYPRE BoomerAMG** (standalone BoomerAMG iterations, relative tolerance 1e-4, maximum 200), one bottom call per V-cycle. BiCGStab is **not** the configured bottom of the masked path. It is used only (a) by the masked path when it sees thin faces that seal several regions (`sealed` case below), and (b) by the **unmasked** composite MLMG (no OBST; AMReX default bottom). The BiCGStab fallback after a failed HYPRE bottom never ran in these runs (`cg_calls` 0 in every masked-HYPRE run).
- Reconciliation of 2.2 vs 20.2 cycles per solve (case B): both are the same solver and the same operator. 20.2 is the mean over composite solves of a run that keeps its fine levels (doc 14's B: `--nofine -1`). 2.2 is the mean over **all** pressure rows of a run whose fine levels vanish after step 5 (the burner vent stays on level 0 by default): 10 composite solves of 22 cycles plus 88 level-0 FFT solves that carry 0 iterations, (5 x 22) / 49 = 2.24 per stage. Section 3.1.
- Ranked drivers of the cycle count (section 5): (1) number of AMR levels and how much of the domain the fine levels cover (14.9 to 20.7 to 29.8-34 cycles); (2) cell aspect ratio (1:1 20.7, 2:1 26, 4:1 diverges); (3) requested tolerance (1.55 cycles per decade); then nothing else measurable. The OBST mask, the coarsened mask/beta, the bottom solver choice and the domain shape do **not** change the cycle count.

## 2. Method and what was counted
- Scratch patch (not committed): `FDSTL_CC_VERB` (>= 3: MLMG iteration history, HYPRE bottom prints iterations per call; prints the end-of-solve line with path, cycles, number of BiCGStab bottom calls and iterations, own residual, tolerance), `FDSTL_PRED_TOL` / `FDSTL_CORR_TOL` (per-stage `tol_rel` and `residual_tol`; the same hook as in doc 15), `FDSTL_CC_BOTTOM` (smoother | bicgstab | hypre for the masked path). Note: the driver only forwards `hypre.*` and `pressure.*` keys, so `pb_mlmg_bottom.solver=...` on the command line is **ignored** at this tip (checked: identical output). The "bottom smoother / bicgstab: no change" row of doc 14 section 7.4 was a no-op: that run's scratch copy had the same filter (confirmed 2026-10-21, see erratum above).
- V-cycles per solve: number of `MLMG: Iteration` lines. Initial residual: `Initial residual (resid0)`, always equal to the rhs norm (zero guess, ratio 1). Final/initial: last fine-level `resid/bnorm`. Convergence factor per cycle: **geo** = (final/initial)^(1/cycles); **asym** = geometric mean over the last 3 cycles.
- **True bottom calls versus fine-level smoothing.** AMReX does not coarsen AMR levels above 0, so a fine AMR level has a single MG level and `mgVcycle_bottom` there is only one smoothing pass (pre-smooth count, 2 sweeps), not a solve. Those passes are **not** counted below: `bottom` columns count only the AMR-level-0 coarsest-grid solve (9x6x6 for the 36x24x24 cases). It runs exactly once per V-cycle, so bottom calls per solve equal the V-cycles per solve (the printed per-call HYPRE and BiCGStab lines confirm it in every run). The TinyProfiler timer `MLMG::mgVcycle_bottom` of doc 14 encloses both.
- Outer Krylov: none on the masked MLMG path (`Krylov::Auto` runs plain V-cycles). It exists on the HYPRE composite path (GMRES) and in the single-level study paths (section 3.3).
- 2 solves per step: predictor (`pred`) and corrector (`corr`). The `pred` count includes the projection solves at initialisation and regrids, so it exceeds the `corr` count by a few.
- One run per row (deterministic; the baseline repeats of docs 14/15 were bit-identical). Cycles are shown as mean/min/max over all solves of the run.

## 3. Census
### 3.1 Cases as run (default tolerance 1e-12, no study flags)
Columns: solves (n pred / n corr); V-cycles pred; V-cycles corr; bottom iterations per bottom call (pred and corr pooled, mean/min/max); final/initial residual (mean); factor geo / asym (pred; corr within 0.01 unless noted); path and bottom.

| Case (input, flags) | n | V-cycles pred | V-cycles corr | bottom its/call | final/initial | factor geo / asym | Path, bottom |
|---|---|---|---|---|---|---|---|
| **B** near-burner 3-level, doc 14 flags: `burner_obst_3level.fds`, `--levels 3 --patch 0 7 0 7 --patch2 0 7 0 3 --nofine -1 --regrid-every 5 --tag-obst --tag-temp 400 --blocking 4 --maxsize 8`, 49 steps (54/180 boxes at start, 12/24 after step 5) | 53/49 | 20.21 / 20 / 22 | 20.18 / 20 / 22 | 3.02 / 2 / 5 | 5.6e-14 | 0.220 / 0.227 | masked MLMG, HYPRE BoomerAMG |
| **B (Role 2 flags)**: same without `--nofine -1` (61/144 boxes at start, **no fine boxes after step 5**) | 5/5 composite (+88 FFT level-0 solves, 0 iterations) | 22.00 / 22 / 22 | 22.00 / 22 / 22 | 2.81 / 2 / 5 | 7.4e-14 | 0.253 / 0.239 | masked MLMG, HYPRE BoomerAMG |
| **S** `run_simple_test_l3.sh` (3 levels, tag-vort auto, regrid 5), steps 1-5 (64 L1 + 1,872 L2 boxes) | 5/5 | 31.60 / 29 / 33 | 34.00 / 33 / 35 | 3.97 / 3 / 4 | 6.4e-14 | 0.382 / 0.420 | masked MLMG, HYPRE BoomerAMG |
| **S** same, steps 6-49 (3-17 L1 boxes, up to 2 levels) | 49/39 | 13.27 / 11 / 15 | 13.38 / 11 / 15 | 3.07 / 3 / 5 | 3.8e-14 | - / 0.104 | masked MLMG, HYPRE BoomerAMG |
| **S** same, whole 49 steps | 54/44 | 14.96 / 11 / 33 | 15.73 / 11 / 35 | 3.2 / 3 / 5 | 3.8e-14 | 0.122 / 0.133 | masked MLMG, HYPRE BoomerAMG |
| S 2-level (`run_simple_test_l2.sh`), 10 steps (steps 1-5 composite) | 5/5 | 30.80 / 28 / 33 | 32.60 / 32 / 33 | 3.97 / 3 / 4 | 7.4e-14 | 0.374 / 0.410 | masked MLMG, HYPRE BoomerAMG |
| S `run_simple_test_l3_obst.sh` (tag-obst, `--nofine -1`, tag-temp 400), 30 steps | 34/30 | 20.12 / 18 / 23 | 20.40 / 19 / 23 | 3.07 / 3 / 5 | 5.6e-14 | 0.218 / 0.227 | masked MLMG, HYPRE BoomerAMG |
| **C** CSIRO `Case_C064_LS` (24x24x8, `--levels 2 --init-patch --patch 10 19 0 23 --blocking 2 --maxsize 8 --regrid-every 2 --tag-phi --nofine 1`), 60 steps | 62/60 | 10.95 / 10 / 11 | 10.97 / 10 / 11 | **7.85 / 6 / 8 (BiCGStab)** | 8.8e-14 | 0.063 / 0.076 | **unmasked MLMG, BiCGStab bottom** (1 call per cycle; 679 calls, 5,337 iterations over the pred solves) |
| vortex ring 48, 2 levels (`run_vortex_adaptive.sh` flags, 1 rank, 10 steps, 552 L1 boxes) | 11/10 | 11.45 / 10 / 12 | 11.50 / 11 / 12 | 4.33 / 3 / 6 (BiCGStab) | 8.8e-14 | 0.071 / 0.088 | unmasked MLMG, BiCGStab bottom |
| `obst_burner_two_level`, 2 levels, 10 steps | 10/10 | 13.80 / 13 / 14 | 13.80 / 13 / 14 | 3.09 / 3 / 5 | 4.8e-14 | 0.107 / 0.104 | masked MLMG, HYPRE BoomerAMG |
| `open_burner_two_level` (burning, tag-vort), 20 steps | 19/15 | 16.42 / 11 / 33 | 18.87 / 11 / 33 | 3.57 / 3 / 5 | 4.0e-14 | 0.144 / 0.155 | masked MLMG, HYPRE BoomerAMG |
| `obst_offgrid_levels`, 2 levels, 6 steps (tag-obst) | 6/6 | 16.83 / 16 / 17 | 16.83 / 16 / 17 | 3.10 / 3 / 5 | 5.4e-14 | 0.162 / 0.146 | masked MLMG, HYPRE BoomerAMG |
| `obst_offgrid_levels`, 3 levels, 6 steps | 6/6 | 31.17 / 31 / 32 | 31.67 / 31 / 32 | 3.03 / 3 / 4 | 8.3e-14 | 0.380 / 0.433 | masked MLMG, HYPRE BoomerAMG |
| `sealed_two_zone` (thin partition, closed), 3 levels, 9 steps (hierarchy gone after step 3) | 4/3 | 25.00 / 25 / 25 | 25.00 / 25 / 25 | **4.06 / 3 / 6 (BiCGStab)** | 7.1e-14 | 0.298 / 0.277 | masked MLMG, **BiCGStab bottom** (thin faces seal regions; 100 calls, 406 iterations) |

Notes. (a) The corrector and predictor behave the same except where noted; the last-cycle factors are identical to 3 digits. (b) Per solve the number of bottom calls equals the cycles. Total bottom calls over a run is therefore n times the cycles of the table. (c) The `n` of the S rows is small for steps 1-5 because the 1,872-box phase lasts 5 steps. (d) `obst_offgrid_levels` and `sealed_two_zone` hierarchies are the static patches of their test scripts and shrink after the first regrid; numbers are for the composite solves that exist. (e) `pyramid_48_4m` (levels 1, projection OFF) and the other `levels 1` scripts use the FDS level-0 path and carry no MLMG solve; not applicable.

Final/initial: the masked path runs to tol/10 (target 1e-13 relative), the unmasked path stops at tol (target 1e-12; the achieved 6e-13 to 9e-14 in the tables above are the max-norm stop on the finest level).

### 3.2 HYPRE composite path (outer Krylov) and the same case through MLMG
Auto selects HYPRE (assembled matrix, GMRES with BoomerAMG) for a masked composite **without** an OPEN vent face and without thin faces. Case `B` with the OPEN vent removed (sealed box, thick burner OBST, 3 static levels, 15 steps): Auto = **HYPRE, GMRES 20 outer iterations on all 30 solves (min 20, max 20)**, own residual 3e-14, tolerance 1e-12. The same input forced to MLMG (`--pressure-backend mlmg --pressure-tol 1e-12`): masked MLMG **20.60 / 20 / 22 (pred), 20.47 / 20 / 22 (corr) V-cycles**, bottom HYPRE BoomerAMG 3.16 its/call, factor 0.224 geo / 0.227 asym. Without `--pressure-tol 1e-12` the study flag lowers the default to 1e-10 (17.5 cycles) and is not comparable.

### 3.3 Single-level paths (vortex ring 48, 10 steps, level 0, 1 rank; outer iterations from the pressure CSV, with `--pressure-tol 1e-12`)
| Script | Path | Outer iterations pred / corr (mean/min/max) |
|---|---|---|
| `run_vortex_fft.sh` | FFT, direct | 0 |
| `run_vortex_mlmg_bicgstab.sh` | unmasked MLMG-preconditioned BiCGStab, BiCGStab bottom | 5.09 / 5 / 6 and 5.10 / 5 / 6 outer iterations (V-cycles per iteration and bottom iterations are not instrumented on this path: **not measured**) |
| `run_vortex_hypre_pcg.sh` | HYPRE PCG with BoomerAMG | 21 / 21 / 21 and 21 / 21 / 21 |

### 3.4 Identification of the 2.2 versus 20.2 discrepancy (case B)
Same operator, same solver, same tolerance (1e-12); the difference is the hierarchy **after step 5**:
- Doc 14 section 7/8 B (20.2): fine levels kept (`--nofine -1`: the burner vent's no-refine buffer is removed, so the OBST/temperature tags keep 12 L1 + 24 L2 boxes). Every one of the 98 solves is a composite MLMG solve of 20-22 cycles (mean 20.2).
- Role 2's B (2.2): `tm/new/simple_test.fds` (the same OBST, vent and burner, no DEVC/SLCF lines) with the flags `--levels 3 --patch 0 7 0 7 --patch2 0 7 0 3 --steps 49 --regrid-every 5 --tag-obst --tag-temp 400 --blocking 4 --maxsize 8` and **no** `--nofine`. The vent stays on level 0 (D-010 no-refine buffer), the tags find nothing at the step-5 regrid and the hierarchy drops to level 0 ("finest 0 boxes L1 0" in every later regrid log line). Steps 1-5 are composite (22 cycles, 10 solves); steps 6-49 are level-0 FFT solves with 0 iterations. The mean over all 98 pressure rows is 220 / 98 = 2.24 cycles per stage (2.24 predictor and 2.24 corrector over 49 rows each). Both inputs, run with the same binary, give the same hierarchy; the two older binaries tried (the Role 2 tip `d09c02b93c` and an older `tip75` build) also dropped the fine levels, so the input file (with or without DEVC/SLCF lines) is not the cause: the flag is.
So 2.2 is the mean over a mostly-level-0 run; the composite solves of that run took 22 cycles, in line with 20.2.

## 4. Controlled variations (case B geometry, 3 static levels `--patch 0 7 0 7 --patch2 0 7 0 3 --nofine -1 --regrid-every 0`, 15 steps, tol 1e-12; 36 L1 / 144 L2 boxes; cycles pred mean/min/max, corr in brackets where it differs by more than 0.3)
| Variation | Path | V-cycles | factor geo / asym | Bottom its/call |
|---|---|---|---|---|
| **Base** (burner OBST on, 3 levels) | masked, HYPRE | **20.73 / 20 / 22** | 0.228 / 0.232 | 3.06 |
| mask nearly off: one-cell OBST in a far corner, burner on the floor | masked, HYPRE | 20.53 / 19 / 22 | 0.226 / 0.238 | 2.48 |
| OBST removed (burner vent moved to the floor) | **unmasked**, BiCGStab | 23.47 / 22 / 25 (tol stop 6e-13, see below) | 0.301 / 0.323 | 21.35 / 20 / 26 |
| OBST aligned to the coarsest MG cell (0.4 m: no half-blocked coarse faces at any MG level) | masked, HYPRE | 21.20 / 20 / 23 | 0.235 / 0.235 | 3.00 / 3 / 3 |
| bottom = BiCGStab (switch) | masked | 20.73 / 20 / 22 | 0.228 / 0.232 | 23.87 / 16 / 29 |
| bottom = smoother only (switch) | masked | **177.3 / 161 / 195** | 0.844 / 0.870 | - |
| 2 levels (patch only, 7% of the fine grid) | masked | 14.93 / 14 / 16 | 0.127 / 0.129 | 3.08 |
| 2 levels, small (2% / 7% / 25% of the domain refined) | masked | 15.07 / 14.93 / 14.87 | 0.128 / 0.127 / 0.126 (asym 0.131 / 0.129 / 0.126) | 3.09 |
| 2 levels, **level 1 covers the whole domain (no C/F interface)** | masked | **29.80 / 27 / 33** | 0.360 / 0.410 | 3.98 |
| same, one-cell OBST | masked | 29.80 / 27 / 33 | 0.360 / 0.410 | 3.99 |
| same, OBST removed | unmasked | 24.60 / 22 / 27 (corr 24.47) | 0.318 / 0.375 | 24.30 / 24 / 29 |
| 3 levels, both fine levels cover the whole domain, OBST removed (3 steps) | unmasked | 33.00 / 33 / 33 | 0.426 / 0.448 | 24.06 |
| 2 levels, OBST removed | unmasked | 22.60 / 21 / 24 | 0.286 / 0.323 | 20.27 |
| cell aspect 2:1 (z extent halved, dz = 0.05, 36x24x24) | masked | 26.18 / 24 / 29 | 0.311 / 0.320 | 3.30 |
| cell aspect 2:1, one-cell OBST | masked | 26.27 / 24 / 29 | 0.311 / 0.326 | 3.54 |
| cell aspect 2:1, OBST removed | unmasked | 33.13 / 30 / 36 | 0.427 / 0.459 | 27.08 |
| **cell aspect 4:1** (z extent 0.6, dz = 0.025) | masked, and unmasked, and one-cell OBST | **diverges**: first solve hit the 200-cycle cap with a residual of 1e16-1e18 (factor 1.21-1.23 per cycle); run aborts | > 1 | - |
| domain aspect 6:2:1 (72x24x12 cells, isotropic cells, vent moved to the new x face) | masked | 20.89 / 20 / 23 | 0.230 / 0.233 | 3.10 |
| tolerance 1e-12 (base) / 1e-10 / 1e-8 (both stages) | masked | 20.73 / 17.60 / 14.60 (corr 20.67 / 17.60 / 14.47) | 0.228 / 0.227 / 0.227 | 3.06 / 3.07 / 3.09 |
| **S** (burning, tag-vort, 10 steps, 60 L1 + 1,854 L2 boxes in steps 1-5), **OBST removed**, burner on the floor | unmasked, BiCGStab | steps 1-5 (12 solves): 27.58 / 25 / 30; steps 6-10 (11 solves): 22.18 / 21 / 23 (with the OBST: steps 1-5 pooled 32.80 / 29 / 35, steps 6-10 level 0 only) | - / 0.36 | 23.4 (steps 1-5), 16.2 (6-10) |
| S static 3 levels (15 steps), OBST on / OBST removed | masked / unmasked | 20.81 (corr 21.33) / 24.50 (corr 24.13) | 0.233 / 0.313 | 2.95 / 21.56 |
| masked B with the burner OPEN vent removed, forced MLMG | masked | 20.60 | 0.224 / 0.227 | 3.16 |
| same, Auto | HYPRE GMRES | 20 outer iterations | - | - |

The aspect-4 run with the first guard stencil (OBST top 0.2 m with patch2 up to 0.2 m): "NotBuilt: a Solid/Known coarse cell lies in the coarse/fine interpolation stencil" (the obstacle crossed the C/F interface); re-run with the patches extended, it diverged as listed.

Notes on reading the table.
- Unmasked rows stop at `tol` (own residual about 6e-13) and masked rows at `tol/10` (about 5e-14); at the asymptotic factor of 0.32 one decade costs 2.0 cycles, so the unmasked 23.5 corresponds to **about 25.5 cycles at the masked tolerance**. The with/without OBST comparisons below use the raw numbers and that adjustment.
- The unmasked operator (plain Poisson, BiCGStab bottom, no OBST) is **not** faster than the masked one for partial patches: asymptotic 0.32 against 0.23, and 0.32 for 2 and for 3 levels where the masked path goes 0.13 and 0.23. It is less sensitive to the fine-level volume (raw 22.6, 23.5, 24.6 cycles for 2 levels, 3 levels, whole-domain level 1; masked 14.9, 20.7, 29.8). A one-cell OBST in a far corner switches the code to the masked operator and returns to the masked numbers (20.5 against 20.7 for the big OBST). So it is the operator/path (masked ABecLap versus the unmasked MLPoisson with FFT-bottom machinery) that differs, not the presence of solid cells. Why the unmasked path is slower was **not isolated** (candidates: smoother or C/F treatment of the unmasked operator); the bottom solver is excluded (the masked path with a BiCGStab bottom gives the same 20.73 cycles as with the HYPRE bottom).

## 5. What drives the cycle count (ranked)
Cycles per solve to 1e-12 (base 20.7 on the burner geometry).
1. **Number of AMR levels and the fine-level volume (+6 to +15 cycles on the masked path; flatter, +1 per step of the series, on the unmasked path).** AMR levels above 0 get only a mini-cycle (2 sweeps) with no coarsening, so smooth error inside a fine patch is left to the level-0 correction. 2 levels, partial patch: **14.9** (factor 0.13); 3 levels: **20.7** (0.23); 2 levels with level 1 over the whole domain (no C/F interface): **29.8** (0.41); S steps 1-5 (3 levels, 626k fine cells): **31.6-34** (0.42); `obst_offgrid` 2 versus 3 levels: **16.8 versus 31.2**. The size of a **partial** patch (2%, 7% or 25% of the domain) changes nothing (15.07, 14.93, 14.87); the C/F interface itself is not the driver (the full-coverage case has none and is the worst). S later phase (small patches): 13.3.
2. **Cell aspect ratio (+5.5 cycles at 2:1, divergence at 4:1).** 2:1: **26.2** (0.32 asymptotic); 4:1: the point GSRB smoother does not converge (factor above 1; also without OBST and with a tiny OBST). Isotropic cells in a long domain (6:2:1): 20.9, no effect. Domains of aspect up to about 2 in the cells are workable; 4 is not without line relaxation or semicoarsening (not available in the MLMG used).
3. **Tolerance (1.55 cycles per decade).** 1e-12: 20.7; 1e-10: 17.6; 1e-8: 14.6. Matches log(10)/log(1/0.228) = 1.57. The factor itself does not depend on the tolerance.
4. **Bottom solve accuracy: no effect** (HYPRE BoomerAMG 3 its/call against BiCGStab 24 its/call, both reaching 1e-4: same 20.73 cycles). A smoother-only bottom costs **177 cycles** (0.84), so a bottom that reduces the coarsest residual (any Krylov or AMG) is required.
5. **OBST mask and coarsened beta: none measurable.** Big OBST 20.73, one-cell OBST 20.53, OBST aligned to the coarsest cells (no half faces) 21.20 (+0.5, not an improvement), burner OBST in the census cases 20.2-22. Removing the OBST switches to the unmasked path, and the sign of the effect depends on the hierarchy: B and S with small static patches need **more** cycles (B 23.5 raw, about 25.5 at the same tolerance, against 20.7; S static 24.5 / 24.1 raw against 20.8 / 21.3), while S with its 1,854-box start needs fewer (steps 1-5: 27.6 raw, about 29.6 at the same tolerance, against 32.8). The OBST is therefore not what sets the count; the hierarchy and the operator path are.
6. Domain/box shape with isotropic cells: none (20.9).

## 6. Against the expectations of doc 13
Doc 13 expects 0.2-0.35 per V(2,2) cycle (15-22 cycles to 1e-10) for the masked operator, and a factor above 0.5 to mean the snapped-mask or alpha-leak rulings are not in force. Measured: 0.23 asymptotic on B and S-static (**20.7 cycles to 1e-12, 17.6 to 1e-10**: inside the expected 15-22), 0.10-0.13 for 2-level patches and S later steps, 0.41-0.43 for the large hierarchies (S steps 1-5, `obst_offgrid` 3 levels, whole-domain level 1), above the expected band but caused by the AMR level structure (the same numbers appear with a one-cell OBST and with no OBST), and a divergence at cell aspect 4. The coarsened mask/beta does not degrade the factor: the OBST-aligned case (no half-blocked coarse faces) is no better than the half-faces case, and the one-cell OBST is equal to the large OBST. No case exceeded 0.5 other than the aspect-4 divergence.

## 7. Not measured / caveats
- A true single-level masked MLMG (no AMR levels at all, same resolution): the `levels 1` path of the driver is the FDS level-0 solver (FFT with the OBST correction, HYPRE or the single-level MLMG-GMRES), not the composite operator, so the 1-level point of the "1 versus 2 versus 3 levels" series is not available; the nearest are 2 levels with a 2% patch (15.1 cycles) and the unmasked rows.
- Why the unmasked composite operator converges slower than the masked one (section 4).
- V-cycles per BiCGSTAB iteration and bottom iterations on the single-level MLMG-BiCGSTAB script (section 3.3).
- Aspect between 2:1 and 4:1 (only 2:1 and 4:1 run); aspect with the HYPRE composite path.
- S with the OBST removed beyond 10 steps (S steps 6-10 without OBST are composite solves, with the OBST they are level-0 solves, so no like-for-like pair exists there); the pyramid and `levels 1` cases.
- Case counts are single runs with 1-6 other jobs on the box; cycle counts are deterministic. An earlier S attempt was killed by the box (memory with six parallel runs) and re-run alone.
- 1 rank only: no distribution effect on the counts was checked.
