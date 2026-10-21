# 14 — Composite pressure solve profile (3-level simple_test)

Status: **MEASURED** on the owner-provided NVIDIA test machine (CPU run: 1 rank, 8 OpenMP threads, GCC 13, Open MPI 5, HYPRE 2.32, AMReX 26.09 built with TinyProfiler), tip `fc21986dda` plus throwaway timers in a scratch copy (no product change, nothing pushed). Estimates are labelled. Runs: 60 steps each, one run per row (no repeats; run-to-run spread about 5%, e.g. three equal-config baselines gave 275, 290, 297 s loop time).

## 1. Case and command
`Verification/Adaptive_Mesh_Refinement/run_simple_test_l3.sh`: simple_test.fds (36x24x24 level 0, burner OBST plus partial OPEN vent), `--two-level-run --levels 3 --regrid-every 5 --tag-vort auto`, driver blocking 4. The CEO case (about 1,900 fine boxes) is the start of this run: steps 1-5 carry **64 L1 + 1,872 L2 boxes** (165k + 461k cells); from the step-5 regrid on the hierarchy drops to 3-100 boxes per fine level (about 25 L1 / 77-100 L2 near step 60, 12k-36k cells).

## 2. A. Exact configuration (as run, default options, no study flags)
| Item | Value |
|---|---|
| Path | composite pressure, `PressureBackend` Auto, `trigger SolveRoutine`; 2 solves per step (predictor, corrector) |
| Selector | masked composite and a mixed face (OPEN vent: `face_neumann`) so Auto is **MLMG** (Auto to HYPRE applies only without face_neumann and thin faces; explicit HYPRE would be NotBuilt with thin faces) |
| Operator | MaskedABecLap, beta 0 on solid faces, area-averaged coarse beta, 3 MG levels on level 0 (max_coarsening_level 2, coarsest 9x6x6), GSRB smoother, plain V-cycles (Krylov Auto = none) |
| Bottom solver | HYPRE BoomerAMG (`bicg_bottom` was 0 in all 122 operator builds: the BiCGStab bottom is chosen only if the bottom has several singular regions AND thin faces; it did not trigger. In the TinyProfiler counts BiCGStab ran in 11% of bottom solves: AMReX fallback when the hypre bottom did not reduce the residual, not the configured bottom) |
| Gauge | D-067 per connected fluid region, one component here |
| Tolerance | `tol_rel` 1e-12 default, V-cycles run to tol/10: stops at **own relative residual about 1e-13** (bnorm); max_iter 200 |
| Residual check | off in default runs (study flags turn it on); D-087 anorm only enters when the check is on |
| Ranks/threads | 1 rank, 8 threads; level 0 is one box |

## 3. B. Time breakdown
Wall per composite solve, split by timers around each phase (`label` = component labelling, `setup` = MaskedMlmg/operator build, `mlmg` = V-cycles, `other` = fold, ghost fills, gauge, average-down, copies):

| Regime (baseline, default build) | Step | Pressure share of step | Per solve | label | setup | MLMG V-cycles | other |
|---|---|---|---|---|---|---|---|
| Steps 1-5, 1,872 L2 boxes (10 solves) | about 52 s | **about 75%** | 19.6 s | **16.5 s (84%)** | 0.43 s (2%) | 1.58 s (8%), 32.7 cycles | 1.1 s (6%) |
| Steps 6-60 (110 solves) | 0.78 s | **about 45%** (0.35 s of 0.78 s) | 0.17 s | 0.072 s (42%) | 0.012 s (7%) | 0.065 s (38%), 13.5 cycles | 0.024 s (14%) |
| Whole 60 steps | 297 s loop | **72.8%** (216 s) | 1.79 s | 80% | 2.6% | 10.7% | 6.2% |

(150-step run before the timers: 390 s loop, pressure 68%, in the CEO's 62-71% range. The 4,572-step owner run with about 260 final boxes: 50%.)

Where the labelling time goes (1,872-box solve, 21-27 rounds): fine-side C/F half of each round 11.4-14.6 s of 12.5-15.9 s; coarse-side 0.29 s; Q2 Known-on-C/F test 0.57 s; same-level sweeps 0.10 s. The fine-side loop calls `BoxArray::contains` on 1,872 boxes for every neighbour of every fine cell, every round. TinyProfiler (whole-run totals) shows the MLMG internals as small: smoothing 1.7%, applyBC 1.6%, `BoxList::complementIn` 1.6% (571k calls from BoxArray searches and ParallelCopy metadata), FillBoundary 1.3%, ParallelCopy 1.1%, `CPC::define` 0.4%, HYPRE bottom setup 0.45% and solve 0.08%; norm/dot reductions below 0.1% (1 rank, so no MPI cost measured).

**Defect found (measured): the workspace is rebuilt on every solve.** `ws_rebuilt = 1` on all 120 solves, no regrid needed. Cause: `impl_matches` compares `W.face_neumann` with the problem's, but `build_impl` returns early for masked problems before storing it (`W.face_neumann = p.face_neumann` sits after the masked return), so any masked composite with an OPEN-vent mixed face never matches. Each solve then rebuilds the workspace, the MaskedMlmg operator and the hypre bottom setup. Throwaway test (store it before the return): setup 0.43 -> 0.046 s per solve at 1,872 boxes and 0.012 -> 0.001 s later; `ws_rebuilt` 0 everywhere; loop 297 -> 290 s (-2%), later steps 0.78 -> 0.72 s (-7%). A fix is a one-line move; no product change was made here.

## 4. C. Convergence
- Zero guess: initial residual = rhs norm (ratio 1.0). V-cycle factor: 0.40 per cycle in the 1,872-box phase (33 cycles to 8e-14), **0.05-0.13** (median 0.12) per cycle later (11-15 cycles). Cycle 1 already cuts 0.55 (first phase) or 0.02-0.16 (later), so the per-cycle cost is spread evenly; the first-phase factor 0.40 is above the 0.35 expectation of doc 13.
- Warm start (`pressure.warm_start=1`, previous-step pressure): initial residual only drops to **0.66** of the rhs norm (median, later steps); iterations 13.5 -> 13.8 (no gain) at 1e-12, 9.3 -> 9.5 at 1e-8. The H field changes too much per step (flow, burner, regrids) for it to help.
- What the divergence needs (measured max |div u - D| over the run, tol_rel as requested; own residual is tol/10):

| tol_rel | cycles per solve (first phase / later) | worst fine-level div | project gate (1e-9) |
|---|---|---|---|
| 1e-12 (default) | 32.7 / 13.5 | 7.3e-11 | PASS |
| 1e-8 | 22.2 / 9.3 | 5.3e-7 | FAIL |
| 1e-6 | 16.9 / 7.2 | 4.0e-5 | FAIL |
| 1e-4 | 11.6 / 5.3 | 3.3e-3 | FAIL |

The divergence error is about 500-1000 times the requested tolerance; the default 1e-12 gives a margin of about 14 against the 1e-9 gate. FDS's own criterion (`VELOCITY_TOLERANCE` = 0.5 x cell size, about 0.05 m/s, not evaluated in refined runs, see `driver/notes/composite-pressure-iteration.md`) is much looser, so the need is set by the project's gate, not by FDS. Tightness needed for the gate: about 1e-11 requested (estimate by interpolation, not run), saving about 1-2 cycles of 13.

## 5. D. Cut options
Pressure per step is the measured quantity; step savings assume the baseline split above (later regime 0.35 of 0.78 s; first phase 39 of 52 s).

| Option | Evidence | Pressure saving | Step saving |
|---|---|---|---|
| Interior-skip in labelling (`vbx.contains(n) || fba.contains(n)`, patch from Role 3, not applied here) | **measured by Role 3**: t_loop 84.2 -> 58.8 s, bit-identical CSV, 90 steps, 3 levels | about 30% of loop | -30% |
| Cache labels (component ids) between regrids; keyed on BoxArray plus class hash (labelling today runs on every solve, 2 per step, regrid every 5) | **estimate**: removes about 90% of label time (1 of 10 solves per regrid period needs it) | 1,872-box phase: -80%; later: -38% | first phase -60%, later -17% |
| Fix workspace reuse (face_neumann stored for masked) | **measured** (throwaway) | setup 100% (2% / 7% of solve) | -2% total, -7% later steps |
| Looser tolerance | **measured**, but fails the 1e-9 gate at 1e-8 and above | 1e-8: MLMG -45% (cycles -31%), 1e-6: -57% | step -4% at 1e-8, -15% at 1e-6 (MLMG is only about 10% of the solve in the first phase, 38% later) |
| Warm start | **measured** | none (0 to +2% cycles) | none |
| Fewer cycles / V(1,1), better bottom (BiCGStab vs BoomerAMG), max_coarsening_level, blocking or box size (fewer larger boxes), `--krylov bicgstab` | **not run** (chain stopped for time); **estimate**: MLMG is 8-38% of the solve, so even halving cycles gives -4% to -19% of the solve; bottom is 0.1-1.5% of MLMG time in TinyProfiler; fewer boxes help labelling roughly linearly in box count (cost is contains-per-cell, about quadratic in box count x rounds) | est. -5 to -20% (cycles); label est. large | est. -2 to -9% |

Order of value: (1) interior-skip, (2) label cache, (3) workspace reuse, then solver tuning. After (1)+(2)+(3) the expected solve in the later regime is about 0.065 s MLMG + 0.02 s other per solve, i.e. pressure about 25% of the step (estimate), where tolerance and cycle count then matter. Items (1)-(3) do not change the solution; tolerance does.

## 6. Not measured / caveats
- 1 rank only (no MPI reduction, no distribution effects); CPU only (GPU laptop not used for this run).
- 60-step windows; the 1,872-box phase lasts 5 steps and then the hierarchy stays small, so full-run fractions depend on window length (60 steps 73%, 150 steps 68%, 4,572 steps 50%).
- The runs with `pressure.warm_start`, tolerances, and the 1e-8 / 1e-6 / 1e-4 rows were launched with the study flags, which add a true-residual check per solve; the timers report solver and phase times only. The tolerance rows 1e-8 and looser end with the run's own divergence-gate failure (exit code 6) after all steps completed.
- Bottom solver, V(1,1), blocking and box-size experiments were queued but not run.

## 7. After cache: ranked next cuts

Status: **MEASURED** on the owner-provided NVIDIA test machine (CPU run, 1 rank, 8 OpenMP threads, GCC 13, Open MPI 5, HYPRE 2.32, AMReX built with TinyProfiler), tip `8c16d56d50` (label cache, workspace reuse, interior-skip) plus throwaway timers and environment switches in a scratch copy (no product change, nothing pushed). One run per row unless stated; run-to-run spread about 5-10% (three equal burner baselines: 12.81, 13.04, 13.11 s pressure). Estimates are labelled.

Cases: (S) `run_simple_test_l3.sh` (3 levels, regrid every 5, tag-vort auto); (B) near-burner 3-level case (`burner_obst_3level.fds`, patch plus patch2 next to the burner, tag-obst/tag-temp, regrid every 5, blocking 4, maxsize 8; 54 L1 / 180 L2 boxes at start, 12 / 24 after step 5). Case B is a reconstruction of the Role 2 case (its exact flags were not recorded; initial box counts differ slightly from Role 2's 61 / 144). Case B ends its gate with a known mass/species-balance FAIL at this tip (div 3.6e-13, unrelated to pressure); the abort was bypassed with a throwaway switch only to get the TinyProfiler print.

### 7.1 Pressure share of t_loop
| Run | t_loop | pressure | share | cycles/solve |
|---|---|---|---|---|
| S 49 steps | 96.0 s | 32.7 s | **34%** | 15.3 |
| S 49 steps, steps 1-5 only (1,872 L2 boxes) | 79.0 s | 25.8 s | 33% | 32.7 |
| S steps 6-49 (by difference) | 17.0 s | 6.98 s | **41%** | 13.3 |
| S 150 steps | 176.2 s | 57.9 s | **33%** | 14.8 (later 14.2) |
| B 49 steps | 25.9 s | 11.9 s | **46%** | 20.2 |
| B to T_END (77 steps; 150 not reachable) | 38.9 s | 18.9 s | **49%** | 20.3 |

So the "about 45%" is case B and the small-hierarchy phase of S; in S the 1,872-box start is dominated by non-pressure work (about 15 s per step against 4 s of pressure).

### 7.2 (1) Labels and stamp
- Labels run only on the first solve after a regrid: 9 of 98 solves in S, 10 of 102 in B (hit rate 91% / 90%). The stamp changes only at regrids; it is cheap (0.03 s over 98 solves in S, 0.004 s in B; 2.5 ms per solve in the 1,872-box phase), so skipping the stamp is not worth it.
- Remaining labelling: S first solve 4.5 s (1,872 boxes, 17% of that phase's pressure), later regrids 0.012-0.038 s each (about 3% of later pressure, 1.2% of later t_loop); B first 0.46 s then 0.06 s per regrid, 1.03 s total = **9% of pressure, 4.0% of t_loop**.
- Finding: in B, 8 of 9 regrids report `changed 0` (identical hierarchy) yet the workspace is dropped: `bind_level`/`unbind_level` bump `layout_gen`, `ensure_comp` rebuilds the runner, `crhs/cphi` and the pressure workspace, so the labels, the masked MLMG (0.012 s) and the BoomerAMG setup are redone. Skipping that rebuild when the layout is unchanged is bit-identical (same labels). Estimated saving in B about 0.6-0.7 s (2.5% of t_loop); none in S (every regrid changes).

### 7.3 (2) V-cycles and cost per cycle
- Factor per cycle (B and S later): cycle 1 about 0.14, then 0.23-0.24 asymptotic (B), S later 0.05-0.13 early, 13-15 cycles. Solve 1 (predictor) and solve 2 (corrector) are the same: 20.2 / 20.2 cycles in B, 14.9 / 15.7 in S.
- MLMG share: B MLMG solve 8.70 s of 11.9 s pressure (73%); S later 5.65 of 6.98 s (81%). TinyProfiler inclusive times, B (S later in brackets), seconds over the run; the groups overlap slightly:

| Phase | B | % pressure (B) | S later | % pressure (S later) |
|---|---|---|---|---|
| smoothing incl. its BC fills (Fsmooth kernel 1.84 / 1.17) | 3.59 | 30% | 2.21 | 32% |
| residuals with C/F corrections (`computeRes*`) | 2.59 | 22% | 1.51 | 22% |
| BoomerAMG bottom | 1.20 | 10% | 1.12 | 16% |
| other MLMG (interp, average-down, norms, prepare) | about 1.3 | 11% | about 0.8 | 11% |
- C/F boundary fills inside these: `applyBC` 2.52 s incl. (FillBoundary 0.91, ParallelCopy 1.02, setVal 0.70); `average_down` 0.48; interp 0.25; reflux 0.73. Reductions (Dot, norminf) are below 0.1 s (1 rank, no MPI cost measured).

### 7.4 (3) Bottom solver
> **Erratum 2026-10-21 (bottom-solver rows).** The `pb_mlmg_bottom.solver=smoother|bicgstab` runs of this section (prof14b `b_botsmooth`, `b_botbicg`, `d_botsm*`, and the earlier `x9_botbicg`) did **not** change the bottom solver. The driver forwards only `hypre.*` and `pressure.*` command-line tokens into ParmParse (`Source/driver/main.cpp` `pressure_flag`, since `25640298df`, 2026-10-08) and calls `amrex::Initialize(argc, argv, false, ...)` (no argv parsing, since `004b593534`), so the `pb_mlmg_bottom.*` keys read in `CompositeSolve.cpp` `build_masked_mlmg` never arrived; the scratch copy used for these runs had the same filter. Evidence: `simple_test_devc.csv` and `simple_test_multi_level.csv` of `b_botsmooth` and `b_botbicg` are byte-identical to the baseline `b0`. Every run used the default HYPRE BoomerAMG bottom. The "bottom smoother / BiCGStab: no change" results are therefore **invalid (not measured)**; the BoomerAMG profile numbers and the persistent-setup estimate are unaffected. Forwarding of `pb_mlmg_bottom.*` was fixed in src (see doc 16 erratum); the comparison needs re-running.

- BoomerAMG bottom is 10% of pressure in B (setup 0.92 s, solve 0.17 s: **the setup is rebuilt at every V-cycle, 5x the solve**) and 16% in S later (setup 0.54, solve 0.15, plus 0.41 s of BiCGStab fallback in 133 solves). Total 4.6% of t_loop in B, 6.6% in S later (est. from the TinyProfiler difference).
- The coarsest grid is 9x6x6 (324 cells). `max_coarsening_level` can only be lowered at runtime (it is already 2, the limit of the 36x24x24 domain): level 1 gave pressure 19.7 s (+65%), level 0 79.8 s (+570%). No larger coarsening is possible.
- ~~Bottom `pb_mlmg_bottom.solver=smoother` or `bicgstab`: cycles unchanged (20.2) and time unchanged within noise (3 paired repeats: baseline 9.37 / 9.51 / 9.59 s MLMG, smoother 9.43 / 9.51 / 9.54 s). No gain from the runtime option.~~ **Invalid (erratum 2026-10-21):** the key was not forwarded, all three variants ran the default BoomerAMG bottom (outputs byte-identical), so the numbers are a repeat of the baseline, not a bottom comparison. Not measured. A persistent BoomerAMG setup would be a code change: est. up to -0.9 s in B (3.6% of t_loop), -0.5 s in S later (3%).

### 7.5 (4) Smoothing
Runtime switch added in the scratch copy only (AMReX has none). B, MLMG seconds / cycles: V(1,1) 8.89 / 24.4, V(2,2) 8.5-9.6 / 20.2, V(3,3) 10.49 / 18.2, V(4,4) 11.11 / 16.4. S 49 steps: V(3,3) 13.6 cycles, loop unchanged. V(2,2) is the optimum; the post-smooth count had no effect on the result (V(2,1) identical to V(2,2): this MLMG path applies the pre count on the way up too). No gain.

### 7.6 (5) Box size
B pressure seconds (t_loop): maxsize 8 (default of this case) 11.8-13.1 (25.5-27.8), maxsize 16 11.15 (17.8), maxsize 32 10.28 (15.3); cycles unchanged. Larger boxes cut pressure only 5-13%, but t_loop by 30-40% (the non-pressure fine-box work), so the share rises to 63-67%: a box-size change helps the step, not the share. S: maxsize 8 is worse (pressure 36.5 s, loop 113 s), maxsize 16 equals the default, blocking 2 worse (39.2 s, 139 s). Blocking 8 is invalid for 36 cells. Box layout changes reduction order only (est.; not bit-compared here).

### 7.7 (6) Second solve, warm start, tolerances
- Initial residual / rhs: solve 1 and 2 both 1.00 from a zero guess. Warm start of the corrector from the predictor pressure: 0.67 of rhs, cycles 20.2 -> 19.4; predictor from the previous corrector: 0.62, 19.1; both: 0.29, **18.3 (-9%)**, MLMG time 8.63 s vs 8.5-9.6 (not detectable); S: 15.3 -> 15.6 (no gain). Reason: asymptotic factor 0.23 dominates, a 3x smaller start saves about 2 cycles of 20.
- The rhs differs for every solve (initial rhs norm 0.833 / 0.809 for predictor / corrector), so there are no repeated solves for an unchanged rhs.
- Looser tolerance on the predictor only (idea: the corrector solve fixes the final divergence): not measured. The throwaway switch did not change the cycle count (20.2 in all runs), so the result is void. Est. -30% cycles on solve 1 from the tolerance table in section 4 (-15% of MLMG, about 5% of t_loop in B); needs a justified-tolerance argument for the predicted state before use (US feeds the corrector fluxes).
- Krylov acceleration (study mode, true-residual check on): `--krylov bicgstab` 5.9 iterations, MLMG 8.11 s; `gmres` 10.2, 8.81 s; plain V-cycles 17.2, 8.13 s. No gain (each Krylov iteration costs about 3 V-cycles).

### 7.8 (7) Other items per solve, B (S later similar)
Per solve 0.1147 s: MLMG 0.084 (73%), label 0.010 (9%), `face_gradient_composite` 0.0089 (7.8%, 3.5% of t_loop), pre-solve (zero solid faces, fold, rhs) 0.0059 (5.1%, 2.3%), H copy-in/out and ghost fills 0.0027 (2.4%, 1.1%), rhs/mean removal/masked MLMG setup 0.0034 (3%), post 0.0004. The workspace and masked MLMG are reused (setup 0.001 s except at regrids).

### 7.9 Ranked cut table (percent of t_loop; B is the 46% case)
| # | Option | Saving in t_loop | Basis | Bit-identity / gates | Simplicity |
|---|---|---|---|---|---|
| 1 | Persistent BoomerAMG bottom setup (reuse across V-cycles and solves; or a cheaper bottom) | 3.6% (B), 3% (S later) | TinyProfiler setup time (upper bound) | identical if matrix unchanged | medium (AMReX hypre object lifetime) |
| 2 | Skip workspace/label/MLMG rebuild on a regrid that leaves the layout unchanged | 2.5% (B), 0 (S) | measured cost per regrid | identical | simple |
| 3 | Predictor-only looser tolerance | about 5% (B) | estimate | not identical; needs tolerance argument | simple, risk to gate argument |
| 4 | Warm start of both solves | about 2-3% (B), 0 (S) | cycles -9% measured, time estimated | not identical (roundoff); tolerance unchanged | simple, keeps state across regrids |
| 5 | `face_gradient_composite` once / only fine-ghost columns | up to 3.5% | estimate (cost measured) | identical if same values | medium |
| 6 | Pre-solve and H copy/ghost trimming | about 1-2% | estimate | identical | medium |
| 7 | Larger boxes (maxsize 16-32) | 5-13% of pressure, 30-40% of t_loop but the share rises | measured | layout change, not bit-compared | trivial, raises share |
| - | V(1,1)/(3,3), Krylov, lower coarsening, warm start in S (bottom smoother: **not measured**, erratum 2026-10-21) | 0 or negative | measured | - | - |

Predicted outcome: items 1-6 sum to about 19% of t_loop in B at best (estimate), taking pressure from 46% to roughly 34% share ((11.9 - 4.9)/(25.9 - 4.9)); 30% is not reached by solver tuning alone. A further 3 s would need about 25% fewer cycles at the same tolerance (a better C/F treatment or smoother; the 0.23 asymptotic factor, against 0.14 in cycle 1, is the target) or a cheaper residual/C/F fill path (items 3-6 of 7.3 are 52% of pressure). In S the post-regrid share is 41% and the same list gives about 33% (estimate).

## 8. Deep dive: ranked optimization table

Status: **MEASURED** unless marked *estimate*. Owner-provided NVIDIA test machine, CPU run, 1 rank, GCC 13, Open MPI 5, HYPRE 2.32, AMReX 26.09 with TinyProfiler. Tip `75308b9640` (persistent BoomerAMG setup, operator reuse, ground-layer corrector fix) plus throwaway environment switches in a private scratch copy of the source and of AMReX (nothing pushed, shared AMReX tree untouched, no product change). Case B (near-burner 3-level, 49 steps) is the main case; case S (`run_simple_test_l3.sh`, 30 steps) is used for tolerance and gate checks. Run-to-run spread 3-7%; repeats are listed where taken. Study flags (`--pressure-csv`) add a true-residual check; the study default tolerance is 1e-10, so every study run here passes `--pressure-tol 1e-12` to match the product default (an earlier set at 1e-10 was discarded).

### 8.1 Baseline at this tip (B, 8 threads unpinned)
t_loop 24.4-26.3 s, pressure 10.4-11.5 s (about 43%), MLMG solve 8.0-8.5 s, **20.2 cycles/solve** (S: 11.4). Divergence at the gate: B 3.6e-13, S 3.7e-11 (gate 1e-9). TinyProfiler (B, 27.6 s): MLMG 8.09 s; smooth 3.77 s inclusive (Fsmooth kernel 1.93), 91.7k `applyBC` calls 2.68 s (29 us each), 66.9k Fsmooth calls (29 us each), FillBoundary 0.99 s, ParallelCopy_nowait 1.07 s, setVal 0.73 s, `computeResWithCrseSolFineCor` 1.60 s, reflux 0.80 s, updateSolBC 0.49 s, updateCorBC 0.42 s, BoxList::complementIn 0.43 s (36k calls), BoomerAMG bottom only 0.22 s. Per-call costs of 10-30 us indicate the solve is overhead/latency bound, not bandwidth bound.

### 8.2 Measured: threading (largest finding)
Pressure and t_loop seconds, B, non-study runs:
| Setting | t_loop | pressure | repeats |
|---|---|---|---|
| 8 threads unpinned (current) | 24.4 / 25.9 / 26.3 | 10.4 / 11.4 / 11.5 | 3 |
| 8 threads pinned to the 8 performance cores | 22.7-23.1 | 9.9 | 3 |
| 4 threads pinned | 18.3 / 19.0 | 7.9 / 8.4 | 2 |
| 2 threads pinned | 19.1 | 8.6 | 2 |
| 1 thread pinned | 18.55 | 8.87 | 3 (spread 0.2%) |
| 1 thread unpinned | 18.6 / 18.7 | 8.9 / 9.0 | 2 |
Fewer threads is faster for this small hierarchy (about 20k cells on level 0, 12-100 boxes of 8^3 on the fine levels): pressure -20 to -25%, t_loop -25 to -30%. The pressure *share* rises to 44-48% because the rest of the step gains more. This is a deployment setting (threads/ranks and pinning), bit-identical in the solve (same operations; reduction order untouched at 1 rank), and it is the largest single seconds-saving found. It does not reduce the share toward 30%.

### 8.3 Measured: cycle count and factor (B unless noted, 1e-12)
- Asymptotic factor is set by the fine AMR levels, which get only a single-level mini cycle (AMReX does not coarsen AMR levels above 0; 2 sweeps down-leg, 2 up-leg). Sweeps at fine levels only (down/up): 2/2 20.2 cycles, solve 8.5 s (baseline); 1/1 25.6, 9.2 s; 2/0 27.8; 0/2 28.0; 4/2 18.2, 8.8 s; 2/4 17.5, 8.5 s; 4/4 16.2, 8.8 s; 8/8 13.2, 10.3 s. Cycles fall (factor 0.23 to 0.15 to 0.10) but time per cycle rises in proportion: **flat optimum, no gain**. S 4/4: 10.4 vs 11.4 cycles, solve 18.7 vs 15.4 s (worse).
- FMG initial guess: `setMaxFmgIter(1)` identical 20.2 cycles; 99: 20.5 cycles, +20% time. No gain. W/F-cycle: no runtime option exists in this AMReX (MLMG has V and FMG only); L0 smoothing is not the limiter (V(3,3) gave -10% cycles at +20% time), so a W-cycle would not help (*estimate*).
- C/F boundary order: `setMaxOrder(3)` changes the discrete C/F stencil: cycles unchanged (17.1 vs 17.2) and the divergence gate fails (1.9) because the face gradient used for the velocity correction no longer matches the operator. Rejected.
- Warm start (earlier measurement): 20.2 to 18.3 cycles, no detectable time gain.
- Krylov: not repeated (section 7.7: 3 V-cycles per Krylov iteration, no gain).

### 8.4 Measured: tolerance and gate (divergence is about residual in max norm, so it scales with the rhs max norm of the case)
| tol_rel (requested; own residual is /10) | B cycles | B div | S cycles | S div | gate 1e-9 |
|---|---|---|---|---|---|
| 1e-12 (default) | 20.2 | 3.6e-13 | 11.4 | 3.7e-11 | pass |
| 3e-12 | - | - | 10.9 | 3.7e-11 | pass |
| 1e-11 | 18.7 | 5.5e-12 | 10.6 | 5.9e-11 | pass |
| 1e-10 | 17.2 | 3.0e-11 | 9.5 | **1.5e-9** | S FAIL |
| 1e-9 | 15.5 | 5.0e-10 | - | - | B pass (margin 2) |
A fixed relative tolerance is too tight for B (margin 2800) and just right for S (margin 27). A stopping rule tied to the gate (absolute max-norm residual of about 1e-10 with margin 10) would give B about -15% cycles and S about -7% (1e-11 row). Outputs for 1e-11 (S) differ from 1e-12 by 1e-15 in the multi-level diagnostics: justified-tolerance material, not bit-identical.
**Predictor-only looser tolerance** (corrector stays 1e-12; gate div unchanged at 3.6e-13 B, 3.7e-11 S, because the corrector decides the final divergence): predictor cycles B 20.2 to 14.2 (1e-8), 11.2 (1e-6), 8.1 (1e-4); MLMG solve 8.5 to 7.2, 6.6, 5.9 s; t_loop -4%, -5%, -8%. S: 1e-8 predictor 7.0 cycles, solve 15.4 to 13.0 s, t_loop -3%; 1e-6: 5.4 cycles, 11.7 s, t_loop -5%. Not bit-identical: S output vs baseline differs by 1.6e-9 (1e-8) and 3.3e-7 (1e-6) in the multi-level diagnostics, with the HRR column showing 0.3% (1e-8) and larger (1e-6) relative differences at small values. 1e-8 is the defensible setting; it needs a justified-tolerance argument for the predicted state (it feeds the corrector fluxes).

### 8.5 Algorithmic alternatives (estimates, with risk)
- Composite L0+L1 then level solve of L2 (Dirichlet from L1), or FAC with sync correction: removes the L2 coupling error only if a sync solve follows; without it the C/F flux mismatch leaves coarse-cell divergence at the interface and breaks the divergence gate and the leaf-only integral gate (D-067/D-091). With a sync solve the cost is about the same as the composite solve. Not recommended; no measurement taken.
- Krylov preconditioning: each iteration costs about 3 V-cycles for a gain below 3x (section 7.7).
- Coarsened fine AMR levels (agglomerated 2-level multigrid on each fine level, upstream feature work): the real lever on the 0.23 factor, since fine-level smoothing cannot reduce patch-wide smooth error. *Estimate* 30-40% fewer cycles at about equal cost per cycle; high effort, unmeasured.

### 8.6 Cost per cycle: avoidable work in AMReX (read from the shared tree, read-only; gains are estimates from the TinyProfiler counts above)
1. `MLCellLinOp.H` `smooth` (about line 1551): the first `applyBC` of a zero-guess smooth is redundant (`MLMG.H` 1806, 1850, 2019 pass `skip_fillboundary=true` only after `setVal(cor,0)`, ghosts are zero and homogeneous BC of zero is zero). 6 of 44 `applyBC` per cycle: **about 1.3% t_loop**, bit-identical. Patch drafted, not measured.
2. `MLLinOp.H` `makeCoarseAmr` (about line 1903): builds a new coarsened BoxArray on every call (2 per cycle), so the ParallelCopy metadata cache misses (CPC::define 3,582 calls) and BoxList::complementIn runs (0.43 s): **up to 1.5% t_loop**, bit-identical. Cache the BoxArray per level in the operator. Patch drafted, not measured.
3. Fuse `setVal(cor,0)` (`MLMG.H` 1804, 1848) into the first sweep: about 1.3% (Role 2: 0.50 s), bit-identical in principle; larger kernel change.
4. updateCorBC/updateSolBC, per-cycle residual test, reflux ghost fill: each follows changed data; no safe skip found.
Sum of upstream-patch items 1-3: about 3-4% t_loop.

### 8.7 GPU path
Measured with the existing single-level study harness (mms, Dirichlet, MLMG with BoomerAMG bottom, managed memory, 1 rank, 100 re-solves, 3 repeats; GPU is the 8 GB RTX 4070 laptop part, CPU is 1 thread on one performance core):
| Cells | Boxes | CPU resolve | GPU resolve | Speedup |
|---|---|---|---|---|
| 32^3 (33k) | 1 | 22.1 ms (24 cycles) | 22.5 ms | 1.0x |
| 32^3 | 8 | 26.7 ms | 25.8 ms | 1.0x |
| 64^3 (262k) | 1 | 205 ms (30 cycles) | 50 ms | 4.1x |
The CSIRO-style crude hierarchy (level 0 about 20k cells, 3 levels, tens of thousands of cells in total, under 10 MB of fields) sits at or below the 32^3 point: **about 1x, no gain** (launch/latency bound, about 0.9 ms per cycle either way). A production-size level 0 (about 2.3M cells, 3 levels, 50-100 MB of fields, fits in 8 GB) is *estimated* at 3-4x on the MLMG solve (doc 05: 2x vs the best 8-rank CPU at 1M cells, 4x vs one thread at 262k). Blockers: the HYPRE bottom (device BoomerAMG latency on 324 unknowns; alternatively smoother or BiCGStab bottom, not yet measured on CPU: the earlier "neutral" result was a no-op run, erratum 2026-10-21), `LoopOnCpu` host loops in the masked operator (`prepareForSolve`, region labelling) and in the per-solve rhs/phi preparation and component labelling (managed-memory page migration each solve), and all FDS Fortran host loops (the rest of the step stays on the CPU). A 3-level GPU composite measurement was not possible with the existing harness (single-level only); the table above is a calibration, not a composite result. Amdahl: MLMG is about 30% of t_loop, so even an infinite GPU speedup caps the t_loop gain at about 30% and leaves the share near 15-20% only if the rest of the step is not also moved.

### 8.8 Ranked table (gain in % of t_loop; B is the main case, S in brackets where measured)
| # | Item | Expected gain | Risk | Effort | Bit-identical / tolerance / gates |
|---|---|---|---|---|---|
| 1 | 1-4 pinned OpenMP threads (or ranks) instead of 8 unpinned | **25-30% of t_loop** (seconds; share rises to 44-48%) | low | trivial (run setting) | identical solve; all gates unchanged |
| 2 | Predictor-only tolerance 1e-8 | 4% B (3% S); 1e-6: 5% | medium (state feeds corrector) | simple | not identical; justified tolerance needed; gate unchanged (corrector decides) |
| 3 | Gate-tied absolute stopping rule (about 1e-10 abs residual, margin 10) | 2.5% B, 1% S (cycles -15% / -7%) | medium | simple | not identical; gate margin kept at 10; needs argument |
| 4 | Skip rebuild on a regrid that leaves the layout unchanged (7.2) | 2.5% B, 0 S | low | simple | identical |
| 5 | AMReX patch: skip first applyBC of zero-guess smooth | 1.3% (estimate) | low | small upstream patch | identical |
| 6 | AMReX patch: cache coarsened BoxArray in makeCoarseAmr | up to 1.5% (estimate) | low | small upstream patch | identical |
| 7 | AMReX patch: fuse cor zeroing into first sweep | about 1.3% (estimate) | medium | medium | identical in principle |
| 8 | Larger fine boxes (maxsize 16-32) | pressure 5-13% (section 7.6) | layout change | trivial | reduction order only; share rises |
| 9 | Coarsened fine AMR levels (agglomerated MG on fine levels) | 5-8% (estimate, 30-40% fewer cycles) | medium | high | not identical; gates unchanged if same operator |
| 10 | `face_gradient_composite`, pre-solve, H copies (7.8) | 3-4% (estimate) | low | medium | identical |
| 11 | GPU MLMG | about 0% on this size; 3-4x MLMG at 2M+ cells (estimate) | high | high | tolerance-level |
| - | Fine-level sweeps, FMG, W-cycle, C/F order 3, Krylov, warm start (bottom variants: not measured, erratum 2026-10-21) | 0 or negative (measured) | - | - | order 3 fails the gate |

### 8.9 Does the combination reach 30%?
Share target: B needs pressure at about 6.4 s against 14.9 s of other work (currently 10.4-11.5 s). Items 2+3+4+5-7+10 sum to about 12-14% of t_loop in B (pressure to about 7.5-8 s; share about 33-35% at 8 threads); in S about 6-8% (share about 30-31% only with the 1e-6 predictor, 32-33% with 1e-8). **The 30% share is reachable only marginally (S with the aggressive predictor) and not robustly; B stops near 33-35% without the coarsened-fine-level work (item 9, -5 to -8%, giving about 30-32%).** Item 1 changes the picture differently: it cuts t_loop by a quarter but raises the share, so share is the wrong target if wall time is the goal. Cheapest to start now: item 1 (setting), item 4, items 5-6 (small upstream patches, bit-identical), then the predictor tolerance study.
