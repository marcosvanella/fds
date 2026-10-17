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
- BoomerAMG bottom is 10% of pressure in B (setup 0.92 s, solve 0.17 s: **the setup is rebuilt at every V-cycle, 5x the solve**) and 16% in S later (setup 0.54, solve 0.15, plus 0.41 s of BiCGStab fallback in 133 solves). Total 4.6% of t_loop in B, 6.6% in S later (est. from the TinyProfiler difference).
- The coarsest grid is 9x6x6 (324 cells). `max_coarsening_level` can only be lowered at runtime (it is already 2, the limit of the 36x24x24 domain): level 1 gave pressure 19.7 s (+65%), level 0 79.8 s (+570%). No larger coarsening is possible.
- Bottom `pb_mlmg_bottom.solver=smoother` or `bicgstab`: cycles unchanged (20.2) and time unchanged within noise (3 paired repeats: baseline 9.37 / 9.51 / 9.59 s MLMG, smoother 9.43 / 9.51 / 9.54 s). **No gain from the runtime option** (the setting may not alter the path that dominates; not investigated). A persistent BoomerAMG setup would be a code change: est. up to -0.9 s in B (3.6% of t_loop), -0.5 s in S later (3%).

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
| - | V(1,1)/(3,3), bottom smoother, Krylov, lower coarsening, warm start in S | 0 or negative | measured | - | - |

Predicted outcome: items 1-6 sum to about 19% of t_loop in B at best (estimate), taking pressure from 46% to roughly 34% share ((11.9 - 4.9)/(25.9 - 4.9)); 30% is not reached by solver tuning alone. A further 3 s would need about 25% fewer cycles at the same tolerance (a better C/F treatment or smoother; the 0.23 asymptotic factor, against 0.14 in cycle 1, is the target) or a cheaper residual/C/F fill path (items 3-6 of 7.3 are 52% of pressure). In S the post-regrid share is 41% and the same list gives about 33% (estimate).
