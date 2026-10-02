# Generator projection: % of modelled time translated and bitwise-tested, weeks 1 to 6

**Everything below is modelled, not measured.** Shares are the inventory's modelled `est_pct_of_total_time` per loop (`docs/inventory/gpu_generator_loop_classes.csv`, source baseline `36975d7`); efforts are estimates; the schedule is a simulation. "Tested" means the generated kernel equals the verbatim upstream loop text bit for bit on the host compiler, not that it has run on a GPU. Weeks count from today: "week k" is the end of week k; "start" is today.

**Assumption change (owner decision, `generator-decisions.md`, decision A): the projection now assumes "bounds change not done".** The 4 to 5 day shared change and the 0.5 day rebase per other engineer are removed (section 7 lists which rows move).

## 1. Starting point and denominators
- **Today (verified by re-running the coverage tool into scratch; docs untouched): 45 loops, 2.479% translated and tested.** The quoted 44 loops, 2.21% is the state before `CHECK_DIVERGENCE` L0363 (0.272): 2.479 - 0.272 = 2.207. The rerun reproduces the tested numbers exactly; 6 loops (0.997%) flip between "classified translatable" and "not translatable" depending on the front-end state of the worktree, which changes no tested number. 54 kernels exist; 8 sit inside the partly translated nests.
- **Headline basis (same as 2.21%): % of total modelled run time** (loop shares sum to 93.298 over 829 loops). 2.479% is on this basis.
- **Non-geometry basis:** 69.444% of total (451 loops). Convert by dividing by 0.69444 (x1.440): today 3.57%. Eligible basis (231 loops, 49.177%): x2.033, today 5.04%. Geometry (378 loops, 23.854%) is outside both.
- Remaining by tier (loops, % of total): T1 wall loops 19 open + 5 partly = 23.5 (8 tested, 0.369); T3 13 open + 1 partly = 8.2; T4 33 open, 7.6 (3 tested, 0.678); T2 1, 0.572; not translatable overall 286 loops, 38.567; partly translated nests 6 loops, 21.628 (upper bound, see 5); retired 207 loops 18.332 (counted done only on the 69.444 "with retired" basis of the coverage doc, not here); host-side 13 loops, 1.935 (never).

## 2. (a) Cumulative % of total modelled time, end of week k (low / base / high)
Base = efforts as in appendix A (2 and 3 engineers: any engineer takes any item; 4 engineers: fixed streams of section 3); low = efforts x1.5, sign-off wait 3 weeks, neighbour-mesh (T4) work never opens, overhead +5 points; high = efforts x0.7, sign-off wait 1 week, T4 opens in week 4, overhead -3 points.
| Engineers | Case | start | wk1 | wk2 | wk3 | wk4 | wk5 | wk6 |
|---|---|---|---|---|---|---|---|---|
| 2 | low | 2.5 | 2.5 | 10.2 | 15.5 | 17.0 | 17.5 | 17.5 |
| 2 | **base** | 2.5 | 10.2 | 15.9 | 17.3 | 17.5 | 26.9 | **29.1** |
| 2 | high | 2.5 | 11.6 | 18.0 | 20.1 | 29.1 | 31.3 | 33.7 |
| 3 | low | 2.5 | 3.9 | 14.0 | 17.0 | 18.6 | 18.8 | 25.0 |
| 3 | **base** | 2.5 | 11.6 | 17.3 | 18.8 | 25.2 | 30.2 | **33.1** |
| 3 | high | 2.5 | 14.7 | 18.7 | 21.2 | 30.8 | 34.0 | 35.1 |
| 4 (fixed split, section 3) | low | 2.5 | 4.3 | 12.5 | 17.2 | 19.0 | 19.8 | 20.8 |
| 4 (fixed split) | **base** | 2.5 | 4.6 | 16.9 | 19.6 | 20.6 | 21.8 | **27.6** |
| 4 (fixed split) | high | 2.5 | 14.1 | 19.2 | 21.5 | 27.7 | 31.5 | 36.3 |
| 4 (reference, any engineer on any item) | base | 2.5 | 7.8 | 18.4 | 20.7 | 26.2 | 31.5 | 32.1 |
Base on the non-geometry basis (x1.440): week 6 = 41.9% (2 eng.), 47.7% (3), 39.7% (4, fixed split); start 3.57%. Loops tested at week 6 (base): 69 / 82 / 85 (start 45); 4-engineer fixed split: 48, 54, 62, 73, 79, 85 at weeks 1 to 6.
**Read this honestly.** (1) The curve is a staircase: 12 loops of at least 1% hold 31.6 of the 44.7 points on the work list. (2) In every case 11.6 to 19.7 points (16.4 in the 4-engineer base) come from completing the 6 partly translated nests (`L0365` 7.679, `L0880` 4.832, `L0401` 3.260, `L0369` 2.399, `L0877` 1.967, `L0381` 1.491). Their share is an upper bound on hot time: the cell sub-nests inside are already translated and tested (8 kernels), so crediting the full nest share on completion moves the metric more than it moves real coverage. (3) Without the nests, base week 6 adds only 6.9 / 10.9 / 8.7 points above start (2 / 3 / 4-fixed engineers), about 0.4 to 0.6 points per engineer-week, which matches the recent measured pace (appendix B). (4) Differences of about 2 points are inside the noise of the greedy order. The fixed 4-engineer split ends below the unconstrained 4 and below 3 engineers because the work is not spread evenly across the streams (section 3), not because of negative scaling.

## 3. The four-engineer split (fixed) and its per-stream load
| Stream | Scope | Tasks | Assigned effort (d) | Modelled pts on its list | Done by week 6 (base, pts) |
|---|---|---|---|---|---|
| 1 GPU Generator Engineer | periodic-step cell loops, reductions and `CYCLE`, function callees, edge nests (accepted-untested, front-end gaps, `L1355`, MMS refs) | 110 | 80.7 | 5.89 | 4.95 |
| 2 Legacy Mapper | the 6 partly translated nests, mass-flux wall nest, pointer, array-constructor and alias machinery (`INF2`, `L0882`, `L0403`, `L0398`, `L0405`, `L0399`, `L0406`, `L0384`, `L0876`) | 15 | 48.5 | 24.33 | 16.98 |
| 3 GPU Wall Loops Engineer | remaining T1 wall loops, T3 hand ports (`WALL_BC` chain), cell-to-wall list, `L1272` | 25 | 87.0 | 8.43 | 1.73 |
| 4 GPU Mesh Data Loops Engineer | T4 loops, T2 (the bounds change is deferred, see section 7) | 31 | 106.0 | 6.01 | 1.50 |
Per-stream effort completed by week 6 at the 78% effective rate (guess overheads below): stream 1 18.5 d, 2 20.0 d, 3 11.0 d, 4 10.5 d. Findings: **stream 2 is the critical path** (24.3 of the 44.7 listed points and 48.5 d against about 23 d of capacity in 6 weeks), so the headline is set by how fast the nests go; stream 4 has almost nothing to do before the neighbour-mesh layout exists (base: week 5), and stream 3 is capped by the serial `WALL_BC` chain. Moving the three `DIVERGENCE_PART_1` nests (`L0365`, `L0369`, `L0381`, 11.6 points, 9 d) from stream 2 to stream 1 changes week 6 to 30.9 (24.0 to 35.9) with a weaker start (2.8 at week 1); a stream-4 engineer lent to stream 3 or 2 before the layout exists would help more than a fifth engineer. Bounds change in stream 4 assumed to go ahead at the start of week 1 (base and high) or week 3 (low); it does not move the percentage.

## 3b. Saturation (week 6, base; low to high in brackets), engineers free to take any item
1 eng. 21.5 (16.1 to 20.7); 2: 29.1 (17.5 to 33.7); 3: 33.1 (25.0 to 35.1); 4: 32.1 (25.8 to 35.4); 5: 32.9; 6: 35.8. **The knee is at 3 engineers when work is free to move; the fixed 4-stream split sits at 27.6.** Why a fourth stops helping: after about week 4 the remaining work is on serial chains or behind gates.
- `GET_SCALAR_FACE_VALUE` explicit-shape rewrite (3 d) -> `L0882`/`L0403` (2 to 2.5 d) -> `L0880`/`L0401` (6 d each): at least 11 working days, then 2 nests in parallel.
- Cell-to-wall list (3 d) -> race family (`L0375`, `L0394`, `L1358`, `L1359`, `L0876`, `L0877`) -> domain-lead sign-off (1 to 3 weeks): credit cannot arrive before the sign-off.
- `WALL_BC` hand port: `L1486` 6 d -> `L1488` 15 d (21 d, over 4 weeks) and `SURFACE_HEAT_TRANSFER` L1485 15 d after `L1486`; one engineer each, cannot be split by loop.
- T4 neighbour-mesh loops (29 on the list, 5.4 points, 98 d): wait for the exchange-buffer layout; not generator work until it exists.
- The shared bounds change is deferred (owner decision); it no longer appears as a work item (section 7).
Coordination overhead (guess, applied as lost time): 2 eng. 8%, 3 eng. 15%, 4 eng. 22% (merge conflicts in `s5gen.py` and the sidecar, full-suite runtime of minutes per run, review waits).

## 4. (b) Order of work that maximises % per week, and what blocks each step (stream in brackets)
Points per engineer-day = modelled share / estimated effort (appendix A).
1. Week 1: `L0365` completion (stream 2) (7.679, 4 d, 1.9 pts/day; `INTERPOLATE1D_UNIFORM` and `TENSOR_DIFFUSIVITY_MODEL` sub-nests remain, read at divg.f90:133-165, 187); `L1355` (1.467, 1.5 d; live-out scalar `DUDX`, a `PRIVATE` fix; stream 1); `L1272` (1.325, 5 d, stream 3). Block: none.
2. Weeks 1 to 2: `L0369` (2.399, 3 d), `L0381` (1.491, 2 d; the `IF (CC_IBM) CALL` hook stays a host call), `L1365`, `L1315`-`L1317`, `L0384`, MMS-function loops, then the 36 accepted-untested loops (0.878 points, 14 d, mechanical). Block: none; the MMS loops are debug-only code.
3. Weeks 2 to 3: shared scalar face-value helper (3 d) and cell-to-wall list (3 d), then `L0882`, `L0403`. Block: pointer-dummy rewrite of `GET_SCALAR_FACE_VALUE`; decision on the wall list (c4).
4. Weeks 3 to 5: `L0880` (4.832, 6 d) and `L0401` (3.260, 6 d) in parallel; `L1272` (1.325, 5 d). Block: step 3.
5. Weeks 3 to 6: race family and `L0877` (1.967 with `L0876`). Block: domain-lead sign-off; exact-sum choice for zone sums (`L0386`, 0.022).
6. Week 5 on: T4 neighbour-mesh loops and `WALL_BC` chain. Block: exchange-buffer layout; chain length.

## 5. (c) Decisions needed from the owner, and the week each gates
| Decision | Needed by | What it gates |
|---|---|---|
| Lower-bound argument change (`ILO..KHI`): **deferred by owner decision; nobody starts it; asked again only when a loop truly needs it** | not scheduled | multi-box-per-mesh runs with lower corner not 1 then need the driver to hand each box as its own array with local origin 1 (section 7); no modelled % depends on it |
| Wall list per box vs per mesh; who writes `UVW_SAVE` (host `MATCH_VELOCITY`, ccib in the deferred bucket) | week 1 | `L0398`, `L0405` (0.30 points) and sub-box tests of wall kernels; no effect on single-box bitwise % |
| Domain-lead sign-offs for the races, reductions and pointers (review table in `docs/upstream-patches/README.md`) | requests open now; answers by week 2 (base) | credit of about 2.8 points (race family 0.804, zone sums 0.022, `L0877` 1.967) |
| Exact-sum choice: **decided** - FDS order by default (GPU terms, serial adds); optional switch for fixed-point sums, off by default, tested separately (`generator-decisions.md`, decision B) | done | zone sums and mass sums in the GPU acceptance; under 0.05 points of modelled time; the switch is GPU Generator Engineer work |
| Upstream patch files (markers, `UNIQUE` assertions, the `BC_INDEX`/`B1_INDEX` check at divg.f90:1581): the owner commits | first batch by week 2 | any rewrite that edits upstream source; none is on the single-box bitwise path |
| Exchange-buffer layout (neighbour-mesh data) | decision by week 3 to open T4 in week 5 | 5.4 points of T4 work |
| GPU validation slots for the Integration Lead | one slot every second week from week 2 (guess) | converts "tested on host" into "run on GPU" (see risk) |

## 6. (d) Main risks to the projection
- **Nest completion is the whole curve.** Effort for the 6 nests (26 d in total) is mostly a guess: only `L0365` was read in detail. If they take twice as long, the base 2-engineer week-6 value falls toward the low band (17.5%).
- **GPU:** 20 of the 54 kernels (rounds 4 to 7) have never run on nvfortran; only 34 were run, and the earlier GPU report shows the host-compiled driver failing 33 of 408 cases on 3 kernels under nvfortran while passing under gfortran. "Tested" is host bitwise; each GPU slot can reopen kernels.
- **Synthetic inputs only:** the bitwise tests use random fields, no real FDS fields.
- **Modelled time is not hot time.** The shares come from a lexical score model with one timer per file. 99 of 181 listed loops are under 0.1% (2.4 points, 120 d). 1.18 points of the list are MMS debug loops (cold in production, 8 loops, 4 of them dead code). Set-up routines (`REASSIGN_WALL_CELLS` 2.037) are left out of the plan because their modelled share is not credible per step.
- **Throughput:** recent rounds on new mechanisms delivered 0.1 to 0.7 points each; the base case needs the batch-round pace (9 to 18 kernels per round) for the cheap items and the nest completions to work as estimated.

## Appendix A: effort assumptions and derivations (efforts are guesses unless marked derived)
| Work type | Loops (share, %) | Effort | Basis |
|---|---|---|---|
| Accepted by the front end, untested | 36 (0.878) | 0.4 d each | derived from rounds of 18, 9, 13 kernels with 1 round = 5 engineer-days (guess) = 0.28 to 0.55 d per kernel |
| Cheap generator gap (live-out scalar, rank, `policy.arrays`, function reference) | 72 (5.617) | 0.8 d (1.0 for MMS refs, 1.5 `L1355`, 2.0 `L0384`) | guess |
| Wall-table T1 loop | 7 (1.75) | 0.5 to 1.5 d; `L1272` 5 d | guess; round 3 built the wall tables with 9 kernels from 5 loops |
| Pointer rewrite, shared helper | 4 + helper (1.07) | helper 3 d, then 2 to 2.5 d each | guess |
| Reduction with tie rule | 1 (0.022 left) | 4 d | rounds 5 and 7: 4 kernels plus two builders |
| Race rewrite needing sign-off | 7 (0.804) + `L0877` | 1 to 4 d plus wait: low 3, base 2, high 1 weeks | guess |
| Partly translated nest completion | 6 (21.628) | 2 to 6 d each (26 d) | guess |
| Callee hand port (`WALL_BC` chain) | 13 T3 (6.846) | 6, 4, 15, 8 d; other T3 2 d | guess |
| Neighbour-mesh / `ONE_D` / ragged table, gated | 30 (6.013) | 3 to 8 d (98 + 8 d) | guess |
| Geometry table | outside the denominator | 0 | out of scope |
| Shared bounds change | 0 | 0 (deferred; was 4.5 d + 0.5 d rebase per other engineer) | owner decision |
Reconciliation (derived): tested 2.479 + work list 44.661 + left out 3.916 (I/O, sequential control and particle loops 1.879; set-up routines 2.037) + other host-side 0.056 = 51.112 = 69.444 - 18.332 retired. The work list holds 181 tasks and 327 engineer-days at base; capacity over 6 weeks is 55 / 77 / 94 days (2 / 3 / 4 engineers after overhead), so the plan covers under 30% of the listed effort.

## Appendix B: throughput evidence (derived from the repository history and the coverage tool)
Rounds: 7 original, 18 (round 2), 9 kernels from 5 loops (round 3), 13 (round 4), 3 (round 5, `CHECK_STABILITY`, 0.125%), 3 (round 6, edge nests, 0.678%), 1 (round 7, `CHECK_DIVERGENCE`, 0.272%). Rounds 5 to 7 added 1.075 points together (0.36 per round); rounds 1 to 4 added 1.404 in total. The base case's new-work-only gain of 0.4 to 0.6 points per engineer-week is on that scale; the nest completions are the part with no precedent.

## 7. Change record: bounds change deferred (assumption "bounds change not done")
- **Removed from the model:** 4.5 d of work in one engineer's week 1 and 0.5 d of rebase for each other engineer (2 engineers: 5.0 d in total; 3: 5.5 d; 4: 6.0 d). Stream 4's assigned effort falls from 110.5 d to 106.0 d (sections 3, 4, 5 and appendix A updated).
- **Consequence stated plainly:** the kernels index one mesh with first cell 1. For a mesh split into several boxes with a lower corner not 1, the driver must hand each box to the kernel as its own array with local origin 1 (box-local arrays, box-local wall tables, sliced 1-D metric arrays, box-local coordinates); the kernels cannot be told the box position. Single-box runs are unaffected, and no modelled percentage depends on the change (it moved no % itself).
- **Week-6 numbers: not re-simulated.** The scheduling simulation behind section 2 was not kept as a script in the repository, so the table in section 2 is **unchanged and still reflects the old assumption**; the figures below are a bounded estimate, not a new run. Freed capacity is 5.0 / 5.5 / 6.0 engineer-days; at the throughput of appendix B (0.4 to 0.6 points per engineer-week, about 0.08 to 0.12 points per day) that is at most +0.4 to +0.7 points at week 6 for the 2- and 3-engineer rows (base 29.1 and 33.1), and it could pull a step of the staircase about a week earlier at weeks 1 to 2. Rows that cannot move: the 4-engineer fixed split at week 6 (27.6 base; stream 2, not stream 4, is the critical path, and stream 4 has no other work before the neighbour-mesh layout exists), and all "low" columns where T4 never opens. Row to re-run when the simulation is rebuilt: 2 and 3 engineers (base and high), weeks 1 to 3.
