# 09 · Device-tier flip-budget gate for the 1-D wall solve port

Owner: AMR Solid Phase Lead · Status: **v0.1, definition for review** · Applies to units U17 and U18 of `06-solve-port-test-design.md` (whole thick-wall and thin-wall solve) on the device tier of gate G15.
Reference revision: FireX `36975d765f`; every `file:line` is from that revision (the working tree is 9 to 12 lines higher). `wall.f90` unless a file is named. Nothing was run to write this page: statements about upstream come from reading the pinned text, and the binomial figures in section 6 come from an exact calculation that the harness will repeat in its own self-tested code.
Rulings and positions used: D-070 (2 ulp per value for libm kernels on the device), D-072 (sub-step cap 10**6, status 301, non-finite input only; back-side snapshot of four fields; flip budget rules), the V&V position of `docs/vv/test-plan.md` section 5.10 (six conditions), D-075 (the driver aborts with the report when a table check refuses).

## 0. Summary

1. A **flip** is a record-call whose discrete decisions differ between host and device. Five decision classes are defined (section 3), each with its place in the solve.
2. The rate is counted only over **libm-exposed, distinct-input** record-calls, separately for real captured calls (P1) and random calls (P2). Adversarial calls (P3) are reported, never counted in a rate.
3. The gate statistic is the exact (Clopper-Pearson) one-sided bound at 95%. Verdict PASS if the upper bound is at most 1e-4, FAIL if the lower bound is above 1e-4, otherwise INCONCLUSIVE. **Minimum N to demonstrate 1e-4 with zero flips: 29,956 calls per population at 95% (46,050 at 99%).** Planned N is 10**5 per population, which tolerates up to 4 flips.
4. **Physical bound:** only partly derivable from the code. Two stopping tolerances give real bounds (the temperature extraction after renoding: 1e-4 K; the oxygen Newton loop: 1e-6 kg/(m² s)), plus the layer-removal thresholds. The sub-step count and most remesh decisions have no code bound; the solver has no error estimator. The T1-class value 1e-10 proposed by V&V is below the solver's own stopping tolerances and cannot be derived. Section 7 replaces it with a two-level rule.
5. **Cap check:** finite inputs never give status 301, no gate case comes within a factor of 10 of the cap, and a poisoned `DT_BC` gives 301 within the cap with byte-identical text on both paths. Reading shows a NaN temperature does **not** in general make the loop run for ever, so the poison matrix is classified from the reference run (section 8).
6. A budget overrun goes to the Solid Phase Lead and the Architect; tolerances are never widened (section 9).

## 1. Definitions

- **Record-call.** One execution of the whole solve (`SOLID_HEAT_TRANSFER`, 1809-3154) for one record, thick-wall entry (`WALL_INDEX`) or thin-wall entry (`THIN_WALL_INDEX`), for one `DT_BC`.
- **Frozen input.** The complete bit image the call reads: the record row (`ONE_D` fields including `TMP(0:NWP+1)`, `MATL_COMP%RHO`, `X`, `DX_OLD`, `N_LAYER_CELLS`, `LAYER_THICKNESS`, `LAYER_THICKNESS_OLD`, `REMESH_NWP`, `SMALLEST_CELL_SIZE`, `DDSUM`, particle fields), the `B1` and `B2` fields read, the gas-side values read by the pyrolysis loop (`B1%ZZ_G`, the gas cell temperature at 3041, `B1%HEAT_TRANS_COEF`), the back-side snapshot, `DT_BC`, `T`, `ICYC`, and the module data (`WALL_INCREMENT`, `CHECK_FO`, `OXPYRO_MODEL`, `TMPMIN`, `TMPMAX`, `I_MAX_TEMP`). Host and device get the same bits. Nothing accumulates across calls: every call starts from its own frozen image (V&V condition 1).
- **Back-side snapshot (D-072).** Four fields only, the ones read at 2165-2169: `TMP_G`, `HEAT_TRANS_COEF`, `Q_RAD_IN` and the sum of `LP_CPUA`. The surface constants of the back side (2174-2185) are case data, not snapshot. A further cross-side read needs a new ruling.
- **Flip-call.** A record-call in which at least one decision trace entry of section 3 differs between host and device. A call with several differing decisions, or a cascade, is **one** flip-call. Per-class counts are reported as well.
- **Exposed call.** A call whose decision chain contains a libm call (static tag, section 4.3). Only exposed calls enter the denominator of the rate.
- **Host** = the verbatim upstream text compiled with `-ffp-contract=off`, serial, with decision-trace stores added by the harness. **Device** = the generated kernel on the GPU (`-gpu=nofma`).

## 2. What the gate compares, in order

For every record-call:
1. Run host and device on the same bit image; both return outputs and the decision trace.
2. **Traces equal** (all five classes): compare every output value under the D-070 rule (2 ulp per value, `bitcmp.py --ulp` semantics; +0 and -0 distance 0; NaN position and infinities identical). A value outside 2 ulp with equal traces is a **failure**, never a flip. The first action on such a failure is to look for an unlisted decision variable (section 3 rule).
3. **Traces differ:** the call is a flip-call. Run the counterfactual check of section 5.3, report it (section 5.2), and count it (section 6).
4. Non-exposed calls are held to the stronger rule: traces equal **and** outputs bitwise (the non-libm rule of D-070). A flip or any value difference in a non-exposed call is a failure, not budgeted.

## 3. Decision variables and where they are

The solve has no tolerance-controlled outer loop. Its discrete outcomes are the following. "Trace" is what the harness records at the site (outcome and, as a hex float, the real-valued argument that was compared).

| Class | Decision | Site (survey revision) | Trace entry | libm in the chain |
|---|---|---|---|---|
| **D1** sub-step count | Divisor per sub-step `MIN(NINT(TIME_STEP_FACTOR*WALL_INCREMENT),MAX(1,NINT(TMP_RATIO)))`, only when `ICYC>WALL_INCREMENT` (guard 2306). `TMP_RATIO` at 2322 from the explicit estimate 2314-2321. Consequences: `DT_BC_SUB` 2324-2325, rebuild of `Q_S` 2327, exit test 2949 and `B1%N_SUBSTEPS` (set to 1 at 2032, incremented 2951). | wall.f90:2306-2337, 2949-2951 | per sub-step: divisor, `TMP_RATIO`; at exit: `N_SUBSTEPS` | `Q_S` from pyrolysis, `R_S**I_GRAD`, back HTC |
| **D2** node and layer structure | Cell disappears 2448; size ratio against `REMESH_RATIO` 2455-2458; temperature check against `NODE_RDT` 2461-2463 and 2467-2469; thin layer or low mass 2508; delamination 2520-2521; thickness drop 2533-2537; remesh trigger 2583-2587; one-cell rules 2610-2613; new layer cell count `GET_N_LAYER_CELLS` (func.f90:5743-5784, `NINT` 5760, loop test 5777); extra shrink 2643-2660; capacity check 2668; zero thickness 2682-2697; final `NWP` and `N_LAYER_CELLS` 2801-2802. | wall.f90:2445-2475, 2502-2680, 2801-2802 | per remesh: `REMESH_LAYER`, `TMP_CHECK`, `CELL_ZERO_CELL` flags, `N_LAYER_CELLS_NEW`, `NWP_NEW`; at exit: `NWP`, `N_LAYER_CELLS(1:N_LAYERS)`, burn-away flag | `REGRID_FACTOR` from pyrolysis, `R_S_NEW` power 2437, `STRETCH_FACTOR**` at 5775 |
| **D3** oxygen Newton iterations | Bounded Newton on `Y_O2_F`: loop `ITER=1,MAX_ITER` (3048), `MAX_ITER=20` when `OXPYRO_MODEL` (3031, 3036), exit when `ABS(M_DOT_ERROR)<M_DOT_ERROR_TOL=1.E-6` (3003, 3126). A non-converged loop ends silently at 3150 and does not write `B2%Y_O2_ITER` (3128), so the **final `ITER` and a converged flag** are traced, not the field. | wall.f90:3031-3150 | per `PERFORM_PYROLYSIS` call: final `ITER`, converged flag, `M_DOT_ERROR` | `EXP`, `**` in `PYROLYSIS`, HTC |
| **D4** temperature extraction iterations | After a regrid, node temperature from enthalpy by Newton: stop when `ABS(TMP-T_NODE)<0.0001` (2784), after 20 iterations average and stop (2785-2788). | wall.f90:2758-2790 | per regridded node: final `ITER`, fallback flag | enthalpy tables are linear interpolation (no libm); inputs come from D2 and pyrolysis |
| **D5** step-function table lookups | Properties looked up by integer temperature **without** interpolation: `K_S(ITMP)`, `C_S(ITMP)` (2075-2081, 2823-2829), `H_R(J,INT(TMP_F))` (3421), `H_R(1,NINT(T_BOIL_EFF))` (3448), `H_R(J,NINT(TMP_S))` (3564-3565). The interpolated lookups (2713-2715, 2768-2775, 3443-3446) are continuous and are not decisions. | as listed | per evaluation: the index; stored as a hash chain per call (the full index list only on replay) | indirect: temperatures carry libm differences |

**Not decisions:** `DT_BC_SUB` itself (follows from D1), the `MIN` with the remaining time (2325) and the exit test (2949). The exit test uses the absolute tolerance `TWENTY_EPSILON_EB` = 20*epsilon (prec.f90:43), so a tail sub-step of size near round-off depends on rounding of `T_BC_SUB`; it can only differ if D1 differs, because the arithmetic that builds `T_BC_SUB` has no libm call.

**Unlisted-decision rule.** If a call has equal traces and a value mismatch above 2 ulp, the gate fails and the investigation either finds a defect or adds a decision site to this table (with a trace entry) by a recorded change. The table is the only place flips are defined.

## 4. The record-call population

### 4.1 Populations

| Pop | Source | Use |
|---|---|---|
| **P1 real** | Capture-and-replay (06 section 5.7, section 11): a patched host build dumps the frozen image of every record-call of `heat_conduction_a`, `energy_budget_solid`, `back_wall_test`, `surf_mass_vent_char_cart_fuel`, then of every other supported input that has a 1-D wall. Every step is captured (not every k-th) so the trajectory is covered. | Gate (rate) |
| **P2 random** | `wall_solve_cases.py` (06 section 6), seeded, draws from the physical ranges of the axis table, no oversampling of thresholds. Each call starts from a state reached by running a few host calls (not from an arbitrary random state), so states are reachable. | Gate (rate) |
| **P3 adversarial** | Constructed ties: `TMP_RATIO` within a few ulp of k+0.5, `RMR` at `REMESH_RATIO`, `M_DOT_ERROR` at 1e-6, temperature at an integer, layer mass at its minimum. | Harness self-test and sensitivity. Flips are reported; **never counted in a rate**. |

### 4.2 Distinct inputs

N counts **distinct frozen images** (SHA-256 of the bit image). Identical records in a uniform wall (a slab of 1000 equal cells) are one call, not 1000. Without this the Bernoulli assumption of section 6 would be false and N could be inflated.

### 4.3 Exposure tag

A call is exposed if any of these holds (the tag is computed in Python from the record features before the run and written to the case file):
- **PYR:** predicted pyrolysis with at least one reaction (`EXP`, `**` in `PYROLYSIS` 3195-3605; D1, D2, D3, D5).
- **GEOM:** `I_GRAD>1` (variable integer power `R_S**I_GRAD` at 2370, 2437, 2709; for `I_GRAD=1` the power is exact on both sides).
- **STRETCH:** a remesh path that calls `GET_N_LAYER_CELLS` with `STRETCH_FACTOR` different from 1 (variable integer power at func.f90:5775 and wall.f90:2648-2656). `SQRT` is correctly rounded and not libm.
- **HTC:** the HTC callee on the path: `VOID` back side (2146), boundary fuel model convection (2235-2256), oxygen model `H_MASS` (3043). The front `HEAT_TRANS_COEF` is an input and not part of this.

The tags are reported with every run. A plain cartesian slab with tabulated properties and no pyrolysis is **not exposed** (everything on its path is `+ - * /` and integer powers), so it can only be bitwise; it adds nothing to the rate. This is deliberate: a pooled rate over all calls could be driven under any budget by adding unexposed calls.

### 4.4 Size

Planned N: at least 10**5 exposed distinct calls in P2. For P1 the number is whatever the supported cases give; it is reported, and section 6 says what follows if it is short. Cost is not the limit on the device (a call is small); the host reference run, the case-file traffic and the counterfactual runs are.

### 4.5 Completeness checks of the frozen image (run once per case family)

- Every field read by the solve is in the image (the read table of 06 section 5.4, checked in both directions).
- **Snapshot closure (D-072):** every field of the back-side record other than the four listed is set to a signalling-NaN pattern; outputs must be bitwise unchanged on host and device. A change means a fifth cross-side read exists and needs a ruling.
- Instrumentation neutrality: the traced build and the untraced build give bitwise equal outputs on the host in all six flag sets.

## 5. Counting and reporting flips

### 5.1 Trace comparison

Traces are sequences per class (section 3). The device returns them in a side array of fixed size (length per class set at build; an overflow counter, nonzero = failure of the harness, not a flip). Equality is exact on outcomes. The real-valued arguments are carried to the report only.

### 5.2 Report format

One JSON object per flip-call in `flips.jsonl` of the run folder, one line of text in the log, and the frozen input in `flips/<pop>_<case>_<call>.case` (16-digit hex bit patterns, one value per line, readable by the Fortran drivers; a Python-side hex-float view is printed by the report tool). Fields:

```
{"v":1,"pop":"P2","case":"rand-0457","call":31822,"input_hash":"<sha256>",
 "input_file":"flips/P2_rand-0457_31822.case","exposure":["PYR","STRETCH"],
 "kernel":"<sidecar kernel id>","device":"<gpu, compiler, flags>","host":"<compiler, flags>",
 "first_class":"D1","site":"wall.f90:2324","eval":3,
 "host":{"outcome":2,"arg":"0x1.8000000000001p+0"},
 "dev":{"outcome":1,"arg":"0x1.7ffffffffffffp+0"},"arg_gap_ulp":2,
 "classes_differing":["D1"],"n_substeps":{"host":7,"dev":6},"nwp":{"host":24,"dev":24},
 "inputs":{"DT_BC":"0x1.999999999999ap-4","TMP_F":"0x1.2c...p+8","TMP_B":"...","TMP_G":"...",
           "HEAT_TRANS_COEF":"...","Q_RAD_IN":"...","back":{"TMP_G":"...","HTC":"...","Q_RAD_IN":"...","LP_CPUA":"..."},
           "DELTA_TMP_MAX":"0x1.9p+4","ICYC":412,"WALL_INCREMENT":1,"NWP":24,"N_LAYERS":2},
 "cf":{"forced":["D1@3"],"levelA_max_ulp":1,"delta":{"TMP_F":{"abs":"<hex>","rel":"<hex>"},"Q_CON_F":{"abs":"<hex>","rel":"<hex>"}},
       "bound":"empirical E_D1 (unsigned)","verdict":"RECORDED"}}
```
(All values above are illustrative.) The summary file `flip_summary.json` has, per population: N distinct exposed calls, N all calls, k, per-class k, per-exposure-tag N and k, per-case k, the upper and lower bounds, the verdict, and the largest `Delta_cf` per class.

### 5.3 Counterfactual host check (the core of the flipped-call rule)

For a flip-call, run the host reference again with the **first differing decision forced to the device outcome** (the harness replaces that one outcome in the verbatim text; all later decisions run free). Repeat while the traces still differ (at most 8 forcings; more is reported as a cascade and fails the call). Then:
- **Level A (port correctness):** the device outputs equal the forced-host outputs within 2 ulp per value. A flipped call is only a branch change if this holds. Failure means the device diverges for another reason (a hidden decision, an FMA, a defect) and is a failure, not a flip.
- **Level B (physical size):** `Delta_cf` = host outputs minus forced-host outputs, per output field, abs and rel. This is the physical difference caused by the flip, measured on the host, independent of the device. It is judged by section 7.

## 6. Statistic and confidence bound

Per population p in {P1, P2}: N_p = exposed distinct-input calls compared, k_p = flip-calls among them, budget b = 1e-4 (hypothesis of V&V condition 4, measured and gated).

- Upper bound: `UCB(k,N)` = the rate q at which P(X<=k) = alpha for X ~ Binomial(N,q) (Beta quantile B(1-alpha; k+1, N-k)). For k = 0: `1 - alpha**(1/N)`, close to `-ln(alpha)/N`.
- Lower bound: `LCB(k,N)` = the rate at which P(X>=k) = alpha (Beta quantile B(alpha; k, N-k+1)); 0 for k = 0.
- alpha = 0.05 (one-sided 95%). A 99% level is an option for phase exit (question V1).
- **Verdict:** PASS if `UCB <= b`. FAIL if `LCB > b`. Otherwise INCONCLUSIVE.

Minimum calls (exact binomial, rate 1e-4):

| flips k | min N at 95% | min N at 99% |
|---|---|---|
| 0 | **29,956** | **46,050** |
| 1 | 47,437 | 66,381 |
| 2 | 62,956 | 84,057 |
| 3 | 77,535 | 100,448 |
| 4 | 91,533 | 116,043 |

At N = 10**5 (planned): upper bounds at 95% are 3.00e-5 (k=0), 4.74e-5, 6.30e-5, 7.75e-5, 9.15e-5 (k=4), and 1.05e-4 for k=5, so up to 4 flips pass at 95% and up to 2 at 99%. FAIL needs the lower bound above 1e-4: k >= 7 at N = 29,957, k >= 16 at N = 10**5, k >= 29 at N = 2*10**5.

A plain "at most 1 per 10**4" reading would call N = 10**4 with one flip a pass. It demonstrates nothing: with N = 10**4 even zero flips gives an upper bound of 3.0e-4, above the budget.

**If N is too small.** The best attainable verdict is INCONCLUSIVE (or FAIL if the flips are many). INCONCLUSIVE does not sign the device tier. Actions, in order: (1) extend sampling (P2 can always be extended; P1 by adding supported inputs with walls); (2) if P1 cannot reach N after all supported inputs, report the shortfall to V&V, who decide whether the gate rests on P2 with P1 reported (not automatic); (3) never convert INCONCLUSIVE into PASS by pooling P2 into P1, by adding unexposed or duplicate calls, or by widening b.

**Assumptions stated.** Calls are treated as independent identically distributed Bernoulli trials. Consecutive steps of one record are correlated, so the bound is slightly optimistic for P1; the report lists k per case, and if all flips sit in one case or one exposure tag the verdict is also computed for that stratum. P1 and P2 are never pooled.

## 7. Physical bound for flipped calls

The question from V&V: does the solver give a bound on the temperature and heat-flux difference after a flip, for example from its stopping tolerances? Answer from reading the code: **partly, for some classes; for the most frequent classes only an empirical bound is possible.**

| Class and site | What the code controls | Derived bound on the difference | Status |
|---|---|---|---|
| D4 extra temperature-extraction iteration (2784) | stop when the Newton step is below 1e-4 K; the other side stops one iteration earlier or later, so the two results differ by about one Newton step | `abs(Delta TMP(I)) <= 1e-4 K` at that node (contractive convergence assumed); heat flux through the face follows from 2918-2926: `abs(Delta Q_CON_F) <= HTCF*DT_BC*max abs(Delta TMP_F)` plus the radiation linearisation term | **Derived.** The non-converged fallback (2785-2788) has no bound |
| D3 extra oxygen iteration (3126) | stop when `abs(M_DOT_ERROR) < 1e-6 kg/(m² s)` (3003) | `abs(Delta M_DOT_O2_PP)` of about 2e-6 kg/(m² s), and `abs(Delta Y_O2_F)` about 2e-6/`H_MASS` | **Derived** (absolute, under contraction). A flip between "converged at iteration 20" and "not converged" has no bound (silent exit at 3150) |
| D2 layer removed for thinness or low mass (2508) | thresholds `MIN_LAYER_THICKNESS` (default: the smaller of 1e-6 m and 10% of the layer, read.f90:8909-8911, 9236) and `MIN_LAYER_MASS` (layer density times 1e-12, read.f90:8937, 9235) | removed mass at most `MIN_LAYER_MASS`, removed thickness at most `MIN_LAYER_THICKNESS` | **Derived** from the thresholds |
| D2 renoding temperature check (2461, 2467) | interpolation error compared with `NODE_RDT` = `RENODE_DELTA_T` of the reacting material, the temperature change that moves the reaction rate by `REAC_RATE_DELTA` = 0.15 (read.f90:7615, 7947-7977; 5000 K if no sensitivity) | node temperature difference at most about `NODE_RDT` (a few K for typical activation energies) | **Derived but weak** (orders above 1e-10) |
| D2 delamination (2520), burn-away zero thickness (2682) | physical events at a threshold | none: a whole layer falls off or not | **No bound**; a rate event |
| D2 size ratio `REMESH_RATIO` = 0.25 (2457), cell count from `GET_N_LAYER_CELLS` (5760, 5777), `NWP` | grid quality heuristics | none: the difference is interpolation and grid error | **Empirical only** |
| D1 sub-step count (2324) | step-size heuristic: the explicit estimate of the change per sub-step against `DELTA_TMP_MAX` (default 25 K, read.f90:9179) with at most `TIME_STEP_FACTOR*WALL_INCREMENT` sub-steps (default factor 10, read.f90:9288). The implicit scheme is Crank-Nicolson (2856-2863), second order, with lagged properties | none. The step size changes by a factor of at most 2 at a flip (divisor k against k+1), and the difference is the difference of two time-truncation errors. A bound of the form `37.5 K * z**2/12` with `z` the largest dimensionless step exists in principle but `z` is not limited by the code (only `CHECK_FO=T` limits it, to about 4, which gives tens of K) | **Empirical only** |
| D5 step-function lookup (2828-2829, 3565) | table jump between neighbouring integer temperatures | the relative jump of `K_S`, `C_S`, `H_R` between `T` and `T+1`, computable per material from the tables | **Computable per flip**; the harness prints it, no temperature bound |

Findings:
1. **The solver has no error estimator.** No tolerance in it bounds the effect of a different sub-step count or a different node count. The sub-step rule limits the explicit estimate of the change per sub-step, which is a step-size heuristic, not an accuracy control.
2. **The T1-class value 1e-10 relative cannot be derived and is wrong for the classes that have a bound.** It is far below the solver's own stopping tolerances (1e-4 K is 3e-7 relative at 300 K; 1e-6 kg/(m² s) is 1e-4 relative at 1e-2 kg/(m² s)). Any D3 or D4 flip would exceed it by orders of magnitude, and a D1 flip almost certainly would. It would be a rule that fails correct flips.

**Proposal: two levels instead of one tolerance.**
- **Level A** (section 5.3, every flip): device equals the host forced to the same decisions, within 2 ulp. This is the real port-correctness test and it is stricter than 1e-10.
- **Level B** (physical size, every flip): `Delta_cf` is compared with a class bound. D3, D4, layer-removal and renoding-check flips use the derived bound of the table. D1, `REMESH_RATIO`, cell-count and D5 flips use an **empirical class bound `E_class`** that is measured before the device gate: a *forced-flip survey* on the host over P1 and P2 (force the alternative decision at the decision site of each call that has a near-threshold decision, at least 1000 calls per class where reachable) gives the distribution of `Delta_cf`; `E_class` = twice the largest observed, signed by V&V and the Solid Phase Lead. Until a class bound is signed, flips of that class are recorded and the class has no PASS.
- Flips of delamination and burn-away are counted in the rate and reported, with no level B bound.

## 8. Cap check in the harness (D-072)

**What the cap is.** A counter `NSUB` in the sub-step loop (2044-2953), cap `10**6`, status 301 in `OD_STATUS`; no `WRITE` or `SHUTDOWN` in the kernel; the host builds the message from the status after the pass (06 section 8). It is a deviation from upstream, which has no cap.

**Reading results that shape the test.**
1. For finite input with `CHECK_FO=F` the count is bounded by about `NINT(TIME_STEP_FACTOR*WALL_INCREMENT)+1` per call (each sub-step is at least `DT_BC` divided by that, except the tail), 11 for the defaults. With `CHECK_FO=T`, `DT_BC_SUB` is limited by `DT_FO` = minimum of `DX_S**2*RHO_C_S/K_S` (2096), so the count can be as large as `DT_BC/DT_FO`. For a steel layer of 1e-6 m cells (`DT_FO` about 8e-8 s) and `DT_BC` = 0.1 s that is about 1.2e6, above the cap, with finite input. So "never triggers on finite input" is **a property of the gate cases, not of the code**, unless the kernel enforces "non-finite only" itself (question A1).
2. A NaN or Inf temperature does **not** in general make the loop run for ever. `DT_BC_SUB` is built from `DT_BC` and integer divisors; `TMP_RATIO` is passed through `MAX` and `NINT`, and with gfortran `MAX(x, NaN)` and `MAXVAL` ignore a NaN (to be confirmed on the reference), `NINT` of a non-finite value gives a large negative integer that `MAX(1,...)` clamps to 1. Then the loop ends after one sub-step with NaN in the outputs. The sure non-terminating input is `DT_BC` = NaN (`T_BC_SUB` becomes NaN and 2949 is never true). `DT_BC` = Inf ends after one sub-step. This corrects the statement in 06 section 2, fact 2, which was conditional on a confirmation by the reference. These intrinsics may behave differently in the device compiler, so which poison hangs can differ between host and device; the harness records it.

**The test (all run under the generator lock, quick tier uses a test-build cap of 1000, device tier uses the real cap once).**
1. **Finite inputs never give 301.** Over every call of P1, P2, P3 (non-poison): `OD_STATUS` = 0, and the largest observed `N_SUBSTEPS` is reported and must be at most cap/10. A gate case within a factor of 10 of the cap fails the run and is redesigned (it would hide the regime in which the cap and the physics meet).
2. **Poison matrix.** For each field class of the frozen image (record temperatures, densities, `B1` gas-side values, the four snapshot fields, `DT_BC`, `T`) and each of NaN, +Inf, -Inf, run the reference with the harness-inserted counter (verbatim text plus the counter and nothing else, 06 section 8 item 6) and the kernel. Expected outcome is taken from the reference, not assumed:
   - reference reaches the cap (certain for `DT_BC` = NaN): the kernel must set status 301 within the cap, leave the record untouched except as the reference does, and leave the canary bands intact;
   - reference terminates: the kernel must also terminate with status 0, same `N_SUBSTEPS`, and outputs equal in NaN positions, infinities and finite values (finite values within 2 ulp; NaN payload not compared).
   A mismatch of category between reference and kernel is a failure. The matrix must contain at least one case of each category (coverage gate); the mandatory hang case is `DT_BC` = NaN.
3. **Text.** The host message for status 301 is built by one template recorded in the sidecar (proposal, same style as 300): format `'(A,I0,A,A)'`, text `ERROR(301): 1-D wall solve did not finish within ` cap ` sub-steps (non-finite input) for ` `TRIM(SF%ID)`. The test compares, byte for byte, three strings: (a) from the status path (`OD_STATUS`, `OD_AUX`, surface id), (b) from the CPU path of the AMR build (the kernel's host fallback, flag set `O2omp_off`, and the harness reference with the cap, both through the mock `SHUTDOWN` with the fixed exit code of 300), (c) a golden string typed in the test file. A change to the template is a visible review event.
4. **Several records.** The lowest failing record index is reported; records that did not fail are bitwise equal to the reference; canaries intact for the failing record and its neighbours.
5. **Finite input at the cap (only if question A1 is answered with the predicate form).** With a test cap of 50 and a finite `CHECK_FO=T` case that needs 200 sub-steps: status 0, outputs bitwise equal to the reference. This shows the cap does not engage on finite input. If the Architect rules the cap fires on count alone, this test is replaced by the assertion of item 1 only.
6. **Mutants:** cap removed (the `DT_BC` = NaN case then hangs; run under a timeout, a timeout counts as detection); status set but loop continues; `OD_AUX` wrong; cap compared with `>` instead of `>=` (exact count test at a test cap); predicate inverted (finite input aborts, non-finite input does not). All must be caught.

## 9. When the budget is exceeded

The device tier is **not signed** when any of these happens: the rate verdict is FAIL for P1 or P2; any flip in a non-exposed call; any Level A failure; any Level B bound exceeded where a bound is signed; any cap assertion of section 8 fails; the trace overflow counter is non-zero. INCONCLUSIVE is also not a sign-off (section 6).

Then the harness writes the run report (section 5.2 summary plus the ten flip-calls with the largest `Delta_cf`, the per-site and per-class counts, per-case and per-tag rates) and the run is reported to the **Solid Phase Lead and the Architect**; V&V logs it on the sign-off list. The host tier, which is bitwise, is unaffected. The tolerance is not widened and the budget is not re-read as a point estimate to make a run pass. The responses open to the Solid Phase Lead and the Architect are: (a) find the decision sites with most flips (the report lists them) and make that chain independent of libm by computing the libm quantity on the host and passing it in (geometry powers, rate constants), which changes no upstream arithmetic; (b) run the exposed records (pyrolysis, non-cartesian) on the host in the whole-wall host pass that D-065 allows, keeping the device for the rest; (c) the Architect re-rules the budget, with the measured rate and `Delta_cf` distribution in front of them.

## 10. Harness self-tests (before any device run)

1. **Flip injection.** A host build of the reference in which every libm call result is perturbed by one ulp (a fixed fraction of calls, seeded; the measured device difference is 1 ulp for about 13% of arguments, test-plan 5.10) must produce flips on the P3 tie cases and none on non-exposed calls; the counter, the report format, the counterfactual check and the statistic are checked on these flips. If the P3 cases do not flip, the case set or the comparator is blind and the gate is not trusted.
2. **Pre-device rate estimate (report only).** The same injected build run over P2 (all exposed calls) gives an estimate of the flip rate before a device exists. It is not a gate: it models the device, it is not the device. A rate near or above 1e-4 here is an early warning.
3. **Margin census (report only).** For every decision site the harness records the gap between the compared argument and its threshold, in ulp. The number of decisions within a few hundred ulp of a threshold, against the observed relative size of libm differences, predicts the flip rate without a huge N. A random 1-ulp perturbation flips a decision with probability of about twice the relative perturbation times the distance scale, which is many orders below 1e-4 for independent decisions; so a measured rate near 1e-4 would point to structural ties (constant boundary conditions, integer temperatures) rather than noise, and the report lists them.
4. **Statistic self-test.** The bound code is checked against the table of section 6 and against a Monte-Carlo coverage test.
5. **Instrumentation neutrality and the snapshot closure** of section 4.5.

## 11. Open questions

| # | Question | Owner |
|---|---|---|
| V1 | Confidence level (95% proposed for each run, 99% at phase exit) and the three-way verdict (PASS, FAIL, INCONCLUSIVE). Is "INCONCLUSIVE does not sign" acceptable? | AMR V&V Lead |
| V2 | The denominator is exposed distinct-input calls only, and non-exposed calls are held bitwise (stronger than D-070 for the libm kernel as a whole). Does this fit condition 4 of test-plan 5.10, which says "per 10**4 record-calls"? | AMR V&V Lead |
| V3 | Replace the T1-class value 1e-10 by the two levels of section 7 (Level A at 2 ulp against the forced host; Level B with derived or measured class bounds). Sign `E_class` after the forced-flip survey. | AMR V&V Lead, with the Solid Phase Lead |
| V4 | Shortfall rule when P1 cannot reach the minimum N: does the gate rest on P2 with P1 reported? | AMR V&V Lead |
| A1 | How is "non-finite input only" (D-072) enforced in the kernel? Proposal: at the cap the kernel tests whether the state and `DT_BC` are all finite; if finite it continues exactly like upstream, if not it sets 301. The alternative (cap on count alone) can abort a legitimate `CHECK_FO=T` thin-layer case (section 8, reading 1). Does the cap also exist on the CPU path of the AMR build, so that both paths print the same text? | AMR Chief Architect |
| A2 | Text and abort policy for status 301 (same as D-075 for 300: the driver aborts with the report). Confirm the template of section 8 item 3. | AMR Chief Architect |
| A3 | Responses (a) and (b) of section 9 when the budget is exceeded: allowed without a new ruling? Delamination and burn-away flips are rate events with no bound: accept? | AMR Chief Architect |
| G1 | Test-build trace feature: stores of outcome and argument at the decision sites of section 3 into a fixed side array, stable site ids in the sidecar, neutrality of the traced build; a returned trace array on the device. | GPU Generator Engineer |
| G2 | Status-301 form of the cap with the finiteness predicate, and the `WRITE` plus `SHUTDOWN` rule of 06 section 8. | GPU Generator Engineer |
| G3 | The `ULP` line format of `check_device_logs.py` for the solve kernel (category `libm`), and device compile flags (`-gpu=nofma`). `MAX`, `MIN`, `MAXVAL` and `NINT` semantics for NaN and Inf in the device compiler: do they match gfortran? (Decides which poisons hang on the device.) | GPU Generator Engineer |
| S1 | Forced-flip survey on the host and the measurement of `E_class` (before the first device run); flip-injection build; the exposure tagger. | AMR Solid Phase Lead |
