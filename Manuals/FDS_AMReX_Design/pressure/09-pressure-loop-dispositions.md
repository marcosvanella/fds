# 09. Dispositions of the pressure loops claimed by the Pressure Solver Lead

Status: complete for the claimed loops. L1206 and L1208 are implemented as host code and tested; the other twelve loops are retired.

Scope: the loops of the Legacy Mapper's ranked list claimed by the AMR Pressure Solver Lead: L1206, L1208, L1210, L1212 to L1214,
L1216 to L1219, L1223 to L1227. Line numbers are those of the survey source (FireX 36975d7). `Source/pres.f90` is unchanged between
36975d7 and the head of the FDS-AMReX reference tree (the merge of FireX 1fb5c300fa): `git diff 36975d7 HEAD -- Source/pres.f90` is empty,
so every cited `pres.f90` line holds at both. AMR-path citations (`Source/driver`, `Source/pressure_backend`) are at the same reference-tree
head. Dispositions: RETIRED (no counterpart in the AMR path; what replaces it is named), HOST CODE, GENERATED KERNEL.

## 1. Summary

| Loop | pres.f90 | What | Disposition | Tests and results |
|---|---|---|---|---|
| L1206 | 729-742 (+744) | CHECK_POISSON residual of the separable equation, `POIS_ERR` | HOST CODE (diagnostic) | bitwise PASS (section 3) |
| L1208 | 768-787 (+792-793) | residual of the inseparable equation, `PRESSURE_ERROR_MAX/_LOC`, controls the pressure iteration | HOST CODE | bitwise PASS (section 3) |
| L1210 | 238-245 | PRHS, cylindrical K,I nest | RETIRED (cylindrical deferred in AMR mode) | none needed |
| L1212, L1213, L1214 | 267-277, 282-292, 297-307 | PRHS stored transposed (IPS 2, 3/6, 5) | RETIRED (FR-037; transposed IPS exists only with stretched grids, IR-002) | none needed |
| L1216, L1217, L1218 | 404-410, 414-420, 424-430 | HP from transposed PRHS (IPS 2, 3/6, 5) | RETIRED (same reasons as L1215 and L1212-L1214) | none needed |
| L1219 | 438-440 | tunnel preconditioner: add H_BAR to HP | RETIRED (charter Q4, D-041) | none needed |
| L1223 to L1227 | 531-537, 544-567, 576-582, 583-589, 591-596 | TUNNEL_POISSON_SOLVER nests | RETIRED (charter Q4, D-041) | none needed |

No loop is a GENERATED KERNEL at this stage. No shared file, no generator sidecar, no file of another owner and no FDS source was edited.

## 2. Evidence per loop

### L1206 (pres.f90:729-742): HOST CODE, diagnostic only

- What it computes: per cell `|RHSS - LHSS|` with `RHSS = -div F - dD/dt` (pres.f90:732-735) and `LHSS` the constant-coefficient
  Laplacian of HP (736-738), into the `WORK8` view (727); `POIS_ERR = MAXVAL(RESIDUAL)` (744). It runs only with `CHECK_POISSON`
  (pres.f90:726; default `.FALSE.`, cons.f90:203).
- Who reads it: `POIS_ERR` is printed to the `.out` file (dump.f90:4467), exchanged for output (main.f90:4421, 4455) and is the gas-phase
  output value of one quantity (dump.f90:8986). It never feeds back into the solution or the iteration control.
- AMR path: the single-level driver calls the unmodified routine on the level-0 meshes after every solve (`fds_p_resid`, fds_step.f90:464-468,
  called at TimeLoop.cpp:1302; `FDS_HOOK_L0_ONLY` keeps it on level 0). Its role as a check of the solve is taken by the backend's own
  residual report (`PressureResult::residual_rel2/relmax/residual_check` against `residual_limit`, PressureIface.H:243-289; enabled in
  the driver with `FDSTL_PBV`, TimeLoop.cpp:836), which checks `b - L H` for the problem actually solved; that is a different, normalised
  quantity and not bitwise `POIS_ERR`. No device kernel is needed (diagnostic, off by default).
- Implemented: `pres_residual_separable` (+ `residual_maxval_maxloc` for `POIS_ERR`) in `Source/pressure_lead/PresResidualLoops.{H,cpp}`,
  so that `POIS_ERR` can also be formed on boxes that the Fortran routine cannot reach (fine levels) when an owner wires it.

### L1208 (pres.f90:768-787): HOST CODE, controls the pressure iteration

- What it computes: per cell `|RHSS - LHSS|` of the inseparable equation, `RHSS` from `F - F_B` (771-774), `LHSS` from
  `div((1/rho) grad p) + lap K` with `P = RHOP*(HP-KRES)` (L1207, 758-764, the backend implementer's loop); then
  `PRESSURE_ERROR_MAX(NM) = MAXVAL`, `PRESSURE_ERROR_MAX_LOC(:,NM) = MAXLOC` (792-793), optional copy to `PP_RESIDUAL` (794).
  It runs while `ITERATE_BAROCLINIC_TERM` is true (750), which `PRESSURE_ITERATION_SCHEME` sets from `BAROCLINIC` (main.f90:1686-1689;
  `BAROCLINIC` defaults to `.TRUE.`, cons.f90:229).
- Not a diagnostic: `MAXVAL(PRESSURE_ERROR_MAX)` switches the baroclinic iteration off once below `PRESSURE_TOLERANCE` (main.f90:1800),
  is half of the exit test together with `VELOCITY_TOLERANCE` (1805-1807) and enters the suspension test (1812-1817); it is also written to
  the velocity error file (1789-1796) and the `.out` file (dump.f90:4444-4448). A wrong value changes the number of pressure iterations and
  therefore the solution.
- AMR path: the single-level driver reproduces the iteration (TimeLoop.cpp:1275-1322): it calls the unmodified Fortran routine on level 0
  (`fds_p_resid`, 1302), reads `PRESSURE_ERROR_MAX` through `fds_p_get_err` (fds_step.f90:477-483) and applies the same three tests
  (TimeLoop.cpp:1312-1319). So today L1208 is host code (the FDS routine itself) on level 0. The backend's residual report does not replace
  it: the backend solves the constant-coefficient problem and reports `b - L H` of that problem, not the variable-density residual.
- Gap for the driver owner: with refinement the composite branch does one solve per stage and does not compute the residual
  (`pressure_scheme_composite`, TimeLoop.cpp:1037-1060, `perr_out = 0`; TimeLoop.H:117 lists the pressure iteration as not built for
  several levels). `VELOCITY_TOLERANCE` and `PRESSURE_TOLERANCE` therefore have no effect in a refined run. When the multi-level iteration is
  built, the Fortran routine cannot run on fine boxes; `pres_residual_inseparable` + `residual_maxval_maxloc` are the box-indexed host code
  for it (the composite maximum over levels and boxes, and whether covered coarse cells count, are the driver owner's decision).
- Device: a GENERATED KERNEL is the later route for a full-GPU step (D-027); the loop list records what it needs (the `WORK8` view with
  `layout = exact`, MAXVAL/MAXLOC after the loop on the host). Not done here.

### L1210 (pres.f90:238-245): RETIRED

- The cylindrical branch of PRHS (`IF (CYLINDRICAL)`, 236). Cylindrical meshes are deferred in AMR mode (D-042; FR-004 rejects them with
  refinement) and the driver aborts on them in every mode it runs: `FdsAmr.cpp:20` ("CYLINDRICAL meshes are not supported (IR-002)") and the
  scope check `chk(IP[14], "CYLINDRICAL")` (TimeLoop.cpp:290). Nothing replaces it in the AMR path; the FDS-only build keeps it.
- If cylindrical meshes come back: host code (a K,I nest with J=1; the generator front end takes perfect K,J,I nests only, s5gen.py:1307,
  the same gap as L1392/L1393).

### L1212, L1213, L1214 (pres.f90:267-307): RETIRED

- They store PRHS transposed for IPS 2, 3/6 and 5. IPS is chosen from the stretched directions `NOC` (init.f90:2358-2364): the transposed
  layouts exist only for meshes stretched in y, z, x+z or y+z. Stretched grids are rejected by the AMReX code in every mode (IR-002, D-030),
  and the driver stops on any `IPS /= 0` (TimeLoop.cpp:320).
- FR-037: the backend solves on the natural box layout and the right-hand side reaches it as PRHS in the IPS=0 layout (`fds_p_get_prhs`,
  fds_step.f90:378; TimeLoop.cpp:842). The natural-layout RHS is L1211 (backend implementer, `pres_compute_rhs_div`). There is no
  transposed store anywhere in the AMR path.

### L1216, L1217, L1218 (pres.f90:404-430): RETIRED

- Copies of the transposed FFT solution back into HP (IPS 2, 3/6, 5). Same reasoning as the backend implementer's retirement of L1215
  (pres.f90:394-400): `pressure_backend` writes phi (H) directly (FR-037) and the driver copies phi to H on the valid cells with a plain loop
  (`solve_level0`, TimeLoop.cpp:855-862; composite branch TimeLoop.cpp:1019-1025). In addition the transposed IPS cannot occur (L1212 above).

### L1219 (pres.f90:438-440) and L1223 to L1227 (pres.f90:531-596): RETIRED

- L1219 adds the 1-D tunnel solution `H_BAR` to HP inside `IF (TUNNEL_PRECONDITIONER)` (436-445). L1223 to L1227 are nests of
  `TUNNEL_POISSON_SOLVER` (505-697): the ULMAT solid zeroing of PRHS (L1223, 530-538), the cross-section sums of PRHS into `TP_CC` (L1224,
  544-567), and the boundary averages `BXS_BAR`/`BXF_BAR` (L1225/L1226 for ULMAT, 576-589; L1227, 591-596). The routine is called only from
  `PRESSURE_ITERATION_SCHEME` when `TUNNEL_PRECONDITIONER` is set (main.f90:1733, 1742).
- Charter Q4 (D-041): `TUNNEL_PRECONDITIONER` is out in AMR mode. The driver never calls the routine (no call site in `Source/driver`);
  its scope check refuses the keyword (TimeLoop.cpp:290), and the backend refuses it too (`HFillOptions::tunnel_preconditioner` gives
  `NotBuilt`, FdsPressureLoops.H). The tunnel cases are covered by MLMG (README Q7 of this directory: obsolete under MLMG; `tunnel_demo`
  stays as a regression and timing case).
- Note for the specification owner: FR-004 and D-041 say the keyword is accepted with a warning and ignored in AMR mode, while D-080 (6) C2
  and the driver refuse it. Either way no tunnel loop runs in the AMR path; the wording should be made consistent.

## 3. Implementation and tests (L1206, L1208)

Code: `Source/pressure_lead/` on the private branch `pressure-lead-loops` of the FDS-AMReX repository (own new directory; README there).
`PresResidualLoops.{H,cpp}` (namespace `pressure_lead::resid`): `pres_residual_separable` (L1206), `pres_residual_inseparable` (L1208),
`residual_maxval_maxloc` (MAXVAL/MAXLOC, pres.f90:744 and 792-793). Fortran-ordered views with explicit lower bounds (an FDS mesh array or
one AMReX box), the cylindrical factors R and RRN kept, the Fortran evaluation order, `-ffp-contract=off`, OpenMP over cells. MAXLOC follows
gfortran: first occurrence in array element order, NaN skipped, all-NaN gives NaN at the first cell, an empty box gives -HUGE at 0.

Tests: `Source/pressure_lead/tests/run_residual_tests.sh` (private build directory; no shared lock needed, about one minute).

| Check | Result |
|---|---|
| Reference | the upstream lines pres.f90:727-744 and 752-793 copied verbatim (pinned by first/last statement and SHA-256; the script stops on drift), including their OpenMP directives, so the OpenMP builds run the reference threaded. `test/make_r2_tests.py` builds references for generated kernels from the generator's argument manifest and does not apply to host code, so the same line-range method is in an own script (`residual_gen.py`). |
| Synthetic cases | 14 cases: sizes 1x1x1 to 9x7x5 and JBAR=1, non-unit R/RRN and non-uniform metrics, Cartesian metrics, ties in every cell, a tie of two cells (first wins), one NaN, all NaN, magnitudes 1e-300 to 1e150. Each case also as four sub-boxes and with shifted index origins (as an AMReX box gives); a view one ghost layer short must be refused. 8238 compared values per run. |
| Flag sets, callee switches, threads | O0, O2, O0omp, O2omp, O2omp_off, O2omp_dpd x dpd, bind x 1, 4, 8 threads (OpenMP sets): 28 runs, 0 differing values; the reference files of every build are byte-identical to the first build. The callee switches and the offload macros cannot change host code (no generated region, no callee directive); they are passed so the matrix is the generator's. |
| Real FDS | a scratch FDS (FireX 36975d7 GNU Release objects, only a scratch copy of pres.f90 with a write-only dump hook recompiled with the same flags and relinked, in a private directory) on two inputs adapted from the backend implementer's dump cases with `CHECK_POISSON=T`: one mesh with fire and an inflow ramp; two MPI ranks with wind, a pressure ramp and `VELOCITY_TOLERANCE=1e-4`, `PRESSURE_TOLERANCE=1e-3`, up to 8 pressure iterations per stage. Steps 1 to 3, predictor and corrector, every pressure iteration: 214 files (118 L1206, 96 L1208; 612,406 values incl. `POIS_ERR`, `PRESSURE_ERROR_MAX`, `PRESSURE_ERROR_MAX_LOC`), 0 differences at O0, O2 and O2omp with 1 and 8 threads. In the two-mesh run `PRESSURE_ERROR_MAX` spans 2e-5 to 40 and falls below the tolerance in 25 of 92 calls, so both sides of the iteration switch were exercised. 8 of the files are committed (`tests/fds_resdump_cases/fds_resdump_cases.tar`) and checked in every run of the matrix. The hook is a test aid on a scratch copy only (not an upstream patch). |
| Mutation check | 14 mutants (one index, operand or sign in each loop, the sum order, the R and RRN factors, the 2.0 constant, last-wins ties, NaN taken as maximum): 14 of 14 caught by the synthetic cases. On the real-FDS files the R/RRN mutants and the tie/NaN mutants are not caught (Cartesian runs have R = RRN = 1, and real residuals have no exact ties or NaN); that is why the synthetic cases carry them. |

Commits (FDS-AMReX repository, private branch `pressure-lead-loops`, not pushed): `2a31d4562d` (code and tests), `8bb1b207d5` (dump
archive as a plain tar, because `*.gz` is ignored by the tree).

## 4. For the register and other owners

- Legacy Mapper (claim register): L1206 and L1208 to `not planned (host)` with the notes of section 2; L1210, L1212 to L1214, L1216 to L1219
  and L1223 to L1227 to `not planned (retired)` with the replacement named in section 2. The register and the generated tables were not
  edited here.
- Driver owner: the multi-level pressure scheme computes no inseparable residual and does not iterate (section 2, L1208); the host functions
  are ready when that is built. Wiring them (or `pres_residual_separable` for `POIS_ERR` on fine levels) is a driver change and was not made.
- Generator owners: nothing requested now. A device kernel for L1208 would need the `WORK8` view (`layout = exact`) and a host-side
  MAXVAL/MAXLOC, as already recorded in the loop list.
- Specification owner: the `TUNNEL_PRECONDITIONER` wording (accept-and-ignore in FR-004/D-041 versus refuse in D-080 (6) C2 and the driver).
