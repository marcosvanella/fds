# 06 · Test design for porting the thick-wall 1-D solve (`SOLID_HEAT_TRANSFER`) to a generator kernel

Owner: AMR Solid Phase Lead · Status: **v0.1, test design for review** · Critical path of `docs/amrex/wall-bc-translation-plan.md` (section 8.3, 24 to 36 working days).
Reference revision: FireX `36975d765f`; every `file:line` is from that revision (the working tree is 9 to 12 lines higher, section 0 of the plan). `wall.f90` unless a file is named. Nothing in `src` or in the generator was changed, and nothing was run to write this page (the shared generator lock is busy): all statements about upstream are from reading the pinned text, and all statements about the existing test style are from reading `wall_sums.py`, `test_wall_sums.py`, `wall_checks.py` and the Wall engineer's `test_wall_*.py`.
Inputs: `wall-bc-translation-plan.md` (sections 3, 4, 8, 10, 18), `05-fine-level-solid-plan.md`, `blocked-loop-families.md` (SP4), `generator-howto.md` (section 4: acceptance rule), rulings D-064, D-065 and D-070.

## 1. What is being ported, and what "done" means

`SOLID_HEAT_TRANSFER(NM,T,DT_BC,PARTICLE_INDEX,WALL_INDEX,CFACE_INDEX,THIN_WALL_INDEX)` is `wall.f90:1809-3154` (1346 lines including its contained `PERFORM_PYROLYSIS`, 2993-3152; the plan's "1186 lines" is the body without the pyrolysis and the header). It advances the 1-D conduction record `ONE_D` of one wall cell over `DT_BC` with an adaptive sub-step, includes pyrolysis (`PYROLYSIS`, 3195-3605, called from `PERFORM_PYROLYSIS`), layer shrinking, renoding and removal, and writes the wall's surface state `B1`/`B2`. The device call sites are pass H2 (thick walls, L1488 line 156) and L1489 (thin walls, line 193) of the plan.

**Scope of this design.** `WALL_INDEX` and `THIN_WALL_INDEX` entries, `HT_DIM=1` (asserted; plan section 5.1). `CFACE_INDEX` and `PARTICLE_INDEX` entries stay on the host (deferred: GEOM, particles). The gas-side pass (near-surface variables, HTC of the front, `ZZ_F`) is not part of the solve; the solve reads its results through `B1`.

**Acceptance (claim protocol, `generator-howto.md` section 4, applied unchanged).** A solve unit is translated when its kernel is in the committed sidecar and is **bitwise equal** to the verbatim upstream text, at 1, 4 and 8 threads, for all six flag sets (`O0`, `O2`, `O0omp`, `O2omp`, `O2omp_off`, `O2omp_dpd`), both callee switches (`dpd`, `bind`), with every mutant caught and every negative check refused. Records are independent of each other in the solve (no shared write, plan section 5.2), so bitwise equality with a serial run is reachable without any ordering device. The tolerance cases are limited to: libm calls on a real device (D-070, section 10), and owner aggregates (section 11).

## 2. Anatomy of the routine (the decomposition input)

Survey lines. "State" is the `ONE_D` record unless stated.

| Block | Lines | What it does | Notes for the port |
|---|---|---|---|
| H | 1809-1850 | header, dummies, ~50 scalars, 25 arrays of extent `NWP_MAX` or `0:NWP_MAX+1` (1832-1836), `ALLOCATABLE` locals `RHO_DOT`, `INT_WGT` (1827) | private arrays of module extent `NWP_MAX`: ~25 x `NWP_MAX` doubles per thread |
| B1 | 1854-1990 | alias set: front (`WALL_INDEX`), thin wall (1901-1932, back-side wall through `MESHES(BACK_MESH)%THIN_WALL`), CFACE/particle (1933-1985, host only); `BACKING` selection, `ISOLATED_THIN_WALL`, `FLAT_PLATE` | `MESHES(...)` lookups: accessor of D-064 (3); `OPTIONAL`/`PRESENT` dummies |
| B2 | 1987-2043 | `LAYER_REMOVED=F`; return if burned away (1991-1999); adjusted densities (2000-2006); zero the flux accumulators (2008-2027); `NWP`, `RHO_S` (2035-2042); `ALLOCATE(RHO_DOT)` (2016-2017) | allocation becomes a private array of extent `(N_MATL, NWP_MAX)` |
| L | 2044-2953 | `SUB_TIMESTEP_LOOP`: no iteration cap, exit at 2947 | the adaptive part |
| L1 | 2046-2064 | grid weights (`GET_WALL_NODE_WEIGHTS`), geometry arrays | callee, func.f90:5860 |
| L2 | 2065-2100 | properties for the Fo check (`CHECK_FO`) | duplicate of L9 |
| L3 | 2101-2190 | front convective flux; back side by `SELECT CASE(BACKING)`: `VOID` (HTC of the back with `BACK_SIDE=.TRUE.`, 2146), `INSULATED`, `EXPOSED` (reads `B1_BACK`, 2163-2190) | back-wall HTC snapshot (D-065 Q2) |
| L4 | 2192-2204 | thickness, radius arrays `R_S` | `**` with real exponent `I_GRAD` |
| L5 | 2206-2223 | pyrolysis: `PERFORM_PYROLYSIS` (predicted) or specified burn rate | contained routine |
| L6 | 2224-2295 | heat sources: user `HEAT_SOURCE` with ramps (2229-2233), boundary fuel model convection (2235-2256, calls HTC), internal radiation two-sweep recursion (2261-2292) | `EVALUATE_RAMP`, HTC callee |
| L7 | 2296-2337 | sub-step control: explicit estimate of `DELTA_TMP`, `TMP_RATIO`, `DT_BC_SUB = DT_BC/MIN(NINT(TIME_STEP_FACTOR*WALL_INCREMENT), MAX(1,NINT(TMP_RATIO)))`, `MIN` with remaining time and `DT_FO`, rebuild `Q_S` if `DT_BC_SUB` changed | **guarded by `ICYC>WALL_INCREMENT`** (module global); `NINT` is a discontinuity |
| L8 | 2339-2357 | store the sub-step's fluxes (`Q_RAD_OUT`, pyrolysis nets) | |
| L9 | 2359-2805 | layer masses and thicknesses (`POINT_LOOP2`, 2373-2430), new coordinates (2432-2475), zero-thickness nodes (2477-2500), layer thickness check (2502-2540), remesh decision (2585-2670), **abort `ERROR(300)`** (2668-2670), zero-thickness wall: `BURNAWAY`, `OBSTRUCTION%MASS=-1.` (2680-2698, shared event, plan row A), regrid with `ALLOCATE(INT_WGT)` and five `INTERPOLATE_WALL_ARRAY` calls (2724-2754) | the largest block (445 lines) |
| L10 | 2815-2848 | properties from material tables `ML%K_S(ITMP)`, `ML%C_S(ITMP)` with `ITMP=MIN(I_MAX_TEMP-1,INT(TMP))`; `K_S` averaging to cell faces; emissivity | temperature dependence is a table lookup by integer temperature |
| L11 | 2850-2908 | assemble the implicit system (`AAS,BBS,CCS,DDS,DDT`), boundary coefficients `RFACF,QDXKF,RFACB,QDXKB`, tridiagonal solve (2900-2908) | sequential recurrence |
| L12 | 2909-2951 | clip to `[TMPMIN,TMPMAX]`, ghost values, `Q_CON_F`/`Q_RAD_OUT` increments, particle production, exit test, `N_SUBSTEPS` | |
| P | 2953-2990 | scale `M_DOT_LAYER_PP`, copy nets to `B1`, `T_IGN`, `INT_FTP`, extinction | |
| C | 2993-3152 | `PERFORM_PYROLYSIS` (contained): per node `PYROLYSIS`, optional bounded Newton on `Y_O2_F` (`O2_LOOP`, `MAX_ITER=20`) | silent fall-through when not converged (3150) |
| Y | 3195-3605 | `PYROLYSIS`: Arrhenius rates (`EXP`, `**`), evaporation (`LOG`, `SQRT`), char/oxidation, particle fluxes | 411 lines, many callees (viscosity, film properties, Rayleigh models) |

Callees by file: `func.f90`: `GET_N_LAYER_CELLS` 5743-5784, `GET_WALL_NODE_COORDINATES` 5800-5838, `GET_WALL_NODE_WEIGHTS` 5860-5994, `GET_INTERPOLATION_WEIGHTS` 6007-6039, `INTERPOLATE_WALL_ARRAY` 6050-6067, `GET_EMISSIVITY` 1979-1999, `EVALUATE_RAMP` 802-857, `HEAT_TRANSFER_COEFFICIENT` 3135-3321, `GET_MASS_FRACTION` 1579-1588, `GET_SPECIFIC_HEAT` 1729-1739.

Three facts about the control flow that shape the tests:
1. **One abort path.** The only `SHUTDOWN` in 1809-3154 is `ERROR(300)` at 2670 (`NWP_NEW > N_CELLS_MAX`). `PYROLYSIS` has none.
2. **Two ways to "not converge".** The sub-step loop has no cap (with a NaN `DT_BC` the sub-step time becomes NaN and the exit test at 2947 is false for ever; a NaN temperature alone ends after one sub-step with NaN outputs, see `09-flip-budget-gate.md`; the NaN cases of section 8 take their expected outcome from the reference run before the watchdog is relied on). The Newton loop `O2_LOOP` simply ends after 20 iterations with the last `Y_O2_F` and without writing `B2%Y_O2_ITER` (3119-3150).
3. **Hidden state.** `ICYC`, `WALL_INCREMENT`, `CHECK_FO`, `OXPYRO_MODEL`, `TMPMIN/TMPMAX`, `I_MAX_TEMP`, `NWP_MAX` are module data, not arguments.

## 3. Testable units in dependency order

Each unit has a verbatim reference range, the inputs and outputs the harness sets and compares, and its equality class: **B** bitwise on the host (all six flag sets), **B/L** bitwise on the host and the libm tolerance on a real device (D-070: 2 ulp per value, applies where `**` with a non-integer exponent, `EXP`, `LOG` or a variable integer power occur), **D** discrete decision (compared exactly; on a device a libm difference can flip it, section 10).

| # | Unit | Reference range | Inputs | Outputs compared | Class |
|---|---|---|---|---|---|
| U0 | Record store: builder, offsets, capacity, canaries (Python, no upstream text) | `type.f90:217-264`, `init.f90:2975-3386` for what a record contains | `SURF`/`MATL` tables, wall list | `OD_OFF`, `OD_N`, `OD_CAP`, MATL offsets, round trip record to `ONE_D` and back | exact integers; reals bitwise on round trip |
| U1 | Grid callees: layer cell counts, node coordinates, node weights, interpolation weights, array interpolation | func.f90:5743-5784, 5800-5838, 5860-5994, 6007-6039, 6050-6067 | layer thicknesses, stretch, cell counts, geometry (cartesian, cylindrical, inner cylindrical, spherical) | node arrays, weights, `INT_WGT`, interpolated arrays | B/L (`STRETCH_FACTOR**(N-1)` is a variable integer power) |
| U2 | Emissivity | func.f90:1979-1999 | `MATL_COMP` densities, `MATL_INDEX`, node index | emissivity | B |
| U3 | Property update (L10) and the Fo copy (L2) | wall.f90:2815-2848, 2065-2100 | `TMP`, `MATL_COMP(N)%RHO`, `RHO_ADJUSTED`, material tables, `PACKING_RATIO`, `BOUNDARY_FUEL_MODEL` | `K_S`, `RHO_C_S`, `RHO_S`, `B1%EMISSIVITY`, `DT_FO` | B (sums in material order; `MINVAL` exact) |
| U4 | Ramp function (shared with the Wall engineer, plan section 4) | func.f90:802-857 | flat ramp tables | value | B; owned elsewhere, only its use inside the solve is tested here |
| U5 | HTC callee for the back side and the boundary fuel model | func.f90:3135-3321 | `SF`, `DTMP`, `BACK_SIDE`, `WALL_INDEX_IN` | HTC, `B2` diagnostics | B/L; owned by the Wall engineer (L1486); the solve tests use it through the accessor of D-064 (3) |
| U6 | Set-up: alias selection, `BACKING`, `ISOLATED_THIN_WALL`, `FLAT_PLATE`, burned-away early return, adjusted densities, accumulator reset | 1854-2043 | record, `B1`, `B2`, `SF`, entry kind | locals and the early-return outcome | B (copies and integer logic) |
| U7 | Sub-step body without pyrolysis, thick, cartesian: grid (L1), front flux and `VOID`/`INSULATED` back (L3), thickness and radius (L4), implicit update (L10 to L12) with `DT_BC_SUB` forced to `DT_BC` | 2044-2064, 2101-2204, 2815-2951 | one layer, one material, constant and table properties | `TMP(0:NWP+1)`, `DELTA_TMP`, `B1%TMP_F`, `TMP_B`, `Q_CON_F`, `Q_RAD_OUT`, `TMP_F_OLD` | B |
| U8 | Sub-step control (L7) with the explicit estimate | 2296-2337 | as U7 plus `DELTA_TMP_MAX`, `TIME_STEP_FACTOR`, `ICYC`, `WALL_INCREMENT`, `CHECK_FO` | `DT_BC_SUB` sequence, number of sub-steps, `N_SUBSTEPS`, all U7 outputs | D for the count; B for the rest |
| U9 | Extended boundary conditions: `EXPOSED` back, `DIRICHLET_BACK`, `TMP_GAS_BACK` ramp, back emissivity variants | 2121-2190, 2870-2896 | back record snapshot (D-065 Q2), ramps | `Q_CON_B`, `Q_RAD_IN_B`, `TMP_B`, back coefficients | B |
| U10 | Heat sources: user `HEAT_SOURCE` with ramp, boundary fuel model convection, internal radiation sweeps | 2224-2295 | `HEAT_SOURCE`, `RAMP_IHS_INDEX`, `INTERNAL_RADIATION`, `KAPPA_S` | `Q_S`, `Q_IR`, `Q_ADD`, `Q_RAD_OUT` | B/L (the fuel-model path calls HTC) |
| U11 | Multi-layer, multi-material, cylindrical and spherical geometry | 2002-2006, 2192-2204, 2815-2848 | 1 to 4 layers, 1 to 3 materials per layer | as U7 | B/L (`R_S**I_GRAD`) |
| U12 | Layer mass and thickness update (`POINT_LOOP2`), new coordinates, zero-thickness nodes | 2359-2500 | `RHO_DOT`, densities, `Q_S` | `MATL_COMP%RHO`, `LAYER_THICKNESS`, `X_S_NEW`, `Q_S` redistribution | B/L |
| U13 | Remesh decision and regrid: layer checks, `GET_N_LAYER_CELLS`, `INT_WGT`, five interpolations, layer removal | 2502-2805 | layers near the thresholds | `NWP`, `N_LAYER_CELLS`, `X`, interpolated arrays, `LAYER_REMOVED` | D for the decision; B/L for values |
| U14 | Pyrolysis, reaction rates and fluxes | 3195-3605 | `RHO_S`, `TMP_S`, material reaction tables, gas-side inputs | `RHO_DOT_OUT`, `M_DOT_G_PPP_*`, `M_DOT_S_PPP`, `Q_DOT_*_PPP`, particle fluxes | B/L |
| U15 | `PERFORM_PYROLYSIS`: point loop, geometry factor sums, bounded Newton loop | 2993-3152 | as U14 plus `OXPYRO_MODEL` | `B1`/`B2` fluxes, `Y_O2_F`, `Y_O2_ITER`, final iteration state | B/L, D for the iteration count |
| U16 | Post-loop (P) | 2953-2990 | accumulators, `T`, `DT_BC` | `M_DOT_LAYER_PP`, `M_DOT_G_PP_*`, `Q_DOT_G_PP`, `T_IGN`, `INT_FTP`, extinction | B |
| U17 | Whole routine, thick wall | 1809-3154 (wall entry) | the full record and all inputs | every field the routine writes (section 5.4) | B (device: B/L) |
| U18 | Thin-wall entry | 1901-1932 and the thin-specific branches (2106, 2895) | thin-wall record, two side walls, back thin wall | as U17 plus the burn-away event | B |
| U19 | Error and watchdog paths | 2668-2670; sub-step cap (section 8) | capacity set too small; NaN temperature | status codes, message text, no write outside the record | exact |

**Dependency order and why.** U0 and U1 to U4 have no dependency on the solve and start first. U5 is another role's deliverable and is mocked until it lands (a bitwise-equal host copy of the verbatim function is used for the reference *and* for the kernel in the mock stage, so the solve tests are not blocked by it). U6 to U9 and U16 are the thin vertical slice (the "core"): one layer, no pyrolysis. U10 and U11 widen the core. U12, U13 (renoding) and U14, U15 (pyrolysis) are independent of each other and of U7 to U11 as units (their inputs are explicit arrays), so they run in a second track and meet the core at U17. U18 and U19 come last.

**Block kernels versus the whole kernel.** The generator translates the whole routine. For bring-up, units U7 to U16 need to run as separate small kernels so a failure points at a block. This assumes the generator can build a kernel from a marked sub-range of the routine with the surrounding state passed as arguments, for test builds only (question T1). If it cannot, the fallback is the verbatim-cut harness of `test_wall_sums.py` (`upstream_blocks`, `probe_source`): the reference side is cut by anchor text; the kernel side is the whole kernel with the other blocks switched off by input flags that make them no-ops (for example no pyrolysis model, remesh disabled). The acceptance gate is always U17/U18 on the whole kernel.

## 4. Ragged record layout: options and what the tests check

The record holds arrays of three kinds: per node `(0:NWP+1)` (`TMP`, `X`, `DELTA_TMP`, `RHO_C_S`, `K_S`, `DX_OLD`, `HEAT_SOURCE`), per layer `(1:N_LAYERS)` (`LAYER_THICKNESS`, `N_LAYER_CELLS`, `STRETCH_FACTOR`, `DDSUM`, ...), per material and node `MATL_COMP(N)%RHO(1:NWP)`, plus scalars. The capacity of a record is `N_CELLS_MAX` of its surface (`type.f90:220`, known at set-up); `NWP` changes during the run inside that capacity (renoding, layer removal), so no reallocation is needed.

| Option | Layout | For | Against |
|---|---|---|---|
| **R1, capacity-offset prefix sums (recommended)** | `OD_OFF(IR)` = start of record IR in each flat array, `OD_CAP(IR)` = capacity, `OD_N(IR)` = current `NWP`; node arrays of record IR occupy `OD_OFF(IR)+0 : OD_OFF(IR)+OD_CAP(IR)+1`; offsets built by an exclusive prefix sum over capacities (the CSR pattern of `wall_checks.csr_cell_walls`); per-material arrays use a second prefix `OM_OFF(IR)` with stride `N_MATL(IR)*(CAP+2)` and element `(N,I)` at `OM_OFF+(N-1)*(CAP+2)+I` | same builder, same checks and same `PTR/LIST` naming as the other wall tables; fixed offsets so renoding never moves data; coalescing is poor but acceptable (one thread per record, a record is a few KB) | memory for unused capacity (2.9 to 8.6 KB per record, `03-g2a-cost-check.md`) |
| R2, uniform stride by surface class | one fixed row length per `SURF` class, `row = class_offset + k*class_stride` | simplest index arithmetic | breaks if one class is large; a second table anyway for the class |
| R3, packed by current `NWP` | prefix sum over current `NWP` | smallest memory | offsets move when `NWP` changes: a record-store rebuild on every renoding step; rejected |

The fine-level plan (`05` section 2.1) uses the same store with one row per record and the face key as sort key; R1 satisfies it with `IR` = record index. Private versus scratch for the solve's local arrays: the 25 locals of extent `NWP_MAX` either become per-thread private arrays (device local memory) or rows of a scratch store with the same R1 offsets (`SC_OFF`); both give identical bits and the tests run both (switch `SCRATCH=private|rows`, a generator policy), because the footprint decides which one a GPU can afford.

**What the layout tests check (all in Python, run at every build, plus the Fortran round trip):**
1. Offsets are the exclusive prefix sum of capacities (any wall order); no two records overlap; the last offset plus capacity equals the allocated length.
2. `0 <= OD_N(IR) <= OD_CAP(IR)` for every record; `OD_CAP` equals the surface's `N_CELLS_MAX`; materials: `OM_OFF` strides equal `N_MATL*(CAP+2)`.
3. Round trip: random `ONE_D` structures (all field kinds, random `NWP`) pack to the flat store and unpack to bitwise-identical structures; elements beyond `OD_N` keep their prior bits (untouched padding).
4. **Canaries:** every record row is padded with a sentinel band of known bit pattern; after any kernel call the bands are intact in every flag set (this is the check that an out-of-capacity store, for example from a wrong `NWP_NEW`, is detected even when no abort happens).
5. Order independence: permuting the record list permutes results and nothing else; records visited in any order and any thread count give the same bits.
6. Rebuild after a refresh gives the same tables as a build from scratch (plan section 11.3), and an incremental append (new obstruction) leaves existing offsets unchanged.
7. Mutants of the builder (off-by-one capacity, `NWP` instead of `NWP+2`, offsets from a non-exclusive prefix sum, material stride from the wrong `N_MATL`) must be caught by items 1 to 4 and by the kernel comparison.

## 5. Reference method: verbatim serial transcription

1. **Text.** The reference is the upstream source text of the range, cut by line range and anchor from `git show 36975d765f:Source/wall.f90` (never retyped), `!$OMP` lines removed, run serially. The same rule as `make_r2_tests.py`, `test_wall_sums.py` (`upstream_blocks`, `probe_source`) and the Wall engineer's `upstream_loop`. An anchor check (first statement normalised) fails the build when the pinned text moved.
2. **Mock module state.** The harness compiles, in a module that stands for the FDS modules: the derived types `BOUNDARY_ONE_D_TYPE`, `BOUNDARY_PROP1_TYPE`, `BOUNDARY_PROP2_TYPE`, `MATL_COMP_TYPE`, `SURFACE_TYPE` (cut verbatim from `type.f90`, member lists reduced only by deleting members the routine never touches, checked by a list of members read and written), `MATERIAL_TYPE` with the temperature tables `K_S(0:I_MAX_TEMP)`, `C_S(0:I_MAX_TEMP)`, the module constants and switches of section 2 fact 3, and the callees compiled verbatim from `func.f90` (U1, U2, U4, U5 in the mock stage). `SHUTDOWN` is a mock that writes the message to a file and stops with a fixed exit code (section 8). `MESHES(...)` lookups are replaced by pointers set by the harness (the accessor of D-064 (3) in the kernel version).
3. **Same inputs.** Case files (Python writes them, Fortran reads them: the `mesh.txt` pattern of the Wall engineer's drivers) set one `ONE_D` record, `B1`, `B2`, `SF` and the globals. The reference routine runs on a copy; the kernel runs on the flat store (R1) from the same case file.
4. **Output set.** Every field the routine writes, found by a read/write table of the pinned text (a Python scan of the routine for assignments to `ONE_D%`, `B1%`, `B2%` and `OBSTRUCTION(...)%MASS`; the solve has no gas-cell write), checked in both directions: a field written by the reference and absent from the compare list is a test failure (list-completeness check), and a field in the list never written by any case is reported (coverage).
5. **Comparison.** Bit patterns (`transfer(x,[0_c_int64_t])`), as the existing drivers do, so a signed zero or a NaN payload counts. The existing drivers also count a case as **vacuous** if the reference left the record unchanged; the solve harness keeps that rule and adds a per-case coverage tag set (section 6).
6. **Build.** gfortran, `-ffp-contract=off`, no `-ffast-math`, `-fno-range-check`, the six flag sets and both callee switches. The reference is compiled with the same flags as the kernel in the same cell, so a host compiler difference never enters.
7. **Capture-and-replay variant** (integration, section 11): the "case file" is a dump of the pre-call state of real records, written by a patched host run, and the reference is the unmodified routine on the dumped state.

## 6. Case generators

A Python generator (`wall_solve_cases.py`, style of `make_meshes` and `make_case`) produces random records with a seed and a list of **coverage tags**. Axes:

| Axis | Values |
|---|---|
| Nodes per layer | 1 (degenerate), 2, 3, 5, 10, 30, up to `N_CELLS_MAX`; total `NWP` 2 to 60; `NWP` at the capacity limit |
| Layers and materials | 1 to 4 layers; 1 to 3 materials per layer (`MATL_COMP`); densities with zero entries (skipped by `<=TWENTY_EPSILON_EB`) |
| Property tables | `K_S(T)` and `C_S(T)` random smooth tables over 0 to 5000 K, constant tables, tables with a step; temperatures at integer boundaries, just below, just above, at `I_MAX_TEMP-1` (clamp) and at the clip limits `TMPMIN`, `TMPMAX` |
| Geometry | cartesian, cylindrical, inner cylindrical, spherical; `I_GRAD` 1, 2, 3; `INNER_RADIUS` 0 and positive |
| Front BC | convective `HEAT_TRANS_COEF` from 0 to 1e4, `TMP_G` above and below `TMP_F`, `Q_RAD_IN` 0 to 1e5, `Q_CONDENSE`, `LP_CPUA` |
| Radiation | `RADIATION` on/off, `EMISSIVITY` specified or computed, `INTERNAL_RADIATION` on/off with `KAPPA_S` ranges, thin-optical and thick-optical |
| Back BC | `VOID` (with and without `DIRICHLET_BACK`, `TMP_GAS_BACK` ramp, specified back emissivity), `INSULATED`, `EXPOSED` (back record from a snapshot) |
| Heat source | none, constant, ramped `HEAT_SOURCE`, boundary fuel model (`PACKING_RATIO`, `BOUNDARY_FUEL_MODEL`) |
| Time | `DT_BC` from 1e-4 to 10 s (`WALL_INCREMENT` 1, 2, 3), `ICYC` below, at and above `WALL_INCREMENT`; `TIME_STEP_FACTOR`; `CHECK_FO` on/off |
| Sub-step triggering | `DELTA_TMP_MAX` small and large so that `TMP_RATIO` takes values below 1, between 1 and 2, 5, 40, and the cap `NINT(TIME_STEP_FACTOR*WALL_INCREMENT)` is hit; flux steps that change `DT_BC_SUB` between sub-steps so the `Q_S` rebuild (2326) runs; an `DT_FO` limit |
| Pyrolysis | none; specified burn rate; predicted: 1 to 3 materials, 1 to 3 reactions each, `A` and `E` spanning inert to fast, `N_S`, `N_T`, `N_O2`, `HEAT_OF_REACTION`, `NU_GAS` to 1 to 3 species, char formation, evaporation (liquid, boiling, `B_NUMBER`), `OXPYRO_MODEL` with `MAX_ITER=20` converging and not converging |
| Layer events | layer shrinking (consumption), layer below `MIN_LAYER_THICKNESS`, below `MIN_LAYER_MASS`, remesh trips (`REMESH_NWP`), nodes going to zero, wall thickness to zero (burn-away), consumable and non-consumable obstruction |
| Entry kind | `WALL_INDEX` thick wall; `THIN_WALL_INDEX` with isolated, one-sided and two-sided back; thin wall with the back thin wall in the same or in another mesh object (accessor) |
| Steady and transient | transient from a cold or hot start; **steady** cases: constant BCs run for 50 to 500 calls until `TMP` stops changing (checks that a long accumulation of rounding steps stays bitwise: the comparison is made at several call counts, not only at the end) |

Each case is a sequence of 1 to N calls on the same record (state carried), because the record evolves and a one-call test cannot reach renoding or burn-away. The harness reports, per tag, how many cases reached it; **a coverage gate** fails the run if a tag has fewer cases than its minimum (default 8; `ERROR(300)` and the burn-away tags have their own fixed minimum 1 and run as dedicated cases). The gate matters because a mutant in a branch no case reaches would otherwise count as "caught" only by luck.

Case counts: quick profile 60 records x 3 calls (one flag set, a minute-scale run on one core); full profile 600 x 5 per flag-set cell. Seeds are fixed in the file and printed.

## 7. Flag sets, callee switches, threads

Per the claim protocol, for every unit U6 to U19 (and for U1 to U3 as far as they are kernels):

- Flag sets: `O0`, `O2`, `O0omp`, `O2omp`, `O2omp_off` (host fallback of the offload branch), `O2omp_dpd` (forced `distribute parallel do`).
- Callee switches: `dpd` and `bind` (`CALLEES`). The solve calls the table callees (grid, interpolation, emissivity, HTC, pyrolysis), so both switches change the generated call form and both must be bitwise.
- Threads: 1 for the non-OpenMP sets; 4 and 8 for the OpenMP sets, with the record loop running over the records: each thread owns whole records. Results must be equal across thread counts and also when the record list is shuffled.
- Compile flags `-ffp-contract=off`. A negative control builds one cell with FMA allowed and requires at least one mismatch on the full profile; if none appears the case set is too weak to see the difference (this keeps the "nofma" rule honest).
- Cost: the full matrix is 6 x 2 builds with 4 or 8 threads. It runs under `flock (local project directory)/src-s5gen/.s5gen.lock`, at most 4 cores, as everything else on this machine. The quick profile (one flag set) is the inner loop; the full matrix is run at the end of each stage, not per edit.

## 8. Host-only error path and the watchdog

**What the host does today.** `ERROR(300)` (2668-2670): `WRITE(MESSAGE,'(A,I5,A,A)') 'ERROR(300): N_LAYER_CELLS_MAX should be at least ',NWP_NEW,' for ',TRIM(SF%ID)` then `SHUTDOWN(MESSAGE,PROCESS_0_ONLY=.FALSE.)`. It is raised after the per-layer cell counts are computed and before `NWP_NEW` is used as an index (2724 on).

**Device form (proposal, needs the Generator Engineer's agreement, question T2).**
- Per record status `OD_STATUS(IR)` (integer, 0 = ok) and `OD_AUX(IR)` (integer, here `NWP_NEW`). At the `IF (NWP_NEW > ONE_D%N_CELLS_MAX)` site the kernel sets `OD_STATUS=300`, `OD_AUX=NWP_NEW` and leaves the record loop body at once (`RETURN` of the record). No `WRITE`, no `SHUTDOWN` in the kernel.
- After the pass the host reduces the status array (max). If non-zero, the host picks the lowest record index with a non-zero status, builds the message with the **same format string** from `OD_AUX` and `SURFACE(...)%ID` of that record's surface, and calls the real `SHUTDOWN`. Records processed after the failing one in the same pass may have been updated; the run aborts anyway.
- The generator rule: a `WRITE` + `SHUTDOWN` pair in a marked range is replaced by a status store plus an exit of the record, with the message template and the argument list recorded in a sidecar field so the host can reproduce the text. Any other `WRITE`, `STOP` or I/O in a marked range is a generator error.
- **Watchdog (new, not in upstream).** The sub-step loop has no iteration cap. On a device a non-terminating record blocks the whole pass. Proposal: a counter `NSUB` in the loop with a large cap (default 10**6, a policy value), status 301 when exceeded. For finite inputs the cap is never reached and the arithmetic is unchanged; for a NaN input upstream hangs and the kernel reports. This is a deliberate deviation and is listed as a decision (T3).

**Tests.**
1. *Reference side.* The mock `SHUTDOWN` writes the message and exits with a fixed code. Each abort case runs in its own process (the harness driver runs reference and kernel as two executables for abort cases); the expected outcome is the message text and the exit code.
2. *Kernel side.* Capacity set below the need (a layer whose `GET_N_LAYER_CELLS` result exceeds `N_CELLS_MAX`): `OD_STATUS=300`, `OD_AUX` equals the `NWP_NEW` printed by the reference, and the host-side message text equals the reference message byte for byte, including the `I5` formatting of `NWP_NEW` and the surface id.
3. *No write outside the record:* the canary bands (section 4 item 4) are intact for the failing record and its neighbours; all other records in the same pass that did not fail are bitwise equal to the reference.
4. *No false positives:* all non-abort cases of every stage have status 0 in all flag sets; `NWP_NEW == N_CELLS_MAX` exactly passes (boundary), `N_CELLS_MAX+1` aborts.
5. *Lowest failing record* is reported when several fail (record order is the list order).
6. *Watchdog:* a NaN in `TMP` and an `Inf` in the flux produce status 301 within the cap, reference run with the cap emulated by the harness (the reference text is the verbatim loop with the same counter inserted by the harness and nothing else); finite cases never set it. A mutant that removes the cap is caught by the NaN case (the test would hang, so the harness runs it under a timeout and treats a timeout as a failure of that mutant's "survival").
7. *Newton non-convergence is not an error:* cases that run 20 iterations without convergence must match the reference bit for bit, including `B2%Y_O2_ITER` staying at its old value (the field is only written on convergence).
8. *Negative generator checks:* a marked range with a `WRITE` that is not part of the recognised `SHUTDOWN` pair, a `STOP`, or a `SHUTDOWN` without a message template is refused with the exit code of the other refusals.

## 9. Mutants and negative checks

Mutants are edits of the generated kernel text (the `KERNEL_MUTANTS` style: `seg_replace(text, kernel, old, new)`, must apply exactly once) and of the table builders (`BUILDER_MUTANTS`). Every mutant must be caught by the unit it targets (and not only by a later one). Minimum list, by unit:

- U1/U12/U13: node weight sign or index shift; geometry factor `R_S(I-1)**I_GRAD - R_S(I)**I_GRAD` terms swapped; stretch exponent `N-1` replaced by `N`; interpolation weight row/column swapped; `INT_WGT` extent; remesh threshold `<` for `<=`, 1.5 factor dropped (2534).
- U3: material loop `CYCLE` condition `<=` to `<`; `ITMP` clamp `I_MAX_TEMP-1` replaced by `I_MAX_TEMP`; `RHO_ADJUSTED` layer index off by one; `VOLSUM` test; `K_S` default 10000 dropped; `PACKING_RATIO` factor dropped.
- U7: tridiagonal forward/backward direction; `DDT` or `AAS` index; `RFACF` coefficient (0.5 or 2 changed); `QDXKF` term dropped; ghost update (`TMP(0)`) dropped; clipping `MIN`/`MAX` swapped; `TMP_F_OLD` assignment moved; `0.5*HTCF*DT_BC_SUB` sign.
- U8: `NINT` replaced by `INT`; `MAX(1,...)` dropped; the `ICYC>WALL_INCREMENT` guard inverted; `DT_FO` dropped from the `MIN`; `Q_S` rebuild condition `/=` to `==`; `T_BC_SUB` update moved after the exit test; exit tolerance `TWENTY_EPSILON_EB` dropped.
- U9: `EXPOSED` reads front HTC instead of back; Dirichlet back ramp index; back emissivity branch order; `Q_LIQUID_B` sign.
- U10: internal radiation sweep direction; `RFLUX_DOWN` formula term; heat-source ramp index; boundary fuel model sign.
- U14/U15: reaction rate exponent (`N_S`, `N_T`, `N_O2`); `EXP` argument sign; species index of `NU_GAS`; `M_DOT_G_PP_ADJUST` versus `ACTUAL`; Newton: `Y_UPPER`/`Y_LOWER` swapped, bisect branch removed, `MAX_ITER` 20 to 19; `B2%Y_O2_ITER` written on non-convergence.
- U16/U17: `M_DOT_LAYER_PP` divisor; `T_IGN` condition; `INT_FTP` exponent; burn-away event dropped or set with the wrong sign; burned-away early return dropped.
- U18: thin-wall back index choice (`_M` before `_P` swapped); isolated flag ignored.
- U19: status not set; status set but the record continues; `OD_AUX` wrong; watchdog cap removed.
- Layout: the builder mutants of section 4 item 7; a kernel that uses `OD_N` where `OD_CAP` is needed for an offset.
- Data: one mutant per input switch that must matter (`CHECK_FO` read as true always; `OXPYRO_MODEL` always off).

**Negative checks** (the generator or the table check must refuse):
- A record whose `OD_N > OD_CAP`; two records with overlapping rows; a surface whose `N_CELLS_MAX` is smaller than `N_CELLS_INI`; a material table shorter than `I_MAX_TEMP+1`.
- A marked range containing an `ALLOCATE` whose extent cannot be bounded from a declared maximum (the private-array rule: extent checked against the upstream `ALLOCATE` text); an `ALLOCATE` with `STAT=` whose status is read (2737 sets `STAT=IZERO` and ignores it: the generator must see the ignore and may then accept).
- A solve called with `HT_DIM>1` (asserted out of AMR mode), or with `CFACE_INDEX`/`PARTICLE_INDEX`.
- A module scalar used by the routine and missing from the policy list (`ICYC`, `WALL_INCREMENT`, `CHECK_FO`, `OXPYRO_MODEL`, `TMPMIN`, `TMPMAX`): generator error naming the symbol.
- A function callee with no flat table; `PRESENT` on a dummy with no binding; an implicit-`SAVE` local (declaration with initial value) not proven assigned before use.
- An atomic or `CRITICAL` in the generated solve (none is allowed; the only shared event, `OBSTRUCTION%MASS=-1.` at 2695, is a slot store, plan row A).

## 10. Tolerance classes and the device-libm question

- **Host, bitwise.** All units on the host in all cells. The routine's sums are `SUM` over small arrays in index order, `MAXVAL`/`MINVAL` (exact), explicit loops; no reduction across records.
- **Device, libm (D-070).** U1, U3, U10 to U15 contain `**` with a real exponent (`R_S**I_GRAD`), variable integer powers (`STRETCH_FACTOR**(N-1)`), `EXP`, `LOG` and `SQRT`. Per D-070 the registry category `libm` applies: per-value tolerance 2 ulp against the host, measured value reported. No own device `pow`.
- **The discontinuity problem.** `NINT(TMP_RATIO)` (2324), the remesh thresholds (2502-2670) and the Newton exit test (`ABS(M_DOT_ERROR)<M_DOT_ERROR_TOL`) turn a 1-ulp difference into a different number of sub-steps, a different node count or one more iteration. The result is then not 2 ulp but a small physical difference, and one record may diverge from the host for the rest of the run. Rule for the device comparison: decisions (sub-step count, `NWP`, layer list, Newton iteration count) are compared first; if they are equal, values are compared with the 2 ulp rule; if a decision differs, the case is reported as a **decision flip** with the input that caused it, not as a failure, and the count of flips per 10**4 record-calls is a metric with a threshold. Proposed threshold: 1 per 10**4 on the test set; to be agreed with the V&V Lead (question T4). On the host there are no flips because reference and kernel share the libm.
- **Owner aggregates and budgets** use the rule of section 11.

## 11. Integration tests with real FDS cases, and the tolerance rule

The cases are those of `05` section 6; run through the capture-and-replay harness (plan section 5.6) first, then as whole runs.

| Case | File | What it exercises | Replay check | Whole-run check |
|---|---|---|---|---|
| `heat_conduction_a` | `Verification/Heat_Transfer/heat_conduction_a.fds` | `SOLID_PHASE_ONLY`, `WALL_INCREMENT=1`, insulated 0.1 m slab, `T_END=2000`, no pyrolysis | records at several steps; sub-step counts equal | wall temperatures vs host FDS: bitwise with the host kernel in the loop; device: libm rule |
| `energy_budget_solid` | `Verification/Energy_Budget/energy_budget_solid.fds` | 1 mm slab, fixed HTC 50, insulated, energy conservation | same | solid energy budget closed to the existing check; run is the FR-022 wall-budget case |
| `back_wall_test` | `Verification/Heat_Transfer/back_wall_test.fds` | 5 mm steel, `EXPOSED` back across four meshes, `INTERPOLATED` boundaries | back snapshot semantics (D-065 Q2): host snapshot mode vs device snapshot mode bitwise; snapshot vs FDS order = algorithm difference, measured | tolerances signed by the V&V Lead (D-065); byte-identical wall temperatures across two layouts (FR-046) |
| `surf_mass_vent_char_cart_fuel` | `Verification/Pyrolysis/surf_mass_vent_char_cart_fuel.fds` | predicted pyrolysis, char, fuel, `WALL_INCREMENT=1`, `DT=0.01`, evaporation-free | `RHO_DOT`, `M_DOT_G_PP_*` per record per step | mass loss and heat release curves vs host FDS |

`couch` and `box_burn_away*` need burn-away (FR-042) and stay cost evidence only. Missing coverage the integration cases do not give (cylindrical or spherical geometry, internal radiation, ramped heat source, boundary fuel model, thin-wall entry, remesh and layer removal, `ERROR(300)`) is covered by the generated cases of section 6 only; this is stated in the report of each stage.

**Tolerance rule (explicit).**
1. **Per record, bitwise** wherever the claim protocol says so: every unit and the whole routine on the host (section 1), including a replay of captured real records: every field of the output list, every call, every flag set. This is the only rule used for acceptance of a kernel.
2. **Device** records: bitwise where no libm call is on the path; otherwise decisions exact and values within 2 ulp (section 10).
3. **Owner aggregates** (area mean `TMP_F`, summed fluxes, obstruction and gas-cell sums built from records): ordered sums of per-record values are bitwise when the order is the ruled one (`wall_sums.py`); aggregates computed by a different but equivalent grouping (a group of `r**2` records against one owner record) use a relative **1e-12** tolerance, as `05` section 5 fixes (a sum of identical values is not always bitwise equal to a multiple).
4. **Budgets** (FR-022): closure to round-off, working threshold 1e-12 relative (`05` section 6, S1 row), unless the V&V Lead sets a looser value for a documented algorithm difference (back wall snapshot, neighbour mass snapshot).
5. **Whole-run comparisons to FDS-only** use the run-level tolerance of the V&V plan, not bit equality, as soon as any documented difference of D-065 enters the case.

## 12. Refinement tests (record fan-out, D-064)

These tests start once the record store (plan `05` WP1/WP2) exists; the solve kernel itself is unchanged (A6), so they test the driver around it. Rulings used: the record-local routine for the `ZZ_F` mass-flux block (D-064 (1), DRAFT), face key `(level, global integer index of the gas-side cell at that level's resolution, IOR in +-1..3)` with sum order `(level,k,j,i,IOR)` (D-064 (2)), accessor of (NM, level) via `POINT_TO_BOX` (D-064 (3)), thin-wall records never deeper than the owner (D-064 (6)).

| # | Test | Check |
|---|---|---|
| F1 | One record per face (Phase 5, depth = owner level) | Record state and `B1` equal to the FDS wall state of a single-mesh equivalent of the same case, bitwise, for the S0 cases (`05` section 6) |
| F2 | Fan-out of one owner to `r**2` records (r = 2 and 4, compounded 2 levels) with the same gas-side inputs | Every record of the group bitwise equal to the others and to the unsplit record advanced alone (identical rows, identical inputs, identical call); owner aggregate within 1e-12 of the unsplit owner value |
| F3 | Deterministic order: group and owner sums in face-key order `(level,k,j,i,IOR)` | Results independent of box layout, rank count, thread count (two layouts, 1 and 4 ranks, 1 and 4 threads, exact-sum switch ON): bitwise on records, bitwise on the ordered owner sums; shuffled order must differ on at least one case (not vacuous) |
| F4 | Per-record HTC from the record's own `TMP_F` | Records with different `TMP_F` in one group get different HTC; mutant: HTC from the owner's `TMP_F` is caught |
| F5 | Record-local mass-flux routine (`CALCULATE_ZZ_F` 1185-1320) at one record per face | Extraction is bit-identical to the in-place block in FDS-only mode (the DRAFT condition of D-064 (1)): replay of captured `B1` before and after, all flag sets |
| F6 | Accessor | A record on a fine mesh object reaches `HEAT_TRANSFER_COEFFICIENT` and the back-side reads without `MESHES(NM)` of a fine number; a mutant that uses `MESHES(NM)` directly aborts as D-056 states and is caught |
| F7 | Split copy at a regrid | Copy is bitwise (A1), `AREA` = parent area / `r**2` exactly, thickness inherited; every record bitwise unchanged across later regrids (A3); a no-change regrid is a no-op |
| F8 | Thin-wall faces | The lateral pass (L1489) runs at owner level; a thin wall whose owner is coarse never gets a deeper record (assert) |
| F9 | Restart by face key | Records written in key order, read back on another rank count: bitwise (A5) |
| F10 | Energy closure | `S1` cases: solid energy change equals net surface flux times area summed over records, to the budget threshold; the area-mean `TMP_F` effect measured (feeds OQ-S5 if it misses) |

Mutants for F2 to F4: group sum in reverse key order; owner area instead of record area; fan-out copies a stale input; group skipped for the last record. All must be caught.

## 13. Staged schedule (working days)

The estimate of `wall-bc-translation-plan.md` section 8.3 (24 to 36 days on the critical path, low confidence) is kept and given a test-design breakdown. Two tracks run in parallel; the critical path is the core track plus the integration of the second track.

**Critical path (core track).**

| Stage | Content | Units | Days |
|---|---|---|---|
| T0 | Harness: cut-and-compile tool for verbatim ranges, mock module, case-file writer/reader, bit compare with the vacuity and coverage gates, field read/write table scanner, flag/thread driver | all | 3-4 |
| T1 | Record store: builder, checks, round trip, canaries; ragged kind in the generator (with the Wall engineer and the Generator Engineer); private/scratch switch | U0 | 3-4 |
| T2 | Core slice: set-up, grid, front flux, `VOID`/`INSULATED` back, implicit update, post-loop; then sub-step control and the explicit estimate; generator: private fixed-size arrays, module scalars, `PRESENT` | U6, U7, U8, U16 | 6-8 |
| T3 | Widening: multi-layer/material, geometries, internal radiation, heat sources and ramps, `EXPOSED` and Dirichlet back, boundary fuel model | U9, U10, U11 | 3-4 |
| T4 | Integrate pyrolysis (units come from the second track) | U14, U15 in U17 | 2-3 |
| T5 | Integrate renoding, layer removal, error path and watchdog | U12, U13, U19 in U17 | 3-4 |
| T6 | Thin-wall entry, back-side snapshot, burn-away slot event | U18 | 2-3 |
| T7 | Full matrix (six sets x two switches x 4/8 threads) on the whole kernel, full mutant list, integration replay on the four cases, tolerance report | U17, U18, section 11 | 2-4 |
| | Sum | | **24-34** |
| | Wait for generator features and for the Wall engineer's HTC callee (not in the stages above; the estimate for features is in the plan: private arrays 1-2, function callees 2-3, inliner 3-4 days, partly parallel) | | 0-2 |
| | **Critical path** | | **24-36** |

**Second track (parallel, off the critical path; one executor each).**

| Stage | Content | Units | Days |
|---|---|---|---|
| B1 | Leaf callees: grid, interpolation, emissivity, property update; with the layout tests | U1, U2, U3 | 4-5 |
| B2 | Pyrolysis: rates, fluxes, then evaporation/char/oxygen, then the Newton loop | U14, U15 | 6-8 |
| B3 | Layer update, remesh decision, regrid, `ERROR(300)` path | U12, U13, U19 | 6-8 |

B1 to B3 total 16 to 21 days of effort, done during T0 to T3 (15 to 20 days) so that the integration stages T4 and T5 start without waiting. If only one executor is available the critical path becomes the sum (about 40 to 55 days); the plan assumes two or three.

**Off the critical path.** Refinement tests F1 to F10 (section 12): 3 to 4 days, after the record store of Role 1 and the driver loop (`05` WP2, WP4), so they do not add to the 24 to 36. The integration replay needs a patched host build to dump records (about 1 day inside T7).

**Milestones and what each gives the project:** after T2 a first device-capable kernel for a plain slab (covers `heat_conduction_a`, `energy_budget_solid`); after T4 the pyrolysis case; after T6 `back_wall_test`; T7 closes the claim.

## 14. Open questions (with owners)

| # | Question | Owner |
|---|---|---|
| T1 | Can the generator build kernels from marked sub-ranges of one routine (block kernels for test builds only) with the surrounding locals as arguments? If not, the fallback is the verbatim-cut harness with switched-off blocks. | GPU Generator Engineer |
| T2 | Is the status-code form of `WRITE` + `SHUTDOWN` (section 8) acceptable as a generator feature, with the message template recorded in the sidecar? | GPU Generator Engineer; the abort policy on the host side: AMR Chief Architect |
| T3 | The sub-step watchdog (cap 10**6, status 301) is a deviation from upstream for non-finite input only. Accept? | AMR Chief Architect |
| T4 | Device libm decision flips (section 10): what rate is acceptable, and is a flip in a sub-step count or node count a failure or a recorded difference? | AMR V&V Lead |
| T5 | Private arrays versus scratch rows for the 25 locals of extent `NWP_MAX` (device local memory budget per thread): which is the default, and what is the largest `NWP_MAX` of the supported cases? | GPU Generator Engineer with the Wall engineer |
| T6 | Layout R1 (capacity offsets) and the per-material stride: confirm, or choose R2. Who builds `OD_*` (plan: Wall engineer builder, Solid Lead review)? | GPU Wall Loops Engineer |
| T7 | Which fields of `B1`/`B2`/`ONE_D` does the routine write (the field scan of section 5.4 gives the list); the Mapper's cross-check against the reads and writes in `WALL_BC` (also needed by `05` WP1). | FDS Legacy Mapper |
| T8 | The accessor of D-064 (3): the Mapper's list of every `MESHES(NM)` lookup on the solid path (the solve has 1880-1883, 1913-1927, 1952-1955, and the HTC callee 3156). | FDS Legacy Mapper |
| T9 | Replay harness: a patched host build that dumps `ONE_D`/`B1`/`B2` before `WALL_BC` (a guarded local patch or a separate tool build). Where does it live? | FDS Legacy Mapper with the Solid Phase Lead |
| T10 | Snapshot of the back side (D-065 Q2): the snapshot is taken at the start of which pass, and is the snapshot record a copy of `B1_BACK` fields only (`TMP_G`, `HEAT_TRANS_COEF`, `Q_RAD_IN`, `LP_CPUA`) or the whole record? The tests assume the four fields read at 2163-2183. | AMR Chief Architect with the Solid Phase Lead |
| T11 | Statistical and decision-flip test design for `MASS_FLUX_VAR` (Q10 of the plan) stays with the Solid Phase Lead; it is not part of the solve (it lives in `CALCULATE_ZZ_F`). | AMR Solid Phase Lead |

## 15. Risks to this design

1. The 24 to 36 days assume the second track in parallel and the generator features arriving on the dates of the plan (section 18.3); with one executor the path is about 40 to 55 days.
2. 771 pointer sites: if the generator cannot rewrite them, the fallback is a hand port, and the verbatim comparison then rests on the replay harness and the generated cases alone (plan section 8.3).
3. The size of the private-array footprint may force the scratch-row form; both are tested, so the cost is a policy switch, not a redesign.
4. Libm decision flips on a device can make records diverge permanently; the rate threshold (T4) must exist before the device gate.
5. Coverage: the generated cases do not reproduce real-case distributions; the four integration cases and the coverage gate are the guard, and the report of each stage lists the tags not reached.
6. The upstream text moves (the solve and `PYROLYSIS` are high-churn): every rebase needs the anchor check, the port-merge check and a rerun of the quick profile, then the full matrix.
