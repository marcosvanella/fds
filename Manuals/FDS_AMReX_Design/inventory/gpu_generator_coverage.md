# GPU kernel generator: coverage after round 3 (cut-cell, cut-face and CFACE loops deferred)

Inputs: `gpu_candidate_loops.csv`, `gpu_loop_blocker_sites.csv`, `gpu_callee_routines.csv`, `gpu_routine_churn.csv`, `gpu_callgraph_survey.md` (FireX `36975d7`). Output: `gpu_generator_loop_classes.csv` (one row per loop; round-2 columns `bucket`, `geometry_reason`, `wall_tier`, `front_end_default`, `front_end_exact`, `tested_kernels`, `status`). Script: `tools/inventory/gpu_generator_coverage.py --round2`. Every time share is the survey's modelled `est_pct_of_total_time` (93.298% over the 829 loops, 0% for `ccib.f90` because its timer is zero), so shares are estimates, not measurements. Target: K2 only (Fortran with OpenMP target). K1 (an AMReX `ParallelFor` backend) is not pursued.

## 1. Scope and denominators
Owner decision: loops keyed on cut cells, cut faces or CFACE are **deferred: geometry refactor** (separate move to structure-of-arrays) and are excluded from targets and from every denominator below.
- **Geometry rule.** A loop is geometry when it is in `ccib.f90`, or its text names `CUT_CELL`, `CUT_FACE`, a `CFACE*`/`CUTCELL*`/`CUTFACE*` name, `CCVAR`, `FCVAR` or a `CC_*` name, or it uses a `CC_*` or `CFACE_TYPE` derived type, or it is tagged as iterating over cut faces and also names one of these. The flag `CC_IBM` alone is not geometry: a guarded call such as `IF (CC_IBM) CALL SET_EXIMRHOZZLIM_3D` (divg.f90:658) keeps loop L0381 (divg.f90:640-660) in scope. `MESH_EXCHANGE` and `pois.f90` loops stay in the denominator as retired or replaced even where their text names these arrays.
- **Deferred bucket: 378 loops, 23.854% of total time.** 308 are in `ccib.f90` (0% time, timer zero). The other 70 (23.854%): `pres.f90` 44 (13.703), `radi.f90` 4 (4.600), `velo.f90` 1 (2.612, `MATCH_VELOCITY`), `part.f90` 8 (1.529), `wall.f90` 1 (0.694, `WALL_BC` CFACE loop), `divg.f90` 4 (0.628), `fire.f90` 5 (0.087), `hvac.f90`, `soot.f90`, `vege.f90` 1 each (about 0). By former subclass (outside `ccib.f90`): module-scalar reductions 3 (14.930), wall/cface layout 48 (4.471), MPI/solver library 3 (2.283), callee 6 (1.472), particles 2, control 5, cut-cell 1, A 2.
- **New denominator: 451 loops, 69.444% of total time** (829 - 378; 93.298 - 23.854). "Translation-eligible" is the part of it that needs kernels: 451 - 207 retired or replaced (`MESH_EXCHANGE` 2 loops 16.536%, `pois.f90` 205 loops 1.796%) - 13 host-side (6 control, 2 I/O, 5 particle; 1.935%) = **231 loops, 49.177%**.

## 2. Headline (non-geometry, 451 loops, 69.444%)
| Class | Loops | Share of loops | Est. time (% of total) | Share of 69.444 |
|---|---|---|---|---|
| A | 128 | 28.4% | 6.308 | 9.1% |
| B | 20 | 4.4% | 1.894 | 2.7% |
| C | 303 | 67.2% | 61.242 | 88.2% |
Geometry bucket (not in the table): A 23, B 18, C 337 loops.

Classes (as in round 1): A = regular flat-array loop; B = needs a primer hygiene change (items 1, 7, 8) or a marker assertion `M`; C = needs flat tables or a layout change, stays on host, or is retired. Change from round 1: 5 loops moved C to B (the fifth, in round 3, is the `B2` wall loop at part.f90:3654, L0906) (the wall loops that the WALL tables now make translatable; item 7), and the 378 geometry loops moved out.

C subclasses (non-geometry): wall/cface layout (now: wall loops over `WALL`/`BOUNDARY_*`/`EXTERNAL_WALL`) 44 loops 22.367%; `MESH_EXCHANGE` 2, 16.536; callee reads components 14, 8.232; probe-found layout blockers 12, 7.174; `pois.f90` 205, 1.796; element pointer 1, 1.325; mesh component 6, 1.297; control flow 6, 0.739; particles 5, 0.670; other layout 7, 0.580; I/O 2, 0.526. No non-geometry loop is a module-scalar reduction or an MPI/solver-library loop (all six are in the geometry bucket), but 52 non-geometry loops contain a reduction pattern (section 5).

## 3. Coverage today (what the generator actually translates and tests)
| Status | Loops | Est. time % |
|---|---|---|
| translated and bitwise-tested (A 11, B 8, C 6) | 25 | 1.227 |
| partly: an inner K,J,I nest is translated and tested, the outer survey loop is not (L0365 divg.f90:128-235, L0369 divg.f90:287-421, L0381 divg.f90:640-660, L0401 divg.f90:995, L0877 mass.f90:868, L0880 mass.f90:65) | 6 | 21.628 |
| front end accepts, no test yet (A 46, B 3) | 49 | 1.055 |
| classified translatable, front end rejects (generator gap; A 71, B 9) | 80 | 6.003 |
| not translatable today (retired 207 and host 13 included) | 291 | 39.531 |
- 34 kernels are generated: 7 from round 1, 18 from round 2 and 9 from round 3. 25 map one-to-one to a survey loop, 8 sit inside the six partly covered loops, and `conductivity` is a whole-array assignment (divg.f90:462) with no loop. Inner-nest time is not separated by the survey model, so the 21.628% of the partly covered loops is an upper bound for what their inner kernels cover.
- **Tested** here means the generated kernel equals the verbatim upstream loop text bit for bit (section 6), not that it ran on a GPU.
- "Classified" is not "translated". Of the 148 A and B loops, 19 are translated and tested, 49 more pass the front end (1.055%), and 80 are rejected for generator gaps. The class-C wall loops add 6 more tested loops (the 25 above are A 11, B 8, C 6).

## 4. What the WALL flat tables unlock
Tables and fields: `docs/adr/drafts/gpu-generator-design.md`, section "WALL flat tables". Tier = what a wall-keyed C loop needs (non-geometry, translation-eligible):
| Tier | Meaning | Loops | Est. time % |
|---|---|---|---|
| T1 | every component read or written is a scalar or a species vector of `WALL`, `BOUNDARY_COORD`, `BOUNDARY_PROP1`, `EXTERNAL_WALL`, `BOUNDARY_PROP2`, `SURFACE`, `SPECIES_MIXTURE` or `CELL`: flat tables are enough | 32 | 23.882 |
| T3 | a callee reads components: the callee must be flattened too | 14 | 8.232 |
| T4 | other derived data: `OMESH`, `BOUNDARY_ONE_D`, `BOUNDARY_THR_D`, `ZONE_MESH`, `EDGE`, thin wall, mesh records | 36 | 8.289 |
| T2 | ragged per-record arrays (`SYNTHETIC_TURBULENCE`, turb.f90) | 1 | 0.572 |
- **T1 by routine** (loops): `DIVERGENCE_PART_1` 5 (outer loops L0365, L0369 included), `MASS_FINITE_DIFFERENCES` 2 (L0880 mass.f90:65-192 4.832%, L0882 mass.f90:224-320 0.302%), `SPECIES_ADVECTION_PART_1_NEW` 2 (L0401 3.260%), `CHECK_MASS_DENSITY` 2 (L0877 mass.f90:868-939 1.967%), `SETTLING_VELOCITY` 1 (1.325%), `COMPUTE_VISCOSITY` 2, `WALL_VELOCITY_NO_GRADH` 3, `ENTHALPY_ADVECTION_NEW` 2, `SPECIES_ADVECTION_PART_2` 2, `CHECK_DIVERGENCE` 1, `NO_FLUX` 1, `CHECK_STABILITY` 1, `COMPUTE_STRAIN_RATE` 1, `CORIOLIS_FORCE` 1, `MERGE_PRESSURE_ZONES` 1, `DIVERGENCE_PART_2` 1, and small `pres.f90`/`part.f90`/`radi.f90` wall loops (5). **Round 3 status of the T1 tier (32 loops; L0906 left it for class B).** Whole-loop translated and bitwise-tested: **6** (L1400 velo.f90:3347-3368, L0372 divg.f90:482-486, and the round-3 loops L1122 pres.f90:4351, L1361 velo.f90:440, L1380 velo.f90:981, L1401 velo.f90:3378). L0906 (part.f90:3654, `B2`) is also translated and tested but is now class B. Partly covered through per-nest markers: 5 (L0365, L0369, L0401, L0877, L0880; for the last three only the cell sub-nests: the `RHO_Z_P` product at mass.f90:70 and divg.f90:998, the `DELTA_RHO_ZZ` zero fill at mass.f90:871 and the `RHO_ZZ` clip assign at mass.f90:931). Of the 29 T1 loops that were blocked before round 3, 5 are now translated and tested (L0906, L1122, L1361, L1380, L1401; 0.062% modelled time) and 3 more are partly covered; 21 stay untranslated. What blocks those 21:
  - reduction or outer-loop `CYCLE`: L0363 (divg.f90:1675), L1348 (velo.f90:3093);
  - race or read-modify-write: L0386 (`USUM(IPZ)`), L1358 (`CELL_COUNTER`), L0375 (accumulation through wall subscripts), the `DELTA_RHO_ZZ` scatter of L0877, and L1359 (velo.f90:355: a `UNIQUE` contract cannot be proven);
  - pointer or callee: L0398 (`U_TEMP=>U_WORK` and a pointer callee), the wall nests of L0880 and L0401 (`FX_P` remap pointer, `Z_TEMP` array constructors), L0405 (`UU` associated outside the routine), L1272 (element pointer with callees);
  - whole-array assignment inside the loop body: L0403, L0882;
  - designators with no table: rank-2 `CONNECTED_ZONES` (L0400), `EXTERNAL_WALL(IW)%NOM` (L1366), the `BOUNDARY_PROP1(WC%BC_INDEX)` alias (L0394), `CELL(IC)%SOLID` cell loops (L0399, L0406), a K,J-only nest (L0876);
  - the rank-4 LOGICAL array `LOG_INTWC` (L1130, L1154) and `WC%B1_INDEX==0` used as a value (L1239).
- **T1 plus the 4 promoted loops.** Four loops moved from C to B in round 2 (a fifth, L0906, in round 3): fire.f90:1957 (L0608, `B1%Q_CONDENSE = 0` over the walls; front end accepts, no test yet), mass.f90:424 (L0858), mass.f90:595 (L0867) and velo.f90:151 (L1353); the last three are tested. They are class B item 7.
- **T3, `WALL_BC` (wall.f90:26-279).** 4 non-geometry loops (L1486 wall.f90:109-119, L1487 :124-137, L1488 :144-185 2.688%, L1489 :192-194 1.370%; 4.819% together) call `ASSIGN_GHOST_VALUE`, `NEAR_SURFACE_GAS_VARIABLES`, `SURFACE_HEAT_TRANSFER` (wall.f90:603), `SOLID_HEAT_TRANSFER` (wall.f90:1809), `CALCULATE_ZZ_F`, `CALCULATE_RHO_F`, `CALCULATE_RHO_D_F`, `WALL_MODEL`, `CALC_DEPOSITION`, `CALC_HVAC_BC` (called from the wall loops at wall.f90:109-194). The WALL tables carry the wall fields; the callees also read `SURFACE`, material and `BOUNDARY_ONE_D` records, so `WALL_BC` is a hand port with generator help (not unlocked by the tables alone). The CFACE loop (wall.f90:201-239, 0.694%) is deferred and the particle loop (wall.f90:247) is host-side.
- **`THERMAL_BC`, `PYROLYSIS`, `SOLID_HEAT_TRANSFER_3D`.** There is no routine named `THERMAL_BC` or `SOLID_HEAT_TRANSFER_3D` at `36975d7`; the thermal switch is `SF%THERMAL_BC_INDEX` (wall.f90:117, 153), and the 3-D conduction exchange is `HT3D_TEMPERATURE_EXCHANGE` (wall.f90:3613; 2 loops, 0.965%, tier T4: needs `BOUNDARY_ONE_D` layer arrays and `BOUNDARY_THR_D`). `SOLID_HEAT_TRANSFER` (3 loops) and `PYROLYSIS`/`PERFORM_PYROLYSIS` (7 loops) iterate over layers and nodes, not mesh space; the survey marks them non-candidates, so they are outside the 829 and outside every share here. Ragged per-layer data (`ONE_D%`) is a second table family (offset and count per wall) that the scalar WALL tables do not provide.
- **`VELOCITY_BC` (velo.f90:1858), `VISCOSITY_BC` (velo.f90:514), `NO_FLUX` (velo.f90:1376), `MATCH_VELOCITY_FLUX`, `PATCH_VELOCITY_FLUX`** are T4: they also read neighbour-mesh arrays (`OMESH(NOM)%U`, `%US`, `%D`, `%KRES`, `%MU`, `%H`, `%HS`) through `EXTERNAL_WALL%IIO_MIN..KKO_MAX`. They wait for the exchange-buffer layout (ADR-001), not for the WALL tables.
- Cross-dependency: `UVW_SAVE`, `U_GHOST`, `V_GHOST`, `W_GHOST` (read by 3 translated wall loops) are written in `ccib.f90` (deferred) and in `velo.f90:2732-2832` (`MATCH_VELOCITY`, whose cut-face loop makes it geometry).

## 5. What remains host-side (realistic list)
| Category | Loops | Est. time % | Notes |
|---|---|---|---|
| retired: `MESH_EXCHANGE` (main.f90) | 2 | 16.536 | replaced by `FillBoundary` (ADR-001) |
| replaced: `pois.f90` pencil loops | 205 | 1.796 | replaced by `amrex::FFT::Poisson` |
| sequential control flow (goto, while, exit, return) | 6 | 0.739 | `INSERT_ALL_PARTICLES`, `CALC_AGGLOMERATION` (soot.f90), `HVAC_BC_IN`, `COMBUSTION_GENERAL_LOAD_BALANCED`, `REMOVE_PARTICLES`, `REMOVE_OLDEST_PARTICLE` |
| particles | 5 | 0.670 | `ParticleContainer`; includes `RADIATION_FVM`/`INTERPOLATE_IL` particle-like loops and `WALL_BC` wall.f90:247 |
| file I/O | 2 | 0.526 | `OPEN_AND_CLOSE`, `WRITE_EWC_TYPE_DIAGNOSTIC` |
| MPI or solver library | 0 here; 3 deferred | 2.283 deferred | `GLMAT_SOLVER`, `GET_BCS_H_MATRIX`, `GET_MATRIX_INDEXES_H` are in the geometry bucket |
| module-scalar reductions | 0 here; 3 deferred | 14.930 deferred | `GET_H_MATRIX`, `GET_MATRIXGRAPH_H_WHLDOM`, `RADIATION_FVM` |
- **Reductions inside translatable loops.** 52 non-geometry loops carry a reduction pattern (`reduction_accum` 36, `reduction_intrinsic` 18, some both). 24 are retired or host-side already (17 `pois.f90`, 2 `MESH_EXCHANGE`, 5 host or particle). Of the remaining 28, 7 are class B whose `SUM`/`MAXLOC` the generator expands in source order and tests bitwise (mass.f90:328-350, 691, divg.f90:1191-1213 and others), and 21 are C that need a `REDUCTION(...)` marker attribute or a hand-written reduction (for example `CHECK_DIVERGENCE` divg.f90:1690-1698, `CHECK_STABILITY` velo.f90:3093-3108). The marker attribute `REDUCTION(op:var)` is parsed and emitted as an OpenMP `reduction` clause; no loop uses it yet, and a floating-point reduction is not order-stable on a GPU, so using it needs the exact-accumulation treatment of ADR-001 (FR-005 (ii)) and is outside the bitwise tests.
- **Not host-side but not wall-table work either:** I/O or allocation inside callees (5 loops, 4.08% callee-only), `alloc_ptr_local` flagged lexically in 172 loops (mostly harmless views).

## 6. Bitwise results (rounds 2 and 3)
The 18 round-2 kernels against the verbatim upstream loop text (module `s5_r2_ref.F90`, built from the real derived types of `type.f90`, the real `cons.f90` parameters, and a flatten step that gathers `WALL`/`BOUNDARY_COORD`/`BOUNDARY_PROP1`/`EXTERNAL_WALL` into the flat tables), 234 cases per flag set (6 box sizes, `NL` 0 or 1, internal walls on or off, species count 1-N), 6 gfortran flag sets, plus the 7 round-1 kernels: all bit-equal, no vacuous case. Round 3 adds 9 kernels (34 in total) and raises the count to 408 cases per flag set; all 6 flag sets and all 12 callee-switch combinations pass with no vacuous case, and the round-1 suite (112 comparisons) is unchanged. 17 deliberate mutants of the generated module (7 earlier, 10 new) are all caught. The generator's negative checks now number 42. All 18 are testable with synthetic inputs; none was untestable. Limits: the wall geometry is synthetic and race-free by construction (no real mesh), the reference runs serially in upstream order, `adv_flux_avg` needs pre-filled `ADV_F*` to avoid a no-op, and no GPU or nvfortran run exists.

## 7. Why the front end still rejects the 163 translation-eligible loops it does not accept
| Reason (first failure) | Loops | Est. time % |
|---|---|---|
| nest is not a perfect K,J,I nest (outer species/zone loop, or wall plus cell mix) | 57 | 31.832 |
| pointer set outside the routine, or alias not in the recognised set (`B2`, `SF`, `OM`, `BR`, `U_TEMP`) | 27 | 7.574 |
| scalar live-out or loop-carried (`IPZ`, `II`, `UVWMAX`) | 12 | 2.511 |
| derived-type designator without a table (`CELL%EDGE_INDEX`, `OMESH`, `ZONE_MESH`) | 13 | 1.840 |
| function reference (callee is a function) | 8 | 1.178 |
| whole-array assignment in the body | 3 | 0.569 |
| rank or subscript contract (rank-2 arrays, non-loop-variable subscripts) | 18 | 0.515 |
| other, `CYCLE`/`EXIT` out of nest 2, no `ALLOCATE` found for `layout=exact` (needs a `policy.arrays` line) 8, wall bounds 3, wall write without `UNIQUE`/`OMP DO` 5 | 25 | 1.000 |
The first row is the big one: the marker route (one marker per inner nest, as done for `DIVERGENCE_PART_1` and `DENSITY`) replaces whole-loop translation, so loop counts understate the work already covered.

## 8. Revised path to about 80% of non-geometry eligible loops and time
Denominator A: the 231 translation-eligible loops (49.177%). Denominator B: all 451 non-geometry loops (69.444%), counting retired or replaced loops as done.
| Step | Adds (loops, time %) | Cumulative eligible loops | Cumulative eligible time | Cumulative of 451 / 69.444 (with retired 207, 18.332) |
|---|---|---|---|---|
| A and B loops (147), all of them | 147, 8.202 | 63.6% | 16.7% | 78.5% / 38.2% |
| + T1 wall loops (WALL tables, 32) | 32, 23.882 | 77.9% | 65.2% | 85.8% / 72.6% |
| + T3 callee loops (14, mainly `WALL_BC`) | 14, 8.232 | 84.0% | 82.0% | 88.9% / 84.5% |
| + T4 and T2 (37; needs `OMESH` exchange buffers, ONE_D tables) | 37, 8.861 | 100% | 100% | 97.1% / 97.2% |
(The 13 host-side loops, 1.935%, never reach 100%.) The revised statement: **about 80% of eligible time needs A, B, T1 and T3 (194 loops, 82.0% of eligible time, 84.5% of the non-geometry denominator once the retired loops are counted); 80% of eligible loops needs A, B and T1 (180 loops, 77.9%) plus about 5 of the T3 loops.** Without the WALL tables the ceiling is A and B: 16.7% of eligible time. Concrete work, in order:
1. Test the 49 loops the front end already accepts (1.055%): mechanical with `make_r2_tests.py`.
2. Close the cheap generator gaps (rows 3 to 7 of section 7: function references 8, `policy.arrays` lines 8, whole-array assignment 3, rank contract 18, live-out scalars via `PRIVATE`): 49 loops, about 5.0 points of total time (mostly A and B).
3. Wall T1 (done in round 3 for 5 loops and for the cell sub-nests of 3 more): `B2`/`SF` aliases, internal-only wall bounds (`N_EXTERNAL_WALL_CELLS+1:`), `UNIQUE` and `IDEMPOTENT` markers, `policy.assoc`, per-nest markers for L0880, L0401, L0877. The 21 T1 loops still open are listed in section 4.
4. T3: flatten the `WALL_BC` callees (hand port, generator for the arithmetic).
5. T4: after the exchange-buffer layout exists.
Priority ordering (translation-eligible routines by time; churn from `gpu_routine_churn.csv`, commits excluding bulk, lines in nests, body-only/surface hunks):
| Routine | Loops (A/B/C) | Time % | Churn | Comment |
|---|---|---|---|---|
| `divg.f90::DIVERGENCE_PART_1` | 20 (11/2/7) | 12.844 | 2, 4, 2/0 | low churn; 9 kernels exist |
| `mass.f90::MASS_FINITE_DIFFERENCES` | 4 (1/1/2) | 5.333 | 1, 0, 0/0 | rank 3 in round 1, now 2: L0880 needs inner markers |
| `wall.f90::WALL_BC` | 4 (0/0/4) | 4.819 | 5, 15, 2/5 | T3; 5 surface hunks, hand port |
| `divg.f90::SPECIES_ADVECTION_PART_1_NEW` | 4 (1/1/2) | 3.626 | 6, 333, 1/31 | high churn: last |
| `mass.f90::CHECK_MASS_DENSITY` | 3 (0/1/2) | 2.666 | 3, 13, 7/0 | cheap: body-only hunks |
| `init.f90::REASSIGN_WALL_CELLS` | 2 (0/0/2) | 2.037 | 0 | T3/T4 |
| `wall.f90::SURFACE_HEAT_TRANSFER` | 1 | 1.832 | 14, 164, 11/2 | T4, highest churn |
| `velo.f90::COMPUTE_VISCOSITY` | 10 (7/1/2) | 1.818 | 2, 113, 3/2 | 4 kernels exist |
| `soot.f90::SETTLING_VELOCITY` | 1 | 1.325 | 0 | T1 |
| `mass.f90::DENSITY` | 16 (12/4/0) | 1.230 | 1, 24, 0/2 | 10 kernels exist |
Round-1 entries that left the list because they are geometry: `GET_H_MATRIX` 5.19%, `GET_MATRIXGRAPH_H_WHLDOM` 5.16%, `RADIATION_FVM` 4.58%, `MATCH_VELOCITY` 2.61%; `WALL_BC` fell from 6.18% to 4.82% (CFACE loop deferred, particle loop host-side).

## 9. Hygiene items (primer §2 numbering and titles)
Primer: `docs/adr/drafts/gpu-friendly-fds-primer.md`. `M` = marker assertion.
| Item | Title | Generator status | Non-geometry loops needing it (lexical) |
|---|---|---|---|
| 1 | Hot loops in own routines with flat arrays as arguments with INTENT and explicit bounds | implemented: `layout = "exact"` takes bounds from the `ALLOCATE` text (`-1:IBP1+1`, `NL:NS`, `0:IBAR`), `exact_check` verifies stencil reach | 29 pointer/allocatable dummies (6 B, 23 C) |
| 2 | No ALLOCATE, I/O, STOP or RANDOM_NUMBER inside loops | checked, refused; private fixed-extent copies of local allocatables | 11 `alloc` + 11 `io` loops, all C |
| 3 | Small PURE leaf functions | not a gate: callees emitted with `declare target` | 8 function references block the front end |
| 4 | Explicit INTENT | derived from use | unchanged from round 1 |
| 5 | Replace a derived-type indirection in a hot loop by a plain array | implemented for `CELL(CELL_INDEX())%SOLID` (mask) and `SPECIES_MIXTURE%MW` (`MW_SPEC`) | 13 designator rejections |
| 6 | Avoid module state in hot loops | module scalars become arguments via `[policy.rename]`/tables | 12 live-out scalar rejections |
| 7 | Split wall loops from cell loops and give each wall loop flat index arrays | implemented: WALL flat tables; a wall race needs `OMP DO` upstream or a `UNIQUE` marker | 4 B (promoted), 32 T1 C (5 translated in round 3) |
| 8 | Avoid pointer re-aiming (`UU => US`) inside a routine | implemented: `X => M%Y` and preamble aliases resolved; pointers set outside the routine are diagnosed | 78 loops (6 B, 72 C) |
Items 2, 3, 4 move no B loop (as in round 1).

## 10. Survey-rule over-reporting, corrected
- `iterates_over = cfaces` was assigned from the loop variable name alone. `func.f90::PACK_CELL` (L0677, func.f90:5183-5201) loops over `ICC=1,CELL_COUNT(NM)` of regular cells and names no `CFACE`/`CUT_` array; it is non-geometry, class C (derived-type component in the body). The geometry rule therefore requires a CC/CFACE name or type in the text as well.
- `dt_alloc_component` is lexical: of 62 non-geometry loops flagged, 5 touch only scalar components per the parser (3.007% time; L0381 divg.f90:640 `CELL%SOLID` and `SPECIES_MIXTURE%RCON`, `WALL_BC` L1486/L1487/L1491, `CHECK_UNSUPPORTED_MESH`); 8 have no parse profile. 30 non-geometry C loops (5.329%) use scalar components only (21 loops, 0.657% in the pure layout subclasses); 33 are T1 and wall tables make them translatable.
- `write_derived_index` flags 138 non-geometry loops; 22 end as A or B (15 A, 7 B) because the "index" is a stencil or a face-offset index. `alloc_ptr_local` flags 172 (68 A, 9 B): harmless views.
- `reduction_accum` is correct for 36 loops (all C) but also fires on counters; the expandable `SUM`/`MAXLOC` forms are 7 class-B loops.
- The 5 `scalar self-update` A-to-C overrides of round 1 (`RESMAX` etc.) are unchanged.

## 11. Spot checks against source (read at `36975d7`)
Classification changes (14):
1. L0608 fire.f90:1957-1960 C to B item 7: `WC=>WALL(IW)`, `B1=>BOUNDARY_PROP1(WC%B1_INDEX)`, `B1%Q_CONDENSE = 0`.
2. L0858 mass.f90:424-436 C to B: `IF (WC%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE`, `SELECT CASE(BC%IOR)`.
3. L0867 mass.f90:595-608 C to B: `EWC%BOUNDARY_TYPE_PREVIOUS`, `UVW_SAVE(IW)`.
4. L1353 velo.f90:151-164 C to B: `U_GHOST(IW)`, `BC%II` writes.
5. L0385 divg.f90:730-753 to geometry: `CCVAR(I,J,K,CC_CGSC) == CC_SOLID` at divg.f90:742.
6. L0391 divg.f90:1525-1543 to geometry: `CC_CGSC` test at divg.f90:1532.
7. L0393 divg.f90:1560-1567 A to geometry: `CCVAR(I,J,K,CC_CGSC)` at divg.f90:1563.
8. L1362 velo.f90:2663-2853 `MATCH_VELOCITY` to geometry: `M2%FCVAR(...,CC_IDCF,...)`/`M2%CUT_FACE(ICF)%ALPHA_CF` at velo.f90:2824-2825.
9. L1490 wall.f90:201-239 to geometry: `CFACE_LOOP: DO ICF=INTERNAL_CFACE_CELLS_LB+1,...` at wall.f90:201.
10. L1128 pres.f90:5038-5260 `GET_H_MATRIX` to geometry: `CCVAR(IIG,JJG,KKG,CC_CGSC)` at pres.f90:5154.
11. L1242 radi.f90:3912-4953 `RADIATION_FVM` to geometry: `CCVAR(I,J,K,CC_CGSC)==CC_SOLID` at radi.f90:4009.
12. L0893 part.f90:1310-1445 to geometry: `CCVAR(II,JJ,KK,CC_CGSC)`, `CC_IDCC` at part.f90:1326-1327.
13. L1149 pres.f90:5863-5869 A to geometry: `CCVAR(I,J,K,CGSC) = IS_SOLID` at pres.f90:5866.
14. Kept non-geometry: L0381 divg.f90:640-660 (only `IF (CC_IBM) CALL ...` at 658) and L0677 func.f90:5183 (no cut-cell name).
Generated constructs (12; line numbers in `generated/s5gen_k2.F90` against upstream):
1. `rho_sum` 452-472 vs mass.f90:687-694: `SUM(ZZ(I,J,K,1:N_TRACKED_SPECIES))` becomes `DOT1` accumulated `M=1..NS` in source order; `CELL(CELL_INDEX(I,J,K))%SOLID` becomes `SOLID(I,J,K)/=0`.
2. `mu_dns` 617-637 vs velo.f90:80-88: `ZZ_GET(1:N_TRACKED_SPECIES) = ZZP(...)` becomes an `M` loop into private `ZZ_GET(NS_MAX)`; `GET_VISCOSITY` receives `MU_RSQMW_Z`, `RSQ_MW_Z`, `I_MAX_TEMP`, `NS` explicitly.
3. `flux_mw_fix` 475-543 vs mass.f90:328-350: `MAXLOC(...,1)` becomes a strict-`>` scan (first maximum wins, as `MAXLOC`); `SPECIES_MIXTURE(N)%MW` becomes `MW_SPEC(N)`; `MAX(0, FX0 - SUM1 - SUM2)` becomes `(FX0 - DOT1) - DOT2`.
4. `kres` 660-678 vs velo.f90:285-294: `0.5*(U2+V2+W2)` becomes `0.5*((U2+V2)+W2)`.
5. `zz_corr` 397-427 vs mass.f90:619-632: species loop `N` outermost, `collapse(4)`, `RHS` private, the right-hand side wrapped left to right.
6. `wall_uvw_interp` 681-706 vs mass.f90:424-436: the aliases `WC=>WALL(IW)`, `BC=>BOUNDARY_COORD(WC%BC_INDEX)` are elided; `WC%BOUNDARY_TYPE` is `W_BOUNDARY_TYPE(IW)`, `BC%IIG` is `BC_IIG(IW)`; `INTERPOLATED_BOUNDARY` stays a named constant (cons.f90:93).
7. `wall_rho_d_dzdn` 777-841 vs divg.f90:192-232: `B1%RHO_D_DZDN_F(N)` becomes `B1_RHO_D_DZDN_F(IW,N)` with dummy `(NWE+NWI,NS)`; `EWC%NIC` is `EW_NIC(IW)` with extent `NWE`; `WC%THIN` is the integer mask `W_THIN`.
8. `wall_rho_d_dzdn` dummies for `RHO_D_DZDX..Z`: `(0:IBAR+1,0:JBAR+1,0:KBAR+1,NL:NS)`, the 4-D lower-bound plumbing.
9. `wall_un_store` 844-872 vs velo.f90:3347-3368: loop bound `N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS_AUX` becomes `NWE + N_INTERNAL_WALL_CELLS_AUX`, `UN_WALLS` extent the same (allocated at velo.f90:3344), `IIG,JJG,KKG,IOR` private.
10. `wall_up_ghost` 738-761 vs velo.f90:151-164: `U_GHOST(NWE)`, `UP(BC_II(IW),BC_JJ(IW),BC_KK(IW))`, `SELECT CASE` over `W_BOUNDARY_TYPE(IW)`.
11. `wall_kp_ghost` 763-775 vs divg.f90:482-486: `KP(BC%II,...) = KP(BC%IIG,...)` becomes the same with `BC_II(IW)`; the `BOUNDARY_LOOP` label is dropped.
12. `wall_uvw_interp` dummies `UU(-1:IBAR+1,0:JBAR+1,0:KBAR+1)`, `VV(0:IBAR+1,-1:JBAR+1,0:KBAR+1)`, `WW(...,-1:KBAR+1)`: bounds from `ALLOCATE(M%U(-1:IBP1,...))` (init.f90, `exact` layout).
Assumption inherited from upstream: `wall_rho_d_dzdn` reads `EW_NIC(IW)` only where `W_BOUNDARY_TYPE` is `OPEN_BOUNDARY` or `INTERPOLATED_BOUNDARY`, which upstream also reads from `EXTERNAL_WALL(IW)`, allocated only for external walls (init.f90:38); the driver must keep that true for internal walls.

## 12. Validation of the classification (unchanged from round 1)
- Sixty-loop random sample (fixed seed): 46 agree, 14 differ (23%): A 18/20, B 8/20, C 20/20. All 206 A and B loops probed: 165 agree, 41 differ (20%). The probe's rules were refined while reading disagreements, so these figures are not an independent accuracy measure.
- Round 2 adds the front-end run on every A and B loop in both layouts and the checks of sections 6 and 11.

## 13. Translator risks
- **Aliasing.** Distinctness comes from the `POINT_TO_MESH` targets; a loop writing one work array and reading another that share storage (`BETAHAT11`/`M11` both aimed at `WORK1`, turb.f90:603 and 626) stays refused; `RHOP`/`UU`/`ZZP` are re-aimed per step (divg.f90:60-76), so the driver must pass the right array (hygiene item 8).
- **Wall races.** The generator requires upstream `!$OMP DO` or a `UNIQUE` marker for any write through wall subscripts; `UNIQUE` is an assertion that distinct walls touch distinct cells, argued geometrically (one face per wall cell; thin pairs processed once through `IF (WC%THIN .AND. BC%IOR<0) CYCLE`, divg.f90:196), tested only on synthetic walls.
- **Layout policy.** `layout = "amrex0"` assumes `0:IBAR+1` boxes; `layout = "exact"` takes bounds from the `ALLOCATE` text; the shim must pass the matching sub-box.
- **nvfortran.** None of the generated kernels has run on a GPU. Known S4/S4d findings (teams loop with callees, private scalars `intent(out)`, `-O2` reassociation, `-Minline`) apply; the wall kernels are new and untried.
- **Syntax drift.** fparser2 parses 33 of 34 sources after one normalisation; `pois.f90:7103` does not (host-only). The anchor check, golden signatures and exit codes (2, 3, 4) catch drift.
