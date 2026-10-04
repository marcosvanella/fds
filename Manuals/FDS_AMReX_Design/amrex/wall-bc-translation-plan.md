# WALL_BC family: translation plan for the source generator (L1486, L1487, L1488, L1489, L1485, L0823)

Status: proposal v0.1 for review by the GPU Wall Loops Engineer, the GPU Generator Engineer, the FDS Legacy Mapper and the AMR Chief Architect.
Author role: AMR Solid Phase Lead. Docs only: nothing in `src`, in the generator or in other roles' documents was changed.

## 0. Revision and citations

- Every line number below is from the survey revision **FireX `36975d765f`** (the revision named in `loop-work-list.md`), read with `git show 36975d765f:Source/<file>`. `wall.f90` unless a file is named.
- The working tree has moved on (`FDS-AMReX` branch, a later head). Its `Source/wall.f90` carries the guarded AMReX hooks (patches 0003/0004/0008): a six-line `#ifdef WITH_AMREX` header, a `USE FDS_AMREX_HOOKS` line and `IF (EXTERNAL_GHOSTS_FILLED) RETURN` in `ASSIGN_GHOST_VALUE`. Line numbers in that tree are **9 higher** at the two `WALL_BC` wall loops and **12 higher** from `SURFACE_HEAT_TRANSFER` on. Anything written against that tree (the notes in `wall_sums.py` and its README, and the line cites in `04-blocked-loop-signoff.md`, see section 17) has to be read with this offset. The generator reads the pinned commit, so sidecar anchors use the survey numbers.
- Modelled shares are the survey-model estimates of `loop_work_list.csv`, not measurements. The measured picture (`03-g2a-cost-check.md`, `05-fine-level-solid-plan.md`): the wall timer is 5-23 % of a step and the 1-D solve inside it is 2-13 %, so the gas-side wall work is the larger device target in three of the four cases measured.

## 1. Summary

| Loop | Lines | Modelled share | What it is | Verdict | Effort (days) |
|---|---|---|---|---|---|
| L1486 | 109-119 | 0.41 % | first wall pass: ghost values, near-surface gas state, heat transfer coefficient (HTC) of thick walls | translate, after generator inlining of the call tree | 3-4 (+3-4 HTC) |
| L1487 | 124-137 | 0.35 % | same, for thin-wall cells; reads two other walls' records | translate, as a separate pass after L1486 | 1.5-2 |
| L1488 | 144-185 | 2.69 % | per-wall dispatcher: ten call sites, three of them with shared writes | translate as a sequence of pass kernels plus the ordered sums of `wall_sums.py`, not as one kernel | 35-52 (mostly the 1-D solve) |
| L1489 | 192-194 | 1.37 % | lateral 1-D solve of thin obstructions (**not** what the register says, see section 17) | translate with the solve; one extra shared event | 1 (after the solve port) |
| L1485 | 888-947 | 1.83 % | coarse-side species flux match of the `INTERPOLATED_BC` case, `EWC%NIC>1` only | **retire** in the AMR route (host abort guard); keep only if the Architect confirms `NIC>1` can occur | 0.5, or 4-5 if kept |
| L0823 | `init.f90` 4940-5054 | 1.73 % | creation/removal of obstructions: re-types wall cells, allocates records | **host only**; the device deliverable is the table-refresh trigger | 1-2 (host/driver) |

Totals: about 74-112 person-days over four roles (section 13). Critical path is the 1-D solve chain (ragged record store, generator support for fixed-size private arrays, the port itself), about 24-36 working days, independent of the other generator features.

## 2. How `WALL_BC` is organised (26-279)

- Return at once if `LEVEL_SET_MODE==1`. `POINT_TO_MESH`. The gas field pointers (`UU`,`VV`,`WW`,`RHOP`,`ZZP`,`PBAR_P`) are aimed by `PREDICTOR`/`CORRECTOR` (about 72-85).
- `CALL_HT_1D` is true when `.NOT.INITIALIZATION_PHASE .AND. CORRECTOR .AND. WALL_COUNTER==WALL_INCREMENT`. Then `DT_BC=T-BC_CLOCK`, `BC_CLOCK=T`, and the 3-D sweep direction cycles (about 88-98). The 1-D solve therefore runs once per `WALL_INCREMENT` steps, in the corrector only.
- Order of the OpenMP loops, each ending in a barrier: L1486 (`WALL_CELL_LOOP_0`), L1487 (thin-wall loop 0), L1488 (`WALL_CELL_LOOP`), L1489 (lateral thin-wall solve), then the CFACE loop (200 on) and the solid-particle loop. CFACE and solid particles are geometry/particle work and are **out of scope** here (stay on the host).
- The barriers matter: L1487 reads records written by L1486, and L1489 reads records written by L1488. The passes must stay separate kernels in the same order.

## 3. Callee tree (T3): pure or inlineable versus stateful

"Pure" means: reads module tables, the gas fields and the caller's own wall records, writes only the caller's own wall records (plus ghost cells, see section 5). "Stateful" means the routine keeps history, writes shared data or has host-only effects.

| Callee (survey lines) | Called from | Class | Notes for the generator |
|---|---|---|---|
| `ASSIGN_GHOST_VALUE` (282-388) | L1486 line 113 | pure, but T4: reads `OMESH(EWC%NOM)` and `MESHES(NOM)%DX/DY/DZ` | does nothing when `EXTERNAL_GHOSTS_FILLED` (patch 0004, default on, required at level > 0). Its loop is L1452 (319-339), a Mesh-side T4 item, not part of this plan |
| `NEAR_SURFACE_GAS_VARIABLES` (402-497) | 116, 129 | pure in the `WALL_INDEX` branch (415-445); the `THIN_WALL` branch (473-495) reads the B1 records of two other walls | `OPTIONAL` dummies and `PRESENT`, three derived-type dummies, `SPACING` (417), `EVALUATE_RAMP` |
| `EVALUATE_RAMP` | 416 and many others | function; pure for plain-table ramps; ramps driven by DEVC/CTRL/external file write `RAMPS()%LAST` | flat ramp tables for plain ramps; the others fall back to the host |
| `HEAT_TRANSFER_COEFFICIENT` (`func.f90` 3135-3321) | 118, 131/133, 742, 762, 781, and inside the 1-D solve (2146) | function, nearly pure: writes `B2%HEAT_TRANSFER_REGIME`, `Z_STAR`, `BLOWING_CORRECTION` (3253-3301) | looks up `MESHES(NM)` (3156), which fails on fine mesh objects (D-056). Two locals are initialised in the declaration (3144) and so are implicitly `SAVE`: they are assigned before use, so they can become plain locals, but the generator must prove it. Contained function `CONSTANT_HTC`. Reads the allocatable `B1%M_DOT_G_PP_ACTUAL` (`ALLOCATED` test at 3287). Calls the four table callees (viscosity, conductivity, specific heat) plus the forced/natural convection models |
| `SURFACE_HEAT_TRANSFER` (603-980) | 154 | per wall, writes only its own B1 and, for `INFLOW_OUTFLOW`, ghost cells (711-718) | `SELECT CASE (SF%THERMAL_BC_INDEX)`; the `CONVECTIVE_FLUX` iteration (776-814) is a loop with `EXIT`; `INTERPOLATED_BC` (821-976) reads `OMESH` (section 7) |
| `SOLID_HEAT_TRANSFER` (1809-3154) | 156, 193 | **stateful**: 1186 lines, 771 pointer sites, owns the 1-D record | section 8.3 |
| `CALCULATE_RHO_D_F` (990-1024) | 160 | pure | species diffusivity table `D_Z` (rank 2), `SIM_MODE` switch |
| `WALL_MODEL` (`turb.f90` 1167-1280) | 165 | pure scalar routine; module constants only; optional outputs are not requested at this call | needs "absent optional argument" support (same feature as `NEAR_SURFACE_GAS_VARIABLES`) |
| `CALC_DEPOSITION` (1645-1737) | 173 | per-wall body plus **shared sums** | `CUNNINGHAM` function, `PBAR`, `PRESSURE_ZONE` cell table |
| `CALC_HVAC_BC` (1745-1796) | 178 | reads HVAC module data (`DUCTNODE`, `DUCT_MF`, `NODE_TMP_EX`, `NODE_ZZ_EX`); writes its own B1 | no device tables exist: host (section 5.4, HVAC) |
| `CALCULATE_ZZ_F` (1035-1460) | 181 | **stateful and shared sums**: 426 lines; updates `Q_IN_SMOOTH`, `QDOTPP_INT`, `K_SUPPRESSION`, ember mass; random draw for `MASS_FLUX_VAR` | section 5.4 |
| `CALCULATE_RHO_F` (1590-1634) | 183 | pure, writes external-only ghost cells (1625-1631) | none beyond table access |

## 4. Data tables

Naming follows `wall_checks.py`/`wall_sums.py`. "WE" is the GPU Wall Loops Engineer, "SL" the AMR Solid Phase Lead.

| Table | Shape | Content | Builder | Owner | State |
|---|---|---|---|---|---|
| `W_*` per-wall flat tables, `BC_IIG/JJG/KKG`, `B1_PRESENT` | (NWALL) | the wall alias fields the loops read | `wall_checks.py` | WE | done / L1239 open |
| `SF_*` per-surface tables incl. `SF_RAMP_INDEX`, `SF_VEL_T`, `W_SURF_INDEX` | (0:N_SURF+reserved) | everything `NEAR_SURFACE_GAS_VARIABLES` reads from `SF` | `wall_checks.py`, probe `test_wall_surf.py` | WE | done for L1486 reads |
| Flat ramp table set | (ramp, point) + per-ramp scalars | plain-table ramps for `EVALUATE_RAMP` | new builder | WE (tables), GPU Generator Engineer (function callee) | new |
| `HTC` extras: `SF` fields for the HTC models, `CP_Z` table | per surface | what `HEAT_TRANSFER_COEFFICIENT` reads | extend `SF_*` | WE | new |
| `THIN_WALL` tables: `TW_WALL_M`, `TW_WALL_P`, `TW_OBST`, `TW_BACK_MESH`, `TW_BACK_INDEX`, `TW_ISOLATED` | (N_THIN_WALL) | thin-wall links | new builder | WE | new |
| `W_LAST_*` ghost-target flag tables (last writer wins) | (NWALL) | one flag per wall and ghost target: "I am the highest wall index that stores here" | `wall_checks.py` last-writer builder | WE | pattern exists (L1402, L1359) |
| `WOBST_START/LIST`, `WCELL_OWNER/START/LIST`, `EW_NBR_OBST` | CSR over walls, per obstruction and per gas cell | wall-to-obstruction and wall-to-gas-cell lists in ascending wall index | `wall_sums.py` | SL | built by a background task, branch `s5-wall-sp`, not yet merged |
| `OD_*` one-dimensional record store (`TMP`, `X`, `RHO_C_S`, `K_S`, `DELTA_TMP`, `HEAT_SOURCE`, layer arrays, `MATL_COMP` arrays) | one record per wall; capacity `N_CELLS_MAX` per record (known at set-up, `type.f90` `BOUNDARY_ONE_D_TYPE`); offsets `OD_OFF(IW)`, counts `OD_N(IW)` | the ragged `ONE_D` allocatables | new builder; layout shared with the record store of `05-fine-level-solid-plan.md` work package 1 | WE (builder), SL (layout review) | new, critical path |
| `B1`/`B2` allocatable species arrays (`ZZ_F`, `ZZ_G`, `RHO_D_F`, `RHO_D_DZDN_F`, `M_DOT_G_PP_*`, `AWM_AEROSOL`, `QDOTPP_INT`, `Q_IN_SMOOTH_INT`, `LP_*`) | (N_SPECIES, NWALL) | rank-2 replacements | new builder | WE | new |
| `EW_NIC`, `EW_NCELL`, `EW_OFF`, overlap ranges | ragged | only if L1485 is kept | T4 plan | WE | not needed if L1485 retires |
| Table-refresh counter | scalar | incremented by `REASSIGN_WALL_CELLS` | upstream patch file (owner commits), host read | FDS Legacy Mapper | proposed in `stage1-generator-answers.md` |

Table checks (host, run after every build and every refresh): `wall_sums.check_all`, `wall_checks.check_all`, plus three new ones that this plan needs: (a) each wall index appears at most once as `WALL_INDEX_M`/`_P` over all thin walls (needed by L1487, section 7.2); (b) ghost-target uniqueness for external walls, and last-writer flags for the rest (section 5.2, row C); (c) record capacity: `OD_N(IW) <= N_CELLS_MAX` of the surface.

## 5. L1488 (144-185): per-wall dispatcher

### 5.1 What it does, line by line

| Lines | Statement | Device note |
|---|---|---|
| 144-148 | loop over all wall cells; `WC`,`BC`,`B1` aliases | wall list (W1 ruling), twin `_wl` entry exists in the generator |
| 149 | `CYCLE` if `NULL_BOUNDARY` | per-wall guard |
| 150-152 | `B2`, `SF` aliases | |
| 153-154 | not thermally thick: `SURFACE_HEAT_TRANSFER` | pass H1 |
| 155-156 | thick and `CALL_HT_1D`: `SOLID_HEAT_TRANSFER(...,SF%HT_DIM*DT_BC,WALL_INDEX=IW)` | pass H2. `HT_DIM>1` (3-D conduction) is rejected in AMR mode (ADR-003, FR-004), so `HT_DIM=1` is asserted |
| 159-161 | `CALCULATE_RHO_D_F` if more than one tracked species and not OPEN/INTERPOLATED | pure |
| 163-167 | `WALL_MODEL` for `SOLID_BOUNDARY` when condensation, deposition or wall output asks for it; writes `B2%U_TAU`, `Y_PLUS` | pure |
| 169-175 | `CALC_DEPOSITION` when deposition, not initialisation, corrector, not solid-phase-only, solid boundary, and (leak path or no node/velocity/volume flow) | shared sums |
| 177-179 | `CALC_HVAC_BC` when `HVAC_SOLVE`, not initialisation, `NODE_INDEX/=0` | host |
| 181 | `CALCULATE_ZZ_F` unless OPEN/INTERPOLATED | stateful, shared sums |
| 183 | `CALCULATE_RHO_F` unless INTERPOLATED | pure, ghost writes |

Serial order inside one wall: heat transfer, `RHO_D_F`, `WALL_MODEL`, deposition, HVAC, `ZZ_F`, `RHO_F`. Every callee reads only values the earlier callees of the same wall wrote, apart from the cross-wall items below.

### 5.2 Writes that race (exact lines) and the rule for each

| # | Write | Lines | Rule |
|---|---|---|---|
| A | `OBSTRUCTION(OBST_INDEX)%MASS`: copy from the neighbour mesh (SET) at 1379, inside `CRITICAL` 1378-1387 (also posts `OMESH(NOM)%REAL_SEND_PKG8`, 1380-1386); subtract at 1390-1391, `CRITICAL` 1389-1392; both corrector only. Third writer: `SOLID_HEAT_TRANSFER` sets `MASS=-1.` for a consumable obstruction at 2695 (branch 2682-2698, not under `CRITICAL`) | | per-wall slots, ordered scan in ascending wall index (SET overwrites the running value, so it is part of the scan; `-1` makes later subtractions start from `-1`). `wall_sums.sum_obst_mass`. No atomics |
| B | `D_SOURCE(IIG,JJG,KKG)` and `M_DOT_PPP(IIG,JJG,KKG,N)`: deposition 1730-1733; thin-obstruction/`LAYER_REMOVED` branch of `ZZ_F` 1402-1423 (`CRITICAL` 1412-1415 species loop; 1420-1422 heat release term) | | per-wall slots, ordered scan per gas cell, wall-major (per wall: deposition species ascending, then `ZZ_F` species terms, then the heat term). `wall_sums.sum_gas_cell` |
| C | Ghost cells `RHOP`,`ZZP`,`RSUM`,`TMP` at `(II,JJ,KK)` and `(II2,JJ2,KK2)`: `INFLOW_OUTFLOW` 711-718; `ZZ_F` 1453-1458; `RHO_F` 1625-1631 | | external walls: one wall per ghost cell by construction, check it (table check b). Internal walls with `INFLOW_OUTFLOW` (open vents, atmospheric profile) write solid cells that several walls can share: gather per target and take the highest wall index (the SP2 rule of `04-blocked-loop-signoff.md`). Do not assert uniqueness |
| D | `B1`/`B2` of the wall itself | | private to the wall; no rule |
| E | `B1%HEAT_TRANS_COEF` of a **back-side wall**, read by the 1-D solve at 2165-2167 (`B1_BACK%TMP_G`, `HEAT_TRANS_COEF`, `Q_RAD_IN`; `B2_BACK%LP_CPUA`) and written in the same loop by `SURFACE_HEAT_TRANSFER` when the back wall is not thick (742, 762, 781) and by `CALC_HVAC_BC` | | upstream is already order-dependent under `SCHEDULE(DYNAMIC)`. See question Q2 |

Particle budgets (`DEPOSIT_PARTICLE_MASS`, 1519-1556) and the 3-D conduction sum (3670-3672) are the same pattern but belong to deferred loops.

### 5.3 Neighbour-mesh reads (T4)

- `CALCULATE_ZZ_F` 1374-1379 reads `MESHES(EWC%NOM)%CELL` and `OBSTRUCTION` to find the obstruction on the other side of an external wall: a **read of the neighbour's value at the start of the loop**. Upstream reads it from a mesh that may be mid-loop; the device version defines it as the value at the start of `WALL_BC` (exchange step before the pass, `EW_NBR_OBST`). Question Q6.
- `HEAT_TRANSFER_COEFFICIENT`, the back-side reads in the solve (1880-1884 wall; 1913-1927 thin wall) and `CALCULATE_ZZ_F` use `MESHES(NM)`. On fine mesh objects these abort by D-056. Question Q3.
- `ASSIGN_GHOST_VALUE` (L1486): skipped in the AMR route. In the compatibility route (`FDSTL_EXTGHOST=0`, level 0 only) it stays on the host or in the Mesh-side T4 package.

### 5.4 Pass plan (the translation unit is the pass, not the loop)

| Pass | Content | Per-wall condition | Shared writes |
|---|---|---|---|
| P1 | L1486 (section 6) | | none |
| P2 | L1487 (section 7) | | none |
| H1 | `SURFACE_HEAT_TRANSFER`, non-thick walls, without the `INTERPOLATED_BC` case | not thick | ghost cells (rule C) |
| H2 | `SOLID_HEAT_TRANSFER`, thick walls, `CALL_HT_1D` only (once per `WALL_INCREMENT` steps) | thick | slot A (the `-1`) |
| G | `RHO_D_F`, `WALL_MODEL`, deposition slots, `ZZ_F`, `RHO_F`, one thread per wall in this order | per the guards above | slots A, B, ghost cells (rule C) |
| S | owner-gather scan: `sum_obst_mass`, `sum_gas_cell` (one thread per owner, serial loop over the CSR segment, no atomics) | | writes `OBSTRUCTION%MASS`, `D_SOURCE`, `M_DOT_PPP` |
| X | `exchange_pkg8` (host) | | `REAL_SEND_PKG8` |
| P3 | L1489 (section 8) | | slot A for thin walls (section 8) |

Splitting by callee is only legal because: (i) cross-wall couplings are exactly the five rows above; (ii) the shared sums are made order-independent of the kernel structure by slots and an ordered scan; (iii) H1 before H2 fixes the back-wall HTC semantics (row E, Q2). The split also means each pass reads values its predecessors finished, which is what a serial run gives per wall.

Pass E and G details that need a decision or special care:

- **H1, `INFLOW_OUTFLOW` (668-721) and the `CONVECTIVE_FLUX` iteration (776-814).** Ghost stores follow rule C; the iteration is a counted loop with `EXIT` and `**ONTH` powers, inlined as written (bitwise discipline: no reassociation, `nofma`).
- **`SURFACE_HEAT_TRANSFER` `INTERPOLATED_BC` case (821-976).** In the AMR route a box interface is not a wall (the driver removes it around the kernels) and a coarse-fine face is handled by the flux hook (D-050, D-061). Pass H1 treats surfaces of this kind as unsupported: a host check aborts the run with a clear message if one is present on a device run. Section 9.
- **`ZZ_F`.** The mass-flux block 1185-1320 runs every stage and uses `T_IGN`, `Q_IN_SMOOTH`, `BURN_DURATION`. The S-Pyro part (1230-1281) uses `SPACING` and ramp tables. `MASS_FLUX_VAR` random draw (1324-1337, `BOX_MULLER`) cannot be bitwise: it uses the counter-based generator of SP-R9 keyed by face and step. Early `RETURN`s when `N_TRACKED_SPECIES==1` (1088-1100) skip the mass consumption: the slot form must carry the same guard. In the fine-level plan this block is run per record (open question OQ-S1 of `05`).
- **HVAC.** `CALC_HVAC_BC` and the HVAC branch of `ZZ_F` (1119-1125) read host-only HVAC tables. First version: if `HVAC_SOLVE` is set, the whole wall pass runs on the host. Q4.
- **Deposition.** Reads `PBAR` and `PRESSURE_ZONE`; `AWM_AEROSOL` updates (1728-1729) are per wall. If `SM%AWM_INDEX` equals `SS%AWM_INDEX` the value is added twice, as upstream does.

### 5.5 Generator features

| Feature | New or existing | Used for |
|---|---|---|
| wall alias set `WC/BC/B1/B2/SF`, wlist twins, `cell_int` tables, `EW_NOM` designator | existing | all wall passes |
| `last_writer` flag tables | existing in the WE branch (L1402, L1359) | rule C |
| `gather=` CSR modes, `unique`, `idempotent`, `private`, `layout="exact"`, `nofma` | existing | S, C |
| Inlining of a subroutine call with `OPTIONAL`/`PRESENT`, derived-type pointer dummies bound to the caller's wall aliases, absent optional arguments | **new** (owner decision: option (b), no upstream edit) | `NEAR_SURFACE_GAS_VARIABLES`, `WALL_MODEL`, `SURFACE_HEAT_TRANSFER`, `CALCULATE_ZZ_F`, `CALC_DEPOSITION` |
| Function callees with a `REAL` result over flat tables | **new** (today only subroutines; four table callees exist) | `EVALUATE_RAMP`, `HEAT_TRANSFER_COEFFICIENT`, `CUNNINGHAM`, `Q_REF_FIT` |
| `SPACING` intrinsic | **new** | 417, 1288 |
| Omit a named statement in a marked range, with a stated reason | **new** (small) | skip line 113 (`ASSIGN_GHOST_VALUE`) in the AMR route |
| Slot store form: a store inside `CRITICAL` becomes a store to `SC_*(IW[,N])` plus an active flag | **new** (requested in the `wall_sums.py` README) | rows A, B |
| Owner-gather kernel template (one thread per owner, serial loop over a CSR segment) | **new** (a variant exists for the gather modes) | pass S |
| Private fixed-size array replacing an `ALLOCATE`d local, with the extent checked against the upstream `ALLOCATE` text | **new** (the `max_extent` policy covers locals with a declared shape only) | the solve (`RHO_DOT`, interpolation weights) |
| Proof that an initialised local (implicit `SAVE`) is assigned before every use, so it can be a plain private | **new** (small) | `HEAT_TRANSFER_COEFFICIENT` |
| Per-wall ragged table kind (`OD_OFF`, `OD_N`, capacity) | **new** | the solve |
| Module `REAL` scalar passed to a callee | **new** (also needed by L1358) | solve, `ZZ_F` |
| Host fallback marking for routines that touch host-only state (HVAC, control-driven ramps) | existing (`host` entries) | HVAC |

### 5.6 Tests for L1488

Common design, section 10. Specific:

- Reference: for each pass, the verbatim survey routine text compiled in a harness with mock module state (the approach of the existing wall probes), run serially. For the whole loop, a capture-and-replay harness: run a small real case on the host to a chosen step, dump the pre-`WALL_BC` state (gas fields, all B1/B2/`ONE_D`, obstructions), run the unmodified loop on the host and the pass sequence on the kernels, and compare every written field bitwise. Cases named in `05-fine-level-solid-plan.md` (`energy_budget_solid`, `surf_mass_vent_char_cart_fuel`, `heat_conduction_a`, `back_wall_test`) plus an atmospheric-profile vent case and a deposition case.
- Order tests (not vacuous): shuffling the wall order in the scan must change at least one `MASS` and one `D_SOURCE` value; the serial-order scan must equal the verbatim loop.
- Mutants: scan order reversed; deposition-all-first instead of wall-major; SET dropped in the scan; burn-away `-1` dropped; ghost last-writer replaced by first-writer; `H1` after `H2`; guard of the early `RETURN` dropped; species index off by one; thin-wall event dropped (section 8).
- Negatives: a build with atomics requested for a `CRITICAL` store must be refused; a ghost target shared by two external walls must fail the table check; a missing table, a missing `PRESENT` binding or a function callee without a flat table must be a generator error.
- Documented difference: cases where the back wall is not thick and has a lower index than the front wall (row E) are not claimed bitwise; the semantics test covers them (Q2).
- `MASS_FLUX_VAR>0`: excluded from the bitwise list; the SP-R9 statistical test (mean, variance, clipping bounds, stream independence across thread counts) applies.

Effort (days): H1 body 4-6; G passes without `ZZ_F` 5.5-6.5; `ZZ_F` 8-12; solve see section 8; assembly and replay harness 7-10 in total (3-4 assembly, 4-6 harness).

## 6. L1486 (109-119) and the near-surface/HTC callees

### 6.1 Line by line

| Lines | Statement |
|---|---|
| 108-109 | `!$OMP DO SCHEDULE(DYNAMIC)`, loop over all wall cells |
| 110-112 | aliases `WC`,`BC`,`B1` |
| 113 | external walls only (`IW<=N_EXTERNAL_WALL_CELLS`): `ASSIGN_GHOST_VALUE(IW,BC,B1)`. Runs **before** the null-boundary test |
| 114 | `CYCLE` if `NULL_BOUNDARY` |
| 115 | `SF` alias |
| 116 | `NEAR_SURFACE_GAS_VARIABLES(T,SF,BC,B1,WALL_INDEX=IW)`: ramp factor, tangential velocity at the gas cell, `TMP_G` (front temperature ramp or cell value), the optional radiation override, `RHO_G`, `ZZ_G` |
| 117-118 | if `CALL_HT_1D` and thermally thick: `B1%HEAT_TRANS_COEF = HEAT_TRANSFER_COEFFICIENT(NM,T,TMP_G-TMP_F,SF,WALL_INDEX_IN=IW)` |

### 6.2 Dependencies, races, T4

- Callees as in section 3. No cross-wall write: `B1` is the wall's own; `HEAT_TRANSFER_COEFFICIENT` writes the wall's own `B2`. No shared sum.
- Ghost writes of `ASSIGN_GHOST_VALUE` are at (II,JJ,KK) and the second layer, one wall per ghost cell (unique); the skip rule makes them irrelevant in the AMR route.
- T4 reads: only inside `ASSIGN_GHOST_VALUE` (skipped) and the `MESHES(NM)` lookup in `HEAT_TRANSFER_COEFFICIENT` (question Q3).

### 6.3 Plan

1. The WE already has the per-surface tables and a probe (`test_wall_surf.py`). The kernel entry waits for the generator inliner (option (b)); the alias binding side stays with the WE.
2. Omit line 113 in the AMR route (feature "omit statement"). In the compatibility route the existing host path is kept.
3. `EVALUATE_RAMP` as a function callee over the flat ramp table; ramps with DEVC/CTRL/external-file drivers fall back to the host (the wall list excludes surfaces that use them, or the whole pass runs on the host).
4. `HEAT_TRANSFER_COEFFICIENT` as a function callee: needs the fine-mesh lookup resolved (Q3), the `SAVE` proof, the allocatable `M_DOT_G_PP_ACTUAL` test replaced by a table flag, and the `B2` stores as flat stores. Sub-tests per `HTC_MODEL` branch (default, log-law, Rayleigh, forced/natural/impact, fixed value, ramped value, blowing correction) and for the back-side call (no `WALL_INDEX`, 3192-3198).

Effort: L1486 entry and tests 3-4 days after the inliner and the ramp function exist; HTC 3-4 days. Owner: WE (tables, entry, tests), GPU Generator Engineer (features).

## 7. L1487 (124-137): thin-wall near-surface pass

### 7.1 Line by line

Runs only if `N_THIN_WALL_CELLS>0` and `CALL_HT_1D`. Per thin wall `ITW`: alias `TW`, `B1`,`BC` of the thin wall; `NEAR_SURFACE_GAS_VARIABLES(...,TW=TW,THIN_WALL_INDEX=ITW)` (129). In the `THIN_WALL` branch (473-495) the gas temperature and incoming radiation of the thin wall are the averages of the two side walls' `B1M%TMP_G` and `B1P%TMP_G`, `Q_RAD_IN` (473-487). Then (130-136): HTC with `WALL_INDEX_IN=TW%WALL_INDEX_M`; else with `_P`; else HTC=0.

### 7.2 Dependencies and races

- It **reads values L1486 wrote** (the side walls' `TMP_G`, `Q_RAD_IN`), so it must be a separate pass after P1, and this is why it cannot be folded into the L1486 kernel.
- It **writes the side wall's `B2`**: `HEAT_TRANSFER_COEFFICIENT` stores `HEAT_TRANSFER_REGIME`, `Z_STAR`, `BLOWING_CORRECTION` into the record of wall `WALL_INDEX_M/_P` (3253-3301), overwriting what P1 wrote there. Race only if one wall index is the side wall of two thin walls. Table check (a) of section 4 makes this explicit; if it fails, the rule is last writer wins by thin-wall index. These three values are diagnostics (output), they do not feed the solve.
- T4: none inside the mesh; with the box view the thin-wall links stay within one box.

### 7.3 Plan and effort

Same features as L1486 plus the `TW` alias set and the thin-wall tables (section 4). Effort 1.5-2 days after L1486. Owner: WE.

## 8. L1489 (192-194): lateral thin-wall solve

### 8.1 What it is

`IF (CALL_HT_1D)`: `DO ITW=1,N_THIN_WALL_CELLS; CALL SOLID_HEAT_TRANSFER(NM,T,3._EB*DT_BC,THIN_WALL_INDEX=ITW)`. The factor 3 is the lateral sweep count. The loop body is one call; it calls neither `CALC_DEPOSITION` nor `CALC_HVAC_BC` and has no gas-cell sums of its own (those entries of the register are wrong, section 17).

### 8.2 Dependencies, races, T4

- Callee: the same `SOLID_HEAT_TRANSFER` as the thick-wall pass, thin-wall branch: back-side thin wall by `MESHES(BACK_MESH)%THIN_WALL(BACK_INDEX)` (1913-1927), reads of `B1_BACK%TMP_G`, `HEAT_TRANS_COEF`, `Q_RAD_IN` of the two back walls (2165-2183). Those values were finished by P1/P2/H1 before this loop starts (barrier), so reading them here is deterministic. The cross-record read is T4-like but within the same box.
- Shared writes: none to gas cells. **One shared event**: a consumable thin obstruction that has burned through sets `OBSTRUCTION(OBST_INDEX)%MASS=-1.` at 2695 via `TW%OBST_INDEX`. In serial order this comes **after** all wall events of L1488 and in ascending thin-wall index. `wall_sums.py` covers only the wall events today. Extension (SL, 0.5-1 day): thin-wall `A` events appended after the wall scan, with a table `TW_OBST` and a test that a burned-through thin wall sets `-1` after a wall subtraction on the same obstruction.
- `HT_DIM>1` returns early (1874, 1909); not reachable in AMR mode (asserted).

### 8.3 The 1-D solve port (the largest piece)

`SOLID_HEAT_TRANSFER` is a stateful routine. Features of the body that decide the plan:

- State is the wall's own `ONE_D` record: ragged allocatables (`BOUNDARY_ONE_D_TYPE`: `X`, `TMP`, `RHO_C_S`, `K_S`, layer arrays, `MATL_COMP(N)%RHO(:)` etc.). Layer removal and renoding (2361-2804) and `REMESH` change the number of cells within the capacity `N_CELLS_MAX`; no reallocation is needed if the record store has that capacity (section 4, `OD_*`).
- A sub-step loop with an adaptive sub-step (2306-2337) and a global `ICYC`; the tridiagonal solve (2856-2908); `PERFORM_PYROLYSIS` (contained, 2993-3152) and `PYROLYSIS` (3195-3605).
- Host-only effects: `SHUTDOWN`/`WRITE` on the error path (2669-2670), `ALLOCATE` of locals (for example `RHO_DOT` at 2016). The port needs a device error flag (per wall status code the host reads after the pass) and the private fixed-size array feature.
- The solve reads gas-side values through `B1` and the back-side through section 5.2, row E. It never writes a back-side record.

Stages (all verified against the unchanged routine on captured records, bitwise):

1. Leaf callees with no state: the helper routines the solve calls (node weights, interpolation, property updates). 5 days.
2. The record store, capacity and the host builder with checks (WE, shared layout with the fine-level plan). 4-6 days.
3. The solve for one layer, one material, no pyrolysis; then multi-layer; then pyrolysis; then layer removal/renoding. 8-12 days with test growth at each step.
4. Thin-wall branch and the back-side reads. 2-3 days.
5. Burn-away events into slot A. 1 day.

Owner: SL for the port and physics review, executed by background tasks; WE for the record store; GPU Generator Engineer for the private-array and ragged-table features. The generator rewrites the text (alias to flat, `ALLOCATE` to private); the routine is not rewritten by hand. A hand port is the fallback if the rewrite cannot cope with the 771 pointer sites, and then the verbatim comparison is by the capture-and-replay harness only.

Effort total for the solve incl. thin-wall branch: 18-28 days (low confidence). L1489 itself, as a pass over the thin-wall list with the extra event: 1 day after the solve.

## 9. L1485 (888-947): `INTERPOLATED_BC` coarse-side species flux

### 9.1 What it does

Inside `CASE (INTERPOLATED_BC)` of `SURFACE_HEAT_TRANSFER`, in the branch for a coarse wall facing several finer cells (`EWC%NIC>1`), with more than one scalar: for each species `N` it looks up `D_Z`; for every overlapping cell of the neighbour mesh (loops over `KKO`,`JJO`,`IIO`, in that order) it computes the neighbour's `RHO_D_DZDN` from `OM%RHO`, `OM%ZZ`, `OM%MU`, the neighbour's cell sizes and centre ramp (`EVALUATE_RAMP(MM%ZC(KKO),I_RAMP_P0_Z)`), and a viscosity call for LES, and sums `ARO*RHO_D_DZDN` into the wall's `B1%RHO_D_DZDN_F(N)`; then a correction (950-953) makes the species sum consistent (`MAXLOC`/`SUM`). The result is the wall-side diffusive flux that `DIVERGENCE_PART_1` moves between the face arrays and the wall record (divg lines 213-226).

### 9.2 Why retire it

- It exists for FDS meshes of different resolution. In the AMR route a same-level box interface is not a wall cell (the driver removes it around the kernels, `level-interface.md`), and a coarse-fine face gets its flux by the interface flux overwrite (D-050, D-061). Patch 0004's note says the same: this branch "needs different-level neighbours, which M2a does not have".
- Highest churn in the wall file (about 1176 changed lines over 15 commits in the survey, `gpu_routine_churn`), reads `OMESH` and `MESHES(NOM)`, needs ragged overlap tables, and its modelled share is not a device cost in the AMR route.

### 9.3 Plan

- Default: no kernel. Table check at build: no wall with `EWC%NIC>1` on a device run; the host aborts with a clear message. 0.5 day (SL with the Legacy Mapper).
- If the Architect says `NIC>1` can occur in the supported inputs: T4 plan with `EW_NIC/EW_NCELL/EW_OFF`, ordered K,J,I sums (bitwise because the order is fixed), ramp and viscosity callees. 4-5 days (WE tables, SL test design). Q1.
- Test if kept: bitwise against the verbatim loop on random overlap geometries with `ARO` below and at 1, 1-4 species, LES on/off; mutants: loop order changed, `MAXLOC` correction dropped.

## 10. Test design shared by all passes

- **Reference:** verbatim survey text (line ranges cited, compiled from the pinned commit), upstream `!$OMP` lines removed, run serially. For callees the verbatim routine is compiled into the harness with mock module state.
- **Matrix:** six flag sets (`O0`, `O2`, `O0omp`, `O2omp`, `O2omp_off`, `O2omp_dpd`); both callee switches (`dpd`, `bind`); 1, 4 and 8 threads for the device-style kernels run on the host; results bitwise equal in all cells, also across thread counts.
- **Inputs:** random box meshes with thick and thin obstructions, several walls per obstruction and per gas cell, external walls with neighbour obstructions, 1-4 species, predictor and corrector, `CALL_HT_1D` on and off, burn-away on and off, deposition on and off.
- **Mutants:** at least one per rule in section 5.2, plus operand, sign and index mutants of each kernel. All must be caught.
- **Negatives:** generator refusals (section 5.6) and table-check failures.
- **Non-vacuous order checks:** shuffled order must differ from serial on at least one case.
- **Excluded from bitwise:** `MASS_FLUX_VAR>0` (SP-R9), cases of row E with a lower-index non-thick back wall.
- **Run discipline:** generator runs and test chains under `flock (local project directory)/src-s5gen/.s5gen.lock`, at most 4 cores.

## 11. L0823 (`init.f90` 4940-5054) and L0822 (4917-4936): `REASSIGN_WALL_CELLS`

### 11.1 What it does

`OBSTRUCTION_LOOP: DO OBST_INDEX=1,N_OBST` skips every obstruction that is not `SCHEDULED_FOR_REMOVAL` or `SCHEDULED_FOR_CREATION`. For the rest it walks the obstruction faces in x, y, z (4958-5023) and calls the contained `GET_BOUNDARY_TYPE` (5061-5204) for each newly exposed or covered wall cell, then `REDEFINE_EDGE` (5213-5226; 5027-5052 computes an edge index and changes nothing else, so it is effectively a no-op). L0822 is the preceding loop over external wall cells with no obstruction (4917-4936), the neighbour-mesh part.

`GET_BOUNDARY_TYPE` re-types the wall cell (`BOUNDARY_TYPE`, `SURF_INDEX`, `OBST_INDEX`), flips `CELL(IC)%SOLID`, sets `B1%T_IGN`, calls `INIT_WALL_CELL` (`init.f90` 2975-3386, allocates the `ONE_D` and property records), calls `SET_DENSITY_AND_MASS_FRACTIONS_AT_WALL`, copies whole `BOUNDARY_ONE_D` derived types between meshes (5188) for burn-away (the exposed surface inherits layers and ignition time), and rewrites `OMESH(...)%WALL_SEND_BUFFER%ITEM_INDEX` (5194-5199).

It is called only from `CREATE_OR_REMOVE_OBSTRUCTIONS` (`main.f90` 1782-1798) when `OBST_CREATED_OR_REMOVED`, after `OPEN_AND_CLOSE` (`init.f90` 4507), `EXCHANGE_GEOMETRY_INFO` and before `GLOBAL_MATRIX_REASSIGN`. It is event-driven (burn-away, control-driven creation/removal); obstruction change under AMR is deferred (FR-042).

### 11.2 Callee tree, races, T4

Impure throughout: allocation, I/O, `SHUTDOWN`, string handling, module global writes, 13 allocatable-pointer sites. Shared writes: the same wall cell can be reached from two obstructions; neighbour-mesh reads of `MESHES(NOM)%WALL/CELL`. None of this is expressible as a wall kernel, and it does not need to be.

### 11.3 Plan

- Keep L0823 and L0822 on the host and reclassify them "host-side set-up/event code" in the register (the projection note `generator-projection.md` already leaves set-up routines out; the modelled 1.73 % per step is not credible because the code runs only on an event).
- The device-relevant deliverable is a **refresh trigger**: after `REASSIGN_WALL_CELLS` the driver rebuilds `WALL_INDEX`, the `W_*` and `SF_*` tables, the thin-wall tables, the `WOBST`/`WCELL` CSRs, the exterior and solid masks and the `OD_*` store; then runs the table checks of section 4. The trigger is a counter incremented by `REASSIGN_WALL_CELLS` that the shim compares every step (`stage1-generator-answers.md`); that is an upstream change, so a patch file for the owner (D-051). Also refresh when the `WALL` array is reallocated (blocks of 1000, `func.f90` 4003-4014).
- Burn-away interplay: the `-1` of row A marks the obstruction for removal (`init.f90` 4542 tests `MASS<TWENTY_EPSILON_EB`); restart (`dump.f90` 3944) and the output of gas-phase residuals (`dump.f90` 9148) read it. The slot scan must therefore finish before those host reads.
- Test: a refresh test with an obstruction removed and one created between two steps: tables after refresh equal tables built from scratch; checks pass; a missed refresh (counter ignored) is a mutant that must be caught by the table-check or by the replay compare.
- Effort: patch file and driver hook 1-2 days (FDS Legacy Mapper with the integration role). A device version of the routine would need allocation and derived-type copies the generator cannot express and gives no benefit; not planned.

## 12. Neighbours that stay as they are

- `ASSIGN_GHOST_VALUE` loop L1452 (319-339): T4, skipped when `EXTERNAL_GHOSTS_FILLED`; only the compatibility route needs it.
- HT3D sweeps L1470/L1471 (3635-3677, 3681-3717): out of scope (3-D conduction is not supported in AMR mode).
- CFACE and solid-particle loops of `WALL_BC`: geometry and particle work, deferred.

## 13. Work split

| Role | Takes | Days (low-high) |
|---|---|---|
| AMR Solid Phase Lead | this plan; `wall_sums.py` and its extension for thin-wall events (0.5-1); the capture-and-replay harness (4-6); pure gas-side passes (`RHO_D_F`, `WALL_MODEL`, `RHO_F`, deposition, HVAC guard: 5.5-6.5); `CALCULATE_ZZ_F` (8-12) incl. the SP-R9 test; the `SURFACE_HEAT_TRANSFER` body (4-6); the 1-D solve port incl. thin-wall branch and L1489 (18-28); L1485 guard (0.5, or 4-5); physics review of every pass; the consumer loops L0394/L0405, to be reconciled with the WE copies; rulings and tests for Q2, Q5 | 43-62 (+5 if L1485 is kept) |
| GPU Wall Loops Engineer | the tables of section 4 (8-12); L1486 entry and tests (3-4); HTC callee (3-4); L1487 (1.5-2); pass assembly with the driver (3-4); L1239; the CSR family and `wall_checks.py` builders; alias binding for the inliner | 18-26 |
| GPU Generator Engineer | inliner for `OPTIONAL`/derived-type pointer dummies (4-6); function callees (2-3); `SPACING` (0.25-0.5); omit-statement (0.5); slot store form and owner-gather template (3-4); private fixed-size arrays (1-2); `SAVE` proof (0.5) | 11-17 |
| FDS Legacy Mapper | claims for L1486-L1489 and the pass split in `loop_claims.csv`; register corrections (section 14); reclassification of L0823/L0822 and the L1485 verdict; the REASSIGN counter patch file and refresh hook with the integration role (1-2); survey versus working-tree line verification; the merge of the WE branch into the generator branch; the S2 package stays as listed in the work list | 1-2 plus coordination |

## 14. Dependency order

1. Merge of the WE branch into the generator branch (freezes lift); table checks (a)-(c); CSR family; `wall_sums.py` merged with its tests.
2. Small generator features: `SPACING`, omit-statement, `SAVE` proof, module `REAL` scalar.
3. Inliner and function callees. Then L1486 with the HTC callee, then L1487.
4. In parallel from step 1: the record store and the private-array feature, then the solve stages 1-5 of section 8.3.
5. Slot store form and owner-gather template, then deposition and `ZZ_F` with the scan; `RHO_D_F`, `WALL_MODEL`, `RHO_F`; `SURFACE_HEAT_TRANSFER` body.
6. Pass assembly (section 5.4) and the replay harness on the named cases; then L1489 with the thin-wall event.
7. L1485 only if Q1 says so. L0823 hook any time (independent).

Interim state worth stating: with the gas-side passes on the device and the solve still on the host, the `B1` fields of every wall cross the host-device boundary once per `WALL_INCREMENT` steps. Correct and testable; the cost is for the GPU plan to judge.

## 15. Top risks

1. The solve port (1186 lines, ragged records, adaptive sub-steps, host-only error path) dominates the effort and its estimate is low-confidence.
2. Generator features (inliner, function callees) are on the critical path of L1486/L1487/L1488 and are not yet designed; the freeze until the branch merge adds waiting time.
3. Order-dependent shared state (obstruction mass with SET, `-1`, subtraction, and the thin-wall event; gas-cell sums) must follow wall-major ascending order; a missed event breaks mass conservation silently.
4. Bitwise claims have limits: back-wall HTC (row E), neighbour-mass read, `MASS_FLUX_VAR`. The plan states these as documented differences; the V&V Lead must accept them.
5. `SURFACE_HEAT_TRANSFER` and `ZZ_F` have high upstream churn; each rebase costs a regeneration plus the port-merge check.
6. Fine-level assumptions: kernels written for "one record per wall" must accept a record list (`05`), and `MESHES(NM)` lookups abort on fine objects.
7. Modelled shares of L0823 and L1485 overstate their device value; counting them as translated work would inflate coverage.

## 16. Open questions

| # | Question | Owner |
|---|---|---|
| Q1 | Can `EWC%NIC>1` occur in any supported AMR input? If not, retire L1485 (section 9). | AMR Chief Architect |
| Q2 | Back-side HTC semantics for the 1-D solve (row E): snapshot at the start of the wall loop (recommended; matches the per-step exchange of FR-046), or "non-thick walls first" (H1 before H2). Upstream is order-dependent. Needs a V&V acceptance for the documented-difference cases. | AMR Chief Architect with AMR V&V Lead |
| Q3 | How do fine mesh objects reach `HEAT_TRANSFER_COEFFICIENT` and the back-side reads (box-view hook or an upstream patch under D-051)? (OQ-S3 of `05`) | AMR Chief Architect with FDS Legacy Mapper |
| Q4 | Is "whole wall pass on the host when `HVAC_SOLVE`" acceptable for the first version? | AMR Chief Architect |
| Q5 | May the mass-flux block of `CALCULATE_ZZ_F` be extracted verbatim as a record-local routine? (OQ-S1 of `05`) | AMR Chief Architect |
| Q6 | The neighbour obstruction mass read (1379): is "value at the start of `WALL_BC`" the contract? | AMR Chief Architect |
| Q7 | Which inliner design for `OPTIONAL`/derived-type dummies and function callees is feasible, with what dates? Can the slot store form be one feature for all three `CRITICAL` patterns? | GPU Generator Engineer |
| Q8 | What do `EW_NCELL`/`EW_OFF` (T4 plan) mean exactly, and does the ragged `ONE_D` store use capacity rows or a prefix-sum layout? Do the WE versions of L0394/L0405 replace the consumer claims of the AMR Solid Phase Lead? | GPU Wall Loops Engineer |
| Q9 | Reclassify L0823/L0822 as host-side; fix the L1489 register text; write the REASSIGN counter patch file; confirm which of L1486-L1489 claims are open. | FDS Legacy Mapper |
| Q10 | Statistical test design for the SP-R9 generator at `MASS_FLUX_VAR`; extension of `wall_sums.py` with thin-wall events (section 8.2). | AMR Solid Phase Lead |

## 17. Corrections to other documents (to be made by their owners)

- `loop_work_list.csv` / `loop-work-list.md`: the L1489 row describes a loop calling `CALC_DEPOSITION`/`CALC_HVAC_BC` with gas-cell sums. At 192-194 it is the lateral thin-wall solve (section 8). Those calls are in L1488 (173, 178). The AMR Solid Phase Lead's own sign-off needs the matching correction.
- `04-blocked-loop-signoff.md` (own document, to be fixed in a separate commit): several `wall.f90` cites are off against the survey. Survey values: obstruction mass copy 1379 and subtraction 1390-1391 (`CRITICAL` blocks 1378-1387 and 1389-1392); `D_SOURCE`/`M_DOT_PPP` 1412-1415 and 1420-1422 (deposition 1730-1733); particle block 1519-1556; back-wall pointers 1880-1884 and 1913-1927; third mass writer 2695.
- `wall_sums.py` and its README quote working-tree line numbers (9-12 higher than the survey); a note on the revision should be added by the author of the files.

## 18. Rulings received

Dated 2026-10-04. Sources: decisions D-064 (fine-level solid plan, `docs/solid/05-fine-level-solid-plan.md`, questions OQ-S1 to OQ-S6) and D-065 (this plan), in `README.md` and in `adr/drafts/ruling-role3-plan.md`, Updates (f) and (g). Answers from the Wall Loops Engineer and the Generator Engineer on Q7 and Q8 are logged here as well. Sections 1 to 17 above are left as written (v0.1); where a ruling changes them, the change is stated here.

### 18.1 D-064 (fine-level solid plan)

| Plan label | Question | Ruling | Consequence for the translation |
|---|---|---|---|
| OQ-S1 (this plan: Q5) | May the mass-flux block of `CALCULATE_ZZ_F` (`wall.f90` 1185-1320) be extracted verbatim as a record-local routine? | Yes, as a guarded `WITH_AMREX` patch, bit-identical in FDS-only mode. Status **DRAFT** until the oneAPI and GNU Debug builds are validated. The Legacy Mapper confirms the block has no dependence on owner loop state beyond its arguments. A6 is extended to cover it. | The `ZZ_F` work (section 5.4, 8 to 12 days) is planned as two parts: the record-local mass-flux routine and the owner remainder. No generator work on the routine starts before the Mapper's confirmation. |
| OQ-S2 | Face key | (level, global integer index of the gas-side cell at that level's resolution, `IOR` in plus or minus 1 to 3). Independent of mesh, box and rank. Sum order (level, k, j, i, `IOR`). The Spec Lead aligns FR-046. | Replaces the Morton-index proposal of `05` section 2.1. Record store sort order and the ordered sums use this key. |
| OQ-S3 (this plan: Q3) | How do fine mesh objects reach `HEAT_TRANSFER_COEFFICIENT` and the back-side reads? | One accessor returns the mesh object of (NM, level) via `POINT_TO_BOX` (D-056). Introduced by a guarded local patch at `HEAT_TRANSFER_COEFFICIENT` (`func.f90` 3156) and at the back-side `MESHES(NM)` lookups. Not an upstream patch. The Legacy Mapper lists every `MESHES(NM)` lookup on the solid path. | The HTC callee (section 6.3, item 4) and the back-side reads of the solve use the accessor; the sites are the ones named in sections 3 and 5.3. Open: the Mapper's list. |
| OQ-S4 | Who builds fine-level OBST wall tables | Role 1 (Data Layout), in Phase 5. Role 3 triggers the rebuild at regrid. | No work for this plan; the refresh trigger of section 11.3 is the same hook. |
| OQ-S5 | Area-mean `TMP_F` against the FR-022 wall budget | Deferred. Bring the WP6 numbers if the budget misses. | None. |
| OQ-S8 (D-064 item 6) | Thin-wall faces | Records are never deeper than their owner; thin-wall records stay at the owner's level. | The lateral thin-wall pass (L1489, section 8) runs at owner level; no fan-out for it. |

D-064 numbers its six items 1 to 6; item 6 answers OQ-S8 of `05`. The radiation-history question (`ILW`, `Q_RAD_IN` belongs to the record or to the owner, OQ-S6 in `05` section 8) and `AREA_ADJUST` (OQ-S7) were not part of D-064 and stay with the Radiation Lead and the Spec Lead.

### 18.2 D-065 (this plan)

| Question | Ruling | Condition and consequence |
|---|---|---|
| Q1 | Retire L1485 (`EWC%NIC>1`) in the AMR route behind a host abort guard. | **Conditional**: the Legacy Mapper and the Solid Phase Lead grep the supported AMR input set and confirm no case has `NIC>1`; inputs that fail are listed as FDS-only and refused (consistent with D-055). Effort stays 0.5 day (section 9.3); the 4 to 5 day T4 plan is not started. |
| Q2 | Back-wall heat transfer coefficient (row E of section 5.2): snapshot of the other side in the AMR route. | FDS-only mode is untouched. The host AMR route offers the same snapshot mode, so the difference to FDS is measured as algorithm, not as port. The V&V Lead signs the tolerances for `back_wall_test` and `heat_conduction_a`. The "non-thick walls first" alternative is dropped. |
| Q4 | Whole wall pass on the host when `HVAC_SOLVE` is on. | Acceptable in Phase 4. HVAC cases are correctness gates only, not GPU performance gates. Listed as a known limitation; reopened if an HVAC case becomes a performance target. |
| Q6 | Neighbour obstruction-mass read (line 1379). | Returns the value at the start of `WALL_BC` of that stage: snapshot before the pass, read-only during it, identical on host and device. The difference from FDS order is recorded like Q2. Burn-away stays deferred (FR-042). |

### 18.3 Answers from other roles

- **Q7 (Generator Engineer), partly answered.** Function callees with a `REAL` result: 2 to 3 days after the current queue. `SPACING` intrinsic: hours. Call-site inlining for `OPTIONAL` and derived-type dummies: 3 to 4 days. The slot-store feature for the three `OMP CRITICAL` patterns (rows A and B of section 5.2) is **not in the generator queue yet**, and its combine order must replay the serial order (ascending wall index, wall-major). Still open: queue position and dates for the slot-store feature and the owner-gather template.
- **Q8 (Wall Loops Engineer), answered.** `EW_NCELL` and `EW_OFF` do not exist. The CSR tables are prefix-sum tables built by `wall_checks.csr_cell_walls`: `W_CSR_GAS_PTR`/`W_CSR_GAS_LIST` (key gas cell) and `W_CSR_WC_PTR`/`W_CSR_WC_LIST` (key wall cell). No ragged `ONE_D` store exists yet; it is built on the same prefix-sum pattern (`OD_OFF`, `OD_N`). L0394 and L0405 were already translated by the Wall engineer (`wall_bc_dp`, `wall_spec_adv2`), so the consumer claims of the Solid Phase Lead on them are **released** (they are reviewed in `docs/solid/07-sp2-sp3-kernel-review.md`). The consumer line in the work split (section 13) no longer applies; the Solid Phase Lead takes L1486 and L1487 instead of L0394 and L0405.

### 18.4 Question status after the rulings

| # | Status | What remains |
|---|---|---|
| Q1 | **Closed, conditional** (D-065) | The input-set grep (Mapper with the Solid Phase Lead). Result decides between "refused inputs listed" and a reopen. |
| Q2 | **Closed** (D-065) | V&V tolerances for `back_wall_test` and `heat_conduction_a`; the host snapshot mode (a switch in the host AMR route). |
| Q3 | **Closed** (D-064 (3)) | The Mapper's list of `MESHES(NM)` lookups on the solid path; the guarded patch (draft). |
| Q4 | **Closed** (D-065) | Known-limitation entry. |
| Q5 | **Closed, DRAFT** (D-064 (1)) | oneAPI and GNU Debug validation; the Mapper's no-hidden-dependence confirmation. |
| Q6 | **Closed** (D-065) | Record the difference from FDS order next to Q2's. |
| Q7 | **Open, partly answered** | Dates for the slot-store feature and the owner-gather template; function callees 2 to 3 days after the current queue; inliner 3 to 4 days; `SPACING` hours. |
| Q8 | **Closed** | Build the ragged `OD_*` store on the prefix-sum pattern (Wall engineer for the builder; layout review in `docs/solid/06-solve-port-test-design.md`). |
| Q9 | **Pending** | FDS Legacy Mapper (reclassify L0823 and L0822 as host-side, fix the L1489 register text, the REASSIGN counter patch file, claim confirmation). |
| Q10 | **Open, mine** | SP-R9 statistical test design; thin-wall events in `wall_sums.py`. |
