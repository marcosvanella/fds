# 05 · Fine-level solid-phase plan: 1-D wall conduction under refinement

Owner: AMR Solid Phase Lead · Status: **v0.1, plan for review**, 2026-10-03 · Reference tree: FireX `36975d765f`; all `file:line` below are from that tree (the working tree has moved).
Tags: **[R]** ruled, with the ID. **[P]** proposal, not decided. Labels `OQ-S1`.. are local to this file, not decision IDs.
Inputs: spec `01` v0.1.2, `02`, `03`, `04`; `adr/drafts/ruling-FR041b-G2.md` (A1-A6); ADR-003; ADR-001 "Wall tables and wall-state ownership" (W1, W2); README decisions D-010, D-033, D-038, D-045, D-047, D-050, D-051, D-053, D-056, D-058, D-061; FR-005, FR-041a/b, FR-045..047; `role3-regrid-transport-plan.md`; `next-work.md`; the GPU wall report.

## At a glance
- Option A is the data design [R]: flat, face-keyed 1-D records at the finest-ever depth; the coarse owner drives its `r²` records with its own gas-side inputs and receives their fluxes in a fixed order.
- `SOLID_HEAT_TRANSFER` stays untouched. One more routine must be split to make that work: the mass-flux block of `CALCULATE_ZZ_F` reads and updates record state every stage [P, OQ-S1].
- Phase 5 is the same design with depth = owner level (one record per face, no split). Splitting and coarsening are FR-041b, after Phase 6.
- One global dt and one driver-owned clock; records carry no clock.
- Gas-side wall work outweighs the 1-D solve in three of four G2a cases, so it is the first GPU target.
- Eight questions need rulings (§8); OBST handling and burn-away are deferred (§9).

## 1. Scope and what stays unchanged
- **In scope:** the wall-side solid phase of a face whose owner sits on any level: `WALL_BC` (wall.f90:26-279), its near-surface gas pass (:108-119, `NEAR_SURFACE_GAS_VARIABLES` :402), the heat transfer coefficient (`HEAT_TRANSFER_COEFFICIENT`, func.f90:3135), `SURFACE_HEAT_TRANSFER` (:603) for non-thick surfaces, the 1-D solve (`SOLID_HEAT_TRANSFER`, :1809-2992, called :156), the mass-flux part of `CALCULATE_ZZ_F` (:1185-1320, called :181), aggregation to the owner, budgets, devices, restart.
- **Unchanged [R, A6 / FR-047 (a)]:** the body of `SOLID_HEAT_TRANSFER`, including adaptive sub-steps (:2306-2336), pyrolysis (:2993, :3195), renoding and layer removal (:2361-2804), the tridiagonal (:2856-2895). `SURFACE_HEAT_TRANSFER` and the FDS formulas of the near-surface pass and of the HTC are unchanged too [P]. A record sees the same kind of inputs as an FDS wall cell; only where they come from changes.
- **Changed:** the driver loop (:143-194 becomes a loop over owner faces plus a loop over records), the wall-table build for fine levels, fan-out to and fan-in from records, the back-side lookup (FR-046), regrid and restart handling.
- **Phasing [R, D-010, roadmap]:** Phase 5 = static hierarchy, FR-041a freeze, exactly one record per face, on the owner level (no split ever happens). FR-041b = after Phase 6: split on first refinement, records outlive owner changes. The data design below is the same in both; Phase 5 uses depth = owner level.
- **Deferred:** OBST creation, removal and burn-away across levels (FR-042), and the items in §9.

## 2. Data
**2.1 Record store [R: option A, A3, A5; P: layout].** Records sit in a flat structure-of-arrays outside the box layout (NFR-048), one row per record, padded to `N_CELLS_MAX` (type.f90:220) per `SURF` class.
- Row content: the solid state of `BOUNDARY_ONE_D` (type.f90:217-264: node depths, `TMP`, per-material densities, layer thicknesses, counters, `PART_MASS`), the solid and ignition part of `BOUNDARY_PROP1` (:292-345: `TMP_F`, `T_IGN`, `Q_IN_SMOOTH(_INT)`, `QDOTPP_INT`, `BURNAWAY`, held `Q_CON_F` and mass fluxes) and of `BOUNDARY_PROP2` (:350-376), plus the frozen setup values: thickness (D-038, A1), `AREA` (record area), `AREA_ADJUST`.
- Not stored per record: coordinates and gas-cell indices; they follow from the face key and the owner.
- Non-thermally-thick surfaces have no 1-D record. They keep one `B1` per owner face, computed by `SURFACE_HEAT_TRANSFER` at the owner's resolution, and are rebuilt after a regrid (FR-041a).
- The exact field split (state / frozen / gas-side input / owner output) is work package WP1. It must cover every accumulator (R-60); a field not on the list is a bug.
- Face key [P]: level-0 column indices of the OBST face, `IOR`, and a Morton sub-index at the record's depth. A group of records is a contiguous key range. The definition is shared with FR-046 (OQ-S2).

**2.2 Face to record(s) [R: ruling §1; P: mapping].** The owner of a face is the level whose valid, uncovered cell is its gas side (FR-045; divg.f90:760). If the finest-ever depth equals the owner level there is one record, identical to an FDS wall cell. Otherwise the owner drives a group of `(r_total)²` records, `r_total` being the ratio from the owner to the finest-ever level (R=2 or 4 per jump, compounded). Fine boxes cover whole coarse cells, so a coarse face is either owned whole or covered whole (assert at build). The mapping is a CSR table owner wall `IW` to first record and count, rebuilt after a regrid on the host (D-047).

**2.3 Relation to W1 and W2 [R].**
- W1 wall tables stay per mesh object (per level and rank), in upstream order with upstream `IW`, and hold **owner faces only**. Records are not wall-table entries. [P] Record kernels get a record list in the W1 style: `DO IRI=1,NRL; IR=RLIST(IRI)`, with `R_OWNER_IW(IR)` giving the owner face. Loop bodies stay verbatim.
- W2 (`UVW_SAVE`, `U_GHOST`..) is read by the velocity side, not by the solid side. Nothing changes for records. Per D-061 the host-staged arrays are uploaded once per stage, and the solid passes stay at the existing `stage_wall_bc` position. Wall fluxes are domain-boundary fluxes, not C/F interface faces (spec `01` §3), so I expect the flux override sets of D-050 and D-061 not to contain them [P]. Role 1's host seam (wall-seam note, `FDSTL_WSEAM`) builds the identity lists `WLIST_EXT/INT`; record lists extend it.

**2.4 What a fine mesh object (D-056, option B) needs from the solid side.** Today a fine box has no wall-table rows (empty lists in the wall-seam note) and no obstructions (D-056). For OBST faces on level L it needs:
1. Wall-table rows for its owner faces (`WALL`, `BOUNDARY_COORD`, `SURF_INDEX`, `B1_INDEX`, `B2_INDEX`, `OD_INDEX`) from a per-level equivalent of `INIT_WALL_CELL` (init.f90:2975-3386), fed by the Legacy Mapper's indexing table and the ADR-003 face masks (OQ-S4).
2. `BOUNDARY_PROP1/2` for the owner face: gas-side inputs (`TMP_G`, `RHO_G`, `ZZ_G`, `U_TANG`, `RDN`, `K_G`, `Q_RAD_IN`) and aggregated outputs (§4). These are the owner's own rows, as in FDS.
3. `BOUNDARY_ONE_D`: not stored in the mesh object. [P] The host route stages one record into a scratch `ONE_D`/`B1`/`B2` slot, calls the unchanged kernel, and copies back; both copies are bitwise (A3). The device route reads the record rows directly.
4. The clock. `WALL_COUNTER` is a module global (cons.f90:588; incremented main.f90:1007), so every mesh object shares it. `BC_CLOCK` is per mesh (mesh.f90:296, init.f90:1178, restart dump.f90:3948) and `DT_BC` is a local (wall.f90:93-94). [P] One driver-owned `BC_CLOCK` and `DT_BC`; each fine object's `BC_CLOCK` copy is set from it before the stage. The step counter `ICYC` read by the kernel (wall.f90:2306) is a global (cons.f90:587).
5. Cross-mesh reads that abort for a fine mesh number: `HEAT_TRANSFER_COEFFICIENT` does `M => MESHES(NM)` (func.f90:3156); the back side does `MESHES(BACK_MESH)%WALL(BACK_INDEX)` (wall.f90:1876-1883, 1911-1927). Both need the box-view route or the FR-046 lookup (OQ-S3).

## 3. Time stepping
- **Single global dt [R, D-050; FR-047 (b)].** The solve runs in the corrector when `WALL_COUNTER==WALL_INCREMENT`, with `DT_BC = T - BC_CLOCK` (wall.f90:89-98). `WALL_INCREMENT` is 2 by default (cons.f90:588) and forced to 1 for surface oxidation (read.f90:7422). ADR-002 Option C as amended by D-050 (one dt, no subcycling) means no per-level clock; a per-level clock returns only if a new ADR reopens subcycling.
- **Records carry no clock [P].** Every record in every group gets the same `DT_BC` in the same step, so a split copy never needs its own phase and stays exact.
- **Order inside a corrector solve step [P]:**
  1. Owner pass: near-surface gas variables for owner faces (:108-119).
  2. Fan-out: copy the owner's gas-side inputs into each record of its group.
  3. Per-record HTC from the record's own `TMP_F` (equivalent of :117-118).
  4. Record loop: `SOLID_HEAT_TRANSFER(NM,T,SF%HT_DIM*DT_BC)` (:156); then the thin-wall lateral pass with `3*DT_BC` (:190-196, see OQ-S8).
  5. Fan-in to the owner (§4).
  6. Owner remainder of `WALL_BC`: `CALCULATE_RHO_D_F`, `CALCULATE_ZZ_F` tail, `CALCULATE_RHO_F`.
- **Every stage, not only solve steps.** The mass-flux block of `CALCULATE_ZZ_F` (:1185-1320) reads and updates record state each predictor and corrector: `T_IGN`, `BURN_DURATION`, `Q_IN_SMOOTH` (:1193-1207), `QDOTPP_INT`, `TAU_LS`, `AREA_ADJUST`. [P] Extract those lines verbatim as a record-local routine (no arithmetic change; bitwise at one record per face), run it per record every stage, then aggregate the fluxes. Without this, ignition and burn-duration differences inside a coarse owner are lost (OQ-S1). The random draw (:1324-1337) moves with it: SP-R9, counter-based, keyed by (face key, step).
- **Held values.** FDS holds the last solve's `Q_CON_F` and mass fluxes between solves. The owner's aggregated copy is refreshed after each solve and held the same way.
- **Initialisation.** `WALL_BC` at main.f90:506 runs with `CALL_HT_1D` false: near-surface pass and BCs only. Records are created at the depth of the t=0 hierarchy (D-058), then split later.
- **Radiation lag [P, as spec `01` §4].** `Q_RAD_IN` is set after `WALL_BC` in the corrector (`COMPUTE_RADIATION`, main.f90:1034; radi.f90:4917). Records take the owner's lagged value (OQ-S6).

## 4. Gas-side coupling under refinement
| Item | Rule |
|---|---|
| Inputs to records | From the owner's own gas cell: `TMP_G`, `RHO_G`, `ZZ_G`, `U_TANG` (wall.f90:402-446), `RDN`, `K_G`, `Q_RAD_IN`. The HTC uses the owner's gas-cell distance `1/RDN` (func.f90:3185, 3218, 3253), so a finer owner has a larger resolved-conduction floor, as FDS on a finer mesh does. Per-record: `TMP_F`, hence the film temperature and the temperature difference in the HTC [P]. |
| Fixed aggregation order | Records of a group are summed in ascending face key; owners then add into gas cells in ascending wall index (04 SP3/SP4). Per-record scratch slots, then ordered sums, no atomics. Depth-0 groups reduce to FDS order. |
| Linear outputs [P] | Area-weighted mean, so owner flux times owner area equals the sum over records: `Q_CON_F`, `Q_RAD_OUT` (radiation uses `OUTRAD_W = BBF*RPI*Q_RAD_OUT`, radi.f90:4238), `M_DOT_G_PP_ADJUST/ACTUAL` per species, `M_DOT_LAYER_PP`, `M_DOT_PART_ACTUAL`, `Q_DOT_G_PP`, `Q_DOT_O2_PP`. The same totals feed `DP` (divg.f90:544, which uses `AREA_ADJUST*Q_CON_F*RDN`; the plain area mean is exact for it while `AREA_ADJUST` is equal inside a group, OQ-S7), `D_SOURCE`/`M_DOT_PPP` (wall.f90:1412-1422), and the obstruction mass sum (:1378-1392, deferred). |
| Nonlinear outputs [P] | Owner `TMP_F` is the area mean (`02`); it feeds `RHO_F`, the blackbody fraction (radi.f90:4229-4235) and the enthalpy of wall species flux (divg.f90:366, 853). The error is second order in the spread of `TMP_F` inside a group. WP6 measures it; OQ-S5 if it breaks the FR-022 wall budget. |
| Budgets, boundary files [R, A4] | Owner aggregated values, owned faces only (FR-045); covered coarse faces count nothing. Zone and global sums over owner faces follow FR-005 (ii) as amended by D-053 (FDS order by default; the exact-sum switch ON for layout tests, FR-005 (iii)). |
| Point wall devices [R, A4] | The finest-ever record containing the point, read directly, no aggregation. Smokeview patches follow D-045. |
| Back side [R: FR-046; P: depth] | Face-key lookup, filled by a per-step exchange. [P] Record to record when the back face has the same depth; otherwise the back owner's aggregated values. |
| Thin OBST gas source | Owner level; `D_SOURCE`, `M_DOT_PPP` from the aggregated flux times owner area (04 SP4). |

## 5. Regrid events
- **No hierarchy change (Phase 5).** The FR-040 checker keeps every stateful face on its level (D-010; static crossing allowed). Stateless faces (inert, `SPECIFIED`) have no record and are rebuilt from `SURF` data (FR-041a).
- **First split [R, A1, A3; FR-041b].** A face with a record at depth d acquires an owner finer than d. The record becomes `r²` rows: solid state copied bit for bit; thickness inherited (A1); `AREA` = parent area / `r²` (exact in binary, since `r²` is a power of two); geometry from the key. `AREA_ADJUST` is inherited at split (OQ-S7). Mass and enthalpy per area are unchanged, so the sum over records of (content × area) is unchanged. Host regrid-time work (D-047); [P] called from Role 3's regrid bracket (`begin_regrid`/`end_regrid`, D-058). At the split instant the group's area-mean `TMP_F` and fluxes equal the parent's to round-off, because the group is `r²` identical rows (a sum of identical values is not always bitwise equal, so the check is a tolerance of 1e-12 relative).
- **Coarsening [R, A3].** No record changes. Only the owner pointer changes: the coarser owner gets the CSR group. Owner values may differ from the former fine owners' at round-off, because they are now sums over a group. After the regrid the owner's gas-side rows are rebuilt from the records by one fan-in plus the near-surface pass, before the first stage.
- **Every later regrid [R, A3].** Every record is bitwise unchanged. Migration between ranks moves rows as raw bytes (not through `PACK_WALL`, func.f90:4482-4530).
- **New wall faces under a regrid [R, A2].** The face set is fixed by geometry (OBST creation is deferred), so regrid adds owners, not faces. A face gets a finer owner only inside the refinable region and up to the level the region permits; tags outside are discarded (IR-008). Outside the region the face keeps one record at its level-0 or static-box resolution.
- **Restart [R, A5].** Replace the per-mesh wall dump (dump.f90:3975-3984, `BC_CLOCK`/`WALL_COUNTER` at :3948, read :4137) by rows written in face-key order with depth, plus (proposal) the owner rows' held values and `BC_CLOCK` and `WALL_COUNTER` once. A restart on another layout or rank count must give bitwise records. The held owner values and the single clock are proposals [P]; the records and face keys are ruled (A5).

## 6. Cases, tests and GPU notes
Cases that use walls and need no burn-away carry the refined tests. `couch` (BURN_AWAY at `couch.fds:48`) and `box_burn_away*` (`box_burn_away1.fds:26-27`) need burn-away, so they stay cost evidence only until FR-042 is ruled.
| Stage | Case | Check |
|---|---|---|
| S0, level 0 only | `Energy_Budget/energy_budget_solid`, `Pyrolysis/surf_mass_vent_char_cart_fuel`, `Heat_Transfer/heat_conduction_a` | Kernel on frozen record and inputs bitwise vs single-mesh FDS (FR-047 (a), a candidate until D-022 is amended); whole run T2; records bitwise across box split, ranks, threads (FR-005 (i), exact-sum switch ON). |
| S1, static 2-level (Phase 5) | Refined variants (A-52) of the three above at 2:1 over the heated face; `back_wall_test` with front and back on different levels, two layouts | One active record per face, zero for covered coarse faces (FR-045 report); solid mass and enthalpy closed against net surface flux and gas gain to round-off (FR-022 walls; <= 1e-12 relative as the working threshold); `back_wall_test` wall temperatures byte-identical across layouts (FR-046). |
| S2, G2 (after Phase 6) | Static 2-level box, heated charring OBST without burn-away: refine, coarsen, refine across its face | Every record bitwise unchanged across each regrid after its first split; solid mass and enthalpy change <= 1e-12 relative; owner surface T continuous to round-off; restart on another rank count bitwise. |
| Needs FR-042 | `box_burn_away1`..`11`, `box_burn_away_2D(_residue)`, `couch` | Deferred (§9). |
| Mode 4 level set | `WUI/level_set_fuel_model_1` | Deferred in detail (§9); the level-0 column rule of spec `01` §2 stays the proposal. |
- **Cost, G2a (`03`).** Worst-case bound at R=2 (3 x share): `couch` 0.16 (wall timer) and 0.06 (isolated solve), gate 25%; `box_burn_away_2D_residue` 0.40 on the solve, over the gate, a 4 s solid-dominated case; no depth cap needed now. Memory: 2.9-8.6 KB per record (assumed node counts); `couch` 39 k records at R=2, 158 k at R=4.
- **Where the time goes.** The isolated solve is about 0-13% of step time (0 to 0.133 across the four cases; `box_burn_away1` is noise), the whole `WALL` timer 5-23%. In three of four cases the gas-side wall work (near-surface pass, HTC, `CALCULATE_ZZ_F`, fan-in) outweighs the 1-D solve (`couch` narrowly: 0.031 against 0.021); only the small charring 2-D case is solve-dominated. It scales with owner faces and with records, so it is the first GPU target.
- **GPU [P].** Fan-out and fan-in are per-record scratch slots and ordered sums over a CSR list (owner to records, gas cell to walls), no atomics (04 SP4); the solve is one thread per record over `RLIST`; the load-balance weight counts records x `NWP` x expected sub-steps (NFR-035), not faces. Regrid-time split and CSR rebuild may stay on the host until Phase 11 (D-047). [P] The owner-to-record CSR is the same table family as the cell-to-wall CSR that the wall-loop work plans for SP1-SP3 (L1358, L1359, L0375); build one. The thin-wall branch of the near-surface pass averages the two side walls' `TMP_G` and `Q_RAD_IN` (wall.f90:473-486, run at :122-137), so it stays a separate pass after the owner pass.

## 7. Work packages (dependency order)
| WP | Content | Owner | Deliverable and check |
|---|---|---|---|
| 1 | Field inventory of record state (R-60) and record schema | Solid Lead | Table of every `ONE_D`/`PROP1`/`PROP2` field with class (state, frozen, input, output); Legacy Mapper cross-check against reads and writes in `WALL_BC`. |
| 2 | Face key, record store, CSR owner to records, depth = owner level | Role 1 Data Layout | Round trip; level-0 state bitwise equal to FDS wall state in the S0 cases; row order independent of layout. |
| 3 | Fine-level OBST wall tables in fine mesh objects; driver clock | Role 1 | Tables equal those of a single-mesh equivalent of a 2-level case; FR-045 active-record report. Needs OQ-S4. |
| 4 | Driver loop of §3 incl. `CALCULATE_ZZ_F` split, fan-out, fan-in, per-record HTC | Solid Lead (spec, review), Role 1 (host code) | Depth-0 bitwise vs FDS on frozen inputs; FR-005 (i) at 1/4 ranks, 1/2 threads. Needs OQ-S1, OQ-S3. |
| 5 | Back-side lookup by face key | Role 1 with the Integration Lead | `back_wall_test` at two layouts byte-identical (FR-046). |
| 6 | Owned-face budgets, wall devices, boundary files, restart by key | Role 1; V&V Lead for tests | A4 and A5 checks; measured spread effect of area-mean `TMP_F` on the wall budget. |
| 7 | Phase 5 acceptance: S1 cases, thresholds | V&V Lead | Case files, baselines, pass thresholds in the V&V plan. |
| 8 | First-split copy and rank migration inside the regrid bracket | Role 3 Regrid (hook), Role 1 (store) | Copies bitwise; records unchanged across later regrids; split is a no-op on a no-change regrid. |
| 9 | G2 run (S2) | V&V Lead, Solid Lead | Pass signals of §6, S2 row. |
| 10 | Device passes: fan-out/fan-in/HTC/near-surface first, solve kernel second | Wall Loops Engineer | FR-005 (i) on GPU; scratch-sum order test as in the SP4 sign-off. |
Order: 1 -> 2 -> 3 -> 4 -> (5, 6) -> 7 = Phase 5 exit. 8 -> 9 only after Phase 6. 10 follows 4, parallel to 5-7, gated by Phase 11.

## 8. Open questions (rulings needed)
- **OQ-S1** Is extracting the mass-flux block of `CALCULATE_ZZ_F` as a record-local routine allowed? A6 names only `SOLID_HEAT_TRANSFER`. The extraction is verbatim, but it changes a second routine. Ruling: AMR Chief Architect.
- **OQ-S2** Face-key definition (shared with FR-046, which names the Solid Lead with the Integration Lead as owners). Ruling: AMR Chief Architect.
- **OQ-S3** Route for `HEAT_TRANSFER_COEFFICIENT` and back-side reads of `MESHES(NM)` on fine mesh objects: box-view hook in the shim, or an upstream patch under D-051. Ruling: Chief Architect, input from the Legacy Mapper.
- **OQ-S4** Who builds fine-level OBST wall tables and when, given D-056 ("no obstructions at fine levels until the geometry work"). Proposal: Role 1, in Phase 5. Ruling: Chief Architect.
- **OQ-S5** If area-mean `TMP_F` breaks the FR-022 wall budget, which closes it: flux-weighted `TMP_F` for the blowing enthalpy, or a looser budget? Decide only if WP6 shows a miss. Ruling: Chief Architect with the V&V Lead.
- **OQ-S6** Does the wall radiation history `ILW` (FR-063) and `Q_RAD_IN` belong to the record or to the owner? Proposal: owner. Ruling: Radiation Lead.
- **OQ-S7** `AREA_ADJUST` at fine levels (init.f90:337-421, D-028 `FDS_AREA` sums): inherit at split, or recompute from fine areas? Proposal: inherit. Ruling: Spec Lead with the Solid Lead.
- **OQ-S8** Thin-wall faces: may they carry records deeper than the owner? Proposal: no. R3-T keeps C/F interfaces off thin walls, so they are not split; the lateral pass runs at owner level. Ruling: Chief Architect.

## 9. Deferred (not planned here)
- **OBST creation, removal, burn-away mass bookkeeping across levels (FR-042, R-58).** To decide: how `OBSTRUCTION%MASS` is held per obstruction when its faces have records on several levels; the per-record contribution (`M_DOT_*_ACTUAL x DT x AREA`, wall.f90:1378-1392, burn-away :2689-2697) and its fixed-order sum (04 SP4); who removes the cells and rebuilds masks and wall tables on each level; what happens to a record group when its OBST goes (drop, or keep for restart); the `-1` mass flag and consumable flags across a level jump; the pressure-matrix rebuild; the T2 criterion for `box_burn_away1` (D-038 revisit clause).
- **3-D heat transfer (HT3D)** under refinement: deferred with GEOM (D-033) and rejected in AMR mode under FR-004 (ADR-003). `HT_DIM>1` surfaces therefore never reach the kernel.
- **Level-set fire spread (SP-R10, mode 4):** level-0 `PHI_LS` column rule versus a per-level solve, record `T_IGN`/`AREA_ADJUST` source for fine faces, simultaneous ignition inside a coarse column (R-61). Needs the Chief Architect's ruling (FR-045 open point).
- **HVAC wall BCs:** owned faces only (FR-034); how the duct mass flux splits across records.
- **Solid particles** (wall.f90:244-274): the particle record travels with the particle, ownership by particle phase; owners: Solid Lead with the Integration Lead as Phase 7 lead (FR-045 open point).
- **CFACE paths** (wall.f90:200-215): belong to the deferred GEOM work (D-033).
