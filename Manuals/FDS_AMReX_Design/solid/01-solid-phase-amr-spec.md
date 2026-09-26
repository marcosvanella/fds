# 01 — Solid Phase & Wall Coupling on AMReX: spec v0.1

Owner: AMR Solid Phase Lead · Status: **v0.1 for review**, 2026-09-26 (outline 2026-09-25) · Proposals only; nothing here is decided unless it cites a D-xxx.
Reference tree: this repository (FireX `36975d765f`), read-only. Line numbers are from that tree.
Binding: uniform grids per level (D-030); GEOM/CC_IBM and HT3D deferred, thin OBSTs required (D-033); FR-006; bitwise per-kernel on frozen inputs, whole runs by tolerance (D-022); refluxing (D-023); Phase 5 stateful-wall freeze (D-010, FR-041a); wall ownership rule (ADR-003 v0.2.1).
Scope: FDS data flow, level mapping, conservation, time stepping, GPU and cost, verification, phase plan, traceability. CFACE paths are listed for completeness; they belong to deferred GEOM. Folded into the project docs in v0.4.17 as FR-045..047, NFR-048, the FR-022 wall part, the FR-041b proposed rule, NFR-035 (walls) and R-55..R-60 (§8).

## 1. Current FDS data flow (FireX)

| Stage | Where | Notes |
|---|---|---|
| Storage | `WALL_TYPE` type.f90:434-457 holds indices into `BOUNDARY_COORD` (189-214), `BOUNDARY_ONE_D` (217-264), `BOUNDARY_PROP1` (292-345), `BOUNDARY_PROP2` (350-376), `BOUNDARY_RADIA` (381-384). `THIN_WALL_TYPE` 487-503; `CFACE_TYPE` 1361-1382. | `ONE_D` = in-depth solid state: node depths `X`, `TMP`, per-material densities (`MATL_COMP`), layer thicknesses, `N_CELLS_MAX`. `B1` = surface state: `TMP_F`, `Q_CON_F`, `Q_RAD_IN`, `M_DOT_G_PP_ADJUST/ACTUAL`, `U_NORMAL(_S)`, `ZZ_F`, `BURNAWAY`, `AREA_ADJUST`. All `ALLOCATABLE` components. |
| Build / rebuild | `INIT_WALL_CELL` init.f90:2975-3386; `ALLOCATE_STORAGE` func.f90:3942-4384; `REASSIGN_WALL_CELLS` init.f90:4898-5228; `CREATE_OR_REMOVE_OBST` init.f90:4852 | Wall cells are per mesh, indexed by `CELL%WALL_INDEX`. |
| Pack / exchange / restart | `PACK_WALL` func.f90:4482-4530; sent main.f90:3467, received :3845; restart dump.f90:3981 | Carries back-side wall data between meshes. |
| Driver | `WALL_BC` wall.f90:26-279, called main.f90:506 (init), :838 (predictor), :1009 (corrector) | `LEVEL_SET_MODE==1` returns at once (:48). |
| 1-D solve cadence | wall.f90:89-98: corrector only, when `WALL_COUNTER==WALL_INCREMENT`; `DT_BC = T - BC_CLOCK` | `WALL_INCREMENT=2` default (cons.f90:588), forced to 1 for surface oxidation (read.f90:7422); counter main.f90:1007-1011. |
| Gas-side near-wall values, h | `NEAR_SURFACE_GAS_VARIABLES` wall.f90:402; `HEAT_TRANSFER_COEFFICIENT` call :117-118 | Loop 0 (:109-119) also fills external ghosts (`ASSIGN_GHOST_VALUE` :282). |
| Thin/specified BCs | `SURFACE_HEAT_TRANSFER` wall.f90:603 (:153-154) | Non-thermally-thick surfaces. |
| 1-D conduction + pyrolysis | `SOLID_HEAT_TRANSFER` wall.f90:1809-2992 (:155-156); thin OBST lateral pass `3*DT_BC` (:190-196) | Sub-step loop :2044, adaptive `DT_BC_SUB` :2323-2325, Fourier limit :2096; `PERFORM_PYROLYSIS` :2993 (called :2213, :2329) and `PYROLYSIS` :3195 (MATL reactions, char/residue, liquid evaporation); layer removal and renoding :2361-2804; burn-away :2689-2697; tridiagonal :2856-2895. |
| Back side | wall.f90:1876-1883, 1911-1927, 2163-2176 | `EXPOSED` backing reads `MESHES(BACK_MESH)%WALL(BACK_INDEX)` directly: a cross-mesh access. |
| Surface mass flux | `CALCULATE_ZZ_F` wall.f90:1035; `M_DOT_G_PP_ADJUST` scaled by `AREA_ADJUST` :1316; blowing `U_NORMAL_S = -UN` :1446 | **Thin OBST or removed layer:** no blowing; mass enters the adjacent gas cell as `D_SOURCE` and `M_DOT_PPP` inside `!$OMP CRITICAL` (:1402-1416). |
| Gas-side sources | divg.f90:213-228 (`RHO_D_DZDN_F` into face diffusive flux); divg.f90:544 (`DP -= AREA_ADJUST*Q_CON_F*RDN - Q_LEAK`); divg.f90:671 (`D_SOURCE`); divg.f90:757-767 (zone `USUM` over solid walls) | Thin walls processed on one side only (divg.f90:196, 334). |
| Mass-flux variation | `CALCULATE_ZZ_F` wall.f90:1324-1337 | `MFT` is the face's own species sum (:1319); `MASS_FLUX_VAR` perturbs it with a per-face `BOX_MULLER` draw and rescales the same face. Local, but the random sequence follows the `SCHEDULE(DYNAMIC)` wall loop (wall.f90:143), so it depends on threads and, in AMR, on layout. No Verification case sets `MASS_FLUX_VAR`. |
| HVAC wall BC | `CALC_HVAC_BC` wall.f90:1745-1797 | Per-face duct mass flux `MFT` (:1761), `TMP_F`, `Q_LEAK`, `U_NORMAL(_S)` or `M_DOT_G_PP_ADJUST`. FR-034 ruling restricts it to owned faces. |
| Solid particles | `WALL_BC` wall.f90:244-274; `DEPOSIT_PARTICLE_MASS` :1468-1581 | Same 1-D solve per particle (`PARTICLE_INDEX`); off-gas deposited into the particle's gas cell in the corrector. |
| Level set to walls | `LEVEL_SET_FIRESPREAD` vege.f90:502-597 | 2-D field `PHI_LS(IIG,JJG)` (type.f90:1075) writes `B2%PHI_LS` and, on arrival, `T_IGN`, `TAU_LS`, `BURN_DURATION` and (coupled fire) `AREA_ADJUST` into the terrain wall face `WALL_INDEX(-3)` (:516-550, 567-592). `CROSSING_DISTANCE` uses `DX`,`DY` (:570-575), so these values depend on resolution. Mass flux then follows the trapezoid at wall.f90:1305-1312. Modes (read.f90:1986-2013): 1 returns before wall BCs (wall.f90:48); 2 and 3 couple wind only; 2 also freezes velocity after ignition (vege.f90:344-347); 4 couples fire. Verification uses modes 1 and 4 only. |
| Radiation in | radi.f90:4917 (`Q_RAD_IN` from `INRAD_W`) after `WALL_BC` in the corrector (main.f90:1034) | Solid sees radiation lagged by one call. |

Key property: given its gas-side inputs, each wall cell's 1-D solve is **local**, except the back-side lookup, the thin-wall volumetric source accumulation, and the random draw of `MASS_FLUX_VAR`.

## 2. Wall cells and solid state on an AMReX hierarchy

**Ownership (ADR-003).** A wall face belongs to the level whose valid, uncovered cell is its gas side (divg.f90:760 test). Exactly one wall record is active per physical face per step. Coarse wall cells under a finer level are inactive (no solve, no source, no output).

**Where state lives.**
- Phase 5 (D-010): per level, on the owning level, built by the per-level `INIT_WALL_CELL` equivalent. A static hierarchy may cross a stateful wall; a regrid may not cut one.
- FR-041b recommendation, **option A, "finest-ever records"** (ADR-003 sub-choice, line 67): key each 1-D record to the OBST face patch at the finest resolution that has ever covered it, stored outside the box layout. When a coarse level owns the face, each of its r² fine records is solved with the coarse face's gas-side inputs; the coarse face's mass, species and heat fluxes are the area-sum of the records. Refine copies nothing (the records exist); first-time refinement splits a coarse record into r² identical copies (exact per-area). Conservation is exact and reversible, and burn-away heterogeneity is kept. Cost: solid work stays at the finest-ever count.
- Option B, remap: coarsen by area-averaging r² profiles. Needs a conservative depth remap (profiles have different node counts after renoding, :2361-2804), per-material mass and enthalpy with T recovered from a nonlinear `RHO_C_S(T)`, plus rules for `T_IGN`, `BURNAWAY`, `INT_FTP`, `PART_MASS`, `A_LP_MPUA`. Lossy on re-refinement. Not recommended.

**Thin OBSTs.** `THIN_WALL` cells and lateral conduction follow the owning level. The R3-T rule (ADR-003 v0.2.1) keeps a C/F interface off thin-wall faces.

**Build and rebuild (Phase 5).** Wall records for a level are built when the level is created, by a per-level equivalent of `INIT_WALL_CELL` (init.f90:2975-3386) over that level's OBST faces. Under the D-010 freeze a regrid may add or remove fine boxes only where no stateful face changes owner; the FR-040 mask checker enforces it. Stateless faces (inert, `SPECIFIED` BCs) are rebuilt from SURF data after any regrid, with nothing to transfer. Restart writes records per owning level with their face keys (replacing dump.f90:3981).

**Special paths inside `WALL_BC` (OQ-4, proposals).**
- Solid particles: the particle's 1-D record travels with the particle and is never a face record, so ownership follows the particle phase. Gas-side inputs and `DEPOSIT_PARTICLE_MASS` use the level whose valid cell contains the particle. The particle lead owns the rule; the solid kernel is unchanged.
- HVAC wall BC: owned faces only (FR-034 ruling). Covered coarse HVAC faces set no `U_NORMAL`, `Q_LEAK` or `M_DOT_G_PP_ADJUST`. Duct flow is split across owned faces by area, as FDS does per face with `NODE_AREA_EX` (:1761), so the vent total is kept.
- Level set, mode 4 (required, D-033): `PHI_LS` stays a 2-D field solved on level 0. An owned terrain face on a finer level takes `T_IGN`, `TAU_LS`, `BURN_DURATION` and `AREA_ADJUST` from the level-0 column that contains it. The fuel energy per area is then the level-0 value (as `AREA_ADJUST` intends, vege.f90:591-593), and ignition within a coarse column is simultaneous. A per-level level-set solve is the alternative; it changes the spread result with resolution and belongs to whoever owns the level-set module. Modes 2 and 3 take the same rule; they have no Verification case, and mode 2's global freeze (vege.f90:344-347) needs nothing per level.

**Back-side coupling.** Replace `MESHES(BACK_MESH)%WALL(BACK_INDEX)` with a lookup through a level- and layout-independent face key, filled by a per-step exchange (front and back may be on different levels or ranks).

**`VARIABLE_THICKNESS`** (charter Q4 item 1). It is a 1-D feature. It only shares noding and grouping code with HT3D (init.f90:1594-1650, 3926-4051; read.f90:841-849). Its thickness comes from the snapped OBST depth, so it depends on the level. **Recommendation: IN**, with thickness computed once at setup from the owning level at t=0 and frozen as part of the solid state. This matches FDS for uniform and static runs and is compatible with option A. `box_burn_away*` also needs FR-042 (burn-away removal).

## 3. Conservation of wall fluxes across coarse/fine

- Wall fluxes are domain-boundary fluxes, not C/F interface fluxes. The owning level writes them into its face fluxes (blowing plus `ZZ_F`, `RHO_D_DZDN_F`) and its `DP`/`D_SOURCE`. Covered coarse wall faces take `average_down_faces` like any covered face, so the D-023 flux register sees them consistently.
- Thin-wall volumetric sources (wall.f90:1411-1414) land in the owning level's gas cell and reach the coarse level through `average_down`.
- Budgets and outputs (HRR, solid mass loss, `OB%MASS` via `M_DOT_LAYER_PP`, zone `USUM` divg.f90:757-767, wall DEVCs) sum owned faces only. Double counting of covered coarse faces is the main risk.
- Transfer (FR-041b): solid mass and enthalpy are conserved to round-off across refine and coarsen (exact under option A).

## 4. Time stepping and subcycling

- Global dt (ADR-002 leaning): keep one global `WALL_COUNTER` and `BC_CLOCK`, so every level solves on the same steps as FDS; this is needed for bitwise kernel checks.
- Subcycling (Phase 6): `BC_CLOCK` and `WALL_COUNTER` become per level; each level's walls advance over that level's accumulated `DT_BC` with time-averaged gas-side inputs. Under option A, fine records under a coarse owner advance with the coarse dt. Solid internal sub-stepping (:2323-2325) already absorbs dt changes.
- Radiation lag (main.f90:1034 after :1009) stays as in FDS.

## 5. GPU, cost and load balance

- The kernel is sequential in depth: variable `NWP`, adaptive per-cell sub-steps, renoding and layer removal, per-material reactions, branching by backing and pyrolysis model. Expect warp divergence and a long tail of cells.
- `ALLOCATABLE` components (type.f90:230-262, 294-302, 352-356) cannot be used in device kernels. Proposed layout: flat SoA padded to `N_CELLS_MAX` (type.f90:220) per SURF class, one thread per wall cell, with the tridiagonal kept in registers or local memory.
- The `!$OMP CRITICAL` accumulation (wall.f90:1412-1415) becomes an atomic on GPU, whose order is nondeterministic and breaks bitwise checks. Use a deterministic gather per gas cell instead (as in D-031).
- Wall cells scale with r² on refined surfaces, and their cost is not proportional to box cell count. The AMReX distribution map needs a per-box weight: wall cells × `NWP` × expected sub-steps. `SCHEDULE(DYNAMIC)` (wall.f90:143) shows FDS already sees this imbalance.
- `MASS_FLUX_VAR`: replace `RANDOM_NUMBER` in AMR mode with a counter-based generator keyed by (face key, step), so each face's draw is independent of threads, boxes and ranks (FR-005 (i)). Refined faces get r² independent draws, as a finer FDS mesh would.
- Per-kernel bitwise candidate: `SOLID_HEAT_TRANSFER` on a frozen wall record and frozen gas-side inputs, compared with single-mesh FDS. Per the FR-047 D-022 check, this is a candidate until the Chief Architect amends D-022; until then FR-005 (i) and whole-run T2 cover it.

## 6. Candidate verification cases (Verification/, real filenames; mesh counts not yet checked)

- Conduction and BCs: `Heat_Transfer/heat_conduction_a`…`_d`, `heat_conduction_kc`, `back_wall`, `back_wall_test`, `back_wall_test_2`, `adiabatic_con_flux`, `adiabatic_net_flux`, `convective_cooling`, `insulated_steel_plate`, `internal_heating`, `SFPE_Case_1`.
- Pyrolysis, char, liquids: `Pyrolysis/pyrolysis_1`, `pyrolysis_2`, `two_step_solid_reaction`, `matl_e_cons_1`…`_9`, `enthalpy`, `shrink_swell`, `cell_burn_away`, `methanol_evaporation`, `water_pool`, `liquid_mixture`, `specified_hrr`.
- Mass conservation at walls: `Pyrolysis/surf_mass_vent_char_cart_fuel`, `surf_mass_vent_nonchar_cart_gas`, `surf_mass_two_species_cart`, `surf_mass_vent_liquid_fuel`.
- Energy budget: `Energy_Budget/energy_budget_solid`, `energy_budget_cold_walls`, `energy_budget_adiabatic_walls`.
- Burn-away and `VARIABLE_THICKNESS` (pending Q4): `Fires/box_burn_away1`…`11`, `box_burn_away_2D_residue`.
- Level set, mode 4: `WUI/level_set_fuel_model_1`, `LS4_ember_ignition`, `LS4_ember_yield`, `ground_vegetation_drag`. Mode 1 cases (`WUI/LS_ellipse_*`, `LS_wind_ramp_*`, `Bova_*`) never call wall BCs and need only the level-0 rule. `Restart/geom_ls_restart_*` uses GEOM and is deferred.
- Excluded (HT3D deferred): `Heat_Transfer/ht3d_*`. Particle-surface cases (`surf_mass_part_*`) follow the particle phase.
- Refined variants to write (A-52, A-46 style): `energy_budget_solid` and `surf_mass_vent_char_cart_fuel` at 2:1 over the burning face; `back_wall_test` with front and back on different levels; a static C/F crossing a stateful wall (D-010 clarification); `level_set_fuel_model_1` with a 2:1 patch over part of the burn area (checks the level-0 rule and HRR against the unrefined run).

## 7. Phase plan (solid and wall work only)

| Phase | Work | Exit check |
|---|---|---|
| 5 | Per-level wall records, ownership mask (FR-045), face key and back-side exchange (FR-046), global `WALL_COUNTER`/`BC_CLOCK` (FR-047 b), owned-face budgets and outputs, HVAC and level-set rules (§2), counter-based `MASS_FLUX_VAR` draw, restart per level, freeze enforced (FR-041a). | FR-022 walls and FR-045 on the A-52 variants; FR-046 on `back_wall_test` at two layouts; FR-047 b trace. |
| 6 | Per-level wall cadence with subcycling; time-averaged gas-side inputs (§4); OQ-3 answered. | FR-022 walls with subcycling. |
| after 6 | FR-041b transfer, option A if G2 supports it; FR-042 burn-away with refinement. | Solid mass and enthalpy change across regrid ≤ 1e-12 relative; `box_burn_away*` refined. |
| 11 | Device port: flat SoA (NFR-048 a), deterministic gather (NFR-048 b), load-balance weight (NFR-035 walls). | FR-005 (i) on GPU; A-measurement of wall share of step time. |

## 8. Traceability, new proposals, open questions

| Spec item | Project ID | State |
|---|---|---|
| SP-R1 single active record | FR-045 | proposed |
| SP-R2 wall budget closure | FR-022 (walls), FR-020 | proposed |
| SP-R3 back-side coupling | FR-046 | proposed |
| SP-R4 kernel parity, cadence | FR-047 | proposed; (a) needs D-022 amendment |
| SP-R5 FR-041b rule | FR-041b proposed rule | open until G2 (OQ-2) |
| SP-R6 device layout | NFR-048 | proposed |
| SP-R7 load-balance weight | NFR-035 (wall part) | proposed |
| SP-R8 `VARIABLE_THICKNESS` frozen at setup | FR-042 note, charter Q4 | pending the project owner's decision (OQ-1) |
| SR-1..SR-6 | R-55..R-60 | open |

New in v0.1 (for the Spec Lead; IDs to assign)
- **SP-R9 Layout-independent `MASS_FLUX_VAR`.** In AMR mode the per-face draw comes from a counter-based generator keyed by (face key, step); results are byte-identical across layouts and rank counts (FR-005 (i)). Verification: a derived case with `MASS_FLUX_VAR` at two layouts. Closes OQ-5 (it is not a cross-cell reduction; the Spec Lead's reading in FR-043 is confirmed).
- **SP-R10 Wall BC ownership for HVAC and level set.** HVAC wall BCs act on owned faces only, with the vent total kept; level-set wall data on finer owned faces come from the level-0 column (§2). Verification: `duct_flow` refined at the vent (FR-043); `level_set_fuel_model_1` 2:1 variant with HRR T2 against unrefined.
- **SR-7 Level-set values depend on resolution.** `BURN_DURATION` and `AREA_ADJUST` use `DX`,`DY` (vege.f90:570-593). A per-level level-set solve would change spread and HRR with refinement; the level-0 rule avoids that but ignites a whole coarse column at once. (M/M)

Open questions
- **OQ-1 (project owner)** `VARIABLE_THICKNESS`: IN with frozen thickness (recommended), or deferred with HT3D? `box_burn_away1` is UNCLEAR for FR-006 until then.
- **OQ-2 (Chief Architect)** FR-041b option A vs B; G2 measures option A. Not raised with the Architect yet (coordinator hold).
- **OQ-3 (with Chief Architect, ADR-002)** Keep `WALL_INCREMENT=2` in AMR mode, or tie it to the level?
- **OQ-4 (Spec Lead to assign)** Owners for the solid-particle rule (particle phase) and the level-set module (level-0 rule vs per-level solve). Proposals in §2.
- **OQ-5** Closed by SP-R9.

## Change log
- v0.1, 2026-09-26: closed the *(unread)* items (`MASS_FLUX_VAR`, HVAC, solid particles, level-set modes); added build/rebuild, special-path rules, phase plan, traceability to v0.4.17 IDs, SP-R9, SP-R10, SR-7.
- Outline, 2026-09-25.
