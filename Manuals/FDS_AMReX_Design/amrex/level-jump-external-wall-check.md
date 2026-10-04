# Level-jump check: does the AMR route ever build an FDS `EXTERNAL_WALL` at a coarse/fine jump? (D-072 (f), D-065 Q1)

## Summary (5 lines)
1. **Stage execution: no.** At a coarse/fine jump no `EXTERNAL_WALL`, no `OMESH` and no wall cell exists; the jump is handled by driver ghost fill and the interface flux overwrite, so the `NIC>1` branches (`wall.f90:897`, `divg.f90:219`) cannot run on a fine level. Level 0 refuses unequal cell sizes, so its mesh faces all have `NIC=1`.
2. **Set-up phase on a multi-resolution input: yes, and then the run aborts.** The unchanged FDS set-up runs first (`main.cpp:97`), builds `EXTERNAL_WALL` with `NIC>1` at every jump (`main.f90:2184`), and only then `build_level0` stops at `FdsAmr.cpp:30`. No stage ever runs on those walls.
3. **Not provable for the designed route.** The code that turns finer `&MESH` entries into levels (`Hierarchy.cpp`, IR-002) is called only by the two-level test harness (`TwoLevelRun.cpp:98`), never by `main.cpp`; the planned wiring takes the mesh list from the FDS set-up, which would have built those walls. No guard on `NIC>1` exists in the driver or in `wall.f90`/`divg.f90`.
4. **Seven inputs:** five in scope (`ns2d_16_int_1to2_refinement`, `obst_coarse_fine_interface`, `dancing_eddies_uglmat_refine`, `race_test_1`, `race_test_4`), one unclear (`duct_flow_uglmat_refine`), one deferred (`geom_stretched_grid`). Only `race_test_1` and `race_test_4` have two or more tracked species. Today all seven abort at level-0 assembly; the four embedded inputs abort there too.
5. **Validation inputs named by the V&V plan** have no coarse/fine face; 325 of 3,132 Validation inputs do (256 with a reaction or species line), none of them in the plan. **Gap for the Architect:** decide how finer `&MESH` lines are kept away from the FDS set-up, and add the `NIC>1` abort guard that D-065 Q1 assumes.

Evidence tree: `Source/` of the FDS-AMReX branch at `991a759f79`; `main.f90` line numbers there are 73 higher than at `36975d7` (`main.f90:2111` there is `main.f90:2184` here). Inputs: `Verification/` and `Validation/` of the same tree. Nothing was run; all findings are from reading code and inputs.

## 1. What the route does at a jump (question (a))

### 1.1 Level 0: the FDS set-up runs unchanged, then the layout is checked
| Step | Code | Effect |
|---|---|---|
| FDS set-up | `Source/driver/main.cpp:97` (`fds_setup(0, ...)`); `Source/main.f90:327` calls `INITIALIZE_MESH_EXCHANGE_1`, which sets `EWC%NIC` at `Source/main.f90:2184`; the routine returns at `main.f90:754-760` | The set-up builds `EXTERNAL_WALL` cells with `NOM>0` and `NIC` for every mesh-to-mesh face of the input, exactly as FDS does. Nothing about levels is known yet. |
| Layout | `Source/driver/main.cpp:99` calls `build_level0` (`FdsSetup.cpp:18-41`), which calls `assemble_level0` (`FdsAmr.cpp:15`) | Aborts "all meshes must have the same cell size (no mesh refinement at level 0)" at `FdsAmr.cpp:30`; aborts if a corner is off the lattice (`FdsAmr.cpp:53`) or meshes overlap (`FdsAmr.cpp:62`); `TRN*` meshes abort at `FdsSetup.cpp:31`, `CYLINDRICAL` at `FdsAmr.cpp:20`. |
| Result | | Level-0 meshes have one cell size, lattice-aligned, disjoint. |

Why every level-0 face has `NIC=1`: `NIC` is the product of the index ranges of the other-mesh cells found by two probe points that sit 0.475 cell widths either side of the wall-cell centre in each tangential direction (`init.f90:3147-3170`, offsets at `init.f90:3161-3163`; range at `main.f90:2184`). With equal cell size and an aligned lattice both probes fall in one cell. A misaligned pair stops with `ERROR(431)` (`init.f90:3185`, `3204`, `3222`). This is a reading of the code, not asserted by any driver test (`grep NIC` finds nothing in `Source/driver`).

What the level-0 `EXTERNAL_WALL` cells do: the driver turns them into open faces in the side data (`fds_mesh_query.f90:57`, `73-75`: code 2 = interface with `NOM>0`), copies box data into `OMESH` for same-level boxes only (`GhostExchange.cpp:93-112`), and runs the unchanged `SURFACE_HEAT_TRANSFER` `INTERPOLATED_BC` case (`wall.f90:833`), which reads `OMESH(EWC%NOM)` (`wall.f90:835-836`). That case contains the `NIC>1` branch (`wall.f90:897`, inside `SINGLE_SPEC_IF` at `wall.f90:861`), which is skipped when `NIC=1`. `DIVERGENCE_PART_1` overwrites the coarse flux only `IF (EWC%NIC>1)` (`divg.f90:219`).

### 1.2 Levels above 0: no wall tables, no `OMESH`
| Fact | Code |
|---|---|
| A fine box is built with zero wall cells: `ALLOCATE(M%WALL(0:0), M%EXTERNAL_WALL(0), ...)`, `N_WALL_CELLS=0`, `N_INTERNAL_WALL_CELLS=0`, `N_EXTERNAL_WALL_CELLS=0` | `Source/driver/fds_fine_level.f90:163-165` (comment at `:111`: "a fine box has no WALL cells ... there is no OMESH neighbour") |
| `INTERPOLATED_MESH` is zero for every fine cell | `fds_fine_level.f90:162` |
| The wall loops that contain the `NIC>1` branches run over `1..N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS`, which is empty on a fine box; the `INTERPOLATED_BC` case needs `EXTERNAL_WALL(WALL_INDEX)` | `divg.f90:201`, `wall.f90:118` and `wall.f90:833-836` |
| Side data of a fine level: a face towards another box of the level or towards the coarser level is open; only the domain edge is a wall | `LevelRegistry.cpp:15-21` |
| A level above 0 must run with `EXTERNAL_GHOSTS_FILLED`; the `OMESH` route is same-level only | `GhostExchange.cpp:121`, `GhostExchange.H:43-46` |
| With `EXTERNAL_GHOSTS_FILLED` the neighbour-mesh branches of the velocity and wall routines are skipped | `velo.f90:528`, `1393`, `1882`, `2665`, `2910`; `wall.f90:312` (`ASSIGN_GHOST_VALUE`) |
| Ghost cells across the jump are filled from the coarser level by the hook installed in `BcStep::exchange` | `GhostExchange.H:40-46`, `GhostExchange.cpp:122-127` |
| Fine face velocity overwrites the coarse faces under it; the interface flux overwrite of the scalar fluxes runs between read-out and cell update; one dt, no subcycling | `TimeLoop.cpp:891-896`, `1019-1037`; `RegistryTransfer.cpp:113`, `172` |
| Fine levels are bound by `bind_level`, which creates the fine box objects | `TimeLoop.cpp:194-232` |

**Answer (a), stage execution: no.** No FDS `EXTERNAL_WALL` exists at a coarse/fine level interface in the stages, and the two `NIC>1` branches have no wall to run on. The same holds for box-to-box faces inside a fine level.

**Answer (a), whole route: not provable.** Two facts stop a flat "never built":
1. The set-up phase builds `NIC>1` walls whenever the FDS input itself carries meshes of different resolution (section 1.1), before the abort at `FdsAmr.cpp:30`. These walls are never used.
2. The step that should make a multi-resolution input into levels is not in the production path. `build_hierarchy_from_meshes` (IR-002: finer `&MESH` becomes a level-1 box; `Source/regrid_transport/Hierarchy.cpp`) is called once in the driver, from the two-level test harness (`TwoLevelRun.cpp:98`), with the level-0 mesh list and an `&AMR_REGION` text generated at `TwoLevelRun.cpp:75-79`. The gate input for that harness is a hand-made single-mesh form of `ns2d_16_int_1to2_refinement` (`Source/driver/tests/cases/ns2d_16_l0.fds`, header comment), not the committed input. The wiring note says the converter takes "the mesh list from the FDS set-up" (`Source/regrid_transport/notes/readf90-cmake-patch-list.md:13-15`). FDS reads every `&MESH` line, so a finer mesh would be an ordinary FDS mesh with ordinary `EXTERNAL_WALL` cells and a mesh number below `NMESHES`, while fine boxes are numbered above `NMESHES` (patch 0007, `Source/driver/patches/0007-mesh-fine-level-boxes.md`). `read.f90` has no `WITH_AMREX` block and no patch in `Source/driver/patches` or `upstream-patches` filters `&MESH` or `&AMR` lines.

### 1.3 Guards
| Guard asked for | Present? |
|---|---|
| Abort when a wall has `NIC>1` (D-065 Q1 "host abort guard") | **No.** No `NIC` reference in `Source/driver`; `wall.f90:897` and `divg.f90:219` carry no `WITH_AMREX` branch. |
| Abort on unequal cell size at level 0 | Yes, `FdsAmr.cpp:30` (message names the rule, not the pair). |
| Abort on overlap / off-lattice corner | Yes, `FdsAmr.cpp:62`, `53`. |
| Input check "finer `&MESH` without an `&AMR` line is an error, AMR mode is never inferred" | In `Hierarchy.cpp:175`, which the production path does not call. |
| Ratio other than 2 or 4 refused | In `Hierarchy.cpp:171`, same status. |
| Overlapping or embedded meshes refused | In `Hierarchy.cpp:226`, same status; in the production path by `FdsAmr.cpp:62`/`53`/`30`. |

## 2. The multi-mesh inputs (question (b))

### 2.1 Reconciliation of the count
The Verification rollup (`inventory/mesh_ratio_cases_rollup.csv`, regenerated from the tree above) lists **seven** inputs with a shared face of tangential ratio 2 or more: `ns2d_16_int_1to2_refinement`, `geom_stretched_grid`, `obst_coarse_fine_interface`, `dancing_eddies_uglmat_refine`, `duct_flow_uglmat_refine`, `race_test_1`, `race_test_4`. This is the "seven" of D-072 (f). The earlier review table has six rows because `race_test_1` and `race_test_4` share one. The four derived 2:1 inputs and the four embedded inputs are **not** among the seven: the derived inputs are V&V copies outside `Verification/`, and the embedded ones have overlapping meshes and therefore no shared faces in the table. Together: 7 + 4 derived + 4 embedded = 15 inputs.

Species count is read from `&REAC` and `&SPEC` lines of the input (not from a run). One `&REAC` makes more than one tracked species; a lone `&SPEC` for water vapour adds one species to air; a lone `&SPEC` for a background or lumped air species adds none. The two inputs with a lone `&SPEC` where this matters were checked against recorded set-up output in the earlier review (one tracked species).

### 2.2 Table
"Reaches branch in FDS" = shared face with ratio 2 or more (so `NIC>1` on the coarse side, `NIC` = product of the tangential ratios) and two or more tracked species. "Today" = what the production path does, by the code in section 1. "Intended" = the ruling in force.

| # | Input | Scope class (`vv/scope_case_list.csv`) | Meshes / ratio | Tracked species | Reaches branch in FDS | Today (code) | Intended AMR-mode status |
|---|---|---|---|---|---|---|---|
| 1 | `Adaptive_Mesh_Refinement/ns2d_16_int_1to2_refinement` | IN, Tier 1 (FR-016 gate) | 13 meshes, 1x2 / 2x1 (2-D) | 1 (`&SPEC` background) | no | set-up, then abort at `FdsAmr.cpp:30` (driver sweep row: `driver/notes/fft-thin-direction-cases.csv:4`) | converted to a level-0 domain plus a level-1 box with an `&AMR MAX_LEVEL=1` line added (`regrid_transport/notes/readf90-cmake-patch-list.md:18-20`); run today only through its single-mesh form with `--two-level-run` |
| 2 | `Pressure_Effects/obst_coarse_fine_interface` | IN, Tier 1 | 2 meshes, 2x2 | 1 (no `&REAC`, no `&SPEC`) | no | set-up, then abort at `FdsAmr.cpp:30` | converted to a refinement level (no `&AMR` line in the input; one has to be added) |
| 3 | `Pressure_Solver/dancing_eddies_uglmat_refine` | IN, Tier 1 | 4 meshes, 1x2 (2-D) | 1 (`&SPEC` lumped air) | no | set-up, then abort (`fft-thin-direction-cases.csv:142`) | converted to a refinement level |
| 4 | `Pressure_Solver/duct_flow_uglmat_refine` | UNCLEAR (thin OBST at the interface; FR-040 R3 open) | 8 meshes, 2x2 | 1 (no `&REAC`, no `&SPEC`) | no | set-up, then abort at `FdsAmr.cpp:30` | undecided until the FR-040 R3 ruling; FDS-only meanwhile |
| 5 | `Thread_Check/race_test_1` | IN (5:1 faces) | 6 meshes, 5x5 on five faces | 2 or more (`&REAC`, no `&SPEC`) | **yes** (`NIC` 25) | set-up builds `NIC=25` walls, then abort at `FdsAmr.cpp:30`; if the hierarchy builder were called, ratio 5 is refused at `Hierarchy.cpp:171` | AMR mode uses the 4:1 copy `vv-runs/inputs/A-35/race_test_1_r4` (`vv/test-plan.md:277`), through the hierarchy; the original stays FDS-only |
| 6 | `Thread_Check/race_test_4` | IN, Tier 2 optional | same layout as 5 | 2 or more | **yes** | as 5 | as 5, copy `race_test_4_r4` |
| 7 | `Complex_Geometry/geom_stretched_grid` | DEFERRED (`&GEOM`, `TRN*`: 31 `TRN*` lines, one `&GEOM`) | 5 meshes, 4x4 | no `&REAC`, no `&SPEC` | no | `TRN*` abort at `FdsSetup.cpp:31` (also before that the level-0 size rule) | FDS-only (deferred scope) |
| 8 | `A-46/level_set_fuel_model_1_2to1` | V&V derived (vv-runs/inputs) | 2 meshes, 2x2, no `&AMR` line | 2 or more (`&REAC` + `&SPEC`) | **yes** (`NIC` 4) | set-up builds `NIC=4` walls, then abort at `FdsAmr.cpp:30` | refinement level from the hierarchy, T2 against its own FDS baseline (`vv/test-plan.md:344-352`); the input needs an `&AMR` line |
| 9 | `A-47/bucket_test_1_2to1` | V&V derived | 5 meshes, 2x2 | 2 (air plus `&SPEC` water vapour) | **yes** | as 8 | as 8 |
| 10 | `A-47/cascadempi_obst_2to1` | V&V derived | 6 meshes, 2x2 | 2 (water vapour) | **yes** | as 8 | as 8 |
| 11 | `A-47/geom_sprk_mass_obst_2to1` | V&V derived | 4 meshes, 2x2 | 2 (water vapour) | **yes** | as 8 | as 8 |
| 12 | `Adaptive_Mesh_Refinement/ns2d_16_emb_1to1_refinement` | OUT (embedded, D-015) | 2 meshes, overlapping, equal size | 1 | no (no shared face) | abort at `FdsAmr.cpp:62` (or `53`) | FDS-only |
| 13 | `Adaptive_Mesh_Refinement/ns2d_16_emb_1to2_refinement` | OUT | 2 meshes, overlapping, 2x1x2 | 1 | no | abort at `FdsAmr.cpp:30` | FDS-only |
| 14 | `Adaptive_Mesh_Refinement/random_meshes` | OUT | 5 meshes, 8 overlapping pairs, 2x1x2 | 1 (no `&REAC`, no `&SPEC`) | no | abort at `FdsAmr.cpp:30` | FDS-only |
| 15 | `Pressure_Solver/dancing_eddies_embed` | OUT | 2 meshes, overlapping, equal size | 1 | no | abort at `FdsAmr.cpp:62` (or `53`) | FDS-only |

Other derived 2:1 inputs (`A-46/duct_flow_2to1`, `fan_test_2to1`, `zone_break_fast_2to1`, `zone_shape_2to1`, `porous_media_r2`) have no `&REAC` and no `&SPEC`, so one tracked species and the branches are not reached; their status is that of rows 8 to 11 apart from the species.

Row notes:
- The driver sweep (`Source/driver/notes/fft-thin-direction-check.md:38`) records the refused-at-level-0 outcome for `ns2d_16_*`, `random_meshes`, `dancing_eddies_embed` and `dancing_eddies_uglmat_refine`; it is a recorded sweep result which I did not re-run. Rows 2, 4, 5, 6, 7, 8 to 11 were not in that sweep; their outcome here follows from the code lines cited.
- Rows 12 and 15 have equal cell sizes: the abort is the overlap or lattice check, not the size check.
- "Refused" in the table means the abort in the production path. No input is converted today; the conversion exists only as a library and a test harness.

### 2.3 What D-072 (f) and D-065 Q1 need
- The condition "the AMR route never builds FDS `EXTERNAL_WALL`s at a coarse/fine jump" **holds for every stage of every fine level** and for level 0 of every input that reaches the time loop (because level 0 is equal-size, `NIC=1`).
- It **does not yet hold as a statement about the route end to end**, because a multi-resolution input is not converted: the FDS set-up sees it whole.
- Retiring L1485 (the species flux match, `wall.f90:897` and following) is safe for any input that reaches a stage. It is not enforced by a guard. Rows 5, 6 and 8 to 11 are the inputs on which baseline FDS runs the branch; they are FDS-baseline-only today and AMR-mode only through the hierarchy once it is wired.

## 3. Validation inputs (question (c))

Scope of the check: the V&V plan states validation against experiments is out of scope and names the validation inputs it uses (`vv/test-plan.md:74`, `628`, `650`, `requirements.md:49`, `485`): the Heskestad flame-height input, its derived coarse and fine copies, the McCaffrey fallback, and the FM_Burner inputs as FR-005 (v) area tests. All are present in the tree (`Validation/`) and the derived copies under `vv-runs/inputs/A-19`. Method: ratio columns from the rollup (parser over every `&MESH`, `MULT` expanded), species from `&REAC`/`&SPEC` lines. Not run.

| Input (plan use) | Meshes | Shared-face ratio | `&REAC`/`&SPEC` | Jump with two or more species? |
|---|---|---|---|---|
| `Heskestad_Flame_Height/Qs=1_RI=10` (NFR-032 level-0 source) | 1 | none | 1 / 0 | no (no face) |
| `A-19/Qs1_RI10_coarse_32x32x80` (level 0) | 1 | none | as above | no |
| `A-19/Qs1_RI10_fine_64x64x160` (uniform-fine reference) | 8, equal size | 1 | as above | no (equal size, `NIC=1`) |
| `McCaffrey_Plume/McCaffrey_14_kW_5` (fallback) | 27 | 1 | 1 / 1 | no |
| `FM_Burner` family A, 2 cm and 1 cm (8 inputs, FR-005 (v) incl. `CH4_2cm`) | 12 or 96 | 1 | 2 / 3 or 2 / 4 | no |
| `FM_Burner` family A, 5 mm (`C2H4`, `C3H6`, `C3H8`, `CH4`; 4 inputs; not named by the plan) | 208 | **2** | 2 / 3 or 2 / 4 | **yes if ever run** (`NIC=2` per tangential direction, up to 4) |
| `FM_Burner` family B, 5 mm (3 inputs; `&GEOM`, out of scope) | 180 | **2** | 1 / 0 | yes, but out of scope |
| `FM_Burner` family B, 2 cm and 1 cm (6 inputs; `&GEOM`) | 96 | 1 | 1 / 0 | no |

Result: **the validation inputs the plan uses have no coarse/fine face**, so they cannot reach the `NIC>1` branches in baseline FDS, and in AMR mode their refinement is made by `&AMR`/`&AMR_REGION` lines on an equal-size level 0 (`vv/test-plan.md:628`), not by differing `&MESH` sizes.

Wider bound, in case the plan's validation set grows: in the whole `Validation/` tree (3,132 inputs) 325 have a shared face of ratio 2 or more, 256 of them with at least one `&REAC` or `&SPEC` line (counts by folder: `USFS_Deep_Fuel_Beds` 111, `USN_Hangars` 26, `NIST_Pool_Fires` 23, `FM_SNL` 19, `Bluff_Body_Flows` 10, `Convection` 10, `NIST_NRC_Corner_Effects` 8, `FM_Burner` 7, `BST_FRS_wood_cribs` 6, `UMD_Line_Burner` 6, `CERTEC_Pool_Fires` 5, `McCaffrey_Plume` 5, `Sandia_Pool_Fires` 4, and 8 folders with 1 to 3; the `McCaffrey_Plume` five are not the fallback input above). The rollup also classes 11 Validation inputs as ratio 3, 11 as non-power-of-2 and 7 as direction-dependent (all outside the 2/4 rule). None of the 3,132 has an overlapping pair. All 325 stay FDS-only unless converted; the `&REAC`/`&SPEC` test is a lower bound on "two or more tracked species" for the `&SPEC` case, so 256 is an upper-bound count of reachable inputs.

Missing from this computer: nothing the plan names. The FM_Burner family B inputs reference geometry files that are absent (`vv/fm_burner_area_check.md` summary), which does not affect the mesh check.

## 4. Gap list for the Architect
1. **Where finer `&MESH` lines go.** Decide that finer `&MESH` entries (and `&AMR*` lines) are removed from what the FDS set-up reads, by a `WITH_AMREX` filter in `read.f90` or by a driver-side input rewrite, before `fds_setup(0)`. Without it, the set-up of every multi-resolution input builds `NIC>1` walls and numbers the finer meshes as level-0 FDS meshes, which clashes with the fine-box numbering above `NMESHES`.
2. **Wire the converter.** `main.cpp` has to call `parse_amr_params` and `build_hierarchy_from_meshes` on the input and give `assemble_level0` only the level-0 meshes; today only the two-level harness does that.
3. **Add the `NIC>1` guard.** D-065 Q1 retires L1485 "behind a host abort guard", but there is none. A guard at the end of `INITIALIZE_MESH_EXCHANGE_1` (abort when any `EWC%NIC>1`, `WITH_AMREX` only) would make the design requirement checkable and costs nothing at run time.
4. **Say how the `_r4` and derived 2:1 inputs get an `&AMR` line.** `Hierarchy.cpp:175` rejects unequal meshes without one, and none of the derived inputs (rows 5 to 11) has one.
5. **Level-0 `NIC=1` test.** An assertion or test that every level-0 `EXTERNAL_WALL` has `NIC=1` after `assemble_level0` would replace the code reading in section 1.1.
6. Documentation: the header comment of `LevelRegistry.H:8-11` still says a level above 0 has no FDS binding; `bind_level` (`TimeLoop.cpp:194-276`) now provides it.
