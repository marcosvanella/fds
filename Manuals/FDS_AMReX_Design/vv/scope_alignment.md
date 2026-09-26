# Scope alignment of the Verification suite (FR-006, D-033; A-41)

Owner: AMR V&V Lead · Draft 2026-09-25, first issue, **updated** after Spec Lead rulings (d) and (e) (D-035, A-46 (4), A-47) · Static text analysis only; nothing was run.
Source: FireX `36975d765f` (this repository, read-only), 941 `Verification/**/*.fds` inputs.
Requirement: **FR-006 "Verification-suite coverage of in-scope features"**, added in spec **v0.4.10** (owner decision **D-033**; V&V action **A-41**). `requirements.md` re-read, header v0.4.14 (FR-006 now records the case list and D-035).
Outputs: `docs/vv/scope_case_list.csv` (one row per input), this file, test-plan §5.8. Rerun with `python3 vv-runs/tools/scope_filter.py` (about 10 s). It re-checks every source citation below and prints `CITE MISMATCH` if a cited line moves.

## 1. Method

1. **Parsing.** `parse_fds_inputs.records()` follows FDS CHECKREAD. A record starts only on a line whose first non-blank character is `&`, runs across continuation lines, and ends at the first `/` outside quotes. Everything after that `/`, and every line not starting with `&`, is a comment. So `CSVF UVWFILE=…` without `&` (csvf_restart_a) and `#&GEOM …` (geom_terrain2) are ignored, as FDS ignores them. The script also:
   - inlines `&CATF OTHER_FILES` (read.f90:248, 284);
   - expands `&MULT` on `&MESH`/`&OBST`/`&HOLE` with FDS's offset formula (read.f90:645-657);
   - takes 2-D directions (one cell on every mesh) out of the resolution and thin tests.
2. **Feature flags and the source keywords behind them** (all re-verified at run time):

| Flag | Rule used | Source (FireX 36975d765f) |
|---|---|---|
| GEOM | any `&GEOM` record | geom.f90:1159 (namelist), 1189 (no GEOM → return) |
| CC_IBM | `&MISC CC_IBM=T`, or any `&GEOM` with T_END−T_BEGIN > 0 | read.f90:1746 (MISC keyword), cons.f90:732 (default F), geom.f90:2696-2698 (GEOM sets CC_IBM=T unless setup-only), ccib.f90:7197 (CC_IBM=T without GEOM is an error) |
| HT3D | `&SURF HT3D=T` (there is no HT3D keyword on OBST or MATL) | read.f90:8042 (SURF keywords HT3D, VARIABLE_THICKNESS), 8462-8465 (HT_DIM=3, SOLID_HEAT_TRANSFER_3D) |
| VARIABLE_THICKNESS (sub-flag) | `&SURF VARIABLE_THICKNESS=T` | read.f90:841, 849 (grouped with HT3D OBSTs) |
| THIN_OBST | any `&OBST` with zero extent in a non-2-D direction; or extent < 0.5·dx of the mesh it is in (FDS's THIN test); or `&MISC THICKEN_OBSTRUCTIONS=T`; or `&OBST THICKEN=T`; or a zero-thickness `&HOLE`. Sub-counts: thin OBSTs touching ≥ 2 meshes, touching meshes of different resolution, and those whose THIN status would change under a ratio-2/4 refinement (0.125·dx ≤ t < 0.5·dx) | read.f90:1765 (THICKEN_OBSTRUCTIONS), 10767 and 10890 (OBST THICKEN), 11315-11317 (`OB%THIN`: snapped to 0 cells and length < 0.5 dx), 11738 (HOLE); FR-040 R1 snap rules |
| HVAC | any `&HVAC`. Sub-flags: node vents (`VENT_ID`/`VENT2_ID` → `&VENT` XB) in ≥ 2 meshes, or in meshes of different resolution | hvac.f90:151, 157 |
| ZONE | explicit `&ZONE`, **or** automatic pressure zones: FireX makes every gas region not reached from an OPEN boundary a new P_ZONE unless `NO_PRESSURE_ZONES` (main.f90:2620, 2652-2682, `N_ZONE = N_ZONE + 1` at 2679). ZONE=1 if there is an explicit ZONE, ≥ 2 automatic zones, a sealed zone next to an OPEN-connected region, or a sealed zone that an OBST/HOLE with DEVC_ID/CTRL_ID can open. A domain that is one sealed zone and nothing else is reported separately (`sealed_domain_single_zone`, background pressure only). Leak paths: `&ZONE LEAK_AREA`, `&SURF LEAK_PATH/LEAK_PATH_ID` | read.f90:13616 (ZONE), 8045 (SURF LEAK_PATH), 1757 (NO_PRESSURE_ZONES); main.f90:2620, 2679 |
| LEVEL_SET | `&MISC LEVEL_SET_MODE>0`. SURF `VEG_LSET_*` without it is an FDS error | read.f90:1754, 1981-1983 (LEVEL_SET_MODE>0 also sets NO_PRESSURE_ZONES), 8343-8345 |
| STRETCHED | any `&TRNX/&TRNY/&TRNZ` | read.f90:1000-1002 |
| MULTI_MESH / MULTI_RES / OVERLAP / RATIO_NOT_2_4 | mesh count after MULT; >1 distinct cell size (2-D direction ignored); `n_overlap_pairs` and `worst_class` from the Legacy Mapper A-33 rollup `docs/inventory/mesh_ratio_cases_rollup.csv` | read.f90:524, 645 |
| HYPRE_DEVICE_RUN | explicit `&PRES HYPRE_DEVICE_RUN` (default T, only active in `WITH_HYPRE_DEVICE` builds). **Ignored for the class** (D-035): noted only | read.f90:10069, cons.f90:569, pres.f90:1177-1178 |
| GPU (rule, D-035 / ruling (d)) | An input is OUT on GPU grounds only if it cannot run on a CPU build, i.e. uses a keyword in `GPU_ONLY_KEYWORDS` (empty at the pin: `HYPRE_DEVICE_RUN` is the only GPU keyword in read.f90 and it is CPU-inert). Ranks per GPU is the job-environment variable `FDS_RANKS_PER_GPU`, not an input keyword. GPU runs of the AMReX code are governed by NFR-043/NFR-047 | main.f90:5108 (`FDS_RANKS_PER_GPU` read from the environment) |
| TUNNEL_PRECONDITIONER | `&PRES TUNNEL_PRECONDITIONER=T` (still open under Q4 per D-033) | read.f90:10068, 10177 |
| RESTART | `&MISC RESTART=T`; producer = same-folder input whose CHID = RESTART_CHID | read.f90:1762 |
| PARTICLES | any `&PART` | read.f90:5933 |
| MISSING_FILES | quoted file values of `OTHER_FILES`, `UVW_FILE`, `CSVF UVWFILE/TMPFILE/SPECFILE`, `EXTERNAL_FILENAME`, `BULK_DENSITY_FILE`, `BINARY_FILE`, `TERRAIN_IMAGE` that do not exist relative to the case folder. Files a sibling run writes are noted rather than counted as missing | read.f90:248, 284, 12963, 17325 |
| VTK | `&DUMP WRITE_FORMAT='VTK'` | read.f90:2382 |
| STORE_SPECIES_FLUX | DEVC/SLCF/BNDF quantity `ADVECTIVE|DIFFUSIVE|TOTAL MASS FLUX X|Y|Z` or `TOTAL MASS FLUX WALL` | read.f90:17011-17014 |
| machine limits | cells ≤ 4 M (16 GB development machine, ~1.2 GB per M cells); FDS ranks needed = max `MPI_PROCESS`+1 ≤ 8 (FDS runs fewer ranks than meshes, read.f90:571) | |

3. **Classification** (first match wins):
   - **DEFERRED:** GEOM, CC_IBM or HT3D, or a restart of a DEFERRED producer.
   - **OUT:**
     - inputs that cannot run on a CPU build (none at the pin; `GPU_Tests/` and explicit `HYPRE_DEVICE_RUN=T` are no longer OUT, D-035 and §2);
     - the four embedded/overlapping cases (decided FDS-baseline-only);
     - TRN-stretched grids (IR-002, D-030);
     - VTK output (no HDF5 build, D-019);
     - missing external files;
     - more than 4 M cells (none).
   - **UNCLEAR:** one-line reason, from the rules in §5.
   - **IN:** everything else.
4. **Gates, tier and tolerance class for IN cases:**
   - Every IN case is gated on "completes": normal exit, no NaN or trap (test-plan §3.4). Setup-only inputs (T_END=0) are gated on completing setup plus the G1 message diff.
   - Cases in the 184-row inventory keep its tier and `tclass` verbatim. The `tclass` is still v0.2 text; D-022 re-derivation is pending (§2).
   - For the rest, a proposal:
     - Verification-Guide dataplot row → T2 with that metric;
     - CHID referenced by a `Utilities/Python/scripts/*.py` → T2 with the script metric (heuristic match);
     - listed in `FDS_Cases.sh` with no FDS metric → T2 against the FDS baseline;
     - not run by FDS (absent or commented) → completes only.
   - Tier **G4** (the full-suite batch, test-plan §5), or **T2** where the case exercises a must-work feature across meshes.
   - `needs_remesh` names the A-35 copies.

## 2. Counts

| Class | Inputs | Main reasons |
|---|---|---|
| **IN** | **705** | 157 are in the 184-row inventory, 548 are new. 127 multi-mesh; 5 multi-resolution (`ns2d_16_int_1to2_refinement`, `obst_coarse_fine_interface`, `dancing_eddies_uglmat_refine`, `race_test_1/_4` via A-35); 134 with particles. Includes 3 of the 4 `*_hypre` inputs and the 9 `GPU_Tests/` inputs (D-035 rule, below) |
| **DEFERRED** | **175** | 128 `&GEOM` (111 in `Complex_Geometry/`), 32 HT3D (`Heat_Transfer/`), 15 `&MISC CC_IBM=T` without GEOM (`Complex_Geometry/*_cc_exp_*`, `saad_CC_explicit_*`). GEOM and HT3D never occur together |
| **OUT** | **26** | 16 TRN-stretched; 4 embedded/overlapping (decided); 4 `VTK/`; 2 missing files (`Controls/ext_heartbeat_std_curve`: `ext_heartbeat_mass_flux.csv` comes from an external driver; `WUI/bulk_density_file`: `../../../cad/…/canopy_foliage.bdf` is not in the tree) |
| **UNCLEAR** | **35** | 12 VARIABLE_THICKNESS, 11 CYLINDRICAL, 7 TUNNEL_PRECONDITIONER (now including `dancing_eddies_ulmat_hypre`), 2 multi-mesh zone leakage, 1 thin OBST at a coarse/fine interface, 1 per-mesh CSVF restart, 1 disjoint multi-resolution meshes (§5) |
| Total | 941 | |

Tolerance-class sources for the 705 IN cases: 157 from the inventory, 379 with a Verification-Guide/script metric (T2), 155 T2 vs FDS baseline (no FDS metric), 14 completes only.

Proposed tiers for the 548 new IN cases: 538 G4 and 10 T2. The 10 T2 promotions are `HVAC/{HVAC_leak_exponent_2, fan_test, fan_test_2, qfan_multi, qfan_test}`, `Heat_Transfer/back_wall_test`, `Pressure_Effects/zone_break_fast_uglmat`, `Pressure_Solver/{duct_flow_ulmat, stairwell}` and `Pyrolysis/methanol_evaporation`.

**Change after ruling (d) (D-035, rerun).** The script now has an explicit GPU rule (`GPU_KEYWORDS_IGNORED_ON_CPU`, `GPU_ONLY_KEYWORDS` in `scope_filter.py`) instead of the two OUT reasons. Class moves:
- `Pressure_Solver/dancing_eddies_uglmat_hypre`, `duct_flow_uglmat_hypre`, `duct_flow_ulmat_hypre`: OUT → **IN** (inventory T2, G6/FR-038 unchanged).
- `Pressure_Solver/dancing_eddies_ulmat_hypre`: OUT → **UNCLEAR**. It also sets `TUNNEL_PRECONDITIONER=T` (input line 10), so it follows Q4, as D-035 says.
- `GPU_Tests/test_gpu` and the 8 `GPU_Tests/HYPRE_GPU_SCALING/test_*` inputs: OUT → **IN**. Checked: none sets `HYPRE_DEVICE_RUN`. They are plain 1 M-cell propane-fire inputs (`SOLVER='UGLMAT'` with `MAX_PRESSURE_ITERATIONS=1`, or `'UGLMAT HYPRE'` for `test_gpu`), 4-32 meshes with no `MPI_PROCESS`, so FDS runs them on ≤ 8 ranks (read.f90:571). The `_RS2/_RS4/_RS8` files are byte-identical to their `NOFRPG`/plain twins apart from `CHID`. The ranks-per-GPU layout comes from `FDS_RANKS_PER_GPU` in the job scripts (`submit_all.sh`; main.f90:5108). They keep their inventory row (T3, `tclass` n/a, no FDS pass criterion), so their gate is **completes only**.
- Cost: the inventory models each at ~2,600 s wall on the development machine (1 M cells, 209-417 steps), about 6.5 h for all 9. Recommendation: run them in G4 (Phase 2 exit, Phase 10) only, and count the four RS twins once for pass/fail (same input).
- This goes beyond D-035, which names only the four `*_hypre` inputs. A-46 (4) expected 696 IN / 35 OUT / 35 UNCLEAR. The 9-case difference is the `GPU_Tests/` inputs, moved by the same rule on the V&V Lead's instruction. The Spec Lead should acknowledge it (or keep them OUT by adding `GPU_Tests/` back as an explicit exception in the script).

Other notes on the IN set:
- `stairwell` needs 10 FDS ranks (`MPI_PROCESS`), so its baseline needs a derived MPI_PROCESS ≤ 7 copy, as for `random_obstructions_fft`.
- Six IN cases set STORE_SPECIES_FLUX: `mass_balance_gas_volume`, `mass_balance_reac`, `mass_balance_reac_2`, `mass_flux_comparison`, `mass_flux_wall_yindex`, `mass_flux_wall_zindex`. They stay out of the D-022 kernel comparison until A-34 closes (test-plan §5.5); the gate column says so.
- **Assumption behind 40 IN cases:** the union of their meshes is not a box (disjoint sub-tests, L-shaped buildings; column `nonbox_domain`). IN assumes AMR mode fills the gaps with solid cells on level 0, which matches FDS's exterior-wall treatment. IR-002 does not say so. If the Integration Lead does not confirm it, these 40 move to UNCLEAR, including inventory cases `simple_duct`, `isentropic`, `hallways`, `device_restart_a/b/base_case` and `1_step_2_step_compare`.

## 3. Coverage of the four must-work features (IN set)

"Multi-mesh" is the subset where the feature crosses a mesh-to-mesh face:
- THIN_OBST: a thin OBST touches ≥ 2 meshes.
- HVAC: node vents lie in ≥ 2 meshes.
- ZONE and LEVEL_SET: the case has pressure zones / level set and its meshes share faces (A-33 `n_shared_faces` > 0). Disjoint sub-test meshes do not count.

"Resolution change" means the feature touches meshes of different cell size.

| Feature | IN | IN multi-mesh | IN at a resolution change | UNCLEAR | DEFERRED | OUT |
|---|---|---|---|---|---|---|
| Thin obstructions | 44 | **10**: `Pressure_Solver/duct_flow`, `duct_flow_ulmat`, `duct_flow_uglmat_pardiso`, `duct_flow_uglmat_hypre`, `duct_flow_ulmat_hypre` (28 thin duct walls across 8 meshes), `random_obstructions_fft` (104), `Pressure_Effects/zone_break_fast_uglmat`, `zone_break_fast_uglmat_hypre`, `Heat_Transfer/back_wall_test`, `Pyrolysis/methanol_evaporation` | **0** | 4 (`duct_flow_uglmat_refine`, `zone_shape_2`, `box_burn_away10/11`) | 8 (HT3D) | 0 |
| HVAC | 45 (36 in `HVAC/`; plus `Pressure_Solver/duct_flow`, `duct_flow_ulmat`, `duct_flow_uglmat_pardiso`, `duct_flow_uglmat_hypre`, `duct_flow_ulmat_hypre`, `Flowfields/jet_fan`, `simple_duct`, `Atmospheric_Effects/stack_effect`, `Aerosols/aerosol_agglomeration_2`) | **5**: `HVAC/fan_test`, `fan_test_2`, `qfan_test`, `qfan_multi`, `HVAC_leak_exponent_2` | **0** | 2 (`duct_flow_uglmat_refine`, `zone_shape_2`) | 2 (`HVAC_geom`, `leak_geom`) | 0 |
| Pressure ZONEs | 54 (3 explicit `&ZONE`: `door_crack`, `leak_test`, `leak_test_2`; 15 with zone opening/closing) | **15**: `Pressure_Effects/zone_shape`, `zone_break_fast_uglmat`, `zone_break_fast_uglmat_hypre`, `Pressure_Solver/ulmat_2zone`, `random_obstructions_fft`, `stairwell`, `Restart/device_restart_a/b/base_case`, `restart_ulmat_a/b`, `HVAC/fan_test`, `fan_test_2`, `qfan_test`, `HVAC_leak_exponent_2` | **0** (partial: `obst_coarse_fine_interface` and `dancing_eddies_uglmat_refine` are one sealed zone spanning a 2:1 interface; `sealed_domain_single_zone`) | 2 (`zone_shape_2`, `HVAC_leak_exponent`) | 8 (e.g. `leak_test_3`, `leak_test_4`, `sphere_leak`, `cascadempi`) | 1 (`divergence_test_1`, TRN) |
| Level-set wildfire | 11 (7 with LEVEL_SET_MODE=1, uncoupled, no gas solve; 4 with mode 4: `LS4_ember_yield`, `LS4_ember_ignition`, `ground_vegetation_drag`, `level_set_fuel_model_1`) | **0** (`LS4_ember_ignition` and `ground_vegetation_drag` have 2 and 3 disjoint meshes; the fire never crosses a mesh face) | **0** | 0 | 5 (`LS_ellipse_*_30deg` and `geom_ls_restart_*` use `&GEOM` terrain) | 0 |

Full per-case lists are in the CSV (filter `cls=IN` and the feature column; `*_multimesh`, `thin_at_res_change`, `hvac_at_res_change`).

Also relevant to refinement: 9 IN cases have sub-cell OBSTs that are THIN on their mesh but would stop being THIN at ratio 2 or 4 (`thin_changes_with_refinement`): `device_restart_a/b/base_case` and the six `Species/hrrpuv_reac_*`. Under FR-040 R1 those OBSTs snap differently per level, and R3 has to grow the fine region around them.

## 4. Gaps

1. **None of the four must-work features has an IN Verification case at a resolution change.** The one input with thin OBSTs across a coarse/fine interface, `Pressure_Solver/duct_flow_uglmat_refine`, is UNCLEAR (§5), and its HVAC vents sit in the coarse meshes only.
2. **Level set has no multi-mesh case and no refinement case.** Every level-set input is single-mesh or uses disjoint meshes. The two 30° ellipse cases and the level-set restart set are DEFERRED (GEOM terrain).
3. **HVAC has no case with vents at different resolutions.** It has five multi-mesh cases (vents in 2 or 8 same-resolution meshes).
4. **Pressure ZONEs:**
   - no multi-zone case at a resolution change;
   - no explicit-`&ZONE` case among the IN multi-mesh cases (the two multi-mesh explicit-ZONE inputs are UNCLEAR);
   - the only zone-breach cases across meshes are `zone_shape` and `zone_break_fast_uglmat[_hypre]`, all same-resolution.
5. **Thin obstructions:** 10 multi-mesh cases (5 are `duct_flow` solver variants of one geometry), none at a resolution change (see 1).
6. Features with no Verification case at all: none of the four. Level set lacks coupled multi-mesh coverage.

Recommendation: FR-006 says the V&V Lead chooses the refined variants (A-41). Derived refined variants under `vv-runs/inputs/A-41/`, in the A-35 style (a remeshed copy, the original kept as FDS baseline), would close gaps 1-5. Proposed first set:
- (a) `duct_flow` → a 2:1 variant with the interface kept off the thin duct walls (or `duct_flow_uglmat_refine` once §5 item 5 is resolved);
- (b) `HVAC/fan_test` with one room refined 2:1 (HVAC vents at a resolution change);
- (c) `Pressure_Effects/zone_shape` with one room refined 2:1 (zone breach across a level jump);
- (d) `WUI/LS4_ember_yield` or `level_set_fuel_model_1` as 2 abutting meshes at 2:1 (level-set front crossing a level jump);
- (e) `zone_break_fast` with a static patch around the breach.

Each needs an FDS baseline run of the variant itself, because AMR results are never compared with the originals (as for A-35).

## 5. UNCLEAR cases (35) and how each gets resolved

Owner's Q4 is pending (A-11). VARIABLE_THICKNESS, CYLINDRICAL and TUNNEL_PRECONDITIONER cases stay UNCLEAR until the project owner answers (Spec Lead ruling).

**Open points with the Chief Architect and the Pressure Solver Lead.** Two UNCLEAR groups do not wait for Q4. They wait for rulings that are open with the Chief Architect and the Pressure Solver Lead:
- **FR-040 R3:** what a *static* FDS multi-resolution input does when a thin OBST's snapped extent differs between levels at the interface: grow the fine region, stop at setup, or accept as given. Affects `duct_flow_uglmat_refine`, and every A-46 refined variant that puts an OBST near a level jump.
- **FR-034 leak-area summation:** how `&ZONE LEAK_AREA` / `SURF LEAK_PATH` leak area is summed when zone walls lie in several boxes and levels. Proposed: sum over uncovered wall faces only, with the exact sum per FR-005 (ii). Affects `zone_shape_2` and `HVAC_leak_exponent`.

V&V reruns the script after each ruling (A-46 (4)).

| Group | Cases | Reason | Resolved by |
|---|---|---|---|
| VARIABLE_THICKNESS (12) | `Fires/box_burn_away1` (T2, FR-042 anchor), `box_burn_away2`-`6`, `8`, `10`, `11`, `box_burn_away_2D`, `box_burn_away_2D_residue`, `Heat_Transfer/ht1d_pile` | The 1-D solid takes its thickness from the snapped OBST and shares the HT3D obstruction-grouping code (read.f90:841-849). D-033 defers "HT3D (3-D solid heat transfer)" and does not name it. With refinement the snapped thickness depends on the level (FR-040 R1) | Project owner/Spec Lead: is VARIABLE_THICKNESS part of deferred HT3D? If not, IN, plus an FR-040 check that the thickness comes from level 0 |
| CYLINDRICAL (11) | `Energy_Budget/test_hrr_2d_cyl` (T2), `Flowfields/cyl_test_1`-`4`, `helium_2d_isothermal` (FR-036 secondary case), `Heat_Transfer/insulated_steel_pipe_2d`, `Radiation/droplet_absorption_cyl`, `plate_view_factor_cyl_{30,60,100}` | D-033: "cylindrical meshes still open under Q4"; IR-002 rejection is an ASSUMPTION | Q4 answer. If rejected: OUT (FDS-only) plus an FR-004 negative test |
| TUNNEL_PRECONDITIONER (7) | `Pressure_Solver/dancing_eddies_default` (T1, FR-002 anchor, FR-033 acceptance case), `dancing_eddies_tight`, `dancing_eddies_ulmat`, `dancing_eddies_ulmat_hypre` (T2, G6/FR-038; D-035), `tunnel_demo`, `Flowfields/divergence_test_3` (T1), `Pressure_Effects/thick_orifice_5cm` | D-033: `TUNNEL_PRECONDITIONER` still open under Q4. It is an FFT-solver option (read.f90:10177-10184; it needs meshes in a single row along x); the AMReX pressure path (FFT::Poisson/MLMG) has no equivalent | Q4 answer. Proposed: IN, with the keyword accepted and ignored in AMR mode (the whole-run T2 comparison already allows a different solver) |
| Zone leakage over several meshes (2) | `Pressure_Effects/zone_shape_2` (T2), `HVAC/HVAC_leak_exponent` | `&ZONE LEAK_AREA` is spread over zone walls in several meshes. With levels, covered coarse wall faces must not add leak area twice; FR-034 does not say how | **Open with the Chief Architect and Pressure Solver Lead**: FR-034 addendum (leak area summed over uncovered wall faces, exact sum per FR-005 (ii)); then IN |
| Thin OBST at a coarse/fine interface (1) | `Pressure_Solver/duct_flow_uglmat_refine` (T2) | 4 thin duct walls cross the z = 3.0 interface between the 16³ and 32³ meshes. FR-040 R3 (proposed) grows the fine region around an OBST whose snap differs between levels, or stops. For a static FDS multi-resolution input the expected outcome (grow, error, or accept) is undefined | **Open with the Chief Architect and Pressure Solver Lead**: FR-040 R3 ruling for static FDS inputs. Recommended: IN with the static hierarchy as given, since it is the suite's only thin-OBST-at-level-jump case (gap 1) |
| Per-mesh CSVF restart (1) | `Restart/csvf_restart_b` | Reads `csvf_restart_{uvw,tmp,spec}_t1_m1.csv` written by `csvf_restart_a`: one file per FDS mesh. Mapping per-mesh CSV fields onto AMReX boxes is not specified (FR-080/IR-001) | Integration Lead; otherwise OUT (FDS-only) |
| Disjoint multi-resolution meshes (1) | `Sprinklers_and_Sprays/porous_media` (T2) | 3 separated meshes at 1:2:4 with no shared faces: not a uniform level 0 (IR-002) and not a nested hierarchy | A-35-style remeshed copy (three uniform single-mesh sub-cases or one uniform mesh); original stays the FDS baseline |

## 6. Changes to the 184-case inventory

The 180 committed-input rows map to 179 inputs (`soborot_superbee_square_wave_128` is on two rows). The 4 derived rows (`VV/…_uglmat`, `VV/restart_test1_continuous`, `VV/test_8mesh_NOFRPG_short`, the determinism row) inherit from their sources, all IN (the `GPU_Tests` source is IN since D-035). Status after this filter: **157 IN, 2 DEFERRED, 9 OUT, 11 UNCLEAR** (first issue: 145 / 2 / 22 / 10).

| Case | Tier | New class | Reason |
|---|---|---|---|
| `Sprinklers_and_Sprays/cascadempi` | T2, anchor (FR-050/051/052) | **DEFERRED** | `&GEOM` → CC_IBM (lines 24-29, GEOM half only). Replacement anchors: §8 |
| `Heat_Transfer/ht3d_energy_conservation_4` | T2-opt | **DEFERRED** | HT3D. A-41 notes it may stay in uniform-mode FR-003 regression; not in the FR-006 list |
| `ns2d_16_emb_1to2_refinement` (T1), `ns2d_16_emb_1to1_refinement` (T2), `random_meshes` (T1), `dancing_eddies_embed` (T2) | | OUT | Embedded/overlapping: FDS-baseline-only (decision in force). The inventory rows stay for the FDS baseline; they are not FR-006 cases |
| `Flowfields/divergence_test_1` | T2 | OUT | TRNX stretched (IR-002). The inventory had it as an ordinary T2 case |
| `Pressure_Solver/{dancing_eddies_uglmat_hypre, duct_flow_uglmat_hypre, duct_flow_ulmat_hypre}` | T2 (G6/FR-038) | IN (unchanged status; was OUT in the first issue) | D-035: `HYPRE_DEVICE_RUN` ignored on CPU builds |
| `Pressure_Solver/dancing_eddies_ulmat_hypre` | T2 (G6/FR-038) | UNCLEAR | D-035 IN rule, but `TUNNEL_PRECONDITIONER=T` (Q4) |
| `GPU_Tests/*` (9) | T3 | IN (unchanged status; was OUT in the first issue) | CPU-runnable, same rule (§2); completes only; beyond D-035, Spec Lead to acknowledge |
| `VTK/*` (4) | T3 | OUT | VTK/HDF5 (D-019). FR-074 still uses them once an HDF5 build exists |
| `Fires/box_burn_away1` | T2, anchor (FR-042) | UNCLEAR | VARIABLE_THICKNESS |
| `Energy_Budget/test_hrr_2d_cyl` | T2 | UNCLEAR | CYLINDRICAL (Q4) |
| `Flowfields/divergence_test_3` (T1), `Pressure_Solver/dancing_eddies_default` (T1, FR-002 anchor, FR-033 acceptance case), `dancing_eddies_tight` (T2), `dancing_eddies_ulmat` (T2), `tunnel_demo` (T2-opt) | | UNCLEAR | TUNNEL_PRECONDITIONER (Q4) |
| `Pressure_Effects/zone_shape_2` | T2 | UNCLEAR | zone leakage over 8 meshes |
| `Pressure_Solver/duct_flow_uglmat_refine` | T2 | UNCLEAR | thin OBST across the 2:1 interface |
| `Sprinklers_and_Sprays/porous_media` | T2 | UNCLEAR | disjoint 1:2:4 meshes |

`race_test_1` (T2) and `race_test_4` (T2-opt) stay IN with `needs_remesh` = A-35 `_r4` copies.

The inventory excluded HVAC and WUI inputs (`case_inventory.md` §7). Of those, 36 `HVAC/` inputs and 11 level-set inputs (`WUI/`) are now IN with proposed gates.

## 7. Recommendations

1. Put this list into the test plan as the FR-006 case list (test-plan §5.8). Keep `case_inventory.csv` as the tiered run list until the refresh job merges the `cls`/`gate` columns. Do not edit it here.
2. Q4 (project owner, pending) moves 30 of the 35 UNCLEAR cases:
   - VARIABLE_THICKNESS vs deferred HT3D (12);
   - CYLINDRICAL (11);
   - TUNNEL_PRECONDITIONER (7, including the FR-002 anchor `dancing_eddies_default`, Tier 1 `divergence_test_3` and `dancing_eddies_ulmat_hypre`).
   FR-040 R3 and FR-034 are open with the Chief Architect and Pressure Solver Lead (3 cases, §5).
3. *(Closed by D-035.)* `HYPRE_DEVICE_RUN=T` inputs are IN with the keyword ignored. The same rule also moves the 9 CPU-runnable `GPU_Tests/` inputs to IN; the Spec Lead should acknowledge this (§2).
4. FR-052 (and the FR-050/051 roles of `cascadempi`): adopt the derived `VV/cascadempi_obst` as anchor (§8, A-47).
5. Write the §4 refined variants (a)-(e) under `vv-runs/inputs/A-41/`, each with its own FDS baseline, before Phase 3 exit. They are the only way to cover the must-work features at a level jump.
6. Confirm the non-box-domain assumption (§2, 40 cases) with the Integration Lead. This is an IR-002 addendum: gaps between meshes become solid level-0 cells.
7. The automatic-zone detection is a static voxel flood fill (heuristic). Check it against the zone count FDS prints in each `.out` during the G1 setup-only sweep, and correct the ZONE column there.

## 8. Particle anchors without GEOM: FR-052, and the FR-050/051 roles of `cascadempi` (ruling (e), A-47)

`cascadempi` is DEFERRED with the rest of `&GEOM` (D-033). FR-052 stays in scope. The inventory also used `cascadempi` as an anchor for FR-050 and FR-051 (`case_inventory.csv`, req `FR-050,FR-051`), and FR-051's verification line names it. So all three need a non-GEOM case.

The requirements (v0.4.14, lines 247-258) ask for:
- **FR-050:** particles move between boxes, levels and ranks, and through regrids, with no loss or duplication. Tested by a particle ledger and an attribute-integrity check; no case is named.
- **FR-051:** particle-gas exchange deposits to and samples from the finest level covering the particle. Named cases: `energy_budget_particles`, `cascadempi` with refinement.
- **FR-052:** particle-wall interaction (deposition, surface motion) on refined levels. Wanted: droplets land on, run over and drip off box tops across several meshes, with a refined-patch boundary cutting a box top and a regrid while droplets are attached.

**The FireX code involved** (all at `36975d765f`; lines checked while writing):
- `type.f90:399`: `LP%WALL_INDEX`, "if liquid droplet has stuck to a wall".
- `part.f90:2400-2402`: a droplet hitting a SOLID wall cell stores `LP%WALL_INDEX`.
- `part.f90:2444`, `2456-2465`: an attached droplet moves down or along the surface at `VERTICAL_VELOCITY` / `HORIZONTAL_VELOCITY`.
- `part.f90:2509-2551`: a droplet leaving the surface (off an edge, under the solid, drip-off) re-looks up or clears `WALL_INDEX`.
- `part.f90:3098-3099`: accumulated water (`A_LP_MPUA`) is stored on the mesh-local wall cell (`B2`).
- `read.f90:6052`: droplets with `SPEC_ID` default to `ADHERE_TO_SOLID`; `6265-6271`, `6461-6463` set the defaults 0.5 and 0.2 m/s.
- `read.f90:8743-8744`, `9164-9165`: SURF `ALLOW_SURFACE_PARTICLES` (default T) and `ALLOW_UNDERSIDE_PARTICLES` (default F).
- `read.f90:17217`: an `AMPUA` DEVC switches on `ACCUMULATE_WATER`.
- Metrics come from `Utilities/Python/FDS_verification_script.py`, reading `Utilities/Python/FDS_verification_dataplot_inputs.csv`; row numbers below are file lines.

**Search.** Among the 705 IN inputs, candidates were those with `&PART`, liquid droplets on walls, or `AMPUA`/`MPUA`/`ADHERE_TO_SOLID`/`HORIZONTAL_VELOCITY`/`ALLOW_*_PARTICLES`. Results:
- IN, droplets on walls: `Sprinklers_and_Sprays/{cascade, bucket_test_1..4, flow_rate, flow_rate_2, water_evaporation_6a, _7}`, `Fires/spray_burner`, `Pyrolysis/ice_cube`, `WUI/hot_rods`.
- DEFERRED: `cascadempi`, `geom_sprk_mass`, `Complex_Geometry/geom_particle_cascade_2`.
- OUT: `VTK/cascade`.
- `Aerosols/*deposition*` is species (aerosol) deposition, not Lagrangian particles, so it does not exercise `WALL_INDEX`.

**No committed IN input has droplets moving over OBST surfaces on more than one mesh.** `bucket_test_1` is the only committed multi-mesh wall-deposition case, and it deposits on the floor only. The multi-mesh OBST case therefore has to be derived. It is easy, because `cascadempi` is really two independent halves: an OBST half and a GEOM half on separate meshes.

### FR-052 (particle-wall interaction)

| Rank | Case | Why | Meshes / cells | Runtime guess (development machine, model ±3×) | Metric | Evidence |
|---|---|---|---|---|---|---|
| **1 (anchor)** | **`VV/cascadempi_obst.fds`** (derived; proposed for `vv-runs/inputs/A-47/`) | `cascadempi` with the GEOM half deleted, so it is the exact test FR-052 was written around. Six OBST boxes (x 0.5-3.5 and 5.5-8.5, three levels). The `OBST-MultiMesh` pair (x 5-7, 7-9) splits the three right-hand box tops at x = 7.0, and nozzle N-2 sits on that face (x = 7.0). Droplets land on, run across (surface motion across a mesh face) and drip off box tops (`ALLOW_UNDERSIDE_PARTICLES=T`). **2:1 variant:** replace the pair with a coarse `IJK=8,16,56` (x 5-7) and a fine `IJK=16,32,112` (x 7-9) mesh. The level jump then cuts every right-hand box top. Box faces lie on the 0.125 m fine grid, so the snapped extents agree on both levels (FR-040 R3 not triggered). Ratio 2 (FR-010). For the "no attached particle changes attachment" check, add an `MPUA` surface-integral DEVC over each right-hand box top (derived-copy addition, T2 vs its own FDS baseline) | 3 meshes (1 + 2 sharing x = 7), 28,672 cells. 2:1 variant: 78,848 cells (fine part can be MULT-split for ranks) | ~1 min on 3 ranks (inventory model for `cascadempi`: 57,344 cells, ~1,334 steps, 58 s on 6 ranks). 2:1 variant: ~3 min with the fine mesh split into 4 (~2,700 steps), ~12 min unsplit | Committed rows 74-75 (`OBST-S-water mass`, `OBST-M-water mass` vs `Sprinklers_and_Sprays/cascadempi.csv` "Expected", Relative Error, area, tol 0.02). The DEVC IDs are kept, so the rows apply with the CHID changed (a V&V-side copy of the two rows). Status T2 | Keep `cascadempi.fds` lines 6, 9, 12 (meshes, MULT), 18 (`ALLOW_UNDERSIDE_PARTICLES=T`), 21-22 (OBST + BOX-MULT), 32 (droplet `PART` with `SPEC_ID`, so it adheres), 40-41 (N-1, N-2 at x = 7.0), 47-55 (OBST AMPUA / vapour DEVCs). Delete 7, 10 (GEOM meshes), 19 (GBOX), 24-29 (`&GEOM`), 42-43 (N-3, N-4), 57-65 (GEOM DEVC/CTRL) |
| 2 (secondary, already R-7) | `VV/geom_sprk_mass_obst.fds` (derived) | Sprinkler over a 2 × 2 m OBST plate with slow surface motion (`HORIZONTAL_VELOCITY=0.05`), so water stays resident. `MPUA` on the plate measures attached water directly. Single mesh after removing the GEOM half. **A 2:1 variant works:** split at x = 1.5 (sprinkler axis and plate centre) with x > 1.5 fine (`IJK=30,60,60`). The plate top then straddles the level jump. The plate is 1 coarse / 2 fine cells thick, so not thin, and on the fine grid | 1 mesh, 27,000 cells. 2:1 variant: 121,500 | ~2-5 min, 1 rank (T_END 50 s, dx 0.1). 2:1 variant: ~5-10 min with the fine part split | Committed row 289 (`m OBST` vs `geom_sprk_mass.csv`, Relative Error, max, tol 0.02); provisional (R-7) | `geom_sprk_mass.fds` lines 3, 10, 17, 20, 26 kept; 4, 18, 21, 27 deleted |
| alt. (committed, no derivation) | `Sprinklers_and_Sprays/bucket_test_1.fds` (IN) | Only committed multi-mesh wall-deposition case. Sprinkler at the point where the four MULT meshes meet. Droplets land and accumulate on the floor (a rectilinear SOLID wall) in all four meshes. Deposition only: no running over or dripping off edges | 4 meshes, 62,500 | ~2-5 min on 4 ranks | Row 58 (`Mass` = floor AMPUA vs `bucket_test_1.csv`, Relative Error, end, tol 0.02) | lines 3-4 (MESH + MULT 2 × 2), 9-10 (`PART SPEC_ID`), 19 (sprinkler XYZ 0,0,4.9), 32 (AMPUA floor integral) |

`Sprinklers_and_Sprays/cascade.fds` (IN, 1 mesh, 14,336 cells, row 73, Relative Error end_1_1 0.02; lines 6, 12, 14-16, 24-25, 38) is the single-mesh form of the same test. It is a good G4 member, but `VV/cascadempi_obst` supersedes it as an anchor.

### FR-050 (particle transport; cascadempi role)

| Rank | Case | Why | Meshes / cells | Runtime | Metric | Evidence |
|---|---|---|---|---|---|---|
| **1** | `Sprinklers_and_Sprays/bucket_test_1.fds` (IN, committed) | Spray starts at the 4-mesh corner, so nearly every droplet crosses one or two mesh faces in flight. About 10⁴-10⁵ droplets, which suits the integer ledger. The floor-mass metric catches lost or duplicated mass. **2:1 variant:** make one or two quadrants fine (`IJK=50,50,50`) so droplets cross level jumps | 4, 62,500 | ~2-5 min, 4 ranks | Row 58 (Relative Error end 0.02) | as above |
| **2** | `VV/cascadempi_obst.fds` (derived, as for FR-052) | Crossing at x = 7 both in flight and while attached, and through the 2:1 jump in the variant. It is the FR-052 anchor too, so one baseline serves both | 3, 28,672 | ~1 min | Rows 74-75 | as above |
| alt. | `Miscellaneous/part_path_ramp_jog.fds` (IN, committed) | One particle on a prescribed `PATH_RAMP` crosses the x = 0.5 mesh face between t ≈ 7 and 10 s. Deterministic, so it suits the per-ID attribute-integrity check. 2:1 variant: mesh 2 at `IJK=22,42,42` | 2, 9,702 | < 1 min | Row 567 (particle X/Z path, Relative Error end 0.001) | lines 3-4 (two meshes split at x = 0.5), 15-25 (ramps), 29 (`INIT … PATH_RAMP`, `N_PARTICLES=1`), 31-32 (`PARTICLE X/Z` DEVCs) |

### FR-051 (two-way coupling; cascadempi role)

`Energy_Budget/energy_budget_particles.fds` stays: it is named in FR-051 and IN (1 mesh, 1,000 cells, ~21 s; rows 209-210, area, 0.01; lines 3, 21, 28, 34-36). It is single-mesh, so a 2:1 split of the sealed box (x > 0.5 fine) is the refinement form. Replacement for the multi-mesh `cascadempi` role:

| Rank | Case | Why | Meshes / cells | Runtime | Metric | Evidence |
|---|---|---|---|---|---|---|
| **1** | `VV/cascadempi_obst.fds` (derived) | `HUMIDITY=0`, so all water vapour comes from droplet evaporation. The committed total (floor liquid + vapour volume integral) is a mass-exchange balance across 2 meshes, and across the 2:1 jump in the variant | 3, 28,672 | ~1 min | Rows 74-75 | `cascadempi.fds` line 16 (`HUMIDITY=0.`), 31-32 (`SPEC WATER VAPOR`, `PART SPEC_ID`), 48, 53 (vapour DEVCs) |
| **2** | `Miscellaneous/part_path_ramp_jog.fds` (IN) | The particle is a helium mass source (`SURF 'mass source'`, `MASS_FLUX=0.01`). As it crosses the mesh face (in the 2:1 variant, the level jump) the release must go to the finest level covering it. `VOL HE` gives the released mass. The committed metric covers only the path; the helium total is T2 vs the FDS baseline | 2, 9,702 | < 1 min | Row 567 + `VOL HE` T2 | lines 11, 27, 29, 33 |

**Recommendation (A-47):**
- Adopt `VV/cascadempi_obst` as the FR-052 anchor, with its 2:1 variant for the static-patch and regrid tests.
- Keep `geom_sprk_mass_obst` as secondary (R-7).
- Use `bucket_test_1` for FR-050.
- Use `energy_budget_particles` plus `cascadempi_obst` for FR-051.
- Each derived copy needs its own FDS baseline (A-35 rule). The committed `cascadempi` baseline covers only the full 6-mesh input.
- The derived inputs are not written yet. They go under `vv-runs/inputs/A-47/` with a README once the Spec Lead accepts the choice.
