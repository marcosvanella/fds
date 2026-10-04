# WP3 generator coverage (in the order of the generator answers) and the WP9b split

Sources read: `amrex/stage1-gpu-spike-plan.md` section 4 (rows WP3 and WP9b), `amrex/stage1-generator-answers.md` Q1 (the table and "Proposed order", items 1 to 7), `amrex/flux-hooks-kernel-review.md` (Finding 3 and "Rulings and notes added after the plan review"), `amrex/blocked-loop-families.md` (DP1), `amrex/loop_work_list.csv` (status, share, blocker per loop), `next-work.md` (Legacy Mapper item 2). Loop ids and shares are those of the loop survey (FireX 36975d7, modelled share of total run time, percent). Status is the work-list status at the generator branch `s5-gen`; the column "after s5-four" differs only where stated (the four wall nests, L0880 and L0882).

## 1. Interpretation used
- **WP3** is the plan's work package "Generator coverage for a periodic step: 39 loops, 7.15 % of modelled time, in the generator's order" (spike plan, section 4 row 3). "The order of the generator answers" is the "Proposed order" list under Q1 of `stage1-generator-answers.md` (items 1 to 7). The 39 loops are exactly the rows of the Q1 table (count checked: 39 ids, 7.153 %).
- **WP9b** is the plan's work package "Split the divergence chain at the DIF hook" (spike plan section 4 row 9b; review Finding 3). It is not a split of a package between owners. The split asked for is therefore read in two ways, both delivered in section 3: (a) the technical split of `DIVERGENCE_PART_1` into part A, the hook, part B; (b) a division of the package into pieces assignable to roles.
- Open question: the next-work item says "WP9b (split the divergence chain at the DIF hook)". If the intent was only to divide WP9b among owners without the technical map, section 3.2 alone is the answer.

## 2. WP3 coverage in generator-answers order
Totals for the 39 Q1 loops (7.153 %): translated at `s5-gen` 29 loops, 1.768 %; at `s5-four` 31 loops, 6.902 % (adds L0880 4.832 and L0882 0.302; the same four nests also make L0403 translated, outside the Q1 table); dead code (test-only) 4 loops, 0.196 %; open 4 loops, 0.055 % (L1391, L1375, L1392, L1393).

Status words: translated = kernel committed (bitwise-tested); translatable-now = front end accepts, nothing else needed; needs-feature = generator feature missing; blocked = waits for a table or decision.

### Item 1. Front end ok: markers and bitwise tests (answers: 0.27 % in "17 loops", 3 work-days)
The answers list 15 ids (L1374, L1394-L1396, L1369-L1371, L0395, L0396, L0881, L1391, L1343-L1346), not 17; their shares add to 0.274 %, so the share is right and the count "17" is a miscount in the answers document (the plan row also says 17).
| Loop | Lines | Share | Class now | Kernel / what is needed |
|---|---|---|---|---|
| L1374 | velo.f90:619-639 | 0.079 | translated | `vflux_vort_tau` |
| L1394, L1395, L1396 | velo.f90:1603-1609, 1613-1619, 1623-1629 | 0.005 each | translated | `vpred_us`, `vpred_vs`, `vpred_ws` |
| L1369, L1370, L1371 | velo.f90:1724-1730, 1734-1740, 1744-1750 | 0.005 each | translated | `vcorr_u`, `vcorr_v`, `vcorr_w` |
| L0395, L0396 | divg.f90:1611-1618, 1621-1630 | 0.012 each | translated | `div2_pred`, `div2_corr` |
| L0881 | mass.f90:201-209 | 0.071 | translated | `rho_rmw_mass` (Species and Combustion Lead; callee `GET_MOLECULAR_WEIGHT`) |
| L1343, L1344, L1345, L1346 | velo.f90:3254-3261, 3267-3275, 3282-3290, 3297-3305 | 0.011 each | translated | `baro_p_rrho`, `baro_fvx`, `baro_fvy`, `baro_fvz` |
| **L1391** | velo.f90:1243-1253 | 0.026 | **translatable-now, no owner** | no generator feature: add the marker and a bitwise test (`VELOCITY_FLUX_CYLINDRICAL`, cylindrical branch) |

### Item 2. `GET_SCALAR_FACE_VALUE` POINTER dummies (answers: 0.275 %, 3 to 4 days)
| Loop | Lines | Share | Class now | Kernel |
|---|---|---|---|---|
| L0652 to L0656 | func.f90:1348-1354, 1358-1368, 1372-1389, 1393-1410, 1414-1435 | 0.005, 0.016, 0.053, 0.053, 0.074 | translated | `gsfv_central`, `gsfv_godunov`, `gsfv_superbee`, `gsfv_minmod`, `gsfv_charm` (whole-field kernels, call-site dispatch) |
| L0657 | func.f90:1439-1453 | 0.074 | translated | `gsfv_mp5` |
Feature that was needed and is done: explicit-shape arguments for the POINTER dummies, element loop for the MP5 whole-array line. Remaining caveat: the MP5 limiter reads a fourth wall element (upstream patches UP-0001 to UP-0003 are tabled by the owner).

### Item 3. `MASS_FINITE_DIFFERENCES` (answers: L0880 4 to 6 days, L0882 1 to 2 days)
| Loop | Lines | Share | Class now (`s5-gen` / after `s5-four`) | Kernels and what is needed |
|---|---|---|---|---|
| L0880 | mass.f90:65-192 | 4.832 | in progress / translated | cell nest `rho_z_p_mass`, the `gsfv_*` call-site dispatch, wall nest `mass_wall_zz`; the pointer aims between them are host code. Feature: wall tables with the species index (done in the four nests) |
| L0882 | mass.f90:224-320 | 0.302 | in progress / translated | wall nest `mass_wall_rmw` (markers lines 224-320 = the whole loop); `spec_wall_zz` and `spec_wall_rmw` are the divg.f90:1021-1084 and 1120-1183 nests of `SPECIES_ADVECTION_PART_1_NEW` (L0401 inner, L0403), not part of L0882 |
| L0883 | mass.f90:328-350 | 0.128 | translated | `flux_mw_fix` |

### Item 4. `CHECK_STABILITY` and `CHECK_DIVERGENCE` (answers: done)
| Loop | Lines | Share | Class now | Kernel |
|---|---|---|---|---|
| L1347, L1348, L1349 | velo.f90:3059-3080, 3093-3108, 3121-3134 | 0.074, 0.014, 0.037 | translated | `cfl_max`, `cfl_wall_max` (device `**ONTH` last-bit difference is an open ruling), `vn_max` |
| L0363 | divg.f90:1675-1711 | 0.272 | translated | `div_extrema` |

### Item 5. Edge-table loops (answers: 0.706 %, design note 2 days plus 4 to 5 days)
| Loop | Lines | Share | Class now | Feature needed |
|---|---|---|---|---|
| L1376, L1377, L1378 | velo.f90:662-714, 720-772, 778-830 | 0.226 each | translated | edge flat tables built (`vflux_fvx`, `vflux_fvy`, `vflux_fvz`) |
| L1392, L1393 | velo.f90:1270-1302, 1306-1337 | 0.014 each | claimed (Generator Engineer), needs-feature | the front end must accept the K,I nest with J fixed (two-loop nest), then the same edge tables (`blocked-loop-families.md`, item 1 of the list before DP1) |
| L1375 | velo.f90:649-653 | 0.001 | claimed (Generator Engineer), needs-feature | callee `EVALUATE_RAMP` (reads derived-type components): callee flatten; the answers list this loop in item "VELOCITY_FLUX 4 of 5 not translatable" |

### Item 6. `DIVERGENCE_PART_2` L0394 (answers: 1 to 2 days, needs the alias-to-table decision)
| Loop | Lines | Share | Class now | Kernel |
|---|---|---|---|---|
| L0394 | divg.f90:1574-1604 | 0.042 | translated | `wall_bc_dp` (Wall Loops Engineer; the `B1` alias is the `B1_INDEX` wall table, patch UP-0004 proposes the matching source change) |

### Item 7. Dead code (answers: leave on the host)
| Loop | Lines | Share | Class now |
|---|---|---|---|
| L1397, L1398 | velo.f90:1648-1656, 1657-1665 | 0.049 each | not planned (test-only; guarded by `PERIODIC_TEST==7 .AND. .FALSE.`, velo.f90:1647) |
| L1372, L1373 | velo.f90:1770-1778, 1779-1787 | 0.049 each | not planned (test-only; guard at velo.f90:1769) |

### Beyond the Q1 table: the `DIVERGENCE_PART_1` nests (the plan's "diffusion chain", WP9b territory)
Owner of all of them: GPU Generator Engineer (`loop-work-list.md` section 4). Loop shares: L0365 7.679, L0369 2.399, L0381 1.491; the small ones add about 0.3.
| Loop | Lines | Share | Class now | Translated sub-nests | What blocks the rest (DP1 table, `blocked-loop-families.md`) |
|---|---|---|---|---|---|
| L0365 | divg.f90:128-235 | 7.679 | in progress | `rho_d_dzd`, `wall_rho_d_dzdn`, `d_z_max`, `rho_d_interp` | `RHO_D = MAX(0,MU)*RSC_T` whole-array forms (lines 122, 124, 146, 148: loop classifier, mine); `IF (CHECK_VN) D_Z_MAX = 0` and `DEL_RHO_D_DEL_Z = 0` (111, 117); `TENSOR_DIFFUSIVITY_MODEL` (187, stays host) |
| L0366 | divg.f90:245-258 | 0.074 | translated | `rho_d_maxloc_fix` | NaN species fractions are not guaranteed to follow the serial `MAXLOC` rule |
| L0369 | divg.f90:287-421 | 2.399 | in progress | `h_rho_d_dzd`, `dp_div_heat`, `del_rho_d_del_z` | `WALL_LOOP_2` (318-386): array sections without lo:hi in the emitter, wall species tables, `MAXLOC(B1%ZZ_F)`, `STORE_SPECIES_FLUX`, `GET_SENSIBLE_ENTHALPY_Z` call; EXIM and output copies (265-283: 4-D whole-array slices) |
| L0381 | divg.f90:640-660 | 1.491 | claimed | `dp_species` (inner) | N loop with `SPECIES_ADVECTION_PART_2` and `IF (CC_IBM) CALL`: host orchestration |
| L0364, L0371, L0373, L0377, L0378, L0379, L0380, L0382, L0383, L0386 | divg.f90:93-95, 463-471, 499-505, 573-582, 591-597, 608-616, 622-629, 668-674, 681-687, 757-767 | 0.099, 0.068, 0.012, 0.025, 0.012, 0.037, 0.025, 0.012, 0.012, 0.022 | claimed (blanket claim on the routine) | none | per loop: L0364 callee flatten (T3); L0379, L0383 generator gap; L0386 table needed; the others have no blocker (front end accepts) |
Line numbers in this subsection are those of the survey (36975d7). In the FDS-AMReX tree the same code sits at different lines (for example the species-sum loop at divg.f90:246-267 and the hook at divg.f90:279; see section 3.1).

## 3. WP9b: split of the divergence chain at the DIF hook

### 3.1 Technical split (FDS-AMReX tree, `Source/divg.f90`, routine `DIVERGENCE_PART_1`)
The hook call is `CALL FDS_HOOK_DIF_FLUX(...)` at divg.f90:279 (under `#ifdef WITH_AMREX`, divg.f90:276-280), between the species-sum fix and the diffusive heat flux. Order of the statements around it:
| Segment | FDS-AMReX lines | Survey id / kernel | Kernel today |
|---|---|---|---|
| **Part A** (before the hook): species fluxes `RHO_D_DZD*`, wall fluxes | divg.f90:137-244 (`DIFFUSIVE_FLUX_LOOP`, includes `TENSOR_DIFFUSIVITY_MODEL` call at 196) | L0365: `rho_d_interp`, `rho_d_dzd`, `d_z_max`, `wall_rho_d_dzdn` | yes, except the host pieces listed in section 2 |
| Part A: species-sum fix | divg.f90:246-270 (`MAXLOC` over species, loop 254-267) | L0366: `rho_d_maxloc_fix` | **yes** (committed; the review of the flux hooks says no kernel, which predates it) |
| Part A: EXIM store | divg.f90:274 (`IF (CC_IBM) CALL SET_EXIMDIFFLX_3D`) | none | out of scope while cut cells are off; stays before the read-out |
| **Hook** | divg.f90:276-280 | none | host control state plus the device gather/scatter of WP9 |
| Part B: flux output copies | divg.f90:282-298 (`IF (STORE_SPECIES_FLUX)`, 4-D whole-array copies) | none | no (host, or a data-movement kernel); reads `RHO_D_DZD*` after the hook |
| Part B: diffusive heat flux and `DEL_RHO_D_DEL_Z` | divg.f90:300-436 (`SPECIES_LOOP`); `WALL_LOOP_2` 333-401; divergence loop 423-431 (the review cites 426-430) | L0369: `h_rho_d_dzd`, `dp_div_heat`, `del_rho_d_del_z` | yes for the cell nests; `WALL_LOOP_2` is host |
| Part B: later | divg.f90:440 onward (specific heat, conduction, `dp_*`) | L0370 to L0386 | partly (`cp_rhg`, `kdtd`, `dp_kdtd`, `dp_species`) |

Consequence for the effort: the two kernels the plan names for WP9b (species-sum fix, `DEL_RHO_D_DEL_Z` divergence) already exist as `rho_d_maxloc_fix` (markers lines 245-258) and `del_rho_d_del_z` (markers lines 408-416) in `s5_markers.toml` at `s5-gen`, with bitwise tests and mutants (`test/make_div1_tests.py:23`, `test/mutation_div1.sh:26-34` including `maxloc_x_last_wins` and `maxloc_z_last_wins`, the last-maximum negative controls). What is left of WP9b is the device launch sequence (A, hook, B), the host pieces still in the chain, the sign-offs and the tests that need a two-level case.

### 3.2 Pieces assignable to roles
Defaults taken from the plan: the new kernels are a **derived copy** of the loops (as `fds_density_split.f90`), not a patch of `divg.f90`; zone-sum updates follow D-053 (per-cell terms on the device, serial addition in FDS order, no atomics).
| Piece | Role | Loop ids / lines | Acceptance test |
|---|---|---|---|
| 9b-1 Chain map and derived-copy statement: assign every statement of `DIVERGENCE_PART_1` from divg.f90:111 to 436 to part A, hook, part B or host-only (section 3.1 is the first draft, to be confirmed against the driver's `run_divergence_part1`) | Legacy Mapper | all of L0365, L0366, L0369 | the Integration Lead confirms each row against `TimeLoop.cpp` (`flux_apply_dif`, `run_divergence_part1`); no statement unassigned |
| 9b-2 Species-sum fix kernel sign-off (D-051 note) and tie test: the kernel re-implements the `MAXLOC` first-maximum rule | Species and Combustion Lead (sign-off), Generator Engineer (kernel and test) | L0366 `rho_d_maxloc_fix` (divg.f90:246-267) | sign-off note filed; test with at least 3 species and exact ties; the mutants `maxloc_x_last_wins` and `maxloc_z_last_wins` must fail it (I have not read whether the current test inputs contain exact ties; check first) |
| 9b-3 Close the NaN limit of `rho_d_maxloc_fix` or document it as accepted | Generator Engineer | L0366 | NaN input case: kernel equals the serial `MAXLOC`, or the limit is recorded in the sign-off |
| 9b-4 Divergence-loop kernel `del_rho_d_del_z` in the split chain, with `DP` re-initialised at the start of part A (no re-run precondition) | Generator Engineer (kernel), Integration Lead (chain sequence) | L0369 sub-nest divg.f90:423-431; `DP` zeroing at the start of the routine | device chain A-hook-B equals the host chain bitwise on `dec1` and `dec4_np4` (tests D-T1, D-T2) |
| 9b-5 Device launch sequence: part A kernels, then the hook (read-out pack or list upload plus scatter, WP9), then part B; remove the double run of `flux_apply_dif` | Integration Lead (kernels), Data Layout implementer (hook module) | L0365, L0366, L0369 | one `DIVERGENCE_PART_1` per stage per level in the device chain; D-T1 (empty override) device = host; D-T2 (no-op override, list = read-out values) through the split chain |
| 9b-6 Zone sums in the chain (`DSUM`, `PSUM`, `USUM` updates) under D-053 | Pressure Solver Lead (rule), Generator Engineer (per-cell terms) | divg.f90 zone-sum lines of `DIVERGENCE_PART_1`; family P1 | per-cell terms bitwise equal; additions serial in FDS order; no atomics (grep of the generated kernels) |
| 9b-7 Host pieces that stay in part B: `STORE_SPECIES_FLUX` copies (divg.f90:282-298) and `WALL_LOOP_2` (333-401): either data-movement kernels or the host path with a stated hand-off cost | Integration Lead | L0369 remainder | decision recorded; if host, the hand-off bytes per stage are measured against the chain time |
| 9b-8 Negative control and timing: scaled override must differ (D-T3); gather read-out equals the same faces of the full-box read-out (D-T4); read-out, upload, scatter and the two syncs at 64^3 and 128^3 (D-T5) | V&V Lead (controls), Integration Lead (timing) | all of the above | D-T3, D-T4, D-T5 of the flux-hooks review, on a two-level case with 1 and 4 boxes |
| 9b-9 Wording: amend D-061 so that the override point is inside part 1 of the chain | Chief Architect | D-061 text | decision log entry |
| 9b-10 Patch 0009 line shift: the hook patch shifts `divg.f90` lines after 267 by 9; any marker or doc that quotes FDS-AMReX lines after 267 must say which tree it uses | Legacy Mapper | divg.f90 references in `flux-hooks-kernel-review.md` and `s5_markers.toml` | a sidecar line check (`port_merge_check`) passes on both trees |

Cannot be tested by the periodic single-box case (no coarse-fine interface, no override): 9b-4, 9b-5 and 9b-8 need a two-level case.

### 3.3 Open questions
1. Is the plan's "3 to 4 work-days" for WP9b still right now that both kernels exist? By the table above the remaining work is 9b-2, 9b-4, 9b-5 and 9b-8; the kernel writing is done.
2. Which role owns the host pieces of part B (9b-7)? I assumed the Integration Lead.
3. The sign-off for the species-sum kernel (9b-2) is not in the register; I found no D-051 note naming `rho_d_maxloc_fix` (search of `docs` for `maxloc` finds only the plan, the review, the families file and the port-status tables).
