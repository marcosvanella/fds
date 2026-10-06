# Loop hand-out for the next round

Source: `amrex/loop_work_list.csv` plus the claim changes recorded with `upstream-1fb5c300-loop-impact.md` (pressure claims, species landings and ruling, radiation updates, the L0630 hold, MP5). Line ranges are at the current FDS-AMReX HEAD `8ad685f435` (FireX 1fb5c300fa merged), mapped from the survey ranges at `36975d7` with `git diff -U0` and checked: the first line of every range holds a `DO` statement at `8ad685f435`. Share is the modelled share of total time from the survey (a model, not a measurement).

## 1. The pool

Open loops are those that are unclaimed, translation-eligible (not host-side or test-only), not on hold and not MP5: **76 loops, 2.186 % of modelled time in total.** Pool rank (`#`) is by modelled share inside this pool. The pool is almost all velocity, turbulence and pressure-side work, and it has **no open loop in the species, mass, combustion or radiation domains, and none in `wall.f90`.** For the Species & Combustion Lead and the Radiation Lead, section 3 therefore offers surveyed loops that are outside the work list and transfer candidates; both need a decision by the Legacy Mapper before they become claimable.

## 2. Hand-out from the open pool

### GPU Mesh Data Loops Engineer (cell loops over mesh data: velo, turb, divg)

| # | Key | file:lines at 8ad685f435 | Routine | Share % | What (read at 8ad685f435) | Blockers |
|---|---|---|---|---|---|---|
| 1 | L1305 | turb.f90:1802-1867 | SYNTHETIC_TURBULENCE | 0.572 | pass 2: per vent (`VENT_APPLY_LOOP`), add each eddy's contribution to the vent-section arrays `VT%U_EDDY/V_EDDY/W_EDDY` (callee `EDDY_VELOCITY_CONTRIB`, turb.f90:1816), scale by `EDDY_C1` and subtract the mean of the normal component (turb.f90:1853-1865) | table needed (vent derived type); outer vent and eddy loops, not a K,J,I nest; the sum over eddies `NE` has to keep its order |
| 3 | L1317 | turb.f90:1033-1049 | TEST_FILTER | 0.284 | midpoint-quadrature test filter: `HAT = sum ORIG(I+L,J+M,K+N)*K3DM(L,M,N)` over a 3x3x3 stencil | generator gap (rank/subscript contract for the inner `L,M,N` stencil loops and the `K3DM` table); same routine as the Legacy Mapper's claimed batch L1315/L1316, so coordinate or fold it in |
| 4 | L1360 | velo.f90:422-444 | COMPUTE_STRAIN_RATE | 0.089 | default-model strain rate `STRAIN_RATE(I,J,K)` from `UU/VV/WW` | generator gap: the front end rejects `DUDX`.. as possibly live after the loop (the PRIVATE-list gap of the Legacy Mapper's L1355) |
| 6 | L1309 | turb.f90:1329-1351 | TENSOR_DIFFUSIVITY_MODEL | 0.068 | x-face nest: `RHO_D_DZDX(:,:,:,N) +=` the tensor-diffusivity term, `D_Z_MAX = MAX(D_Z_MAX,KX)` | generator gap (live-out scalars, MAX update); the routine is called from DIVERGENCE_PART_1 (L0365, claimed by the GPU Generator Engineer, call not translated), so agree the split with that owner |
| 7 | L1310 | turb.f90:1353-1375 | TENSOR_DIFFUSIVITY_MODEL | 0.068 | y-face nest, as L1309 (`RHO_D_DZDY`, `KY`) | as L1309 |
| 8 | L1311 | turb.f90:1377-1399 | TENSOR_DIFFUSIVITY_MODEL | 0.068 | z-face nest, as L1309 (`RHO_D_DZDZ`, `KZ`) | as L1309 |

### AMR Pressure Solver Lead (pressure and the pressure side of velo)

The new pressure claims (L1206, L1208, L1210, L1212-L1214, L1216-L1219, L1223-L1227) are all near 0 % of modelled time. The only open pressure-domain loop with real share is L1365. The remaining open `pres.f90` loops are below 0.01 % and blocked (section 4), so the list continues with velo.f90 force terms that feed the Poisson right-hand side through FVX/FVY/FVZ. There are no open divg.f90 loops: DIVERGENCE_PART_1 is under a routine claim.

| # | Key | file:lines at 8ad685f435 | Routine | Share % | What (read at 8ad685f435) | Blockers |
|---|---|---|---|---|---|---|
| 2 | L1365 | velo.f90:1419-1474 | NO_FLUX | 0.347 | per obstruction, three face nests drive FVX/FVY/FVZ toward zero inside solids using the HP gradient (velo.f90:1423-1472) | generator gap: outer loop over obstructions is not a K,J,I nest (host launch per box as in pvf_*); faces shared by overlapping obstructions get the same bits, so IDEMPOTENT fits (work-list entry) |
| 14 | L1391 | velo.f90:1255-1265 | VELOCITY_FLUX_CYLINDRICAL | 0.026 | cylindrical vorticity `OMY` and stress `TXZ` | none (front end accepts; no bitwise test yet); cylindrical meshes are deferred in AMR mode (charter Q4, D-042), so the value is mainly for uniform-mode parity with the L1210 claim |
| 29 | L1385 | velo.f90:908-915 | DIRECT_FORCE | 0.011 | `FVX -= RRHO*FVEC(1)*TIME_RAMP_FACTOR*SIN_THETA` | generator gap: constant-subscript rank-1 `FVEC(1)` as a by-value scalar; ramp factor from the host (velo.f90:899-905) |
| 30 | L1386 | velo.f90:929-936 | DIRECT_FORCE | 0.011 | `FVY` update with `FVEC(2)`, `COS_THETA` | as L1385 |
| 31 | L1387 | velo.f90:950-957 | DIRECT_FORCE | 0.011 | `FVZ` update with `FVEC(3)` | as L1385 |

### GPU Wall Loops Engineer (gas-side wall loops)

There is no open loop in `wall.f90`: L1486-L1489, L1485 and L1452 are claimed by the AMR Solid Phase Lead, and L1470/L1471 are closed. Of these, L1486 (`WALL_CELL_LOOP_0`, wall.f90:118-128, near-surface gas values) is purely gas-side; a split with the Solid Phase Lead is a possible option, but this file changes no claim. The open wall-index loops elsewhere are:

| # | Key | file:lines at 8ad685f435 | Routine | Share % | What (read at 8ad685f435) | Blockers |
|---|---|---|---|---|---|---|
| 23 | L0911 | part.f90:4643-4763 | PARTICLE_MOMENTUM_TRANSFER | 0.015 | external-wall loop over interpolated boundaries (corrector, NMESHES>1) with nests over the neighbour-mesh index box `EWC%IIO_MIN..KKO_MAX` | neighbour mesh (T4): `OM => OMESH(NOM)`, `M2 => MESHES(NOM)` (part.f90:4655-4656) |
| 32 | L0890 | part.f90:605-607 | INSERT_VENT_PARTICLES | 0.008 | wall loop calling `PARTICLE_FACE_INSERT(WALL_INDEX=IW)` | callee flatten (T3); particle creation, which may end up host-side |
| 43 | L1204 | pres.f90:1931-1994 | ULMAT_SOLVE_ZONE | 0.003 | external-wall loop setting the HP ghost value per wall for Neumann or Dirichlet pressure BCs, filtered by `ZONE_MESH(..)%CONNECTED_ZONE_PARENT` | inventory tags it neighbour mesh (T4), but the loop text reads only `EXTERNAL_WALL(IW)%PRESSURE_BC_TYPE`, `BOUNDARY_COORD`, a zone table and own-mesh arrays, so re-check the tag; per-wall ghost store like `wall_hs_bt` |

### AMR Solid Phase Lead (wall solid and pyrolysis, solid-phase set-up, level set)

| # | Key | file:lines at 8ad685f435 | Routine | Share % | What (read at 8ad685f435) | Blockers |
|---|---|---|---|---|---|---|
| 17 | L0817 | init.f90:4889-4900 | CREATE_OR_REMOVE_OBST | 0.017 | reset `RHO/RHOS` (and `TMP` when `SOLID_PHASE_ONLY`), `ZZ/ZZS` to ambient and zero `DEL_RHO_D_DEL_Z` when the obstruction touches the mesh boundary (init.f90:4897), in the cells of a created or removed obstruction | generator gap (derived-type section `SPECIES_MIXTURE(1:N)%ZZ0`); event-driven (called at init.f90:4709, 4715); runs right after BLOCK_CELL (init.f90:4885), so it is paired with the L0630 hold (impact report section 2) |
| 76 | L1340 | vege.f90:763-857 | LEVEL_SET_SPREAD_RATE | 0.000 | level-set spread-rate J,I nest (`FLUX_ILOOP`) | generator gap: 2-D J,I nest (front end takes K,J,I); level-set wildfire must work with refinement (charter Q4) |
| 73 | L1334 | vege.f90:885-894 | LEVEL_SET_ADVECT_FLUX | 0.000 | x-face level-set flux `F_X` via `SCALAR_FACE_VALUE_LS` on a 4-point stencil | 2-D nest; function callee; pointer `PHI_LS_P` |
| 74 | L1335 | vege.f90:896-905 | LEVEL_SET_ADVECT_FLUX | 0.000 | y-face flux `F_Y`, as L1334 | as L1334 |
| 75 | L1336 | vege.f90:907-914 | LEVEL_SET_ADVECT_FLUX | 0.000 | `FLUX_LS_P` from the face fluxes, `CYCLE` unless `SURFACE(LS_SURF_INDEX(I,J))%VEG_LSET_SPREAD` | 2-D nest; per-surface flag gather |

Share is about 0.017 % in total. The Solid Phase Lead's main work stays the claimed WALL_BC set (L1488 2.688 %, L1485, L1489, L1486, L1487).

## 3. Owners with no open loop on the work list

### AMR Species & Combustion Lead (species, combustion, mass, soot, chem/fire cell loops)

All species, mass, combustion and soot loops on the work list are translated or claimed. L1272, L0866, L0875 and L0397 are now translated (s5-gen `fe8017c109`, `0b1284e5d7`); L0876/L0877 are hand-written, owned by this lead with the Legacy Mapper as reviewer; L0596/L1273 are blocked on data layout. The candidates below are therefore surveyed loops (`inventory/gpu_generator_loop_classes.csv` keys) that the work list leaves out by class, plus two Legacy Mapper claims. Each needs a Legacy Mapper decision (reclassify, or release the claim) before it can be claimed.

| Key | file:lines at 8ad685f435 | Routine | Share % | Status | What (read at 8ad685f435) | Blockers |
|---|---|---|---|---|---|---|
| L1270 | soot.f90:280-389 | CALC_AGGLOMERATION | 0.336 | off-list: survey class host control (sequential control flow) | per species, per gas cell: Brownian, gravitational and shear kernel `PHI`, then a per-cell `DO WHILE` sub-stepped agglomeration (soot.f90:356-383) and the write-back of the bins to `ZZ` | per-cell `DO WHILE` (needs a loop-shape feature or a hand kernel); callees `GET_VISCOSITY`, `CUNNINGHAM`; private arrays sized by `MAXVAL(N_PARTICLE_BINS)`; rank-2/3 bin tables; `PBAR(K,PRESSURE_ZONE)` gather, which the generator now has (`0b1284e5d7`) |
| L0398 | divg.f90:850-954 | ENTHALPY_ADVECTION_NEW | 0.260 | claimed by Legacy Mapper (transfer candidate) | wall loop with accumulate into `U_DOT_DEL_RHO_H_S` (work-list entry) | pointer scratch (S2) and the CSR cell-to-wall gather (S3) |
| L0609 | fire.f90:1963-2224 | CONDENSATION_EVAPORATION | 0.087 | off-list: geometry-deferred | species loop over condensable species: gas-phase cell nest (fire.f90:1974-2080) and solid-surface wall loop (fire.f90:2084-2223) | geometry bucket only through the `CC_IBM` / `CCVAR` test (fire.f90:1978-1980); species-mixture and species tables; wall part writes `B1%AWM_AEROSOL` |
| L0597 | fire.f90:123-137 | COMBUSTION_GENERAL_LOAD_BALANCED | 0.002 | off-list: survey class host control | marks chemically active cells and appends them to `CHEM_ACTIVE_CELLS` (counter, fire.f90:131-134) | stream compaction in source order; callee `CHECK_CHEMICALLY_ACTIVE_STATE`; `RETURN` on `STOP_STATUS` (fire.f90:129) |
| L0401 | divg.f90:1010-1103 | SPECIES_ADVECTION_PART_1_NEW | 3.260 | claimed by Legacy Mapper (transfer candidate) | cell nest (translated) and wall nest with pointer scratch | the remaining blocker was the MP5 fourth wall element (UP-0001/UP-0002); with the ruling "MP5 not implemented in the generator", what is left is a non-MP5 sign-off of the committed `rho_z_p_divg` and `spec_wall_zz` |

### AMR Radiation Lead (radi.f90, init.f90, dump.f90; radiation-adjacent and wall/table-heavy)

The work list has four radiation loops, and none is open: L1239 is translated, L1243 is translated, L1245 stays host-side and L1248 is in progress. L1242 and the stacked UII/QR/QR_W kernels are excluded as asked. The survey holds only 11 `radi.f90` loops, and the radiation set-up loops are not in it, so they have no keys. `dump.f90` has no surveyed loop. The open `init.f90` loop (L0817) is listed under the Solid Phase Lead. The candidates are therefore surveyed loops outside the work list. They have survey keys but no work-list rank, and each needs a Legacy Mapper reclassification (particle loops are host-side by class, `tools/inventory/loop_work_list.py:26`; geometry rows are outside every denominator, `:376`).

| Key | file:lines at 8ad685f435 | Routine | Share % | Status | What (read at 8ad685f435) | Blockers |
|---|---|---|---|---|---|---|
| L1238 | radi.f90:5105-5144 | ADD_VOLUMETRIC_HEAT_SOURCE | 0.024 | off-list: geometry-deferred (host-control class) | per `INIT` with `HRRPUV>0`: ramp factor (`EVALUATE_RAMP`, radi.f90:5109), then a K,J,I nest over 0:KBP1 adding to `Q` and, in mode 1, to `KFST4_GAS` (radi.f90:5122-5129) | outer `N_INIT` loop and `SELECT CASE(IN%SHAPE)` (host launch per INIT, like pvf_*); the cut-cell block (radi.f90:5132-5140) is the geometry part and would stay out |
| L1244 | radi.f90:4993-5038 | RADIATION_FVM | 0.001 | off-list: particle (host by class) | incoming radiative flux on solid particles: masked `SUM` of `BR%BAND(IBND)%ILW` over the view cone, external flux ramp | particle storage layout; ragged `BR_ILW` table (as L1243/L1248); callees `EVALUATE_RAMP`, `CALC_PARTICLE_TEMP_ORIENTATION` |
| L1240 | radi.f90:3894-3901 | RADIATION_FVM | 0.000 | off-list: particle (host by class) | zero `B1%Q_RAD_IN` of particles where `TMP_GAS_FRONT<=0`, the particle twin of L1239 (`rad_wall_qin_zero`) | particle storage layout; particle `B1` table |
| L1250 | radi.f90:3679-3695 | INTERPOLATE_IL | 0.000 | off-list: particle (host by class) | particle twin of L1248: interpolate `ILW` between angles for oriented solid particles | as L1248 (ragged `BR_ILW`, keep the `ILW_OLD` copy); do it with L1248's host version |
| L1241 | radi.f90:3902-3908 | RADIATION_FVM | 0.000 | off-list: geometry-deferred | CFACE twin of L1239 | GEOM/CC_IBM deferred (charter Q4, D-033) |
| L1249 | radi.f90:3665-3676 | INTERPOLATE_IL | 0.000 | off-list: geometry-deferred | CFACE twin of L1248 | GEOM/CC_IBM deferred (charter Q4, D-033) |

## 4. Open pool loops not handed out this round

| # | Key | file:lines at 8ad685f435 | Routine | Share % | Class | Blocker |
|---|---|---|---|---|---|---|
| 5 | L1121 | pres.f90:4012-4026 | CHECK_UNSUPPORTED_MESH | 0.086 | blocked | blocked |
| 9 | L1312 | turb.f90:1420-1440 | TENSOR_DIFFUSIVITY_MODEL | 0.063 | needs-feature | generator gap |
| 10 | L1313 | turb.f90:1442-1462 | TENSOR_DIFFUSIVITY_MODEL | 0.063 | needs-feature | generator gap |
| 11 | L1314 | turb.f90:1464-1485 | TENSOR_DIFFUSIVITY_MODEL | 0.063 | translatable-now | none |
| 12 | L1275 | turb.f90:825-837 | CALC_VARDEN_LEONARD_TERM | 0.037 | needs-feature | generator gap |
| 13 | L1323 | turb.f90:639-651 | VARDEN_DYNSMAG | 0.037 | needs-feature | generator gap |
| 15 | L1325 | turb.f90:696-716 | VARDEN_DYNSMAG | 0.021 | needs-feature | generator gap |
| 16 | L1354 | velo.f90:196-205 | COMPUTE_VISCOSITY | 0.021 | translatable-now | none |
| 18 | L0677 | func.f90:5248-5266 | PACK_CELL | 0.016 | blocked | neighbour mesh (T4) |
| 19 | L1379 | velo.f90:982-990 | CORIOLIS_FORCE | 0.016 | translatable-now | none |
| 20 | L1381 | velo.f90:1004-1012 | CORIOLIS_FORCE | 0.016 | needs-feature | generator gap |
| 21 | L1382 | velo.f90:1018-1026 | CORIOLIS_FORCE | 0.016 | needs-feature | generator gap |
| 22 | L1383 | velo.f90:1032-1040 | CORIOLIS_FORCE | 0.016 | needs-feature | generator gap |
| 24 | L1276 | turb.f90:193-200 | COMPRESSION_WAVE | 0.011 | translatable-now | none |
| 25 | L1277 | turb.f90:201-208 | COMPRESSION_WAVE | 0.011 | translatable-now | none |
| 26 | L1278 | turb.f90:209-216 | COMPRESSION_WAVE | 0.011 | translatable-now | none |
| 27 | L1279 | turb.f90:217-224 | COMPRESSION_WAVE | 0.011 | translatable-now | none |
| 28 | L1324 | turb.f90:672-684 | VARDEN_DYNSMAG | 0.011 | needs-feature | generator gap |
| 33 | L1280 | turb.f90:227-233 | COMPRESSION_WAVE | 0.005 | translatable-now | none |
| 34 | L1281 | turb.f90:234-240 | COMPRESSION_WAVE | 0.005 | translatable-now | none |
| 35 | L1282 | turb.f90:241-247 | COMPRESSION_WAVE | 0.005 | translatable-now | none |
| 36 | L1283 | turb.f90:248-254 | COMPRESSION_WAVE | 0.005 | translatable-now | none |
| 37 | L1284 | turb.f90:256-262 | COMPRESSION_WAVE | 0.005 | translatable-now | none |
| 38 | L1285 | turb.f90:263-269 | COMPRESSION_WAVE | 0.005 | translatable-now | none |
| 39 | L1286 | turb.f90:270-276 | COMPRESSION_WAVE | 0.005 | translatable-now | none |
| 40 | L1287 | turb.f90:277-283 | COMPRESSION_WAVE | 0.005 | translatable-now | none |
| 41 | L1322 | turb.f90:576-587 | VARDEN_DYNSMAG | 0.005 | needs-feature | generator gap |
| 42 | L1351 | velo.f90:124-130 | COMPUTE_VISCOSITY | 0.005 | translatable-now | none |
| 44 | L1198 | pres.f90:1837-1846 | ULMAT_SOLVE_ZONE | 0.002 | blocked | neighbour mesh (T4) |
| 45 | L0680 | func.f90:2845-2847 | COMPUTE_WIND_COMPONENTS | 0.001 | blocked | callee flatten (T3) |
| 46 | L1130 | pres.f90:5580-5601 | GET_H_REGFACES | 0.001 | needs-feature | table needed |
| 47 | L1154 | pres.f90:2031-2047 | ULMAT_GET_H_REGFACES | 0.001 | needs-feature | table needed |
| 48 | L1190 | pres.f90:1686-1694 | ULMAT_SOLVE_ZONE | 0.001 | blocked | neighbour mesh (T4) |
| 49 | L1192 | pres.f90:1713-1720 | ULMAT_SOLVE_ZONE | 0.001 | blocked | neighbour mesh (T4) |
| 50 | L1194 | pres.f90:1782-1790 | ULMAT_SOLVE_ZONE | 0.001 | blocked | neighbour mesh (T4) |
| 51 | L1196 | pres.f90:1808-1815 | ULMAT_SOLVE_ZONE | 0.001 | blocked | neighbour mesh (T4) |
| 52 | L1318 | turb.f90:1060-1065 | TEST_FILTER | 0.001 | needs-feature | generator gap |
| 53 | L1319 | turb.f90:1069-1074 | TEST_FILTER | 0.001 | needs-feature | generator gap |
| 54 | L1320 | turb.f90:1078-1083 | TEST_FILTER | 0.001 | needs-feature | generator gap |
| 55 | L0884 | part.f90:4606-4617 | CLIP_PARTICLE_DRAG | 0.000 | translatable-now | none |
| 56 | L1203 | pres.f90:1923-1925 | ULMAT_SOLVE_ZONE | 0.000 | needs-feature | generator gap |
| 57 | L1289 | turb.f90:888-890 | FILL_EDGES | 0.000 | needs-feature | generator gap |
| 58 | L1290 | turb.f90:893-895 | FILL_EDGES | 0.000 | needs-feature | generator gap |
| 59 | L1291 | turb.f90:898-900 | FILL_EDGES | 0.000 | needs-feature | generator gap |
| 60 | L1292 | turb.f90:903-905 | FILL_EDGES | 0.000 | needs-feature | generator gap |
| 61 | L1293 | turb.f90:910-912 | FILL_EDGES | 0.000 | needs-feature | generator gap |
| 62 | L1294 | turb.f90:915-917 | FILL_EDGES | 0.000 | needs-feature | generator gap |
| 63 | L1295 | turb.f90:920-922 | FILL_EDGES | 0.000 | needs-feature | generator gap |
| 64 | L1296 | turb.f90:925-927 | FILL_EDGES | 0.000 | needs-feature | generator gap |
| 65 | L1297 | turb.f90:932-934 | FILL_EDGES | 0.000 | needs-feature | generator gap |
| 66 | L1298 | turb.f90:937-939 | FILL_EDGES | 0.000 | needs-feature | generator gap |
| 67 | L1299 | turb.f90:942-944 | FILL_EDGES | 0.000 | needs-feature | generator gap |
| 68 | L1300 | turb.f90:947-949 | FILL_EDGES | 0.000 | needs-feature | generator gap |
| 69 | L1330 | vege.f90:618-623 | GET_BOUNDARY_VALUES | 0.000 | blocked | callee flatten (T3) |
| 70 | L1331 | vege.f90:624-629 | GET_BOUNDARY_VALUES | 0.000 | blocked | callee flatten (T3) |
| 71 | L1332 | vege.f90:630-639 | GET_BOUNDARY_VALUES | 0.000 | blocked | callee flatten (T3) |
| 72 | L1333 | vege.f90:678-699 | FILL_BOUNDARY_VALUES | 0.000 | blocked | neighbour mesh (T4) |

Handed out from the pool: 19 loops, 1.598 % of the pool's 2.186 %. The rest is mostly COMPRESSION_WAVE, FILL_EDGES, TEST_FILTER edge loops and ULMAT set-up loops at or below 0.04 % each.

## 5. Work list regeneration

The work list was not regenerated. `tools/inventory/loop_work_list.py` takes its ranges from the survey CSV at `36975d7` and attaches kernels by matching marker `lines` (also `36975d7`, or `bee11f0329` for radiation) against them. `--src-ref` only changes the `DO`-line check. Pointing the script at `8ad685f435` therefore needs the survey ranges and the s5-gen markers re-pinned together; the impact report section 5 lists the marker ranges, and section 6 below lists the loop ranges. `--check` was already failing before these edits, because s5-gen has new kernels the list does not show yet (L0397, L0866, L0875; L1272 as a four-kernel loop), and it fails now with the claim changes too. A dry run into a scratch copy, with the generator change described next, gives 78 translated loops (11.879 %) and one loop on hold.

`tools/inventory/loop_work_list.py` was edited so that the next regeneration shows the new claim statuses. The changes:
- `hold-upstream` and `closed-owner-ruling` are passed through as statuses with their own work classes ("on hold (upstream)", "closed (owner ruling)"), and those loops get no proposed owner.
- L1272 is a composite of `settle_zero`, `settle_vel`, `settle_wall_zero` and `settle_div`.

## 6. All work-list loops: survey range and range at 8ad685f435

| Key | Survey (36975d7) | At 8ad685f435 | Moved |
|---|---|---|---|
| L0363 | divg.f90:1675-1711 | divg.f90:1690-1726 | +15 |
| L0364 | divg.f90:93-95 | divg.f90:102-104 | +9 |
| L0365 | divg.f90:128-235 | divg.f90:137-244 | +9 |
| L0366 | divg.f90:245-258 | divg.f90:254-267 | +9 |
| L0369 | divg.f90:287-421 | divg.f90:302-436 | +15 |
| L0370 | divg.f90:435-443 | divg.f90:450-458 | +15 |
| L0371 | divg.f90:463-471 | divg.f90:478-486 | +15 |
| L0372 | divg.f90:482-486 | divg.f90:497-501 | +15 |
| L0373 | divg.f90:499-505 | divg.f90:514-520 | +15 |
| L0374 | divg.f90:512-523 | divg.f90:527-538 | +15 |
| L0375 | divg.f90:532-554 | divg.f90:547-569 | +15 |
| L0376 | divg.f90:561-570 | divg.f90:576-585 | +15 |
| L0377 | divg.f90:573-582 | divg.f90:588-597 | +15 |
| L0378 | divg.f90:591-597 | divg.f90:606-612 | +15 |
| L0379 | divg.f90:608-616 | divg.f90:623-631 | +15 |
| L0380 | divg.f90:622-629 | divg.f90:637-644 | +15 |
| L0381 | divg.f90:640-660 | divg.f90:655-675 | +15 |
| L0382 | divg.f90:668-674 | divg.f90:683-689 | +15 |
| L0383 | divg.f90:681-687 | divg.f90:696-702 | +15 |
| L0384 | divg.f90:696-716 | divg.f90:711-731 | +15 |
| L0386 | divg.f90:757-767 | divg.f90:772-782 | +15 |
| L0394 | divg.f90:1574-1604 | divg.f90:1589-1619 | +15 |
| L0395 | divg.f90:1611-1618 | divg.f90:1626-1633 | +15 |
| L0396 | divg.f90:1621-1630 | divg.f90:1636-1645 | +15 |
| L0397 | divg.f90:814-822 | divg.f90:829-837 | +15 |
| L0398 | divg.f90:835-939 | divg.f90:850-954 | +15 |
| L0399 | divg.f90:946-966 | divg.f90:961-981 | +15 |
| L0400 | divg.f90:1301-1319 | divg.f90:1316-1334 | +15 |
| L0401 | divg.f90:995-1088 | divg.f90:1010-1103 | +15 |
| L0402 | divg.f90:1097-1105 | divg.f90:1112-1120 | +15 |
| L0403 | divg.f90:1120-1183 | divg.f90:1135-1198 | +15 |
| L0404 | divg.f90:1191-1213 | divg.f90:1206-1228 | +15 |
| L0405 | divg.f90:1237-1266 | divg.f90:1252-1281 | +15 |
| L0406 | divg.f90:1269-1289 | divg.f90:1284-1304 | +15 |
| L0596 | fire.f90:1902-1918 | fire.f90:1902-1918 | no |
| L0608 | fire.f90:1957-1961 | fire.f90:1957-1961 | no |
| L0630 | func.f90:5509-5524 | func.f90:5581-5627 | +72 |
| L0652 | func.f90:1348-1354 | func.f90:1348-1354 | no |
| L0653 | func.f90:1358-1368 | func.f90:1358-1368 | no |
| L0654 | func.f90:1372-1389 | func.f90:1372-1389 | no |
| L0655 | func.f90:1393-1410 | func.f90:1393-1410 | no |
| L0656 | func.f90:1414-1435 | func.f90:1414-1435 | no |
| L0657 | func.f90:1439-1453 | func.f90:1439-1453 | no |
| L0677 | func.f90:5183-5201 | func.f90:5248-5266 | +65 |
| L0680 | func.f90:2845-2847 | func.f90:2845-2847 | no |
| L0817 | init.f90:4878-4889 | init.f90:4889-4900 | +11 |
| L0822 | init.f90:4917-4936 | init.f90:4928-4947 | +11 |
| L0823 | init.f90:4940-5054 | init.f90:4951-5065 | +11 |
| L0858 | mass.f90:424-436 | mass.f90:430-442 | +6 |
| L0859 | mass.f90:442-455 | mass.f90:448-461 | +6 |
| L0860 | mass.f90:462-472 | mass.f90:468-478 | +6 |
| L0861 | mass.f90:484-495 | mass.f90:490-501 | +6 |
| L0862 | mass.f90:505-512 | mass.f90:511-518 | +6 |
| L0863 | mass.f90:527-534 | mass.f90:533-540 | +6 |
| L0865 | mass.f90:556-564 | mass.f90:562-570 | +6 |
| L0866 | mass.f90:570-577 | mass.f90:576-583 | +6 |
| L0867 | mass.f90:595-608 | mass.f90:601-614 | +6 |
| L0868 | mass.f90:619-632 | mass.f90:625-638 | +6 |
| L0869 | mass.f90:650-660 | mass.f90:656-666 | +6 |
| L0870 | mass.f90:666-677 | mass.f90:672-683 | +6 |
| L0871 | mass.f90:687-694 | mass.f90:693-700 | +6 |
| L0872 | mass.f90:709-716 | mass.f90:715-722 | +6 |
| L0874 | mass.f90:738-746 | mass.f90:744-752 | +6 |
| L0875 | mass.f90:752-759 | mass.f90:758-765 | +6 |
| L0876 | mass.f90:799-849 | mass.f90:805-855 | +6 |
| L0877 | mass.f90:868-939 | mass.f90:874-945 | +6 |
| L0878 | mass.f90:947-961 | mass.f90:953-967 | +6 |
| L0879 | mass.f90:980-987 | mass.f90:986-993 | +6 |
| L0880 | mass.f90:65-192 | mass.f90:71-198 | +6 |
| L0881 | mass.f90:201-209 | mass.f90:207-215 | +6 |
| L0882 | mass.f90:224-320 | mass.f90:230-326 | +6 |
| L0883 | mass.f90:328-350 | mass.f90:334-356 | +6 |
| L0884 | part.f90:4606-4617 | part.f90:4606-4617 | no |
| L0890 | part.f90:605-607 | part.f90:605-607 | no |
| L0906 | part.f90:3654-3659 | part.f90:3654-3659 | no |
| L0911 | part.f90:4643-4763 | part.f90:4643-4763 | no |
| L1121 | pres.f90:4012-4026 | pres.f90:4012-4026 | no |
| L1122 | pres.f90:4351-4356 | pres.f90:4351-4356 | no |
| L1130 | pres.f90:5580-5601 | pres.f90:5580-5601 | no |
| L1154 | pres.f90:2031-2047 | pres.f90:2031-2047 | no |
| L1190 | pres.f90:1686-1694 | pres.f90:1686-1694 | no |
| L1192 | pres.f90:1713-1720 | pres.f90:1713-1720 | no |
| L1194 | pres.f90:1782-1790 | pres.f90:1782-1790 | no |
| L1196 | pres.f90:1808-1815 | pres.f90:1808-1815 | no |
| L1198 | pres.f90:1837-1846 | pres.f90:1837-1846 | no |
| L1203 | pres.f90:1923-1925 | pres.f90:1923-1925 | no |
| L1204 | pres.f90:1931-1994 | pres.f90:1931-1994 | no |
| L1206 | pres.f90:729-742 | pres.f90:729-742 | no |
| L1207 | pres.f90:758-764 | pres.f90:758-764 | no |
| L1208 | pres.f90:768-787 | pres.f90:768-787 | no |
| L1209 | pres.f90:65-228 | pres.f90:65-228 | no |
| L1210 | pres.f90:238-245 | pres.f90:238-245 | no |
| L1211 | pres.f90:250-260 | pres.f90:250-260 | no |
| L1212 | pres.f90:267-277 | pres.f90:267-277 | no |
| L1213 | pres.f90:282-292 | pres.f90:282-292 | no |
| L1214 | pres.f90:297-307 | pres.f90:297-307 | no |
| L1215 | pres.f90:394-400 | pres.f90:394-400 | no |
| L1216 | pres.f90:404-410 | pres.f90:404-410 | no |
| L1217 | pres.f90:414-420 | pres.f90:414-420 | no |
| L1218 | pres.f90:424-430 | pres.f90:424-430 | no |
| L1219 | pres.f90:438-440 | pres.f90:438-440 | no |
| L1220 | pres.f90:450-462 | pres.f90:450-462 | no |
| L1221 | pres.f90:466-477 | pres.f90:466-477 | no |
| L1222 | pres.f90:481-492 | pres.f90:481-492 | no |
| L1223 | pres.f90:531-537 | pres.f90:531-537 | no |
| L1224 | pres.f90:544-567 | pres.f90:544-567 | no |
| L1225 | pres.f90:576-582 | pres.f90:576-582 | no |
| L1226 | pres.f90:583-589 | pres.f90:583-589 | no |
| L1227 | pres.f90:591-596 | pres.f90:591-596 | no |
| L1239 | radi.f90:3886-3892 | radi.f90:3887-3893 | +1 |
| L1243 | radi.f90:4961-4970 | radi.f90:4965-4974 | +4 |
| L1245 | radi.f90:5044-5059 | radi.f90:5048-5063 | +4 |
| L1248 | radi.f90:3651-3662 | radi.f90:3651-3662 | no |
| L1272 | soot.f90:54-175 | soot.f90:54-175 | no |
| L1273 | soot.f90:412-472 | soot.f90:412-472 | no |
| L1275 | turb.f90:819-831 | turb.f90:825-837 | +6 |
| L1276 | turb.f90:187-194 | turb.f90:193-200 | +6 |
| L1277 | turb.f90:195-202 | turb.f90:201-208 | +6 |
| L1278 | turb.f90:203-210 | turb.f90:209-216 | +6 |
| L1279 | turb.f90:211-218 | turb.f90:217-224 | +6 |
| L1280 | turb.f90:221-227 | turb.f90:227-233 | +6 |
| L1281 | turb.f90:228-234 | turb.f90:234-240 | +6 |
| L1282 | turb.f90:235-241 | turb.f90:241-247 | +6 |
| L1283 | turb.f90:242-248 | turb.f90:248-254 | +6 |
| L1284 | turb.f90:250-256 | turb.f90:256-262 | +6 |
| L1285 | turb.f90:257-263 | turb.f90:263-269 | +6 |
| L1286 | turb.f90:264-270 | turb.f90:270-276 | +6 |
| L1287 | turb.f90:271-277 | turb.f90:277-283 | +6 |
| L1289 | turb.f90:882-884 | turb.f90:888-890 | +6 |
| L1290 | turb.f90:887-889 | turb.f90:893-895 | +6 |
| L1291 | turb.f90:892-894 | turb.f90:898-900 | +6 |
| L1292 | turb.f90:897-899 | turb.f90:903-905 | +6 |
| L1293 | turb.f90:904-906 | turb.f90:910-912 | +6 |
| L1294 | turb.f90:909-911 | turb.f90:915-917 | +6 |
| L1295 | turb.f90:914-916 | turb.f90:920-922 | +6 |
| L1296 | turb.f90:919-921 | turb.f90:925-927 | +6 |
| L1297 | turb.f90:926-928 | turb.f90:932-934 | +6 |
| L1298 | turb.f90:931-933 | turb.f90:937-939 | +6 |
| L1299 | turb.f90:936-938 | turb.f90:942-944 | +6 |
| L1300 | turb.f90:941-943 | turb.f90:947-949 | +6 |
| L1305 | turb.f90:1796-1861 | turb.f90:1802-1867 | +6 |
| L1309 | turb.f90:1323-1345 | turb.f90:1329-1351 | +6 |
| L1310 | turb.f90:1347-1369 | turb.f90:1353-1375 | +6 |
| L1311 | turb.f90:1371-1393 | turb.f90:1377-1399 | +6 |
| L1312 | turb.f90:1414-1434 | turb.f90:1420-1440 | +6 |
| L1313 | turb.f90:1436-1456 | turb.f90:1442-1462 | +6 |
| L1314 | turb.f90:1458-1479 | turb.f90:1464-1485 | +6 |
| L1315 | turb.f90:979-995 | turb.f90:985-1001 | +6 |
| L1316 | turb.f90:1003-1019 | turb.f90:1009-1025 | +6 |
| L1317 | turb.f90:1027-1043 | turb.f90:1033-1049 | +6 |
| L1318 | turb.f90:1054-1059 | turb.f90:1060-1065 | +6 |
| L1319 | turb.f90:1063-1068 | turb.f90:1069-1074 | +6 |
| L1320 | turb.f90:1072-1077 | turb.f90:1078-1083 | +6 |
| L1321 | turb.f90:501-534 | turb.f90:507-540 | +6 |
| L1322 | turb.f90:570-581 | turb.f90:576-587 | +6 |
| L1323 | turb.f90:633-645 | turb.f90:639-651 | +6 |
| L1324 | turb.f90:666-678 | turb.f90:672-684 | +6 |
| L1325 | turb.f90:690-710 | turb.f90:696-716 | +6 |
| L1330 | vege.f90:618-623 | vege.f90:618-623 | no |
| L1331 | vege.f90:624-629 | vege.f90:624-629 | no |
| L1332 | vege.f90:630-639 | vege.f90:630-639 | no |
| L1333 | vege.f90:678-699 | vege.f90:678-699 | no |
| L1334 | vege.f90:885-894 | vege.f90:885-894 | no |
| L1335 | vege.f90:896-905 | vege.f90:896-905 | no |
| L1336 | vege.f90:907-914 | vege.f90:907-914 | no |
| L1340 | vege.f90:763-857 | vege.f90:763-857 | no |
| L1343 | velo.f90:3254-3261 | velo.f90:3288-3295 | +34 |
| L1344 | velo.f90:3267-3275 | velo.f90:3301-3309 | +34 |
| L1345 | velo.f90:3282-3290 | velo.f90:3316-3324 | +34 |
| L1346 | velo.f90:3297-3305 | velo.f90:3331-3339 | +34 |
| L1347 | velo.f90:3059-3080 | velo.f90:3093-3114 | +34 |
| L1348 | velo.f90:3093-3108 | velo.f90:3127-3142 | +34 |
| L1349 | velo.f90:3121-3134 | velo.f90:3155-3168 | +34 |
| L1350 | velo.f90:80-88 | velo.f90:89-97 | +9 |
| L1351 | velo.f90:115-121 | velo.f90:124-130 | +9 |
| L1352 | velo.f90:137-145 | velo.f90:146-154 | +9 |
| L1353 | velo.f90:151-164 | velo.f90:160-173 | +9 |
| L1354 | velo.f90:187-196 | velo.f90:196-205 | +9 |
| L1355 | velo.f90:204-252 | velo.f90:213-261 | +9 |
| L1356 | velo.f90:256-279 | velo.f90:265-288 | +9 |
| L1357 | velo.f90:285-294 | velo.f90:294-303 | +9 |
| L1358 | velo.f90:306-351 | velo.f90:315-360 | +9 |
| L1359 | velo.f90:355-363 | velo.f90:364-372 | +9 |
| L1360 | velo.f90:413-435 | velo.f90:422-444 | +9 |
| L1361 | velo.f90:440-483 | velo.f90:449-492 | +9 |
| L1363 | velo.f90:2891-3015 | velo.f90:2921-3045 | +30 |
| L1364 | velo.f90:1376-1400 | velo.f90:1388-1415 | +12 |
| L1365 | velo.f90:1404-1459 | velo.f90:1419-1474 | +15 |
| L1366 | velo.f90:1463-1559 | velo.f90:1478-1574 | +15 |
| L1367 | velo.f90:1858-1897 | velo.f90:1877-1919 | +19 |
| L1369 | velo.f90:1724-1730 | velo.f90:1743-1749 | +19 |
| L1370 | velo.f90:1734-1740 | velo.f90:1753-1759 | +19 |
| L1371 | velo.f90:1744-1750 | velo.f90:1763-1769 | +19 |
| L1372 | velo.f90:1770-1778 | velo.f90:1789-1797 | +19 |
| L1373 | velo.f90:1779-1787 | velo.f90:1798-1806 | +19 |
| L1374 | velo.f90:619-639 | velo.f90:631-651 | +12 |
| L1375 | velo.f90:649-653 | velo.f90:661-665 | +12 |
| L1376 | velo.f90:662-714 | velo.f90:674-726 | +12 |
| L1377 | velo.f90:720-772 | velo.f90:732-784 | +12 |
| L1378 | velo.f90:778-830 | velo.f90:790-842 | +12 |
| L1379 | velo.f90:970-978 | velo.f90:982-990 | +12 |
| L1380 | velo.f90:981-987 | velo.f90:993-999 | +12 |
| L1381 | velo.f90:992-1000 | velo.f90:1004-1012 | +12 |
| L1382 | velo.f90:1006-1014 | velo.f90:1018-1026 | +12 |
| L1383 | velo.f90:1020-1028 | velo.f90:1032-1040 | +12 |
| L1385 | velo.f90:896-903 | velo.f90:908-915 | +12 |
| L1386 | velo.f90:917-924 | velo.f90:929-936 | +12 |
| L1387 | velo.f90:938-945 | velo.f90:950-957 | +12 |
| L1388 | velo.f90:1040-1046 | velo.f90:1052-1058 | +12 |
| L1389 | velo.f90:1048-1054 | velo.f90:1060-1066 | +12 |
| L1390 | velo.f90:1081-1198 | velo.f90:1093-1210 | +12 |
| L1391 | velo.f90:1243-1253 | velo.f90:1255-1265 | +12 |
| L1392 | velo.f90:1270-1302 | velo.f90:1282-1314 | +12 |
| L1393 | velo.f90:1306-1337 | velo.f90:1318-1349 | +12 |
| L1394 | velo.f90:1603-1609 | velo.f90:1622-1628 | +19 |
| L1395 | velo.f90:1613-1619 | velo.f90:1632-1638 | +19 |
| L1396 | velo.f90:1623-1629 | velo.f90:1642-1648 | +19 |
| L1397 | velo.f90:1648-1656 | velo.f90:1667-1675 | +19 |
| L1398 | velo.f90:1657-1665 | velo.f90:1676-1684 | +19 |
| L1399 | velo.f90:514-550 | velo.f90:523-562 | +9 |
| L1400 | velo.f90:3347-3368 | velo.f90:3381-3402 | +34 |
| L1401 | velo.f90:3378-3408 | velo.f90:3412-3442 | +34 |
| L1402 | velo.f90:3414-3450 | velo.f90:3448-3484 | +34 |
| L1452 | wall.f90:319-339 | wall.f90:331-351 | +12 |
| L1470 | wall.f90:3635-3677 | wall.f90:3647-3689 | +12 |
| L1471 | wall.f90:3681-3717 | wall.f90:3693-3729 | +12 |
| L1485 | wall.f90:888-947 | wall.f90:900-959 | +12 |
| L1486 | wall.f90:109-119 | wall.f90:118-128 | +9 |
| L1487 | wall.f90:124-137 | wall.f90:133-146 | +9 |
| L1488 | wall.f90:144-185 | wall.f90:153-194 | +9 |
| L1489 | wall.f90:192-194 | wall.f90:201-203 | +9 |
