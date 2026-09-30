# S4d map: DIVERGENCE_PART_1 (FDS `divg.f90:22-785`)

Status: static map only (read of the FireX reference tree `src/Source/divg.f90`, `func.f90`, `cons.f90`, `type.f90` plus the Legacy Mapper
survey `docs/inventory/gpu_callgraph_survey.md`, `gpu_candidate_loops.csv`, `gpu_callee_routines.csv`). Nothing was compiled or run for this
map. Line counts are counted from the source; time shares are the survey's model scores (no per-loop timing exists).
The scope check (section 6) stops the spike here.

## 1. Routine shape
`DIVERGENCE_PART_1(T,DT,NM)`: 764 source lines, ~50 loop nests (cell, face, wall-cell, CFACE), ~25 `CALL`s. It is a sequence of stages that
each write one 3D field and feed the next through the `WORK1..9` / `SWORK1..3` scratch fields, so the stages are ordered and the
wall-cell (boundary) loops sit between the cell loops and modify the same arrays (e.g. `RHO_D_DZDX` is written by the cell nest, patched
by the wall loop, then read by the next cell nest).

Entry state (no dummy carries a field): `POINT_TO_MESH(NM)` sets ~100 module pointers (`MESH_POINTERS`) into `MESHES(NM)`; the routine then
aliases `DP,RHOP,UU,VV,WW,ZZP,PBAR_P` to the predictor or corrector set (module-level `POINTER` variables `UU,VV,WW,RHOP,ZZP`) and
points `RTRM,CP,KP,KDTD*,RHO_D,H_RHO_D_DZD*,R_H_G,U_DOT_DEL_*` at `WORK1..9` / `SWORK1..3` (reused under different names in different stages).

## 2. Data touched
- Mesh field pointers (allocatable components of `MESHES(NM)`): `TMP,RHO,RHOS,ZZ,ZZS,U,V,W,US,VS,WS,D,DS,Q,QR,MU,RSUM,D_SOURCE,D_Z_MAX,
  DEL_RHO_D_DEL_Z,DIF_FX/FY/FZ(+S),WORK1..9,SWORK1..3,WORK_PAD,PRESSURE_ZONE,INTERPOLATED_MESH,CCVAR`; 1D metrics `R,RRN,RC,XC,ZC,DX,DY,DZ,
  RDX,RDY,RDZ,RDXN,RDYN,RDZN,RHO_0`; 2D `PBAR,PBAR_S,R_PBAR(K,IPZ)`; zone sums `DSUM,PSUM,USUM`.
- Global tables (module `allocatable` in `cons.f90`): `D_Z(0:I_MAX_TEMP,N)`, `H_SENS_Z(0:I_MAX_TEMP,N)`, `CP_Z`, `K_RSQMW_Z`, `RSQ_MW_Z`
  (I_MAX_TEMP = 5000), scalars/flags `SIM_MODE,PREDICTOR,CORRECTOR,CHECK_VN,CC_IBM,CYLINDRICAL,TENSOR_DIFFUSIVITY,STORE_SPECIES_FLUX,
  CONSTANT_SPECIFIC_HEAT_RATIO,N_TOTAL_SCALARS,N_TRACKED_SPECIES,PERIODIC_TEST,STRATIFICATION,GVEC,MU_DNS,RSC_T,RPR_T,CPOPR,GM1OG`.
- Derived types (all through pointers, in the wall/boundary loops and in three cell nests): `WALL(IW)` (`WALL_TYPE`: `BOUNDARY_TYPE`,
  `BC_INDEX`, `B1_INDEX`, `THIN`), `BOUNDARY_COORD(:)` (`II,JJ,KK,IIG,JJG,KKG,IOR`), `BOUNDARY_PROP1(:)` (`ZZ_F(:)`, `RHO_D_F(:)`, `RHO_D_DZDN_F(:)`
  allocatable components, `U_NORMAL(_S)`, `TMP_F`, `RDN`, `AREA`, `K_G`, `Q_CON_F`, `Q_LEAK`), `EXTERNAL_WALL(IW)%NIC`, `CFACE(:)`,
  `CELL(CELL_INDEX(I,J,K))%SOLID` (2-level indirection in two hot cell nests), `SPECIES_MIXTURE(N)%RCON,%SC_T_USER,%DEPOSITING`.
- Local allocatables / assumed-shape: `ZZ_GET(1:N_TRACKED_SPECIES)` (allocated per OpenMP thread, rank-1), `D_Z_N(0:I_MAX_TEMP)` (a 5001-double
  automatic array copied per species, passed to `INTERPOLATE1D_UNIFORM` with an assumed-shape lower-bound dummy), `RHO_D_DZDN_GET(1:N_TRACKED_SPECIES)`.
- Reduction-like updates: `DSUM/PSUM/USUM(IPZ)` (indexed accumulation over cells/walls/cfaces; order-sensitive), `MAXLOC`/`SUM` over species inside cell nests.

## 3. Cell nests (structured, device-friendly in principle), by survey cost
Survey model score (rank_score) and estimated share; "call" = calls a source routine. Line numbers are `divg.f90`.
| nest (lines) | what | call tree / touches | score |
|---|---|---|---|
| 128-235 (SURVEY L0365; whole diffusion block, 4 deep incl. species loop) | RHO_D, D_Z interpolation, RHO_D_DZD*, wall patch, D_Z_MAX | `INTERPOLATE1D_UNIFORM` (assumed-shape, uses `LBOUND/UBOUND`), `TENSOR_DIFFUSIVITY_MODEL` (198 lines, allocatable pointers, only if tensor diffusivity); wall loop uses WALL/BC/B1/EWC/SM; 27 module vars, 17 mesh pointers | 621.8 (7.68 % of run) |
| 287-421 (L0369; species loop) | H_RHO_D_DZD* (3 x `GET_SENSIBLE_ENTHALPY_Z` per cell), wall correction `WALL_LOOP_2`, two flux-divergence nests | `GET_SENSIBLE_ENTHALPY_Z(N,T,H)` (11 lines, module `H_SENS_Z` allocatable, `I_MAX_TEMP`); wall loop is derived-type heavy (BC, B1, WC, `SUM/MAXLOC` of `ZZ_F`) | 194.2 (2.40 %) |
| 640-660 (L0381; species advection stage) | u.grad(rho Z_n) then `DP` update with `GET_SENSIBLE_ENTHALPY_Z`, `CELL%SOLID`, `SM%RCON` | `SPECIES_ADVECTION_PART_1_NEW` (247 lines), `SPECIES_ADVECTION_PART_2` (72), `GET_SCALAR_FACE_VALUE` (158, pointer dummies, limiter SELECT), `SET_EXIMRHOZZLIM_3D` (CC only) | 120.8 (1.49 %) |
| 696-716 (L0384; MMS) | manufactured-solution source | `VD2D_MMS_Z_SRC`, `GET_SENSIBLE_ENTHALPY_Z`; only `PERIODIC_TEST==7` | 53.2 |
| 730-753 (L0385; pressure-zone sums) | `DSUM/PSUM` accumulation | `ADD_CUTCELL_PSUM/ADD_LINKEDCELL_PSUM` (CC only); order-sensitive scatter to `IPZ` | 25.8 |
| 435-444 (L0370) | CP and R_H_G per cell | `GET_SPECIFIC_HEAT(ZZ_GET,CP,T)` (module `CP_Z`, `DOT_PRODUCT` over species, per-thread `ZZ_GET`) | 5.8 |
| 463-471 (L0371) | thermal conductivity per cell | `GET_CONDUCTIVITY` (module `K_RSQMW_Z`, `RSQ_MW_Z`), `CELL(...)%SOLID` | 5.5 |
| 512-522 (L0374) | k dT/dx,y,z face fluxes | pure arithmetic, 7 arrays | 6.0 |
| 245-258 (L0366) | flux-sum correction with `MAXLOC`/`SUM` over species | arithmetic over `ZZP(...,1:N)`, species loop inside the cell loop | 6.0 |
| 499-505, 561-582, 591-629, 668-687 (L0373, L0376-L0380, L0382, L0383) | D_Z_MAX, div k grad T (Cartesian and cylindrical), `DP -= U_DOT_DEL_RHO_H_S`, RTRM scaling, source add, stratification | pure arithmetic, 2-9 arrays each | 1-3 each |
| 171-181 (inside L0365) | RHO_D * grad Z (face fluxes) | pure arithmetic, `ZZP(...,N)` 4D | part of L0365 |

Also called from the routine (outside the loops above): `ENTHALPY_ADVECTION_NEW` (182 lines: cell nest calling `GET_SENSIBLE_ENTHALPY(ZZ_GET,...)`
over a padded (-1:KBP1+1) range, `GET_SCALAR_FACE_VALUE` x3, a wall loop, and a CFACE-free correction; `TARGET` local 4x4x4 work arrays and
pointers into them), `MERGE_PRESSURE_ZONES` (host), and the CC_IBM (cut-cell) chain `CC_VELOCITY_FLUX`, `CFACE_PREDICT_NORMAL_VELOCITY`,
`SET_EXIM*`, `CC_DIVERGENCE_PART_1` (~9 kL of `ccib.f90` in the survey; out of scope when `CC_IBM` is off).

## 4. Wall / boundary loops (derived-type, irregular)
`WALL_LOOP3` (93, calls `PREDICT_NORMAL_VELOCITY`, 101 lines: `WALL`,`VENTS`,`SURFACE`, ramp evaluation, string-free but pointer chasing),
`WALL_LOOP` (192-233, 2 nested `SELECT CASE`), `WALL_LOOP_2` (318-380), `BOUNDARY_LOOP` (482), `CORRECTION_LOOP` (532-554, updates
`B1%K_G`, `DP` scatter, `KDTD*` zero), `WALL_LOOP4` (757), `CFACE_LOOP` (770). About 200 source lines in this routine, all indexed by wall
cell with data-dependent scatter into the cell fields (`RHO_D_DZDX(BC%IIG-1,...)`), `CYCLE` on boundary type, and per-thread privates.
Each is a small share of run time (0.4-3.4 by score) but they are interleaved between the cell nests.

## 5. Blockers seen for K1 and K2
1. Derived types with allocatable / pointer components (`BOUNDARY_PROP1%ZZ_F`, `RHO_D_DZDN_F`, `WALL`, `BOUNDARY_COORD`, `CELL`, `SPECIES_MIXTURE`)
   in every wall loop and in three cell nests. K2 (`declare target` + flat explicit-shape dummies) and K1 (plain structs) both need them
   flattened to struct-of-arrays (about 20 integer/real arrays, plus `ZZ_F(N)`/`RHO_D_F(N)` as 2D arrays).
2. Module state: `H_SENS_Z, CP_Z, K_RSQMW_Z, RSQ_MW_Z, D_Z` are module allocatables used inside the callees; the callees must take them as explicit-shape
   arguments (K2) or captured pointers (K1). `GET_SENSIBLE_ENTHALPY_Z`, `GET_SPECIFIC_HEAT`, `GET_CONDUCTIVITY` are 10-11 lines each and only
   read those tables: the easy callees.
3. Assumed-shape / pointer dummies and `LBOUND/UBOUND` (`INTERPOLATE1D_UNIFORM`, `GET_SCALAR_FACE_VALUE`, `SPECIES_ADVECTION_PART_2`) and per-thread
   `ALLOCATE(ZZ_GET)` inside `!$OMP PARALLEL` (not allowed in a device region; needs fixed-size local or an explicit scratch array).
4. Pointer aliasing of scratch fields (`WORK1..9`, `SWORK1..3` reused under 4+ names) and stage ordering with host-side boundary patches between nests:
   a resident-data port needs either the wall loops on the device (blocker 1) or a per-stage device/host sync of the touched fields.
5. Indexed accumulation `DSUM/PSUM/USUM(IPZ)`: order-sensitive scatter; bitwise parity needs a fixed-order two-phase gather (the S4 D-031 pattern) or
   stays on the host.
6. `TENSOR_DIFFUSIVITY_MODEL` (198 lines, alloc pointers), CC_IBM chain, `MERGE_PRESSURE_ZONES`, MMS: outside a lean spike.

## 6. Scope check
- Whole `DIVERGENCE_PART_1` on device (cell + wall loops + `ENTHALPY_ADVECTION_NEW` + `SPECIES_ADVECTION_PART_1_NEW/2` + `GET_SCALAR_FACE_VALUE` +
  `PREDICT_NORMAL_VELOCITY`, without CC_IBM and tensor diffusivity): routine 764 lines + callees 182 + 247 + 72 + 158 + 101 + ~40 (small function
  callees, `VD2D`, `INTERPOLATE1D_UNIFORM`) = about 1 560 Fortran lines to touch, ~2 100 with `TENSOR_DIFFUSIVITY_MODEL`, before the
  struct-of-arrays plumbing. That is well above the ~600-line limit and needs derived types and allocatable module state flattened across
  about 12 routines and ~20 loop nests (about 30 including the wall loops).
- **Decision: stop after the map** (no K1/K2 port, no harness, no timing). Estimate for the coordinator:
  | option | nests | callee routines | lines to change (Fortran ref + K1 + K2 ports) | blocks |
  |---|---|---|---|---|
  | full (no CC, no tensor) | ~30 | ~12 | ~1 560 F + ~1 500 K1 + ~1 600 K2 | derived types (blockers 1, 3, 4, 5) |
  | cell-only subset (wall patches stay host-side, fields resident, boundary values passed as flat arrays) | 12-14 | 3 small (`GET_SENSIBLE_ENTHALPY_Z`, `GET_SPECIFIC_HEAT`, `GET_CONDUCTIVITY`) | ~250 F + ~250 K1 + ~300 K2 + harness ~400 | blocker 2 only (tables as arguments); `CELL%SOLID` as a flat mask |
- The cell-only subset is the recommended lean version if S4d is to continue. Candidate nests: 298-311 (3 calls to `GET_SENSIBLE_ENTHALPY_Z` per cell,
  the callee test), 171-181 and 512-522 (pure arithmetic face fluxes), 561-569 (`DP` update), 435-444 (`GET_SPECIFIC_HEAT` with a per-cell
  species dot product; needs a fixed local array instead of `ALLOCATE(ZZ_GET)`), 646-656 (call plus `CELL%SOLID` mask plus `RCON`).
- Input for a bitwise reference: the S4 FDS-derived frozen fields (`fds_sb_r_1`, `fds_g8_r_1`) carry rho, Z and TMP only (no `U,V,W`, `Q`, `MU`,
  the `H_SENS_Z` and `CP_Z` tables or boundary state), so a harness would use synthetic fields with realistic ranges (tables built from
  a polynomial) unless a fuller FDS restart dump is produced.

## 7. What a fuller FDS restart dump would need for a real-input test of the seven ported nests
Status: the S4d kernels (findings section 12) ran on synthetic fields only. To replace them by FDS data, FDS would have to write, for ONE mesh at ONE step
(the predictor entry state of `DIVERGENCE_PART_1`, after `POINT_TO_MESH`), the arrays below, plus the reference outputs of the same nests for the bitwise check.
Sizes are for a 128^3 mesh, ghost-padded (0:IBAR+1)^3 = 130^3 = 2.197 M cells, 17.6 MB per real(8) 3D array, N = number of tracked species (4 in the harness; 1 species = 1 array).
The whole FDS mesh-array dump would be dominated by the 4D species arrays: (N+1) x 17.6 MB per 4D array.

| nest (divg.f90) | needs at entry (FDS name) | produced for comparison | size at 128^3 (N = 4) |
|---|---|---|---|
| 171-181 | `ZZP(0:IBAR+1,...,1:N)` (predictor `ZZS`-alias or corrector `ZZ`, whichever `ZZP` points to), `RHO_D` (WORK array holding RHO_D just before the nest), `RDXN,RDYN,RDZN` | `RHO_D_DZDX/Y/Z(...,1:N)` right after the nest (before the wall patch) | ZZP 70 MB + RHO_D 18 MB + outputs 3 x 70 MB |
| 298-311 | `TMP`, `RHO_D_DZDX/Y/Z(:,:,:,N)` AFTER the host wall patch (i.e. the state the nest actually reads), `H_SENS_Z(0:I_MAX_TEMP,1:N)` | `H_RHO_D_DZDX/Y/Z` before the wall correction at line 318 | 3 x 70 MB + TMP 18 MB + table 160 kB + outputs 3 x 18 MB |
| 435-444 | `ZZP`, `TMP`, `CP_Z(0:I_MAX_TEMP,1:N)` (only when `.NOT.CONSTANT_SPECIFIC_HEAT_RATIO`) | `CP`, `R_H_G` (WORK9) | outputs 2 x 18 MB |
| 463-471 | `ZZP`, `TMP`, `K_RSQMW_Z(0:I_MAX_TEMP,1:N)`, `RSQ_MW_Z(1:N)`, `CELL(CELL_INDEX(I,J,K))%SOLID` as an integer field (i.e. `CELL_INDEX` + `CELL(:)%SOLID`, flattened by FDS), `SIM_MODE` | `KP` (WORK4) | 1 mask 8.8 MB (int4) + output 18 MB |
| 512-522 | `TMP`, `KP` (as after 463-471 and after the LES/DNS `KP` modifications at lines 476-495), `RDXN,RDYN,RDZN` | `KDTDX/Y/Z` before the wall/`CORRECTION_LOOP` patch | outputs 3 x 18 MB |
| 561-569 | `KDTDX/Y/Z` AFTER the host wall patch and `CORRECTION_LOOP` (line 532-554), `RDX,RDY,RDZ`, `Q`, `QR`, `DP` on entry (`D_SOURCE`/`DP` state before the nest) | `DP` after the nest | 3 x 18 MB + Q, QR, DP 3 x 18 MB |
| 646-656 | `TMP`, `RSUM`, `R_H_G` (from 435-444), `DEL_RHO_D_DEL_Z(:,:,:,1:N)` AFTER the wall/flux-divergence nests, `U_DOT_DEL_RHO_Z` (per species, computed by `SPECIES_ADVECTION_*` before each species iteration), `RHOP`, `CELL%SOLID`, `SM%RCON` per species, `H_SENS_Z`, and `DP` on entry | `DP` after each species pass (N snapshots or one final) | 70 MB + about 6 x 18 MB in + DP |
Everything else the nests read is geometry: `RDX,RDY,RDZ,RDXN,RDYN,RDZN` (1D, a few kB, or `DX(I)` etc.) and `IBAR,JBAR,KBAR,N_TRACKED_SPECIES,I_MAX_TEMP`.

Minimum FDS output: one binary dump (raw Fortran-order arrays plus a small header of the integers above) of the 4D fields ZZP, RHO_D_DZD* (after patch), DEL_RHO_D_DEL_Z (after patch) and the 3D fields RHO_D, TMP,
KP, KDTD* (after patch), Q, QR, DP (entry), RSUM, R_H_G, U_DOT_DEL_RHO_Z (per species), RHOP, plus CP_Z, H_SENS_Z, K_RSQMW_Z, RSQ_MW_Z, the SOLID mask and the reference outputs listed above, taken with a
tracked-species count N equal to the case (N = 4 gives about 0.8 GB for everything listed; N = 2 about 0.5 GB). Taking the snapshots at the nest boundaries needs a small debug hook in `DIVERGENCE_PART_1` (WRITE of the arrays between the nests), which the
current FDS does not have; the frozen S4 inputs cannot serve as substitutes. A one-step dump on a small mesh (e.g. 32^3, about 15 MB at N = 4) would be enough for the bitwise test; the 128^3 size is only needed for real-data timing (table-lookup access pattern).
Which nests are chained in the real routine: 171-181 -> wall patch -> 298-311 -> wall correction; 512-522 -> patch -> 561-569; 646-656 reads `DEL_RHO_D_DEL_Z` written after the flux-divergence nests (lines 350-421).
