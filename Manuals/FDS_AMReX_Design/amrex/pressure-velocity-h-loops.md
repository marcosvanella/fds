# Pressure and H loops: candidate list for the AMR Pressure Solver Lead

Owner of the list: FDS Legacy Mapper. Reader: AMR Pressure Solver Lead. Companion of `docs/amrex/loop-work-list.md` (ranking, claim protocol, coverage) and `docs/amrex/generator-howto.md` (how to add a loop and check it).

The tables between the `GENERATED` markers are written by `tools/inventory/loop_work_list.py` from the survey (`docs/inventory/gpu_generator_loop_classes.csv`, FireX 36975d7 line numbers), the claim register (`docs/amrex/loop_claims.csv`) and the committed generator sidecar. Line numbers were checked against the survey source. Do not edit the tables by hand; send corrections to the Legacy Mapper.

## 1. What is covered

H is the total-pressure-like variable of the FDS projection step, `H = |u|^2/2 + p~/rho` (`HP` is the value on the predictor/corrector time level, `KRES` the kinetic-energy part). The loops that produce, solve for, correct with, or check H are:

| Stage | Where | Loops |
|---|---|---|
| Right-hand side of the Poisson problem | `pres.f90` PRESSURE_SOLVER_COMPUTE_RHS | L1209 (boundary arrays), L1210 to L1214 (RHS nests, 3-D and transposed layouts) |
| FFT driver | `pres.f90` PRESSURE_SOLVER_FFT | L1215 to L1218 (PRHS to HP copies), L1219 (tunnel H_BAR), L1220 to L1222 (boundary fills of H) |
| Residual checks | `pres.f90` PRESSURE_SOLVER_CHECK_RESIDUALS | L1206 to L1208 |
| Tunnel and ULMAT solvers | `pres.f90` TUNNEL_POISSON_SOLVER, ULMAT_* | L1223 to L1227, L1190 to L1204 (inventory entries only) |
| Velocity from H | `velo.f90` VELOCITY_PREDICTOR / VELOCITY_CORRECTOR, BAROCLINIC_CORRECTION | L1394 to L1396, L1369 to L1371, L1343 to L1346 (translated), L1397/L1398 and L1372/L1373 (dead code) |
| H and velocity at boundaries | `velo.f90` NO_FLUX, WALL_VELOCITY_NO_GRADH | L1364, L1365, L1366, L1400 to L1402 |
| Forcing terms in the momentum flux | `velo.f90` VELOCITY_FLUX, DIRECT_FORCE, CORIOLIS_FORCE, PATCH_VELOCITY_FLUX | L1374 to L1390 |
| Pressure-zone divergence | `divg.f90` DIVERGENCE_PART_2, MERGE_PRESSURE_ZONES, CHECK_DIVERGENCE | L0363, L0394 to L0396, L0400 |

Section D lists the `pres.f90` routines whose loops are geometry-deferred (cut-cell matrices, CFACE, CCVAR) and are outside this list.

## 2. Classification

| Class | Meaning |
|---|---|
| translatable-now | The generator front end accepts the loop as written. Missing: a sidecar entry and the bitwise test. |
| needs-feature | Translatable once a named generator feature exists; the feature is named in the table. |
| blocked | Waits for something outside the loop (neighbour-mesh data, callee flattening, CSR table). |
| claimed | Another owner holds it; shown in the Owner column. |
| translated | Kernel in the committed sidecar; the kernel name is in the table. |
| test-only | Executed only for a manufactured-solution test or guarded by dead code. Not planned. |

The `Read` column: `R` the loop text was read for this list, `I` inventory entry only (read before claiming).

## 3. Notes for the Pressure Solver Lead

- **The FFT path may disappear.** The Poisson solve itself (`pois.f90`, 205 loops) is retired in the survey: FR-037 keeps one FFT per box, and the AMReX FFT Poisson solver may replace it. The loops around the solve (L1210 to L1222) exist only if the FDS driver is kept. Check `docs/pressure/00-fds-pressure-baseline.md` and `docs/pressure/01-amr-mapping-spec.md` before investing in L1212 to L1218 (the transposed IPS layouts); if the AMReX solver takes the PRHS and HP arrays directly, those loops become layout conversions or go away.
- **ULMAT, GLMAT and the matrix routines** (section D, plus L1190 to L1204) are solver internals or geometry-deferred and are left to AMReX.
- **L1205 COMPUTE_VELOCITY_ERROR** is classed geometry-deferred only because its body is guarded by `IF (CC_IBM)` (`pres.f90`, near line 862) and it reads neighbour-mesh data. It is a diagnostic; it is not on this list.
- **Dependencies on other owners.** L1209 and L1364 need neighbour-mesh data (`DX_OTHER`, exchange buffers): GPU Mesh Data Loops Engineer. L1402 and L0375 are with the GPU Wall Loops Engineer. L0394 is with the AMR Solid Phase Lead as a CSR consumer. L1375, L1392, L1393 are with the GPU Generator Engineer (EVALUATE_RAMP callee, edge tables).
- **Dead and test-only code.** L1397/L1398 and L1372/L1373 sit behind `PERIODIC_TEST==7 .AND. .FALSE.` (`velo.f90:1647`, `velo.f90:1769`); L1388 and L1389 only run when `PERIODIC_TEST==7` (`velo.f90:842`). Their modelled shares (0.049 to 0.325 %) are not real; skip them.
- **NO_FLUX L1365 and L1366.** L1365 (`velo.f90:1404-1459`) writes the same bits to a face from every obstruction that contains it, because the value depends on the face only; this fits the IDEMPOTENT marker, with `pvf_fvx` as the nearest kernel. L1366 (`velo.f90:1463-1559`) writes faces through wall subscripts; a thin pair of walls (IOR=+1 at one cell, IOR=-1 at the next) can name the same face, so the last wall index would win. It needs a uniqueness proof or a gather before it is translated; it is claimed by the GPU Mesh Data Loops Engineer.

## 4. Suggested order

1. **Accepted by the front end today:** L1211 (RHS, `pres.f90:250-260`), L1207 (P = RHOP*(HP-KRES), `pres.f90:758-764`, same formula as the translated baroclinic kernel `baro_p_rrho`), L1215 (PRHS to HP copy), L1379 (Coriolis cell-centred velocities), L1391 (cylindrical vorticity/stress). Each needs a sidecar entry and a bitwise test.
2. **Constant-subscript forcing (one feature):** L1385 to L1387 and L1381 to L1383 index `FVEC(n)` and `OVEC(n)` with constants; the feature is "rank-1 array subscript is a constant, not a loop variable". One piece of generator work unlocks six loops.
3. **Layout contract (one feature):** L1212 to L1214 and L1216 to L1218 store through permuted subscripts; one layout feature unlocks six loops, but only if the FFT driver is kept (see the first note).
4. **Nest shape (one feature):** L1210 and L1220 to L1222 are not K,J,I nests; the same feature unlocks L1392 and L1393 (GPU Generator Engineer).
5. **L1365** after step 1 (closest kernel `pvf_fvx`).
6. **L1366** after the uniqueness proof.
7. **L1209 and L1364** after the neighbour-mesh table exists.

Share of the whole list: the pressure and H loops are small in the model (the largest solver-side loop, L1209, is 0.008 %; the matrix and Poisson internals carry the weight and are not in scope). Prioritise by dependency (steps 2 to 4 unlock several loops each), not by share.

<!-- GENERATED-BEGIN:pvh -->
Counts over sections A to C (77 loops): blocked 12, claimed 12, needs-feature 20, test-only 6, translatable-now 2, translated 25.

### A. pres.f90, translation-eligible, non-geometry

| Loop | file:lines | Routine | What | Share % | Class | Feature or blocker | Closest kernel | Owner | Read |
|---|---|---|---|---|---|---|---|---|---|
| L1209 | pres.f90:65-228 | PRESSURE_SOLVER_COMPUTE_RHS | Poisson boundary arrays BXS..BZF from wall data (Neumann, Dirichlet, interpolated, open) | 0.008 | claimed | MESHES(NOM)%DX(EWC%IIO_MIN), VENTS%PRESSURE_RAMP_INDEX, EVALUATE_RAMP, synthetic-eddy arrays | wall_us_pred (wall gather); no neighbour-mesh model yet | AMR Pressure Backend Implementer | R |
| L1210 | pres.f90:238-245 | PRESSURE_SOLVER_COMPUTE_RHS | Poisson RHS, cylindrical K,I nest | 0.000 | needs-feature | front end takes K,J,I nests only (same as L1392/L1393) | vflux_fvx (flux differences) | open — available | R |
| L1211 | pres.f90:250-260 | PRESSURE_SOLVER_COMPUTE_RHS | Poisson RHS PRHS (3-D, IPS 1/4/7) | 0.002 | claimed | front end accepts; no bitwise test committed | dp_kdtd (divg.f90:561-570, same flux-difference form) | AMR Pressure Backend Implementer | R |
| L1212 | pres.f90:267-277 | PRESSURE_SOLVER_COMPUTE_RHS | Poisson RHS PRHS transposed (IPS 2: PRHS(J,I,K)) | 0.002 | needs-feature | layout contract: output subscript is not I+/-c | dp_kdtd with a permuted store | open — available | R |
| L1213 | pres.f90:282-292 | PRESSURE_SOLVER_COMPUTE_RHS | Poisson RHS PRHS transposed (IPS 3,6: PRHS(K,J,I)) | 0.002 | needs-feature | layout contract | dp_kdtd with a permuted store | open — available | R |
| L1214 | pres.f90:297-307 | PRESSURE_SOLVER_COMPUTE_RHS | Poisson RHS PRHS transposed (IPS 5: PRHS(I,K,J)) | 0.002 | needs-feature | layout contract | dp_kdtd with a permuted store | open — available | R |
| L1215 | pres.f90:394-400 | PRESSURE_SOLVER_FFT | copy PRHS to HP after the FFT solve (IPS 1,4,7) | 0.000 | claimed | front end accepts; whether the loop survives depends on the solver driver (FR-037) | vcorr_u (plain cell copy form) | AMR Pressure Backend Implementer | R |
| L1216 | pres.f90:404-410 | PRESSURE_SOLVER_FFT | copy PRHS to HP (IPS 2: PRHS(J,I,K)) | 0.000 | needs-feature | layout contract | as L1215 | open — available | R |
| L1217 | pres.f90:414-420 | PRESSURE_SOLVER_FFT | copy PRHS to HP (IPS 3,6) | 0.000 | needs-feature | layout contract | as L1215 | open — available | R |
| L1218 | pres.f90:424-430 | PRESSURE_SOLVER_FFT | copy PRHS to HP (IPS 5) | 0.000 | needs-feature | layout contract | as L1215 | open — available | R |
| L1219 | pres.f90:438-440 | PRESSURE_SOLVER_FFT | tunnel preconditioner: add H_BAR to HP (loop over I, section assignment) | 0.000 | needs-feature | inside !$OMP MASTER; only with TUNNEL_PRECONDITIONER | none | open — available | R |
| L1220 | pres.f90:450-462 | PRESSURE_SOLVER_FFT | H boundary fill in x (LBC/MBC/NBC code tests) | 0.000 | claimed | front end takes K,J,I nests only | wall_hs_bt (boundary store) | AMR Pressure Backend Implementer | R |
| L1221 | pres.f90:466-477 | PRESSURE_SOLVER_FFT | H boundary fill in y | 0.000 | claimed | front end takes K,J,I nests only | as L1220 | AMR Pressure Backend Implementer | R |
| L1222 | pres.f90:481-492 | PRESSURE_SOLVER_FFT | H boundary fill in z | 0.000 | claimed | front end takes K,J,I nests only | as L1220 | AMR Pressure Backend Implementer | R |
| L1223 | pres.f90:531-537 | TUNNEL_POISSON_SOLVER | (inventory only: cell/face loop; wall gather) | 0.000 | blocked | needs neighbour-mesh or ragged per-wall data (inventory) | wall_up_ghost, wall_us_pred | open — available | I |
| L1224 | pres.f90:544-567 | TUNNEL_POISSON_SOLVER | (inventory only: cell/face loop; wall gather; non-perfect nest (outer loop)) | 0.007 | blocked | needs neighbour-mesh or ragged per-wall data (inventory) | wall_up_ghost, wall_us_pred | open — available | I |
| L1225 | pres.f90:576-582 | TUNNEL_POISSON_SOLVER | (inventory only: cell/face loop; wall gather; non-perfect nest (outer loop)) | 0.000 | blocked | needs neighbour-mesh or ragged per-wall data (inventory) | wall_up_ghost, wall_us_pred | open — available | I |
| L1226 | pres.f90:583-589 | TUNNEL_POISSON_SOLVER | (inventory only: cell/face loop; wall gather; non-perfect nest (outer loop)) | 0.000 | blocked | needs neighbour-mesh or ragged per-wall data (inventory) | wall_up_ghost, wall_us_pred | open — available | I |
| L1227 | pres.f90:591-596 | TUNNEL_POISSON_SOLVER | (inventory only: cell/face loop; wall gather; non-perfect nest (outer loop)) | 0.000 | blocked | needs neighbour-mesh or ragged per-wall data (inventory) | wall_up_ghost, wall_us_pred | open — available | I |
| L1206 | pres.f90:729-742 | PRESSURE_SOLVER_CHECK_RESIDUALS | CHECK_POISSON residual (diagnostic) | 0.001 | needs-feature | layout = exact: WORK8(1:IBAR,..) has no ALLOCATE line; diagnostic only | dp_kdtd | open — available | R |
| L1207 | pres.f90:758-764 | PRESSURE_SOLVER_CHECK_RESIDUALS | P = RHOP*(HP-KRES) over the full box | 0.000 | claimed | front end accepts; no bitwise test committed | baro_p_rrho (velo.f90:3254-3261, same formula) | AMR Pressure Backend Implementer | R |
| L1208 | pres.f90:768-787 | PRESSURE_SOLVER_CHECK_RESIDUALS | inseparable Poisson residual (ITERATE_BAROCLINIC_TERM) | 0.001 | needs-feature | layout = exact: WORK8 view; MAXVAL/MAXLOC after the loop stay on the host | baro_fvx, cfl_max (reduction) | open — available | R |
| L1190 | pres.f90:1686-1694 | ULMAT_SOLVE_ZONE | (inventory only: cell/face loop; wall gather; derived-type designator table; ) | 0.001 | blocked | needs neighbour-mesh or ragged per-wall data (inventory) | wall_up_ghost, wall_us_pred | open — available | I |
| L1192 | pres.f90:1713-1720 | ULMAT_SOLVE_ZONE | (inventory only: cell/face loop; wall gather; zone table) | 0.001 | blocked | needs neighbour-mesh or ragged per-wall data (inventory) | wall_up_ghost, wall_us_pred | open — available | I |
| L1194 | pres.f90:1782-1790 | ULMAT_SOLVE_ZONE | (inventory only: cell/face loop; wall gather; derived-type designator table; ) | 0.001 | blocked | needs neighbour-mesh or ragged per-wall data (inventory) | wall_up_ghost, wall_us_pred | open — available | I |
| L1196 | pres.f90:1808-1815 | ULMAT_SOLVE_ZONE | (inventory only: cell/face loop; wall gather; zone table) | 0.001 | blocked | needs neighbour-mesh or ragged per-wall data (inventory) | wall_up_ghost, wall_us_pred | open — available | I |
| L1198 | pres.f90:1837-1846 | ULMAT_SOLVE_ZONE | (inventory only: cell/face loop; wall gather; derived-type designator table; ) | 0.002 | blocked | needs neighbour-mesh or ragged per-wall data (inventory) | wall_up_ghost, wall_us_pred | open — available | I |
| L1203 | pres.f90:1923-1925 | ULMAT_SOLVE_ZONE | (inventory only: cell/face loop; non-perfect nest (outer loop)) | 0.000 | needs-feature | loop nest ['I'] at line 1923 is not a perfect K,J,I nest, nor an index loop (N) around one (collapse contract) (inventory) | zz_corr (species N outermost) | open — available | I |
| L1204 | pres.f90:1931-1994 | ULMAT_SOLVE_ZONE | (inventory only: wall loop; wall gather; derived-type designator table; zone ) | 0.003 | blocked | needs neighbour-mesh or ragged per-wall data (inventory) | wall_up_ghost, wall_us_pred | open — available | I |
| L1154 | pres.f90:2031-2047 | ULMAT_GET_H_REGFACES | LOG_INTWC marks, ULMAT variant (set-up) | 0.001 | needs-feature | P4; host acceptable | wall_up_ghost | open — available | I |
| L1121 | pres.f90:4012-4026 | CHECK_UNSUPPORTED_MESH | Dirichlet-set check over meshes (set-up check) | 0.086 | blocked | mesh loop with process test; host | none | open — available | R |
| L1122 | pres.f90:4351-4356 | COPY_CCVAR_IN_HS | HS from WALL boundary type for IS_WALLT (CC_IBM helper) | 0.001 | translated | translated | kernel wall_hs_bt | - | R |
| L1130 | pres.f90:5580-5601 | GET_H_REGFACES | LOG_INTWC marks for internal solid faces (set-up) | 0.001 | needs-feature | P4: rank-4 LOGICAL array has no table kind; host acceptable | wall_up_ghost | open — available | R |

### B. velo.f90, H assembly, no-flux, velocity update and forcing

| Loop | file:lines | Routine | What | Share % | Class | Feature or blocker | Closest kernel | Owner | Read |
|---|---|---|---|---|---|---|---|---|---|
| L1374 | velo.f90:619-639 | VELOCITY_FLUX | vorticity and stress tensor (VELOCITY_FLUX) | 0.079 | translated | translated | kernel vflux_vort_tau | - | R |
| L1375 | velo.f90:649-653 | VELOCITY_FLUX | gravity components GX/GY/GZ by EVALUATE_RAMP | 0.001 | claimed | callee is a function of a ramp table (O1) | vflux_fvx host_inputs GX/GY/GZ | GPU Generator Engineer | R |
| L1376 | velo.f90:662-714 | VELOCITY_FLUX | (inventory only: cell/face loop; wall gather; edge tables) | 0.226 | translated | translated | kernel vflux_fvx | GPU Generator Engineer | I |
| L1377 | velo.f90:720-772 | VELOCITY_FLUX | (inventory only: cell/face loop; wall gather; edge tables) | 0.226 | translated | translated | kernel vflux_fvy | GPU Generator Engineer | I |
| L1378 | velo.f90:778-830 | VELOCITY_FLUX | (inventory only: cell/face loop; wall gather; edge tables) | 0.226 | translated | translated | kernel vflux_fvz | GPU Generator Engineer | I |
| L1385 | velo.f90:896-903 | DIRECT_FORCE | DIRECT_FORCE: FVX -= RRHO*FVEC(1)*ramp*SIN_THETA | 0.011 | needs-feature | rank-1 array subscript is a constant | vcorr_u | open — available | R |
| L1386 | velo.f90:917-924 | DIRECT_FORCE | DIRECT_FORCE: FVY update | 0.011 | needs-feature | as L1385 | vcorr_v | open — available | R |
| L1387 | velo.f90:938-945 | DIRECT_FORCE | DIRECT_FORCE: FVZ update | 0.011 | needs-feature | as L1385 | vcorr_w | open — available | R |
| L1379 | velo.f90:970-978 | CORIOLIS_FORCE | CORIOLIS_FORCE: cell-centred velocities UP,VP,WP | 0.016 | translatable-now | front end accepts; no bitwise test committed | up_deardorff, wall_coriolis_ghost | open — available | R |
| L1380 | velo.f90:981-987 | CORIOLIS_FORCE | CORIOLIS_FORCE: ghost velocities at external walls | 0.005 | translated | translated | kernel wall_coriolis_ghost | - | R |
| L1381 | velo.f90:992-1000 | CORIOLIS_FORCE | CORIOLIS_FORCE: FVX += 2(OVEC(2)WBAR-OVEC(3)VBAR) | 0.016 | needs-feature | rank-1 array subscript is a constant, not a loop variable | vcorr_u (face-loop FVX update) | open — available | R |
| L1382 | velo.f90:1006-1014 | CORIOLIS_FORCE | CORIOLIS_FORCE: FVY update | 0.016 | needs-feature | as L1381 | vcorr_v | open — available | R |
| L1383 | velo.f90:1020-1028 | CORIOLIS_FORCE | CORIOLIS_FORCE: FVZ update | 0.016 | needs-feature | as L1381 | vcorr_w | open — available | R |
| L1388 | velo.f90:1040-1046 | MMS_VELOCITY_FLUX | manufactured-solution force in FVX (PERIODIC_TEST==7 only) | 0.315 | test-only | called only if PERIODIC_TEST==7 (velo.f90:842) | none | open — available | R |
| L1389 | velo.f90:1048-1054 | MMS_VELOCITY_FLUX | manufactured-solution force in FVZ (PERIODIC_TEST==7 only) | 0.325 | test-only | called only if PERIODIC_TEST==7 (velo.f90:842) | none | open — available | R |
| L1390 | velo.f90:1081-1198 | PATCH_VELOCITY_FLUX | PATCH_VELOCITY_FLUX: velocity patch forcing, three face nests | 1.270 | translated | translated | kernel pvf_fvx;pvf_fvy;pvf_fvz;pvf_fvx;pvf_fvy;pvf_fvz | GPU Mesh Data Loops Engineer | R |
| L1391 | velo.f90:1243-1253 | VELOCITY_FLUX_CYLINDRICAL | cylindrical vorticity and stress (OMY, TXZ) | 0.026 | translatable-now | front end accepts; no bitwise test committed | vflux_vort_tau (velo.f90:619-639) | open — available | R |
| L1392 | velo.f90:1270-1302 | VELOCITY_FLUX_CYLINDRICAL | cylindrical edge nest (K,I with J fixed) | 0.014 | claimed | front end takes K,J,I nests only (O1, deferred) | vflux_fvx | GPU Generator Engineer | I |
| L1393 | velo.f90:1306-1337 | VELOCITY_FLUX_CYLINDRICAL | cylindrical edge nest (K,I with J fixed) | 0.014 | claimed | as L1392 | vflux_fvx | GPU Generator Engineer | I |
| L1364 | velo.f90:1376-1400 | NO_FLUX | NO_FLUX: HP ghost value from the neighbour mesh (average over the overlap) | 0.426 | claimed | needs the exchange-buffer layout (ADR-001); the sum keeps its K,J,I order | none (new); pvf_* for the host rule | GPU Mesh Data Loops Engineer | R |
| L1365 | velo.f90:1404-1459 | NO_FLUX | NO_FLUX: FVX/FVY/FVZ driven to zero inside obstructions (three face nests per obstruction) | 0.347 | needs-feature | outer loop over obstructions is not a K,J,I nest; overlapping obstructions write the same bits to a face (value depends on the face only), so IDEMPOTENT fits | pvf_fvx (box and SOLID mask, host launches one kernel per box), vpred_us (HP gradient form) | open — available | R |
| L1366 | velo.f90:1463-1559 | NO_FLUX | NO_FLUX: FVX/FVY/FVZ at walls toward the specified normal velocity | 0.056 | claimed | two wall records can name one face (a thin pair writes FVX(II) from IOR=1 and from IOR=-1 of the next cell): the last wall index would win; the owner must prove or gather | wall_us_pred, wall_un_store (velo.f90:3347-3408) | GPU Mesh Data Loops Engineer | R |
| L1394 | velo.f90:1603-1609 | VELOCITY_PREDICTOR | (inventory only: cell/face loop) | 0.005 | translated | translated | kernel vpred_us | - | I |
| L1395 | velo.f90:1613-1619 | VELOCITY_PREDICTOR | (inventory only: cell/face loop) | 0.005 | translated | translated | kernel vpred_vs | - | I |
| L1396 | velo.f90:1623-1629 | VELOCITY_PREDICTOR | (inventory only: cell/face loop) | 0.005 | translated | translated | kernel vpred_ws | - | I |
| L1397 | velo.f90:1648-1656 | VELOCITY_PREDICTOR | manufactured-solution U in VELOCITY_PREDICTOR (dead code) | 0.049 | test-only | guard is PERIODIC_TEST==7 .AND. .FALSE. (velo.f90:1647): never executes | none | open — available | R |
| L1398 | velo.f90:1657-1665 | VELOCITY_PREDICTOR | manufactured-solution W in VELOCITY_PREDICTOR (dead code) | 0.049 | test-only | guard is PERIODIC_TEST==7 .AND. .FALSE. (velo.f90:1647): never executes | none | open — available | R |
| L1369 | velo.f90:1724-1730 | VELOCITY_CORRECTOR | (inventory only: cell/face loop) | 0.005 | translated | translated | kernel vcorr_u | - | I |
| L1370 | velo.f90:1734-1740 | VELOCITY_CORRECTOR | (inventory only: cell/face loop) | 0.005 | translated | translated | kernel vcorr_v | - | I |
| L1371 | velo.f90:1744-1750 | VELOCITY_CORRECTOR | (inventory only: cell/face loop) | 0.005 | translated | translated | kernel vcorr_w | - | I |
| L1372 | velo.f90:1770-1778 | VELOCITY_CORRECTOR | manufactured-solution U in VELOCITY_CORRECTOR (dead code) | 0.049 | test-only | guard is PERIODIC_TEST==7 .AND. .FALSE. (velo.f90:1769): never executes | none | open — available | R |
| L1373 | velo.f90:1779-1787 | VELOCITY_CORRECTOR | manufactured-solution W in VELOCITY_CORRECTOR (dead code) | 0.049 | test-only | guard is PERIODIC_TEST==7 .AND. .FALSE. (velo.f90:1769): never executes | none | open — available | R |
| L1343 | velo.f90:3254-3261 | BAROCLINIC_CORRECTION | (inventory only: cell/face loop) | 0.011 | translated | translated | kernel baro_p_rrho | - | I |
| L1344 | velo.f90:3267-3275 | BAROCLINIC_CORRECTION | (inventory only: cell/face loop) | 0.011 | translated | translated | kernel baro_fvx | - | I |
| L1345 | velo.f90:3282-3290 | BAROCLINIC_CORRECTION | (inventory only: cell/face loop) | 0.011 | translated | translated | kernel baro_fvy | - | I |
| L1346 | velo.f90:3297-3305 | BAROCLINIC_CORRECTION | (inventory only: cell/face loop) | 0.011 | translated | translated | kernel baro_fvz | - | I |
| L1400 | velo.f90:3347-3368 | WALL_VELOCITY_NO_GRADH | (inventory only: wall loop; wall gather; live-out scalars (PRIVATE list)) | 0.014 | translated | translated | kernel wall_un_store | - | I |
| L1401 | velo.f90:3378-3408 | WALL_VELOCITY_NO_GRADH | (inventory only: wall loop; wall gather; live-out scalars (PRIVATE list)) | 0.021 | translated | translated | kernel wall_us_pred | - | I |
| L1402 | velo.f90:3414-3450 | WALL_VELOCITY_NO_GRADH | WALL_VELOCITY_NO_GRADH: U/V/W at the wall face from the stored normal velocity | 0.022 | translated | translated | kernel wall_uvw_nograd | GPU Wall Loops Engineer | I |

### C. divg.f90, pressure-related (DIVERGENCE_PART_2, zones, CHECK_DIVERGENCE)

| Loop | file:lines | Routine | What | Share % | Class | Feature or blocker | Closest kernel | Owner | Read |
|---|---|---|---|---|---|---|---|---|---|
| L0400 | divg.f90:1301-1319 | MERGE_PRESSURE_ZONES | MERGE_PRESSURE_ZONES: CONNECTED_ZONES flags | 0.037 | translated | translated | kernel wall_connect_zones | - | I |
| L0394 | divg.f90:1574-1604 | DIVERGENCE_PART_2 | DIVERGENCE_PART_2: DP correction at solid cells from walls | 0.042 | translated | translated | kernel wall_bc_dp | GPU Wall Loops Engineer | I |
| L0395 | divg.f90:1611-1618 | DIVERGENCE_PART_2 | DIVERGENCE_PART_2: DIV predictor | 0.012 | translated | translated | kernel div2_pred | - | R |
| L0396 | divg.f90:1621-1630 | DIVERGENCE_PART_2 | DIVERGENCE_PART_2: DIV corrector | 0.012 | translated | translated | kernel div2_corr | - | R |
| L0363 | divg.f90:1675-1711 | CHECK_DIVERGENCE | CHECK_DIVERGENCE: residual and divergence extrema with locations | 0.272 | translated | translated | kernel div_extrema | - | R |

### D. pres.f90 loops outside the work list: geometry-deferred (cut-cell, CFACE, CCVAR)

| Routine | Loop ids | Share % |
|---|---|---|
| COPY_H_OMESH_TO_MESH | L1123, L1124 | 0.224 |
| GET_BCS_H_MATRIX | L1127 | 0.388 |
| GET_H_MATRIX | L1128 | 5.193 |
| GET_MATRIXGRAPH_H_WHLDOM | L1137, L1138, L1139, L1140 | 5.162 |
| GET_MATRIX_INDEXES_H | L1142 | 0.574 |
| GLMAT_SOLVER | L1143 | 1.321 |
| GLMAT_SOLVER_SETUP | L1144, L1145 | 0.134 |
| PRESSURE_SOLVER_CHECK_RESIDUALS_U | L1147, L1148 | 0.177 |
| SET_CCVAR_CGSC_H | L1149 | 0.000 |
| ULMAT_BCS_H_MATRIX | L1153 | 0.014 |
| ULMAT_H_MATRIX | L1161, L1162, L1163 | 0.062 |
| ULMAT_MATRIXGRAPH_H | L1169, L1170, L1171, L1172, L1173, L1174, L1175, L1176, L1177, L1178, L1179 | 0.155 |
| ULMAT_SOLVER_SETUP | L1182, L1183 | 0.072 |
| ULMAT_SOLVE_ZONE | L1186, L1187, L1188, L1189, L1191, L1193, L1195, L1197, L1199, L1201, L1202 | 0.015 |
| COMPUTE_VELOCITY_ERROR | L1205 | 0.212 |
<!-- GENERATED-END:pvh -->

## 5. Regenerating

`python3 tools/inventory/loop_work_list.py` rewrites the tables above (and `loop-work-list.md`, `loop_work_list.csv`). `--check` exits 1 if a committed file is stale.
