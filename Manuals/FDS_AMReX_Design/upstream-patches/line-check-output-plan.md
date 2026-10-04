# Line-number check of the output plan and the output patches (driver 0011 and 0012, UP-0011)

Scope: `Source/driver/notes/output-plan.md` (blob `04ae981082`, the same text in the working copy and in `git show HEAD:...`) and the patches it leads to. Reference tree: FDS-AMReX HEAD `2c639da09b`. `Source/main.f90`, `dump.f90` and `func.f90` have the same blobs at HEAD as at the commit that added the plan (`04ae981082`), so every "current HEAD" line number of the plan is checked against the file it was written from. `hvac.f90` is named in the plan's method text only and carries no citation. Method as in `line-check-0010-0011.md`: `git show <rev>:<path>` copies in a scratch folder, `grep -n` and `sed -n` for each citation, `patch -p1 --dry-run` for each hunk. Nothing was built.

## What exists and what does not
- The plan names its two output steps "patch 0010" (step A, level-0 dump loop in `FDS_SETUP(MODE=3)`) and "patch 0011" (step B, output-mesh set). Under the numbering scheme (A-63, `README.md` Process rules) they are driver patches **0011** and **0012**. **Neither patch file exists** in `Source/driver/patches/` (0001 to 0009 only), in `upstream-patches/` (0010 is the level-0 guard; the others are `UP-NNNN`), or on any branch. No source contains `FDS_HOOK_SET_MESH_DUMPS`, `OUT_MESH` or `N_OUT_MESHES`. So there is no hunk to dry-run for 0011 and 0012; this check covers the plan, gives Role 1 the exact anchors, and flags what the hunks must replicate.
- `UP-0011` (pressure RHS dump) is **not** an output patch of this plan. It is re-checked here only because the request names it: it still applies (section 4).
- The plan's own wording of the numbers (lines 20 to 29, 45, 46, 50 of the plan: "patch 0010", "patch 0011", "numbers from 0010") predates A-63 and should read 0011 and 0012. Patches 0006 and 0007 are cited correctly.

## 1. Anchors Role 1 needs (FDS-AMReX HEAD)
Step A (driver 0011): run the per-mesh dump loop in `FDS_SETUP(MODE=3)`.
| What | File:line (HEAD) | Content / note |
|---|---|---|
| MODE=3 block | main.f90:120-142 | `IF (MODE==3) THEN` (120) ... `RETURN` (141), `ENDIF` (142); inside `#ifdef WITH_AMREX` (119-147) |
| Insertion point | main.f90:131-132 | between `CALL UPDATE_CONTROLS(T,DT,CTRL_STOP_STATUS,.FALSE.)` (131) and `CALL DUMP_GLOBAL_OUTPUTS` (132): the same place as in `MAIN_LOOP` |
| Reference order in MAIN_LOOP | main.f90:1209 `UPDATE_CONTROLS`; 1212-1215 (`WITH_HDF5` only) `INITIALIZE_VTKHDF_FILES`, `MPI_BARRIER`; 1216-1221 `N_WRITTEN=0` and the five `WROTE_*` resets; 1222-1225 `DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX` / `CALL DUMP_MESH_OUTPUTS(T,DT,NM,.FALSE.)` / `N_WRITTEN=N_WRITTEN+1` / `ENDDO`; 1226-1235 (`WITH_HDF5` only) the `FAKEWRITE` padding loop; 1238 `DUMP_GLOBAL_OUTPUTS` | the new block copies 1216-1225 (and 1226-1235 inside `#ifdef WITH_HDF5` if VTK is wanted) |
| Variables | `N_WRITTEN` is declared in `FDS_SETUP` (main.f90:114); `WROTE_SL3D`, `_SL2D`, `_SMOKE3D`, `_BNDF`, `_PART` are module variables (cons.f90:290-294); `MESHES_PER_PROCESS` main.f90:113 | no new declaration needed except the flag |
| Not replicated today | main.f90:1210 `IF (CTRL_STOP_STATUS) STOP_STATUS = CTRL_STOP` has no counterpart in MODE=3 (131) | not part of step A, noted |
| Flag, Fortran side | fds_amrex_hooks.f90:29 (`STEP_OUTPUTS_PATCHED`, put the new `LOGICAL, PUBLIC, SAVE` after it), 39 (`PUBLIC ::` list), 71-75 (`FDS_HOOK_SET_STEP`, put `FDS_HOOK_SET_MESH_DUMPS` after 75); main.f90:15 (`USE FDS_AMREX_HOOKS, ONLY: OUT_T,OUT_DT,OUT_ICYC,STEP_OUTPUTS_PATCHED`) | add the flag to the `USE` list |
| Driver side | TimeLoop.cpp:52-53 (declarations `fds_hook_step_outputs`, `fds_hook_set_step`), 339 (`fds_outputs = fds_hook_step_outputs() != 0 ...`), 1644 (`void outputs(const StepRecord& r)`), 1648 `if (fds_outputs)`, 1651 `fds_hook_set_step(...)`, 1653 `fds_setup(3, "", &dtd)` | set the new flag once before the loop (near 339) or before 1653 |
| Writers called by `DUMP_MESH_OUTPUTS` | dump.f90:77-287; `DUMP_PART` 169, `DUMP_ISOF` 177, `DUMP_SMOKE3D` 185, `DUMP_SLCF` 193/201/217, `DUMP_BNDF` 209, `DUMP_PROF` 227, `DUMP_UVW` 237/240, `DUMP_TMP` 250, `DUMP_SPEC` 259, `DUMP_MMS` 270, `DUMP_ROTCUBE_MMS` 276, `SANDIA_OUT` 282 | writers run only if `WRITE_SMV` (167) except the CSV dumps from 226 on |

Step B (driver 0012): output-mesh set. The edits the plan lists, with the real counts.
| Site | File:line (HEAD) | Note |
|---|---|---|
| File-name and unit tables | dump.f90:405-451 | **35** `ALLOCATE` lines sized by `NMESHES` or `NMESHES+1` (20 `FN_*`, 15 `LU_*`), not "about 25"; `FN_PROF`/`LU_PROF` (394-395) are sized by `N_PROF`, not by mesh |
| `ASSIGN_FILE_NAMES` mesh loop | dump.f90:453 (`MESH_LOOP: DO NM=1,NMESHES`), 455 skips meshes of other ranks | |
| `WRITE_SMOKEVIEW_FILE` | dump.f90:2056, 2058, 2060 (`XLEVEL/YLEVEL/ZLEVEL(0:2*NMESHES)`), loops 2066, 2202, mesh loop 2457-2679 (`GRID` 2499, `PDIM` 2505, `TRNX` 2511, `OBST` 2541), write 2683-2723 (MPI-IO branch 2685-2692, gather branch 2698-2723) | called from main.f90:282 and 285 |
| `ADD_EXTERIOR_VENTS` | dump.f90:2831 (loop 2843); called at dump.f90:1817 | |
| `.smv` file-listing strings | produced in `INITIALIZE_MESH_DUMPS` (dump.f90:963), not in `WRITE_SMOKEVIEW_FILE`: `ISOG` 1114, `SMOKF3D` 1170, `SLCC/SLCT/SLCF` 1260-1264, `BNDC/BNDF` 1471-1473, `PRT5` 1506; `NMESHES` at 1519 and 1532 | flushed by `WRITE_STRINGS` (main.f90:4231; rank 0 writes `MESHES(1)` 4266-4267, then `DO NOM=2,NMESHES` 4269) |
| `WRITE_STRINGS` | main.f90:4231-4342 (send loop `LOWER_MESH_INDEX..UPPER_MESH_INDEX` 4245; reset 4337); called in MODE=3 at main.f90:135 and in MAIN_LOOP at 1242 | **missing from the plan's inventory**; depends on the mesh count |
| Output clocks and counters | func.f90:432 (`INITIALIZE_OUTPUT_CLOCKS`), 449-474 (the `SET_OUTPUT_CLOCK` calls), 496 (`SET_OUTPUT_CLOCK`); per-mesh counters are allocated `COUNTER(LOWER_MESH_INDEX:UPPER_MESH_INDEX)` at func.f90:536-538; no `NMESHES` in 432-550 | the counters are sized by the **rank's mesh range**, not by `NMESHES`; an output-mesh number outside that range is out of bounds |
| Per-mesh dump initialisation | main.f90:529 (`INITIALIZE_OUTPUT_CLOCKS`), 534-536 (`INITIALIZE_MESH_DUMPS(NM)` for the rank's meshes, skipped if `TGA_SURF_INDEX>=1`), 539 `EXCHANGE_NPATCH_INFO`, 542 `EXCHANGE_NSLICE_INFO`, 668 `EXCHANGE_NOBST_INFO` (only if `WRITE_VTK_GEOM`, 666) | `EXCHANGE_NOBST_INFO` at 296 is commented out |
| Global dump initialisation | main.f90:620-621 (`INITIALIZE_GLOBAL_DUMPS`, rank 0), dump.f90:627 | |
| Set-up order | main.f90:232 `ASSIGN_FILE_NAMES`, 282/285 `WRITE_SMOKEVIEW_FILE`, 327 `INITIALIZE_MESH_EXCHANGE_1`, 529, 535, 539, 542, 621 | the output-mesh count must exist before 232 |

## 2. Citation check of the plan (67 line-number citations, all match)
"File:line" column of section 2 (63 numbers, ranges counted once), plus the line numbers in other columns and text (405-451, 2457-2459, 2688-2722, main.f90:1223). Every routine start was read with `grep -nE "^ *SUBROUTINE <name>"`; every loop line was read with `sed`.
| Citation in the plan | Found at HEAD | Result |
|---|---|---|
| `ASSIGN_FILE_NAMES` 292, loop 453; tables 405-451 | 292; 453 `MESH_LOOP: DO NM=1,NMESHES`; 405 `ALLOCATE(FN_XYZ(NMESHES))`, 451 `ALLOCATE(FN_PART_VTK(...,NMESHES+1))` | match |
| `WRITE_SMOKEVIEW_FILE` 1766, loops 2066-2202, 2457, 2457-2459, 2688-2722 | 1766; 2066 `DO NM=1,2*NMESHES`, 2202 `YLOOP3: DO NY=1,2*NMESHES`; 2457-2459 mesh loop and rank skip; 2688-2722 MPI-IO and gather write | match (the range 2688-2722 spans both branches; the MPI-IO branch starts at 2685) |
| `ADD_EXTERIOR_VENTS` 2831 (2843) | 2831; 2843 `MESH_LOOP: DO NM=1,NMESHES` | match |
| `WRITE_STL_FILE` 1575 (1615, 1627) | 1575; both `DO NM=1,NMESHES` | match |
| `INITIALIZE_DIAGNOSTIC_FILE` 3037 (3085, 3091, 3782) | 3037; the three `NMESHES` loops | match |
| `EXCHANGE_NSLICE_INFO`, `_NPATCH_`, `_NOBST_` main.f90:4872, 5007, 4794 | 4872, 5007, 4794 | match (description corrected, section 3 item 1) |
| `INITIALIZE_OUTPUT_CLOCKS`, `SET_OUTPUT_CLOCK` func.f90:432, 496 | 432, 496 | match (description corrected, item 2) |
| `INITIALIZE_MESH_DUMPS` 963, `INITIALIZE_GLOBAL_DUMPS` 627 | 963, 627 | match |
| `MPI_INITIALIZATION_CHORES` main.f90:1334 | 1334 | match |
| `DUMP_MESH_OUTPUTS` dump.f90:77; call main.f90:1223 | 77; 1223 | match |
| `DUMP_SLCF` 6746, `DUMP_SLICE_GEOM` 6261, `_DATA` 6399, `DUMP_CFACES_GEOM` 6326, `DUMP_SLCF_VTK` 7323 | all five | match |
| `DUMP_BNDF` 11950, `_VTKHDF` 12177 | both | match |
| `DUMP_SMOKE3D` 5052, `GET_SMOKE3D_QQ` 5117, `_VTKHDF` 5175 | all three | match |
| `DUMP_ISOF` 4860; "calls in 77-287" | 4860; calls at 169-282 (list in section 1) | match |
| `DUMP_PART` 4496, `_VTKHDF` 4622 (loops 4654, 4743) | 4496, 4622; 4654 `DO NMNM=1,NMESHES`, 4743 `MESH_LOOP_HDF: DO NMNM=1,NMESHES` | match |
| `DUMP_PROF` 11362 | 11362 | match (description corrected, item 6) |
| `DUMP_RESTART` 3871 (4017), `READ_RESTART` 4058 (4250) | 3871, 4058; 4017 and 4250 `OTHER_MESH_LOOP: DO NOM=1,NMESHES` | match; the restart file names are `FN_RESTART(NMESHES)` (dump.f90:442) |
| `UPDATE_DEVICES_1/2` 7772, 8338; `DUMP_DEVICES` 11155; `UPDATE_HRR` 11637; `UPDATE_MASS` 11873; `DUMP_HRR` 11830; `DUMP_MASS` 11929 | all seven | match |
| `WRITE_DIAGNOSTICS` 4315 (4386, 4461), `EXCHANGE_DIAGNOSTICS` main.f90:4347 | 4315; 4386, 4461 `DO NM=1,NMESHES`; 4347 | match |
| `DUMP_GEOM` 12538, `DUMP_HVAC` 11535, `DUMP_CONTROLS` 11334 | all three | match |
| `INITIALIZE_BACK_WALL_EXCHANGE` 2410, `ZONE_BOUNDARY_EXCHANGE` 2907, `INITIALIZE_PRESSURE_ZONES` 2679, `INITIALIZE_MESH_EXCHANGE_1` 2132 | all four | match |

Name and structure claims (14), all found: `FDS_SETUP(MODE=3)` end-of-step block (main.f90:120-142: `UPDATE_GLOBAL_OUTPUTS` per rank mesh 127-129, `EXCHANGE_GLOBAL_OUTPUTS` 130, `UPDATE_CONTROLS` 131, `DUMP_GLOBAL_OUTPUTS` 132, `WRITE_STRINGS` 135, `WRITE_DIAGNOSTICS` 138); `DUMP_MESH_OUTPUTS` and `DUMP_RESTART` not called in MODE=3 (comment main.f90:123); the driver calls `fds_setup(3, ...)` after every step (TimeLoop.cpp:1653); `FINE_LEVEL` array (mesh.f90:516, type 511-515) and `POINT_TO_BOX` (mesh.f90:962); `BOX_VIEW` pointer fields RHO, ZZ, TMP, U, V, W, H, D, MU, Q, KRES, RSUM (fds_amrex_hooks.f90:31-37); `BUILD_FINE_BOX` (driver/fds_fine_level.f90:2) and `FDS_FINE_B_SET_VIEW` (driver/fds_box_obj.f90:60); `exact_sum_hierarchy` (driver/ExactSum.H:41); `MESHES(NM)` fields `KRES` (mesh.f90:29), `CELL` (226), `WALL` (298), `N_SLCF` (321), `N_UNIQUE_SLCF` (325); `PROCESS(NM)` (cons.f90:405); ADR-004 D1 to D8 and spike S-A (`adr/ADR-004-smokeview-output.md:34-73, 107`); D-028 and D-054 (`README.md` decision table); FR-072 and FR-074 (`requirements.md`); `_uvw_t<n>_m<NM>.csv`, `_tmp_t..`, `_spec_t..` file names carry the mesh number (dump.f90:239, 249, 258).

Not found as stated: `OutputPlan`, `OUT_MESH`, `OUT_POINT_TO`, `FDS_HOOK_SET_MESH_DUMPS`, `AMR_LEVEL`: they are proposals of the plan (no source defines them), as the plan says.

## 3. Corrections to the plan text (no line number is wrong)
1. `EXCHANGE_NSLICE_INFO`, `EXCHANGE_NPATCH_INFO` and `EXCHANGE_NOBST_INFO` (plan section 2, "for the `.smv`") exchange the data of the **VTK** output (comments main.f90:4792, 4870, 5005). They are called unconditionally at 539 and 542, and 668 only with `WRITE_VTK_GEOM`. They do not feed the `.smv`. They still index by `NMESHES` (e.g. 4874 `USE ... NMESHES`, `MESHES(1)%N_UNIQUE_SLCF`), so they follow the output mesh count only if VTK output is wanted.
2. `SLCF_COUNTER(NM)` and the other per-mesh counters (plan: "one counter per mesh") are allocated `(LOWER_MESH_INDEX:UPPER_MESH_INDEX)` (func.f90:536-538) for the calls with `MESH_SPECIFIC_COUNTER=.TRUE.` (449-474), so they exist only for the compute meshes of the rank. Output meshes numbered outside that range need their own counter range or the clocks must be indexed by a rank-local output index.
3. `IBLANK` is not a per-mesh item: `WRITE_SMOKEVIEW_FILE` writes one global `IBLANK` parameter (dump.f90:1962-1970). Per mesh it writes `GRID`, `PDIM`, `TRNX`, `OBST` (2499, 2505, 2511, 2541). `CELL` is a `MESH_TYPE` field (mesh.f90:226), `IBLANK` is not.
4. `M%N_BNDF` and `M%N_SMOKE3D` are not `MESH_TYPE` fields: `N_BNDF` and `N_SMOKE3D` are global (data.f90:15). `M%N_SLCF` (mesh.f90:321) is a mesh field.
5. The "about 25 allocation lines" are 35 lines in dump.f90:405-451, plus 3 in `WRITE_SMOKEVIEW_FILE` (2056-2060).
6. `DUMP_PROF` (plan: "wall profile files per mesh") writes one file per profile `N` (`FN_PROF(N_PROF)` dump.f90:395, name `_prof_<N>.csv` 400), and each call acts only on profiles whose `PF%MESH` equals the mesh (11393). The key is the compute mesh of the wall, so it does not depend on output meshes. Answer to the plan's question: it does not belong in the output-mesh list.
7. `SANDIA_OUT` is not in `dump.f90`: it is `sandia_out` at turb.f90:2412 (public list turb.f90:24, call dump.f90:282, `USE TURBULENCE` 80).
8. `INITIALIZE_MESH_DUMPS` (plan: "opens and headers the per-mesh files") also creates the `.smv` listing strings (section 1, step B). `INITIALIZE_GLOBAL_DUMPS` (627) writes the global files on rank 0 and is not per mesh.
9. Inventory entries to add: `WRITE_STRINGS` (main.f90:4231, loops over `NMESHES` at 4269); `INITIALIZE_VTKHDF_FILES` (vtkf.f90:1327-1357) which calls `INITIALIZE_VTKHDF_SLCF` (1375; `NMESHES` loops 1395, 1424), `INITIALIZE_VTKHDF_SMOKE3D`, `INITIALIZE_VTKHDF_BNDF_STEP1/2` and `INITIALIZE_VTKHDF_PART` (1516; 1534, 1559); the `NMESHES` loops of the VTK-HDF writers (vtkf.f90:844, 866, 947, 970, 1055, 1079, 1164, 1200) and `WRITE_VTKHDF_GEOM_FILE` (1619; `2*NMESHES` at 1640, 1654-1656, 1677-1678). All are `WITH_HDF5` only. Answer to the plan's question: the VTK-HDF loops belong in the list.
10. Cut-cell files: `DUMP_CFACES_GEOM` (6326) and `WRITE_CFACES` (6310) are reached only with `CC_IBM`; the driver already refuses `CC_IBM` (TimeLoop.cpp:288 `chk(IP[7], "CC_IBM")`), so they stay out of the list. The same check (290) refuses Lagrangian particles for the whole run (`chk(IP[17] > 0, ...)`), so the plan's "particles refused at fine boxes" is stricter today: no particles at all in the driver.
11. Patch numbers: see "What exists" above.

## 4. UP-0011 (pressure RHS dump), re-checked on the current tree
`patch -p1 --dry-run`, files exported with `git show`:
| Base | main.f90 hunks (-1650,6; -1669,6) | pres.f90 hunks (-9,7; -1076,6) |
|---|---|---|
| FDS-AMReX HEAD `2c639da09b` (main.f90 blob `2e1327bbe8`, pres.f90 `58dd2780c9`) | both at 1723 and 1748 (offset +73) | clean |
| `bee11f0329` | both at 1723 and 1748 (offset +73) | clean |
| `36975d765f` | clean | clean |
| master `ce1f659cd4` | both succeed (offset -131) | hunk 2 at offset -4; hunk 1 clean |
Applied for real in the scratch copy: calls at patched main.f90:1729 (`PRESSURE_SOLVER_DUMP(NM,T,DT,.FALSE.)`, after the `PRESSURE_SOLVER_COMPUTE_RHS` loop, original 1723, and before `SELECT CASE(PRES_FLAG)` original 1728) and 1754 (`.TRUE.`, after the `END SELECT` original 1743, before the residual comment original 1745); `PRESSURE_ITERATION_SCHEME` starts at original 1674, `PRESSURE_ITERATION_LOOP` is the `DO` label at 1696; pres.f90 PUBLIC list line 12, routine 1090-1244 (155 lines), `POISSON_FACE_TYPE` 1251, `END MODULE PRES` 1267 (original 1079); the patch adds 188 lines to pres.f90 (6114 to 6302), the note says "about 190"; no added line is longer than 132 characters. The signoff note's citations all hold; the 0009 label is gone (`README.md:10` of the frozen case says UP-0011). Driver-side citations of the note moved: `FDSTL_PDUMP` is now TimeLoop.cpp:836 (was 834 in the earlier check) and the file name is built at 843 (was 841), after an unrelated driver edit; `pb_fds_frozen` is m4_modes.cpp:1-2. The level-0 guard patch 0010 also still applies clean on HEAD (no offset).

## 5. Not verifiable without a build
- That the level-0 arrays of the writers alias the MultiFabs for every field `DUMP_SLCF`, `DUMP_BNDF`, `DUMP_SMOKE3D` and `DUMP_ISOF` read (plan section 3 step A). `BOX_VIEW` (fds_amrex_hooks.f90:31-37) lists 20 fields; fields the writers read beyond these (for example the quantities of `GET_QUANTITY`, wall arrays, `OMESH`) were not traced.
- That `DUMP_MESH_OUTPUTS` with the driver's T, DT and ICYC produces files whose headers equal an FDS run (plan section 5), and the T2-tolerance statement.
- That the `.smv` lists files the driver never writes today (plan section 1): the listing strings are created at set-up (dump.f90:1114-1506 via `INITIALIZE_MESH_DUMPS`), which supports it, but no run was made.
- Compile of the added `FDS_HOOK_SET_MESH_DUMPS` flag and of the MODE=3 block; the `-DWITH_HDF5` build of the `FAKEWRITE` loop.
