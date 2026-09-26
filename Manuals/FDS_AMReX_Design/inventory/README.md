# FDS mesh-data inventory (for the AMReX / AMR port)

**Reference base: FireX**
- Path: `Source`
- Branch `AMReX`, HEAD `36975d765fcead401e14b094a04f910ac42eab8a` (committed 2026-09-24)
- 34 `.f90` files, 180,840 lines, analysed read-only.
  - The task brief said 35 files; the checkout has 34 (`ls *.f90 | wc -l`).

Every `file:line` in this directory refers to that commit, unless it is explicitly marked as the old base.
- The old base is `(local FDS master checkout)/Source` at `ce1f659`. It is used only in `base_delta.md` and in the line remap (`tools/inventory/remap.py`).
- Generators are in `(local project tools directory)/inventory/*.py`. They use the Python stdlib only; each takes a few seconds.
- Classifications that come from heuristics are marked *uncertain*.

> **Status: Section 4 complete (final report).**
> - Checkpoint 2 (accepted): Section 1 (`mesh_fields.csv`), Item A, Item B, FR-016, prototype-P1 scoping, `base_delta.md`.
> - Section 4 (this revision): `routine_field_access.csv`, `cross_mesh_access.csv`, `mesh_exchange.csv`, `exchange_call_sites.csv`, `uniform_grid_assumptions.csv`; the refinement-ratio answer (under the coarse/fine rule); the A-20 upstream slips (under "Snapping and uniform-grid assumptions"); `global_reductions.csv` (every MPI reduction and cross-mesh accumulation, for the fixed-point accumulation requirement); two 40-row Section 4 spot checks.
> - Correction to checkpoint 2: at a level jump the **second** ghost layer is a zero-gradient copy of the first, not neighbour data (wall.f90:351-386). This is fixed below and in `kernel_footprint_mass.md`.

## Files

All CSVs are generated on FireX (the `FDS_SRC` default). Set `FDS_SRC=(local FDS master checkout)/Source FDS_INV_OUT=<dir>` to regenerate them against the old base.

| file | rows | content |
|---|---|---|
| `mesh_fields.csv` | 1,328 | Every MESH_TYPE component (444), every component of the 37 derived types reachable from it (881), and 3 pointer-only rows |
| `allocation_sites.csv` | 2,820 | Every `ALLOCATE` item (2,536) and `MOVE_ALLOC(TO=)` (284), with its owning type resolved (362 sites / 283 distinct MESH_TYPE fields) |
| `pressure_fields_access.csv` | 876 | Item A: per-routine access to the pressure-solver fields, for MESH_TYPE and OMESH_TYPE |
| `point_to_mesh_calls.csv` | 214 | Every live `CALL POINT_TO_MESH` (163 routines) |
| `point_to_mesh_pointers.csv` | 405 | One row per MESH_POINTERS pointer: its target, remap line, and usage counts |
| `module_globals.csv` | 1,551 | Module and PROGRAM-level state outside MESH_TYPE used by mesh kernels, with a `blocks_pure_kernel` flag |
| `mesh_id_dependencies.csv` | 166 | State whose meaning depends on the mesh **number** (breaks if box ids change on regrid) |
| `lbound_extent_dependencies.csv` | 451 | Code that depends on the 0/-1 lower bounds or on IBAR/JBAR/KBAR being the owned extent |
| `obst_wall_cface_indexing.csv` | 22 | Curated: how CELL/WALL/OBST/CFACE/THIN_WALL/BOUNDARY_* storage indexes cells |
| `interface_averaging_sites.csv` | 164 | FR-016: every coarse/fine (EWC `IIO/JJO/KKO_MIN:MAX`) loop, min-corner pick, count normalisation and range store |
| `kernel_footprint_mass.csv` + `kernel_footprint_mass.md` | 337 | Prototype P1: the complete data footprint of `mass.f90`, including ghost widths, and what a per-box shim must provide |
| `base_delta.md` | – | What changed from ce1f659 to FireX, affecting the inventory |
| `routine_field_access.csv` | 12,667 | Section 4: every routine × every MESH_TYPE (11,746 rows) and OMESH_TYPE (921 rows) field it touches, with arguments resolved through the callee's INTENT, plus a runtime/setup `phase` |
| `cross_mesh_access.csv` | 2,109 | Section 4: every reference to another mesh's storage (`MESHES(i)` with i ≠ own, `OMESH(j)`, the neighbour's `OMESH(NM)`), classified by kind and access |
| `mesh_exchange.csv` | 85 | Section 4: every CODE-guarded block of POST_RECEIVES, MESH_EXCHANGE and MESH_CC_EXCHANGE, with its stage (post_recv / pack_send / start_wait / unpack / post_process), fields, buffers, MPI requests and purpose |
| `exchange_call_sites.csv` | 89 | Section 4: every call of MESH_EXCHANGE, POST_RECEIVES, MESH_CC_EXCHANGE and the other 25 `*EXCHANGE*` routines, with code, guard, time-loop stage and neighbouring calls |
| `uniform_grid_assumptions.csv` | 231 | Section 4: single-resolution and uniform-grid assumptions (first-cell widths, ghost spacing copy, snapping, global min cell size, mesh-number sequence, min-corner picks), plus curated rows for A-20, a unit mix, the implicit ratio limit and solver restrictions |
| `global_reductions.csv` | 259 | Section 4: every MPI collective (119), every serial accumulation over all meshes (42) and every per-process pre-sum over own meshes into a global (98), with MPI op, order dependence, what it feeds, a physics/dt/bookkeeping/output flag and the pre-sum locations |
| `upstream_issue_candidates.md` | – | Ready-to-file drafts for firemodels/fds: #1 shape-area adjust decomposition dependence (run-confirmed); #2–#4 thin-OBST collapse slips from A-20 (found by reading) |
| `mesh_ratio_cases.csv` | 143,708 | Section 4 (f): one row per shared mesh face in every FireX Verification/Validation input (MULT-expanded, CATF-inlined): pair, axis, tangential coarse/fine ratios, integer_tiling, class |
| `mesh_ratio_cases_rollup.csv` | 4,073 | One row per input: mesh count, shared faces, worst class, max ratio, class counts, embedded/overlapping pairs |

## How the tables are generated

Regeneration order (caches are keyed by an md5 of `FDS_SRC`: `/tmp/fds_globals_<tag>.json`, `/tmp/fds_routines_<tag>.json`, `/tmp/fds_ptr_usage_<tag>.json`):
1. `alloc_scan`
2. `mesh_fields`
3. `point_to_mesh_calls`
4. `routine_index`
5. `globals_scan`
6. `point_to_mesh_pointers`
7. `field_access --fields ...` (writes pressure_fields_access)
8. `mesh_id_deps`
9. `lbound_extent_scan`
10. `interface_sites`
11. `kernel_footprint mass`
12. `obst_wall_cface_indexing`
13. `phase` (runtime/setup reachability cache `/tmp/fds_phase_<tag>.json`; imported by the Section 4 scripts)
14. `routine_field_access` (runs `field_access` over all MESH_TYPE + OMESH_TYPE components; must run after step 7)
15. `cross_mesh`
16. `mesh_exchange` (writes `mesh_exchange.csv` and `exchange_call_sites.csv`)
17. `uniform_grid_scan` (imports the pick_min rows of `interface_averaging_sites.csv`, so it must run after step 10)
18. `readme_tables` (prints the count and Item A tables in this README)
19. `spotcheck_sec4 [seed]` (prints 40 random Section 4 rows next to the cited source lines)
20. `global_reductions` (raw scan, cache `/tmp/fds_globred_<tag>.json`, ~15 s)
21. `global_reductions_csv` (curated feeds/flags/exclusions; writes `global_reductions.csv` and prints the counts and exclusions)
22. `mesh_ratio_cases` (parses `$FDS_ROOT/{Verification,Validation}/**/*.fds`, default `FDS_ROOT=(repo root)`; writes `mesh_ratio_cases.csv` and `mesh_ratio_cases_rollup.csv`; ~20 s)

Library modules:
- **`fortran_lex.py`** streams logical statements. It:
  - strips comments, respecting strings;
  - upper-cases code;
  - joins `&` continuations and splits on `;`;
  - strips labels and construct names;
  - maps back to physical lines.

  Preprocessor lines are skipped, so **both `#ifdef` branches are scanned**.
- **`fortran_scope.py`** has two parts:
  - `collect_globals()` collects the derived types, module and PROGRAM declarations, and the USE graph.
  - `Walker` resolves names through the frames: associate → local → host → USE (handling ONLY lists and renames).
- **`designator.py`** resolves `A(i)%B(j)%C` chains.
- **`fcommon.py`** holds shared helpers.
- **`routine_index.py`** builds the routine index and call graph.
- **`argres.py`** resolves dummy INTENT for actual arguments.
- **`remap.py`** maps old-base line citations to FireX.
- **`spotcheck_fields.py [seed] [n]`** prints random `mesh_fields` rows next to the cited source lines.

## `mesh_fields.csv` (Section 1)

Columns:
- `field`
- `parent_type`
- `fortran_type`
- `rank`
- `declared_bounds`: deferred, plus the bounds used at allocation
- `decl_file`, `decl_line`
- `alloc_sites`: `file:line[bounds]`
- `category`
- `notes`: curated facts with file:line, the `doc:` Doxygen comment, flags, and `#ifdef` alternatives

Categories: `grid_geometry`, `cell_centred_solution`, `face_centred_flux`, `pressure_solver`, `radiation`, `scratch_work`, `ghost_boundary`, `mpi_exchange`, `particle`, `wall_obst_cutcell`, `other`.
- MESH_TYPE rows are categorised from a curated list.
- Nested rows default to their parent type's category.

| category | MESH_TYPE (444) | nested (881 in 37 types) | pointer-only (3) | total (1,328) |
|---|---|---|---|---|
| wall_obst_cutcell | 102 | 548 |  | 650 |
| mpi_exchange | 4 | 165 |  | 169 |
| other | 99 | 59 | 3 | 161 |
| pressure_solver | 50 | 29 |  | 79 |
| grid_geometry | 71 |  |  | 71 |
| ghost_boundary | 8 | 42 |  | 50 |
| particle | 12 | 35 |  | 47 |
| scratch_work | 35 |  |  | 35 |
| face_centred_flux | 27 |  |  | 27 |
| cell_centred_solution | 25 |  |  | 25 |
| radiation | 11 | 3 |  | 14 |

Largest nested types: OMESH_TYPE (151), VENTS_TYPE (76), CC_CUTFACE_TYPE (65), CC_CUTCELL_TYPE (58), OBSTRUCTION_TYPE (57), BOUNDARY_PROP1_TYPE (47), BOUNDARY_ONE_D_TYPE (39), LAGRANGIAN_PARTICLE_TYPE (35).

Differences from ce1f659:
- +13 MESH_TYPE fields for output and VTK (mesh.f90:159, 323-329).
- +1 SLICE_TYPE field (`SLCF_NAME`, type.f90:1603).
- +4 ZONE_SOLVE_TYPE RS fields (type.f90:1816-1818).

See `base_delta.md`.

## Notable observations (Section 1, FireX lines)

1. **Access model is a global pointer swap.**
   - `MESHES` is at mesh.f90:356 (allocated at read.f90:579).
   - The `MESH_POINTERS` module (mesh.f90:361) declares the module pointers. `POINT_TO_MESH(NM)` (mesh.f90:505-916) re-points 402 of them.
   - There are 214 live calls in 163 routines. Only one is in main.f90 (main.f90:734); main.f90 mostly uses `M=>MESHES(NM)`.
   - Two calls hard-code mesh 1: `POINT_TO_MESH(1)` at vtkf.f90:1875 and 2482.
   - Several routines assume the caller has already done the swap: comments at ccib.f90:370, 510, 666 and 2602.
   - So the "current mesh" is hidden global state, which blocks re-entrant or tile-parallel kernels.
2. **Dead pointers and fields.**
   - `MTR`, `MSR` and `WEM` are declared in MESH_POINTERS (mesh.f90:370) but are not MESH_TYPE components and are never pointed or used.
   - 42 MESH_TYPE components have no pointer alias.
   - `SAVE2` (mesh.f90:88), `PWORK4` (mesh.f90:121) and `PHI_W` (mesh.f90:352) are never allocated; they are only aliased, at mesh.f90:574, 648 and 904.
   - The ChkMemErr label `'PHI_W'` at vege.f90:222 actually belongs to `PHI_WS`.
   - OMESH `IFEP_R_1`/`IFEP_R_2` (type.f90:1071) are never allocated. They are read only inside `DO IFEP=1,M2%NFEP_R(1|2)` loops (e.g. ccib.f90:4228, 4326). `NFEP_R` defaults to 0 (type.f90:1058), and only `NFEP_R(3:5)` are ever incremented (ccib.f90:3827, 3864, 3908). So these are dead paths (*uncertain*: static analysis only).
3. **Ghost widths differ by field.**
   - Two ghost layers `(-1:IBP1+1)`: `TMP`, `RHO`, `RHOS`, `ZZ`, `ZZS`, `WORK_PAD`, "for flux limiters" (init.f90:522-529).
   - Face velocities `U/US`, `V/VS`, `W/WS` have an extra ghost layer only in their own normal direction (init.f90:533-535).
   - Everything else is `0:IBP1`.
   - Face fluxes `FVX/FVY/FVZ` use cell-shaped `0:IBP1` storage (init.f90:543-545).
   - `kernel_footprint_mass.md` confirms that mass.f90 really indexes the second ghost layer.
4. **Poisson work arrays are 1-based with no ghosts.**
   - `PRHS(ITRN,JTRN,KTRN)` is at init.f90:2445.
   - The BC arrays `BXS..BZF` are at init.f90:2448-2453.
   - Both are re-allocated with the same shape by the ULMAT setup (pres.f90:1280-1303).
5. **Background state is a 1-D z-column per mesh.** `P_0(0:KBP1)` (init.f90:449), `RHO_0` and `TMP_0` share that shape, and `PBAR` is (K, zone).
6. **Neighbour data is replicated.**
   - `OMESH(1:NMESHES)` is allocated in every mesh (main.f90:2077).
   - OMESH holds copies of neighbour fields over the receive window (main.f90:2203-2264).
   - Off-process `MESHES(NOM)` is partly allocated:
     - `PRESSURE_ZONE` at main.f90:1355, filled by MPI_BCAST at main.f90:2853;
     - `WALL` at main.f90:2497;
     - `CELL_INDEX` at func.f90:5357.
7. **Halo fill happens in two stages.**
   - First, MESH_EXCHANGE (main.f90:3117-3975) unpacks into **OMESH**. That is true even on the same process: `RHOP2 => M2%RHO`, direct copy at main.f90:3339-3341.
   - Then per-mesh routines (ASSIGN_GHOST_VALUE wall.f90:282-388, COPY_H_OMESH_TO_MESH pres.f90:4054-4208, NO_FLUX velo.f90:1394) average or inject OMESH values into the mesh's own ghosts over the EWC ranges (see FR-016).
   - `TMP` is not exchanged.
8. **Layout.**
   - `WALL(0:N_WALL_CELLS_DIM)` lists external walls first (main.f90:2104-2108, init.f90:74-107).
   - `WALL` grows by re-allocation (func.f90:4011).
   - `CELL_INDEX` is a sparse IJK → cell map (read.f90:11622); see the obst/wall/cface section.

## Item A (Pressure Solver Lead)

`pressure_fields_access.csv` columns:
- `field`, `parent_type`, `file`, `routine` (module::host::routine)
- `access`
- `via`: `ptr` for a MESH_POINTERS pointer, the base name (`M`, `OM`, `M2`, `MESHES`), or `alias:P` for a local pointer chain such as `DPVOL => DP => D` (ccib.f90:8379)
- `n_refs`, `first_line`, `lines`
- `callee`, `arg_evidence`

`access` values:
- `read`
- `write`
- `ptr_target`
- `arg_read` / `arg_readwrite`: the actual argument's dummy INTENT is resolved by `argres.py`
- `alloc` / `dealloc`
- `inquiry`

Access counts on FireX: read 410, write 191, ptr_target 153, arg_read 53, alloc 38, arg_readwrite 13, inquiry 11, dealloc 7.

MESH_TYPE summary (generated by `readme_tables.py`; the OMESH copies are in the CSV):

| field | decl | ALLOCATE site [bounds] | writer routines (first write site) | # reader routines | # routines passing it as an actual argument (of which readwrite) |
|---|---|---|---|---|---|
| H | mesh.f90:27 | init.f90:556[0:IBP1,0:JBP1,0:KBP1] | CC_H_INTERP (ccib.f90:7102), READ_RESTART (dump.f90:4096), INITIALIZE_MESH_VARIABLES_1 (init.f90:790), COPY_H_OMESH_TO_MESH (pres.f90:4130), GLMAT_SOLVER (pres.f90:3577), PRESSURE_SOLVER_FFT (pres.f90:397), ULMAT_SOLVE_ZONE (pres.f90:1892), NS_ANALYTICAL_SOLUTION (turb.f90:135), SHUNN_MMS_3 (turb.f90:3402), NO_FLUX (velo.f90:1399) | 30 | 4 (3) |
| HS | mesh.f90:28 | init.f90:557[0:IBP1,0:JBP1,0:KBP1] | CC_H_INTERP (ccib.f90:7102), COPY_CC_UNKH_TO_HS (ccib.f90:22944), READ_RESTART (dump.f90:4101), INITIALIZE_MESH_VARIABLES_1 (init.f90:791), COPY_CCVAR_IN_HS (pres.f90:4350), COPY_HS_IN_CCVAR (pres.f90:4268), COPY_H_OMESH_TO_MESH (pres.f90:4196), GLMAT_SOLVER (pres.f90:3577), PRESSURE_SOLVER_FFT (pres.f90:397), RESTORE_HS_AFTER_SETUP (pres.f90:3281), ULMAT_SOLVE_ZONE (pres.f90:1892), NS_ANALYTICAL_SOLUTION (turb.f90:136), SHUNN_MMS_3 (turb.f90:3403), NO_FLUX (velo.f90:1399) | 23 | 5 (3) |
| FVX | mesh.f90:30 | init.f90:543[0:IBP1,0:JBP1,0:KBP1] | CC_BAROCLINIC_CORRECTION (ccib.f90:1677), CC_MATCH_VELOCITY_FLUX (ccib.f90:5216), CC_VELOCITY_FLUX (ccib.f90:15060), GET_LINKED_FV (ccib.f90:15315), ROTATED_CUBE_VELOCITY_FLUX (ccib.f90:6162), INITIALIZE_MESH_VARIABLES_1 (init.f90:786), PARTICLE_MOMENTUM_TRANSFER (part.f90:4767), BAROCLINIC_CORRECTION (velo.f90:3234), CORIOLIS_FORCE (velo.f90:997), DIRECT_FORCE (velo.f90:900), MATCH_VELOCITY_FLUX (velo.f90:2951), MMS_VELOCITY_FLUX (velo.f90:1043), NO_FLUX (velo.f90:1419), PATCH_VELOCITY_FLUX (velo.f90:1124), VELOCITY_FLUX (velo.f90:711), VELOCITY_FLUX_CYLINDRICAL (velo.f90:1300) | 25 | 0 (0) |
| FVY | mesh.f90:31 | init.f90:544[0:IBP1,0:JBP1,0:KBP1] | CC_BAROCLINIC_CORRECTION (ccib.f90:1678), CC_MATCH_VELOCITY_FLUX (ccib.f90:5268), CC_VELOCITY_FLUX (ccib.f90:15061), GET_LINKED_FV (ccib.f90:15316), INITIALIZE_MESH_VARIABLES_1 (init.f90:787), PARTICLE_MOMENTUM_TRANSFER (part.f90:4768), BAROCLINIC_CORRECTION (velo.f90:3235), CORIOLIS_FORCE (velo.f90:1011), DIRECT_FORCE (velo.f90:921), MATCH_VELOCITY_FLUX (velo.f90:2975), NO_FLUX (velo.f90:1436), PATCH_VELOCITY_FLUX (velo.f90:1158), VELOCITY_FLUX (velo.f90:769) | 23 | 0 (0) |
| FVZ | mesh.f90:32 | init.f90:545[0:IBP1,0:JBP1,0:KBP1] | CC_BAROCLINIC_CORRECTION (ccib.f90:1679), CC_MATCH_VELOCITY_FLUX (ccib.f90:5321), CC_VELOCITY_FLUX (ccib.f90:15062), GET_LINKED_FV (ccib.f90:15317), ROTATED_CUBE_VELOCITY_FLUX (ccib.f90:6208), INITIALIZE_MESH_VARIABLES_1 (init.f90:788), PARTICLE_MOMENTUM_TRANSFER (part.f90:4769), BAROCLINIC_CORRECTION (velo.f90:3236), CORIOLIS_FORCE (velo.f90:1025), DIRECT_FORCE (velo.f90:942), MATCH_VELOCITY_FLUX (velo.f90:2999), MMS_VELOCITY_FLUX (velo.f90:1051), NO_FLUX (velo.f90:1453), PATCH_VELOCITY_FLUX (velo.f90:1191), VELOCITY_FLUX (velo.f90:827), VELOCITY_FLUX_CYLINDRICAL (velo.f90:1335) | 25 | 0 (0) |
| D | mesh.f90:25 | init.f90:560[0:IBP1,0:JBP1,0:KBP1] | CCENTHALPY_ADVECTION (ccib.f90:9542), CC_CONDUCTION_HEAT_FLUX (ccib.f90:10670), CC_DIFFUSIVE_HEAT_FLUXES (ccib.f90:10249), CC_DIVERGENCE_PART_1 (ccib.f90:8387), DIVERGENCE_PART_1 (divg.f90:84), DIVERGENCE_PART_2 (divg.f90:1540), READ_RESTART (dump.f90:4095), INITIALIZE_MESH_VARIABLES_1 (init.f90:793), COMPRESSION_WAVE (turb.f90:214), VISCOSITY_BC (velo.f90:548) | 16 | 2 (1) |
| DS | mesh.f90:26 | init.f90:561[0:IBP1,0:JBP1,0:KBP1] | DIVERGENCE_PART_1 (divg.f90:84), DIVERGENCE_PART_2 (divg.f90:1540), READ_RESTART (dump.f90:4100), INITIALIZE_MESH_VARIABLES_1 (init.f90:794), COMPRESSION_WAVE (turb.f90:215), VISCOSITY_BC (velo.f90:546) | 12 | 1 (1) |
| DDDT | mesh.f90:24 | init.f90:559[0:IBP1,0:JBP1,0:KBP1] | DIVERGENCE_PART_2 (divg.f90:1619), INITIALIZE_MESH_VARIABLES_1 (init.f90:792) | 4 | 0 (0) |
| KRES | mesh.f90:29 | init.f90:558[0:IBP1,0:JBP1,0:KBP1] | CC_COMPUTE_KRES (ccib.f90:3728), INITIALIZE_MESH_VARIABLES_1 (init.f90:789), NS_ANALYTICAL_SOLUTION (turb.f90:148), COMPUTE_VISCOSITY (velo.f90:291), VISCOSITY_BC (velo.f90:544) | 23 | 0 (0) |
| RHO | mesh.f90:39 | init.f90:525[-1:IBP1+1,-1:JBP1+1,-1:KBP1+1] | CC_CHECK_MASS_DENSITY (ccib.f90:12574), CC_DENSITY (ccib.f90:11857), GET_RHOZZ_CC_3D (ccib.f90:13309), INIT_CUTCELL_DATA (ccib.f90:7788), READ_RESTART (dump.f90:4104), CREATE_OR_REMOVE_OBST (init.f90:4882), INITIALIZE_MESH_VARIABLES_1 (init.f90:525), SPEC_INIT (init.f90:5518), TMP_INIT (init.f90:5448), CHECK_MASS_DENSITY (mass.f90:854), DENSITY (mass.f90:691), SAAD_MMS_1 (turb.f90:3362), SHUNN_MMS_3 (turb.f90:3401), ASSIGN_GHOST_VALUE (wall.f90:343), CALCULATE_RHO_F (wall.f90:1627), SURFACE_HEAT_TRANSFER (wall.f90:717) | 68 | 7 (0) |
| RHOS | mesh.f90:40 | init.f90:526[-1:IBP1+1,-1:JBP1+1,-1:KBP1+1] | CC_DENSITY (ccib.f90:11796), GET_RHOZZ_CC_3D (ccib.f90:13250), INIT_CUTCELL_DATA (ccib.f90:7789), CREATE_OR_REMOVE_OBST (init.f90:4881), INITIALIZE_MESH_VARIABLES_1 (init.f90:526), SPEC_INIT (init.f90:5519), TMP_INIT (init.f90:5449), CHECK_MASS_DENSITY (mass.f90:854), DENSITY (mass.f90:509), ASSIGN_GHOST_VALUE (wall.f90:343), CALCULATE_RHO_F (wall.f90:1627), SURFACE_HEAT_TRANSFER (wall.f90:717) | 39 | 6 (0) |
| PRHS | mesh.f90:89 | init.f90:2445[ITRN,JTRN,KTRN];  pres.f90:1287[ITRN,JTRN,KTRN] | INITIALIZE_POISSON_SOLVER (init.f90:2458), PRESSURE_SOLVER_COMPUTE_RHS (pres.f90:243), TUNNEL_POISSON_SOLVER (pres.f90:534), ULMAT_SOLVER_SETUP (pres.f90:1297) | 4 | 1 (1) |
| BXS | mesh.f90:90 | init.f90:2448[JDIM,KBP1];  pres.f90:1290[JDIM,KBP1] | GET_PRES_CFACE_BCS (ccib.f90:2359), INITIALIZE_POISSON_SOLVER (init.f90:2459), PRESSURE_SOLVER_COMPUTE_RHS (pres.f90:81), PRESSURE_SOLVER_FFT (pres.f90:441), TUNNEL_POISSON_SOLVER (pres.f90:604), ULMAT_SOLVER_SETUP (pres.f90:1298), ULMAT_SOLVE_ZONE (pres.f90:1926) | 5 | 1 (0) |
| BXF | mesh.f90:90 | init.f90:2449[JDIM,KBP1];  pres.f90:1291[JDIM,KBP1] | GET_PRES_CFACE_BCS (ccib.f90:2389), INITIALIZE_POISSON_SOLVER (init.f90:2460), PRESSURE_SOLVER_COMPUTE_RHS (pres.f90:83), PRESSURE_SOLVER_FFT (pres.f90:442), TUNNEL_POISSON_SOLVER (pres.f90:605), ULMAT_SOLVER_SETUP (pres.f90:1299), ULMAT_SOLVE_ZONE (pres.f90:1927) | 5 | 1 (0) |
| BYS | mesh.f90:90 | init.f90:2450[IBP1,KBP1];  pres.f90:1292[IBP1,KBP1] | GET_PRES_CFACE_BCS (ccib.f90:2419), INITIALIZE_POISSON_SOLVER (init.f90:2461), PRESSURE_SOLVER_COMPUTE_RHS (pres.f90:85), ULMAT_SOLVER_SETUP (pres.f90:1300) | 4 | 1 (0) |
| BYF | mesh.f90:90 | init.f90:2451[IBP1,KBP1];  pres.f90:1293[IBP1,KBP1] | GET_PRES_CFACE_BCS (ccib.f90:2449), INITIALIZE_POISSON_SOLVER (init.f90:2462), PRESSURE_SOLVER_COMPUTE_RHS (pres.f90:87), ULMAT_SOLVER_SETUP (pres.f90:1301) | 4 | 1 (0) |
| BZS | mesh.f90:90 | init.f90:2452[IBP1,JDIM];  pres.f90:1294[IBP1,JDIM] | GET_PRES_CFACE_BCS (ccib.f90:2479), INITIALIZE_POISSON_SOLVER (init.f90:2463), PRESSURE_SOLVER_COMPUTE_RHS (pres.f90:89), ULMAT_SOLVER_SETUP (pres.f90:1302) | 4 | 1 (0) |
| BZF | mesh.f90:90 | init.f90:2453[IBP1,JDIM];  pres.f90:1295[IBP1,JDIM] | GET_PRES_CFACE_BCS (ccib.f90:2509), INITIALIZE_POISSON_SOLVER (init.f90:2464), PRESSURE_SOLVER_COMPUTE_RHS (pres.f90:91), ULMAT_SOLVER_SETUP (pres.f90:1303) | 4 | 1 (0) |
| P_0 | mesh.f90:94 | init.f90:449[0:M%KBP1] | INITIALIZE_ATMOSPHERE (init.f90:475) | 13 | 0 (0) |
| PRESSURE_ZONE | mesh.f90:107 | init.f90:586[0:IBP1,0:JBP1,0:KBP1];  main.f90:1355[0:M4%IBP1,0:M4%JBP1,0:M4%KBP1] | ASSIGN_PRESSURE_ZONE (func.f90:5586), INITIALIZE_MESH_VARIABLES_1 (init.f90:586), INITIALIZE_PRESSURE_ZONES (main.f90:2623), MPI_INITIALIZATION_CHORES (main.f90:1356), ZONE_BOUNDARY_EXCHANGE (main.f90:2879) | 41 | 1 (1) |

**Exchange codes that carry these fields** (MESH_EXCHANGE main.f90:3117-3975; POST_RECEIVES main.f90:2945-3111):
- **CODE 5** (PKG7, REQ7) carries `FVX/FVY/FVZ` at the interface face and `H` (predictor) or `HS` (corrector) for the two cells either side.
  - Pack at main.f90:3352-3395; unpack into OMESH at main.f90:3740-3770.
- **Inside `PRESSURE_ITERATION_SCHEME`** (main.f90:1601-1745), CODE 5 is called at:
  1. main.f90:1617 (CC_IBM);
  2. main.f90:1636 (after BAROCLINIC_CORRECTION);
  3. main.f90:1663 (after GLMAT_SOLVER at main.f90:1662, followed by COPY_H_OMESH_TO_MESH at main.f90:1664);
  4. main.f90:1687 (before COMPUTE_VELOCITY_ERROR).

  The FFT solver runs per mesh at main.f90:1659 and ULMAT at main.f90:1668. Errors are combined by MPI_ALLGATHERV at main.f90:1702.
- **CODE 3** sends `HS, US, VS, WS` (main.f90:909) and **CODE 6** sends `H, U, V, W` (main.f90:1105).
  - Shared pack at main.f90:3423-3456; unpack at main.f90:3791-3811.
- **CODE 1/4** (PKG1) carry `D/DS`, `KRES`, `RHOS/RHO`, `ZZS/ZZ`, `MU` and `Q`, two cells deep.
  - Pack at main.f90:3301-3345; unpack at main.f90:3704-3733.
- **Never exchanged:** `DDDT`, `PRHS`, `BXS..BZF` and `P_0`.
- **PRESSURE_ZONE** is sent by MPI_BCAST at main.f90:2853.
- **Global reductions:** `DSUM/PSUM/USUM` are MPI_ALLREDUCEd in EXCHANGE_DIVERGENCE_INFO (main.f90:2038-2040).
- **FireX GPU/RS additions:**
  - GLMAT gathers to an RS master, solves, then scatters (pres.f90:3414-3491).
  - HYPRE device migration (pres.f90:1753-1770, 3011-3044, 3455-3472).
  - DEFINE_RS_COMM_INFO (main.f90:5096-5146, called at main.f90:153).

  See `base_delta.md`.

## Item B

### `module_globals.csv`

Columns: `module`, `variable`, `decl_file`, `decl_line`, `kind`, `written_by`, `read_by_sample`, `blocks_pure_kernel`, `notes`, `fortran_type`, `dims`, `n_writers`, `n_runtime_writers`, `n_kernel_users`, `mesh_id_indexed`.

Definitions:
- **Mesh kernel:** a routine that touches MESH_TYPE data directly (576 routines).
- **Module global:** any name that resolves to a MODULE, or to the PROGRAM FDS specification part (host-associated into main.f90 contained routines). MESH_POINTERS names are excluded.
- **Runtime:** reachable in the call graph from calls inside `MAIN_LOOP` (main.f90:695-1213). That gives 695 runtime routines, against 1,007 setup routines.
- **`blocks_pure_kernel`:**
  - `yes` = a POINTER global, or written by any runtime routine;
  - `maybe` = passed as an actual argument whose INTENT could not be resolved, in a runtime routine;
  - `no` = PARAMETER constants, or written only at setup. Setup-only globals still need a device/const mirror.

| kind | yes | maybe | no |
|---|---|---|---|
| runtime_scalar (751) | 112 | 3 | 636 |
| runtime_array (395) | 145 | | 250 |
| pointer (69) | 69 | | |
| constant (336) | | | 336 |

Examples:
- `DT`/`T` (PROGRAM FDS, main.f90:65; written in MAIN_LOOP at main.f90:715 and 933).
- `PREDICTOR`/`CORRECTOR` (cons.f90:192-193; written at main.f90:748-749).
- `ICYC` (cons.f90:587).
- `T_USED` (cons.f90:719; 82 runtime writers).
- `DSUM` (cons.f90:648; written at divg.f90:739 and main.f90:2038).
- `REQ1`/`N_REQ1` (main.f90:85/87; written in MESH_EXCHANGE main.f90:3634, 3195).
- `MESHES` itself (mesh.f90:356; 17 runtime writers).

**Flag (a): mesh-number dependencies (`mesh_id_dependencies.csv`)**

| kind | rows | what it means |
|---|---|---|
| module_array_by_mesh | 68 | Module arrays dimensioned by NMESHES, e.g. `NUNKZ_LOC` (ccib.f90:66, allocated at ccib.f90:22969) |
| component_array_by_mesh | 10 | Type components dimensioned by NMESHES, e.g. INITIALIZATION_TYPE `PARTICLE_INSERT_CLOCK` (type.f90:1685, read.f90:13343) |
| component_stores_mesh_id | 19 | Components holding a mesh number, e.g. DEVICE_TYPE `MESH` (devc.f90:87), BOUNDARY_ONE_D `BACK_MESH` (type.f90:226) |
| hardcoded_mesh_1 | 24 | 100 `MESHES(1)%` occurrences use mesh 1 as a global holder for slice/VTK output state: `N_UNIQUE_SLCF` (mesh.f90:325; 19 sites, e.g. main.f90:4814), `N_SLCF_O` (mesh.f90:323). Also `POINT_TO_MESH(1)` at vtkf.f90:1875 and 2482 |
| local_array_by_mesh | 45 | Local arrays sized by NMESHES, e.g. `REQ0` (ccib.f90:18942) |

Also:
- `PROCESS(NM)` (cons.f90:405) and `LOWER/UPPER_MESH_INDEX` (cons.f90:402) fix the mesh-to-rank map at READ_MESH (read.f90:712-716).
- `OMESH(NOM)` and `EWC%NOM` are indexed by mesh number.

**Flag (b): lower bounds and IBAR-as-extent (`lbound_extent_dependencies.csv`)**

Detail rows:

| category | rows |
|---|---|
| pointer_bounds_remap | 58 |
| deferred_shape_dummy_keeps_bounds | 48 |
| local_array_mesh_extent | 20 |
| first_element_actual | 12 |
| assumed_shape_dummy_explicit_lb | 11 |
| explicit_shape_dummy_mesh_extent | 5 |

Aggregated per routine:
- loop_owned_extent: 738 statements in 101 routines;
- loop_ghost_extent: 414 in 65;
- ghost_index_literal: 666 in 113;
- coord_to_index: 61 in 18.

Implication: MESH_POINTERS must be re-pointed at FAB data with FDS lower bounds (e.g. `RHO(-1:,-1:,-1:) => fab`), and IBAR/JBAR/KBAR must equal the box extent.

### `point_to_mesh_pointers.csv`

- 402 remaps in POINT_TO_MESH (mesh.f90:513-914).
- `MTR`, `MSR` and `WEM` are unused.
- 42 MESH_TYPE components have no alias.
- Columns give each pointer's target, its remap line and its usage counts (from `globals_scan`).

### `obst_wall_cface_indexing.csv` (22 curated rows)

- **CELL_INDEX is sparse.**
  - It is filled by READ_OBST (read.f90:11622-11688).
  - `IC>0` only for the 2-cell shell (0, 1, IBAR, IBP1) and for cells around OBSTs. Every other cell maps to `CELL(0)`, the shared plain-gas cell.
  - The exterior ring is SOLID by BLOCK_CELL (read.f90:11697-11702); OBST interiors are set at read.f90:11708. The ghost is reopened at interpolated boundaries (init.f90:3281).
  - CELL_TYPE is at type.f90:2175-2187.
- **WALL.**
  - External walls come first: `N_EXTERNAL_WALL_CELLS` formula at read.f90:700; numbering order ±1, ±2, ±3 at init.f90:74-107.
  - Internal walls are appended at init.f90:195-330.
  - `WALL(0)` is a null wall (init.f90:48, 51-52).
  - `CELL%WALL_INDEX(-IOR)` is set at init.f90:3298.
- **Storage slots.**
  - `WALL%BC/OD/TD/B1/B2/BR_INDEX` (type.f90:436-441) are slot indices from ALLOCATE_STORAGE (func.f90:3942-4384; WALL branch 4003-4034; slot logic 4110-4134).
  - Slots come from `NEXT_AVAILABLE_*_SLOT` / `*_OCCUPANCY` (mesh.f90:311-318), so in general `BC_INDEX ≠ IW`.
- **Other types.**
  - BOUNDARY_COORD: type.f90:189-212.
  - EXTERNAL_WALL: type.f90:462-482.
  - THIN_WALL: type.f90:487-500; built at init.f90:3417, 4373.
  - OBSTRUCTION `I1..K2`: type.f90:1116-1121, set at read.f90:11157-11169.
  - CFACE: type.f90:1361-1382; built at geom.f90:12147.
  - `CUT_FACE%CFACE_INDEX`: type.f90:1305.
  - CCVAR/FCVAR use `NGUARD=5` (geom.f90:33).
  - `FIND_WALL_INDEX` (func.f90:5449-5490) does a cross-mesh search (HT3D setup, init.f90:1951).
- **AMR consequence:** every AMReX box needs its own sparse CELL/WALL build (READ_OBST + INIT_WALL_CELL equivalents), with external walls on every box face.

## FR-016: coarse/fine interfaces (`interface_averaging_sites.csv`, 164 rows)

Columns: `file`, `line`, `end_line`, `routine`, `phase`, `kind`, `operation`, `direction`, `range_source`, `loop_vars`, `fields_read`, `fields_written`, `accumulators`, `weights`, `notes`.

Row breakdown:
- **Kinds:**
  - `range_loop`: 119 in total — average 75, scatter_or_copy 14, area_or_count_sum 10, other 18, sum 1, pick_min_or_max 1;
  - `pick_min_corner`: 13, plus 9 inline;
  - `count_normalisation`: 13;
  - `range_extent`: 6 (the IMIN/IMAX… receive bounding box, main.f90:2118-2128);
  - `range_span`: 3 (AREA_RATIO, main.f90:2132-2136);
  - `range_store`: 1.

  The Section 4 spot check found that the 6 range_extent and 3 range_span rows had been classified as `pick_min_corner_inline` in checkpoint 2. They use the full MIN..MAX range, so the classifier was fixed. The file grew from 161 to 164 rows because the three `_MAX+1` extent lines are now captured.
- **Phases:** runtime 70, runtime+setup 61, setup 33.
- **Files:** ccib 59, pres 28, velo 27, geom 14, main 12, part 9, init 5, wall 5, fire 2, vege 2, turb 1.

**The coarse/fine rule (verified on FireX):**
1. **Range definition.** INIT_WALL_CELL (init.f90:3136-3180) places two probe points at ±0.475 of this cell's face width, offset by `MESH_SEPARATION_DISTANCE` normal to the face.
   - SEARCH_OTHER_MESHES returns the neighbour cell for each probe, and MIN/MAX of the two gives `EWC%IIO/JJO/KKO_MIN:MAX`.
   - If the neighbour is finer, the range spans several cells; if it is equal or coarser, a single cell (FireX comment at pres.f90:5419-5422).
   - The range is stored at init.f90:3317-3326.
   - `EWC%NIC` is the per-wall receive count (main.f90:2111).
2. **Alignment.** ERROR 431 (init.f90:3184-3226) requires fine cells to tile the coarse face exactly, and forbids mixed refinement across the two tangential directions.
3. **`AREA_RATIO`** (main.f90:2131-2138) = own face area / area spanned by the range. It is ≈1 when the neighbour is equal or finer and <1 when it is coarser, and it is tested as `>0.9` meaning "this side is not the fine side".
4. **Runtime: one code path serves both directions.**
   - ASSIGN_GHOST_VALUE (wall.f90:319-339) uses `ARO = MIN(1, A_other/A_this)`. The coarse side gets the area-weighted mean of the fine cells; the fine side gets injection of the single coarse cell. The first ghost layer is always filled from the neighbour. The **second** ghost layer gets neighbour data only when the tangential face areas match (|AREA_RATIO−1|<0.01, wall.f90:356) and none of the four cells spanning the boundary is solid (wall.f90:361-371). Otherwise, which includes **every level jump**, it is a zero-gradient copy of the first ghost (wall.f90:382-386). The limiter stencil is therefore first-order at coarse/fine interfaces. *(Correction: checkpoint 2 said both layers were filled from the neighbour.)*
   - MATCH_VELOCITY (velo.f90:2724-2803) saves `UVW_SAVE` = own U, then sets U(0) = ½(own + area-mean of the other side).
   - Arithmetic means divided by `N_INT_CELLS`:
     - NO_FLUX (velo.f90:1391-1399);
     - VISCOSITY_BC (velo.f90:526-541);
     - VELOCITY_BC (velo.f90:1873-1888);
     - COMBUSTION_BC (fire.f90:1909-1917);
     - vege.f90:678-705;
     - COMPUTE_VELOCITY_ERROR (pres.f90:924-1038);
     - COPY_H_OMESH_TO_MESH (pres.f90:4109-4126, 4175-4192).
   - Min-corner pick (the first fine cell only):
     - wall.f90:357 (2nd-order ghost, only when |AREA_RATIO−1|<0.01, wall.f90:356), wall.f90:839, 1371;
     - turb.f90:2284;
     - init.f90:1092, 1413, 4924, 5164;
     - pres.f90:124-139 (DX_OTHER);
     - ccib.f90:22887-22889;
     - main.f90:2890.

     These are candidate accuracy defects at refinement interfaces.
   - geom.f90:8988 runs only when the range is larger than one cell (coarse reads fine).
5. **No reflux.** DENSITY restores the pre-average `UVW_SAVE` at interpolated faces (mass.f90:424-436, 595-608). The coarse and fine face mass fluxes are therefore not matched, and conservation across a level jump is not guaranteed.

### Refinement ratio: what FDS itself allows (Section 4, question a)

Question: the team capped AMR refinement ratios at 4 because AMReX MLMG's InterpBndryData asserts ratio ≤ 4. Does FDS limit the cell-size ratio between adjacent meshes?

**Answer:** FDS enforces **no explicit maximum**. It requires **integer nesting per coarse face**, but **not a power of two**, and **not the same ratio in each direction**. The 4:1 figure is User Guide advice only, so the team's cap of 4 is stricter than FDS and matches that advice.

1. **No maximum in the code.**
   - read.f90 has no ratio check. Its MESH/TRN errors cover line syntax, 2-D/cylindrical consistency, MULT_ID, MPI_PROCESS assignment, degenerate XB and non-monotonic transformations, never neighbour sizes: ERROR(110)-(122) at read.f90:541-816 and (124)-(126) at read.f90:1182-1376.
   - The only mesh-to-mesh size test is ERROR(431) in INIT_WALL_CELL (init.f90:3184-3226), and it has no ratio bound.
   - The only "4" is User Guide advice: "Mesh refinement ratios of more than 4:1 to should be avoided if possible" (sic; `Manuals/FDS_User_Guide/FDS_User_Guide.tex:1127`).
   - The UG's other 4:1 (FDS_User_Guide.tex:9097) is about **cell aspect ratio** and the CFL norm, not refinement.
2. **Implicit hard limit: ratio < 40** (*derived from the code, not tested by a run*).
   - For each exterior face, INIT_WALL_CELL places two probes at ±0.475 of the **own** tangential face width (init.f90:3161-3163), i.e. 0.025·Δ_own inside each edge, and offsets them by MESH_SEPARATION_DISTANCE normal to the face (init.f90:3164ff).
   - The neighbour cells found by the probes give IIO_MIN:IIO_MAX.
   - If the neighbour's first/last fine cell is narrower than 0.025·Δ_coarse (ratio ≥ 40), the probes skip it. The span X(IIO_MAX)−X(IIO_MIN−1) is then shorter than the coarse face, and the alignment test raises ERROR(431) (init.f90:3194-3208).
   - At exactly 40 the probe lies on a fine-cell face, and which cell it lands in depends on rounding.
3. **Integer nesting per coarse face, not a power of two.**
   - When the neighbour is finer or equal in a tangential direction (`MM%DX(IIO_MIN)<=M%DX(I)`), the fine cells found must tile the coarse face exactly: |span − Δ_coarse| / Δ_fine ≤ ALIGNMENT_TOLERANCE (init.f90:3194-3202).
   - ALIGNMENT_TOLERANCE defaults to 0.001 (cons.f90:598) and can be set on MISC (read.f90:1745).
   - Each coarse face may see fine cells from **one** mesh only: NOM_CHECK, ERROR(431) at init.f90:3184-3189.
   - Any integer ratio (2, 3, 5, …) passes. With stretched grids (TRN*), the real requirement is that faces coincide for each coarse cell, so the ratio need not even be constant along an interface.
   - The UG figure captions agree: "an integral number of fine cells abutting each coarse cell and each coarse cell only sees fine cells from one mesh" (FDS_User_Guide.tex:1085; the non-integral case is rejected at :1102).
   - Doc note: the UG writes the rejection criterion with "<" (FDS_User_Guide.tex:1120), but the code rejects when the value is **>** the tolerance (init.f90:3195). The UG sign is a typo.
   - The ALIGNMENT_TOLERANCE docstring "Maximum ratio of sizes of abutting grid cells" (cons.f90:598) is misleading: the value is a relative tiling tolerance, not a ratio cap.
4. **Directions.**
   - The two tangential directions may have **different** ratios. The only restriction is that one side may not be finer in one tangential direction and coarser in the other (1 % tolerance, init.f90:3210-3226). So 2:1 in y with 4:1 in z passes, and so does 1:1 in y with 2:1 in z.
   - The **normal** direction is not compared at all. Normal widths never enter ERROR(431), and the normal probe offset is MESH_SEPARATION_DISTANCE = min(1 mm, 0.05·CHARACTERISTIC_CELL_SIZE) (read.f90:916), where CHARACTERISTIC_CELL_SIZE is the global minimum per-mesh cell size (read.f90:796).
   - An isotropic AMReX ratio therefore satisfies FDS. An anisotropic `IntVect` ratio also does, provided no face is finer in one tangential direction and coarser in the other.
5. **Code paths downstream of the check are ratio-generic.** All EWC consumers loop over IIO_MIN:IIO_MAX (see the rule above).
   - The FireX GLMAT graph comments mention only 2:1 and 4:1 (pres.f90:5419-5422), but the loop is generic.
   - EWC_TYPE reconciliation uses AREA_RATIO thresholds of 0.9 (pres.f90:3722-3820). Against a coarser neighbour, AREA_RATIO is 1/r per coarser tangential direction, so the 0.9 cut misclassifies only 1<r<1.11. Such non-nesting ratios are already rejected by ERROR(431) on the coarse side (*derived*).
6. **Related topology limit (GLMAT only, *uncertain* effect).**
   - CHECK_UNSUPPORTED_MESH (pres.f90:3907-4048, called at pres.f90:3709) allocates `MESH_GRAPH(1:6,NMESHES)` (pres.f90:3940). It fills it with `COUNT=COUNT+1` over every CONNECTED_MESH, with no bound check (pres.f90:3942-3947), and then walks NMLOC=1,6 (pres.f90:3967).
   - A mesh connected to more than 6 others writes out of bounds. A coarse box abutting many fine boxes hits this easily.
7. **What a level jump actually costs.**
   - The second ghost layer is first-order at every level jump (wall.f90:351-386; see item 4 of the rule above).
   - The ghost-cell width is the **own** first-cell width, not the neighbour's: DX(0)=DX(1) (read.f90:1306, 1351, 1387).
   - Both points hold for any ratio ≠ 1, so they do not argue for ratio 2 over 4.

**Bottom line for the AMR design.** A cap of 4 violates nothing in FDS. Ratios 2 and 4, isotropic or not, are within what FDS already accepts between multi-meshes. Neither FDS nor our port needs powers of two, but AMReX does. The binding constraints are: fine cells tile each coarse face, each coarse face sees exactly one fine box (automatically true for AMReX boxes, which are coarsenable by the ratio), and no face is finer in one tangential direction and coarser in the other.

## Prototype P1 scoping: `kernel_footprint_mass.md` / `.csv`

- **Scope.** The per-routine footprint of every routine in mass.f90, plus the per-step callers and ghost producers:
  - callers: MASS_FINITE_DIFFERENCES at main.f90:767/948, DENSITY at main.f90:780/949;
  - ghost producers: MESH_EXCHANGE(1/4) at main.f90:790/959 and WALL_BC → ASSIGN_GHOST_VALUE.
- **Headlines.**
  - RHO(S) and ZZ(S) are read over `-1:IBAR+2` (mass.f90:70-76, 201-209), so **2 ghost layers are required**. This follows from the SUPERBEE default (cons.f90:710).
  - The ZZ ghost ring is *written* by the M_DOT_PPP add (mass.f90:477-479, 640-641).
  - There is no direct cross-mesh access; the only `MESHES(NM)%` reference is to M_DOT_PPP (mass.f90:477, 639).
  - CHECK_MASS_DENSITY depends on WALL_INDEX (mass.f90:831-836).
- **Shim requirements.** The md ends with an 11-item list of what the shim must provide per box.

## Section 4 (a): `routine_field_access.csv` (12,667 rows)

Columns:
- `field`, `parent_type` (MESH_TYPE or OMESH_TYPE), `category` (from `mesh_fields.csv`), `file`, `routine`;
- `phase`: runtime, setup, runtime+setup, or unreached, from `phase.py` (MAIN_LOOP is main.f90:695-1213);
- `access`: read, write, write_component, ptr_target, ptr_lhs, alloc, dealloc, inquiry, arg_read, arg_write, arg_readwrite;
- `via`: how the field was reached, e.g. `M:` for M%X, `ptr:` for a MESH_POINTERS name, `alias:` for a local pointer, `MESHES:`, `M2:`, `OM:`;
- `n_refs`, `first_line`, `lines`;
- `callee` and `arg_evidence`, filled for `arg_*` rows.

There is one row per (field, routine, access, via) combination.

**Coverage.**
- 11,746 MESH_TYPE rows and 921 OMESH_TYPE rows (133 OMESH fields).
- 582 routines.
- **All 444 MESH_TYPE components** are touched by at least one routine.
- Access counts: read 7,194; ptr_target 1,973; write 1,485; write_component 472; arg_read 461; alloc 431; inquiry 281; arg_write 120; arg_readwrite 114; dealloc 97; ptr_lhs 39.
- Phases: runtime+setup 5,468; setup 4,253; runtime 2,886; unreached 58.

**Argument resolution.** Every actual argument that is a mesh field is resolved through the callee's dummy argument by `argres.py`. **No row is left `arg_unknown`.** The evidence (`arg_evidence`) for the 695 arg rows is:
- explicit `INTENT(IN/OUT/INOUT)`: 415;
- MPI buffer rule (send buffer = read, receive buffer = write, MPI_IN_PLACE reduction = readwrite): 118;
- no INTENT but the dummy is only read in the callee body: 97;
- dummy passed on to a further callee and resolved there: 36;
- no INTENT and the dummy is assigned in the body: 29.

**Findings.**
- **271 MESH_TYPE fields are written by runtime-reachable code.** 164 more are written only during setup (geometry, index maps, boundary descriptors), which is the natural "rebuild on regrid" set.
  - The most widely written fields at runtime are CUT_CELL (30 routines), BOUNDARY_PROP1 (27), CUT_FACE (22), WORK1 (16), U/V/W/ZZ/FVX/FVZ/LAGRANGIAN_PARTICLE (15 each), WORK2 (14) and BOUNDARY_COORD (13).
  - The WORKn scratch arrays are shared by unrelated routines, so a per-box kernel cannot assume they are private to it.
- **9 components are dead:** declared but never written. DXMAX, DYMAX, DZMAX (mesh.f90:170ff), LBC2/MBC2/NBC2 (mesh.f90:92), PHI_W, PWORK4 and SAVE2. They can be dropped from the AMR box type.
- **How fields are reached:** through MESH_POINTERS names (`ptr:` 5,201 rows), `M%` (2,562), local aliases (2,392), and directly through `MESHES(...)%` (1,219). The MESH_POINTERS indirection is the dominant access path, which is why a per-box POINT_TO_MESH equivalent is the cheapest port route (see Item B).

## Section 4 (b): `cross_mesh_access.csv` (2,109 rows)

Columns: `file`, `routine`, `phase`, `exchange_routine` (yes if the routine is one of the exchange routines), `kind`, `chain`, `via`, `field`, `field_type`, `access`, `n_refs`, `first_line`, `lines`.

The table has every reference to storage that is not the current mesh's own MESH_TYPE. References are canonicalised to `MESHES(i)[%OMESH(j)]`, including MESH_POINTERS `OMESH` and local `TYPE(MESH_TYPE/OMESH_TYPE), POINTER` bindings. The kinds are:
- `omesh_copy` (791): this mesh's OMESH(NOM) copy of neighbour data.
- `all_meshes` (681): MESHES(i) inside `DO i=1,NMESHES`. The loop is often filtered, e.g. by a neighbour or process test (main.f90:2135; ccib.f90:4165-4171).
- `other_mesh` (474): MESHES(NOM) for a specific other mesh.
- `neighbour_omesh` (94): the **neighbour's** OMESH(NM), written directly on the same process.
- `literal_mesh` (69): e.g. MESHES(1).

Access counts: read 1,031; write 282; arg 224 (not INTENT-resolved here; see `routine_field_access`); alloc 160; ptr_target 156; ptr_bind 136; inquiry 80; dealloc 24; write_component 16. 582 rows are inside exchange routines.

**Runtime:** 826 rows in 69 routines (omesh_copy 442, other_mesh 203, all_meshes 113, neighbour_omesh 48, literal 20). Outside the exchange routines, the heaviest users are:
- CC_COMPUTE_VELOCITY_ERROR and WRITE_DIAGNOSTICS (30 rows each);
- CC_RESCALE_OMESH_STAGGERED_TO_CUTFACE (24), CC_MATCH_VELOCITY (22), GET_H_GUARD_CUTCELL (19);
- SEARCH_OTHER_MESHES (17);
- COMPUTE_VELOCITY_ERROR and VELOCITY_BC (16 each);
- GET_H_MATRIX_CC and GET_LINKED_VELOCITIES (15 each);
- SURFACE_HEAT_TRANSFER (14), MATCH_VELOCITY (13), ASSIGN_GHOST_VALUE (11).

**Findings.**
1. **Runtime writes into another mesh's MESH_TYPE (not OMESH).** On MPI runs these land in the *local replica* of the neighbour (partially allocated MESHES(NOM)). On a single process they land in the neighbour's real storage.
   - MESH_EXCHANGE CODE 18: `M4 => MESHES(NOM)`, `M4%OBSTRUCTION(OBST_INDEX)%MASS = …` (main.f90:3960-3964).
   - EXCHANGE_BACK_CFACE_DATA: `MESHES(NOM)%BOUNDARY_PROP1(B1)%TMP_G/HEAT_TRANS_COEF` (ccib.f90:358-359).
   - MESH_CC_EXCHANGE: writes the neighbour's CUT_FACE/RCEDGE/IBEDGE interpolation storage through `M1=>MESHES(NOM)` (ccib.f90:4171, 4235, 4354, 4371), and MESHES(NM)%CUT_CELL in the receive loop (ccib.f90:4776).
   - COPY_CC_MUNKH_TO_UNKH (called at pres.f90:4260): writes MESHES(NOM)%CUT_CELL(ICC)%UNKH and MESHES(NOM)%CCVAR(…,CC_UNKH) over the EWC range (ccib.f90:22893-22906). Its comment says "ghost cut-cell storage for mesh NOM". When NOM is on the same process this overwrites the real neighbour's arrays (*uncertain* whether the values then always agree).
   - EXCHANGE_DIAGNOSTICS writes CFL, DIVMX, … for every mesh (main.f90:4343-4365), and WRITE_DIAGNOSTICS resets MESHES(NM)%DT_RESTRICT_STORE for all meshes (dump.f90:4472). These make MESHES(:) a global diagnostics replica on rank 0.
   - EXCHANGE_GLOBAL_OUTPUTS writes MESHES(DV%MESH)%N_STRINGS/STRING (main.f90:4574-4577).
   - EXCHANGE_GEOMETRY_INFO allocates MESHES(NOM)%CELL_INDEX and CELL for every off-process neighbour and broadcasts the **whole** neighbour CELL_INDEX/CELL over MPI_COMM_NEIGHBORS (func.f90:5352-5368). It runs at setup and again whenever obstructions are created or removed (main.f90:1789-1790).
2. **Same-process shortcut.** In MESH_EXCHANGE's pack stage, a same-process neighbour's OMESH(NM) is written directly (`neighbour_omesh` rows: main.f90:3238-3240, 3342-3344, 3390-3392, 3415-3417, 3508, 3535, 3575-3581, 3590, 3600; ccib.f90:4683). The MPI unpack blocks are then skipped (`process_scope=mpi_only` in `mesh_exchange.csv`). A shim must therefore implement both paths or force one.
3. **Kernels read neighbour geometry and cell data, not just halo values.** Examples: ASSIGN_GHOST_VALUE reads MM%CELL_INDEX/CELL (wall.f90:361-371); AREA_RATIO reads M2%X over the EWC range (main.f90:2135); PRESSURE_SOLVER_COMPUTE_RHS reads DX_OTHER (pres.f90:124-139); SEARCH_OTHER_MESHES reads M2%XS/CELLSI (func.f90:5425-5428); others are SURFACE_HEAT_TRANSFER, MATCH_VELOCITY(_FLUX), COMPUTE_VELOCITY_ERROR, GET_H_MATRIX and the CC_* routines.
   - For every neighbour box, the AMR shim must expose: box extents and coordinates (X/Y/Z, DX/RDXN), CELL_INDEX + CELL(SOLID), and for CC_IBM also CCVAR/FCVAR/CUT_CELL/CUT_FACE indices.

## Section 4 (c): `mesh_exchange.csv` (85 rows) and `exchange_call_sites.csv` (89 rows)

**`mesh_exchange.csv` columns:** `routine`, `file`, `codes`, `stage`, `start_line`, `end_line`, `process_scope` (both / mpi_only / same_process_only / n_mpi>1), `condition`, `mesh_fields`, `omesh_fields`, `buffers`, `mpi_calls`, `purpose`, `purpose_src`.
- There is one row per CODE-guarded outermost IF block in POST_RECEIVES (main.f90:2945ff), MESH_EXCHANGE (main.f90:3117ff) and MESH_CC_EXCHANGE (ccib.f90:3946-5144).
- The stage comes from the source markers "Start the communications" (main.f90:3618) and the receive loop (main.f90:3686), and from the ccib equivalents (ccib.f90:4725, ~4750, and TINTP at ~4991).
- `purpose` is curated, with a citation in `purpose_src`.

**`exchange_call_sites.csv` columns:** `file`, `line`, `caller`, `phase`, `callee`, `code_or_args`, `guard`, `time_loop_stage` (setup / predictor / corrector / pressure_iteration, or blank outside MAIN_LOOP stages), `prev_call`, `next_call`.
- It covers 28 callees: MESH_EXCHANGE 37 calls, POST_RECEIVES 9, EXCHANGE_GEOMETRY_INFO 5, MESH_CC_EXCHANGE 4, and 24 other `*EXCHANGE*` routines with 1-3 calls each.

**MESH_EXCHANGE codes:** 0-11 and 14-20 are used. **Codes 12 and 13 are never called and have no block.**

| code | content | MPI request | call sites (main.f90) |
|---|---|---|---|
| 0 | allocate buffers, persistent requests | – | setup 309-310 |
| 1 / 4 | RHOS/ZZS/D (+MU, KRES, Q) two cells deep, predictor / corrector | REQ1 | 790 / 959; setup 399-400 |
| 2 | radiation IL | REQ5 | 1040, 1112; setup 510 |
| 3 | HS, US/VS/WS (end of predictor) | REQ3 | 909; GLOBAL_MATRIX_REASSIGN 1822-1828 |
| 5 | FVX/FVY/FVZ + H/HS inside the pressure iteration | REQ7 | 1617, 1636, 1663, 1687 |
| 6 | H, U/V/W (end of corrector) + back-wall data | REQ3, REQ6 | 1014, 1105; setup 401, 2559 |
| 7 / 11 | particle orphan counts / particle buffers | REQ2 / – | 801, 1000 / 803, 1002 |
| 8, 9, 19, 10 | back-wall setup: counts, index lists, buffer sizes, persistent REQ6 | REQ6 | INITIALIZE_BACK_WALL_EXCHANGE 2369-2553 |
| 14 | level-set values | REQ14 | 794; setup 429, 433 |
| 15-18 | OBST mass loss and new OBST mass (PKG8) | REQ15 | 1074-1079 |
| 20 | particle drag FVX_D/FVY_D/FVZ_D | REQ4 | 991 |

**MESH_CC_EXCHANGE** is called from MESH_EXCHANGE when CC_IBM (main.f90:3135), and at setup from INIT_CUTCELL_DATA (ccib.f90:7975-7977).
- It acts on codes 1, 3, 4, 5 and 6.
- It returns immediately for 0, 2 and >6, and for 3 when CALL_FROM_GLMAT_SETUP (ccib.f90:3976-3979).
- It uses REQ11/REQ112/REQ12/REQ13 (ccib.f90:4726-4747).
- Rows: pack_send 16, start_wait 4, unpack 10, post_process 1 (ccib.f90:4991-4997): FILL_GCCUTCELL_SPECIES after codes 1/4, and CC_H_INTERP + CC_RHO0W_INTERP after codes 3/6. ELSEIF branches are merged into the row's `codes`, and single-line early-RETURN guards (ccib.f90:3976, 3979) are labelled as such.

**Time-loop sequence** (from `exchange_call_sites.csv`):
- Predictor: 1 (790), 14 (794), 7 (801), 11 (803), 3 (909).
- Corrector: 4 (959), 20 (991), 7 (1000), 11 (1002), 6 (1014), 2 (1040), 15-18 (1074-1079), 6 (1105), 2 (1112).
- Pressure iteration: 5.
- EXCHANGE_DIVERGENCE_INFO is called at main.f90:594, 845 and 1065; HT3D_TEMPERATURE_EXCHANGE at 1016.

**AMR implications.**
- Codes 1/4, 3/6, 5 and 20 are halo fills with a fixed per-wall footprint over the EWC range. They map onto FillPatch/FillBoundary plus a coarse/fine averaging step.
- Codes 7/11 (particles), 15-18 (OBST mass), 2 (radiation) and the back-wall set (8/9/19/10/6) exchange **object lists** keyed by mesh number and item index. They need AMReX particle redistribution or custom neighbour lists that are rebuilt on every regrid (see `mesh_id_dependencies.csv`).

## Section 4 (d): Snapping and uniform-grid assumptions (`uniform_grid_assumptions.csv`, 231 rows)

Columns: `category`, `file`, `line`, `routine`, `phase`, `expression`, `stretch_aware` (yes / no / n/a), `curated`, `amr_impact`, `notes`.

Automated categories. Only mesh objects count: the component parent must be MESHES or a TYPE(MESH_TYPE) variable.

| category | rows | meaning |
|---|---|---|
| `snap_computational` | 103 (69 stretch-aware) | index from physical coordinate, e.g. `NINT(GINV(X-XS,1,NM)*RDXI)` or `…*DXI` |
| `coord_lookup_table` | 34 | CELLSI/CELLSJ/CELLSK lookups (stretch-aware; e.g. func.f90:5425-5428, 6775-6777) |
| `first_cell_width` | 31 | DX(1)/DX(IBAR) taken as the mesh spacing |
| `interface_min_corner` | 22 | imported pick_min rows of `interface_averaging_sites.csv` |
| `representative_size` | 13 | CELL_SIZE, DXMIN, … as a single per-mesh length (read.f90:791-796, 1280-1374) |
| `global_min_size` | 9 | uses of the global CHARACTERISTIC_CELL_SIZE |
| `ghost_spacing_copy` | 6 | DX(0)=DX(1), RDX(0)=1/DX(1) and the y/z equivalents (read.f90:1306-1389) |
| `mesh_sequence` | 6 | neighbour = mesh NM±1 (tunnel solver) |
| curated | 8 | see below |

Phases: setup 157, runtime 46, runtime+setup 23, unreached 5.

**A-20: two upstream slips in the thin-OBST collapse (READ_INPUT::READ_OBST, setup).** Both were re-verified on FireX and are listed as `upstream_slip` rows.
- **read.f90:11193**: `IF(GINV(XB3-XS,2,NM)-REAL(OB%J1,EB) < REAL(OB%J2,EB) - GINV(XB4-YS,2,NM))` uses **XS** where **YS** is intended. The guard two lines above uses `XB3-YS` (read.f90:11192).
- **read.f90:11210**: `IF(GINV(XB5-ZS,3,NM)-REAL(OB%I1,EB) < REAL(OB%I2,EB) - GINV(XB6-ZS,3,NM))` compares **OB%I1/OB%I2** in the z branch where **OB%K1/OB%K2** is intended. It then assigns K2=K1 or K1=K2 (read.f90:11211-11213). If the x branch already collapsed I1=I2, the z decision depends on x data.
- **Both slips are in master:** ce1f659 read.f90:11060 and :11077, same text. FireX did not introduce them.
- **Bit-parity requirement.** Level 0 of the AMReX port **must reproduce both slips exactly** (same XS and same I1/I2 operands), or thin-OBST placement changes and bit-parity with FDS master is lost. Fix them only behind a flag, upstream, and never silently in the port.
- **AMR hazard.** The y decision depends on the mesh origin **XS**. If a refined box evaluates it with the box's own XS, the same OBST can collapse to different faces on different boxes, levels or regrids. The level ≥1 path must evaluate the collapse with a defined, box-independent origin (*design decision needed*).
- **Related `unit_mix` row (read.f90:11176, same pattern in all three branches: 11176, 11193, 11210).**
  - `GINV(XB1-XS,1,NM) - REAL(OB%I1,EB)` subtracts an integer cell index from a computational coordinate in metres. OB%I1 = NINT(GINV(XB1-XS,1,NM)*RDXI) (read.f90:11157), with RDXI = 1/DXI (read.f90:780).
  - The comparison is dimensionally correct only for 1 m cells. For dx < 1 m, thin OBSTs collapse to the lower face almost regardless of where they sit inside the cell, and the choice is resolution-dependent.
  - The guard threshold `0.25_EB/RDXI` (read.f90:11175) is in metres and is correct.
  - This is a third upstream defect (*derived by reading, not run*). It carries the same bit-parity rule for level 0.

**Other findings.**
1. **Ghost spacing is the own first-cell width**, DX(0)=DX(1) (read.f90:1306, 1351, 1387), not the neighbour's. Code that uses RDXN(0)/DX(0) at an interpolated face therefore sees a symmetric spacing even at a level jump. The exception is PRESSURE_SOLVER_COMPUTE_RHS, which fetches DX_OTHER from the neighbour's min-corner cell (pres.f90:124-139).
2. **Global resolution scalars.** The global minimum CHARACTERISTIC_CELL_SIZE (read.f90:796) sets:
   - the initial DT (read.f90:187);
   - MESH_SEPARATION_DISTANCE and NEIGHBOR_SEPARATION_DISTANCE (read.f90:916-917);
   - the default VELOCITY_TOLERANCE and PRESSURE_TOLERANCE (read.f90:10161-10162);
   - HRRPUV cutoffs (vege.f90:1267, dump.f90:2265).

   With AMR these scale with the finest level present at setup, and they do not change when refinement is added later (*design decision*).
3. **Snapping is per-mesh-lattice.** The 103 snap sites snap to the lattice of the mesh doing the snapping. They are stretch-aware where GINV is used, which is 69 of 103; the other 34 multiply by DXI, which is equivalent on uniform meshes. Coarse and fine boxes therefore snap the same object to different faces. Examples are OBST/VENT/DEVC/PATCH bounds (e.g. velo.f90:1098-1103, 1131-1136, 1165-1170). Levels must agree on a snapping rule, e.g. snap on the coarsest level and refine (*design decision*).
4. **Solver restrictions.**
   - ERROR(425): a mesh can be stretched in at most two directions (init.f90:2366).
   - The tunnel solver assumes x-consecutive mesh numbers: I_OFFSET (main.f90:1316) and MESHES(NM±1)%DX (pres.f90:541-542, 686-690), with a YS/YF equality check (read.f90:10182). It is incompatible with AMR box numbering.
5. **First-cell-width sites** are mostly output: SANDIA_OUT (turb.f90:2420-2435), SMOKE3D VTKHDF (dump.f90:5197) and SPECTRAL_OUTPUT, whose comment reads "obviously, assumes uniform grid spacing" (turb.f90:2357; unreached here, NM>1 returns at turb.f90:2354). The runtime numerics exception is the pressure RHS above.

## Section 4 (e): Global reductions (`global_reductions.csv`, 259 rows)

**Purpose.** This section is the complete list of reductions that exact fixed-point accumulation must cover. It lists every MPI reduction and every cross-mesh accumulation that feeds physics or time-step control, plus the bookkeeping and output-only ones, each flagged so they can be told apart.

**Columns.**
- `kind`, `file`, `line`, `routine`, `phase`, `call`
- `variables`, `mpi_op`, `datatype`
- `order_dependent`: `yes` means a floating-point SUM, whose result depends on summation order and so on the decomposition/rank count. `no` means MIN/MAX/MINLOC/MAXLOC, integer/logical ops, or gathers.
- `feeds`
- `physics_or_output_only`: `physics`, `physics_via_devc` (reaches CTRL/setpoints), `dt_control`, `bookkeeping` or `output_only`
- `feeds_src`: the line where the reduced value is consumed
- `per_process_presum`: where each process first sums over its own meshes
- `notes`

**Kinds.**
- `mpi_collective` (119): ALLREDUCE 82, GATHERV 19, ALLGATHERV 10, REDUCE 6, EXSCAN 1, ALLGATHER 1.
- `serial_loop_all_meshes` (42): accumulation in a loop over 1..NMESHES / N_MPI_PROCESSES / N_SUBDEVICES.
- `per_process_presum` (98): accumulation into a module/PROGRAM global, either inside an own-mesh loop or inside a routine called once per own mesh.

**Flags.**
- All rows: physics 104, output_only 74, bookkeeping 62, dt_control 11, physics_via_devc 8.
- MPI rows only: output_only 40, bookkeeping 37, physics 34, dt_control 5, physics_via_devc 3.
- 134 rows are order-dependent in total.

### Order-dependent sums that feed physics: 91 rows, 13 quantities

The 91 rows are 18 MPI + 68 per-process pre-sum + 5 serial.

| quantity (feeds) | MPI reduction | per-process pre-sum over own meshes | phase |
|---|---|---|---|
| DSUM/PSUM/USUM → D_PBAR_DT (zone background-pressure rise, divergence) | main.f90:2038-2040 | divg.f90:739, 750, 765-766, 776-777, 1506; plus 745/747 through the dummies of ADD_CUTCELL_PSUM (ccib.f90:739) and ADD_LINKEDCELL_PSUM (ccib.f90:719) | runtime |
| RAD_Q_SUM/KFST4_SUM → RTE_SOURCE_CORRECTION_FACTOR (used main.f90:1767) | main.f90:1760-1761 | radi.f90:4040-4041, 4076-4077, 4120-4121 | runtime |
| NODE_PROPERTIES → HVAC node boundary conditions | main.f90:4979 | hvac.f90:2415-2533 (39 rows; HVAC_BC_IN, run once per mesh from the main.f90:815-822 loop) | runtime |
| TC_ARRAY → DEVC spatial statistics → CTRL/setpoints | main.f90:4479 (MPI_SUM under CASE(1), main.f90:4442) | subdevice loop main.f90:4456-4459 | runtime |
| SUM_FH → GLMAT RHS compatibility (mean removed per zone) | pres.f90:3396, 3403 | pres.f90:3382-3383, 3392-3393; 3402 (SUM intrinsic) | runtime |
| SUM_GAUGE → GLMAT pressure gauge shift | pres.f90:3546 | pres.f90:3509-3510, 3527-3528, 3540-3541 | runtime |
| SUM_XH → GLMAT zero-mean shift | pres.f90:3554 | pres.f90:3553 (SUM intrinsic) | runtime |
| VENT_MEAN_SUM/CNT → synthetic eddy method (SEM) mean inflow | turb.f90:1889-1890 | turb.f90:1877-1884 | runtime |
| VENT_TOTAL_AREA → B1%AREA_ADJUST of split vents | main.f90:5059 | read.f90:12864 (serial) | setup |
| FDS_AREA_GEOM → GEOM VOLUME_FLOW U_NORMAL_0 | geom.f90:12310 | geom.f90:12306 | setup |
| GEOM_AREA_SURF_OLD/NEW → CFACE AREA_ADJUST | geom.f90:14619-14620 | geom.f90:14345, 14401 | setup |
| P_ZONE volume (REAL_BUFFER) → pressure-zone volumes | main.f90:2731 | main.f90:2699 | setup |
| MULTIPLIER%FDS_AREA → B1%AREA_ADJUST of shaped multiplied OBSTs | **none** | init.f90:366-387 (12 rows) | setup |

**Notes.**
- The pre-sums on pres.f90/turb.f90 local variables are not rows of their own; they appear in the `per_process_presum` column of the MPI row.
- The ULMAT sums (pres.f90:1684-1876) cover a single mesh and are not global.
- **Implication:** for bit-reproducibility across box layouts, the fixed-point accumulation must cover both stages, i.e. the per-process pre-sum and the MPI SUM. Making the MPI SUM exact is not enough.
  - The order within one mesh (the cell loop) also changes when boxes are regridded. The per-cell contributions must therefore be accumulated exactly too, not only the per-mesh partials.

### Upstream defect, run-confirmed: MULTIPLIER%FDS_AREA

- ADJUST_OBST_SHAPE_AREA (init.f90:339-432) is called once per mesh (init.f90:917).
- It accumulates MR%FDS_AREA in the global MULTIPLIER array (init.f90:366-387), one entry per MULT, never reset between meshes and never MPI-reduced.
- It immediately uses that array for B1%AREA_ADJUST (init.f90:405-426).
- The routine is byte-identical in master ce1f659.

**Confirmed 2026-09-25** with a MULT-voxel `SHAPE='CYLINDER'` (R=H=0.1 m; top and side emit different species at 0.01 kg/m²/s). The inputs are in `(local project tools directory)/inventory/runs/fds_area/`. Effective area ratios versus exact (top / side):

| variant | top | side |
|---|---|---|
| 1 mesh | 1.000 | 1.000 |
| 2 meshes, 1 process | 0.750 | 1.500 |
| 2 meshes, 2 processes | 1.000 | 2.000 |
| off-centre split (x=0.14), 1 process | 0.871 | 1.700 |
| two cylinders sharing one MULT, 1 mesh | 0.500 | 0.500 |

- FireX and master binaries give identical numbers, and a face-count model of the mechanism reproduces them to 8 digits.
- The ready-to-file draft is `upstream_issue_candidates.md` #1.
- For the AMR port, the fix needs a global (exact, fixed-point) shape-area sum keyed by OBST line.

### Order-independent controls (MIN/MAX)

- Time step: DT = MINVAL(DT_NEW) (main.f90:715, 741, 891; gathers at main.f90:738, 884).
- CFL-norm switch: MAXVAL(MAX_CELL_ASPECT_RATIO) (main.f90:161-163).
- Pressure-iteration exit tests: MAXVAL (main.f90:1724-1740; gather at main.f90:1702).
- Global CHARACTERISTIC_CELL_SIZE: MIN (read.f90:796).
- DEVC MIN/MAX statistics: TC2_ARRAY with MINLOC/MAXLOC (main.f90:4481). On ties the lowest rank wins, so the result depends on the decomposition but not on floating-point rounding.

**Output-only sums.** Q_DOT, M_DOT, Q_DOT_SUM, M_DOT_SUM and MASS_DT (MPI_REDUCE, main.f90:4602-4641) feed only the .csv outputs.

### Exclusions (printed by `global_reductions_csv`)

- integer/logical accumulations: 416 (exact in any order)
- timers: 145
- linear-system row assembly: 44 (single-mesh rows; *uncertain* at box faces)
- coordinate nudges off a mesh face: 18
- input-record, property, CONTROL or HVAC-network loops that are not over meshes: 29
- per-mesh storage and slots: 14
- constant increments: 1

### Limitations specific to this table

- **INTENT(INOUT) dummies.** Accumulation through them is not scanned in general. The known cases were added by grep: ADD_CUTCELL_PSUM, ADD_LINKEDCELL_PSUM, ADD_Q_DOT_CUTCELLS, HVAC NODE_* pointers, and the GEOM_AREA_SURF dummies.
- **IF-guarded MIN/MAX selection** (`IF (x<m) m=x`) is only partly caught.
- **Reductions inside external solvers are not listed.** These are the MKL cluster sparse solver for GLMAT (pres.f90:3435, 3439, 4663), HYPRE PCG (pres.f90:1760, 3461) and PARDISO (pres.f90:1742, 2915, 2939). Their internal dot products and norms are order-dependent by process count, so the GLMAT/ULMAT pressure is not bit-reproducible across decompositions regardless of the fixed-point sums above.
- **Gathers are marked `no`.** The tunnel TDMA system (pres.f90:633-647) and the GLMAT RHS (pres.f90:3414, 3417) are then solved serially in rank order: deterministic for a fixed decomposition, but tied to the box layout (*uncertain*).

## Section 4 (f): Refinement ratios in the FireX V&V inputs (`mesh_ratio_cases.csv`, `mesh_ratio_cases_rollup.csv`)

**What the parse does.** The script is `tools/inventory/mesh_ratio_cases.py`: stdlib only, no FDS runs, about 20 s single-process. It covers all 4,073 `.fds` files under FireX `Verification/` (941) and `Validation/` (3,132).
- **CATF:** `&CATF OTHER_FILES` are inlined.
- **MULT expansion:** MESH lines with a MULT are expanded as READ_MESH does (read.f90:633-658). The MULT defaults, SKIP ranges and N_LOWER/N_UPPER sequential mode follow read.f90:1643-1728. XB is reordered as in CHECK_XB (func.f90:327-337).
- **Stretching:** TRNX/Y/Z applicability follows READ_TRAN (read.f90:1032-1062).
- **Shared faces** use FDS's own tolerance:
  - |XF_a − XS_b| < MESH_SEPARATION_DISTANCE = MIN(1e-3, 0.05·CHARACTERISTIC_CELL_SIZE) (read.f90:916, 790-796);
  - the tangential overlap must exceed that distance in both directions.
- **Ratios:** `ratio_t1`/`ratio_t2` are coarse/fine per tangential axis.
- **`integer_tiling`** means |r − round(r)| ≤ 0.001 and both overlap bounds lie on both grids within 0.001·Δfine. This is the per-coarse-face ALIGNMENT_TOLERANCE test (init.f90:3194-3202; cons.f90:598).
- **Cross-check:** mesh counts agree with the V&V Lead's `vv/verification_case_survey.csv` for **all 941 Verification cases (0 mismatches)**.

**Classes.** Each shared face gets one class, applied in this order:
1. `stretched`: a TRN* in a tangential direction on either mesh. Ratios are then in computational (uniform) coordinates and `integer_tiling` is `n/a`.
2. `non_tiling`
3. `direction_dependent`: the two tangential ratios differ.
4. `uniform`, `supported_2_4`, `ratio3` or `other_nonpow2`, by the common ratio.

In 2-D inputs (IJK(2)=1 on every mesh; read.f90:703-707), the single y cell is ignored for the class (409 faces).

**Rollup.** The worst class per case uses the order uniform < supported_2_4 < ratio3 < other_nonpow2 < direction_dependent < stretched < non_tiling. Cases without a face class get one of these labels instead:
- `single_mesh`
- `no_shared_faces`: all meshes separated
- `embedded_or_overlapping`: no shared faces but overlapping volumes
- `no_mesh`: CATF include fragments
- `parse_error`: two `Build_Input_Files/*Template.fds` with `param_*` placeholders

The column `n_overlap_pairs` counts embedded or overlapping mesh pairs in every case.

### Counts

**Shared faces (143,708 rows):**

| class | faces |
|---|---|
| uniform | 121,886 |
| supported_2_4 | 12,772 (4:4 9,177; 2:2 3,586; 2-D 1:2 9) |
| stretched | 7,453 |
| ratio3 | 961 |
| other_nonpow2 | 378 |
| direction_dependent | 258 |
| non_tiling | **0** |

**Cases (4,073), worst class per case:**

| worst class | total | Verification | Validation |
|---|---|---|---|
| single_mesh | 2,010 | 748 | 1,262 |
| uniform | 1,438 | 140 | 1,298 |
| supported_2_4 | 301 | 5 | 296 |
| no_shared_faces | 208 | 44 | 164 |
| stretched | 59 | 0 | 59 |
| other_nonpow2 | 13 | 2 | 11 |
| ratio3 | 11 | 0 | 11 |
| direction_dependent | 7 | 0 | 7 |
| embedded_or_overlapping | 2 | 2 | 0 |
| no_mesh | 22 | 0 | 22 |
| parse_error | 2 | 0 | 2 |

- 2,039 cases have more than one mesh.
- The largest tangential ratio anywhere is **8** (BST_FRS_6, meshes 12/16, 8:4).
- No input has a power of two above 4 as a same-in-both-directions ratio.
- 16 single-mesh Verification inputs use TRN stretching; see the rollup `notes`.

### Every case outside `uniform` / `supported_2_4`

**ratio3 (11):**
- Validation/Convection: `impinging_jet_Re_1e5_{coarse,medium,fine}`, `impinging_jet_Re_4e5_{coarse,medium,fine}` (35 meshes, 3:3).
- Validation/McCaffrey_Plume: `McCaffrey_{14,22,33,45,57}_kW_45` (158 meshes, 3:3).

**other_nonpow2 (13), all 5:5:**
- Verification/Thread_Check: `race_test_1`, `race_test_4` (6 meshes, 0.05 m | 0.01 m).
- Validation/Aalto_Woods: `Roomcorner_M12`, `Roomcorner_modified` (4 meshes).
- Validation/TAMU_Jet_Fires: `TAMU_T16`, `TAMU_T17`, `TAMU_T18`, `TAMU_T19a` … `TAMU_T19f` (126 meshes).

**direction_dependent (7):**
- Validation/BST_FRS_wood_cribs: `BST_FRS_6` (4:2 and **8:4**; also 2:2 and 4:4 faces).
- Validation/NIST_Pool_Fires: `NIST_Methanol_1m_pan_1cm_grid`, `NIST_Methanol_1m_pan_1cm_grid_predicted` (3:4).
- Validation/UMD_Line_Burner: `methane_dx_1p25cm`, `propane_dx_1p25cm` (2:3), and `methane_dx_p625cm`, `propane_dx_p625cm` (4:6).

**stretched (59):**
- Validation/Crown_Fires: all 32 `{Heil,Pike1,RedF,UNC}_{pre,post}_{2,4,9,13}ms` (72 meshes).
- Validation/FHWA_Tunnel: `IFAB-07`, `-08`, `-09`, `-10`, `-11`, `-13`, `-14`, `-15`, `-19`, `-22`, `-24` (12 meshes).
- Validation/LNG_Dispersion: `Burro3`, `Burro7`, `Burro8`, `Burro9`, `Coyote3`, `Coyote5`, `Coyote6`, `Falcon1`, `Falcon3`, `Falcon4`, `MaplinSands27`, `MaplinSands34`, `MaplinSands35`.
- Validation/NIST_USFS_Camp_Swift: `CS_BB1_0p10dx` (1,024 meshes), `CS_BB1_0p25dx` (400).
- Validation/Restivo_Experiment: `Restivo`.

**embedded_or_overlapping (2):** Verification/Adaptive_Mesh_Refinement `ns2d_16_emb_1to1_refinement` (embedded, 1:1) and `ns2d_16_emb_1to2_refinement` (embedded, 2:2 in x/z).

**Embedded meshes inside otherwise `uniform` cases (2):** `Adaptive_Mesh_Refinement/random_meshes` (8 overlapping pairs, including 2:2) and `Pressure_Solver/dancing_eddies_embed` (1 pair, 2:2).

### Cases the V&V anchors may use

These were checked against `docs/vv/case_inventory.csv`, `test-plan.md`, `requirements.md` and `pressure/01-amr-mapping-spec.md` by grepping the stem of every case outside uniform/supported.

- **Heskestad** (`Validation/Heskestad_Flame_Height`, all 57 inputs): `single_mesh`. No interfaces, so there is no ratio issue. The NFR-032 level-0 input and the A-19 fine reference live in `vv-runs/inputs/A-19/`, outside FireX, and were not parsed.
- **shunn3**: `single_mesh` (all N and CFL variants, plus the 12 `Complex_Geometry/shunn3_*_cc_exp_*`); `shunn3_4mesh_{32,64,128,256,512}` are `uniform` (4 faces each).
- **FM_Burner** (21 inputs):
  - the 1 cm and 2 cm inputs are `uniform`;
  - all 7 `*_5mm` inputs are `supported_2_4`, with 2:2 (128–144 faces). The m2 block of 0.005 m meets m1 at 0.01 m.
- **Flagged from the V&V lists:**
  - **`Thread_Check/race_test_1`** (Tier 2) and **`race_test_4`** (Tier 2-optional): these are the NFR-011 race/deadlock cases, listed in requirements.md:47 and :302. Their 5:1 interfaces are **outside a 2/4 cap**, so the AMReX code could not run them with the same mesh layout under that cap. Baseline-FDS runs are unaffected.
  - **`ns2d_16_emb_1to2_refinement`** (Tier 1, anchor=yes) and **`ns2d_16_emb_1to1_refinement`** (Tier 2): embedded, overlapping meshes, not face-sharing. FDS one-way couples embedded meshes. requirements.md:111 already treats emb_1to2 as a non-FR-016 case.
  - **`random_meshes`** and **`dancing_eddies_embed`** (requirements.md:42, refinement analogues): uniform at shared faces but containing embedded 2:2 meshes.
  - The other refinement analogues are `supported_2_4`: `ns2d_16_int_1to2_refinement` (FR-016 primary, 8 faces 2:2 in 2-D), `dancing_eddies_uglmat_refine` (2:2) and `duct_flow_uglmat_refine` (4 faces 2:2).
  - No other case in the V&V lists falls outside uniform/supported.

### Limitations

- Periodic wrap-around neighbours are not paired.
- Stretched faces get no tiling test.
- Only `.fds` files are parsed. CATF fragments are inlined into the including case but are also listed on their own as `no_mesh` or `single_mesh`.
- Namelist parsing is textual. Array sections like `XB(2)=` and repeat counts `n*v` are handled; anything else is a `parse_error`.
- Mesh counts are verified against the survey for Verification only.

## Validation (FireX)

- **mesh_fields.**
  - 444 MESH_TYPE components (mesh.f90:16-354). This matches an independent naive regex count.
  - 32 random rows (spot check seed 20260925, re-run on FireX) agree with the cited declaration and allocation lines: **32/32**.
- **allocation_sites.**
  - An `rg` for `ALLOCATE(M%X(` / `ALLOCATE(MESHES(n)%X(` finds 415 distinct lines. **All 415 have a scanner record.**
  - A continuation-line ALLOCATE attribution bug was fixed in `alloc_scan.py`.
  - turb.f90:37-46 is commented out (`!ALLOCATE`), so it is correctly skipped.
- **pressure_fields_access.** 34 random rows (seed 20260925) were checked against the source.
  - All `first_line` citations contain the field.
  - All access classes are consistent, e.g. arg_readwrite PRHS at pres.f90:349 and OM%FVX at ccib.f90:15813.
  - An earlier 40-row check on the old base found and fixed one DEALLOCATE-behind-IF misclassification.
- **kernel_footprint_mass.** All 1,159 line citations were checked: each cited line contains the field (0 misses).
- **obst_wall_cface_indexing.** All curated citations were verified line by line and corrected where they had drifted.
- **interface_averaging_sites.** The coarse/fine rule sites listed below were read by hand. The per-row `operation` classes are heuristic (*uncertain*), except where the `notes` column carries curated semantics.
- **Section 4 spot check 1** (`spotcheck_sec4 20260925`): 40 random rows from routine_field_access, cross_mesh_access, mesh_exchange, exchange_call_sites and uniform_grid_assumptions, each checked against the cited source lines.
  - Result: **39/40**. The miss was a READ input list at dump.f90:4253, which cross_mesh had classified as a write.
  - After the fix in `cross_mesh.py`: 40/40.
- **Section 4 spot check 2** (`spotcheck_sec4 20260926`), 40 new rows: **39/40**.
  - The miss was a uniform_grid row at main.f90:2136, imported as `interface_min_corner` but actually the AREA_RATIO span.
  - Fixed in `interface_sites.py` with the new kinds `range_span` (main.f90:2132-2136) and `range_extent` (main.f90:2118-2128). interface_averaging_sites and uniform_grid_assumptions were regenerated (164 and 231 rows).
- **global_reductions.** Each of the 18 order-dependent physics MPI rows and its pre-sum lines were read by hand. TC2_ARRAY (main.f90:4481) was corrected from "order-dependent" to MINLOC/MAXLOC.
  - The completeness of the MPI rows was checked against `rg "CALL MPI_(ALL)?REDUCE|MPI_(ALL)?GATHERV?|MPI_EXSCAN"` over FireX. All 119 source call sites are in the table. Three rows cite the first line of a continued statement (`IF (...) &`): geom.f90:12310, main.f90:884 and pres.f90:4689, where the CALL itself is on the next line.
  - The per-process pre-sum rows are heuristic (*uncertain*); see the limitations in Section 4 (e).
- **README citations.** All were re-checked on FireX for this checkpoint. Old-base citations were translated with `remap.py`; the table is in `base_delta.md`.

## Known limitations of the static analysis

- Both `#ifdef` branches are scanned. For duplicate declarations, the first is kept and the others are listed in `notes`.
- Allocation-on-assignment is not detected.
- The MESH_POINTERS binding is assumed flow-insensitively: a MESH_POINTERS name means "the current mesh".
- Local pointer aliases use a textual reaching-definition heuristic. Loop-carried aliases are not modelled.
- Runtime/setup reachability comes from a static call graph. Calls through procedure pointers and generic interfaces may be missed (*uncertain*).
- Section 4 tables inherit the flow-insensitive pointer binding. A local pointer rebound in a loop is attributed to every target it is bound to.
- "All meshes" loops are recognised from their bounds (1/0/2 .. NMESHES / N_MPI_PROCESSES). A loop over a filtered list, or one with a CYCLE on ownership, is classified by its bounds only (*uncertain*).
- `cross_mesh_access.csv` rows of access `arg` are not INTENT-resolved; `routine_field_access.csv` is.
- The implicit refinement-ratio limit (r < 40) is derived from init.f90:3161-3208 and was not confirmed by a run.
- global_reductions: accumulation through INTENT(INOUT) dummies (only the known cases), partial MIN/MAX detection, and no external-solver internals. See Section 4 (e).
- `kernel_footprint` index extents are symbolic, from the enclosing DO bounds and linear subscripts (`I±c`, IBAR/IBP1 terms). Wall-coordinate subscripts are reported as `wall:<expr>`, unevaluable subscripts as `?`, and whole-array references as `whole`.

## Remaining (optional)

1. Run a two-mesh case with refinement ratio > 40 to confirm the implicit ERROR(431) limit, and a case with more than 6 connected meshes to see the effect of the MESH_GRAPH bound in CHECK_UNSUPPORTED_MESH (pres.f90:3940-3967).
2. (Done 2026-09-25: the MULTIPLIER%FDS_AREA defect is confirmed by runs; see Section 4 (e) and `upstream_issue_candidates.md`.)
3. Do a further spot check of `global_reductions.csv` per_process_presum rows and of `routine_field_access.csv` arg rows.

## Checkpoint 2026-09-25 (low-spend pause until about Sep 29)

Status: the inventory tables and A-33 (mesh_ratio_cases) are complete. Nothing is running.

Findings added tonight (FireX 36975d7):
- The pressure `SOLVER` string is parsed by `DEFINE_PRES_METHOD` (`func.f90:7143`, UGLMAT cases at `7154-7159`, unknown string gives ERROR(371) at `7182`). It is called from `READ_PRES` (`read.f90:10073-10086`), which reads every `&PRES` line. The block at `read.f90:10090-10150` is commented out; an earlier citation of it was wrong.
- `Adaptive_Mesh_Refinement/ns2d_16_int_1to2_refinement.fds` is fully periodic (`PERIODIC_TEST=1` and periodic vents, lines 20 and 26-29) and exists only at N=16. Its existing `&PRES` is at line 22.
- The Smokeview mesh list (`GRID`/`PDIM`/`TRNX`/`OBST`) is written once at setup (`main.f90:216-221`, `dump.f90:2457-2541`). After that the `.smv` file is only appended to (`main.f90:255`, `651`), so FDS never changes the mesh set during a run. This matters for FR-072.

Open items to resume:
- A-34: check whether slice output reaches `JJ=0` for `ADV_FX` (`dump.f90:9850-9898`), and check the CC_IBM `CC_REGFACE` lists (`ccib.f90:8023-8034`).
- Add the list of fused face nests to `kernel_footprint_mass.md` (`mass.f90:328-341`, `463-468`).
- Hand-check `race_test_4` for remeshing under D-030.
- UGLMAT with overlapping meshes (the `emb_*` cases) is not verified.
- FR-072 support: list every output routine that loops over meshes, if the Architect or Integration Lead asks for it.
