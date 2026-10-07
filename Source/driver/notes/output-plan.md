# Output plan for the driver: FR-072 / ADR-004 (Role 1 side)

Status: **Step A (0011) and Step B tables (0012) APPLIED** (`c0d3573c7d`; Architect + Legacy Mapper approval). **Step C landed** in `OutputPlan.{H,cpp}`: static layout + fill + D6 `composite_mass_exact` (unit check `test_output_plan` / `run_output_plan_check.sh`). `bind_ready()` stays false until `N_OUT_MESHES > NMESHES`. Inventory membership locked to Legacy Mapper FR-072 confirm vs `docs/upstream-patches/line-check-output-plan.md` (no further Mapper ask). Prefer option (a) `OUT_MESH` / `N_OUT_MESHES`; do **not** expand 0012 into remaining mesh loops until `N_OUT_MESHES > NMESHES` is validated (Mapper: tables-only until then). Basis: ADR-004 D1–D8; line numbers = current HEAD after 0011/0012 (dump ~+11, main ~+14 from the pre-apply pins in the Mapper line-check).

## 1. What the driver writes today
`TimeLoop` calls `FDS_SETUP(MODE=3)` after every step (patch 0006): `UPDATE_GLOBAL_OUTPUTS` for the meshes of the rank, `EXCHANGE_GLOBAL_OUTPUTS`, `UPDATE_CONTROLS`, optional per-mesh `DUMP_MESH_OUTPUTS` when `OUT_MESH_DUMPS` (0011, default off), `DUMP_GLOBAL_OUTPUTS`
(`_devc.csv`, `_hrr.csv`, `_mass.csv`, ...), `WRITE_STRINGS`, `WRITE_DIAGNOSTICS`. The set-up (`.smv`) is FDS's own. With `FDSTL_MESH_DUMPS` the driver also writes slice/boundary/… files in uniform mode (0011 check). Everything below concerns level 0 first (uniform mode, ADR-004 D1) and then the refined levels.

## 2. Inventory of the routines that loop over meshes or are written per mesh
Kind: S = set-up, P = per output time per mesh, G = global, R = restart.
**OutputPlan membership** = §2 rows that feed the output-mesh path; **REFUSE** = not in OutputPlan (driver refuses the physics, or Mapper closed as refuse).

| Routine | File:line (HEAD) | Kind | Mesh dependence | ADR-004 / OutputPlan | Driver action |
|---|---|---|---|---|---|
| `ASSIGN_FILE_NAMES` | dump.f90:303 (loop 464) | S | `FN_*`/`LU_*` sized by `NMESHES_OUT` | tables by output count | 0012 APPLIED |
| `WRITE_SMOKEVIEW_FILE` | dump.f90:1777 | S | `GRID/PDIM/TRNX/OBST` per mesh; one global `IBLANK` | rank 0, output-mesh order (D1) | 0012 APPLIED |
| `ADD_EXTERIOR_VENTS` | dump.f90:2842 | S | exterior `VENT`s per mesh | follows output mesh set | with 0012 expand (later) |
| `WRITE_STL_FILE` | dump.f90:1575 | S | geometry | unchanged (level 0) | none |
| `INITIALIZE_DIAGNOSTIC_FILE` | dump.f90:3037 | S | `.out` header per mesh | level-0 diagnostics | none |
| `EXCHANGE_NSLICE_INFO`, `EXCHANGE_NPATCH_INFO`, `EXCHANGE_NOBST_INFO` | main.f90:4920, 5055, 4842 | S | MPI counts for **VTK** | output meshes | with 0012 expand |
| `INITIALIZE_OUTPUT_CLOCKS`, `SET_OUTPUT_CLOCK` | func.f90:432, 496 | S | counters on rank mesh range | output meshes | with 0012 expand |
| `INITIALIZE_MESH_DUMPS(NM)`, `INITIALIZE_GLOBAL_DUMPS` | dump.f90:974, 638 | S | opens/headers per-mesh and global | output meshes | with 0012 expand |
| `MPI_INITIALIZATION_CHORES` | main.f90:1351 | S | compute MPI tables | compute meshes | none |
| `DUMP_MESH_OUTPUTS` | dump.f90:77 | P | per mesh `NM` | writers on filled output mesh (D8) | 0011 APPLIED |
| `DUMP_SLCF` (+ `DUMP_SLICE_GEOM`, `DUMP_SLICE_GEOM_DATA`), `DUMP_SLCF_VTK` | dump.f90:6757, 6272, 6410, 7334 | P | slice files; `AMR_LEVEL` (D2) to add | output mesh | none for level 0 |
| `DUMP_CFACES_GEOM` | dump.f90:6337 | P | cut-cell faces (`CC_IBM`) | **REFUSE** (CC_IBM refused, TimeLoop.cpp:287); `WRITE_CFACES` dump.f90:6321 not listed for OutputPlan | none |
| `DUMP_BNDF` (+ `_VTKHDF`) | dump.f90:11966, 12193 | P | boundary from `WALL` | D5 | with output-mesh builder |
| `DUMP_SMOKE3D`, `GET_SMOKE3D_QQ`, `_VTKHDF` | dump.f90:5063, 5128, 5186 | P | soot/HRRPUV | output mesh | none for level 0 |
| `DUMP_ISOF`; `_uvw_t*` / `_tmp_t*` / `_spec_t*` CSV, `DUMP_MMS` | dump.f90:4871; calls in 77–287 | P | per mesh | output mesh (verification; L0) | none for level 0 |
| `SANDIA_OUT` | turb.f90:2412 (`sandia_out`; PUBLIC; call dump.f90:293) | P | Sandia verification | output mesh (L0) | none for level 0 |
| `DUMP_PART`, `DUMP_PART_VTKHDF` | dump.f90:4507, 4633 | P | particles | **REFUSE**: Mapper refuses **all** Lagrangian at the driver (stronger than “refused at fine boxes”; TimeLoop.cpp:289 `IP[17]`) | none |
| `DUMP_PROF` | dump.f90:11378 | P | one file per profile; acts on `PF%MESH==NM` | level 0 / action **none** (not an output-mesh item) | none |
| `DUMP_RESTART`, `READ_RESTART` | dump.f90:3871, 4058 | R | one restart per mesh | out of ADR-004 | report |
| `UPDATE_DEVICES_1/2`, `DUMP_DEVICES`, `UPDATE_HRR`, `UPDATE_MASS`, `DUMP_HRR`, `DUMP_MASS` | dump.f90:7783, 8349, 11155, 11653, 11889, 11846, 11945 | G | per-mesh then MPI sums | D6: composite uncovered cells, ExactSum | driver `composite_mass_exact` / OutputPlan |
| `WRITE_DIAGNOSTICS`, `EXCHANGE_DIAGNOSTICS` | dump.f90:4315; main.f90:4395 | G | CFL/div extrema | L0 + fine extrema later | with fine stepping |
| `WRITE_STRINGS` | main.f90:4279 (send; `DO NOM=2,NMESHES`) | S/P | `.smv` string exchange | output mesh set | with 0012 expand |
| `INITIALIZE_VTKHDF_FILES` (+ `INITIALIZE_VTKHDF_SLCF` etc.) | vtkf.f90:1327 | S | VTK-HDF init; `NMESHES` loops | output meshes when VTK on | with 0012 + `WITH_HDF5` |
| `DUMP_GEOM`, `DUMP_HVAC`, `DUMP_CONTROLS` | dump.f90:12538, 11535, 11334 | G | global | unchanged | none |

Per-mesh state the writers read: `MESHES(NM)` fields, `WALL`, `CELL`/`IBLANK`, module arrays indexed by `NM` (`*_COUNTER`, `FN_*/LU_*`, `PROCESS(NM)`, `M%N_SLCF`; `N_BNDF`/`N_SMOKE3D` are global). Anything indexed by mesh number must exist for output mesh numbers once `N_OUT_MESHES > NMESHES` is validated.

## 3. Design (0011/0012 APPLIED behind `WITH_AMREX`; Step C in driver)
**Step A (0011, APPLIED):** `FDS_SETUP(MODE=3)` runs the per-mesh dump loop when `OUT_MESH_DUMPS` (default off). Uniform-mode check: `tests/run_mesh_dumps_check.sh`.

**Step B (0012, APPLIED, option (a) first part):** `NMESHES_OUT = MAX(NMESHES,N_OUT_MESHES)`; tables sized accordingly. Default `N_OUT_MESHES=0` → size `NMESHES`. **Do not** expand remaining `NMESHES` loops until `N_OUT_MESHES > NMESHES` is validated.

**Step C (driver, `OutputPlan`):** (1) static layout from the registry — level-0 `&MESH` outside the refinable region, tiles at `L_out` inside, ≤ `max_grid_size`, blocking-aligned; writer ranks = contiguous blocks balanced by cell count via `DistributionMapping(Vector<int>)`; (2) fill: `ParallelCopy` COPY from level 0 + `average_down` from finer (host); `OverrideSync` on face data before face fill; (3) `AMR_LEVEL` MultiFab; (4) bind to Fortran output mesh objects — **gated** on `N_OUT_MESHES > NMESHES` (`bind_ready()` false until then; reuse `FDS_FINE_B_SET_VIEW` / `BUILD_FINE_BOX` device when opened); (5) derived slices from FDS on the output mesh; (6) D6 via `composite_mass_exact` (`exact_sum_hierarchy` / uncovered product).

## 4. Order of work
1. ~~Step A~~ done (0011).
2. ~~Step C layout + fill + writer DM + D6~~ done (`OutputPlan`, `composite_mass_exact`, `test_output_plan`).
3. Rank-count byte identity spike S-A on frozen fields from `OutputPlan::fill` (optional TwoLevelRun `FDSTL_OUTPLAN` hook deferred while Role 3 owns that file).
4. Expand 0012 mesh loops / bind OUT_MESH — **after** `N_OUT_MESHES > NMESHES` validation (`bind_ready` gated false).
5. VTK (FR-074): `WITH_HDF5`, not tested here.
6. Stage 2 background writers (NFR-047): after writer snapshot refactor.

## 5. Tests
- Uniform mode (0011): `run_mesh_dumps_check.sh`.
- Rank-count byte identity (S-A): frozen-field raw dumps from `OutputPlan::fill` on 1/2/4 ranks when run; Smokeview `.sf`/`.bf` wait on bind gate.
- Composite sums: `composite_mass_exact` in `test_output_plan` (1 and 4 ranks via `run_output_plan_check.sh`) equals uncovered FP composite mass within ExactSum tolerance; must **not** equal the sum that includes covered coarse cells (double-count control).

## 6. Questions (Mapper closed; Architect only)
- ~~Legacy Mapper inventory / DUMP_PROF / DUMP_CFACES / VTKHDF~~ — **closed** (FR-072 confirm): DUMP_PROF keep (L0, action none); DUMP_CFACES_GEOM/WRITE_CFACES **REFUSE**; INITIALIZE_VTKHDF_FILES match (0012 + WITH_HDF5); particles **all** Lagrangian refused at driver; membership OK with pin refresh.
- Architect: restart of a hierarchy in Phase 9 scope? (ADR-004 silent.) Option (a) preferred and in use for 0012.
