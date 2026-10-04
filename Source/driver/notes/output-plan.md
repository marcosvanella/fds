# Output plan for the driver: FR-072 / ADR-004 (Role 1 side)

Status: plan and inventory, written from the source. The Legacy Mapper's list of mesh-looping output routines was requested through the Chief Architect's queue; it was not available to this
run, so the inventory below was built by scanning `Source/dump.f90`, `main.f90`, `hvac.f90` and `func.f90` for loops over `NMESHES` and for the per-mesh writers. The Legacy Mapper should
check it against theirs (section 6). Basis: ADR-004 D1-D8 (static disjoint output meshes at `L_out`, writer layout independent of the compute layout, stage 1 = synchronous unchanged
`dump.f90` writers), patch 0006 (`FDS_SETUP(MODE=3)`, end-of-step global outputs). Line numbers are those of the current HEAD.

## 1. What the driver writes today
`TimeLoop` calls `FDS_SETUP(MODE=3)` after every step (patch 0006): `UPDATE_GLOBAL_OUTPUTS` for the meshes of the rank, `EXCHANGE_GLOBAL_OUTPUTS`, `UPDATE_CONTROLS`, `DUMP_GLOBAL_OUTPUTS`
(`_devc.csv`, `_hrr.csv`, `_mass.csv`, ...), `WRITE_STRINGS`, `WRITE_DIAGNOSTICS`. The set-up (`WRITE_SMOKEVIEW_FILE`, the `.smv`) is FDS's own. Not called: `DUMP_MESH_OUTPUTS` (slice, boundary, 3-D smoke,
isosurface, particle, Plot3D, profile, VTK) and `DUMP_RESTART`. In a run that has `&SLCF` or `&BNDF` lines the `.smv` therefore lists files that the driver never writes (a Smokeview load error).
Everything below concerns level 0 first (uniform mode, ADR-004 D1 "the output mesh set is exactly the `&MESH` set") and then the refined levels.

## 2. Inventory of the routines that loop over meshes or are written per mesh
Kind: S = set-up only (rank 0 or all ranks, before the loop), P = per output time, per mesh (called as `DUMP_MESH_OUTPUTS(T,DT,NM,.FALSE.)` for the meshes of the rank, main.f90:1223),
G = global (reductions over meshes and ranks), R = restart.

| Routine | File:line | Kind | Mesh dependence | ADR-004 treatment | Driver action |
|---|---|---|---|---|---|
| `ASSIGN_FILE_NAMES` | dump.f90:292 (loop 453) | S | file-name and unit tables `FN_*(.,NMESHES)`, `LU_*` sized by `NMESHES` (405-451) | tables must be sized by the OUTPUT mesh count | patch 0011 (output mesh count) |
| `WRITE_SMOKEVIEW_FILE` | dump.f90:1766 (loops 2066-2202, 2457) | S | writes `GRID/PDIM/TRNX/OBST/IBLANK` per mesh of the rank, rank-order concatenation (2457-2459, 2688-2722) | rank 0 writes the whole `.smv` in output-mesh order (D1) | patch 0011; level 0 unchanged |
| `ADD_EXTERIOR_VENTS` | dump.f90:2831 (2843) | S | adds exterior `VENT`s per mesh for the `.smv` | follows the output mesh set | with 0011 |
| `WRITE_STL_FILE` | dump.f90:1575 (1615, 1627) | S | geometry output, mesh loop | unchanged (geometry, level 0) | none |
| `INITIALIZE_DIAGNOSTIC_FILE` | dump.f90:3037 (3085, 3091, 3782) | S | `.out` header per mesh | level-0 meshes (diagnostics are not output meshes) | none |
| `EXCHANGE_NSLICE_INFO`, `EXCHANGE_NPATCH_INFO`, `EXCHANGE_NOBST_INFO` | main.f90:4872, 5007, 4794 | S | MPI exchange of per-mesh slice, boundary-patch and OBST counts for the `.smv` | output meshes | patch 0011 |
| `INITIALIZE_OUTPUT_CLOCKS`, `SET_OUTPUT_CLOCK` | func.f90:432, 496 | S | `SLCF_COUNTER(NM)`, `BNDF_COUNTER`, ... one counter per mesh | output meshes | patch 0011 |
| `INITIALIZE_MESH_DUMPS(NM)`, `INITIALIZE_GLOBAL_DUMPS` | dump.f90:963, 627 | S | opens and headers the per-mesh files | output meshes | patch 0011 |
| `MPI_INITIALIZATION_CHORES` | main.f90:1334 | S | per-mesh MPI tables (`PROCESS`, counts) | compute meshes; only the output exchange routines above change | none |
| `DUMP_MESH_OUTPUTS` | dump.f90:77 | P | per mesh `NM`: clocks `*_COUNTER(NM)`, calls the writers below | the writers run on a filled output mesh (D8 stage 1) | patch 0010 (call it in mode 3) |
| `DUMP_SLCF` (+ `DUMP_SLICE_GEOM`, `DUMP_SLICE_GEOM_DATA`, `DUMP_CFACES_GEOM`), `DUMP_SLCF_VTK` | dump.f90:6746, 6261, 6399, 6326, 7323 | P | slice files per mesh; `AMR_LEVEL` slice quantity (D2) to add | output mesh | none for level 0 |
| `DUMP_BNDF` (+ `_VTKHDF`) | dump.f90:11950, 12177 | P | boundary files per mesh from `WALL` records | D5: patches of the output mesh from the owning wall faces | with the output-mesh builder |
| `DUMP_SMOKE3D`, `GET_SMOKE3D_QQ`, `_VTKHDF` | dump.f90:5052, 5117, 5175 | P | per mesh, soot/HRRPUV fields | output mesh | none for level 0 |
| `DUMP_ISOF`; the `_uvw_t*`, `_tmp_t*`, `_spec_t*` CSV dumps (`DUMP_UVW`, `DUMP_TMP`, `DUMP_SPEC`), `DUMP_MMS`, `SANDIA_OUT` | dump.f90:4860; calls in 77-287 | P | per mesh, file names carry `NM` | output mesh (verification dumps; level 0 only) | none for level 0 |
| `DUMP_PART`, `DUMP_PART_VTKHDF` | dump.f90:4496, 4622 (loops 4654, 4743) | P | particles per mesh | D4: particle goes to the output mesh that contains it | later (particles refused at fine boxes) |
| `DUMP_PROF` | dump.f90:11362 | P | wall profile files per mesh | level 0 | none |
| `DUMP_RESTART`, `READ_RESTART` | dump.f90:3871 (4017), 4058 (4250) | R | one restart file per mesh, `MESHES(NM)` fields | restart of a hierarchy is a separate item (not in ADR-004) | out of scope, report |
| `UPDATE_DEVICES_1/2`, `DUMP_DEVICES`, `UPDATE_HRR`, `UPDATE_MASS`, `DUMP_HRR`, `DUMP_MASS` | dump.f90:7772, 8338, 11155, 11637, 11873, 11830, 11929 | G | per-mesh sums then MPI sums (`EXCHANGE_GLOBAL_OUTPUTS`) | D6: composite grid, uncovered cells, D-028 exact sums (`ExactSum.H`), ties to lowest global key | driver side, section 4 |
| `WRITE_DIAGNOSTICS`, `EXCHANGE_DIAGNOSTICS` | dump.f90:4315 (4386, 4461), main.f90:4347 | G | CFL/divergence extrema per mesh (`MESHES(NM)%CFL`...) | level 0 plus level > 0 extrema (needs the fine boxes' `CFL`, `DIVMX`) | with the fine-level stepping |
| `DUMP_GEOM`, `DUMP_HVAC`, `DUMP_CONTROLS` | dump.f90:12538, 11535, 11334 | G | global | unchanged | none |
| `INITIALIZE_BACK_WALL_EXCHANGE`, `ZONE_BOUNDARY_EXCHANGE`, `INITIALIZE_PRESSURE_ZONES`, `INITIALIZE_MESH_EXCHANGE_1` | main.f90:2410, 2907, 2679, 2132 | not output | mesh loops of the exchange set-up | not output routines; listed because they size by `NMESHES` | none |

Per-mesh state that the writers read: `MESHES(NM)` fields (`RHO`, `ZZ`, `TMP`, `U`, `V`, `W`, `H`, `D`, `MU`, `Q`, `KRES`, `RSUM`...), the `WALL` records, `CELL`/`IBLANK`, plus module arrays indexed by `NM`:
`*_COUNTER(NM)`, `*_CLOCK`, `FN_*/LU_*(.,NM)`, `PROCESS(NM)`, `M%N_SLCF`, `M%N_BNDF`, `M%N_SMOKE3D`... Anything indexed by the mesh number must exist for the output mesh numbers (section 3).

## 3. Design, with the patches it needs (numbers from 0010; all DRAFT, behind `WITH_AMREX`)
**Step A, level 0, uniform mode (patch 0010, small):** extend `FDS_SETUP(MODE=3)` so that it also runs the per-mesh dump loop of `MAIN_LOOP` (`DUMP_MESH_OUTPUTS(T,DT,NM,.FALSE.)` for `NM=LOWER_MESH_INDEX..UPPER_MESH_INDEX`, between `UPDATE_CONTROLS` and `DUMP_GLOBAL_OUTPUTS`), switched
by a driver flag (`FDS_HOOK_SET_MESH_DUMPS`, default off, so patch 0006's behaviour is unchanged). The fields the writers read are the level-0 arrays, which alias the MultiFabs, so no copy is needed. ADR-004 D1 requires uniform mode to
reproduce the baseline files; the check is the file set, the headers byte for byte, and the values within the T2 tolerance of the fields (the driver's fields are not bitwise those of FDS after full steps, D-054), see section 5.

**Step B, output meshes (patch 0011):** the writers index everything by the mesh number, and the output mesh set differs from the compute mesh set (D1). Two options for the Architect:
(a) a second module array `OUT_MESH(:)` of `MESH_TYPE` (the same device as the fine-level array `FINE_LEVEL` of patch 0007) and an `OUT_POINT_TO(NM)`; the dump routines call `POINT_TO_BOX`-style lookups; the file tables and counters are sized by `N_OUT_MESHES`.
(b) temporarily re-pointing `MESHES`/`NMESHES` at the output meshes while the writers run (ADR-004 D8 stage 1 wording). This touches every writer's use of `NMESHES` in one place (a module variable), but any routine that also reads compute-mesh state (e.g. `PROCESS`, MPI maps) is then wrong during the swap.
Role 1 prefers (a): it does not change what `NMESHES` means for the compute side. Estimated patch size: `dump.f90` table sizes (about 25 allocation lines), `ASSIGN_FILE_NAMES`, `WRITE_SMOKEVIEW_FILE` loops, the three `EXCHANGE_N*_INFO` routines, `INITIALIZE_OUTPUT_CLOCKS`, `INITIALIZE_MESH_DUMPS`. Each edit is a one-token change from `NMESHES` to the output count when `WITH_AMREX`, so the `OFF` source is unchanged.

**Step C, driver side (no patch):** (1) `OutputPlan`: output mesh layout from the registry (level-0 `&MESH` blocks outside the refinable region, tiles at `L_out` inside it, blocking-factor aligned, at most `max_grid_size`), writer ranks as contiguous blocks balanced by cell count, `DistributionMapping(Vector<int>)`; (2) fill: per output time an output MultiFab per output mesh by `ParallelCopy` COPY mode from level 0 (piecewise-constant injection for a coarser source) and `average_down` from finer levels (host arithmetic, exact scale factors; face data after the overrides are in sync); (3) the `AMR_LEVEL` field; (4) bind the filled arrays to the output mesh object (the same `FDS_FINE_B_SET_VIEW` device used for fine boxes; `BUILD_FINE_BOX` builds the metrics and `CELL`/`IBLANK` tables); (5) derived slice quantities come from the unchanged FDS functions on the output mesh, which needs `H`, `D`, `MU`, `KRES`, `Q`, `RSUM` on the output mesh as well as the primitives: either transferred like the primitives or recomputed by the kernels on the output mesh; (6) D6 global outputs from the composite grid with `exact_sum_hierarchy` (exists in `ExactSum.H`) and the MINLOC/MAXLOC tie rule.

## 4. Order of work and what each step needs
1. Step A with its test (host only, this repository, level 0). Result: `.sf`, `.bf`, `.s3d`, `.iso`, `.prt5`, `.q` files of the driver in uniform mode.
2. Output-mesh object builder (shared with the fine boxes) and fill for a static two-level case, writer layout, rank-count byte identity (spike S-A of ADR-004) on 1, 2, 4 ranks: needs the level-1 fields to exist (the two-level run of the queue) and patch 0011.
3. D6 global outputs on the composite grid (driver only, independent of 0011): `UPDATE_HRR`/`UPDATE_MASS` sum over level-0 meshes today; with a covered coarse region they double count the covered cells, so the per-mesh sums must use the covered mask.
4. VTK (FR-074) through the same output meshes: needs `WITH_HDF5`, not tested here.
5. Stage 2 background writers (NFR-047): after the snapshot refactor of the writers (not Role 1 only).

## 5. Tests
- Uniform mode (step A): a level-0 case with `&SLCF`, `&BNDF`, `&ISOF`: file list and sizes equal to an FDS reference run, headers byte-equal, values within the T2 tolerance of the field comparison; `.smv` lists only files that exist.
- Rank-count byte identity of the driver's own output (1, 2, 4 ranks, same BoxArray): `cmp`-identical files (ADR-004 S-A); this part is exact because the driver's fields are bitwise rank independent (decomposition checks).
- Composite sums: HRR and mass on a two-level run with a covered region equal the same run's level-0-only sums plus the fine contribution, and do not count covered cells twice.

## 6. Questions
- Legacy Mapper: compare section 2 with your list; in particular whether `DUMP_PROF`, `DUMP_CFACES_GEOM`/`WRITE_CFACES` (cut cells, refused in AMR mode) and the HDF5 `INITIALIZE_VTKHDF_FILES` loops belong in it.
- Architect: option (a) or (b) of step B; whether restart of a hierarchy is in scope for Phase 9 (ADR-004 does not mention it).
