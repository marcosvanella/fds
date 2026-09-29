# A-46: vents on gas-gap mesh faces (non-box level 0)

Follow-up to `docs/adr/drafts/ruling-nonbox-level0.md` (ruling N2): list every `&VENT` that FDS attaches to a gas-gap mesh face in the 40 non-box IN cases, split into `OPEN` and others. Analysis only; no simulation beyond three FDS setup-only (T_END=0) checks on 1 rank.

- Script: `vv-runs/tools/gap_face_vents.py` (re-run with `python3 vv-runs/tools/gap_face_vents.py`)
- Per-vent CSV: `docs/vv/gap_face_vents.csv`, one row per FDS vent entry (one input `&VENT` per mesh and MULT copy) with overlap on a gas-gap face
- Per-case CSV: `docs/vv/gap_face_summary.csv`, all 40 cases, including cases with no gap-face vents

## Scope
Rows of `scope_case_list.csv` with `cls = IN` and `nonbox_domain = 1`: **40 rows**, which matches the ruling. All 40 inputs are under `src/Verification/`. None uses `&CATF`, `DB=` vents or TRNX/Y/Z stretching.

## Method
1. **Meshes / level 0.** `&MESH` records are expanded with `&MULT` using `scope_filter.expand_xb`, the same parser that computed `nonbox_domain` (fill = sum of mesh volumes / bounding-box volume). **All 40 cases have one cell size across all their meshes** (`cell_sizes` column; `n_res = 1` in the scope list). No mesh is finer or coarser than another, so every mesh is level 0 and IR-002's multi-resolution question does not arise here. Two cases (`check_kappa`, `species_props`) are 2-D (one cell in y).
2. **Gas-gap face regions.** For each face of each mesh: skip it if it lies on the bounding-box boundary (exterior). Otherwise take the face rectangle and remove the parts where another mesh holds the cells directly across the face. What remains is the gap-face region. The areas are computed exactly by rectangle algebra on a compressed grid. `gap_face_area` in the summary is this geometric area. It does not subtract OBSTs.
3. **Vents, snapped as FDS does.** Each `&VENT` is processed per mesh with the `READ_VENT` logic (`read.f90:12221-12433`):
   - `PBX/PBY/PBZ` and `MB=` become mesh-relative planes (XB = that mesh's bounds, `read.f90:12228-12252`).
   - `MULT_ID` copies are expanded.
   - The XB is clipped to the mesh, and the indices come from `NINT((x-XS)/DX)` (Fortran round-half-away) clamped to `[0,IBAR]`.
   - An entry is rejected if it lies half a cell or more outside the mesh, collapses to a line after snapping, or is fully blocked by non-removable OBSTs in the first cell layer (the BLOCKED test).
   - An accepted entry with `I1==I2` in `{0, IBAR}` (and likewise for J and K) sits on a mesh face. Its snapped rectangle is intersected with that face's gap region.
   - The CSV records the overlap area, the fraction of the snapped entry that overlaps, a `partial` flag (fraction < 1, meaning the rest of the entry borders a neighbour mesh), and `obst_frac` (the share of the overlap backed by an OBST in the first cell layer).
   - A near-miss pass also flags XB vents whose plane lies within one cell of a gap face that FDS does not attach there. There were none.
4. **Class:** `OPEN`, `MIRROR`, `PERIODIC` or `other`. The SURF_ID is always recorded.

## Assumptions
- Uniform grids, so GINV reduces to `(x-XS)/DX`. This holds for all 40 cases.
- `&HOLE` is ignored in the OBST-coverage check. Only `velocity_bc_test` has holes, and none of them touches the first cell layer behind its gap-face vents.
- OBSTs with `DEVC_ID`/`CTRL_ID`/`REMOVABLE` count as removable for the BLOCKED test, following FDS.
- `DB=` vents lie on the bounding box by construction and are never gap-face vents. None occur in these cases.
- A gap face gets FDS's default surface. 15 cases set `&SURF DEFAULT=T` (listed in the summary `notes`), so their gap walls get that surface, not INERT.

## Counts
| | input `&VENT` records | FDS vent entries (per mesh) |
|---|---|---|
| OPEN | 20 | 31 |
| MIRROR | 5 | 14 |
| other | 18 | 35 |
| all non-OPEN | 23 | 49 |

- 21 of the 40 cases have at least one gap-face vent.
- 11 cases have OPEN vents on gap faces.
- 10 cases have only non-OPEN vents there: `part_attenuation`, `condensation_3`, `particle_anisotropic_radi`, `screen_drag_1`, `screen_drag_2`, `sphere_drag_1`, `ground_vegetation_drag`, `ground_vegetation_load`, `ground_vegetation_radi`, `vegetation_absorb`.
- `bi_dir` has both an OPEN vent and non-OPEN vents (VEL surfaces `1p0` and `10p0` via PBX).
- No PERIODIC vent sits on a gap face.
- 19 cases have gap faces but no vents on them. Their gap walls are plain default-surface walls.

### Cases with OPEN vents on gap faces
| case | vent lines (spec) | entries |
|---|---|---|
| Controls/bi_dir | 23 (MB=XMAX) | 2 |
| Flowfields/velocity_bc_test | 42, 43, 44, 45 (PBZ) | 4 |
| HVAC/qfan_multi | 52 (MB=XMIN), 56 (MB=XMAX) | 2 |
| Heat_Transfer/back_wall_test | 45 (MB=XMIN), 46 (MB=XMAX), 48 (MB=YMAX) | 5 |
| Pressure_Solver/hallways | 15 (XB, mesh 5 XMIN face x=3) | 1 |
| Pyrolysis/shrink_swell | 25 (MB=XMIN), 26 (MB=XMAX) | 10 |
| Radiation/hot_spheres | 8 (MB=XMIN), 9 (MB=XMAX) | 2 |
| Restart/device_restart_a | 29 (XB, mesh 4 XMAX face x=3.2) | 1 |
| Restart/device_restart_b | 30 (same vent) | 1 |
| Restart/device_restart_base_case | 28 (same vent) | 1 |
| WUI/LS4_ember_ignition | 19 (MB=XMIN), 20 (MB=XMAX) | 2 |

## Ruling cross-references
- **N3 UGLMAT comparison cases.**
  - `hallways`: one OPEN vent covers the whole gap face at x=3 (y 3-4, z 3-4) of mesh 5.
  - `device_restart_a`: one OPEN vent, 7.6 m x 2.8 m, on the x=3.2 face of mesh 4.
  - **`simple_duct` has no OPEN vent at all.** Its three vents are two HVAC vents and one LEAK vent at interior planes of mesh 3 (z=1.0, 1.1, 1.55). They attach to obstructions, not mesh faces. In the UGLMAT comparison, `simple_duct` exercises only sealed gap walls, not the OPEN Dirichlet placement.
- **N6 `stairwell`.** No gap-face vent of any kind. Its OPEN vent (line 451, y=-7.0) and its `Extract` vent (line 455, z=31.0) both lie on the bounding-box boundary. On vents, the N6 fallback condition "no gap-face OPEN vent" is met.
- **N2 list** (`device_restart_a/b/base_case`, `hallways`): confirmed. The analysis finds 7 more OPEN cases that N2 does not name.

## Edge cases
- **Half-cell snapping (`device_restart_a/b/base_case`).** The vent is given at x=3.1. FDS snaps it to the mesh-4 XMAX face at x=3.2 (DX=0.4) and snaps z 0-2.6 up to 0-2.8. The setup-only run confirms this: the `.smv` shows indices `24 24 14 33 0 7`, IOR -1.
- **Gap part fully behind an OBST (`back_wall_test` line 48, mesh 1).** The MB=YMAX entry on mesh 1 overlaps the gap face only on a 0.05 m x 1 m strip (x 0.50-0.55). That strip is entirely backed by OBST line 39, so no gas is present there. The rest of the entry borders mesh 4. The same vent's entries on meshes 2 and 3 are real gas-gap OPEN faces. This is the only `partial = 1` entry.
- **Partly obstructed OPEN gap faces.**
  - `velocity_bc_test`: PBZ=2.6 and 7.4 are 20% behind the duct-wall OBSTs, and PBZ=7.6 is 1%.
  - `qfan_multi`: 5%.
  - `back_wall_test` mesh 3 YMAX: 5%.
  - In these areas FDS keeps the solid wall, and only the gas part acts as OPEN.
- **Neighbour-covered parts of MB/PB vents.** MB and PB vents also produce entries on faces that border another mesh, or on bounding-box faces. Only the gap-region overlap is listed. FDS lists these other entries in the `.smv` too (for example `back_wall_test` mesh 4 YMIN), but they do not act on interpolated boundaries.
- **Anisotropic cells.** `qfan_multi` uses 0.05 x 0.1 x 0.1 cells, the same in every mesh. The resolution is uniform across meshes, just not isotropic.

## FDS setup-only cross-check
Setup-only runs (T_END=0, 1 rank, `OMP_NUM_THREADS=1`) in `vv-runs/A-46-gapvents/setup/` covered `hallways`, `device_restart_a` and `back_wall_test`. All three stopped normally ("Set-up only") with no warnings. The per-mesh `VENT` blocks in the `.smv` match the script's mesh, face, IOR and snapped extents for every gap-face entry in these three cases:
- `hallways`: mesh 5, `0 0 0 16 0 16`, surface OPEN.
- `device_restart_a`: see above.
- `back_wall_test`: the MB entries on meshes 1-4, including the partial YMAX entry on mesh 1 (`0 21 20 20 0 20`).

The `.smv` does not show whether a face borders a gap or a neighbour. That comes from the geometric analysis.
