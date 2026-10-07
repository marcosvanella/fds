# Patch 0011 (DRAFT): `main.f90`, per-mesh dump loop in `FDS_SETUP(MODE=3)`

Status: DRAFT, needs validation by the Architect and the Legacy Mapper before it is applied. Apply after 0006 (it edits the `MODE==3` block that 0006 adds). Touches only
`Source/main.f90` (one `USE` list entry, one comment, 14 added lines, all inside the `#ifdef WITH_AMREX` branch of `FDS_SETUP`). Step A of `notes/output-plan.md` (A-63 number 0011).

## FireX re-anchor (merge `8ad685f435`, tip after FireX `1fb5c300`)
`patch -p1 --dry-run` applies clean on the merged tree (no offset). Insertion points:
- `USE` list: main.f90:15 (`USE FDS_AMREX_HOOKS, ONLY: ...`, add `OUT_MESH_DUMPS`)
- Insert: main.f90:131-132, between `CALL UPDATE_CONTROLS(...)` and `CALL DUMP_GLOBAL_OUTPUTS` inside `MODE==3`
- MAIN_LOOP reference dump loop: main.f90:1219-1228 (`N_WRITTEN=0`, five `WROTE_*`, `DUMP_MESH_OUTPUTS` at 1226; was 1216-1225 before the merge)
Sources: `docs/upstream-patches/line-check-firex-1fb5c300.md`, `line-check-output-plan.md`.

## What it does
- Insertion point: the `MODE==3` block, between `CALL UPDATE_CONTROLS(T,DT,CTRL_STOP_STATUS,.FALSE.)` and `CALL DUMP_GLOBAL_OUTPUTS`; `OUT_MESH_DUMPS` is added to the
  `USE FDS_AMREX_HOOKS` list of `FDS_SETUP`.
- When `OUT_MESH_DUMPS` is true it runs the per-mesh dump loop of `MAIN_LOOP` (reference sequence: `N_WRITTEN=0`, the five `WROTE_*` resets, then
  `DUMP_MESH_OUTPUTS(T,DT,NM,.FALSE.)` for `NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX`).
- `OUT_MESH_DUMPS` is a flag of `fds_amrex_hooks.f90` (default false, so the behaviour of 0006 is unchanged), set by `fds_hook_set_mesh_dumps`; the driver calls it when
  `FDSTL_MESH_DUMPS` is set and the FDS writers are available.
- Not replicated: the VTK-HDF padding loop and `INITIALIZE_VTKHDF_FILES` of `MAIN_LOOP`, and `IF (CTRL_STOP_STATUS) STOP_STATUS = CTRL_STOP`.

## Evidence (re-checked after FireX)
- `patch -p1 --dry-run` clean on merge tip `8ad685f435` / current HEAD (no offset).
- Preprocessed `main.f90` with `WITH_AMREX` undefined identical to the unpatched file.
- `USE_AMREX=OFF` bitwise (`check_off_bitwise.sh` on a `git archive` with 0011+0012): `shunn3_32` 1 rank PASS (16 files), `shunn3_4mesh_32` 4 ranks PASS (47 files).
- `USE_AMREX=ON` with the patch applied: `tests/run_mesh_dumps_check.sh` PASS — `shunn3_32` np1: 24 slice/boundary files, sizes match, 18 bitwise equal to the baseline; `shunn3_4mesh_32` np4: 80 files, sizes match (content not gated on multi-mesh; driver fields are not bitwise after a full run).

## Unverified
- VTK-HDF padding loop and `INITIALIZE_VTKHDF_FILES` not replicated.
- Behaviour with `WITH_HDF5` / `VTK_HDF` not covered.


## Unverified
- The mode-3 loop has not been exercised: `tests/run_mesh_dumps_check.sh` (slice and boundary files from `fds_amr --run` with `FDSTL_MESH_DUMPS=1` against the baseline
  `.sf` files, 1 and 4 ranks) was written but not run.
- Behaviour with `WITH_HDF5` and `VTK_HDF` output is not covered.
