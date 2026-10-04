# Patch 0011 (DRAFT): `main.f90`, per-mesh dump loop in `FDS_SETUP(MODE=3)`

Status: DRAFT, needs validation by the Build Chiefs and the Architect before it is applied. Apply after 0006 (it edits the `MODE==3` block that 0006 adds). Touches only
`Source/main.f90` (one `USE` list entry, one comment, 14 added lines, all inside the `#ifdef WITH_AMREX` branch of `FDS_SETUP`). Step A of `notes/output-plan.md`.

## What it does
- Insertion point: the `MODE==3` block, between `CALL UPDATE_CONTROLS(T,DT,CTRL_STOP_STATUS,.FALSE.)` and `CALL DUMP_GLOBAL_OUTPUTS`; `OUT_MESH_DUMPS` is added to the
  `USE FDS_AMREX_HOOKS` list of `FDS_SETUP`.
- When `OUT_MESH_DUMPS` is true it runs the per-mesh dump loop of `MAIN_LOOP` (reference sequence: `N_WRITTEN=0`, the five `WROTE_*` resets, then
  `DUMP_MESH_OUTPUTS(T,DT,NM,.FALSE.)` for `NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX`).
- `OUT_MESH_DUMPS` is a flag of `fds_amrex_hooks.f90` (default false, so the behaviour of 0006 is unchanged), set by `fds_hook_set_mesh_dumps`; the driver calls it when
  `FDSTL_MESH_DUMPS` is set and the FDS writers are available.
- Not replicated: the VTK-HDF padding loop and `INITIALIZE_VTKHDF_FILES` of `MAIN_LOOP`, and `IF (CTRL_STOP_STATUS) STOP_STATUS = CTRL_STOP`.

## Evidence
- `patch -p1 --dry-run` passes on the committed tree.
- Preprocessed `main.f90` with `WITH_AMREX` undefined (with and without `WITH_HDF5`) is identical to the unpatched file.
- `USE_AMREX=OFF`, `git archive` of HEAD with 0011 and 0012 applied, `tests/check_off_bitwise.sh`: `shunn3_32` 1 rank PASS (16 files bitwise identical to the baseline),
  `shunn3_4mesh_32` 4 ranks PASS (47 files).
- `USE_AMREX=ON`: the patched tree (0011 and 0012) compiled and linked, and the driver ran with it.

## Unverified
- The mode-3 loop has not been exercised: `tests/run_mesh_dumps_check.sh` (slice and boundary files from `fds_amr --run` with `FDSTL_MESH_DUMPS=1` against the baseline
  `.sf` files, 1 and 4 ranks) was written but not run.
- Behaviour with `WITH_HDF5` and `VTK_HDF` output is not covered.
