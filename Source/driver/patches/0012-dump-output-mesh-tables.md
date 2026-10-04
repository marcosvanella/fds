# Patch 0012 (DRAFT): `dump.f90`, output-mesh tables sized by the output mesh count

Status: DRAFT, needs validation by the Build Chiefs and the Architect before it is applied. Independent of 0011. Touches only `Source/dump.f90` (a file-scope macro,
a `USE` line under `#ifdef WITH_AMREX`, and `NMESHES` replaced by `NMESHES_OUT` in the table allocations). Step B of `notes/output-plan.md`, option (a), first part only.

## What it does
- `#define NMESHES_OUT MAX(NMESHES,N_OUT_MESHES)` under `WITH_AMREX`, `#define NMESHES_OUT NMESHES` otherwise. `N_OUT_MESHES` is a variable of `fds_amrex_hooks.f90`
  (default 0, so the tables keep the size `NMESHES` unless the driver asks), set by `fds_hook_set_out_meshes` before `fds_setup(0)` (test switch `FDSTL_OUT_MESHES` in `main.cpp`).
- Insertion points in `dump.f90`: the file-scope macro before `MODULE DUMP`; the `USE` line after `USE COMP_FUNCTIONS`; the allocations of `ASSIGN_FILE_NAMES` (`FN_*`/`LU_*`
  from `FN_XYZ` to the VTK file names, 35 lines); the `XLEVEL/YLEVEL/ZLEVEL` allocations in `WRITE_SMOKEVIEW_FILE`.

## Evidence
- `patch -p1 --dry-run` passes on the committed tree.
- Preprocessed `dump.f90` with `WITH_AMREX` undefined (with and without `WITH_HDF5`) is identical to the unpatched file apart from comments and blank lines.
- `USE_AMREX=OFF`, `git archive` of HEAD with 0011 and 0012 applied, `tests/check_off_bitwise.sh`: `shunn3_32` 1 rank PASS (16 files), `shunn3_4mesh_32` 4 ranks PASS (47 files).
- `USE_AMREX=ON`: the patched tree compiled and linked.

## Unverified / not in this patch
- No run with `N_OUT_MESHES` greater than `NMESHES` has been made.
- The loops over the output meshes still run over the compute meshes: `ASSIGN_FILE_NAMES` (mesh loop), `WRITE_SMOKEVIEW_FILE` loops, `ADD_EXTERIOR_VENTS`, `WRITE_STRINGS`,
  the `EXCHANGE_N*_INFO` routines (VTK only), `INITIALIZE_MESH_DUMPS` and the output-clock counters (`func.f90`, sized by the rank's mesh range). They need an output-mesh to
  rank map and are not done here.
- The edit of `notes/output-plan.md` (numbering 0011/0012, the description errors listed in the Legacy Mapper's line check) is not done.
