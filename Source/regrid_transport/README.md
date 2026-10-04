# regrid_transport (Role 3): regrid and multi-level transport, Phase 3

Stub (first commit by the Chief Architect, 2026-10-02). Milestones R0 to R6 are in `notes/plan.md`.

| path | content |
|---|---|
| `OWNERS.md` | ownership and edit rules |
| `notes/` | plan and design notes |
| `tests/` | test scripts, run at 1 and 4 ranks, one thread |
| `CMakeLists.txt` | library `fds_regrid_transport_core` (parser) behind the interface target `fds_regrid_transport`; CTest registration; not yet added to the FDS build |
| `AmrInput.H/.cpp` | `&AMR` and `&AMR_REGION` parser and range checks (R0), no AMReX needed |
| `IBox.H`, `Hierarchy.H/.cpp` | IR-002 grouping of `&MESH` by cell size, mesh-pair ratio check, blocking-factor and nesting checks, static hierarchy, dump, refinable region (R1), no AMReX needed |
| `InputConverter.H/.cpp`, `ConvertInputTool.cpp` | D-076 input converter: finer `&MESH` lines to a hierarchy, level-0-only input text for the unchanged FDS set-up (`notes/input-converter.md`; driver patch `notes/driver-patch-main-input-converter.patch`), no AMReX needed |
| `RegridInterface.H` | proposal (declarations only) of the interface to the driver registry, for Role 1 review |
| `notes/flux-override-interface.md` | flux read-out and override input requested from the driver (interface flux overwrite) |
| `tests/` | `test_amr_input.cpp`, `test_hierarchy.cpp`, `mesh_text.H` (test-only mesh reader), `cases/` (mesh lines of two Verification inputs) |
| `notes/readf90-cmake-patch-list.md` | patch list for `read.f90` and top-level CMake (Architect) |

## R2a (AMReX part)
- `RegridAmrCore.H/.cpp`: `AmrCore` subclass; `init_static()` installs the grids of `Hierarchy` and calls `LevelListener::make_level` (level 0 is adopted by the driver, so it is not notified by default).
- `LevelOps.H/.cpp`: `average_down_cells` (fine to coarse volume average of covered cells), `fill_cf_ghosts_pc` (coarse-fine ghost cells, piecewise-constant, layers 1 and 2 hold the same coarse value).
- `DriverAdapter.H/.cpp`: `make_cf_ghost_hook(LevelRegistry&)` for `TimeLoop::set_cf_ghost_hook`, `average_down_registry`. Cell-centred scalars only; H/HS and face velocities are skipped.
- `tests/test_amrcore.cpp` (ctest `regrid_transport_amrcore`, also run with `mpirun -np 2/4`). The CMake part is built when AMReX is found; the installed AMReX needs Fortran enabled and a HYPRE prefix
  (`-DHYPRE_ROOT=<prefix> -DCMAKE_PREFIX_PATH=<prefix>`), as for the driver build.


## Ghost rules, wiring, checks
- `LevelOps`/`DriverAdapter`: FDS coarse-fine ghost rules (see `notes/fr016-ghost-check.md`), exposed as the per-level hook `make_cf_ghost_hook`.
- `TimeLoopWiring`: `install_cf_ghost_hooks(loop, finest)` and `average_down_hierarchy(loop)` for the driver's TimeLoop (compile-checked as `fds_rt_wiring_check`).
- `FluxOverrideOps`: face lists and area-sum of fine fluxes onto coarse faces (`notes/flux-override-interface.md`); tested without the flux hooks.
- `tests/fr016_ghost_check.cpp`: baseline comparison, enabled with `-DRT_FR016_DUMP=<prefix>`.
