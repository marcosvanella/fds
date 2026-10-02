# R4 step 4: RegridAmrCore wiring (what exists, what is missing)

## What exists (tested in `test_regrid_core`, 1, 2 and 4 ranks)
- `RegridAmrCore::ErrorEst` runs a user tag function (`set_tag_function`, normally built from `TagOps.H`), force-tags the footprints of static finer `&MESH` boxes, clips tags to `Level::taggable` (once-per-level message, every discarded tag counted in `RegridStats`) and counts the tags.
- AMReX calls `ErrorEst` with `ngrow = 0`, buffers with `N_ERROR_BUF`, coarsens by the blocking factor and then calls the virtual `ManualTagsPlacement`. That hook clips the buffered blocks to the refinable region (counted as `buffer_blocks_discarded`). `MakeBaseGrids` is not virtual, so `init_from_tags` repeats AMReX's loop (`MakeNewGrids` + `MakeNewLevelFromScratch`) until `MAX_LEVEL`.
- `MakeNewLevelFromScratch`, `RemakeLevel`, `ClearLevel` forward to a `LevelListener` (layout) and a `LevelDataTransfer` (data): `fill_initial_level` (t = 0, evaluate the input directly on the new level, repeating tag/make/fill), `fill_new_level` / `fill_remade_level` (later regrids, interpolate from the coarser level, copy old fine data where it overlaps), `hierarchy_done(initial)` (one average-down, ghost refill, derived fields).
- `init_from_tags` (t = 0) and `regrid_dynamic` (between full steps, with `begin_regrid` / `end_regrid` around it) drive these. Box-to-rank assignment can be supplied (`DmFn`); grids and data are identical for a different assignment (FR-015 check, hash compared at 1, 2, 4 ranks).
- Test (mock store with rho*Z, 3 components; 32^3 periodic, 3 levels, a blob moving 1.5 level-1 cells per regrid): composite change per regrid 2e-16 (limit 1e-12), zero clips, the feature stays on the finest level, no cell of the finest level outside the region, hysteresis active. Negative controls: no tag buffer loses the feature, stale coarse data under the fine patches breaks the budget.
- Already in Role 1's `LevelRegistry`: `LevelListener` implementation, `retired_level(l)`, `release_retired(l)`, covered mask.

## Region and nesting
Levels below `MAX_LEVEL` may exceed the refinable region by the proper-nesting cover of the next finer level (`N_PROPER`, snapped to the blocking factor). Only the finest level is strictly inside the region. The "region" check of the test uses that rule.

## Missing (not written, nothing faked)
Role 1:
1. `begin_regrid` / `end_regrid` in `LevelRegistry` (defaults do nothing today; `release_retired` must be called after the transfer).
2. Fine-level FDS mesh objects (per-level `Fields` on any BoxArray at run time) and the "initial fields on level L" entry point (the t = 0 evaluation of the input).
3. An adapter from `LevelRegistry` fields to `LevelDataTransfer`: it must map the FDS index map (faces offset by one, ZZ as mass fraction including passive scalars) to the cell and face operators of `CellTransfer.H` / `FaceTransfer.H`, call `prolong_conserved` for rho and rho*Z (and rebuild rho as the species sum), `prolong_faces_normal_linear` for velocity, and fill coarse ghost layers first.
4. Thermo hooks to recompute TMP and the derived fields on new levels, and `run_divergence_part1` / `max_divergence_error` to print `max|div u - D|` after a regrid.
Role 2: a rebuild entry point for the pressure backend after a hierarchy change.

## Not verified
Level-0 boxes that are not blocking-factor multiples under `MakeNewGrids` nesting; clusters from the grid-efficiency step that exceed a non-box refinable region (only a single box region was tested); behaviour at walls and OBST; GPU runs.
