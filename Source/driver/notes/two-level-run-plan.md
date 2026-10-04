# Two-level end-to-end run of periodic `ns2d_16` (1-to-2): state, gap list, timing status

Status: NOT run. All numbers requested for this milestone (composite mass change, max |div u - D|, level-0 + level-1 CPU time) are therefore absent; nothing below is a physics result.
When the run exists every number it produces is labelled pre-validation until the Intel Build Chief reports a pass for patches 0007-0009 (applied as DRAFT).
This run had no way to message other agents, so that status was not obtained.

## What exists (and is tested, level 0 / unit level)
- Flux hooks on level 0 (`FluxStages`, patch 0009): ADV and DIF read-out, DIF override applied between `apply_flux_divergence` and `run_divergence_part1`, `max_divergence_error` = max |div u - D|
  (5e-15 against max |D| = 0.65 on the periodic Taylor-Green checks, `tests/run_div_error_check.sh`). Role 3's `FluxOverrideOps` builds the face lists; the ordering follows D-061.
- Registry: `make_level/remake_level/clear_level`, `begin_regrid/end_regrid`, `fill_initial_level`, `covered_mask`; `RegistryTransfer` (S11.1) implements Role 3's `LevelDataTransfer`
  (conservative prolongation, conservative average-down, bitwise copy over the overlap); unit test passes at 1 and 4 ranks.
- Draft fine-box builder (`BUILD_FINE_BOX`, patch 0007 `FDS_FINE_B`) and the fine-ready `fds_p_*` wrappers.

## What is missing for the run (blockers, in dependency order)
1. **Composite pressure solve.** `pressure_backend` returns `NotBuilt` for more than one level, for covered cells and for masked cells (Role 2 milestone 2). Without it a real two-level step cannot enforce div u = D on the composite grid.
2. **Level-1 binding in `TimeLoop`.** `TimeLoop` holds a `LevelCtx` for level 0 only; `select(lev>0)` aborts. Needed for level 1: fine box objects from the registry, views on the level-1 fields, a `BcStep`
   with the `ext_ghost` exchange, a `SideData`, fine-box domain-edge and coarse-fine interface walls (for the periodic case only the interface walls), the `INTERPOLATED_MESH` mask and fine-level pressure stages.
3. **Case layout.** `Verification/Adaptive_Mesh_Refinement/ns2d_16_int_1to2_refinement.fds` is 12 coarse ring meshes plus one fine mesh; the registry run needs level 0 as the whole 16 x 1 x 16 periodic box with a
   level-1 patch over the central 8 x 8 coarse cells (Role 3's `tests/cases/ns2d_16_int_1to2.fds` adds `&AMR MAX_LEVEL=1`). A one-mesh level-0 case does not exist yet. The case has constant density and one species, so
   the composite mass change is trivially zero there; a variant with a density or species gradient is needed for a meaningful conservation number.
4. **Validation of 0007-0009** on the oneAPI build (Intel Build Chief): required before any physics claim.

## Feasible step that does not need 1 or 2 (proposed next; not started)
A prescribed-velocity two-level transport test in `driver_unit_tests` (AMReX only): registry with one refined patch, an exactly divergence-free stream-function velocity, density and species with a gradient, subcycled
level-1 advance, fluxes recorded on the coarse-fine faces, `FluxOverrideOps` lists, override applied, `RegistryTransfer::hierarchy_done`. Reports the composite mass and species mass change (should be at round-off
level) and the discrete |div u - D| of the composite grid. It tests the override wiring and the transfer, not the FDS operators on level 1.

## Timing (milestone 4)
- A level-0 + level-1 timing needs items 1 and 2 above; not available.
- Existing level-0 baseline (host CPU, one core, shared machine): `notes/wp2-periodic-baseline.md` (128^3 periodic: 5.8-6.5 s per step for 2.1e6 cells, pressure about one third, driver RSS 2.54 GB).
  Those are not idle-machine numbers (NFR-030).
- Cost items of the AMR layer that can be timed before the binding exists: `RegistryTransfer` fill and average-down per regrid, flux-list build and apply per coarse step. Not timed yet.
- GPU needs: none for the host runs described here. The GPU items are the wall-seam kernels (`notes/wall-seam-design.md`) and the transfer kernels (the Role 3 operators already use `ParallelFor`). Laptop runs go through the project coordinator.
