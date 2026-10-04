# Interface flux overwrite: how a stage is wired (Role 3 side, R2b)

Plain words: a *face flux* is the amount of a quantity crossing one cell face per unit area and time. At the edge of a fine patch the coarse faces take the
area-weighted mean of the fine face fluxes that tile them (the "overwrite"); with one time step on all levels this is all that is needed, no flux register.

## Code
- `FluxOverrideOps` builds the lists for one coarse level from the fine level's face fluxes (periodic wrap, shared faces listed in every box that holds them).
- `FluxStageRunner::set_overrides(FluxAccess&, kind)` does it for every level, finest first, and ALWAYS calls `set_flux_override` on every level (the finest gets an
  empty set, which clears the lists of the previous stage). `overwrite=false` clears everything: the negative control and the FDS interface behaviour of the FR-016 baseline.
- `GhostShare` counts the covered coarse cells that serve two or more ghost faces (D-059); the number must be printed at each regrid (the regrid wiring owns that call).
- `tests/mock_transport.H` is a small conservative transport code with the stage structure of one FDS SSPRK2 step, used to test all of the above without the driver's
  fine-level binding (test `test_flux_stage`, 1 and 4 ranks): composite mass conserved to ~1e-17 with the overwrite ON, 1e-3..1e-4 drift with it OFF; uniform state stays
  uniform; three levels; 3-D patch with a lug; patch on periodic edges; full-domain level 1 bitwise equal to the single-level run; override equals the independent mean and
  keeps the species sum zero.

## Stage order the driver loop must follow (per stage, predictor and corrector)
1. `stage_viscosity` (COMPUTE_VISCOSITY + MASS_FINITE_DIFFERENCES) on all levels
2. `flux_readout_adv` on all levels            -> `runner.set_overrides(fs, Adv)`
3. `stage_density` on all levels (consumes the ADV lists)
4. `stage_exchange(code 1 or 4)` and `stage_boundary` on all levels, coarse first (the coarse-fine ghost hook runs here), `stage_velocity_flux`, `stage_wall_bc`
5. `flux_readout_dif` on all levels (= DIVERGENCE_PART_1)   -> `runner.set_overrides(fs, Dif)`
6. `flux_apply_dif` on all levels (re-run where a list exists)
7. average-down of the cell fields (`average_down_hierarchy`)

## Finding for Role 1 (FluxStages as committed)
`FluxStages::apply_flux_divergence(level)` runs `stage_density` and then immediately `flux_readout_dif(level)` (DIVERGENCE_PART_1). In FDS order DIVERGENCE_PART_1 comes
AFTER the exchange of RHOS/ZZS (steps 4), and with several levels, finest-first, a fine level's DIF read-out would run before the coarse level has its updated density and
before the coarse-fine ghost hook has filled the fine ghost cells of that stage. The DIF read-out of a level must therefore not be part of `apply_flux_divergence`; it
belongs to step 5 above. Either split `apply_flux_divergence` (density only) and move the DIF read-out into `run_divergence_part1` / a new phase, or the stage loop calls the
primitive `TimeLoop` entry points directly as listed. Role 3's runner only needs `stage_flux` and `set_flux_override` from `FluxAccess` and works with either.

## Requirements on the velocities
The advective override value is face value x the single interface-face velocity. The coarse interface-face velocity must equal the mean of the fine faces that tile it
(prescribed fields: build the coarse field by average-down of the fine face velocities, or from a stream function sampled at the common nodes, as the mock test does).
