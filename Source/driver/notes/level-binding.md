# Level binding (S12): TimeLoop::bind_level and the two-level transport test

Status: fine-level physics on level > 0 is allowed on the ifx side (patches 0005-0009 passed there); GNU Debug confirmation is pending, so every result below is labelled
pre-validation. The DT_NEW(*) amendment of patch 0007 (velo.f90) is applied by the Architect, not here.

## API (TimeLoop.H)
- `bind_level(lev)`: builds the FDS mesh objects of the boxes of level `lev` (`fds_fine_level_create`, mesh numbers `NM0 + box`), binds the registry fields as views
  (`fds_fine_b_set_view`, 22 field codes), makes the per-level `BcStep` (`ext_ghost = true`: the fine box has no wall cells, its ghost layer is filled by the exchange and
  by the coarse-fine fill) and the clip scratch, pushes a `LevelCtx`. Levels are bound in order; only the top bound level can be rebound (after a regrid that remade it) or
  unbound (`unbind_level`). `level_bound(lev)` asks. It aborts with a message when the tree lacks patch 0007.
- `fds_fine_level_create` gives a new box the single pressure zone of the level-0 meshes of its rank (a case with several zones is refused, code 3) and copies the zone-wide
  background pressure state (`PBAR`, `PBAR_S`, `R_PBAR`, `D_PBAR_DT`, `D_PBAR_DT_S`, `U_LEAK`; PBAR uniform in K only without STRATIFICATION).
- Stage functions (`stage_density(lev, ...)`, `stage_divergence1(lev)`, `flux_readout_*`, ...) take the level; `advance()` loops over the bound levels but stops when more than one
  is bound, because the composite pressure is not wired into the driver (see below).
- `stage_state(predictor, first_pass)`, `stage_init_divergence()`: set the global FDS state before a stage by hand (used by the check and the transport test).

## Protocol for a new level
A new level has none of the box-own arrays that the previous step left behind. `DEL_RHO_D_DEL_Z` (read by the density stage of the next step as `DEL_RHO_D_DEL_Z__0`),
`D`, `KRES`, `MU` are produced by one corrector-form `DIVERGENCE_PART_1` (`stage_state(false, false)`, then `stage_divergence1(lev)`): it works on RHO and ZZ, the state of the new
time level. The predictor form would read RHOS and ZZS, which a new level does not have (the first version of the check did this and got DEL_RHO_D_DEL_Z = 0 on level 1).

## FluxStages order (D-061, notes/flux-stage-wiring.md of Role 3)
`FluxStages::apply_flux_divergence` does the density stage only (ADV read-out, ADV overrides, density). The DIF read-out is a separate phase, `readout_dif(level)`, that runs
after exchange/boundary/vflux/wall_bc, followed by the DIF overrides and `run_divergence_part1`:

1. ADV read-out, ADV overrides, density, finest level first;
2. exchange, boundary, velocity flux, wall BC;
3. DIF read-out (`readout_dif`), DIF overrides, `run_divergence_part1`, finest level first.

The single-level path is bitwise unchanged (regression below).

## What `--level-bind-check` proves
Level 1 is a ratio-1 copy of level 0 (same boxes, state copied into the registry fields). The first-stage positions run on both levels and every field is compared.
- `ns2d_16_l0` (uniform state), 1 rank: 32 fields equal, including the ADV and DIF flux arrays.
- `dec1` (1 rank) and `dec4_np4` (4 ranks), Shunn MMS with non-uniform species: 31 fields equal; one KNOWN GAP (below).
- Not bitwise, with the reason: FVY (hidden direction, 1e-17 against FVX of order 1) and D after the first DIVERGENCE_PART_1 (2e-15) are rounding noise; FVX/FVZ in the faces
  next to the periodic seam differ by up to 2e-10 in ns2d (level 0 takes edge vorticity/stress from its EDGE objects, the fine box has none and uses the interior stencil).
- KNOWN GAP, open: the predictor `DIVERGENCE_PART_1` of a level with non-uniform species gives a DS that differs from level 0 by up to 0.4 % (probe: `L1CHK_PROBE=1`) although
  U V W RHOS ZZS TMP RSUM KRES MU, PBAR_S, D_PBAR_DT, the pressure zone and the coordinates are bitwise equal. The corrector form (D) agrees to 2e-15. Suspects: a box-own array
  read only by the predictor path (conduction, species advection face values at the cell halo -1 / IBP1+1, which are outside the registry ghost width of 1). Not resolved.

## Gaps of the binding (not bound yet)
Fine-box domain-edge and interface wall cells, the `INTERPOLATED_MESH` mask, the fine-level pressure stages (composite solve of Role 2 is not wired into `solve_poisson_impl`),
the `RegistryTransfer::derive` hook for EOS and `run_divergence_part1`, and a composite D_PBAR_DT over the levels of a zone.

## Single-level regression
`tests/run_level_bind_check.sh <build> <reference fds_amr> 10`: `ns2d_16_l0` (1 rank), `dec1` (1 rank), `dec4_np4` (4 ranks), 10 steps with `--exact-zone-sums`: all final-field
files and the step log bitwise equal to the executable built before bind_level (Role 3's `bd/fds_amr`).

## Two-level prescribed-velocity transport test (tests/two_level_transport.H, test_units)
See README S12.2. Numbers (pre-validation): composite rho change 2.2e-15 (1 rank) / 4.4e-16 (4 ranks), rho*Z1 2.2e-15 / 6.3e-16, rho*Z2 2.2e-15 / 3.4e-16, across two regrids;
overwrite OFF control: 7.2e-4 / 2.2e-3 / 2.3e-3. Interface velocity equals the mean of the fine faces exactly. Species average-down of rho*Z: 1.1e-16; a linear Z average would give 1.6e-2
(the test fails for it).
