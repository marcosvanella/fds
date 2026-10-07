# Three or more levels (D-053)

Extends `--two-level-run` (S14) to `max_level >= 2` (`--levels N` with N = number of levels = max_level+1).

## What is in place

| Piece | Status |
|---|---|
| `LevelRegistry` / `make_level` for lev >= 2 | yes (unit test `test_levels` covers level 2) |
| `TimeLoop::bind_level` / `unbind_level` in order | yes (mesh offset via `fds_fine_level_create`; FINE_LEVEL stacks) |
| `RegistryTransfer::derive` parent stage arrays (D, DS, RSUM, MU, KRES, H, HS) | yes, any level >= 1 from level-1 |
| Nested tagging (`patch` on L0, `patch2` on L>=1) | yes in `TwoLevelRun.cpp` |
| CF ghost hooks on all fine levels | yes (`install_cf_ghost_hooks(loop, finest)`) |
| Composite pressure over all bound levels | uses Role 2 `PressureIface` / `solve_composite` with `p.levels.size() == nl`; Role 2 already tests three levels (`pb_comp_three_levels_*`) |
| Flux overwrite / face average-down | already loops `nl` / finest-first |
| `--exact-zone-sums` | accepted on `--two-level-run` (zone sums remain level-0 only; known gap) |
| `--bind-only` | hierarchy + bind, no advance (gate when a full step is not wanted) |

## CLI

```
fds_amr <case.fds> --two-level-run --levels 3 [--patch i0 i1 k0 k1] [--patch2 i0 i1 k0 k1] \
  [--steps N] [--bind-only] [--exact-zone-sums] [--no-overwrite] ...
```

Defaults: `--patch 4 11 4 11`, `--patch2 6 9 6 9` (nested inside patch), ratio 2, blocking 4.

## Test

`tests/run_three_level_check.sh <build> [work] [steps]`: bind-only PASS, then a short advance. See also `driver_unit_tests` (`test_levels` three-level registry block).

## Role 2 dependency

The driver does not edit `Source/pressure_backend`. Full advance needs the composite MLMG path over three levels (already exercised in Role 2's harness). If a driver advance fails with `NotBuilt` or a composite status other than Ok/NotConverged, that is a Role 2 blocker for the driver path; bind+registry still land.

## Known gaps (unchanged)

- Zone sums are level 0 only.
- `PRESSURE_ITERATIONS` / L1208 residual iteration not applied on the multi-level branch (`notes/composite-pressure-iteration.md`).
- Same-level fine seam FV averaging not generally verified.

## Results (`ns2d_16_l0`, ratio 2, patch 4..11 / patch2 6..9, blocking 4)

| Check | Result |
|---|---|
| `--levels 3 --bind-only` | PASS: registry 3 levels, TimeLoop bound 3; L1 256 cells (1 box), L2 576 cells (4 boxes; tag buffer grew the nest) |
| `--levels 3 --steps 10 --exact-zone-sums` | PASS: mass change 1.57e-15; max\|div u − D\| 5.68e-13 (L0 2.5e-14, L1 2.6e-13, L2 5.7e-13); composite MLMG 25 solves |
| `driver_unit_tests` `test_levels` | PASS (includes make_level(2)) |

Role 2 composite at 3 levels: no blocker (driver path Ok).
