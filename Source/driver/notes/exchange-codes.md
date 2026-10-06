# MESH_EXCHANGE codes after the FireX merge (driver mapping)

Source: FireX `1fb5c300fa` merged at `8ad685f435`; `docs/next-work.md` r6 Role 1 item 1; `docs/inventory/upstream_merge_reports/2026-10-06.md`.

## Mapping the driver uses (`ExchangeFields.H`)

| Code | When (FDS) | Fields the driver fills | Path |
|---|---|---|---|
| 1 | after predictor density | RHOS, ZZS, MU, KRES, D, TMP, RSUM | driver ghost fill + `BcStep::after_exchange` |
| 4 | after corrector density | RHO, ZZ, MU, KRES, DS, TMP, RSUM | driver |
| 3 | after predictor pressure | US, VS, WS, HS | driver |
| 6 | after corrector (and init) | U, V, W, H | driver (velocity/pressure ghosts only) |
| 5 | inside pressure iteration | FVX, FVY, FVZ, H or HS | driver (`exchange_om` for OMESH) |
| 66 | `INITIALIZE_BACK_WALL_EXCHANGE` (once) | full `PACK_WALL` / `PACK_THIN_WALL` | **FDS host only** |
| 67 | every step after `MESH_EXCHANGE(6)`; HT3D | `PACK_WALL_LITE` | **FDS host only** |

After FireX, code 6 no longer carries back-wall data. That moved to 66 (full pack, set-up) and 67 (lite pack, every step). The driver's code 6 therefore fills only U, V, W, H ghosts, which was already the mapping and needs no field change.

## Rule for 66 and 67

Never route 66 or 67 through the driver ghost fill (`ghost_exchange`, `BcStep::exchange`, `BcStep::after_exchange`, `TimeLoop` stage exchange). `exchange_fields(66)` and `exchange_fields(67)` return empty; `after_exchange` returns immediately for any code other than 1, 3, 4, 6. The time loop never calls exchange 66 or 67. FDS continues to run those codes on the host (`MAIN_LOOP`, `INITIALIZE_BACK_WALL_EXCHANGE`).

## Test

`test_exchange_codes` in `tests/test_units.cpp`: code 6 returns exactly {U,V,W,H}; codes 66 and 67 return empty; `TimeLoop.cpp` has no call site of exchange(66) or exchange(67).
