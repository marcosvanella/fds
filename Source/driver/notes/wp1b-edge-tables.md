# WP1b edge tables and WALL_UVW_INTERP (for the Integration Lead)

## Flat edge-table layout (run; same as stage-1 section 7.13b)
FDS holds derived types `CELL(CELL_INDEX(I,J,K))%EDGE_INDEX(1:12)`, `EDGE(IE)%OMEGA(-2:2)`, `EDGE(IE)%TAU(-2:2)`. The driver flattens them for the dump (and later the WP8 shim) as:

| name | shape | type in the dump file | notes |
|---|---|---|---|
| `CE_EDGE_INDEX` | `(0:IBP1,0:JBP1,0:KBP1,12)` | float64 in `scr_*.bin` (integers cast to REAL; harness writes `.i32`) | 0 where `CELL_INDEX=0` |
| `ED_OMEGA` | `(0:NEDGE,-2:2)` | float64 | `NEDGE = UBOUND(MESHES(NM)%EDGE)` (`EDGE` is allocated `0:N_EDGES_DIM`) |
| `ED_TAU` | `(0:NEDGE,-2:2)` | float64 | same |

Transfer: `fds_k_xfer` names `CE_EDGE_INDEX` / `ED_OMEGA` / `ED_TAU` (`fds_kernels.f90`). Dump: existing `stage_scratch` under tag `c_vfluxE`, right after corrector `VELOCITY_FLUX`, together with `WORK1..WORK6` (stress/vorticity), `DS`, `FVX/FVY/FVZ`. Switch: `FDSTL_STAGES=<icyc>` (same as the other scratch dumps). Companion raw fields for the harness: `FDSTL_STAGEG=<icyc>` → `raw_c_prevflux_{US,VS,WS,MU,RHOS}_b*.bin`.

Older diagnostic `FDSTL_EDGES` / `fds_p_edge_dump` still writes a packed `(NE,10)` OMEGA+TAU blob at the predictor wall-BC point; it is not the 7.13b layout. Prefer `c_vfluxE`.

Wall-table dump (`FDSTL_WALLTAB`, `FDS_P_WALL_TABLE_DUMP`) already lives in the tree (S13.1). Still owed for a full WP1b (not this item): corrector-pass fluxes beyond the brackets, `RHO_D`, `Q/QR`, `DELTA_RHO_ZZ`, previous-step `DEL_RHO_D_DEL_Z`.

## Dump commands
```
source Source/driver/tests/env.sh
# edge tables + WORK/DS/FV* at corrector VELOCITY_FLUX of cycle 1; also raw c_prevflux and static metrics
FDSTL_STAGES=1 FDSTL_STAGEG=1 mpirun -np 1 <build>/fds_amr <case>.fds --run --steps 1 --outdir o --chid t
# files: o/scr_c_vfluxE_{CE_EDGE_INDEX,ED_OMEGA,ED_TAU,WORK1..6,DS,FVX,FVY,FVZ}_b0.bin
#        o/raw_c_prevflux_{US,VS,WS,MU,RHOS,...}_b0.bin , o/scr_static_*_b0.bin
Source/driver/tests/run_edge_tables_dump_check.sh <build-dir>   # smoke: files exist, ranks/bounds match the layout
```
Harness prep (unchanged): `prototypes/s4_cuda_mass/stage1/vflux_edge/prep_vflux.py <outdir> <vfe-dir> <box>`.

## What replaces `WALL_UVW_INTERP` in the two-level / AMR path
FDS `DENSITY` restores the normal face velocity of every `INTERPOLATED_BOUNDARY` wall from `UVW_SAVE` (`mass.f90:424-436`; generated kernel `wall_uvw_interp`; corrector twin `wall_uvw_interp_corr` at 595-608). In the driver that loop still runs at **level 0** inside `DENSITY_PRE_CLIP` (`fds_density_split.f90:86-94` predictor, `:204-212` corrector). `UVW_SAVE` is filled by `fds_p_save_uvw` (`fds_step.f90`), called from `BcStep::after_exchange` (`GhostExchange.cpp`) on MESH_EXCHANGE codes 3 and 6 when `EXTERNAL_GHOSTS_FILLED` is set (because `MATCH_VELOCITY` returns at once).

At a **fine level (level > 0)** a box has no external wall cells (`fds_fine_level.f90` allocates empty `UVW_SAVE` / `WALL`; `NWL=0`), so `wall_uvw_interp` is a no-op. Face velocities that FDS would have restored come from:

| face kind | routine | call site |
|---|---|---|
| coarse-fine interface | `composite_fill_fine_face_ghosts` (Role 3 cf ghost hooks via `install_cf_ghost_hooks`) | `TimeLoop.cpp` `s_exchange` for codes 3 and 6 (`US/VS/WS` or `U/V/W`); also the scalar cf hook inside `BcStep::exchange` |
| domain-edge (periodic / Neumann) | `mirror_domain_edges` | `GhostExchange.cpp` after the exchange on a level > 0 |

Level-0 box-interface walls with `FDSTL_IFACE`: `TimeLoop::iface` → `fds_p_iface_walls` turns them into `NULL_BOUNDARY` (or restores them) around the density/flux stages, so the UVW_SAVE restore does not act on those faces; the single-mesh-split path uses `fds_k_neutralize` to set `UVW_SAVE` = own face velocity for the same reason.
