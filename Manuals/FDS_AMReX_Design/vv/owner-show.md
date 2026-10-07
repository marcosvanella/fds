# Owner show-and-tell (V&V)

Owner made **adaptive regrid** primary. Frozen **static R=2** cold-flow show script
is the **runnable fallback** on disk (not a git repo) at:

`vv-runs/owner-show/` — `README.md`, `run_show.sh`, `inputs/` (symlinks to driver test cases).

- **Fallback (runnable):** `./run_show.sh` → `ns2d_16_l0` + `blob_16_l0` `--two-level-run`. PASS-class tip: `c0d3573c7d`.
- **Primary / stretch:** Role 1 mid-run regrid + `--plotfile` (RHO/U/V/W) at `b644989940`; Role 3 TAG_SPEC_* at `6d24385755`. Stub targets `--regrid-every N --plotfile DIR` once Role 3 confirms blob e2e. `./run_show.sh adaptive` exits 2 (prints flags + tip SHA). Static R=2 fallback unchanged.
- **D-089 (Architect):** `FREEZE_VELOCITY` allowed for cold-flow adaptive AMR demos (blob show / P3-B).
- Pack excludes ratio-4 blob; still **R=2 only**.
- Pack helper: `Source/driver/tests/run_plotfile_pack.sh`.

**Viz (show gate):** native AMReX multilevel plotfile (`WriteMultiLevelPlotfile`) —
artifact directory `Header` + `Level_*/Cell_H`, fields **RHO, U, V, W**; CLI
`--plotfile` default-off. Integration smoke at `prototypes/s4_cuda_mass/stage1/plotfile_smoke/`
(`check_plotfile.sh`, `paraview_smoke.py`; docs `amrex/plotfile-smoke.md`).

**Concrete plt (Role 1 pack, `check_plotfile.sh` PASS, tip `b644989940`):**

- `(local project directory)/scratch/role1-plotpack/ns2d_static/plt_ns2d_static_final`
- `(local project directory)/scratch/role1-plotpack/blob_static/plt_blob_static_final`
- `(local project directory)/scratch/role1-plotpack/ns2d_regrid/plt_ns2d_regrid_final` (also `_s5`, `_s10` mid-regrid dumps)

**FDS VTKHDF is not the show gate** (see `vv/show_and_tell_notes.md`).

## Owner movies (on disk)

Rebuild: `prototypes/s4_cuda_mass/stage1/plotfile_smoke/make_movies.sh`.
Overlays: dark-grey cell mesh at actual `dx` + AMR box outlines.

| File under `vv-runs/owner-show/movies/` | Content |
|---|---|
| `blob_adapt_rho_boxes.mp4` / `.gif` | adaptive blob s5→s10→s15→s20→final; RHO + mesh + moving L1 boxes |
| `blob_adapt_vel_boxes.mp4` / `.gif` | same series; \|u\| + mesh + moving L1 boxes |
| `ns2d_regrid_vel_boxes.mp4` / `.gif` | ns2d regrid \|u\| + mesh + boxes |
| `blob_static_rho_boxes.mp4` / `.gif` | blob_static RHO + mesh + boxes (hold) |

Adaptive inputs: `vv-runs/owner-show/blob_adapt_regrid5/plt_blob_adapt_{s5,s10,s15,s20,final}`.
Index: `vv-runs/owner-show/movies/MOVIES.md`.
