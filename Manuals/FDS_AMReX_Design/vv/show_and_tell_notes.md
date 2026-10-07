# Show-and-tell notes (V&V)

Owner: AMR V&V Lead. Short owner-facing notes; not a gate.

## Owner show-and-tell pack

Architect pivot: **adaptive 3-D vortex ring** is the pack primary under
`vv-runs/owner-show/`. Tracer blob adaptive is **interim regression only** (Role 3
tip `370c4412cf` TAG_SPEC path) — **held**, not the final primary. Static R=2
`--regrid-every 0` remains fallback. Prefer **FREEZE_VELOCITY off** for vortex
(D-089 lift `8c8381b946` allowed if needed).

**TAG_VORT tip landed:** `50c4f0588f` (ancestors D-090 `e3309bec16`, Role 3
`a8a5b37994`). Design: `Source/regrid_transport/notes/r4-tagging-regrid-design.md`.
Do **not** mix R4 combustion / species item-7 tagging into the vortex pack.

| Mode | Command | Status |
|---|---|---|
| Primary | `./run_show.sh` / `vortex` / `adaptive` | **PASS** tip `50c4f0588f` (`--tag-vort auto`) |
| Interim | `./run_show.sh blob` | TAG_SPEC blob smoke; prior PASS `blob_adapt_regrid5/` |
| Fallback | `./run_show.sh static` | ns2d + blob static |

Binary: `vv-runs/owner-show/build/ompi_gnu_release_tip/` → `ompi_gnu_release_tip_50c4f0588f/fds_amr`.
Evidence: `vv-runs/owner-show/RESULT.md`.

## Viz path (show gate)

Critical-path viz is **native AMReX multilevel plotfile** via
`WriteMultiLevelPlotfile` (static, mid-regrid, and denser `--plot-every` dumps):

| Item | Detail |
|---|---|
| Primary case | `vortex_ring_16_l0` under `--two-level-run` + `--tag-vort auto` |
| Interim | `blob_16_l0` adaptive TAG_SPEC |
| Tip (TAG_VORT) | `50c4f0588f` |
| Tip (write base) | Role 1 `b644989940` |
| Artifact | Plotfile directory: `Header` + `Level_*/Cell_H`; **RHO, U, V, W** (+ OMEGA_Y on vortex) |
| CLI | `--plotfile [DIR]` / `--plot-every N` (pack default 2) |
| Pack helper | `Source/driver/tests/run_plotfile_pack.sh` |
| ParaView smoke | `prototypes/s4_cuda_mass/stage1/plotfile_smoke/` |
| Movies | dark-grey mesh + AMR boxes → `vv-runs/owner-show/movies/` |

**Concrete plt paths:**

| Case | Path |
|---|---|
| vortex adaptive (PRIMARY) | `vv-runs/owner-show/out/vortex_adapt/plt_vortex_adapt_*` |
| vortex movie (ω_y preferred) | `vv-runs/owner-show/movies/vortex_adapt_omegay_xz_boxes.mp4` / `.gif` |
| vortex movie (|u| panel) | `vv-runs/owner-show/movies/vortex_adapt_vel_boxes.mp4` / `.gif` |
| ns2d static | `(local project directory)/scratch/role1-plotpack/ns2d_static/plt_ns2d_static_final` |
| blob static | `(local project directory)/scratch/role1-plotpack/blob_static/plt_blob_static_final` |
| ns2d regrid | `(local project directory)/scratch/role1-plotpack/ns2d_regrid/plt_ns2d_regrid_final` (+ `_s5`, `_s10`) |
| blob adaptive (interim) | `vv-runs/owner-show/blob_adapt_regrid5/plt_blob_adapt_*` |

**FDS VTKHDF is not on the critical path** / not the show gate. A-14 remains background.

Related: FR-073 (AMReX plotfiles), FR-074 / A-14 (FireX VTKHDF + HDF5 build),
`docs/amrex/plotfile-smoke.md`.
