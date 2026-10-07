# Owner show-and-tell (V&V)

Architect pivot: **3-D vortex ring** is the pack **primary** (adaptive regrid).
Tracer blob adaptive is **interim regression only** (held — not final primary).
Frozen **static R=2** (`--regrid-every 0`) remains the fallback. On-disk pack (not a
git repo):

`vv-runs/owner-show/` — `README.md`, `RESULT.md`, `run_show.sh`, `inputs/`.

- **Primary (PASS):** `./run_show.sh` / `adaptive` / `vortex` →
  `inputs/vortex_ring_16_l0.fds` + `vortex_ring_16_uvw.csv`; tip **`50c4f0588f`**
  (`--tag-vort auto --regrid-every 5 --plot-every 2`). Binary:
  `vv-runs/owner-show/build/ompi_gnu_release_tip/` → `…_tip_50c4f0588f/fds_amr`.
  Prefer **FREEZE_VELOCITY off**. No R4 combustion / `TAG_HRRPUV` / species item-7.
- **Interim:** `./run_show.sh blob` → Role 3 tip `370c4412cf` TAG_SPEC TRACER/0.05/0.8,
  `--regrid-every 5 --plotfile --plot-every 2`. Prior PASS: `blob_adapt_regrid5/`.
- **Fallback:** `./run_show.sh static` → `ns2d_16_l0` + `blob_16_l0`. PASS-class tip: `c0d3573c7d`.
- Ancestors: D-090 `e3309bec16`; Role 3 thr×1000 `a8a5b37994`; Role 1 plotfile `b644989940`;
  D-089 FREEZE lift `8c8381b946` (allowed; vortex prefers OFF).
- Pack helper: `Source/driver/tests/run_plotfile_pack.sh`.

**Viz (show gate):** native AMReX multilevel plotfile (`WriteMultiLevelPlotfile`) —
`Header` + `Level_*/Cell_H`; fields **RHO, U, V, W** (+ OMEGA_Y on vortex tip dumps);
`--plotfile` / `--plot-every`. Integration smoke: `prototypes/s4_cuda_mass/stage1/plotfile_smoke/`.
Movies: dark-grey mesh + AMR boxes (`make_movies.sh` / `render_vel_boxes` / `frames/render_omegay_xz.py` → `vv-runs/owner-show/movies/`).

**Primary vortex plt (tip `50c4f0588f`, PASS):**

- `vv-runs/owner-show/out/vortex_adapt/plt_vortex_adapt_{s2…s20,final}`
- Movies: `vv-runs/owner-show/movies/vortex_adapt_omegay_xz_boxes.mp4` / `.gif` (preferred mid-y xz ω_y); also `vortex_adapt_vel_boxes.mp4` / `.gif`
- Evidence: `vv-runs/owner-show/RESULT.md`, `logs/run_vortex_e2e.log`

**Concrete plt (Role 1 pack, tip `b644989940`):**

- `(local project directory)/scratch/role1-plotpack/ns2d_static/plt_ns2d_static_final`
- `(local project directory)/scratch/role1-plotpack/blob_static/plt_blob_static_final`
- `(local project directory)/scratch/role1-plotpack/ns2d_regrid/plt_ns2d_regrid_final`

**Interim adaptive blob dumps:** sparse `vv-runs/owner-show/blob_adapt_regrid5/plt_blob_adapt_*`; dense `vv-runs/owner-show/adaptive_blob/plt_adaptive_blob_*` (~25 dumps; `--steps 40 --plot-every 2 --regrid-every 5`). Dense velocity movie: `movies/blob_adapt_dense_vel_boxes.mp4` / `.gif`.

**FDS VTKHDF is not the show gate** (see `vv/show_and_tell_notes.md`).
