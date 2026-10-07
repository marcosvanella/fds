# Show-and-tell notes (V&V)

Owner: AMR V&V Lead. Short owner-facing notes; not a gate.

## Owner show-and-tell pack

Owner made **adaptive regrid** the primary narrative; static R=2 `--two-level-run`
(`ns2d_16_l0` / `blob_16_l0`) remains the **runnable fallback** in
`vv-runs/owner-show/`. Role 1 tip `b644989940`: mid-run regrid + `--plotfile`
(RHO/U/V/W); pack helper `Source/driver/tests/run_plotfile_pack.sh`. Role 3
TAG_SPEC_* namelist at `6d24385755`; **blob e2e confirm still pending** before
enabling adaptive for real. Stretch stub targets `--regrid-every N --plotfile DIR`.
Still **R=2 only**. **D-089 (Architect):** `FREEZE_VELOCITY` allowed for cold-flow adaptive AMR demos (blob show / P3-B).

## Viz path (show gate)

Critical-path viz is **native AMReX multilevel plotfile** via
`WriteMultiLevelPlotfile` on the static R=2 two-level cases (and mid-regrid dumps):

| Item | Detail |
|---|---|
| Cases | `ns2d_16_l0`, `blob_16_l0` under `--two-level-run` |
| Tip | Role 1 `b644989940` (write landed) |
| Artifact | Plotfile directory: `Header` + `Level_*/Cell_H` (+ siblings); **RHO, U, V, W** |
| CLI | `--plotfile [DIR]` — **default off** |
| Pack helper | `Source/driver/tests/run_plotfile_pack.sh` |
| ParaView smoke | `prototypes/s4_cuda_mass/stage1/plotfile_smoke/` (`check_plotfile.sh`, `paraview_smoke.py`) |

**Concrete plt paths** (Role 1 pack, `check_plotfile.sh` PASS):

| Case | Path |
|---|---|
| ns2d static | `(local project directory)/scratch/role1-plotpack/ns2d_static/plt_ns2d_static_final` |
| blob static | `(local project directory)/scratch/role1-plotpack/blob_static/plt_blob_static_final` |
| ns2d regrid | `(local project directory)/scratch/role1-plotpack/ns2d_regrid/plt_ns2d_regrid_final` (+ `_s5`, `_s10`) |

**FDS VTKHDF is not on the critical path** / not the show gate. A-14 (GNU HDF5 /
FireX VTKHDF build) remains background and does not block the pack.

Related: FR-073 (AMReX plotfiles), FR-074 / A-14 (FireX VTKHDF + HDF5 build),
`docs/amrex/plotfile-smoke.md`.
