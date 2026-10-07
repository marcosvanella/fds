# Vortex-ring owner-show (TAG_VORT adaptive)

Primary adaptive demo for `vv-runs/owner-show/`: a self-advecting 3-D vortex ring
refined by undivided vorticity magnitude (`TAG_VORT`).

## Current owner pack (Re÷100 + floor INERT) — active

| Item | Value |
|------|-------|
| Grid | **48³** as **4 meshes** (2×2×1 of 24×24×48), np=4 |
| Case | `Source/driver/tests/cases/vortex_ring_48_4m.fds` (+ `_uvw_m{1..4}.csv`) |
| Also | `vortex_ring_48_8m.fds` (2×2×2 of 24³) |
| Pack outdir | `/workspace/fds-amr/vv-runs/owner-show/vortex_ring_re100_wall/` |
| ν | **VISCOSITY=DIFFUSIVITY=1.0E-2** → **Re_advect≈25** (baseline≈2460 @ μ=1e-4) |
| Γ, cz | **−0.8**, **0.60** (downward −z toward floor) |
| R, a | **R=0.30**, **a=0.20** (a/dx≈9.6 on 48³; Architect) |
| BC | **PBZ=0 INERT** (floor); **PBX/PBY PERIODIC**; **PBZ=L MIRROR** (top). Expect driver `periodic=(110)`. **No OPEN, no WALL**. |
| IC loader | `PERIODIC_TEST=2` + UVW/CSVF (loader only; not a BC claim) |

Generator defaults emit this retarget:

```bash
make_vortex_ring_case.py 48 cases --4m --8m   # defaults: --re-scale 100 --downward --floor-inert
# old periodic baseline: ... --baseline-periodic
```

## Archived PASS (pre-retarget)

| Item | Value |
|------|-------|
| Pack | `vv-runs/owner-show/vortex_ring/` + `archive_re_periodic_pass/` |
| Case | `vortex_ring_32_4m.fds` — PERIODIC all sides, ν=1e-4, Γ=+0.8, cz=0.40 |
| Result | **TWO-LEVEL RUN PASS** steps 2000; plts left in place |

## Cases (paths)

| Grid | Driver test path | Pack path |
|------|------------------|-----------|
| **48³ 4-mesh (np=4)** | `.../vortex_ring_48_4m.fds` | **overnight primary** → `vortex_ring_re100_wall/` |
| 48³ 8-mesh | `.../vortex_ring_48_8m.fds` | optional |
| 32³ 4-mesh (archived) | `.../vortex_ring_32_4m.fds` | `vortex_ring/` (periodic PASS) |
| 16³ / 32³ l0 | `.../vortex_ring_*_l0.fds` | inputs/ |

### Multi-mesh UVW (important)

FDS `UVW_INIT` with `PERIODIC_TEST=2` reads **each** mesh from the **start** of
`UVW_FILE` (mesh-sized). A global N³ CSV does **not** tile. Multi-mesh cases emit
per-mesh slices and `&CSVF` so `UVW_RESTART` overwrites each mesh from its slice.

## Driver CLI (overnight / owner viz)

```bash
mpirun -np 4 fds_amr vortex_ring_48_4m.fds --two-level-run \
  --regrid-every 5 --tag-vort auto \
  --maxsize 8 --blocking 4 \
  --plotfile DIR --plot-every 2 \
  --steps 2000 --outdir … --chid vortex_ring
```

| Knob | Recommended |
|------|-------------|
| `--tag-vort` | `auto` (= **0.45 × max\|VORT\|_L0**) |
| `TAG_KEEP` | **0.85** |
| `--maxsize` / `--blocking` | **8** / **4** |
| `--regrid-every` / `--plot-every` | **5** / **2** |
| np | **= NMESHES** (4→`_4m`) |

Plotfile fields: **RHO, U, V, W, OMEGA_Y**. `VORT` tagging-only.

## Architect notes

- OPEN vents refused for this pack (TimeLoop / inhomogeneous Poisson) — deferred.
- Floor is **INERT** (not WALL). Sides PERIODIC; top **MIRROR** (z non-periodic → `periodic=(110)`).
- Do not scale Γ for Re; keep |Γ|=0.8, **R=0.30**, **a=0.20** (a/dx≈9.6 on 48³; Architect).
