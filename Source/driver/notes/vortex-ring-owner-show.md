# Vortex-ring owner-show (TAG_VORT adaptive)

Primary adaptive demo for `vv-runs/owner-show/`: a self-advecting 3-D vortex ring
on a periodic cube, refined by undivided vorticity magnitude (`TAG_VORT`).

## Cases

| Grid | Driver test path | Pack path |
|------|------------------|-----------|
| 16³ | `Source/driver/tests/cases/vortex_ring_16_l0.fds` (+ `_uvw.csv`) | `vv-runs/owner-show/inputs/` |
| 32³ 1-mesh | `.../vortex_ring_32_l0.fds` | same |
| **32³ 4-mesh (np=4)** | `.../vortex_ring_32_4m.fds` + `_4m_uvw_m{1..4}.csv` | **overnight primary** |
| 32³ 8-mesh (np=8) | `.../vortex_ring_32_8m.fds` + `_8m_uvw_m{1..8}.csv` | pop-os / denser ranks |

Generator: `Source/driver/tests/make_vortex_ring_case.py [N] [outdir] [--4m] [--8m] [--sx S --sy S --sz S]`.

IC: Rosenhead-regularised Biot-Savart circular filament (x–y plane, self-advects in +z).
`FREEZE_VELOCITY` off. No dye / `TAG_SPEC_*`.

### Multi-mesh UVW (important)

FDS `UVW_INIT` with `PERIODIC_TEST=2` reads **each** mesh from the **start** of
`UVW_FILE` (mesh-sized). A global N³ CSV does **not** tile. Multi-mesh cases emit
per-mesh slices and `&CSVF UVWFILE='..._uvw_m1.csv'` so `UVW_RESTART` overwrites
each mesh from its slice after the first read. `UVW_FILE` is the m1 slice (avoids EOF).

## Driver CLI (overnight / owner viz)

```bash
# Box or laptop — np must equal NMESHES (FDS ERROR 112 otherwise)
mpirun -np 4 fds_amr vortex_ring_32_4m.fds --two-level-run \
  --regrid-every 5 --tag-vort auto \
  --maxsize 8 --blocking 4 \
  --plotfile DIR --plot-every 2 \
  --steps 2000 --outdir … --chid vortex_ring
```

### Recommended TAG / grid (tight L1, not fat)

| Knob | Recommended | Notes |
|------|-------------|-------|
| `--tag-vort` | `auto` (= **0.45 × max\|VORT\|_L0**) | Raise further (e.g. fixed 0.5–0.6×max) if L1 still fat |
| `TAG_KEEP` | **0.85** (emitted on TAG_VORT path) | Was 0.8; 0.85–0.9 hugs ring core |
| `--maxsize` | **8** (try 4–8) | Smaller → more / finer L1 boxes |
| `--blocking` | **4** | Must divide maxsize |
| `--regrid-every` | 4–5 | |
| `--plot-every` | **2** | Dense frames |
| R | **2** | `--two-level-run` default |
| np | **= NMESHES** | 4→`_4m`, 8→`_8m`; do not oversubscribe past mesh count |

Plotfile fields: **RHO, U, V, W, OMEGA_Y** (OMEGA_Y = undivided ω_y = ΔzU − ΔxW).
`VORT` (\|ω\|) remains tagging-only via `Fields::ensure("VORT")` — not dumped.

## Still open

- Owner movies: dark-grey mesh + AMR boxes, prefer ω_y on mid-y xz plane.
- Threshold fine-tune on 32³ if auto still fat/thin after first overnight.
