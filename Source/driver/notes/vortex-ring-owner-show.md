# Vortex-ring owner-show (TAG_VORT adaptive)

Primary adaptive demo for `vv-runs/owner-show/`: a self-advecting 3-D vortex ring
on a periodic cube, refined by undivided vorticity magnitude (`TAG_VORT`).

## Cases

| Grid | Driver test path | Pack path |
|------|------------------|-----------|
| 16³ | `Source/driver/tests/cases/vortex_ring_16_l0.fds` (+ `_uvw.csv`) | `vv-runs/owner-show/inputs/vortex_ring_16_l0.fds` |
| 32³ | `Source/driver/tests/cases/vortex_ring_32_l0.fds` (+ `_uvw.csv`) | (regenerate with the script) |

Generator: `Source/driver/tests/make_vortex_ring_case.py [N] [outdir]`.

IC: Rosenhead-regularised Biot-Savart circular filament (x–y plane, self-advects in +z).
`FREEZE_VELOCITY` off. No dye / `TAG_SPEC_*`.

## Driver CLI

```bash
fds_amr vortex_ring_16_l0.fds --two-level-run \
  --regrid-every 5 --tag-vort auto \
  --plotfile DIR --plot-every 2 \
  --steps 20 --outdir … --chid vortex_adapt
```

- `--tag-vort auto` → `TAG_VORT = 0.35 * max|VORT|` on level 0 after the first undivided `|ω|` fill.
- `--tag-vort <thr>` → fixed undivided threshold.
- Without `--tag-vort`, `--regrid-every>0` keeps the blob `TAG_SPEC_ID='TRACER'` emit.

TwoLevelRun fills cell-centred undivided `|ω|` into Fields `VORT` (`ensure` only — not in
`kDefault`) before ErrorEst / each mid-regrid, and wires `TagFields.vort`.

## Still open (long movie pack)

- Rebuild tip binary with this Role 1 tip; short smoke on 16³.
- Owner-show `run_show.sh vortex` path + dense plotfile movie (dark-grey mesh + AMR boxes).
- Threshold / regrid-interval tuning for a clean ring track on 16³ and 32³.
