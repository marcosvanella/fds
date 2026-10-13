#!/usr/bin/env bash
# Vortex ring 48^3, 4 meshes, adaptive: 2 levels, regrid every 5 steps on vorticity, 8^3 boxes, plotfiles every 5 steps.
. "$(dirname "$0")/_common.sh"
run_case vortex_ring_48_4m.fds vr48_amr --levels 2 --regrid-every 5 --tag-vort auto --blocking 8 --maxsize 8 --plotfile plt --plot-every 5 "$@"
