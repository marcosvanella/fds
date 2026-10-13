#!/usr/bin/env bash
# Vortex ring 48^3, 4 meshes (2x2x1), uniform grid (level 0 only), FFT Poisson solver (the default).
. "$(dirname "$0")/_common.sh"
run_case vortex_ring_48_4m.fds vr48_fft --levels 1 "$@"
