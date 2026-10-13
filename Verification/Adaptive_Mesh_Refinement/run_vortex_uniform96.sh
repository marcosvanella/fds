#!/usr/bin/env bash
# Vortex ring 96^3 uniform fine grid, 4 meshes (2x2x1, 48x48x96 each), level 0 only, FFT. Pass solver flags to switch backend.
. "$(dirname "$0")/_common.sh"
run_case vortex_ring_96_2x2x1.fds vr96_uniform --levels 1 "$@"
