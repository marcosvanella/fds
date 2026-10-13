#!/usr/bin/env bash
# Pyramid (vortex ring over a 3-step OBST pyramid), 48^3, 4 meshes, level 0 only (P0 baseline).
. "$(dirname "$0")/_common.sh"
run_case pyramid_48_4m.fds pyramid_l0 --levels 1 --projection OFF --log-every 10 "$@"
