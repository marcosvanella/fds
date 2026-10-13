#!/usr/bin/env bash
# Vortex ring 48^3, 4 meshes, level 0 only, HYPRE PCG with BoomerAMG preconditioner.
. "$(dirname "$0")/_common.sh"
run_case vortex_ring_48_4m.fds vr48_hypre --levels 1 --pressure-backend hypre --krylov cg "$@"
