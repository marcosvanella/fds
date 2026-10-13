#!/usr/bin/env bash
# Vortex ring 48^3, 4 meshes, level 0 only, AMReX MLMG with BiCGStab.
. "$(dirname "$0")/_common.sh"
run_case vortex_ring_48_4m.fds vr48_mlmg --levels 1 --pressure-backend mlmg --krylov bicgstab "$@"
