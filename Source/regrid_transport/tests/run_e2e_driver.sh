#!/bin/bash
# Driver-level end-to-end checks of the multi-level transport (prepared; each one SKIPs with its reason until the driver can run a level > 0).
# usage: run_e2e_driver.sh <driver-build> [work]       Exit 0 all run ones pass, 1 a check failed, 77 everything skipped.
# Prerequisite (Role 1): a TimeLoop entry that binds a level > 0 (fine-level FDS boxes through BUILD_FINE_BOX, Fields of the registry, flux arrays) and a way to run
# stages on it. Detected here by the presence of `bind_level` in Source/driver/TimeLoop.H.
# Checks (cases: tests/cases/ns2d_16_int_1to2.fds and variants; run on 1 and 4 ranks):
#   E1 full-domain level 1 equals the uniform-fine run (field by field, round-off)
#   E2 uniform state stays uniform across the interface
#   E3 static two-level ns2d_16_int_1to2: composite mass and species conserved to round-off over a few steps with the overwrite ON; negative control OFF drifts
#   E4 FR-016 ghost check through a real TimeLoop stage (compare the hook's ghosts with the instrumented-FDS dumps, RT_FR016_DUMP)
#   E5 D-059: print the number of covered cells that serve two or more faces; the 4 corner-zone KRES cells of FR-016 are the accepted exception; corner blob tests
set -u
BLD=${1:-}; WORK=${2:-/tmp/rt_e2e_driver}
HERE=$(cd "$(dirname "$0")" && pwd); TL=$HERE/../../driver/TimeLoop.H
skip() { echo "SKIP $1: $2"; }
if [ -z "$BLD" ] || [ ! -x "$BLD/fds_amr" ]; then echo "SKIP all: driver executable fds_amr not found (argument 1)"; exit 77; fi
if ! grep -q "bind_level" "$TL"; then
  for t in "E1 full-domain level 1 vs uniform fine" "E2 uniform state across the interface" "E3 static two-level conservation + negative control" "E4 FR-016 ghost check through a TimeLoop stage" "E5 D-059 shared ghost cells and corner blobs"; do
    skip "$t" "the driver cannot bind a level > 0 yet (no bind_level in Source/driver/TimeLoop.H); waiting on Role 1. The same physics is covered by test_flux_stage on the mock transport."
  done
  exit 77
fi
echo "bind_level found in TimeLoop.H, but the driver-level cases are not written yet: extend this script"; exit 77
