#!/bin/bash
# run_sealed_zones_check.sh <build-dir> [<run-dir>]: two sealed isothermal compartments (cases/sealed_two_zone.fds, 2 &ZONEs) split by a partition (x 0.78-0.81) thinner than a
# level-0 cell, 3 levels forced at step 0 over the partition (--init-patch), derefined by the regrid after step 3 (pressure zones inherited by every level, D7/D8; derefinement
# keeps T, Z, background pressure and divergence, D9). Passes when the DEREFINE line shows relative changes of T, Z, p0 <= 1e-12, max|div u - D| <= 1e-9 on every step, no
# abort, and the run gate passes. Needs the per-cell zone inheritance of fds_fine_level.f90 (PATCH-D7); without it the run aborts (fds_fine_level_create code 3).
set -u
BLD=$1; RUN=${2:-$BLD/run/sealedzones}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" && cd "$RUN" || exit 1
cp "$HERE/cases/sealed_two_zone.fds" .
FDSTL_ZONES=1 OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} "$BLD/fds_amr" sealed_two_zone.fds --two-level-run --levels 3 --init-patch --patch 4 11 0 11 --patch2 6 9 3 8 --tag-temp 5000 --tag-obst-snap 0 --steps 9 --regrid-every 3 --blocking 4 --maxsize 8 --outdir . --chid sz --log-every 1 > run.log 2>&1
bad=$(grep -a "MULTI-LEVEL step" run.log | sed -E 's/.*max\|div u - D\| ([^ ]+).*/\1/' | awk '$1+0>1e-9{n++} END{print n+0}')
if [ "$bad" = 0 ] && grep -q "RESULT derefinement invariance (1 derefining regrids).*PASS" run.log && grep -q "RESULT gate PASS" run.log; then echo "PASS sealed_zones: $(grep -m1 'DEREFINE after' run.log | cut -c1-200)"; exit 0; fi
echo "FAIL sealed_zones (see $RUN/run.log; steps with div>1e-9: $bad)"; grep -a -m3 -E "DEREFINE|RESULT|abort|failed" run.log | cut -c1-250; exit 1
