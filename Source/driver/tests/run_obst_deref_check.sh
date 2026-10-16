#!/bin/bash
# run_obst_deref_check.sh <build-dir> [<run-dir>]: derefinement over an OBST that is off the level-0 grid (cases/obst_offgrid_levels.fds). Three levels are forced at step 0
# (--init-patch, nested patches over the OBST), no tag criterion keeps them (--tag-temp 5000, --tag-obst-snap 0), so the regrid after step 3 removes both levels at once and
# the OBST is back on its level-0 snap. Passes when the DEREFINE line shows relative changes of T, Z and background pressure <= 1e-12 (leaf field) and max|div u - D| <= 1e-9,
# and the run gate passes. At 0482a047b8: all changes 0, div 2.4e-14.
set -u
BLD=$1; RUN=${2:-$BLD/run/obstderef}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" && cd "$RUN" || exit 1
cp "$HERE/cases/obst_offgrid_levels.fds" .
OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} "$BLD/fds_amr" obst_offgrid_levels.fds --two-level-run --levels 3 --init-patch --patch 0 11 0 7 --patch2 0 5 0 3 --tag-temp 5000 --tag-obst-snap 0 --steps 9 --regrid-every 3 --blocking 4 --maxsize 8 --outdir . --chid od --log-every 1 > run.log 2>&1
if grep -q "RESULT derefinement invariance (1 derefining regrids).*PASS" run.log && grep -q "RESULT gate PASS" run.log; then echo "PASS obst_deref: $(grep -m1 'DEREFINE after' run.log | cut -c1-200)"; exit 0; fi
echo "FAIL obst_deref (see $RUN/run.log)"; grep -m2 -E "DEREFINE|RESULT gate" run.log | cut -c1-250; exit 1
