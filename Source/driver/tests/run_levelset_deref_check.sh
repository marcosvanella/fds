#!/bin/bash
# run_levelset_deref_check.sh <build-dir> [<run-dir>]: derefinement of a level-set run with a moving, non-uniform flow. The CSIRO C064 crude case is run 2 levels (--init-patch, --nofine 1) for 24 steps,
# no tag criterion keeps the fine level (--tag-temp 5000), so the regrid after step 12 removes it. Passes when the DEREFINE line shows relative changes of T, RHO, Z and background pressure
# <= 1e-12 on the leaf field (the covered level-0 cells the removal uncovers must hold the fine average, not their own stale coarse update: regrid_dynamic syncs fine -> coarse
# first), max|div u - D| <= 1e-9, and the run gate passes. Measured with that sync switched off: T 5.1e-4, RHO 6.7e-4, Z up to O(1), mass gate 22x the limit.
set -u
BLD=$1; RUN=${2:-$BLD/run/levelset_deref}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" && cd "$RUN" || exit 1
cp "$HERE/cases/levelset_tag_phi_ground.fds" c.fds
sed -i 's/^&HEAD CHID=[^,]*,/\&HEAD CHID='"'c'"',/' c.fds
sed -i '1a &AMR MAX_LEVEL=0, BLOCKING_FACTOR=2, MAX_GRID_SIZE=12 /' c.fds
OMP_NUM_THREADS=2 "$BLD/fds_amr" c.fds --two-level-run --levels 2 --init-patch --patch 5 9 0 11 --blocking 2 --maxsize 8 --regrid-every 12 --tag-temp 5000 --nofine 1 --steps 24 --outdir . --chid c > run.log 2>&1
if grep -q "RESULT derefinement invariance (1 derefining regrids).*PASS" run.log && grep -q "RESULT gate PASS" run.log; then echo "PASS levelset_deref: $(grep -m1 'DEREFINE after' run.log | cut -c1-200)"; exit 0; fi
echo "FAIL levelset_deref (see $RUN/run.log)"; grep -m2 -E "DEREFINE|RESULT gate" run.log | cut -c1-250; exit 1
