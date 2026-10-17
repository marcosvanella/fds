#!/bin/bash
# run_levelset_regrid_check.sh <build-dir> [<run-dir>]: 2-level --tag-phi level-set run (CSIRO C064 crude, 60 s, --nofine 1) with a regrid every 2 steps and a plotfile after each
# regrid, so fine boxes are created and the layout moves many times. Passes when (1) the composite mass and species gate passes (the fuel the level-0 stateful vents and the fine
# walls above them add is counted once: leaf walls only), and (2) check_levelset_regrid.py finds on every plot: PHI constant in k on both levels, the fine PHI equal to the level-0
# PHI of the parent cell (fill_phi: PHI is only used by the tag and the plot, both refresh it first, so a per-step refresh is not needed), and layout events happened.
# Not covered (measured, see scratch notes): --nofine 0 (fine k=0 over the vents) fails the gate; box removal inside a level never occurs in this case (tag-phi only grows).
set -u
BLD=$1; RUN=${2:-$BLD/run/levelset_regrid}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" && cd "$RUN" || exit 1
cp "$HERE/cases/levelset_tag_phi_ground.fds" c.fds
sed -i 's/^&HEAD CHID=[^,]*,/\&HEAD CHID='"'c'"',/' c.fds
sed -i '1a &AMR MAX_LEVEL=0, BLOCKING_FACTOR=2, MAX_GRID_SIZE=12 /' c.fds
rm -rf plt
OMP_NUM_THREADS=2 "$BLD/fds_amr" c.fds --two-level-run --levels 2 --init-patch --patch 5 9 0 11 --blocking 2 --maxsize 8 --regrid-every 2 --tag-phi --nofine 1 --plotfile plt --plot-every 2 --outdir . --chid c > run.log 2>&1
grep -aq "gate PASS" run.log || { echo "FAIL levelset_regrid (gate, see $RUN/run.log)"; exit 1; }
python3 "$HERE/check_levelset_regrid.py" "$RUN"
