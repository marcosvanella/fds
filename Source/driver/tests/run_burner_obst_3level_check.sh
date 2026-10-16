#!/bin/bash
# run_burner_obst_3level_check.sh <build-dir> [<run-dir>]: 3-level MASS_FLUX burner on an OBST (cases/burner_obst_3level.fds = scratch/m1/b3 simple_test) with L2 touching the floor next to the burner
# (--patch2 0 7 0 3): the L1 box touches the sides of the burner OBST (coarse SOLID cells across its outline). 15 steps, ~1 min. Passes when the run gate passes (composite mass and species balance
# vs the counted burner inflow) and no fine-box solid wall carries a velocity. Fails at 960e69850c (mass residual 1e-5 of initial, the high-side faces of a nodal U/V/W box were taken a cell
# too far out in composite_zero_solid_faces, so 8..32 faces next to the OBST stayed open) and passes with the enclosedCells fix.
set -u
BLD=$1; RUN=${2:-$BLD/run/burnerobst3}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" && cd "$RUN" || exit 1
cp "$HERE/cases/burner_obst_3level.fds" simple_test.fds
OMP_NUM_THREADS=${OMP_NUM_THREADS:-2} "$BLD/fds_amr" simple_test.fds --two-level-run --levels 3 --patch 0 7 0 7 --patch2 0 7 0 3 --steps 15 --blocking 4 --maxsize 8 --outdir . --chid c > run.log 2>&1
if grep -q "RESULT gate PASS" run.log; then echo "PASS burner_obst_3level: $(grep -a -m1 'RESULT steps' run.log | sed -E 's/.*(balance residual[^(]*).*/\1/')"; exit 0; fi
echo "FAIL burner_obst_3level (see $RUN/run.log)"; grep -a -m3 "RESULT" run.log | cut -c1-300; exit 1
