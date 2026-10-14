#!/bin/bash
# run_fine_burner_check.sh <build-dir> [<run-dir>]: simple_test two levels with the burner OBST on the fine level (--tag-obst 1 --nofine -1) next to the
# OPEN door face (face_neumann): masked composite MLMG with mixed domain faces, fine OBST wall velocity (MASS_FLUX) restored after the solid-face zeroing.
# Before: abort at step 1 (mixed face built on unmasked MLMG only); then HRR 0 (burner inflow zeroed on the fine solid faces).
# Pass: run gate PASS (fine div <= 1e-9, mass and species balance), backend MLMG, HRR > 0 kW.
set -u
BLD=$1; RUN=${2:-$BLD/run/fineburner}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" && cd "$RUN" || exit 1
sed 's/T_END=5.0/T_END=0.5/' "$HERE/../../../Verification/Adaptive_Mesh_Refinement/simple_test.fds" > simple_test.fds
OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} "$BLD/fds_amr" simple_test.fds --two-level-run --tag-obst 1 --nofine -1 --steps 200 > run.log 2>&1
gate=$(grep -c "RESULT gate PASS" run.log); mlmg=$(grep -c "composite pressure: .*backend MLMG" run.log)
hrr=$(awk '/RESULT composite HRR final/{print $6}' run.log)
if [ "$gate" -eq 1 ] && [ "$mlmg" -eq 1 ] && awk -v h="${hrr:-0}" 'BEGIN{exit !(h + 0 > 1)}'; then echo "PASS fine_burner: gate PASS, MLMG, HRR $hrr kW"; exit 0; fi
echo "FAIL fine_burner: gate PASS $gate, MLMG $mlmg, HRR '${hrr:-}' kW (see $RUN/run.log)"; exit 1
