#!/bin/bash
# run_leaf_devc_check.sh <build-dir> [<run-dir>]: two-level simple_test with gas statistic devices over the whole domain (cases/leaf_devc_two_level.fds).
# Leaf rule (no double counting): the VOLUME statistic is the analytic gas volume 3.6*2.4*2.4 - 0.4*0.4*0.2 = 20.704 m3 in every row, and the VOLUME INTEGRAL of DENSITY equals
# the _mass.csv total of the same time (both leaf). The fine level must exist in the run (more than one level in the log).
set -u
BLD=$1; RUN=${2:-$BLD/run/leafdevc}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" && cd "$RUN" || exit 1
OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} "$BLD/fds_amr" "$HERE/cases/leaf_devc_two_level.fds" --two-level-run --steps 60 > run.log 2>&1
python3 - <<'P' || { echo "FAIL leaf_devc (see $RUN)"; exit 1; }
import csv,sys
d=list(csv.reader(open('leaf_devc_devc.csv')))[2:]; m=list(csv.reader(open('leaf_devc_mass.csv')))[2:]
mt={round(float(r[0]),4):float(r[1]) for r in m}
worst_v=max(abs(float(r[1])-20.704) for r in d); worst_m=max(abs(float(r[2])-mt[round(float(r[0]),4)]) for r in d if round(float(r[0]),4) in mt)
nsum=max(float(r[5]) for r in d)/float(d[0][5])
print(f"rows {len(d)}, max|VOLUME-20.704| {worst_v:.2e}, max|devc mass - mass.csv| {worst_m:.2e}, SUM growth (fine cells present) {nsum:.3f}")
sys.exit(0 if worst_v<1e-9 and worst_m<1e-6 and nsum>1.2 else 1)
P
echo "PASS leaf_devc"
