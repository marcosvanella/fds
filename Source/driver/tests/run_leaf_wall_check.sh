#!/bin/bash
# run_leaf_wall_check.sh <build-dir> [<run-dir>]: simple_test with the burner OBST on the fine level (hidden fine walls on the burner top and the floor) and wall SURFACE AREA devices.
# Leaf rule: the SURFACE AREA of the whole floor is the analytic 3.6*2.4 - 0.4*0.4 = 8.48 m2 and of the burner top 0.16 m2 in every row (parent walls replaced by the hidden fine walls,
# no double count, no gap), while the convective-heat-flux statistics of the fine walls enter (log shows the fine level).
set -u
BLD=$1; RUN=${2:-$BLD/run/leafwall}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" && cd "$RUN" || exit 1
OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} "$BLD/fds_amr" "$HERE/cases/leaf_wall_two_level.fds" --two-level-run --tag-obst 1 --nofine -1 --patch 0 7 0 7 --steps 60 > run.log 2>&1
python3 - <<'P' || { echo "FAIL leaf_wall (see $RUN)"; exit 1; }
import csv,sys
rows=list(csv.reader(open('leaf_wall_devc.csv'))); h=rows[1]; d=rows[2:]
a=h.index('A_FLOOR'); b=h.index('A_BURN')
wf=max(abs(float(r[a])-8.48) for r in d); wb=max(abs(float(r[b])-0.16) for r in d)
print(f"rows {len(d)}, max|A_FLOOR-8.48| {wf:.2e}, max|A_BURN-0.16| {wb:.2e}")
sys.exit(0 if wf<1e-12 and wb<1e-12 and len(d)>=3 else 1)
P
grep -q "levels 0..1" run.log || { echo "FAIL leaf_wall: no fine level"; exit 1; }
echo "PASS leaf_wall"
