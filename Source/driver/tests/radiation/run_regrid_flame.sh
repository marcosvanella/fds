#!/bin/bash
# Radiation across regrids with a flame (simple_test, radiation on, --tag-temp so fine boxes appear and move while it burns).
# usage: run_regrid_flame.sh FDS_AMR_BINARY OUTDIR [LEVELS=2] [T_END=1.0] [TSI=1]  (env TAGT = --tag-temp K, default 330)   (1 rank, niced; needs FDS_RAD_LEAF support in radi.f90)
set -e
A=$1; D=$2; LEV=${3:-2}; TE=${4:-1.0}; TSI=${5:-1}; HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$D"; cd "$D"
cat > in.fds <<EOF
&HEAD CHID='in', TITLE='simple_test radiation regrid' /
&MESH IJK=36,24,24, XB=0.0,3.6,0.0,2.4,0.0,2.4 /
&TIME T_END=$TE /
&REAC FUEL='PROPANE', SOOT_YIELD=0.01 /
&SURF ID='BURNER', HRRPUA=1000. /
&OBST XB=0.0,0.4,1.0,1.4,0.0,0.2 /
&VENT XB=0.0,0.4,1.0,1.4,0.2,0.2, SURF_ID='BURNER' /
&VENT XB=3.6,3.6,0.8,1.6,0.0,2.0, SURF_ID='OPEN' /
&RADI TIME_STEP_INCREMENT=$TSI /
&AMR BLOCKING_FACTOR=4 /
&TAIL /
EOF
export FDS_RAD_LEAF=1 FDS_RAD_REG=0.0,3.6,0.0,2.4,0.0,2.4 OMP_NUM_THREADS=1
PATCH="--patch 0 15 2 17"; [ "$LEV" = 3 ] && PATCH="$PATCH --patch2 0 9 3 12"   # initial boxes round the burner (no tag exists at t=0); later regrids follow --tag-temp
nice "$A" in.fds --two-level-run --levels "$LEV" --steps 1000000 --regrid-every 3 --tag-temp ${TAGT:-330} $PATCH --outdir . --chid in > run.log 2>&1 || echo "driver exit $? (see run.log)"
python3 "$HERE/regrid_flame_check.py" run.log
