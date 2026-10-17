#!/bin/bash
# Regression: 3-level simple_test with radiation ON and dynamic regrid (1 rank). Before the fix the run segfaulted in rad_src_fill (RESHAPE of the saved wall
# incident flux RAD_SRC(NM)%WQ, whose shape belongs to the box NM of the previous hierarchy) at the first radiation source fill after the regrid of step 25.
# usage: run_regrid_l3_radiation_check.sh FDS_AMR_BINARY OUTDIR [STEPS=32]    (~2-3 min)
# exit 0 = PASS (run exits 0, finishes the steps, a regrid happened, no Fortran runtime error / segfault in run.log), 1 = FAIL
A=$(readlink -f "$1"); D=$2; N=${3:-32}; HERE=$(cd "$(dirname "$0")" && pwd)
[ -x "$A" ] || { echo "usage: $0 FDS_AMR_BINARY OUTDIR [STEPS]"; exit 2; }
source "$HERE/../env.sh" >/dev/null 2>&1; export OMP_NUM_THREADS=1
rm -rf "$D"; mkdir -p "$D"; cd "$D"
cat > in.fds <<'EOT'
&HEAD CHID='in', TITLE='3-level simple_test, radiation on, regrid' /
&MESH IJK=36,24,24, XB=0.0,3.6,0.0,2.4,0.0,2.4 /
&TIME T_END=5.0 /
&RADI RADIATION=.TRUE. /
&REAC FUEL='PROPANE', SOOT_YIELD=0.01 /
&SURF ID='BURNER', HRRPUA=1000., COLOR='RED' /
&OBST XB=0.0,0.4,1.0,1.4,0.0,0.2 /
&VENT XB=0.0,0.4,1.0,1.4,0.2,0.2, SURF_ID='BURNER' /
&VENT XB=3.6,3.6,0.8,1.6,0.0,2.0, SURF_ID='OPEN' /
&AMR BLOCKING_FACTOR=2 /
&TAIL /
EOT
"$A" in.fds --two-level-run --levels 3 --regrid-every 5 --tag-vort auto --steps "$N" --outdir . --chid in > run.log 2>&1; rc=$?
NR=$(grep -ac "regrid_dynamic" run.log)
if [ $rc -eq 0 ] && [ "$NR" -ge 1 ] && ! grep -aqi "fortran runtime error\|segmentation\|SIGSEGV" run.log; then echo "PASS (rc=0, $NR regrids, $N steps)"; exit 0; fi
echo "FAIL (rc=$rc, regrids=$NR)"; exit 1
