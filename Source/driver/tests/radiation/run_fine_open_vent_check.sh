#!/bin/bash
# Regression: radiation_fvm crashed (segfault / unset-pointer read) with radiation ON when a fine box sits on the OPEN vent or fine boxes live on more ranks than level-0 meshes
# (3-level simple_test, 4 level-0 meshes, static patches over the OPEN vent at x=3.6, 10 steps, 1/2/4 ranks). Causes fixed in radi.f90: VENTS(WC%VENT_INDEX) of a fine box (WALL_LOOP1),
# CELL_ILW(0,.) write of a hidden OBST wall, and RADF_CLOCK(RADF_COUNTER(NM)) read for NM>NMESHES with no RADF (a non-short-circuit .AND.).
# usage: run_fine_open_vent_check.sh FDS_AMR_BINARY OUTDIR [RANKS="1 2 4"]   (~10 s per rank count; a build with -fcheck=bounds,pointer on radi.f90 catches the OOB reads even when they do not crash)
# exit 0 = PASS (every rank count finishes with exit 0, no Fortran runtime error / segfault in run.log), 1 = FAIL
A=$1; D=$2; RANKS=${3:-"1 2 4"}; HERE=$(cd "$(dirname "$0")" && pwd); A=$(readlink -f "$A")
[ -x "$A" ] || { echo "usage: $0 FDS_AMR_BINARY OUTDIR [RANKS]"; exit 2; }
source "$HERE/../env.sh" >/dev/null 2>&1; export OMP_NUM_THREADS=1
fail=0
for NP in $RANKS; do
   W="$D/np$NP"; rm -rf "$W"; mkdir -p "$W"; cd "$W"
   {
   echo "&HEAD CHID='simple_test', TITLE='fine box over the OPEN vent' /"
   k=0; for XB in "0.0,1.8,0.0,1.2" "0.0,1.8,1.2,2.4" "1.8,3.6,0.0,1.2" "1.8,3.6,1.2,2.4"; do
      echo "&MESH IJK=18,12,24, MPI_PROCESS=$((k*NP/4)), XB=$XB,0.0,2.4 /"; k=$((k+1)); done
   cat <<'EOT'
&TIME T_END=2.0 /
&REAC FUEL='PROPANE', SOOT_YIELD=0.01 /
&SURF ID='BURNER', HRRPUA=1000., COLOR='RED' /
&OBST XB=0.0,0.4,1.0,1.4,0.0,0.2 /
&VENT XB=0.0,0.4,1.0,1.4,0.2,0.2, SURF_ID='BURNER' /
&VENT XB=3.6,3.6,0.8,1.6,0.0,2.0, SURF_ID='OPEN' /
&AMR BLOCKING_FACTOR=4 /
&TAIL /
EOT
   } > simple_test.fds
   timeout 600 mpirun --oversubscribe --allow-run-as-root -np $NP "$A" simple_test.fds --two-level-run --levels 3 --steps 10 \
      --patch 32 35 0 7 --patch2 34 35 0 3 --outdir . --chid simple_test > run.log 2>&1
   rc=$?
   if [ $rc -ne 0 ] || grep -a -q -i "Fortran runtime error\|SIGSEGV\|Segmentation\|^Segfault" run.log; then
      echo "np=$NP FAIL (exit $rc): $(grep -a -m1 -i 'Fortran runtime error\|Segfault' run.log)"; fail=1
   else echo "np=$NP PASS"; fi
done
exit $fail
