#!/bin/bash
# W2a (PV-10) particle bed level check, setup only (--bind-only): a gas-only cut-down C064 (x -10..70, y -20..20, z 0..20 at 1 m) with the
# C064 bed XB passed by --bed (the &INIT PART_ID path is the same list once particle cases pass the driver scope, W0/W7).
# (a) 2 levels: L_ins 1 (0.5 m) covers the bed; (b) 3 levels: L_ins 2 (0.25 m); (c) 1 level: no bed line (no-op);
# (d) a heat-flux VENT inside the bed kept on level 0 (D-010, default --nofine): the bed level cannot cover the bed, the run stops;
# (e) the C064 IGNITION burner (HRRPUA+RAMP, TMP_FRONT, SPREAD_RATE/XYZ; D-010 exempt) inside the bed: 2 levels with the default --nofine
#     build the bed level and run (rc 0);
# (f) like-dx (PV-13): a small burner case whose level 1 (0.5 m, over the --bed) covers the whole burner, against stock 1 level at 0.5 m: BURNING RATE and HRRPUA
#     surface integrals and the HRR-file MLR_PROPANE agree to 1e-12 relative at every output (per-face SPREAD_RATE timing and AREA_ADJUST; leaf faces only: a counted
#     covered level-0 face or a skipped leaf face breaks it). The 1 m 1-level run differs (asserted); its burner-MLR L1 difference to 0.5 m is printed as a convergence number.
# usage: run_bed_level_check.sh <build-dir (absolute)> [work-dir]
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
B=$1; W=${2:-$PWD/bedlev}; rc=0; rm -rf "$W"; mkdir -p "$W"
cat > "$W/cut.fds" <<'FDS'
&HEAD CHID='cut', TITLE='cut-down C064 domain, gas only' /
&AMR MAX_LEVEL=0, BLOCKING_FACTOR=2, MAX_GRID_SIZE=32 /
&MESH IJK=80,40,20, XB=-10.0,70.0,-20.0,20.0,0.0,20.0 /
&TIME T_END=1. /
&VENT MB='XMIN', SURF_ID='OPEN' /
&VENT MB='XMAX', SURF_ID='OPEN' /
&VENT MB='YMIN', SURF_ID='OPEN' /
&VENT MB='YMAX', SURF_ID='OPEN' /
&VENT MB='ZMAX', SURF_ID='OPEN' /
&TAIL /
FDS
sed -e 's/^&TAIL/\&SURF ID='"'B'"', NET_HEAT_FLUX=50. \/\n\&VENT XB=10.,20.,-5.,5.,0.,0., SURF_ID='"'B'"' \/\n\&TAIL/' "$W/cut.fds" > "$W/cutv.fds"
run() {  # levels expected-line [case nofine]
  (cd "$W" && timeout 900 mpirun --bind-to none --oversubscribe -np 1 "$B/fds_amr" "${3:-cut}.fds" --two-level-run --levels "$1" --bind-only --blocking 2 --maxsize 32 --nofine "${4:-2}" \
     --regrid-every 0 --bed 0 100 -50 50 0 0.21 --outdir . --chid "l$1${3:-}" > "l$1${3:-}.log" 2>&1)
  if [ -z "$2" ]; then grep -q "L_ins" "$W/l$1${3:-}.log" && { echo "FAIL levels $1: unexpected bed line"; rc=1; } || echo "PASS levels $1: no bed line (no-op)"
  elif [ -n "${3:-}" ]; then grep -q "$2" "$W/l$1$3.log" && ! grep -q "BIND-ONLY PASS" "$W/l$1$3.log" && echo "PASS levels $1 $3: stops ($2)" || { echo "FAIL levels $1 $3: no stop"; rc=1; }
  else grep -q "$2" "$W/l$1.log" && grep -q "BIND-ONLY PASS" "$W/l$1.log" && echo "PASS levels $1: $2" || { echo "FAIL levels $1"; grep -aE "L_ins|LEVEL|Abort|error" "$W/l$1${3:-}.log" | head; rc=1; }; fi
}
run 2 "particle bed level L_ins 1, dx 0.5 0.5 0.5"
run 3 "particle bed level L_ins 2, dx 0.25 0.25 0.25"
run 1 ""
run 2 "does NOT cover the whole bed" cutv
sed -e 's/^&TIME T_END=1. \//\&TIME T_END=1., DT=0.05, LOCK_TIME_STEP=T \/\n\&DUMP DT_DEVC=0.05, SIG_FIGS=16, SIG_FIGS_EXP=3 \/\n\&REAC FUEL='"'PROPANE'"' \//' \
    -e 's/^&TAIL/\&RAMP ID='"'ig'"', T=0., F=0. \/\n\&RAMP ID='"'ig'"', T=1., F=1. \/\n\&SURF ID='"'IGNITION'"', HRRPUA=1500., TMP_FRONT=100., RAMP_Q='"'ig'"' \/\n\&VENT XB=0.0,1.6,-20.,20.,0.,0., SURF_ID='"'IGNITION'"', XYZ=0,0,0, SPREAD_RATE=1. \/\n\&DEVC ID='"'MLRB'"', QUANTITY='"'BURNING RATE'"', SPATIAL_STATISTIC='"'SURFACE INTEGRAL'"', SURF_ID='"'IGNITION'"', XB=0,1.6,-20,20,0,0 \/\n\&TAIL/' "$W/cut.fds" > "$W/cutb.fds"
for l in 1 2; do (cd "$W" && timeout 900 mpirun --bind-to none --oversubscribe -np 1 "$B/fds_amr" cutb.fds --two-level-run --levels $l --steps 20 --blocking 2 --maxsize 32 --regrid-every 0 \
  --bed 0 100 -50 50 0 0.21 --outdir . --chid b$l > b$l.log 2>&1 && mv cut_devc.csv b$l.csv) || { echo "FAIL burner levels $l: rc"; rc=1; }; done
grep -q "L_ins 1, dx 0.5 0.5 0.5$" "$W/b2.log" && [ -s "$W/b2.csv" ] && echo "PASS burner: bed level built at default --nofine, run rc 0" || { echo "FAIL burner"; rc=1; }
cat > "$W/lk.fds" <<'FDS'
&HEAD CHID='lk' /
&AMR MAX_LEVEL=0, BLOCKING_FACTOR=2, MAX_GRID_SIZE=32 /
&MESH IJK=8,8,4, XB=-2.0,6.0,-4.0,4.0,0.0,4.0 /
&TIME T_END=2., DT=0.05, LOCK_TIME_STEP=T /
&DUMP DT_DEVC=0.05, DT_HRR=0.05, SIG_FIGS=16, SIG_FIGS_EXP=3 /
&REAC FUEL='PROPANE' /
&RAMP ID='ig', T=0., F=0. /
&RAMP ID='ig', T=1., F=1. /
&SURF ID='IGNITION', HRRPUA=1500., TMP_FRONT=100., RAMP_Q='ig' /
&VENT XB=0.0,1.6,-2.,2.,0.,0., SURF_ID='IGNITION', XYZ=0,0,0, SPREAD_RATE=1. /
&VENT MB='XMIN', SURF_ID='OPEN' /
&VENT MB='XMAX', SURF_ID='OPEN' /
&VENT MB='YMIN', SURF_ID='OPEN' /
&VENT MB='YMAX', SURF_ID='OPEN' /
&VENT MB='ZMAX', SURF_ID='OPEN' /
&DEVC ID='MLRB', QUANTITY='BURNING RATE', SPATIAL_STATISTIC='SURFACE INTEGRAL', SURF_ID='IGNITION', XB=-2,6,-4,4,0,0 /
&DEVC ID='HRRB', QUANTITY='HRRPUA', SPATIAL_STATISTIC='SURFACE INTEGRAL', SURF_ID='IGNITION', XB=-2,6,-4,4,0,0 /
&TAIL /
FDS
sed 's/IJK=8,8,4/IJK=16,16,8/' "$W/lk.fds" > "$W/lkf.fds"
for r in "f1 lkf 1" "c1 lk 1" "c2 lk 2"; do set -- $r; (cd "$W" && timeout 900 mpirun --bind-to none --oversubscribe -np 1 "$B/fds_amr" $2.fds --two-level-run --levels $3 --steps 30 --blocking 2 --maxsize 32 --regrid-every 0 \
  --bed 0 2 -2 2 0 0.21 --outdir . --chid $1 > $1.log 2>&1 && cut -d, -f1-3 lk_devc.csv > $1.dv && cut -d, -f1,13 lk_hrr.csv > $1.hr) || { echo "FAIL like-dx run $1"; rc=1; }; done
cmpr() { paste -d, "$W/$1" "$W/$2" | awk -F, -v n=$3 'NR>2{for(c=2;c<=n;c++){x=$c; y=$(c+n); d=x-y; if(d<0)d=-d; if(x!=0)z++; r=x!=0?d/(x<0?-x:x):d; if(r>m)m=r}} END{printf "%.3e %d", m, z}'; }
read md zd <<< "$(cmpr f1.dv c2.dv 3)"; read mh zh <<< "$(cmpr f1.hr c2.hr 2)"; mc=$(paste -d, "$W/f1.dv" "$W/c1.dv" | awk -F, 'NR>2{d=$2-$5; s+=d<0?-d:d; t+=$2} END{printf "%.3e", s/t}')
awk -v a=$md -v b=$mh -v c=$mc -v z=$zd 'BEGIN{exit !(a<=1e-12 && b<=1e-12 && c>1e-3 && z>=20)}' \
  && echo "PASS like-dx: leaf burner MLR/HRRPUA max rel $md, HRR-file MLR max rel $mh ($zd nonzero); 1 m vs 0.5 m burner MLR (convergence, not a gate) sum|d|/sum $mc" \
  || { echo "FAIL like-dx: $md $mh $mc $zd"; rc=1; }
exit $rc
