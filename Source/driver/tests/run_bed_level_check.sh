#!/bin/bash
# W2a (PV-10) particle bed level check, setup only (--bind-only): a gas-only cut-down C064 (x -10..70, y -20..20, z 0..20 at 1 m) with the
# C064 bed XB passed by --bed (the &INIT PART_ID path is the same list once particle cases pass the driver scope, W0/W7).
# (a) 2 levels: L_ins 1 (0.5 m) covers the bed; (b) 3 levels: L_ins 2 (0.25 m); (c) 1 level: no bed line (no-op);
# (d) a heat-flux VENT inside the bed kept on level 0 (D-010, --nofine 1): the bed level cannot cover the bed, the run stops.
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
  (cd "$W" && timeout 900 mpirun --bind-to none --oversubscribe -np 1 "$B/fds_amr" "${3:-cut}.fds" --two-level-run --levels "$1" --bind-only --blocking 2 --maxsize 32 --nofine "${4:-0}" \
     --regrid-every 0 --bed 0 100 -50 50 0 0.21 --outdir . --chid "l$1${3:-}" > "l$1${3:-}.log" 2>&1)
  if [ -z "$2" ]; then grep -q "L_ins" "$W/l$1${3:-}.log" && { echo "FAIL levels $1: unexpected bed line"; rc=1; } || echo "PASS levels $1: no bed line (no-op)"
  elif [ -n "${3:-}" ]; then grep -q "$2" "$W/l$1$3.log" && ! grep -q "BIND-ONLY PASS" "$W/l$1$3.log" && echo "PASS levels $1 $3: stops ($2)" || { echo "FAIL levels $1 $3: no stop"; rc=1; }
  else grep -q "$2" "$W/l$1.log" && grep -q "BIND-ONLY PASS" "$W/l$1.log" && echo "PASS levels $1: $2" || { echo "FAIL levels $1"; grep -aE "L_ins|LEVEL|Abort|error" "$W/l$1${3:-}.log" | head; rc=1; }; fi
}
run 2 "particle bed level L_ins 1, dx 0.5 0.5 0.5"
run 3 "particle bed level L_ins 2, dx 0.25 0.25 0.25"
run 1 ""
run 2 "does NOT cover the whole bed" cutv 1
exit $rc
