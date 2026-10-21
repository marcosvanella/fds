#!/bin/bash
# W2a (PV-10) particle bed level check, setup only (--bind-only): a gas-only cut-down C064 (x -10..70, y -20..20, z 0..20 at 1 m) with the
# C064 bed XB passed by --bed (the &INIT PART_ID path is the same list once particle cases pass the driver scope, W0/W7).
# (a) 2 levels: L_ins 1 (0.5 m) covers the bed; (b) 3 levels: L_ins 2 (0.25 m); (c) 1 level: no bed line (no-op).
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
run() {  # levels expected-line
  (cd "$W" && timeout 900 mpirun --bind-to none --oversubscribe -np 1 "$B/fds_amr" cut.fds --two-level-run --levels "$1" --bind-only --blocking 2 --maxsize 32 --nofine 0 \
     --regrid-every 0 --bed 0 100 -50 50 0 0.21 --outdir . --chid "l$1" > "l$1.log" 2>&1)
  if [ -z "$2" ]; then grep -q "L_ins" "$W/l$1.log" && { echo "FAIL levels $1: unexpected bed line"; rc=1; } || echo "PASS levels $1: no bed line (no-op)"
  else grep -q "$2" "$W/l$1.log" && grep -q "BIND-ONLY PASS" "$W/l$1.log" && echo "PASS levels $1: $2" || { echo "FAIL levels $1"; grep -aE "L_ins|LEVEL|Abort|error" "$W/l$1.log" | head; rc=1; }; fi
}
run 2 "particle bed level L_ins 1, dx 0.5 0.5 0.5"
run 3 "particle bed level L_ins 2, dx 0.25 0.25 0.25"
run 1 ""
exit $rc
