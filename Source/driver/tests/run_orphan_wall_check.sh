#!/bin/bash
# run_orphan_wall_check.sh <fds_amr binary> [<run-dir>]: fine SOLID wall faces without a level-0 wall (no solid takes their convective flux; ruling 2, counted, not fixed).
# CSIRO C064 level-set crude (cases/levelset_tag_phi_ground.fds), 2 levels --ls-fine, 20 steps: the count and the energy (stdout RESULT line and the .out line) must both be 0.
set -u
BIN=$1; RUN=${2:-$PWD/orphan_wall_run}; SRC=$(cd "$(dirname "$0")" && pwd)/cases/levelset_tag_phi_ground.fds
mkdir -p "$RUN" && cd "$RUN" || exit 1
python3 -c "import sys; s=open(sys.argv[1]).read(); open('c.fds','w').write(s if '&AMR' in s else s.replace('&TAIL', '&AMR BLOCKING_FACTOR=2 /\n&TAIL'))" "$SRC"   # no &AMR line in the case: the default blocking factor 8 does not divide 12x12x4
OMP_NUM_THREADS=1 FDSTL_NO_MESH_DUMPS=1 timeout 900 "$BIN" c.fds --two-level-run --levels 2 --init-patch --patch 5 9 0 11 --blocking 2 --maxsize 8 --regrid-every 2 --tag-phi --ls-fine \
    --steps 20 > run.log 2>&1 || { echo "FAIL orphan_wall: run rc (see $RUN/run.log)"; exit 1; }
rm -f ./*.s3d ./*.s3d.sz ./*.sf ./*.bf ./*.prt5
L=$(grep -a "RESULT fine walls without a level-0 wall" run.log); O=$(grep -a "AMR fine walls without a level-0 wall" ./*.out 2>/dev/null)
echo "$L"; echo "out: $O"; grep -a "RESULT WARNING fine walls remapped" run.log
N=$(echo "$L" | sed 's/.*wall: \([^ ]*\) .*/\1/'); E=$(echo "$L" | sed 's/.*energy \([^ ]*\) J.*/\1/')
[ -n "$O" ] || { echo "FAIL orphan_wall: no .out line"; exit 1; }
python3 -c "import sys; sys.exit(0 if float('$N') == 0 and float('$E') == 0 else 1)" && echo "PASS orphan_wall" || { echo "FAIL orphan_wall: count $N energy $E J"; exit 1; }
