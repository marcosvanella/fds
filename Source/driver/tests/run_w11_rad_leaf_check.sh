#!/bin/bash
# run_w11_rad_leaf_check.sh <fds_amr> <c10 cut-down C064 input (4 meshes at 1 m, &AMR line, vegetation_model.txt beside it)> [work-dir] [ranks] [steps]
# P0 W11: radiation particle terms on the leaf. c10 at 2 levels tagged on temperature, regrid every 5 steps, radiation on, FDS_RAD_LEAF=1 (radi.f90 RADW11 line per
# level-0 mesh and radiation update): over the cells covered by level 1 the direct particle terms KAPPA_PART and KFST4_PART sum to exactly 0 and the averaged-down
# emission (RAD_SRC component 2) is present. PASS: covered cells in every line, both sums 0 in every line, emission > 0 in every line, TWO-LEVEL gate PASS.
# 1 level: identical to stock is V1 (the RADW11 line is written only with FDS_RAD_LEAF=1).
set -u
BIN=$1; IN=$2; W=${3:-$PWD/w11_rad_leaf}; NP=${4:-1}; N=${5:-12}
rm -rf "$W"; mkdir -p "$W"; cp "$IN" "$W/c10.fds"; cp "$(dirname "$IN")/vegetation_model.txt" "$W/"
cd "$W" && FDS_RAD_LEAF=1 OMP_NUM_THREADS=1 timeout 1800 mpirun --bind-to none --oversubscribe -np "$NP" "$BIN" c10.fds --two-level-run --levels 2 --tag-temp 320 \
    --regrid-every 5 --steps "$N" > log 2>&1; rc=$?
rm -f *.s3d* *.sf *.bf *.prt5* *.restart
grep -a "^RADW11" log | awk -v rc=$rc -v np=$NP -v g="$(grep -ac 'RESULT gate PASS' log)" '{n++; if($4==0)nc++; if($5!=0||$6!=0)nz++; if(!($7>0))ne++}
  END{ok=(rc==0 && g>0 && n>0 && nc==0 && nz==0 && ne==0); printf "%s w11_rad_leaf np%s: %d RADW11 lines, without covered cells %d, nonzero direct particle terms in covered cells %d, without averaged emission %d, rc %d, gate PASS %d\n", ok?"PASS":"FAIL", np, n, nc, nz, ne, rc, g; exit !ok}'
