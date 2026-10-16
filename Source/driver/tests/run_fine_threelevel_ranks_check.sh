#!/bin/bash
# run_fine_threelevel_ranks_check.sh <build-dir> [<run-dir>]: simple_test split into 4 meshes, three levels (--levels 3), at 1, 2 and 4 MPI ranks, for
#   vent:   fine boxes over the OPEN door region (fine walls whose level-0 parent mesh is on another rank: BND_STAGE_RATE exchange)
#   burner: fine boxes next to the burner (burner stays on L0, D-010)
# Pass: every run exits 0 with "gate PASS", and net boundary inflow and through-flow agree across the rank counts (relative 1e-6; the through-flow used to be rank-0 only).
# Known open: --tag-obst 1 --nofine -1 (burner on the fine level) at 3 levels fails the species balance (notes/fine-obst-walls.md R22).
set -u
BLD=$1; RUN=${2:-$BLD/run/threelevel_ranks}; HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" && cd "$RUN" || exit 1
python3 - "$HERE/../../../Verification/Adaptive_Mesh_Refinement/simple_test.fds" <<'P'
import sys,re
s=open(sys.argv[1]).read().replace('T_END=5.0','T_END=0.6')
m=re.search(r"&MESH[^/]*/",s).group(0)
four="\n".join("&MESH IJK=18,12,24, XB=%s,%s,%s,%s,0.0,2.4 /"%(x0,x1,y0,y1) for x0,x1 in((0.0,1.8),(1.8,3.6)) for y0,y1 in((0.0,1.2),(1.2,2.4)))
s=s.replace(m,four+"\n&RADI RADIATION=.FALSE. /")
open("simple_test.fds","w").write(s)
P
fail=0
run() { # name args
  local name=$1; shift; local ref=""
  for np in 1 2 4; do
    d=$RUN/${name}_$np; mkdir -p $d; python3 -c "import re,sys; n=[0]; np_=$np
def f(m):
    k=n[0]*np_//4; n[0]+=1; return m.group(0).replace(' /',', MPI_PROCESS=%d /'%k)
open('$d/simple_test.fds','w').write(re.sub(r'&MESH[^/]*/',f,open('simple_test.fds').read()))"; ( cd $d && OMP_NUM_THREADS=1 mpirun --oversubscribe -np $np "$BLD/fds_amr" simple_test.fds --two-level-run --levels 3 "$@" --steps 12 > run.log 2>&1 ); rc=$?
    g=$(grep -c "RESULT gate PASS" $d/run.log)
    v=$(grep -a "RESULT steps" $d/run.log | sed 's/.*boundary inflow \([-0-9.e+]*\) (through-flow \([-0-9.e+]*\)).*/\1 \2/')
    echo "$name np=$np rc=$rc gate=$g inflow/through: $v"
    [ $rc -eq 0 ] && [ $g -eq 1 ] || fail=1
    if [ -z "$ref" ]; then ref="$v"; else
      awk -v a="$ref" -v b="$v" 'BEGIN{split(a,x," ");split(b,y," ");for(i=1;i<=2;i++){d=x[i]-y[i];if(d<0)d=-d;m=(x[i]<0?-x[i]:x[i]);if(d>1e-6*m+1e-14)exit 1}}' || { echo "  MISMATCH vs np=1"; fail=1; }
    fi
  done
}
run vent --patch 4 15 4 19 --patch2 6 13 6 15
run burner --patch 0 15 0 7 --patch2 10 13 0 3
[ $fail -eq 0 ] && echo "PASS threelevel_ranks" || echo "FAIL threelevel_ranks (see $RUN)"
exit $fail
