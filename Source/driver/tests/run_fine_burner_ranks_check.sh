#!/bin/bash
# run_fine_burner_ranks_check.sh <build-dir> [<run-dir>]: burner on the fine level (--nofine -1), simple_test split into 4 meshes on 1 and 4 MPI ranks, 2 levels; at 4 ranks the
#   rank boundaries cut the fine level. The fine layout must not depend on the rank count (AMReX refine_grid_layout chopped it per rank).
# Pass: every run exits 0 with "gate PASS", and the leaf fuel MLR (c_multi_level.csv mlr_leaf) and the HRR-file Q_CONV, Q_COND, MLR_PROPANE agree across the rank counts (relative 1e-11).
set -u
BLD=$(realpath "$1"); RUN=${2:-$BLD/run/fine_burner_ranks}
IN=$(cd "$(dirname "$0")" && pwd)/../../../Verification/Adaptive_Mesh_Refinement/simple_test.fds
mkdir -p "$RUN" && cd "$RUN" || exit 1
python3 - "$IN" <<'P'
import sys,re
s=open(sys.argv[1]).read().replace('T_END=5.0','T_END=0.6')
m=re.search(r"&MESH[^/]*/",s).group(0)
four="\n".join("&MESH IJK=18,12,24, XB=%s,%s,%s,%s,0.0,2.4 /"%(x0,x1,y0,y1) for x0,x1 in((0.0,1.8),(1.8,3.6)) for y0,y1 in((0.0,1.2),(1.2,2.4)))
s=s.replace(m,four+"\n&RADI RADIATION=.FALSE. /\n&DUMP SIG_FIGS=15, SIG_FIGS_EXP=3 /")
open("simple_test.fds","w").write(s)
P
fail=0
run() { # name args
  local name=$1; shift; local ref=""
  for np in 1 4; do
    d=$RUN/${name}_$np; mkdir -p $d; python3 -c "import re; n=[0]; np_=$np
def f(m):
    k=n[0]*np_//4; n[0]+=1; return m.group(0).replace(' /',', MPI_PROCESS=%d /'%k)
open('$d/simple_test.fds','w').write(re.sub(r'&MESH[^/]*/',f,open('simple_test.fds').read()))"; ( cd $d && OMP_NUM_THREADS=1 mpirun --oversubscribe -np $np "$BLD/fds_amr" simple_test.fds --two-level-run "$@" --nofine -1 --steps 12 > run.log 2>&1 ); rc=$?
    g=$(grep -c "RESULT gate PASS" $d/run.log)
    v=$(python3 -c "
import csv,glob
r=list(csv.DictReader(open(glob.glob('$d/*multi_level.csv')[0])));h=list(csv.reader(open('$d/simple_test_hrr.csv')))
print(r[-1]['mlr_leaf'].strip(),*[h[-1][h[1].index(k)].strip() for k in ('Q_CONV','Q_COND','MLR_PROPANE')])" 2>/dev/null)
    echo "$name np=$np rc=$rc gate=$g mlr_leaf Q_CONV Q_COND MLR_PROPANE: $v"
    [ $rc -eq 0 ] && [ $g -eq 1 ] && [ -n "$v" ] || fail=1
    if [ -z "$ref" ]; then ref="$v"; else
      awk -v a="$ref" -v b="$v" 'BEGIN{n=split(a,x," ");if(split(b,y," ")!=n||n<2)exit 1;for(i=1;i<=n;i++){d=x[i]-y[i];if(d<0)d=-d;m=(x[i]<0?-x[i]:x[i]);if(d>1e-11*m+1e-14)exit 1}}' || { echo "  MISMATCH vs np=1"; fail=1; }
    fi
  done
}
run l2 --levels 2 --patch 0 7 0 7
[ $fail -eq 0 ] && echo "PASS fine_burner_ranks" || echo "FAIL fine_burner_ranks (see $RUN)"
exit $fail
