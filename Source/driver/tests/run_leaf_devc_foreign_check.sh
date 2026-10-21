#!/bin/bash
# run_leaf_devc_foreign_check.sh <build-dir> [<run-dir>]: cases/leaf_devc_foreign.fds (CSIRO cut-down, 4 meshes, two levels) at 1, 2 and 4 MPI ranks with
# FDS_AMR_DEBUG_FOREIGN_BOXES set: every fine box that touches no domain face goes to the rank after its mesh's, so the gas devices in those boxes take the
# leaf devc path over another rank's mesh (dump.f90 FOREIGN pass, main.f90 LEAF_SUM_IN_RANK_ORDER).
# Pass: every run exits 0, boxes were moved at np 2 and 4, and every devc entry (c_cat_devc.csv, CATF) agrees with np 1 to 1e-11 relative (1e-13 floor).
set -u
BLD=$1; RUN=${2:-$BLD/run/leafdevc_foreign}; HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" || exit 1
fail=0
for np in 1 2 4; do
  d=$RUN/np$np; rm -rf $d; mkdir -p $d; cp "$HERE/cases/leaf_devc_foreign_vegetation_model.txt" $d/vegetation_model.txt
  python3 -c "import re; n=[0]
def f(m):
    k=n[0]*$np//4; n[0]+=1; return m.group(0).replace(' /',', MPI_PROCESS=%d /'%k)
open('$d/c.fds','w').write(re.sub(r'&MESH[^/]*/',f,open('$HERE/cases/leaf_devc_foreign.fds').read()))"
  ( cd $d && FDS_AMR_DEBUG_FOREIGN_BOXES=1 OMP_NUM_THREADS=1 mpirun --oversubscribe -np $np "$BLD/fds_amr" c.fds --two-level-run --levels 2 --blocking 2 --maxsize 16 \
      --regrid-every 0 --nofine 0 --bed 0 25 -15 15 0 0.21 --outdir . --chid c > run.log 2>&1 ); rc=$?
  mv=$(grep -a "boxes on a non-parent rank" $d/run.log | tail -1 | sed 's/.*: \([0-9]*\) of.*/\1/')
  echo "np=$np rc=$rc moved boxes: ${mv:-0}"
  [ $rc -eq 0 ] || fail=1
  [ $np -gt 1 ] && [ "${mv:-0}" -eq 0 ] && fail=1
done
python3 - "$RUN" <<'P' || fail=1
import sys,numpy as np
R=sys.argv[1]; h=[x.strip() for x in open(f'{R}/np1/c_cat_devc.csv').read().split('\n')[1].split(',')]
r=np.genfromtxt(f'{R}/np1/c_cat_devc.csv',delimiter=',',skip_header=2); bad=0
for n in (2,4):
    b=np.genfromtxt(f'{R}/np{n}/c_cat_devc.csv',delimiter=',',skip_header=2)
    if b.shape!=r.shape: print('np',n,'shape'); bad=1; continue
    d=np.abs(b-r); rel=d/np.maximum(np.abs(r),1e-13); rel[d<=1e-13]=0; m=rel.max(0); j=int(m.argmax())
    print(f'devc np1 vs np{n}: max rel {m.max():.2e} ({h[j]})'); bad|=int(m.max()>1e-11)
sys.exit(int(bad))
P
[ $fail -eq 0 ] && echo "PASS leaf_devc_foreign" || echo "FAIL leaf_devc_foreign (see $RUN)"
exit $fail
