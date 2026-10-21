#!/bin/bash
# run_fire_depth_leaf_check.sh <build-dir> [<run-dir>]: FIRE DEPTH on leaf cells (D-084, Architect ruling 2026-10-21). cases/leaf_devc_foreign.fds (CSIRO cut-down,
# 4 meshes at 1 m, level 1 at 0.5 m over the bed x 0..25, y -15..15, z 0..0.5) to T=1.5 s. FIRE DEPTH lines inside the level-1 region at z 0.1: DA at y=-0.3 (inside
# the 0.5 m row y -0.5..0, which is half of a 1 m level-0 row) and DB at y=-1.0 (on a level-0 and level-1 face: the scan points go to the row above it, y -1..-0.5, as VB).
# Reference: the leaf VOLUME statistic of HRRPUV over the same 0.5 m cell row with QUANTITY_RANGE 40 kW/m3 (= MIN(200,20/0.5), leaf threshold), divided by the
# row cross-section 0.25 m2: the leaf length above threshold. A covered 1 m level-0 cell counted by FIRE DEPTH (1 m steps, 20 kW/m3) breaks the match.
# Run at np 1, 2, 4 with FDS_AMR_DEBUG_FOREIGN_BOXES (fine boxes on a non-parent rank: FOREIGN pass). Pass: DA == VA/0.25 and DB == VB/0.25 to 1e-9 in every row,
# a nonzero depth in some row, and np 2/4 equal to np 1.
set -u
BLD=$1; RUN=${2:-$BLD/run/firedepth_leaf}; HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" || exit 1
fail=0
for np in 1 2 4; do
  d=$RUN/np$np; rm -rf $d; mkdir -p $d; cp "$HERE/cases/leaf_devc_foreign_vegetation_model.txt" $d/vegetation_model.txt
  python3 -c "import re; n=[0]
def f(m):
    k=n[0]*$np//4; n[0]+=1; return m.group(0).replace(' /',', MPI_PROCESS=%d /'%k)
s=re.sub(r'&MESH[^/]*/',f,open('$HERE/cases/leaf_devc_foreign.fds').read())
s=s.replace('T_END=0.5','T_END=1.5').replace('&TAIL',\"&DEVC ID='DA', QUANTITY='FIRE DEPTH', XB=0,25,-0.3,-0.3,0.1,0.1 /\n&DEVC ID='VA', QUANTITY='HRRPUV', SPATIAL_STATISTIC='VOLUME', QUANTITY_RANGE=40.,1.E9, XB=0,25,-0.5,0,0,0.5 /\n&DEVC ID='DB', QUANTITY='FIRE DEPTH', XB=0,25,-1.0,-1.0,0.1,0.1 /\n&DEVC ID='VB', QUANTITY='HRRPUV', SPATIAL_STATISTIC='VOLUME', QUANTITY_RANGE=40.,1.E9, XB=0,25,-1.0,-0.5,0,0.5 /\n&TAIL\")
open('$d/c.fds','w').write(s)"
  ( cd $d && FDS_AMR_DEBUG_FOREIGN_BOXES=1 OMP_NUM_THREADS=1 mpirun --oversubscribe -np $np "$BLD/fds_amr" c.fds --two-level-run --levels 2 --blocking 2 --maxsize 16 \
      --regrid-every 0 --nofine 0 --bed 0 25 -15 15 0 0.21 --outdir . --chid c > run.log 2>&1 ); rc=$?
  echo "np=$np rc=$rc"; [ $rc -eq 0 ] || fail=1
  python3 - $d/c_cat_devc.csv $np <<'P' || fail=1
import sys,numpy as np
f=sys.argv[1]; h=[x.strip().strip('"') for x in open(f).read().split('\n')[1].split(',')]
a=np.genfromtxt(f,delimiter=',',skip_header=2); a=a[a[:,0]>0.05]
c=lambda k:a[:,h.index(k)]
ea=np.abs(c('DA')-c('VA')/0.25).max(); eb=np.abs(c('DB')-c('VB')/0.25).max(); mx=max(c('DA').max(),c('DB').max())
print(f"np={sys.argv[2]} rows {len(a)} DA {np.round(c('DA'),3).tolist()} DB {np.round(c('DB'),3).tolist()} max|D-V/0.25| {ea:.1e} {eb:.1e}")
sys.exit(int(not (len(a)>0 and ea<1e-9 and eb<1e-9 and mx>0)))
P
done
python3 - "$RUN" <<'P' || fail=1
import sys,numpy as np
R=sys.argv[1]; r=np.genfromtxt(f'{R}/np1/c_cat_devc.csv',delimiter=',',skip_header=2); bad=0
for n in (2,4):
    b=np.genfromtxt(f'{R}/np{n}/c_cat_devc.csv',delimiter=',',skip_header=2); d=np.abs(b-r); rel=(d/np.maximum(np.abs(r),1e-13)); rel[d<=1e-13]=0
    print(f'devc np1 vs np{n}: max rel {rel.max():.2e}'); bad|=int(rel.max()>1e-11)
sys.exit(bad)
P
[ $fail -eq 0 ] && echo "PASS fire_depth_leaf" || echo "FAIL fire_depth_leaf (see $RUN)"
exit $fail
