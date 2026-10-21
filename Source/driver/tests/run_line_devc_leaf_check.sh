#!/bin/bash
# run_line_devc_leaf_check.sh <build-dir> [<run-dir>]: TRANSMISSION / PATH OBSCURATION and LAYER HEIGHT / UPPER / LOWER TEMPERATURE on leaf cells (dump.f90
# LEAF_LINE_DEVICE, D-084; FIRE DEPTH has run_fire_depth_leaf_check.sh). A covered level-0 cell counted by either device breaks the check.
# (a) cases/leaf_devc_foreign.fds (CSIRO cut-down, 1 m, level 1 at 0.5 m over the bed z 0..0.5) to 1.5 s at np 1, 2, 4 with FDS_AMR_DEBUG_FOREIGN_BOXES:
#     PO (PATH OBSCURATION) and TR (TRANSMISSION) along x 0..25 at y=-0.3, z=0.1 against the leaf VOLUME INTEGRAL of soot density RS over the 0.5 m row
#     y -0.5..0, z 0..0.5: -ln(1-PO/100) = -ln(TR/100)*25 = 8700*RS/0.25 (soot mass extinction 8700 m2/kg) to 1e-8, soot present, np 2/4 equal to np 1.
# (b) 16x16x8 m box at 1 m, level 1 (0.5 m) over the full height of x,y 2..10 (--bed 2 10 2 10 0 8), a 2 MW HRRPUV block x 6.5..7.5, y 6..7, z 0..1 under a closed
#     ceiling (XMIN open) whose plume and ceiling jet form a hot layer over the column cell x,y 6..6.5, 2 s (0.01 s steps) at np 1 and 2 (foreign boxes): LAYER HEIGHT and LOWER TEMPERATURE at (6.2,6.2) equal the two-layer
#     reduction of the 16 leaf point temperatures of that 0.5 m column (to 1e-6 m and 1e-9 K); reading the covered 1 m column (averaged down) breaks it.
set -u
BLD=$1; RUN=${2:-$BLD/run/line_devc_leaf}; HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" || exit 1
fail=0
for np in 1 2 4; do
  d=$RUN/np$np; rm -rf $d; mkdir -p $d; cp "$HERE/cases/leaf_devc_foreign_vegetation_model.txt" $d/vegetation_model.txt
  python3 -c "import re; n=[0]
def f(m):
    k=n[0]*$np//4; n[0]+=1; return m.group(0).replace(' /',', MPI_PROCESS=%d /'%k)
s=re.sub(r'&MESH[^/]*/',f,open('$HERE/cases/leaf_devc_foreign.fds').read())
s=s.replace('T_END=0.5','T_END=1.5').replace('&TAIL',\"&DEVC ID='PO', QUANTITY='PATH OBSCURATION', XB=0,25,-0.3,-0.3,0.1,0.1, TIME_AVERAGED=F /\n&DEVC ID='TR', QUANTITY='TRANSMISSION', XB=0,25,-0.3,-0.3,0.1,0.1, TIME_AVERAGED=F /\n&DEVC ID='RS', QUANTITY='DENSITY', SPEC_ID='SOOT', SPATIAL_STATISTIC='VOLUME INTEGRAL', XB=0,25,-0.5,0,0,0.5, TIME_AVERAGED=F /\n&TAIL\")
open('$d/c.fds','w').write(s)"
  ( cd $d && FDS_AMR_DEBUG_FOREIGN_BOXES=1 OMP_NUM_THREADS=1 mpirun --oversubscribe -np $np "$BLD/fds_amr" c.fds --two-level-run --levels 2 --blocking 2 --maxsize 16 \
      --regrid-every 0 --nofine 0 --bed 0 25 -15 15 0 0.21 --outdir . --chid c > run.log 2>&1 ); rc=$?
  echo "np=$np rc=$rc"; [ $rc -eq 0 ] || fail=1
  python3 - $d/c_cat_devc.csv $np <<'P' || fail=1
import sys,numpy as np
f=sys.argv[1]; h=[x.strip().strip('"') for x in open(f).read().split('\n')[1].split(',')]
a=np.genfromtxt(f,delimiter=',',skip_header=2); a=a[a[:,0]>0.05]; c=lambda k:a[:,h.index(k)]
ref=8700.*c('RS')/0.25; po=-np.log(1-c('PO')/100.); tr=-np.log(c('TR')/100.)*25.
s=np.maximum(np.abs(ref),1e-30); e1=(np.abs(po-ref)/s).max(); e2=(np.abs(tr-ref)/s).max()
print(f"np={sys.argv[2]} rows {len(a)} extinction integral max {ref.max():.3e}; max rel |PO-ref| {e1:.1e}, |TR-ref| {e2:.1e}")
sys.exit(int(not (len(a)>0 and ref.max()>1e-6 and e1<1e-8 and e2<1e-8)))
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
d=$RUN/layer; rm -rf $d; mkdir -p $d
python3 - $d/l.fds <<'P'
s="""&HEAD CHID='l' /
&AMR MAX_LEVEL=0, BLOCKING_FACTOR=2, MAX_GRID_SIZE=32 /
&MESH IJK=8,16,8, XB=0.,8.,0.,16.,0.,8., MPI_PROCESS=0 /
&MESH IJK=8,16,8, XB=8.,16.,0.,16.,0.,8., MPI_PROCESS=1 /
&TIME T_END=2.0, DT=0.01, LOCK_TIME_STEP=T /
&DUMP DT_DEVC=0.01, SIG_FIGS=16, SIG_FIGS_EXP=3 /
&RADI RADIATION=F /
&INIT XB=6.5,7.5,6.0,7.0,0.0,1.0, HRRPUV=2000. /
&VENT MB='XMIN', SURF_ID='OPEN' /
&DEVC ID='LH', QUANTITY='LAYER HEIGHT', XB=6.2,6.2,6.2,6.2,0.,8., TIME_AVERAGED=F /
&DEVC ID='LT', QUANTITY='LOWER TEMPERATURE', XB=6.2,6.2,6.2,6.2,0.,8., TIME_AVERAGED=F /
"""
for k in range(16): s+="&DEVC ID='P%d', QUANTITY='TEMPERATURE', XYZ=6.2,6.2,%.2f, TIME_AVERAGED=F /\n"%(k,0.25+0.5*k)
open(__import__('sys').argv[1],'w').write(s+"&TAIL /\n")
P
for np in 1 2; do
  sed -i "s/XB=8.,16.,0.,16.,0.,8., MPI_PROCESS=[01]/XB=8.,16.,0.,16.,0.,8., MPI_PROCESS=$((np-1))/" $d/l.fds
  ( cd $d && FDS_AMR_DEBUG_FOREIGN_BOXES=1 OMP_NUM_THREADS=1 mpirun --oversubscribe -np $np "$BLD/fds_amr" l.fds --two-level-run --levels 2 --blocking 2 --maxsize 8 \
      --regrid-every 0 --nofine 0 --bed 2 10 2 10 0 8 --outdir . --chid l > run$np.log 2>&1 && mv l_devc.csv l$np.csv ); rc=$?
  echo "layer np=$np rc=$rc"; [ $rc -eq 0 ] || { fail=1; continue; }
  python3 - $d/l$np.csv $np <<'P' || fail=1
# reference: the two-layer reduction (GET_LAYER_HEIGHT_INTEGRALS, UPDATE_DEVICES_2) of the 16 leaf point temperatures of the column cell (0.5 m each)
import sys,numpy as np
f=sys.argv[1]; h=[x.strip().strip('"') for x in open(f).read().split('\n')[1].split(',')]
a=np.genfromtxt(f,delimiter=',',skip_header=2); a=a[a[:,0]>0.01]; c=lambda k:a[:,h.index(k)]
T=np.stack([c('P%d'%k) for k in range(16)],1)+273.15; dz=0.5; H=8.
zf=np.arange(16)*dz; LH=[]; zp=0.
for r in range(len(T)):
    t=T[r]; I1=(t*dz).sum(); I2=(dz/t).sum(); TL=t[0]
    I3=sum(t[k]*dz if zf[k]>=zp else (t[k]*(zf[k]+dz-zp) if zf[k]+dz>zp else 0.) for k in range(16))
    TU=TL if abs(I3)<=1e-12 else I3/(H-zp)
    zi=TL*(I1*I2-H*H)/(I1+I2*TL**2-2*TL*H); LH.append(H if TU-TL<1. else min(zi,H)); zp=zi
LH=np.array(LH)[1:]; TLs=T[1:,0]; a=a[1:]  # row 1 on: the previous Z_INT (UPPER T, I_3) is known
e1=np.abs(c('LH')-LH).max(); e2=np.abs(c('LT')+273.15-TLs).max(); spread=(T.max(1)-T.min(1)).max()
print(f"layer np={sys.argv[2]} rows {len(a)} LH {np.round(c('LH'),3).tolist()} ref {np.round(LH,3).tolist()}; max|LH-ref| {e1:.1e} m, max|LT-ref| {e2:.1e} K, T spread {spread:.0f} K")
sys.exit(int(not (len(a)>0 and e1<1e-6 and e2<1e-9 and spread>20 and (c('LH')<7.99).any())))
P
done
[ $fail -eq 0 ] && echo "PASS line_devc_leaf" || echo "FAIL line_devc_leaf (see $RUN)"
exit $fail
