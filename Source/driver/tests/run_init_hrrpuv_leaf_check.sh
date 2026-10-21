#!/bin/bash
# run_init_hrrpuv_leaf_check.sh <build-dir> [<run-dir>]: &INIT HRRPUV on fine boxes (radi.f90 ADD_VOLUMETRIC_HEAT_SOURCE; IN%VOLUME_ADJUST is sized NMESHES, a fine
# box NM>NMESHES read past it). Gas-only box 16x16x8 m at 1 m with level 1 (0.5 m) over --bed 2 10 2 10 0 0.21, an HRRPUV block x 3.2..6.7, y 3.3..5.6, z 0..1.3
# (not aligned to either grid) inside level 1, RADIATION off (mode 0) and on (mode 1). Pass: every HRR-file row from t>0 of the 2-level run equals the 1-level run
# and the analytic HRRPUV*block volume (3.5*2.3*1.3*500 kW) to 1e-9 relative. Build <build-dir> with -fcheck=bounds to catch the out-of-bounds read itself.
set -u
BLD=$1; RUN=${2:-$BLD/run/init_hrrpuv}; mkdir -p "$RUN" || exit 1; cd "$RUN" || exit 1
fail=0
for rad in F T; do
cat > h$rad.fds <<FDS
&HEAD CHID='h$rad' /
&AMR MAX_LEVEL=0, BLOCKING_FACTOR=2, MAX_GRID_SIZE=32 /
&MESH IJK=16,16,8, XB=0.,16.,0.,16.,0.,8. /
&TIME T_END=0.5, DT=0.05, LOCK_TIME_STEP=T /
&DUMP DT_HRR=0.05, DT_DEVC=0.05, SIG_FIGS=16, SIG_FIGS_EXP=3 /
&RADI RADIATION=$rad /
&INIT XB=3.2,6.7,3.3,5.6,0.0,1.3, HRRPUV=500. /
&VENT MB='XMIN', SURF_ID='OPEN' /
&VENT MB='XMAX', SURF_ID='OPEN' /
&VENT MB='YMIN', SURF_ID='OPEN' /
&VENT MB='YMAX', SURF_ID='OPEN' /
&VENT MB='ZMAX', SURF_ID='OPEN' /
&TAIL /
FDS
  for l in 1 2; do
    OMP_NUM_THREADS=1 timeout 900 "$BLD/fds_amr" h$rad.fds --two-level-run --levels $l --blocking 2 --maxsize 32 --regrid-every 0 --nofine 0 --bed 2 10 2 10 0 0.21 \
      --outdir . --chid h${rad}$l > h${rad}$l.log 2>&1 || { echo "FAIL rad=$rad levels $l: rc (see $RUN/h${rad}$l.log)"; grep -a -m3 -i "bound\|error" h${rad}$l.log; fail=1; continue; }
    mv h${rad}_hrr.csv h${rad}$l.hrr
  done
  python3 - h${rad}1.hrr h${rad}2.hrr $rad <<'P' || fail=1
import sys,numpy as np
a=np.genfromtxt(sys.argv[1],delimiter=',',skip_header=2); b=np.genfromtxt(sys.argv[2],delimiter=',',skip_header=2)
a=a[a[:,0]>0]; b=b[b[:,0]>0]; ref=3.5*2.3*1.3*500.
e12=np.abs(b[:,1]-a[:,1]).max()/ref; e2=np.abs(b[:,1]-ref).max()/ref
print(f"rad={sys.argv[3]} rows {len(a)}/{len(b)} HRR 1L {a[-1,1]:.6f} 2L {b[-1,1]:.6f} analytic {ref:.6f}; max rel 2L-1L {e12:.1e}, 2L-analytic {e2:.1e}")
sys.exit(int(not (len(a)==len(b)>0 and e12<1e-9 and e2<1e-9)))
P
done
[ $fail -eq 0 ] && echo "PASS init_hrrpuv_leaf" || echo "FAIL init_hrrpuv_leaf (see $RUN)"
exit $fail
