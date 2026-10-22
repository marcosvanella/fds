#!/bin/bash
# run_open_step_check.sh <fds_amr binary> [<run-dir>] [ground=tmpf|thick]: cases/open_step.fds (rad_init_csiro without the bed: &WIND L profile, TMPA 32 C, thick DIRT ground, OPEN sides and top;
# radiation and noise off, DT 0.02 locked) as 2 levels with one level-1 box over the whole domain and as 1 level at the level-1 cell size, 10 steps.
# Stock: the OPEN ghost cells take TMP_F/ZZ_F/RHO_F (wall.f90 INFLOW_OUTFLOW) and the velocity ghosts VELOCITY_BC of the open wall at every stage.
# Pass: the level-1 fields (FDSTL_LEVEL_CSV: RHO, TMP, U, V, W, H, D, ZZ) cell by cell and every gas DEVC (mass, mean T, velocity statistics, points near the open sides) of
# steps 1-10 within 1e-10 of the 1L ones, relative to max(max|1L| over the run, 1e-12 x the quantity's natural scale: 36.6 kg, 305 K, 1 m/s; fields: their max|1L| over the run; D: max|U| RDX of the step, see below).
# ground=tmpf (default, gated): the GROUND is a fixed-temperature surface (TMP_FRONT=TMPA, no solid), so the gas sees the same wall in 1L and 2L.
# ground=thick (report only, always exit 0): the 0.2 m DIRT ground; its solid is the level-0 one (L0Q, by ruling), step 1 D of the ground cells differs by -1.27e-5 1/s.
set -u
BIN=$1; RUN=${2:-$PWD/open_step_run}; GROUND=${3:-tmpf}; HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" && cd "$RUN" || exit 1
for L in 1 2; do
    rm -rf r$L; mkdir r$L; cp "$HERE/cases/leaf_devc_foreign_vegetation_model.txt" r$L/vegetation_model.txt
    G=s/x/x/; [ "$GROUND" = tmpf ] && G="s/^&SURF ID='GROUND'.*/\&SURF ID='GROUND', DEFAULT=T, TMP_FRONT=32. \//"
    if [ $L = 1 ]; then sed -e 's/IJK=8,8,4/IJK=16,16,8/' -e "$G" "$HERE/cases/open_step.fds" > r$L/c.fds; else sed -e "$G" "$HERE/cases/open_step.fds" > r$L/c.fds; fi
    (cd r$L && FDSTL_LEVEL_CSV=. OMP_NUM_THREADS=1 timeout 600 "$BIN" c.fds --two-level-run --levels $L --patch 0 7 0 3 --maxsize ${MAXSIZE:-16} --nofine -1 --steps 10 --outdir . > log 2>&1) \
        || { echo "FAIL open_step: $L-level run rc (see $RUN/r$L/log)"; exit 1; }
    rm -f r$L/*.s3d* r$L/*.sf r$L/*.bf r$L/*.prt5
done
python3 - <<'P'; rc=$?
import csv,sys
a=list(csv.reader(open('r1/open_step_cat_devc.csv')));b=list(csv.reader(open('r2/open_step_cat_devc.csv')))
h=[x.strip() for x in a[1]]
A=[[float(v) for v in r] for r in a[3:]];B=[[float(v) for v in r] for r in b[3:]]   # steps 1..10 (the t=0 row of the 2-level run is level 0)
if len(A)<10 or len(B)<10: print("  fewer than 10 steps"); sys.exit(1)
NAT={'M':36.6,'TM':305.,'T1':305.,'UM':1.,'UMAX':1.,'WMAX':1.,'VMIN':1.,'U1':1.,'U2':1.}   # natural scale of each quantity (kg, K, m/s)
sc=[max(max(abs(r[c]) for r in A),1e-12*NAT.get(h[c],1.)) for c in range(len(h))]   # absolute floor: 1e-12 of the natural scale
w=(0.,'',0.)
for x,y in zip(A,B):
    for c in range(1,len(h)):
        e=abs(x[c]-y[c])/sc[c]
        if e>w[0]: w=(e,h[c],x[0])
print("  DEVC worst |2L-1L|/max|1L| %.3e (%s at t=%.2f)"%w)
import glob
def rd(L,s):
    d={}
    for f in glob.glob('r%d/level_s%d_r*.csv'%(L,s)):
        for n,c,i,j,k,v in csv.reader(open(f)): d[(n,int(c),int(i),int(j),int(k))]=float(v)
    return d
F1=[rd(1,s) for s in range(1,11)];F2=[rd(2,s) for s in range(1,11)]
fs={}
for F in F1:
    for key,v in F.items(): fs[key[0]]=max(fs.get(key[0],0.),abs(v))
# D is gated on a term-scaled absolute norm, |dD| <= 1e-10 max|U| RDX (max|U| of the 1L step over the domain, RDX of the level-1 cell): D is a near-cancelling
# sum of terms of about U RDX (|D| ~ 5e-5 1/s against terms ~ 18 1/s here), so relative to max|D| its round-off is large; nudging TMPA by 1e-10 K in the 1L run alone
# moved D by 3.1e-11 of max|D| (TMP 3.3e-13), and the 2L-1L D difference is a few ulps of the terms (4e-14 1/s at step 2). Every other field: 1e-10 of its max|1L|.
import re
m=re.search(r"&MESH IJK=(\d+),\d+,\d+, XB=([-\d.]+),([-\d.]+)",open('r1/c.fds').read()); RDX=int(m.group(1))/(float(m.group(3))-float(m.group(2)))
wf=(0.,'',0,None);wd=(0.,0,None)
for s,(a,b) in enumerate(zip(F1,F2),1):
    if not a or set(a)!=set(b): print("  level-1 fields of step %d missing or of another layout"%s); sys.exit(1)
    dsc=1e-10*max(abs(v) for key,v in a.items() if key[0]=='U')*RDX
    for key,v in a.items():
        if key[0]=='D':
            e=abs(v-b[key])/dsc*1e-10
            if e>wd[0]: wd=(e,s,key[1:])
            continue
        e=abs(v-b[key])/max(fs[key[0]],1e-300)
        if e>wf[0]: wf=(e,key[0],s,key[1:])
print("  fields worst |2L-1L|/max|1L| %.3e (%s step %d, comp/i/j/k %s)"%wf)
print("  D worst |2L-1L|/(max|U| RDX) %.3e (step %d, comp/i/j/k %s)"%wd)
sys.exit(0 if max(w[0],wf[0],wd[0])<=1e-10 else 1)
P
[ -n "${KEEP_LEVEL_CSV:-}" ] || rm -f r1/level_s*_r*.csv r2/level_s*_r*.csv
[ "$GROUND" = thick ] && { echo "REPORT open_step thick ground (L0Q, not gated): rc $rc"; exit 0; }
[ $rc = 0 ] && echo "PASS open_step" || { echo "FAIL open_step (see $RUN)"; exit 1; }
