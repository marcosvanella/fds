#!/bin/bash
# run_nearwall_mu_check.sh <build-dir> [<run-dir>]: near-wall viscosity over the burner, two-level leaf vs one level at the same dx (0.5 m).
# cases/nearwall_mu_2l.fds (C064 cut-down c10, gas only, 1 m level 0, --bed over the whole domain: level 1 at 0.5 m everywhere) and cases/nearwall_mu_1l108.fds (the same case
# on one level at 0.5 m with the 108-box layout and rank map of that level 1), T_END 2 s, 4 ranks. Burner first-cell VISCOSITY devices MU1..MU3 (z 0.25 m) and the HRR-file Q_DIFF.
# Before the fix the fine ground ghost was the stress-matched value clipped to -U (no slip): WALE near-wall MU about 4x low, burner Q_DIFF -6.5 vs -24 kJ over [0,2] s.
# The ground-momentum defect behind the earlier expected failure: fine boxes took a uniform U_WIND = U0 (never the profile at their own ZC) and their t=0 velocity
# was prolonged from level 0, so the near-ground U on level 1 was ~30% high. Fixed: profile at the box ZC (setup and every VELOCITY_FLUX) and the stock IC on the fine level.
# Pass: both runs exit 0; time mean over [1,2] s of MU1+MU2+MU3 within 35% and Q_DIFF integral over [0,2] s within 25% (2L vs 1L).
BLD=$1; RUN=${2:-$BLD/run/nearwall_mu}; HERE=$(cd "$(dirname "$0")" && pwd)
rm -rf "$RUN"; mkdir -p "$RUN/2l" "$RUN/1l"
V="$HERE/../../../Validation/CSIRO_Grassland_Fires/FDS_Input_Files/vegetation_model.txt"
cp "$HERE/cases/nearwall_mu_2l.fds" "$RUN/2l/c.fds"; cp "$HERE/cases/nearwall_mu_1l108.fds" "$RUN/1l/c.fds"; cp "$V" "$RUN/2l/"; cp "$V" "$RUN/1l/"
(cd "$RUN/2l" && FDSTL_NO_MESH_DUMPS=1 OMP_NUM_THREADS=1 mpirun --oversubscribe -np 4 "$BLD/fds_amr" c.fds --two-level-run --bed -10 30 -20 20 0 20 > run.log 2>&1); r2=$?
(cd "$RUN/1l" && OMP_NUM_THREADS=1 mpirun --oversubscribe -np 4 "$BLD/fds_amr" c.fds --run > run.log 2>&1); r1=$?
echo "rc 2L $r2 1L $r1"; [ $r2 -eq 0 ] && [ $r1 -eq 0 ] || { echo "FAIL nearwall_mu"; exit 1; }
rm -f "$RUN"/*/*.s3d* "$RUN"/*/*.sf "$RUN"/*/*.bf "$RUN"/*/*.prt5*
python3 - "$RUN" <<'EOP'
import sys,numpy as np
def L(d,f):
    p=f"{sys.argv[1]}/{d}/c_cat_{f}.csv"; h=[x.strip().strip('"') for x in open(p).read().split('\n')[1].split(',')]
    return h,np.genfromtxt(p,delimiter=',',skip_header=2)
r={}
for d in ['2l','1l']:
    h,a=L(d,'devc'); m=(a[:,0]>=1)&(a[:,0]<=2); mu=sum(a[m,h.index(k)] for k in ('MU1','MU2','MU3')).mean()
    h,b=L(d,'hrr'); qd=np.trapezoid(b[:,h.index('Q_DIFF')],b[:,0]); r[d]=(mu,qd)
emu=r['2l'][0]/r['1l'][0]-1; eqd=r['2l'][1]/r['1l'][1]-1
print(f"burner first-cell MU sum mean [1,2] s: 2L {r['2l'][0]:.4e} 1L {r['1l'][0]:.4e} ({emu:+.1%}); Q_DIFF [0,2] s: 2L {r['2l'][1]:.2f} 1L {r['1l'][1]:.2f} kJ ({eqd:+.1%})")
ok=abs(emu)<0.35 and abs(eqd)<0.25
print("PASS nearwall_mu" if ok else "FAIL nearwall_mu"); sys.exit(0 if ok else 1)
EOP
