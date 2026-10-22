#!/bin/bash
# run_dt_parity_check.sh <build-dir> [<run-dir>]: time step of the whole-domain two-level run vs one level at the same dx, steps 1-10 within 1e-8 (relative).
# Cases: cases/nearwall_mu_2l.fds (level 1 at 0.5 m over the whole domain) and cases/nearwall_mu_1l108.fds (one level at 0.5 m, 108 meshes), T_END cut to 0.6 s, 4 ranks.
# Starting DT: read.f90 takes CFL_MAX*CHARACTERISTIC_CELL_SIZE/VEL_CHAR with the smallest cell of all meshes, so the two-level run starts from the finest level's (both 0.1087 s).
# EXPECTED FAILURE (open). Level 1 now starts from the stock t=0 state (run_t0_state_check.sh): worst rel 4.49e-5 (was 1.59e-3). With NOISE=.FALSE. in both
# cases 3.57e-6. The step-1 dt is noise sensitive at 1e-4: reordering the 108 MESH lines of the 1-level case into the level-1 box order (the INITIAL_NOISE seed is
# the mesh number) moves the 1-level dt by 2.0e-4, and the 2-level run then differs by 1.60e-4. Remaining suspects: the fine boxes have no wall heat flux in
# DIVERGENCE_PART_1 (hidden walls), and a noise-sensitive difference after t=0 not yet located.
# The tolerance is not to be relaxed; remove this mark only when it passes.
BLD=$1; RUN=${2:-$BLD/run/dt_parity}; HERE=$(cd "$(dirname "$0")" && pwd)
rm -rf "$RUN"; mkdir -p "$RUN/2l" "$RUN/1l"
V="$HERE/../../../Validation/CSIRO_Grassland_Fires/FDS_Input_Files/vegetation_model.txt"
for d in 2l 1l; do cp "$V" "$RUN/$d/"; done
sed 's/T_END=[0-9.]*/T_END=0.6/' "$HERE/cases/nearwall_mu_2l.fds" > "$RUN/2l/c.fds"; sed 's/T_END=[0-9.]*/T_END=0.6/' "$HERE/cases/nearwall_mu_1l108.fds" > "$RUN/1l/c.fds"
(cd "$RUN/2l" && FDSTL_NO_MESH_DUMPS=1 OMP_NUM_THREADS=1 mpirun --oversubscribe -np 4 "$BLD/fds_amr" c.fds --two-level-run --bed -10 30 -20 20 0 20 > run.log 2>&1); r2=$?
(cd "$RUN/1l" && FDSTL_NO_MESH_DUMPS=1 OMP_NUM_THREADS=1 mpirun --oversubscribe -np 4 "$BLD/fds_amr" c.fds --run > run.log 2>&1); r1=$?
rm -f "$RUN"/*/*.s3d* "$RUN"/*/*.sf "$RUN"/*/*.bf "$RUN"/*/*.prt5* "$RUN"/*/run_final_*.bin
echo "rc 2L $r2 1L $r1"; [ $r2 -eq 0 ] && [ $r1 -eq 0 ] || { echo "FAIL dt_parity"; exit 1; }
python3 - "$RUN" <<'EOP'
import sys,re,csv
R=sys.argv[1]
a=[float(m.group(1)) for m in re.finditer(r'TWO-LEVEL step \d+ t \S+ dt (\S+)',open(R+'/2l/run.log').read())][:10]
b=[float(r[2]) for r in list(csv.reader(open(R+'/1l/run_driver_steps.csv')))[1:11]]
w=0.0
for i,(x,y) in enumerate(zip(a,b)): e=abs(x-y)/y; w=max(w,e); print(f"step {i+1}: dt 2L {x:.10g} 1L {y:.10g} rel {e:.2e}")
ok=len(a)==10 and len(b)==10 and w<1e-8
print(f"worst rel {w:.2e}"); print("PASS dt_parity" if ok else "FAIL dt_parity"); sys.exit(0 if ok else 1)
EOP
