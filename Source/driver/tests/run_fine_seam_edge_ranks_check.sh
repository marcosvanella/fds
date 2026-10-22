#!/bin/bash
# run_fine_seam_edge_ranks_check.sh <fds_amr binary> [<run-dir>]: fine-level boxes on 1 rank vs 2 ranks must give bit-identical results.
# cases/fine_seam_edge_ranks.fds: CSIRO C064 level-set crude as 4 L0 meshes (2x2; MPI_PROCESS 0,0,1,1 at 2 ranks), 2 levels, --ls-fine --nofine 0, regrid every 2 steps,
# 30 s (the lsfine_multirank_gate L2_nf0_m4 rows). The fine boxes have hidden OPEN domain walls; the seams between boxes cross the domain faces.
# Guards 2c12ae8401: a regrid handed the old layout's vent faces (vent_old) to the new layout rank-locally and with the OLD box index, so fine_wall_fv drove the
# face of an unrelated box (or none, when that index was on another rank): the composite RHS of the first step after a layout change (step 59) depended on the
# rank count. Also covers the domain edges of hidden OPEN walls at box seams (b1a806cfdc, Integration's rule: edges on any domain face).
# Pass: same level-1 layout at every regrid (TWO-LEVEL mid-regrid lines) and c_hrr.csv, c_devc.csv, c_multi_level.csv identical byte for byte, except the multi_level
# columns mass, mass_rel_change, inflow (rank-partial sums: rank tolerance, lsfine_multirank_gate).
set -u
BIN=$1; RUN=${2:-$PWD/fine_seam_edge_ranks_run}; HERE=$(cd "$(dirname "$0")" && pwd); MPIRUN=${MPIRUN:-mpirun}
TF="--two-level-run --levels 2 --init-patch --patch 5 9 0 11 --blocking 2 --maxsize 8 --regrid-every 2 --tag-phi --ls-fine --plotfile pf --plot-every 1000000 --nofine 0"
for r in 1 2; do
    d=$RUN/r$r; rm -rf "$d"; mkdir -p "$d"
    if [ $r = 1 ]; then sed 's/, MPI_PROCESS=[0-9]*//' "$HERE/cases/fine_seam_edge_ranks.fds" > "$d/c.fds"; else cp "$HERE/cases/fine_seam_edge_ranks.fds" "$d/c.fds"; fi
    (cd "$d" && OMP_NUM_THREADS=1 timeout 900 $MPIRUN -np $r "$BIN" c.fds $TF --outdir . --chid c > run.log 2>&1) || { echo "FAIL fine_seam_edge_ranks: $r-rank run rc (see $d/run.log)"; exit 1; }
    rm -rf "$d"/pf* "$d"/*.s3d* "$d"/*.sf "$d"/*.bf "$d"/*.prt5
done
python3 - "$RUN/r1" "$RUN/r2" <<'P'
import csv,os,re,sys
a,b=sys.argv[1:3]; bad=[]
lay=lambda d:[re.sub(r' t .*?changed','changed',re.sub(r'leaf fuel.*','',l)) for l in open(os.path.join(d,'run.log'),errors='replace') if 'mid-regrid' in l]
if lay(a)!=lay(b): bad.append('level-1 layouts differ')
for f in ('c_hrr.csv','c_devc.csv','c_multi_level.csv'):
    A=list(csv.reader(open(os.path.join(a,f)))); B=list(csv.reader(open(os.path.join(b,f))))
    h=A[0] if f=='c_multi_level.csv' else A[1]
    if len(A)!=len(B): bad.append(f+': row count'); continue
    for i,(x,y) in enumerate(zip(A,B)):
        d=[h[j] for j,(p,q) in enumerate(zip(x,y)) if p!=q and h[j] not in ('mass','mass_rel_change','inflow')]
        if d: bad.append(f'{f}: first diff row {i} {d[0]}: {x[h.index(d[0])]} vs {y[h.index(d[0])]}'); break
print('\n'.join('  '+s for s in bad) if bad else f'  {len(lay(a))} regrids, csv identical')
sys.exit(1 if bad else 0)
P
[ $? -eq 0 ] && echo "PASS fine_seam_edge_ranks" || { echo "FAIL fine_seam_edge_ranks (see $RUN)"; exit 1; }
