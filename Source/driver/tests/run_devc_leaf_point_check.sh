#!/bin/bash
# run_devc_leaf_point_check.sh <fds_amr binary> [<run-dir>]: gas DEVCs on leaf cells (D-084). cases/devc_leaf_point.fds (two level-0 meshes, one level-1 patch over the whole
# domain) as 2 levels on 1 rank (--maxsize 16) and on 2 ranks (--maxsize 8, boxes on both ranks), and as 1 level at the level-1 cell size, 3 steps. FDSTL_LEVEL_CSV writes the
# finest level's fields after each step; every DEVC must equal the stock definition evaluated on those leaf fields to 1e-12: VOLUME INTEGRAL / MEAN / MAX / MIN over the
# leaf cells of the extent, and a point device at the leaf cell (TEMPERATURE) or the nearest leaf face (U-VELOCITY, read.f90 CELL_FACE rule).
# Old code: U1 (a point 0.1 m from XMIN) read the level-0 face (1.7e-2), and the VOLUME INTEGRAL counted the fine layer at the mesh seam twice (6e-2).
set -u
BIN=$1; RUN=${2:-$PWD/devc_leaf_point_run}; HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" && cd "$RUN" || exit 1
rc=0
for C in "2 1 16" "2 2 8" "1 2 8"; do
    set -- $C; L=$1; NP=$2; MS=$3; D=l${L}np$NP
    rm -rf $D; mkdir $D; cp "$HERE/cases/leaf_devc_foreign_vegetation_model.txt" $D/vegetation_model.txt
    if [ $L = 1 ]; then sed 's/IJK=4,8,4/IJK=8,16,8/' "$HERE/cases/devc_leaf_point.fds" > $D/c.fds; else cp "$HERE/cases/devc_leaf_point.fds" $D/c.fds; fi
    (cd $D && FDSTL_LEVEL_CSV=. OMP_NUM_THREADS=1 timeout 600 mpirun -np $NP "$BIN" c.fds --two-level-run --levels $L --patch 0 7 0 3 --maxsize $MS --nofine -1 --steps 3 --outdir . > log 2>&1) \
        || { echo "FAIL devc_leaf_point: $D run rc (see $RUN/$D/log)"; exit 1; }
    rm -f $D/*.s3d* $D/*.sf $D/*.bf $D/*.prt5
    python3 - $D <<'P' || rc=1
import csv,sys,glob
D=sys.argv[1]
a=list(csv.reader(open(D+'/devc_leaf_point_cat_devc.csv')));h=[x.strip() for x in a[1]]
dx=0.25;NI,NJ,NK=16,16,8;worst=(0.,'',0)
for s in (1,2,3):
    F={}
    for f in glob.glob(D+'/level_s%d_r*.csv'%s):
        for n,c,i,j,k,v in csv.reader(open(f)): F[(n,int(c),int(i),int(j),int(k))]=float(v)
    row=[float(x) for x in a[2+s]]
    cells=[(i,j,k) for i in range(NI) for j in range(NJ) for k in range(NK)]
    U=[F[('U',0,i+1,j,k)] for i,j,k in cells];V=[F[('V',0,i,j+1,k)] for i,j,k in cells];W=[F[('W',0,i,j,k+1)] for i,j,k in cells]
    E={'M':sum(F[('RHO',0)+c] for c in cells)*dx**3,'TM':sum(F[('TMP',0)+c] for c in cells)/len(cells)-273.15,'UM':sum(U)/len(U),'UMAX':max(U),'WMAX':max(W),'VMIN':min(V),
       'T1':F[('TMP',0,0,8,0)]-273.15,'U1':F[('U',0,round(0.1/dx),8,7)],'U2':F[('U',0,round(3.9/dx),0,2)]}
    for n in E:
        e=abs(row[h.index(n)]-E[n])/max(abs(E[n]),1e-300)
        if e>=worst[0]: worst=(e,n,s)
print("  %s: worst |DEVC - leaf field value| / |value| %.2e (%s, step %d)"%((D,)+worst)); sys.exit(0 if worst[0]<=1e-12 else 1)
P
    rm -f $D/level_s*_r*.csv
done
[ $rc = 0 ] && echo "PASS devc_leaf_point" || { echo "FAIL devc_leaf_point (see $RUN)"; exit 1; }
