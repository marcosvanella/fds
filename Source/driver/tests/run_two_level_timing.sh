#!/bin/bash
# Timing of the two-level driver against the single-level driver (notes/two-level-timing.md). Per case (ns2d, blob) and coarse size N (16, 32 ...), at 1 rank (one mesh) and 4 ranks
# (four meshes): single-level coarse (--run), two-level (fine patch over the middle half), REF two-level (fine patch = whole domain) and single-level uniform fine (--run on 2N).
# Stage seconds from <chid>_driver_profile.txt (FDSTL_PROFILE), best of REPS runs by t_loop. Writes timing.csv in the work dir.
# usage: run_two_level_timing.sh <build-dir> [work-dir] [sizes="16 32"] [steps=20] [reps=3]
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
B=$1; W=${2:-$PWD/timing}; SIZES=${3:-"16 32"}; ST=${4:-20}; REPS=${5:-3}
mkdir -p "$W"; cd "$W"; export FDSTL_PROFILE=1
mkcase() { # base N meshes(1|4) -> file name; the l0 case with N x 1 x N cells, as one mesh or as 2 x 2 meshes, T_END large (the step count decides)
  local b=$1 n=$2 m=$3 f=${1}_${2}_${3}m.fds
  python3 - "$HERE/cases/${b}_16_l0.fds" "$n" "$m" "$f" <<'PY'
import re,sys
src,n,m,out=sys.argv[1],int(sys.argv[2]),int(sys.argv[3]),sys.argv[4]
s=open(src).read()
mm=re.search(r"&MESH IJK=16,1,16, XB=([^/]*)/",s); x0,x1,y0,y1,z0,z1=[float(v) for v in mm.group(1).split(',') if v.strip()]
if m==1: mesh="&MESH IJK=%d,1,%d, XB=%r,%r,%r,%r,%r,%r /"%(n,n,x0,x1,y0,y1,z0,z1)
else:
    xh,zh=(x0+x1)/2,(z0+z1)/2; h=n//2
    mesh="\n".join("&MESH IJK=%d,1,%d, XB=%r,%r,%r,%r,%r,%r, MPI_PROCESS=%d /"%(h,h,a,b,y0,y1,c,d,k) for k,(c,d,a,b) in enumerate([(z0,zh,x0,xh),(z0,zh,xh,x1),(zh,z1,x0,xh),(zh,z1,xh,x1)]))
s=s.replace(mm.group(0),mesh); s=re.sub(r"T_END=[0-9.]+","T_END=1000.0",s); s=re.sub(r"CHID='[^']*'","CHID='%s'"%out[:-4],s)
open(out,"w").write(s)
PY
  echo $f
}
one() { # tag np case args... : best of REPS by t_loop; appends tag,np,ncoarse,nfine,profile values to timing.csv
  local tag=$1 np=$2 c=$3; shift 3; local best=""
  for r in $(seq $REPS); do
    rm -rf out_$tag; mkdir -p out_$tag
    timeout 3000 mpirun --bind-to none --oversubscribe -np $np "$B/fds_amr" $c "$@" --steps $ST --outdir out_$tag --chid $tag > $tag.log 2>&1
    local p=out_$tag/${tag}_driver_profile.txt; [ -f $p ] || { echo "FAIL $tag (no profile, see $tag.log)"; return 1; }
    local tl; tl=$(awk '$1=="t_loop"{print $2}' $p)
    if [ -z "$best" ] || python3 -c "import sys; sys.exit(0 if $tl < $best else 1)"; then best=$tl; cp $p best_$tag.txt; cp out_$tag/${tag}_driver_perf.csv best_${tag}_perf.csv; fi
  done
  echo "$tag,$np,$(paste -sd, <(awk '{print $2}' best_$tag.txt)),$(tail -1 best_${tag}_perf.csv | cut -d, -f6)" >> timing.csv
  echo "$tag np$np t_loop $best s"
}
echo "tag,np,t_loop,t_pressure,t_poisson_solve,t_fds_outputs,ghost_fill_boundary,fill_omesh,match_periodic_faces,fds_bc_routines,avg_down,cf_ghost_fill,flux_overwrite,steps" > timing.csv
for b in ns2d blob; do
  for n in $SIZES; do
    q=$((n / 4)); h=$((3 * n / 4 - 1)); e=$((n - 1))
    for np in 1 4; do
      m=1; [ $np = 4 ] && m=4
      c=$(mkcase $b $n $m); f=$(mkcase $b $((2 * n)) $m)
      one ${b}${n}_np${np}_single $np $c --run
      one ${b}${n}_np${np}_two $np $c --two-level-run --patch $q $h $q $h
      one ${b}${n}_np${np}_ref $np $c --two-level-run --patch 0 $e 0 $e
      one ${b}${n}_np${np}_fine $np $f --run
    done
  done
done
