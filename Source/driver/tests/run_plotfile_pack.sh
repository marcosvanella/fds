#!/bin/bash
# Show-and-tell plotfile pack: static R=2 --two-level-run on ns2d_16_l0 / blob_16_l0; optional --levels 3, --regrid-every, --plot-every.
# Fields: RHO,U,V,W. After write, points Integration smoke: /workspace/fds-amr/prototypes/s4_cuda_mass/stage1/plotfile_smoke/check_plotfile.sh
# usage: run_plotfile_pack.sh <build-dir> [work-dir] [steps=10]
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
SMOKE="/workspace/fds-amr/prototypes/s4_cuda_mass/stage1/plotfile_smoke"
B=$1; W=${2:-$PWD/plotpack}; ST=${3:-10}; rm -rf "$W"; mkdir -p "$W"; rc=0
run1() {
  local tag=$1 c=$2; shift 2
  mkdir -p "$W/$tag"; cp "$HERE/cases/$c.fds" "$W/$tag/"
  (cd "$W/$tag" && timeout 900 mpirun --bind-to none --oversubscribe -np 1 "$B/fds_amr" "$c.fds" --two-level-run --steps "$ST" --patch 4 11 4 11 --blocking 4 --maxsize 16 \
      --outdir . --chid "$tag" --plotfile . "$@" > run.log 2>&1) || { echo "FAIL $tag"; grep -E "Abort|error|STOP|PLOTFILE|mid-regrid" "$W/$tag/run.log" | head -30; return 1; }
  grep -q "PLOTFILE wrote" "$W/$tag/run.log" || { echo "FAIL $tag: no PLOTFILE line"; return 1; }
  local plt; plt=$(grep "PLOTFILE wrote" "$W/$tag/run.log" | tail -1 | awk '{print $3}')
  # resolve relative path written under case dir
  [ -d "$plt" ] || plt="$W/$tag/$(basename "$plt")"
  if [ -x "$SMOKE/check_plotfile.sh" ]; then
    "$SMOKE/check_plotfile.sh" "$plt" || { echo "FAIL $tag: Integration check_plotfile"; return 1; }
  else
    [ -f "$plt/Header" ] || { echo "FAIL $tag: no Header under $plt"; return 1; }
  fi
  echo "PASS $tag -> $plt"
}
run1 ns2d_static ns2d_16_l0 || rc=1
run1 blob_static blob_16_l0 || rc=1
run1 ns2d_l3 ns2d_16_l0 --levels 3 --patch2 6 9 6 9 || rc=1
run1 ns2d_regrid ns2d_16_l0 --regrid-every 5 || rc=1
run1 ns2d_cadence ns2d_16_l0 --plot-every 2 || rc=1
echo "Pack under $W. Integration: $SMOKE/check_plotfile.sh <plt>; $SMOKE/paraview_smoke.py <plt>"
exit $rc
