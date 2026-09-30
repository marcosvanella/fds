#!/bin/bash
# S3 kernel checks (USE_AMREX=ON build, threads=1): the wrapped, unmodified FDS kernels against frozen-input dumps
# written by the scratch reference generator (tests/refdump, unpatched FDS + write-only hooks).
# usage: run_kernelcheck.sh <build-dir> <ref-runs-dir> [work-dir]
#   <ref-runs-dir>/<case>/ holds <case fds file> and ref.dump (4-mesh: ref.dump, ref.dump.2, .3, .4)
#   cases: shunn3_32, csmag_32, shunn3_32_clip, shunn3_4mesh_32__1mesh, shunn3_4mesh_32 (native, 4 ranks + window of the 1-mesh dump, with and without the clip active
#   (shunn3_4mesh_32__1mesh_clip): the 4-box gather clip against the single-mesh clip)
# Every kernel tag must print BITWISE-OK; any DIFFER makes the script fail.
# ghost modes: dump = ghosts from the dump; full+bc / face+bc = level ghost fill (full or face-neighbour-only) with the boundary-condition-step ghost values
# (corner strips, FDS WALL_BC) taken from the dump. Plain --ghost=full (no +bc) is not gated: it differs by ~1e-16 in 3 elements next to the domain corners (README).
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
BLD=$1; REF=$2; WORK=${3:-$BLD/run/kernelcheck}
mkdir -p "$WORK"; rc=0
run() { # label np case fds dump extra...
  local label=$1 np=$2 case=$3 fds=$4 dump=$5; shift 5
  local d="$WORK/$label"; rm -rf "$d"; mkdir -p "$d"; cp "$REF/$case/$fds" "$d/"; cp "$REF/$case"/*.csv "$d/" 2>/dev/null
  (cd "$d" && timeout 900 mpirun --bind-to none -np $np "$BLD/fds_amr" "$fds" --kernelcheck "$dump" "$@" > out.log 2>&1)
  local nok nbad; nok=$(grep -ac "BITWISE-OK" "$d/out.log"); nbad=$(grep -ac "DIFFER" "$d/out.log")
  echo "[$label np=$np] BITWISE-OK tags: $nok  DIFFER tags: $nbad"
  grep -a "DIFFER" "$d/out.log" | head -5
  [ "$nok" -gt 0 ] && [ "$nbad" = 0 ] || rc=1
}
for g in dump full+bc face+bc; do
  run "shunn3_32.$g"            1 shunn3_32 shunn3_32.fds "$REF/shunn3_32/ref.dump" --ghost=$g
  run "csmag_32.$g"             1 csmag_32 csmag_32.fds "$REF/csmag_32/ref.dump" --ghost=$g
  run "shunn3_32_clip.$g"       1 shunn3_32_clip shunn3_32.fds "$REF/shunn3_32_clip/ref.dump" --ghost=$g
  run "1mesh.$g"                1 shunn3_4mesh_32__1mesh shunn3_4mesh_32.fds "$REF/shunn3_4mesh_32__1mesh/ref.dump" --ghost=$g
  run "4mesh_native.$g"         4 shunn3_4mesh_32 shunn3_4mesh_32.fds "$REF/shunn3_4mesh_32/ref.dump" --ghost=$g
  run "4mesh_window.$g"         4 shunn3_4mesh_32 shunn3_4mesh_32.fds "$REF/shunn3_4mesh_32__1mesh/ref.dump" --window --ghost=$g
  run "4mesh_window_clip.$g"    4 shunn3_4mesh_32 shunn3_4mesh_32.fds "$REF/shunn3_4mesh_32__1mesh_clip/ref.dump" --window --ghost=$g
done
[ $rc = 0 ] && echo "ALL KERNEL CHECKS BITWISE" || echo "KERNEL CHECKS: SOME DIFFER (see out.log files)"
exit $rc
