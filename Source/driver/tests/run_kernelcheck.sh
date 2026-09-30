#!/bin/bash
# S3 kernel checks (USE_AMREX=ON build, threads=1): the wrapped, unmodified FDS kernels against frozen-input dumps
# written by the scratch reference generator (tests/refdump, unpatched FDS + write-only hooks).
# usage: run_kernelcheck.sh <build-dir> <ref-runs-dir> [work-dir]
#   <ref-runs-dir>/<case>/ holds <case fds file> and ref.dump (4-mesh: ref.dump, ref.dump.2, .3, .4)
#   cases: shunn3_32, csmag_32, shunn3_32_clip, shunn3_4mesh_32__1mesh, shunn3_4mesh_32 (native, 4 ranks + window of the 1-mesh dump, with and without the clip active
#   (shunn3_4mesh_32__1mesh_clip): the 4-box gather clip against the single-mesh clip)
# Every kernel tag must print BITWISE-OK; any DIFFER makes the script fail.
# ghost modes: dump = ghosts from the dump; full+bc / face+bc = level ghost fill (full or face-neighbour-only) with the boundary-condition-step ghost values
# (corner strips, FDS WALL_BC) taken from the dump.
# full / face (plain, S4): no boundary-condition value from the dump; the driver runs FDS's own VISCOSITY_BC/VELOCITY_BC on the frozen state (GhostExchange).
# A frozen snapshot cannot give the exact pre-boundary-step state (the boundary-face values are rewritten between the boundary step and the dump point), so
# a handful of elements differ in the last bits (README, "Plain ghost modes"; with the dump's boundary-face strips, the +strips runs, all tags are bitwise): the plain gate is (1) every BCCHAIN tag (the real chain, from the pre-boundary
# state) BITWISE-OK, (2) every other tag BITWISE-OK or with at most PLAIN_MAX_PPM parts per million of its elements differing, (3) nothing else differs.
set -u
PLAIN_MAX_PPM=${PLAIN_MAX_PPM:-100}
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
  [ "$nok" -gt 0 ] && { [ "$nbad" = 0 ] || [ "$STRICT" = 0 ]; } || rc=1
}
STRICT=1
plain() { # like run, plus the plain-mode gate on the per-tag counts
  local label=$1 np=$2
  STRICT=0; run "$@"; STRICT=1
  local d="$WORK/$label" res
  res=$(grep -aE "^  [A-Z0-9_]+ +(BITWISE-OK|DIFFER)" "$d/out.log" | sed 's/:/ /; s/,//g' | awk -v ppm="$PLAIN_MAX_PPM" '
    { tag=$1; st=$2; el=$8; bd=$10; tot+=bd
      if (st=="DIFFER" && (tag ~ /^BCCHAIN/ || bd*1000000 > ppm*el)) { bad=1; printf "  plain gate FAIL: %s %d of %d\n", tag, bd, el } }
    END { if (bad) print "FAIL"; else printf "PASS (%d bit differences in total, all within tolerance)\n", tot }')
  echo "[$label np=$np] plain gate: $res"
  case "$res" in *FAIL*) rc=1;; esac
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
for g in full face; do
  plain "shunn3_32.$g"            1 shunn3_32 shunn3_32.fds "$REF/shunn3_32/ref.dump" --ghost=$g
  plain "csmag_32.$g"             1 csmag_32 csmag_32.fds "$REF/csmag_32/ref.dump" --ghost=$g
  plain "shunn3_32_clip.$g"       1 shunn3_32_clip shunn3_32.fds "$REF/shunn3_32_clip/ref.dump" --ghost=$g
  plain "1mesh.$g"                1 shunn3_4mesh_32__1mesh shunn3_4mesh_32.fds "$REF/shunn3_4mesh_32__1mesh/ref.dump" --ghost=$g
  plain "4mesh_native.$g"         4 shunn3_4mesh_32 shunn3_4mesh_32.fds "$REF/shunn3_4mesh_32/ref.dump" --ghost=$g
  plain "4mesh_window.$g"         4 shunn3_4mesh_32 shunn3_4mesh_32.fds "$REF/shunn3_4mesh_32__1mesh/ref.dump" --window --ghost=$g
  plain "4mesh_window_clip.$g"    4 shunn3_4mesh_32 shunn3_4mesh_32.fds "$REF/shunn3_4mesh_32__1mesh_clip/ref.dump" --window --ghost=$g
done
# plain+strips (S6): the plain modes with the boundary-face strips of U/V/W/US/VS/WS taken from the dump after the replay (FDSKC_STRIPDUMP=1). The strips are the only
# part of the state that a frozen snapshot cannot give exactly (see README "Plain ghost modes"); with them every tag is required BITWISE (ppm gate 0).
SAVE_PPM=$PLAIN_MAX_PPM; PLAIN_MAX_PPM=0
for g in full face; do
  export FDSKC_STRIPDUMP=1
  plain "shunn3_32.$g+strips"            1 shunn3_32 shunn3_32.fds "$REF/shunn3_32/ref.dump" --ghost=$g
  plain "csmag_32.$g+strips"             1 csmag_32 csmag_32.fds "$REF/csmag_32/ref.dump" --ghost=$g
  plain "shunn3_32_clip.$g+strips"       1 shunn3_32_clip shunn3_32.fds "$REF/shunn3_32_clip/ref.dump" --ghost=$g
  plain "1mesh.$g+strips"                1 shunn3_4mesh_32__1mesh shunn3_4mesh_32.fds "$REF/shunn3_4mesh_32__1mesh/ref.dump" --ghost=$g
  plain "4mesh_native.$g+strips"         4 shunn3_4mesh_32 shunn3_4mesh_32.fds "$REF/shunn3_4mesh_32/ref.dump" --ghost=$g
  plain "4mesh_window.$g+strips"         4 shunn3_4mesh_32 shunn3_4mesh_32.fds "$REF/shunn3_4mesh_32__1mesh/ref.dump" --window --ghost=$g
  plain "4mesh_window_clip.$g+strips"    4 shunn3_4mesh_32 shunn3_4mesh_32.fds "$REF/shunn3_4mesh_32__1mesh_clip/ref.dump" --window --ghost=$g
  unset FDSKC_STRIPDUMP
done
PLAIN_MAX_PPM=$SAVE_PPM
[ $rc = 0 ] && echo "KERNEL CHECKS PASS: dump, full+bc, face+bc BITWISE; plain full, face: BCCHAIN bitwise, other tags within the PLAIN_MAX_PPM gate; plain+strips (exact boundary-face strips): all tags BITWISE" || echo "KERNEL CHECKS: SOME DIFFER (see out.log files)"
exit $rc
