#!/bin/bash
# D-056 option B, item (1): kernels on a fine-level box. Needs a build of the scratch tree with patches 0007 and 0008 and -DFDS_AMR_FINE_B_DRAFT=ON (see patches/0008-kernels-point-to-box.md).
# `--run --fine-b-shadow`: every kernel call on the level-0 box (single-mesh case that covers the whole domain) is repeated on a fine-level copy FINE_LEVEL(1)%BOX(1) (mesh number NMESHES+1,
# reached through POINT_TO_BOX) and the results are compared bit by bit; `--fine-b-build` does the same with a fine box that BUILD_FINE_BOX made from scratch (metrics, tables, arrays) and
# only the state arrays copied. The run itself is the unchanged level-0 run (its final fields are compared with the plain run: shadow must not change them).
# usage: run_fine_b_shadow_check.sh <build-dir> [steps] [work-dir]
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
BLD=$1; STEPS=${2:-20}; WORK=${3:-$BLD/run/fineshadow}; rm -rf "$WORK"; mkdir -p "$WORK"; rc=0
python3 "$HERE/make_periodic_case.py" 16 "$WORK/gen" 0.4 > /dev/null || exit 1
for CASE in dec1 tg16; do
  if [ -f "$HERE/cases/$CASE.fds" ]; then cp "$HERE/cases/$CASE.fds" "$WORK/"; else cp "$WORK/gen/$CASE.fds" "$WORK/gen/${CASE}_uvw.csv" "$WORK/"; fi
  MODES="plain --fine-b-shadow"; [ "$CASE" = tg16 ] && MODES="plain --fine-b-shadow --fine-b-build"   # dec1 has a one-cell y direction: every cell is next to a y face, there is no core to compare in build mode
  for MODE in $MODES; do
    D="$WORK/${CASE}_${MODE#--}"; mkdir -p "$D"; cd "$D"; cp "$WORK/$CASE".* . ; cp "$WORK"/${CASE}_uvw.csv . 2>/dev/null
    OPT=""; [ "$MODE" != plain ] && OPT=$MODE
    timeout 1500 mpirun --bind-to none --oversubscribe -np 1 "$BLD/fds_amr" $CASE.fds --run --outdir . --chid $CASE --steps $STEPS --exact-zone-sums $OPT > run.out 2> run.err; r=$?
    if [ "$MODE" = plain ]; then
      if [ $r -eq 0 ]; then echo "PASS $CASE plain run"; else echo "FAIL $CASE plain run (exit $r)"; rc=1; fi
    else
      if [ $r -eq 0 ] && grep -aq "FINE-B SHADOW PASS" $CASE.out; then echo "PASS $CASE $MODE: $(grep -ac "FINE-B SHADOW " $CASE.out) report lines, all kernels bitwise"; grep -a "FINE-B SHADOW " $CASE.out | head -16
      else echo "FAIL $CASE $MODE (exit $r)"; grep -a "FINE-B SHADOW\|MISMATCH" $CASE.out run.err | head -30; rc=1; fi
    fi
  done
  # the shadow must not change the level-0 results: final-field dumps and manifests of the shadowed runs are byte-identical to those of the plain run
  for MODE in $MODES; do [ "$MODE" = plain ] && continue; M=${MODE#--}; ok=1
    for f in "$WORK/${CASE}_plain/${CASE}"_final_*; do cmp -s "$f" "$WORK/${CASE}_$M/$(basename $f)" || ok=0; done
    if [ $ok = 1 ]; then echo "PASS $CASE level-0 final fields byte-identical with $MODE"; else echo "FAIL $CASE level-0 final fields differ with $MODE"; rc=1; fi
  done
done
exit $rc
