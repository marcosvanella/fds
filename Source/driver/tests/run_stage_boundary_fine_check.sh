#!/bin/bash
# stage_boundary(1,3) and stage_boundary(1,6) on a bound fine level (the boundary routines FDS runs after the exchanges of the predictor and the corrector: fds_p_save_uvw, fill_omesh,
# VELOCITY_BC). Mode `fds_amr ns2d_16_l0.fds --two-level-run --stage-boundary-test [legacy]`. Three runs, all with FDSTL_SKIPAFT=0 (nothing skipped on the fine level):
#   1. fixed code: both calls return and level 0 U/V/W (ghost layers included) are bit-for-bit unchanged                                  -> must PASS
#   2. `legacy`: fds_p_save_uvw with the guard it had before S14.1 (level-0 mesh numbers only): the call on a fine box must stop the run   -> must abort with the guard message
#   3. FDSTL_LEGACY_BCSTEP_OFFSET=1: the BcStep of the fine level calls the FDS routines with box index + 1 (the offset 0 it had before S14.4), i.e. on level-0 mesh objects:
#      VELOCITY_BC overwrites the periodic ghost layers of level 0 U and W                                                               -> the check must FAIL (reproducer, notes/fine-velocity-bc.md)
# usage: run_stage_boundary_fine_check.sh <build-dir> [work-dir]
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
BLD=$1; WORK=${2:-$BLD/run/stageboundary}; rm -rf "$WORK"; mkdir -p "$WORK"; cp "$HERE/cases/ns2d_16_l0.fds" "$WORK/"; rc=0
cd "$WORK"
run() { FDSTL_SKIPAFT=0 timeout 600 mpirun --bind-to none --oversubscribe -np 1 "$BLD/fds_amr" ns2d_16_l0.fds --two-level-run --outdir . --chid sb "$@" > "$OUT.out" 2> "$OUT.err"; }
OUT=fixed; run --stage-boundary-test; r=$?
if [ $r -eq 0 ] && grep -q "STAGE-BOUNDARY-FINE: returned, level 0 U/V/W max change 0 (PASS)" fixed.out; then echo "PASS stage_boundary(1,3|6) on the fine level returns, level 0 untouched"; else echo "FAIL fixed code (exit $r)"; rc=1; fi
OUT=legacy; run --stage-boundary-test legacy; r=$?
if [ $r -ne 0 ] && grep -q "fds_amr ERROR in fds_p_save_uvw: mesh number" legacy.err && ! grep -q "STAGE-BOUNDARY-FINE: returned" legacy.out; then echo "PASS control: the pre-S14.1 guard aborts on the fine box (exit $r)"; else echo "FAIL legacy-guard control (exit $r)"; rc=1; fi
OUT=offset0; FDSTL_LEGACY_BCSTEP_OFFSET=1 run --stage-boundary-test; r=$?
if [ $r -ne 0 ] && grep -q "(FAIL)" offset0.out; then echo "PASS control: with mesh offset 0 the fine-level call overwrites level 0 U/W ghosts (reproducer): $(grep -a 'U/V/W max change' offset0.out | sed 's/.*change //')"; else echo "FAIL offset-0 reproducer (exit $r)"; rc=1; fi
exit $rc
