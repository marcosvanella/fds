#!/bin/bash
# Single-level equivalence (plan section 6): an input with `&AMR MAX_LEVEL=0` must give BYTE-IDENTICAL results to the same input without any &AMR line, on the three M2a
# cases (shunn3_32 on 1 rank, shunn3_4mesh_32 on 1 and 4 ranks, csmag_32 on 1 rank). Compared: every <chid>_final_*.bin and <chid>_driver_steps.csv.
# Today the driver does not read &AMR at all (FDS ignores the namelist), so this checks that adding the line changes nothing; once the AMR front end of the driver is in, the
# same script is the real max_level=0 test (nothing in it depends on how the line is consumed).
# usage: run_single_level_equiv.sh <driver-build> [steps] [work]    Exit 0 pass, 1 fail, 77 skipped (driver binary or baseline inputs missing).
set -u
BLD=${1:-}; STEPS=${2:-6}; WORK=${3:-/tmp/rt_single_level_equiv}
HERE=$(cd "$(dirname "$0")" && pwd)
DRVTESTS=$HERE/../../driver/tests
if [ -z "$BLD" ] || [ ! -x "$BLD/fds_amr" ]; then echo "SKIP: driver executable fds_amr not found (give the driver build directory as argument 1)"; exit 77; fi
source "$DRVTESTS/env.sh" > /dev/null 2>&1
if [ ! -d "$BASELINE/shunn3_32" ]; then echo "SKIP: baseline inputs not found under BASELINE=$BASELINE"; exit 77; fi
rm -rf "$WORK"; mkdir -p "$WORK"; rc=0
run() { # dir case np amrline(0/1)
  local D="$WORK/$1"; mkdir -p "$D"; cp "$BASELINE/$2/$2.fds" "$D/"; for c in "$BASELINE/$2"/*.csv; do [ -f "$c" ] && case $c in *_hrr.csv|*_mass.csv|*_cpu.csv|*_steps.csv|*_devc.csv|*_mms.csv|*_uvw_t*_m*.csv) ;; *) cp "$c" "$D/";; esac; done
  if [ "$4" = 1 ]; then sed -i 's|^&TAIL|\&AMR MAX_LEVEL=0 /\n\&TAIL|' "$D/$2.fds"; grep -q '^&AMR MAX_LEVEL=0' "$D/$2.fds" || { echo "FAIL: could not insert the &AMR line into $2.fds"; rc=1; }; fi
  (cd "$D" && timeout 900 mpirun --bind-to none --oversubscribe -np $3 "$BLD/fds_amr" $2.fds --run --steps $STEPS --outdir . --chid $2 --quiet --exact-zone-sums > run.out 2> run.err) || { echo "FAIL: run $1 exited with an error"; tail -3 "$D/run.err"; rc=1; }
}
for T in "shunn3_32 1" "shunn3_4mesh_32 1" "shunn3_4mesh_32 4" "csmag_32 1"; do
  set -- $T
  run "$1_np$2_plain" $1 $2 0; run "$1_np$2_amr0" $1 $2 1
  n=0; bad=0
  for f in "$WORK/$1_np$2_plain/$1"_final_*.bin "$WORK/$1_np$2_plain/$1_driver_steps.csv"; do
    [ -f "$f" ] || continue; n=$((n+1)); b=$(basename "$f")
    cmp -s "$f" "$WORK/$1_np$2_amr0/$b" || { bad=$((bad+1)); echo "  differs: $b"; }
  done
  if [ $n -gt 0 ] && [ $bad = 0 ]; then echo "PASS $1 np=$2: $n files byte-identical with and without '&AMR MAX_LEVEL=0'"; else echo "FAIL $1 np=$2: $bad of $n files differ (or none found)"; rc=1; fi
done
exit $rc
