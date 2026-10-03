#!/bin/bash
# Interface flux hooks (notes/flux-hooks-design.md), host-side tests T1/T2/T3 and a negative control. Needs a build whose FDS tree carries patch 0009 (divg.f90) for the DIF part;
# without it only the ADV part is exercised and the script says so.
#   FDSTL_FLUXCHK=1  every density / divergence-1 call is preceded by a read-out and runs with an EMPTY override list (T1: result bitwise the unhooked one; the DIF pass runs
#                    DIVERGENCE_PART_1 twice, T3: the re-run is idempotent once the zone sums are restored)
#   FDSTL_FLUXCHK=2  the override lists carry the read-out value of EVERY face of every box (T2, no-op override: bitwise the unhooked result; tests the face index conversion and the
#                    second copy of the density loop; with 4 ranks the box-box faces are listed by both boxes)
#   FDSTL_FLUXCHK=3  as 2 with the values scaled by 1.001: the result must DIFFER (the hooks are effective), for ADV only, DIF only and both.
# usage: run_flux_hook_check.sh <build-dir> [steps] [work-dir]
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
BLD=$1; STEPS=${2:-10}; WORK=${3:-$BLD/run/fluxhook}; rm -rf "$WORK"; mkdir -p "$WORK"; rc=0
run() {  # name case np chk kinds
  local D="$WORK/$1"; mkdir -p "$D"; cd "$D"; cp "$HERE/cases/$2.fds" .
  FDSTL_FLUXCHK=$4 FDSTL_FLUXCHK_KINDS=$5 timeout 900 mpirun --bind-to none --oversubscribe -np $3 "$BLD/fds_amr" $2.fds --run --outdir . --chid $2 --steps $STEPS --exact-zone-sums > run.out 2> run.err || { echo "FAIL $1: run exit $?"; tail -3 run.err; rc=1; }
}
same() { local ok=1; for f in "$WORK/$1/$3"_final_*.bin; do cmp -s "$f" "$WORK/$2/$(basename $f)" || ok=0; done; cmp -s "$WORK/$1/$3_driver_steps.csv" "$WORK/$2/$3_driver_steps.csv" || ok=0; echo $ok; }
for T in "dec1 1" "dec4_np4 4"; do
  set -- $T; C=$1; NP=$2
  run ${C}_plain $C $NP 0 3
  for CHK in 1 2; do
    run ${C}_chk$CHK $C $NP $CHK 3
    if [ "$(same ${C}_plain ${C}_chk$CHK $C)" = 1 ]; then echo "PASS $C np=$NP FDSTL_FLUXCHK=$CHK: final fields and step log bitwise equal to the unhooked run ($(grep -a '^FLUXCHK' $WORK/${C}_chk$CHK/run.out | head -1 | cut -c1-150))"
    else echo "FAIL $C np=$NP FDSTL_FLUXCHK=$CHK: differs from the unhooked run"; rc=1; fi
  done
  for K in 1 2 3; do
    run ${C}_ctl$K $C $NP 3 $K
    if [ "$(same ${C}_plain ${C}_ctl$K $C)" = 0 ]; then echo "PASS $C np=$NP negative control (scaled override, kinds mask $K): result differs, the hooks act"
    else echo "NOTE $C np=$NP negative control kinds mask $K: result did NOT differ"; [ $K = 2 ] && echo "     (DIF part needs patch 0009 in the FDS tree)"; [ $K != 2 ] && rc=1; fi
  done
done
exit $rc
