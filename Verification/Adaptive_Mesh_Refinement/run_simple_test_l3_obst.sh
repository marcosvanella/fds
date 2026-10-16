#!/usr/bin/env bash
# Fast 3-level OBST/wall smoke on the 36x24x24 simple_test (burner on an OBST at level 0): levels refined at the OBST and the hot plume (--tag-obst, --tag-temp),
# regrid every 5 steps, default 30 steps (about 1 min on one core). Passes when every step reaches the end and max|div u - D| <= 1e-9 on every level.
# Usage: FDS_AMR=.../fds_amr ./run_simple_test_l3_obst.sh   (STEPS=n to change the length; the end-of-run mass gate is not part of this smoke check)
: "${NP:=1}"
: "${STEPS:=30}"
. "$(dirname "$0")/_common.sh"
run_case simple_test.fds simple_l3_obst --levels 3 --regrid-every 5 --tag-obst --nofine -1 --tag-temp 400 --blocking 4 --maxsize 8 --log-every 1 "$@" || true
log=${RUN_DIR:-$HERE/out_simple_l3_obst}/simple_l3_obst.log
n=$(grep -c "^MULTI-LEVEL step " "$log")
bad=$(awk '/^MULTI-LEVEL step /{for(i=1;i<=NF;i++) if($i=="D|"){ if ($(i+1)+0 > 1e-9) c++ }} END{print c+0}' "$log")
if [ "$n" -eq "$STEPS" ] && [ "$bad" -eq 0 ]; then echo "PASS simple_test_l3_obst: $n steps, max|div u - D| <= 1e-9"; exit 0; fi
echo "FAIL simple_test_l3_obst: $n of $STEPS steps, $bad steps with div > 1e-9 ($log)"; exit 1
