#!/bin/bash
# S2 driver tests (USE_AMREX=ON build, threads=1, gfortran + Open MPI).
# usage: [FDS_INVENTORY_CSV=<mesh_fields.csv>] run_driver_tests.sh <build-dir> [--thread-sweep]
#   0. (optional, when FDS_INVENTORY_CSV is set) the field table against the inventory bounds
#   1. driver_unit_tests on 1 and 4 ranks (IR-005 index round trip, ghosts, fill, side data, registry, IR-007 skeleton)
#   2. fds_amr --selftest on shunn3_32 (1 rank) and shunn3_4mesh_32 (4 ranks): FDS bounds vs table, alias read/write, side data from FDS cells
#   3. the side-data hash of the 1-mesh and the 4-mesh runs must be equal (same 32^3 domain, layout independent)
#   4. run_stage_boundary_fine_check.sh (stage_boundary(1,3|6) on a bound fine level, with its two negative controls)
#   5. (optional, DRIVER_DENSITY_CHECK=1; about 25 min on one shared core) run_two_level_density_check.sh: variable-density two-level run, flux overwrite ON vs OFF, 1 and 4 ranks
# FDS ends the set-up run with STOP (exit 0), so the FDS-linked checks are judged from their PASS/FAIL lines.
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
BLD=$1; SW=${2:-}
rc=0
INV=${FDS_INVENTORY_CSV:-}
if [ -n "$INV" ]; then python3 "$HERE/check_inventory.py" "$BLD/driver_unit_tests" "$INV" || rc=1; fi
for np in 1 4; do
  out=$(mpirun --bind-to none -np $np "$BLD/driver_unit_tests" $SW 2>&1)
  echo "$out" | grep -E "^(PASS|FAIL|ALL PASS|SOME FAILED)" | sed "s/^/[unit np=$np] /"
  echo "$out" | grep -q "^ALL PASS" || { echo "[unit np=$np] FAIL"; rc=1; }
done
h=()
for spec in "shunn3_32 1" "shunn3_4mesh_32 4"; do
  set -- $spec; CASE=$1; NP=$2
  RUN=$BLD/run/selftest_$CASE; rm -rf "$RUN"; mkdir -p "$RUN"; cp "$BASELINE/$CASE/$CASE.fds" "$RUN/"
  (cd "$RUN" && timeout 600 mpirun --bind-to none -np $NP "$BLD/fds_amr" "$CASE.fds" --selftest > stdout.txt 2> stderr.txt)
  grep -E "^(PASS|FAIL|SELFTEST)|INFO|SIDEDATA_HASH" "$RUN/stdout.txt" | sed "s/^/[$CASE np=$NP] /"
  grep -q "^SELFTEST PASS" "$RUN/stdout.txt" && ! grep -q "^FAIL" "$RUN/stdout.txt" || { echo "[$CASE np=$NP] FAIL"; rc=1; }
  h+=("$(grep -o 'SIDEDATA_HASH [0-9a-f]*' "$RUN/stdout.txt" | head -1)")
done
if [ -n "${h[0]}" ] && [ "${h[0]}" = "${h[1]}" ]; then echo "PASS side-data hash equal for 1 mesh and 4 meshes: ${h[0]}"; else echo "FAIL side-data hash differs: '${h[0]}' vs '${h[1]}'"; rc=1; fi
bash "$HERE/run_stage_boundary_fine_check.sh" "$BLD" "$BLD/run/stageboundary" || rc=1
if [ "${DRIVER_DENSITY_CHECK:-0}" = 1 ]; then bash "$HERE/run_two_level_density_check.sh" "$BLD" "$BLD/run/density" || rc=1; fi
[ $rc = 0 ] && echo "ALL DRIVER TESTS PASS" || echo "DRIVER TESTS FAILED"
exit $rc
