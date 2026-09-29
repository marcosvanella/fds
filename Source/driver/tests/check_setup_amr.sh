#!/bin/bash
# S1 smoke test: fds_amr (USE_AMREX=ON build) runs the FDS set-up and builds the level-0 BoxArray.
# usage: check_setup_amr.sh <build-dir-with-fds_amr> <case> <np> <expected-box-count>
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
BLD=$1; CASE=$2; NP=$3; NB=$4
RUN=$BLD/run/$CASE; rm -rf "$RUN"; mkdir -p "$RUN"; cp "$BASELINE/$CASE/$CASE.fds" "$RUN/"
(cd "$RUN" && timeout 300 mpirun --bind-to none -np "$NP" "$BLD/fds_amr" "$CASE.fds" > stdout.txt 2> stderr.txt) || { echo "FAIL: run"; exit 1; }
grep -q "level 0: $NB box(es)" "$RUN/stdout.txt" || { echo "FAIL: expected $NB boxes"; exit 1; }
grep -q "STOP: Set-up only" "$RUN/stderr.txt" || { echo "FAIL: no clean FDS stop"; exit 1; }
echo "PASS $CASE np=$NP: level-0 BoxArray with $NB box(es) built after the unchanged FDS set-up"
