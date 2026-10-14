#!/usr/bin/env bash
# simple_test fire case, level 0 only (uniform). Single mesh -> NP=1.
# Long gate run: STEPS=900 (covers T_END=5). Short smoke: default STEPS=10.
: "${NP:=1}"
. "$(dirname "$0")/_common.sh"
run_case simple_test.fds simple_l0 --levels 1 "$@"
