#!/usr/bin/env bash
# simple_test fire case, 2-level adaptive (regrid every 5, tag-vort auto).
# Long gate run (r5_wbc): STEPS=900 covers T_END=5. Short smoke: default STEPS=10.
: "${NP:=1}"
. "$(dirname "$0")/_common.sh"
run_case simple_test.fds simple_l2 --levels 2 --regrid-every 5 --tag-vort auto "$@"
