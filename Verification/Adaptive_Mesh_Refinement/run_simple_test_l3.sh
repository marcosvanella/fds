#!/usr/bin/env bash
# simple_test fire case, 3-level adaptive (regrid every 5, tag-vort auto).
# Long gate run (r3lev): STEPS=2000. Short smoke: default STEPS=10.
: "${NP:=1}"
. "$(dirname "$0")/_common.sh"
run_case simple_test.fds simple_l3 --levels 3 --regrid-every 5 --tag-vort auto "$@"
