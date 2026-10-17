#!/bin/bash
# Flat WUI level-set cases (Verification/WUI): stock FDS vs fds_amr unrefined (--run), optional refined (--two-level-run --init-patch ... --tag-phi).
# usage: run_wui_levelset_compare.sh STOCK_FDS FDS_AMR OUTDIR [CASES]   (default cases: the flat ones below; REFINED=1 also runs the refined variant)
# prints, per case: stock vs amr final HRR (kW), final value of every device of the case, max relative difference; exit 0 if every unrefined case agrees to TOL (default 1e-6).
S=$(readlink -f "$1"); A=$(readlink -f "$2"); D=$3; CASES=${4:-"level_set_fuel_model_1 LS_ellipse_0ms_00deg LS_ellipse_0ms_30deg LS_ellipse_5ms_00deg LS_ellipse_5ms_30deg LS_wind_ramp_lin LS_wind_ramp_quad"}
HERE=$(cd "$(dirname "$0")" && pwd); V=$(cd "$HERE/../../../Verification/WUI" && pwd); TOL=${TOL:-1e-6}
source "$HERE/env.sh" >/dev/null 2>&1; export OMP_NUM_THREADS=1
mkdir -p "$D"; fail=0
for c in $CASES; do
  for k in stock amr; do w="$D/$c/$k"; rm -rf "$w"; mkdir -p "$w"; cp "$V/$c".* "$V"/*.csv "$w"/ 2>/dev/null; done
  case $c in LS*) for k in amr; do sed -i 's#^&TAIL#\&AMR BLOCKING_FACTOR=1 /\n\&TAIL#' "$D/$c/$k/$c.fds"; done ;; esac   # 50x50x5 cells: the default blocking factor 8 does not divide the level-0 domain
  (cd "$D/$c/stock" && timeout ${TMO:-3000} "$S" $c.fds > run.log 2>&1; echo $? > rc)
  (cd "$D/$c/amr" && timeout ${TMO:-3000} "$A" $c.fds --run > run.log 2>&1; echo $? > rc)
  if [ -n "$REFINED" ]; then
    w="$D/$c/ref"; rm -rf "$w"; mkdir -p "$w"; cp "$V/$c".* "$V"/*.csv "$w"/ 2>/dev/null
    case $c in LS*) sed -i 's#^&TAIL#\&AMR BLOCKING_FACTOR=1 /\n\&TAIL#' "$w/$c.fds" ;; esac
    (cd "$w" && timeout ${TMO:-3000} "$A" $c.fds --two-level-run --init-patch --patch ${PATCH:-4 12 4 12} --regrid-every 2 --tag-phi --blocking 4 --maxsize 8 > run.log 2>&1; echo $? > rc)
  fi
  python3 "$HERE/wui_levelset_compare.py" "$D/$c" "$c" || fail=1
done
exit $fail
