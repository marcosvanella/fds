#!/bin/bash
# FireX exchange-code check (r6 Role 1 item 1): code 6 = U/V/W/H only; codes 66/67 never routed through the driver ghost fill.
# usage: run_exchange_code_check.sh [Source/driver dir]
set -u
HERE=$(cd "$(dirname "$0")" && pwd); D=${1:-$(cd "$HERE/.." && pwd)}; rc=0
# 1. Mapping: ExchangeFields.H case 6 returns U,V,W,H; 66/67 fall through to default empty.
grep -q 'case 6: return {"U", "V", "W", "H"}' "$D/ExchangeFields.H" || { echo "FAIL ExchangeFields.H case 6"; rc=1; }
! grep -E 'case 66:|case 67:' "$D/ExchangeFields.H" || { echo "FAIL ExchangeFields.H must not list 66/67 as driver codes"; rc=1; }
# 2. TimeLoop never calls exchange 66 or 67 (they stay in FDS MAIN_LOOP / INITIALIZE_BACK_WALL_EXCHANGE).
if grep -En 's_exchange\(66|s_exchange\(67|s_after\(66|s_after\(67|after_exchange\(66|after_exchange\(67|[^a-z_]exchange\(66|[^a-z_]exchange\(67' "$D/TimeLoop.cpp"; then
  echo "FAIL TimeLoop.cpp must not call exchange/after for 66 or 67"; rc=1
else
  echo "PASS TimeLoop.cpp has no 66/67 exchange call"
fi
grep -Eq 's_exchange\(6|exchange\(6' "$D/TimeLoop.cpp" || { echo "FAIL TimeLoop.cpp must still call code 6"; rc=1; }
# 3. after_exchange only handles 1, 3, 4, 6 (so 66/67 return immediately).
grep -q 'code != 1 && code != 3 && code != 4 && code != 6' "$D/GhostExchange.cpp" || { echo "FAIL GhostExchange.cpp after_exchange gate"; rc=1; }
# 4. FDS still owns 66/67 call sites (MAIN_LOOP / INITIALIZE_BACK_WALL_EXCHANGE).
# FireX c2581655bb reverted the lightweight WALL exchange: MAIN_LOOP calls 66 again and 67 has no call site (code 67 remains a valid CODE value).
grep -q 'MESH_EXCHANGE(66)' "$D/../main.f90" && echo "PASS FDS main.f90 still calls MESH_EXCHANGE(66)" || { echo "FAIL FDS lost 66 call sites"; rc=1; }
[ $rc = 0 ] && echo "PASS exchange-code mapping (6 = U/V/W/H; 66/67 host-only)"
exit $rc
