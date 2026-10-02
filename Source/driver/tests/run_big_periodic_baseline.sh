#!/bin/bash
# WP2 (stage-1 GPU spike plan): periodic N^3 CPU baseline for the driver and for the USE_AMREX=OFF FDS on the same input (N >= 128).
# usage: run_big_periodic_baseline.sh <driver-build> <off-build-or-refbin-dir containing fds> <work> [N=128] [steps=10]
# The input is generated (tests/make_periodic_case.py, about 56 MB of UVW csv at N=128, not committed): the same T_END (= `steps` constant steps of the case,
# dt = 9.38020e-3 at N=128, so T_END is set to steps*9.38e-3 - epsilon) makes both codes run the same number of steps, FDS stops at T_END.
# Output: <work>/baseline.csv with one row per variant: wall_s, peak RSS of the largest rank (kB), loop time, pressure time and f_pres.
#   driver : loop = t_loop_s of <chid>_driver_perf.csv, f_pres = t_pressure/t_loop (the driver's own profile, FDSTL_PROFILE=1)
#   off    : loop = MAIN-loop total of <chid>_cpu.csv (sum of the subroutine columns without set-up is not separated by FDS: Total T_USED is reported, f_pres = PRES / Total)
# Valid only on an idle machine (load < 1, >= 6 GB free at N=128; the 128^3 driver peaks at 2.5 GB RSS). A memory watchdog kills a run when MemAvailable < 350 MB.
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
B=$1; OFFB=$2; W=$3; N=${4:-128}; STEPS=${5:-10}
mkdir -p "$W"; OUT=$W/baseline.csv
DT=$(python3 -c "print('%.6f' % (${STEPS}*9.38e-3 - 1e-6))")   # 9.38e-3 is the first-step dt of the case at N=128 (IJK-independent to 4 digits for the 0.56549 m box only at N=128)
[ "$N" != 128 ] && DT=$(python3 -c "print('%.6f' % (${STEPS}*9.38e-3*128/$N - 1e-6))")   # dt scales with dx (advective limit) as a first guess; steps are then approximate
[ -f "$W/tg$N/tg${N}_uvw.csv" ] || python3 "$HERE/make_periodic_case.py" $N "$W/tg$N" $DT
echo "variant,N,steps,load1,free_gb,wall_s,max_rss_kb_one_rank,exit,loop_s,pressure_s,f_pres,t_end" > "$OUT"
run() {   # variant, binary, args...
  v=$1; bin=$2; shift 2
  D=$W/run_$v; rm -rf "$D"; mkdir -p "$D"; cp "$W/tg$N/tg$N.fds" "$W/tg$N/tg${N}_uvw.csv" "$D/"
  sed -i "s/T_END=[0-9.]*/T_END=$DT/; s/T=0.67/T=$DT/" "$D/tg$N.fds"
  load=$(cut -d' ' -f1 /proc/loadavg); free=$(awk '/MemAvailable/{printf "%.1f",$2/1048576}' /proc/meminfo)
  ( while sleep 2; do a=$(awk '/MemAvailable/{print $2}' /proc/meminfo); if [ "$a" -lt 358400 ]; then echo "watchdog: MemAvailable $a kB, killing" > "$D/watchdog.txt"; pkill -f "$bin" ; break; fi; done ) &
  WD=$!
  res=$(FDSTL_PROFILE=1 python3 "$HERE/perf_one.py" "$D" 1 "$bin" "$@")
  kill $WD 2>/dev/null; wait $WD 2>/dev/null
  echo "$v,$N,$STEPS,$load,$free,$res" > "$D/row.txt"
}
run driver "$B/fds_amr" tg$N.fds --run --outdir . --chid tg$N
p=$(sed -n 2p "$W/run_driver/tg${N}_driver_perf.csv" 2>/dev/null)
echo "$(cat $W/run_driver/row.txt),$(echo $p | cut -d, -f7),$(echo $p | cut -d, -f8),$(echo $p | cut -d, -f11),$DT" >> "$OUT"
if [ -x "$OFFB/fds" ]; then
  run off "$OFFB/fds" tg$N.fds
  c=$(sed -n 2p "$W/run_off/tg${N}_cpu.csv" 2>/dev/null)
  tot=$(echo $c | awk -F, '{print $NF}'); pr=$(echo $c | awk -F, '{print $6}')
  fp=$(python3 -c "print('%.4f' % (float($pr)/float($tot)))" 2>/dev/null)
  echo "$(cat $W/run_off/row.txt),$tot,$pr,$fp,$DT" >> "$OUT"
fi
cat "$OUT"
