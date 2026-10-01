#!/bin/bash
# NFR-030 / A-22 / NFR-031 measurement for the M2a cases (measured, not required): wall time, pressure share f_pres, peak memory.
# usage: run_perf.sh <driver-build> <unpatched-FDS-build or -> <work> [reps]
# For each case/ranks it runs the driver (fds_amr --run) and, when given, the USE_AMREX=OFF fds of the same source tree, REPS times each (default 3), on the same machine.
# Valid only on an idle machine (load < 1, >= 8 GB free, NFR-030 / test-plan 9.2): the script records the load before every run and flags a busy machine; a run on
# a busy machine is invalid, not failed. Output: <work>/perf.csv.
set -u
HERE=$(cd "$(dirname "$0")" && pwd); source "$HERE/env.sh" > /dev/null 2>&1
B=$1; OFFB=$2; W=$3; REPS=${4:-3}
CASES="shunn3_32:1 csmag_32:1 shunn3_4mesh_32:1 shunn3_4mesh_32:4"
CSDIR=${CSMAG_DIR:-/workspace/fds-amr/scratch/role1-s3-work/ref_runs/csmag_32}
mkdir -p "$W"; OUT=$W/perf.csv
echo "case,np,variant,rep,load1,free_gb,wall_s,max_rss_kb_one_rank,exit,driver_loop_s,driver_pressure_s,driver_solve_s,f_pres,f_solve,field_bytes_sum" > "$OUT"
for cn in $CASES; do
  c=${cn%%:*}; n=${cn##*:}
  if [ "$c" = csmag_32 ]; then SRC=$CSDIR; else SRC=$BASELINE/$c; fi
  for variant in driver off; do
    [ "$variant" = off ] && [ "$OFFB" = "-" ] && continue
    for r in $(seq 1 $REPS); do
      D=$W/$c.np$n.$variant.$r; rm -rf "$D"; mkdir -p "$D"; cp "$SRC/$c.fds" "$D/"
      for f in "$SRC"/*uvw*.csv; do [ -f "$f" ] && cp "$f" "$D/"; done
      load=$(cut -d' ' -f1 /proc/loadavg); free=$(awk '/MemAvailable/{printf "%.1f",$2/1048576}' /proc/meminfo)
      if [ $variant = driver ]; then
        res=$(python3 "$HERE/perf_one.py" "$D" $n "$B/fds_amr" $c.fds --run --outdir . --chid $c --quiet)
        p=$(sed -n 2p "$D/${c}_driver_perf.csv" 2>/dev/null)
        loop=$(echo $p | cut -d, -f7); pres=$(echo $p | cut -d, -f8); sol=$(echo $p | cut -d, -f9); fp=$(echo $p | cut -d, -f11); fs=$(echo $p | cut -d, -f12); fb=$(echo $p | cut -d, -f16)
      else
        res=$(python3 "$HERE/perf_one.py" "$D" $n "$OFFB/fds" $c.fds)
        loop=; pres=; sol=; fp=; fs=; fb=
      fi
      echo "$c,$n,$variant,$r,$load,$free,$res,$loop,$pres,$sol,$fp,$fs,$fb" >> "$OUT"
    done
  done
done
cat "$OUT"
