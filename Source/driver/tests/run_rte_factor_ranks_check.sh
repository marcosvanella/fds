#!/bin/bash
# run_rte_factor_ranks_check.sh <fds_amr binary> [<run-dir>]: cases/rte_factor_ranks.fds (wide HRRPUA burner, radiation on), level 1 pinned by
#   --patch 2 29 0 7 --maxsize 8 (128 boxes, static), at 1, 2 and 4 MPI ranks. The RTE source correction factor (fds_kernels.f90 FDS_K_RTE_FACTOR) sums
#   the flame cells of many boxes on every rank; it enters Q_RADI and HRR. Pass: gate PASS at every rank count and the DEVC file (point temperatures in
#   the fine boxes, no reduction) bitwise identical to np=1; the HRR file is compared for information only (its columns are rank sums). Also checks that Q_RADI is nonzero (radiation and flame cells present).
set -u
BIN=$1; RUN=${2:-$PWD/rte_factor_ranks}; HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$RUN" || exit 1
fail=0
for np in 1 2 4; do
  d=$RUN/np$np; rm -rf "$d"; mkdir -p "$d"
  python3 - "$HERE/cases/rte_factor_ranks.fds" "$d/c.fds" $np <<'P'
import re,sys
n=[0]; k=int(sys.argv[3])
def f(m):
    r=n[0]*k//4; n[0]+=1; return re.sub(r'MPI_PROCESS=\d','MPI_PROCESS=%d'%r,m.group(0))
open(sys.argv[2],'w').write(re.sub(r'&MESH[^/]*/',f,open(sys.argv[1]).read()))
P
  ( cd "$d" && OMP_NUM_THREADS=1 mpirun --oversubscribe -np $np "$BIN" c.fds --two-level-run --patch 2 29 0 7 --maxsize 8 --nofine -1 > run.log 2>&1 ); rc=$?
  g=$(tr -d '\000' < "$d/run.log" | grep -ac "RESULT gate PASS")
  echo "np=$np rc=$rc gate=$g"
  [ $rc -eq 0 ] && [ "$g" -eq 1 ] || fail=1
  rm -f "$d"/*.s3d "$d"/*.s3d.sz "$d"/*.smv
  if [ $np -gt 1 ]; then
    if cmp -s "$RUN/np1/rte_factor_ranks_devc.csv" "$d/rte_factor_ranks_devc.csv"; then echo "  devc bitwise = np1"; else echo "  devc DIFFERS from np1"; fail=1; fi
    cmp -s "$RUN/np1/rte_factor_ranks_hrr.csv" "$d/rte_factor_ranks_hrr.csv" && echo "  hrr bitwise = np1 (info)" || echo "  hrr differs from np1 (info: the HRR columns are rank sums)"
  fi
done
python3 - "$RUN/np1/rte_factor_ranks_hrr.csv" <<'P' || fail=1
import csv,sys
r=list(csv.reader(open(sys.argv[1])))[2:]; q=max(abs(float(x[3])) for x in r)
print("max |Q_RADI| %.6g kW"%q); sys.exit(0 if q>1.0 else 1)
P
[ $fail -eq 0 ] && echo "PASS rte_factor_ranks" || echo "FAIL rte_factor_ranks (see $RUN)"
exit $fail
