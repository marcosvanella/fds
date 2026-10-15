#!/bin/bash
# run_l0_bitcheck.sh <new-bin> <old-bin> <abs-path-case.fds>: --run 20 on both, compare every csv/bin byte for byte (timing files excluded); run in a scratch dir
set -u
for t in new old; do rm -rf l0_$t; mkdir l0_$t; done
(cd l0_new && OMP_NUM_THREADS=1 $1 $3 --run 20 > run.log 2>&1)
(cd l0_old && OMP_NUM_THREADS=1 $2 $3 --run 20 > run.log 2>&1)
n=0; bad=0
for f in l0_new/*.csv l0_new/*.bin; do b=$(basename $f); case $b in *perf.csv|*_cpu.csv|*_steps.csv) continue;; esac; n=$((n+1)); cmp -s $f l0_old/$b || { echo "DIFF $b"; bad=$((bad+1)); }; done
echo "L0 --run 20: $n csv compared, $bad differ"
