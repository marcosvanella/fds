#!/bin/bash
# S9.5 (A)(ii): FDS pressure codes and driver boundary strings for every verification input with a one-cell direction (`fds_amr case.fds --pressure-bc`, set-up only, no time loop).
# usage: run_pressure_bc_sweep.sh <fds_amr> <case-list> <work-dir>      (<case-list>: one path relative to Verification/ per line; Verification root = $VERIF or <repo>/Verification)
# Output: <work-dir>/pbc_lines.txt (one `PBC ...` line per case, or `NOLINE <case>: <first error>`; resumable); the table made from it is notes/fft-thin-direction-cases.csv.
# resumable sweep: appends to pbc/pbc_lines.txt; cases already present are skipped; relative ../../Utilities paths resolve through a mirror root
source "$(dirname "$0")/env.sh" >/dev/null 2>&1
HERE=$(cd "$(dirname "$0")" && pwd); B=$1; LIST=$2; W=$3; V=${VERIF:-$(cd "$HERE/../../../Verification" && pwd)}; UT=$(cd "$V/../Utilities" && pwd); OUT=$W/pbc_lines.txt; mkdir -p $W; touch $OUT
sed -i '/^finished$/d' $OUT
while read f; do
  grep -q " $f[ :]" $OUT && continue
  g=$(dirname $f); r=$W/root; d=$r/Verification/$g; rm -rf $r; mkdir -p $d; ln -s "$UT" $r/Utilities
  for x in $V/$g/*; do [ -f "$x" ] && ln -s "$x" $d/; done
  rm -f $d/$(basename $f); cp $V/$f $d/
  (cd $d && timeout 90 mpirun --bind-to none --oversubscribe -np 1 $B $(basename $f) --pressure-bc > out.txt 2> err.txt < /dev/null)
  l=$(grep "^PBC " $d/out.txt | head -1)
  if [ -z "$l" ]; then echo "NOLINE $f: $(grep -h -i -m1 'error\|abort\|MPI_PROCESS' $d/out.txt $d/err.txt | head -1 | cut -c1-150)" >> $OUT; else echo "$l" | sed "s#PBC [^ ]* #PBC $f #" >> $OUT; fi
done < "$LIST"
echo finished >> $OUT
