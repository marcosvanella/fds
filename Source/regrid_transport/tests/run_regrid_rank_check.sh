#!/bin/bash
# FR-015: the grids and the data of the dynamic regrid test must be identical at 1, 2 and 4 ranks. Usage: run_regrid_rank_check.sh <test_regrid_core> <mpiexec> <numproc flag>
exe="$1"; mpiexec="${2:-mpirun}"; flag="${3:--np}"
h1=$("$mpiexec" "$flag" 1 "$exe" 2>&1 | grep '^HASH')
h2=$("$mpiexec" "$flag" 2 "$exe" 2>&1 | grep '^HASH')
h4=$("$mpiexec" "$flag" 4 "$exe" 2>&1 | grep '^HASH')
echo "1 rank : $h1"; echo "2 ranks: $h2"; echo "4 ranks: $h4"
[ -n "$h1" ] && [ "$h1" = "$h2" ] && [ "$h1" = "$h4" ]
