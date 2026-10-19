#!/bin/bash
# run_ls_fine_drive_check.sh <build-dir> [<run-dir>]: --ls-fine (the LS grid drives the fire: level-0 walls carry the LS cell means, fine hidden walls the LS flux under them) on the CSIRO C064 level-set crude case.
# (1) one level, 24x24x8, 60 s: the HRR over 20-60 s stays within 2% of the original level set (the LS grid is bit-identical to it, the wall flux differs by rounding only, so only turbulent scatter remains).
# (2) two levels (--init-patch 5 9, --regrid-every 2 --tag-phi --nofine 0, 60 s): TWO-LEVEL RUN PASS (fine div, mass and species balance gates), fire ignited.
set -u
BLD=$1; RUN=${2:-$BLD/run/ls_fine_drive}
HERE=$(cd "$(dirname "$0")" && pwd)
rc=0
mkdir -p "$RUN"
for v in orig drive; do
  D=$RUN/one_$v; rm -rf "$D"; mkdir -p "$D"; cd "$D" || exit 1
  cp "$HERE/cases/levelset_tag_phi_ground.fds" c.fds
  sed -i "s/^&HEAD CHID=[^,]*,/\&HEAD CHID='c',/;s/IJK=12,12,4/IJK=24,24,8/" c.fds
  sed -i '1a &AMR MAX_LEVEL=0, BLOCKING_FACTOR=2, MAX_GRID_SIZE=24 /' c.fds
  FL=""; [ $v = drive ] && FL=--ls-fine
  OMP_NUM_THREADS=2 "$BLD/fds_amr" c.fds --run $FL --outdir . --chid c > run.log 2>&1
done
python3 - "$RUN" <<'PY' || rc=1
import sys, numpy as np
m = []
for v in ('orig', 'drive'):
    a = np.genfromtxt(f'{sys.argv[1]}/one_{v}/c_hrr.csv', delimiter=',', skip_header=2)
    m.append(a[a[:, 0] >= 20, 1].mean())
r = abs(m[1] - m[0]) / m[0]
print(f"{'PASS' if r < 0.02 and m[0] > 0 else 'FAIL'} one level: mean HRR 20-60 s orig {m[0]:.4g} kW, drive {m[1]:.4g} kW, diff {100*r:.2f}%")
sys.exit(0 if r < 0.02 and m[0] > 0 else 1)
PY
D=$RUN/two; rm -rf "$D"; mkdir -p "$D"; cd "$D" || exit 1
cp "$HERE/cases/levelset_tag_phi_ground.fds" c.fds
sed -i "s/^&HEAD CHID=[^,]*,/\&HEAD CHID='c',/" c.fds
sed -i '1a &AMR MAX_LEVEL=1, BLOCKING_FACTOR=2, MAX_GRID_SIZE=24 /' c.fds
OMP_NUM_THREADS=2 "$BLD/fds_amr" c.fds --two-level-run --levels 2 --init-patch --patch 5 9 0 11 --blocking 2 --maxsize 8 --regrid-every 2 --tag-phi --nofine 0 --ls-fine --outdir . --chid c > run.log 2>&1
if grep -aq "TWO-LEVEL RUN PASS" run.log; then echo "PASS two levels: $(grep -a 'RESULT gate PASS' run.log | cut -c1-160)"; else echo "FAIL two levels (see $D/run.log)"; rc=1; fi
exit $rc
