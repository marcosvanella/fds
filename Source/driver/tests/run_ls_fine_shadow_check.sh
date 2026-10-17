#!/bin/bash
# run_ls_fine_shadow_check.sh <build-dir> [<run-dir>]: --run --ls-fine-shadow on the CSIRO C064 level-set crude case, one level, at 12x12x4 (10 m cells) and 24x24x8 (5 m cells).
# The level set of the global LS grid module (fds_ls_grid.f90, a copy of the vege.f90 level-set numerics; grid = the level-0 mesh) is advanced beside the original LEVEL_SET_FIRESPREAD
# in every stage and compared bit for bit (--ls-fine-shadow: wind sampled from the level-0 mesh in Fortran; --ls-fine-shadow-gather: wind gathered from the leaf cells by LsFine.cpp): PHI (all values incl. ghost cells), T_IGN, TAU_LS, BURN_DURATION, AREA_ADJUST of every ground cell. Fails on any difference or if nothing ignited.
set -u
BLD=$1; RUN=${2:-$BLD/run/ls_fine_shadow}
HERE=$(cd "$(dirname "$0")" && pwd)
rc=0
for FLAG in --ls-fine-shadow --ls-fine-shadow-gather; do
for IJK in 12,12,4 24,24,8; do
D=$RUN/${FLAG#--}_ijk_${IJK//,/_}
rm -rf "$D"; mkdir -p "$D" && cd "$D" || exit 1
cp "$HERE/cases/levelset_tag_phi_ground.fds" c.fds
sed -i 's/^&HEAD CHID=[^,]*,/\&HEAD CHID='"'c'"',/' c.fds
sed -i 's/IJK=12,12,4/IJK='"$IJK"'/' c.fds
sed -i '1a &AMR MAX_LEVEL=0, BLOCKING_FACTOR=2, MAX_GRID_SIZE=24 /' c.fds
OMP_NUM_THREADS=2 "$BLD/fds_amr" c.fds --run $FLAG --outdir . --chid c > run.log 2>&1
if grep -aq "RUN COMPLETE\|completed successfully" run.log && [ -f c_ls_shadow.txt ] && grep -q "PHI values differing 0, cells with differing T_IGN/TAU/BURN_DURATION/AREA_ADJUST 0," c_ls_shadow.txt && ! grep -q "ignited cells 0$" c_ls_shadow.txt; then
  echo "PASS ${FLAG#--} IJK=$IJK: $(sed 's/^LS shadow: stages compared with the original, //' c_ls_shadow.txt)"
else echo "FAIL ${FLAG#--} IJK=$IJK (see $D/run.log, $D/c_ls_shadow.txt): $(cat c_ls_shadow.txt 2>/dev/null)"; rc=1; fi
done
done
exit $rc
