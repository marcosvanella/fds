# Shared helper for the FDS-AMReX run_*.sh scripts (sourced, not run).
# Environment:
#   FDS_AMR   path to the fds_amr executable (required)
#   NP        MPI ranks (default 4; the inputs here have 4 meshes, MPI_PROCESS 0..3)
#   STEPS     time steps (default 10; a short smoke run)
#   MPIRUN    MPI launcher (default "mpirun --bind-to none --oversubscribe" for Open MPI;
#             with Intel MPI use MPIRUN=mpiexec, without --oversubscribe)
#   RUN_DIR   where outputs go (default ./out_<chid> next to this script)
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
CASES=$REPO/Source/driver/tests/cases
: "${FDS_AMR:?set FDS_AMR to the fds_amr executable}"
NP=${NP:-4}
STEPS=${STEPS:-10}
MPIRUN=${MPIRUN:-mpirun --bind-to none --oversubscribe}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

# run_case INPUT.fds CHID [driver flags ...]
run_case() {
  local input=$1 chid=$2; shift 2
  local dir=${RUN_DIR:-$HERE/out_$chid}
  mkdir -p "$dir"
  cp "$HERE/$input" "$dir/"
  # initial-velocity CSVs: the 48^3 cases use the committed slices, the 96^3 case is generated on demand
  case $input in
    vortex_ring_48_4m.fds|pyramid_48_4m.fds) cp "$CASES"/vortex_ring_48_4m_uvw_m?.csv "$dir/" ;;
    vortex_ring_96_2x2x1.fds)
      [ -f "$dir/vortex_ring_96_2x2x1_uvw_m4.csv" ] || \
        python3 "$REPO/Source/driver/tests/make_vortex_ring_case.py" 96 "$dir" --sx 2 --sy 2 --sz 1 >/dev/null
      cp "$HERE/$input" "$dir/"; rm -f "$dir/vortex_ring_96_uvw.csv" ;;
  esac
  cd "$dir"
  echo "+ $MPIRUN -np $NP $FDS_AMR $input --two-level-run --steps $STEPS --outdir . --chid $chid $*"
  $MPIRUN -np "$NP" "$FDS_AMR" "$input" --two-level-run --steps "$STEPS" --outdir . --chid "$chid" "$@" 2>&1 | tee "$chid.log"
}
