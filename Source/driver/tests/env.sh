# Source this file to get the GNU 14.2 + Open MPI 5.0.7 + AMReX environment used for the driver (see ../README.md).
# Not a test itself. Override BASELINE/ROOT for another machine.
source /workspace/gnu_ompi/env_gnu_ompi.sh
export OMP_NUM_THREADS=1
export HWLOC_LIBXML=0
export HYPRE_FIREX=/workspace/firemodels-gnu/libs/hypre/63331f19
export BASELINE=${BASELINE:-/workspace/fds-amr/vv-runs/baseline/gnu_ompi_firex-36975d7}
# Common configure options of the reference (FireX) build, plus fixed version strings so that binaries are comparable.
export FDS_CMAKE_COMMON="-DCMAKE_BUILD_TYPE=Release -DUSE_SYSTEM_HYPRE=ON -DCMAKE_PREFIX_PATH=$HYPRE_FIREX \
-DUSE_SYSTEM_SUNDIALS=ON -DSUNDIALS_DIR=$SUNDIALS_DIR -DBUILD_DATE=fixed -DBUILD_DATE_XLF=fixed -DGIT_DATE=fixed \
-DGIT_BRANCH=FDS-AMReX -DGIT_HASH=FDS-6.11.1-1244-g36975d765f -DGIT_DIRTY="
