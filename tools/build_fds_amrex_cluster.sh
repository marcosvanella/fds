#!/usr/bin/env bash
# build_fds_amrex_cluster.sh - build FDS-AMReX (fds_amr) with Intel oneAPI (icx/icpx/ifx) + Intel MPI on a cluster.
#
# USAGE
#   module load <intel compilers> <intel mpi> <mkl> cmake git   # or:  source /opt/intel/oneapi/setvars.sh
#   export FIREMODELS=/path/to/firemodels                        # dir that holds the third-party clones and libs/
#   tools/build_fds_amrex_cluster.sh [--smoke] [--clean-fds] [--help]
#
#   --smoke      after the build run the 10-step simple_test (1 rank) as a check of the binary
#   --clean-fds  remove the FDS build directory first (the dependencies are kept)
#   -h|--help    print this header
#
# WHAT IT DOES (every step is skipped when its result is already there, so a second run only re-checks and runs make)
#   1. checks the toolchain: mpiicx, mpiicpx, mpiifx, icx, ifx, mpiexec, cmake >= 3.24, git, make, pkg-config; MKL (MKLROOT, the mkl module) is checked but optional
#   2. HYPRE   63331f19c72b2ec363824907c0575a6edf7bbf3c (v2.32.0-24)  clone $FIREMODELS/hypre   -> $LIBS/hypre/63331f19
#   3. SUNDIALS v7.5.0 (MPI, OpenMP, Fortran modules)                 clone $FIREMODELS/sundials -> $LIBS/sundials/v7.5.0
#   4. FFTW    3.3.10 (double+float, OpenMP) - only if pkg-config finds no fftw3/fftw3f  -> $LIBS/fftw/3.3.10
#   5. AMReX   99ddfda6722b63f7e128ee6d6822ff29046bb7a9 (26.09, development)  clone $FIREMODELS/amrex -> $LIBS/amrex/99ddfda6
#      built against the SAME HYPRE as FDS (the FDS CMake aborts on a mismatch)
#   6. FDS-AMReX: cmake -DUSE_AMREX=ON ... ; make fds  ->  $FDS_BUILD_DIR/fds_amr
#   Clones are made at the $FIREMODELS level (next to fds/) only when missing. An existing clone is never modified, except a
#   detached checkout of the pinned commit when its working tree is clean. No network is used when the clones and the FFTW tarball exist.
#   The FDS-AMReX source is the repository this script lives in (tools/..).
#
# RECIPE: AMReX version/options/driver command/dependencies come from the AMReX Integration Lead's notes and README_FDS_AMReX.md
#   (sections 1-2). Intel flags: D-079, Fortran/C/C++ "-O2 -fp-model=precise" (NOT -O3: ifx -O3 -fiopenmp miscompiles Source/pois.f90).
#
# ENVIRONMENT (all optional except FIREMODELS)
#   FIREMODELS     base directory (default: parent of the repository, the firemodels convention). Created if missing.
#   LIBS           install root of the libraries (default $FIREMODELS/libs/impi_intel)
#   FDS_BUILD_DIR  FDS build directory (default <repo>/Build/cmakeb/impi_intel_amrex_rel)
#   NJOBS          parallel make jobs (default min(nproc,8))
#   BUILD_NICE     niceness for configure/build steps (default 10; 0 disables)
#   FDS_USE_MKL    OFF (default, as the recipe says: static MKL exports FFTW symbols that break the AMReX FFT pressure solve) or ON
#   REQUIRE_MKL    1 = fail when MKLROOT is not set (always required when FDS_USE_MKL=ON); default: only a warning
#   ONEAPI_SETVARS path of oneAPI setvars.sh, sourced when mpiifx is not on PATH (default /opt/intel/oneapi/setvars.sh)
#   INTEL_MODULES  space separated module names to "module load" when the tools are not found (e.g. "intel-oneapi-compilers intel-oneapi-mpi intel-oneapi-mkl")
#   FFTW_PREFIX    use this FFTW install (needs lib/pkgconfig/fftw3.pc and fftw3f.pc) instead of searching/building
#   GIT_HYPRE_URL, GIT_AMREX_URL, GIT_SUNDIALS_URL, FFTW_URL   source locations (defaults: GitHub, fftw.org)
#   FDS_AMR_REPO   repository to build (default: the one containing this script)
#
# Run it on a login node only if the site allows compiling there; otherwise run it inside an interactive/batch job.
# Run the result with OMP_NUM_THREADS=1..4 and I_MPI_PIN=0 or Intel MPI pinning; do not pass --oversubscribe to Intel mpiexec.

set -euo pipefail

SCRIPT_VERSION=1
HYPRE_SHA=63331f19c72b2ec363824907c0575a6edf7bbf3c
AMREX_SHA=99ddfda6722b63f7e128ee6d6822ff29046bb7a9
SUNDIALS_TAG=v7.5.0
FFTW_VER=3.3.10
# D-079 flags for Intel
INTEL_FFLAGS="-O2 -fp-model=precise"
INTEL_CFLAGS="-O2 -fp-model=precise -DNDEBUG"

SMOKE=0; CLEAN_FDS=0
for a in "$@"; do
  case $a in
    --smoke) SMOKE=1 ;;
    --clean-fds) CLEAN_FDS=1 ;;
    -h|--help) sed -n '2,/^set -euo/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown argument: $a (see --help)" >&2; exit 2 ;;
  esac
done

log()  { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
die()  { printf '[%s] ERROR: %s\n' "$(date +%H:%M:%S)" "$*" >&2; exit 1; }
T0=$(date +%s)

# ---------------------------------------------------------------- locations
SELF=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=${FDS_AMR_REPO:-$(cd "$SELF/.." && pwd)}
[ -f "$REPO/CMakeLists.txt" ] && [ -d "$REPO/Source/driver" ] || die "$REPO is not an FDS-AMReX checkout (no Source/driver)"
FIREMODELS=${FIREMODELS:-$(cd "$REPO/.." && pwd)}
mkdir -p "$FIREMODELS"; FIREMODELS=$(cd "$FIREMODELS" && pwd)
LIBS=${LIBS:-$FIREMODELS/libs/impi_intel}
FDS_BUILD_DIR=${FDS_BUILD_DIR:-$REPO/Build/cmakeb/impi_intel_amrex_rel}
FDS_USE_MKL=${FDS_USE_MKL:-OFF}
BUILD_NICE=${BUILD_NICE:-10}
GIT_HYPRE_URL=${GIT_HYPRE_URL:-https://github.com/hypre-space/hypre.git}
GIT_AMREX_URL=${GIT_AMREX_URL:-https://github.com/AMReX-Codes/amrex.git}
GIT_SUNDIALS_URL=${GIT_SUNDIALS_URL:-https://github.com/LLNL/sundials.git}
FFTW_URL=${FFTW_URL:-https://www.fftw.org/fftw-$FFTW_VER.tar.gz}
if [ -z "${NJOBS:-}" ]; then NJOBS=$(nproc 2>/dev/null || echo 4); [ "$NJOBS" -gt 8 ] && NJOBS=8; fi
HYPRE_PREFIX=$LIBS/hypre/${HYPRE_SHA:0:8}
SUNDIALS_PREFIX=$LIBS/sundials/$SUNDIALS_TAG
AMREX_PREFIX=$LIBS/amrex/${AMREX_SHA:0:8}
FFTW_PREFIX_BUILT=$LIBS/fftw/$FFTW_VER
SRC_CACHE=$FIREMODELS/libs/src
WORK=$LIBS/build          # out-of-source build trees of the dependencies
mkdir -p "$LIBS" "$SRC_CACHE" "$WORK"

nicer() { if [ "$BUILD_NICE" != 0 ] && command -v nice >/dev/null; then nice -n "$BUILD_NICE" "$@"; else "$@"; fi; }

# ---------------------------------------------------------------- toolchain / required modules
have_tools() { command -v mpiicx >/dev/null && command -v mpiicpx >/dev/null && command -v mpiifx >/dev/null && command -v icx >/dev/null && command -v ifx >/dev/null && command -v mpiexec >/dev/null; }
if ! have_tools; then
  if [ -n "${INTEL_MODULES:-}" ] && { type module >/dev/null 2>&1 || [ -f /etc/profile.d/modules.sh ] || [ -f /usr/share/lmod/lmod/init/bash ]; }; then
    type module >/dev/null 2>&1 || . /etc/profile.d/modules.sh 2>/dev/null || . /usr/share/lmod/lmod/init/bash 2>/dev/null || true
    for m in $INTEL_MODULES; do log "module load $m"; module load "$m" || die "module load $m failed"; done
  fi
fi
if ! have_tools; then
  SV=${ONEAPI_SETVARS:-/opt/intel/oneapi/setvars.sh}
  if [ -f "$SV" ]; then log "sourcing $SV"; set +u; source "$SV" --force >/dev/null 2>&1 || true; set -u; fi
fi
missing=""
for t in mpiicx mpiicpx mpiifx icx ifx mpiexec cmake git make; do command -v "$t" >/dev/null || missing="$missing $t"; done
if [ -n "$missing" ]; then
  die "required tools not found:$missing. Load the Intel compiler (icx/ifx), Intel MPI (mpiicx/mpiicpx/mpiifx/mpiexec), cmake and git modules, e.g. 'module load intel-oneapi-compilers intel-oneapi-mpi intel-oneapi-mkl cmake', or source oneAPI setvars.sh (set INTEL_MODULES or ONEAPI_SETVARS to let this script do it)."
fi
# MKL: NOT needed for this build (FDS_USE_MKL=OFF keeps the real FFTW in charge of the AMReX FFT solver). It is only a hard
# requirement when FDS_USE_MKL=ON or REQUIRE_MKL=1; otherwise a missing MKLROOT is only a warning.
if [ -n "${MKLROOT:-}" ] && [ -d "$MKLROOT" ]; then log "MKL found ($MKLROOT); not linked unless FDS_USE_MKL=ON"
elif [ "${FDS_USE_MKL}" = ON ] || [ "${REQUIRE_MKL:-0}" = 1 ]; then
  die "MKL is required (FDS_USE_MKL=ON or REQUIRE_MKL=1) but MKLROOT is not set: load the mkl module (e.g. module load intel-oneapi-mkl) or source setvars.sh"
else log "note: MKLROOT not set (MKL module not loaded); fine, MKL is not used with FDS_USE_MKL=OFF"
fi
CMV=$(cmake --version | sed -n '1s/[^0-9]*\([0-9.]*\).*/\1/p')
[ "$(printf '%s\n3.24\n' "$CMV" | sort -V | head -1)" = 3.24 ] || die "cmake $CMV is too old (need >= 3.24)"
case "$FDS_USE_MKL" in ON|OFF) ;; *) die "FDS_USE_MKL must be ON or OFF" ;; esac
[ "$FDS_USE_MKL" = ON ] && log "WARNING: FDS_USE_MKL=ON: static MKL FFTW symbols can break the AMReX FFT pressure solve (README_FDS_AMReX.md 2.2)"

export CC=mpiicx CXX=mpiicpx FC=mpiifx
export I_MPI_PIN=${I_MPI_PIN:-0} HWLOC_LIBXML=0
FP_ICX=$(icx --version | head -1); FP_IFX=$(ifx --version | head -1); FP_MPI=$(mpiexec --version 2>&1 | head -1)
log "FIREMODELS=$FIREMODELS  LIBS=$LIBS  repo=$REPO"
log "toolchain: $FP_ICX | $FP_IFX | $FP_MPI | cmake $CMV | NJOBS=$NJOBS"
TOOLKEY="v$SCRIPT_VERSION|$FP_ICX|$FP_IFX|$FP_MPI|$INTEL_FFLAGS|$INTEL_CFLAGS"

# ---------------------------------------------------------------- helpers
# stamp_ok PREFIX KEY: true when PREFIX/.fds_amrex_stamp equals KEY
stamp_ok() { [ -f "$1/.fds_amrex_stamp" ] && [ "$(cat "$1/.fds_amrex_stamp")" = "$2" ]; }
stamp_set() { printf '%s' "$2" > "$1/.fds_amrex_stamp"; }

# ensure_clone DIR URL SHA_OR_TAG [full]: clone when missing; check out the pin when the tree is clean
ensure_clone() {
  local dir=$1 url=$2 pin=$3
  if [ ! -d "$dir/.git" ] && [ ! -f "$dir/.git" ]; then
    [ -e "$dir" ] && die "$dir exists and is not a git clone"
    log "cloning $url -> $dir"
    if [ "${4:-}" = tag ]; then git clone --depth 1 --branch "$pin" "$url" "$dir"
    else git clone --quiet --filter=blob:none "$url" "$dir" 2>/dev/null || { rm -rf "$dir"; git clone --quiet "$url" "$dir"; }
    fi
  fi
  local want have
  if [ "${4:-}" = tag ]; then
    want=$(git -C "$dir" rev-parse "$pin^{commit}" 2>/dev/null || true)
    if [ -z "$want" ]; then git -C "$dir" fetch --quiet --depth 1 origin "refs/tags/$pin:refs/tags/$pin"; want=$(git -C "$dir" rev-parse "$pin^{commit}"); fi
  else
    want=$pin
    git -C "$dir" cat-file -e "$pin^{commit}" 2>/dev/null || git -C "$dir" fetch --quiet origin "$pin" 2>/dev/null || git -C "$dir" fetch --quiet origin
    git -C "$dir" cat-file -e "$pin^{commit}" 2>/dev/null || die "commit $pin not found in $dir (remote $url)"
  fi
  have=$(git -C "$dir" rev-parse HEAD)
  if [ "$have" != "$want" ]; then
    [ -z "$(git -C "$dir" status --porcelain --untracked-files=no)" ] || die "$dir has local changes and is not at $pin; refusing to change it"
    log "checking out $pin in $dir (detached)"
    git -C "$dir" -c advice.detachedHead=false checkout --quiet "$want"
  fi
}

cmake_cfg_build_install() { # SRC BUILD ARGS...
  local src=$1 bld=$2; shift 2
  rm -rf "$bld"
  nicer cmake -S "$src" -B "$bld" "$@" > "$bld.configure.log" 2>&1 || { tail -30 "$bld.configure.log"; die "configure failed ($bld.configure.log)"; }
  nicer cmake --build "$bld" -j "$NJOBS" > "$bld.build.log" 2>&1 || { tail -40 "$bld.build.log"; die "build failed ($bld.build.log)"; }
  cmake --install "$bld" > "$bld.install.log" 2>&1 || { tail -30 "$bld.install.log"; die "install failed ($bld.install.log)"; }
}
CFLAGS_REL="$INTEL_CFLAGS"
COMMON=(-DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=mpiicx -DCMAKE_CXX_COMPILER=mpiicpx -DCMAKE_Fortran_COMPILER=mpiifx
        "-DCMAKE_C_FLAGS_RELEASE=$INTEL_CFLAGS" "-DCMAKE_CXX_FLAGS_RELEASE=$INTEL_CFLAGS" "-DCMAKE_Fortran_FLAGS_RELEASE=$INTEL_FFLAGS")

# ---------------------------------------------------------------- 2. HYPRE
HYPRE_KEY="hypre|$HYPRE_SHA|$TOOLKEY"
if stamp_ok "$HYPRE_PREFIX" "$HYPRE_KEY"; then log "HYPRE: up to date ($HYPRE_PREFIX)"
else
  log "HYPRE: building $HYPRE_SHA"
  ensure_clone "$FIREMODELS/hypre" "$GIT_HYPRE_URL" "$HYPRE_SHA"
  mkdir -p "$HYPRE_PREFIX"
  cmake_cfg_build_install "$FIREMODELS/hypre/src" "$WORK/hypre" "${COMMON[@]}" -DCMAKE_C_STANDARD=99 \
    -DCMAKE_INSTALL_PREFIX="$HYPRE_PREFIX" -DCMAKE_POSITION_INDEPENDENT_CODE=OFF -DBUILD_SHARED_LIBS=OFF \
    -DHYPRE_ENABLE_MPI=ON -DHYPRE_ENABLE_OPENMP=OFF -DHYPRE_ENABLE_BIGINT=OFF -DHYPRE_ENABLE_MIXEDINT=OFF \
    -DHYPRE_ENABLE_SINGLE=OFF -DHYPRE_ENABLE_LONG_DOUBLE=OFF -DHYPRE_ENABLE_COMPLEX=OFF \
    -DHYPRE_ENABLE_HYPRE_BLAS=ON -DHYPRE_ENABLE_HYPRE_LAPACK=ON -DHYPRE_ENABLE_PERSISTENT_COMM=OFF \
    -DHYPRE_ENABLE_FEI=OFF -DHYPRE_ENABLE_HOPSCOTCH=OFF -DHYPRE_ENABLE_LTO=OFF -DHYPRE_ENABLE_PRINT_ERRORS=OFF \
    -DHYPRE_ENABLE_TIMING=OFF -DHYPRE_ENABLE_STRICT_CHECKING=OFF -DHYPRE_ENABLE_HOST_MEMORY=ON \
    -DHYPRE_BUILD_EXAMPLES=OFF -DHYPRE_BUILD_TESTS=OFF \
    -DHYPRE_ENABLE_CUDA=OFF -DHYPRE_ENABLE_HIP=OFF -DHYPRE_ENABLE_SYCL=OFF -DHYPRE_ENABLE_GPU_AWARE_MPI=OFF \
    -DHYPRE_ENABLE_UNIFIED_MEMORY=OFF -DHYPRE_ENABLE_UMPIRE=OFF -DHYPRE_ENABLE_SUPERLU=OFF -DHYPRE_ENABLE_DSUPERLU=OFF \
    -DHYPRE_ENABLE_MAGMA=OFF -DHYPRE_ENABLE_CALIPER=OFF
  stamp_set "$HYPRE_PREFIX" "$HYPRE_KEY"
fi
HYPRE_LIBDIR=$HYPRE_PREFIX/lib; [ -f "$HYPRE_LIBDIR/libHYPRE.a" ] || HYPRE_LIBDIR=$HYPRE_PREFIX/lib64
[ -f "$HYPRE_LIBDIR/libHYPRE.a" ] || die "libHYPRE.a not found under $HYPRE_PREFIX"

# ---------------------------------------------------------------- 3. SUNDIALS
SUN_KEY="sundials|$SUNDIALS_TAG|$TOOLKEY"
if stamp_ok "$SUNDIALS_PREFIX" "$SUN_KEY"; then log "SUNDIALS: up to date ($SUNDIALS_PREFIX)"
else
  log "SUNDIALS: building $SUNDIALS_TAG"
  ensure_clone "$FIREMODELS/sundials" "$GIT_SUNDIALS_URL" "$SUNDIALS_TAG" tag
  mkdir -p "$SUNDIALS_PREFIX"
  # same options as the FetchContent block of the top-level CMakeLists.txt (MPI, OpenMP, Fortran modules, shared)
  cmake_cfg_build_install "$FIREMODELS/sundials" "$WORK/sundials" "${COMMON[@]}" -DCMAKE_INSTALL_PREFIX="$SUNDIALS_PREFIX" \
    -DENABLE_MPI=ON -DENABLE_OPENMP=ON -DBUILD_FORTRAN_MODULE_INTERFACE=ON -DENABLE_HYPRE=OFF \
    -DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_ENABLE_CXX=OFF -DEXAMPLES_ENABLE_F2003=OFF -DEXAMPLES_INSTALL=OFF \
    -DBUILD_SHARED_LIBS=ON -DBUILD_STATIC_LIBS=ON -DBUILD_TESTING=OFF
  stamp_set "$SUNDIALS_PREFIX" "$SUN_KEY"
fi
SUNDIALS_CMAKE=$SUNDIALS_PREFIX/lib/cmake/sundials; [ -d "$SUNDIALS_CMAKE" ] || SUNDIALS_CMAKE=$SUNDIALS_PREFIX/lib64/cmake/sundials
[ -d "$SUNDIALS_CMAKE" ] || die "SUNDIALS cmake config not found under $SUNDIALS_PREFIX"

# ---------------------------------------------------------------- 4. FFTW (AMReX FFT, the default level-0 pressure solver)
fftw_pc_dir=""
if [ -n "${FFTW_PREFIX:-}" ]; then
  fftw_pc_dir=$FFTW_PREFIX/lib/pkgconfig; [ -f "$fftw_pc_dir/fftw3.pc" ] || fftw_pc_dir=$FFTW_PREFIX/lib64/pkgconfig
  [ -f "$fftw_pc_dir/fftw3.pc" ] && [ -f "$fftw_pc_dir/fftw3f.pc" ] || die "FFTW_PREFIX=$FFTW_PREFIX has no fftw3.pc/fftw3f.pc"
elif command -v pkg-config >/dev/null && pkg-config --exists fftw3 fftw3f; then
  log "FFTW: system $(pkg-config --modversion fftw3)"
else
  FFTW_KEY="fftw|$FFTW_VER|$TOOLKEY"
  fftw_pc_dir=$FFTW_PREFIX_BUILT/lib/pkgconfig
  if stamp_ok "$FFTW_PREFIX_BUILT" "$FFTW_KEY"; then log "FFTW: up to date ($FFTW_PREFIX_BUILT)"
  else
    log "FFTW: building $FFTW_VER"
    tb=$SRC_CACHE/fftw-$FFTW_VER.tar.gz
    if [ ! -s "$tb" ]; then
      log "downloading $FFTW_URL"
      if command -v curl >/dev/null; then curl -fsSL -o "$tb.part" "$FFTW_URL"; else wget -q -O "$tb.part" "$FFTW_URL"; fi
      mv "$tb.part" "$tb"
    fi
    rm -rf "$WORK/fftw-$FFTW_VER"; mkdir -p "$WORK" "$FFTW_PREFIX_BUILT"
    tar -xzf "$tb" -C "$WORK"
    for prec in double float; do
      pdir=$WORK/fftw-$FFTW_VER-$prec; rm -rf "$pdir"; cp -r "$WORK/fftw-$FFTW_VER" "$pdir"
      extra=""; [ "$prec" = float ] && extra="--enable-float"
      ( cd "$pdir" && CC=icx CFLAGS="-O2 -fp-model=precise" nicer ./configure --prefix="$FFTW_PREFIX_BUILT" --enable-shared --enable-static \
          --enable-threads --enable-openmp $extra > "$pdir.configure.log" 2>&1 && nicer make -j "$NJOBS" > "$pdir.build.log" 2>&1 \
          && make install > "$pdir.install.log" 2>&1 ) || { tail -30 "$pdir".*.log; die "FFTW $prec build failed"; }
    done
    rm -rf "$WORK/fftw-$FFTW_VER"
    stamp_set "$FFTW_PREFIX_BUILT" "$FFTW_KEY"
  fi
fi
[ -n "$fftw_pc_dir" ] && export PKG_CONFIG_PATH="$fftw_pc_dir${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
[ -n "$fftw_pc_dir" ] && export LD_LIBRARY_PATH="${fftw_pc_dir%/pkgconfig}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
command -v pkg-config >/dev/null || die "pkg-config is required to let AMReX find FFTW"

# ---------------------------------------------------------------- 5. AMReX (same HYPRE as FDS)
AMREX_KEY="amrex|$AMREX_SHA|hypre=$HYPRE_SHA|fftw=${fftw_pc_dir:-system}|$TOOLKEY"
if stamp_ok "$AMREX_PREFIX" "$AMREX_KEY"; then log "AMReX: up to date ($AMREX_PREFIX)"
else
  log "AMReX: building $AMREX_SHA"
  ensure_clone "$FIREMODELS/amrex" "$GIT_AMREX_URL" "$AMREX_SHA"
  mkdir -p "$AMREX_PREFIX"
  cmake_cfg_build_install "$FIREMODELS/amrex" "$WORK/amrex" "${COMMON[@]}" -DCMAKE_INSTALL_PREFIX="$AMREX_PREFIX" \
    -DAMReX_SPACEDIM=3 -DAMReX_PRECISION=DOUBLE -DAMReX_MPI=ON -DAMReX_OMP=ON -DAMReX_GPU_BACKEND=NONE \
    -DAMReX_FORTRAN=ON -DAMReX_FORTRAN_INTERFACES=ON -DAMReX_EB=ON -DAMReX_PARTICLES=ON \
    -DAMReX_LINEAR_SOLVERS=ON -DAMReX_LINEAR_SOLVERS_EM=ON -DAMReX_LINEAR_SOLVERS_INCFLO=ON -DAMReX_AMRLEVEL=ON -DAMReX_PROBINIT=ON \
    -DAMReX_FFT=ON -DAMReX_FFTW_SEARCH=PKGCONFIG \
    -DAMReX_HYPRE=ON -DHYPRE_ROOT="$HYPRE_PREFIX" \
    -DAMReX_PIC=OFF -DAMReX_IPO=OFF -DAMReX_BUILD_SHARED_LIBS=OFF -DAMReX_INSTALL=ON \
    -DAMReX_ENABLE_TESTS=OFF -DBUILD_TESTING=OFF
  stamp_set "$AMREX_PREFIX" "$AMREX_KEY"
fi
AMREX_CMAKE=$AMREX_PREFIX/lib/cmake/AMReX; [ -d "$AMREX_CMAKE" ] || AMREX_CMAKE=$AMREX_PREFIX/lib64/cmake/AMReX
[ -d "$AMREX_CMAKE" ] || die "AMReX cmake config not found under $AMREX_PREFIX"
export LD_LIBRARY_PATH="$SUNDIALS_PREFIX/lib:$SUNDIALS_PREFIX/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# ---------------------------------------------------------------- 6. FDS-AMReX
[ "$CLEAN_FDS" = 1 ] && { log "removing $FDS_BUILD_DIR"; rm -rf "$FDS_BUILD_DIR"; }
mkdir -p "$(dirname "$FDS_BUILD_DIR")"
# Configure only when the build directory is new or its options changed: the FDS CMake embeds the configure time (BUILD_DATE) in
# every Fortran compile line, so a re-configure on each run would recompile everything. (The embedded build date is therefore the date of the configure.)
CFG_KEY="$TOOLKEY|$AMREX_PREFIX|$HYPRE_PREFIX|$SUNDIALS_PREFIX|$FDS_USE_MKL|$REPO"
if [ -f "$FDS_BUILD_DIR/CMakeCache.txt" ] && [ "$(cat "$FDS_BUILD_DIR/.fds_amrex_cfgkey" 2>/dev/null)" = "$CFG_KEY" ]; then
  log "FDS-AMReX: build directory already configured ($FDS_BUILD_DIR)"
else
  log "FDS-AMReX: configure $FDS_BUILD_DIR"
  nicer cmake -S "$REPO" -B "$FDS_BUILD_DIR" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER=mpiicx -DCMAKE_CXX_COMPILER=mpiicpx -DCMAKE_Fortran_COMPILER=mpiifx \
    "-DCMAKE_Fortran_FLAGS_RELEASE=$INTEL_FFLAGS" "-DCMAKE_C_FLAGS_RELEASE=$INTEL_CFLAGS" "-DCMAKE_CXX_FLAGS_RELEASE=$INTEL_CFLAGS" \
    -DAMREX_INSTALL_DIR="$AMREX_PREFIX" -DAMReX_DIR="$AMREX_CMAKE" -DFDS_USE_MKL="$FDS_USE_MKL" \
    -DHYPRE_DIR="$HYPRE_LIBDIR/cmake/HYPRE" -DHYPRE_INCLUDE_DIRS="$HYPRE_PREFIX/include" -DHYPRE_LIBRARIES="$HYPRE_LIBDIR/libHYPRE.a" \
    -DSUNDIALS_DIR="$SUNDIALS_CMAKE" \
    -DUSE_AMREX=ON -DUSE_HYPRE=ON -DUSE_OPENMP=ON -DUSE_SUNDIALS=ON -DUSE_SYSTEM_HYPRE=ON -DUSE_SYSTEM_SUNDIALS=ON \
    -DBUILD_TESTING=OFF > "$FDS_BUILD_DIR.configure.log" 2>&1 || { tail -40 "$FDS_BUILD_DIR.configure.log"; die "FDS configure failed ($FDS_BUILD_DIR.configure.log)"; }
  printf '%s' "$CFG_KEY" > "$FDS_BUILD_DIR/.fds_amrex_cfgkey"
fi
log "FDS-AMReX: make fds -j $NJOBS"
nicer cmake --build "$FDS_BUILD_DIR" --target fds -j "$NJOBS" > "$FDS_BUILD_DIR.build.log" 2>&1 || { tail -40 "$FDS_BUILD_DIR.build.log"; die "FDS build failed ($FDS_BUILD_DIR.build.log)"; }
FDS_AMR=$FDS_BUILD_DIR/fds_amr
[ -x "$FDS_AMR" ] || die "$FDS_AMR was not produced"
log "built $FDS_AMR ($(du -h "$FDS_AMR" | cut -f1))"

# ---------------------------------------------------------------- optional smoke test
if [ "$SMOKE" = 1 ]; then
  SM=${SMOKE_DIR:-$FDS_BUILD_DIR/smoke}
  rm -rf "$SM"; mkdir -p "$SM"
  log "smoke: simple_test, level 0, 1 rank, 10 steps (OMP_NUM_THREADS=1)"
  ( export FDS_AMR NP=1 STEPS=10 MPIRUN=mpiexec RUN_DIR="$SM" OMP_NUM_THREADS=1
    bash "$REPO/Verification/Adaptive_Mesh_Refinement/run_simple_test_l0.sh" ) > "$SM.log" 2>&1 || { tail -30 "$SM.log"; die "smoke run failed ($SM.log)"; }
  tail -5 "$SM.log"
  ls "$SM"/simple_l0* >/dev/null 2>&1 || die "smoke produced no output files in $SM"
  log "smoke: OK (output in $SM)"
fi

log "DONE in $(( $(date +%s) - T0 )) s. Run with: OMP_NUM_THREADS=1 mpiexec -n <ranks> $FDS_AMR case.fds --levels ... (see README_FDS_AMReX.md)"
