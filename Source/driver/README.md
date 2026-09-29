# Source/driver: FDS-AMReX driver (Role 1, milestone M2a)

C++ AmrCore driver plus the Fortran `bind(C)` glue (ADR-001 Option A/C; C++ only for driver and glue, NFR-049). Owner: Role 1
(see `OWNERS.md`). Edits to existing FDS files are delivered as numbered patches in `patches/` and applied by the Chief Architect;
everything they add is inside `#ifdef WITH_AMREX`, so `USE_AMREX=OFF` is unchanged (IR-006).

## Status (S0 + S1)
| Step | State |
|---|---|
| S0 toolchain | done, section "Environment" |
| S1 skeleton | `main.cpp`, `FdsAmr.{H,cpp}`, `fds_mesh_query.f90`, `CMakeLists.txt`; patches 0001, 0002; the skeleton runs the unchanged FDS set-up and builds the level-0 `BoxArray` for `shunn3_32` (1 box) and `shunn3_4mesh_32` (4 boxes, 4 ranks) |
| S2 onward | not started (fields, shim, ghost fill, time loop, pressure hook-up, output) |

## Kernel-facing file rules (M2a)
Every file that touches kernels states in its header: (a) passive scalars (`N_TOTAL_SCALARS` beyond the tracked species) are handled by
`Fields.cpp` (S2; it sizes the `RHO_ZZ` component list from `N_TOTAL_SCALARS`); (b) only uniform Cartesian metrics are used, `R(I)` and
`RRN(I)` are dropped, and `CYLINDRICAL` or `TRNX/TRNY/TRNZ` meshes are rejected (IR-002). No K1/K2 kernels in M2a.

## Environment (S0)
Reference toolchain: Debian 13, GCC/gfortran 14.2.0, Open MPI 5.0.7, CMake 3.31, AMReX 26.09 at `/workspace/amrex-install`
(MPI, OpenMP, FFT on; double precision; 3D; no GPU), FireX HYPRE `63331f19` (reports itself as 2.32.0) at
`/workspace/firemodels-gnu/libs/hypre/63331f19`, SUNDIALS 7.5.0 at `/workspace/firemodels-gnu/libs/sundials/v7.5.0`.

```bash
# 1. Only if gfortran/mpifort are missing from PATH (after a machine reset the apt layer is gone; /workspace persists):
bash /workspace/gnu_ompi/install_gnu_ompi.sh --source=cache     # idempotent, pinned versions, self-check must end with "VERIFICATION PASSED"
# 2. Every shell:
source Source/driver/tests/env.sh     # sources /workspace/gnu_ompi/env_gnu_ompi.sh, sets OMP_NUM_THREADS=1, HWLOC_LIBXML=0, BASELINE, FDS_CMAKE_COMMON
```
`env_gnu_ompi.sh` sets `CC=mpicc CXX=mpic++ FC=mpifort` and prints a warning if gfortran is not 14.2.0 or Open MPI is not 5.0.7.
Run with `mpirun --bind-to none -np N <exe> case.fds` and one OpenMP thread (threads=1 for the gate).

Configure and build, always out of tree, with the patches applied on a copy of the tree (not in the reference tree):
```bash
source Source/driver/tests/env.sh
cmake -S <patched-tree> -B <build-dir> $FDS_CMAKE_COMMON                 # USE_AMREX=OFF: stand-alone FDS, executable <build-dir>/fds
cmake -S <patched-tree> -B <build-dir> $FDS_CMAKE_COMMON -DUSE_AMREX=ON  # AMReX driver, executable <build-dir>/fds_amr
cmake --build <build-dir> -j4
```
`FDS_CMAKE_COMMON` = the reference options (`Release`, GNU `-O3`, OpenMP on) with the offline HYPRE/SUNDIALS copies
(`-DUSE_SYSTEM_HYPRE=ON -DUSE_SYSTEM_SUNDIALS=ON`) and fixed date/version strings. The reference binary was built with the FetchContent
copy of the same HYPRE commit; a scratch rebuild of the unpatched tree with these options reproduces the baseline output bitwise (see 0001 note).
The tree is compiled once per option set: the Fortran build takes about 4 minutes on 4 cores.

### Baseline (FireX `36975d765f`)
- Binary: `/workspace/fds-amr/build/firex-36975d7/ompi_gnu_rel/fds` (documented sha256 `9d3b5991...8050ec`). **The directory
  `/workspace/fds-amr/build/firex-36975d7/` (and `build/env/`) does not exist at present**, so this binary is unavailable and was not re-run. The checked-in output of that binary is still available and is the comparison reference.
- Baseline outputs: `/workspace/fds-amr/vv-runs/baseline/gnu_ompi_firex-36975d7/<case>/` (`manifest.json` per case, `SHA256SUMS` at the top).
  `shunn3_32`, `shunn3_4mesh_32` (+ `__1mesh`, `__sf17`, `__1mesh_sf17`) are present.
- **`csmag_32` has no baseline** there (no directory with that name; the input is `Verification/Turbulence/csmag_32.fds`). It must be
  captured by the V&V lead before any `csmag_32` bitwise claim.

## Layout of this directory
| File | Purpose |
|---|---|
| `main.cpp` | `MPI_Init_thread(FUNNELED)`, `amrex::Initialize`, `fds_setup(0)`, level-0 layout, `amrex::Finalize`, `fds_setup(mode>0)` (FDS end: `MPI_Finalize`, `STOP`) |
| `FdsAmr.{H,cpp}` | `build_level0()`: one box per FDS `&MESH`, owner rank = FDS `PROCESS(NM)`, `Geometry` with periodicity from the FDS PERIODIC vents |
| `fds_mesh_query.f90` | read-only `bind(C)` queries of the FDS mesh data (module `FDS_MESH_QUERY`) |
| `CMakeLists.txt` | included from the top-level file when `USE_AMREX=ON` (patch 0002) |
| `patches/` | numbered patches against the reference tree, each with a `.md` note |
| `tests/` | `env.sh`, `check_off_bitwise.sh` (IR-006), `check_setup_amr.sh` (S1) |

## Fortran/C++ boundary (IR-005)
Update this table with every new crossing. "Owner" is who allocates and frees; "Lifetime" is how long the crossing object is valid.
Fortran arrays are column-major with FDS index bounds; AMReX FABs are column-major with box bounds (the +1 face-index offset of staggered
arrays is handled in S2).

| Crossing (C symbol) | Direction | Type and layout | Owner | Lifetime | Status |
|---|---|---|---|---|---|
| `fds_setup(int mode, const char* fname, double* dt_out)` | C++ -> Fortran | `mode` by value (C_INT); `fname` NUL-terminated (C_CHAR array); `dt_out` scalar C_DOUBLE, may be null only for mode>0 | caller | call | S1; patch 0001 |
| `fds_get_nmeshes()` | Fortran -> C++ | `int` | Fortran | call | S1 |
| `fds_get_mesh(nm, ijk[3], xb[6], &rank, &nonuniform)` | Fortran -> C++ | `nm` 1-based FDS mesh index (value); out: 3 `int`, 6 `double` (`XS,XF,YS,YF,ZS,ZF`), 2 `int` | caller provides storage | call | S1 |
| `fds_get_domain(periodic[3], &cyl, &n_tracked, &n_total, &nranks)` | Fortran -> C++ | 3 `int` flags + 4 `int` | caller provides storage | call | S1 |
| FAB memory bound to `MESHES(NM)%RHO, RHOS, ZZ, U, V, W, ...` (W2 alias, `p1_alias.c`) | C++ -> Fortran | descriptor of an ALLOCATABLE component made to describe FAB memory, FDS lower bounds | AMReX (MultiFab); the alias must be released before Fortran can touch the component | one kernel call per box | planned (S2/S3), isolated in one file |
| Per-box side data (`CELL`, `EXTERNAL_WALL`, face mask, metrics) | Fortran-owned | FDS derived types | Fortran | run | planned (S3) |
| `PressureIface.H` (`amrex::FFT::Poisson` call-through) | C++ -> Role 2 | `MultiFab` RHS in, `H` out | driver | step | planned (S6); signature agreed with Role 2, recorded here |

Rules of the boundary: no Fortran `STOP` or `MPI_Finalize` while AMReX is alive (patch 0001 defers `END_FDS` until after
`amrex::Finalize`); the C++ main owns MPI initialization (FUNNELED); FDS keeps `MPI_COMM_WORLD`.
