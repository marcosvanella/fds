# Source/driver: FDS-AMReX driver (Role 1, milestone M2a)

C++ AmrCore driver plus the Fortran `bind(C)` glue (ADR-001 Option A/C; C++ only for driver and glue, NFR-049). Owner: Role 1
(see `OWNERS.md`). Edits to existing FDS files are delivered as numbered patches in `patches/` and applied by the Chief Architect;
everything they add is inside `#ifdef WITH_AMREX`, so `USE_AMREX=OFF` is unchanged (IR-006).

## Status (S0 to S2)
| Step | State |
|---|---|
| S0 toolchain | done, section "Environment" |
| S1 skeleton | `main.cpp`, `FdsAmr.{H,cpp}`, `fds_mesh_query.f90`, `CMakeLists.txt`; patches 0001, 0002; the skeleton runs the unchanged FDS set-up and builds the level-0 `BoxArray` for `shunn3_32` (1 box) and `shunn3_4mesh_32` (4 boxes, 4 ranks) |
| S2 fields | `Fields.{H,cpp}` MultiFab registry (D-031 ghost widths), +1 face offset, `SideData.{H,cpp}` per-box side data, `fds_box_shim.f90` + `fds_alias.c` (FAB memory seen by Fortran with FDS bounds), IR-005 round-trip tests, IR-007 tile-race skeleton; results in "Tests (S2)" |
| S3 onward | not started (kernel shim, stage kernels, time loop, pressure hook-up, output) |

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
| `FdsAmr.{H,cpp}` | `assemble_level0()` (pure C++, no FDS needed): one box per mesh description, owner rank, `Geometry` with periodicity; the checks (uniform dx, lattice, no overlap) |
| `FdsSetup.{H,cpp}` | `build_level0()` (asks the FDS set-up for the meshes) and `fds_cell_walls()` (per-cell SOLID and wall codes) |
| `Fields.{H,cpp}` | field table, index maps FDS <-> AMReX, `FdsView`, the `Fields` registry, `fill_ghosts()` |
| `SideData.{H,cpp}` | per-box side data rebuilt from FDS cell data: SOLID, face mask (valid+2), source flag; layout-independent hash |
| `fds_box_shim.f90`, `fds_alias.c` | Fortran side of the FAB alias (bind, release, bounds, element access, explicit-shape pass test) |
| `fds_mesh_query.f90` | read-only `bind(C)` queries of the FDS mesh data (module `FDS_MESH_QUERY`) |
| `CMakeLists.txt` | included from the top-level file when `USE_AMREX=ON` (patch 0002) |
| `patches/` | numbered patches against the reference tree, each with a `.md` note |
| `tests/` | `env.sh`, `check_off_bitwise.sh` (IR-006), `check_setup_amr.sh` (S1); S2: `test_units.cpp` (`driver_unit_tests`, no FDS), `selftest_fds.cpp` (`fds_amr case.fds --selftest`), `check.H`, `check_inventory.py`, `run_driver_tests.sh` |

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
| `fds_get_cell_walls(nm, n[3], flags[0:6][k][j][i])` | Fortran -> C++ | `int` array, first index fastest: 0 = SOLID, 1..6 = wall code of faces -x,+x,-y,+y,-z,+z (0 none, 1 wall, 2 mesh-to-mesh interface) | caller provides storage | call | S2 |
| `fds_shim_bind(nm, name, base, lb[4], ext[4], stride[4])`, `fds_shim_release(nm, name)` | C++ -> Fortran | makes `MESHES(NM)%<name>` (3-D, or 4-D for ZZ/ZZS) describe FAB memory with the FDS lower bounds and element strides (W2 alias, `fds_alias.c`); bind deallocates the FDS array first | AMReX (MultiFab); release before Fortran can touch the component | one kernel call per box | S2, tested |
| `fds_shim_bounds`, `fds_shim_access`, `fds_shim_pass_test` | C++ -> Fortran | query bounds, get/set one element with FDS indices, pass to an explicit-shape dummy and report the address | caller | call | S2, test support |
| Per-box side data (`CELL`, `EXTERNAL_WALL`, metrics) | Fortran-owned | FDS derived types | Fortran | run | planned (S3); the face mask is C++ (`SideData`) |
| Per-box side data (`CELL`, `EXTERNAL_WALL`, face mask, metrics) | Fortran-owned | FDS derived types | Fortran | run | planned (S3) |
| `PressureIface.H` (`amrex::FFT::Poisson` call-through) | C++ -> Role 2 | `MultiFab` RHS in, `H` out | driver | step | planned (S6); signature agreed with Role 2, recorded here |

Rules of the boundary: no Fortran `STOP` or `MPI_Finalize` while AMReX is alive (patch 0001 defers `END_FDS` until after
`amrex::Finalize`); the C++ main owns MPI initialization (FUNNELED); FDS keeps `MPI_COMM_WORLD`.

## Field registry (S2)
- Table: `Fields.cpp` `field_table()` (25 arrays). Cell-centred ng=2 (TMP, ZZ, ZZS), ng=1 for the other cell arrays; U/V/W and US/VS/WS nodal in one direction, ng=1;
  RHO/RHOS ng=3; ZZ/ZZS carry `N_TOTAL_SCALARS` components (this is where passive scalars are handled). FX/FY/FZ, ADV_*, DIF_*, SWORK*
  are per-box scratch and are refused by the registry. `tests/check_inventory.py` compares the table with `mesh_fields.csv`
  (`FDS_INVENTORY_CSV=<path> tests/run_driver_tests.sh`): 25 of 25 declared bounds match.
- Index maps: AMReX cell index = box lower corner + FDS index - 1; for a face array the nodal direction gets +1 more (FDS `U(I)` is the
  face at the high side of FDS cell I, AMReX node `lo+I`; FDS `U(-1)` and `U(0)` are the ghost and the low boundary face). `to_amrex`/`to_fds`
  in `Fields.H` are the single place for it; the IR-005 test round-trips every element of every array and compares with the FDS window.
- RHO/RHOS: native ng=3, FDS allocates 2. The FDS window is a strided sub-block of the FAB (`FdsView` carries strides). gfortran 14 honours
  element access through such a strided alias descriptor (test `strided`, INFO line). Passing it to an explicit-shape dummy (how the kernels take arrays) makes gfortran
  build a temporary, and in the S2 test the write did not reach the FAB (INFO line): S3 must pass the full native FAB with native bounds (or a ng=2 copy) to kernels that take RHO/RHOS by explicit shape.
- Ghost fill: `Fields::fill_ghosts` does a full fill (faces, edges, corners, periodic images). A face-only fill (FillBoundary with `cross=true`)
  is NOT offered: in AMReX 26.09 the `cross` flag is applied only to messages between ranks, not to same-rank copies, so the result
  depends on the box-to-rank map (seen in the S2 test on 1 rank versus 4 ranks). OPEN QUESTION for the Architect: the D-031 redundant
  valid+1 clip reads RHO/RHOS/ZZ edge and corner ghosts; full fill provides them, but the FDS MESH_EXCHANGE does not, so
  bitwise agreement with FDS may need the clip restricted to face ghosts (decide at S3, K-level DENSITY).
- Side data: `SideData` (iMultiFab, 8 components, ng=2): 0 SOLID, 1..6 face mask (1 = wall, 0 = open), 7 source flag (1 only inside the domain).
  Mesh-to-mesh interface faces (wall cell in FDS) are open in the mask; wall faces at a domain edge, periodic edges included, stay closed
  (single-mesh FDS has a wall cell there as well). The rebuild is a full recompute from the FDS cell data (`rebuild` = construct again).

## Tests (S2)
```
source Source/driver/tests/env.sh
# out-of-tree build: cmake -S <tree with patches 0001+0002> -B <build> $FDS_CMAKE_COMMON -DUSE_AMREX=ON && cmake --build <build> -j4
FDS_INVENTORY_CSV=<path>/mesh_fields.csv Source/driver/tests/run_driver_tests.sh <build> [--thread-sweep]
```
`driver_unit_tests` (no FDS): IR-005 round trip FDS -> AMReX -> FDS for every registered array on 4 layouts, ghost widths and staggering,
lower-bound remap, window versus native FAB, ghost fill, side data, registry; IR-007 skeleton: a cell 7-point kernel and an x-face kernel over
tile sizes and layouts, bitwise equal and every element written once. `--thread-sweep` repeats it with 2, 4, 8 OpenMP threads (not part of
the threads=1 gate). The skeleton uses stand-in kernels; the real stage kernels replace them from S3.
`fds_amr case.fds --selftest` (FDS-linked): bounds of the FDS allocations equal the table, Fortran read/write through the alias in both
directions, explicit-shape pass, side data from FDS cell data (hash equal for 1 mesh and 4 meshes of the same domain).
