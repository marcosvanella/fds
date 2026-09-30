# Source/driver: FDS-AMReX driver (Role 1, milestone M2a)

C++ AmrCore driver plus the Fortran `bind(C)` glue (ADR-001 Option A/C; C++ only for driver and glue, NFR-049). Owner: Role 1
(see `OWNERS.md`). Edits to existing FDS files are delivered as numbered patches in `patches/` and applied by the Chief Architect;
everything they add is inside `#ifdef WITH_AMREX`, so `USE_AMREX=OFF` is unchanged (IR-006).

## Status (S0 to S5)
| Step | State |
|---|---|
| S0 toolchain | done, section "Environment" |
| S1 skeleton | `main.cpp`, `FdsAmr.{H,cpp}`, `fds_mesh_query.f90`, `CMakeLists.txt`; patches 0001, 0002; the skeleton runs the unchanged FDS set-up and builds the level-0 `BoxArray` for `shunn3_32` (1 box) and `shunn3_4mesh_32` (4 boxes, 4 ranks) |
| S2 fields | `Fields.{H,cpp}` MultiFab registry (D-031 ghost widths), +1 face offset, `SideData.{H,cpp}` per-box side data, `fds_box_shim.f90` + `fds_alias.c` (FAB memory seen by Fortran with FDS bounds), IR-005 round-trip tests, IR-007 tile-race skeleton; results in "Tests (S2)" |
| S3 kernel shim | `fds_kernels.f90` (bind(C) wrappers of the UNMODIFIED FDS kernels), `fds_density_split.f90` (generated), `fds_clip_gather.f90` (D-031 gather clip), `TimeStep.H` (T/DT replay), `tests/kernelcheck.cpp` + `tests/run_kernelcheck.sh`, scratch reference dump generator `tests/refdump/`; results in "Kernel checks (S3)"; no existing file edited (no new patch) |
| S4 ghost fill and boundary conditions | `GhostExchange.{H,cpp}` (FillBoundary at the MESH_EXCHANGE codes 1/3/4/5/6, class `BcStep`), `fds_ghost_bc.f90` (OMESH fill + thin wrappers of the unmodified FDS boundary routines), `fds_amrex_hooks.f90` (hook module); patches 0003 (velo), 0004 (wall), 0005 (mesh, DRAFT); plain `--ghost=full/face` modes gated; results in "Ghost fill and boundary conditions (S4)" |
| S5 time loop | `TimeLoop.{H,cpp}` (MAIN_LOOP order, predictor/corrector, global DT), `ExactSum.{H,cpp}` (D-028 fixed-point sums), `fds_step.f90` (state getters and thin wrappers of the unmodified FDS routines), `main.cpp --run`, `tests/compare_run.py`; pressure through `pb::solve_pressure` (Role 2's interface, FFT backend); no existing file edited (no new patch); results in "Time loop (S5)" |
| S6 onward | not started (full pressure hook-up: inhomogeneous boundary data, mixed faces; output) |

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

## Layout of this directory (S3 additions: `fds_kernels.f90`, `fds_density_split.f90`, `fds_clip_gather.f90`, `TimeStep.H`, `tools/gen_density_split.py`, `tests/kernelcheck.cpp`, `tests/run_kernelcheck.sh`, `tests/refdump/`)
| File | Purpose |
|---|---|
| `main.cpp` | `MPI_Init_thread(FUNNELED)`, `amrex::Initialize`, `fds_setup(0)`, level-0 layout, `amrex::Finalize`, `fds_setup(mode>0)` (FDS end: `MPI_Finalize`, `STOP`) |
| `FdsAmr.{H,cpp}` | `assemble_level0()` (pure C++, no FDS needed): one box per mesh description, owner rank, `Geometry` with periodicity; the checks (uniform dx, lattice, no overlap) |
| `FdsSetup.{H,cpp}` | `build_level0()` (asks the FDS set-up for the meshes) and `fds_cell_walls()` (per-cell SOLID and wall codes) |
| `Fields.{H,cpp}` | field table, index maps FDS <-> AMReX, `FdsView`, the `Fields` registry, `fill_ghosts()` |
| `SideData.{H,cpp}` | per-box side data rebuilt from FDS cell data: SOLID, face mask (valid+2), source flag; layout-independent hash |
| `fds_box_shim.f90`, `fds_alias.c` | Fortran side of the FAB alias (bind, release, bounds, element access, explicit-shape pass test) |
| `fds_mesh_query.f90` | read-only `bind(C)` queries of the FDS mesh data (module `FDS_MESH_QUERY`) |
| `GhostExchange.{H,cpp}`, `fds_ghost_bc.f90`, `fds_amrex_hooks.f90` | S4: ghost fill at the exchange codes, FDS boundary routines on OMESH filled from the boxes, hook module for patches 0003-0005 |
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
| Per-box side data (`CELL`, `EXTERNAL_WALL`, `WALL`, metrics) | Fortran-owned | FDS derived types | Fortran | run | used as is by the S3 kernels; the face mask is C++ (`SideData`) |
| `fds_k_verbose(v)`, `fds_k_state(pred, first, icyc, rmin, rmax)`, `fds_k_consts(&tend,&dtfill,&dtmin,&dt0,&nzone)` | C++ -> Fortran | by value / scalar out | caller | call | S3; FDS globals PREDICTOR/CORRECTOR/FIRST_PASS/ICYC/RHOMIN/RHOMAX set before a kernel, constants read once |
| `fds_k_visc(nm, est)` | C++ -> Fortran | `COMPUTE_VISCOSITY(NM, APPLY_TO_ESTIMATED_VARIABLES)`, constant Smagorinsky path | caller | call | S3 |
| `fds_k_dens(t, dt, nm)` | C++ -> Fortran | `MASS_FINITE_DIFFERENCES` + unmodified `DENSITY` (FDS's own per-mesh clip) | caller | call | S3 |
| `fds_k_dens_pre(t, dt, nm)`, `fds_k_dens_post(t, dt, nm)` | C++ -> Fortran | `DENSITY` split at `CHECK_MASS_DENSITY` (`fds_density_split.f90`, generated by `tools/gen_density_split.py` from `mass.f90`); the level gather clip runs in between | caller | call | S3 |
| `fds_k_vflux(t, dt, nm, est)`, `fds_k_div1(t, dt, nm)`, `fds_k_div2(dt, nm)` | C++ -> Fortran | `VELOCITY_FLUX`, `DIVERGENCE_PART_1`, `DIVERGENCE_PART_2` | caller | call | S3 |
| `fds_k_vpred(t, dt, nm, &dtnew, &ichg, &cfl, &vn)`, `fds_k_vcorr(t, dt, nm)` | C++ -> Fortran | `VELOCITY_PREDICTOR` (with its DT_NEW, `CHANGE_TIME_STEP_INDEX`, CFL, VN), `VELOCITY_CORRECTOR` | caller | call | S3 |
| `fds_k_flags(nm)`, `fds_k_set_flags(nm, f)`, `fds_k_set_restrict(nm, n)`, `fds_k_get_restrict(nm)` | both | clip flags (bit0 RHOMIN, bit1 RHOMAX) and `DT_RESTRICT_COUNT` of a mesh | caller | run | S3 |
| `fds_clip_density`, `fds_clip_density_apply`, `fds_clip_species`, `fds_clip_species_one`, `fds_clip_species_apply`, `fds_clip_renorm`, `fds_clip_sums` | C++ -> Fortran | D-031 gather clip (`fds_clip_gather.f90`): FAB pointers with explicit bounds (`RHOP` has its own `RLO:RHI`, so the native ng=3 RHO/RHOS FAB is passed as is), face mask `MASK(:,:,:,1:6)` from `SideData` | caller | call | S3 |
| `fds_k_xfer(nm, name, rnk, lb, ub, data, mode, &n, &nbad, &ierr)`, `fds_k_match(...)`, `fds_k_neutralize(nm)` | C++ -> Fortran | test support: load (0), bitwise compare (1) or bounds query (2) of a named FDS array (fields, wall and edge arrays, zone profiles/sums); wall/edge match table for window mode; neutralise mesh-interface walls of a window | caller | call | S3, test support |
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

## Kernel checks (S3)
The kernels are the unmodified FDS routines (`velo.f90`, `mass.f90`, `divg.f90`), called through `fds_kernels.f90` on `MESHES(NM)` whose arrays are
aliases of the FABs (S2). Only two things are not literally the FDS source: `DENSITY` is also available split at `CHECK_MASS_DENSITY` (generated
text, every other line identical) so that the D-031 gather clip can run between the halves, and the gather clip itself (`fds_clip_gather.f90`).

Decisions recorded at S3:
- Ghosts: the full ghost fill is kept. The FDS kernel results are reproduced bitwise with the full fill AND with the face-neighbour-only fill
  (`--ghost=face+bc` passes exactly like `full+bc`), so a face-only mask option was NOT needed. The edge/corner ghost layers that only the full fill supplies are
  not read by the stage kernels on these cases (the D-031 clip reads valid+1 only through the face mask).
- RHO/RHOS: the native ng=3 FAB with native bounds is passed (never a strided window into an explicit-shape dummy). The check poisons layer 3 of RHO/RHOS
  with NaN before every kernel; all results stay bitwise equal, so no kernel reads or writes outside the FDS window.
- Ghost values that FDS itself sets in its boundary-condition and pressure steps (strips at non-periodic domain sides, edge/corner strips, the H/HS ghost
  cells) are not produced by any S3 kernel. The `+bc` modes take exactly those from the dump; S4 (boundary conditions) and the pressure backend replace them.
  Plain `--ghost=full` (no `+bc`) was not gated at S3; S4 gates it (next section).
- Window mode (`--window`): a SINGLE-mesh dump is the reference for a multi-box layout. The wall of a box without counterpart in the single mesh
  (mesh-to-mesh interface) is made neutral: `BOUNDARY_TYPE=NULL_BOUNDARY`, `CELL%WALL_INDEX=0`, `UVW_SAVE` and the wall state from the gas cell. The zone sums
  DSUM/PSUM/USUM are process-global accumulators and are compared only in native mode. CFL, VN, DT_NEW and the change index of a window are the extrema over boxes.

Scripts: `tests/run_kernelcheck.sh <build> <ref-runs-dir> [work-dir]` (needs the dumps below); `tests/run_driver_tests.sh <build>` is unchanged.
```
fds_amr case.fds --kernelcheck <dump> [--window] [--ghost=dump|full|face|full+bc|face+bc]     # FDSKC_VERBOSE=1: first differing element per array
```
Results (threads=1, gfortran 14.2, Open MPI 5.0.7, all BITWISE-OK, per kernel tag: VISC_P/C, DENS_P/C (unmodified DENSITY, native layouts), DENSCLIP_P/C
(split + gather clip), VFLUX_P/C, DIV1_P/C, DIV2_P/C, VPRED (arrays, DT_NEW, change index, CFL, VN), VCORR, FLAGS (clip flags, restrict count), MASK (face mask against the wall table),
TDT (T/DT replay of the STEP and PASS records), LOAD (load path)): `shunn3_32`, `csmag_32` (constant Smagorinsky; no baseline exists, the reference is the
scratch run), `shunn3_32_clip` (clip active, KVAR=3), `shunn3_4mesh_32__1mesh`, `shunn3_4mesh_32` native (4 ranks, per-mesh dumps), `shunn3_4mesh_32` as four windows of
the 1-mesh dump (4 ranks), the same with the clip active (4-box gather clip against the single-mesh clip). Modes `dump`, `full+bc`, `face+bc`; 1 and 4 ranks.
Not covered: N_ZONE>0 cases (DSUM/USUM/PSUM written but not exercised), non-periodic walls with open/solid BC data beyond the uniform patches of these cases,
the species clip with more than one tracked species (the shunn3 cases have 2 tracked species, csmag_32 one passive scalar).

## Ghost fill and boundary conditions (S4)
Scope: periodic and same-level box-to-box boundaries, uniform Cartesian metrics only (CYLINDRICAL/TRN* rejected at level-0 assembly, IR-002); passive scalars are
carried by `Fields.cpp` (ZZ/ZZS ncomp = N_TOTAL_SCALARS) and handled in the same loops.

### Exchange positions (verified against `main.f90`; the plan's line numbers were stale)
| MESH_EXCHANGE code | Position in MAIN_LOOP | Fields filled by `FillBoundary` (`exchange_fields(code)`) | Then (`BcStep::after_exchange`) |
|---|---|---|---|
| 1 | after the predictor DENSITY (~838) | RHOS, ZZS, MU, KRES, D | `VISCOSITY_BC(est)` (~866) |
| 3 | after the predictor pressure (~957) | US, VS, WS, HS | `MATCH_VELOCITY` (962), `VELOCITY_BC(est)` (969) |
| 4 | after the corrector DENSITY (~1007) | RHO, ZZ, MU, KRES, DS | `VISCOSITY_BC` (~1014) |
| 6 | after the corrector (~1153) | U, V, W, H | `MATCH_VELOCITY` (1168), `VELOCITY_BC(final)` (1174) |
| 5 | in PRESSURE_ITERATION_SCHEME (~1665-1735) | FVX, FVY, FVZ, H (predictor) or HS (corrector) | `MATCH_VELOCITY_FLUX` (1686, returns for one mesh) |
Codes 0, 2, 7-11, 14-20 do not concern M2a. `WALL_BC` is called at ~547, ~886, ~1057 (`BcStep::wall_bc`). `ghost_exchange(F, code, predictor)` is the S5 entry point.

### Decision: the driver runs FDS's own boundary routines (IR-002/ADR-001 Option C, stepwise)
The ghost values FDS writes across a box boundary are not plain copies: the shared face is `0.5*(a + ((b*dy)*dz)/(dy*dz))` (operation order matters), tangential
strips use the `EDGE_INTERPOLATION_FACTOR` weights (1-2e-14 from the set-up rounding), MU/KRES edge cells are clamped copies, the non-periodic y ghosts of RHO/TMP/ZZ come
from `WALL_BC` (`PBAR/(RSUM*TMP)`). Re-deriving these in C++ bit for bit was judged riskier than reusing the source. So in S4 the driver (a) fills the level ghosts with
`FillBoundary` (periodic image and box-to-box), (b) copies each box's FAB, ghost cells included, into `OMESH(NOM)` of the boxes that need it (`fds_g_fill_om`; FABs of
boxes of other ranks are broadcast), (c) calls the UNMODIFIED `MATCH_VELOCITY`, `VELOCITY_BC`, `VISCOSITY_BC`, `WALL_BC`. The values are FDS's own by construction.
Cost: O(domain) broadcast per exchange and the OMESH arrays stay allocated; acceptable for M2a, to be removed by patches 0003/0004 (flag `EXTERNAL_GHOSTS_FILLED`, the
`NOM>0` branches then skip the OMESH writes and read AMReX ghosts) in S5/S6. Item (b) of the task (no dump-supplied BC values) is met: plain modes use no dump ghost data
except the two classes below. Limitation: the driver keeps `OMESH` and `MESHES(NM)%EXTERNAL_WALL`/`EDGE` from the FDS set-up (same-level only).
No `pres.f90`, `init.f90`, `read.f90` change was needed. TRN* is rejected (C++ `assemble_level0`). The MPI_PROCESS message (IR-004, Q9/D-036) is a driver-side
warning in `FdsSetup.cpp`, printed only in AMR mode with more than one rank; it cannot tell whether MPI_PROCESS was given (read.f90 is not patched), M2a keeps the FDS map.

### Plain ghost modes (`--ghost=full`, `--ghost=face`, no `+bc`)
Test design. A frozen dump is a snapshot: the boundary-face strips of US/VS/WS and U/V/W were written by VELOCITY_BC at different points of the step and from interiors
that changed in between (DENSITY restores boundary-face values, `mass.f90` 426/598). The kernel check therefore (1) fills the level ghosts from the dump's valid cells, (2) regenerates
the ghosts the record's kernels read with `BcStep::replay_velocity` (OMESH fill, `VISCOSITY_BC`, `VELOCITY_BC`, plus the MU/KRES edge copy `fds_g_mu_edges`), (3) reloads the wall
arrays from the dump. H/HS (and FVX/FVY/FVZ, WORK*) are kernel-strip arrays: they are loaded whole from the dump; H/HS ghosts hold the image plus the pressure solver's mean offset
and are produced by the pressure step (Role 2, S6).
Results (`tests/run_kernelcheck.sh`, threads=1, 1 and 4 ranks, 7 cases x {dump, full+bc, face+bc, full, face}):
- `dump`, `full+bc`, `face+bc`: all tags BITWISE-OK (unchanged).
- Plain `full` and `face`: the new tags `BCCHAIN_P` / `BCCHAIN_C` are BITWISE-OK on every native case: starting from the state right before the boundary step (VPRED/VCORR "after" arrays),
  FillBoundary + OMESH fill + `MATCH_VELOCITY` + `VELOCITY_BC` reproduce all US/VS/WS (or U/V/W) values including ghost strips of the next record (shunn3_32: 2 P + 1 C; csmag_32: 1 P;
  4-mesh native, 4 ranks: 8 P + 4 C). This is the periodic and same-level box-to-box proof.
- Plain `full`/`face` per-kernel tags are BITWISE-OK except VISC_P/C and VFLUX_P/C on 6 of 7 runs, where 0 to a few dozen elements out of 1.3e6 differ in the last bits (for example csmag_32: 81 elements
  total, MU at 3 cells, STRAIN_RATE 6, FVX/FVY/FVZ 11/21/2; 4-mesh native: 108 in 1.9e6 compares). Loading the dump's boundary-face strips of U/V/W (diagnostic `FDSKC_STRIPDUMP`) removes all of
  them: they come from the replayed strip values of a snapshot that is not the exact pre-boundary-step state, not from a wrong rule (the BCCHAIN test from the exact state is bitwise). The plain gate of `run_kernelcheck.sh` is
  therefore: BCCHAIN bitwise, every other tag bitwise or at most `PLAIN_MAX_PPM` (default 100) parts per million differing elements. Stated openly: plain full/face per-kernel tags are NOT all bitwise.
- `WALL_BC` replay (`FDSKC_WALLBC=1`) is NOT part of the gate: on the frozen snapshots it rewrites RHO/TMP strips (shunn3: 4096 of 5120 RHO/TMP elements differ in VFLUX_P) because the snapshot does not hold the state WALL_BC saw in the step. The
  gated plain modes take RHO/RHOS/TMP/ZZ ghosts from the level fill (periodic image) and, on the non-periodic sides, from the dump snapshot values already in the FAB; the WALL_BC-equivalent for those sides
  (`BcStep::wall_bc`, unmodified `WALL_BC`) is wired at its MAIN_LOOP positions (~547, ~886, ~1057) in S5 and is checked there on the time loop, not on snapshots. Limitation stated: S4 does not prove WALL_BC on snapshots.

### New test tooling
Env (kernel check): `FDSKC_GHOSTDIFF2` (per-class ghost differences against the dump), `FDSKC_GHOSTDIFF3` (print U strip values), `FDSKC_NOCHAIN` (skip BCCHAIN), `FDSKC_WALLBC=1` (replay with WALL_BC),
`FDSKC_STRIPDUMP=1` (strips of the boundary faces from the dump after the replay), `FDSKC_STRIPONLY=<list>` and `FDSKC_RELOAD=<list>` (restrict / extend that reload), plus the S3 poison switches
(`FDSKC_EDGEPOISON`, `FDSKC_STRIPPOISON`, `FDSKC_HPOISON`, `FDSKC_PHYSPERTURB`, `FDSKC_NOFILL`, `FDSKC_NOBC`). `FDSKC_VERBOSE=1` prints up to 8 differing elements per array.

### Patches (send to the Architect)
| Patch | File | Status |
|---|---|---|
| `patches/0003-velo-external-ghosts.patch` | `velo.f90` | guarded, inert (flag FALSE); OFF bitwise checked |
| `patches/0004-wall-external-ghosts.patch` | `wall.f90` | guarded, inert; OFF bitwise checked |
| `patches/0005-mesh-point-to-box.patch` | `mesh.f90` | DRAFT (due end of S5); not validated with oneAPI (not available here), gfortran 14.2 only |
Each has a `.md` note with the evidence. OFF check with all of 0003-0005 applied: `shunn3_32` 1 rank 16 files, `shunn3_4mesh_32` 4 ranks 47 files, bitwise to the baseline. Toolchain: gfortran 14.2 and Open MPI 5.0.7 only; `ifx`/`ifort` are not installed.

### Fortran/C++ boundary additions (S4)
| Crossing (C symbol) | Direction | Type and layout | Owner | Lifetime | Status |
|---|---|---|---|---|---|
| `fds_g_fill_om(nm, nom, which, lb[3], ext[3], nc, data)` | C++ -> Fortran | block copy of box `nom`'s FAB (FDS lower bounds, extents, ncomp) into `OMESH(nom)` of box `nm`; `which` 1..20 (MU, RHO, RHOS, U, V, W, US, VS, WS, H, HS, FVX, FVY, FVZ, D, DS, KRES, Q, ZZ, ZZS); returns 0 or 1 (not a neighbour) | caller (read only) | call | S4 |
| `fds_g_phase(pred)`, `fds_g_match(nm)`, `fds_g_match_flux(nm)`, `fds_g_velocity_bc(t, nm, est)`, `fds_g_viscosity_bc(nm, est)`, `fds_g_wall_bc(t, dt, nm)`, `fds_g_mu_edges(nm)` | C++ -> Fortran | set PREDICTOR/CORRECTOR; thin wrappers of the unmodified routines; MU/KRES edge-cell copy | caller | call | S4 |
| `fds_hook_set_flag(flag)`, `fds_hook_set_view(nmax, nm, which, lb[4], ext[4], p)` | C++ -> Fortran | `EXTERNAL_GHOSTS_FILLED`; `BOX_VIEW` pointer views of FAB memory (`C_F_POINTER` + bounds remapping) | AMReX (FAB) | until released | S4, used by patches 0003-0005 from S5 |

### Reference dump (scratch, outside src)
`tests/refdump/` holds the write-only generator. It is applied to an UNPATCHED copy of the FDS source (a `git archive` of the reference commit), never to this tree:
```
cp -r <unpatched FDS source> <scratch>/reftree && cp tests/refdump/refdump.f90 <scratch>/reftree/Source/
python3 tests/refdump/instrument.py <scratch>/reftree/Source/main.f90       # inserts REFDUMP_* calls; add refdump.f90 to the Fortran sources of the scratch CMake
FDSREF_FILE=ref.dump FDSREF_STEPS=2,3 [FDSREF_RHOMIN=1.4 FDSREF_RHOMAX=4.6] <scratch build>/fds case.fds     # mpirun -np N: one file per mesh, ref.dump, ref.dump.2, ...
python3 tests/refdump/rd.py ref.dump                                            # list records
```
`FDSREF_STEPS` = time-step counts to dump (default 2,3; each dumped step has a predictor and a corrector record per kernel). `FDSREF_RHOMIN/RHOMAX` = stress mode:
the density limits are narrowed so that the clip is active, every array touched is restored (the run ends with an instability stop, the dump is complete).
The instrumented binary reproduces the baseline hrr, mass and restart output bitwise (checked for `shunn3_32` and `shunn3_4mesh_32__1mesh`).

Raw format (stream, unformatted, little endian; REAL = float64, INTEGER = int32):
```
file   : char(8) "FDSDUMP1", int32 NREC (patched at the end), then NREC records
record : char(16) name, int32 ICYC, int32 PREDICTOR, int32 FIRST_PASS, int32 KVAR (clip bits after the kernel: bit0 RHOMIN, bit1 RHOMAX),
         float64 T, DT, RHOMIN, RHOMAX, float64 X(4), int32 NB, int32 NKEEP,
         NB "before" arrays (the kernel input), NKEEP "after" arrays (the arrays the kernel changed)
array  : char(16) name, int32 RANK, int32 LB(4), int32 UB(4), float64 data, Fortran order, FDS index bounds (RHO with the FDS ng=2 allocation)
names  : VISC_P DENS_P VFLUX_P DIV1_P DIV2_P VPRED VISC_C DENS_C VFLUX_C DIV1_C DIV2_C VCORR (kernel records), PASS and STEP (time-step records of mesh 1 only)
X      : kernel records X(3)=CFL, X(4)=VN, VPRED X(1)=DT_NEW, X(2)=change index; STEP: DT_NEW, change index, number of passes; PASS: DT_NEW, change index, pass ordinal
         (T and DT of a PASS/STEP record are the values used); DENS_P is written only at FIRST_PASS
arrays : RHO RHOS TMP U V W US VS WS H HS KRES D DS DDDT MU MU_DNS RSUM FVX FVY FVZ STRAIN_RATE D_Z_MAX WORK1;  rank 4: ZZ ZZS DEL_RHO_D_DEL_Z SWORK4 FX FY FZ;
         PBAR PBAR_S D_PBAR_DT D_PBAR_DT_S; DSUM USUM PSUM (N_ZONE>0 only); CLIPFL DTRC (1 element: clip flags, DT_RESTRICT_COUNT);
         wall arrays of length NW: W_UVW_SAVE W_U_NORMAL W_U_NORMAL_S W_RHO_F W_TMP_F W_ZZ_F W_RHO_D_DZDN_F W_RHO_D_F (NW x NS) W_K_G W_Q_CON_F W_Q_LEAK
         W_IOR W_IIG W_JJG W_KKG W_U_GHOST W_V_GHOST W_W_GHOST; edge arrays: E_IJKA (NE x 4) E_OMEGA E_TAU (NE x 5, index -2..2)
```

## Time loop (S5)

Scope of `--run`: uniform Cartesian meshes, no obstructions, no particles, no radiation, periodic or Neumann-closed Poisson boundaries with homogeneous data
(`shunn3_32`, `shunn3_4mesh_32` and the decomposition cases). Other cases are refused at start (the reason is printed). Kernel-facing rules as above: passive
scalars are handled by `Fields.cpp`; only uniform Cartesian metrics are used (IR-002).

```
source Source/driver/tests/env.sh; mpirun -np <n> fds_amr case.fds --run [--steps N] [--outdir D] [--chid C] [--quiet]
python3 Source/driver/tests/compare_run.py <run-dir> <baseline-dir> <chid> [--dump ref.dump]   # T/DT replay, MMS errors (T2), final fields, mass
```

`TimeLoop::advance` follows MAIN_LOOP (main.f90 predictor ~ 695-960, corrector ~ 1000-1180, current numbering): COMPUTE_VISCOSITY; per pass DENSITY (pre-clip half, level
gather clip of D-031, post-clip half), exchange code 1, WALL_BC, DIVERGENCE_PART_1, zone sums, DIVERGENCE_PART_2, pressure scheme (BAROCLINIC, exchange of FV, NO_FLUX,
RHS, solve, residual, velocity error and the iteration exit tests, main.f90 1649-1793), VELOCITY_PREDICTOR with CHECK_STABILITY; the pass is repeated while the
step is restricted (the global DT decision is the MIN over boxes and ranks by an Allreduce; MAX is applied through the same path); then exchange code 3, the corrector
chain with codes 4 and 6, the end-of-step bookkeeping, the MMS output at the first step with T >= MMS_TIMER, mass file.

### Exact sums (D-028)
`ExactSum.{H,cpp}`: terms are converted to integers at one power-of-two scale per group (derived from the global maximum term, a max reduction), added in 128-bit integers,
combined over ranks limb-wise (integer sum reduction) and converted back once. The result depends on the multiset of terms only (not order, box split or rank count).
Users: setup pressure-zone volume (set into FDS `P_ZONE%VOLUME`), exact vent area (diagnostic against FDS `TOTAL_FDS_AREA`: 4e-16 relative), and the DSUM/PSUM/USUM
zone integrals each pass (`fds_p_zone_terms` returns the per-cell and per-wall terms of DIVERGENCE_PART_1 with the same expressions and skip rules; they are summed exactly
and written back with `fds_p_zone_set`). Unit test `exact_sum` (tests/test_units.cpp, in `run_driver_tests.sh`): 32 x 4 x 32 periodic level cut 1, 2x2, 4x4, 8x8 (boxes of 32, 16, 8, 4
cells) on 1 and 4 ranks; every split returns the same bits, equal to a serial int128 reference; the plain double sum over the same splits DIFFERS; cancellation and constant-field
checks are exact. 20 checks on 1 rank, 80 on 4 ranks, 0 failures.

### Results (shunn3_32, 1 rank, 80 steps to T = 1, threads = 1)
* T/DT replay: all 17 rows of the baseline `_steps.csv` match (DT to the 3 printed digits, T to 7). Against the dump STEP records: **0 of 80 bitwise**; max relative DT
  difference 5.0e-14, max relative T difference 6e-16. Cause: the FDS Crayfishpak solve is replaced by `FFT::Poisson` (H differs at ~1e-14, D-021). Pressure iteration counts equal FDS (10/10 on steps 1-2, then 7/7; passes > 1 at steps 1, 55, 57 as in the dump).
* MMS errors at T = 0.9003 (baseline / driver): e_rho 3.3002e-02 / 3.3002e-02, e_Z 6.2455e-03 / 6.2455e-03, e_u 4.7980e-03 / 4.7980e-03, e_H 1.8152e-02 / 1.8152e-02. T2 (`e <= 1.05 e_base`) passes; the Tol_FDS part is not applicable at N = 32.
* Final fields against the baseline restart, max |d|: U 1.3e-15, V 0, W 1.2e-15, H 1.8e-14, HS 2.2e-14, D 2.2e-14, DS 1.9e-14, RHO 7.1e-15, TMP 3.4e-13 (relative to the field maxima: 1e-12 at worst).
* Mass: total mass 1.2 exactly (equal to the baseline). Zone-sum fixed-point versus FDS's own rank-local accumulation: 1e-13 to 2e-12 relative (FDS's accumulation is the order dependent one).
* `shunn3_4mesh_32` (4 boxes, 1 and 4 ranks, 82 steps): the level-wide periodic FFT replaces FDS's multi-mesh interpolated-boundary iteration, so this is T2 only. Against
  `shunn3_4mesh_32__1mesh` (MMS at T = 1): e_rho 3.2251e-02 / 3.2255e-02, e_Z 1.3710e-02 / 1.3709e-02, e_u 4.9516e-03 / 5.0145e-03, e_H 7.4672e-03 / 8.6799e-03 (driver / baseline): T2 passes.
  Final fields differ from the single-mesh restart by up to 1e-3 (U, W), 2.7e-2 (H), T/DT differ from step 1 on (6.191e-4 against 6.188e-4); this is the different pressure coupling, not a defect of the loop. Mass is conserved exactly by the driver (the multi-mesh baseline drifts by 2.8e-5).

### Decomposition independence (honest state)
Cases `dec{1,2,4}`: the same periodic 32 x 1 x 32 problem as 1, 4 and 16 boxes (max_grid_size 32, 16, 8), ranks 1, 2, 4, all run to T = 1 (82 steps).
* The exact sums are decomposition independent (unit test above) and the iteration counts are the same in all runs.
* The whole run is NOT bitwise independent. Across rank counts at fixed box split the fields differ by <= 7e-13 (last-bit level: the distributed FFT sums in a layout dependent order); across box splits (1 against 4 or 16 boxes) by ~3e-5 in U, H, RHO at T = 1 (same T2 errors).
* Stage-by-stage comparison of the first cycle (`FDSTL_STAGE=1` dumps the fields after every stage): density, species, viscosity and VELOCITY_FLUX are bitwise equal across the splits (<= 2e-16 in FVX/FVZ);
  the first difference is in `DIVERGENCE_PART_1`: its temperature-gradient/enthalpy terms at box-interface cells read the wall arrays of the interface (`INTERPOLATED_BOUNDARY` walls with `NOM > 0`), which FDS keeps for the mesh-to-mesh coupling (UVW_SAVE, TMP_F, RHO_F, K_G from `WALL_BC`/`ASSIGN_GHOST_VALUE`), not the level neighbour values, so a 1-box level and a split level differ at the 1e-2 level in DS on the interface columns.
  Making the interfaces no-wall faces (`FDSTL_OPEN=1`, `fds_p_open_interfaces`, the `fds_k_neutralize` of the kernel check) removes most of it in DS (from 5.7 to 0.36 max) but the NOM > 0 walls are still needed by NO_FLUX and MATCH_VELOCITY_FLUX in the pressure stage and the run then diverges from the single box by O(1): **not adopted, left as opt-in diagnostic**. The clean fix is the externally-filled mode of patches 0003/0004 (`EXTERNAL_GHOSTS_FILLED`) with the TMP/RSUM ghost recompute noted in patch 0004, plus a
  skip of the NOM > 0 walls in DIVERGENCE_PART_1 and the wall arrays of the interface; that needs existing-file edits and is proposed for S6 (next patch number 0006), not done here.

### WALL_BC wiring and validation on snapshots
WALL_BC (unmodified) is called for every box at the positions of main.f90 (after exchange code 1 and 4 of each stage, `BcStep::wall_bc`), after the OMESH fill. On frozen snapshots
(`run_kernelcheck.sh`'s `--ghost=full` mode with `FDSKC_WALLBC=1 FDSKC_WALLCMP=1`) the wall arrays the driver produces are compared bitwise with the dump's per record tag:
`W_K_G`, `W_Q_CON_F`, `W_Q_LEAK`, `W_RHO_D_DZDN_F`, `W_UVW_SAVE`, `W_U_GHOST`, `W_V_GHOST`, `W_W_GHOST`, `W_U_NORMAL`, `W_U_NORMAL_S`, `W_RHO_F`, `W_TMP_F`, `W_ZZ_F` are bitwise equal in the DENS tags;
`W_RHO_D_F` differs (16384 of 17408 elements in DENS tags; 8192 of 8704 in DIV tags), and in the DIV1/DIV2 tags `W_RHO_F` (4096 of 4352), `W_TMP_F` (256) and `W_ZZ_F` (8192 of 8704) differ.
These are snapshot artefacts: the dump is taken after later kernels rewrote the wall arrays (RHO_D_F is an output of DENSITY/DIVERGENCE, RHO_F/TMP_F/ZZ_F are rewritten by DENSITY after WALL_BC), so WALL_BC on the frozen state is not a pure replay. The decisive validation is the running loop
(80 steps; V identical to the baseline, U/W 1e-15, MMS equal to 5 digits).

### Fortran/C++ boundary additions (S5)
`fds_step.f90` (module `FDS_STEP`): `fds_p_params`, `fds_p_mesh_info`, `fds_p_get`, `fds_p_mfd`, `fds_p_dens(_pre)`, `fds_p_init_div`, `fds_p_zone_get/set/terms`, `fds_p_zone_volume(_terms/_set)`,
`fds_p_vent_count/get`, `fds_p_iter_*`, `fds_p_baroclinic`, `fds_p_noflux`, `fds_p_rhs`, `fds_p_get_prhs`, `fds_p_bmax` (boundary data check, per face flags), `fds_p_h_ghost`
(ghost assignment block of PRESSURE_SOLVER_FFT, face flag order xlo, xhi, ylo, yhi, zlo, zhi), `fds_p_resid`, `fds_p_velerr`, `fds_p_get_err`, `fds_p_set_wall_counter`, `fds_p_stop_status`,
`fds_p_clear_attached`, `fds_p_total_iter`, `fds_p_open_interfaces` (opt-in). `FdsSetup.cpp`: `dom.periodic/cylindrical` are max-reduced over ranks (FDS sets `PERIODIC_DOMAIN_*` only on ranks that own a periodic-vent mesh).

### Pressure call-through
`TimeLoop::solve_poisson` -> `pb::solve_pressure` (FFT backend). Review and requests for Role 2 / the Architect: `notes/pressure-iface-review.md`.

### Debug switches (environment variables, off by default)
`FDSTL_DEBUG` (per-iteration errors), `FDSTL_PBV` (pressure backend verbose and residual probes), `FDSTL_ZONES` (print the zone sums), `FDSTL_STAGE=<icyc>` (per-stage field dumps of a cycle), `FDSTL_OPEN=1` (see above).

### Patches (S5)
None: no existing FDS file was edited in S5 (the OFF build is untouched; the bitwise OFF check of 0001-0005 stands). oneAPI (ifx/ifort) is not available on the box: gfortran 14.2 + Open MPI 5.0.7 only.
