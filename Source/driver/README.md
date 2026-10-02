# Source/driver: FDS-AMReX driver (Role 1, milestone M2a)

C++ AmrCore driver plus the Fortran `bind(C)` glue (ADR-001 Option A/C; C++ only for driver and glue, NFR-049). Owner: Role 1
(see `OWNERS.md`). Edits to existing FDS files are delivered as numbered patches in `patches/` and applied by the Chief Architect;
everything they add is inside `#ifdef WITH_AMREX`, so `USE_AMREX=OFF` is unchanged (IR-006).

## Status (S0 to S9)
| Step | State |
|---|---|
| S0 toolchain | done, section "Environment" |
| S1 skeleton | `main.cpp`, `FdsAmr.{H,cpp}`, `fds_mesh_query.f90`, `CMakeLists.txt`; patches 0001, 0002; the skeleton runs the unchanged FDS set-up and builds the level-0 `BoxArray` for `shunn3_32` (1 box) and `shunn3_4mesh_32` (4 boxes, 4 ranks) |
| S2 fields | `Fields.{H,cpp}` MultiFab registry (D-031 ghost widths), +1 face offset, `SideData.{H,cpp}` per-box side data, `fds_box_shim.f90` + `fds_alias.c` (FAB memory seen by Fortran with FDS bounds), IR-005 round-trip tests, IR-007 tile-race skeleton; results in "Tests (S2)" |
| S3 kernel shim | `fds_kernels.f90` (bind(C) wrappers of the UNMODIFIED FDS kernels), `fds_density_split.f90` (generated), `fds_clip_gather.f90` (D-031 gather clip), `TimeStep.H` (T/DT replay), `tests/kernelcheck.cpp` + `tests/run_kernelcheck.sh`, scratch reference dump generator `tests/refdump/`; results in "Kernel checks (S3)"; no existing file edited (no new patch) |
| S4 ghost fill and boundary conditions | `GhostExchange.{H,cpp}` (FillBoundary at the MESH_EXCHANGE codes 1/3/4/5/6, class `BcStep`), `fds_ghost_bc.f90` (OMESH fill + thin wrappers of the unmodified FDS boundary routines), `fds_amrex_hooks.f90` (hook module); patches 0003 (velo), 0004 (wall), 0005 (mesh, DRAFT); plain `--ghost=full/face` modes gated; results in "Ghost fill and boundary conditions (S4)" |
| S5 time loop | `TimeLoop.{H,cpp}` (MAIN_LOOP order, predictor/corrector, global DT), `ExactSum.{H,cpp}` (D-028 fixed-point sums), `fds_step.f90` (state getters and thin wrappers of the unmodified FDS routines), `main.cpp --run`, `tests/compare_run.py`; pressure through `pb::solve_pressure` (Role 2's interface, FFT backend); no existing file edited (no new patch); results in "Time loop (S5)" |
| S6 decomposition independence, pressure hook-up | box-interface walls handled in the driver (`TimeLoop::iface`, `fds_p_iface_walls`), the `EXTERNAL_GHOSTS_FILLED` path (patches 0003/0004) is the default route (`fds_p_save_uvw`), TMP/RSUM exchanged with the level ghosts; `tests/run_decomp_check.sh` (cases `tests/cases/dec*.fds`); `FDS_HOOK_SET_VIEW`/`POINT_TO_BOX` unit test in `--selftest`; plain kernel-check modes regenerated bitwise (`+strips`); `FDSTL_PDUMP` hand-over for Role 2's one-solve check; no new patch (no 0006 needed); results in "Decomposition independence (S6)" |
| S6b csmag_32 | root cause of the csmag_32 H discrepancy, three driver fixes (periodic face match, domain-edge MU/KRES), `tests/run_csmag_check.sh`; section "csmag_32 investigation" |
| S7 output writers, measurements, gate report | FDS's own writers driven per step through `fds_setup(mode=3)` (patch 0006): `CHID_devc.csv`, `_hrr.csv`, `_mass.csv`, `_steps.csv`, `_cpu.csv`, `.out`; field dump `<chid>_final_<FIELD>.bin` + `_final_manifest.txt`; `<chid>_driver_perf.csv` (loop time, pressure share f_pres, peak RSS); `tests/run_outputs_check.sh`, `tests/run_perf.sh`, `tests/perf_one.py`; `notes/m2a-gate-report.md` (D-006 gate table) |
| S8 M2a closure | OFF bitwise re-check on the committed tree, the three remaining gate items against the official baselines (`tests/run_m2a_baselines.sh`, `tests/m2a_compare.py`), regression cases for the three periodic-only bugs (`tests/run_periodic_regression.sh`), diagnostics (`FDSTL_PROFILE`, `FDSTL_STAGES`, `FDSTL_SKIP_FIX`), 128^3 periodic case generator (`tests/make_periodic_case.py`); section "S8"; answers to Role 3 and the Integration Lead in `notes/interface-answers.md` |
| S9 per-level interface, WP2 baseline, D-053 | `LevelRegistry.{H,cpp}`, `Level` (alias `Level0`), per-level `TimeLoop` stage entry points, global dt over levels, uncovered-cell exact sums, coarse-fine ghost hook, per-level `SideData` rebuild (`notes/level-interface.md`); zone sums in FDS order by default with `--exact-zone-sums` as the GPU option; `tests/run_big_periodic_baseline.sh` and `notes/wp2-periodic-baseline.md` |
| after S9 | flux hooks (task 3) and the wall-state ownership seam (task 4): not started in this entry; inhomogeneous Poisson data and mixed faces: see `notes/pressure-iface-review.md` |

### Scope of the bitwise claims (architect ruling)
Bitwise equality is claimed ONLY for (a) the `USE_AMREX=OFF` build against the baseline (`tests/check_off_bitwise.sh`) and (b) the frozen-input kernel tests (`tests/run_kernelcheck.sh`: every wrapped FDS kernel, dump/full+bc/face+bc, and the exact-strip plain runs). Full time steps with the FFT pressure backend are T2-only: 0 of 80 STEP records of `shunn3_32` are bitwise (FFT::Poisson replaces the Crayfishpak solve, H differs at ~1e-14), and the ~1e-3 difference of the multi-mesh case against its single-mesh twin is the different pressure coupling. Obstruction masks stay NotBuilt: cases with OBST are refused.

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
**Intel oneAPI builds** (ifx 2026.1, Intel MPI; validation of patches 0005 and 0006 by the Intel Build Chief, whose own README is outside this repository). Facts to follow, as reported to the project: (1) configure with `-DCMAKE_CXX_FLAGS=-fp-model=precise`, otherwise the unit test `tile_race` fails one check (the C++ floating-point contraction differs); (2) use `mpiicpx` as the C++ compiler wrapper (with `mpiifx` for Fortran); (3) do not pass `--oversubscribe` to `mpirun`/`mpiexec` under Intel MPI (the test scripts in `tests/` pass it for Open MPI and must be adapted); (4) `-ffree-line-length-none` is set for the driver Fortran files only when the Fortran compiler is GNU (generator expression in `CMakeLists.txt`). Not tested on this box (no oneAPI here).

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
| `patches/0007-mesh-fine-level-boxes.patch` | `mesh.f90` | DRAFT (D-056 option B, applies on top of 0005): second module array of fine-level mesh objects by level, `POINT_TO_MESH_OBJECT`, abort on fine mesh numbers; guarded; not oneAPI validated; OFF bitwise checked |
| `patches/0008-kernels-point-to-box.patch` | `velo.f90`, `mass.f90`, `divg.f90`, `wall.f90`, `turb.f90` | DRAFT with 0007 (D-056 option B): `#define POINT_TO_MESH POINT_TO_BOX` under `WITH_AMREX` before `MODULE` in each file (38 call sites); not oneAPI validated; OFF bitwise checked |
| `patches/0006-main-step-outputs.patch` | `main.f90` | guarded (`WITH_AMREX`), new in S7: `FDS_SETUP(MODE=3)` runs the end-of-step FDS output sequence; OFF bitwise checked with 0001-0006 (`shunn3_32` 16 files, `shunn3_4mesh_32` 47 files); not validated with oneAPI |
Each has a `.md` note with the evidence. OFF check with all of 0003-0006 applied: `shunn3_32` 1 rank 16 files, `shunn3_4mesh_32` 4 ranks 47 files, bitwise to the baseline. Toolchain: gfortran 14.2 and Open MPI 5.0.7 only; `ifx`/`ifort` are not installed.

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
Users: setup pressure-zone volume (set into FDS `P_ZONE%VOLUME`), exact vent area (diagnostic against FDS `TOTAL_FDS_AREA`: 4e-16 relative), and, **only with `--exact-zone-sums`
(env `FDSTL_EXACT_ZONES=1`)**, the DSUM/PSUM/USUM zone integrals each pass (`fds_p_zone_terms` returns the per-cell and per-wall terms of DIVERGENCE_PART_1 with the same expressions and skip rules; they are summed exactly
and written back with `fds_p_zone_set`). **D-053 (from S9): the default keeps FDS's summation order for USUM/DSUM/PSUM** (the kernels accumulate per rank in box order, the driver adds the MPI_SUM
reduction that `main.f90` does, and writes the reduced values back); the exact sum is a switchable option, off by default. Consequences, measured: the default changes the final fields of the
M2a cases by at most 4e-13 relative to the exact mode (`dec1`, `dec2`, `dec4`); exact mode reproduces the previous driver bitwise (final fields and `_driver_steps.csv` byte-identical for dec1, dec2, dec4/4 ranks);
the decomposition check runs its bitwise first-step stages with the exact option and checks the default separately (`tests/run_decomp_check.sh`); `zone_rel` in the step log is computed only in exact mode or with `FDSTL_ZONE_DIAG=1`
(`FDSTL_ZONES=1` prints the zone sums in either mode). Set-up volumes and areas stay D-028 sums in both modes. Unit test `exact_sum` (tests/test_units.cpp, in `run_driver_tests.sh`): 32 x 4 x 32 periodic level cut 1, 2x2, 4x4, 8x8 (boxes of 32, 16, 8, 4
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

### Decomposition independence (S6)
Cases `tests/cases/dec{1,2,4}.fds` (+ `dec{2,4}_np4.fds` with `MPI_PROCESS`): the same periodic 32 x 1 x 32 problem as 1, 4 and 16 boxes; `tests/run_decomp_check.sh <build>` runs them at 1 and 4 ranks.

Cause (found by stage-by-stage comparison of cycle 1, `FDSTL_STAGE=1`): FDS couples boxes through `INTERPOLATED_BOUNDARY` walls with `NOM > 0`; those walls (a) make `DIVERGENCE_PART_1`,
`VELOCITY_FLUX` (`ENTHALPY_ADVECTION_NEW`) and `VELOCITY_BC` add wall-face terms from `B1%RHO_F`, `TMP_F`, `ZZ_F`, `UVW_SAVE` and the interface edge data (`EDGE%OMEGA/TAU`) where a single mesh has plain
interior faces; (b) fill the TMP/RSUM ghost layer across the interface from the OMESH average (`ASSIGN_GHOST_VALUE`), which the level ghost fill does not do; (c) FDS's own set-up ends with
MATCH_VELOCITY / VISCOSITY_BC / VELOCITY_BC / an initial DIVERGENCE_PART_1 on the multi-mesh state, which leaves interface data that a 1-box level does not have.

Fix, driver only (no FDS file edited, no patch 0006):
* `fds_p_iface_walls` (`fds_step.f90`) turns the box-interface walls into `NULL_BOUNDARY` (and `CELL%WALL_INDEX = 0`) around `VELOCITY_FLUX`, `DIVERGENCE_PART_1` and `VELOCITY_BC`, and restores them afterwards
  (`TimeLoop::iface`; the periodic images at the domain edge stay walls, they are walls of the single mesh too). `FDSTL_IFACE=0` disables it (S5 behaviour), `=2` also converts the domain-edge images.
* The initial state is redone with the interfaces as regular faces (`visc`, exchange codes 4 and 6 with their boundary step, initial `DIVERGENCE_PART_1`).
* TMP and RSUM are exchanged with the level ghost fill at the codes 1 and 4.
* `EXTERNAL_GHOSTS_FILLED` (patches 0003/0004, applied in the tree) is now set by the driver (`fds_hook_set_flag`, default on; `FDSTL_EXTGHOST=0` returns to the OMESH-average route). FDS then skips the OMESH writes
  of ASSIGN_GHOST_VALUE, VISCOSITY_BC, VELOCITY_BC, NO_FLUX, the H ghosts and MATCH_VELOCITY(_FLUX); what MATCH_VELOCITY leaves for later routines (`UVW_SAVE`, `BOUNDARY_TYPE_PREVIOUS`) is set by
  `fds_p_save_uvw` from the box's own face value (the two boxes share the face: the average with the neighbour is the value itself). Result: the 1-box run with the flag on equals the flag-off run to 1e-14 (last bits),
  and the split runs equal the 1-box run as below.

Results (first cycle, both passes, `FDSTL_STAGE=1`; 1 box against 4 boxes (1 and 4 ranks) and 16 boxes (1 and 4 ranks)):
* Bitwise equal (byte-identical): RHOS, ZZS (density), FVX, FVZ, MU (VELOCITY_FLUX), MU, KRES, TMP, RSUM (after DIVERGENCE_PART_1), for pass 1 and pass 2, for every split and rank count tested.
* DS, D, DDDT (DIVERGENCE_PART_1/2): bitwise equal for the 4-box splits; for the 16-box split one cell of 1024 differs in the last bit (4e-16 in DS, 1.3e-15 in DDDT: a summation-order effect at a box corner).
* Pressure-dependent stages (H, HS, US, WS, U, W, and pass 2 D/DS/DDDT which read the pass-1 velocities): differ at 1e-16 (velocities) to 1e-14 (H) and in a few percent of the cells: the FFT sums in a layout dependent order (the distributed transform changes with the box split and rank count). This is the irreducible
  part without a decomposition independent (fixed-order) Poisson solve.
* Whole run (82 steps, T = 1): final U W H HS D DS RHO TMP ZZ of 4 and 16 boxes at 1 and 4 ranks equal the 1-box run to <= 5.7e-13 relative (before S6: 3e-5 absolute, 1e-3 relative at T = 1).
* The S5 results stand: `shunn3_32` (1 rank, 80 steps) T/DT 17/17 rows, MMS errors equal to 5 digits (T2 pass), fields <= 1.4e-14 absolute to the baseline restart, mass exact; `shunn3_4mesh_32` (1 and 4 ranks): T2 pass
  against `shunn3_4mesh_32__1mesh`, and now within 3e-13 of a 1-box run of the same domain (the different pressure coupling of the multi-mesh baseline remains the reason it differs from the baseline).

### Plain kernel-check modes, regenerated (S6)
The 0 to 108 last-bit differences of the plain `--ghost=full/face` modes come from the boundary-face strips of U/V/W/US/VS/WS, which a frozen snapshot cannot hold in the exact pre-boundary-step state (DENSITY restores the face values
after the boundary step). `run_kernelcheck.sh` now also runs every plain case with those strips taken from the dump after the replay (`FDSKC_STRIPDUMP=1`, the `*+strips` runs) and requires ALL tags BITWISE (ppm gate 0):
7 cases x {full, face}, 1 and 4 ranks, 0 bit differences in every tag (the gate of the plain runs without strips stays at `PLAIN_MAX_PPM=100`; the worst observed case there, `4mesh_native` with 482312 compared elements per tag, has VISC_P/C 32 (66 ppm), VFLUX_P 18 (37 ppm) and VFLUX_C 26 (54 ppm) differing elements, 108 in total over the four tags, each tag below 100 ppm). The 100 ppm allowance therefore only covers the snapshot strips, not a kernel defect: the same kernels are bitwise with the exact strips.

### csmag_32 (3-D, six periodic faces): investigation and fix (S6b)
Symptom (S6): H differed from a scratch reference dump by ~1e-2 after cycle 1 although `lap(phi) = PRHS` held to 1e-12 for the driver. Findings, in order:
1. The reference is the scratch instrumented run of the UNPATCHED FDS on `csmag_32.fds`, not an official baseline (none exists; the V&V lead has to capture it). Its input has no `&PRES FISHPAK_BC`, so FDS (FFT solver, velocity tolerance 8.8e-3) treats the six periodic vents of the single mesh as mesh-to-self interpolated boundaries: Dirichlet-coupled pressure problem, 5 pressure iterations in step 1 (pressure error ~1e-11, velocity error 0.59e-2). That H is not a periodic Poisson solution: `lap(H_ref) = PRHS` is 3.5 % off, its mean is -3e-3 (gauge shift) and the remainder (shape, up to 6.6e-3) sits next to the domain faces. The driver's level solve is the periodic FFT solve. This remaining difference is a LEGITIMATE SOLVER DIFFERENCE (the same class as the 4-mesh case): T2 only, and the plain `csmag_32` has no T2 verdict until the official baseline exists.
2. Hidden behind it were three DRIVER bugs, visible only in a fully periodic 3-D case (shunn3 is 2-D with one cell in y and periodic in x, z only, which does not exercise them). They were found by comparing against a second reference dump of the same input with `&PRES FISHPAK_BC=0,0,0` (FDS's own periodic Crayfishpak solve, `tests/cases/csmag_32_fishpak.fds`):
   * The two copies of a flow face at opposite periodic domain sides (`U(0)`/`U(IBAR)` etc.) were not averaged: FDS's `MATCH_VELOCITY` does `0.5*(a + ((b*dA1)*dA2)/(dA1*dA2))` for PERIODIC/INTERPOLATED wall faces, patch 0003 makes it return under `EXTERNAL_GHOSTS_FILLED`, and the driver did the match only for box interfaces (same face in AMReX). Now `BcStep::match_periodic_faces` does it (layout independent, `ParallelCopy`), then refills the ghosts and OMESH, with `UVW_SAVE` taken before the match.
   * MU and KRES in the edge and corner cells of the domain: FDS ends COMPUTE_VISCOSITY with clamped copies (`MU(0,:,0) = MU(1,:,1)` ...); the driver's later periodic ghost fill replaced them by periodic images. New `fds_g_mu_edges_dom` restores the clamped copies on domain edges only (interface edges stay with the ghost fill); the corner cells take the adjacent interior cell because FDS writes them from z-ghost cells that hold the mirror of the gas cell there.
   * Diagnostics added (off by default): `FDSTL_EDGES=<icyc>` (EDGE%OMEGA/TAU dump of pass 1), extra `stage()`/`stage_raw()` points `pred_match`, `c_vflux`, `c_prevflux`, `p1_pre_dens`, `c_pre`.
3. Result against the FISHPAK reference dump (`tests/run_csmag_check.sh`, gate 1e-12 relative, floor 1e-12): step 1 and 2 FVX, H, HS, DS, corrector U/W agree to <= 1.3e-14 absolute (FVX 7e-15, H 7e-16, U/W 2e-16). T/DT of both steps equal to the printed digits. Before the fix the same comparison gave FVX 4e-2, H 1e-2, DS 2e-12 (step 2). The kernel, decomposition and T2 checks are unchanged (below). No existing FDS file was touched: no patch 0006 (patch 0003's note is extended: the driver's shared-face match now also covers the periodic domain faces).
4. What T2 says: the plain `csmag_32` driver run is not bitwise to the plain reference (by design, the coupled-Dirichlet iteration is replaced by the periodic solve: H ~1e-2, U/W ~3e-3 after cycle 1), T/DT equal. No T2 error measure exists for this case (no analytical solution, no official baseline): verdict NOT YET.

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

## Output writers and measurements (S7)
- **FDS-format files.** With patch 0006 the driver calls `fds_setup(3)` at the end of every step (after setting T, DT, ICYC with `fds_hook_set_step`). That runs FDS's own `UPDATE_GLOBAL_OUTPUTS`, `EXCHANGE_GLOBAL_OUTPUTS`, `UPDATE_CONTROLS`, `DUMP_GLOBAL_OUTPUTS`, `WRITE_DIAGNOSTICS` on the box data, so `CHID_devc.csv`, `_hrr.csv`, `_mass.csv`, `_steps.csv`, `_cpu.csv` and the per-step blocks of the `.out` file are FDS's own writers, not a re-implementation. `fds_p_zero_dot` clears `Q_DOT`/`M_DOT` after `T = T + DT` as MAIN_LOOP does (without it the HRR file accumulates). Slice, boundary, particle and restart dumps are not written. `FDSTL_NO_FDS_OUTPUTS=1` switches the call off; a tree without patch 0006 runs without these files (`fds_hook_step_outputs()` returns 0).
- **Result against the baseline (`tests/run_outputs_check.sh <build> <work>`).** `shunn3_32`: `_mass.csv` byte-identical to the baseline; `_hrr.csv` same rows and times, columns within 8e-10 of the column scale (the full step is T2, not bitwise); `.out` pressure-iteration lines equal; `_steps.csv` row count equal. `shunn3_4mesh_32`: `_mass.csv` np1 = np4 bitwise, total mass constant; `_hrr.csv` np1 vs np4 within 9e-10 of scale (the file differs from the 4-mesh baseline from the second row on, because T/DT differ, T2-only). `csmag_32`: `_devc.csv` with the `KE` device (5 rows, t = 0 value equal to the scratch FDS reference; KE at T_END differs by 1e-5 relative from the scratch reference, which is not an official baseline).
- **Field dump.** `<chid>_final_<FIELD>.bin` (U V W H HS US VS WS D DS RHO TMP, ZZ1..n; float64; valid cells of the level in FDS order) and `<chid>_final_manifest.txt` (sizes, T, step, field list). Repeated runs give identical files (all three cases, 1 and 4 ranks).
- **Measurements.** `<chid>_driver_perf.csv`: ranks, boxes, grid, steps, loop time, time in `pressure_scheme` (f_pres = that over the loop time), time in `solve_poisson`, time in the FDS writers, peak and end RSS (`VmHWM`, `VmRSS`, summed and maximum over ranks), field bytes. `tests/run_perf.sh <driver> <OFF fds> <work> [reps]` repeats the three cases against the `USE_AMREX=OFF` build of the same tree with the load and free memory logged (NFR-030 idle-machine rule). Numbers and caveats: `notes/m2a-gate-report.md`, section 5.
- **Gate report.** `notes/m2a-gate-report.md`: every D-006 in/out item and gate item marked PASS, PARTIAL or NOT YET with its evidence; the items that need other people (official `csmag_32` baseline, one frozen solve vs GLMAT from the pressure role, GLMAT 4-mesh baseline of the 32 case, oneAPI) are listed.

### Debug switches (environment variables, off by default)
`FDSTL_DEBUG` (per-iteration errors), `FDSTL_PBV` (pressure backend verbose and residual probes), `FDSTL_ZONES` (print the zone sums), `FDSTL_STAGE=<icyc>` (per-stage field dumps of a cycle), `FDSTL_IFACE=0|1|2` (box-interface handling, default 1), `FDSTL_EXTGHOST=0|1` (EXTERNAL_GHOSTS_FILLED route, default 1), `FDSTL_STAGEG=<icyc>` (raw FAB dumps with ghosts), `FDSTL_WALLS=<icyc>` (print the wall arrays of pass 1), `FDSTL_EDGES=<icyc>` (write `edges_p1_b<box>.bin`: EDGE OMEGA then TAU, 5 columns each), `FDSTL_NO_FDS_OUTPUTS=1` (do not call `fds_setup(3)`), `FDSTL_PDUMP=<icyc>` (write `pdump_<icyc>_{P,C}_{rhs,phi}.bin` and `pdump_dx.txt`: PRHS and the level solution on the level's cell layout, I fastest, for a one-solve comparison against another Poisson solver).

### Patches (S5, S6)
S7: patch 0006 (`main.f90`, guarded). None in S5 or S6: no existing FDS file was edited in S5 or S6 (S6 needed no patch 0006: the interface handling is done in the driver around unmodified routines, the EXTERNAL_GHOSTS_FILLED route uses the applied 0003/0004; 0005 POINT_TO_BOX stays a DRAFT until a oneAPI validation and the alias route keeps working, it is used for the kernels) (the OFF build is untouched; the bitwise OFF check of 0001-0005 stands). oneAPI (ifx/ifort) is not available on the box: gfortran 14.2 + Open MPI 5.0.7 only.

## S8: M2a closure, regression cases, diagnostics
- **OFF check on the committed tree** (`tests/check_off_bitwise.sh`, repository HEAD with patches 0001 to 0006, no `WITH_AMREX`): `shunn3_32`, `shunn3_4mesh_32`, `shunn3_4mesh_32__glmat`, `csmag_32`, `csmag_32__fishpak_bc000` are bitwise identical to the baselines (details in `patches/0006-main-step-outputs.md`).
- **Gate items** (`tests/run_m2a_baselines.sh <driver-build> <work>`; comparisons in `tests/m2a_compare.py`): results and the T2 definition in `notes/m2a-gate-report.md`.
- **Regression cases for the periodic-only bugs** (`tests/run_periodic_regression.sh <driver-build> <work>`, 1 rank on the official single-box baseline and 4 ranks on 2x2x1 boxes from `tests/make_periodic_case.py`): `periodic_face_match`, `mu_edge_corner`, `kres_edge_corner`. Positive leg: the invariant of `tests/check_periodic_raw.py` holds in the `FDSTL_STAGEG` dump of cycle 2. Negative leg: the same run with the fix switched off (`FDSTL_SKIP_FIX=match|mu|kres`) must violate it, so the test provably sees the bug.
- **`tests/make_periodic_case.py <N> <outdir> [t_end] [sx sy sz]`**: a fully periodic 3-D LES case of N^3 cells that the driver accepts (csmag_32 recipe, Taylor-Green start, optional box split). Example: 128^3, 10 steps, 1 rank, from `/workspace/fds-amr/scratch/role1-s8-work/big/tg128` (regenerate with the script; the 56 MB csv is not committed).
- **Diagnostic environment variables (all off by default, no effect on results):** `FDSTL_PROFILE=1` writes `<chid>_driver_profile.txt` (seconds in `FillBoundary`, `fill_omesh`, the periodic face match, the FDS boundary routines, next to the loop, pressure, solve and writer times); `FDSTL_STAGES=<icyc>` dumps the Fortran-side arrays that are not registered (`FX FY FZ`, `WORK1..9`, `SWORK1..4`, `DEL_RHO_D_DEL_Z` of pass 1) and, once, the mesh metrics `RDX..RDZN R RRN` and the species tables (`MU_RSQMW_Z K_RSQMW_Z CP_Z H_SENS_Z RSQ_MW_Z MWR_Z MW`) as `scr_<tag>_<name>_b<box>.bin` (int32 header[16] = rank, lb[4], ub[4]; then float64, first index fastest), read through `fds_k_xfer` mode 3; **`scr_static_SOLID_b<box>.bin`** (S9.3) is `CELL%SOLID` of every box (SideData component 0, which `fds_cell_walls` fills from the FDS cells), same layout: header rank 3, FDS index range lb = 0,0,0 to ub = IBP1,JBP1,KBP1 (valid cells plus one ghost layer, a ghost holds SideData's fill: neighbour-box value or periodic image, 0 at a closed domain edge, not FDS's own `CELL(0,..)%SOLID`), float64 0/1, I fastest, one file per box (the other dumps of `static` are mesh 1 only); `FDSTL_SKIP_FIX=match,mu,kres` switches the named periodic fix off (fault injection for the regression test only).

## S9: per-level interface, WP2 baseline, D-052/D-053
- **Per-level interface** (`notes/level-interface.md`: what exists, what is global, what a level > 0 needs, date estimates). `Level0` is now an alias of `Level`; the stage statements of `advance()` moved verbatim into `TimeLoop::Impl::s_*` bodies run through `for_levels` (level 0 is the only bound level); `LevelRegistry` implements Role 3's `fdsrt::LevelListener` (`Source/regrid_transport/RegridInterface.H`, include path added in `CMakeLists.txt`). Unit test `levels` (1 and 4 ranks) and `selftest_fds` (cf-hook request) cover the new code.
- **D-053 zone sums**: FDS summation order is the default, `--exact-zone-sums` / `FDSTL_EXACT_ZONES=1` the option (section "Exact sums"); separate test of the option `tests/run_zone_sum_check.sh <driver-build>` (exact sums of dec4, dec2/4 ranks, dec4/4 ranks bitwise equal to the 1-box run; default sums equal the FDS-order sums; the two differ by rounding only); `ILO..KHI` lower-bound kernel arguments are deferred until a loop needs them (D-052), nothing to do here.
- **WP2 periodic 128^3 CPU baseline**: `tests/run_big_periodic_baseline.sh <driver-build> <off-build> <work> [N] [steps]`, numbers and caveats in `notes/wp2-periodic-baseline.md`.
- **S9.3 `CELL%SOLID` dump**: `scr_static_SOLID_b<box>.bin` under `FDSTL_STAGES` (section "Debug switches"), test `tests/run_solid_dump_check.sh <driver-build>` (an obstruction over four boxes).
- **S9.4 flux hooks**: design note `notes/flux-hooks-design.md` (ADV in the generated density copy, DIF through a guarded patch 0008 in `divg.f90`, three phases, empty-set and no-op tests). No code yet.
- **S9.5 thin directions of the FFT backend**: `PressureBcMap.H` (FDS pressure code to `pb::BC` per direction, the one-cell rules), mode `fds_amr case.fds --pressure-bc`, sweep `tests/run_pressure_bc_sweep.sh`, notes and per-case table `notes/fft-thin-direction-check.md`, `notes/fft-thin-direction-cases.csv` (for Role 2, D-057). Behaviour change: a Dirichlet face on a one-cell x or z direction is refused with a message (before: silently solved without its term); a fully open TWO_D box gives the ignored y direction `DD`.
- **S9.6 D-056 option B, DRAFT**: `patches/0007-mesh-fine-level-boxes.patch` (+ `.md`): fine-level mesh objects in a second module array `FINE_LEVEL(:)` indexed by level, reached through `POINT_TO_BOX`; `POINT_TO_MESH` aborts on a fine mesh number. Driver side (always built): `FDS_HOOK_FINE_GUARD` at the entry of every `fds_k_*` wrapper that takes a mesh number (aborts on a number above NMESHES unless `fds_hook_set_fine_ready(1)`, never called), mode `--fine-guard-test`, `tests/run_fine_guard_check.sh`. Draft driver file `draft/fds_fine_mesh_b.f90` and `-DFDS_AMR_FINE_B_DRAFT=ON` (needs 0007 applied; off by default; modes `--fine-b-selftest`, `--fine-b-abort N`). Not oneAPI validated, not used for physics.
- **S10.1 fine-level kernels (D-056 option B, DRAFT)**: patch 0008 routes the kernel files to `POINT_TO_BOX`; driver `fds_box_obj.f90` (`BOX_OBJ(NM)`, `FDS_FINE_B_SET_VIEW`), shadow hook `FDS_HOOK_SHADOW`, the `fds_p_*` wrappers fine-ready or `FDS_HOOK_L0_ONLY`, `draft/fds_fine_box_b.f90` (`BUILD_FINE_BOX`, `--run --fine-b-shadow|--fine-b-build`), `tests/run_fine_b_shadow_check.sh`. Details and results: `patches/0007-...md` section "Fine-level kernels" and `patches/0008-...md`. Pressure stages (RHS, residual, velocity error, `H_GHOST`) and domain-edge walls of fine boxes are open.

## Layout interface (stable for Role 3)
These headers are the interface that the hierarchy code may include. Stable means: names, meanings, index maps and ghost widths below do not change without a note in this README and a heads-up to the consumers; additions are allowed. Everything else (all `fds_*.f90`, `fds_alias.c`, the `bind(C)` names, `TimeLoop::Impl`, the scratch dump files) is internal and may change.

| Header | What it guarantees |
|---|---|
| `Fields.H` | `FieldSpec`/`field_table()`: the registry of the FDS arrays as one `MultiFab` each on the level-0 `BoxArray` (names `RHO RHOS TMP ZZ ZZS U V W US VS WS H HS KRES D DS DDDT MU RSUM FVX FVY FVZ`, see `default_field_names()`). Ghost widths (D-031): `RHO`, `RHOS` 3; `TMP`, `ZZ`, `ZZS` 2; all others 1. Staggering: `U V W US VS WS` are nodal in their own direction. `ZZ`/`ZZS` carry `N_TOTAL_SCALARS` components (component c is FDS scalar c+1). Index maps `to_amrex`/`to_fds`: valid cell FDS I = a - lo + 1; face FDS `U(I)` (the high face of cell I) is AMReX face `lo + I` (an AMReX face index is the LOW face of the cell of that index), so FDS I = a - lo in a nodal direction; FAB bounds equal the FDS allocation bounds (`fds_bounds`). `Fields::operator[]`, `has`, `spec`, `names`, `fill_ghosts` (FillBoundary incl. periodic images), `bytes`. Per-scalar scratch (`FX FY FZ`, `WORK*`, `SWORK*`, `ADV_FX`, `DIF_FX`) is NOT registered (`is_per_box_scratch`). |
| `FdsAmr.H` | `Level` (alias `Level0`: geometry, `BoxArray`, `DistributionMapping`, per-mesh info, periodicity, `dx`, **`level`, `ref_ratio_from_parent`, `fds_mesh_offset`, `fds_bound()`**) and `assemble_level0`: one box per FDS `&MESH`, box i is FDS mesh i+1, owner rank = FDS process; no re-boxing; uniform Cartesian cells only. `make_layout_level()` builds a layout-only level (no FDS mesh) from a regrid layout. |
| `LevelRegistry.H` (new, S9) | `LevelRegistry : fdsrt::LevelListener`: per level a `Level`, `Fields`, `SideData`, covered-by-finer mask; `make_level` / `remake_level` (old objects readable via `retired_*` until `release_retired`) / `clear_level`; `rebuild_side_data(level)`; `covered_mask(level)`; `fds_bound(level)` (false for every level > 0 until fine-level FDS mesh objects exist). Level 0 is adopted from `TimeLoop::registry()`. |
| `SideData.H` | `iMultiFab` with ng=2 on the level-0 layout: comp 0 SOLID, comps 1..6 face masks (zero only between two boxes of one level), comp 7 source-in-domain flag; `hash()` is independent of layout and rank count. |
| `GhostExchange.H` | `BcStep::exchange(code, predictor)` = `ghost_exchange` + the **coarse-fine ghost hook** (`CfGhostRequest`, `CfGhostHook`, S9); `ghost_exchange(F, code)` fills the ghost layers listed by `exchange_fields(code)` at the MESH_EXCHANGE positions 1, 3, 4, 5, 6; `BcStep::after_exchange` runs the boundary-condition step of FDS for the periodic and box-to-box case (incl. the periodic face match and the domain edge/corner values of `MU`, `KRES`); physical-boundary ghosts of the fill are not touched. |
| `TimeLoop.H` | `TimeLoop::advance()` is one MAIN_LOOP iteration with a global `dt` (MIN over all levels, boxes and ranks) and the same `dt` in predictor and corrector; `set_pressure_hook` replaces the pressure solve of a stage; `StepRecord`, `fields()`, `bcstep()`; **S9: `registry()`, `num_levels()`, the per-level stage entry points `stage_viscosity/density/exchange/boundary/velocity_flux/wall_bc/divergence1/divergence2/velocity_update/velocity_correct(level, ...)`, `global_dt()`, `set_cf_ghost_hook()`**; a level without an FDS binding aborts with a message. `RunOptions::exact_zone_sums` (D-053). |
| `FdsSetup.H` | `build_level0()` and `fds_cell_walls()`: layout and cell wall data from the FDS set-up. |
| `TimeStep.H` | header-only T/DT rules of FDS on plain doubles. |
| `ExactSum.H` | `exact_group_sums`, `exact_sum`, `exact_sum_product`: sums that do not depend on box split, order or rank count (D-028); **S9: `exact_sum_uncovered`, `exact_sum_product_uncovered`, `exact_sum_hierarchy` (cells covered by a finer level contribute nothing; one scale for the hierarchy)**. |
| `Source/pressure_backend/PressureIface.H` | not mine (Role 2); the driver is its only caller (`TimeLoop::solve_poisson`). |

Guarantees that hold across all of them: the numbers a kernel sees in a box are FDS's own (unmodified kernels, same operation order); with `USE_AMREX=OFF` nothing changes (`tests/check_off_bitwise.sh`); decomposition changes alter only the reductions (`tests/run_decomp_check.sh`). The multi-level generalisation (S9) is in: the level index, per-level stage functions, the global dt, the uncovered-cell sums, the coarse-fine ghost hook and the per-level `SideData` rebuild exist and are described in `notes/level-interface.md` (what a level > 0 still needs: fine-level FDS `MESH_TYPE` objects, section 3 of that note; flux read-out and override: next task). The single-level results are unchanged (byte-identical final fields in exact-zone-sum mode, `tests/run_driver_tests.sh` test `levels`).
