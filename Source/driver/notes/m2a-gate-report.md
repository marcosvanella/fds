# M2a gate report (Role 1, driver at the S7 commit)

Gate wording: D-006, `Manuals/FDS_AMReX_Design/spec-responses.md` (a), and `roadmap.md` (Demo gate M2a). Toolchain: gfortran 14.2, Open MPI 5.0.7, AMReX 26.09,
`OMP_NUM_THREADS=1`, 1 and 4 ranks. Reference: FireX `36975d765f` baseline (`vv-runs/baseline/gnu_ompi_firex-36975d7`).
Rulings applied (Architect): bitwise applies only to `USE_AMREX=OFF` and to frozen-input kernel tests. Full steps with the FFT backend are T2-only (0 of 80 STEP
records of `shunn3_32` are bitwise; the ~1e-3 difference between the 4-mesh and the 1-mesh/driver coupling is accepted as T2-only). Obstruction masks stay NotBuilt
(OBST cases are refused). T2 here means the metric of `tests/compare_run.py`: T/DT sequence to the printed digits, MMS error of the driver <= 1.05 x baseline error,
total mass, final fields against the baseline restart file.

Legend: PASS = evidence in this repository, reproducible with the named script. NOT YET = cannot be decided by Role 1 (needs another person or a missing input).
FAIL = none.

## 1. Case table (what each case must show)

| # | Gate item | Case | Status | Evidence |
|---|---|---|---|---|
| 1.1 | Bitwise explicit-kernel checks on frozen input vs baseline (mass/species, viscosity, velocity flux, divergence, predictor/corrector) | `shunn3_32` | PASS | `tests/run_kernelcheck.sh`: every wrapped kernel bitwise in the dump, full+bc and face+bc modes; plain modes within `PLAIN_MAX_PPM`, all tags bitwise with `+strips` |
| 1.2 | same | `csmag_32` | PASS (scratch reference) | same script, `csmag_32` 19 tags bitwise in face+strips. The frozen input comes from a scratch run of the UNPATCHED FDS on the plain input, not an official baseline (see 4.1) |
| 1.3 | Bitwise kernel checks vs a derived single-mesh copy (A-24) | `shunn3_4mesh_32` | PASS | `run_kernelcheck.sh` 4mesh window/clip modes at 4 ranks: all tags bitwise (16 of 16 in `face+strips`), derived single-mesh copy `shunn3_4mesh_32__1mesh` |
| 1.4 | Whole run T2 vs baseline (CSV and field dump at `T_END`) | `shunn3_32` | PASS (T2) | 17 of 17 T/DT rows equal; MMS errors at T=0.9003 equal to 5 digits (e_rho 3.3002e-2, e_Z 6.2455e-3, e_u 4.7980e-3, e_H 1.8152e-2); mass exact; final fields vs baseline restart U 1.3e-15, W 1.2e-15, H 1.8e-14, HS 2.2e-14, D 2.2e-14, DS 1.9e-14, RHO 7.1e-15, TMP 3.4e-13; pressure-iteration counts equal FDS; `_mass.csv` bitwise equal, `_hrr.csv` within 8e-10 of the column scale (`tests/run_outputs_check.sh`). Not bitwise: 0 of 80 STEP records (accepted, T2-only) |
| 1.5 | Whole run T2 | `csmag_32` | NOT YET | The plain input has no official baseline, so there is no T2 reference and no analytical error measure. Against the scratch plain reference H differs by up to 9.5e-3 (gauge 2.9e-3) and U/W by ~3e-3 after cycle 1: FDS solves a Dirichlet-coupled problem there (the six periodic vents of a single mesh are mesh-to-self interpolated boundaries; 5 iterations in step 1), the driver solves the periodic problem. Against FDS's own periodic solve of the same input (`&PRES FISHPAK_BC=0,0,0`, `tests/run_csmag_check.sh`) the driver agrees to <= 1.3e-14 absolute for FVX, H, HS, DS, U, W in steps 1 and 2, and T/DT are equal. KE device value at `T_END` differs from the scratch reference by 1.0e-5 relative. Needs the V&V lead to capture the official baseline |
| 1.6 | Full steps T2 vs the derived single-mesh copy | `shunn3_4mesh_32` | PASS (T2) | `tests/compare_run.py` vs `shunn3_4mesh_32__1mesh`: e_rho 3.2253e-2 vs 3.2255e-2, e_Z 1.3710e-2 vs 1.3709e-2, e_u 4.953e-3 vs 5.0145e-3, e_H 7.466e-3 vs 8.680e-3, all within 1.05 x; at 1 and 4 ranks. T/DT differ from step 1 (T 6.189e-4 vs 6.188e-4): accepted as T2-only |
| 1.7 | Full steps T2 vs the 4-mesh baseline (FFT) | `shunn3_4mesh_32` | PASS (T2) | same metric against `shunn3_4mesh_32`; driver mass exact (baseline drifts by 2.8e-5); final fields differ by ~1e-3 from the baseline (the FDS multi-mesh interpolated pressure coupling versus one level solve: accepted as T2-only) |
| 1.8 | Full steps T2 vs the 4-mesh baseline with `SOLVER='GLMAT'` | `shunn3_4mesh_32` | NOT YET | No GLMAT 4-mesh baseline of the 32 case exists in the baseline set (only `shunn3_4mesh_128__sf17_glmat`). Needs the V&V lead to run it |
| 1.9 | `H` from one frozen solve within eps_H vs GLMAT (FR-002) | `shunn3_4mesh_32` | NOT YET | Needs the pressure role (Role 2). Hand-over provided: `FDSTL_PDUMP=<icyc>` writes PRHS and the level solution (`pdump_<icyc>_{P,C}_{rhs,phi}.bin`, `pdump_dx.txt`) |
| 1.10 | Decomposition: 1 mesh vs 4 meshes vs 4 ranks | `shunn3_4mesh_32` | PASS | final fields 1-mesh vs 4-mesh vs np1 vs np4 within ~1e-14 (TMP 4e-13); `_mass.csv` of np1 and np4 bitwise equal; `_hrr.csv` np1 vs np4 within 9e-10 of scale |

## 2. In the demo (D-006 list)

| # | Item | Status | Evidence |
|---|---|---|---|
| 2.1 | `USE_AMREX` CMake option, gfortran 14.2 + Open MPI 5 only | PASS | patch 0002; builds with the commands in `README.md` |
| 2.2 | One-level AmrCore driver, one box per `&MESH` | PASS (with a note) | one `BoxArray` box per mesh (`FdsAmr.cpp`), `max_grid_size` = mesh size. The driver uses `amrex::Geometry`/`BoxArray`/`MultiFab` for one level; it does not instantiate the `AmrCore` class, because the demo has no refinement. The Architect may want this stated in the gate |
| 2.3 | Shim binds all fields the full step touches | PASS | field registry (`Fields.cpp`), alias and `FDS_HOOK_SET_VIEW`/`POINT_TO_BOX` (patch 0005 is a DRAFT: gfortran only) |
| 2.4 | Whole-step kernels through the shim: mass/species, viscosity and constant Smagorinsky, velocity flux, divergence, predictor/corrector, stability/dt, `UVW_FILE` init | PASS | `TimeLoop.cpp`, `fds_kernels.f90`; kernel checks (1.1 to 1.3); `csmag_32` runs with `UVW_FILE` and constant Smagorinsky |
| 2.5 | `FillBoundary` ghost fill, periodic and box-box | PASS | `GhostExchange.cpp`; decomposition check `tests/run_decomp_check.sh` (1, 4, 16 boxes at 1 and 4 ranks): 9 of 9 reduction-free stages bitwise, final fields <= 7e-13 relative. S6b added the periodic domain-face match and the domain-edge MU/KRES copies (found by `csmag_32`) |
| 2.6 | Pressure outside the shim: FFT fast path for one box | PASS | `pb::solve_pressure` (FFT backend), residual `lap(phi) = PRHS` to 1e-12 on `csmag_32` |
| 2.7 | Pressure, multi-box | PASS (FFT) | the level-wide FFT covers the multi-box periodic cases; the composite MLMG is not used (the roadmap text says multi-box is covered by the level FFT here). Homogeneous Poisson data only; mixed open/closed faces and obstruction masks return NotBuilt and OBST cases are refused (ruling) |
| 2.8 | Output: `CHID_devc.csv`, `CHID_hrr.csv`, `CHID_mass.csv`, `.out`, `_steps.csv`, `_cpu.csv` | PASS | written by FDS's own routines through `fds_setup(mode=3)` (patch 0006). `tests/run_outputs_check.sh` PASS: `shunn3_32` `_mass.csv` bitwise, `_hrr.csv` columns within 8e-10 of scale, `.out` pressure iterations equal, `_steps.csv` row count equal; `csmag_32` `_devc.csv` has the `KE` device (5 rows, t=0 equal to the reference). Without patch 0006 the driver still runs and writes only the `_driver_*` files |
| 2.9 | Field dump for comparison | PASS | `<chid>_final_<FIELD>.bin` (U V W H HS US VS WS D DS RHO TMP ZZ1..n, float64, valid cells of the level, FDS order) plus `<chid>_final_manifest.txt`; two identical runs give identical files for all three cases |
| 2.10 | Global CFL/dt selection | PASS | `TimeStep.H` replay and global MIN over boxes and ranks; T/DT sequence equal to the baseline for `shunn3_32` |
| 2.11 | IR-006: with `USE_AMREX=OFF` the three cases are T0 vs baseline | PARTIAL | patches 0001 to 0006 applied, `tests/check_off_bitwise.sh`: `shunn3_32` 1 rank 16 of 16 files and `shunn3_4mesh_32` 4 ranks 47 of 47 files bitwise. `csmag_32` has no OFF baseline (no official baseline) so its OFF check is NOT YET |
| 2.12 | Runs at 1 rank and 4 ranks, `OMP_NUM_THREADS=1` | PASS | `shunn3_32` 1, `csmag_32` 1, `shunn3_4mesh_32` 1 and 4 ranks (`tests/run_perf.sh`, `tests/run_outputs_check.sh`). 4 ranks on `csmag_32` was not run (single box) |
| 2.13 | NFR-030 measured (not required), f_pres (A-22) | measured | section 5 |
| 2.14 | Roadmap: FR-005 (i) and (iv) on the same cases | PARTIAL | (i): reduction-free stages byte-identical across box split and rank count (2.5). (iv): two identical runs per case give byte-identical field dumps at a fixed rank count and one thread. Not tested by Role 1: thread-count independence (only 1 thread built and run), three repeats, 2 and 8 ranks. Verification is by the V&V lead |
| 2.15 | Periodic only; no walls, combustion, radiation, particles | PASS | `check_scope()` in `TimeLoop.cpp` refuses cases outside the scope (OBST, reactions, radiation, particles, HVAC, other pressure solvers, ...) |

## 3. Out of the demo (must stay refused or absent)

OBST, VENT other than PERIODIC, OPEN boundaries and wall cells (refused, ruling); combustion, radiation, particles, HVAC, more than one pressure zone, restart
(refused or not written); Smokeview and VTK checks (no slice/boundary/particle files are written by the driver); refinement and regridding, OpenMP > 1 thread, oneAPI
build, `MPI_PROCESS` re-mapping (not done). The `.smv` and `.out` skeleton FDS writes at set-up and end is present, unvalidated for Smokeview.

## 4. Items that need other people

| # | Item | Owner | Status |
|---|---|---|---|
| 4.1 | Official `csmag_32` baseline (FireX run of `Turbulence/csmag_32.fds`); then T2 and the OFF bitwise check of 2.11 for that case | V&V lead | NOT YET |
| 4.2 | One frozen solve `H` within eps_H vs GLMAT on `shunn3_4mesh_32` (FR-002) | Role 2 (pressure) | NOT YET (hand-over `FDSTL_PDUMP`) |
| 4.3 | `SOLVER='GLMAT'` 4-mesh baseline of the 32 case | V&V lead | NOT YET |
| 4.4 | oneAPI validation of patches 0005 (DRAFT) and 0006, and of the whole driver (gfortran only here) | Architect | NOT YET; the demo is gfortran-only by D-006 |
| 4.5 | Acceptance of the `csmag_32` pressure difference as a legitimate solver difference (periodic solve vs coupled-Dirichlet iteration of a single mesh with six self-periodic vents) | Architect | for decision |

## 5. NFR-030 wall time, f_pres (A-22), NFR-031 memory note (measured, not required)

Method: `tests/run_perf.sh`, 3 runs per case and variant, median; driver (`fds_amr --run`) against the `USE_AMREX=OFF` `fds` built from the same source tree on the same machine, 1 thread. Load before every run
was 0.58 to 0.72 (1-minute average, below the NFR-030 limit of 1) and about 9 GB were free, with the 8-core box otherwise idle, so the runs meet the idle-machine rule; the run set is short (0.6 to 1.9 s) and 4 ranks share
8 cores. Values are indicative, not an NFR-030 verdict: the NFR cases are `openmp_test64a.fds` and an 8-rank anchor, not the M2a cases.

| Case | Ranks | Driver wall (s) | FDS OFF wall (s) | Ratio | f_pres (driver) | Poisson solve share | Peak RSS driver / FDS (MB per rank) | Ratio |
|---|---|---|---|---|---|---|---|---|
| `shunn3_32` | 1 | 1.86 | 1.52 | 1.23 | 0.59 | 0.15 | 53.6 / 37.4 | 1.43 |
| `csmag_32` | 1 | 0.60 | 0.62 | 0.97 | 0.32 | 0.13 | 102.1 / 78.6 | 1.30 |
| `shunn3_4mesh_32` | 1 | 1.29 | 1.50 | 0.86 | 0.23 | 0.04 | 56.2 / 39.9 | 1.41 |
| `shunn3_4mesh_32` | 4 | 0.97 | 0.65 | 1.49 | 0.39 | 0.05 | 50.2 / 34.3 | 1.46 |

- Wall time is the whole process (set-up, time loop, end), so set-up (0.3 to 0.4 s of FDS initialisation, from wall minus loop time) is a large part of the short runs. The time-loop-only driver time is in `<chid>_driver_perf.csv`.
- f_pres = driver time in `pressure_scheme` (RHS, boundary terms, solve, velocity-error evaluation, ghost fills, all iterations) divided by the time-loop time. The FFT solve itself is the "Poisson solve share" column. FDS's own `_cpu.csv` PRES share of the baseline is 0.09 (`shunn3_32`) and 0.02 (`shunn3_4mesh_32`) of the total time; it covers only FDS's solver call, so it does not compare directly with f_pres.
  The two definitions differ; the high f_pres is dominated by the pressure iteration loop (baroclinic term, no-flux, RHS, velocity-error checks and ghost exchanges run per iteration, and these cases take several iterations per solve).
- NFR-030 (<= 1.25 x baseline, proposed): 3 of 4 rows are within it, the 4-rank 4-mesh row (1.49) is over it. The 4-rank result is dominated by MPI start-up and 4 ranks sharing 8 cores for a 1 s run. Not an NFR-030 verdict (see above).
- NFR-031 (peak RSS <= 1.3 x baseline, proposed, Phase 2+): the driver peaks at 1.30 to 1.46 x the FDS OFF build per rank (largest descendant `maxrss` of the whole process). The driver holds AMReX, the field MultiFabs with ghost cells (7.1 MB on `csmag_32`, 0.9 to 1.1 MB on the 2-D cases) and, in addition, the FDS set-up state that the S2 alias does not replace (the driver copies the set-up values into the FABs and re-points the descriptors; whether the original arrays are released
  was not measured). The breakdown of the excess was not measured; the numbers are a note, not a pass. To be re-measured with a per-allocation count.
- Reproduce: `tests/run_perf.sh <driver-build> <OFF-fds-build> <work> 3`.

## 6. Summary

| Group | PASS | PARTIAL | NOT YET | FAIL |
|---|---|---|---|---|
| Case table (10 items) | 7 (1.1 to 1.4, 1.6, 1.7, 1.10) | 0 | 3 (1.5, 1.8, 1.9) | 0 |
| In the demo (15 items) | 12 | 2 (2.11 `csmag_32` OFF check, 2.14) | 0 | 0 |
| Measured only | 2.13 | | | |
| Other people (5 items) | 0 | 0 | 4 (4.1 to 4.4) | 4.5 for decision |

Reproduce: `tests/run_driver_tests.sh`, `run_kernelcheck.sh`, `run_decomp_check.sh`, `run_csmag_check.sh`, `run_outputs_check.sh`, `run_perf.sh`, `check_off_bitwise.sh` (commands in `README.md`).
