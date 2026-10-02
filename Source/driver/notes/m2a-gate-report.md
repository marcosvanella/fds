# M2a gate report (Role 1, driver at the S8 commit)

Gate wording: D-006, `Manuals/FDS_AMReX_Design/spec-responses.md` (a), and `roadmap.md` (Demo gate M2a). Toolchain: gfortran 14.2, Open MPI 5.0.7, AMReX 26.09,
`OMP_NUM_THREADS=1`, 1 and 4 ranks. Reference: FireX `36975d765f` baseline (`vv-runs/baseline/gnu_ompi_firex-36975d7`; the M2a directories `csmag_32*` and `shunn3_4mesh_32__glmat` come from the rebuilt Release reference binary `vv-runs/refbin/gnu_ompi_firex-36975d7/fds`).
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
| 1.2 | same | `csmag_32` | PASS | same script, `csmag_32` 19 tags bitwise in face+strips. The frozen input comes from an instrumented build of the same FDS source run on the official `csmag_32.fds` apart from `T_END` (0.02 instead of 0.67; the first cycles do not depend on it) |
| 1.3 | Bitwise kernel checks vs a derived single-mesh copy (A-24) | `shunn3_4mesh_32` | PASS | `run_kernelcheck.sh` 4mesh window/clip modes at 4 ranks: all tags bitwise (16 of 16 in `face+strips`), derived single-mesh copy `shunn3_4mesh_32__1mesh` |
| 1.4 | Whole run T2 vs baseline (CSV and field dump at `T_END`) | `shunn3_32` | PASS (T2) | 17 of 17 T/DT rows equal; MMS errors at T=0.9003 equal to 5 digits (e_rho 3.3002e-2, e_Z 6.2455e-3, e_u 4.7980e-3, e_H 1.8152e-2); mass exact; final fields vs baseline restart U 1.3e-15, W 1.2e-15, H 1.8e-14, HS 2.2e-14, D 2.2e-14, DS 1.9e-14, RHO 7.1e-15, TMP 3.4e-13; pressure-iteration counts equal FDS; `_mass.csv` bitwise equal, `_hrr.csv` within 8e-10 of the column scale (`tests/run_outputs_check.sh`). Not bitwise: 0 of 80 STEP records (accepted, T2-only) |
| 1.5 | Whole run T2 | `csmag_32` | PASS (periodic variant: strict; plain input: T2 form, see note) | Official baselines `csmag_32` and `csmag_32__fishpak_bc000` (28 steps, T_END 0.67), `tests/run_m2a_baselines.sh`. **Periodic variant** (`&PRES FISHPAK_BC=0,0,0`, the same discretisation as the driver): T and DT of all listed steps equal to the printed digits; the KE device series (29 rows) equal to the 8 printed digits (relative difference 0, final KE 7.8083125E-3), and to 1.7e-15 on the `__sf17` copy; final fields against a rerun of the reference binary (its devc is byte-identical to the baseline) U/V/W/US/VS/WS <= 4.6e-15, H 6.0e-16, HS 8.0e-16, D/DS <= 2.0e-14 (absolute, fields zero to round-off), RHO 1.6e-15, TMP 4.0e-13. **Plain input** (default `SOLVER=FFT`, which FDS runs as a coupled-Dirichlet iteration for the six self-periodic vents of one mesh): the driver differs from the baseline in final KE by 1.5e-4 relative (7.8094950E-3 baseline, driver 7.8083125E-3, the same value as the periodic-variant baseline), U/V/W by 3.4e-3 to 4.6e-3 and H by 1.5e-2 of the field maximum. T2 form (e_driver <= 1.05 e_base, `m2a_compare.py csmag-t2`) with e_base the deviation of FDS's own periodic-variant run from the plain run: ratio 1.0000 for KE and every field. This is a PASS of the stated form but it is close to tautological (the driver reproduces the periodic-variant solution to round-off), so the honest statement is: the driver is no further from the plain baseline than FDS's own periodic solver setting is. The requirements give no numerical bar for `csmag_32` (no Verification Guide metric, chaotic-case margin not calibrated), so V&V should confirm that this form is the acceptance form |
| 1.6 | Full steps T2 vs the derived single-mesh copy | `shunn3_4mesh_32` | PASS (T2) | `tests/compare_run.py` vs `shunn3_4mesh_32__1mesh`: e_rho 3.2253e-2 vs 3.2255e-2, e_Z 1.3710e-2 vs 1.3709e-2, e_u 4.953e-3 vs 5.0145e-3, e_H 7.466e-3 vs 8.680e-3, all within 1.05 x; at 1 and 4 ranks. T/DT differ from step 1 (T 6.189e-4 vs 6.188e-4): accepted as T2-only |
| 1.7 | Full steps T2 vs the 4-mesh baseline (FFT) | `shunn3_4mesh_32` | PASS (T2) | same metric against `shunn3_4mesh_32`; driver mass exact (baseline drifts by 2.8e-5); final fields differ by ~1e-3 from the baseline (the FDS multi-mesh interpolated pressure coupling versus one level solve: accepted as T2-only) |
| 1.8 | Full steps T2 vs the 4-mesh baseline with `SOLVER='GLMAT'` | `shunn3_4mesh_32` | PASS (T2) | `shunn3_4mesh_32__glmat` (4 ranks, 82 steps). The driver cannot read `SOLVER='GLMAT'` (refused: its pressure solve is the level FFT), so it runs the default-solver input and is compared with the GLMAT baseline: MMS error of mesh 1 (the only mesh that `_mms.csv` holds) e_rho 3.2253e-2 vs 3.2250e-2, e_Z 1.3710e-2 vs 1.3710e-2, e_u 4.9528e-3 vs 4.9520e-3, e_H 7.4652e-3 vs 7.4719e-3, all within 1.05 x; final fields against the four assembled restart files U/W 3.8e-5, D 6.3e-5, RHO 7.9e-6, TMP 1.5e-8 (relative to the field maximum), H/HS 1.3e-3 (a pressure field, gauge and tolerance dependent); total mass exact (GLMAT baseline drift 0). At 1 and 4 ranks identical. Informational: T and DT differ from the GLMAT baseline by up to 7.6e-4 and 4.1e-3 relative over the listed steps (T2-only, ruling). The default (FFT) 4-mesh baseline differs more (U/W 1.6e-3, see 1.7) |
| 1.9 | `H` from one frozen solve within eps_H vs GLMAT (FR-002) | `shunn3_4mesh_32` | PASS (by Role 2, other owner) | Role 2 reports the one-solve check passed: 4e-12 against eps_H = 1e-8. Not re-run by Role 1; the hand-over hook is `FDSTL_PDUMP=<icyc>` |
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
| 2.11 | IR-006: with `USE_AMREX=OFF` the three cases are T0 vs baseline | PASS | re-checked on the COMMITTED tree (repository HEAD with patches 0001 to 0006 applied, no `WITH_AMREX`, `git archive` + out-of-tree build), `tests/check_off_bitwise.sh`: `shunn3_32` 1 rank 16 of 16 files, `shunn3_4mesh_32` 4 ranks 47 of 47, `shunn3_4mesh_32__glmat` 4 ranks 47 of 47, `csmag_32` 1 rank 21 of 21, `csmag_32__fishpak_bc000` 1 rank 21 of 21 files bitwise identical to the official baselines (reference binary rebuilt with the same FireX 36975d7 source). Same result as before the merge of the newer FireX commit |
| 2.12 | Runs at 1 rank and 4 ranks, `OMP_NUM_THREADS=1` | PASS | `shunn3_32` 1, `csmag_32` 1, `shunn3_4mesh_32` 1 and 4 ranks (`tests/run_perf.sh`, `tests/run_outputs_check.sh`). 4 ranks on `csmag_32` was not run (single box) |
| 2.13 | NFR-030 measured (not required), f_pres (A-22) | measured | section 5 |
| 2.14 | Roadmap: FR-005 (i) and (iv) on the same cases | PARTIAL | (i): reduction-free stages byte-identical across box split and rank count (2.5). (iv): two identical runs per case give byte-identical field dumps at a fixed rank count and one thread. Not tested by Role 1: thread-count independence (only 1 thread built and run), three repeats, 2 and 8 ranks. Verification is by the V&V lead |
| 2.15 | Periodic only; no walls, combustion, radiation, particles | PASS | `check_scope()` in `TimeLoop.cpp` refuses cases outside the scope (OBST, reactions, radiation, particles, HVAC, other pressure solvers, ...) |
| 2.16 | Regression cases for the three periodic-only driver bugs | PASS | `tests/run_periodic_regression.sh`: `periodic_face_match` (domain-face velocity match), `mu_edge_corner` (MU in domain edge/corner cells), `kres_edge_corner` (KRES likewise). Each: positive leg (fix on: invariant holds, and on the single-box baseline the whole run matches the reference restart to the csmag gates) and negative leg (env `FDSTL_SKIP_FIX=match|mu|kres` switches the fix off: the invariant is violated, 720/392/392 cells on 1 rank, 128/408/408 on 4 ranks with 2x2x1 boxes; with the fix off the whole-run comparison also fails for match and mu, while KRES in edge cells reaches no field of this case, so its stage invariant is the test). Both at 1 and 4 ranks. |

## 3. Out of the demo (must stay refused or absent)

OBST, VENT other than PERIODIC, OPEN boundaries and wall cells (refused, ruling); combustion, radiation, particles, HVAC, more than one pressure zone, restart
(refused or not written); Smokeview and VTK checks (no slice/boundary/particle files are written by the driver); refinement and regridding, OpenMP > 1 thread, oneAPI
build, `MPI_PROCESS` re-mapping (not done). The `.smv` and `.out` skeleton FDS writes at set-up and end is present, unvalidated for Smokeview.

## 4. Items that need other people

| # | Item | Owner | Status |
|---|---|---|---|
| 4.1 | Official `csmag_32` baseline | V&V lead | DONE (used in 1.5 and 2.11) |
| 4.2 | One frozen solve `H` within eps_H vs GLMAT on `shunn3_4mesh_32` (FR-002) | Role 2 (pressure) | PASSED by Role 2 (4e-12 vs 1e-8), not re-run here |
| 4.3 | `SOLVER='GLMAT'` 4-mesh baseline of the 32 case | V&V lead | DONE (used in 1.8) |
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
- NFR-030 (<= 1.25 x baseline, proposed): 3 of 4 rows are within it, the 4-rank 4-mesh row (1.49) is over it. Profile hint (opt-in `FDSTL_PROFILE=1`, `<chid>_driver_profile.txt`, `shunn3_4mesh_32`, 4 ranks, median of 3, loop 0.79 s on a machine with load about 2, so indicative only): the extra time sits in the OMESH fill (`BcStep::fill_omesh`: every registered field of every box is broadcast to all ranks and block-copied into FDS's OMESH arrays at each boundary step and pressure iteration), 0.50 s of the 0.79 s loop; ghost fills (`FillBoundary`) 0.013 s, periodic face match 0.009 s, FDS boundary routines 0.017 s, pressure solve incl. FFT plan rebuild per call 0.04 s, FDS output writers 0.008 s. The same OMESH fill is 0.37 to 0.42 s of the 1.2 to 1.3 s loop with 1 rank (and 16 s of 58 s for a 128^3 single box), so it is the first thing to remove (broadcast only the strips a neighbour reads; patches 0003/0004 are the hook). Not an NFR-030 verdict (see above).
- NFR-031 (peak RSS <= 1.3 x baseline, proposed, Phase 2+): the driver peaks at 1.30 to 1.46 x the FDS OFF build per rank (largest descendant `maxrss` of the whole process). The driver holds AMReX, the field MultiFabs with ghost cells (7.1 MB on `csmag_32`, 0.9 to 1.1 MB on the 2-D cases) and, in addition, the FDS set-up state that the S2 alias does not replace (the driver copies the set-up values into the FABs and re-points the descriptors; whether the original arrays are released
  was not measured). The breakdown of the excess was not measured; the numbers are a note, not a pass. To be re-measured with a per-allocation count.
- Reproduce: `tests/run_perf.sh <driver-build> <OFF-fds-build> <work> 3`.

## 6. Summary

| Group | PASS | PARTIAL | NOT YET | FAIL |
|---|---|---|---|---|
| Case table (10 items) | 10 (1.1 to 1.10; 1.9 by Role 2; 1.5 plain input is a T2-form pass, see its note) | 0 | 0 | 0 |
| In the demo (16 items) | 14 (incl. 2.11 for all five baselines, 2.16) | 1 (2.14: thread-count independence and repeat runs not tested here) | 0 | 0 |
| Measured only | 2.13 | | | |
| Other people (5 items) | 4.1, 4.2, 4.3 done | 0 | 4.4 oneAPI | 4.5 for decision |

Not covered by Role 1: oneAPI (gfortran only, D-006), FR-005 beyond the checks listed in 2.5 and 2.14, the NFR-030 verdict (measured only).

Reproduce: `tests/run_driver_tests.sh`, `run_kernelcheck.sh`, `run_decomp_check.sh`, `run_csmag_check.sh`, `run_outputs_check.sh`, `run_m2a_baselines.sh`, `run_periodic_regression.sh`, `run_perf.sh`, `check_off_bitwise.sh` (commands in `README.md`).
