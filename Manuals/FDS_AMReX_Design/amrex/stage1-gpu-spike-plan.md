# Stage-1 GPU spike: plan and first piece

Scope: stage 1 of `docs/adr/drafts/gpu-staging-plan.md`: generated K2 cell kernels (branch `s5-gen`, HEAD `6c87023310`, `S5_CALLEE_DPD`) on real FDS fields, data resident on the GPU, device FFT pressure solve for a periodic case, end-to-end step time against the CPU driver, bitwise or tolerance check against the CPU driver. Labels: **run** = executed; **compiled** = built only; **read** = from source or documents. GPU: cc 8.9 test-machine GPU (8 GB), nvfortran 26.9, `-O2 -mp=gpu -gpu=cc89,nofma -Minline`, nvcc `--fmad=false`, AMReX CUDA build without fast math.

## 1. What exists today (survey of `Source/driver`, read-only; read)
- **Fields:** `Fields.cpp` keeps 25 FDS arrays as one `MultiFab` each on the level-0 `BoxArray` (box = FDS mesh; ghost widths per D-031). The unmodified FDS Fortran sees FAB memory through `fds_shim_bind`/`fds_alias.c` (array descriptor aimed at the FAB), i.e. through host pointers.
- **Kernels:** `fds_kernels.f90` wraps whole unmodified FDS routines (`fds_k_visc`, `fds_p_dens_pre/post` with the D-031 gather clip between, `fds_k_vflux`, `fds_k_div1/div2`, `fds_k_vpred/vcorr`); `TimeLoop.cpp` calls them per box in MAIN_LOOP order. No generated, K1 or K2 kernel is called anywhere. The wall chain (`WALL_BC`, `MATCH_VELOCITY`, `VELOCITY_BC`, `VISCOSITY_BC`) is unmodified host code on `OMESH` copies filled from the FABs (`fds_g_fill_om`).
- **Periodic case:** periodic ghosts by `FillBoundary`; the two copies of a periodic face are averaged by `BcStep::match_periodic_faces`. FDS still holds periodic faces as `INTERPOLATED` wall cells, so the host wall routines still run: no case is truly wall-free in the FDS data model.
- **Pressure:** `solve_poisson` -> `pb::solve_pressure` (Role 2, `Source/pressure_backend`) -> `FFT::Poisson<MultiFab>`. The RHS goes Fortran -> host buffer -> MultiFab by host loops. Against FDS it is T2 only (H differs at about 1e-14).
- **Device readiness:** `ParallelFor` count is 0 in `Source/driver/*.cpp`, `pressure_backend/*.cpp` and `ExactSum.H`: every `Array4` loop (RHS fill, mean removal, residual, exact sums, ghost fills) is a host loop and would fault on device FABs.
- **AMReX options:** the existing test machine CUDA install has CUDA on, **FFT off, MPI off, OpenMP off**. The FFT component has cuFFT branches but needs `AMReX_FFT=ON` (new install). Runtime knobs: `amrex.the_arena_is_managed` (all FABs unified; ADR-001 rejects it for production, fine for a spike), `amrex.use_gpu_aware_mpi`; MLMG also runs on GPU.
- **CPU baseline (Role 1 S7 report, 1 rank, 1 thread):** `shunn3_32` 1.86 s for 80 steps, f_pres 0.59; `csmag_32` 0.60 s. No real case of 128^3 or more exists.
- **Generated kernels:** 34, about 1.2% of modelled time. A periodic step also needs the `VELOCITY_FLUX` nests (4 of 5 loops "not translatable today"), the `MASS_FINITE_DIFFERENCES` flux nest ("not translatable"), `VELOCITY_PREDICTOR/CORRECTOR` (front end accepts, untested), `DIVERGENCE_PART_2`, `CHECK_STABILITY` (reduction) (`docs/inventory/gpu_generator_loop_classes.csv`).

## 2. Gap to GPU-resident stepping
(1) FABs sit on the host for the Fortran shim, kernels need device pointers: ruling (a) below, device arena with explicit staging at every host routine. (2) The driver calls whole FDS routines; a generated kernel replaces one loop inside one, so it needs a hook patch in the routine or driver stage code with host remainders (more hand-offs). (3) 34 kernels are not a step (section 1). (4) Driver and backend C++ are host-only. (5) The old test machine AMReX install had no FFT or MPI (now replaced by a separate FFT and MPI install, section 6.3); the `__int128` exact sums are unchecked on the device.

## 3. Decisions (Architect rulings, folded in)
- **(a) Memory:** field data on the device arena in the CUDA driver build; managed memory only as a debug option. K2 uses `is_device_ptr`; explicit once-per-step uploads with checksums.
- **(b) Kernel bounds (mandatory):** generated kernels take box lo/hi offsets as arguments, no fixed `IBAR/JBAR/KBAR`; single-box output stays bitwise unchanged (owner: generator).
- **(c) Reductions:** `CHECK_STABILITY` min/max are order independent as they are. Sums (zone `DSUM/PSUM/USUM`, mass): fixed-order per-box tree on the device, box partials combined in global box-index order (gather, not `MPI_Allreduce`); not bitwise against the FDS CPU path.
- **(d) Acceptance:** as in section 5, measured against the CPU run of the same driver and layout; mass drift per step at most max(1e-13 relative, 2x the CPU drift). GPU flags rule NFR-051 applies (K2 `nofma`, nvcc `--fmad=false`, fast math off).
- **D-050 (later, not now):** single global dt and interface flux overwrite (two-phase flux / override / divergence split; note `Source/regrid_transport/notes/flux-override-interface.md`); to be reviewed against this plan after stage 1.

## 4. Work packages
| WP | What | Owner | Size |
|---|---|---|---|
| 0 | Kernel-level real-field harness (section 7.1) | AMReX Integration Lead | done |
| 1 | All non-wall cell kernels on real fields: run with synthetic tables (section 7.2); real tables, work arrays and the wall kernels still need dump data | Integration Lead; dump: Role 1 | first pass done; 1-2 d left |
| 2 | Periodic real case of 128^3 or more and its CPU driver baseline | Role 1 | 1 d |
| 3 | Generator coverage for a periodic step (flux nests, predictor/corrector, `DIVERGENCE_PART_2`, `CHECK_STABILITY`, `BAROCLINIC_CORRECTION`); box lo/hi offset arguments (ruling b) | Legacy Mapper | 5-8 d |
| 4 | GPU driver build (new AMReX CUDA install with MPI and FFT), device arena (ruling a), shim for device FABs, kernel seam; the install is done (section 7.3) | Role 1 (driver), Integration Lead (kernel calls) | 3-5 d |
| 5 | Device pressure: `FFT::Poisson` on a periodic box is run (section 7.4); RHS fill and common layer as `ParallelFor`, even/odd boundaries, MultiFab layouts of the driver | Role 2 (backend), Role 1 (call site) | 2-3 d |
| 6 | Device ghost fill and periodic face match, `dt` min/max and zone sums as device reductions (ruling c) | Role 1 | 2-3 d |
| 7 | End-to-end runs and measurements (section 6) | Integration Lead | 1-2 d |

Estimate: 15-25 work-days, 2-3 calendar weeks over three owners; GPU time a few hours on the test machine, no cloud cost. Main uncertainty: WP3. WP4-7 can start with the device arena and host remainders to measure the hand-offs; the first honest end-to-end number needs part of WP3.

## 5. Acceptance criteria
- **A1 kernels, bitwise:** each generated cell kernel, device vs the same file built for the host (gfortran `-O2 -ffp-contract=off`), on real frozen fields: 0 differing elements per output array (+0/-0 counts). Against the FDS result where the dump holds it: bitwise for pure loops; about 1e-15 relative where the FDS value includes a source term (stated per kernel).
- **A2 stage chain, bitwise:** pressure frozen (`FDSTL_PDUMP` fed to both), every `FDSTL_STAGE` array of a step equal between GPU and CPU driver.
- **A3 full periodic step, tolerance** (FFT sums and the box-tree sums depend on layout), against the CPU run of the same driver and layout: fields within 1e-12 relative, T/DT rows equal to the printed digits, MMS errors within 1.05x of the CPU driver (`shunn3_32`), mass drift per step at most max(1e-13 relative, 2x the CPU drift).
- **A4:** GPU run repeats bitwise.

## 6. Measurements (when WP4-7 exist; run)
Per step: wall time split (kernels, FFT, `FillBoundary`, host routines, hand-offs), host/device copy count and bytes, launches; against the CPU driver (1 rank, 1 and 4 threads) on `shunn3_32` and the larger case; device memory per cell; size sweep 32^3 to 192^3 for the crossover.

## 7. Results so far (`prototypes/s4_cuda_mass/stage1/`)
### 7.1 WP0: ten generated K2 kernels on real fields, kernel level
Inputs: reference-dump records of unmodified FDS (`shunn3_32`, 32x1x32, 2 species; `csmag_32`, 32^3, 1 species), second step. Fields go into one-box MultiFabs on the device arena and are passed as device pointers to 10 generated kernels (`zzs_pred, rhos_sum, rsum_pred, adv_flux_store, rho_z_p_mass, zz_corr, rho_sum, adv_flux_avg, kres, up_deardorff`) in FDS order. `SOLID` zero, `R = RRN = 1`, metrics rebuilt as in `read.f90`. Larger sizes are `csmag_32` repeated periodically 4x (128^3) and 6x (192^3): real values, not a solved flow state, kernel level only. The other kernels of the 34 need property tables or work arrays that the dumps do not hold.
- **Bitwise (run):** device vs host, 15 output arrays per case, 0 differing elements of 100848 (`shunn3_32`), 626320 (`csmag_32`), 33.5 million (128^3), 110.7 million (192^3). Host 4-thread build vs host serial: 0 differences in all four; two gfortran versions on two machines: 0 differences at 32^3.
- **Against FDS (run):** `KRES` bit-identical on all interior cells (1024 of 1024; 32768 of 32768). `RHOS` and `RHO` of `shunn3_32` agree to 2.2e-16 relative (784 and 795 of 1024 bit-identical; the manufactured source of `PERIODIC_TEST=7` cancels only to rounding). Not comparable: `ZZS` of `shunn3_32` (source term) and the density kernels of `csmag_32` (one species: FDS takes another path).
- **Kernel time (run;** median of 50 device-synchronised calls, 20 at 192^3; host 5 calls): 128^3 ten-kernel chain 4.9 ms on the device, 32.1 ms on one host thread, 27.8 ms on four (6.5x, 5.6x); per kernel 5.0x to 8.4x vs one thread (128^3), 4.5x to 7.7x (192^3). At 32^3 launch and synchronisation (about 5 us per call) limit the device: `csmag_32` 3.8x to 13.5x, `shunn3_32` (1024 cells) 0.3x to 1.6x, i.e. no gain, as the Architect expected. Tables: `results/test machine/res_gpu/timetable.md`. Not an end-to-end number: no pressure, boundary conditions or hand-offs.

### 7.2 WP1: the other 13 non-wall kernels (23 of the 34; the 11 wall kernels are not run)
Same harness and protocol; four more stages on the same dump records: `cp_rhg, conductivity, kdtd, dp_kdtd, rho_d_dzd, h_rho_d_dzd, dp_species` (divergence part 1 chain), `mu_dns`, `flux_mw_fix, flux_mw_fix_zz, rho_z_p_divg` (mass), `rho_zz_clip_assign, delta_rho_zz_zero` (density clip). Real fields: `ZZ`, `TMP`, `RHO`, `RSUM`, `DEL_RHO_D_DEL_Z`, `FX/FY/FZ` (cropped to the kernel bounds). **Synthetic inputs (the dumps do not hold them):** the six property tables (`CP_Z, K_RSQMW_Z, RSQ_MW_Z, H_SENS_Z, MU_RSQMW_Z, MW_SPEC`, 5001 temperature entries, smooth positive functions), the work arrays `RHO_D, Q, QR, U_DOT_DEL_RHO_Z, DELTA_RHO_ZZ` (scaled copies of real fields), the scalars `RCON` and `RHO_ZZ_MIN`; `RDXN/RDYN/RDZN` equal `RDX/RDY/RDZ` (uniform grid); `SOLID` zero. No FDS truth exists for these stages.
- **Bitwise (run):** device vs host serial, 23 output arrays per case, 0 differing elements of 144288 (`shunn3_32`), 1154520 (`csmag_32`), 63.9 million (128^3), 212.2 million (192^3). Host 4-thread vs serial: 0 differences in all four; two gfortran versions on two machines: 0 of 1.30 million at 32^3.
- **Coverage limit (run):** `csmag_32` has one species with uniform `ZZ`, so its species-gradient outputs (`RHO_D_DZD*`, `H_RHO_D_DZD*`) are zero and `DELTA_RHO_ZZ` is zero after the zeroing kernel; only `shunn3_32` exercises them (and it is 2-D, 1024 cells).
- **Kernel time (run;** 18 calls, medians of 50 device-synchronised calls, 20 at 192^3, host 5 and 3 calls): sum over the calls, device vs host 1 thread vs host 4 threads: 128^3 7.2 ms vs 77.1 ms vs 32.4 ms (10.7x, 4.5x); 192^3 24.1 ms vs 269.3 ms vs 134.7 ms (11.2x, 5.6x); `csmag_32` 0.16 ms vs 2.10 ms vs 0.59 ms; `shunn3_32` 0.115 ms vs 0.135 ms vs 0.065 ms (1.2x, 0.6x: no gain at that size). Per kernel at 128^3: 5.6x to 30.8x against one thread (the table-lookup kernels `conductivity`, `mu_dns`, `cp_rhg` are the largest). Tables: `results/timetable_wp1.md`. Kernel level only, no hand-offs.
- **Open item:** the harness executable returns exit code 139 after printing "AMReX finalized" (all output files are complete before; seen in a short diagnostic job, cause not investigated; the exit status of the WP0 jobs was not recorded). Needs a look before the driver build (RDC and the nvfortran runtime at teardown are suspects).

### 7.3 New AMReX CUDA install (compiled)
Separate prefix on the test machine, same AMReX source commit (99ddfda) copied unmodified; the existing install is untouched. CMake options: `AMReX_GPU_BACKEND=CUDA`, `CMAKE_CUDA_ARCHITECTURES=89`, `AMReX_FFT=ON`, `AMReX_MPI=ON` (Intel MPI), `AMReX_OMP=OFF`, `AMReX_FORTRAN=OFF`, `AMReX_LINEAR_SOLVERS=OFF`, `AMReX_CUDA_FASTMATH=OFF`, `AMReX_FASTMATH=OFF`, `AMReX_GPU_RDC=ON`, `CMAKE_CUDA_FLAGS=--fmad=false`, Release, double precision, 3-D, 12 build jobs. It links `cufft`; the run needs the cuFFT library directory on `LD_LIBRARY_PATH`. Configuration and log: `results/laptop_wp1_fft/amrex_fft_install_config.txt`, `build_amrex_fft.sh/.log`. OpenMP is off: if the driver needs it, a further install is required.

### 7.4 Device FFT Poisson solve on a periodic box (run)
`fft/fft_poisson.cpp`: `FFT::Poisson<MultiFab>`, periodic in all directions, single box, `dx = 1/32`, right-hand side from exact integer and power-of-two arithmetic (bitwise identical on both builds), data on the device arena, 20 to 50 solves per size in one job. Reference: the same source built with the CPU AMReX install (FFTW) on a different host, 1 and 4 threads (times are not same-machine), and an independent spectral solve (numpy, second-order periodic Laplacian eigenvalues). The solution is defined up to a constant, so all comparisons remove the mean.
| Grid | device solve (ms) | host 1 thread (ms) | host 4 threads (ms) | device vs 1 thread | device vs 4 threads |
|---|---|---|---|---|---|
| 32x1x32 | 0.053 | 0.021 | not run | 0.4x | - |
| 128^3 | 2.37 | 74.4 | 19.5 | 31x | 8.2x |
| 192^3 | 8.87 | 299.9 | 114.5 | 34x | 12.9x |
| 256^3 | 19.7 | 777.7 | 308.4 | 39x | 15.6x |
Solver construction (plan set-up, once) 53 to 78 ms on the device; the first solve adds about 0.15 ms. Times are solves on device-resident data: no right-hand-side assembly, no host/device copies.
- **Accuracy (run):** device vs spectral reference 1.1e-15 to 3.3e-15 of max|phi| (host: 1.0e-15 to 3.2e-15); residual of the discrete Laplacian at most 6.4e-15 for a right-hand side of size 0.7. Device vs host: not bitwise (414 of 1024, 1.82 of 2.10 million and 6.38 of 7.08 million elements differ) but within 1.3e-15, 7.1e-16 and 1.0e-15 of max|phi| (32x1x32, 128^3, 192^3). Host 1 thread vs 4 threads: bitwise equal. Tolerance for A3: 1e-12 relative, far above these.
- **Repeat (run):** all solves of one job are bitwise equal on the device (51, 31, 21 and 11 solves; solution compare plus an integer checksum); a second process was not run.
- **Not covered:** even/odd (Neumann/Dirichlet) boundaries on the device, the driver's MultiFab layout, multi-box, and the cost of moving the right-hand side and `phi` to and from the host.

## 8. Risks
Host/device hand-offs (host wall routines read and write the same FABs every stage); unified-memory page migration (managed memory is a debug option only, ruling a); per-step gather of wall tables; device FFT through AMReX FFT/cuFFT (periodic box verified, section 7.4; even/odd boundaries unverified; plan set-up costs 50 to 80 ms once); no-FMA and fast-math-off in every translation unit (K2 `nofma`, nvcc `--fmad=false`, `AMReX_CUDA_FASTMATH=OFF`, host `-ffp-contract=off`); K2 kernels are not idempotent and not callable from several threads on one box (OpenMP tiling needs disjoint boxes); nvfortran reassociates unparenthesised sums, so the host reference stays gfortran; `__int128` fixed-point sums on the device unchecked (ruling c replaces them by a fixed-order box tree); until ruling b is implemented the kernel bounds are tied to one box per mesh; harness exit code 139 at teardown (section 7.2).

## 9. Open questions
Answered by the rulings: device arena (a), box lo/hi offsets in the generator (b), reductions (c), acceptance reference (d); the AMReX install (section 7.3) is built.

**Role 1 (driver)**
1. Where in `TimeLoop` can a generated kernel replace a loop of an FDS routine (hook patch like 0003/0004), or should driver stage code replace whole routines? With the device arena, which FABs are uploaded once per step and which host routines stay (checksums at each hand-off)?
2. Can the reference dump (or `FDSTL_STAGEG`) also write `CELL%SOLID`, the species tables (`MU_RSQMW_Z, K_RSQMW_Z, CP_Z, H_SENS_Z, RSQ_MW_Z, MW`, `I_MAX_TEMP`), `WORK1..9`/`SWORK1..4` and the metrics (`RDX..RDZN, R, RRN`)? WP1 ran on synthetic tables; real ones would give FDS truth for the divergence and viscosity kernels. Is there, or can you make, a periodic 3-D case of 128^3 or more that the driver accepts (the `csmag_32` initial field is a 32^3 file), with a CPU baseline?
3. In a fully periodic step, which host routines remain (`WALL_BC` at its three positions, `MATCH_VELOCITY`, `VELOCITY_BC`, `VISCOSITY_BC`, the `OMESH` fill), and can a device-only periodic path be bitwise to the BCCHAIN test?
4. Does the driver need OpenMP in the CUDA build? The new install has OpenMP off.

**Legacy Mapper (generator)**
1. Order, blockers and time for the loops a periodic step needs beyond the 34: `VELOCITY_FLUX` nests, `MASS_FINITE_DIFFERENCES` flux nest, `VELOCITY_PREDICTOR/CORRECTOR`, `DIVERGENCE_PART_2`, `CHECK_STABILITY`, `BAROCLINIC_CORRECTION`.
2. Ruling (b): lo/hi offsets as arguments; when can the harness take the first such build (single-box output must stay bitwise unchanged; WP0 and WP1 outputs are the regression set)?
3. Who fills and refreshes the flat tables (`MU_RSQMW_Z`, `MW_SPEC`, wall tables) on the device, and when?
