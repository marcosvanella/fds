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
(1) FABs sit on the host for the Fortran shim, kernels need device pointers: managed memory (spike) or device arena with staging at every host routine. (2) The driver calls whole FDS routines; a generated kernel replaces one loop inside one, so it needs a hook patch in the routine or driver stage code with host remainders (more hand-offs). (3) 34 kernels are not a step (section 1). (4) Driver and backend C++ are host-only. (5) AMReX install has no FFT or MPI; the `__int128` exact sums are unchecked on the device.

## 3. Work packages
| WP | What | Owner | Size |
|---|---|---|---|
| 0 | Kernel-level real-field harness (section 6) | AMReX Integration Lead | done |
| 1 | All non-wall cell kernels on real fields: needs missing inputs in the dump; chain the `DIVERGENCE_PART_1` nests | Integration Lead; dump: Role 1 | 1-2 d |
| 2 | Periodic real case of 128^3 or more and its CPU driver baseline | Role 1 | 1 d |
| 3 | Generator coverage for a periodic step (flux nests, predictor/corrector, `DIVERGENCE_PART_2`, `CHECK_STABILITY`, `BAROCLINIC_CORRECTION`); per-box bounds | Legacy Mapper | 5-8 d |
| 4 | GPU driver build (new AMReX CUDA install with MPI and FFT), arena choice, shim for device FABs, kernel seam | Role 1 (driver), Integration Lead (kernel calls) | 3-5 d |
| 5 | Device pressure: RHS fill and common layer as `ParallelFor`, `FFT::Poisson` on device MultiFabs | Role 2 (backend), Role 1 (call site) | 2-3 d |
| 6 | Device ghost fill and periodic face match, `dt` and zone sums as device reductions | Role 1 | 2-3 d |
| 7 | End-to-end runs and measurements (section 5) | Integration Lead | 1-2 d |

Estimate: 15-25 work-days, 2-3 calendar weeks over three owners; GPU time a few hours on the test machine, no cloud cost. Main uncertainty: WP3. WP4-7 can start with managed memory and host remainders to measure the hand-offs; the first honest end-to-end number needs part of WP3.

## 4. Acceptance criteria
- **A1 kernels, bitwise:** each generated cell kernel, device vs the same file built for the host (gfortran `-O2 -ffp-contract=off`), on real frozen fields: 0 differing elements per output array (+0/-0 counts). Against the FDS result where the dump holds it: bitwise for pure loops; about 1e-15 relative where the FDS value includes a source term (stated per kernel).
- **A2 stage chain, bitwise:** pressure frozen (`FDSTL_PDUMP` fed to both), every `FDSTL_STAGE` array of a step equal between GPU and CPU driver.
- **A3 full periodic step, tolerance** (FFT sums depend on layout): fields within 1e-12 relative of the CPU driver, T/DT rows equal to the printed digits, MMS errors within 1.05x of the CPU driver (`shunn3_32`), total mass exact (fixed-point sums).
- **A4:** GPU run repeats bitwise.

## 5. Measurements (when WP4-7 exist; run)
Per step: wall time split (kernels, FFT, `FillBoundary`, host routines, hand-offs), host/device copy count and bytes, launches; against the CPU driver (1 rank, 1 and 4 threads) on `shunn3_32` and the larger case; device memory per cell; size sweep 32^3 to 192^3 for the crossover.

## 6. First piece: generated K2 kernels on real fields, kernel level (`prototypes/s4_cuda_mass/stage1/`)
Inputs: reference-dump records of unmodified FDS (`shunn3_32`, 32x1x32, 2 species; `csmag_32`, 32^3, 1 species), second step. Fields go into one-box MultiFabs on the device arena and are passed as device pointers to 10 generated kernels (`zzs_pred, rhos_sum, rsum_pred, adv_flux_store, rho_z_p_mass, zz_corr, rho_sum, adv_flux_avg, kres, up_deardorff`) in FDS order. `SOLID` zero, `R = RRN = 1`, metrics rebuilt as in `read.f90`. Larger sizes are `csmag_32` repeated periodically 4x (128^3) and 6x (192^3): real values, not a solved flow state, kernel level only. The other kernels of the 34 need property tables or work arrays that the dumps do not hold.
- **Bitwise (run):** device vs host, 15 output arrays per case, 0 differing elements of 100848 (`shunn3_32`), 626320 (`csmag_32`), 33.5 million (128^3), 110.7 million (192^3). Host 4-thread build vs host serial: 0 differences in all four; two gfortran versions on two machines: 0 differences at 32^3.
- **Against FDS (run):** `KRES` bit-identical on all interior cells (1024 of 1024; 32768 of 32768). `RHOS` and `RHO` of `shunn3_32` agree to 2.2e-16 relative (784 and 795 of 1024 bit-identical; the manufactured source of `PERIODIC_TEST=7` cancels only to rounding). Not comparable: `ZZS` of `shunn3_32` (source term) and the density kernels of `csmag_32` (one species: FDS takes another path).
- **Kernel time (run;** median of 50 device-synchronised calls, 20 at 192^3; host 5 calls): 128^3 ten-kernel chain 4.9 ms on the device, 32.1 ms on one host thread, 27.8 ms on four (6.5x, 5.6x); per kernel 5.0x to 8.4x vs one thread (128^3), 4.5x to 7.7x (192^3). At 32^3 launch and synchronisation (about 5 us per call) limit the device: `csmag_32` 3.8x to 13.5x, `shunn3_32` (1024 cells) 0.3x to 1.6x, i.e. no gain, as the Architect expected. Tables: `results/test machine/res_gpu/timetable.md`. Not an end-to-end number: no pressure, boundary conditions or hand-offs.
- **Not done:** device FFT (blocker: new AMReX install with FFT on, needs approval), anything inside the driver, multi-box layouts, wall kernels, `S5_CALLEE_BIND` build, device repeat-run check.

## 7. Risks
Host/device hand-offs (host wall routines read and write the same FABs every stage); unified-memory page migration if managed memory is used; per-step gather of wall tables; device FFT through AMReX FFT/cuFFT (new install; even/odd boundaries and plan set-up cost on the device unverified); no-FMA and fast-math-off in every translation unit (K2 `nofma`, nvcc `--fmad=false`, `AMReX_CUDA_FASTMATH=OFF`, host `-ffp-contract=off`); K2 kernels are not idempotent and not callable from several threads on one box (OpenMP tiling needs disjoint boxes); nvfortran reassociates unparenthesised sums, so the host reference stays gfortran; `__int128` fixed-point sums on the device unchecked; kernel bounds are tied to one box per mesh.

## 8. Open questions
**Role 1 (driver)**
1. For a CUDA driver build, may the spike use `amrex.the_arena_is_managed=1` so `fds_shim_bind` and the host routines keep working beside device kernels, or do you want a device arena with staging? Where in `TimeLoop` can a kernel replace a loop of an FDS routine (hook patch like 0003/0004), or should driver stage code replace whole routines?
2. Can the reference dump (or `FDSTL_STAGEG`) also write `CELL%SOLID`, the species tables (`MU_RSQMW_Z, K_RSQMW_Z, CP_Z, H_SENS_Z, RSQ_MW_Z, MW`, `I_MAX_TEMP`), `WORK1..9`/`SWORK1..4` and the metrics (`RDX..RDZN, R, RRN`)? Is there, or can you make, a periodic 3-D case of 128^3 or more that the driver accepts (the `csmag_32` initial field is a 32^3 file), with a CPU baseline?
3. In a fully periodic step, which host routines remain (`WALL_BC` at its three positions, `MATCH_VELOCITY`, `VELOCITY_BC`, `VISCOSITY_BC`, the `OMESH` fill), and can a device-only periodic path be bitwise to the BCCHAIN test?
4. May the driver be built against a new AMReX CUDA install with MPI and FFT on (separate prefix, existing one untouched)? I will not build it without approval.

**Legacy Mapper (generator)**
1. Order, blockers and time for the loops a periodic step needs beyond the 34: `VELOCITY_FLUX` nests, `MASS_FINITE_DIFFERENCES` flux nest, `VELOCITY_PREDICTOR/CORRECTOR`, `DIVERGENCE_PART_2`, `CHECK_STABILITY`, `BAROCLINIC_CORRECTION`.
2. The dummies have fixed FDS bounds from `IBAR/JBAR/KBAR`; for several boxes per mesh (AMR, `max_grid_size`) can the generator take lower bounds or an offset as arguments?
3. Can `CHECK_STABILITY` (CFL, VN, DT_NEW) and the zone sums come out as kernels with an order-independent result (max/min, fixed-point)?
4. Who fills and refreshes the flat tables (`MU_RSQMW_Z`, `MW_SPEC`, wall tables) on the device?
