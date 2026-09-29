# 05. Masked single-level pressure solve: MLMG vs assembled HYPRE (A-56). PARTIAL

> **PARTIAL. Stopped early by a budget pause.** The numbers below are real measurements. The comparison is not finished. See "Not yet measured".
>
> **Update:** the section "1M-unknown CPU head-to-head" adds a single-level CPU comparison at about 1M unknowns and a recommended default backend.

## Setup

- Harness: `scratch/pressure-signoff/mlmg_mask/src/h2h.cpp`, driven by `h2h_matrix.sh`. Raw log: `h2h_matrix.log`. Pinned HYPRE v2.32.0-24 and AMReX 99ddfda, CPU only, `OMP_NUM_THREADS=1`, at most 4 MPI ranks.
- Operator: a masked 7-point operator. OPEN domain faces are Dirichlet (a diagonal term of 2/dx²). There is one pinned cell per sealed component. The mean is removed per sealed component (ruling-nonbox-level0 §5).
- **RHS: synthetic, smooth, zero mean per sealed component.** It is not FDS's own RHS, so iteration counts on real flow fields are not measured.
- Backends:
  - **M**: MLMG with the required settings (a>0 in masked cells, HYPRE bottom on level 0, per-component pins, no `setNSolve`), using the default HYPRE bottom options.
  - **Mf**: same as M, but the HYPRE bottom uses FDS's UGLMAT BoomerAMG options (PCG + BoomerAMG, coarsen 8, relax 18, 1 sweep, strong threshold 0.25, interp 6, trunc 0).
  - **H**: the whole level assembled into one HYPRE matrix and solved with PCG + BoomerAMG, using the same FDS options.
- Metrics:
  - `setup` is the operator/matrix build plus the solver or preconditioner setup. It is paid again at every regrid.
  - `re-solve` is the median time of a solve that reuses the setup (pressure is solved twice per step between regrids).
  - `true rel. L2` is the residual computed independently with one common operator.
  - `mem` is peak RSS summed over ranks (getrusage).
  - The hallways numbers are the median of 3 repeats. The stairwell numbers are a single repeat.

## Results

### hallways (`Pressure_Solver/hallways.fds`, dx = 0.0625), tolerance 1e-12

| Variant | Ranks | Iterations | Setup (s) | Re-solve (s) | True rel. L2 | Mem (MB) |
|---|---|---|---|---|---|---|
| open, M | 1 | 20 | 0.048 | 0.270 | 9.5e-13 | 119 |
| open, Mf | 1 | 20 | 0.046 | 0.276 | 9.5e-13 | 119 |
| open, H | 1 | 24 | 0.220 | 0.091 | 2.7e-12 | 76 |
| open, M | 4 | 20 | 0.015 | 0.079 | 9.5e-13 | 225 |
| open, Mf | 4 | 20 | 0.016 | 0.079 | 9.5e-13 | 225 |
| open, H | 4 | 23 | 0.092 | 0.038 | 2.6e-12 | 185 |
| sealed, M | 1 | 8 (cap) | 0.368 | 1.54 | 5.7e-13 | 360 |
| sealed, Mf | 1 | 8 (cap) | 0.480 | 0.948 | 5.7e-13 | 360 |
| sealed, H | 1 | 23 | 0.216 | 0.088 | 3.0e-12 | 76 |
| sealed, M | 4 | 8 (cap) | 0.106 | 0.484 | 5.7e-13 | 475 |
| sealed, Mf | 4 | 8 (cap) | 0.153 | 0.286 | 5.7e-13 | 475 |
| sealed, H | 4 | 24 | 0.091 | 0.037 | 3.0e-12 | 185 |

In the sealed variants MLMG hit the 8-iteration cap and reported "not converged". The true residual was nonetheless 5.7e-13. This matches D-032: MLMG's convergence flag is not trusted on singular problems.

### stairwell, mesh union 192x176x552 (18.65 M cells), tolerance 1e-10, one repeat

This is the `ba=union` form: the union of the meshes, with no mask and no sealed component (the OPEN vent at `stairwell.fds:451` makes it non-singular).

| Backend | Ranks | Iterations | Setup (s) | Re-solve (s) | True rel. L2 | Mem (MB) |
|---|---|---|---|---|---|---|
| M | 4 | 3 | 1.82 | 10.4 | 2.9e-11 | 1441 |
| Mf | 4 | 3 | 1.77 | 1.92 | 4.4e-11 | 1386 |
| H | 4 | 27 | 1.71 | 1.15 | 1.6e-10 | 982 |
| M | 1 | 3 | 2.33 | 37.3 | 2.8e-11 | 1170 |

## Preliminary observations (not a recommendation)

1. **The assembled HYPRE solve (H) was fastest to re-solve and used the least memory in every measured case.**
   - On the stairwell at 4 ranks it re-solved about 1.7 times faster than Mf, with about 30% less memory.
   - On sealed hallways it re-solved 8 to 13 times faster than Mf.
2. **Setup cost is close on the stairwell** (1.7 to 1.8 s for all three).
   - On open hallways H's setup is about 5 times slower than MLMG's, although it is still under 0.25 s.
   - On sealed hallways H's setup is faster than MLMG's, because MLMG's level-0 coarsening is limited by the pin (§5 of the ruling).
3. **In MLMG, the HYPRE bottom options matter a lot on large or sealed problems.** Using FDS's BoomerAMG options (Mf) cut the stairwell re-solve from 10.4 s to 1.9 s. If MLMG is chosen, Mf should be the configuration.
4. **H stops short of the requested tolerance on the true residual,** reaching 2.6 to 3.0e-12 for a 1e-12 request and 1.6e-10 for a 1e-10 request. HYPRE PCG's stopping test is not the same norm as the common true-residual check. A fair comparison needs either a matched true-residual stop or a tolerance rescaled for H.
5. These observations apply only to a single level. MLMG's advantage, the composite multilevel solve with coarse-fine consistency, is not exercised by a single-level test. H would need a separate composite assembly for more than one level.

## 1M-unknown CPU head-to-head

### Method

- **Machine:** an owner-provided NVIDIA test machine. CPU: Intel Core i9-13900HX, 24 physical cores (8 performance cores and 16 efficiency cores, SMT off), 62 GB RAM. It is a test machine, so the CPU throttles under sustained load (see "Timing caveats"). No GPU was used.
- **Toolchain:**
  - gcc 13.3.0, called through the Intel MPI 2021.18 wrappers `mpicc`/`mpicxx` (MPI 4.1); cmake 3.28.
  - HYPRE v2.32.0: a shallow clone of the tag, built with `configure --with-MPI --disable-fortran`, `CFLAGS=-O2`, CPU only.
  - AMReX 99ddfda: CMake Release (`-O3`), 3D, `AMReX_MPI=ON`, `AMReX_OMP=OFF`, `AMReX_GPU_BACKEND=NONE`, `AMReX_HYPRE=ON` against that HYPRE.
  - Note that HYPRE is exactly the v2.32.0 tag this time. The earlier sections used v2.32.0-24.
- **Harness:** `h2h.cpp` from `scratch/pressure-signoff/mlmg_mask/`, extended only where these cases needed it. The extended copy and the scripts are in `scratch/pressure-signoff/pressure_1M/`. The raw per-run logs stay in the run folder `pressure_1M/` on the test machine. The additions are:
  - `case=mms` and `case=mms_neumann`, with the error computed against the exact solution;
  - `refine=` for hallways;
  - `warmup=`, one untimed full pass (setup plus all solves) before the timed pass;
  - `htol=` and `hrecompute=` for HYPRE PCG's stopping test.
- **Runs:**
  - `OMP_NUM_THREADS=1`, with Intel MPI `I_MPI_PIN=1` and one rank per physical core.
  - 1 and 8 ranks use performance cores only. 16 ranks use the 8 performance cores plus 8 efficiency cores.
  - Every configuration was run 3 times. The tables show medians.
  - Each run times, after one warm-up pass: setup, the first solve, and the median of 5 re-solves that reuse the setup. The initial guess is zero for every solve.
- **Backends:**
  - **Mf:** MLMG with a HYPRE bottom solver using FDS's BoomerAMG options (`fdsb_opts.txt`), with `hypre.recompute_preconditioner=0`. This is the recommended MLMG configuration.
  - **H:** the whole level assembled into one ParCSR matrix and solved with PCG + BoomerAMG, using FDS's UGLMAT options.
  - Plain MLMG with the default bottom was only probed (see below).
- **Stopping, matched on the true residual:**
  - Every backend had to reach a **true** relative L2 residual of at most 1e-10, computed with one common operator. I tightened the stop tolerance where that was needed; I did not report results at a nominal tolerance.
  - **H:** PCG at `tol=1e-10` (two-norm) met the target in all 9 H configurations (at most 9.9e-11). PCG's own final relative residual agreed with the true residual to 3 digits, and `HYPRE_PCGSetRecomputeResidual` changed nothing, so no tightening was needed. The 1.6e-10 gap seen earlier on the stairwell did not reproduce in these cases.
  - **MLMG:** MLMG's stopping test is a max norm. On mms, `tol=1e-10` stopped after 19 iterations at a true relative L2 of 1.31e-10, so mms uses `tol=4e-11` (20 iterations, 4.2e-11). On hallways `tol=1e-10` already gave 4.2e-11.
  - Because iterations are whole V-cycles or PCG steps, both backends end somewhat below 1e-10.
- **Setup definition:**
  - H: numbering, IJ assembly and the BoomerAMG setup (`HYPRE_ParCSRPCGSetup`).
  - Mf: the operator and coefficient build, plus MLMG's lazy setup inside the first `solve()` (level preparation, and HYPRE IJ assembly and BoomerAMG setup for the bottom). The lazy part is measured as first solve minus the median re-solve.
  - So the Mf "first solve" includes that lazy setup, while the H first solve does not.
- **Memory:** peak RSS from getrusage, summed over ranks and as the maximum on one rank. It includes MPI buffers and the harness's replicated cell-class arrays.

### Cases

- **(a) mms:**
  - Unit cube, 100^3 cells, dx = 0.01, 1,000,000 unknowns.
  - Exact solution u = 0.37 + sin(πx)sin(πy)sin(πz). Dirichlet u = 0.37 on all six faces, applied as the same 2/dx² face term the harness uses for OPEN faces.
  - The RHS is the analytic Laplacian, -3π² sin sin sin.
  - `max_grid_size=20` gives 125 boxes of 20^3. The suggested size of 25 or 32 was not used: either one produces 25^3 boxes, which MLMG cannot coarsen at all.
- **(b) hallways:**
  - `Pressure_Solver/hallways.fds` geometry and mask, refined 2.5 times: 320x160x160, dx = 0.025, 1,152,000 gas unknowns.
  - It has internal masked cells, Neumann walls and the OPEN (Dirichlet) vent slab.
  - **Synthetic, smooth, zero-mean RHS** (the harness's own RHS; not an FDS RHS).
  - Two BoxArray layouts, both with `max_grid_size=32`:
    - **hallC** (`ba=cover`): 250 boxes covering the domain, 8.19M cells, most of them masked.
    - **hallD** (`ba=drop`): the 84 boxes that contain gas, 2.75M cells.
  - H solves only the 1,152,000 gas unknowns in both layouts. The layout changes only how H's rows are spread across ranks.

### Results (median of 3 runs; memory is the peak RSS sum / max per rank)

**(a) mms, 100^3, Dirichlet on all faces.** For both backends the error against the exact solution is max 8.22e-5 and L2 (RMS) 2.91e-5, identical to 4 digits. The algebraic error is negligible next to the discretisation error.

| Backend | Ranks | Iterations | Setup (s) | First solve (s) | Re-solve (s) | True rel. L2 | Mem sum / max (MB) |
|---|---|---|---|---|---|---|---|
| Mf | 1 | 20 | 0.060 | 1.01 | 0.976 | 4.2e-11 | 237 / 237 |
| Mf | 8 | 20 | 0.028 | 0.301 | 0.280 | 4.2e-11 | 695 / 92 |
| Mf | 16 | 20 | 0.016 | 0.243 | 0.231 | 4.2e-11 | 1283 / 84 |
| H | 1 | 18 | 2.17 | 0.936 | 0.936 | 4.4e-11 | 497 / 497 |
| H | 8 | 18 | 1.25 | 0.374 | 0.376 | 9.9e-11 | 1030 / 136 |
| H | 16 | 18 | 1.46 | 0.569 | 0.531 | 3.5e-11 | 1683 / 112 |

**(b) hallways refined 2.5x, 320x160x160, 1,152,000 gas unknowns, synthetic zero-mean RHS**

| Layout | Backend | Ranks | Iterations | Setup (s) | First solve (s) | Re-solve (s) | True rel. L2 | Mem sum / max (MB) |
|---|---|---|---|---|---|---|---|---|
| hallD (84 boxes) | Mf | 1 | 15 | 0.112 | 1.30 | 1.26 | 4.2e-11 | 527 / 527 |
| hallD | Mf | 8 | 15 | 0.025 | 0.457 | 0.449 | 4.2e-11 | 1110 / 146 |
| hallD | Mf | 16 | 15 | 0.052 | 0.944 | 0.933 | 4.2e-11 | 1792 / 123 |
| hallD | H | 1 | 24 | 2.50 | 1.40 | 1.40 | 9.6e-11 | 624 / 624 |
| hallD | H | 8 | 26 | 1.72 | 0.740 | 0.742 | 4.1e-11 | 1244 / 200 |
| hallD | H | 16 | 23 | 1.56 | 0.751 | 0.752 | 8.9e-11 | 1969 / 148 |
| hallC (250 boxes) | Mf | 1 | 15 | 0.321 | 3.21 | 3.10 | 4.2e-11 | 1379 / 1379 |
| hallC | Mf | 8 | 15 | 0.062 | 1.14 | 1.13 | 4.2e-11 | 2021 / 264 |
| hallC | Mf | 16 | 15 | 0.127 | 2.03 | 1.84 | 4.2e-11 | 2845 / 188 |
| hallC | H | 1 | 24 | 2.52 | 1.41 | 1.41 | 9.6e-11 | 741 / 741 |
| hallC | H | 8 | 23 | 2.80 | 0.920 | 0.919 | 9.5e-11 | 1334 / 230 |
| hallC | H | 16 | 24 | 2.71 | 1.30 | 1.27 | 5.6e-11 | 2038 / 198 |

- Every run converged, with the same iteration count on every re-solve.
- H's iteration count depends slightly on the partition, because PMIS coarsening is partition-dependent.
- On hallC, H is slower than on hallD at 8 and 16 ranks. The cover layout balances all 8.19M cells, so the gas rows are unevenly spread across ranks.

**Extra probes (8 ranks, a single run each, machine idle beforehand, so not directly comparable with the medians above):**

- **Plain MLMG with the default HYPRE bottom (M):**
  - mms re-solve 0.257 s, against 0.183 s for Mf in the same state (both 19 iterations).
  - hallC re-solve 1.04 s, against 1.00 s for Mf.
  - Mf stays the preferred MLMG configuration.
- **Pure-Neumann mms** (u = cos πx cos πy cos πz, one pinned cell, error taken after removing the mean difference; `ba=cover`, one repeat):

| Backend | Iterations | Setup (s) | First solve (s) | Re-solve (s) | True rel. L2 | Error max / L2 | Mem sum (MB) |
|---|---|---|---|---|---|---|---|
| Mf | 3 | 1.34 | 2.32 | 0.987 | 2.6e-12 | 8.22e-5 / 2.91e-5 | 1285 |
| H | 24 | 1.64 | 0.582 | 0.557 | 8.8e-11 | 8.22e-5 / 2.91e-5 | 1025 |

With the pin, MLMG built **only one MG level** (`mg_levels=1`). The pinned cell stops level-0 coarsening (ruling-nonbox-level0 §5). The HYPRE "bottom" solve is therefore the whole 1M-cell level, which is the same weakness as sealed hallways in the earlier section.

### Timing caveats

- **Throttling.** The test machine's clocks drop under sustained all-core load, and setup times are the most sensitive. Three measurements of mms at 8 ranks, all single runs:
  - Mf setup: 0.016 s idle, 0.028 s matrix median, 0.156 s straight after a heavy 16-rank run.
  - H setup: 0.48 s idle, 1.25 s matrix median, 2.05 s straight after a heavy 16-rank run.
  - First solves: Mf 0.21 / 0.30 / 0.73 s and H 0.25 / 0.37 / 0.52 s, in the same order.
  - The matrix interleaves all configurations, so the comparisons within it are made in similar machine states. The absolute numbers carry roughly a ±2x machine-state uncertainty.
  - **The ordering never changed: H setup was 8 to 90 times Mf setup in every state and configuration.** The low end, 8 times, is hallC at 1 rank.
- **Discarded first pass.** A first full pass of the matrix overlapped a CUDA compile run by another job on the same machine. It was discarded, and the tables come from a rerun on an otherwise idle machine. The discarded logs are kept in the run folder.
- **Scaling beyond 8 ranks is not meaningful on this machine.**
  - From 1 to 8 ranks, re-solves sped up 2.7 to 3.5 times for Mf and 1.5 to 2.5 times for H.
  - From 8 to 16 ranks, both backends got slower or stayed flat. The extra ranks sit on efficiency cores, which slow down under the static partition, and the two memory channels saturate.
  - H's setup barely scales at all (2.2 s → 1.3 to 1.5 s on mms).

### Recommendation: default backend for the solver-agnostic pressure interface

**Make MLMG with the FDS BoomerAMG bottom (Mf) the default, on a BoxArray with all-solid boxes dropped (the hallD form). Keep the assembled HYPRE PCG+BoomerAMG (H) as a selectable alternate backend behind the same interface.**

1. **Setup cost decides it.** The interface is rebuilt at every regrid, and Mf's setup is 0.016 to 0.11 s against 1.25 to 2.8 s for H at about 1M unknowns (0.32 s for Mf on the fully covered hallC at 1 rank). That is 8 to 90 times cheaper, and the ordering held in every machine state.
2. **Re-solve cost:**
   - On mms, Mf is faster at 8 and 16 ranks (0.28 vs 0.38 s; 0.23 vs 0.53 s) and about equal at 1 rank.
   - On hallD, Mf is faster at 1 and 8 ranks (1.26 vs 1.40 s; 0.45 vs 0.74 s). H is 20% faster at 16 ranks (0.75 vs 0.93 s). Break-even there is about 8 solves, or 4 time steps at two pressure solves per step.
   - H wins clearly only when MLMG has to smooth a fully covered, mostly masked level (hallC: 0.92 vs 1.13 s at 8 ranks, 1.41 vs 3.10 s at 1 rank). The interface should therefore never hand MLMG all-solid boxes.
3. **Memory is comparable.** Mf uses less on mms and hallD at 1 and 8 ranks (for example 527 vs 624 MB at 1 rank on hallD) and more on hallC (1379 vs 741 MB at 1 rank). Nothing here comes close to a constraint.
4. **Multilevel.** MLMG solves composite multilevel problems natively, with coarse-fine consistency. H would need a separate composite assembly (coarse-fine interface stencils, and a global numbering rebuilt at every regrid), which makes its setup-cost disadvantage worse.
5. **Known exception: sealed or pure-Neumann domains.** With a pin, MLMG cannot coarsen level 0, and H re-solves about 1.8 times faster (0.56 vs 0.99 s on the pure-Neumann mms). Until the pin handling is replaced (for example by MLMG's native singular-solve path with the per-component mean constraint), the interface should allow H for sealed components.

### GPU-readiness notes (analysis only; nothing was run on the GPU)

**MLMG (Mf) on GPU**

- **Builds:**
  - AMReX with `AMReX_GPU_BACKEND=CUDA` and `AMReX_CUDA_ARCH=89` (the test GPU is an Ada-generation test machine part with 8 GB), with `AMReX_HYPRE=ON`.
  - HYPRE built `--with-cuda --with-gpu-arch=89`, using the same CUDA toolkit and host compiler as AMReX. Add `--enable-unified-memory` if any hypre option in use lacks GPU support, or if host code touches hypre data. That is the conservative choice for AMReX's hypre interface until the device path has been checked.
  - CUDA-aware MPI is optional.
- **Harness changes:**
  - The coefficient, mask and RHS construction is written as host `BoxIterator` loops and must become `amrex::ParallelFor` kernels.
  - The replicated cell-class arrays must become distributed iMultiFabs.
  - MLMG's smoothing, restriction and interpolation already run on the device.
- **The bottom solve:**
  - The HYPRE bottom handles only the coarsest level, which is small and runs GPUs inefficiently. AMReX's native CG/BiCGStab bottom is a fallback that avoids hypre on the device.
  - In the pinned (sealed) case the "bottom" becomes the whole level (`mg_levels=1` above). The GPU performance of hypre then decides the result, which is one more reason to fix the pin.
- **Precision on this GPU.** A consumer GPU runs FP64 at 1/64 of its FP32 rate.
  - Stencil application and SpMV are bandwidth-bound, so this matters less for them.
  - AMG setup (SpGEMM) and coarse solves are more exposed.
  - At about 1M unknowns the footprint (about 0.25 to 0.5 GB per rank on CPU) fits comfortably in 8 GB.

**Assembled HYPRE (H) on GPU**

- **Build and runtime:**
  - HYPRE `--with-cuda --with-gpu-arch=89`, optionally `--enable-unified-memory`. In v2.32 unified memory is required only when a selected option lacks GPU support. It is a safety net while porting.
  - Runtime calls: `HYPRE_Initialize`, `HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE)`, `HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE)`, `HYPRE_SetSpGemmUseVendor(0)`, `HYPRE_SetUseGpuRand(1)`, and a device memory pool (hypre's own or Umpire).
- **Device-side assembly:**
  - The global numbering becomes a device scan over the unknown mask.
  - The row, column and value arrays are filled in one `ParallelFor` per box and passed to `HYPRE_IJMatrixSetValues`/`AddToValues` as device pointers in large chunks.
  - The RHS and solution vectors are loaded the same way.
  - The harness's host loop, which does one `SetValues` per box from host vectors, cannot be carried over. A composite multilevel assembly would need the same treatment for the coarse-fine rows.
- **BoomerAMG options:**
  - FDS's CPU choices are all on hypre's list of GPU-enabled options: PMIS coarsening (8), extended+i interpolation (6), l1-Jacobi relaxation (18) with 1 sweep, strong threshold 0.25.
  - A GPU setup still differs from the FDS CPU settings in the following ways:
    - relaxation order must be 0 (set it explicitly; C/F order is supported only for relaxation 7 and 18);
    - `KeepTranspose=1` and `RAP2=0`;
    - interpolation truncation (`PMaxElmts` 4 or so) instead of FDS's `trunc_factor=0`, to bound operator complexity within 8 GB;
    - usually one level of aggressive coarsening, with multipass (4 or 8) or second-stage (5, 6 or 7) interpolation;
    - the interpolation variant hypre recommends on GPU: the matrix-based extended+i form (hypre's GPU default) or interpolation 14 or 18 in place of 6;
    - relaxation choices restricted to 7, 18, 3/4/6, 11/12 or 16. Two-stage Gauss-Seidel (11/12) and Chebyshev (16) are the usual GPU alternatives to l1-Jacobi.
  - These change iteration counts, so the fairness rule (matched true residual) has to be re-established on the GPU.

## Solve-only rerun (pinned, CPU+GPU, with FFT)

FDS sets up the pressure matrix once per run and then re-solves it hundreds of thousands of times. This rerun therefore measures the **solve phase only**; setup appears only as a side note.

### Method

- **Machine and toolchain:** the same owner-provided NVIDIA test machine and CPU toolchain as the section above.
  - gcc 13.3.0 through the Intel MPI 2021.18 wrappers.
  - HYPRE v2.32.0 (`-O2`).
  - AMReX 99ddfda, rebuilt with `AMReX_FFT=ON` against FFTW 3.3.10. FFTW was built from source in the run folder, in double and single precision, because the system had no double-precision FFTW.
- **Harness changes** (`scratch/pressure-signoff/pressure_1M/src/h2h.cpp`):
  - A new **FFT** backend using `amrex::FFT::Poisson` with `Boundary::odd` on all faces. This is the cell-centred sine transform, which places a zero value on the face.
  - A solve-only timing loop:
    - setup once, then one untimed warm-up solve, then N timed solves;
    - every solve starts from x = 0, reset with a device-side copy;
    - the device is synchronized, and all ranks meet at a barrier, around each timed solve;
    - the harness reports the median and 90th percentile (p90) per solve.
  - What each backend's timed solve contains:
    - **H:** resetting x to 0 plus `HYPRE_ParCSRPCGSolve`, with b fixed. Gathering the solution back into the MultiFab is done once afterwards and is not timed.
    - **Mf and FFT:** the solve writes straight into the MultiFab.
- **Pinning:** ranks were pinned to the **performance cores only**. `/sys/devices/cpu_core/cpus` lists logical CPUs 0,2,4,…,14; the odd CPUs are offline and CPUs 16–31 are efficiency cores. The runs used `I_MPI_PIN=1`, `I_MPI_PIN_PROCESSOR_LIST=0,2,4,6,8,10,12,14` and `OMP_NUM_THREADS=1`, at 1, 4 and 8 ranks. Each configuration ran once, with N = 100 to 300 timed solves.
- **Stopping:** each backend stops at a true relative L2 residual of at most 1e-10, checked with the common operator. `tol=1e-10` was enough for both Mf and H here: Mf reached 2.6e-11 on mms, where its max-norm test now sits well below the L2 target, and H stayed at or below 9.6e-11. FFT is direct; its iteration count is reported as 1, and its true residual is about 3.5e-13.
- **Cases:**
  1. **mms, 100^3, homogeneous Dirichlet.** FFT::Poisson supports only homogeneous Dirichlet (odd), even (Neumann) and periodic boundaries, so all three backends use u = sin(πx)sin(πy)sin(πz) with u = 0 on every face (`hb=0`). The RHS is the analytic -3π²u. On this case the FFT operator and the harness operator (the 2/dx² face term) are the same discrete operator. 125 boxes of 20^3.
  2. **hallD:** hallways refined 2.5 times, open (OPEN vent), all-gas-free boxes dropped. 84 boxes, 1,152,000 unknowns, synthetic zero-mean RHS.
  3. **sealD:** the same geometry without the vent (`hallways_sealed`). One pinned cell, RHS mean removed, 84 boxes.
- **H CPU options:** the FDS UGLMAT set, as above.

### CPU results: solve only (one run per configuration; median / p90 over N timed solves)

| Case | Backend | Ranks | Iterations | Median solve (s) | p90 (s) | N | Warm-up solve (s) | True rel. L2 | Error max / L2 | Setup (s), side note |
|---|---|---|---|---|---|---|---|---|---|---|
| mms | FFT | 1 | 1 | 0.0226 | 0.0229 | 300 | 0.025 | 3.5e-13 | 8.22e-5 / 2.91e-5 | 0.003 |
| mms | FFT | 4 | 1 | 0.0086 | 0.0098 | 300 | 0.011 | 3.5e-13 | 8.22e-5 / 2.91e-5 | 0.001 |
| mms | FFT | 8 | 1 | 0.0046 | 0.0048 | 300 | 0.006 | 3.5e-13 | 8.22e-5 / 2.91e-5 | 0.001 |
| mms | Mf | 1 | 23 | 1.06 | 1.10 | 200 | 1.14 | 2.6e-11 | 8.22e-5 / 2.91e-5 | ~0.1 |
| mms | Mf | 4 | 23 | 0.402 | 0.704 | 200 | 0.366 | 2.6e-11 | 8.22e-5 / 2.91e-5 | ~0.02 |
| mms | Mf | 8 | 23 | 0.779 | 0.782 | 200 | 0.231 | 2.6e-11 | 8.22e-5 / 2.91e-5 | ~0.02 |
| mms | H | 1 | 21 | 1.05 | 1.06 | 200 | 1.06 | 6.9e-11 | 8.22e-5 / 2.91e-5 | 2.2 |
| mms | H | 4 | 22 | 0.429 | 0.466 | 200 | 0.405 | 2.9e-11 | 8.22e-5 / 2.91e-5 | 0.85 |
| mms | H | 8 | 22 | 0.436 | 0.467 | 200 | 0.646 | 3.7e-11 | 8.22e-5 / 2.91e-5 | 1.5 |
| hallD | Mf | 1 | 15 | 1.27 | 1.27 | 100 | 1.31 | 4.2e-11 | - | 0.16 |
| hallD | Mf | 4 | 15 | 0.738 | 0.771 | 100 | 0.477 | 4.2e-11 | - | ~0.03 |
| hallD | Mf | 8 | 15 | 0.853 | 0.946 | 100 | 0.374 | 4.2e-11 | - | ~0.03 |
| hallD | H | 1 | 24 | 1.37 | 1.37 | 100 | 1.37 | 9.6e-11 | - | 2.5 |
| hallD | H | 4 | 24 | 0.584 | 0.590 | 100 | 0.554 | 3.4e-11 | - | 1.0 |
| hallD | H | 8 | 26 | 0.734 | 0.780 | 100 | 0.769 | 4.1e-11 | - | 1.6 |
| sealD | Mf | 1 | 4 | 4.78 | 4.79 | 100 | 7.75 | 3.6e-12 | - | 3.1 |
| sealD | Mf | 4 | 4 | 2.26 | 2.27 | 100 | 3.40 | 3.6e-12 | - | 1.2 |
| sealD | Mf | 8 | 4 | 2.33 | 2.43 | 100 | 2.57 | 3.6e-12 | - | 0.26 |
| sealD | H | 1 | 25 | 1.43 | 1.43 | 100 | 1.43 | 8.1e-11 | - | 2.5 |
| sealD | H | 4 | 26 | 0.633 | 0.639 | 100 | 0.606 | 6.5e-11 | - | 1.0 |
| sealD | H | 8 | 27 | 0.776 | 0.808 | 100 | 0.743 | 4.3e-11 | - | 0.76 |

- For MLMG, "setup" is derived as the operator build plus the warm-up solve minus the median solve. Wherever the long series throttled, that difference goes negative, so those entries show only the operator-build part ("~").
- In sealD, MLMG again built only **one MG level** (`mg_levels=1`, the pin blocks coarsening). Each "iteration" is therefore a full HYPRE bottom solve on the whole level.

**Throttling. Pinning to performance cores did *not* remove it.**

- Pinning removed the efficiency-core stragglers seen in the earlier section.
- Long solve series at 4 and 8 ranks, however, run into the test machine's sustained package-power and thermal limit. The CPU package read about 90 °C during the runs.
- The clearest evidence is the warm-up solve against the median:
  - Mf at 8 ranks took 0.23 s for the warm-up solve on mms, but its median over 200 solves was 0.78 s (hallD: 0.37 s against 0.85 s).
  - Mf at 4 ranks on mms is bimodal (median 0.40 s, p90 0.70 s).
  - H ran straight after Mf, already throttled, and stayed flat (mms at 8 ranks: median 0.44 s, p90 0.47 s).
- The 1-rank numbers are stable: warm-up, median and p90 agree to within 5%.
- So on this machine the 8-rank medians are **sustained, power-limited** numbers. That is arguably the regime FDS would run in, but it is not a clean strong-scaling measurement. 4 ranks is the fastest CPU configuration for Mf and H on every case except FFT.

### GPU results

GPU runs were still pending when these CPU results were written (see below).

## Not yet measured

Done in the 1M-unknown section: the 1e-10 hallways comparison (now at 1.15M unknowns), first-solve times, per-rank memory, the stop matched on the true residual, the GPU-readiness notes, and a recommended default backend (single level, CPU).

Still open:

- **The stairwell with the mask and all-gap boxes dropped (`ba=drop`).** This is the case A-56 asked for, and it still has not run. The drop layout is the configuration the recommendation relies on.
- Stairwell repeats 2 and 3, the stairwell at a 1e-12 tolerance, and Mf and H on the stairwell at 1 rank.
- The HYPRE PCG gap between its recursive and true residuals seen on the stairwell (1.6e-10 for a 1e-10 request). It did not reproduce on mms or hallways with the v2.32.0 tag, and the cause was not investigated.
- Iteration counts with an RHS built the way FDS builds it (hallways here used the synthetic zero-mean RHS).
- **Composite multilevel solves.** Everything here is single level.
- Any GPU run of either backend.
- Scaling past 8 ranks on a node with uniform cores. On the test machine, ranks 9 to 16 run on efficiency cores.
- Timings on a thermally stable machine. Test machine throttling moves setup times by up to about 4 times, and by more for MLMG setups that take under 0.1 s (see "Timing caveats").
- The pure-Neumann mms at 1 and 16 ranks and with repeats (one run at 8 ranks only), plus a fix for MLMG's pin-limited level-0 coarsening on sealed domains.
- Plain MLMG with the default bottom (M) in the full matrix (probed only).
- The hallways manufactured-solution error (not applicable to the synthetic RHS) and a grid-convergence study for mms (one grid only).

## How to resume

Rerun the remaining blocks of `h2h_matrix.sh` (from the `stairU` 1-rank Mf/H runs onward, then `stairD`) on at most 4 cores. The harness asserts that `ba=union` has no sealed component. Then add the GPU notes and the recommendation.
