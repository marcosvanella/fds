# S4 findings: FDS mass kernel as K1 (C++ `ParallelFor`) and K2 (Fortran OpenMP `target`), CUDA compile + CPU T1

Spike S4 of ADR-001. The development machine has no GPU and no CUDA driver, so everything below is **compile-only for the device**
plus CPU runs of the same kernel source. Under ADR-001 v0.4 (see §1) this is exactly **S4a**, the part that gates
acceptance. **S4b** (the device run) is handed off in §8.

| item | status | run or only compiled |
|---|---|---|
| CUDA AMReX 99ddfda (26.09), sm_80, CUDA 13.3 | built and installed, no compile warnings | compiled only |
| K1 for sm_80 (nvcc 13.3, `--fmad=false`) | builds and links; 12 K1 kernels, 16-98 regs, no spills | compiled only (no device) |
| K2 for sm_80 (nvfortran 26.9 `-mp=gpu -gpu=cc80,nofma`) | builds and links against CUDA AMReX; 10 kernels, 38-110 regs, no spills, no local memory | compiled only |
| K1, CPU AMReX build, 1/2/4 ranks | 69/69 checks pass; 60/60 comparisons T1 **and T0 bitwise** | **run** |
| K2 host fallback (gfortran `-fopenmp`), 1/2/4 ranks | 69/69 checks pass; 60/60 comparisons T1 and T0 bitwise | **run** |
| FMA-contraction proxy (CPU, `-march=native -ffp-contract=fast`) | T0 lost everywhere; **T1 fails** in the two-level case (1.7e-3) | **run** |
| run on an NVIDIA GPU (S4b) | not possible here; script + references packed | not run |
| S4b, K1 on a cc 8.9 test-machine GPU | 29/29 comparisons T1 **and T0 bitwise**; clip counts 978 / 2520 / 154354 / 7970 | **run** (see §10) |
| S4b, K2 on the same GPU, after the §10 fix | 29/29 comparisons T1 **and T0 bitwise**; same clip counts; 144 K2 sync points | **run** (see §10) |

## 1. What was read and how it was interpreted

- **ADR version drift.** The task named ADR-001 v0.3.4. The file was v0.3.8 at the start and v0.4 by the end
  (`docs/adr/ADR-001-driver-architecture.md:5-6`). v0.4 splits S4 (`:22`): **S4a** = SDK installed; K1 and K2 written;
  both compiled for CUDA (D-029); K1 CPU build and K2 host fallback match the shimmed kernel at T1; effort per variant
  and the device-data mechanism the compiler accepts recorded. **S4b** = both variants at T1 on the device, plus K2
  stream ordering and sync cost (R-39); S4b does not gate acceptance. Spike row: `:451`.
- **Kernel route.** K1 = C++ `amrex::ParallelFor` (goes first). K2 = Fortran with `!$omp target teams loop` +
  `collapse`, device data via `has_device_addr` (fallback `is_device_ptr` + `c_f_pointer`), `map` only for host
  scalars, compiled with nvfortran (`:451`, D-029 `docs/README.md:67`). Both kernels are built. K1 is single-source
  (CPU and CUDA). D-043/NFR-049 (`requirements.md:487`): C++ physics is allowed only if the review picks K1, so K1
  here is a review variant, not a decision.
- **D-031 clip as implemented.** Pass order (a)-(g) with two host reductions per stage and level (`ADR-001:287-299`).
  One-species branch: one reduction (`:295`). Redundant ghost form: RHOP ng=3, RHO_ZZ ng=2, one pre-clip
  `FillBoundary`, density clip+apply on valid+1 (`:301-307`). Two-phase target form: density terms on valid+2, species
  terms on valid+1, then a gather in K,J,I order (`:335`). Coarse-side mask (`:316`). Flags only from valid,
  uncovered cells (`:319-322`). The new covered-coarse-cell acceptance check (`:41`, `:335`).
- **T1** (`requirements.md:25`): |x − x_ref| ≤ 1e-10·max(‖x_ref‖∞ over the field, 1e-30) for every kernel output
  field on frozen input, with `OMP_NUM_THREADS=1`. T0 = bitwise (`:24`). Scope: explicit kernels on frozen input
  against single-mesh FDS (D-022, `:29`). Fields checked: RHO, ZZ(1), ZZ(2), TMP after N steps (P1 raw format).
- **References** (P1 `check_clip.sh` sections 1-7, all valid for S4, see §6). (1) P1 gather outputs, which
  `check_clip.sh` shows are bitwise equal to the single-mesh FDS clip (`runs/clip/check_clip_output.txt`).
  (2) The interior of the single-mesh FDS restarts `fds_ref/restarted_sb/sb_r_1.restart` and
  `restarted_g8/g8_r_1.restart`. (3) The P1 two-level output `tl_g16_np1_L0/L1`. That output predates the
  v0.3.3/v0.3.4 coarse-side rules, so it is compared with `s4.coarse_mask=0` (§6).
  `prototypes/s4_cuda_mass/make_refs.py` copies them to `refs/` (11 MB, `SHA256SUMS`); the FDS restarts are
  converted to raw.
- **Needs a GPU:** only S4b (device T1, stream/sync timing). Nothing in S4a does.

## 2. CUDA AMReX build (step 1)

Command (`prototypes/s4_cuda_mass/amrex_cuda_cmake.sh`, run with `nice`, `-j3`; log `logs/amrex_cuda_build.log`):
```
. (local NVIDIA HPC SDK directory)/env_nvhpc.sh          # nvcc V13.3.73, g++ 14.2
cmake -S (local AMReX checkout) -B (local workspace)/amrex-build-cuda -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=(local workspace)/amrex-install-cuda -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ \
  -DCMAKE_CUDA_COMPILER=$(command -v nvcc) -DCMAKE_CUDA_HOST_COMPILER=g++ \
  -DAMReX_GPU_BACKEND=CUDA -DAMReX_CUDA_ARCH=8.0 -DAMReX_SPACEDIM=3 -DAMReX_PRECISION=DOUBLE \
  -DAMReX_MPI=OFF -DAMReX_OMP=OFF -DAMReX_FORTRAN=OFF -DAMReX_FORTRAN_INTERFACES=OFF \
  -DAMReX_LINEAR_SOLVERS=OFF -DAMReX_EB=OFF -DAMReX_PARTICLES=OFF -DAMReX_FFT=OFF -DAMReX_AMRLEVEL=OFF \
  -DAMReX_TINY_PROFILE=OFF -DAMReX_ENABLE_TESTS=OFF -DAMReX_BUILD_TUTORIALS=OFF -DAMReX_CUDA_FASTMATH=OFF \
  -DAMReX_CUDA_ERROR_CROSS_EXECUTION_SPACE_CALL=ON -DAMReX_CUDA_ERROR_CAPTURE_THIS=ON
cmake --build (local workspace)/amrex-build-cuda -j3 && cmake --install (local workspace)/amrex-build-cuda
```
- **Result:** all 91 CUDA translation units built with no compiler warnings; `libamrex_3d.a` (53 MB) installed.
  AMReX 99ddfda has no problem with CUDA 13.3, so no workaround was needed. CUDA 12.9's nvcc was not tried: it
  fails on glibc 2.41 (`docs/build/nvhpc-sdk.md`), so it is not an alternative here.
- `AMReX_CUDA_ARCH` is deprecated and maps to `CMAKE_CUDA_ARCHITECTURES` (CMake warning,
  `Tools/CMake/AMReXCUDAArchs.cmake:157`). `AMReX_BUILD_TUTORIALS` is unused. Neither matters.
- **`AMReX_CUDA_FASTMATH` defaults to ON** (`Tools/CMake/AMReXCUDAOptions.cmake:204`). It appends
  `--use_fast_math` to AMReX's CUDA flags (`AMReXParallelBackends.cmake:145-147`), which are PUBLIC and so reach every
  downstream TU (`:203-205`): approximate division and sqrt, FTZ, FMA. **It must be OFF for T0/T1.** The install
  records `AMReX_FASTMATH OFF` (`lib/cmake/AMReX/AMReXConfig.cmake:151`).
- Fortran interfaces, MPI and OMP are off: K1 and K2 need none of them. K2's Fortran is compiled in the S4 project,
  not in AMReX. On a GPU machine this gives one process per GPU, which matches NFR-043's single-GPU test machine.

## 3. K1/K2 implementation (branch `s4-cuda-mass`, `amrex/s4_mass/`)

- **`s4_mass_k1.H`** (350 lines, 279 code lines) holds the physics as `ParallelFor` bodies with FDS names and the
  FDS operation order. SCALAR_FACE_VALUE `:51-111` (CENTRAL, GODUNOV, SUPERBEE, MINMOD, CHARM; MP5 not ported).
  Kernels: RHO_Z_P `:118`, FACE_VALUES `:134`, MW_CORRECTION `:148`, DENSITY_UPDATE `:168`, SUM_SPECIES `:194`,
  POST_CLIP `:208`, CLIP_TERMS `:231` (phase 1: 7 terms plus a flag byte), CLIP_GATHER `:282` (phase 2: k-1, j-1,
  i-1, self, i+1, j+1, k+1; no atomics), DENSITY_APPLY `:303`, SPECIES_ONE `:313`, SPECIES_APPLY `:319`,
  RENORM `:330`. No globals, module state, recursion, `printf` or virtual calls. Every input is an argument (IR-007).
- **`s4_mass_k2.F90`** (473 lines, 372 code lines) is the same 10 kernels (+ `scalar_face_value`, `declare target`,
  `:19`) as BIND(C) subroutines over explicit-shape dummies in AMReX index space. Directive shape (`:87-89`):
  `!$omp target teams loop collapse(3) has_device_addr(...) map(to: MWR_Z) private(...)`. Wrappers are in
  `s4_mass_k2.H` (107 lines). No field data is mapped. The target regions are synchronous (no `nowait`).
- **`s4_driver.cpp`** (644 lines) is a static 1-2 level `AmrCore` in native layout (no shim). It reads P1's parameter
  names and writes P1's raw format.
  - Setup: host init in pinned memory, then copied, with no contraction in any build (`:33-58`). The mask is built
    natively (`:222-262`, coarse side `:253-256`).
  - Clip: `count_flags` over valid, uncovered cells (`:350`); clip pass order `:364-479`; pre-clip fill `:371-372`.
  - Test hooks: poison `:373`, covered or uncovered injection `:382-398`.
  - Other: `k2sync` `:287`; average-down `:536` (AMReX `amrex_avgdown`: adds and one multiply,
    `AMReX_MultiFabUtil_3D_C.H:364-370`).
- **Choices where the ADR is open.**
  - (a) The "host `ParallelAllReduce::Or`" of `:287-297` does not exist for arrays in AMReX 26.09. Only a scalar
    `bool&` overload exists (`AMReX_ParallelReduce.H:236`). The packed OR is `ParallelAllReduce::Max` on 0/1 ints
    (`s4_driver.cpp:414`, `:442`). That is one reduction each, so at most 2 host syncs per stage and level, as
    required (`:297`).
  - (b) The pre-clip fill is a same-level `FillBoundary` on every level, as in P1. Coarse/fine ghosts are never clip
    sources and lie behind wall faces.
  - (c) R = RRN = 1 (uniform Cartesian grid, D-030). PBAR is a scalar (N_ZONE = 0, no gravity), as in P1.
  - (d) `s4.coarse_mask=0` reproduces A-37/P1 for reference comparison. The default is 1 (v0.3.3/v0.3.4 rules).

## 4. Device compile results (step 2)

K1: `prototypes/s4_cuda_mass/build_cuda.sh` (cmake against `(local workspace)/amrex-install-cuda`; CMake adds
`--fmad=false -Xptxas=-v -Xcompiler=-ffp-contract=off`, `CMakeLists.txt:28-34`). The compile line (log
`logs/build_cuda.log`) is `nvcc -ccbin=g++ -O3 -std=c++20 --generate-code=arch=compute_80,code=[compute_80,sm_80]
--fmad=false -Xptxas=-v ... -maxrregcount=255 ... -rdc=true`. `cuobjdump --list-elf` gives `s4_mass.1.sm_80.cubin`.
Per-kernel figures (`logs/ptxas_table_fmad_false.txt`, from `ptxas_table.py`; SASS scan of
`cuobjdump -sass`):

| K1 kernel | regs | stack B | spill st/ld | K2 kernel (`cuobjdump -res-usage`) | regs | stack/local |
|---|---|---|---|---|---|---|
| RHO_Z_P | 56 | 0 | 0/0 | rho_z_p | 94 | 0/0 |
| FACE_VALUES | 50 | 0 | 0/0 | face_values | 64 | 0/0 |
| MW_CORRECTION | 48 | 0 | 0/0 | mw_correction | 60 | 0/0 |
| DENSITY_UPDATE | 98 | 0 | 0/0 | density_update (+sum) | 108 | 0/0 |
| SUM_SPECIES | 52 | 0 | 0/0 | (fused above) | | |
| POST_CLIP | 72 | 0 | 0/0 | post_clip | 110 | 0/0 |
| CLIP_TERMS | 54 | 0 | 0/0 | clip_terms | 84 | 0/0 |
| CLIP_GATHER | 36 | 0 | 0/0 | clip_gather | 38 | 0/0 |
| CLIP_DENSITY_APPLY | 16 | 0 | 0/0 | clip_density_apply | 38 | 0/0 |
| CLIP_SPECIES_ONE | 16 | 0 | 0/0 | (NS=1 only) | | |
| CLIP_SPECIES_APPLY | 20 | **16** | 0/0 | clip_species_apply | 42 | 0/0 |
| CLIP_RENORM | 52 | 0 | 0/0 | clip_renorm | 72 | 0/0 |

- **K1 local memory.** No spills anywhere (58 entries, max 98 registers). The only local memory is the 16 B stack
  frame of CLIP_SPECIES_APPLY: 2 `STL`, no `LDL`. It comes from `std::min(RHOP, std::max(RHO_ZZ_MIN, ...))`
  returning a reference to a local (`s4_mass_k1.H:324-325`). This is harmless; writing it with ternaries would
  remove it.
- **K1 FMA evidence.** With `--fmad=false`, the PTX has 0 `fma.rn.f64`, 116 `mul.rn.f64` and 132 `add.rn.f64`.
  With `--fmad=true` (`build-cuda-fmad`), it has 31 `fma.rn.f64` plus 248 contractible `mul/add.f64`. The SASS DFMAs
  left under `--fmad=false` (e.g. CLIP_TERMS 56 = 8 × 7) are the IEEE `div.rn.f64` Newton sequences (1 MUFU.RCP64H
  plus 7 DFMA per division), which are correctly rounded. DENSITY_UPDATE, SUM, GATHER and the apply kernels have 0
  DFMA. With `--fmad=true`, DENSITY_UPDATE gains 27, RHO_Z_P 29 and POST_CLIP 29 contractions.
- **K2 compile** (`logs/nvfortran_k2_minfo.log`): all 10 kernels are "Generating NVIDIA GPU code"; every inner
  species loop is "run sequentially", so the summation order is kept. The PTX with `nofma` has 0 `fma.rn.f64`
  (45 without it). Two findings:
  - (i) **nvfortran 26.9 rejects `has_device_addr`** (`NVFORTRAN-S-0034-Syntax error at or near identifier
    has_device_addr`). The ADR fallback is used under `#if defined(__NVCOMPILER)`: `is_device_ptr` on the
    explicit-shape dummies, which nvfortran accepts (OpenMP 5.1 treats non-`c_ptr` items as `has_device_addr`).
    `c_f_pointer` *inside* a target region fails (`NVFORTRAN-S-1058 ... pgf90_c_f_ptr_i8`). On the host it is fine,
    but it is not needed because BIND(C) array dummies already receive the device address. gfortran 14 accepts both
    clauses.
  - (ii) Without `-Minline`, `face_values` is mapped to teams only (blockidx, 1 thread per team) because of the
    `scalar_face_value` call. With `-Minline` it gets teams+threads(128). This affects performance only.
- **K2 link** (`logs/build_cuda_k2.log`): it needs one device link done by the NVHPC driver.
  - Set `CUDA_RESOLVE_DEVICE_SYMBOLS OFF`, use Fortran as the linker language, link with
    `-mp=gpu -gpu=cc80 -Mnomain -cuda -c++libs`, and add `CUDA::cudart` (`CMakeLists.txt:53-64`).
  - With a separate nvcc device link, `__fatbinwrap_*` or `__cudaRegisterLinkedBinary_*` are undefined.
  - Keep CMake's NVHPC default `-fast` off (it implies host FMA and flush-to-zero; forced to `-O2`, `:53`).
  - The linked binary has one sm_80 cubin with the nvcc and `nvkernel_s4k2_*` kernels.
- **Device-compile blockers:** none in K1 or K2. The P1 shim, the Fortran module state and `clip_gather.f90` are
  not used. Running a binary here aborts in `AMReX_GpuDevice.cpp:270` ("CUDA driver version is insufficient"),
  which is expected with no driver.

## 5. CPU vs GPU code-path differences and required flags

- **FMA contraction** is the only difference found that changes results.
  - Proxy (CPU, kernels contracted, host init not: `S4_CPU_FMA_PROXY`, `logs/run_cpu_t1_fmaproxy.txt`): T0 is lost
    in every case. Single level stays within T1 (max 1.4e-14), but clip decisions change: rhomax clips 2591 cells vs
    2520, rhomin 140911 vs 154354.
  - **Two-level blob case fails T1:** level 0 T1 error 1.67e-3 (ZZ(2); RHO 1.1e-4), level 1 1.55e-3. The SUPERBEE
    limiter branches and clip thresholds turn ulp changes into O(1e-3) differences within 12 steps.
  - **The GPU run must use `--fmad=false` (nvcc) and `-gpu=nofma` (nvfortran) for T1, not only for T0.** It must also
    use `AMReX_CUDA_FASTMATH=OFF`, and no `-use_fast_math`. The host side is g++ with `-ffp-contract=off`.
- **Math functions.** `amrex::min/max` and `amrex::Math::abs/copysign` are `std::` on both paths
  (`AMReX_Algorithm.H:25-66`, `AMReX_Math.H:48`). The kernels use only + − × ÷, min/max, abs and copysign. Division
  is `div.rn.f64`, correctly rounded, like the CPU. No `exp`/`pow`/`sin` runs on the device; those run only in the
  host initialisation.
- **Reductions.** Results depend only on integer flag counts (`ParReduce` sums of 0/1, order-independent) and the
  0/1 `Max` reduction. `MultiFab::sum` feeds diagnostic prints only: the clip mass lines can differ in the last digits
  between backends while the fields are bitwise equal. Example: the tracer rel. change prints −1.738e-15 in S4 vs
  −1.521e-15 in P1 on bitwise-identical output.
- **Other.** No atomics (gather form). No device `printf`. Tiling is ignored on the GPU. The `-rdc` default of AMReX
  (`AMReX_GPU_RDC`) does not affect arithmetic.

## 6. CPU-side T1 results (step 3), all **run**

Commands: `prototypes/s4_cuda_mass/build_cpu.sh` (AMReX `(local AMReX install)`, P1's cmake flags,
`-ffp-contract=off`); `./run_cpu_t1.sh` (K1); `./build_cpu.sh $PWD/build-cpu-k2 -DS4_K2=ON &&
./run_cpu_t1.sh $PWD/build-cpu-k2 s4.kernel=k2`. `OMP_NUM_THREADS=1`, `mpirun -np 1|2|4`.
Outputs: `logs/run_cpu_t1_k1.txt`, `logs/run_cpu_t1_k2.txt`. Case definitions: `t1_cases.sh`. Every comparison
prints the per-field T1 error and T0 status. Every one gave **max T1 error 0 and byte-identical files** for both K1
and K2:

| section (P1 `check_clip.sh` origin) | configurations (max_grid_size / ranks) | reference | result |
|---|---|---|---|
| 1 SUPERBEE, FDS initial field, 8 steps (§1) | 32/1, 8/1, 16/2, 8/2, 16/4, 8/4, fill2 8/2 | P1 = single-mesh FDS clip | 7/7 T0; 978 species-clipped cells = P1 |
| 2 FDS restart, FDS dt sequence (§2) | 32/1, 8/1, 8/2, 8/4 | FDS `sb_r_1.restart` | 4/4 T0 against FDS itself |
| 3 rhomax=1.85, rhomin=1.333 (§4) | 32/1, 8/1, 8/2, 8/4, fill2 8/4 | P1 | 10/10 T0; 2520 and 154354 density-clipped cells = P1 |
| 4 poisoned ghosts (§5) | GODUNOV, SUPERBEE, rhomax; redundant and fill2; 8/2 | FDS `g8_r_1`, `sb_r_1`; P1 | 5/5 T0; no flag set from ghosts |
| 5 two-level, `coarse_mask=0` (§7) | 16/1, 8/1, 8/2, 8/4, fill2 8/2; L0 and L1 | P1 `tl_g16_np1` | 10/10 T0 |
| 6 two-level, v0.3.3/v0.3.4 rules | 8/1, 8/2, 16/4, 8/4, fill2 8/2 vs 16/1 | own 16/1 | 10/10 T0; 7970 level-1 clipped cells |
| 7 covered-coarse-cell-only clip (new, ADR `:41`) | 16/1, 8/4 | un-injected run | flag not set; L0+L1 T0 |
| 7b coarse-side mask (ADR `:316`) | 16/1, 8/4 | own | uncovered vs level clip-mass change 3e-17 (1e-12 tol) |
| 8 IR-007 tile check | `s4.tiling=1`, 2 and 3 threads | P1 / own | 8/8 T0 |

Notes:
- **Section 6.** In this case level 0 never clips, so the v0.3.4 rules give output identical to P1 (T0).
- **Section 7.** The tracer RHO_ZZ is set to −1e-3 in covered cell (16,16,16) at step 3. The control with P1 rules
  sets the flag and changes level 0 (T1 FAIL), which shows the check is sensitive.
- **Section 7b.** The injection is at the uncovered cell (7,16,16). With P1 rules, 5.1e-9 (1/6 of the exchange) is
  taken from the covered cell and lost at average-down.
- **Totals.** 60 comparisons and 69 checks for K1 and for K2. They cover every gather-clip case of P1
  `check_clip.sh` §1-7, plus the new covered-cell and coarse-side checks.
- **Not applicable:** P1 `check_clip.sh` §8 (regressions of the `p1.clip=fds` shim) does not apply to S4.
- **K2 sync points:** 144 over 8 steps on one level, i.e. **9 per stage and level** (+1 when the density apply
  runs). They are placed around each K2 group (`s4_driver.cpp:287`). The minimum is 1 before each K2 group that
  follows AMReX-stream work, because a target region without `nowait` returns completed.

## 7. Effort per variant (S4a record, Phase 11 estimate input)

- **K1:** 279 code lines for 12 kernels. It is a line-by-line port (FDS names, `Array4` indexing). The work was in
  the D-031 two-phase clip and the index map.
- **K2:** 372 code lines (with the duplicated directives) plus 96 lines of C++ wrappers for the same kernels.
  Porting K1 to K2 was mechanical. The extra work was bounds plumbing (every array passes pointer + lo/hi + ncomp),
  `private` lists (needed for correctness under `loop`; D-029 mentions only `map`/`collapse`), the nvfortran
  `has_device_addr` gap, `-Minline`, and the mixed nvcc/nvfortran link (§4).
- Both variants passed every CPU check on the first full run after the harness fixes. No kernel bug was found by
  T1.

## 8. Handoff: what is left for the NVIDIA test machine (S4b)

- **Minimum hardware and software.**
  - One NVIDIA GPU with compute capability ≥ 8.0. sm_86/89 mobile parts run the sm_80 cubin. For cc ≥ 9.0, build with
    `CUDA_ARCH` set; PTX JIT is also safe because the PTX carries explicit `.rn` operations.
  - Driver ≥ 580.65.06 (CUDA 13.x).
  - ≥ 2 GB free GPU memory: the arena is capped at 1 GiB via `amrex.the_arena_init_size=1073741824`; the cases use
    under 100 MB.
  - NVHPC 26.9 (nvcc 13.3; nvfortran for K2), g++ ≤ 14, CMake ≥ 3.24, python3 + numpy, an AMReX clone at
    99ddfda6722b. FP64 speed does not matter.
- **Package:** `prototypes/s4_cuda_mass/s4_gpu_handoff.tar.gz` (3.2 MB, from `make_handoff.sh`) contains the
  source, the scripts, and `refs/` with `SHA256SUMS`.
- **Command:** `tar xzf s4_gpu_handoff.tar.gz && cd s4_gpu_handoff/s4_cuda_mass &&
  AMREX_SRC=<amrex@99ddfda> AMREX_PFX=<prefix> ./run_gpu_t1.sh --k2 [--fma-control]`. It runs these steps:
  - (0) Preflight: driver, compute capability, nvcc, numpy, reference checksums.
  - (1) Builds AMReX CUDA with FASTMATH OFF (the §2 command, with `CUDA_ARCH` from nvidia-smi).
  - (2) Builds K1 (`--fmad=false`), prints the ELF list and ptxas table, and runs sections 1-7b at 1 GPU with
    max_grid_size 32/16/8, plus a run-to-run repeat.
  - (3) Builds and runs K2 (`-gpu=ccXX,nofma -Minline`, `is_device_ptr`) and prints its sync count.
  - (4) Optional: the `--fmad=true` control, informational only.
  - The harness was exercised here with `--cpu-selftest=<CPU binary>`: 29 comparisons, all PASS.
- **Expected and pass/fail.**
  - **Pass (S4b):** every "T1 PASS" line and exit code 0, for K1 and K2.
  - **Expected in addition:** T0 YES in every comparison (as on the CPU); clip counts 978, 2520, 154354 and 7970;
    the run-to-run repeat bitwise equal.
  - **If T1 passes but T0 does not:** report it, and check the compile lines for `--fmad=false` / `nofma` and for
    FASTMATH.
  - **If K2 fails to start or reads garbage:** `is_device_ptr` did not carry AMReX's device pointers. This is the
    ADR's "Overturns if" (`:451`) and drops K2 (R-39).
- **Still to measure on the device (R-39):**
  - K2 ordering against AMReX's stream. The current code syncs; test removing the post-K2 syncs.
  - The per-step sync cost (time with and without `k2sync`).
  - Kernel times for K1 vs K2 (`timings` line). An `nsys` profile is optional.

## 9. Blockers, ambiguities, repository state

- **Branch name.** The task and ADR say branch `AMReX`. this repository is on **`FDS-AMReX`** at
  `1c8f150b1b` (docs commit above `36975d765f`). That worktree was not changed.
- **Worktree metadata.** `git worktree add` stores the new branch and worktree metadata in the main git dir
  `(local FDS master checkout)` (the repo's common dir); this cannot be avoided with the given command. Nothing was pushed.
- **ADR v0.3.4 → v0.4 drift** (§1). **D-029 wording** allows only `map`/`collapse` clauses; K2 also needs
  `private` and the device-address clause. **ADR `:287`** names `ParallelAllReduce::Or`, which AMReX provides
  only for a scalar (§3).
- **P1's two-level reference predates the v0.3.3/v0.3.4 rules.** It is reproduced with `s4.coarse_mask=0`; the new
  rules are checked by sections 6, 7 and 7b.
- **Not ported:** MP5, PBAR columns (N_ZONE > 0), R ≠ 1. The optional `ifx` K2 compile was not tried.

## 10. S4b: K2 device failure, root cause and fix

**First device run of K2** (before the fix), **run**: every valid cell clipped in every stage (species clip count
1048576 vs 978), total rho mass +306 %, tracer mass x11.8, max T1 error 2.65e2. K1 passed 29/29 T0 on the same GPU.
Box-split and repeat checks passed, so the result was deterministic but wrong.

**Isolation.** A standalone harness (`prototypes/s4_cuda_mass/k2_repro/harness.c`, `build.sh`) calls the real
`s4_mass_k2.F90` kernels on `cudaMalloc` buffers (as the AMReX CUDA arena) in the driver's stage order and hashes
every array after every kernel, device build vs nvfortran host build. No implicit maps of field data were seen
(`NV_ACC_NOTIFY=3`: only the 24-byte `map(to:)` vectors move). `is_device_ptr` on the explicit-shape dummies, the
bounds and the value/reference scalars were all correct; the first diverging kernel was `s4k2_clip_terms`.

**Root cause 1 (the gross failure): nvfortran 26.9 offload defect.** `s4_mass_k2.F90` at `0bbb0c7cb5`, `s4k2_clip_terms`:
`:292` `!$omp private(d, flag, QMIN, QMAX, ...)` together with `:300-301` `QMIN = QMIN_IN; QMAX = QMAX_IN` (a
loop-invariant copy of a `value` dummy at the top of the collapsed loop body). `-Minfo=mp` reports
`288, Generating implicit private(qmin)` at the **target** line instead of the loop line, and the device threads
then compare against a wrong bound: every cell is flagged "clipped high" (`CF = CLIPPED|ACTIVE`), the density and
species applies run everywhere, and mass is created. The CPU paths (gfortran, nvfortran host) are correct.
Minimal reproducer: `prototypes/s4_cuda_mass/k2_repro/mini/` (`mini.F90`, `mini_main.c`, `build_mini.sh`), **run**
on the device: the `private(QMIN, QMAX)` form gives 4096/4096 elements wrong, the direct form 0 wrong.
Dropping `QMIN, QMAX` from `private` also avoids it (tested, **run**), but leaves a formally shared scalar.

**Root cause 2 (last-bit differences, found after fix 1): reassociation.** With fix 1, T1 passed but T0 did not
(~6e-16, ~6000 cells). nvfortran host and device results were equal, so the harness was run with gfortran vs
nvfortran host (box, **run**): the first difference is `s4k2_density_update` (`:207-214` at `0bbb0c7cb5`), `RHS = - DEL + (..)*RDX + (..)*RDY +
(..)*RDZ` (unparenthesised). Fortran allows a compiler to reassociate that sum; nvfortran -O2 does, gfortran and
the C++ K1 do not.

**Fix** (`amrex/s4_mass/s4_mass_k2.F90`, numerics, tolerances and `nofma`/fast-math-off flags unchanged):
- `s4k2_clip_terms`: local `QMIN` removed, `QMIN_IN` read directly; `QMAX = merge(RHOP(i,j,k), QMAX_IN, SPECIES /= 0)`
  (cell-dependent, cannot be hoisted); `QMIN, QMAX` no longer a private copy of an invariant. Same values, same
  comparisons.
- `s4k2_density_update`: explicit parentheses giving the FDS/K1 left-to-right order for `RHS` and for the
  corrector `0.5*((RHO*ZZ + RHOS*ZZS) - DT*RHS)`.

**Results after the fix:**
- Harness: device == nvfortran host == gfortran host, all arrays, all kernels incl. the corrector (**run**).
- Full `run_gpu_t1.sh --k2` on a cc 8.9 test-machine GPU (1 GPU, 1 process), **run**: K1 29/29 comparisons T1 and T0
  bitwise, 36 checks passed; **K2 29/29 comparisons T1 and T0 bitwise, 36 checks passed**; clip counts 978, 2520,
  154354, 7970 for both (= P1); run-to-run repeat bitwise; exit 0. K2 sync points 144 per case (the earlier 160
  counted density applies that the broken flags switched on).
- K2 host fallback (gfortran, CPU AMReX, 1/2/4 ranks) with the fixed source: 69/69 checks, 60/60 T0 (**run**).
- Timings from the same run, single sample, informational only (SUPERBEE g32, 8 steps): loop 6.1 ms (K1) vs
  8.3 ms (K2); mass kernels 2.6 ms vs 4.8 ms; clip 2.9 ms vs 2.9 ms. The K2 sync cost was not isolated.

**Left open:** report the privatisation defect to NVIDIA with `k2_repro/mini`; decide whether a K2 coding rule
("no private copies of loop-invariant values in offload loops; parenthesise sums whose order matters") goes into
the ADR; the R-39 items of §8 (removing post-K2 syncs, sync cost, profile) were not measured.
