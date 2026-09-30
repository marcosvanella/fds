# S4 K2: one Fortran source for CPU and GPU builds

Question: can `amrex/s4_mass/s4_mass_k2.F90` (OpenMP `target` directives) serve both offload and CPU builds?
Answer: **yes**, with a 3-macro include (`s4_omp.inc`); committed locally on `s4-cuda-mass` (not pushed).
"run" = executed, "compiled" = only built. Box = 8 shared cores (max 4 used), gfortran 14.2, OpenMPI 5.
Test machine = cc 8.9 test-machine GPU, 24 cores, nvfortran 26.9, ifx 2026.1.1, Intel MPI. Prototype code, logs, scripts:
`prototypes/s4_cuda_mass/single_src/` (`test machine/` holds the test machine copies). Standalone: no AMReX rebuild needed.

## 1. Method
- `bench.F90` (no AMReX) drives the real kernels in the `k2_repro/harness.c` stage order (23 stages, 3 species), N=16,
  and hashes every array after every stage. Reference: gfortran build with directives as comments (serial). Sign of
  zero is canonicalised in the hash (see 5.1). `harness.c` (414 per-array hashes, N=8) was also run for gfortran.
- Timing: N=96 box, min of 9 reps (us), OMP_NUM_THREADS=1 and 4, bound cores. Kernels: face_values+mw_correction
  (one direction), clip_terms, density_update, rho_z_p. Shared machine: expect a few tens of percent noise.
- Call modes: 0 direct from serial code; 1 every thread of an `!$omp parallel` calls the kernel on the whole box
  (FDS-style region); 2 inside `!$omp single`; 3 every thread calls it on its own k-slab (tile/MFIter style).
- Flags: exact arithmetic (`-ffp-contract=off`; nvfortran `-Mnofma`, `-gpu=cc89,nofma`; ifx `-fp-model=strict -no-fma`).

## 2. Compiler x mode matrix (all **run**; bitwise = per-stage hashes equal to the serial reference, modes 0/2/3)
| Build (new single source unless noted) | Compiles | Directive executes as | Threads (4 requested) | Bitwise | face+mw us, 1 thr / 4 thr |
|---|---|---|---|---|---|
| gfortran, no -fopenmp (old and new) | yes | comments | 1 | yes (harness 414/414) | 14.7k / - |
| gfortran -fopenmp, OLD K2 (`target teams loop`, `has_device_addr`) | yes | host fallback | 4 | yes | 18.1k / 6.5k |
| same + `-foffload=disable` | yes | host fallback (identical) | 4 | yes | same |
| gfortran -fopenmp, new (`parallel do collapse(2)`) | yes | host threads | 4 | yes | 16.0k / 4.3k |
| gfortran -fopenmp -DS4_OFFLOAD, new | yes | host fallback | 4 | yes | not timed |
| gfortran, OLD K2 with `is_device_ptr` for all | yes | host fallback | 4 | yes | - |
| nvfortran, no -mp | yes | comments | 1 | yes | 14.7k / - |
| nvfortran `-mp` or `-mp=multicore`, new | yes | host threads | 4 | yes | 15.1k / 3.8k |
| nvfortran `-mp=multicore -DS4_OFFLOAD` (`target teams loop`, `is_device_ptr`, `map(to:)`, host arrays) | yes | host multicore | 4 | yes | 14.1k / 3.5k |
| nvfortran `-mp=gpu` (no macro needed), device memory | yes | GPU | n/a | yes | 1.2k (GPU) |
| ifx, no -qopenmp | yes | comments | 1 | yes | 23.5k / - |
| ifx `-qopenmp`, new | yes | host threads | 4 | yes | 23.3k / 5.9k |
| ifx `-qopenmp -DS4_OFFLOAD` | yes | host | 4 | yes | - / 7.9k |
| ifx `-fiopenmp -fopenmp-targets=spir64 -DS4_OFFLOAD` | yes | host (no Intel GPU; on-device iterations 0) | 4 | yes | - / 7.6k |
GPU case (device memory via `omp_target_alloc`): old and new K2 give the same GPU times (face+mw 1198, clip_terms 395,
density_update 830, rho_z_p 327 us). ifx offload to a real Intel GPU: **not tested** (none). gfortran nvptx/amdgcn
offload: **not tested** (not installed). nvfortran `-mp=gpu` with host arrays: not tested (K2 needs device memory).
`is_device_ptr` and `map(to:)` on host arrays in host fallback: fine on gfortran, ifx and nvfortran multicore (no
compiler message; results equal). Error messages seen are only those in section 3.

Host timing, N=96, min us, 1 / 4 threads (kernels: clip_terms, density_update, rho_z_p; gfortran, **run**):
| gfortran | clip_terms | density_update | rho_z_p |
|---|---|---|---|
| serial (no -fopenmp) | 5.0k | 11.4k | 3.2k |
| OLD `target teams loop collapse(3)` | 7.9k / 2.3k | 14.3k / 5.6k | 3.5k / 1.6k |
| new `parallel do collapse(2)` | 5.8k / 1.5k | 13.7k / 3.2k | 3.3k / 0.85k |
| `parallel do collapse(3)` | 8.5k / 1.9k | 21.7k / 5.1k | 7.0k / 1.8k |
nvfortran multicore new, 1 / 4 thr: clip 3.7k / 1.7k, density 9.6k / 3.4k, rho_z_p 1.9k / 0.85k (serial 3.8k, 7.8k, 2.0k).
`collapse(3)` on a host loop costs +40..190 % at 1 thread (blocks vectorisation of the i loop): host uses `collapse(2)`.

## 3. Threading coexistence (**run**; probe.F90 counts executions per iteration and OS threads, data per thread)
- `target teams loop` on the host is multithreaded: 4 threads for gfortran (host fallback), nvfortran multicore, ifx.
  Speedup at 4 threads: about 2.3x (gfortran) versus 3.5x for `parallel do collapse(2)`; nvfortran and ifx about equal.
- Called from inside an `!$omp parallel` region (mode 1): every construct except an orphaned `!$omp do` runs the
  loop once per thread (4.00 executions per iteration). Non-idempotent kernels then race: mode 1 is wrong (hash differs)
  for every shape on all three compilers at 4 threads (correct at 1 thread). gfortran adds oversubscription for
  the target shapes: 16 OS threads (each caller starts its own 4-thread team); `parallel do` stays at 4 (nested
  region inactive); nvfortran and ifx use 4.
- Inside `!$omp single` (mode 2): target shapes use the full team, `parallel do` runs on 1 thread. Correct everywhere.
- Per-thread disjoint boxes (mode 3, tile style): correct for all shapes and compilers. gfortran timing at 4 threads
  (face+mw, clip, density): `target teams loop` 6.4k, 1.7k, 3.8k versus `parallel do collapse(2)` 3.9k, 1.5k, 3.4k.
- Orphaned `!$omp do collapse(3)` (workshares over the caller's team) is correct in mode 1 (bitwise equal) and would
  match the FDS `!$OMP DO` style, but it is not valid under `single`/mode 3 (hangs on gfortran, wrong under mode 3), so it
  is not a general choice. `!$omp loop` needs `bind()` (gfortran: "'bind' clause not specified on a 'loop' construct not
  nested inside another OpenMP construct"); `loop bind(parallel)` ran on 1 thread outside a region (gfortran, ifx) and
  hung under `single`; plain `!$omp loop collapse(3)` gives nvfortran multicore "NVFORTRAN-F-0000-Internal compiler
  error. expected region markers to always match -1 (k2_ss.F90: 135)". `!$omp loop` cannot call omp_* routines (nvfortran
  gpu "NVFORTRAN-S-1224-LOOP construct may not contain calls to the OpenMP Runtime API"; gfortran similar).
- `metadirective`: gfortran 14.2 "Error: Unclassifiable OpenMP directive". nvfortran 26.9 and ifx accept it. With
  `when(device={kind(gpu)}: ...)` nvfortran `-mp=gpu` picked the host branch; with `target_device={kind(gpu)}` it launched
  a CUDA kernel (mini test only, not K2; 1 launch seen). Not usable as the single mechanism (gfortran), so not adopted.
- MPI (no MPI calls in the bench; ranks x threads with binding, face+mw us, N=64, **run**). Box, OpenMPI, ranks x threads
  1x4 / 2x2 / 4x1: serial 4.5k / 4.4k / 4.5k versus new `parallel do` 1.15k / 2.3k / 4.5k. Test machine, Intel MPI, nvfortran
  multicore: 2x2 2.1k vs serial 4.4k; 4x1 4.4k vs 4.8k. No oversubscription when cores = ranks x threads. Default
  (unbound) rank placement not characterised. Real MPI plus GPU (several ranks on one GPU): not tested.
- FDS style (`grep` in `Source/mass.f90`, divg, velo): `!$OMP PARALLEL PRIVATE(...)` ... `!$OMP DO` ... `!$OMP END PARALLEL` and
  `!$OMP PARALLEL DO ... SCHEDULE(STATIC)` inside a routine, kernels called from serial code, so mode 0 is the FDS-like case.

## 4. Recommended pattern (verified, **run** on gfortran, nvfortran, ifx; committed as `s4_omp.inc` + `s4_mass_k2.F90`)
Every kernel keeps the literal sentinel; only the directive tail and the clause pieces are macros:
```
!$omp S4_LOOP S4_DEV((RZP, RHOP, ZZP)) S4_MAPTO(MWR_Z) private(DOT, MW_G, N)
```
`s4_omp.inc` (included once, after the S4K2_COLL / S4K2_CL block):
- `#if defined(S4_OFFLOAD) || defined(__NVCOMPILER_OPENMP_GPU)` offload: `S4_LOOP` = `target teams loop collapse(S4K2_COLL) S4K2_CL`,
  `S4_MAPTO(x)` = `map(to: x)`, `S4_DEV(l)` = `is_device_ptr l` if `__NVCOMPILER` else `has_device_addr l`.
- else host: `S4_LOOP` = `parallel do collapse(2)`, `S4_MAPTO`, `S4_DEV` empty (the `map(to:)` vectors are read-only shared).
Lines that differ by compiler (all inside `s4_omp.inc`): (1) the offload switch: nvfortran sets `__NVCOMPILER_OPENMP_GPU`
itself with `-mp=gpu`; gfortran (nvptx/amdgcn) and ifx (spir64) need `-DS4_OFFLOAD` because `_OPENMP` is defined for host-only
builds too; (2) `is_device_ptr` (nvfortran rejects `has_device_addr`) vs `has_device_addr`; (3) `collapse(3)` vs `collapse(2)`.
The 10 duplicated `#if/#else` blocks of the old file are gone (K2: 16 lines added, 75 removed). The `S4K2_*` tuning macros
(NOWAIT, THREADS, TEAMS, COLL) keep working on the offload branch (nvfortran compile with `-DS4K2_NOWAIT` checked, **compiled**).
Preprocessor limits (**run** on all three):
- `-cpp` (gfortran; `-fpp` is rejected: "unrecognized command-line option '-fpp'; did you mean '-cpp'"), `-Mpreprocess` or the `.F90` suffix
  (nvfortran), `-fpp` (ifx). `.F90` suffix already turns it on for the other two.
- A macro that expands to the whole `!$omp ...` line works on gfortran and nvfortran (directive active, team size 4 seen) but fails on
  ifx (`error #5082: Syntax error, found END-OF-STATEMENT when expecting one of: ( % . = =>`). Clause-only tail macros work everywhere,
  also after a continuation `&`, and `S4_DEV((a, b))` with a doubled parenthesis keeps the commas inside one macro argument.
- Fallback if a compiler ever rejects a tail macro: the `#ifdef`-wrapped duplicate directive (the old layout), one pair per kernel.
- `s4_omp.inc` must sit next to the source or be on `-I`; nvfortran found it beside the source with no `-I`.

## 5. Verification of the committed change
- CPU T1 on the development machine, K2 host build (gfortran -fopenmp, 1/2/4 ranks, tiled 2/3 thread checks): **run**, 60 comparisons, 69 checks
  passed, "S4 T1: ALL PASS" (same numbers as before the change; log `t1_k2_single_src.log`).
- gfortran harness.c, 6 flag combinations plus 3 alternative shapes: 414/414 hashes equal. nvfortran and ifx: table above.
- GPU: `-mp=gpu` build of the committed source gives the same GPU results and times as the old source in the standalone bench
  (not re-run through the AMReX CUDA T1; the device-side directive text is unchanged from the old `__NVCOMPILER` branch).

## 6. Findings and open problems
1. nvfortran -O2 vs gfortran: `clip_gather` `DELTA` differs only in the sign of zero (`0 - T` folded to `-T`). Values equal; `-Kieee`
   removes it. A bytewise T0 compare against K1 could flag it; the earlier GPU T0 compare passed.
2. Host `parallel do` in K2 is safe under the AMReX tiled MFIter loop only because nested regions are inactive by default; if nested
   parallelism is enabled (OMP_MAX_ACTIVE_LEVELS>1) threads multiply. Not tested.
3. A GPU build must not use the host shape: `-mp=gpu` with `parallel do` on device pointers produced no output (crash), so the switch must
   follow the actual build (auto for nvfortran, `S4_OFFLOAD` otherwise).
4. Calling K2 on the same box from several threads is wrong for any shape (the kernels are not idempotent); use disjoint boxes/tiles.
5. Harness only: nvfortran -O2 miscompiled the bench driver in mode 3 (bench built -O1); ifx loses results when a pointer-valued
   function is an actual argument; neither touches K2.
6. Not tested: ifx on an Intel GPU, gfortran offload, ifx/nvfortran through AMReX, Windows, more than 4 threads on the development machine.

## 7. Closing GPU T1 run through AMReX CUDA (single-source K2)
Result (**run**, cc 8.9 test-machine GPU, 1 GPU, 1 process, worktree commit cf05615691 shipped as a tarball, AMReX CUDA build reused, only the S4
driver and kernels rebuilt, `run_gpu_t1.sh --k2 --skip-amrex-build`, executed once):
- K1 (nvcc, `--fmad=false`): 29 comparisons, 29 bitwise (T0), 36 checks passed, "S4 T1: ALL PASS".
- K2 (`s4_mass_k2.F90` including `s4_omp.inc`, nvfortran 26.9 `-O2 -mp=gpu -gpu=cc89,nofma -Minline -Minfo=mp`, link
  `-mp=gpu -gpu=cc89 -Mnomain -cuda -c++libs`): 29 comparisons, 29 bitwise (T0), 36 checks passed, "S4 T1: ALL PASS". The offload branch
  of the macros was the active one (device code generated for the target regions; 144 K2 sync points in the decoupled case).
- Final line: "S4 GPU T1: PASS (K1 and K2)". The two 'FAIL' lines in each suite are the informational P1-rules control (`s4.coarse_mask=0`).
- This closes the open item in section 5 (single-source K2 not yet through the AMReX CUDA T1). No source change was needed.
- Logs: `prototypes/s4_cuda_mass/single_src/test machine/t1close/` (`t1close_run.log`, `build_gpu_k1.log`, `build_gpu_k2.log`).
