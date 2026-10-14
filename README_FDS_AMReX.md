# FDS-AMReX branch: build and run guide

This branch runs FDS on top of [AMReX](https://amrex-codes.github.io/amrex/) block-structured adaptive mesh refinement (AMR).
The `fds_amr` executable reads a normal FDS input file. Driver flags after the input file control levels, regridding and the pressure solver.
This page covers dependencies, CPU and GPU builds, the GPU kernel generator, and three runnable cases. Inputs and launch scripts are in
[`Verification/Adaptive_Mesh_Refinement/`](Verification/Adaptive_Mesh_Refinement/). Design documents are in
[`Manuals/FDS_AMReX_Design/`](Manuals/FDS_AMReX_Design/README.md). Driver internals are in [`Source/driver/README.md`](Source/driver/README.md).

Status labels used below: **verified** = the command ran on the Linux reference box (GNU/OpenMPI, CPU, np 4, a few steps).
**not verified here** = taken from the project's build records, not re-run for this page.

---

## 1. Dependencies

| Component | Version / pin | Notes |
|---|---|---|
| Compilers (GNU) | GCC/gfortran 14.2.0 | Debian 13 packages |
| MPI (GNU) | Open MPI 5.0.7 | `mpicc`, `mpicxx`, `mpifort` |
| Compilers (Intel) | oneAPI `icx`/`icpx`/`ifx` with Intel MPI (`mpiicx`, `mpiicpx`, `mpiifx`) | not verified here |
| CMake | 3.31 | |
| AMReX | `development` at commit `99ddfda6722b63f7e128ee6d6822ff29046bb7a9` (reports itself as 26.09) | 3-D, double, MPI + OpenMP, FFT on, HYPRE on, EB on, Fortran interfaces on, static, CPU (`AMReX_GPU_BACKEND=NONE`) |
| HYPRE | FireX pin `63331f19c72b2ec363824907c0575a6edf7bbf3c` = `v2.32.0-24-g63331f19c` | CPU, MPI, static. The FireX CMake labels it "3.0.0" but the code is 2.32.0+24: ask for `find_package(HYPRE 2.32.0)` |
| FFT | FFTW 3.3.10 (`libfftw3-dev`, `libfftw3-mpi-dev`) | used by AMReX `FFT::Poisson`, the default level-0 pressure solver |
| SUNDIALS | 7.5.0 | as in FireX (`USE_SUNDIALS=ON`) |
| Python 3 + numpy | any recent | to generate initial-velocity CSVs; `fparser` 0.2.5 for the GPU generator |

**One HYPRE only.** AMReX and FDS must link the *same* HYPRE build (project risk R-30). The CMake checks this
(`Source/driver/CMakeLists.txt` compares the HYPRE version string inside `libamrex_3d.a` with the one FDS links).

### 1.1 HYPRE (as built for FireX)
Source: `https://github.com/hypre-space/hypre.git`, commit `63331f19c7` (the `GIT_TAG` in the top-level `CMakeLists.txt`). CMake source dir is `src/`.
```bash
git clone https://github.com/hypre-space/hypre.git && git -C hypre checkout 63331f19c72b2ec363824907c0575a6edf7bbf3c
cmake -S hypre/src -B hypre-build -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DHYPRE_ENABLE_MPI=ON \
      -DCMAKE_C_COMPILER=mpicc -DCMAKE_INSTALL_PREFIX=$HYPRE_ROOT
cmake --build hypre-build -j8 && cmake --install hypre-build
```
HYPRE is C only. Intel: use `-DCMAKE_C_COMPILER=mpiicx`.

### 1.2 AMReX
```bash
git clone https://github.com/AMReX-Codes/amrex.git && git -C amrex checkout 99ddfda6722b63f7e128ee6d6822ff29046bb7a9
cmake -S amrex -B amrex-build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$AMREX_ROOT \
      -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_Fortran_COMPILER=mpifort \
      -DAMReX_SPACEDIM=3 -DAMReX_PRECISION=DOUBLE -DAMReX_MPI=ON -DAMReX_OMP=ON -DAMReX_FORTRAN=ON -DAMReX_FORTRAN_INTERFACES=ON \
      -DAMReX_FFT=ON -DAMReX_HYPRE=ON -DHYPRE_ROOT=$HYPRE_ROOT -DAMReX_EB=ON -DAMReX_LINEAR_SOLVERS=ON \
      -DAMReX_PARTICLES=ON -DAMReX_GPU_BACKEND=NONE -DAMReX_BUILD_SHARED_LIBS=OFF
cmake --build amrex-build -j8 && cmake --install amrex-build
```
(Options as recorded in the reference install's `CMakeCache.txt`. Unlisted options keep their defaults.)

---

## 2. Building FDS-AMReX (CPU)

### 2.1 Patches that are not committed yet
Some fixes outside `Source/driver/` are kept as patches in [`Source/driver/patches/`](Source/driver/patches/). Every current build applies them on a **copy**
of the tree (e.g. a `git worktree`), never in the shared checkout:

| Patch | File | State at the current tip |
|---|---|---|
| `0013-mesh-dddt-pointer.patch` | `Source/mesh.f90` | **apply** (not committed) |
| `0014-init-zero-closed-wall.patch` | `Source/init.f90` | **apply** (not committed) |
| `0015-fine-wall-bc-wire.patch` | `Source/init.f90` | **apply** (not committed) |
| `0016-regrid-old-fine-face-period.patch` | | already committed, do not apply |
| `0017-fine-wall-noflux.patch` | | already committed, do not apply |
| `0001`–`0012` | | already in the branch history |

Check a patch with `git apply --check <patch>`: if it applies it is still needed. If `git apply --check -R <patch>` succeeds, it is already in.

```bash
git worktree add --detach ../fds-amr-build HEAD
cd ../fds-amr-build
for p in 0013 0014 0015; do git apply Source/driver/patches/$p-*.patch; done
```

### 2.2 GNU + Open MPI (verified)
```bash
cmake -S . -B build-gnu -DCMAKE_BUILD_TYPE=Release -DUSE_AMREX=ON \
      -DAMReX_DIR=$AMREX_ROOT/lib/cmake/AMReX \
      -DUSE_SYSTEM_HYPRE=ON -DCMAKE_PREFIX_PATH=$HYPRE_ROOT \
      -DUSE_SYSTEM_SUNDIALS=ON -DSUNDIALS_DIR=$SUNDIALS_ROOT/lib/cmake/sundials \
      -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_Fortran_COMPILER=mpifort \
      -DBUILD_TESTING=OFF
cmake --build build-gnu -j8          # executable: build-gnu/fds_amr
```
This matches the reference builds (Release `-O3`, `USE_OPENMP=ON`, `USE_HYPRE=ON`, `USE_SUNDIALS=ON`). Use `-DBUILD_TESTING=ON` for the ctest unit tests.
Run with one OpenMP thread: `export OMP_NUM_THREADS=1`.

**MKL is off by default** (`FDS_USE_MKL=OFF`). Static MKL exports FFTW symbols that override the real
`libfftw3` and break the AMReX FFT pressure solve. Add `-DFDS_USE_MKL=ON` only if you need MKL
(e.g. the PARDISO/cluster-sparse solver) and do not use the AMReX FFT solver.

### 2.3 Intel oneAPI (icx/ifx) + Intel MPI (not verified here)
Same line with the Intel wrappers. Build HYPRE and AMReX with the same wrappers:
```bash
source <oneAPI>/setvars.sh
cmake -S . -B build-intel -DCMAKE_BUILD_TYPE=Release -DUSE_AMREX=ON \
      -DAMReX_DIR=$AMREX_INTEL_ROOT/lib/cmake/AMReX \
      -DUSE_SYSTEM_HYPRE=ON -DCMAKE_PREFIX_PATH=$HYPRE_INTEL_ROOT \
      -DUSE_SYSTEM_SUNDIALS=ON -DSUNDIALS_DIR=$SUNDIALS_INTEL_ROOT/lib/cmake/sundials \
      -DCMAKE_C_COMPILER=mpiicx -DCMAKE_CXX_COMPILER=mpiicpx -DCMAKE_Fortran_COMPILER=mpiifx \
      -DCMAKE_CXX_FLAGS=-fp-model=precise -DBUILD_TESTING=OFF
cmake --build build-intel -j8
```
Facts from the Intel validation (see `Source/driver/README.md`): use `-fp-model=precise` for C++ (otherwise unit test `tile_race` fails) and
`mpiicpx`/`mpiifx`. Do **not** pass `--oversubscribe` to Intel `mpiexec` (with the scripts: `MPIRUN=mpiexec`). The plain-FDS pyramid reference (§6) was run with an Intel
build on an NVIDIA test machine. FDS-AMReX builds on that machine are recorded as stopping at step 1 of the pyramid case. This has not been diagnosed; the box binary runs it.

---

## 3. GPU build (NVHPC, OpenMP target, "K2")

Decision ([ADR-001](Manuals/FDS_AMReX_Design/adr/ADR-001-driver-architecture.md), accepted): GPU kernels are **K2** = Fortran with OpenMP `target`
offload, compiled with `nvfortran`. **K1** (restricted C++ `amrex::ParallelFor`) is only a per-kernel fallback.

**What exists today:** the K2 kernels come from a generator (`s5gen.py`, below). They are checked bit for bit against the upstream loops on the host and,
with side-by-side drivers, on an NVIDIA GPU. **They are not wired into `fds_amr` yet.** The FDS-AMReX CMake has no GPU option, so there is no
end-to-end GPU `fds_amr` build. A GPU build covers only the loops the generator translates (§3.6). Everything else stays on the host.
An AMReX CUDA install exists for experiments (`AMReX_GPU_BACKEND=CUDA`). It is not used by `fds_amr`.

Compiler flags for K2 (generator `gpu/CMakeLists.txt`, ADR-001 rule 4): compile `-mp=gpu -gpu=ccXX,nofma -Minline -Minfo=mp`, link
`-mp=gpu -gpu=ccXX -Mnomain -cuda -c++libs`. Never `-fast`. The CPU reference is compiled `-Mnofma -Kieee` (matches gfortran `-ffp-contract=off`).

### 3.1 s5gen.py: where and inputs
Branch **`s5-gen`**. Directory `amrex/s4_mass/s5_gen/` on that branch (check it out in a separate worktree). Main docs there: `README.md`, `gpu/README_device.md`,
`callsite/README_callsite.md`.
- **Input:** the upstream FDS Fortran `Source/*.f90`, read-only. By default it is read through git at the upstream commit named in `s5_markers.toml`.
  With `--src-dir <dir>` it reads a synced directory instead. The generator resolves module symbols across the whole directory (parser: `fparser` 0.2.5).
- **Markers** (which loops become kernels):
  - `s5_markers.toml`: a sidecar file, so upstream stays unmodified. `[policy.*]` sections (wall tables, arrays, assoc) plus one `[[kernel]]` per loop:
    file, routine, line range, `anchor` (first statement, catches drift) and `args` (pinned argument order).
  - Inline tags in upstream source (plain comments for other compilers), read with `--inline <file.f90>`:
    ```
    !$GPU KERNEL <name> [HOST] [UNIQUE] [PRIVATE(a,b)] [REDUCTION(op:v)] [COLLAPSE(n)] [LAYOUT(exact)] [ARGS(a,b,...)]
       ... one K,J,I (or wall) loop nest ...
    !$GPU END KERNEL <name>
    ```
    `HOST` keeps a loop on the CPU. `UNIQUE` asserts that distinct walls write distinct cells. An inline tag and a sidecar entry generate identical text.
  - `markers/golden_signatures.json`: the golden signature of every kernel.

### 3.2 Run it
```bash
cd amrex/s4_mass/s5_gen                       # on branch s5-gen
python3 -m venv venv && venv/bin/pip install fparser==0.2.5
venv/bin/python s5gen.py                      # upstream via git at the pinned commit -> generated/
venv/bin/python s5gen.py --src-dir <dir with *.f90>
venv/bin/python s5gen.py --inline mass.f90    # also scan inline !$GPU KERNEL tags
venv/bin/python s5gen.py --update-golden      # only after a reviewed, intended kernel change
```
Exit codes: 0 ok, 2 a loop cannot be translated or its argument set differs from the pinned one, 3 golden signature differs, 4 parse error. Each failure names the kernel and the reason.

### 3.3 Outputs (`generated/`)
`s5gen_k2.F90` (the K2 kernels), `s5gen_k2.H` (C/C++ interface), `s5gen_args.json` (argument manifest: name, role, dims, intent),
`s5gen_report.md` (per kernel: source lines, directive, pinned arguments, derived intents, private scalars, storage of each array).

### 3.4 Bitwise test
```bash
./build_and_test.sh all                          # flag sets O0, O2, O0omp, O2omp, O2omp_off, O2omp_dpd; both callee directives (dpd, bind)
venv/bin/python test/test_generator_checks.py
```
Each generated kernel runs against the verbatim upstream loop (extracted by line range) on randomised inputs at several box sizes. Every output array must
match bit for bit. Use `-O0 -fopenmp` (`O0omp`) to catch a dropped `private` clause. On the GPU, `gpu/dev_stage.sh` → `dev_build.sh`
(`NVHPC_ROOT=…  GPU_CC=ccNN`) → `dev_run.sh NX NY NZ TAG` prints one hash per output array for host and device builds and compares them.
Overlapping test runs take a shared lock via `test/with_s5lock.sh`.

### 3.5 How a kernel gets wired in
1. Add a `[[kernel]]` entry to `s5_markers.toml` (or an inline `!$GPU KERNEL` tag), regenerate, review `s5gen_report.md`, then `--update-golden`.
2. Build: `integration/s5gen.cmake` provides `s5gen_add_kernels(<target> <upstream-Source-dir>)`, which runs the generator at build time and adds
   `s5gen_k2.F90` to the target. The build fails on parse error, drift or signature change. `integration/s5gen.mk` does the same for make.
3. At the call site, replace the upstream loop with a call to the generated kernel, passing the arguments in the pinned order from `s5gen_args.json`.
   Wall loops need the flat WALL tables gathered first. `callsite/s5_gsfv_dispatch.F90` shows a call-site shim for `GET_SCALAR_FACE_VALUE`.
4. Daily upstream merges: `port_merge_check.py OLD NEW` sorts each hunk as body-only, surface-changing or outside-kernel.

**Gap:** no `fds_amr` call site uses a generated kernel yet. The repo does not show who owns the first wiring step into the FDS-AMReX build or which
routine goes first. (A question to the FDS Legacy Mapper on this point could not be sent from this session. Still open.)

### 3.6 Loops translated today
139 kernels in `generated/s5gen_report.md` (upstream `36975d7`), from `divg.f90`, `func.f90`, `velo.f90`, `mass.f90` and the wall routines:
`DIVERGENCE_PART_1` cell loops and callees (`rho_d_dzd`, `kdtd`, `dp_species`, …), species advection / mass (`zzs_pred`, `rho_zz_*`, `tdm_*`,
`settle_*`), viscosity models (`mu_dns`, `visc_wale`, `visc_vreman`, `mu_deardorff`, `mu_consmag`), velocity predictor/corrector (`vpred_*`, `vcorr_*`),
`VELOCITY_FLUX` edge and baroclinic nests (`vflux_*`, `baro_*`, `cw3_*`, `cw4p_*`, `cw4c_*`), `CHECK_STABILITY`/`CHECK_DIVERGENCE` reductions
(`cfl_max`, `cfl_wall_max`, `vn_max`, `div_extrema`), wall loops (`wall_*`, `wall_*_wl`), the six `GET_SCALAR_FACE_VALUE` limiters (`gsfv_*`) and
`PATCH_VELOCITY_FLUX` (`pvf_*`). Not translated (examples from the README there): loops that call `EVALUATE_RAMP`, K,I nests with fixed J,
neighbour-mesh loops, and the solid-phase `L1486` nest (blocked on a callee that uses a module variable). See `STATUS.txt` on that branch.

---

## 4. Running: the vortex ring on a uniform grid (several meshes)

Case: a 3-D vortex ring (Rosenhead/Biot-Savart initial velocity, Γ = −0.8, R = 0.3, core 0.2, moving down) in the unit cube. INERT floor, MIRROR top,
PERIODIC sides, ν = 1e-2. Inputs:
- `vortex_ring_48_4m.fds`: 48³ cells, 4 meshes of 24×24×48 (x/y split 2×2×1). Initial velocity from `Source/driver/tests/cases/vortex_ring_48_4m_uvw_m{1..4}.csv`.
- `vortex_ring_96_2x2x1.fds`: the uniform **fine** variant, 96³ cells, 4 meshes of 48×48×96. Its CSVs (~38 MB) are generated by the script on first run
  (`Source/driver/tests/make_vortex_ring_case.py 96 <dir> --sx 2 --sy 2 --sz 1`).
- `vortex_ring_48_4x4x1.fds`: the same ring at Re 800 (ν = 1e-3, `T_END` 5), 48³ cells on 16 meshes (4×4×1, 4 per rank, np 4); the input of the long run. Its CSVs come from
  `make_vortex_ring_case.py 48 <dir> --sx 4 --sy 4 --sz 1` (generated by `_common.sh` on first run).

The `.fds` files are plain FDS input: the pressure solver and the adaptivity are chosen only through `fds_amr` flags (the one exception is the plain-FDS
reference `pyramid_96_4x4x1.fds`, see section 5).

**Split meshes in x and y only.** Cases split in z (e.g. `*_8m.fds`, 2×2×2) are known to be broken.

All runs use `fds_amr <input> --two-level-run --levels 1 ...`. `--levels 1` means "level 0 only", i.e. the uniform grid, through the same time loop,
logs and plotfiles as the adaptive runs. The scripts copy the input and CSVs into `out_<chid>/` and run there:
```bash
cd Verification/Adaptive_Mesh_Refinement
export FDS_AMR=/path/to/build-gnu/fds_amr   # NP=4 STEPS=10 by default; MPIRUN="mpirun --bind-to none --oversubscribe"
./run_vortex_fft.sh
```

### 4.1 The pressure solver
Pick the backend with `--pressure-backend auto|mlmg|hypre` and the Krylov method with `--krylov none|cg|gmres|bicgstab`
(`Source/driver/main.cpp`, `TimeLoop.cpp`). Without these flags, level 0 uses the direct FFT solver. `--pressure-tol T` sets the relative tolerance.
`--pressure-csv FILE` names the per-solve CSV (default `<chid>_pressure.csv` when any solver flag is given). Expert keys `hypre.<key>=<value>` /
`pressure.<key>=<value>` go to AMReX ParmParse (e.g. `hypre.strong_threshold=0.5`).

| Solver | Flags | Script | Status |
|---|---|---|---|
| Parallel FFT Poisson (AMReX `FFT::Poisson`, FFTW) | *(none)* | `run_vortex_fft.sh` | verified |
| BiCGStab + MLMG (AMReX geometric multigrid preconditioner) | `--pressure-backend mlmg --krylov bicgstab` | `run_vortex_mlmg_bicgstab.sh` | verified |
| HYPRE PCG + BoomerAMG | `--pressure-backend hypre --krylov cg` | `run_vortex_hypre_pcg.sh` | verified |

Exact commands (from inside the run directory, with the CSVs next to the input):
```bash
mpirun --bind-to none --oversubscribe -np 4 $FDS_AMR vortex_ring_48_4m.fds --two-level-run --steps 10 --outdir . --chid vr48_fft   --levels 1
mpirun --bind-to none --oversubscribe -np 4 $FDS_AMR vortex_ring_48_4m.fds --two-level-run --steps 10 --outdir . --chid vr48_mlmg  --levels 1 --pressure-backend mlmg  --krylov bicgstab
mpirun --bind-to none --oversubscribe -np 4 $FDS_AMR vortex_ring_48_4m.fds --two-level-run --steps 10 --outdir . --chid vr48_hypre --levels 1 --pressure-backend hypre --krylov cg
mpirun --bind-to none --oversubscribe -np 4 $FDS_AMR vortex_ring_96_2x2x1.fds --two-level-run --steps 10 --outdir . --chid vr96_uniform --levels 1   # + any solver flags
```
Notes: FFT is direct, so `--krylov` with FFT is an error. With obstacles (masked cells), `auto` picks HYPRE. `--krylov cg` needs a symmetric matrix.
A good run ends with `MULTI-LEVEL RUN PASS`. Per-step data is in `<chid>_multi_level.csv`.

---

## 5. Adaptivity

### 5.1 AMReX words in plain language
- **Level:** one resolution of the whole problem. Level 0 covers the domain. Level 1 is finer (here 2× in each direction) and covers only part of it, and so on.
- **Box:** a rectangular block of cells, given by its lower and upper cell indices on a level.
- **BoxArray:** the list of boxes that make up one level. Boxes do not overlap. They are spread over MPI ranks (the *DistributionMapping*).
- **MultiFab:** a field (e.g. U, RHO) stored on a BoxArray: one array per box, each with some ghost cells, possibly several components.
- **Ghost cells:** a layer of extra cells around each box. They are filled from neighbouring boxes, from the coarser level, or from boundary conditions, so a stencil can read past the box edge.
- **Tagging:** marking the cells that need more resolution (here: high vorticity).
- **Regrid:** building new finer-level boxes around the tagged cells and moving the data onto them (copy where the old fine level existed, interpolate from coarse elsewhere).
- **Blocking factor:** every fine box has sizes and corners that are multiples of this number of cells. It keeps the boxes aligned and coarsenable.
- **Max grid size:** the largest box edge allowed. Bigger regions are chopped into several boxes.
- **Average-down:** replace coarse cells that lie under fine cells with the average of those fine cells, so both levels agree.
- **Coarse-fine interface:** the boundary between a fine level and the coarse cells around it. Fluxes there are corrected ("refluxing"/flux overwrite) so mass is conserved, and fine ghost cells there are interpolated from the coarse level.

### 5.2 Flags (`--two-level-run`)
| Flag | What it does |
|---|---|
| `--two-level-run` | Run the multi-level (AMR) time loop. Must be the first flag after the input file. |
| `--levels N` | Number of levels (`max_level = N-1`, default 2). `--levels 1` = level 0 only. |
| `--steps N` | Number of time steps. The loop stops on steps, not on `T_END`. |
| `--regrid-every N` | Re-tag and regrid every N steps (0 = keep the initial fine boxes fixed). |
| `--tag-vort auto\|<thr>` | Tag cells by undivided vorticity |ω|·Δx. `auto` sets the threshold to 0.45 × the current maximum on level 0, so it follows the flow. A number gives a fixed threshold (the input headers print a hint). |
| `--tag-obst` | **In progress (P1), not in the committed tip:** tag cells near obstacles. |
| `--blocking B` / `--maxsize M` | Blocking factor and max grid size of fine boxes. Use **equal** values, e.g. `--blocking 8 --maxsize 8` (8³ cubes). Never use blocking < maxsize. |
| `--ratio R` | Refinement ratio (default 2). |
| `--plotfile [DIR]` | Write AMReX plotfiles (`plt_<chid>_s<step>`, fields RHO, U, V, W, OMEGA_Y) under DIR (default `--outdir`). Readable in ParaView, VisIt, yt, Amrvis. |
| `--plot-every N` | Plotfile every N steps (plus a `_final` one). |
| `--projection ON\|AUTO\|OFF` | Divergence projection after the initial grid is built (POST_REGRID_PROJECTION). `AUTO` (default) keeps the result only if divergence went down enough. `ON` always keeps it. `OFF` skips it. |
| `--outdir D`, `--chid C`, `--log-every N` | Output directory, case name of outputs, log cadence. |

### 5.3 Adaptive vortex ring (verified)
```bash
./run_vortex_adaptive.sh
# = mpirun --bind-to none --oversubscribe -np 4 $FDS_AMR vortex_ring_48_4m.fds --two-level-run --steps 10 --outdir . --chid vr48_amr \
#     --levels 2 --regrid-every 5 --tag-vort auto --blocking 8 --maxsize 8 --plotfile plt --plot-every 5
```
The solver flags of §4.1 can be added. Plotfiles are a few MB each, so keep `--plot-every` sparse on long runs.

---

## 6. The pyramid case (vortex ring over obstacles)

The same ring (ν = 1e-3, `T_END` 6) moves down onto a three-step `&OBST` pyramid. Inputs in `Verification/Adaptive_Mesh_Refinement/`:
- `pyramid_48_4m.fds`: FDS-AMReX input, 48³, 4 meshes (2×2×1). It uses the same initial-velocity CSVs as `vortex_ring_48_4m.fds`.
- `pyramid_96_4x4x1.fds`: plain-FDS 96³ reference, 16 meshes (4×4×1). It needs `vortex_ring_96_4x4x1`-style CSVs from
  `make_vortex_ring_case.py 96 <dir> --sx 4 --sy 4 --sz 1 --re-scale 10` (file names as in the input). The file already carries the
  reference pressure solver `&PRES SOLVER='UGLMAT', VELOCITY_TOLERANCE=1.E-6, MAX_PRESSURE_ITERATIONS=100 /` (the default FFT solver hits the
  100-iteration cap every step on this case). Run it with plain FDS on 16 ranks (or 8 ranks with two meshes each after editing `MPI_PROCESS`). Not run for this page.

Level-0-only FDS-AMReX baseline (verified):
```bash
./run_pyramid_l0.sh
# = mpirun --bind-to none --oversubscribe -np 4 $FDS_AMR pyramid_48_4m.fds --two-level-run --steps 10 --outdir . --chid pyramid_l0 \
#     --levels 1 --projection OFF --log-every 10
```
`--projection OFF` because the initial projection does not know about solids yet. The level-0 solve is one pressure solve per stage, with no FDS
pressure iteration on the obstacle normal velocity. For the full baseline run use `STEPS=760 ./run_pyramid_l0.sh --plotfile plt --plot-every 5` (≈ T_END).

---

## 7. `simple_test` fire case (species / combustion)

Propane burner in a small open compartment. Input: `Verification/Adaptive_Mesh_Refinement/simple_test.fds`
(same physics as the Role-3 species gate runs). Single `&MESH` (36×24×24), so the run scripts default to `NP=1`.

**Physics.** Species transport and mixing-controlled combustion (`&REAC FUEL='PROPANE'`), radiation off
(`&RADI RADIATION=.FALSE.`), one `OPEN` vent, burner `SURF` (`HRRPUA=1000`) on a small `&OBST` face that stays on
level 0 (D-010: stateful-wall coverage freeze; this burner has no `MATL_ID`, so the freeze does not bind it,
but the burner geometry is on L0). `&AMR BLOCKING_FACTOR=2` is present; `--two-level-run` ignores `&AMR MAX_LEVEL`
(use `--levels`).

**Patches.** Needs the prime WALL_BC + flux-match gate + fine-level rebind path (Role-3 patches 0004 / 0005 / 0006).
On current HEAD those land as `patches/0023-prime-wall-bc.patch`, `patches/0024-fluxmatch-gate.patch`,
`patches/0025-fine-level-rebind.patch` (commits `93e321629b`, `d3699aa628`, `ff24969ae1`). No extra patch apply
is required when building from that tip.

**Scripts** (source `_common.sh`; short smoke `STEPS=10` by default):
```bash
cd Verification/Adaptive_Mesh_Refinement
export FDS_AMR=/path/to/fds_amr
./run_simple_test_l0.sh          # --levels 1
./run_simple_test_l2.sh          # --levels 2 --regrid-every 5 --tag-vort auto
./run_simple_test_l3.sh          # --levels 3 --regrid-every 5 --tag-vort auto
```
Exact long-run flags (pressure: default FFT on the single unmasked mesh at L0; composite MLMG once fine levels exist —
same as the gate binaries, no `--pressure-backend` override):
```bash
# level 0 only (uniform), to T_END=5
$FDS_AMR simple_test.fds --two-level-run --levels 1 --steps 900
# 2-level adaptive (r5_wbc)
$FDS_AMR simple_test.fds --two-level-run --levels 2 --steps 900 --regrid-every 5 --tag-vort auto
# 3-level adaptive (r3lev)
$FDS_AMR simple_test.fds --two-level-run --levels 3 --steps 2000 --regrid-every 5 --tag-vort auto
```

**Optional `--tag-temp [T]`** (default off; `T` defaults to 400 K): also tag level cells with `TMP > T` (as `TAG_TMP_ABOVE = T - 293.15`, OR with `--tag-vort`), so the hot plume and ceiling layer stay refined. Adaptive runs only (`--regrid-every > 0`). The burner and the D-010 no-refine buffer (`--nofine`, default 2 level-0 cells) stay on level 0.
Example: `./run_simple_test_l2.sh --tag-temp`.

**Expected gate output (measured).** Limit for mass/species is `1e-6 × (initial mass + through-flow)` on an open domain;
ratios are `|balance residual| / limit` (PASS ≤ 1). Divergence gate: worst fine-level `|div u − D|` ≤ 1e-9 and a
count of steps with fine div > max(2·L0, 1e-10). Flux-match stale-register count must be 0.

| Case | Script | Steps / stop | Wall time | Mass ratio | Z1 / Z2 / Z3 ratios | Worst \|div\| (level, step) | FLUXMATCH | Gate |
|---|---|---|---|---|---|---|---|---|
| L0 only | `run_simple_test_l0.sh` | 451 steps, `T_END=5` | ~11 min | 0.222 | 0.390 / 0.598 / 1.05e-11 | 1.70e-13 (L0) | n/a (no C/F; stale 0) | PASS |
| 2-level | `run_simple_test_l2.sh` | 869 steps, `T_END=5` | ~5 min | 0.226 | 0.397 / 0.609 / 3.4e-11 | 1.38e-10 (fine; step 352) | stale 0 (when 0005 present) | FAIL only on fine-div step count (29) |
| 3-level | `run_simple_test_l3.sh` | 2000 steps (`t=4.709`) | ~80 min | 0.226 | 1.496 / 0.552 / 1.075 | 0.217 (L1, step 1980); L0 3.61e-11; L2 7.72e-10 | ADV worst Z1/Z2/Z3/rho 0.769/0.0125/0.126/0.646; DIF 0.0394/0.00180/0.0394/1.6e-16; stale 0 | FAIL (species Z1/Z3 and L1 div spikes; 1335 fine-div steps) |

Level-0 numbers are from a matching `fds_amr` binary with `--levels 1` (451 steps to `T_END=5`). 2-level and 3-level numbers are from the Role-3 species gate logs (`r5_wbc`, `r3lev`) on the binary with 0004+0005+0006.
The 3-level run keeps mass at the same 0.226 ratio as 2-level and passes the flux-register stale gate (0), but
regresses on species (Z1/Z3 > 1) and shows late L1 divergence spikes after regrids.

---

## 8. Where to start reviewing

- `Source/driver/` (C++ driver: `main.cpp` flags, `TimeLoop.cpp`, `TwoLevelRun.cpp`), `Source/pressure_backend/` (FFT/MLMG/HYPRE behind one interface),
  `Source/regrid_transport/` (regrid, post-regrid projection).
- Design: [`Manuals/FDS_AMReX_Design/README.md`](Manuals/FDS_AMReX_Design/README.md), ADR-001 (architecture, GPU K2), ADR-002 (time stepping), `pressure/`.
