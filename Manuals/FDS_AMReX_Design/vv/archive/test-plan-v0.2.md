# FDS-AMR Verification & Validation Test Plan (DRAFT v0.2)

Owner: AMR V&V Lead · Status: **draft for review by the project owner**, 2026-09-25 · Nothing in this plan is approved yet.

| Item | Value |
|---|---|
| Refactor base / acceptance baseline source | **FireX `36975d765f`** (`36975d765fcead401e14b094a04f910ac42eab8a`), local branch `AMReX`, this repository (read-only). Recorded with `git -C (repo root) rev-parse --short HEAD`. `git describe` = `FDS-6.11.1-1244-g36975d765f`. Working tree clean. |
| Existing binaries | GNU/OpenMPI and Intel/Intel MPI builds of master **`ce1f659`** under `(local FDS master checkout)/Build/…`. **Smoke tests only** (toolchain, MPI, comparison tooling). They can't start on the development machine today (missing MPI/oneAPI runtimes). Never used as references. |
| Requirements read | `docs/requirements.md` **v0.2**, revised by the Spec & Program Lead during this work (baseline re-pinned to FireX, FR-038 and FR-074 added). This plan follows v0.2. requirements.md was **not** edited. |
| Companion files | `environment.md` (machine, binaries, smoke results, Smokeview, toolchain facts), `case_inventory.md/.csv` (every candidate run with tier, estimate and class), `verification_case_survey.csv` (all 941 FireX inputs), tools in `(local V&V run directory)/tools/` |

Conventions: **T0/T1/T2/T3** are the tolerance classes of requirements.md §2.2. **Tier 1/2/3** are case-set sizes (run frequency) and are unrelated to T-classes. Runtimes are model estimates (±3×) until measured. "Derived" inputs (`VV/…`) are V&V copies of committed inputs with one documented change. They don't exist yet and are created under `vv-runs/`, never in the worktree.

---

## 1. Purpose and scope

The refactor moves FDS `MESH` data onto AMReX grids, first in *uniform mode* (zero refinement levels) and then with static and dynamic refinement. V&V has to show, phase by phase (roadmap.md), that:
1. in uniform mode, results match baseline FDS to the class each requirement names (FR-001/002/003/037/038/061, IR-001/006);
2. with refinement, the code conserves (FR-020..024), reproduces FDS's existing coarse-to-fine coupling where that is the stated reference (FR-016), and is more accurate than the coarse grid and approaches the uniform-fine grid (FR-014, NFR-032);
3. restart, determinism, parallelism, output and performance requirements hold (FR-015/070-074/080, NFR-010..034).

This plan covers verification (code against analytic solutions, conservation, and baseline FDS). Validation against experiments is out of scope until compute beyond the 8-core development machine is available (charter Q6).

## 2. Tolerance classes as V&V will apply them

| Class | requirements.md definition (v0.2) | How V&V evaluates it | Tool |
|---|---|---|---|
| **T0 bitwise** | Byte-identical `CHID_devc/hrr/mass.csv`, timing/version columns excluded, same ranks, `OMP_NUM_THREADS=1` | Same, on V&V copies of the input with `&DUMP SIG_FIGS=17` added so the CSV text round-trips doubles exactly (see D-1). Timing columns are removed by name pattern plus a per-case list (e.g. DEVC `cpu` in `dancing_eddies_*`). | `compare_csv.py --class T0` |
| **T1 tight** | Every column, every output time: \|x−x_base\| ≤ 1e-10·max(‖x_base‖∞,1e-30) | Same, per column, no interpolation (output times must match). Also on SIG_FIGS=17 copies, because 8 significant figures can't resolve 1e-10. | `compare_csv.py --class T1` (`--t1-factor`) |
| **T2 verification** | Meets the FDS verification acceptance metric; error metric no worse than baseline by >5% relative (proposed) | V&V proposal in §9 (R-1): FDS criterion **and** `e ≤ (1+m)·e_base + 0.1·Tol_FDS`, with m = 5% by default and a calibrated m for chaotic cases. Plot-only rows get a V&V numeric metric (observed order, L2 error). | dataplot metrics re-implemented in numpy (no pandas/matplotlib on the development machine), `compare_csv.py` for series |
| **T3 physical** | Case-specific, multi-level vs uniform-fine | Proposed values in §9 (R-2), calibrated at Phase 4 | per-case scripts (to write) |

Which requirement uses which class, and which class each case is judged at: §7 and the `tclass` column of `case_inventory.csv`.

## 3. Baselines, reference results and the FireX build request

### 3.1 What exists today
- **Smoke-test binaries (ce1f659).** `fds_ompi_gnu_linux` (GNU Fortran 14.2.0, OpenMPI soname .40, OpenMP, HYPRE 3.0.0, SUNDIALS 7.5.0, no MKL) and `fds_impi_intel_linux` (ifx 2026.1.1, Intel MPI 2021.18, MKL static, no OpenMP). Details are in environment.md §2. They are used **only** to prove the toolchain, MPI launch and comparison tooling (G0 below).
- **Smoke status (corrected 2026-09-25):**
  - GNU: **PASS** (version print). From a clean shell, `ldd` resolves `libmpi_usempif08.so.40` to `/lib/x86_64-linux-gnu` (Debian `libopenmpi40`/`openmpi-bin` 5.0.7), and `/usr/bin/mpirun -np 1 fds_ompi_gnu_linux` prints the version banner (GCC 14.2.0, Open MPI 5.0.7, HYPRE 3.0.0, SUNDIALS 7.5.0). The earlier FAIL came from a non-clean shell environment and is withdrawn.
  - Intel: FAIL. `libimf.so` exists nowhere on disk and there is no `/opt/intel/oneapi/setvars.sh` or `mpiexec`. The oneAPI install is gone (the development machine rebooted Sep 23). Reinstalling oneAPI is a prerequisite for the Intel FireX build.

### 3.2 Where reference results must come from
**Reference (baseline) results for every acceptance test must come from FDS built from FireX `36975d765f`. That build does not exist yet.** requirements.md §2.1 (v0.2) defines baseline as a CMake, out-of-tree FireX build with the HYPRE GPU options OFF, using the same compiler, MPI, flags, ranks and threads as the build under test. It also says results from the earlier `ce1f659` baseline are not acceptance references. The `ce1f659` binaries can't serve even as a proxy: FireX changes radiation, init, pressure (GPU HYPRE paths, resource sets) and output code (`case_inventory.md` §1a). **A build request will go to the GNU Build Chief and the Intel Build Chief.** V&V does not build (A-08 records the build; it does not perform it). No reference capture happens until both builds exist and pass G0.

### 3.3 Build request spec (to the GNU and Intel Build Chiefs)

| # | Item | GNU build | Intel build |
|---|---|---|---|
| 1 | Source | this repository at `36975d765fcead401e14b094a04f910ac42eab8a` (branch `AMReX`). Check `git rev-parse HEAD` and an empty `git status --short` before and after. **The worktree is read-only**: nothing may be written into it. | same |
| 2 | Build system (primary, = baseline per requirements §2.1) | CMake ≥ 3.24, **out of tree**: `cmake -S (repo root) -B (local build directory)/firex-36975d7/ompi_gnu_rel -DCMAKE_BUILD_TYPE=Release`. Take the compiler wrappers from preset `ompi_gnu` (`CC=mpicc CXX=mpic++ FC=mpifort`). **Don't use the preset `binaryDir`** (`./Build/cmakeb/<preset>`): it is inside the worktree. | same, preset `impi_intel` wrappers (`mpiicx`/`mpiifx`), build dir `…/impi_intel_rel` |
| 3 | Compiler / MPI | Same as the existing binary: GNU Fortran 14.2.0 (Debian 14.2.0-19) + OpenMPI 5.0.x (soname .40), as NFR-020 requires | Same as the existing binary: ifx 2026.1.1 + Intel MPI 2021.18, with MKL visible to CMake (`find_package(MKL CONFIG)` → `WITH_MKL`) so the `*_pardiso` cases work as they do in the existing Intel binary |
| 4 | Options | `USE_HYPRE=ON`, `USE_HYPRE_NVIDIA/AMDGPU/INTELGPU=OFF`, `USE_SUNDIALS=ON`, `USE_OPENMP=ON` (CMake default). Record any deviation. | same. Note: the existing Intel makefile binary had **no** OpenMP; CMake enables it by default. Keep the CMake default and record it. |
| 5 | HYPRE | FireX CMake fetches hypre commit `63331f19` (requires ≥ 2.32.0) over the network (A-06/NFR-022). Alternative: `USE_SYSTEM_HYPRE=ON` with the existing v3.0.0 libs (`(local GNU third-party library tree)/libs/hypre`). Build Chief chooses; **record which**, because it sets the HYPRE version for every HYPRE reference (FR-016, FR-038). | same (`(local Intel third-party library tree)/libs/hypre`) |
| 6 | SUNDIALS | v7.5.0 (fetched, or the existing libs); record | same |
| 7 | HDF5 / VTK | Not available in CMake (no option). Leave out of the primary build; record "HDF5: no". | same |
| 8 | Provenance (A-08) | `fds -V` must print a non-empty revision (CMake takes `git describe`; if building from a copy, pass `-DGIT_HASH=…`). Deliver: `CMakeCache.txt`, verbose compile lines (actual flags), `mpifort --version` / `mpifort -showme`, HYPRE/SUNDIALS versions, `ldd` output, sha256 of the binary, and the environment script that makes it run **on the development machine** | same, plus `mpiifx -show`, MKL version |
| 9 | Runnable on the development machine, reboot-proof | OpenMPI 5.0.7 runtime is present (Debian packages). Write down the environment setup (a script plus notes) so it survives the next reboot. | **Prerequisite:** reinstall oneAPI (ifx 2026.1.1, Intel MPI 2021.18, MKL) under `/opt/intel/oneapi` with `setvars.sh`; it is missing since the Sep 23 reboot. Write down the setup so it survives the next reboot. |
| 10 | Debug variant (requested, lower priority) | Same configuration, `CMAKE_BUILD_TYPE=Debug` (preset `ompi_gnu_db`), for A-09 NaN/bounds checks on the three AMR cases | `impi_intel_db` |
| 11 | Optional: HDF5 variant for FR-074 / A-14 | Makefile target **`ompi_gnu_linux`** (`Build/ompi_gnu_linux/make_fds.sh`) with HDF5 1.14.5 (`build_thirdparty_libs.sh` builds it only if `$FIREMODELS/hdf5` is a clone at tag `hdf5_1.14.5`; otherwise it silently omits HDF5). Makefile builds run **inside** `Build/<target>`, so they must use an exported copy of the commit (e.g. `git archive 36975d765f`) outside the worktree. | target **`impi_intel_linux`** (no OpenMP; `-O2 -ipo`) |
| 12 | Later (Phase 3): instrumented baseline for FR-016(a) | Same commit plus a V&V ghost-cell dump patch that doesn't change arithmetic. Accepted only if it is T0 against item 2 on all Tier 1 CSVs. Spec to follow with the Pressure and Integration Leads. | same |

Build names for reference: makefile targets `ompi_gnu_linux`, `impi_intel_linux` (+ `_db`, `_dv`, `impi_intel_linux_openmp`); CMake presets `ompi_gnu_rel`, `impi_intel_rel` (+ `_db`, `_dv`).

**V&V acceptance of a delivered build (G0):** `fds -V` shows `36975d765f`. `ns2d_16` runs on 1 rank and `obst_activation_default` on 4 ranks, each finishing normally in < 2 min. A setup-only sweep (`T_END=0` copies) over all 941 inputs records which ones stop at setup; that list becomes the IR-001 reference.

### 3.4 Reference capture procedure (after the builds pass G0)
1. Order: A-09 runs first (three AMR cases + the UGLMAT copy, GNU and Intel, plus the Debug build), then Tier 1, then Tier 2.
2. Every run: `OMP_NUM_THREADS=1` unless the case is a threading case. Load < 2 and free RAM ≥ ranks×0.5 GB + 1 GB checked before start. `mpirun --bind-to none` (shared machine). Watchdog per §9 (R-10).
3. Every Tier 1 case runs twice: once as committed, and once as a SIG_FIGS=17 copy for T0/T1 use. Chaotic cases (fire/LES) and the FR-002 cases also get the extra variants their class needs (GLMAT copy, alternate rank count).
4. Stored under `(local V&V run directory)/baseline/<toolchain>/<case>/` with a manifest: input sha256, binary sha256, `fds -V`, env, ranks/threads, launcher line, wall time, exit status, load/free-RAM at start. Outputs are read-only after capture.
5. T2 margin calibration ensemble (for chaotic cases): GNU vs Intel and firebot rank count vs one alternate rank count, giving σ for R-1.
6. Measured wall times replace the model estimates, and cases are re-tiered.

## 4. Environment and execution rules
- Development machine: 8 cores (1 socket, AVX-512), 16 GB RAM with **no swap**, shared with other workloads (2.4–6.4 GB free observed). At most 8 ranks per case and at most 8 ranks in total across concurrent cases. No case above 12 GB (NFR-031).
- Python 3.13 + numpy only (no pandas, scipy or matplotlib). FDS's `Utilities/Python` plotting scripts can't run here as-is, so V&V re-implements the needed metrics in numpy (§8).
- Timing tests (NFR-030/033/034) are valid only on an otherwise idle machine: load < 1, ≥ 8 GB free, no other V&V jobs. Otherwise the run is flagged invalid rather than failed.
- Smokeview: not installed. Proposed reference version SMV 6.11.2 (§9, R-9).

## 5. Test gates

| Gate | When | Content | Class / criterion |
|---|---|---|---|
| **G0** Toolchain smoke | Each new binary (incl. ce1f659 once runtimes exist) | `fds -V`; `ns2d_16` (1 rank); `obst_activation_default` (4 ranks); `compare_csv.py` self-test (identical, perturbed and re-timed synthetic CSVs) | Runs / exit 0. No T-class (ce1f659 output is never a reference). |
| **G1** Input compatibility (IR-001) | Phase 2, then each milestone | Setup-only (`T_END=0`) copies of all 941 inputs | Same set of setup stops and messages as baseline (message diff) |
| **G2** Tier 1 regression | Every merge to `AMReX` (≈ 8 min model, budget ≤ 30 min) | 40 runs (§6) | T1 single-mesh (T0 stretch, Q1); T2 vs FFT baseline plus T1 vs GLMAT copy for same-resolution multi-mesh (FR-002); T1 for the first 100 steps plus T2 for chaotic cases; T0/T1 for restart; T0 for determinism |
| **G3** Tier 2 | Nightly while code changes; every phase exit | 123 runs (≈ 79 min sequential, ≈ 30–40 min packed on 8 cores) + optional 8 | T2 (FDS criterion + margin, R-1) |
| **G4** Suite non-regression (FR-003) | Phase 2 exit, Phase 10 | All firebot-listed cases with ≤ 8 ranks (870 of 884) | T2. Exceptions recorded and approved by the project owner |
| **G5** Build-option equivalence (IR-006) | Phase 2 onward | FR-001 anchors with `USE_AMREX=OFF` | **T0** |
| **G6** FireX HYPRE options (FR-038) | Phase 2 | The 4 `Pressure_Solver/*_hypre.fds` at 1 and 4 ranks, `FDS_RANKS_PER_GPU` unset and =2 (the variable changes the matrix partition on CPU too: RS masters gather and solve, pres.f90:3411-3440) | T1 vs baseline at the same setting |
| **G7** Static two-level (FR-016) | Phase 3 (a–c), Phase 4 (d) | `ns2d_16_int_1to2_refinement` + `VV/…_uglmat` (+ `emb_1to2` for (a) only, D-4) | (a) T0 on first-step fine-side ghosts; (c) Pressure Lead rule (H minus mean, rel. L2 < 10× solver tol; FireX HYPRE tol 1e-12, A-16); (d) T2 vs UGLMAT-HYPRE baseline |
| **G8** Conservation (FR-020/021/022/024) | Phase 3/4 | `species_conservation_1..4`, `energy_budget_*`, `simple_duct`, mass-balance cases with a refined patch | Requirement numbers (1e-12/step, 1e-10 cumulative; FR-022 with floor, R-5) |
| **G9** Accuracy (FR-014, NFR-032) | Phase 4, 10 | `ns2d_{8..64}` + static patch, `ns2d_16_*_refinement`, plume case (R-8) | T3 (R-2), order ≥ 1.8, factor R-4 |
| **G10** Determinism / restart (FR-015, FR-080/081) | Phase 3, 9 | repeat runs, hierarchy dumps; `restart_test1a/b` + continuous copy, `device_restart_*`, `restart_ulmat_*` | T0 (hierarchy identical); restart per R-11 |
| **G11** Parallel / threads (NFR-010/011/012) | Phase 2+ | anchors at 1/2/4/8 ranks where meshes allow (D-8); `race_test_1/4` | T2 across ranks in uniform mode, T3 in AMR mode; watchdog (R-10) |
| **G12** Performance / memory (NFR-030/031/033/034) | Phase 2 (measure), 10 (meet) | `openmp_test64a`; one 8-rank multi-mesh anchor; scaling case (R-12) | Timing ratios, median of 3, idle machine only |
| **G13** Output (FR-071/072/074) | Phase 9 | header diff; Smokeview load; VTK (needs HDF5 build) | Exact header match (I); D for Smokeview/ParaView |

Negative tests (FR-004, FR-044, IR-002): one input per excluded feature, expecting `SETUP_STOP` and the named `ERROR(nnn)`. Written in Phase 3 as the supported feature set is defined.

## 6. Case selection (details and per-run estimates in `case_inventory.md`)

**Tier 1 (40 runs, ≈ 8 min sequential model estimate).**
- Interfaces/pressure: `obst_activation_{default,ulmat}`, `divergence_test_2/3`, `dancing_eddies_{1mesh,default,uglmat_refine}`, `obst_coarse_fine_interface`, `lapse_rate`.
- Refinement analogues: `random_meshes`, **`ns2d_16_int_1to2_refinement` + derived `VV/…_uglmat`**, **`ns2d_16_emb_1to2_refinement`**.
- Interface invariance: `soborot_superbee_square_wave_128{,_1mesh}`, `shunn3_4mesh_128`.
- Convergence: `ns2d_{8..64}{,_nupt1}`, `saad_512_cfl_*`, `shunn3_{32,64,128}`.
- Conservation: `energy_budget_tmix`, `species_conservation_1/2`.
- Fire: `1_step_2_step_compare`.
- Restart: `restart_test1a/b` + derived continuous run.
- Determinism pair.

**Tier 2 (123 runs, ≈ 79 min sequential / ≈ 30–40 min packed)** adds the pressure-solver variants, the remaining refinement analogues (`dancing_eddies_embed`, `duct_flow_uglmat_refine`, `ns2d_16_emb_1to1_refinement`, `porous_media`, `race_test_1`), the full MMS / scalar / limiter convergence series, multi-mesh symmetry and volume-flow cases, pressure zones, mass-balance cases, fires, restarts, and every anchor case not in Tier 1 (`energy_budget_particles`, `zone_break_fast`, `zone_shape`, `cascadempi`, `cannon_ball`, the two radiation anchors, `openmp_test64a`). **Optional (8):** `tunnel_demo`, `energy_budget_adiabatic_walls`, `race_test_4`, `shunn3_512`, `ht3d_energy_conservation_4`, `VV/test_8mesh_NOFRPG_short` (HYPRE cost), and the Intel-only `*_uglmat_pardiso` pair.

**Adaptive_Mesh_Refinement folder (scope item; A-09).** `ns2d_16_int_1to2_refinement` (13 meshes / 448 cells, fine 16×16 patch abutting 12 coarse 4×4 meshes, ratio 2) and `ns2d_16_emb_1to2_refinement` (2 meshes / 512 cells, fine mesh embedded in the coarse) are **Tier 1**. `ns2d_16_emb_1to1_refinement` (2 / 320, same-resolution control) is **Tier 2**. All are 2-D, periodic, T_END = 30 s, with `VELOCITY_TOLERANCE=1e-6` and up to 100 pressure iterations. Estimates: ~10–30 s each; the 13-mesh case is 0.5–5 min on 1 rank, depending on iteration counts. Memory is negligible.

*Why commented out:* all three lines (FDS_Cases.sh 4-6) come from commit `042aa2624d` (R. McDermott, 2024-11-08), *"updates to amr test cases, demonstrating errors in both embedded mesh strategy and interpolated refinement"*. That commit also created the `int` case. They were added deliberately disabled, as demonstrations of known error. They are not broken or slow. There is no dataplot row or Guide section, and nothing has changed since. A-09's "pass" therefore needs a V&V definition: runs to T_END with no NaN on both toolchains, and the analytic-error values are recorded as the baseline values (`case_inventory.md` §4a).

**FireX-only folders.**
- **`GPU_Tests`**: 9 inputs, all 1 M cells, propane fire, UGLMAT/HYPRE, plus HPC job scripts. The physics runs on CPU (`HYPRE_DEVICE_RUN` is compiled out without `WITH_HYPRE_DEVICE`). The job scripts don't apply. The 16/32-mesh variants need more than 8 ranks. Full runs take 20–45 min each, so Tier 3, except one 20-step copy (T2-opt). **Highly relevant to the HYPRE pressure work:** the only ≥1 M-cell multi-mesh UGLMAT/HYPRE inputs, and the natural source for the FR-016(c)/FR-031/A-16 tolerance, NFR-033 scaling and FR-038 resource-set checks.
- **`VTK`**: 4 CPU inputs demonstrating the VTKHDF writer. Their data output needs an HDF5 build, which CMake can't produce. They turn off `.smv` output and have no pass/fail criterion, so they are Demonstration only. Not relevant to the solver, but relevant to FR-072/074 and Q5.

Full assessment: `case_inventory.md` §6.

**Excluded, with reasons:** `case_inventory.md` §7 (GEOM, >8 ranks, hour-plus runs, local physics sub-models covered by G4).

## 7. Requirement → test traceability (class used)

| Req | Test / cases | Class | Phase |
|---|---|---|---|
| FR-001 | single-mesh anchors (ns2d_16/32, shunn3_32, species_conservation_1, energy_budget_*, cannon_ball, radiation anchors, openmp_test64a, restart_test1*) | T1 (T0 stretch); chaotic ones T1 over 100 steps + T2 (D-2) | 2 |
| FR-002 | multi-mesh same-res anchors (shunn3_4mesh_32, layer_4mesh, dancing_eddies_default, symmetry_test_mpi, duct_flow, …) | T1 vs GLMAT copy (see D-3) or T2 vs 1-mesh; T2 vs FFT baseline | 2 |
| FR-003 | G4 | T2 | 2, 10 |
| FR-004, FR-044, IR-002 | negative tests | SETUP_STOP + message | 3 |
| FR-010/013/015 | hierarchy dumps, debug assertions, repeat runs | exact box lists (T0-equivalent) | 3 |
| FR-011 | per-criterion tagging unit cases | exact cell sets | 3 |
| FR-012 | regrid mass/species deltas | 1e-12 rel (req.) | 3 |
| FR-014 | ns2d series + patch, ns2d_16_*_refinement | T3 + factor R-4 + order ≥ 1.8 | 4 |
| FR-016 | G7 | T0 (a), Pressure-Lead rule (c), T2 (d) | 3–4 |
| FR-020/021/024 | species_conservation_1..4, simple_duct, mass_balance_* with patch | req. numbers | 3 |
| FR-022 | energy_budget_tmix/particles/dns_100/adiabatic_walls with patch | R-5 | 4–5 |
| FR-023, FR-031, FR-032 | divergence_test_1..3, obst_activation, CHECK_POISSON output | req. numbers (Pressure Lead) | 4 |
| FR-030/033/035 | duct_flow, tunnel_demo, dancing_eddies_default, ns2d periodic; pressure_iteration3d_default, random_obstructions_fft | T2 | 4–5 |
| FR-034 | zone_break_fast, zone_shape (+ zone_shape_2, zone_break_fast_uglmat_hypre) | T2 | 4/6 |
| FR-036 | shunn3_128 (PRESSURE_TOLERANCE 1e-6, analytic H) + helium_2d_isothermal (R-6) | T2 | 4 |
| FR-037 | FR-001 cases with FFT selected | D + T1 | 2 |
| FR-038 | G6 | T1 | 2 |
| FR-040/041/042/043 | mask checker; burning-OBST regrid case; obst_activation_default, box_burn_away1; random_meshes, duct_flow | req. / T2 | 5 |
| FR-050/051/052 | particle ledger; energy_budget_particles, cascadempi; FR-052 case R-7 | exact counts / T2 | 7 |
| FR-060/061 | plate_view_factor_cart_30, radiating_polygon_square_20 | T2 (060); T1 (061 via FR-001) | 2, 8 |
| FR-070 | all Tier 1 DEVC outputs; regrid-jump check | T1 uniform; T3 jump | 9 |
| FR-071 | header diff on all Tier 1/2 CSVs | exact (I) | 9 |
| FR-072 | Smokeview load of Tier 1 cases | D, SMV 6.11.2 (R-9) | 9 |
| FR-074 | VTK/*.fds (HDF5 build, A-14) | identical VTKHDF (D-7) | 9 |
| FR-080/081 | restart_test1a/b + continuous, device_restart_*, restart_ulmat_*, csvf_restart_a | R-11; T2 (081) | 9 |
| IR-001 | G1 | message diff | 2 |
| IR-004 | random_meshes, layer_4mesh (MPI_PROCESS inputs) | behaviour + warning text | 3 |
| IR-006 | G5 | T0 | 2 |
| NFR-010/011/012 | G11 | T2 uniform / T3 AMR; watchdog R-10 | 2+ |
| NFR-030/031 | openmp_test64a + 8-rank multi-mesh anchor (duct_flow_uglmat_refine proposed) | ratio, median of 3 | 2, 10 |
| NFR-032 | plume case R-8 | T3 + ≤ 50% wall | 10 |
| NFR-033 | scaling case R-12 | ≥ 60% efficiency (proposed) | 10 |
| NFR-040 | harness over the anchor set using T0–T3 | D | 2 |

## 8. Tooling (delivered and planned)
- **Delivered:**
  - `parse_fds_inputs.py`: static survey. Follows FDS `CHECKREAD` rules and expands MULT including SKIP ranges.
  - `make_inventory.py` and `render_inventory_md.py`: inventory, cost model, class mapping.
  - `compare_csv.py`: CSV comparison with `--class T0|T1` presets. Implements §2 T0/T1 exactly: timing-column exclusion patterns, T1 column-max scaling, `--exact`, per-column tolerances, interpolation for T2-style series, JSON reports. Tested on synthetic CSVs only (`vv-runs/tools/testdata`), since no real FDS output exists yet.
- **Planned (Phase 2, NFR-040):**
  - `run_case.sh` harness: pre-flight load/RAM check, watchdog, manifest.
  - numpy re-implementations of the dataplot metrics (`end`, `max`, `mean`, `area`, `end_1_n`, `all`, `tolerance`, `slope`) and of the analytic-error scripts needed by Tier 1/2 (`ns2d.py` RMS error, `shunn_mms.py`, `saad_mms_temporal_error.py`, `soborot_mass_transport.py`, `mass_balance*.py`).
  - A T3 statistics script (time-window mean/RMS at probes).
  - A hierarchy-dump comparator once the dump format exists.

---

## 9. Responses to requirements.md TBD(V&V)

All numbers below are **proposals, to be calibrated** against the FireX baseline once it exists (A-08/A-10). None is a commitment until the project owner approves.

### 9.1 Answers

| # | TBD(V&V) item (requirements.md v0.2) | V&V answer |
|---|---|---|
| **R-0** | §2.3: does each anchor case run on the 8-core development machine in reasonable time? | **Yes, all 37 anchor inputs fit**: max 8 ranks, < 0.4 GB each. Total ≈ 53 min sequential model estimate (±3×). Most take under 1 min. Longest *(est.)*: `energy_budget_adiabatic_walls` ≈ 14 min (1 mesh, 1 rank, ~5.3 k steps); `tunnel_demo` ≈ 9 min (8 ranks, whole box); `race_test_1` and `race_test_4` ≈ 5 min each; `openmp_test64a` ≈ 5 min (×3 for NFR-030 medians); `random_obstructions_fft` ≈ 3 min. **Rank caveats:** `random_obstructions_fft` has 16 meshes (firebot `-p 16`), so run it on 1 rank or as an 8-rank `MPI_PROCESS` copy. `ns2d_16_int_1to2_refinement` has 13 meshes: 1 rank, or `-p 13` oversubscribed (acceptable for 448 cells). `race_test_4` needs 1 rank × 4 threads. **Status caveats:** the two `ns2d_16_*_1to2` anchors are commented out of FDS_Cases.sh on purpose (A-09, §6). Pair references not in the anchor list are needed: `device_restart_base_case` (dataplot reference for `device_restart_a/b`) and a derived continuous run for `restart_test1a/b`. Per-case table: `case_inventory.csv` (`anchor=yes`). |
| **R-1** | §2.2 T2 margin (5% relative, proposed) | Keep 5% as the default, with two additions: **(i) absolute floor:** pass if `e ≤ (1+m)·e_base + 0.1·Tol_FDS` (Tol_FDS = the dataplot Error_Tolerance or the script limit), so near-zero metrics (e.g. divergence max ≈ 1e-15 against a 1e-13 tolerance) aren't judged on round-off. **(ii) chaotic cases** (fires/LES): `m = max(5%, 3σ)`, where σ is the relative spread of the metric over the baseline ensemble {GNU, Intel} × {firebot ranks, alternate ranks} (§3.4 step 5). The FDS criterion itself (`e ≤ Tol_FDS`) must always hold. **(iii) plot-only / convergent-series rows** (`N/A`, `Convergent Series`; e.g. the ns2d, pulsating, vort2d series): V&V metric is observed order `p ≥ min(p_base, nominal) − 0.1` and finest-grid L2 error per (i). |
| **R-2** | §2.2 T3 values | Proposed T3 for multi-level vs uniform-fine (and AMR-mode rank comparisons). **Analytic cases:** refined-region L2 error ≤ 2× uniform-fine (R-4) and global L2 ≤ coarse. **Time-averaged plume quantities** (window ≥ last 50% of run, ≥ 10 puffing periods): centreline ΔT and w within **10%** of uniform-fine at each probe height inside the refined region. **Integrated HRR / burner MLR:** mean within **2%** (1% where HRR is prescribed). **Compartment layer height:** within **5% or one fine cell**; layer temperatures within **5%**. **Turbulent point probes:** mean within **10%**, RMS within **25%**, never instantaneous values. **Discrimination rule (always applies):** \|AMR − fine\| ≤ **0.5·\|coarse − fine\|** for the same quantity; if the coarse and fine runs don't differ enough for this to mean anything, the case isn't a T3 case. **Calibration:** T3 = max(proposal, 2× measured GNU-vs-Intel spread of the uniform-fine run), capped by the discrimination rule. |
| **R-3** | FR-003: full suite or the subset that fits | **Both, at different frequencies.** Tier 1+2 (163 runs, ≈ 90 min sequential / ≈ 40 min packed) at every phase exit. Full firebot-listed suite restricted to ≤ 8 ranks (**870 of 884** active cases) at Phase 2 exit and Phase 10: crude model ≈ 15–30 h sequential, i.e. an overnight-to-weekend batch on 8 cores with the development machine to itself. Seven Heat_Transfer cases dominate (`ht3d_ibeam` alone ≈ 9 h *est.*). The 14 cases needing 9–64 ranks become recorded exceptions (the project owner approval, per FR-003) unless extra compute appears (Q6). `GPU_Tests` and `VTK` aren't in the firebot suite. |
| **R-4** | FR-014: factor on uniform-fine error | **2.0**: refined-region L2 error ≤ 2 × uniform-fine L2 error in the same region. For 2nd order and ratio 2, coarse ≈ 4× fine, so 2× means at least half the log-distance gained. Global L2 ≤ coarse, as written. **Note:** baseline FDS is expected to fail "refined ≤ coarse" on `ns2d_16_*_refinement` (commit 042aa2624d documents this), so FR-014 is demonstrated on the AMR code only; baseline numbers are context (A-09). |
| **R-5** | FR-022: closure ≤ baseline closure × 1.05 | Agree, with an absolute floor: `|closure_AMR| ≤ max(1.05·|closure_base|, 0.005·max_t|Q_TOTAL|)` over the comparison window. Rationale: baseline closure can be close to zero (tmix, dns_100), where a multiplicative rule is ill-conditioned. |
| **R-6** | FR-036: baroclinic case | **`Scalar_Analytical_Solution/shunn3_128.fds`** (variable-density MMS; sets `PRESSURE_TOLERANCE=1e-6`, so `ITERATE_BAROCLINIC_TERM` really iterates; analytic H error available; 16 k cells, ~20 s *est.*) as primary. Secondary: `Flowfields/helium_2d_isothermal.fds` (strong density ratio, 1 mesh, 10 k cells, ~1 min *est.*, dataplot abs 0.01). Gate: T2 plus pressure-iteration counts no higher than baseline +1. |
| **R-7** | FR-052: particle-wall case | Primary **`Sprinklers_and_Sprays/cascadempi.fds`** (droplets land on, run over and drip off box tops across 6 meshes; 57 k cells, ~1 min) with a static, later dynamic, refined patch whose boundary cuts a box top, regridding while droplets are attached. Secondary: a derived **OBST-only copy of `geom_sprk_mass.fds`** (drop the GEOM half), whose accumulated-surface-water output measures wall attachment directly. |
| **R-8** | NFR-032: plume case | **Derived from `Fires/circular_burner.fds`** (8 meshes, 5 cm, propane burner, T_END 20 s, 64 k cells, ~1 min *est.*), with centreline T/w probe columns added. Uniform-fine = 2.5 cm (512 k cells, ≈ 0.6 GB, ≈ 30–40 min at 8 ranks *est.*); coarse = as committed; AMR = 5 cm base + 2.5 cm over the burner/plume. Fits the development machine and uses T3 R-2 plume metrics. Alternative if a HYPRE-heavy case is wanted: `GPU_Tests/test_8mesh_NOFRPG` (8 mm uniform-fine = 1 M cells), but its uniform-fine reference costs hours. |
| **R-9** | FR-072: Smokeview version | Smokeview is **not installed on the development machine**. FireX `36975d765f` = `FDS-6.11.1` + 1244 commits. The matching release bundle is **FDS-6.11.1_SMV-6.11.2** (released 2026-07-10; standalone SMV then at 6.11.1). Proposal: pin **SMV 6.11.2**, falling back to the SMV nightly of the same date if FireX-only output doesn't load. FireX still defaults to `SMOKE3D_VERSION=1`. Installation needs someone with authority; the installed version goes in environment.md §8. |
| **R-10** | NFR-011: wall-clock watchdog | Per case: `timeout = max(10 min, 3 × measured baseline wall time at the same ranks/threads)`, then SIGTERM, a 60 s grace period, and SIGKILL. **Stall detector:** fail as a hang if neither `CHID.out` nor `CHID_devc.csv` has grown for `max(5 min, 20 × baseline median wall time per output interval)`. Hangs are recorded as FAIL with the last `.out` lines and rank states. Before baselines exist: 10 min for Tier 1 cases, 60 min for Tier 2. |
| **R-11** | (related, FR-080 "T0 proposed; T1 minimum") | V&V wants to measure first. FDS's own `device_restart` check allows **5 °C** absolute (dataplot), which hints that baseline restart may not be bitwise. Proposal: uniform-mode restart must be **no worse than baseline's own restart-vs-continuous difference**, and T0 if baseline achieves T0. Measure on `restart_test1a/b` + continuous copy and `device_restart_*` at capture. |
| **R-12** | NFR-033: scaling case | **`GPU_Tests/HYPRE_GPU_SCALING/test_8mesh_NOFRPG.fds`**, shortened to 20 steps, 8 fixed meshes of 125 k cells mapped by `MPI_PROCESS` copies onto 1/2/4/8 ranks. The decomposition stays fixed and only the rank mapping changes: a clean strong-scaling test with the HYPRE solve inside (≈ 12 min at 1 rank, ≈ 2 min at 8 *est.*). Secondary, non-reacting: `duct_flow_uglmat_refine` (147 k cells, 8 meshes). The 60% target is meaningful only on an idle machine (§4). |

### 9.2 Disagreements and clarifications (tolerances and cases)

| # | Item | Concern | V&V proposal |
|---|---|---|---|
| **D-1** | §2.2 T0 definition | FDS writes CSVs with `SIG_FIGS=8` by default (read.f90:2374). "Byte-identical CSV" at 8 digits isn't bitwise: differences up to ~5e-9 relative can vanish in rounding, while a 1-ulp difference near a rounding boundary flips a digit. At 8 digits T0 is neither stronger nor weaker than T1 (1e-10). | Evaluate T0 (and T1) on V&V copies with `&DUMP SIG_FIGS=17` (output only; numerics unchanged). Keep an explicit per-case list of timing columns to exclude (e.g. DEVC `QUANTITY='CPU TIME'`, ID `cpu`). |
| **D-2** | FR-001 "all single-mesh anchors meet T1" | Several anchors are chaotic LES/fires (`energy_budget_adiabatic_walls`, `restart_test1a/b`, `openmp_test64a`). Any change in floating-point operation order grows exponentially there. Full-run T1 is then effectively T0 and forces arithmetic-order preservation. That is a Q1 decision, not a tolerance detail. | For chaotic cases: T1 over the first 100 steps + T2 over the full run, unless Q1 chooses bitwise-preserving refactoring. Keep full-run T1 for laminar/DNS/analytic anchors. |
| **D-3** | FR-002 "T1 against baseline GLMAT" | The composite AMReX solve (FR-031 proposes 1e-10 relative) and GLMAT (HYPRE PCG, `HYPRE_SOLVER_TOL=1e-12`, or PARDISO) agree only to solver tolerance per solve. Over hundreds of steps the difference easily exceeds 1e-10 of a column's max. | T1 on the first step (or first 10 steps) against GLMAT, then `|Δ| ≤ max(1e-10, 100·tol_solver)·‖x_base‖∞` over the run, or T2. Needs the Pressure Lead's FR-031 number. |
| **D-4** | FR-016 case list includes `ns2d_16_emb_1to2_refinement` | That input is **not abutting**: the fine mesh overlaps the coarse one. FDS embedded meshes are one-way coupled ("the larger mesh receives no information from the mesh embedded within", User Guide ≈ line 1112), and UGLMAT "is used for non-overlapping meshes" (User Guide 9455). So the FR-016 reference (baseline UGLMAT-HYPRE) can't be run on it, and (b)–(d) don't apply. | FR-016 = `ns2d_16_int_1to2_refinement` for (a)–(d). Use `emb_1to2` for FR-016(a) at most and for FR-014, with `emb_1to1` as the control. |
| **D-5** | FR-016 reference availability | The committed `int_1to2` input uses FFT + iteration, not UGLMAT. The UGLMAT-HYPRE reference is a V&V-derived copy (`SOLVER='UGLMAT HYPRE'`). pres.f90 handles refinement and periodic seams (comment at pres.f90:3245-3246), but whether it runs on this 13-mesh periodic 2:1 case is unknown. FR-016(a) also needs ghost-cell dumps that baseline FDS doesn't write. | Add "A-09b: run the UGLMAT copy on baseline". Request the instrumented-baseline build variant (§3.3 item 12), accepted only if T0 against the plain baseline. |
| **D-6** | §2.3 note "The GPU cases cannot run on the development machine (D-004)" | Only the GPU offload can't run. With a CPU build, `HYPRE_DEVICE_RUN` is compiled out and every `GPU_Tests` input runs on CPU (the 16/32-mesh ones need >8 ranks or an MPI_PROCESS copy). They are the suite's only 1 M-cell UGLMAT/HYPRE multi-mesh inputs. | Reword to "GPU execution can't be tested; the inputs run on CPU". Use them for NFR-033 (R-12), FR-038 (resource sets) and A-16 (HYPRE tolerance). |
| **D-7** | §2.1 baseline = CMake vs FR-074 "identical VTKHDF … in a build with HDF5" | FireX CMake has **no HDF5 option**. Only the makefile wires `WITH_HDF5`, and makefile targets differ in flags (Intel `-O2 -ipo`, no OpenMP) from CMake Release. So the §2.1 baseline can't produce VTK data. | A-14 decides. Either a secondary makefile/HDF5 baseline, compared like-for-like only with a makefile-built test binary (§3.3 item 11), or a CMake HDF5 option added in the `AMReX` branch, with the baseline built from that branch's option-off state (which conflicts with "baseline = unmodified FireX"). |
| **D-8** | NFR-010 / NFR-011 "all anchor cases at 8 ranks", "1, 2, 4 and 8 ranks" | In uniform mode FDS assigns whole meshes to ranks, so a 1-mesh input can't use more than 1 rank. 15 of the 37 anchor inputs (including the named 1-mesh partners) have 1 mesh, 10 have 2–4, and 30 of 37 have fewer than 8. | "at min(8, n_meshes) ranks". Where a rank sweep is really wanted, use V&V MULT-split copies, which change the problem and so compare at T2 against the 1-mesh original. |
| **D-9** | §2.2 T2 "5% relative" | Ill-conditioned when e_base ≈ 0 and too tight for chaotic metrics, where run-to-run spread from rank count or compiler alone can exceed 5%. | R-1 (floor + calibrated σ). |
| **D-10** | FR-022 "× (1 + 0.05)" | Same near-zero problem. | R-5 floor. |
| **D-11** | FR-080 uniform T0 | Possibly stricter than baseline itself (R-11). | Measure, then set. |
| **D-12** | NFR-030/033 timing | The development machine is shared: other jobs used 9.5–13 GB of RAM during this session. A median of 3 runs doesn't protect against that. | Valid only on an idle machine (§4). Otherwise the run is invalid, not failed. |
| **D-13** | Anchor list completeness | `device_restart_a/b` need `device_restart_base_case`; `restart_test1a/b` need a continuous companion; `layer_4mesh`, `shunn3_4mesh_32` and `dancing_eddies_default` need their 1-mesh partners (already named). | Add `device_restart_base_case` to §2.3; V&V supplies the derived continuous run. |

## 10. Open questions and actions for V&V
- **Blocking:**
  - (1) MPI/oneAPI runtimes on the development machine (environment.md §7).
  - (2) FireX builds per §3.3, from the GNU and Intel Build Chiefs.
  - (3) Smokeview install (R-9).
- **Needs a decision from others:**
  - Q1 (T0 vs T1 → D-2).
  - Q5 and A-14 (VTK/HDF5 → D-7).
  - A-16 and FR-031 (solver tolerance → D-3, G7).
  - Pressure Lead's FR-032 number.
  - Q9 (MPI_PROCESS in AMR mode).
- **V&V next steps once builds land:** G0 → G1 sweep → A-09 (+A-09b) → Tier 1 capture on both toolchains → Tier 2 → calibration ensemble (R-1/R-2) → replace model estimates with measurements and re-tier.

## 11. Files
- `docs/vv/test-plan.md` (this file)
- `docs/vv/environment.md`
- `docs/vv/case_inventory.md` and `case_inventory.csv`
- `docs/vv/verification_case_survey.csv` (FireX) and `verification_case_survey_ce1f659.csv` (v0.1)
- `vv-runs/tools/{parse_fds_inputs.py, make_inventory.py, render_inventory_md.py, compare_csv.py, testdata/}`
- `vv-runs/smoke/{run_smoke.sh, inputs/, gnu/, intel/}`
