# FDS-AMR Verification & Validation Test Plan (DRAFT v0.4.5)

Owner: AMR V&V Lead · Status: **draft for review by the project owner**, 2026-09-25 · Nothing in this plan is approved yet. v0.3 brings the plan in line with `requirements.md` v0.3.1: v0.3, plus the bit-parity scope D-022/D-023 and FR-005, published while this revision was being written. It does not re-argue items the Spec Lead has decided (§9). The v0.2 text, including the full §9 rationale, is kept at `docs/vv/archive/test-plan-v0.2.md`.

**v0.4.5 (2026-09-25)** applies the Spec Lead rulings to §5.8 (requirements header read at v0.4.14). The section changed materially:
- GPU rule (D-035): CPU-inert GPU keywords no longer make a case OUT. The three `*_hypre` inputs without TUNNEL_PRECONDITIONER and the 9 CPU-runnable `GPU_Tests/` inputs are IN; `dancing_eddies_ulmat_hypre` is UNCLEAR (Q4).
- New counts: 705 IN / 175 DEFERRED / 26 OUT / 35 UNCLEAR. The `GPU_Tests/` part goes beyond D-035's four inputs, so the Spec Lead should acknowledge it.
- FR-040 R3 and the FR-034 leak-area summation are recorded as open points with the Chief Architect and Pressure Solver Lead.
- VARIABLE_THICKNESS, CYLINDRICAL and TUNNEL_PRECONDITIONER stay UNCLEAR pending Q4.
- New non-GEOM particle anchors proposed: FR-052 `VV/cascadempi_obst` (A-47), plus FR-050/051 replacements for `cascadempi` (`scope_alignment.md` §8).
- v0.4.4 is archived at `docs/vv/archive/test-plan_v0.4.4.md`.

**v0.4.4 (2026-09-25)** adds the FR-006 scope suite. FR-006 ("Verification-suite coverage of in-scope features") entered in requirements v0.4.10 under owner decision D-033 and V&V action A-41; `requirements.md` header read at v0.4.11. Changes:
- New §5.8 "Scope suite (requirement FR-006 from spec v0.4.10)". Its case list is `docs/vv/scope_case_list.csv` class IN. The section also sets pass criteria, how UNCLEAR cases get resolved, and that DEFERRED cases (GEOM/CC_IBM/HT3D) stay FDS-only baselines.
- New files: `docs/vv/scope_case_list.csv`, `docs/vv/scope_alignment.md`, `vv-runs/tools/scope_filter.py`.
- No other section changes. §6 tiers and `case_inventory.md/.csv` are unchanged; the status changes for inventory cases are listed in `scope_alignment.md` §6 for the next inventory refresh.
- v0.4.3 is archived at `docs/vv/archive/test-plan_v0.4.3_pre-scope.md`.

**v0.4.3 (2026-09-25)** aligns the plan with the requirements v0.4.3 items (the file header now reads v0.4.4, last read; v0.4.4 changes only the kernel-language items K2/IR-007/NFR-045, and none of them changes a test here). Changes:
- FR-016(a) covers the first ghost layer only (§5.2).
- FR-005 (ii): exact summation starts at the per-box sum. New FR-005 (v) test for setup areas and volumes (§5.3).
- FR-010 negative tests: ratio rejection, blocking-factor report, MLMG coarsening warning (§5.4).
- D-022 frozen-input kernel comparison excludes `STORE_SPECIES_FLUX` and CC_IBM cases until A-34 closes (§5.5).
- NFR-011 uses the A-35 derived race_test inputs (§5.6).
- S1/NFR-030 timing rules (§5.7).
- A-35 and the ratio-inventory cross-check (§6, §10).

v0.3 is archived at `docs/vv/archive/test-plan-v0.3.md`.

| Item | Value |
|---|---|
| Refactor base / acceptance baseline source | **FireX `36975d765f`** (`36975d765fcead401e14b094a04f910ac42eab8a`), local branch `AMReX`, this repository (read-only). Recorded with `git -C (repo root) rev-parse --short HEAD`. `git describe` = `FDS-6.11.1-1244-g36975d765f`. Working tree clean. |
| FireX reference builds | Build requests went to the **GNU and Intel Build Chiefs on 2026-09-25, GNU first** (spec §3.3). **GNU: delivered** (Release + Debug; §3.5). **Intel: not yet built**; blocked on the oneAPI reinstall (A-18, R-32). |
| Existing binaries | GNU/OpenMPI and Intel/Intel MPI builds of master **`ce1f659`** under `(local FDS master checkout)/Build/…`. **Smoke tests only** (toolchain, MPI, comparison tooling). GNU starts from a clean shell and passes the version print; Intel cannot start because oneAPI is missing (§3.1). Never used as references. |
| Requirements read | **Now: `docs/requirements.md` v0.4.4** (header; changelog in `docs/README.md`, last read). This plan implements the v0.4.3 items: D-028 (exact setup areas/volumes, FR-005 (v)), D-030 (legacy non-{2,4} inputs stay FDS-only), A-33 closed, A-34 narrowed, A-35, FR-010 blocking-factor and MLMG-coarsening checks, FR-016(a) first ghost layer, NFR-011 A-35 inputs, NFR-030 timing rules (1)-(6). requirements.md was **not** edited. *Previous read (v0.3):* `docs/requirements.md` **v0.3.1** (Spec & Program Lead, 2026-09-25; file still being edited, last read). v0.3.1 adds the bit-parity scope (D-022: T0/T1 against baseline only for explicit kernels on frozen input vs single-mesh FDS; whole runs T2), refluxing and tolerance-only multi-level comparisons (D-023), FR-005 (decomposition/rank/thread independence), A-24/A-25 and R-36. At the time of writing, `docs/README.md` (decision/action log, changelog) and `risks.md` still show v0.3 and do not yet define D-022, D-023, A-24, A-25 or R-36. The requirements text is authoritative. requirements.md was **not** edited. |
| Companion files | `environment.md`, `case_inventory.md/.csv`, `verification_case_survey.csv`: still v0.2 content (see §11); tools in `(local V&V run directory)/tools/`; A-19 inputs in `(local V&V run directory)/inputs/A-19/`; A-35 race_test remesh in `(local V&V run directory)/inputs/A-35/` (README); FM_Burner area check `docs/vv/fm_burner_area_check.md`. The Legacy Mapper inventory is at `docs/inventory/` (`mesh_ratio_cases*.csv`, `global_reductions.csv`, `upstream_issue_candidates.md`) |

Conventions: **T0/T1/T2/T3** are the tolerance classes of requirements.md §2.2. **eps_H** is the single-solve pressure tolerance of requirements.md §2.1a. **Tier 1/2/3** are case-set sizes (run frequency) and are unrelated to T-classes. Runtimes are model estimates (±3×) until measured. "Derived" inputs (`VV/…`) are V&V copies of committed inputs with a documented change. They are created under `vv-runs/`, never in the worktree. Item IDs "R-n"/"D-n" without a leading zero are this plan's §9 items. "D-0nn" are decision-log IDs in `docs/README.md`.

---

## 1. Purpose and scope

The refactor moves FDS `MESH` data onto AMReX grids, first in *uniform mode* (zero refinement levels) and then with static and dynamic refinement. V&V has to show, phase by phase (roadmap.md), that:
1. in uniform mode, results match baseline FDS to the class each requirement names (FR-001/002/003/037/038/061, IR-001/006), results do not depend on box split, rank or thread count beyond FR-005's limits, and the per-step projection-solver selection behaves as specified (FR-039);
2. with refinement, the code conserves (FR-020..024), reproduces FDS's existing coarse-to-fine coupling where that is the stated reference (FR-016), and is more accurate than the coarse grid and approaches the uniform-fine grid (FR-014, NFR-032);
3. restart, determinism, parallelism, output and performance requirements hold (FR-015/070-074/080, NFR-010..034).

This plan covers verification (code against analytic solutions, conservation, and baseline FDS). Validation against experiments is out of scope until compute beyond the 8-core development machine is available (charter Q6). The one validation input used, Heskestad `Qs=1_RI=10` for NFR-032, serves as a code-to-code AMR-vs-uniform-fine comparison, not as validation.

## 2. Tolerance classes as V&V will apply them

| Class | requirements.md v0.3 definition (short) | How V&V evaluates it | Tool |
|---|---|---|---|
| **T0 bitwise** | Byte-identical output, used only for (1) **kernel-level parity with single-mesh FDS**: output arrays of an explicit kernel (mass, velocity predictor, divergence) on frozen input (D-022), and (2) **self-comparisons** (restart vs continuous, IR-006, run-to-run reproducibility, FR-005, A-09b): byte-identical `CHID_devc/hrr/mass.csv`, timing/version columns excluded, same ranks, `OMP_NUM_THREADS=1`, on derived `&DUMP SIG_FIGS=17` copies (D-014) | (1) array-by-array comparison of kernel dumps; (2) CSV comparison, timing columns removed by name pattern plus a per-case list (e.g. DEVC `cpu` in `dancing_eddies_*`). **Never** applied to whole runs of the AMR code against FDS, or to multi-mesh FDS as a kernel reference. | kernel-dump comparator (planned); `compare_csv.py --class T0` |
| **T1 tight** | Same uses as T0: every kernel output field, or every CSV column at every output time: \|x−x_base\| ≤ 1e-10·max(‖x_base‖∞,1e-30) (proposed), on the same SIG_FIGS=17 copies | As defined; per column/field, no interpolation (output times must match) | `compare_csv.py --class T1` (`--t1-factor`) |
| **T2 verification** | FDS criterion e ≤ Tol_FDS always holds, and e ≤ (1+m)·e_base + 0.1·Tol_FDS with m = 5%, or max(5%, 3σ) for chaotic cases; plot-only/convergent series: observed order ≥ min(p_base, nominal) − 0.1 plus finest-grid L2 rule (provisional, D-013) | As defined. σ comes from the calibration ensemble (§3.4 step 5). | dataplot metrics re-implemented in numpy (no pandas/matplotlib on the development machine), `compare_csv.py` for series |
| **T3 physical** | Time-averaged plume ΔT and w 10%, HRR/MLR 2% (1% prescribed), layer height 5% or one fine cell, layer T 5%, probe mean 10%/RMS 25%; discrimination rule \|AMR − fine\| ≤ 0.5·\|coarse − fine\| always; calibration max(value, 2× GNU-vs-Intel spread) (provisional, D-013) | Window ≥ last 50% and ≥ 10 puffing periods; per-case scripts | T3 statistics script (planned) |
| **eps_H** (not a T-class) | `H` agreement for **same-discretisation, single solve on frozen input**: eps_H = max(1e-8, 2.4e-12·N²), N = largest cell count per direction, mean removed, relative L2 (D-012). Not valid for multi-step outputs (use T2/T3) or for coarse-fine comparisons against UGLMAT (FR-016(c)). | Used for FR-002 (composite single level vs baseline GLMAT), FR-037 (`amrex::FFT::Poisson` vs baseline FFT) and FR-039 (FFT vs MLMG). Needs a frozen-input solve harness in the code under test (§5.1). | frozen-input H comparator (planned) |

**Bit-parity scope (requirements §2.2, D-022/D-023), as V&V applies it:** whole runs of the AMReX code, single-level included, are judged at **T2** against FDS after the first pressure solve (Crayfishpak is replaced). Kernel parity is against **single-mesh** FDS only (committed, or a derived single-mesh copy per A-24). Multi-level runs are never compared bitwise with multi-mesh FDS. eps_H still applies to single solves on frozen input. Kernel-level parity needs frozen-input kernel dumps from both codes: from baseline FDS this is an instrumented build like A-09b (§3.6). **Until A-34 closes, the frozen-input kernel comparison excludes every case that sets `STORE_SPECIES_FLUX` or runs CC_IBM** (requirements §2.2 note, IR-007). The cases this removes are listed in §5.5. Whole-run T2 comparisons of those cases are unaffected.

Which requirement uses which class, and which class each case is judged at: §7 and the `tclass` column of `case_inventory.csv` (the inventory's `tclass` still reflects v0.2 and must be re-derived for D-022).

## 3. Baselines, reference results and the FireX build request

### 3.1 What exists today
- **Smoke-test binaries (ce1f659).** `fds_ompi_gnu_linux` (GNU Fortran 14.2.0, OpenMPI soname .40, OpenMP, HYPRE 3.0.0, SUNDIALS 7.5.0, no MKL) and `fds_impi_intel_linux` (ifx 2026.1.1, Intel MPI 2021.18, MKL static, no OpenMP). Details are in environment.md §2. They are used **only** to prove the toolchain, MPI launch and comparison tooling (G0 below).
- **Smoke status (corrected 2026-09-25):**
  - GNU: **PASS** (version print). From a clean shell, `ldd` resolves `libmpi_usempif08.so.40` to `/lib/x86_64-linux-gnu` (Debian `libopenmpi40`/`openmpi-bin` 5.0.7), and `/usr/bin/mpirun -np 1 fds_ompi_gnu_linux` prints the version banner (GCC 14.2.0, Open MPI 5.0.7, HYPRE 3.0.0, SUNDIALS 7.5.0). The earlier FAIL came from a non-clean shell environment and is withdrawn.
  - Intel: FAIL. `libimf.so` exists nowhere on disk and there is no `/opt/intel/oneapi/setvars.sh` or `mpiexec`. The oneAPI install is gone (the development machine rebooted Sep 23). Reinstalling oneAPI is a prerequisite for the Intel FireX build.

### 3.2 Where reference results must come from
**Reference (baseline) results for every acceptance test must come from FDS built from FireX `36975d765f`** (requirements.md §2.1: CMake, out of tree, HYPRE GPU options OFF, same compiler, MPI, flags, ranks and threads as the build under test; `ce1f659` results are not acceptance references). The `ce1f659` binaries can't serve even as a proxy: FireX changes radiation, init, pressure (GPU HYPRE paths, resource sets) and output code (`case_inventory.md` §1a), and it links a different HYPRE (§3.5). V&V does not build: A-08 records the build, it does not perform it.

**Build requests** went to the GNU Build Chief and the Intel Build Chief on **2026-09-25**, GNU first, with the spec in §3.3. Status:
- **GNU: delivered 2026-09-25** (§3.5). G0 and the Tier 1 capture are being run by another V&V worker, which records the results in `baseline_status.md`. This plan does not duplicate that record.
- **Intel: pending.** It needs oneAPI reinstalled first (A-18). Until then, every acceptance item that names both toolchains (NFR-020, the T2 σ ensemble, T3 calibration, A-09 on both toolchains) is GNU-only and recorded as incomplete.

### 3.3 Build request spec (to the GNU and Intel Build Chiefs; sent)

| # | Item | GNU build | Intel build |
|---|---|---|---|
| 1 | Source | this repository at `36975d765fcead401e14b094a04f910ac42eab8a` (branch `AMReX`). Check `git rev-parse HEAD` and an empty `git status --short` before and after. **The worktree is read-only**: nothing may be written into it. | same |
| 2 | Build system (primary, = baseline per requirements §2.1) | CMake ≥ 3.24, **out of tree**: `cmake -S (repo root) -B (local build directory)/firex-36975d7/ompi_gnu_rel -DCMAKE_BUILD_TYPE=Release`. Take the compiler wrappers from preset `ompi_gnu` (`CC=mpicc CXX=mpic++ FC=mpifort`). **Don't use the preset `binaryDir`** (`./Build/cmakeb/<preset>`): it is inside the worktree. | same, preset `impi_intel` wrappers (`mpiicx`/`mpiifx`), build dir `…/impi_intel_rel` |
| 3 | Compiler / MPI | Same as the existing binary: GNU Fortran 14.2.0 (Debian 14.2.0-19) + OpenMPI 5.0.x (soname .40), as NFR-020 requires | Same as the existing binary: ifx 2026.1.1 + Intel MPI 2021.18, with MKL visible to CMake (`find_package(MKL CONFIG)` → `WITH_MKL`) so the `*_pardiso` cases work as they do in the existing Intel binary |
| 4 | Options | `USE_HYPRE=ON`, `USE_HYPRE_NVIDIA/AMDGPU/INTELGPU=OFF`, `USE_SUNDIALS=ON`, `USE_OPENMP=ON` (CMake default). Record any deviation. | same. Note: the existing Intel makefile binary had **no** OpenMP; CMake enables it by default. Keep the CMake default and record it. |
| 5 | HYPRE | FireX CMake fetches hypre commit `63331f19` (requires ≥ 2.32.0). **Outcome (GNU):** fetched pin used (`USE_SYSTEM_HYPRE=OFF`). `63331f19c7` = `v2.32.0-24`, an unreleased commit between 2.32 and 2.33, **not 3.0.0**. The system v3.0.0 libs were rejected by `find_package(HYPRE 2.32.0)` (SameMajorVersion). See §3.5. | **Must use the same fetched commit `63331f19c7`**, so both toolchains share one HYPRE for FR-016/FR-038 references. Record its `git describe`. |
| 6 | SUNDIALS | v7.5.0 (fetched, or the existing libs); record. GNU used the system 7.5.0 libs. | same |
| 7 | HDF5 / VTK | Not available in CMake (no option). Leave out of the primary build; record "HDF5: no". **No HDF5 build exists yet**; it waits on D-019 (A-14, Q5). | same |
| 8 | Provenance (A-08) | `fds -V` must print a non-empty revision (CMake takes `git describe`; if building from a copy, pass `-DGIT_HASH=…`). Deliver: `CMakeCache.txt`, verbose compile lines (actual flags), `mpifort --version` / `mpifort -showme`, HYPRE/SUNDIALS versions, `ldd` output, sha256 of the binary, and the environment script that makes it run **on the development machine** | same, plus `mpiifx -show`, MKL version |
| 9 | Runnable on the development machine, reboot-proof | OpenMPI 5.0.7 runtime is present (Debian packages). Write down the environment setup (a script plus notes) so it survives the next reboot. | **Prerequisite:** reinstall oneAPI (ifx 2026.1.1, Intel MPI 2021.18, MKL) under `/opt/intel/oneapi` with `setvars.sh`; it is missing since the Sep 23 reboot. Write down the setup so it survives the next reboot. |
| 10 | Debug variant (requested, lower priority) | Same configuration, `CMAKE_BUILD_TYPE=Debug` (preset `ompi_gnu_db`), for A-09 NaN/bounds checks on the three AMR cases. **Delivered** (§3.5). | `impi_intel_db` |
| 11 | Optional: HDF5 variant for FR-074 / A-14 | Makefile target **`ompi_gnu_linux`** (`Build/ompi_gnu_linux/make_fds.sh`) with HDF5 1.14.5 (`build_thirdparty_libs.sh` builds it only if `$FIREMODELS/hdf5` is a clone at tag `hdf5_1.14.5`; otherwise it silently omits HDF5). Makefile builds run **inside** `Build/<target>`, so they must use an exported copy of the commit (e.g. `git archive 36975d765f`) outside the worktree. On hold pending D-019. | target **`impi_intel_linux`** (no OpenMP; `-O2 -ipo`) |
| 12 | Later (Phase 3): instrumented baseline for FR-016(a) | **A-09b**, placeholder spec in §3.6. | same |

Build names for reference: makefile targets `ompi_gnu_linux`, `impi_intel_linux` (+ `_db`, `_dv`, `impi_intel_linux_openmp`); CMake presets `ompi_gnu_rel`, `impi_intel_rel` (+ `_db`, `_dv`).

**V&V acceptance of a delivered build (G0):** `fds -V` shows `36975d765f`. `ns2d_16` runs on 1 rank and `obst_activation_default` on 4 ranks, each finishing normally in < 2 min. A setup-only sweep (`T_END=0` copies) over all 941 inputs records which ones stop at setup; that list becomes the IR-001 reference.

### 3.4 Reference capture procedure (after a build passes G0)
1. Order (requirements A-09, R-34): run `ns2d_16_int_1to2_refinement` and its derived `SOLVER='UGLMAT HYPRE'` copy on the FireX build **first**, on the Release and the Debug build of each toolchain, before FR-016 is treated as a gate. Then the two `emb` cases (A-09), Tier 1 and Tier 2. The int_1to2 pair is in Tier 1, so its GNU Release result comes with the GNU Tier 1 capture. The GNU Debug runs follow, and Intel follows A-18. **FR-016 does not gate anything until A-09 shows that baseline runs the case** (pass = runs to T_END with no NaN; analytic-error values recorded as baseline values). If the UGLMAT copy fails, FR-016 goes back to the Spec Lead and Pressure Lead for a replacement case (R-34).
2. Every run: `OMP_NUM_THREADS=1` unless the case is a threading case. Load < 2 and free RAM ≥ ranks×0.5 GB + 1 GB checked before start. `mpirun --bind-to none` (shared machine). Environment from the Build Chief's script (GNU: `(local build directory)/env/gnu_ompi_env.sh`). Watchdog per NFR-011 (R-10).
3. Every Tier 1 case runs twice: once as committed, and once as a SIG_FIGS=17 copy (D-014) for the T0/T1 self-comparisons that remain under D-022 (restart vs continuous, IR-006, reproducibility, A-09b acceptance). Chaotic cases (fire/LES) and the FR-002 cases also get the extra variants their class needs (GLMAT copy, alternate rank count).
4. Stored under `(local V&V run directory)/baseline/<toolchain>/<case>/` with a manifest: input sha256, binary sha256, `fds -V`, env, ranks/threads, launcher line, wall time, exit status, load/free-RAM at start. **The HYPRE version in the manifest comes from the build's `PROVENANCE.md`, never from the FDS run header** (the header prints a hard-coded "3.0.0"; §3.5). Outputs are read-only after capture.
5. T2 margin calibration ensemble (for chaotic cases): GNU vs Intel and firebot rank count vs one alternate rank count, giving σ for T2. Incomplete until the Intel build exists.
6. Measured wall times replace the model estimates, and cases are re-tiered.

### 3.5 FireX GNU build as delivered (GNU Build Chief, 2026-09-25)

| Item | Value |
|---|---|
| Release binary | `(local build directory)/firex-36975d7/ompi_gnu_rel/fds`, sha256 `9d3b5991447a70b8cd1d067e7db6fec8506479b32fd62961339f8d387e8050ec` (checked by V&V with `sha256sum`) |
| Revision | `FDS-6.11.1-1244-g36975d765f-AMReX` |
| Toolchain | GCC 14.2.0, Open MPI 5.0.7, OpenMP on, no MKL |
| Libraries | SUNDIALS 7.5.0 (system libs, `(local GNU third-party library tree)/libs/sundials/v7.5.0`, shared); HYPRE fetched commit `63331f19c7` (= `v2.32.0-24`, static); HDF5 not used |
| Fortran flags | `-O3 -cpp -std=f2018 -frecursive -ffpe-summary=none -fall-intrinsics -fopenmp` |
| Environment script | `(local build directory)/env/gnu_ompi_env.sh` (sets `HWLOC_LIBXML=0`) |
| Provenance | `(local build directory)/firex-36975d7/ompi_gnu_rel/provenance/PROVENANCE.md` (CMake cache, configure/verbose build logs, flags, ldd, third-party hashes) |
| Debug binary | `(local build directory)/firex-36975d7/ompi_gnu_db/fds`, built with `-O0 -ggdb -finit-real=snan -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace` (for A-09 NaN/bounds checks) |
| HDF5 / VTK variant | none yet; waits on D-019 (A-14, Q5) |

**HYPRE note (risk to FR-016/FR-038 references; also NFR-021, R-30):**
- FireX's CMake pin `63331f19c7` is `v2.32.0-24`, an unreleased commit between 2.32 and 2.33. FireX requires HYPRE 2.x ≥ 2.32, so the development machine's system HYPRE 3.0.0 was rejected.
- FDS's run header still prints "Hypre library version: 3.0.0", because CMakeLists hard-codes that string. **Manifests take the HYPRE version from provenance, never from the header.**
- Consequences:
  - The FireX reference uses a different HYPRE from the old `ce1f659` smoke binary (a real 3.0.0). Nothing HYPRE-related may be compared across the two.
  - The **Intel build must use the same fetched commit**, or the GNU and Intel HYPRE references differ.
  - requirements.md §2.1/NFR-021 and risks R-30 describe `63331f19` as "labelled 3.0.0" and name the development machine's 3.0.0 as the HYPRE. The linked code is 2.32.0-24. Reported to the Spec Lead; not edited here.
  - If the AMR build later links a different HYPRE (e.g. one HYPRE for FDS and AMReX MLMG, R-30), FR-016 (UGLMAT-HYPRE reference) and FR-038 (T2 vs baseline HYPRE runs) are no longer like-for-like. Either keep `63331f19c7` in the build under test, or rebuild and re-capture every HYPRE baseline on the new version and record the difference (R-09). Trigger: any HYPRE change, or baseline `H` shifting beyond eps_H on a single frozen solve.

### 3.6 A-09b: instrumented baseline build (placeholder spec)
- **Purpose:** FR-016(a) needs fine-side ghost values of ρ, ρY_α, temperature and the ghost-averaged properties (`MU`, `KRES`, `D`/`DS`) on the first step of `ns2d_16_int_1to2_refinement`, on both sides of the coarse-fine interface. Baseline FDS doesn't write them. Coarse-fine (ratio-2) interface data also serves FR-016(b)/(c) diagnostics.
- **Build:** FireX `36975d765f` plus a V&V ghost-cell / coarse-fine dump patch that doesn't change arithmetic (writes only, after the ghost fill; no reordering of operations). Applied to a **scratch copy** exported outside the read-only worktree (e.g. `git archive 36975d765f` into `(local build directory)/scratch/…`), built out of tree with the same CMake configuration, compiler, flags and HYPRE commit as §3.5. **Never committed** to branch `AMReX` or anywhere else (D-005). The patch file and its sha256 are stored with the build provenance.
- **Owner:** AMR V&V Lead (with the GNU Build Chief), per requirements A-09b; due before the Phase 3 unit checks.
- **Acceptance:** the instrumented build is accepted only if it is **T0 against the uninstrumented build** of the same toolchain (§3.5) on **all Tier 1 CSVs** (SIG_FIGS=17 copies, same ranks, `OMP_NUM_THREADS=1`). Any T0 failure rejects the patch.
- **Open (detailed spec needed):**
  - which arrays and which code points (after which `MESH_EXCHANGE` code / ghost-averaging site in `docs/inventory/interface_averaging_sites.csv`);
  - the dump format and the matching dump in the AMR code;
  - whether H and the post-projection face velocities are dumped for FR-016(c)/FR-032 diagnostics.
- Possible extension (v0.3.1, D-022): the same scratch-copy mechanism could also dump the explicit-kernel outputs (mass, velocity predictor, divergence) on frozen input that the FR-001/FR-002 kernel-parity checks need from single-mesh FDS. To be decided with the spec.
- These need input from the **AMR Pressure Solver Lead** and the **AMReX Integration Lead**. Not yet requested.

## 4. Environment and execution rules
- Development machine: 8 cores (1 socket, AVX-512), 16 GB RAM with **no swap**, shared with other workloads (2.4–6.4 GB free observed). At most 8 ranks per case and at most 8 ranks in total across concurrent cases. No case above 12 GB (NFR-031).
- Python 3.13 + numpy only (no pandas, scipy or matplotlib). FDS's `Utilities/Python` plotting scripts can't run here as-is, so V&V re-implements the needed metrics in numpy (§8).
- Timing tests (NFR-030/033/034, and the NFR-032 wall-clock ratio) are valid only on an otherwise idle machine: load < 1, ≥ 8 GB free, no other V&V jobs. Otherwise the run is flagged invalid rather than failed (D-013).
- Smokeview: not installed. Reference version SMV 6.11.2 (FR-072).
- FireX runs use the Build Chief's environment script (GNU: `(local build directory)/env/gnu_ompi_env.sh`).

## 5. Test gates

| Gate | When | Content | Class / criterion |
|---|---|---|---|
| **G0** Toolchain smoke | Each new binary (incl. ce1f659 once runtimes exist) | `fds -V`; `ns2d_16` (1 rank); `obst_activation_default` (4 ranks); `compare_csv.py` self-test (identical, perturbed and re-timed synthetic CSVs) | Runs / exit 0. No T-class (ce1f659 output is never a reference). GNU FireX build: in progress by the capture worker (`baseline_status.md`). |
| **G1** Input compatibility (IR-001) | Phase 2, then each milestone | Setup-only (`T_END=0`) copies of all 941 inputs | Same set of setup stops and messages as baseline (message diff) |
| **G2** Tier 1 regression | Every merge to `AMReX` (≈ 8 min model, budget ≤ 30 min) | 40 runs (§6) | Whole runs vs FDS: **T2** (single-mesh and same-resolution multi-mesh; FR-002 also T2 vs the single-mesh equivalent and vs multi-mesh FFT and GLMAT baselines) (D-022). Kernel parity (FR-001/FR-002): explicit kernels on frozen input at **T1, T0 per Q1**, vs single-mesh FDS (A-24 copies for multi-mesh anchors). Restart per FR-080 (G10). Determinism/reproducibility T0 (FR-005(iv)). |
| **G3** Tier 2 | Nightly while code changes; every phase exit | 123 runs (≈ 79 min sequential, ≈ 30–40 min packed on 8 cores) + optional 8 | T2 |
| **G4** Suite non-regression (FR-003) | Phase 2 exit, Phase 10 | All firebot-listed cases with ≤ 8 ranks (870 of 884) | T2. The 14 cases needing 9–64 ranks are exceptions recorded and approved by the project owner |
| **G5** Build-option equivalence (IR-006) | Phase 2 onward | FR-001 anchors with `USE_AMREX=OFF` | **T0** |
| **G6** FireX HYPRE options (FR-038) | Phase 2 | The 4 `Pressure_Solver/*_hypre.fds` at 1 and 4 ranks, `FDS_RANKS_PER_GPU` unset and =2 (the variable changes the matrix partition on CPU too: RS masters gather and solve, pres.f90:3411-3440) | **T2** vs baseline at the same setting (v0.3.1: whole runs, multi-mesh inputs, D-022); same HYPRE commit as baseline (§3.5) |
| **G7** Static two-level (FR-016) | Phase 3 (a–c), Phase 4 (d); **gating only after A-09 shows baseline runs the case (R-34)** | `ns2d_16_int_1to2_refinement` + `VV/…_uglmat` (`SOLVER='UGLMAT HYPRE'`) only (D-015). `ns2d_16_emb_1to2_refinement` is informational and non-gating ((a)-type comparison and FR-014 only). | (a) T0 on first-step fine-side ghosts, **first ghost layer only** (§5.2), vs the A-09b instrumented baseline (component-level only); multi-step and end-state comparisons vs `int_1to2` FDS are T2/T3, never bitwise (D-023); (b) coarse side not matched, FR-020/021 instead; (c) FR-032 normal-velocity mismatch at machine zero (target), `H` error against the exact `ns2d` solution no worse than UGLMAT-HYPRE's, observed order ≥ 1.8. **eps_H does not apply here** (not same-discretisation; replaces the v0.2 "10× solver tol" rule; A-16 closed); (d) T2 vs UGLMAT-HYPRE baseline |
| **G8** Conservation (FR-020/021/022/024) | Phase 3/4 | `species_conservation_1..4`, `energy_budget_*`, `simple_duct`, mass-balance cases with a refined patch | Requirement numbers (1e-12/step, 1e-10 cumulative; FR-022 with absolute floor) |
| **G9** Accuracy (FR-014, NFR-032) | Phase 4, 10 | `ns2d_{8..64}` + static patch, `ns2d_16_*_refinement`; NFR-032 plume: Heskestad `Qs=1_RI=10` (D-008) on the trimmed 32×32×80 level 0, with the A-19 64×64×160 uniform-fine reference (§6) | FR-014: refined ≤ coarse, refined-region L2 ≤ 2.0× uniform-fine, order ≥ 1.8. NFR-032: T3 (Lf within one fine cell 0.057 m, HRR 1%, centreline ΔT and w 10%, discrimination rule) at ≤ 50% of uniform-fine wall time |
| **G10** Determinism / restart (FR-015, FR-080/081) | Phase 3, 9 | repeat runs, hierarchy dumps; `restart_test1a/b` + continuous copy, `device_restart_*` (+ `device_restart_base_case`), `restart_ulmat_*` | Hierarchy identical. FR-080 uniform mode: T0 proposed, **T1 minimum**. Only if baseline's own restart-vs-continuous difference misses T1 (measured at capture, A-08) does the criterion become "no worse than baseline's own restart-vs-continuous difference". AMR mode: T1 with an identical hierarchy after restart. FR-081: T2 |
| **G11** Parallel / threads / decomposition (NFR-010/011/012, FR-005) | Phase 2+ | Baseline FDS: {1,2,4,8} ranks ≤ mesh count (D-016), MULT-split copies for wider sweeps (T2 vs original). AMReX code: every anchor at 1, 2, 4 and 8 ranks, two `max_grid_size` values, two thread counts; a pressure-zone case (`zone_break_fast`) for FR-005(ii); FR-005(v) setup-area/volume cases (§5.3); 3 repeats per configuration; `race_test_1/4` (AMR mode: A-35 `_r4` copies; FDS baselines: committed originals; §5.6) | FR-005: reduction-free explicit stages byte-identical across box split/ranks/threads (i); zone integrals and gauge exact, with exact accumulation from the per-cell/per-box sum up, not only at the Allreduce (ii); `H` within eps_H across box split/ranks (iii); bitwise run-to-run at fixed ranks/threads/layout, `OMP_DYNAMIC=false` (iv); setup areas/volumes byte-identical across box split/ranks/threads (v). Whole runs T2 across ranks (uniform), T3 (AMR); NFR-011 watchdog |
| **G12** Performance / memory (NFR-030/031/033/034) | Phase 2 (measure), 10 (meet) | `openmp_test64a`; one 8-rank multi-mesh anchor; NFR-033 `GPU_Tests/HYPRE_GPU_SCALING/test_8mesh_NOFRPG` 20-step copy | Timing ratios, median of 3, idle machine only; timing-run rules (1)-(6) of NFR-030 (§5.7) |
| **G13** Output (FR-071/072/074) | Phase 9 | header diff; Smokeview load; VTK (needs HDF5 build, D-019) | Exact header match (I); D for Smokeview/ParaView |
| **G14** Pressure single-solve checks (FR-002 `H`, FR-037, FR-039, FR-005(iii)) | Phase 2 (prototype P2 for FR-039), then Phase 4 | Frozen-input solves (§5.1) | **eps_H** (§2) on a single solve after gauge fixing; switch timing exact. Multi-step outputs of the same runs are judged T2, never eps_H |

Negative tests (FR-004, FR-044, IR-002): one input per excluded feature, expecting `SETUP_STOP` and the named `ERROR(nnn)`. Written in Phase 3 as the supported feature set is defined. FR-010 refinement-ratio and blocking-factor negative tests: §5.4.

### 5.1 FR-039 (and FR-037/FR-002 `H`) proposed test

FR-039: the driver uses `amrex::FFT::Poisson` while the hierarchy is single-level and composite MLMG once a refined level exists. Both paths use the same pressure gauge, and MLMG uses `setMaxOrder(2)`. Owner: AMR Chief Architect / AMReX Integration Lead. Acceptance is part of the prototype P2 check (roadmap), then Phase 4. V&V supplies the cases, the frozen-input comparator and the switch check.

| Sub-test | Method | Cases (single mesh, one BC type per face, so FFT-able) | Criterion |
|---|---|---|---|
| FR-039(a) FFT vs MLMG, same input | Dump the Poisson RHS and BCs at steps 1 and ~N/2 of a uniform run. Solve once with `FFT::Poisson` and once with MLMG (maxorder 2, same gauge applied after the solve), then compare H (mean removed, relative L2) | `ns2d_16` (N=16), `ns2d_32` (32), `shunn3_32` (32; periodic), `csmag_32` (32; 3-D periodic), `dancing_eddies_1mesh` (N=300; inflow/open/wall faces + OBSTs, same all-cell operator in both). Optional: `shunn3_128` (128) | eps_H: 1e-8 for N ≤ 64; 3.9e-8 (N=128); 2.2e-7 (N=300) |
| FR-037 FFT fast path vs baseline | Same frozen input, `FFT::Poisson` vs baseline `pois.f90` FFT (FR-001 cases with one BC type per face) | as (a) | eps_H; plus D (FR-001 cases run with the FFT path selected) |
| FR-002 `H` vs baseline GLMAT | Same frozen input, composite single-level solve vs baseline GLMAT | `shunn3_4mesh_32` (M2a), then the other FR-002 anchors | eps_H |
| FR-039(b) switch | A derived `ns2d_32` copy with a static refined patch added at step n1 and removed at step n2 (needs the `&AMR` namelist, IR-003). Check the per-step solver log. | derived `VV/ns2d_32_switch` | Solver = FFT exactly on single-level steps and MLMG exactly on multi-level steps. Full-run outputs T2 vs baseline |
| FR-039(c) order at maxorder 2 (P2 open item) | Refinement study with maxorder 2 and 3 on the `ns2d` series + static patch | `ns2d_{16,32,64}` + patch | Observed order ≥ 1.8 at the coarse-fine interface (FR-014). If order 2 fails: order 3 on multi-level hierarchies, and the eps_H switch check then applies only to the order-2 single-level comparison |

Needs from the code under test: a frozen-input solve hook (RHS/BC dump + replay) and a per-step solver-selection log line. These go into the P2 interface. Not yet agreed with the Integration Lead.

### 5.2 FR-016(a) ghost-layer scope (v0.4.3)

The bitwise (T0) ghost-value check at a level jump covers **only the first fine-side ghost layer** of the transported and diffused fields (ρ, ρY_α, T, and the ghost-averaged `MU`, `KRES`, `D`/`DS`) on step 1, compared with the A-09b instrumented dump of `int_1to2`. The second ghost layer is not compared, either bitwise or against a separate tolerance.

**Rationale:**
- At a level jump FDS fills the second layer as a zero-gradient copy of the first (wall.f90:351-386; first order). The P1 fill-patch (`FillPatchTwoLevels` + `PCInterp`) is not required to reproduce that, so a bitwise check there would test an FDS artefact rather than a requirement.
- The second layer enters only through the wider advection stencils next to the interface, i.e. through the interface fluxes. Those are covered by two checks: the conservation/refluxing checks (G8, FR-020/021/024: round-off with refluxing, D-023) and the whole-run T2 comparison of `int_1to2` against its UGLMAT-HYPRE baseline (G7 (d)).
- A second-layer error large enough to matter shows up there. One too small to show up there isn't worth a dedicated check.
- Confirmed by the V&V Lead and recorded in requirements FR-016(a).

### 5.3 FR-005 (ii) and (v): exact sums from the per-box level (v0.4.3, D-028)

- **(ii) amended:** the exact (order-independent, fixed-point) accumulation of `DSUM`/`PSUM`/`USUM` and the pressure gauge starts at the per-cell/per-box sum, not only at the `MPI_Allreduce`. 68 of the 91 order-dependent rows in `inventory/global_reductions.csv` are per-process loops before the MPI call, and regridding reorders cells within a box.
  - Test (unchanged case, stronger check): on `zone_break_fast` (and on `zone_shape` for multiple zones), `D_PBAR_DT` and the gauge constant, dumped per step, are **byte-identical** across 1/2/4/8 ranks, two `max_grid_size` values and 1/2/4 threads. In AMR mode they also stay identical before and after a regrid that only reorders boxes.
- **(v) new, setup-time areas and volumes:** the driver computes `FDS_AREA`/`AREA_ADJUST` and every setup-time area and volume sum in `global_reductions.csv` over the whole domain, with the exact sum.

| Sub-test | Input | Configurations | Criterion |
|---|---|---|---|
| (v-a) cylinder spanning boxes | Legacy Mapper MULT-voxel cylinder reproducer (R=H=0.1 m; variants a–e of upstream issue candidate #1: 1 mesh, 2 meshes/1 process, 2 meshes/2 processes, moved split, two cylinders sharing one MULT) as V&V copies | AMR code: level 0 split into 1, 2 and 4 boxes (two `max_grid_size`), 1/2/4/8 ranks, 1/2/4 threads | `AREA_ADJUST` per face class (top/side/bottom) and the resulting fuel MLR **byte-identical** across all configurations (T0), and equal to single-mesh FDS (see note) |
| (v-b) real case | `FM_15cm_Burner_CH4_2cm.fds` (MULT-voxel cylinder, top split over 4 meshes; §6 note and `fm_burner_area_check.md`) and its single-mesh copy `vv-runs/inputs/fm_burner_area/CH4_2cm_1mesh.fds` | as (v-a) | fuel MLR / (MASS_FLUX·πR²·tanh t) = 1.00000 in every configuration. FDS as committed gives 0.52083 on 1 process and 1.00000 with 1 mesh per process, so FDS multi-mesh is **not** a reference here. |
| (v-c) volumes | an anchor with pressure zones and OBSTs (`zone_shape`) plus `random_meshes` | as (v-a) | zone volumes and total gas volume byte-identical across configurations, equal to single-mesh FDS |

- How the value is observed: FDS reports no adjusted area. Setup-only FDS never computes it (T_END=0 stops before `ADJUST_OBST_SHAPE_AREA`, init.f90:917; main.f90:238-248). The AMR code needs a setup diagnostic that prints `FDS_AREA`/`AREA_ADJUST` per shaped OBST/MULT and the zone volumes, **also in setup-only mode**. Until it exists, the fuel MLR of a few time steps is the observable, as in the FM_Burner check.
- **Note / clarification requested:** "reproduces single-mesh FDS" can hold bitwise only if FDS's sequential floating-point sum over one mesh happens to equal the exact sum rounded once. For equal voxel faces that is usually, but not always, the case. V&V applies **T0 across AMR configurations** and **T1 (1e-10 relative) against single-mesh FDS** unless the Spec Lead decides otherwise.

### 5.4 FR-010 negative and input-check tests (v0.4.3)

AMR mode accepts only refinement ratios 2 and 4, the same in every direction (D-030). Expected behaviour: a setup stop at input time with an error naming the mesh pair and the ratio. The blocking-factor check reports alongside it, and the MLMG coarsening check is a WARNING only.

| Test | Input | Expected |
|---|---|---|
| N-1 3:1 | `VV/ratio3_int.fds`: copy of `ns2d_16_int_1to2_refinement` with the fine patch `IJK=24,1,24` (3:1 against the 4×4 coarse meshes) | Rejected: error names the fine mesh and each coarse neighbour and "ratio 3" |
| N-2 direction-dependent | `VV/ratio_dirdep_int.fds`: same copy with the fine patch `IJK=16,1,32` (2:1 in x, 4:1 in z; each ratio is in {2,4} but they differ) | Rejected: error names the pair and both ratios |
| N-3 5:1 legacy | committed `Thread_Check/race_test_1.fds` | Rejected: meshes 1–3, 2–3, 3–4, 3–5, 3–6 at ratio 5 |
| N-4 legacy 3:1 and direction-dependent (secondary) | one each from the D-030 list, setup-only: `impinging_jet_*` (3:1) and `BST_FRS_6` (8:4) | Rejected with the named pairs |
| P-1 2-D positive control | committed `ns2d_16_int_1to2_refinement` | **Accepted**. Its single y cell (IJK(2)=1 on every mesh) must not count as a direction ratio of 1 (the inventory ignores y in 2-D inputs; inventory README §4 (f)). Otherwise the FR-016 gating case would be rejected. |
| P-2 A-35 remesh | `vv-runs/inputs/A-35/race_test_1_r4.fds` | Accepted for ratios (all five interfaces 4:1). With the default `blocking_factor` 8 it stops with the named blocking-factor error (level-0 domain 34×18×32; 34 and 18 not divisible by 8). With `amr.blocking_factor = 2 8` it starts. |
| W-1 MLMG coarsening warning | as P-2 | WARNING printed, coarsest MLMG grid **17×9×16** reported; the run continues |

The derived N-1/N-2 copies are created under `vv-runs/inputs/FR-010/` in Phase 3, when the `&AMR` input exists (IR-003). They are Verification T (setup stop + message text).

### 5.5 D-022 frozen-input kernel comparison: cases excluded until A-34 closes

- **Rule:** exclude every case that sets `STORE_SPECIES_FLUX` or runs CC_IBM.
- **What triggers them:**
  - FDS sets `STORE_SPECIES_FLUX` internally when a DEVC/SLCF/BNDF asks for `ADVECTIVE|DIFFUSIVE|TOTAL MASS FLUX X|Y|Z` or `TOTAL MASS FLUX WALL` (read.f90:17011-17014).
  - Any `&GEOM` turns CC_IBM on for a time-stepping run (geom.f90:2698), as does `&MISC CC_IBM=T`.
- **Scan** (non-comment lines of the 180 committed inputs behind our 184 inventory rows; the 4 derived rows inherit from sources with none of these):

| Case | Tier / anchor | Trigger | Was it a kernel-comparison candidate? |
|---|---|---|---|
| `Sprinklers_and_Sprays/cascadempi.fds` | T2, anchor (FR-050/051/052) | `&GEOM` (lines 24-26) → CC_IBM | yes: FR-002 kernel vs a single-mesh/GLMAT copy (inventory `tclass`) |
| `Species/mass_flux_wall_yindex.fds` | T2 | `TOTAL MASS FLUX WALL` | yes (v0.2 `tclass` "T1 (T0 stretch)", single mesh) |
| `Species/mass_flux_wall_zindex.fds` | T2 | `TOTAL MASS FLUX WALL` | yes (as above) |
| `Species/mass_balance_gas_volume.fds` | T2 | `TOTAL MASS FLUX X/Y/Z` | yes (as above) |
| `Species/mass_balance_reac.fds` | T2 (8 meshes) | `TOTAL MASS FLUX WALL` | yes ("T1 first 100 steps"; would need an A-24 copy) |

- None of the FR-001 single-mesh anchors listed in §7 and none of the A-24 copies (`shunn3_4mesh_32`, `dancing_eddies_default`, `symmetry_test_mpi`) is affected. Their whole runs stay T2.
- Re-add these cases when A-34 closes: either `STORE_SPECIES_FLUX` output is shown not to read the J=0/K=0 rows, or those rows are restored.

### 5.6 NFR-011 race tests (A-35, D-030)

- **AMR mode** uses the A-35 derived inputs: `vv-runs/inputs/A-35/race_test_1_r4.fds` and `race_test_4_r4.fds`, with mesh 3 at `IJK=24,24,16` (dx 0.0125) so all five interfaces are 4:1.
  - Level 0: 34×18×32 at dx 0.05, `amr.blocking_factor = 2 8` (ratio-4 level 1) or `2`.
  - Pass: same metric as upstream (race_test_4 vs race_test_1 TMP/VEL at t=4.49–4.51 s, abs. error 0.01), no hang (watchdog).
  - The AMR reference values come from FDS runs of the `_r4` copies (1 and 4 threads), not from the originals.
- **FDS baselines** stay the committed `race_test_1.fds`/`race_test_4.fds`.
- The remesh raises the burner top and vent plane from 0.120 to 0.125 m, halves the burning-cell count (32 → 16; total fuel flow preserved by AREA_ADJUST) and moves the probe cell. So `_r4` results are never compared with the original baselines.
- Details: `vv-runs/inputs/A-35/README.md`.

### 5.7 S1 / NFR-030 timing-run rules (v0.4.3)

Every AMReX timing run (NFR-030, ADR-002 S1 / D-007) records and reports the following. V&V reads "acceptance requires all of" as: a run missing any item is invalid, not failed, like a run on a busy machine.
1. The pressure solver used (`amrex::FFT::Poisson` or MLMG), from the per-step solver log (G14).
2. Pressure-solve time separately from total time, from the TinyProfiler regions `FFT::Poisson::solve` and `MLMG::solve()`, inclusive and exclusive, per rank.
3. One level-0 (uniform) run forced onto MLMG, for a same-grid FFT-vs-MLMG comparison, so that AMR cost ratios don't mix the solver switch with the added cells.
4. Profiled builds link a separate AMReX install with `AMReX_TINY_PROFILE=ON`; the production-like install (`(local AMReX install)`) stays OFF. Headline NFR-030/D-007 ratios come from the clean install and the breakdown from the profiled one, unless the profiled build's overhead is measured to be under 2% on the same case (median of 3, idle machine).
5. The FR-010 MLMG coarsening WARNING, if any, and the coarsest MLMG grid size.
6. The rules apply to S1 once the AMReX driver exists, not to the A-26 FDS fine-reference runs.

### 5.8 Scope suite (requirement FR-006 from spec v0.4.10)

FR-006 (D-033, A-41) requires the AMReX code to run the FDS Verification inputs that cover in-scope features. The scope is:
- deferred, not dropped (FR-044): GEOM/CC_IBM and HT3D;
- must work with refinement: thin obstructions, HVAC, pressure ZONEs, level-set wildfire.

This section turns FR-006 into a run list and pass rules. The analysis behind it is `scope_alignment.md`.

**Case list.** The FR-006 case list is every row of `docs/vv/scope_case_list.csv` with `cls = IN`. That is 705 of the 941 `Verification/**/*.fds` inputs at FireX `36975d765f` (v0.4.5 rerun). The other classes:

| cls | Count | Meaning |
|---|---|---|
| IN | 705 | All features in scope and runnable on a CPU build here |
| DEFERRED | 175 | Uses GEOM, CC_IBM (MISC or implied by GEOM) or HT3D |
| OUT | 26 | TRN stretching (IR-002/D-030); the four decided embedded/overlapping cases; VTK (D-019); missing external files |
| UNCLEAR | 35 | Needs a decision; reason in the `reason` column |

GPU rule (D-035):
- A case is OUT on GPU grounds only if it cannot run on a CPU build (a keyword in the script's `GPU_ONLY_KEYWORDS`; none at this pin).
- `&PRES HYPRE_DEVICE_RUN` is ignored on CPU builds (pres.f90:1177-1178).
- The ranks-per-GPU layout is the `FDS_RANKS_PER_GPU` job variable, not an input keyword.
- GPU runs of the AMReX code are verified under NFR-043/NFR-047, not here.
- The 9 `GPU_Tests/` inputs are IN at completes-only, run in G4 only (about 2,600 s each, inventory model). The four `_RSn` twins differ from their partners only in CHID.

Other rules for the list:
- `scope_case_list.csv` is generated. It is never edited by hand. Rerun `python3 vv-runs/tools/scope_filter.py` (static parse, about 10 s, no FDS run). The script re-checks its FireX source citations and reports any that no longer match.
- The list goes with the source commit. A new FireX base means a rerun and a diff of the `cls` column, recorded in §10.
- Derived inputs replace originals only where `needs_remesh` names a copy (at present `race_test_1`/`_4` → A-35 `_r4`, AMR mode only). The FDS baseline is always the committed original. The one exception is `stairwell`, whose baseline needs an MPI_PROCESS ≤ 7 derived copy (as `random_obstructions_fft` already has).

**Pass criteria for an IN case.**
1. *Completes:* normal exit at `T_END` with no NaN, trap, `SETUP_STOP` or watchdog stop (§3.4). Setup-only inputs (`T_END=0`) complete setup, and their G1 message diff is clean.
2. *Assigned gate:* the `tier`, `tolerance_class` and `gate` columns of the CSV.
   - Cases in the 184-row inventory keep their inventory tier and class (G2/G3 rules above).
   - Other IN cases are G4 (FR-003, T2 vs the FDS baseline), using the committed Verification-Guide metric where one exists. Ten must-work multi-mesh cases are promoted to T2 in G3.
   - Cases FDS itself does not run (absent from or commented out of `FDS_Cases.sh`) are completes-only.
   - G7 (FR-016), G11 (FR-005) and G13 (output) apply where §5 already assigns them.
   - `STORE_SPECIES_FLUX` cases take no part in the D-022 kernel comparison until A-34 closes (§5.5); their whole-run T2 is unaffected.
3. *Baseline failures:* an IN case that fails on the FireX baseline is recorded, not counted against the AMReX code, and reported to the project owner.
4. *Exceptions:* any other exception (a case dropped, a class loosened) needs the project owner's approval and is recorded in §10.

Milestones:
- Phase 2: every IN case in uniform mode.
- From Phase 3: the refined variants (below).
- Phase 10: complete.

**Refined variants (A-41).** None of the four must-work features has an IN Verification case at a resolution change, and level set has no multi-mesh case (`scope_alignment.md` §4). V&V therefore adds derived refined variants under `vv-runs/inputs/A-41/`, A-35-style: the original stays the FDS baseline, and each variant gets its own FDS baseline. First set:
- thin duct walls across a 2:1 interface (`duct_flow`);
- HVAC vents at a resolution change (`fan_test`);
- a zone breach across a level jump (`zone_shape`, `zone_break_fast`);
- a level-set front crossing a 2:1 interface (`LS4_ember_yield` or `level_set_fuel_model_1`).

These are T2 vs their own FDS baseline and gate from Phase 3.

**Resolving UNCLEAR cases.** Each UNCLEAR row names one reason group. Each group is closed by a single recorded decision (decision-log ID), after which the script rule is updated and rerun, and the cases move to IN, OUT or DEFERRED. Until then they run in the baseline only and do not gate.

| Group | Cases | Decided by |
|---|---|---|
| VARIABLE_THICKNESS vs deferred HT3D | 12 | Project owner / Spec Lead |
| CYLINDRICAL | 11 | Q4 |
| TUNNEL_PRECONDITIONER (includes `dancing_eddies_ulmat_hypre`) | 7 | Q4 |
| ZONE leakage over several meshes | 2 | FR-034 leak-area summation (open with the Chief Architect and Pressure Solver Lead) |
| Thin OBST at a static coarse/fine interface | 1 | FR-040 R3 ruling (open with the Chief Architect and Pressure Solver Lead) |
| Per-mesh CSVF restart files | 1 | Integration Lead |
| Disjoint multi-resolution meshes | 1 | A-35-style remesh |

The three Q4 groups (VARIABLE_THICKNESS, CYLINDRICAL, TUNNEL_PRECONDITIONER; 30 cases) stay UNCLEAR until the project owner answers Q4.

The 40 IN cases with a non-box mesh union (`nonbox_domain`) are IN on the assumption that inter-mesh gaps become solid level-0 cells. If the Integration Lead does not confirm this, they move to UNCLEAR.

**Particle anchors without GEOM (A-47, proposed).** `cascadempi` is DEFERRED, but FR-050/051/052 stay in scope. Proposed anchors (details, metrics and file:line evidence in `scope_alignment.md` §8):
- **FR-052:** derived `VV/cascadempi_obst` (the GEOM half of `cascadempi` removed). It has 3 meshes and 28,672 cells, and box tops are split by the x = 7 mesh face. Its 2:1 variant puts the level jump through the box tops. Metric: committed dataplot rows `OBST-S/OBST-M-water mass`, T2. Secondary: `VV/geom_sprk_mass_obst` (R-7).
- **FR-050:** `bucket_test_1` (4 meshes, spray from the 4-mesh corner), then `VV/cascadempi_obst`.
- **FR-051:** `energy_budget_particles` (unchanged), plus `VV/cascadempi_obst` (evaporation mass balance across meshes).

Each derived copy gets its own FDS baseline and gates from Phase 7.

**DEFERRED cases** stay FDS-only baselines (captured as usual, used for G1 and as FR-004/FR-044 negative-test sources) until the scope changes. With refinement they must stop at setup (FR-004). A-41 notes that `ht3d_energy_conservation_4` may stay in the uniform-mode FR-003 regression; that is outside FR-006. When a deferred feature returns to scope, its cases are reclassified by rerunning the script with the rule removed.

## 6. Case selection (details and per-run estimates in `case_inventory.md`)

**Tier 1 (40 runs, ≈ 8 min sequential model estimate).**
- Interfaces/pressure: `obst_activation_{default,ulmat}`, `divergence_test_2/3`, `dancing_eddies_{1mesh,default,uglmat_refine}`, `obst_coarse_fine_interface`, `lapse_rate`.
- Refinement analogues: `random_meshes`, **`ns2d_16_int_1to2_refinement` + derived `VV/…_uglmat`** (FR-016 gating case, after A-09), **`ns2d_16_emb_1to2_refinement`** (informational, non-gating; FR-014).
- Interface invariance: `soborot_superbee_square_wave_128{,_1mesh}`, `shunn3_4mesh_128`.
- Convergence: `ns2d_{8..64}{,_nupt1}`, `saad_512_cfl_*`, `shunn3_{32,64,128}`.
- Conservation: `energy_budget_tmix`, `species_conservation_1/2`.
- Fire: `1_step_2_step_compare`.
- Restart: `restart_test1a/b` + derived continuous run.
- Determinism pair.

**Tier 2 (123 runs, ≈ 79 min sequential / ≈ 30–40 min packed)** adds the pressure-solver variants, the remaining refinement analogues (`dancing_eddies_embed`, `duct_flow_uglmat_refine`, `ns2d_16_emb_1to1_refinement`, `porous_media`, `race_test_1` (AMR mode: A-35 `_r4` copy)), the full MMS / scalar / limiter convergence series, multi-mesh symmetry and volume-flow cases, pressure zones, mass-balance cases, fires, restarts, and every anchor case not in Tier 1 (`energy_budget_particles`, `zone_break_fast`, `zone_shape`, `cascadempi`, `cannon_ball`, the two radiation anchors, `openmp_test64a`, `device_restart_base_case`). **Optional (8):** `tunnel_demo`, `energy_budget_adiabatic_walls`, `race_test_4` (AMR mode: A-35 `_r4` copy), `shunn3_512`, `ht3d_energy_conservation_4`, `VV/test_8mesh_NOFRPG_short` (HYPRE cost, NFR-033), and the Intel-only `*_uglmat_pardiso` pair.

**Adaptive_Mesh_Refinement folder (scope item; A-09).** `ns2d_16_int_1to2_refinement` (13 meshes / 448 cells, fine 16×16 patch abutting 12 coarse 4×4 meshes, ratio 2) and `ns2d_16_emb_1to2_refinement` (2 meshes / 512 cells, fine mesh embedded in the coarse) are **Tier 1**. `ns2d_16_emb_1to1_refinement` (2 / 320, same-resolution control) is **Tier 2**. All are 2-D, periodic, T_END = 30 s, with `VELOCITY_TOLERANCE=1e-6` and up to 100 pressure iterations. Estimates: ~10–30 s each; the 13-mesh case is 0.5–5 min on 1 rank, depending on iteration counts. Memory is negligible.

*Why commented out:* all three lines (FDS_Cases.sh 4-6) come from commit `042aa2624d` (R. McDermott, 2024-11-08), *"updates to amr test cases, demonstrating errors in both embedded mesh strategy and interpolated refinement"*. That commit also created the `int` case. They were added deliberately disabled, as demonstrations of known error. They are not broken or slow. There is no dataplot row or Guide section, and nothing has changed since. A-09's "pass" is the V&V definition in §3.4 step 1 (`case_inventory.md` §4a).

**NFR-032 plume (Tier 3, Phase 1 S1 and Phase 10; D-008, A-19).**
- Level 0 / coarse (decided 2026-09-25 by Spec Lead, Integration Lead and Chief Architect): the committed `Validation/Heskestad_Flame_Height/FDS_Input_Files/Qs=1_RI=10.fds` (33×33×80) trimmed by one cell on the high x/y side to **32×32×80**, `XB=-1.8,1.690909,-1.8,1.690909,-0.45,8.59` (dx = dy = 3.6/33 and dz = 0.113 m unchanged; AMR level 0 `blocking_factor=8`, `max_grid_size = 16 16 8` per the D-024 amendment: 40 boxes, 5 per rank on 8 ranks). The trim keeps the burner at 9×9 cells; an evenly centred 32-cell grid would snap it to 10×10 (+23 % area, changes Q*). V&V input with the probes: `vv-runs/inputs/A-19/Qs1_RI10_coarse_32x32x80.fds` (81,920 cells), also the ADR-002 S1 input (`CFL_FILE`).
- Uniform-fine reference: `vv-runs/inputs/A-19/Qs1_RI10_fine_64x64x160.fds`, same XB, 655,360 cells in 8 meshes of 64×64×20 (81,920 each) stacked in z. Estimates (**not measured**): ≈ 2.3 h on 8 ranks (range 0.5–4.5 h), ≈ 1.6 GB total (range 1–2.5 GB); coarse ≈ 50 min on 1 rank.
- Burner check (setup-only FDS runs, GNU Release, 2026-09-25): committed 33×33×80 and new 32×32×80 both snap the burner to 9×9 cells (`.smv` OBST indices `12 21 12 21 3 4`), the fine input to 18×18×2 (`24 42 24 42 6 8`); footprint ±0.49091 m, 0.963967 m² on all three, so the area-adjusted HRR (1512.7 kW nominal) and Q* are unchanged. Details, changes, probe verification (`verify_A19_probes.py`) and assumptions: that folder's `README.md`. Superseded v1 (66×66×160 / 33×33×80) is in `vv-runs/inputs/A-19/superseded/`.
- Requirements v0.3.2 records this as D-024. Minor text points for the Spec Lead (A-19 README): `max_grid_size=16` gives 20 level-0 boxes, not 32 (now moot: the D-024 amendment sets `max_grid_size = 16 16 8`, 40 boxes); a 10×10 snap would change burner area and Q*, not total HRR (FDS area-adjusts HRRPUA); memory ≈ 1.4 GB there vs ≈ 1.6 GB here (both estimates).
- The V&V alternative `Fires/circular_burner.fds` (R-8) stays as a fallback only (D-020 open); requirements' own fallback is `McCaffrey_14_kW_5.fds`.

**FireX-only folders.**
- **`GPU_Tests`**: 9 inputs, all 1 M cells, propane fire, UGLMAT/HYPRE, plus HPC job scripts. The physics runs on CPU (`HYPRE_DEVICE_RUN` is compiled out without `WITH_HYPRE_DEVICE`). The job scripts don't apply. The 16/32-mesh variants need more than 8 ranks. Full runs take 20–45 min each, so Tier 3, except one 20-step copy (T2-opt) that serves NFR-033 and FR-038.
- **`VTK`**: 4 CPU inputs demonstrating the VTKHDF writer. Their data output needs an HDF5 build (D-019). They turn off `.smv` output and have no pass/fail criterion, so they are Demonstration only. Relevant to FR-072/074 and Q5.

**Refinement-ratio cross-check (A-33 inventory vs our lists; v0.4.3, D-030).** Sources: `docs/inventory/mesh_ratio_cases.csv` and `mesh_ratio_cases_rollup.csv`, crossed with the 184 rows of `case_inventory.csv` (read-only). 180 are in the rollup; the 4 derived rows are not, and they inherit from their sources.
- Worst class per case: single_mesh 99, uniform 70, supported_2_4 4, no_shared_faces 3, embedded_or_overlapping 2, other_nonpow2 2.
- **Only two of our cases use a ratio outside {1,2,4}:** `Thread_Check/race_test_1` (T2) and `race_test_4` (T2-opt), each with five 5:1 faces. They need the A-35 remesh for AMR mode (done, §5.6); the originals stay the FDS baselines.
- None of our cases has a 3:1, direction-dependent or TRN-stretched interface.
- The embedded/overlapping cases have no shared-face ratio and stay FDS-baseline-only for the multi-level comparison (D-015; FR-014/FR-043 use them as AMR-code checks, not as ratio inputs):
  - `ns2d_16_emb_1to2_refinement` (T1), `ns2d_16_emb_1to1_refinement` (T2);
  - `random_meshes` (T1; worst class uniform, but 8 overlapping pairs at 1:1 and 2:1);
  - `dancing_eddies_embed` (T2; 2:1 overlap).
- Supported 2:1 cases need no action: `ns2d_16_int_1to2_refinement` (2-D, y ignored), `dancing_eddies_uglmat_refine`, `obst_coarse_fine_interface`, `duct_flow_uglmat_refine`.
- Of the 43 requirements §2.3 anchors, only race_test_1/4 are outside {2,4}.

**Shaped-OBST area defect (FM_Burner; FR-005 (v)).** FDS's `AREA_ADJUST` for CYLINDER OBSTs split across meshes depends on the mesh-to-process map (upstream issue candidate #1).
- The 12 OBST-cylinder FM_Burner validation inputs deliver 0.52083 × the specified fuel when all meshes are on 1 process, and exactly 1.00000 with 1 mesh per process (the upstream layout) or on a single mesh.
- None of them is in our lists. They are FR-005 (v) test inputs (§5.3), not FDS baselines, unless run as a single-mesh copy or 1 mesh per process.
- Details: `fm_burner_area_check.md`.

Full assessment: `case_inventory.md` §6. **Excluded, with reasons:** `case_inventory.md` §7 (GEOM, >8 ranks, hour-plus runs, local physics sub-models covered by G4).

## 7. Requirement → test traceability (class used)

| Req | Test / cases | Class | Phase |
|---|---|---|---|
| FR-001 | single-mesh anchors (ns2d_16/32, shunn3_32, species_conservation_1, energy_budget_*, cannon_ball, radiation anchors, openmp_test64a, restart_test1*) | Explicit kernels on frozen input T1 (T0 per Q1); whole runs T2 (D-022) | 2 |
| FR-002 | multi-mesh same-res anchors (shunn3_4mesh_32 first at M2a, layer_4mesh, dancing_eddies_default, symmetry_test_mpi, duct_flow, …) | Kernels on frozen input T1 (T0 per Q1) vs an equivalent single-mesh input (committed `layer_1mesh`; derived copies for shunn3_4mesh_32, dancing_eddies_default, symmetry_test_mpi, A-24); whole runs T2 vs single-mesh and vs multi-mesh FFT and GLMAT baselines; `H` single frozen solve vs GLMAT within **eps_H** (G14) | 2 |
| FR-003 | G4 | T2 | 2, 10 |
| FR-005 | G11 (1/2/4/8 ranks, two `max_grid_size`, two thread counts, 3 repeats; `zone_break_fast`/`zone_shape` for (ii); §5.3 cylinder reproducer, FM_Burner CH4 2cm + single-mesh copy, `zone_shape`/`random_meshes` volumes for (v)) | (i) byte-identical explicit stages; (ii) exact zone sums/gauge from the per-box sum; (iii) eps_H; (iv) bitwise repeat; (v) setup areas/volumes T0 across configurations, T1 vs single-mesh FDS (§5.3) | 2 |
| FR-004, FR-044, IR-002 | negative tests | SETUP_STOP + message | 3 |
| FR-010/013/015 | hierarchy dumps, debug assertions, repeat runs; ratios {2,4} only (MLMG limit ≤ 4; D-030); §5.4 negative tests N-1..N-4 (3:1, direction-dependent, 5:1, legacy), positive controls P-1/P-2, blocking-factor report, MLMG coarsening warning W-1 | exact box lists (T0-equivalent); setup stop + named pair/ratio; warning text + coarsest grid | 3 |
| FR-011 | per-criterion tagging unit cases | exact cell sets | 3 |
| FR-012 | regrid mass/species deltas | 1e-12 rel (req.) | 3 |
| FR-014 | ns2d series + patch, ns2d_16_*_refinement (AMR code only; baseline numbers are context) | refined ≤ coarse, factor 2.0, order ≥ 1.8 | 4 |
| FR-016 | G7 (`int_1to2` + UGLMAT copy only, after A-09; `emb_1to2` informational) | T0 (a, component-level, first ghost layer only, via A-09b; §5.2); (c) FR-032 machine zero + exact-solution error ≤ UGLMAT + order ≥ 1.8, no eps_H; T2 (d); multi-step never bitwise (D-023); HYPRE commit as §3.5 | 3–4 |
| FR-020/021/024 | species_conservation_1..4, simple_duct, mass_balance_* with patch | req. numbers (round-off with refluxing, D-023) | 3 |
| FR-022 | energy_budget_tmix/particles/dns_100/adiabatic_walls with patch | max(1.05·\|closure_base\|, 0.005·max\|Q_TOTAL\|) | 4–5 |
| FR-023, FR-031, FR-032 | divergence_test_1..3, obst_activation, CHECK_POISSON output | req. numbers; FR-031 tolerance TBD(Pressure Lead). eps_H is not a residual criterion | 4 |
| FR-030/033/035 | duct_flow, tunnel_demo, dancing_eddies_default, ns2d periodic; pressure_iteration3d_default, random_obstructions_fft | T2 | 4–5 |
| FR-034 | zone_break_fast, zone_shape (+ zone_shape_2, zone_break_fast_uglmat_hypre) | T2 | 4/6 |
| FR-036 | shunn3_128 (PRESSURE_TOLERANCE 1e-6, analytic H) + helium_2d_isothermal | T2 + iterations ≤ baseline + 1 | 4 |
| FR-037 | FR-001 cases with the FFT path selected; single frozen solve vs baseline FFT (G14, §5.1) | D + eps_H | 2 |
| FR-038 | G6 | T2 (v0.3.1), same HYPRE commit (§3.5) | 2 |
| **FR-039** | G14 / §5.1: FFT vs MLMG single solve on ns2d_16/32, shunn3_32, csmag_32, dancing_eddies_1mesh; switch run `VV/ns2d_32_switch`; maxorder 2 vs 3 study | eps_H (single solve, after gauge); exact switch steps; order ≥ 1.8 | 2 (P2), 4 |
| FR-040/041/042/043 | mask checker; burning-OBST regrid case; obst_activation_default, box_burn_away1; random_meshes, duct_flow | req. / T2 | 5 |
| FR-050/051/052 | particle ledger; energy_budget_particles, cascadempi; FR-052 cascadempi + OBST-only geom_sprk_mass copy | exact counts / T2 | 7 |
| FR-060/061 | plate_view_factor_cart_30, radiating_polygon_square_20 | T2 (060); 061 via FR-001/002: kernels T1, whole runs T2 | 2, 8 |
| FR-070 | all Tier 1 DEVC outputs; regrid-jump check | T2 uniform (D-022); T3 jump | 9 |
| FR-071 | header diff on all Tier 1/2 CSVs | exact (I) | 9 |
| FR-072 | Smokeview load of Tier 1 cases | D, SMV 6.11.2 | 9 |
| FR-074 | VTK/*.fds (HDF5 build; D-019 open) | identical VTKHDF | 9 |
| FR-080/081 | restart_test1a/b + continuous, device_restart_* (+ base case), restart_ulmat_*, csvf_restart_a | FR-080: T0 proposed, T1 minimum; baseline-difference fallback only if baseline misses T1; AMR T1 + identical hierarchy. T2 (081) | 9 |
| IR-001 | G1 | message diff | 2 |
| IR-004 | random_meshes, layer_4mesh (MPI_PROCESS inputs) | behaviour + warning text | 3 |
| IR-006 | G5 | T0 | 2 |
| NFR-010/011/012 | G11 (NFR-011 AMR mode on the A-35 race_test `_r4` copies, FDS baselines on the originals; §5.6) | explicit stages bitwise and `H` within eps_H across ranks (FR-005); whole runs T2 uniform / T3 AMR; watchdog | 2+ |
| NFR-020 | both toolchains build and pass FR-001 | blocked (A-18); GNU FireX baseline exists | 2 |
| NFR-021 | HYPRE/AMReX pins recorded; FireX HYPRE = `63331f19c7` (v2.32.0-24) per provenance | I | 1 |
| NFR-030/031 | openmp_test64a + 8-rank multi-mesh anchor; A-19 memory; timing rules (1)-(6) incl. solver, TinyProfiler pressure time, forced-MLMG level-0 run, separate profiled install (§5.7) | ratio, median of 3, idle machine; clean-install headline | 2, 10 |
| NFR-032 | Heskestad Qs=1_RI=10, level 0 32×32×80 trimmed + A-19 64×64×160 fine reference (D-008) | T3 + ≤ 50% wall | 10 |
| NFR-033 | test_8mesh_NOFRPG 20-step, MPI_PROCESS 1/2/4/8 ranks; secondary duct_flow_uglmat_refine | ≥ 60% efficiency (proposed) | 10 |
| NFR-040 | harness over the anchor set using T0–T3 | D | 2 |

## 8. Tooling (delivered and planned)
- **Delivered:**
  - `parse_fds_inputs.py`: static survey. Follows FDS `CHECKREAD` rules and expands MULT including SKIP ranges.
  - `make_inventory.py` and `render_inventory_md.py`: inventory, cost model, class mapping.
  - `compare_csv.py`: CSV comparison with `--class T0|T1` presets. Implements §2 T0/T1 exactly: timing-column exclusion patterns, T1 column-max scaling, `--exact`, per-column tolerances, interpolation for T2-style series, JSON reports. Tested on synthetic CSVs (`vv-runs/tools/testdata`); first real FDS output comes from the GNU capture.
  - `vv-runs/inputs/A-19/make_A19_inputs.py`: generator for the NFR-032 fine reference (64×64×160) and level-0/coarse input (32×32×80); `verify_A19_probes.py`: independent probe/device/burner-snap checker.
- **Planned (Phase 2, NFR-040):**
  - `run_case.sh` harness: pre-flight load/RAM check, watchdog, manifest (HYPRE version from provenance).
  - numpy re-implementations of the dataplot metrics (`end`, `max`, `mean`, `area`, `end_1_n`, `all`, `tolerance`, `slope`) and of the analytic-error scripts needed by Tier 1/2 (`ns2d.py` RMS error, `shunn_mms.py`, `saad_mms_temporal_error.py`, `soborot_mass_transport.py`, `mass_balance*.py`).
  - A T3 statistics script (time-window mean/RMS at probes; A-19 probe naming `T_z…`/`W_z…`).
  - A frozen-input `H` comparator (mean removal, relative L2, eps_H(N)) and a solver-selection log checker (G14).
  - A hierarchy-dump comparator once the dump format exists.

---

## 9. Status of the V&V items folded into requirements v0.3

The items below are this plan's v0.2 answers (R-n) and disagreements (D-n). Their rationale is in `archive/test-plan-v0.2.md` §9. Status is as recorded in requirements.md v0.3 and the `docs/README.md` decision log. v0.3.1 (D-022) changes the status of D-2 and D-3 (below). All accepted numbers are **provisional until calibrated** against the measured FireX baseline (A-08/A-10) and approved by the project owner.

### 9.1 Answers (R-n)

| # | Topic | v0.3 status | Where in v0.3 | Remaining V&V action |
|---|---|---|---|---|
| R-0 | Anchor runtimes, rank caveats | accepted | §2.3 runtime check; D-013 | Replace estimates with measured times |
| R-1 | T2 margin (floor, 3σ) | accepted | §2.2 T2; D-013 | σ ensemble (needs Intel) |
| R-2 | T3 values, discrimination rule | accepted | §2.2 T3, NFR-032; D-013 | Calibrate at Phase 4 |
| R-3 | FR-003 at two frequencies | accepted | FR-003; D-013 | — |
| R-4 | FR-014 factor 2.0 | accepted | FR-014; D-013 | — |
| R-5 | FR-022 floor | accepted | FR-022; D-013 | — |
| R-6 | FR-036 case shunn3_128 (+ helium_2d_isothermal) | accepted | FR-036; D-013 | — |
| R-7 | FR-052 case cascadempi (+ OBST-only geom_sprk_mass copy) | accepted | FR-052; D-013 | Create derived copy |
| R-8 | NFR-032 case circular_burner | **open, not adopted** | NFR-032 keeps Heskestad `Qs=1_RI=10` (D-008); **D-020** | A-19 inputs delivered (§6); alternative kept until A-19 / ADR-002 S1 report |
| R-9 | Smokeview 6.11.2 | accepted | FR-072; D-013 | Install needs authority over the development machine |
| R-10 | NFR-011 watchdog | accepted | NFR-011; D-013 | Implement in `run_case.sh` |
| R-11 | FR-080 restart rule | **adjusted** | FR-080: T0 proposed, **T1 minimum**; "no worse than baseline's own restart difference" only if baseline misses T1; D-013 | Measure baseline restart-vs-continuous at capture (G10) |
| R-12 | NFR-033 scaling case | accepted | NFR-033; D-013 | Create 20-step MPI_PROCESS copies |

### 9.2 Disagreements and clarifications (D-n)

| # | Topic | v0.3 status | Where in v0.3 | Remaining V&V action |
|---|---|---|---|---|
| D-1 | T0/T1 on `SIG_FIGS=17` derived copies | accepted | §2.2 T0/T1; **D-014** | Copies under `vv-runs/` |
| D-2 | Chaotic FR-001 anchors: T1 100 steps + T2 | **open in the decision log (D-017, Q1)**; **superseded in v0.3.1 text**: FR-001 whole runs are T2 for every anchor, and kernels are T1/T0 per Q1 (D-022) | FR-001; D-017 (log not yet updated) | Spec Lead to close D-017 or restate what Q1 still decides (kernel T0 vs T1) |
| D-3 | FR-002: T1 first 10 steps, then solver-tolerance bound | **open in the decision log (D-018)**; **resolved in v0.3.1 text**: full-run T1 vs GLMAT dropped ("resolves V&V's objection"), whole runs T2, kernels vs single-mesh T1/T0, `H` eps_H (FR-002, D-022) | FR-002; D-018 (log not yet updated) | Spec Lead to close D-018; FR-031 tolerance still needed for FR-023/FR-032 |
| D-4 | FR-016 = `int_1to2` only; `emb_1to2` informational | accepted | FR-016; **D-015** | — |
| D-5 | UGLMAT reference copy; instrumented baseline | accepted | FR-016, A-09 (UGLMAT copy run first), **A-09b** (instrumented build), R-34; D-015 | §3.4 step 1, §3.6 |
| D-6 | GPU_Tests inputs run on CPU | accepted | §2.3 note; D-013 | — |
| D-7 | VTKHDF baseline needs HDF5 build | **open** | **D-019**, depends on A-14/Q5; one option (HDF5 option in `AMReX` branch) conflicts with D-005 | No HDF5 build until decided |
| D-8 | Ranks = min(8, meshes) | accepted | §2.3, NFR-010/011; **D-016** | — |
| D-9 | T2 absolute floor / chaotic margin | accepted | §2.2 T2; D-013 | — |
| D-10 | FR-022 floor | accepted | FR-022; D-013 | — |
| D-11 | FR-080 T0 vs baseline | **adjusted** | as R-11 | as R-11 |
| D-12 | Timing valid only on idle machine | accepted | NFR-030/033; D-013 | — |
| D-13 | `device_restart_base_case` anchor; continuous restart companion | accepted | §2.3; D-013 | Create continuous run |

## 10. Open questions and actions for V&V
- **Blocking:**
  - (1) oneAPI runtime on the development machine (A-18, R-32). The OpenMPI runtime is fine (§3.1).
  - (2) Intel FireX build per §3.3 (Intel Build Chief), with the **same HYPRE commit `63331f19c7`** as GNU.
  - (3) Smokeview install (FR-072).
- **Needs a decision from others:**
  - Q1: now only kernel-level T0 vs T1 (FR-001/002). D-017/D-018 are superseded by D-022 in the v0.3.1 text; the decision log still lists them open.
  - Pressure Lead's FR-031 tolerance → FR-023, FR-032.
  - Kernel-parity dumps from baseline FDS (FR-001/FR-002 under D-022): extend the A-09b instrumented build or specify a separate one (Integration Lead, Pressure Lead, Legacy Mapper). A-24 derived single-mesh copies: owner not stated in the requirements; V&V can make them.
  - Q5 and A-14 → D-019 (VTK/HDF5 baseline).
  - D-020: Spec Lead, after the A-19/S1 report.
  - Q9 (MPI_PROCESS in AMR mode).
  - Spec Lead: HYPRE label in requirements §2.1/NFR-021 and R-30 ("3.0.0" vs the linked v2.32.0-24), and whether the AMR build keeps `63331f19c7` (§3.5).
  - Spec Lead/project owner: whether A-09 on the GNU Release + Debug builds is enough to start treating FR-016 as a gate while the Intel build is blocked.
  - ~~Integration Lead / Chief Architect: the D-008 level-0 grid 33×33×80 is odd in x/y~~ **Resolved 2026-09-25 (D-024, requirements v0.3.2)**: level 0 = 32×32×80 trimmed on the high side, `blocking_factor=8`, `max_grid_size = 16 16 8` (D-024 amendment, 40 boxes; §6).
  - Integration Lead: frozen-input solve hook and solver log for G14 (§5.1).
  - (v0.4.3) Spec Lead: FR-005 (v) says the exact setup sums "reproduce single-mesh FDS". V&V applies T0 across AMR configurations and T1 against single-mesh FDS, because FDS's sequential float sum need not equal the exact sum rounded once (§5.3). Confirm or tighten.
  - (v0.4.3) Integration Lead: a setup diagnostic that prints `FDS_AREA`/`AREA_ADJUST` per shaped OBST/MULT and the zone volumes, also in setup-only mode. FDS computes neither in setup-only runs (§5.3).
  - (v0.4.3) Integration Lead / input converter: the FR-010 ratio check must ignore the single y cell of 2-D inputs, or `ns2d_16_int_1to2_refinement` (FR-016) is rejected (§5.4 P-1).
  - (v0.4.3) Chief Architect, for the record: `amr.blocking_factor = 2 8` is valid only with a ratio-4 level 1 (AMReX needs bf0·ratio ≥ bf1: 2·2 < 8 fails). A ratio-2 hierarchy on the same level 0 needs `2 4` or `2` (§5.6; A-35 README).
- **V&V actions (from the requirements action log):**

| ID | Action | Status |
|---|---|---|
| A-08 | Record the FireX builds; capture baselines | GNU build recorded (§3.5); G0 + Tier 1 capture running (`baseline_status.md`); Intel blocked (A-18) |
| A-09 | Run `int_1to2` + UGLMAT copy first (Release + Debug, both toolchains), then the `emb` cases; FR-016 gates only after (R-34) | open; GNU Release/Debug binaries available |
| **A-09b** | Instrumented ghost-cell / coarse-fine dump build from a scratch copy, never committed; accepted only if T0 vs uninstrumented build on all Tier 1 CSVs (with GNU Build Chief) | open; placeholder spec §3.6; detailed spec needs Pressure Solver Lead + AMReX Integration Lead input |
| A-10 | Calibrate provisional numbers | open until measured baseline |
| A-14 | VTK/HDF5 decision support | open (D-019) |
| **A-19** | NFR-032 uniform-fine reference `Qs=1_RI=10` (with Pressure Solver Lead); ADR-002 S1 with `CFL_FILE`; uniform-fine runtime/memory | **v2 inputs generated and setup-checked (T_END=0), no time-stepping run**: `vv-runs/inputs/A-19/Qs1_RI10_coarse_32x32x80.fds` (level 0 / S1, 81,920 cells) and `Qs1_RI10_fine_64x64x160.fds` (655,360 cells, 8 × 64×64×20), trimmed XB `-1.8,1.690909,-1.8,1.690909,-0.45,8.59`; burner 9×9 / 18×18 cells, 0.963967 m² on committed, coarse and fine (§6). Estimates: fine ≈ 2.3 h on 8 ranks, ≈ 1.6 GB. v1 66×66×160 in `superseded/`. S1 and measurements wait for run authorisation and the A-26 exclusive window (D-024; path and grid now match requirements v0.3.2) |
| **A-35** | Remesh `race_test_1` (and hand-check/remesh `race_test_4`) for NFR-011 in AMR mode (D-030) | **inputs done 2026-09-25**, not yet run in AMR mode (needs M3). Ruling confirmed: mesh 3 → `IJK=24,24,16`, all 5 interfaces 4:1, faces tile. race_test_4 has identical meshes, so it gets the identical remesh (option A of 4). Level 0 34×18×32, `blocking_factor 2 8` or `2`. Setup-checked (28,656 cells vs 37,440). Side effects: burner top/vent 0.120 → 0.125 m, 32 → 16 burning cells with fuel flow preserved, probe cell moved. `vv-runs/inputs/A-35/README.md`. Next: FDS runs of the `_r4` copies at 1 and 4 threads as the AMR-mode reference (short runs only until run authorisation). |
| A-33 (input) | Legacy Mapper ratio scan, cross-checked against our lists | closed by the Legacy Mapper. V&V cross-check done (§6): only race_test_1/4 are outside {2,4}. |
| FM_Burner (R/FR-005 v) | Area-adjust check of the FM_Burner inputs | done: `fm_burner_area_check.md`. 0.52083 × fuel with all meshes on 1 process, 1.00000 with 1 mesh per process or single mesh. Not in our lists; used as an FR-005 (v) input. |

- **Also from v0.3.1:** A-24 (derived single-mesh copies of `shunn3_4mesh_32`, `dancing_eddies_default`, `symmetry_test_mpi`) and FR-005 verification (V&V is verifier for (i)-(iv)); A-25 (global-reduction list) is not V&V's.
- **V&V next steps:** GNU G0 → G1 sweep → A-09 (GNU Release + Debug) → Tier 1 capture (GNU) → Tier 2 → Intel once A-18 is done → calibration ensemble → replace model estimates with measurements and re-tier.

## 11. Files
- `docs/vv/test-plan.md` (this file, v0.4.3); `docs/vv/archive/test-plan-v0.3.md` (previous version); `docs/vv/archive/test-plan-v0.2.md` (full §9 rationale)
- `docs/vv/fm_burner_area_check.md` (FM_Burner / FR-005 (v) area check); work files in `vv-runs/inputs/fm_burner_area/`
- `vv-runs/inputs/A-35/{race_test_1_r4.fds, race_test_4_r4.fds, check_mesh_interfaces.py, rt*_check.txt, README.md, setup_check/}`
- `docs/vv/environment.md` (v0.2 content; still says no FireX binary exists and HYPRE 3.0.0. To be updated from §3.5)
- `docs/vv/case_inventory.md` and `case_inventory.csv` (v0.2 content; e.g. `emb_1to2` still marked "FR-016 disputed", no A-19 rows)
- `docs/vv/verification_case_survey.csv` (FireX) and `verification_case_survey_ce1f659.csv` (v0.1)
- `vv-runs/tools/{parse_fds_inputs.py, make_inventory.py, render_inventory_md.py, compare_csv.py, testdata/}`
- `vv-runs/smoke/{run_smoke.sh, inputs/, gnu/, intel/}`
- `vv-runs/inputs/A-19/{Qs1_RI10_fine_64x64x160.fds, Qs1_RI10_coarse_32x32x80.fds, make_A19_inputs.py, verify_A19_probes.py, README.md, setup_check/, superseded/}` (superseded v1 66×66×160 / 33×33×80 files in `superseded/`)
- Baseline status: `baseline_status.md` (maintained by the capture worker, not by this plan)
