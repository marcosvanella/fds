# FDS-AMR V&V — Environment record

*Snapshot note: sections 4 (runtime environment search), 5 (smoke tests) and 7 (runtime actions) of the working document described the state of the shared development machine on 2026-09-25 and are omitted here. The toolchain and build facts needed to reproduce results are kept.*

Status: **DRAFT v0.2**, captured 2026-09-25 and updated for the FireX re-base. Everything below is **measured on the development machine** unless marked *(assumed)* or *(inferred)*.

## 1. Source

### 1a. Refactor base and acceptance baseline: FireX `36975d765f` (NEW)

| Item | Value |
|---|---|
| Worktree | this repository (read-only for V&V) |
| HEAD | `36975d765f` = `36975d765fcead401e14b094a04f910ac42eab8a` (`git -C (repo root) rev-parse --short HEAD`), local branch `AMReX` |
| Commit | "Merge pull request #16596 from cxp484/FireX", 2026-09-24 15:38:07 -0400. `git describe --tags` = `FDS-6.11.1-1244-g36975d765f` |
| Remote | `origin https://github.com/firemodels/fds.git` |
| Relation to ce1f659 | `ce1f659` is an ancestor. `git diff --stat ce1f659 36975d765f`: 72 files, about +11.5k/−1.0k lines (GPU HYPRE, VTK/HDF5 writer, radiation/init changes, 13 new Verification inputs, 4 inputs gained `HYPRE_DEVICE_RUN=T`). Details: `case_inventory.md` §1a |
| Working tree | Clean (`git status --short` printed nothing) |
| Build status | **No FireX binary exists.** Nothing was built (that is outside V&V scope). "Baseline FDS" for acceptance tests is a FireX 36975d765f build; see test-plan.md §3 for the request to the Build Chiefs |

### 1b. Existing binaries' source: master `ce1f659` (smoke tests only)

| Item | Value |
|---|---|
| Repo | `(local FDS master checkout)` (read-only for V&V: nothing built, run or written there) |
| HEAD | `ce1f659` = `ce1f659cd419aa4f7bb6f60c77e3b009ebbfe824` (`git -C (local FDS master checkout) rev-parse --short HEAD`) |
| Commit | "Merge pull request #16530 from mcgratta/master", 2026-08-31 17:46:24 -0400 |
| Working tree | Not clean, and it was already like this before V&V touched anything: `git status` shows 5 deleted files under `Manuals/Sphinx_Online_Docs/build/...` and the two binaries below as untracked. These changes don't affect the Fortran sources. |

## 2. Binaries (built from ce1f659; used ONLY for toolchain/MPI/tooling smoke tests, never for references)

| | GNU / OpenMPI | Intel / Intel MPI |
|---|---|---|
| Path | `(local FDS master checkout)/Build/ompi_gnu_linux/fds_ompi_gnu_linux` | `(local FDS master checkout)/Build/impi_intel_linux/fds_impi_intel_linux` |
| sha256 | `6d60b02c0242dab3b942aeeadf9899613d65262fe9787fdd475db440ab8f3804` | `325148b9e17925aad59b47779f6a394550f9c42e0a59e3a3dfc1a3c642ad57a5` |
| Size / mtime | 12,704,088 B / 2026-09-23 | 149,717,944 B / 2026-09-23 |
| Compiler (from `.comment` and strings) | GNU Fortran 14.2.0 (Debian 14.2.0-19) | Intel Fortran (ifx) 2026.1.1 ("Intel(R) Fortran 26.0-1156"); icx 2026.1.0/2026.1.1 also present |
| Make flags (`Build/makefile` target) | `-O3 -std=f2018 -frecursive -ffpe-summary=none -fall-intrinsics` + `-fopenmp` (OpenMP enabled, libgomp linked) | `-O2 -ipo -no-wrap-margin -DUSE_IFPORT` (no OpenMP) |
| MPI | OpenMPI, built with `/usr/bin/mpicc` and `/usr/bin/mpifort` (from `firemodels-gnu/sundials/BUILDDIR/CMakeCache.txt`). NEEDED libs: `libmpi.so.40`, `libmpi_usempif08.so.40`, which means OpenMPI 4.x/5.x ABI. Exact version unknown *(inferred)*. | Intel MPI 2021.18: RUNPATH `/opt/intel/oneapi/mpi/2021.18/lib`, built with `mpiicx`/`mpiifx` (from CMakeCache). NEEDED: `libmpi.so.12`, `libmpifort.so.12`, `libimf.so` |
| Third-party (static) | HYPRE v3.0.0 and SUNDIALS v7.5.0 (`(local GNU third-party library tree)/libs`). **No MKL** (0 `pardiso`/`mkl_` symbols). | HYPRE v3.0.0 and SUNDIALS v7.5.0 (`(local Intel third-party library tree)/libs`). **MKL statically linked** (~25,900 MKL/PARDISO symbols). |
| Embedded revision | Empty. The `Revision`/`Revision Date` strings are blank, so the binary can't identify its own commit. Provenance comes only from the file location and mtime. | Same, empty |

Consequences for testing:
- **GNU build, no MKL:** `SOLVER='UGLMAT'`, `'GLMAT'` and `'ULMAT'` fall back to HYPRE (`func.f90` `DEFINE_PRES_METHOD`, `#ifndef WITH_MKL`). Cases that request `... PARDISO` explicitly (`dancing_eddies_uglmat_pardiso`, `duct_flow_uglmat_pardiso`) are only meaningful on the Intel binary. On GNU they are expected to fail or misbehave; this is not yet verified.
- **GNU build is OpenMP-enabled.** Always export `OMP_NUM_THREADS=1` unless a case needs threads on purpose (e.g. `race_test_4`).

## 3. Development machine

| Item | Value |
|---|---|
| OS | Debian GNU/Linux 13 (trixie), kernel 6.12.94+, glibc 2.41 |
| CPU | 8 × "Intel(R) Xeon(R) Processor" (1 socket, 1 thread/core, 1 NUMA node, AVX2 + AVX-512F) |
| RAM | 16 GB total, **no swap**. Shared with other workloads: 9.5–13 GB was in use and **2.4–6.4 GB was available** |
| Python | 3.13.5, numpy 2.2.4. **No pandas, scipy or matplotlib**, so FDS's own `Utilities/Python` plotting scripts can't run here as-is |

## 6. Intended run commands (to use once runtimes exist; NOT yet validated)

```bash
# common
export OMP_NUM_THREADS=1
cat /proc/loadavg; free -m          # require: load < 2 and > (ranks x 0.5 GB + 1 GB) available

# GNU / OpenMPI (needs the OpenMPI runtime that provides libmpi.so.40 + libmpi_usempif08.so.40, e.g. Debian 13 openmpi-bin/libopenmpi-dev)
FDS=(local FDS master checkout)/Build/ompi_gnu_linux/fds_ompi_gnu_linux
$FDS -V                                               # version
$FDS case.fds                                         # 1 rank (singleton MPI init; assumed to work with OpenMPI, else use mpirun -np 1)
mpirun --bind-to none -np N $FDS case.fds             # N ranks (N <= meshes; N==meshes unless MPI_PROCESS is set)
# (--bind-to none because the development machine is shared; revisit binding for timing studies)

# Intel / Intel MPI (needs oneAPI compiler runtime 2026.1 (libimf) + Intel MPI 2021.18 under /opt/intel/oneapi)
source /opt/intel/oneapi/setvars.sh
FDS=(local FDS master checkout)/Build/impi_intel_linux/fds_impi_intel_linux
$FDS -V
mpiexec -n N $FDS case.fds
```

## 8. Smokeview (FR-072 version question)

- **No Smokeview on the development machine**: `which smokeview smv` finds nothing, and neither does a filesystem search outside the source trees.
- FireX `36975d765f` describes as `FDS-6.11.1-1244-g…`. The matching public release bundle is **FDS-6.11.1_SMV-6.11.2** (GitHub release, published 2026-07-10). The standalone Smokeview release at that time was SMV-6.11.1 (2026-05-27). Source: firemodels release pages / pages.nist.gov/fds-smv/downloads.html, checked 2026-09-25. FireX still writes Smoke3D format version 1 by default (`SMOKE3D_VERSION=1`, cons.f90:178); a test flag for format 2 was added 2026-09-09 (`e169e63d44`).
- Proposal: pin **SMV 6.11.2** (from the FDS-6.11.1 bundle) as the FR-072 reference viewer. If FireX-only output fails to load in it, fall back to the Smokeview nightly test bundle built from the same date. The installed version must be recorded here once someone with authority installs it (nothing was installed).

## 9. Toolchain facts needed for the FireX build request (read from Build, not built)

- Makefile targets: `ompi_gnu_linux` (gfortran, `-O3 … -fopenmp`, output `fds_ompi_gnu_linux`) and `impi_intel_linux` (ifx, `-O2 -ipo`, no OpenMP, output `fds_impi_intel_linux`). Debug variants: `ompi_gnu_linux_db`, `impi_intel_linux_db`. Each is driven by `Build/<target>/make_fds.sh`, which sources `Build/Scripts/build_thirdparty_libs.sh` and then runs `make -f ../makefile <target>`.
- Third-party versions in the FireX makefile: HYPRE `v3.0.0`, SUNDIALS `v7.5.0` (both the same as the existing binaries), HDF5 `v1.14.5` (**new, optional**: built only if `$FIREMODELS/hdf5` is a clone at tag `hdf5_1.14.5`; otherwise the build silently goes ahead without `WITH_HDF5` and VTK output writes nothing). The existing third-party trees `(local GNU third-party library tree)/libs` and `(local Intel third-party library tree)/libs` contain hypre and sundials, but **no hdf5**.
- CMake (`CMakeLists.txt`): HYPRE fetched at commit `63331f19` (needs ≥ 2.32.0; this differs from the makefile's v3.0.0). No HDF5 option. New GPU options `USE_HYPRE_{NVIDIA,AMDGPU,INTELGPU}` default OFF. The makefile embeds `GITHASH_PP`/`GITDATE_PP`/`BUILDDATE_PP` from git; the existing binaries have these fields empty, which suggests they were built outside a git worktree. The FireX build should fill them (A-08 provenance).
