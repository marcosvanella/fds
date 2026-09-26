# ADR-001: Driver architecture — C++ AMReX driver with FDS Fortran kernels vs FDS as a Fortran program on AMReX F_Interfaces

| Field | Value |
|---|---|
| Status | Proposed. **Driver choice (Option A, C++ AmrCore driver) is owner-confirmed via D-027**; kernel style (K1/K2) open until the P1 readability review (NFR-044); K2 redefined as Fortran + OpenMP `target` offload (v0.3) |
| Version | v0.3.7 (2026-09-25): owner decisions: NVIDIA is the only GPU target, AMD out of scope; second-compiler check optional (ifx only); AMR mode uniform grids per level, stretched cases FDS-only; see "v0.3.7 changes". v0.3.6 (2026-09-25): pressure section reduced to pointers to ADR-002 v0.2 (maxorder 2 confirmed by P2); HYPRE per D-026; see "v0.3.6 changes". v0.3.5 (2026-09-25): D-031 pass order: one-species case, gate and sync count stated; see "v0.3.5 changes". v0.3.4 (2026-09-25): D-031 rulings: two-phase density terms on valid+2; clip flags only from uncovered valid cells; see "v0.3.4 changes". v0.3.3 (2026-09-25): D-031 pass order with two host OR reductions, coarse-side mask, A-37 closed (47 checks pass), two-phase target form; see "v0.3.3 changes". v0.3.2 (2026-09-25): D-031 ghost-depth ruling (redundant density clip over valid+1, one pre-clip `FillBoundary`); see "v0.3.2 changes". v0.3.1 (2026-09-25): species/density clipping rewritten as a layout-independent gather (D-031). v0.3 (2026-09-25): owner decision changes K2 from OpenACC to Fortran + OpenMP `target` offload, with a no-copy device-data rule; see "v0.3 changes". v0.2 (2026-09-25): owner decisions D-027 and NVIDIA-only GPU target |
| Type | Full ADR |
| Date | 2026-09-25 |
| Owner | AMR Chief Architect (for the project owner) |
| Deciders | Project owner; AMReX Integration Lead; FDS Legacy Mapper; AMR Spec & Program Lead |
| **Code base** | **FireX `36975d765f` (2026-09-24), local branch `AMReX`, this repository (read-only).** Acceptance "baseline FDS" is the same commit (Spec Lead). |
| Other evidence | AMReX `99ddfda` (shallow, CHANGES.md head 26.09); IAMR `0a63b79`, incflo `878827d`, PeleLMeX `5c21556`, ERF `78707d3`, amrex-tutorials `1a73f32`. FDS master `ce1f659` (FireX's merge-base) is quoted only where teammates cited it. |
| Depends on | ADR-003, but only partially (see "Pressure path" below); linked in the decision log, R-03 and FR-030 |

All counts were run with `rg`/`wc` on `Source` on 2026-09-25. "Refs" are textual occurrences. Teammate inventory numbers (`docs/inventory/`, `docs/amrex/mapping.md`, `docs/pressure/`) were produced against `ce1f659`. FireX changed 22 source files (+8,621/−847 lines, `git diff ce1f659 --stat -- Source`), mostly `vtkf.f90` (new), `dump.f90`, `main.f90`, `pres.f90`, `read.f90` and `radi.f90`, so their line numbers in those files have moved. Where both are quoted, the FireX figure comes first.

Cross-references: risks R-02, R-03, R-08, R-11, R-12, R-18, R-21, R-22, R-23, R-26, R-29, R-30, R-31, R-36, R-38, R-39; requirements FR-001/002, FR-005, FR-010, FR-030, FR-037, FR-039, FR-050..052, FR-072/073, IR-002, IR-004, IR-005, IR-006, IR-007, NFR-012, NFR-030/031, NFR-043, NFR-044; owner questions Q2 (answered by D-027), Q5, Q8, Q9, Q11; decisions D-012, D-021, D-027 (supersedes D-004), D-028, D-029, D-031; assumption A-37; roadmap Phase 11 / M11.

## v0.3.7 changes (owner decisions, 2026-09-25)
- **NVIDIA is the only GPU target; AMD is out of scope, not deferred.** "Decision needed: AMD deferred or dropped?" is closed as decided. AMD/HIP/rocFFT planning text is removed or marked "out of scope per owner decision 2026-09-25".
- The second-compiler kernel compile check becomes **optional and non-gating**: `ifx` is kept (already installed), `amdflang` is dropped.
- K1 is still described as portable in principle, but portability is no longer a decision criterion for K1 vs K2.
- **AMR mode uses only uniform grids on each level;** the 59 stretched-grid cases stay FDS-only (confirms D-030). Charter Q11 (c) (x/y-stretched meshes on the GPU) and the "stretched grids as hard requirement" question are closed; details in ADR-002 v0.2.1.

## v0.3.6 changes (2026-09-25)
- "Pressure path and global reductions": solver choice, per-step selection, gauge/eps_H, MLMG order, ratio cap and global reductions now point to ADR-002 v0.2 instead of duplicating it; the open maxorder fallback is replaced by "maxorder 2 confirmed by P2". Only the GPU-specific points stay here.
- HYPRE per D-026: FireX pins v2.32.0-24 `63331f19c`; AMReX built against it at `(local AMReX install built against HYPRE 2.32)`; CPU-only; used only as MLMG's bottom solver. Pressure Lead's P2 open item marked done.

## v0.3.5 changes (2026-09-25)
- D-031 pass order: `N_TRACKED_SPECIES==1` case (one host reduction); renormalisation gate, species-stage `RHOP` ceiling in ghosts and the at-most-2-syncs rule stated explicitly, with prototype line references. The prototype currently issues 2 + NS separate reductions; production packs them.

## v0.3.4 changes (2026-09-25)
- Ruling: two-phase density terms on valid+2 accepted; the input ghost depth (ng=3 for `RHOP`) still suffices.
- Ruling: clip flags are set only from valid cells not covered by a finer level (per-level fine-covered mask, like the ghost exclusion); clipping covered cells is allowed but may be skipped. Flags stay OR-reduced per level. New acceptance check for the two-phase build: a two-level case where only a covered coarse cell clips must not renormalise the uncovered cells.

## v0.3.3 changes (2026-09-25)
- D-031 pass order corrected: two host OR reductions (density flags, then the per-species flag vector), not one OR over all flags.
- Coarse-side mask (Legacy Mapper): coarse cells next to covered cells neither gather from nor push into them.
- A-37 closed (Integration Lead, `docs/amrex/p1-findings.md` §13, `prototypes/p1_mass_shim/check_clip.sh`, 47 checks pass): bitwise across box splits and rank counts, bitwise vs single-mesh FDS, both ghost variants bitwise equal, flags from valid cells only, two-level runs bitwise within each level.
- Cost recorded (unoptimised gather +25-40 % loop time); target form is the §13 two-phase version, accepted if it passes the same 47 checks. Native layout: density MultiFab ng=3, species ng=2; the ng=3 `RHOP` temporary is shim-only.

## v0.3.2 changes (2026-09-25)
- D-031 ghost-depth ruling: the density clip+apply runs redundantly over the grown tile (valid+1), so a single pre-clip `FillBoundary` (ng=3 for `RHOP`, ng=2 for `RHO_ZZ`) replaces the second exchange; the face mask is defined on valid+2. The second `FillBoundary` becomes a rejected alternative. Stencil checked against FireX `mass.f90:799-922`; the ghost counts stand. Added three conditions the stencil implies: `SOLID` on valid+2, metrics on valid+3, clip flags from valid cells only.

## v0.3.1 changes (2026-09-25)
- New subsection "Species/density clipping (D-031)" under the kernel interface: P1 finding (Integration Lead, `docs/amrex/p1-findings.md` §6.2) that `CHECK_MASS_DENSITY` is box-layout-dependent when clipping triggers; decision: gather rewrite with bitwise parity to single-mesh FDS, no exemption or tolerance; explicit face-mask argument; host-side flag reduction between kernel passes.
- New kernel-interface rule (both styles): domain-wide reductions happen on the host between kernel passes, never inside a kernel.
- Line numbers from the Legacy Mapper, checked against FireX `mass.f90`; the final-step `DT` nudge is cited at FireX `main.f90:719` (line 702 in the handed-over note is `DIAGNOSTICS = .FALSE.`; `ce1f659` has the nudge at `main.f90:622`).

## v0.3 changes (owner decision, 2026-09-25)
- **K2 changes from OpenACC to Fortran with simple OpenMP `target` offload.** Owner rationale: FDS already uses only simple OpenMP constructs (356 `!$OMP` lines, 0 `!$ACC`; see "Threading" and "FireX GPU work"), so staying in OpenMP helps readability for FDS's Fortran-only developers and keeps portability (OpenMP offload is supported by nvfortran, AMD amdflang/Cray ftn and Intel ifx). CUDA Fortran may be used for profiled hot loops.
- K2 is restricted to `!$omp target teams loop` plus `collapse`; data via `has_device_addr`; `map` only for host scalars. Kernels must also compile with a second offload compiler (amdflang or ifx) as a portability check.
- **No-copy device-data rule** for K2 (new subsection "Kernel interface and device data"); AMReX managed memory rejected for production.
- CUDA Fortran is allowed only for hot loops with a profiled gap on real NVIDIA hardware; each keeps its OpenMP version as fallback and correctness reference.
- "K2 locks the kernel layer to NVIDIA" is withdrawn: with OpenMP, only the optional CUDA Fortran kernels are NVIDIA-only. AMD stays deferred; whether deferred or dropped is a new owner question.
- OpenACC moves to "Rejected alternatives". The kernel-interface rule "no `!$omp` inside device kernels" is reworded for K2; the Spec & Program Lead owns the requirement text (IR-007, NFR-044) and has been asked to update it. S4/S5 now build the OpenMP-offload Fortran variant (still blocked on A-31).

## v0.2 changes (owner decisions, 2026-09-25)
- **Full time step on the GPU; only I/O may stay on the host** (D-027, supersedes D-004; NFR-043; charter O7). On the development machine NFR-043 is verifiable only by a real GPU-backend compile plus a host-fallback run (D-029, R-38); on-device acceptance waits for test hardware (roadmap Phase 11, M11).
- **GPU target: NVIDIA only for now; AMD is deferred and not a current requirement** (owner answer to charter Q11 (a); Q11 (b) test hardware and (c) x/y-stretched meshes stay open). The spec docs still list Q11 as open and name ifx + SYCL as a compile route (NFR-043, D-029, A-30); updating them is the Spec & Program Lead's call.
- **FDS developers are Fortran-only, and C++ maintainability is a stated concern** (D-027; NFR-044; charter O8).
- Consequences here: Option B rejected definitively; the `POINT_TO_BOX` shim becomes a CPU-only stepping stone with mandatory kernel extraction (R-26 now governs timing only); new sections "Kernel implementation style" (K1 vs K2, co-equal, decided by the P1 readability review), "Layering for maintainability" and "Pressure path and global reductions"; spike plan and owner questions updated.

## Context

### What the driver must host (FireX)
- **Size.** 34 `.f90` files (35 entries in `Source/` including `README.md`), **180,840 lines**. Largest files: `geom` 27,735; `ccib` 24,046; `read` 17,466; `dump` 13,276; `prop` 8,935; `rcal` 7,775; `pois` 7,474; `func` 7,221; `pres` 6,114; `init` 5,531; `radi` 5,323; `hvac` 5,171; `main` 5,149; `part` 5,015; `vtkf` 4,470 (new in FireX).
- **Global per-mesh state.** `TYPE MESH_TYPE` spans `mesh.f90:16-354`. All meshes sit in one global `TYPE (MESH_TYPE), SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: MESHES` (`mesh.f90:356`). The Mapper counts 431 `MESH_TYPE` fields and 1,314 components across the 37 reachable types (`docs/inventory/mesh_fields.csv`, ce1f659). FireX adds 8 slice/VTK fields.
- **Pointer remapping.**
  - `POINT_TO_MESH(NM)` (`mesh.f90:505-916`, module `MESH_POINTERS` at `mesh.f90:361`) makes **402** `=> M%` associations. There were 393 at ce1f659, which matches the Mapper's "393 of 396 module pointers".
  - It is called at **214** sites: 216 `rg` matches minus 2 comments (`main.f90:4950`, `ccib.f90:2602`), across 19 files. By file: ccib 54, pres 33, geom 30, turb 15, read 11, velo 11, dump 10, fire 9, vtkf 9, init 7, part 5, soot 4, vege 4, divg 3, radi 3, mass 2, wall 2, hvac 1, main 1.
  - At ce1f659 it was 200 sites in 152 routines (`docs/inventory/point_to_mesh_calls.csv`). My enclosing-routine scan on FireX gives ~163 routines (heuristic). The Mapper's `docs/inventory/base_delta.md` independently confirms 200 → 214 sites and 152 → 163 routines, and flags two new literal `CALL POINT_TO_MESH(1)` calls (`vtkf.f90:1875, 2482`) that use mesh 1 as a holder for global VTK slice metadata. The shim must preserve that, or regrid box renumbering will break it.
  - Some routines assume the call already happened, e.g. `! Assumes POINT_TO_MESH(NM) has been called.` (`ccib.f90:370`, `:1071`, `:1093`).
- **Direct global access.** `MESHES(NM)%` has **3,037** refs (2,926 at ce1f659). Of these, geom has 1,944 and ccib 542, i.e. 82%; then main 132, pres 93, dump 84, part 50. `MESHES(<idx>)%` has 3,560 refs, `MESHES(NOM)` 382, `=> MESHES(` 530.
- **Kernel shape.** Kernels take a mesh number, not arrays. Example: `DENSITY(T,DT,NM)` (`mass.f90:365`) calls `POINT_TO_MESH(NM)` (`:397`), sets external-wall ghosts in a loop (`:424`), then loops `DO K=1,KBAR / J / I` (`:443-445`) inside `!$OMP DO`. Whole-mesh loop counts: `DO I=1,IBAR` 217, `DO K=1,KBAR` 210.
- **Ghost indexing.** Arrays are local and 0-based.
  - Two ghost layers for `TMP, RHO, RHOS, ZZ, ZZS, WORK_PAD` (`init.f90:524-529`).
  - Face velocities get one extra layer on the low side of their own direction, e.g. `U(-1:IBP1,0:JBP1,0:KBP1)` (`init.f90:533-535`).
  - Most other fields have one layer, e.g. `FVX(0:IBP1,...)` (`init.f90:543-545`) and `H(0:IBP1,...)` (`init.f90:556`).
  - Token counts: `0:IBP1` 206, `IBP1` 438, `SIZE/LBOUND/UBOUND(` 375.
  - Metrics are 1-D per-mesh arrays (`DX(` 574, `RDX(` 129, `RDXN(` 90). Stretched grids are supported (`NAMELIST /TRNX/`, `read.f90:1000`); AMReX has uniform spacing per level.
- **Threading.** 356 `!$OMP` lines and 66 `!$OMP PARALLEL` regions, all inside kernels.
- **MPI.**
  - Several meshes may share a rank, but only as contiguous blocks (`LOWER/UPPER_MESH_INDEX`, `ERROR(117)`, `read.f90:711-722`).
  - 436 `CALL MPI_*` statements (main 214, ccib 94, pres 31, geom 28, vtkf 27, dump 14).
  - The halo exchange is hand-written: `POST_RECEIVES` (`main.f90:2945-3111`) and `MESH_EXCHANGE(CODE)` (`main.f90:3117-3975`, about 860 lines). Neighbour structure is fixed at setup (`INITIALIZE_MESH_EXCHANGE_1`, `main.f90:2059-2301`). Each mesh holds `OMESH(NMESHES)` (`main.f90:2077`). Counts: `NMESHES` 499, `OMESH` 730.
- **Two-stage ghost fill** (Mapper, spot-checked).
  - `MESH_EXCHANGE` first unpacks into OMESH copies, even on the same process. Kernels then copy or average OMESH into their ghosts: `ASSIGN_GHOST_VALUE` (`wall.f90:282`ff), `COPY_H_OMESH_TO_MESH` (`pres.f90:4054`), `NO_FLUX` (`velo.f90:1348`).
  - Kernel-side cross-mesh reads: `OMESH(NOM)%` 249, `EWC%NOM` 84, `DO IIO=EWC%IIO_MIN` loops 93 (ccib 43, velo 22, part 9, pres 8, geom 6).
- **Static coarse-fine interfaces already exist.** `EXTERNAL_WALL(IW)` holds neighbour index ranges `IIO/JJO/KKO_MIN:MAX`, and kernels average over them into the ghost cell.
  - The same loop copies when the neighbour is coarser and averages when it is finer: `velo.f90:514-547` for `MU, KRES, D/DS`; `wall.f90:321-330` with area weights `ARO`; `vege.f90:678`.
  - Exception: `turb.f90:2283` says `! assumes no refinement`.
  - AMReX equivalents (Integration Lead; verified): fine side `FillPatchTwoLevels` + `PCInterp` (`AMReX_Interpolater.H:420`), coarse side `average_down` + flux register (FR-016).
  - So kernels already tolerate a neighbour of different resolution **through ghost cells**, which favours reusing them on per-box data with ghosts filled by AMReX.
- **Global syncs** any driver must reproduce:
  - dt = `MINVAL(DT_NEW)` (`main.f90:715, 738-741, 885-891`);
  - zone integrals `DSUM/PSUM/USUM` via `MPI_ALLREDUCE` (`main.f90:2038-2040`);
  - HVAC network on rank 0 (`main.f90:829`).
- **I/O.**
  - Smokeview `.smv` is written once (`WRITE_SMOKEVIEW_FILE`, `dump.f90:1766-3019`), with one `GRID` block per mesh (`:2499`) and `OBST` blocks (`:2541`).
  - Files are keyed by mesh number (`CHID_<NM>_<N>.bf` `:524`; `CHID_<NM>.prt5` `:552`). Restart is per mesh (`DUMP_RESTART`, `:3871`).
  - **FireX adds VTK output.** `vtkf.f90` (module `VTK_FDS_INTERFACE`) writes VTK/VTKHDF `UnstructuredGrid` data (`vtkf.f90:758`) and a ParaView state file. It is selected with `&DUMP WRITE_FORMAT='SMV'|'VTK'|'BOTH'` (`read.f90:2382, 2429-2442`). A non-Smokeview output path is therefore already upstream-sanctioned.
- **FireX GPU work (item c).**
  - The GPU model is **library offload of the pressure linear solve only**. HYPRE is built with its own CUDA/HIP/SYCL backend (`CMakeLists.txt:16-18, 180-185`: `USE_HYPRE_NVIDIA/AMDGPU/INTELGPU` → `HYPRE_ENABLE_CUDA/HIP/SYCL`; makefile adds `-DWITH_HYPRE_DEVICE`, `Build/makefile:124-126`).
  - Fortran calls `HYPRE_SETMEMORYLOCATION(HYPRE_MEMORY_DEVICE)`/`HYPRE_SETEXECUTIONPOLICY(HYPRE_EXEC_DEVICE)` (`pres.f90:1177`ff). It migrates the IJ matrix and vectors host↔device around each solve (`HYPRE_IJVECTORMIGRATE`, `pres.f90:1753-1767, 3454-3469, 4818-4853`).
  - Ranks are grouped into "resource sets", one per GPU (`FDS_RANKS_PER_GPU` env var, `MPI_COMM_RS`, `main.f90:5105-5126`). Unknowns are gathered to the RS master, which alone solves (`pres.f90:3411-3418`).
  - There are **zero** `!$OMP TARGET`/`!$ACC` directives in `Source/`, so every physics kernel stays on the CPU.
  - Reading: FDS developers are moving toward GPU **through libraries**, not by porting Fortran kernels. That aligns with an AMReX driver whose MLMG (or HYPRE-in-AMReX) runs on the device. Under D-027 it is not enough: every physics kernel must also run on the device, which FireX gives no precedent for.

### AMReX F_Interfaces (`(local AMReX checkout)/Src/F_Interfaces`) — verified contents
- 64 files, 10,140 lines, in 5 directories: `Base`, `AmrCore`, `LinearSolvers`, `Particle`, `Octree`.
- **AmrCore:** init/regrid/callbacks (`AMReX_amrcore_mod.F90:11-23`); FillPatch single/two-level plus face variants (`AMReX_fillpatch_mod.F90:14-23`); `FluxRegister` and `FlashFluxRegister`; tagging.
- **MultiFab** accepts per-direction ghosts and nodality: `amrex_fi_new_multifab(mf,ba,dm,nc,ng,nodal)` with `ng(3), nodal(3)` (`Base/AMReX_multifab_mod.F90:177-182`). So FDS's per-field-group ghost widths and staggered faces are expressible.
- **Linear solvers:** only `amrex_poisson` and `amrex_abeclaplacian`, both cell-centred and multi-level. MLMG bottom solvers include HYPRE and PETSc (`AMReX_multigrid_mod.F90:9-14`). `set_acoeffs`/`set_bcoeffs` are exposed (`AMReX_abeclaplacian_mod.F90:15-16`), so a variable-coefficient ∇·(β∇H) projection can be driven from Fortran. There are no nodal, EB or tensor operators, and no overset-mask binding (`rg overset` finds nothing). The pressure spec's option E-2 relies on that mask (`docs/pressure/01-amr-mapping-spec.md` §E). There is no `MacProjector` either, which lives in AMReX-Hydro, C++.
- **Particles:**
  - One fixed `bind(C)` type `amrex_particle` with `pos(3), vel(3), id, cpu` (`AMReX_particlecontainer_mod.F90:18-23`), built as `AmrParticleContainer<NSTRUCTREAL=BL_SPACEDIM, NSTRUCTINT=0>` (`AMReX_particlecontainer_fi.cpp:8-13`).
  - Only add/get/count (per MFIter or grid), `redistribute` and `write` are available (`:28-38`). There are no runtime real/int components.
  - FDS `LAGRANGIAN_PARTICLE_TYPE` (`type.f90:389-429`) carries many attributes plus per-particle `BOUNDARY_ONE_D` surface storage. That data would need a side array re-synchronised on every redistribute, or would have to go through C++. (Teammate input, verified line-by-line against the clone: verdict partial, leaning poor.) The Integration Lead's `docs/amrex/driver-options.md` (DRAFT, row 6) agrees: "partial → effectively missing for FDS". Exposure: 160 of 941 verification cases use particles in the FireX survey (`docs/vv/verification_case_survey.csv`, `part=True`). driver-options.md quotes 157 of 928, which is the `ce1f659` survey.
- **Absent:** EB (no `EB2`/`EBFArray` symbols), nodal/MAC projection, EB redistribution, `AmrLevel`.
- **GPU — own source check.** Nothing in `F_Interfaces` is GPU-aware except the C++ internals of `FlashFluxRegister` (`Gpu::DeviceVector`, `AMReX_FlashFluxRegister.H:112-113`) and one OpenMP guard in the octree. `amrex_fi_multifab_dataptr*` hands Fortran a raw `Real*` into FAB memory (`Base/AMReX_multifab_fi.cpp:50, 63`) with no host/device staging, so on a GPU build Fortran would dereference device memory. The docs agree: "The Fortran interface of AMReX does not currently have GPU support. AMReX recommends porting Fortran code to C++ when coding for GPUs." (`Docs/sphinx_documentation/source/GPU.rst:105-106`). AMReX's own test suite builds the Fortran-interface tests only when `AMReX_GPU_BACKEND STREQUAL NONE`, with the comment "The Fortran interface tests do not work on GPU" (`Tests/CMakeLists.txt:168-169`). CMake does not block the library combination, so it is unsupported rather than forbidden. (This check was done independently. The Integration Lead's `driver-options.md` row 7, which landed afterwards, reaches the same verdict, "missing", with the same citations.) driver-options.md adds two Fortran-route limits that ADR-001 did not list. First, the refinement ratio is a scalar per level, and FillPatch and flux registers are isotropic (`AMReX_fillpatch_fi.cpp:232, 282`; `AMReX_fluxregister_fi.cpp:11`, all building `IntVect(rr,rr,rr)`). Second, there is no Fortran checkpoint-Header helper. Both favour Option A.
- **Build:** `AMReX_FORTRAN_INTERFACES` defaults OFF (`Tools/CMake/AMReXOptions.cmake:289-290`).
- **Maintenance.** The clone is shallow, so there is no `git log`. In `CHANGES.md`, the last *feature* work is 24.10 (average-down functions #4124, nvfortran fix #4115), then 23.10 (face FillPatch, #3541–#3553), 21.04 (#1793) and 18.08 (particles). The only 26.x item is a bugfix (#5604, 26.09, line 52). Reading: kept compiling, grows only on user demand (e.g. FLASH-X → `FlashFluxRegister`). Evidence is thin; confirm with AMReX developers.
- A working Fortran-driven subcycling example exists: `amrex-tutorials/ExampleCodes/FortranInterface/Advection_F` (1,608 lines).

### GPU offload routes in AMReX @ `99ddfda` (verified for v0.2)
- AMReX GPU backends are `NONE|SYCL|CUDA|HIP` only (`Tools/CMake/AMReXOptions.cmake:124-125`). OpenMP `target` and OpenACC are not backends, and AMReX's CMake has no offload option (`rg` over `Tools/CMake`). The GNU-make docs allow `USE_ACC=TRUE` for PGI, Cray and GNU (`GPU.rst:137`) and say OpenMP offload is supported only with IBM compilers (`GPU.rst:140`).
- Pragma kernels on AMReX memory are documented: a C++ `MFIter` loop passes `BL_TO_FORTRAN_BOX/ANYD` to a Fortran routine, which marks the FAB pointer `deviceptr` (OpenACC) or `is_device_ptr` (OpenMP `target`) (`GPU.rst:1457-1530`). The next section notes that CUDA/HIP launches are asynchronous (`GPU.rst:1536`ff), so pragma regions on the compiler's own queue need explicit ordering against AMReX's stream (R-39).
- AMReX CI builds a CUDA AMReX with the NVIDIA HPC SDK (`nvc`/`nvc++`/`nvfortran`, job `tests-nvhpc-nvcc`, `.github/workflows/cuda.yml:188-252`, Fortran compiler at `:242`) and a HIP AMReX with ROCm `flang` (`.github/workflows/hip.yml:21, 66`; AMD/HIP out of scope per owner decision 2026-09-25). Neither job compiles OpenACC/OpenMP-target code, so Fortran offload on AMReX memory is not tested upstream.
- Toolchains on the development machine (checked 2026-09-25): gfortran; Intel oneAPI 2026.1 at `/opt/intel/oneapi` with `ifx` 2026.1.1 and `icpx` (ifx offload targets Intel GPUs, not NVIDIA; it serves as K2's optional, non-gating second-compiler check, not tried); **no** `nvfortran`, `nvc++` or `nvcc` (NVIDIA HPC SDK install is A-31, open); no GPU (R-38).

### Ecosystem precedent
- IAMR, incflo and PeleLMeX contain 0 Fortran files. Their kernels are `amrex::ParallelFor` lambdas (`rg -c ParallelFor` over `Source/`: PeleLMeX 262, incflo 153; ERF 1,478).
- ERF (C++ AmrCore) calls legacy WRF Fortran microphysics through `BIND(C)` (`ERF_module_mp_morr_two_moment_isohelper.F90:29`):
  - a whole-FAB bridge with tiling disabled on the Fortran path (`ERF_AdvanceMorrison.cpp:204-208`);
  - pinned host memory on GPU builds (`:264-268`);
  - a **C++ port as default**, with the Fortran path kept as the reference answer (`use_morr_cpp_answer = true`, `:188`).

## Decision drivers (ranked) and how each option scores
The pressure path is **not** decisive on its own, so each driver is scored independently.

| # | Driver | A: C++ driver + box kernels | B: Fortran on F_Interfaces | Evidence |
|---|---|---|---|---|
| 1 | Preserve validated physics (minimal kernel rewrite) | same | same | The work is set by FDS global state, not driver language (see shim analysis below). Under D-027 every time-step kernel is rewritten for the device anyway, validated against its shimmed original |
| 2 | **Pressure path** | good (`MacProjector`, `MLABecLaplacian` + overset mask, EB) | **adequate** for cell-centred IBM (E-1) or variable-β; **not** for EB (E-3) or overset-mask masking (E-2) without new bindings | `AMReX_abeclaplacian_mod.F90:15-16`; pressure spec §E, §G.1 |
| 3 | **Particles** (FR-050..052) | good (`ParticleContainer` with runtime SoA comps) | partial→poor: fixed pos/vel struct plus an FDS side array synced on every redistribute | `AMReX_particlecontainer_mod.F90:18-38`, `type.f90:389-429` |
| 4 | **Geometry** (ADR-003) | EB available if chosen | EB unavailable without writing bindings | no EB symbols in F_Interfaces |
| 5 | **GPU: full time step on device** (D-027, NFR-043; NVIDIA) | required and reachable: kernels ported in style K1 or K2; FFT::Poisson/MLMG on device | none; AMReX recommends C++. **Disqualifying under D-027** | `GPU.rst:105-106`; `multifab_fi.cpp:50,63`; `Tests/CMakeLists.txt:168-169`; FireX `pres.f90:1177` |
| 6 | **Maintenance and access to new AMReX features** | full, same day | only what someone binds. Last feature work 24.10; option OFF by default | CHANGES.md; `AMReXOptions.cmake:289` |
| 7 | Team skills / single language (FDS developers Fortran-only, D-027) | worse (two languages); mitigated by the layering below and the kernel style (K1/K2) chosen by the P1 readability review (NFR-044) | better, but moot: B cannot meet driver 5 | owner statement; review pending |
| 8 | I/O (Smokeview/VTK/plotfile) | neutral: either option keeps `dump.f90`/`vtkf.f90` in Fortran | neutral | `read.f90:2429-2442` |

B wins only on driver 7. It is adequate on 1, 2 (within limits) and 8, and loses on 3, 4, 5 and 6. Under D-027, driver 5 alone rejects B. The pressure path alone favours C++ only if ADR-003 picks EB or masking via overset mask. That is why this ADR depends on ADR-003 only partially: drivers 3, 5 and 6 favour A regardless of ADR-003. Driver 7 is now addressed inside Option A, by layering and kernel style, not by the driver language.

## Options

### Option A — C++ AmrCore driver; FDS physics as box kernels (owner-confirmed, D-027)
C++ owns the `AmrCore` subclass, time loop, regrid, FillPatch/FillBoundary, flux registers, MLMG/`MacProjector`/`FFT::Poisson`, `ParticleContainer` and checkpointing. During the CPU migration, Fortran kernels are `BIND(C)` routines taking `(lo, hi, array views, dx, ...)`. The end state is device kernels in the style chosen under "Kernel implementation style" (K1 or K2).
- **Pros:** see drivers 2–6.
  - The ~860-line `MESH_EXCHANGE`, `POST_RECEIVES`, OMESH and most of the 436 MPI calls retire. The two-stage OMESH ghost fill collapses to one `FillBoundary`/`FillPatch` per MultiFab.
  - Direct ERF precedent.
- **Cons:**
  - Mixed-language build and debugging.
  - Kernels must take explicit arrays and `lo:hi` bounds (214 `POINT_TO_MESH` sites; 3,037 `MESHES(NM)%` refs, mostly in out-of-scope geom/ccib). Under D-027 this is required for every kernel, not eventual.
  - The driver layer needs real C++ skill (see "Layering for maintainability").

### Option B — FDS stays a Fortran program on AMReX F_Interfaces — rejected definitively (D-027)
The AMReX Fortran interface has no GPU support: "The Fortran interface of AMReX does not currently have GPU support. AMReX recommends porting Fortran code to C++ when coding for GPUs." (`GPU.rst:105-106`); its tests are built only without a GPU backend (`Tests/CMakeLists.txt:168-169`); `amrex_fi_multifab_dataptr*` hands out raw FAB pointers with no staging (`Base/AMReX_multifab_fi.cpp:50, 63`). A full time step on the device (NFR-043) is therefore impossible on this route. The analysis below is kept as the record.
- **Pros:**
  - One language.
  - AmrCore, FillPatch (incl. faces), FluxRegister and cell-centred multi-level MLMG with HYPRE bottom are all available. This is enough for a non-EB, particle-free, CPU-only prototype with E-1/variable-β pressure. Advection_F shows the pattern.
- **Cons:**
  - Particles partial→poor.
  - No EB, no overset mask, no `MacProjector`, no GPU.
  - Maintenance-mode layer, so every gap becomes our binding code.
  - Still needs the same kernel refactor as A.

### Option C — Hybrid: A with a transitional per-box `POINT_TO_BOX` shim (recommended migration mechanism)
Under the C++ driver, a Fortran shim binds FDS's module pointers to FAB memory for one box. It uses `C_F_POINTER` plus F2003 bounds remapping (`U(-1:,0:,0:) => view`, zero-copy). It synthesises per-box scalars and metrics (`IBAR..`, `X, XC, DX, RDX, RDXN`) and calls the existing, unmodified kernel. Kernels migrate one by one to explicit-argument device kernels in the chosen style (K1 or K2), each behind a legacy-vs-new comparison.

**Status under D-027: CPU-only migration stepping stone, not an end state.** Module pointers are host descriptors in process-global state (`REAL(EB), POINTER` at `mesh.f90:367-376`), so shimmed kernels cannot run on the device. Extracting every kernel onto the device is mandatory (D-027 (b); roadmap Phase 11). The shim's value is that each kernel is validated against FDS on AMReX data before it is rewritten.

**Scope decisions for the shim:**
1. **The pressure path is excluded entirely.** `PRESSURE_SOLVER_*` and `PRESSURE_ITERATION_SCHEME`'s solve are replaced wholesale by one composite MLMG/`MacProjector` call on MultiFabs, for these reasons:
   - The `pois.f90` FFT (`H3CZSS`, `pois.f90:187`) needs a whole rectangular mesh.
   - The iteration's mesh-interface repair has no purpose under a composite solve.
   - Exchange CODE 5 (`main.f90:1617, 1636, 1663, 1687`) and `COPY_H_OMESH_TO_MESH` (`pres.f90:4054`) drop out.
   - The solid and baroclinic reasons for iterating remain as the pressure spec defines them (§D; ADR-002/003).
   - The uniform level-0 solver `amrex::FFT::Poisson` (FR-037, D-021) is the only exception, and it lives behind the composite-solver interface (see "Pressure path and global reductions").
2. **Cost record (Integration Lead, verified where possible):**
   - Bounds remapping is zero-copy (meets IR-005).
   - Module pointers are process-global, so one box at a time per rank, no OpenMP tiling over boxes (FDS inner `!$OMP` stays; NFR-012, R-21), and no GPU until kernels are extracted. On a GPU build, shimmed kernels need their FABs in a pinned or managed arena (driver-options.md §5), so they are a transition cost, not a configuration to ship.
   - Per-box side data (`WALL`, `EXTERNAL_WALL`, `CELL_INDEX`/`CELL`, `X/RDX` metrics) must be rebuilt on every regrid.
   - Mesh-number-indexed state (`OMESH(NOM)`, `EXTERNAL_WALL%NOM`, `NMESHES`-sized arrays: 499 `NMESHES` refs) breaks when box ids change.
3. **Exit condition (R-26), timing only.** Whether kernels leave the shim is decided (D-027: all of them, onto the device). R-26 now decides only **when**: extraction is brought forward if **regrid rebuild exceeds 10% of step time**, or if **no kernel-extraction date is set by M4**; otherwise it completes in Phase 11. Extraction can start as soon as the kernel style is recorded (NFR-044).

**Feasibility and pitfalls, stated plainly:**

| Pitfall | Evidence | Assessment |
|---|---|---|
| Lower-bound remap | FDS box-local 0-based; FABs use global indices | Solvable, zero-copy; per-box 1-D metrics synthesised (uniform per level) |
| Face-index convention | FDS `U(I)` = *high* face of cell I, real faces `0:IBAR` (`init.f90:533`); AMReX x-nodal i = *low* face | +1 offset in the normal direction; the most likely off-by-one source; unit-test every staggered array |
| Ghost widths differ per field and per direction | 2 / 1 / asymmetric (`init.f90:524-545`) | Expressible per MultiFab (`ng(3)`, `nodal(3)`, `multifab_mod.F90:177-182`). Group fields by width. 375 `SIZE/LBOUND/UBOUND(` refs need an audit for extent-derived loop bounds (Mapper `module_globals.csv`, pending) |
| OMESH cross-mesh reads | 249 `OMESH(NOM)%`, 84 `EWC%NOM`, 93 `IIO_MIN` loops; kernels *write* ghosts from OMESH (`velo.f90:514-547`) | **Main edit surface.** With AMReX filling ghosts first, the `NOM>0` branches would overwrite them. Mark box-boundary wall cells "filled externally" (edits ~84 sites) or fake a per-box OMESH. Not zero-edit |
| Direct `MESHES(NM)%` access | 3,037 refs, 82% geom/ccib | Those routines cannot use the shim. That is acceptable only while GEOM/CC_IBM stays out of scope (FR-044, R-07). pres/main/dump/part (359 refs) need review, but pres leaves the shim anyway (decision 1) |
| Box-sized unstructured data | `N_EXTERNAL_WALL_CELLS = 2*IBAR*JBAR+...` (`read.f90:700`); `INIT_WALL_CELL` (`init.f90:2975-3386`) | Per-box wall store owned by the driver, rebuilt on regrid (ADR-003, FR-041); cost bounded by the R-26 exit condition |
| Call-order preconditions | `ccib.f90:370/1071/1093` | Enter the shim at the same call-tree level as `POINT_TO_MESH` today |
| Thread safety | global module pointers | One box at a time; tile = box; keep `max_grid_size` ≈ legacy mesh size |
| GPU | host Fortran, global pointers | **Incompatible by construction.** CPU-only migration stage; extraction mandatory (D-027) |
| Regrid invalidation | FAB memory moves | Re-bind (~400 pointer assignments) on every entry; never cache |

Net: feasible for the core gas-phase kernels (mass, velo, divg, turb, fire, radi, wall BCs) with honest edits at cross-mesh sites. Not zero-edit, and not viable for GEOM/CC_IBM.

### Option D — "each AMReX box is an FDS mesh" (dynamic MESHES) — rejected
The static topology is wired in everywhere (`NMESHES` 499 refs, `OMESH(NMESHES)`, neighbours fixed in `INITIALIZE_MESH_EXCHANGE_1`, contiguous rank blocks). Per-mesh init is heavy (`init.f90` 5,531 lines). There is no AMR time interpolation or refluxing.

### Option E — full rewrite up front — rejected for phase 1
~181k lines, and all V&V would have to be re-earned at once. Under D-027 every time-step kernel is rewritten anyway (in K1 or K2), but incrementally, kernel by kernel, behind legacy comparisons; input parsing, setup and output are not rewritten (see layering).

## Kernel implementation style (D-027, NFR-044)
Both candidates need the same rewrite of each kernel: explicit array arguments and `lo:hi` bounds instead of `POINT_TO_MESH` module pointers, no `MESHES(NM)%`/`OMESH` access, no derived-type records inside the loop. They differ in the language of the loop body and the toolchain. Both target NVIDIA, the only GPU target; AMD is out of scope per owner decision 2026-09-25.

### K1 — restricted "Fortran-style C++" `ParallelFor` bodies
- Each kernel is an `amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {...})` over `Array4` views: plain `(i,j,k)` loops, the same index order and global lower bounds FDS developers already read. Physics code uses no templates, classes, inheritance or operator overloading; a written style guide fixes the allowed subset and naming (FDS variable names kept).
- Builds for CUDA through AMReX's supported backend. The same source is portable in principle (AMReX CI also builds HIP, `.github/workflows/hip.yml`), but portability is no longer a criterion (AMD out of scope per owner decision 2026-09-25). PeleLMeX, ERF and incflo use this pattern throughout (counts under "Ecosystem precedent").
- Toolchain: AMReX's supported CUDA path (nvcc or NVIDIA HPC SDK `nvc++`). The device-code rules are enforced by the compiler (`AMReX_CUDA_ERROR_CROSS_EXECUTION_SPACE_CALL`, `AMReX_CUDA_ERROR_CAPTURE_THIS`, as in `cuda.yml:244-245`).
- Cost: FDS's Fortran developers read and review C++ syntax (lambdas, `amrex::Real`, 0-based component index); the kernel no longer matches the FDS Fortran source line by line, so each port is checked against the shimmed Fortran kernel (T1 on frozen input).

### K2 — Fortran kernels with simple OpenMP `target` offload on AMReX device memory (v0.3, owner decision)
- Each kernel stays Fortran. The C++ `MFIter` loop passes box bounds and FAB device addresses through `BIND(C)` (`GPU.rst:1457-1530` documents the pattern); the loop nest carries one `!$omp target teams loop` directive.
- **Allowed subset:** `!$omp target teams loop` plus `collapse`; data via `has_device_addr` (see "Kernel interface and device data"); `map` only for small host scalars/constants. Not allowed: `distribute parallel do`, nested or combined constructs beyond this one, `declare target` on module data, or other complex constructs. The style guide shared with K1 fixes naming (FDS variable names kept).
- **Toolchain:** `nvfortran` for the kernels, with AMReX built for CUDA by the same NVIDIA HPC SDK (`nvc++`, as in `cuda.yml:188-252`). The SDK is free and compiles without a GPU (A-31).
- **Second-compiler check (optional, non-gating; v0.3.7):** kernel code may also be compiled with Intel `ifx` 2026.1.1, which is already on the development machine (not tried). It is a source-level check only and gates nothing. `amdflang` is dropped (AMD out of scope per owner decision 2026-09-25).
- **CUDA Fortran (optional, NVIDIA-only):** allowed only for hot loops where profiling on real NVIDIA hardware shows a meaningful gap against the OpenMP version. Each such kernel keeps its OpenMP version as fallback and correctness reference (T1 against it). None exist yet; none can be justified before GPU hardware (R-38, Q11 (b)).
- Caveats:
  - Offload is not an AMReX backend (`AMReXOptions.cmake:124-125`); AMReX is still built for CUDA and the offload compiler must share its device runtime (R-39). Upstream CI builds AMReX with `nvfortran` but compiles no offload code, so the pairing is untested upstream. AMReX's GNU-make docs list OpenMP offload as supported only with IBM compilers (`GPU.rst:140`), so the docs give no nvfortran precedent.
  - AMReX kernels run asynchronously on AMReX's stream, `target` regions on the OpenMP runtime's queue; every C++/Fortran kernel boundary needs a synchronization or the OpenMP queue wired onto AMReX's stream (R-39). The cost per boundary is unmeasured.
  - `POINT_TO_MESH` module pointers (`mesh.f90:367-376`) are host descriptors and cannot enter `target` regions as they are, so kernels need explicit array arguments: the same rewrite as K1.
  - **Portability:** not a criterion. Only the optional CUDA Fortran kernels are NVIDIA-specific, and NVIDIA is the only GPU target; AMD/HIP pairing questions are out of scope per owner decision 2026-09-25.
- Status here: `nvfortran` is not installed (A-31 open), so K2 cannot be compiled for the target yet.

| | K1 restricted C++ | K2 Fortran + OpenMP `target teams loop` |
|---|---|---|
| Language FDS developers read | restricted C++ | Fortran, same OpenMP family FDS already uses |
| NVIDIA toolchain | AMReX CUDA (nvcc or nvc++) | nvfortran + nvc++ (NVIDIA HPC SDK), one toolchain |
| Second compiler | not needed | optional, non-gating `ifx` compile (kernel source only) |
| NVIDIA-only code | none | optional CUDA Fortran hot loops only (OpenMP version kept) |
| AMD | out of scope per owner decision 2026-09-25 (source portable in principle) | out of scope per owner decision 2026-09-25 |
| AMReX upstream support | backend, CI-built, used by PeleLMeX/ERF/incflo | documented pattern; not a backend; not CI-tested |
| Stream/sync | native | sync or stream wiring at every boundary (R-39) |
| Device data | `Array4` captured by value | `has_device_addr` on explicit-shape dummies (fallback `is_device_ptr` + `c_f_pointer`); no `map` of field data |
| Kernel rewrite (explicit args) | required | required |
| Compilable on the development machine now | CUDA: no (A-31); a SYCL device compile with the installed `icpx` is possible in principle (not tried) | NVIDIA: no (A-31); optional ifx compile possible in principle (not tried) |

### Kernel interface and device data for K2 (IR-007)
- **No-copy rule.** AMReX-owned arrays already live in device arena memory, so `map(tofrom:)` on them at kernel entry would copy again or be wrong. Kernel entry wrappers pass AMReX `Array4` data as explicit-shape dummy arrays and declare them `has_device_addr` (OpenMP 5.1) on the `target teams loop` construct. If the first nvfortran compile (A-31) shows `has_device_addr` unsupported for Fortran dummy arrays, fall back to `is_device_ptr` on `type(c_ptr)` arguments with `c_f_pointer` inside the target region. The mechanism is confirmed by that first compile.
- `map` clauses are allowed only for small host scalars/constants, never for AMReX-owned field data.
- **AMReX managed memory** (`amrex.the_arena_is_managed=1`) is rejected for production (page-migration cost; it hides missing-device-data bugs) and allowed only as a debugging aid.
- **Threading rule (reworded for v0.3):** no host-threading OpenMP (`parallel do`) inside device kernels; `target` directives appear only in kernel files. AMReX handles host threading and tiling from the driver. This replaces "no `!$omp` inside device kernels" (IR-007), which would forbid K2's own directives.
- The Spec & Program Lead owns the IR-007/NFR-044 text (currently "OpenACC Fortran", `deviceptr`, "no `!$omp` threading inside the kernel") and has been asked to update it; IR-007 carries the no-copy rule.

### Species/density clipping (D-031)
Applies to both kernel styles and to the CPU shim path.

**Problem** (P1 finding, Integration Lead, `docs/amrex/p1-findings.md` §6.2; line numbers from the Legacy Mapper, FireX `mass.f90`). `CHECK_MASS_DENSITY` gives box-layout-dependent results whenever clipping triggers:
- Species: neighbour masses gated by `WALL_INDEX` (`mass.f90:907-912`), `CONST` (`:915`), scatter into `DELTA_RHO_ZZ` (`:916-922`), interior-only apply (`:931-937`, skipped per species unless `CLIP_RHO_ZZ(N)`, `:927`), early return (`:943`), renormalisation (`:947-961`).
- Density: the same pattern, neighbour masses at `:831-836`, scatter at `:840-846`, apply at `:853-854`; `CLIP_RHOMIN`/`CLIP_RHOMAX` set at `:815`/`:819`.
- Cause: every mesh-boundary face in FDS is a wall cell, including `INTERPOLATED_BOUNDARY` (`init.f90:76-107, 3253, 3298`), so `WALL_INDEX≠0` and `MASS_N=0` across the interface (`mass.f90:907-912`). A clipped cell next to a box edge spreads over fewer neighbours and gets a different `CONST` (`:915`). Nothing is lost into ghost cells; the redistribution stays within the box (clipping itself is not exactly conservative, p1-findings §6.2).
- Second cause: the renormalisation runs only if something clipped in that box (per-box flags, `:943`), so an ulp-level change reaches every cell of a box that clipped anywhere.
- Multi-mesh FDS behaves the same way (the code is per mesh).

**Decision: gather rewrite; no exemption or tolerance.** The byte-identical rule for explicit stages (T0 kernel parity, D-022; IR-007) stands.
- Each cell gathers the contributions it would receive, in FDS's K,J,I scatter order: k-1, j-1, i-1, self, i+1, j+1, k+1. This reproduces the floating-point summation order of single-mesh FDS bitwise.
- Ghost data: one pre-clip `FillBoundary` only; see "Ghost depth" below.
- Pass order per stage and level (v0.3.3; p1-findings §13.3), with two host-side `ParallelAllReduce::Or` reductions:
  - (a) density clip + gather over valid+1 (see "Ghost depth"); `CLIP_RHOMIN`/`CLIP_RHOMAX` set from valid, uncovered cells only (`:799-849`; see "Clip flags");
  - (b) **OR #1** on the host over the density flags `CLIP_RHOMIN`, `CLIP_RHOMAX`;
  - (c) density apply gated by the reduced density flags (`:853-854`);
  - (d) species gathers for all species (independent of each other: species N reads only `RHO_ZZ(:,N)` and the clipped `RHOP`), per-species flags `CLIP_RHO_ZZ(N)` set from valid, uncovered cells (`:870-925`);
  - (e) **OR #2** on the host over the per-species flag vector `CLIP_RHO_ZZ(1:NS)` (one packed reduction);
  - (f) per-species apply (`:927`, `:931-937`) gated by the reduced species flags;
  - (g) renormalisation (`:947-961`, replacing the per-box return at `:943`) gated by the reduced density flags OR'd with every reduced species flag (as `:943`; P1 `p1_driver.cpp:409`, `rmin || rmax || anyz`).
  - **One species** (`N_TRACKED_SPECIES==1`): after (c), FDS copies `RHO_ZZ=RHOP` and returns (FireX `mass.f90:858-860`; P1 `clip_gather.f90:163`), gated by the reduced density flag. Steps (d)-(g) do not run, so there is one host reduction, not two.
  - **Species-stage `RHOP` ceiling in ghosts** (`RHO_ZZ_MAX`, `:887`): it comes either from the redundant valid+1 density apply (the target) or from a second 2-layer `FillBoundary` of `RHOP` after (c) (the `fill2` cross-check).
  - **Host syncs:** production packs `CLIP_RHOMIN`/`CLIP_RHOMAX` into one reduction and all species flags into one, so there are at most 2 host syncs per stage and level whatever the species count. The P1 prototype issues them separately (2 + NS `ParallelAllReduce::Or` calls, `p1_driver.cpp:368-369, 404`).
  - Gating by the reduced flags is needed for bitwise parity, not tidiness: an unconditional apply turns a −0.0 ρZ into +0.0, and renormalising when single-mesh FDS would not changes ulps (§13.3-13.4).
- A gather writes only its own cell, so there are no atomics on the GPU.

**Ghost depth (v0.3.2 ruling).** The density clip+apply runs redundantly over the grown tile (valid+1), so the species gather sees clipped `RHOP` in 1 ghost layer without a second `FillBoundary`.
- Stencil (checked against FireX `mass.f90`): the species clip at a cell reads `RHOP` only at that cell (`RHO_ZZ_MAX = RHOP(I,J,K)`, `:887`; the neighbour masses at `:907-912` clamp with that same cell's `RHO_ZZ_MAX`) and `RHO_ZZ` at the cell and its 6 neighbours. The density clip is analogous: `RHOP` at the cell and its 6 neighbours (`:809-836`). Both read `SOLID` at the clipping cell only (`:811`, `:885`); `VC` uses `DX/DY/DZ` of the neighbours (`:801-807`, `:822-828`).
- Hence: the species gather at cell i uses clipped `RHOP` at i-1..i+1 and `RHO_ZZ` out to i±2; the density clip over valid+1 gathers from clips at valid+2, which read pre-clip `RHOP` out to valid+3.
- **Single pre-clip `FillBoundary`: ng=3 for `RHOP` (`RHOS` in the predictor, `RHO` in the corrector, `:789-795`), ng=2 for `RHO_ZZ`.** `RHO_ZZ` is not changed by the density stage (the one-species branch at `:858-861` returns first).
- The face mask and `SOLID` must be defined on valid+2; cell metrics on valid+3 (trivial under uniform spacing per level).
- Clip flags are set from valid cells only (confirmed by the A-37 ghost-poisoning test), and on coarse levels only from uncovered ones (see "Clip flags"); redundant ghost-cell clips never set a flag (otherwise a coarse-fine or physical-boundary ghost could trigger a renormalisation single-mesh FDS would not run).
- Bitwise identical to the alternative (a second `FillBoundary` between density apply and species clip): the ghost computation uses the same operations, in the same order, with the same mask as the owning box's valid computation (confirmed by A-37 in every case tested). The redundant form saves one exchange plus a device sync per call; on the CPU development machine the measured difference is small (§13.5), so the saving matters for many ranks and on the GPU.
- **Native layout:** the ng=3 `RHOP` temporary exists only because the shim's FDS-bounds arrays (`RHO/RHOS/ZZ/ZZS`, ng=2) cannot hold a third ghost layer (§13.3). In the native layout the density MultiFab is allocated with ng=3 and the species MultiFab with ng=2.
- Coarse-fine and physical-boundary ghosts are irrelevant: the mask is nonzero on those faces, so no contribution crosses them in either direction.
- Rejected alternative: a second `FillBoundary` after the density apply. A-37 showed it bitwise equal to the redundant form; it stays in the prototype as a cross-check (`p1.clip_ghost=fill2`), not as the target form.

**Face mask.** An explicit kernel argument, separate from `WALL_INDEX`:
- Defined on valid+2 (see "Ghost depth"). Built from the `WALL_INDEX` single-mesh FDS would compute: domain-boundary faces, including OPEN vents and periodic faces, are nonzero; exposed OBST faces only, both sides of thin obstructions.
- Only faces between boxes on the same level are zeroed.
- Faces at coarse-fine boundaries stay nonzero: no clip redistribution across levels, because fine ghost cells are interpolated and a cross-level transfer would be neither conservative nor refluxed.
- **Coarse side too** (Legacy Mapper, v0.3.3; see also "Clip flags" for covered cells): on the coarse level, the face between an uncovered cell and a covered cell (one lying under the fine level) is nonzero. The coarse cell neither gathers from nor pushes into the covered cell, because average-down overwrites covered cells and any clip mass moved there would disappear. This matches FDS, where both sides of a mesh interface have wall cells (`init.f90:76-107`).
- Layout independence therefore holds for box splits within a level, not across changes of level structure.

**Clip flags (v0.3.4 ruling).** Clip flags are set only from valid cells that are not covered by a finer level. This uses the per-level fine-covered mask, in the same way ghost cells are excluded.
- Rationale: average-down overwrites covered cells, so a clip there must not switch on renormalisation for the uncovered cells of the level.
- Clipping covered cells is allowed (their values are discarded) but may be skipped as an optimisation. The coarse-side mask already keeps uncovered cells from interacting with covered ones.
- Flags stay OR-reduced per level, which matches FDS's per-mesh behaviour.

**Kernel-interface rule (both styles):** domain-wide reductions happen on the host between kernel passes, never inside a kernel.

**Parity note.** Bitwise tests must replay FDS's exact time-step sequence (the final step is nudged to hit `T_END`, FireX `main.f90:719`). For clipped cases the reference is single-mesh FDS, not multi-mesh FDS.

**Result: A-37 closed** (Integration Lead, `docs/amrex/p1-findings.md` §13; `prototypes/p1_mass_shim/check_clip.sh`, 47 checks pass):
- The gather clip is bitwise across `max_grid_size` 32 vs 8 and 1 vs 2 ranks, and bitwise vs single-mesh FDS with 1 box and with 64 boxes on 2 ranks (where FDS's own clip differs by 5.2e-5 in RHO).
- The density clip was exercised with forced limits (`p1.rhomax=1.85`, `p1.rhomin=1.333`); bitwise.
- Both ghost variants (`fill2`, `redundant`) are bitwise equal.
- A ghost-poisoning test (`p1.clip_poison=1`) confirms that flags come from valid cells only.
- Two-level runs (static patch) are bitwise across box splits within each level.

**Cost and target form.** The unoptimised gather adds 25-40 % to loop time, because each neighbour's clip is recomputed 7 times per target (per pass) through a non-inlined call (§13.5). The target form is the §13.5 two-phase version: phase 1 computes each cell's clip terms once (a per-cell "clipped" byte plus the self term `CONST*SUM_MASS_N/VC(0)` and the 6 neighbour terms `CONST*MASS_N(d)/VC(d)`, same expressions and order as FDS); phase 2 gathers them in the K,J,I order above, only where a neighbour clipped. The terms are computed over valid+1 for the species gather; under the redundant ghost form the density terms are computed on valid+2, because the density gather covers valid+1 (ruling, v0.3.4). The input ghost depth still suffices: density terms on valid+2 read `RHOP` out to valid+3, which ng=3 already gives. §13.5 calls this form bitwise-safe but it was not built in A-37. It is accepted provided it passes the same 47 checks plus one new check: a two-level case where only a covered coarse cell clips must not renormalise the uncovered cells of that level ("Clip flags"). The same form suits the GPU better than recomputation (§13.6).

### Leaning
None. **K1 and K2 are co-equal candidates; the P1 readability review decides** (NFR-044, D-027 (c), D-029). The trade-off to weigh is readability for Fortran-only developers and continuity with FDS's existing OpenMP (favour K2) against upstream support and stream simplicity (favour K1). Portability is not a criterion: NVIDIA is the only GPU target and AMD is out of scope per owner decision 2026-09-25.

### P1 readability review record (NFR-044)
The mass kernel (`MASS_FINITE_DIFFERENCES`, `mass.f90:20`, and `DENSITY`, `mass.f90:365`) is written in K1 and K2 and reviewed by FDS Fortran developers. The K2 variant is OpenMP-offload Fortran (`nvfortran`, `target teams loop`); an OpenACC variant is not built (v0.3). Per D-029 the review happens only after both variants compile for a real GPU backend (nvfortran/nvc++ + CUDA after A-31); the K1 variant goes first. The K2 variant may also get the optional, non-gating `ifx` compile. Both variants also run on the host (K1 CPU build, K2 host fallback) against the shimmed Fortran kernel.

| Item | Value |
|---|---|
| Status | **pending** (blocked on A-31) |
| Reviewers | TBD (FDS Fortran developers) |
| Variants and commits | TBD |
| Device compile (backend, compiler versions) | TBD (K2: nvfortran OpenMP offload; `has_device_addr` vs fallback recorded here) |
| K2 optional `ifx` compile (non-gating) | TBD |
| Host result vs shimmed kernel | TBD |
| Outcome (K1 / K2) and reasons | TBD |
| Owner sign-off (if K2: any CUDA Fortran kernels) | TBD |

## Layering for maintainability
- **(a) Driver, regrid, pressure: C++.** `AmrCore` subclass, time loop, FillPatch/flux registers, regrid and side-data rebuild orchestration, `FFT::Poisson`/MLMG, particles container, checkpoint. Small and stable once written; needs real C++ skill, owned by the AMReX-side team.
- **(b) Physics kernels: restricted style (K1 or K2).** Where FDS developers work day to day. One kernel = one `ParallelFor` body or one Fortran loop nest under `!$omp target teams loop`, with explicit arguments.
- **(c) Input parsing, setup and output: may stay Fortran on the host.** `read.f90`, setup in `init.f90`, `dump.f90`/`vtkf.f90`. This matches NFR-043's I/O exception and keeps the largest FDS-specific code unchanged. Open: whether regrid-time side-data rebuild (`WALL`, `CELL_INDEX`, `EXTERNAL_WALL`) counts as part of the time step under NFR-043 (Spec & Program Lead).

## Pressure path and global reductions (room decisions, 2026-09-25; details in ADR-002 v0.2)
- **Solver choice, per-step selection, shared gauge and eps_H agreement:** see ADR-002 v0.2, "Accepted decisions / Pressure solver" (`amrex::FFT::Poisson` replaces porting `pois.f90`, D-021; FR-037, FR-039, D-012).
- **MLMG order:** `setMaxOrder(2)`, confirmed by P2; the FR-039 fallback to order 3 is not needed. Evidence and the other P2 rulings (mean removal, `average_down_faces`, HYPRE as bottom solver only): ADR-002 v0.2.
- **Refinement ratio ≤ 4, supported {2,4}:** see ADR-002 v0.2, "Refinement ratios" (FR-010, D-030).
- **Global scalar reductions:** exact fixed-point sums from per-box sums, computed domain-wide; see ADR-002 v0.2, "Decomposition requirements" (FR-005 (ii), (v); D-028; R-36). On GPU builds the per-box accumulation runs on the device.
- **On the GPU (NFR-043):** `FFT::Poisson` via cuFFT for single-level uniform runs; composite MLMG on the device otherwise. `PoissonHybrid` (z-stretched only; device branch at `AMReX_FFT_Poisson.H:715-780`) is no longer needed: stretched-grid cases stay FDS-only (owner decision 2026-09-25, D-030), so R-29/A-29 are moot. HYPRE is CPU-only in our builds, including the FireX-pinned v2.32.0-24 `63331f19c` (`HYPRE_USING_CUDA`/`HYPRE_USING_GPU` undefined, `(local GNU third-party library tree)/libs/hypre/63331f19/include/HYPRE_config.h:93, 144`; D-026), so MLMG's native bottom solver is the device default (a GPU HYPRE would need its own build per backend; R-30). Stretched meshes need no GPU pressure path: AMR mode uses uniform grids on each level and stretched cases stay FDS-only, which closes charter Q11 (c) (ADR-002 v0.2.1).

## Recommendation
**Adopt Option A (C++ AmrCore driver), owner-confirmed via D-027. Migrate through the Option C `POINT_TO_BOX` shim as a CPU-only stage with the three scope decisions above, then extract every kernel onto the device in style K1 or K2, chosen by the P1 readability review. Reject B definitively.**

The driver decision rests on four points, the first now decisive on its own:
1. GPU: the full time step must run on the device (D-027, NFR-043), and F_Interfaces has no GPU support (`GPU.rst:105-106`).
2. Particles: the fixed `amrex_particle` struct cannot carry FDS particle state.
3. Maintenance: the F_Interfaces layer is maintained but not developed, so every new AMReX feature would need our bindings.
4. ADR-003's EB option needs C++.

The pressure path is not among them: F_Interfaces can drive a cell-centred variable-β projection, though `FFT::Poisson` would have needed a C++ wrapper (FR-030).

Confidence: **high** for rejecting B (owner requirement plus upstream documentation); **medium** for the shim's edit count and regrid cost, pending S1/S3; kernel style **open** until the P1 review. What would change it: nothing short of withdrawing D-027 reopens B; the K1/K2 choice turns on the review, on the first nvfortran compile (device-data mechanism, R-39).

## Rejected alternatives and why
- **B:** rejected definitively (D-027): no GPU support in F_Interfaces. It also loses on particles, maintenance and EB, and wins only on single language. The kernel refactor is unchanged.
- **D:** static mesh topology; no AMR time interpolation or reflux.
- **E:** V&V and schedule risk; D-027's rewrite happens incrementally instead.
- **Shim as end state:** incompatible with the device by construction (D-027 (b)).
- **OpenACC for K2** (v0.2 candidate; replaced by owner decision, 2026-09-25): the most mature Fortran offload model on nvfortran, and for simple loops its performance is expected to be similar to OpenMP `loop`. Rejected because it is NVIDIA-only in practice, adds a second directive dialect alongside FDS's existing OpenMP, and is less portable. No performance numbers have been measured for FDS kernels in either model.
- **Full OpenMP offload feature set** (`distribute parallel do`, complex/nested constructs): harder to read for FDS developers and more compiler variance; K2 is limited to `target teams loop` plus `collapse`.
- **Second `FillBoundary` in D-031 clipping** (between density apply and species clip): one extra exchange plus a device sync per call for a bitwise-identical result (A-37); replaced by the redundant density clip over valid+1; kept only as a prototype cross-check.
- **Single OR over all clip flags in D-031:** wrong order; the density apply must be gated by the reduced density flags before the species gathers read the clipped `RHOP`, so two reductions are needed.
- **`map` of AMReX field data / AMReX managed memory in production:** `map` copies data already in the device arena (breaks IR-005's no-per-step-copy rule); managed memory costs page migration and hides missing-device-data bugs. Managed memory stays a debugging aid only.

## Consequences and risks
- **+** AMReX features without binding debt; ~860 lines of exchange and most MPI calls retire; one route to the device for driver, pressure and kernels.
- **−** Two languages; transitional shim code; every time-step kernel is rewritten (Phase 11 / M11; schedule TBD until S4).
- **R-26 (shim becomes permanent / too costly):** now a timing rule only; extraction itself is required.
- **R-38 (no GPU on the development machine):** NFR-043 is verified here by a real GPU-backend compile and host-fallback runs only (D-029); device correctness and the no-transfer rule wait for hardware.
- **R-39 (offload toolchain shares AMReX's device runtime, unproven):** applies to K2 only (nvfortran OpenMP offload, and any CUDA Fortran kernels); tested in the P1 review once nvfortran is installed (A-31), together with the `has_device_addr` mechanism.
- **R-31 (rank model):** FireX's `FDS_RANKS_PER_GPU` gather-to-master does not carry over to MLMG on device.
- **K2 NVIDIA-specific code:** CUDA Fortran kernels, if any, keep their OpenMP fallback. Portability is not a criterion (AMD out of scope per owner decision 2026-09-25).
- **C++ skill:** layer (a) needs C++ developers who stay with the project; layers (b)/(c) are sized for Fortran developers.
- **R-02:** this ADR chooses "adapter first, progressive refactor".
- **R-21/NFR-012:** FDS inner OpenMP on shimmed kernels.
- **Unstructured per-box state under regrid** (WALL/CFACE/1-D conduction/particles) — *highest technical risk*. Phase 2–4 use static or wall-avoiding regrids until ADR-003 lands.
- **R-23:** stretched grids (`read.f90:1000`) are rejected in AMR mode (IR-002); AMR mode uses only uniform grids on each level and the 59 stretched cases stay FDS-only (owner decision 2026-09-25, D-030).
- **R-08/R-18 (I/O and restart):** the FireX VTK path (`vtkf.f90`) and AMReX plotfiles are candidates for refined output; Smokeview needs static `GRID`s. **Q5 is an owner decision.** Output stays on the host (NFR-043 exception).
- **Teammate line numbers drift:** inventory and pressure docs cite ce1f659. Re-anchor them on FireX before M1.

## Open questions for teammates
- **FDS Legacy Mapper:**
  - Re-base the inventory on FireX `36975d765f`.
  - Classify the ~163 `POINT_TO_MESH` callers as shim-able (module pointers only), OMESH-reading, or `MESHES(NM)%`-direct.
  - `module_globals.csv` has landed: 1,551 module-level variables, 326 flagged `blocks_pure_kernel=yes` and 3 `maybe` (top modules: GLOBAL_CONSTANTS 63, FDS 45, CC_SCALARS 42, OUTPUT_CLOCKS 27, GLOBMAT_SOLVER 24). Next: tag which of the 326 are written inside the Phase-1 kernel set, so the shim's per-box save/restore list can be sized (R-26), and which are read inside kernels (they become kernel arguments or device constants under K1/K2).
- **AMReX Integration Lead:**
  - Install the NVIDIA HPC SDK (A-31) so both P1 variants compile for CUDA.
  - For K2: document how an `nvfortran` OpenMP `target` region is ordered against AMReX's CUDA stream, and measure the per-boundary cost (R-39). In the first compile (A-31), confirm `has_device_addr` on explicit-shape Fortran dummies, or switch to the `is_device_ptr` + `c_f_pointer` fallback.
  - Confirm pointer re-binding rules after `FillBoundary`/regrid.
  - D-031: build the §13.5 two-phase clip (terms once, then gather; density terms on valid+2), rerun the 47 `check_clip.sh` checks, and add the covered-coarse-cell check (only a covered coarse cell clips; uncovered cells not renormalised). A-37 is closed for the unoptimised form.
  - Confirm that the per-box fixed-point accumulation (FR-005 (ii)) can run in device code.
  - F_Interfaces maintenance status is moot for the driver (D-027).
- **AMR Spec & Program Lead:**
  - Update charter Q11, NFR-043, D-029, R-39 and A-30 for the NVIDIA-only target; (v0.3.7) AMD is out of scope, not deferred (A-32 moot); Q11 (c) closed (stretched cases FDS-only); NFR-045/A-36 become optional and non-gating, `ifx` only; R-29/A-29 moot. Q11 (b) stays open.
  - (Asked, v0.3) Update IR-007, NFR-044, R-38 and R-39 from OpenACC to OpenMP-offload Fortran: reworded threading rule (no host-threading `parallel do` inside device kernels; `target` directives only in kernel files), the no-copy device-data rule (`has_device_addr`; `map` only for host scalars), the second-compiler portability check, and CUDA Fortran only for profiled hot loops with OpenMP fallback.
  - Does regrid-time side-data rebuild count as part of the time step under NFR-043?
  - Minimum feature set for the Phase 2 demo; Smokeview as a hard requirement (Q5)? Stretched grids: *answered* (owner, 2026-09-25): not in AMR mode; stretched cases stay FDS-only.
- **AMR Pressure Solver Lead:** confirm the pressure path sits wholly outside the shim. *Done* (P2, ADR-002 v0.2): maxorder 2 vs 3 study and FFT/MLMG eps_H check (FR-039).
- **AMR V&V Lead:** S1 tolerance class T1 vs FireX baseline, with `OMP_NUM_THREADS=1`; T1 check of each K1/K2 port against its shimmed kernel.
- **GNU/Intel build chiefs:**
  - Mixed C++/Fortran/MPI link on gfortran 14.2 + OpenMPI 5 and oneAPI (NFR-020); add an NVIDIA HPC SDK + CUDA configuration (compile-only here).
  - Optional, non-gating `ifx` compile of K2 kernel files (amdflang dropped).
  - `C_F_POINTER` + bounds-remap idioms on all compilers.
  - AMReX OpenMP vs FDS `-fopenmp`.
  - HYPRE per D-026: FireX pins v2.32.0-24 `63331f19c` (FireX CMake mislabels it "3.0.0"), and AMReX is built against it at `(local AMReX install built against HYPRE 2.32)` (R-30). HYPRE is CPU-only and used only as MLMG's bottom solver; a CUDA HYPRE build only if HYPRE is ever wanted on device.

## Spike plan
| Spike | Scope | Pass | Overturns if |
|---|---|---|---|
| **S1 Shim feasibility** (2 wk) | Minimal C++ AmrCore (from `Advection_AmrCore`), single level, periodic, no walls. Shim binds `RHO, ZZ, U, V, W, TMP`; runs `DENSITY` and `COMPUTE_VISCOSITY` built from the FireX sources (outside the read-only tree). | T1 vs FireX baseline; edit count per kernel recorded; overhead ≤ 10% | pervasive kernel edits needed, or overhead > 30% (then extract kernels directly, skipping the shim) |
| ~~S2 F_Interfaces counter-spike~~ | Cancelled by D-027: its reopen condition required GPU out of scope. | — | — |
| **S3 Regrid wall-state rebuild** (2 wk, with ADR-003) | Rebuild `WALL`/`BOUNDARY_ONE_D` for a box with one OBST after a synthetic regrid | conserved exactly; rebuild < 10% of step time at a 10-step interval (R-26) | full re-init needed per regrid |
| **S4 GPU cost probe = P1 readability variants** (1–2 wk) | Port the S1 mass kernel to K1 (first) and K2 (OpenMP-offload Fortran, `target teams loop`; no OpenACC variant); compile both for CUDA (D-029, after A-31); K2 optionally also with `ifx` (non-gating); run K1 CPU build and K2 host fallback | both match the shimmed kernel at T1; effort per variant recorded (Phase 11 estimate); sync points per step counted for K2; device-data mechanism (`has_device_addr` or fallback) recorded | K2 cannot share AMReX's device memory (neither `has_device_addr` nor `is_device_ptr`) or stream with nvfortran (drops K2; R-39) |
| **S5 P1 readability review** (NFR-044) | FDS Fortran developers review both variants (K2 = OpenMP-offload Fortran; blocked on A-31) | outcome recorded in "P1 readability review record" above | — (decides K1 vs K2) |

## Decision needed from the project owner
1. **Q8: are C++ components acceptable in the code base?** Implied yes by D-027 (C++ driver decided); direct confirmation still requested (A-30).
2. **Q2** is answered by D-027. **Q11 (a)** answered: NVIDIA is the only GPU target; AMD out of scope (owner, 2026-09-25). **Q11 (c)** closed: AMR mode uses uniform grids on each level and stretched cases stay FDS-only (owner, 2026-09-25). Still open: **Q11 (b)** GPU test hardware.
3. ~~AMD: deferred or dropped?~~ **Decided** (owner, 2026-09-25): AMD is out of scope; the second-compiler check is optional and non-gating (`ifx` only).
4. **Q5: are refined-level outputs Smokeview-native, or VTK (FireX `vtkf.f90`) / AMReX plotfiles?**

## References
FireX `Source/{mesh,init,mass,main,read,dump,vtkf,pres,velo,wall,turb,vege,type,ccib,pois}.f90`, `CMakeLists.txt`, `Build/makefile` @ `36975d765f`. AMReX `Src/F_Interfaces/**`, `Src/FFT/AMReX_FFT_Poisson.H`, `Src/LinearSolvers/MLMG/{AMReX_MLLinOp,AMReX_MLCellLinOp}.H`, `Src/Boundary/AMReX_InterpBndryData.H`, `Docs/sphinx_documentation/source/{Fortran_Chapter,GPU}.rst`, `Tools/CMake/AMReXOptions.cmake`, `Tests/CMakeLists.txt`, `.github/workflows/{cuda,hip}.yml`, `CHANGES.md` @ `99ddfda`. HYPRE v2.32.0-24 `63331f19c` `HYPRE_config.h` (`(local GNU third-party library tree)/libs/hypre/63331f19`; D-026). ERF `Source/Microphysics/Morrison/*`; PeleLMeX, incflo, ERF `Source/`. Tutorials `FortranInterface/Advection_F`, `Amr/Advection_AmrCore`. Teammate docs: `docs/{risks,requirements,roadmap,charter,README,spec-responses}.md`, `docs/inventory/*`, `docs/amrex/{mapping,driver-options}.md`, `docs/pressure/0{0,1}-*.md`.
