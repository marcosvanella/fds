# Driver options: C++ AmrCore driver + FDS Fortran kernels vs FDS on AMReX F_Interfaces

| | |
|---|---|
| Author | AMReX Integration Lead |
| Status | **DRAFT** (first draft, for ADR-001) |
| Date | 2026-09-25 |
| FDS base | FireX-based local branch `AMReX` @ `36975d765f` (`36975d765fcead401e14b094a04f910ac42eab8a`), this repository (read-only) |
| AMReX base | `(local AMReX checkout)` @ `99ddfda` (read-only) |
| Companion docs | `docs/amrex/mapping.md` (concept mapping); `docs/adr/ADR-001-driver-architecture.md` (decision); `docs/pressure/00-fds-pressure-baseline.md` and `01-amr-mapping-spec.md` (pressure); `docs/inventory/*.csv`; `docs/vv/verification_case_survey.csv` |

Unless marked otherwise, AMReX paths are relative to `amrex/Src/F_Interfaces/` and FDS paths to `fds-amr/src/Source/`.
Nothing here was built or run. Every verdict comes from reading source and documentation.

ADR-001 asks the Integration Lead to confirm or dispute its F_Interfaces verdict. **Answer: confirmed.**
- The table below agrees with ADR-001 §"AMReX F_Interfaces — verified contents" and with its decision-matrix rows 2, 4 and 5.
- There is one addition, in the MLMG row: the overset mask can probably be emulated through `set_acoeffs`/`set_bcoeffs`. This is unverified.

---

## 1. Fortran-interface feature table (for ADR-001)

Verdicts:
- **works**: usable for FDS as-is.
- **partial**: usable with a noted restriction, or needs a small new binding.
- **missing**: no Fortran binding. You would have to write C++ glue or go through C++.

The F_Interfaces tree has about 6.8k lines of Fortran module code in `Base/`, `AmrCore/`, `LinearSolvers/`, `Particle/` and `Octree/`, plus thin `*_fi.cpp` glue.

| # | Feature | Verdict | What exists (AMReX citation) | Limits relevant to FDS |
|---|---|---|---|---|
| 1 | **AmrCore** (hierarchy, regrid, level callbacks) | **works** | `amrex_amrcore_module`: publics `AmrCore/AMReX_amrcore_mod.F90:11-23`; `bind(c)` callbacks `make_level`/`clear_level`/`error_est` `:26-47`; `amrex_amrcore_init` `:178-192`; `init_from_scratch`/`regrid`/`post_regrid` `:261-291`; explicit-grid setters `amrex_set_boxarray/distromap/geometry/finest_level` `:14-15,238-242`. The C++ side `FAmrCore` subclasses `AmrCore` with function pointers (`AmrCore/AMReX_FAmrCore.H`). | (i) Configured only through ParmParse (`amr.*`, `geometry.*`), because `FAmrCore` is default-constructed. FDS `&MESH` input has to be translated to ParmParse (`Base/AMReX_parmparse_mod.F90` is available). (ii) The ref ratio is a **scalar per level** (`amrex_ref_ratio(:)` `:58,186-187`; glue uses `MaxRefRatio`, `AMReX_amrcore_fi.cpp` ~20-26), whereas C++ `AmrMesh` holds `Vector<IntVect>` (`AmrCore/AMReX_AmrMesh.H:30,156-162`). FDS multi-resolution meshes with anisotropic ratios (`mapping.md` §3) cannot be expressed. (iii) CPU only. |
| 2a | **FillPatch, single level** | **works** | `amrex_fillpatch` generic `AmrCore/AMReX_fillpatch_mod.F90:14-18`; `amrex_fillpatch_single` `:421-456` → `FillPatchSingleLevel` (`AMReX_fillpatch_fi.cpp:200-285`) | The physical-BC callback `amrex_physbc_proc(mf,scomp,ncomp,time,geom)` (`Base/AMReX_physbc_mod.F90:24-29`) goes through `FPhysBC`, which ignores `nghost`/`bccomp` (`Base/AMReX_FPhysBC.cpp`) and runs on the host. FDS open/wall ghost logic must go inside this callback, one whole MultiFab at a time. |
| 2b | **FillPatch, two level** (cell and face) | **works** | `amrex_fillpatch_two` `:458-538`; `amrex_fillpatch_two_faces` `:540-673`; `amrex_fillcoarsepatch_default` `:676-727` and `_faces` `:730-825`; `amrex_fillpatcher` class with fill/fillcf and RK3/RK4 time interpolation `:220-419`; pre/post interpolation hooks `amrex_interp_hook_proc`/`_arr_proc` `:42-71` | Ratio is isotropic only: glue builds `IntVect(rr)` / `IntVect(nghost)` (`AMReX_fillpatch_fi.cpp:96`). No EB-aware interpolation. |
| 2c | **Interpolaters** | **works** | `AmrCore/AMReX_interpolater_mod.F90`: `pc`(0), `node_bilinear`(1), `cell_bilinear`(2), `quadratic`(3), `lincc`(4), `cell_cons`(5), `protected`(6), `quartic`(7), `face_divfree`(8), `face_linear`(9); mapped at `AMReX_fillpatch_fi.cpp:9-20` | You cannot register a user-defined `Interpolater` (e.g. an FDS-specific conservative scheme for scalars with a positivity clip) from Fortran. The hooks in 2b are the only customisation. `face_divfree` is available for velocity. |
| 3 | **Flux registers** | **works** (restricted) | `amrex_fluxregister`: fineadd, crseinit, crseadd, setval, reflux, overwrite (`AmrCore/AMReX_fluxregister_mod.F90:14-26`); `amrex_flash_fluxregister` (`AmrCore/AMReX_flash_fluxregister_mod.F90`) | Isotropic `IntVect(rr)` (`AMReX_fluxregister_fi.cpp:11`). `fineadd_1fab` is `RunOn::Cpu` (`:45`). Reflux uses a plain Geometry volume with no EB volume fractions (`:97`). Enough for FDS Phase 1: one global dt, no subcycling, so the flux register is needed only for conservative synchronisation of scalars at coarse/fine faces (FR-016). |
| 4 | **MLMG** (operators, BCs, bottom solvers) | **partial** | Operators: only `amrex_poisson` and `amrex_abeclaplacian`, both **cell-centred** and multi-level via `geom(0:)`, `ba(0:)`, `dm(0:)` (`LinearSolvers/AMReX_poisson_mod.F90:53`, `AMReX_abeclaplacian_mod.F90:76`). `amrex_abeclaplacian` exposes `set_scalars`, `set_acoeffs` and `set_bcoeffs` (`AMReX_abeclaplacian_mod.F90:14-16`). `set_bcoeffs` takes `bcoef(amrex_spacedim)` **face** MultiFabs (`:163-173`) → C++ `setBCoeffs` (`AMReX_abeclaplacian_fi.cpp:47`). **So a variable-coefficient, face-centred β = 1/ρ projection ∇·(β∇p) = rhs works from Fortran.** `amrex_linop`: `set_maxorder`, `set_domain_bc(lobc,hibc)`, `set_coarse_fine_bc`, `set_level_bc` (`AMReX_linop_mod.F90:12-19`). `amrex_multigrid`: bottom solvers smoother/bicgstab/cg/**hypre**/petsc (`AMReX_multigrid_mod.F90:9-14`; mapping `AMReX_multigrid_fi.cpp:82-95`); `solve`, `get_grad_solution`, `get_fluxes`, `comp_residual` (`:24-35`). BC constants: dirichlet, neumann, reflect_odd, inhomog_neumann, robin, symmetry, periodic, … (`Src/Boundary/AMReX_lo_bctypes_mod.F90`; enum `Src/Boundary/AMReX_LO_BCTYPES.H:27-39`). | **Missing:** `MLNodeLaplacian` (nodal projection), `MLEBABecLap` / all EB operators, `MLTensorOp`, the overset-mask constructor (`Src/LinearSolvers/MLMG/AMReX_MLABecLaplacian.H:41-45`; no `overset` symbol in F_Interfaces), and `MacProjector` (AMReX-Hydro, C++). **Partial:** Robin: `set_level_bc` passes no robin a/b/f, whereas C++ `setLevelBC` has them (`Src/LinearSolvers/MLMG/AMReX_MLLinOp.H:286-297`). FDS open boundaries are Dirichlet on H (`pres.f90:186-191`) and walls are Neumann, so Robin is not needed. **Overset emulation (unverified idea):** zero `bcoef` on faces touching solid cells, plus `acoef=1` in solid cells (0 in gas), turns solids into decoupled "known" cells. That would be the E-2 option (`docs/pressure/01-amr-mapping-spec.md` §E) without the mask binding. The per-zone null-space handling (spec §C.4) is still needed either way. |
| 5 | **Embedded boundaries** | **missing** | Zero EB symbols in F_Interfaces. `Src/F_Interfaces/CMakeLists.txt` adds no EB sources. | This blocks ADR-003 option EB / pressure option E-3 on the Fortran route. |
| 6 | **Particles** | **partial → effectively missing for FDS** | A single fixed particle type `amrex_particle` {pos, vel, id, cpu} (`Particle/AMReX_particlecontainer_mod.F90:18-23`) instantiated as `AmrParticleContainer<AMREX_FI_NSTRUCTREAL=BL_SPACEDIM, NSTRUCTINT=0>` (`AMReX_particlecontainer_fi.cpp:8-13`); redistribute, write (checkpoint), add/get/num particles (`mod.F90:28-35`) | FDS `LAGRANGIAN_PARTICLE_TYPE` carries a large attribute set plus a pointer to the `BOUNDARY_ONE_D` record (`type.f90:389-429`). None of that fits. You would need new C++ glue per attribute layout, or keep FDS's own particle storage and redistribute by hand. 157 of 928 verification cases use particles (`docs/vv/verification_case_survey.csv`). |
| 7 | **GPU** | **missing** | "The Fortran interface of AMReX does not currently have GPU support. AMReX recommends porting Fortran code to C++" (`Docs/sphinx_documentation/source/GPU.rst:105-106`). `amrex_fi_multifab_dataptr*` hands Fortran a raw pointer into FAB memory with no staging (`Base/AMReX_multifab_fi.cpp:50,63`). AMReX builds its Fortran-interface tests only when `AMReX_GPU_BACKEND STREQUAL NONE` ("do not work on GPU", `Tests/CMakeLists.txt:168-169`). | CMake does not forbid the combination (`Tools/CMake/AMReXOptions.cmake:289-291`), so it is unsupported rather than blocked. Under GPU you would need managed memory, and the host callbacks (physbc, `fineadd_1fab`) would stay on the CPU. |
| 8 | **AmrLevel** (`Amr` class, `StateData`, `Derive`) | **missing** (not needed) | `Docs/sphinx_documentation/source/Fortran_Chapter.rst:10-13`: the Fortran interface covers everything **except AmrLevel and particles**. `AmrCore/AMReX_amr_mod.F90` only aggregates modules. | Irrelevant if we follow the AmrCore pattern on either route. AmrLevel gives built-in subcycling/checkpoint, which Phase 1 (single global dt) does not need. |
| 9 | **Tagging** | **works** | `amrex_tagboxarray` with `dataPtr` (`AmrCore/AMReX_tagbox_mod.F90:11-15`), passed to the `error_est` callback (`AMReX_amrcore_mod.F90:40-47`) | Tag values are `character(c_char)` per cell. The FDS criteria (gradients, HRRPUV, near-OBST buffer) are ordinary Fortran loops. |
| 10a | **Plotfile** | **works** | `amrex_write_plotfile`, `amrex_write_hdf5plotfile` (`Base/AMReX_plotfile_mod.F90:15-17`). HDF5 only under `AMREX_USE_HDF5` (`Base/AMReX_plotfile_fi.cpp:5,32`), i.e. `AMReX_HDF5=ON` (`Tools/CMake/AMReXOptions.cmake:366`). | Plotfiles are not Smokeview input. The output path is its own work item (`mapping.md` §9). |
| 10b | **Checkpoint/restart** | **partial** | `amrex_multifab_write`/`read` via VisMF (`Base/AMReX_vismf_fi.cpp`); particle `write` (row 6) | There is no Fortran helper to write/read the checkpoint Header (BoxArray per level, time, step). `rg` finds no write/read in `boxarray_mod`, `distromap_mod` or `geometry_mod`, so the header I/O and the BoxArray rebuild on restart have to be hand-rolled (small glue). |
| — | MultiFab (cell, face, nodal) and MFIter | **works** | `amrex_multifab_build(…, nc, ng, nodal)` (C binding takes `ng(3)`, `nodal(3)`: `Base/AMReX_multifab_mod.F90:177-182`). fill_boundary, parallel_copy, sum_boundary, override_sync, average_sync `:56-62`. Owner mask `:22`. `amrex_mfiter` tilebox/nodaltilebox/growntilebox `:140-154`. `dataPtr` returns a pointer with **global** lower bounds `dp(bx%lo(1):,…)` `:758-786`. Average down cell/node/faces (`Base/AMReX_multifabutil_mod.F90:11-14`). | This is what makes the POINT_TO_MESH shim possible on either route (§4). |
| — | ParmParse, init with a user communicator | **works** | `Base/AMReX_parmparse_mod.F90`. `amrex_init(comm, arg_parmparse, proc_parmparse)` duplicates the communicator (`Base/AMReX_init_mod.F90:20-62`; `Base/AMReX_parallel_mod.F90:65-88`). | AMReX's Fortran MPI glue uses legacy `use mpi` (`Base/AMReX_fi_mpi_mod.F90:2`). FDS uses `USE MPI_F08`. They coexist: pass `MPI_COMM_WORLD%MPI_VAL`. |
| — | Octree AMR, Runge–Kutta helpers | **works** | `Octree/AMReX_octree_mod.F90`; rk2/rk3/rk4 `Base/AMReX_rungekutta_mod.F90:10-13` | Not needed for FDS (FDS uses its own predictor–corrector). |

**Maintenance signal** (from ADR-001; I re-checked the `CHANGES.md` lines):
- The last feature work on the F_Interfaces layer was in 24.10 (#4124, #4115). The face FillPatch bindings came in 23.10 (#3541–#3553). 26.09 has a bugfix (#5604).
- Reading: the layer is **maintained but not developed**. Any new AMReX capability we want (EB, overset mask, nodal projection, particles with attributes, GPU) means writing our own bindings. **Unverified:** this has not been confirmed with the AMReX developers.

### 1.1 What the table means for the pressure path (corrected)

The team plan replaces `PRESSURE_SOLVER`/`PRESSURE_ITERATION_SCHEME` (`main.f90:1601-1745`) wholesale with one composite MLMG solve at one global dt. For that plan:

- **H-form, constant coefficient** (FDS FFT/GLMAT analogue, E-1 IBM forcing): `amrex_poisson`. **Works from Fortran.**
- **Variable coefficient, β = 1/ρ on faces** (a true density-weighted projection): `amrex_abeclaplacian` + `set_scalars(0,1)` + `set_bcoeffs(beta_faces)`. **Works from Fortran.** Variable-coefficient projections do **not** force C++.
- **Solids removed from the operator** (E-2, UGLMAT analogue): the overset mask is not bound. The likely Fortran alternative is β=0 on solid faces plus acoef in solid cells (row 4). A small C++ binding for the mask constructor (~50 lines) is the other way. Either way this is **not** a forcing issue.
- **EB cut-cell operator** (E-3; ADR-003 EB): **forces C++.**
- **Nodal (cell-vertex) projection**: **forces C++**. FDS pressure is cell-centred with face velocities (MAC layout, `mapping.md` §7), so nodal projection is not required.
- **HYPRE**: available as MLMG bottom solver on both routes (`AMReX_multigrid_fi.cpp:82-95`).
- **GPU for the solver**: only on the C++ route, as MLMG on device. FireX shows the shape of this today (below).

So the pressure path pushes toward C++ **only if ADR-003 chooses embedded boundaries** (or if we later want the solver on the GPU through AMReX).

---

## 2. The two options

### Option (a): C++ AmrCore driver calling FDS Fortran kernels via `bind(C)` on per-box pointers

**Division of work.**
- A C++ class derived from `amrex::AmrCore` owns:
  - the hierarchy, regrid, `FillPatch`/`FillBoundary`, flux registers and average-down;
  - the MLMG pressure solve (or `MacProjector`);
  - particles, checkpoint/plotfile, and the time loop (a translation of `MAIN_LOOP`, `main.f90:695-1213`, and its predictor/corrector halves at `:748` and `:928`).
- FDS physics stays in Fortran. Each kernel becomes `BIND(C)` and receives `(lo, hi, array pointers + bounds, dx, …)` for one box or tile inside an `MFIter` loop.
- This is AMReX's documented pattern for Fortran kernels: `Docs/sphinx_documentation/source/Basics.rst:2724-2779` (sec:basics:fortran; `BL_TO_FORTRAN_BOX/ANYD`).
- ERF uses the same pattern for WRF microphysics (per ADR-001: `ERF_AdvanceMorrison.cpp:204-208`, whole-FAB with tiling off on the Fortran path; pinned arena under GPU `:264-268`; not re-checked by me).

**Pros**
- Every AMReX capability is reachable:
  - EB (`MLEBABecLap`, EB flux registers);
  - `MacProjector`, the overset mask and nodal operators;
  - `ParticleContainer<…>` with arbitrary real/int attribute counts (compile-time SoA/AoS), e.g. `AmrParticles` (`Src/AmrCore/AMReX_AmrParticles.H`);
  - anisotropic `IntVect` ref ratios (`AMReX_AmrMesh.H:30,156-162`);
  - full checkpoint header I/O as in the AmrCore tutorials;
  - MLMG and later kernels on the GPU.
- GPU is incremental. Hot kernels are ported one at a time to `amrex::ParallelFor` while the rest stay host-Fortran (needs a pinned/managed arena for those FABs, like ERF). The solver goes to the device first, which is the direction FireX already takes with HYPRE (`pres.f90:1177-1184,4431-4438`; §5).
- It tracks AMReX development, so no binding debt accumulates.
- OpenMP comes from `MFIter` tiling if the kernel accepts a tile box (`lo`/`hi` inside a larger FAB). FDS kernels written as `DO K=1,KBAR` loops over a whole mesh need `lo/hi` loop bounds to use tiles. Until then run `TilingIfNotGPU()` off and keep the existing `!$OMP` inside kernels (349 `!$OMP` lines in `Source/*.f90`, my count; ADR-001 counts 359).

**Cons**
- There are two languages at the top level. FDS's `MAIN_LOOP`, `CHANGE_TIME_STEP_LOOP` (`main.f90:774-897`), the dump schedule and `INITIALIZE_*` sequencing have to be re-expressed in C++ or called as coarse Fortran "phase" routines. The second is feasible as a transition: C++ calls `FDS_PHASE_X(level)` and the Fortran side loops over boxes through the shim.
- `bind(C)` interfaces need interoperable argument types. FDS derived types (`WALL_TYPE`, `BOUNDARY_ONE_D_TYPE`, `LAGRANGIAN_PARTICLE_TYPE`) cannot cross the boundary. They stay Fortran-owned and are referenced by box index (see §4).
- The Fortran module state (`MESH_VARIABLES`, `GLOBAL_CONSTANTS`, `READ_INPUT` results) stays global in Fortran. C++ has to call a Fortran `READ_DATA` and query sizes through C getters. This works but is plumbing-heavy.
- The build becomes mixed-language (§6).

### Option (b): FDS stays a Fortran program on `Src/F_Interfaces`

**Division of work.**
- FDS keeps `PROGRAM FDS` and `MAIN_LOOP`.
- It calls `amrex_init`, `amrex_amrcore_init` and `amrex_init_from_scratch`, and it provides `make_level`/`error_est` callbacks in Fortran.
- It stores fields in `amrex_multifab` and loops `amrex_mfiter`. It gets per-box Fortran pointers with global bounds (`dataPtr`, `AMReX_multifab_mod.F90:758-786`) and calls the existing kernels through the shim (§4).
- Pressure goes through `amrex_poisson` or `amrex_abeclaplacian` + `amrex_multigrid`.
- The tutorial `ref/amrex-tutorials/ExampleCodes/FortranInterface/Advection_F` shows the full pattern (subcycling + reflux in `Source/evolve_mod.F90:68-137`).

**Pros**
- One language for FDS developers, who are all Fortran developers. `PROGRAM FDS` and the input/output flow stay recognisable.
- For a **non-EB, CPU-only, particle-free or particle-light** Phase 1 everything needed exists:
  - AmrCore with explicit or tagged grids;
  - cell and face FillPatch with `face_divfree`;
  - flux registers;
  - cell-centred composite MLMG with variable face β and a HYPRE bottom solver;
  - plotfile output.
- The pressure path does **not** force C++ (§1.1).
- Smaller build change: add the AMReX library and `enable_language(CXX)` only for linking, since the F_Interfaces glue is compiled inside AMReX.

**Cons**
- **Hard walls:**
  - no EB (ADR-003 EB option is closed);
  - no GPU, officially (`GPU.rst:105-106`);
  - no attribute-bearing particles;
  - no nodal projection or overset mask without writing bindings;
  - isotropic ref ratio only.
- **Binding debt:** every missing item means writing C++ glue plus a Fortran module inside our tree against a layer AMReX maintains but does not develop. At that point we are writing C++ anyway, just in a worse place.
- **Performance traps:**
  - the physical-BC callback is whole-MultiFab, host-only, with `nghost` ignored (`Base/AMReX_FPhysBC.cpp`);
  - `fineadd_1fab` is CPU (`AMReX_fluxregister_fi.cpp:45`).
  Both are fine on CPU.
- AMReX's own tests exclude this path on GPU (`Tests/CMakeLists.txt:168-169`). The ecosystem of reference codes on it is thin: the Fortran tutorials, and `FlashFluxRegister`, whose name suggests a FLASH-X user (not checked).

### Side-by-side

| Criterion | (a) C++ driver + Fortran kernels | (b) Fortran on F_Interfaces |
|---|---|---|
| Phase-1 scope (non-EB, CPU, one global dt, composite MLMG) | yes | **yes** |
| Variable-β (1/ρ) projection | yes (`MLABecLaplacian`, `MacProjector`) | **yes** (`set_bcoeffs`, `AMReX_abeclaplacian_mod.F90:14-16,163-173`) |
| Solids removed from operator (E-2) | yes (overset mask) | probably, via β=0 / acoef emulation (unverified), or ~50-line binding |
| EB (ADR-003 EB, E-3) | yes | **no**: forces C++ |
| Nodal projection | yes | no: forces C++ (not needed for FDS's MAC layout) |
| Particles with FDS attributes | yes (templated container) | no (pos/vel only) |
| GPU (solver) | yes (MLMG/HYPRE on device) | no (officially unsupported) |
| GPU (kernels) | incremental C++ port | no |
| Anisotropic ref ratio | yes | no |
| Checkpoint | full | MultiFab only; header hand-rolled |
| Developer familiarity | mixed | high |
| Binding maintenance | none (native API) | ours, for every missing feature |
| Migration shim (§4) usable | yes | yes |
| Reversibility | the shim + Fortran kernels are reusable from (b) | the kernels + shim are reusable from (a) |

The last row matters. The expensive work is the same on both routes:
- converting ~200 `POINT_TO_MESH` call sites (216 `CALL POINT_TO_MESH` in `Source/*.f90`, my count);
- removing 382 `MESHES(NOM)` cross-mesh references;
- restructuring the WALL/OBST data model.

So the driver choice is cheaper to revisit than it looks, **if** the kernel interface is designed as "box + pointers" from day one.

---

## 3. Confirm/dispute ADR-001

- **Confirmed:** ADR-001's F_Interfaces inventory (lines 64-74), its decision-matrix rows 2, 4 and 5, and its statement that the decision "does not rest on the pressure path" (line 166).
- **Confirmed:** ADR-001's GPU evidence (`multifab_fi.cpp:50,63`, `Tests/CMakeLists.txt:168-169`). I re-read both.
- **Addition:** E-2 (overset-mask-like masking) is not strictly a binding gap. It can probably be expressed through `acoef`/`bcoef` (row 4). This is unverified, and prototype P2 below tests it.
- **Addition:** the ref ratio is scalar per level on the Fortran route. FDS multi-resolution setups with different refinement per direction (`mapping.md` §3) are C++-only.
- **Recount note:** my counts differ slightly from ADR-001 because of regex differences, not substance:

| Count | Mine | ADR-001 |
|---|---:|---:|
| `CALL POINT_TO_MESH` | 216 | 200 |
| `=>` lines in `POINT_TO_MESH` (`mesh.f90:505-916`) | ~403 | 368 |
| `MESHES(NM)%` | 3,037 | 2,926 |
| `MESHES(NOM)` | 382 | 381 |
| `!$OMP` | 349 | 359 |
| `CALL MPI_` | 436 | 369 |

---

## 4. The POINT_TO_MESH shim (transition mechanism, usable from either option)

**Idea.** Make "one FDS mesh" mean "one AMReX box" for the duration of a kernel call. The Fortran pointers in `MESH_POINTERS` (`mesh.f90:361`; `REAL(EB), POINTER` declarations of `U,V,W,…,RHO,…` at `:367-376`) are pointed at FAB data with Fortran 2003 **bounds remapping**, so the unmodified kernel body still sees FDS-local indices. `POINT_TO_MESH` (`mesh.f90:505-916`) is the natural place for this: it already does ~400 `U=>M%U`-style associations (`mesh.f90:513`, `:534`, …; the full list is in `docs/inventory/point_to_mesh_pointers.csv`).

**Feasibility: yes, with three preconditions.**

1. **Index mapping.** `amrex_multifab%dataptr(mfi)` returns a contiguous pointer whose lower bounds are the **global** indices of the grown FAB box (`Base/AMReX_multifab_mod.F90:758-786`). Remapping only the lower bounds gives FDS local indices. For a box with valid cells `lo..hi`:
   - cell-centred, ng=2 (e.g. `RHO(-1:IBP1+1,…)`, `init.f90:525`): `RHO(-1:,-1:,-1:) => dp` gives FDS `I` ↔ global cell `lo+I-1`;
   - x-face, nodal=(1,0,0), ng=1 (e.g. `U(-1:IBP1,0:JBP1,0:KBP1)`, `init.f90:533`): `U(-1:,0:,0:) => dp` gives FDS face `I` ↔ global face `lo+I` (FDS face `I` is the forward face I+½);
   - `IBAR = hi(1)-lo(1)+1` and so on, set per box.

   **Precondition:** every MultiFab's ghost width must equal the FDS allocation exactly (RHO/ZZ/WORK_PAD ng=2 at `init.f90:525-529`; U/V/W ng=1 with the extra staggered ghost at `:533-536`; H ng=1 at `:556`). If they differ, the pointer extents are wrong. `docs/inventory/mesh_fields.csv` (bounds column) and `lbound_extent_dependencies.csv` provide the per-array table. Multi-component fields such as `ZZ(:,:,:,N)` map to `ncomp` MultiFabs, with the 4-D `dataptr` rank matching.
2. **ALLOCATABLE → POINTER in MESH_TYPE.** All 195 array components of `MESH_TYPE` are `ALLOCATABLE` (`mesh.f90:16-354`, e.g. `U` at `:18`, `RHO` at `:39`), and 3,037 references go through `MESHES(NM)%…`/`M%…` without passing through the `MESH_POINTERS` pointers. For those references to see FAB data, the components that migrate into MultiFabs must become `POINTER, CONTIGUOUS`. Then `ALLOCATED()` → `ASSOCIATED()`, and `DEALLOCATE`/`MOVE_ALLOC` sites change. This is a mechanical but wide edit (see `docs/inventory/allocation_sites.csv`).
3. **No cross-mesh access inside kernels.** The 382 `MESHES(NOM)` references (OMESH buffers, `MESH_EXCHANGE` `main.f90:3117-3975`, wall-to-neighbour lookups in `init.f90`) must be replaced. Coarse/fine and neighbour data then comes only through ghost cells filled by `FillPatch`/`FillBoundary` before the kernel runs.

**Costs and limits (these are the "shim tax").**
- **One box at a time per rank.** The module-level `MESH_POINTERS` are global state, so the shim cannot run two boxes concurrently on one rank:
  - no `MFIter` OpenMP tiling (`!$omp parallel` around the `MFIter` loop would race on the pointers);
  - OpenMP parallelism has to come from **inside** the kernels, i.e. the existing `!$OMP` directives (349 lines), exactly as FDS does per mesh today.
  - A future tile-parallel version needs the kernels to take arrays as dummy arguments (the Option-(a) end state), not module pointers.
- **No GPU** on the shim path. Module-pointer kernels are host Fortran. GPU requires a C++ `ParallelFor` port of a kernel, with its data passed as `Array4`, not module pointers.
- **Per-box side data is rebuilt on every regrid.** FDS precomputes mesh-shaped side structures, all of which assume a fixed mesh:
  - `WALL`/`WALL_INDEX`, `EXTERNAL_WALL`, `CELL_INDEX`/`CELL` (`mesh.f90:225-226`; `CELL_TYPE` `type.f90:2175-2187`), `OBSTRUCTION` subsets;
  - coordinate arrays `X`,`XC`,`RDX`,`RDXN`,… (stretched grids; 17 stretched verification cases);
  - `OMESH` buffers;
  - `INIT_WALL_CELL` `init.f90:2975-3386`, `OPEN_AND_CLOSE` `:4507`, `REASSIGN_WALL_CELLS` `:4898`.

  They must be rebuilt per box in `make_level`/`remake_level` (Option b callbacks `AMReX_amrcore_mod.F90:26-47`; Option a `MakeNewLevelFrom*`/`RemakeLevel`). The expensive part is the **1-D wall/solid-phase state** (`BOUNDARY_ONE_D`, surface temperatures and profiles). It lives in wall records, not in fields, so it cannot be FillPatched. It needs a regrid-stable wall ID and explicit re-homing, conserving the stored energy and mass, whenever a wall cell changes box or level. This is the hardest part of the whole migration and is route-independent (`mapping.md` §5; prototype P3).
- **Setup amortisation.** Per-box setup (pointer association, side data lookup) costs O(boxes). Use a large `amr.max_grid_size` (e.g. 64–128) and `amr.blocking_factor` ≥ 8 so there are few, large boxes. FDS meshes are typically 32³–64³ per rank, so this matches FDS decomposition habits.
- **Scalar mesh metadata.** `IBAR/JBAR/KBAR`, `IBP1…`, `XS/XF`, `DX`, `NM`-keyed lookups and `CELL_COUNT` offsets are per-mesh module scalars. The shim sets them per box, and anything that indexes global arrays by `NM` (restart files, device output, `MESH_ID` dependencies: `docs/inventory/mesh_id_dependencies.csv`, 166 rows) must switch to (level, box) keys.
- **Index-offset constants.** Staggered arrays have asymmetric extra ghosts (U is `-1:IBP1` in x but `0:JBP1` in y). Each face MultiFab therefore needs its own nodal flag and ng, and code that uses `LBOUND`/`UBOUND` or hard-coded `0`/`IBP1` limits (`lbound_extent_dependencies.csv`) must be checked against the box.

**Verdict.** The shim is a sound **transition** mechanism for CPU. It lets the kernels be validated bit-for-bit against FDS before any kernel is rewritten (P1). It is not an end state: it forecloses tiling and GPU and inherits the global-state design. Plan the kernel-by-kernel exit from the shim, to explicit dummy-argument `bind(C)` kernels, as part of Phase 2.

---

## 5. GPU implications

- **FireX already runs "host Fortran kernels + GPU linear solver".**
  - `DEFINE_RS_COMM_INFO` (`main.f90:5096-5146`) groups ranks into resource sets per GPU (`FDS_RANKS_PER_GPU`).
  - `GLMAT_SOLVER` (`pres.f90:3299`) gathers the RHS to the resource-set master (`:3411-3419`). Only masters solve (`:3422`). HYPRE vectors are migrated to and from the device (`:3453-3458`, `:3469-3473`), and the result is scattered back (`:3485-3492`).
  - The device policy is set at `pres.f90:1177-1184,4431-4438`, and `HYPRE_DEVICE_RUN` defaults to `.TRUE.` (`cons.f90:569`).
  - Build flags: `USE_HYPRE_NVIDIA/AMDGPU/INTELGPU` (`CMakeLists.txt:16-18`) turn on `HYPRE_ENABLE_CUDA/HIP/SYCL` (`:181-185`) and `WITH_HYPRE_DEVICE` (`:207`).
- **Option (a)** keeps and improves this shape. MLMG itself runs on the device (AMReX GPU build), with HYPRE-on-device as the bottom solver. Host kernels read/write MultiFabs in a managed or pinned arena (as ERF does). Hot kernels (advection `VELOCITY_FLUX` `velo.f90:563`, `MASS_FINITE_DIFFERENCES` `mass.f90:20`, `DENSITY` `:365`) are then ported to `ParallelFor` in priority order, removing host↔device traffic step by step.
  - **Caveat:** a GPU AMReX build makes **every** MultiFab default to device memory. Host-Fortran kernels need `The_Pinned_Arena()`/managed arenas for their FABs. Until kernels move, the solve pays a copy per pressure iteration. That is the same trade-off FireX already makes with HYPRE today, so it is not a regression.
- **Option (b)** has no GPU route. `GPU.rst:105-106`, AMReX tests are CPU-only (`Tests/CMakeLists.txt:168-169`), the physbc callbacks run on the host, and `fineadd_1fab` is CPU. Even the solver-only GPU pattern would require AMReX to be built for GPU with Fortran reading device pointers from `dataptr` (`Base/AMReX_multifab_fi.cpp:50,63`), which is unsupported.
- **Resource sets vs AMReX.** AMReX's GPU model is one rank per GPU (or `amrex.the_arena_is_managed` with oversubscription). FireX's `FDS_RANKS_PER_GPU` gather-to-master pattern does not carry over to MLMG. With MLMG on device, every rank owns device data. Hosts with more CPU cores than GPUs would need MPS or fewer, fatter ranks with OpenMP. **Unverified:** the performance of that trade.

---

## 6. CMake integration

**Current FDS build** (FireX `CMakeLists.txt`):
- Fortran-only (`project(… LANGUAGES Fortran)` `:8-12`, `enable_language(Fortran)` `:13`, cmake ≥ 3.24 `:1`).
- Options `USE_HYPRE` (`:15`), `USE_HYPRE_{NVIDIA,AMDGPU,INTELGPU}` (`:16-18`), `USE_SYSTEM_HYPRE` (`:19`), `USE_SUNDIALS` (`:21`), `USE_OPENMP` (`:24`).
- HYPRE is fetched and built in-tree by default (`FetchContent` `:157-192`, pinned `GIT_TAG` `:164`, `HYPRE_ENABLE_FMANGLE CAPS` `:171`) or found as `HYPRE 2.32.0` (`:196`) → `HYPRE::HYPRE` (`:203`).
- MPI via `MPI::MPI_Fortran` (`:101-102`). Optional MKL (`:111-119`).
- `Source/vtkf.f90` is compiled (`:62`), but HDF5 is **not** wired in CMake. `WITH_HDF5` appears only in `Build/makefile:135`.

**AMReX side** (`Tools/CMake/AMReXOptions.cmake`):
- `AMReX_SPACEDIM` (`:25`, default 3), `AMReX_FORTRAN` default **OFF** (`:97`), `AMReX_PRECISION` DOUBLE (`:112`), `AMReX_GPU_BACKEND` NONE|SYCL|CUDA|HIP (`:124-125`), `AMReX_MPI` ON (`:265`), `AMReX_OMP` OFF (`:276`), `AMReX_EB` OFF (`:286`), `AMReX_FORTRAN_INTERFACES` dependent on `AMReX_FORTRAN`, default OFF (`:289-291`), `AMReX_LINEAR_SOLVERS` ON (`:293`), `AMReX_FFT` OFF (`:306`), `AMReX_PARTICLES` ON (`:312`), `AMReX_HYPRE` OFF (`:356`), `AMReX_HDF5` OFF (`:366`).
- F_Interfaces sources are added when enabled (`Src/CMakeLists.txt:135-137`). LinearSolvers/Particle bindings are added only if those components are on (`Src/F_Interfaces/CMakeLists.txt`). Turning on F interfaces adds the MPI Fortran component (`Tools/CMake/AMReXParallelBackends.cmake:18-28`).
- An installed AMReX exports `AMReX::amrex_3d` via `find_package(AMReX)`. The config file errors if AMReX was built with Fortran and the consumer did not enable Fortran (`Tools/CMake/AMReXConfig.cmake.in:181-185`), and it finds MPI Fortran (`:196-201`).
- `AMReX_HYPRE` uses AMReX's own `FindHYPRE.cmake` and links a target named `HYPRE` (`Tools/CMake/AMReXThirdPartyLibraries.cmake:147-159`). With CUDA it also links cuSPARSE/cuRAND (`:149-155`).

**Option (a) build sketch**
```cmake
project(fds VERSION 6.9.1 LANGUAGES Fortran CXX)   # add CXX (and CUDA/HIP when GPU)
set(AMReX_SPACEDIM 3) ; set(AMReX_MPI ON) ; set(AMReX_OMP ${USE_OPENMP})
set(AMReX_EB OFF)         # ON only if ADR-003 picks EB
set(AMReX_HYPRE ${USE_HYPRE}) ; set(AMReX_FORTRAN OFF)   # FDS Fortran does not need AMReX Fortran modules
set(AMReX_GPU_BACKEND NONE)  # CUDA/HIP/SYCL later
add_subdirectory(external/amrex)   # or find_package(AMReX 26.09 REQUIRED)
add_executable(fds Source/driver/main.cpp Source/driver/FdsAmr.cpp ${FDS_FORTRAN_SOURCES})
target_link_libraries(fds PRIVATE AMReX::amrex_3d MPI::MPI_Fortran MPI::MPI_CXX HYPRE::HYPRE)
set_target_properties(fds PROPERTIES LINKER_LANGUAGE CXX)   # C++ main; add Fortran runtime libs
```
- **Linker language.** With a C++ `main`, CMake needs the Fortran runtime libraries on the link line. CMake handles this with `CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES` when both languages are enabled. With GPU, the CUDA/HIP linker must drive the link (AMReX's `setup_target_for_cuda_compilation`).
- **HYPRE must be one build shared by FDS's Fortran calls and AMReX.** FDS's fetched HYPRE (`FMANGLE CAPS`) has to satisfy AMReX's `FindHYPRE` (version ≥ 2.20, target `HYPRE`). Options: point AMReX at the FDS-built HYPRE install, or build HYPRE once externally with `USE_SYSTEM_HYPRE=ON`. **Unverified:** whether the FetchContent `HYPRE` target satisfies AMReX's `find_package(HYPRE)` inside the same configure. The likely answer is "no", which would require building HYPRE out-of-tree first. Also: HYPRE device builds must match `AMReX_GPU_BACKEND`, and integer width (`HYPRE_BigInt`) must be consistent.
- **Precision/dimension.** FDS is double (`EB`) and 3-D. Use `AMReX_PRECISION=DOUBLE`, `AMReX_SPACEDIM=3`. For FDS 2-D cases (287 in the V&V survey), keep 3-D with one cell in y (as FDS does today). A separate `amrex_2d` build is not needed.
- `AMReX_FFT` can stay OFF (MLMG only). Turn it ON only if a spike compares AMReX's FFT Poisson on uniform single-level grids.

**Option (b) build sketch.** The same, but with `AMReX_FORTRAN=ON` and `AMReX_FORTRAN_INTERFACES=ON`. The Fortran module files (`amrex_*_module.mod`) are then produced by AMReX's compiler, so **AMReX and FDS must use the same Fortran compiler and version** (`.mod` files are not portable). `CXX` is still needed in `project()` for linking the glue, and the program stays `PROGRAM FDS`. The GNUmake route in the tutorials uses `USE_F_INTERFACES=TRUE`. The CMake example gate is `AMReX_FORTRAN_INTERFACES` (`ref/amrex-tutorials/ExampleCodes/CMakeLists.txt:247-249`).

**Both options.**
- Keep `Build/makefile` working for FDS-only builds during transition. Put AMReX behind `option(USE_AMREX … OFF)` with `#ifdef WITH_AMREX` guards, as with `WITH_HYPRE`.
- Wire HDF5 into CMake if `vtkf.f90` output is kept (`mapping.md` §9). It is currently makefile-only.

---

## 7. Recommendation (**DRAFT**, for ADR-001; not a decision)

> **DRAFT recommendation: Option (a), a C++ `AmrCore` driver calling FDS Fortran kernels via `bind(C)` on per-box pointers, with the POINT_TO_MESH shim (§4) as the Phase-1 migration mechanism.** This agrees with ADR-001 (Option A, with Option C "mesh-view shim" as the migration path).

Rationale, stated honestly:
1. **The pressure path is not the reason.** The team plan (one composite MLMG solve replacing `PRESSURE_SOLVER`/`PRESSURE_ITERATION_SCHEME`, one global dt, no subcycling) can be driven from Fortran. `amrex_poisson` covers constant-coefficient H, and `amrex_abeclaplacian` with `set_bcoeffs` covers a face-centred β = 1/ρ, ∇·(β∇p) solve (`AMReX_abeclaplacian_mod.F90:14-16,163-173`). Only **EB operators and nodal projection** force C++, so the pressure path pushes toward C++ **only if ADR-003 chooses embedded boundaries**.
2. **What does push toward C++:**
   - GPU: none on F_Interfaces, and FireX is already heading toward device solvers (§5);
   - particles with FDS attributes: 157/928 V&V cases;
   - EB, if ADR-003 picks it;
   - anisotropic ref ratios;
   - the maintained-not-developed status of F_Interfaces, which turns each of the above into our own binding debt.
3. **Option (b) remains viable** for a non-EB, CPU-only Phase 1, and the shared work (shim, `MESHES(NOM)` removal, wall data model) is reusable. If the project owner rules GPU, EB and AMReX-managed particles permanently out of scope, ADR-001's Spike S2 overturn condition applies and (b) should be reopened.
4. **Kernel-interface rule for both routes:** new or refactored kernels take arrays and `lo/hi` as dummy arguments (`bind(C)`-compatible). They should not rely on `MESH_POINTERS` module state, so the shim can be retired kernel by kernel and the driver choice stays reversible.

---

## 8. Proof-of-concept prototypes (ranked; **not built**)

### P1: shim + one kernel on a static two-level C++ AmrCore (highest priority)
- **Scope:**
  - A minimal C++ `AmrCore` subclass with two levels, static fine patch, ratio 2, CPU, MPI.
  - Fortran `bind(C)` wrappers `fds_box_begin(lev, lo, hi, dx, ptrs…)`/`fds_box_end` that do the §4 bounds remapping for `RHO`, `ZZ`, `U/V/W`, `D` and the coordinate arrays.
  - Calls unmodified `MASS_FINITE_DIFFERENCES` (`mass.f90:20`) and `DENSITY` (`mass.f90:365`) per box, with `FillPatchTwoLevels` (cell_cons for scalars, face_divfree for velocity) providing ghosts and `average_down` + flux register for synchronisation.
  - No pressure solve (prescribed divergence-free velocity), no walls except periodic/open.
- **Success criteria:**
  1. Single level, single box, same grid as an FDS mesh: `RHO`, `ZZ` after one step match FDS **bit-for-bit** (same compiler flags).
  2. Single level, N boxes: matches the one-box result to round-off. This validates ghost widths and face index mapping (FDS `I` ↔ `lo+I` for faces, `lo+I-1` for cells).
  3. Two levels: total species mass is conserved to round-off with reflux on. A passive blob crosses the coarse/fine interface without visible artefacts.
  4. The per-box shim overhead is measured (target < 5 % of kernel time at `max_grid_size` 64).
- **Answers:** whether the shim works without editing kernel bodies; the ALLOCATABLE→POINTER fallout; real ghost-width requirements.

### P2: composite MLMG pressure solve, driven from both C++ and Fortran
- **Scope:**
  - Build the FDS pressure RHS on MultiFabs from a snapshot of FDS fields.
  - Solve with MLMG three ways:
    - (i) `MLPoisson` in C++ (H-form, E-1);
    - (ii) `MLABecLaplacian` in C++ with face β = 1/ρ;
    - (iii) the same two solves via `amrex_poisson`/`amrex_abeclaplacian` from Fortran (this doubles as ADR-001's Spike S2 for the pressure path).
  - Also test the E-2 emulation (β=0 on solid faces + `acoef` in solids, §1 row 4) against the C++ overset-mask constructor.
  - Bottom solvers: native and HYPRE.
- **Success criteria:**
  1. Single level, uniform grid: the solution matches FDS FFT (`pois.f90` via `PRESSURE_SOLVER`) to solver tolerance on at least three standard cases (open, periodic, closed box with the compatibility condition). The Fortran-route solution is identical to the C++-route solution.
  2. Two levels: the composite solution is continuous across the coarse/fine interface, and the projected velocity divergence is ≤ `VELOCITY_TOLERANCE` without the FDS iteration loop (compare `docs/pressure/01-amr-mapping-spec.md` §A.1).
  3. E-2 emulation reproduces the overset-mask solution. Disconnected zones are handled (spec §C.4).
  4. Time to solution is recorded against FFT and against FireX GLMAT/HYPRE on the same case. There is no pass threshold for this first measurement, but it feeds ADR-001/roadmap.
- **Answers:** whether the Fortran route really covers Phase-1 pressure (confirming §1.1); MLMG cost vs FFT.

### P3: wall/OBST data model under regrid
- **Scope:**
  - Walls on a two-level hierarchy with one regrid that moves a fine patch across an `OBST`.
  - Build `CELL_INDEX`/`WALL` per box in `MakeNewLevelFromScratch`/`RemakeLevel` from a regrid-stable geometric wall ID (face coordinates + orientation).
  - Re-home `BOUNDARY_ONE_D` state (1-D temperature profile) when a wall cell changes level: split on refinement, energy-weighted merge on coarsening.
- **Success criteria:**
  1. Wall and solid energy are conserved to round-off across regrid.
  2. Surface temperature fields are continuous before and after regrid.
  3. The rebuild cost per regrid is measured relative to one time step.
- **Answers:** the largest route-independent risk (`mapping.md` §5 and top risks).
- **Alternative P3 (if the output team prefers):** an output adapter writing AMReX hierarchies to FireX VTKHDF (`vtkf.f90`) or plotfile + a Smokeview-compatible slice path (`mapping.md` §9).

---

## 9. Open items and unverified statements

- The maintenance status of F_Interfaces ("maintained, not developed") is inferred from `CHANGES.md`. It has not been confirmed with the AMReX developers.
- The E-2 emulation through `acoef`/`bcoef` is reasoned, not tested (P2).
- MLMG vs FFT/GLMAT performance on FDS cases: unmeasured (P2).
- Whether FDS's FetchContent HYPRE can satisfy AMReX's `FindHYPRE` in one configure: untested.
- Whether a FireX CMake build compiles `vtkf.f90` without `WITH_HDF5`: not built.
- ERF line citations are taken from ADR-001 and were not re-checked.
- `docs/inventory/*.csv` line numbers: spot-checked against FireX (`mesh.f90:18` for `U` matches). They are not fully re-validated.
- The GPU resource-set vs one-rank-per-GPU performance trade-off (§5): unmeasured.
