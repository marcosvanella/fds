# FDS → AMReX concept mapping (first draft)

Owner: AMReX Integration Lead. Status: **DRAFT for discussion**, 2026-09-25.

**Base commits**
- **FDS: FireX-based local branch `AMReX` @ `36975d765f`** (this repository, read-only). It is roughly 8.6k lines ahead of master `ce1f659` (`git diff --stat ce1f659 HEAD -- Source`: 22 files, +8,621/−847), mainly in `pres.f90`, `main.f90`, `dump.f90`, `read.f90`, `radi.f90` and the new `vtkf.f90`. All FDS citations below were rechecked against this tree.
- **AMReX: `(local AMReX checkout)` @ `99ddfda`** (shared clone, read-only).
- Reference apps: `(local reference-code clones)/{incflo,PeleLMeX,IAMR,ERF,amrex-tutorials}`.

FDS paths are relative to `fds-amr/src/Source/`. AMReX paths are relative to `amrex/Src/` unless written in full. Line numbers refer to the commits above.

Cross-references (not repeated here):
- Field-by-field inventory of `MESH_TYPE` (431 fields; 1,314 across all derived types): `docs/inventory/mesh_fields.csv`. `POINT_TO_MESH` call sites: `docs/inventory/point_to_mesh_calls.csv`. Note: `mesh.f90` gained 8 lines in FireX (VTK slice bookkeeping at `mesh.f90:159,323-329`), so any inventory line numbers generated against `ce1f659` are off by 0–8 lines after `mesh.f90:158`.
- Pressure formulation and the iteration loop in detail: `docs/pressure/00-fds-pressure-baseline.md`. §4 below only covers the AMReX mapping.
- Driver architecture decision: `docs/adr/ADR-001-driver-architecture.md`. The Fortran-interface feature table and the recommendation are in `docs/amrex/driver-options.md`.
- V&V case statistics quoted below come from `docs/vv/verification_case_survey.csv` (928 cases).

Fit scale: **good** means a direct 1:1 AMReX object exists. **partial** means AMReX has the mechanism but FDS semantics need adaptation code or an algorithm change. **poor** means no AMReX equivalent, or the FDS design is structurally at odds with AMReX.

## Summary table

| # | FDS concept | AMReX counterpart | Fit |
|---|---|---|---|
| 1 | `MESH_TYPE` / `MESHES(NM)` | `Box` + per-level `BoxArray`/`DistributionMapping`/`Geometry`; field arrays → `MultiFab` components; per-mesh lists → `LayoutData`/per-box containers | partial (fields good; unstructured per-mesh lists and stretched grids poor) |
| 2 | Ghost cells, `OMESH`, `MESH_EXCHANGE(CODE)` | `FabArray::FillBoundary`, `FillPatchSingleLevel/TwoLevels`, `ParallelCopy`, `SumBoundary` | good for same-level halos; partial for FDS's coarse/fine interpolated boundaries |
| 3 | Multi-mesh with different resolutions | `AmrCore` levels, `ref_ratio`, `FluxRegister`/`YAFluxRegister`, `average_down` | partial (FDS is non-nested multi-block, not nested AMR) |
| 4 | Pressure: per-mesh FFT + iterations, ULMAT, UGLMAT/HYPRE | `MLMG` + `MLPoisson` (cell-centred), `MLABecLaplacian`, `MLNodeLaplacian`, HYPRE bottom/`HypreMLABecLap`, `FFT::Poisson[Hybrid]` | good (numerically better), but changes the FDS algorithm |
| 5 | OBST/WALL cells, GEOM/CC_IBM | EB2 (`EB2::Build`, `IF_Box`, STL `IndexSpace_STL`), `EBFArrayBoxFactory`, `EBFluxRegister`, state redistribution | poor→partial (OBST = staircase, EB = cut cell; GEOM ≈ EB conceptually but data models differ) |
| 6 | Lagrangian particles (`part.f90`) | `ParticleContainer` (SoA + runtime comps), `Redistribute`, `ParticleToMesh`/`MeshToParticle` | partial (attached 1-D solid storage has no counterpart) |
| 7 | Staggered U,V,W | face-centred `MultiFab` with nodal `IndexType` in one direction; face FillPatch (`face_linear`/`face_divfree`) | good (index offset needs care) |
| 8 | Predictor/corrector, global `DT` | `AmrCore` without subcycling (incflo/PeleLMeX pattern). `Amr`/`AmrLevel` subcycling optional | good if no subcycling. Poor with subcycling (pressure/divergence constraint) |
| 9 | Smokeview output, devices | `WriteMultiLevelPlotfile`, particle `Checkpoint`/`WritePlotFile`, point sampling | poor for Smokeview compatibility; good for AMReX-native I/O |

---

## 1. `MESH_TYPE` / `MESHES(NM)` vs `Box`/`BoxArray`/`DistributionMapping`/`Geometry`/`MultiFab`/`FArrayBox`

**FDS evidence**
- `TYPE MESH_TYPE` spans `mesh.f90:16-354`. The global array `TYPE (MESH_TYPE), SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: MESHES` is at `mesh.f90:356`. `POINT_TO_MESH(NM)` (`mesh.f90:505-916`) aliases module pointers to one mesh (`M=>MESHES(NM)`, `mesh.f90:512`).
- A mesh bundles four different kinds of data:
  - **Structured fields**: `U,V,W,US,...,RHO,TMP,H,HS,D,DS,MU,...` (`mesh.f90:18-62`), `ZZ(:,:,:,N)` (`mesh.f90:64-85`), integer maps `CELL_INDEX`, `PRESSURE_ZONE`, `MUNKH` (`mesh.f90:107-109,225`).
  - **1-D grid metrics**: `X,Y,Z,XC,DX,RDX,RDXN,HX,...` (`mesh.f90:185-214`), with cell counts `IBAR/JBAR/KBAR` (`mesh.f90:147-149`). They allow non-uniform (stretched) spacing via `TRNX/TRNY/TRNZ` (`TRNX_ID...` `mesh.f90:141`).
  - **Unstructured per-mesh lists**: `WALL`, `EXTERNAL_WALL`, `THIN_WALL`, `OBSTRUCTION`, `VENTS`, `CELL`, `CFACE`, `CUT_CELL`, `LAGRANGIAN_PARTICLE`, `BOUNDARY_*` storage (`mesh.f90:219-319`).
  - **Neighbour/exchange/output state**: `OMESH(NMESHES)`, `NEIGHBORING_MESH`, `SLICE`, `PATCH` (`mesh.f90:156-157,301,322,335-336`), plus FireX VTK slice bookkeeping (`N_SLCF_VTK`, `UNIQUE_SLICE_NAMES`, ..., `mesh.f90:323-329`).
- Per the inventory, of the 431 `MESH_TYPE` fields, 117 are rank-3/4 arrays. Categories: 102 wall/obst/cut-cell, 71 grid geometry, 50 pressure solver, 35 scratch, 27 face-centred, 25 cell-centred solution, 12 particle (`docs/inventory/mesh_fields.csv`).
- Allocation is per mesh with local 0-based indexing and 1–2 ghost layers: `RHO/TMP(-1:IBP1+1,...)` (`init.f90:524-525`), `ZZ(-1:IBP1+1,...,N_TOTAL_SCALARS)` (`init.f90:527`), `U(-1:IBP1,0:JBP1,0:KBP1)` (`init.f90:533`), `H/KRES/D/DS/MU(0:IBP1,...)` (`init.f90:556-562`).
- MPI decomposition assigns whole meshes to ranks: `DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX` (e.g. `main.f90:205`) and `SNODE = PROCESS(NM)` (`main.f90:3142`). The mesh count and layout are fixed at input time (`&MESH` lines).

**AMReX equivalent (evidence)**
- `Box` (`Base/AMReX_Box.H`) is an index-space rectangle in a *global* index space per level. `BoxArray` (`Base/AMReX_BoxArray.H:680`) is the union of boxes on a level. `DistributionMapping` (`Base/AMReX_DistributionMapping.H:50`) maps box → rank. `Geometry` (`Base/AMReX_Geometry.H:82`) holds the physical domain, periodicity and a **single uniform `dx` per direction per level** (`Base/AMReX_CoordSys.H:79-94`, `CellSize()`). It supports `cartesian`, `RZ` and `SPHERICAL` coordinates (`Base/AMReX_CoordSys.H:35`).
- `MultiFab` (`Base/AMReX_MultiFab.H:37`) is a distributed collection of `FArrayBox` (`Base/AMReX_FArrayBox.H:234`), one per box, with `ncomp` components and `ngrow` ghosts. `iMultiFab` (`Base/AMReX_iMultiFab.H`) is the integer version. Arbitrary per-box objects go in `LayoutData<T>` (`Base/AMReX_LayoutData.H`).
- Levels are owned by `AmrMesh`/`AmrCore` (`AmrCore/AMReX_AmrMesh.H:82`, `AmrCore/AMReX_AmrCore.H:29`).

**Proposed mapping**
- One FDS 3-D field `M%X(:,:,:)` becomes one component of a level `MultiFab` (cell-centred, or face-centred for staggered quantities, §7). `M%ZZ(:,:,:,1:N)` becomes one `MultiFab` with `N` components. Keep FDS's ghost widths: `ngrow=2` for `RHO,TMP,ZZ`, `ngrow=1` for `H,D,MU`.
- An FDS mesh becomes an FDS *input region*, not a runtime object. `&MESH` lines only define the level-0 domain and, optionally, static refined regions (see §3). At runtime the unit of work is an AMReX box/tile. `NM` loops become `MFIter` loops.
- Integer maps (`CELL_INDEX`, `PRESSURE_ZONE`, solid masks) become `iMultiFab`. Unstructured lists (`WALL`, `CFACE`, `CUT_CELL`, `OBSTRUCTION` fragments, particles) become per-box containers (`LayoutData<std::vector<...>>`, or a Fortran array indexed by local box index), rebuilt on every regrid.
- Grid metrics (`DX(I)`, `RDXN(I)`...) collapse to scalars `dx(lev)` on uniform AMReX levels. Kernels that index `DX(I)` either receive constant arrays of the right extent (cheapest migration) or are rewritten to scalars.
- Index convention: AMReX `dataPtr` returns Fortran pointers with *global* bounds (`F_Interfaces/Base/AMReX_multifab_mod.F90:758-771`, `dp(bx%lo(1):,...)`). FDS kernels expect `1..IBAR` with a 0/-1 ghost. There are two options: (a) pass `lo,hi` and loop globally, or (b) re-point with Fortran lower-bound remapping (`RHO(-1:,-1:,-1:) => dp` for ng=2 cells, `U(-1:,0:,0:) => dp` for ng=1 x-faces). That keeps kernel bodies unchanged but makes `I` box-local: cell `I` ↔ global `lo+I-1`, face `I` ↔ global face `lo+I` (§7; `driver-options.md` §4).

**Fit: partial.** Structured fields map cleanly. The fit is poor for (i) stretched grids, (ii) the large body of per-mesh unstructured state that FDS builds once at initialization, and (iii) global `MESHES(NOM)` access to other meshes.

**Risks**
- R1.1 **Stretched grids** (`NAMELIST /TRNX/ /TRNY/ /TRNZ/`, `read.f90:1000-1002`): AMReX levels are uniform. 17/928 V&V cases use stretching (`verification_case_survey.csv`, column `stretched`). Mitigations: drop stretching and replace it with static refinement, or map it through a coordinate transform (ERF-style terrain metrics are a precedent, but that is a large effort). `FFT::PoissonHybrid` does allow a non-uniform *last* dimension (`FFT/AMReX_FFT_Poisson.H:158-160,227-229`), but nothing else in AMReX does.
- R1.2 **Direct cross-mesh access**: `MESHES(NOM)%...` appears 381 times per ADR-001's count, e.g. `M2%DY(JJO)` in `MATCH_VELOCITY` (`velo.f90:2691`) and `MESHES(NOM)%DX(EWC%IIO_MIN)` in `PRESSURE_SOLVER_COMPUTE_RHS` (`pres.f90:123`). None of these exist in an AMReX box model, where neighbour data is only reachable through ghost cells. Each one must be replaced by ghost-cell reads or a precomputed per-box table.
- R1.3 **Regrid invalidates per-mesh setup**: WALL/EXTERNAL_WALL/OMESH/CUT_CELL tables are built at init (`INITIALIZE_MESH_EXCHANGE_1` `main.f90:2059`, `INIT_WALL_CELL` `init.f90:2975-3386`, geometry setup in `geom.f90`). With dynamic AMR they must be rebuilt inside `MakeNewLevel*/RemakeLevel`, and their *state* (wall temperatures, 1-D solid profiles in `BOUNDARY_ONE_D`) must be interpolated or transferred. This is the single largest refactor item outside the pressure solver.
- R1.4 **Box granularity**: AMReX will cut level-0 into many boxes (`max_grid_size` default 32/64 in 3-D, `blocking_factor` 8: `AmrCore/AMReX_AmrMesh.H:30-40`). FDS per-mesh overheads (wall-cell loops, FFT set-up, output files per mesh) scale with the number of boxes.
- R1.5 **2-D / cylindrical**: FDS `TWO_D`/`CYLINDRICAL` are runtime flags (`pres.f90:321`). AMReX dimension is compile time (`AMReX_SPACEDIM`, `Tools/CMake/AMReXOptions.cmake:25`). 287/928 V&V cases are 2-D. Either run them as 3-D with one cell in y, or build a 2-D variant. RZ is supported by `CoordSys::RZ`.

## 2. Ghost cells and `MESH_EXCHANGE` vs `FillBoundary` / `FillPatch`

**FDS evidence**
- Neighbour discovery happens once: `INITIALIZE_MESH_EXCHANGE_1(NM)` (`main.f90:2059-2301`) allocates `OMESH(NMESHES)` (`main.f90:2077`). It walks `EXTERNAL_WALL` cells whose `NOM` is the other mesh to build receive index lists `IIO_R/JJO_R/KKO_R/IOR_R` (`main.f90:2143-2164`) and area ratios for non-matching resolutions (`main.f90:2131-2138`). Embedded meshes are detected geometrically (`main.f90:2174-2180`). `OMESH_TYPE` is at `type.f90:1033-1078`.
- The exchange is `MESH_EXCHANGE(CODE)` (`main.f90:3117-3975`) plus `POST_RECEIVES(CODE)` (`main.f90:2945-3111`). The codes are semantic, not just per-field. For example:
  - 1/4: `RHO,ZZ` after predictor/corrector (`main.f90:3301`)
  - 3/6: `H/HS` and velocities (`main.f90:3423`)
  - 5: `FVX/FVY/FVZ` face fluxes during pressure iterations (`main.f90:3352`)
  - 2: radiation (`main.f90:3486`)
  - 11: particles (`main.f90:3525`), with 7 as the same-rank variant (`main.f90:3514`)
  - 20: particle drag (`main.f90:3399`)
  - 8/9/10/19: wall data (`main.f90:3228-3293`)
  - 15–18: OBST mass (`main.f90:3587-3610,3952`)
  - 14: level set/terrain (`main.f90:3547`; FireX calls it under `TERRAIN_CASE`, `main.f90:794`)
  - CC_IBM adds its own `MESH_CC_EXCHANGE(CODE)` (`main.f90:3135`).
- Received data lands in **OMESH copies** of neighbour arrays (e.g. `OM%US`, `OM%U` used in `MATCH_VELOCITY`, `velo.f90:2686-2707`), *not* directly in the ghost cells. Ghost values are then set per external wall cell by physics code, e.g. the interpolated-boundary branch of the wall routines (`wall.f90:353-375`, `SECOND_ORDER_INTERPOLATED_BOUNDARY`) and pressure BCs (`pres.f90:117-130`).
- The time loop calls exchanges at fixed points: `main.f90:790,794,801,803,909,959,991,1000,1002,1014,1040,1074-1079,1105,1112`.

**AMReX equivalent (evidence)**
- Same-level halo fill with periodicity: `FabArray::FillBoundary` (`Base/AMReX_FabArray.H:1392-1399`). Non-conforming copies: `ParallelCopy`. Accumulating ghost contributions (particle deposition, drag): `SumBoundary` (exposed to Fortran too, `F_Interfaces/Base/AMReX_multifab_mod.F90:57-61`).
- Coarse/fine + physical BC ghost fill with time interpolation: `FillPatchSingleLevel` / `FillPatchTwoLevels` (`AmrCore/AMReX_FillPatchUtil.H:95,126,174,228,286`). There is a face-array overload at `:286` and a reusable `FillPatcher` (`AmrCore/AMReX_FillPatcher.H`). Interpolaters are in `AmrCore/AMReX_Interpolater.H` and `AMReX_MFInterpolater.H`.

**Proposed mapping**
- Codes 1, 3, 4, 6, 5 and 14 (pure field halos) become `FillPatch` of the relevant `MultiFab` group (`FillBoundary` on a single level). FDS's two-step "OMESH copy then wall-cell ghost assignment" becomes one step: FillPatch fills ghosts directly, and FDS's interpolated-boundary logic in `wall.f90`/`velo.f90` is removed for same-level neighbours.
- Code 11 (particles) → `ParticleContainer::Redistribute`. Code 20 (drag) → deposit then `SumBoundary`. Code 2 (radiation) → FillBoundary of intensity/`UIID` per angle band. Codes 8/9/19 (wall back-side data across meshes) and 15–18 (OBST mass across meshes) have **no AMReX analogue**. They exchange unstructured per-object data and need custom MPI (or `ParallelDescriptor` helpers) keyed by global OBST/wall IDs.
- The physical BC fill (`amrex_physbc_proc` / `PhysBCFunct`) must reproduce FDS's external-boundary ghost logic for OPEN/MIRROR/SOLID vents, currently spread across `wall.f90` and `velo.f90:VELOCITY_BC` (`velo.f90:1799`).

**Fit: good** for same-level halos, where AMReX removes roughly 800 lines of pack/unpack in `main.f90`. **Partial** for coarse/fine interfaces (see §3) and for unstructured exchanges (wall/OBST).

**Risks**
- R2.1 FDS ghost values are **physics-dependent** (e.g. `UVW_SAVE` and the `BOUNDARY_TYPE_PREVIOUS` bookkeeping in `mass.f90:424-432,598`). A generic FillPatch will not reproduce FDS bit-for-bit, so V&V baselines will shift at interfaces.
- R2.2 Exchange placement is tangled with algorithm steps (e.g. `MESH_EXCHANGE(5)` inside `PRESSURE_ITERATION_SCHEME`, `main.f90:1617,1636,1663,1687`). Each call site needs an explicit "which MultiFabs, how many ghosts, what time" decision.
- R2.3 `FPhysBC` ignores the `nghost` and `bccomp` arguments and calls the user routine for the whole `MultiFab` (`F_Interfaces/Base/AMReX_FPhysBC.cpp`). So a Fortran physical-BC routine must handle every ghost region itself.

## 3. FDS multi-mesh with different resolutions vs AmrCore levels, refinement ratios, flux registers

**FDS evidence**
- FDS supports abutting and embedded meshes of different resolution, but as **peer blocks at a single "level"**, all advanced with the same `DT` (`DT = MINVAL(DT_NEW)`, `main.f90:715,738-741`).
- Alignment is enforced: `ERROR(431): MESH ... is out of alignment` (`init.f90:3185-3222`). Cell faces must coincide within `ALIGNMENT_TOLERANCE`, which in practice means integer ratios.
- The coarse/fine coupling is **FDS-specific conservative matching**:
  - `MATCH_VELOCITY` (`velo.f90:2618-2860`) forces normal face velocities at `INTERPOLATED_BOUNDARY` wall cells to agree using area-weighted sums over the other mesh's faces (`DA_OTHER`, `velo.f90:2682-2714`).
  - `MATCH_VELOCITY_FLUX` (`velo.f90:2863`) does the same for momentum fluxes `FVX..`.
  - The mass transport then uses those matched velocities (`UVW_SAVE`, `mass.f90:424-432`).
  - The Poisson BC at an interpolated boundary is a distance-weighted Dirichlet average of the two sides' `H` (`pres.f90:115-130`).
  - Obstructions inside a finer mesh are hidden in the coarse one (`read.f90:11459`).
- Usage: only 12/928 V&V cases have more than one resolution; 181 are multi-mesh (`verification_case_survey.csv`).

**AMReX equivalent (evidence)**
- `AmrCore` (`AmrCore/AMReX_AmrCore.H:29`) manages properly nested levels. `ref_ratio` is a `Vector<IntVect>` (anisotropic allowed in C++, `AmrCore/AMReX_AmrMesh.H:30,156-162`). The Fortran wrapper exposes only a scalar per level (`amrex_fi_get_ref_ratio` uses `MaxRefRatio`, `F_Interfaces/AmrCore/AMReX_amrcore_fi.cpp`).
- Tagging uses `ErrorEst` / `TagBoxArray` (`AmrCore/AMReX_TagBox.H`, `AMReX_ErrorList.H`). Static refinement can be expressed with tagging-by-region or by setting grids directly.
- Refluxing: `FluxRegister` (`AmrCore/AMReX_FluxRegister.H:22`: `CrseInit`/`FineAdd`/`Reflux`), `YAFluxRegister` (`Boundary/AMReX_YAFluxRegister.H:32`), `EBFluxRegister` (`EB/AMReX_EBFluxRegister.H:62`). Averaging down: `average_down`, `average_down_faces` (`Base/AMReX_MultiFabUtil.H`).

**Proposed mapping**
- **Phase 1 (static refinement)**: express each FDS fine mesh as a level-1 (or level-2) static region with `ref_ratio=2` or `4`. FDS's area-weighted `MATCH_VELOCITY` corresponds to `average_down_faces` of the fine normal velocity onto coarse faces at the C/F boundary, followed by a composite projection. FDS's scalar consistency via matched velocities corresponds to `FluxRegister` refluxing of `rho*Z` and `rho` fluxes (FDS face fluxes `FX/FY/FZ` and `ADV_F*`/`DIF_F*`, `mesh.f90:68-85`).
- **Phase 2 (dynamic AMR)**: tag on HRR (`Q`), temperature, vorticity, and species gradients. Regrid every N steps.
- Non-nested FDS layouts (e.g. a fine mesh that spans a coarse mesh edge, or side-by-side meshes of different `dx` whose union is not a nested hierarchy) must be converted. The coarse level must cover the whole domain, and fine patches must be nested with `blocking_factor` granularity. **This is an input-compatibility break.**

**Fit: partial.** The mechanisms (nesting, interpolation, refluxing) are standard and better defined than FDS's. But FDS semantics (peer meshes, arbitrary integer ratios per face, anisotropic ratios such as refining only z) do not map 1:1.

**Risks**
- R3.1 FDS allows anisotropic ratios (e.g. finer `DZ` only; see the checks at `init.f90:3212-3218`). AMReX C++ supports `IntVect` ratios; the Fortran interface does not.
- R3.2 Conservation on EB + AMR needs `EBFluxRegister` plus redistribution. This is C++ only (see driver-options.md).
- R3.3 `n_error_buf`, `blocking_factor` and `grid_eff` will produce different fine-patch footprints than hand-placed FDS meshes, which is a V&V reproducibility risk.

## 4. Pressure solver: FFT/ULMAT/UGLMAT(+HYPRE) with pressure iterations vs MLMG

(Algorithmic detail lives in `docs/pressure/00-fds-pressure-baseline.md`; this section only covers the AMReX mapping.)

**FDS evidence (FireX tree)**
- FDS solves a **constant-coefficient** Poisson equation for `H = p̃/ρ + |u|²/2`: `PRHS = TRM1+TRM2+TRM3+TRM4` with `TRM4 = -DDDT` (`pres.f90:256-257`). The baroclinic term goes to the RHS (`BAROCLINIC_CORRECTION`, `velo.f90:3216`) and is iterated.
- **Default FFT path**: `PRESSURE_SOLVER_FFT(NM)` (`pres.f90:318`) calls FISHPACK-derived `H3CZSS/H2CZSS/H2CYSS/H3CSSS` (`pres.f90:320`, `pois.f90:187` `H3CZSS`). It is one direct solve **per mesh**, and stretching is allowed in up to two directions (`IPS` transposes, `init.f90:2356+`). At an `INTERPOLATED_BOUNDARY` the BC is a Dirichlet value averaged with the neighbour's `H` from the previous iterate (`pres.f90:115-130`). The result is a block-Jacobi/Schwarz iteration between meshes.
- **Pressure iterations**: `PRESSURE_ITERATION_SCHEME` (`main.f90:1601-1745`). Each iteration exchanges `FVX..` and `H` (`MESH_EXCHANGE(5)`, `main.f90:1617,1636,1663,1687`), selects FFT/GLMAT/UGLMAT/ULMAT (`main.f90:1659-1668`), and iterates until `VELOCITY_TOLERANCE`/`PRESSURE_TOLERANCE` or `MAX_PRESSURE_ITERATIONS` (`cons.f90:547-559`). The velocity error comes from interpolated mesh interfaces *and* from the direct-forcing treatment of OBSTs in the FFT solve (`COMPUTE_VELOCITY_ERROR`, `pres.f90:802`).
- **ULMAT** (`MODULE LOCMAT_SOLVER`, `pres.f90:1084`; `ULMAT_SOLVER_SETUP` `pres.f90:1131`; `ULMAT_SOLVE_ZONE` `pres.f90:1463`): a per-mesh unstructured matrix over gas cells only. Solved with MKL PARDISO or HYPRE PCG+BoomerAMG on `MPI_COMM_SELF` (`pres.f90:2962,3028,1760`).
- **GLMAT/UGLMAT** (`MODULE GLOBMAT_SOLVER`, `pres.f90:3129`; `GLMAT_SOLVER` `pres.f90:3299`; `GLMAT_SOLVER_SETUP` `pres.f90:3662`): one global matrix across meshes, per pressure zone (`ZSL_COMM(IPZ)%COMM`, `pres.f90:4736`). "U" means gas cells only, i.e. obstructions are excluded exactly. Solver choice via `DEFINE_PRES_METHOD` (`func.f90:7143-7219`). V&V usage: 87/928 cases select a matrix solver (40 UGLMAT, 24 ULMAT, 14 ULMAT HYPRE, ...; the rest are FFT).
- **FireX-specific (new in this base)**: a GPU "resource set" path.
  - `DEFINE_RS_COMM_INFO` (`main.f90:5096-5146`) splits `MPI_COMM_WORLD` into `MPI_COMM_RS` groups of `FDS_RANKS_PER_GPU` ranks (env var, `main.f90:5108-5124`).
  - In `GLMAT_SOLVER` all ranks of a resource set `MPI_GATHERV` their RHS to the RS master (`pres.f90:3411-3419`). Only masters solve (`MASTER_IF`, `pres.f90:3422`). HYPRE vectors are migrated to device memory (`HYPRE_IJVECTORMIGRATE(...,HYPRE_MEMORY_DEVICE)`, `pres.f90:3453-3458`), PCG runs, the solution migrates back (`pres.f90:3469-3473`) and is `MPI_SCATTERV`ed (`pres.f90:3485-3492`).
  - Device policy is selected with `HYPRE_SETMEMORYLOCATION/HYPRE_SETEXECUTIONPOLICY` (`pres.f90:1177-1184,4431-4438`) under `WITH_HYPRE_DEVICE`, toggled at runtime by `HYPRE_DEVICE_RUN` (`cons.f90:569`).
  - CMake options `USE_HYPRE_NVIDIA/AMDGPU/INTELGPU` enable HYPRE CUDA/HIP/SYCL and define `WITH_HYPRE_DEVICE` (`fds-amr/src/CMakeLists.txt:16-18,154-208`). Resource-set bookkeeping was added to `ZONE_SOLVE_TYPE` (`NUNKH_LOCAL_RS`, `TOT_NNZ_H_RS`, `NUNKH_LOC_RS`, `UNKH_IND_RS`; `type.f90:1813-1818`).
  - Net effect: FireX already runs *only the linear solve* on the GPU while all FDS kernels stay host Fortran. That is exactly the "host Fortran kernels + GPU linear solver" split AMReX would give us (see driver-options.md §GPU).

**AMReX equivalent (evidence)**
- `MLMG` (`LinearSolvers/MLMG/AMReX_MLMG.H:26`) is a composite multi-level geometric multigrid. Cell-centred operators:
  - `MLPoisson` (`AMReX_MLPoisson.H:30`): constant coefficient; FDS's `∇²H` maps directly here.
  - `MLABecLaplacian` (`AMReX_MLABecLaplacian.H:20`): `αa − β∇·(b∇)` with face `b`; needed if we move to `∇·(1/ρ ∇p)` form.
  - EB variants `MLEBABecLap` (`AMReX_MLEBABecLap.H`, `setEBDirichlet` `:149-187`).
- Nodal operator: `MLNodeLaplacian` (`AMReX_MLNodeLaplacian.H:31`). It is used for approximate nodal projections by incflo (`ref/incflo/src/projection/incflo_apply_nodal_projection.cpp:184-188`) and PeleLMeX (`ref/PeleLMeX/Source/PeleLMeX_Projection.cpp:771-787`). Their `NodalProjector`/`MacProjector` live in AMReX-Hydro, **not** in `(local AMReX checkout)`.
- Domain BCs: `LinOpBCType` includes Dirichlet, Neumann, inhomogNeumann, Robin, Periodic (`Boundary/AMReX_LO_BCTYPES.H:27-39`). Robin coefficients come through `setLevelBC(amrlev, levelbcdata, robin_a, robin_b, robin_f)` (`AMReX_MLLinOp.H:286-297`).
- Bottom solvers: `BottomSolver::{smoother,bicgstab,cg,bicgcg,cgbicg,hypre,petsc,custom}` (`AMReX_MLLinOp.H:40-42`). The HYPRE interfaces are in `Extern/HYPRE` (`AMReX_HypreABecLap3.H`, `AMReX_HypreMLABecLap.H:34`, `AMReX_HypreNodeLap.H`, `AMReX_HypreIJIface.H`). The build option `AMReX_HYPRE` defaults to OFF (`Tools/CMake/AMReXOptions.cmake:356`).
- Global FFT solvers: `FFT::Poisson` (periodic/Dirichlet/Neumann, `FFT/AMReX_FFT_Poisson.H:27-80`) and `FFT::PoissonHybrid` (non-uniform spacing in the **last** dimension, `:158-237`). They are single-level, whole-domain and C++ only. Build option `AMReX_FFT` defaults to OFF (`AMReX_FFT`, `Tools/CMake/AMReXOptions.cmake:306`).

**Proposed mapping (consistent with the team plan: one composite MLMG solve, global dt, no subcycling in Phase 1)**
- Replace `PRESSURE_ITERATION_SCHEME` + `PRESSURE_SOLVER_*` wholesale with **one composite cell-centred MLMG solve per predictor/corrector stage** over all levels. This is a MAC-type projection on FDS's staggered grid, so MLPoisson/MLABecLaplacian is the correct family. MLNodeLaplacian would force a collocated/nodal redesign and is *not* recommended.
  - Keep FDS's H-form with `MLPoisson` for Phase 1 (smallest change to `VELOCITY_PREDICTOR/CORRECTOR`, `velo.f90:1572,1685`).
  - Optionally move later to `MLABecLaplacian` with `b = 1/ρ` on faces (the variable-density form used by incflo/PeleLMeX), which removes the baroclinic iteration.
- Interfaces between mesh blocks disappear (a composite solve has no inter-mesh Dirichlet lag), so the "velocity error at mesh interfaces" motivation for iterations goes away.
- OBSTs: either (a) EB with covered cells (MLEB*; C++ only), or (b) no EB, with solid faces given `b=0` in MLABecLaplacian so that the operator is exactly UGLMAT-like, gas-only. Option (b) reproduces UGLMAT semantics without EB and **works from the Fortran interface**: `amrex_abeclaplacian` exposes `set_scalars/set_acoeffs/set_bcoeffs` (`F_Interfaces/LinearSolvers/AMReX_abeclaplacian_mod.F90:14-16`) with face β MultiFabs, so a variable-coefficient β = 1/ρ projection also works from Fortran. Only EB operators and nodal projection force C++ (see `driver-options.md` §1.1). Velocity iterations for OBSTs then become unnecessary. Pressure zones become separate solves or a masked single solve; `CHECK_UNSUPPORTED_MESH` (`pres.f90:3907`) documents the Dirichlet/Neumann-disconnected-domain issue that remains.
- Open boundaries → Dirichlet on H (FDS `BXS = P_EXTERNAL/RHO + KRES`, `pres.f90:188-190` region); walls → Neumann (`BXS = HX(0)*(-FVX + DUNDT)`, `pres.f90:81`) as inhomogeneous Neumann through `setLevelBC`.
- HYPRE: use MLMG with `BottomSolver::hypre` for the coarsest level, or `HypreMLABecLap`. The FireX GPU resource-set gather/scatter is not needed: AMReX MLMG runs natively on GPU (C++ build), and its HYPRE bottom solve uses device memory when AMReX is built with a GPU backend.

**Fit: good** for the numerics; MLMG is exactly what FDS's UGLMAT approximates. It is a deliberate **algorithm change**: results at mesh interfaces and near OBSTs will differ from FFT+iterations, and V&V must re-baseline.

**Risks**
- R4.1 **Stretched-grid FFT cases** have no MLMG equivalent (uniform levels). `FFT::PoissonHybrid` only covers z-stretching on a single level.
- R4.2 **Performance** of MLMG vs per-mesh FISHPACK FFT for the 90% of cases that use FFT. MLMG is iterative (typically 5–15 V-cycles to 1e-10 relative). The FFT solve is direct. We need a benchmark (driver-options.md Prototype P2).
- R4.3 **Pressure zones / `PBAR` / `D_PBAR_DT`** (`mesh.f90:94-109`), plus HVAC coupling and leakage, are global per-zone integrals. MLMG solves one linear system; zone solvability constraints need per-zone RHS mean removal as in `GLMAT_SOLVER` (`pres.f90:3395-3406`).
- R4.4 **Tunnel preconditioner** (`TUNNEL_POISSON_SOLVER`, `pres.f90:505`) has no counterpart. MLMG agglomeration/consolidation probably makes it unnecessary, but this is untested.
- R4.5 The FireX HYPRE-GPU investment (`pres.f90` +582 lines) is **superseded**, not reused, by the MLMG route. Only its build knowledge (HYPRE 2.32+, `HYPRE::HYPRE` target, `CMakeLists.txt:196-208`) carries over.

## 5. OBST/WALL cells and GEOM complex geometry vs AMReX EB

**FDS evidence**
- **OBST** blocks are snapped to cell faces. `CELL_TYPE%SOLID/OBST_INDEX/WALL_INDEX(-3:3)` (`type.f90:2175-2187`) marks solid cells, and `CELL_INDEX(I,J,K)` (`mesh.f90:225`) indexes into `CELL`. Every gas/solid face carries a `WALL_TYPE` (`type.f90:434-457`) pointing into `BOUNDARY_COORD/ONE_D/PROP1/PROP2/RADIA` storage, which holds the 1-D solid heat conduction and pyrolysis state per wall cell. `THIN_WALL` covers zero-thickness obstructions (`mesh.f90:298`).
- OBSTs can be **consumable/removable** (`type.f90:1140-1141`). `CREATE_OR_REMOVE_OBSTRUCTIONS` (`main.f90:1782-`) calls `OPEN_AND_CLOSE` (`init.f90:4507`), a global `MPI_ALLREDUCE` (`main.f90:1788`), `REASSIGN_WALL_CELLS` (`init.f90:4898-5228`) and matrix re-setup (`GLOBAL_MATRIX_REASSIGN`, `main.f90:1806`). So geometry changes at runtime.
- **GEOM / CC_IBM** (`CC_IBM`, `cons.f90:732`): triangulated surfaces are intersected with the Cartesian grid in `MODULE COMPLEX_GEOMETRY` (`geom.f90:3`, 27.7k lines) and handled numerically in `MODULE CC_SCALARS` (`ccib.f90:7`, 24k lines). Cut cells and faces are `CC_CUTCELL_TYPE`/`CC_CUTFACE_TYPE` (`type.f90:1391-1453,1279`). Boundary faces are `CFACE_TYPE` (`type.f90:1361`). Small cut cells are **linked** to neighbours (`IJK_LINK`, `type.f90:1400`; `FINEST_LINK_LEV`, `mesh.f90:242`). There is exactly one cut cell per Cartesian cell (`NCELL` "Now fixed at 1 by blocking", `type.f90:1393`). CC_IBM has its own exchange (`MESH_CC_EXCHANGE`, `main.f90:3135`) and prefers UGLMAT.

**AMReX equivalent (evidence)**
- EB2 builds implicit-function or STL geometry once per level hierarchy: `EB2::Build` (`EB/AMReX_EB2.H:167-260`), `IF_Box`/`IF_Union` (`EB/AMReX_EB2_IF_Box.H`, `AMReX_EB2_IF_Union.H`), STL (`EB/AMReX_EB2_IndexSpace_STL.H`, `STLtools` `EB/AMReX_EB_STL_utils.H:25`). Data: `EBFArrayBoxFactory`, `EBCellFlag` (single-valued cells), volume/area fractions and centroids. Conservation across C/F: `EBFluxRegister` (`EB/AMReX_EBFluxRegister.H:62`). Small-cell stability: state redistribution (`EB/AMReX_EB_StateRedistribute.cpp`, `AMReX_EB_Redistribution.H`). EB solvers: `MLEBABecLap`, `MLNodeLaplacian_eb.cpp`. Output: `EB_WriteMultiLevelPlotfile` (`Base/AMReX_PlotFileUtil.H:140`).
- Build flag `AMReX_EB` (default OFF, `Tools/CMake/AMReXOptions.cmake:286`). **No EB bindings in `Src/F_Interfaces`** (no `EB2`/`EBFArrayBox`/`EBCellFlag` symbols; `rg` count 0).

**Proposed mapping**
- **OBST (Phase 1): no EB.** Represent OBSTs as an `iMultiFab` solid mask plus per-direction face masks, rebuilt at regrid from the global OBST list. FDS OBSTs are grid-aligned *on the level they were defined for*; on a coarser level they may cut cells. Rule: OBST extents are snapped to the finest level that covers them, and coarse levels treat partially covered cells by volume fraction (or conservatively "solid if >50%"). Each gas/solid face keeps a per-box `WALL` record (list in `LayoutData`). The 1-D solid state lives with the OBST surface patch, *not* the box, so it survives regrids. This is the key data-model change.
- **GEOM (Phase 2+)**: AMReX EB is the natural target (single-valued cut cells, linking ≈ state redistribution, `CFACE` ≈ EB boundary face with `bcent/bnorm`). But FDS's geometry engine (`geom.f90`) and AMReX's EB2 are two independent intersectors. We either adopt EB2 (STL input) and drop `geom.f90` intersection, or feed FDS geometry into EB2 through a custom `GeometryShop`. Either way `ccib.f90` numerics must be ported to EB data structures. This is the subject of ADR-003.

**Fit: poor→partial.** OBST staircases need no EB but do need a new wall-cell data model. GEOM ≈ EB in concept but is a rewrite in practice. Dynamic (burn-away) geometry is poorly supported by EB, which is built once and costly to rebuild.

**Risks**
- R5.1 **Burn-away/removable OBSTs with EB**: EB2 is built at initialization. Changing it at runtime means rebuilding the index space and factories on all levels. This is why Phase 1 keeps OBST as masks (cheap to change).
- R5.2 **Wall-cell state across regrids**: `BOUNDARY_ONE_D` profiles (e.g. `mesh.f90:302-307`) must be conserved when a wall face moves from level ℓ to ℓ+1 (1 coarse face → r² fine faces). This needs a remap policy (copy/average) and energy conservation checks.
- R5.3 The choice of EB forces the C++ route (no Fortran EB API). See driver-options.md.

## 6. Lagrangian particles (`part.f90`) vs `ParticleContainer`

**FDS evidence** (`part.f90` is unchanged in FireX)
- `LAGRANGIAN_PARTICLE_TYPE` (`type.f90:389-429`): 14 integer indices (including `WALL_INDEX`, `CFACE_INDEX`, `DUCT_INDEX`, `BC_INDEX/OD_INDEX/B1/B2/BR_INDEX` into per-mesh `BOUNDARY_*` storage), 4 logicals and 17 reals. Particles are stored per mesh in `LAGRANGIAN_PARTICLE(:)` with count `NLP` (`mesh.f90:302,319`). Particles can carry a **variable-length 1-D interior heat-conduction state** through `OD_INDEX`.
- Routines:
  - insertion `INSERT_ALL_PARTICLES` (`part.f90:118`)
  - `MOVE_PARTICLES` (`part.f90:1824`) / `MOVE_IN_GAS` (`part.f90:2594`)
  - two-way coupling `PARTICLE_MASS_ENERGY_TRANSFER` (`part.f90:3485`) and `PARTICLE_MOMENTUM_TRANSFER` (`part.f90:4626`), which deposit into `FVX_D` etc.
  - migration via `REMOVE_PARTICLES` → `ADD_TO_PARTICLE_SEND_BUFFER` (`part.f90:4778,4854`), packed by `PACK_PARTICLE` (`func.f90:4399-4467`) and exchanged by `MESH_EXCHANGE(11)`/`(7)`
  - drag halo sum `MESH_EXCHANGE(20)`
- 157/928 V&V cases use particles.

**AMReX equivalent (evidence)**
- `ParticleContainer_impl` (`Particle/AMReX_ParticleContainer.H:148`) supports compile-time AoS/SoA components *and* runtime components (`AddRealComp`/`AddIntComp`, `:1316-1349`). `Redistribute(lev_min,lev_max,nGrow)` (`:514-531`) handles migration and level changes. `Checkpoint` (`:787-815`) handles I/O. `AmrParticleContainer_impl` (`AmrCore/AMReX_AmrParticles.H:294`). Deposition/interpolation helpers are in `Particle/AMReX_ParticleMesh.H`.
- Fortran interface: one fixed type `amrex_particle` = `pos(3), vel(3), id, cpu` (`F_Interfaces/Particle/AMReX_particlecontainer_mod.F90:18-23`), hard-wired `AmrParticleContainer<BL_SPACEDIM,0>` (`AMReX_particlecontainer_fi.cpp:8-13`).

**Proposed mapping**
- Put scalar attributes (the 17 reals and the integer tags) in SoA runtime components.
- Store references to wall/CFACE/duct as **global IDs**, not per-mesh indices. Per-mesh indices become invalid after `Redistribute`/regrid.
- Variable-length interior state (`OD_INDEX`): either a fixed maximum number of runtime comps per particle class (pad), or a side table keyed by particle `id/cpu` that migrates with the particle. The latter needs custom communication.
- Deposition (`FVX_D`, `M_DOT_PPP`, `QR_W`) → deposit then `SumBoundary`, then average-down for AMR.

**Fit: partial.** Transport, migration and I/O map well; the attached solid-phase storage does not. **Poor from Fortran** because the Fortran particle API has no attributes.

**Risks**
- R6.1 **Variable-length particle state**, as above.
- R6.2 Particles stuck to walls (`WALL_INDEX`, `CFACE_INDEX`) need the wall data model from §5 to be regrid-stable.
- R6.3 Level-crossing particles with subcycling. Not an issue in Phase 1 (no subcycling).

## 7. Staggered U,V,W vs face-centred MultiFabs (IndexType)

**FDS evidence**: `U(-1:IBP1,0:JBP1,0:KBP1)`, `V(0:IBP1,-1:JBP1,0:KBP1)`, `W(0:IBP1,0:JBP1,-1:KBP1)` (`init.f90:533-535`, same for `US..`). `U(I,J,K)` is on the **forward** x-face of cell `I` (`mesh.f90:12-14`). Fluxes `FVX..` (`mesh.f90:30-38`) and species face fluxes `FX..`/`ADV_F*`/`DIF_F*` (`mesh.f90:68-85`) are face quantities too.

**AMReX equivalent**: `IndexType` (`Base/AMReX_IndexType.H:35,241`). A MultiFab built on `convert(ba, IntVect::TheDimensionVector(d))` is nodal in direction `d`. From Fortran: `amrex_multifab_build(mf, ba, dm, nc, ng, nodal)` with the `nodal(3)` flag (`F_Interfaces/Base/AMReX_multifab_mod.F90:177-182`). Face FillPatch uses `face_linear_interp`/`face_divfree_interp` (`F_Interfaces/AmrCore/AMReX_interpolater_mod.F90`, ids 8/9) and `FillPatchTwoLevels(Array<MF*,DIM>...)` (`AmrCore/AMReX_FillPatchUtil.H:286`). Averaging: `average_down_faces` (`Base/AMReX_MultiFabUtil.H`).

**Index convention**: AMReX face `i` of a cell-box `[lo,hi]` is the *low* face of cell `i`, so faces run `lo..hi+1`. FDS `U(I)` is the *high* face of cell `I`. With FDS cell `I` ↔ AMReX cell `lo+I-1`, **FDS `U(I)` ↔ AMReX face `lo+I`**: FDS `U(0)` is the low domain face `lo`, `U(IBAR)` is `hi+1`. FDS's extra ghost `U(-1)` and `U(IBP1)` need `ngrow=1` in x on the face MultiFab. `amrex_multifab%dataptr` returns a pointer whose lower bounds are the global indices of the grown FAB box (`Base/AMReX_multifab_mod.F90:758-786`). Remapping only the lower bounds, as in `U(-1:,0:,0:) => dp` for a face MultiFab with ng=1, gives exactly FDS local index `I` ↔ AMReX global face `lo+I`. Likewise `RHO(-1:,-1:,-1:) => dp` (ng=2) gives cell `I` ↔ `lo+I-1`. So the offset is a single constant per direction per box (details in `driver-options.md` §4). It must be written once in the shim and unit-tested.

**Fit: good.** Staggered MAC grids are first-class in AMReX (incflo/IAMR use face-centred `umac`).

**Risks**
- R7.1 An off-by-one in the face mapping is silent and catastrophic. Prototype P1 must include a divergence-free round-trip test.
- R7.2 Owner ambiguity on shared faces between boxes. AMReX FillBoundary treats the face on a box boundary as valid in both boxes, so writes must be consistent (`OverrideSync`/`override_sync` is exposed in Fortran, `AMReX_multifab_mod.F90:60`).

## 8. Predictor/corrector time stepping vs AmrCore (subcycling or not)

**FDS evidence**: `MAIN_LOOP` (`main.f90:695-1213`). There is a single global `DT = MINVAL(DT_NEW)` (`main.f90:715,738-741`) with CFL/VN from `CHECK_STABILITY` (`velo.f90:3028`). The predictor starts at `main.f90:748` (mass/density `767-780`, divergence `840-850`, pressure `855`, `VELOCITY_PREDICTOR` `863`) inside `CHANGE_TIME_STEP_LOOP` (`main.f90:774-897`), which *repeats the predictor with a smaller DT* if stability is violated. The corrector starts at `main.f90:928` (density `948-949`, combustion `973`, particles `987`, radiation `1034`, divergence `1056-1085`, pressure `1090`, `VELOCITY_CORRECTOR` `1095`).

**AMReX equivalent**: `AmrCore` leaves time integration to the application. incflo (`ref/incflo/src/incflo.H:28`) and PeleLMeX (`ref/PeleLMeX/Source/PeleLMeX.H:67`) both derive from `AmrCore` and advance **all levels with one dt, no subcycling**, with one composite projection per stage. `Amr`/`AmrLevel` provide subcycling (`Amr/AMReX_Amr.H:123-146`, `subcyclingMode`, `nCycle`) with sync projections as in IAMR. The Fortran `Advection_F` tutorial shows Fortran subcycling with reflux (`ref/amrex-tutorials/ExampleCodes/FortranInterface/Advection_F/Source/evolve_mod.F90:68-137`). `FillPatcher` supports RK-stage fills (`F_Interfaces/AmrCore/AMReX_fillpatch_mod.F90:375-419`).

**Proposed mapping**: Phase 1 is **no subcycling**. One global dt, and the predictor and corrector each do: FillPatch → per-level kernels over `MFIter` → reflux/average-down → one composite MLMG solve. FDS's `CHANGE_TIME_STEP_LOOP` stays at the top level (restart the predictor on all levels). Subcycling is deferred: for low-Mach variable-density flow it requires synchronization projections (IAMR-style), which is a major research effort.

**Fit: good** (no subcycling). **Poor** with subcycling.

**Risks**
- R8.1 Without subcycling, the finest level's CFL sets dt everywhere. With r=2 per level and 2–3 levels this costs 4–8× in coarse-level work. That is acceptable for Phase 1, but the ADR-002 framing should record it.
- R8.2 Step rejection (`CHANGE_TIME_STEP_LOOP`) requires keeping old-time state on all levels (AMReX apps usually keep `old`/`new` MultiFabs, which also FillPatch needs).

## 9. Smokeview output/devices (and FireX VTK-HDF) vs plotfiles

**FDS evidence (FireX tree)**
- Smokeview: the `.smv` index file is written once by `WRITE_SMOKEVIEW_FILE` (`dump.f90:1766-3019`), with one `GRID` block per mesh (`dump.f90:2499`). Data files are keyed by mesh number: slices `CHID_<NM>_<N>.sf` (`dump.f90:501`), boundary files `.bf` (`dump.f90:524`), particles `CHID_<NM>.prt5` (`dump.f90:552`), and SMOKE3D, isosurfaces and PLOT3D similarly. `DUMP_MESH_OUTPUTS(T,DT,NM,FAKEWRITE)` (`dump.f90:77-287`) drives per-mesh output.
- Devices: `UPDATE_DEVICES_1(T,DT,NM)` (`dump.f90:7772-8330`) samples at points located per mesh; `DUMP_DEVICES` (`dump.f90:11155`) runs on rank 0. Restart is per mesh: `DUMP_RESTART(T,DT,NM)` (`dump.f90:3871-4050`).
- **FireX adds `vtkf.f90`** (`MODULE VTK_FDS_INTERFACE`, `vtkf.f90:3`, 4,470 lines):
  - Writes VTKHDF (HDF5) files through **parallel HDF5/MPI-IO** (`H5PSET_FAPL_MPIO_F(...MPI_COMM_WORLD...)`, `vtkf.f90:596,695,742,776`) with VTK type `UnstructuredGrid` (`vtkf.f90:758`). Covers gas-phase/slice/SMOKE3D/boundary/particle/geometry (`BUILD_VTK_*`, `WRITE_VTKHDF_*`, `INITIALIZE_VTKHDF_*`, `vtkf.f90:37-3793`), plus a ParaView state file (`WRITE_PARAVIEW_STATE_FILE`, `vtkf.f90:3992-4438`).
  - Driven by `DUMP_*_VTKHDF` in dump.f90 (`dump.f90:123-160`, `DUMP_PART_VTKHDF` `4622`, `DUMP_SMOKE3D_VTKHDF` `5175`, `DUMP_BNDF_VTKHDF` `12177`). Enabled by `VTK_HDF` (`cons.f90:283`, default `.TRUE.`) and the `WRITE_VTK`/`DT_VTK` inputs (`read.f90:2347-2410`), compiled only `#ifdef WITH_HDF5` (`vtkf.f90:10-12`).
  - Because HDF5 writes are collective, ranks with fewer meshes make **fake writes** so every rank participates (`DUMP_MESH_OUTPUTS(...,LOWER_MESH_INDEX,.TRUE.)` loop, `main.f90:1150-1160`).
  - Note: `WITH_HDF5` is set only by the legacy makefile (`Build/makefile:135`). The FireX `CMakeLists.txt` lists `vtkf.f90` (`:62`) but never finds HDF5 or defines `WITH_HDF5`, so a CMake build gets the stubs.

**AMReX equivalent (evidence)**
- `WriteMultiLevelPlotfile` (`Base/AMReX_PlotFileUtil.H:87`) writes a native multi-level plotfile readable by VisIt/ParaView/yt/amrvis. `EB_WriteMultiLevelPlotfile` (`:140`) and an HDF5 variant `WriteMultiLevelPlotfileHDF5` exist; the latter requires `AMReX_HDF5` (default OFF, `Tools/CMake/AMReXOptions.cmake:366`; `#ifdef AMREX_USE_HDF5` in `F_Interfaces/Base/AMReX_plotfile_fi.cpp:5,32`). Particles: `Checkpoint`/`WritePlotFile`. Raw MultiFab I/O: `VisMF` (`amrex_multifab_write/read`, `F_Interfaces/Base/AMReX_vismf_fi.cpp`). AMReX has **no VTKHDF writer** (only SENSEI adaptors mention AMR VTK types).

**Proposed mapping**
- Checkpoint/restart → AMReX-native (per-level `VisMF` + BoxArray headers + particle `Checkpoint`), replacing per-mesh `DUMP_RESTART`. Wall/OBST 1-D state needs its own checkpoint section.
- Field output → plotfiles for 3-D fields (replacing PLOT3D/SMOKE3D-as-data). **Smokeview compatibility** requires an adapter: either (a) an AMReX-level → Smokeview writer that emits one "mesh" per AMReX box per output time (Smokeview supports many meshes, but not boxes that change between frames), or (b) Smokeview learns plotfiles. Option (a) breaks when regridding changes box layouts, which argues for (b) or for VTK.
- **FireX VTKHDF is the more promising bridge than Smokeview.** VTKHDF has an `OverlappingAMR` type (per the VTK file-format documentation; not verifiable in these trees) that maps naturally onto AMReX levels/boxes. `vtkf.f90`'s parallel-HDF5 machinery could be retargeted from "one UnstructuredGrid piece per FDS mesh" to "one AMR block per AMReX box". The collective-write/fake-write pattern maps onto AMReX's rank-local box lists. This needs an HDF5-enabled build on both sides (`AMReX_HDF5=ON` if we also use AMReX's HDF5 plotfiles; FDS CMake `WITH_HDF5` must be wired, which is a known gap).
- Devices: sample by `(x,y,z)` → finest level containing the point → owning box. Use a per-step list of devices per local box, rebuilt at regrid, with an `MPI_Reduce` to rank 0. Same for HRR/mass integrals (`DUMP_HRR`), which become composite sums with fine-covered masks.

**Fit: poor** for Smokeview file compatibility. **Good** for AMReX-native plotfile/checkpoint. **Partial** for re-using FireX VTKHDF.

**Risks**
- R9.1 Users depend on Smokeview. Losing it is a product decision for the project owner.
- R9.2 Output volume. AMR changes box layouts, so time-series viewers must handle changing topology (plotfile readers do; Smokeview does not).
- R9.3 Device/HRR integrals must exclude fine-covered coarse cells (`amrex_imultifab_build_owner_mask` and masks exist; `AMReX_multifab_mod.F90:22`).

## 10. Other couplings that any AMR design must carry (brief)

- **Radiation** (`COMPUTE_RADIATION` `radi.f90:3744`, `RADIATION_FVM` `radi.f90:3772`; FireX changed `radi.f90` by 213 lines): FVM discrete-ordinates sweeps across meshes, with intensities exchanged by `MESH_EXCHANGE(2)` (`main.f90:1040,1112`) and angle increments amortized over steps. Sweeping across an AMR hierarchy has no AMReX primitive. The likely plan is to solve radiation on one level (coarse or a designated level) and interpolate `QR`. Poor fit; needs its own design.
- **HVAC** (`hvac.f90`) and **pressure zones**: global networks/integrals. Keep them as-is on rank 0 with reductions over boxes.
- **Level set / terrain** (`MESH_EXCHANGE(14)` under `TERRAIN_CASE`, `main.f90:794`): 2-D fields; `MultiFab` with one z-cell or a separate 2-D BoxArray.

## Top mapping risks (ranked)

1. **Wall/OBST data model across regrids** (§1 R1.3, §5 R5.2): FDS's per-mesh, init-time `WALL`/`BOUNDARY_*` tables have to become regrid-stable, box-independent objects. This is the biggest refactor outside pressure.
2. **Pressure algorithm change** (§4): composite MLMG replaces FFT+iterations. It is numerically better, but requires V&V re-baselining and a performance check against FISHPACK (P2).
3. **Stretched grids / anisotropic refinement / 2-D-as-runtime-flag** (§1 R1.1, R1.5; §3 R3.1): input-compatibility breaks affecting 17 + 287 V&V cases.
4. **Smokeview compatibility** (§9): needs a product decision. FireX VTKHDF is a better bridge than native Smokeview.
5. **Particles with attached solid state** (§6).
6. **GEOM/CC_IBM ↔ EB**: two independent geometry engines. This decision belongs to ADR-003.
