# FR-072 / FR-076: Smokeview output from an AMR run, AMReX-side notes (input to ADR-004)

| Field | Value |
|---|---|
| Status | Design notes (hypothesis), AMReX Integration side. Nothing here is decided. It answers the architecture owner's five questions on approach (A), which is already chosen for ADR-004. |
| Sources | AMReX `99ddfda` (`Src/...`). FDS `Source/...` at FireX `36975d765f`: the worktree HEAD is one docs-only commit above it, and `git diff 36975d765f HEAD -- Source` is empty. Project docs `docs/...`. |
| Builds on | `docs/amrex/multi-gpu-mpi-spec.md` §5 (AsyncOut, background thread, VisMF async). That material is referenced here, not repeated. |
| Verification | Source reading only. Nothing was built or run. Items marked **unverified** have not been checked in source or by test. |

## 0. Requirement basis (what the docs actually say)

- **FR-072** (`docs/requirements.md:320-324`): in AMR mode, `.smv`, slice, boundary, 3-D smoke and particle files carry refined-level data. Open points 1-4 are at `:321`. Plot3D is open point (2). **FR-074** (`:327`) adds VTK. **FR-073** (`:325`) makes AMReX plotfiles optional.
- **FR-076** (`requirements.md:334`): Smokeview output must be "identical in content" at any rank count and after redistribution. Each file must be written once and completely, under an ownership-independent mapping. The test allows "equal at the output format's precision where a gather changes summation order". FR-076 covers Smokeview-format output only. AMReX plotfiles and checkpoints are outside it.
- **D-036** (`docs/README.md:74`): AMReX load balancing overrides `MPI_PROCESS`. **R-48** (`risks.md:56`): the risk is files that are missing or inconsistent after boxes move. **A-40** (`README.md:130`) and **A-48** (`README.md:138`) own this design.
- **D-027** (`README.md:65`): output is written on the host and is the only exception to the GPU time step.
- **D-028** (`README.md:66`) and **FR-005 (ii)** (`requirements.md:78`) define the project's exact, order-independent sum: fixed-point int128 accumulation starting at the per-cell/per-box level.
- **Correction:** the task brief cites "D-025's background writers", but **D-025 is "FR-016 gating"** (`README.md:63`). The background-writer basis is actually NFR-047's cores-per-GPU clause and its A-43 constraints (`requirements.md:465-466`), R-46 (`risks.md:54`), A-43/A-44 (`README.md:133-134`) and spec §5. This note designs to those.

## 1. Q1: fixed output meshes and a writer rank fixed at setup

**Answer: confirmed. No AMReX obstacle was found.**
- `DistributionMapping(const Vector<int>&)` and `(Vector<int>&&)` (`Src/Base/AMReX_DistributionMapping.H:103,109`) only store the map (`AMReX_DistributionMapping.cpp:311-321`). They do no validation or balancing.
- The only constraints I found:
  - one entry per box: `BL_ASSERT(dm.ProcessorMap().size() == bxs.size())` (`AMReX_FabArrayBase.cpp:199`), which fires in debug builds only;
  - every entry must be a rank of the current AMReX communicator, because communication translates ranks with `ParallelContext::global_to_local_rank` on `CommunicatorSub()` (`AMReX_FabArrayCommI.H:901-907, 974-989`).
- **The output BoxArray need not be disjoint.**
  - `CPC::define` intersects each destination box with the source BoxArray independently (`AMReX_FabArrayBase.cpp:405-442`). Disjointness only switches the thread-safe flags (`:444-460`).
  - Because the output BoxArray is not an `AmrMesh` level, `max_grid_size`, `blocking_factor` and `ChopGrids` never touch it.
- **Still recommended** (design rules, not AMReX requirements):
  - output patches on one level should be disjoint, so each cell has exactly one file value (Smokeview behaviour with overlapping meshes is **unverified**);
  - each level-L patch should be aligned to the ratio of level L-1, and `coarsen(patch_L)` should lie inside the union of level L-1 output patches, so the piecewise-constant fill (§2) always has a filled coarse source.
- **Writer assignment:** use `writer(m) = f(m, NProcs)`, fixed at setup and never derived from the compute DM. It may differ between rank counts, because each file is written whole by one rank.
  - **Caveat from FDS:** the setup `.smv` concatenates each rank's mesh blocks **in rank order**. The mesh loop keeps only owned meshes (`Source/dump.f90:2457-2459`). The blocks are then joined by `MPI_EXSCAN` offsets (`dump.f90:2688-2690`) or by a rank-ordered `MPI_GATHERV` to rank 0 (`:2706-2722`).
  - So round-robin writers would reorder the `.smv` when the rank count changes. Either rank 0 writes the whole mesh section in output-mesh order, or writers own **contiguous, increasing** blocks of output meshes. Run-time `.smv` strings are already written in mesh order (rank 0 walks `NOM=2..NMESHES`, `main.f90:4158-4200`).

## 2. Q2: fill by ParallelCopy plus average_down, and determinism

**Answer: confirmed for cell-centred data from disjoint sources, with three hazards (H1-H3).**

1. **COPY is the default.** `ParallelCopy(src, period, op = FabArrayBase::COPY)` (`AMReX_FabArray.H:971-974`; enum `COPY=0, ADD=1` at `AMReX_FabArrayBase.H:411`). ADD is used only by `ParallelAdd` (`FabArray.H:968-970`) or an explicit `op`. The default `snghost = 0` restricts sources to valid regions (`FabArray.H:962-967`).
2. **How the copy runs.** `CPC::define` creates one tag per (source box, destination box) intersection (`FabArrayBase.cpp:359-442`) and sorts the send and receive lists (`:463-471`). `ParallelCopy_nowait` does the local tags immediately (`AMReX_FabArrayCommI.H:653-666`, via `PC_local_cpu`, whose COPY branch is a plain `copy` at `AMReX_PCI.H:27-29`). `ParallelCopy_finish` unpacks remote buffers afterwards (`FabArrayCommI.H:705-740`). The 1-rank/1-box shortcut is a straight element copy (`:454-488`).
3. **Disjoint sources give bitwise layout independence.** When the source boxes are disjoint (one level's cell-centred valid regions), every destination cell is written exactly once, with a verbatim copy of the owning source value. The result therefore cannot depend on the source DM, on the writer DM, or on local-versus-remote order. No arithmetic is involved.
4. **H1, overlapping sources (face-centred or nodal data, periodic images).** A cell covered by two sources is written twice:
   - local tags go first and remote unpacking later (item 2), so the value that survives depends on which source box shares a rank with the destination, i.e. on the DM;
   - on GPU, COPY with non-thread-safe tags falls back to atomic stores (`AMReX_PCI.H:108-111`; `AMReX_FBI.H:1234-1237`), so the order is not defined at all.

   Face velocities on a face shared by two boxes are exactly this case unless both copies are bitwise equal. **Fix:** call `OverrideSync` on face data before output. It overrides every duplicate with one value (`FabArray.H:1488-1496`), taken from the lowest grid number (`OwnerMask`, `AMReX_MultiFab.H:789`). Alternatively, build cell-centred output fields on the compute side with a per-cell kernel before the copy.
5. **`average_down` is pure per-cell arithmetic.**
   - The 3-D path (`AMReX_MultiFabUtil.H:725-847`) computes in place when `coarsen(fineBA) == crseBA` with the same DM (`:742`). Otherwise it averages into a temporary on the fine DM and ParallelCopies it (COPY) into the coarse MultiFab (`:847`).
   - The kernel `amrex_avgdown` (`AMReX_MultiFabUtil_3D_C.H:377-395`) sums the fine cells in a fixed `kref/jref/iref` loop and multiplies by `1/(rx*ry*rz)`. For the supported ratios {2,4} (FR-010) that factor is 1/8 or 1/64, an exact power of two. The result depends only on the fine values, not on the layout.
   - The geometry-taking overload forwards to this path in 3-D (`AMReX_MultiFabUtil.cpp:379-380`).
   - **Do not use `sum_fine_to_coarse`:** it ends in a `ParallelCopy(..., FabArrayBase::ADD)` (`MultiFabUtil.cpp:487-488`).
6. **Proposed fill for the output MultiFab `O_L` of level L** (writer DM, fixed program order; each `ParallelCopy` is blocking, `FabArrayCommI.H:336-338`):
   - (a) `ParallelCopy` `O_{L-1}` into a temporary on `coarsen(outBA_L, r)` with the writer DM, then inject piecewise-constant per cell into `O_L` (a local kernel with no reductions);
   - (b) `ParallelCopy` the level-L composite valid data on top (COPY). If the driver has not already synced covered cells, first `average_down` level L+1 into a temporary copy of level L (item 5).

   `O_0` is filled from level 0 only. AMReX's `InterpFromCoarseLevel` with `pc_interp` (`AMReX_Interpolater.H:420, 955`) would do (a), but it needs BC functors. Its determinism was not traced (**unverified**), so (a) is preferred. Where no level-L data exists, the file carries the injected coarse value. A flag value instead is an owner choice (open question 2).
7. **H2, the compute BoxArray depends on the rank count by default.** `refine_grid_layout = true` (`Src/AmrCore/AMReX_AmrMesh.H:52`) makes `ChopGrids` split level 0 and new fine grids until there are at least `NProcs` boxes (`AMReX_AmrMesh.cpp:491-540, 569, 954`). The output path stays layout-independent. The composite solution itself, however, differs across rank counts wherever the physics is not box-split-independent (FR-005 (iii), pressure to eps_H). `OwnerMask` also follows grid numbers. **Consequence:** test output-path bit-identity on a frozen field, and whole-run FR-076 separately, using `amr.refine_grid_layout=0` or the FR-076 precision clause.
8. **H3, order-dependent small-data sums.**
   - HRR and mass: `Q_DOT_SUM`/`M_DOT_SUM` use `MPI_REDUCE(MPI_SUM)` (`main.f90:4602-4607`), `MASS_DT` likewise (`:4638-4641`).
   - Devices: per-mesh subdevice sums, then `MPI_ALLREDUCE` (`main.f90:4431-4435, 4479`).
   - AMReX's `ParallelDescriptor::ReduceRealSum` is a plain `MPI_Allreduce(MPI_SUM)` (`AMReX_ParallelDescriptor.H:1321, 1351-1352`).

   Per-box floating-point partial sums change with the box split and the rank order. **Fix:** reuse the FR-005 (ii)/D-028 exact fixed-point accumulation, from per-cell values to an integer MPI sum, for `hrr.csv`, mass and integral devices. The cost is negligible, and the bytes become identical. FR-005 (ii) lists "device statistics feeding controls" but not `hrr.csv` (open question 3). MIN/MAX devices are order-independent. Tie handling in `MINLOC`/`MAXLOC` was not checked (**unverified**).
9. **Host placement.** `O_L` can live in `The_Pinned_Arena()`, where a GPU-build ParallelCopy unpacks with device kernels into pinned memory. Alternatively it can live on the device and be snapshotted to host the way `VisMF::AsyncWriteDoit` does (`AMReX_VisMF.cpp:2471-2494`). Which is faster is a measurement for A-44 (**unverified**).

## 3. File layout and owner rule

- **FDS naming:** all files are named per mesh in `ASSIGN_FILE_NAMES` (`dump.f90:292`), whose loop keeps only meshes this process owns (`:453-455`).
  - `CHID_NM.xyz` (`:464`);
  - `CHID_NM_N.sf` plus `.bnd`/`.rle` (`:501-504`);
  - `CHID_NM_N.s3d` (`:488-490`);
  - `CHID_NM_N.bf` (`:524-525`);
  - `CHID_NM.prt5` (`:552-553`).
- **Per-mesh writers:** `main.f90:623-627` calls `DUMP_MESH_OUTPUTS(T,DT,NM,...)` for the rank's own meshes. That routine starts with `POINT_TO_MESH(NM)` (`dump.f90:89`) and dispatches (`:167-227`) to:
  - `DUMP_PART` (`:4496`);
  - `DUMP_ISOF` (`:4860`);
  - `DUMP_SMOKE3D` (`:5052`);
  - `DUMP_SLCF` with `IFRMT` 0 = slice, 2 = 3-D slice, 1 = Plot3D (`:6746`);
  - `DUMP_BNDF` (`:11950`);
  - `DUMP_PROF` (`:11362`).
- **Other `POINT_TO_MESH` callers:** `UPDATE_GLOBAL_OUTPUTS` (`:63`), the VTKHDF writers (`:4749, 5191, 7466, 12250, 12488`), `WRITE_CFACES` (`:6318`) and restart (`:3897, 4090`).
- **Global files:** `DUMP_HRR` (`:11830`) and `DUMP_DEVICES` (`:11155`) run on rank 0 after the reductions in H3. The `.smv` mesh section is `WRITE_SMOKEVIEW_FILE` (`:1766`, loop `:2457`), called once at `main.f90:216-221` and later only appended to (`main.f90:255`; strings at `:651`). All of these lines were checked.
- **Owner rule:** every file is a function of exactly one output mesh m (or of the global state, for `.smv`/CSV). Only `writer(m)` writes it, in FDS i,j,k order from `O_L`. Its bytes then depend only on m's data, whatever the compute DM. `.smv` mesh blocks are written in output-mesh order (§1 caveat).
- **What a snapshot writer needs as arguments** instead of `POINT_TO_MESH`:
  - mesh extents and coordinates (`IBAR/JBAR/KBAR`, `X/Y/Z`, `RDX`...);
  - the output fields with the ghost layers the writer interpolates from;
  - for BNDF, the wall/patch arrays (`WALL`, `CELL_INDEX`, boundary patches) mapped onto the output mesh;
  - the particle list;
  - the per-mesh output state (clocks and counters, `LU_*`, `FN_*`).

  **Stage 1 (synchronous):** the writer rank fills `MESHES(m)` from `O_L` and calls the unmodified `DUMP_*` on the main thread. **Stage 2:** refactor to explicit arguments so the writer can run in a background job.
- **BNDF on refined levels:** this needs an owned-wall gather onto output meshes (FR-045 one active record, R-55). It is not designed here (open question 4).
- **AMReX `VisMF`/plotfiles are not rank-independent.** The file number depends on the writing rank and `nfiles` (`AMReX_NFiles.H:153-163`), and FABs are written by their owners. They are therefore unsuitable as the Smokeview path, but fine for FR-073 plotfiles and checkpoints, which FR-076 does not cover.

## 4. Q3: writers on host cores beside the GPU rank

**AsyncOut** (facts beyond spec §5.1):
- one `std::thread` per rank with a FIFO queue (`AMReX_BackgroundThread.cpp:5-8, 19-50`);
- `Finish` drains the queue (`:52-60`), and `Submit` takes any `std::function` (`AMReX_AsyncOut.cpp:95-103`);
- `Submit` dereferences the thread unconditionally, so callers must guard with `UseAsyncOut()` (`:70`);
- the `MPI_THREAD_MULTIPLE` abort fires whenever `async_out && nfiles < NProcs` (`:35-44`), even if our jobs make no MPI calls, so set `amrex.async_out_nfiles >= ranks` (the default of 64 is clamped to `NProcs`, `:28-32`);
- I saw no queue-depth limit (`BackgroundThread.cpp:38-50`), so snapshot memory must be bounded by calling `Finish()` before the next output or by our own counter;
- Fortran I/O on a non-main thread is **unverified** (spec §5.1).

**Separate writer ranks:**
- **ParallelCopy cannot reach ranks outside the AMReX communicator.** Both DMs are interpreted on `ParallelContext::CommunicatorSub()`, with rank translation (`FabArrayCommI.H:901-907, 974-989`). `amrex::Initialize(MPI_Comm, ...)` (`AMReX.H:109-113, 187-188`) and `ParallelContext::push/pop` (`AMReX_ParallelContext.H:105-111`) only change which communicator that is. This confirms the owner's belief.
- **AMReX does have `amrex::MPMD::Copier`** (`AMReX_MPMD.H:26-52`). It copies FabArray data between two AMReX programs over `MPI_COMM_WORLD`:
  - tags come from box intersections and are sorted (`AMReX_MPMD.cpp:250-277`), and unpacking uses COPY (`MPMD.H:180-186`), so it is deterministic for disjoint sources, like ParallelCopy;
  - limits: exactly two programs, identified by `MPI_APPNUM`, argc or executable hash (`MPMD.cpp:57-93`); the second program's ranks must follow the first's (`:148-156`); construction is collective; `send` blocks in `Waitall` until the writer receives (`MPMD.H:113-118`).

**Recommendation: AsyncOut inside the compute job.**
- The main thread does the collective ParallelCopy into the writer DM (§2). The writer rank takes a host snapshot of `O_L` and submits one job per output mesh. The job writes only its own files, with no MPI and no module-state access.
- This works in route (a) (writers spread over the N ranks per GPU) and in route (b) (the thread runs on a spare host core). It needs no second executable and no MPMD launch.
- Until Stage 2 and the Fortran thread check are done, run the Stage 1 writer synchronously.
- Keep the MPMD writer program as the fallback if A-44 shows compute ranks stalling on snapshot or write. Its cost is a separate writer executable (or a distinct argc) and blocking sends.

## 5. Q4: particles

- **The AMReX id is rank-dependent.** `idcpu` packs a 39-bit id and a 24-bit creating-rank field (`Src/Particle/AMReX_Particle.H:40-53, 57-88`). The id comes from a per-process static counter (`the_next_id`, `:421, 591`), incremented under `omp atomic capture` (`:595-605`). Both parts therefore depend on the creating rank and on the thread schedule. **Nothing layout-independent is built in**, which confirms the owner's expectation.
- **User int components exist:** compile-time `NArrayInt`/`NStructInt` (`AMReX_ParticleContainer.H:139-158`) and run-time `AddIntComp` (`:1349, 1377`).
- **FDS tags today:** `LP%TAG` (`type.f90:396`, a default INTEGER) is assigned `PARTICLE_TAG = PARTICLE_TAG + NMESHES` at insertion (`part.f90:413-414, 777-778, 1275-1276, 1373-1374`), and the counter is seeded with the mesh number (`init.f90:923`). Tags are unique only because meshes are fixed. They go to `.prt5` in storage order (`dump.f90:4546-4553`).
- **Proposal:**
  - store a 64-bit key in two int components: (source id [INIT/device/surface class], source-local element such as a point index or a global (level, i, j, k, ior) wall face, per-element insertion count);
  - each step, assign the 32-bit `.prt5` TAG as base + position in the key order of that step's new particles, using an integer `MPI_Exscan` of per-source counts in fixed source order (exact, order-independent);
  - at output, `copyParticles(pc, local=false)` into a container defined on the level-0 output BoxArray and writer DM, which redistributes (`AMReX_ParticleContainer.H:703-713, 247`);
  - sort by TAG per output mesh before calling the FDS `.prt5` writer.
- **Uniqueness caveat:** 32-bit TAG uniqueness over long runs is the same limit FDS has. Whether Smokeview uses TAG to track particles is **unverified**.

## 6. Q5: conflicts with `multi-gpu-mpi-spec.md`

1. **Conflict.** §5.3 says route (a) parallelises host work "with no code change: each rank handles only its own boxes, as FDS does per mesh today" (`:234`), and the default hypothesis (`:243`) keeps "FDS writers unchanged". In AMR mode, per-rank own-box writing is exactly what FR-076/R-48 forbid. **Amend:** write per writer-owned output mesh (§1-§3), with Stage 1 as the minimal change.
2. **Correction.** §5.1 "no facility for dedicated I/O ranks" (`:207`) and §8 item 11 "AMReX offers neither" (`:340`) are partly wrong: `amrex::MPMD::Copier` exists, with the two-program limits in §4.
3. **Consistent, add a note.** `async_out_nfiles` and `MPI_THREAD_MULTIPLE` (`:203, 239, 283`) match the source. Add that the abort applies even though the FDS jobs make no MPI calls, and that `Submit` needs the `UseAsyncOut()` guard.
4. **Consistent.** The `POINT_TO_MESH` and thread-safety rows (`:235, 339`) agree with this note and motivate Stages 1 and 2.
5. **Gap.** The §5.3 reproducibility row cites only FR-005 (iii) (`:240`). The output path must also meet FR-076 (H1-H3).
6. **Consistent.** Snapshot memory (§8 item 12, `:341`) also applies to `O_L` snapshots.

## 7. Approach (B), for the record

Writing refined boxes as Smokeview meshes would need a new mesh list at every regrid, which FDS never writes (§3). The box numbering is not stable: `ChopGrids` depends on `NProcs` (H2), and regridding rebuilds the BoxArrays. So file names and mesh indices would change with regrids and with the rank count. Smokeview support for a changing mesh set is **unverified** (FR-072 open point 4). The fixed per-level output patches of §1 give refined-resolution data with a static mesh list, which is why (A) holds up from the AMReX side.

## 8. Verification tests (proposals, not run)

1. **Output path on a frozen field:** analytic cell- and face-centred data on two levels; compare file bytes at 1, 2 and 4 ranks. Repeat with the DM strategy set to `ROUNDROBIN`, `KNAPSACK` and `SFC` (`AMReX_DistributionMapping.H:58, 188`), and with a random permutation passed through `DistributionMapping(Vector<int>)`. Every file and the `.smv` must match bitwise.
2. **H1 probe:** make shared-face copies differ deliberately. Without `OverrideSync` the bytes should change with the DM; with it they should not.
3. **Two-level run with at least one regrid and a forced redistribution between outputs:** use `amr.refine_grid_layout=0` for the bitwise form. Otherwise apply the FR-076 precision clause.
4. **Async on/off:** `amrex.async_out=0/1` with `async_out_nfiles >= ranks`; the bytes must be identical.
5. **Particles:** the same insertion history at 1, 2 and 4 ranks must give identical `.prt5` after the TAG sort.
6. **HRR/devices with exact sums:** `hrr.csv` and `devc.csv` must match bitwise across ranks. A float-sum control run should show the difference.

## 9. Open questions for the architecture owner

1. Should writers own contiguous blocks of output meshes, or should rank 0 write the whole `.smv` mesh section (§1 caveat)?
2. Where no level-L data exists, should the file carry the injected coarse value or a flag value?
3. Should FR-005 (ii) exact sums be extended to `hrr.csv`, mass and integral devices, so FR-076 holds bitwise and not only at format precision?
4. BNDF on refined output patches: which wall records feed a fixed output mesh (FR-045 ownership)? Is Plot3D in scope (FR-072 open point 2)?
5. Is `amr.refine_grid_layout=0` acceptable for FR-076 acceptance runs, or must the tests pass with the default?
6. Is the `.prt5` TAG scheme in §5 acceptable, and is the 32-bit TAG limit tolerable?
