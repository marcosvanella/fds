# ADR-004: Smokeview-format output on the AMR hierarchy (FR-072, FR-076)

| Field | Value |
|---|---|
| Status | **Accepted** (owner answered §8). The §6 items are verification follow-ups, not open decisions. |
| Type | Full ADR |
| Date | 2026-09-26 |
| Owner | AMR Chief Architect (with AMReX Integration Lead) |
| Deciders | Project owner (§8), AMR Chief Architect, AMReX Integration Lead |
| Evidence base | FireX `36975d765f` (Source/ unchanged on branch `FDS-AMReX`); AMReX `99ddfda`; Smokeview tag `SMV-6.11.2` (`27cbe596`). AMReX side: `amrex/fr072-output-amrex-notes.md`. |

Supersedes the A-40 "test both candidates" plan. Option (b) is now rejected from source evidence (§4), so only option (a) needs to be tested (§7).

## 1. Context

- FR-072 requires AMR runs to write `.smv`, slice, boundary, 3-D smoke and particle files that carry refined-level data. Level-0-only output is rejected by the owner decision. Output is written on the host (D-027). FR-074 requires VTK as well.
- FR-076 requires Smokeview output whose content does not depend on rank count or on boxes moving between ranks (D-036).
- The background-writer basis is NFR-047's cores-per-GPU clause (requirements.md:465-466), R-46 and A-43/A-44. D-025 is FR-016 gating and not relevant here.
- FDS declares the Smokeview mesh set once at setup, and a mesh block is never added later (FR-072 rationale; `WRITE_SMOKEVIEW_FILE` dump.f90:1766, mesh loop at 2457-2541; append-only afterwards, main.f90:255). Slice and boundary files bind to a fixed mesh and fixed index bounds (dump.f90:1215, 1260-1264, 1471-1473). Plot3D binds to a per-mesh grid file (dump.f90:1056, 6842).
- **Smokeview 6.11.2 does not resolve overlapping meshes of different resolution.** Point lookup takes the first listed mesh, whatever its resolution (`smv_geometry.c:710-735`). Embedded-mesh masking works only when a mesh lies wholly inside one other mesh (`smv_geometry.c:1310-1437`), and it is used only for slices behind the optional `.ini` switch `SKIPEMBEDSLICE`, default off (`readsmv.c:4307`, `IOslice.c:5125-5303`). The cell-centred slice mask looks mis-indexed (node fill at `smv_geometry.c:1430`, cell read at `IOslice.c:4165`). Boundary patches from both meshes are drawn (`IOboundary.c:1375-1389, 3666`). Smoke3d overlap zeroing runs only at `.smv` read (`IOsmoke.c:2890-2963`, `readsmv.c:2870`). Plot3D has no handling (`IOplot3d.c:1309`). Particles are drawn at absolute positions, unclipped (`IOpart.c:229-230, 595-608`). The "hide overlaps" key is dead code (`smokeviewvars.h:1867`, `callbacks.c:1918`). Each file has its own time list (`update.c:662-675`), and there is no missing-data marker within a frame (`shared/readslice.c:89`).
- FDS output reductions are plain MPI sums (main.f90:4479, 4602-4607, 4638-4641). AMReX `ReduceRealSum` is the same (ParallelDescriptor.H:1351-1352).

## 2. Decision drivers

1. Owner decision: refined-level data must appear in Smokeview-format files.
2. FR-076: output content is independent of rank ownership, rank count and redistribution.
3. The output must display correctly in unmodified Smokeview 6.11.2 (no Smokeview changes; §1 overlap facts).
4. Keep `dump.f90` writers unchanged at first (D-027 host output; readability for FDS developers), then move them to background writers (NFR-047).
5. Output size and write time stay bounded.
6. Uniform mode (max level 0) stays byte-identical to baseline (FR-072).

## 3. Decision

**D1. Static, disjoint output meshes.** At setup the domain is split into output meshes that **do not overlap**:
- Outside the *refinable region*, output meshes are the level-0 `&MESH` blocks, split into boxes where they border the region.
- Inside the refinable region, output meshes are tiles at the *output level* `L_out`. The tiles are aligned to the blocking factor, and each is at most `max_grid_size` at that level.
- The refinable region is the union of the finer `&MESH` static boxes (IR-002) and any refinable boxes the user declares (spec delta S2). Where it is empty, AMR mode writes level-0 meshes only and warns once. `L_out` defaults to max level and can be capped by the user (S2).
- In uniform mode the output mesh set is exactly the `&MESH` set, so baseline output is unchanged.
- Rank 0 writes the whole `.smv` at setup, including every mesh block, `OBST` and `IBLANK`, in output-mesh order. Per-step entries stay rank-0 appends as today. The `.smv` does not depend on rank count; this addresses the rank-order concatenation in dump.f90:2457-2459, 2688-2690, 2706-2722.
- `OBST` blocks for an output mesh are snapped with the D-009 rule at that mesh's cell size.

**D2. Fill rule (answers FR-072 open point 1).** At each output time every output cell holds composite data at the output mesh's resolution:
- A finer level contributes through `average_down`, which is plain per-cell arithmetic with exact scale factors.
- A coarser level contributes through piecewise-constant injection (a per-cell copy, no arithmetic).
- Face and node data call `OverrideSync` before any copy (hazard H1 in the AMReX notes).
- A new output quantity, `AMR_LEVEL` (the level whose valid data filled the cell), is available for slices, so users can see the true resolution (S3).
- Derived slice quantities are computed on the output mesh from the filled primitive fields by the unchanged FDS output functions. Integrals and devices never come from output meshes (D6).

**D3. Writer layout independent of compute layout.** Each output mesh has a fixed writer rank, set once at setup. Writer ranks own contiguous, increasing blocks of output meshes, balanced by cell count. The layout uses `DistributionMapping(const Vector<int>&)` (AMReX_DistributionMapping.H:103), independent of the compute DistributionMapping and of load balancing (D-036).
- A `ParallelCopy` in COPY mode (FabArray.H:971-974) fills an output MultiFab. For non-overlapping cell-centred sources every destination cell receives exactly one copy, so file content is byte-identical across rank counts and redistributions **for a given solution**.
- Ghost cells of output meshes are filled by the same FillPatch path as the compute levels. Output meshes are grouped by resolution for this.
- Solution differences across box splits are governed by FR-005, not by this ADR (hazard H2).

**D4. Particles.** Each particle is written to exactly one output-mesh file: the mesh that contains its position, lower faces inclusive. Within the file, particles are sorted by a layout-independent `TAG`.
- `TAG` comes from a key of (source, element, insertion count) held in particle int components (ParticleContainer.H:1349).
- The key is numbered with an integer `MPI_Exscan` in fixed source order.
- Particles move to the writer layout with `copyParticles`.
- This replaces FDS's mesh-seeded tag (part.f90:413-414, init.f90:923) in AMR mode.

**D5. Boundary files.** Each output mesh's boundary patches are its wall faces at the output resolution, and a patch value is taken from the owning wall face (FR-045):
- where the owner is coarser, the owner's value is copied;
- where finer records exist (FR-041b Option A, `adr/drafts/ruling-FR041b-G2.md`), the output value is their area-weighted mean in the fixed face-key order.
Coarse faces under a finer owner are never written twice, because D1 output meshes do not overlap.

**D6. Global outputs and devices.** `CHID_hrr.csv`, `CHID_mass.csv` and `CHID_devc.csv` are computed on the composite grid (uncovered cells and owned faces) with the D-028 exact fixed-point sum, and written by rank 0 (FR-070, FR-071, H3). `sum_fine_to_coarse` is not used, because it ends in an ADD (MultiFabUtil.cpp:487-488). MINLOC and MAXLOC ties go to the lowest global cell key: level-0-equivalent (i,j,k), then level.

**D7. Scope of file types.**
- 3-D smoke (both `SMOKE3D_VERSION` layouts) and Plot3D use the D1-D3 path. Plot3D is **IN**, which answers FR-072 open point 2.
- In AMR mode `&DUMP WRITE_FORMAT` defaults to `'BOTH'`, following the owner decision that AMR runs write VTK too. It stays honoured if the user sets it (open point 3).
- VTK (FR-074) writes the same output meshes through the same snapshot in Phase 9. A native overlapping-AMR VTK layout is an optional later extension.
- Radiation outputs follow this ADR (FR-077). RADF files are written per output mesh.

**D8. Background writers (NFR-047).** Staged:
- **Stage 1 (Phase 9 entry):** synchronous writing. The writer rank points `MESHES(m)` at the filled output mesh m and calls the unchanged `dump.f90` writers.
- **Stage 2:** the writers are refactored to explicit arguments (no `MESHES` pointers, no MPI). Each writer rank takes a host snapshot and submits one AMReX `AsyncOut` job per output mesh. That gives one background thread per rank with a FIFO queue (BackgroundThread.cpp:19-60).
- Stage 2 settings:
  - `async_out_nfiles` ≥ ranks, to avoid the `MPI_THREAD_MULTIPLE` abort (AsyncOut.cpp:35-44);
  - calls are guarded by `UseAsyncOut()`;
  - snapshot memory is capped at one output interval's data, with `Finish` called if exceeded.
- The fallback is the `amrex::MPMD::Copier` two-program route (AMReX_MPMD.H:26-52). A separate writer communicator inside one program is rejected, because `ParallelCopy` needs both layouts on one communicator.

## 4. Rejected alternatives

1. **Option (b): refined boxes written as Smokeview meshes per regrid.** The `.smv` mesh set is fixed at setup and only appended to afterwards (§1). A changing mesh set would need a new `.smv` per regrid or Smokeview changes. It would also make file names and mesh indices depend on the regrid history.
2. **Overlapping per-level output meshes (level-0 meshes whole, fine meshes on top).** In Smokeview 6.11.2 this works only for slices, only with `SKIPEMBEDSLICE`, only for fine meshes wholly inside one coarse mesh, and with a suspected cell-centred mask bug. Boundary, 3-D smoke, Plot3D and particles would display twice (§1).
3. **Blanking covered coarse cells per frame.** There is no missing-data marker in slice or boundary frames, and `IBLANK` is static per mesh (§1). A sentinel value would corrupt colour bars.
4. **Each rank writes its own boxes (the old multi-GPU default).** It violates FR-076 and R-48, and the AMReX spec was corrected accordingly.
5. **Level-0-only Smokeview output plus full-resolution VTK.** Rejected by the owner decision (FR-072).
6. **Output at full refined resolution wherever refinement currently exists, without a declared region.** Mesh blocks cannot follow refinement (alternative 1). Declaring the whole domain as refinable gives this result at the cost of writing `L_out` everywhere, which the user can choose.

## 5. Consequences and risks

- (+) Unmodified Smokeview displays every file type correctly, with no overlap artefacts. The `.smv` and all files are layout-independent. Stage 1 reuses `dump.f90` unchanged.
- (−) Inside the refinable region, output is always written at `L_out`, even where the solution is coarser at that moment. Size grows with (refinable volume) × r^(3·L_out). This is mitigated by capping `L_out` and by the S4 size estimate printed at setup.
- (−) Derived slice quantities are computed from averaged primitives. They are for display and are not conservative. Integrals come from D6.
- (−) Output meshes that split a level-0 `&MESH` renumber Smokeview meshes relative to the input. Devices and CSV files are unaffected.
- Risks to add: R-6x, Fortran I/O on a non-main thread (gfortran, nvfortran) is unverified (spike S-B). R-6y, output size inside large refinable regions.

## 6. Open items (owners)

- Determinism of `InterpFromCoarseLevel` / `pc_interp` for the D2 coarse injection (AMReX Integration Lead): trace it and confirm it is a pure copy.
- Whether Smokeview uses particle `TAG` for tracking (AMReX Integration Lead, during the Smokeview check).
- Boundary patches at an output-resolution boundary that crosses a thin OBST (AMR Solid Phase Lead): check them against R3-T.

## 7. Spike plan (replaces A-40)

- **S-A, FR-076 byte identity (CPU).** Run a two-level case (static fine box, plus one regrid inside the output window) with a fixed BoxArray (`amr.refine_grid_layout=0`, fixed `max_grid_size`) at 1, 2 and 4 ranks, plus one run with a forced redistribution between output times. **Pass:** `.smv`, `.sf`, `.bf`, `.s3d`, `.q` and `.prt5` files are `cmp`-identical across all runs.
- **S-B, background writer thread safety.** A Fortran unformatted-write job runs on the `AsyncOut` thread under gfortran (CPU build) and nvfortran (GPU test machine). **Pass:** the files are identical to the synchronous writer's in 100 repeats, and ThreadSanitizer is clean on gfortran. If it fails, Stage 2 uses the MPMD fallback.
- **S-C, Smokeview display.** Needs Smokeview 6.11.2 installed where the owner permits. Load S-A's output. **Pass:** zero load errors; slices, boundary, smoke and particles appear once each in the refined region at `L_out`; the `AMR_LEVEL` slice shows the level pattern.

## 8. Owner decision (answered)

The question was whether to accept refined-level data at a fixed output level `L_out` over refinable regions declared at setup, with level-0 resolution elsewhere.

**Answer: accepted.** Smokeview cannot handle a mesh count that changes within one simulation. Output that follows refinement is therefore not a goal, and no Smokeview changes are planned. D1-D8 stand as written, and S2 is no longer pending.

## 9. Spec deltas (for the Spec Lead)

- S1. FR-072: close open points 1-3 per D2 and D7. Point 4 is moot, because the mesh set never changes. Verification adds the `AMR_LEVEL` slice check and S-C.
- S2. New interface requirement: AMR-mode inputs for the refinable region (boxes) and the output level cap `L_out` (default max level). Namelist names are left to the Spec Lead and Integration Lead.
- S3. New output quantity `AMR_LEVEL` (slice, AMR mode only).
- S4. A setup report of output mesh count and bytes per output interval.
- S5. FR-076: verification per S-A (fixed BoxArray, varying rank count and redistribution; `cmp`-identical). The note that solution differences across box splits fall under FR-005.
- S6. FR-070 and FR-071: D-028 exact sums, and the D6 MINLOC/MAXLOC tie rule.
- S7. Particle `TAG` rule (D4) in AMR mode.
- S8. `WRITE_FORMAT` default `'BOTH'` in AMR mode (D7).
- S9. FR-077: RADF files per output mesh.
