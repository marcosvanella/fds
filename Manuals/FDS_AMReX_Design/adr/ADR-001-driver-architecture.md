# ADR-001: Driver architecture — C++ AMReX driver with FDS Fortran kernels vs FDS as a Fortran program on AMReX F_Interfaces

| Field | Value |
|---|---|
| Status | **Accepted (v0.7, 2026-10-02).** Owner sign-off 2026-10-02 (relayed by the project coordinator): K2 (Fortran with OpenMP `target` offload, compiled with `nvfortran`) is the default kernel style; K1 (restricted C++ `ParallelFor`) is allowed only as a per-kernel fallback. Reason given by the owner: tests show similar GPU performance at a much lower maintenance cost. The driver choice (Option A, C++ AmrCore driver) is owner-confirmed via D-027; C++ is limited to the driver and AMReX glue (D-043); the regrid-time rebuild may run on the host until Phase 11 (D-047, NFR-043); S5 reviewers were named in D-048. A1, A2 (S4a) and A3 (S5 review and sign-off) are done; S4b ran on an NVIDIA GPU. Sign-off gaps: `signoff-gaps.md` |
| Version | v0.7 (2026-10-02, **Accepted**): K2 default and K1 per-kernel fallback (owner sign-off); K2 rule 7 from NVIDIA's reply; wall-table rulings W1/W2; merge log. v0.6 draft (2026-09-29): new rule "GPU kernel data layout" (owner-approved), survey numbers, and the K1/K2 evidence summary; the K1/K2 decision stays open for the owner. v0.5 (2026-09-29): A1 done (NVHPC 26.9, CUDA 13.3); S4a done; S4b run on an NVIDIA GPU (cc 8.9) with K1 and K2 bitwise equal to the CPU; new "K2 coding rules" (no private copy of loop-invariant values, parenthesised order-sensitive sums, widened D-029 clause list) and a "GPU build flags" requirement; `ParallelAllReduce::Or` wording fixed; S5 package written; see "v0.5 changes". v0.4 (2026-09-26): owner decisions recorded: regrid-time side-data rebuild host-allowed until Phase 11 (D-047, A4, NFR-043); S5 reviewers named (D-048); S4 split into a compile-only part (gates acceptance) and a GPU-run part (does not); acceptance checklist; pressure and geometry text aligned with ADR-002 v1.1 and ADR-003 v1.1; see "v0.4 changes". v0.3.8 (2026-09-26): Q8 aligned with D-043; radiation coupling D-039 and the FR-062 ruling in the kernel interface; regrid-time rebuild recommendation (for owner confirmation); Q5 answered by ADR-004; masked level 0 needs the C++ MLMG path; see "v0.3.8 changes". v0.3.7 (2026-09-25): owner decisions: NVIDIA is the only GPU target, AMD out of scope; second-compiler check optional (ifx only); AMR mode uniform grids per level, stretched cases FDS-only; see "v0.3.7 changes". v0.3.6 (2026-09-25): pressure section reduced to pointers to ADR-002 v0.2 (maxorder 2 confirmed by P2); HYPRE per D-026; see "v0.3.6 changes". v0.3.5 (2026-09-25): D-031 pass order: one-species case, gate and sync count stated; see "v0.3.5 changes". v0.3.4 (2026-09-25): D-031 rulings: two-phase density terms on valid+2; clip flags only from uncovered valid cells; see "v0.3.4 changes". v0.3.3 (2026-09-25): D-031 pass order with two host OR reductions, coarse-side mask, A-37 closed (47 checks pass), two-phase target form; see "v0.3.3 changes". v0.3.2 (2026-09-25): D-031 ghost-depth ruling (redundant density clip over valid+1, one pre-clip `FillBoundary`); see "v0.3.2 changes". v0.3.1 (2026-09-25): species/density clipping rewritten as a layout-independent gather (D-031). v0.3 (2026-09-25): owner decision changes K2 from OpenACC to Fortran + OpenMP `target` offload, with a no-copy device-data rule; see "v0.3 changes". v0.2 (2026-09-25): owner decisions D-027 and NVIDIA-only GPU target |
| Type | Full ADR |
| Date | 2026-09-25 |
| Owner | AMR Chief Architect (for the project owner) |
| Deciders | Project owner; AMReX Integration Lead; FDS Legacy Mapper; AMR Spec & Program Lead |
| **Code base** | **FireX `36975d765f` (2026-09-24), local branch `AMReX`, this repository (read-only).** Acceptance "baseline FDS" is the same commit (Spec Lead). |
| Other evidence | AMReX `99ddfda` (shallow, CHANGES.md head 26.09); IAMR `0a63b79`, incflo `878827d`, PeleLMeX `5c21556`, ERF `78707d3`, amrex-tutorials `1a73f32`. FDS master `ce1f659` (FireX's merge-base) is quoted only where teammates cited it. |
| Depends on | ADR-003, but only partially (see "Pressure path" below); linked in the decision log, R-03 and FR-030 |

All counts were run with `rg`/`wc` on `Source` on 2026-09-25. "Refs" are textual occurrences. Teammate inventory numbers (`docs/inventory/`, `docs/amrex/mapping.md`, `docs/pressure/`) were produced against `ce1f659`. FireX changed 22 source files (+8,621/−847 lines, `git diff ce1f659 --stat -- Source`), mostly `vtkf.f90` (new), `dump.f90`, `main.f90`, `pres.f90`, `read.f90` and `radi.f90`, so their line numbers in those files have moved. Where both are quoted, the FireX figure comes first.

Cross-references: risks R-02, R-03, R-08, R-11, R-12, R-18, R-21, R-22, R-23, R-26, R-29, R-30, R-31, R-36, R-38, R-39; requirements FR-001/002, FR-005, FR-010, FR-030, FR-037, FR-039, FR-050..052, FR-072/073, IR-002, IR-004, IR-005, IR-006, IR-007, NFR-012, NFR-030/031, NFR-043, NFR-044; owner questions Q2 (answered by D-027), Q5, Q8, Q9, Q11; decisions D-012, D-021, D-027 (supersedes D-004), D-028, D-029, D-031; assumption A-37; roadmap Phase 11 / M11.

## v0.7 changes (2026-10-02): K1/K2 decided, K2 rule 7, wall rulings
- **Owner decision 2026-10-02 (relayed by the project coordinator): K2 is the default kernel style, K1 only as a per-kernel fallback.** Reason: similar GPU performance at much lower maintenance cost. Status set to **Accepted**; the P1 readability review record, "Leaning", "Layering" (d), "Rejected alternatives" and "Decision needed" item 6 are updated. The reviewer conditions on K2 become requirements: a GPU CI check for the hidden rules, a lint rejecting `private` items the loop only reads, a script flagging unparenthesised multi-term sums, generated or templated bounds plumbing, and kernel header notes for passive scalars and dropped cylindrical terms. CUDA Fortran stays as decided 2026-09-25 (profiled hot loops only, OpenMP fallback each); none is approved yet.
- **K2 rule 7 rewritten** from NVIDIA's reply (default `distribute parallel do` for kernels with callees; guarded `bind(teams,parallel)` option; no "initialise before call" workaround). Rules 3, 4 and 6 aligned.
- **New section "Wall tables and wall-state ownership"** (rulings W1 and W2) and an **upstream merge log**.
- New draft `drafts/gpu-staging-plan.md`.

## v0.6 changes (2026-09-29, draft): GPU data-layout rule and K1/K2 evidence summary
- **New rule, owner-approved: "GPU kernel data layout" (IR-007 addition)**, in its own section after "K2 coding rules and GPU build flags". GPU kernels take flat arrays and scalars only. The rule is backed by the call-graph survey numbers (`docs/inventory/gpu_callgraph_survey.md`).
- **K1/K2 evidence summary** added to "Leaning" (performance, merge burden, single source). The evidence is complete except one closing GPU T1 run. The coordinator recommends K2 as the default with K1 as a per-kernel fallback, conditional on that run and on the S5 reviews. **The owner has not decided, and this draft records no decision.** The S5 reviewers still split (Species & Combustion Lead: K1; FDS Legacy Mapper: K2 with CUDA Fortran for profiled hot loops).
- **New section "Upstream FireX merge protocol for ported kernels"** (owner request): naming rule, rename script, kernel-to-upstream mapping table with a daily merge check, handling of edits that add or change a field a kernel reads (bitwise gate), and the host-only fallback. **K2 coding rule 7 (kernel shape, S4d):** `target teams distribute parallel do` for kernels with callees or local arrays; the allowed-subset text and the rejected-alternatives line are updated to match.
- Status stays "Proposed"; acceptance checklist item 3 (S5 plus owner sign-off) is still open.

## v0.5 changes (2026-09-29): A1 and S4a done, S4b run, K2 coding rules, S5 package
Evidence: `docs/amrex/s4-cuda-mass-findings.md` (§4-§10), kernels in `src-s4/amrex/s4_mass/` at `b51f4361b3`.
- **A1 done.** NVIDIA HPC SDK 26.9 (CUDA 13.3 default; `docs/build/nvhpc-sdk.md`) is installed. It provides `nvfortran` and `nvc++`, so the items that waited on A-31 are closed.
- **A2 (S4a) done.** K1 and K2 compile for CUDA (sm_80, no spills; ptxas and `-Minfo` tables in findings §4). On the CPU (K1 build, K2 host fallback) both are byte-identical to P1/FDS: 60 comparisons and 69 checks at 1, 2 and 4 ranks. Not ported: MP5, `N_ZONE>0`, refinement ratio other than 1. Effort (findings §7): K1 279 code lines for 12 kernels; K2 372 lines plus 96 lines of C++ wrappers, a mechanical port from K1 whose extra work was bounds plumbing, `private` lists, the `is_device_ptr` clause and the mixed nvcc/nvfortran link.
- **S4b ran** (it is not a gate; it was done early) on an NVIDIA GPU with compute capability 8.9, one GPU and one process. After the two fixes below, K1 and K2 each pass 29/29 comparisons at T1 and T0 (bitwise) and 36 checks, with clip counts equal to P1 (978, 2520, 154354, 7970) and a bitwise run-to-run repeat. `is_device_ptr` carries AMReX's device pointers to nvfortran, so S4's "Overturns if" condition (K2 cannot share AMReX's device memory) is **not** met and K2 stays a candidate. Single-sample timings, informational only (SUPERBEE, 32³ box, 8 steps): loop 6.1 ms (K1) and 8.3 ms (K2). K2 has 144 host synchronisation points per case (9 per stage and level). Removing the post-K2 syncs, their cost and a profile (R-39) are not measured yet.
- **Two K2 failures found on the device**, both invisible on the CPU paths:
  1. An nvfortran 26.9 offload defect: a `private` copy of a loop-invariant value (`QMIN = QMIN_IN` at the top of the collapsed loop) is handled wrongly, every cell is flagged as clipped and mass is created (total density +306 %).
  2. nvfortran reassociates an unparenthesised sum at `-O2` (gfortran and the C++ K1 do not), which gave last-bit differences (T0 lost, T1 kept).
  Both are fixed in the kernel source without changing the numerics. K1 showed neither. This is evidence for the S5 reviewers to weigh, not a decision.
- **D-029 clause list widened** (see "K2 coding rules"): `map`, `collapse`, `private` and `is_device_ptr`. nvfortran 26.9 rejects `has_device_addr` (syntax error) and `c_f_pointer` inside a target region fails to compile, so the `is_device_ptr` route on explicit-shape dummies, which this ADR named as the fallback, is the rule for nvfortran; gfortran keeps `has_device_addr`. The Spec & Program Lead owns the D-029/IR-007 wording.
- **GPU build flags are a requirement** (see "GPU build flags"). With FMA contraction on, the two-level case fails T1 (about 1.6e-3).
- **`ParallelAllReduce::Or` wording** in "Species/density clipping": in AMReX 26.09 it takes a single value, so the packed reductions use `Max` over 0/1 integers.
- P1's two-level reference predates the v0.3.3/v0.3.4 rules; S4 reproduces it with `s4.coarse_mask=0`, and the covered-cell checks (findings §6, sections 6, 7, 7b) cover the new rules.
- **S5 answers so far (2 of 3): split, K1 vs K2 + CUDA Fortran; recorded in the review record. Both ports need a header note on passive scalars and on the dropped cylindrical terms.**
- **S5 package written** (`drafts/s5-readability-review-package.md`) for the three reviewers. The Chief Architect accepts or amends this ADR when the answers are in.
- Acceptance checklist status: 1 (A1) done, 2 (A2) done, 3 (S5 review plus owner sign-off) open, 4 (Accepted) waits on 3.

## v0.4 changes (2026-09-26): owner decisions; acceptance checklist
- **Regrid-time side-data rebuild (A4, owner decision D-047, spec v0.4.28, NFR-043).** The rebuild of per-box side data at regrid (`WALL`, `CELL_INDEX`, `EXTERNAL_WALL`, wall records, particle bookkeeping) may run on the host (CPU) until Phase 11; it must be device-capable by Phase 11 (M11). It is not counted in the device time step before then. The v0.3.8 recommendation is adopted unchanged, and "Layering" (c) and "Decision needed" item 5 are closed.
- **S5 reviewers (owner decision D-048, spec v0.4.28; NFR-044 and roadmap P1 name them):** the project owner, the AMR Species & Combustion Lead and the FDS Legacy Mapper. The owner signs off the outcome. "Decision needed" item 6 is reduced to that sign-off.
- **S4 split.** *S4a (compile-only, gates acceptance):* NVIDIA HPC SDK installed (A-31); the mass kernel written as K1 and K2; both compiled for CUDA (nvc++/nvfortran, D-029); K1 CPU build and K2 host fallback match the shimmed kernel at T1; effort per variant and the device-data mechanism accepted by the compiler (`has_device_addr` or the `is_device_ptr` + `c_f_pointer` fallback) recorded. *S4b (GPU run, does not gate acceptance):* on the owner-provided NVIDIA test machine, both variants match at T1 on the device, and K2's ordering against AMReX's CUDA stream and its sync cost per step are measured (R-39). S4b runs before Phase 11 kernel extraction starts. If the review picks K2 and S4b then meets S4's "Overturns if" condition, only the kernel style reopens, and K1 is the fallback.
- **Aligned with ADR-002 v1.1 and ADR-003 v1.1 (D-046):** the pressure solve sits behind a solver-agnostic interface (MLMG, or assembled-matrix HYPRE PCG + BoomerAMG), both in layer (a); a device HYPRE build becomes a GPU-phase work item. EB is one of two complex-geometry candidates, and it still needs C++ if chosen.
- **Acceptance checklist (exactly what is left):**
  1. A1: NVIDIA HPC SDK installed (compile-only; AMReX Integration Lead with a build chief).
  2. A2: S4a passes as defined above (AMReX Integration Lead; running now, compile-only).
  3. A3: S5 review of both variants by the three named reviewers; outcome, reasons and any CUDA Fortran kernels recorded in the "P1 readability review record"; owner sign-off.
  4. The Chief Architect records the outcome, sets Status to Accepted and marks the rejected kernel style in "Rejected alternatives". Writing only.
  Nothing else gates acceptance: no owner question is open for ADR-001; GPU test hardware (Q11 (b)), S4b and device timing (R-39) are implementation items.

## v0.3.8 changes (2026-09-26)
- **Q8 answered (D-043):** C++ is used for the driver and the AMReX glue only (time loop, FillPatch and flux registers, regrid, pressure solvers, particle container, checkpoint). Physics stays in Fortran unless the NFR-044 readability review picks K1. "Decision needed" item 1 is closed; A-30 is answered.
- **Radiation coupling (D-039, FR-062 ruling, `drafts/rulings-IR008-FR062.md` §2.1):** the default and only required mode is a fully parallel per-box lagged sweep. It adds one kernel-interface rule (below): the per-box sweep kernel reads box-face intensities from the previous exchange only (double-buffered), and the face exchange is a host-side step between passes, the same for same-rank and cross-rank faces. No ordering between boxes; `RADIATION_ITERATIONS` gives K passes per step.
- **Regrid-time side-data rebuild (NFR-043), recommendation for owner confirmation:** host execution is allowed through the CPU phases; the rebuild must be device-capable by Phase 11 (M11). Rationale: it runs once per regrid interval, not every step, and S3 bounds its cost (< 10 % of step time at a 10-step interval, R-26); keeping it on the host until then avoids porting `init.f90` wall setup before the kernels it serves. This replaces the open point in "Layering" (c).
- **Q5 answered by ADR-004 (accepted):** Smokeview-format output on static output meshes plus VTK; output stays on the host under NFR-043's I/O exception. "Decision needed" item 4 is closed.
- **Masked (non-box) level 0** (`drafts/ruling-nonbox-level0.md`) uses single-level masked MLMG with an overset mask, which is C++-only (layer (a)). Consistent with D-043.

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
- **GPU target: NVIDIA only; AMD is out of scope** (owner answer to charter Q11 (a), 2026-09-25; v0.3.7). Q11 (c) is closed (stretched cases FDS-only); Q11 (b) test hardware stays open.
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
- Toolchains on the development machine (checked 2026-09-25): gfortran; Intel oneAPI 2026.1 at `/opt/intel/oneapi` with `ifx` 2026.1.1 and `icpx` (ifx offload targets Intel GPUs, not NVIDIA; it serves as K2's optional, non-gating second-compiler check, not tried); NVIDIA HPC SDK 26.9 with `nvfortran`, `nvc++` and CUDA 13.3 installed on 2026-09-26 (A1 done, v0.5; compile-only here); no GPU on the development machine (R-38).

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
- **Data-movement kernels (D-066):** infrastructure kernels in the C++ driver layer that only move or check data (pack, gather, scatter, checksum) and contain no physics may be written as C++ `ParallelFor` (K1) without a recorded fallback reason; every physics kernel stays K2.
- Cost: FDS's Fortran developers read and review C++ syntax (lambdas, `amrex::Real`, 0-based component index); the kernel no longer matches the FDS Fortran source line by line, so each port is checked against the shimmed Fortran kernel (T1 on frozen input).

### K2 — Fortran kernels with simple OpenMP `target` offload on AMReX device memory (v0.3, owner decision)
- Each kernel stays Fortran. The C++ `MFIter` loop passes box bounds and FAB device addresses through `BIND(C)` (`GPU.rst:1457-1530` documents the pattern); the loop nest carries one `!$omp target teams loop` directive.
- **Allowed subset (v0.5):** `!$omp target teams loop` plus `collapse`, `private` and the device-address clause (v0.6: `target teams distribute parallel do` instead of `teams loop` for a kernel that calls a routine or uses a local array, see the kernel-shape rule) (`is_device_ptr` on nvfortran, `has_device_addr` on gfortran; see "Kernel interface and device data" and "K2 coding rules"); `map` only for small host scalars/constants. Not allowed: nested or combined constructs beyond the one directive, `declare target` on module data, or other complex constructs. The style guide shared with K1 fixes naming (FDS variable names kept).
- **Toolchain:** `nvfortran` for the kernels, with AMReX built for CUDA by the same NVIDIA HPC SDK (`nvc++`, as in `cuda.yml:188-252`). The SDK is free and compiles without a GPU (A-31).
- **Second-compiler check (optional, non-gating; v0.3.7):** kernel code may also be compiled with Intel `ifx` 2026.1.1, which is already on the development machine (not tried). It is a source-level check only and gates nothing. `amdflang` is dropped (AMD out of scope per owner decision 2026-09-25).
- **CUDA Fortran (optional, NVIDIA-only):** allowed only for hot loops where profiling on real NVIDIA hardware shows a meaningful gap against the OpenMP version. Each such kernel keeps its OpenMP version as fallback and correctness reference (T1 against it). None exist yet; none can be justified before GPU hardware (R-38, Q11 (b)).
- Caveats:
  - Offload is not an AMReX backend (`AMReXOptions.cmake:124-125`); AMReX is still built for CUDA and the offload compiler must share its device runtime (R-39). Upstream CI builds AMReX with `nvfortran` but compiles no offload code, so the pairing is untested upstream. AMReX's GNU-make docs list OpenMP offload as supported only with IBM compilers (`GPU.rst:140`), so the docs give no nvfortran precedent.
  - AMReX kernels run asynchronously on AMReX's stream, `target` regions on the OpenMP runtime's queue; every C++/Fortran kernel boundary needs a synchronization or the OpenMP queue wired onto AMReX's stream (R-39). The cost per boundary is unmeasured.
  - `POINT_TO_MESH` module pointers (`mesh.f90:367-376`) are host descriptors and cannot enter `target` regions as they are, so kernels need explicit array arguments: the same rewrite as K1.
  - **Portability:** not a criterion. Only the optional CUDA Fortran kernels are NVIDIA-specific, and NVIDIA is the only GPU target; AMD/HIP pairing questions are out of scope per owner decision 2026-09-25.
- Status (v0.5): `nvfortran` 26.9 is installed; K2 compiles for CUDA and matches at T1 and T0 on the CPU and on an NVIDIA GPU after two source fixes (see "K2 coding rules").

| | K1 restricted C++ | K2 Fortran + OpenMP `target teams loop` |
|---|---|---|
| Language FDS developers read | restricted C++ | Fortran, same OpenMP family FDS already uses |
| NVIDIA toolchain | AMReX CUDA (nvcc or nvc++) | nvfortran + nvc++ (NVIDIA HPC SDK), one toolchain |
| Second compiler | not needed | optional, non-gating `ifx` compile (kernel source only) |
| NVIDIA-only code | none | optional CUDA Fortran hot loops only (OpenMP version kept) |
| AMD | out of scope per owner decision 2026-09-25 (source portable in principle) | out of scope per owner decision 2026-09-25 |
| AMReX upstream support | backend, CI-built, used by PeleLMeX/ERF/incflo | documented pattern; not a backend; not CI-tested |
| Stream/sync | native | sync or stream wiring at every boundary (R-39) |
| Device data | `Array4` captured by value | `is_device_ptr` on explicit-shape dummies (nvfortran; `has_device_addr` on gfortran); no `map` of field data |
| Kernel rewrite (explicit args) | required | required |
| Compiles for CUDA and runs on an NVIDIA GPU (v0.5) | yes; T1 and T0 bitwise | yes after two source fixes; T1 and T0 bitwise |

### Kernel interface and device data for K2 (IR-007)
- **No-copy rule.** AMReX-owned arrays already live in device arena memory, so `map(tofrom:)` on them at kernel entry would copy again or be wrong. Kernel entry wrappers pass AMReX `Array4` data as explicit-shape dummy arrays and declare them `has_device_addr` (OpenMP 5.1) on the `target teams loop` construct. If the first nvfortran compile (A-31) shows `has_device_addr` unsupported for Fortran dummy arrays, fall back to `is_device_ptr` on `type(c_ptr)` arguments with `c_f_pointer` inside the target region. The mechanism is confirmed by that first compile.
- `map` clauses are allowed only for small host scalars/constants, never for AMReX-owned field data.
- **AMReX managed memory** (`amrex.the_arena_is_managed=1`) is rejected for production (page-migration cost; it hides missing-device-data bugs) and allowed only as a debugging aid.
- **Threading rule (reworded for v0.3):** no host-threading OpenMP (`parallel do`) inside device kernels; `target` directives appear only in kernel files. AMReX handles host threading and tiling from the driver. This replaces "no `!$omp` inside device kernels" (IR-007), which would forbid K2's own directives.
- The Spec & Program Lead owns the IR-007/NFR-044 text (currently "OpenACC Fortran", `deviceptr`, "no `!$omp` threading inside the kernel") and has been asked to update it; IR-007 carries the no-copy rule.

### K2 coding rules and GPU build flags (v0.5, from S4b; requirement text owned by the Spec & Program Lead)
These apply to every K2 kernel file. Evidence: `docs/amrex/s4-cuda-mass-findings.md` §5 and §10.
1. **No private copy of a loop-invariant value in an offload loop.** A scalar dummy, or any expression that does not depend on the loop indices, is read directly inside the loop. It is never copied into a local first (`QMIN = QMIN_IN` at the loop top) and never put on the `private` list. `private` is for temporaries that each iteration assigns from cell data or indices. Reason: nvfortran 26.9 handles such a copy wrongly at the teams level and the device threads use a wrong value (every cell clipped, mass created). The CPU paths do not show it.
2. **Parenthesise every sum whose order matters.** A sum that must reproduce the FDS or K1 order (T0) is written with explicit parentheses in that left-to-right order, for example `RHS = ((A + B) + C) + D` and `0.5*((X + Y) - DT*RHS)`. Fortran lets a compiler reassociate an unparenthesised sum, and nvfortran `-O2` does; it honours parentheses. Species loops stay sequential `do N` loops (`-Minfo` reports "run sequentially").
3. **Clause list (widened D-029; rule 7 for the directive).** `!$omp target teams loop collapse(3)` for kernels without callees, `!$omp target teams distribute parallel do collapse(3)` for kernels with callees (rule 7), `private(...)` for loop-local temporaries only, and the device-address clause: `is_device_ptr(...)` on the explicit-shape dummies for nvfortran, `has_device_addr(...)` for gfortran, selected under `#if defined(__NVCOMPILER)`. `map` only for small host scalars. No `c_f_pointer` inside a target region (nvfortran fails to compile it; BIND(C) array dummies already carry the device address). **D-068:** `reduction(max:)` and `reduction(min:)` are added to the clause list (exact and order independent; a kernel using them counts as ported only after a device run). `reduction(+)` is forbidden except through the zone-sum order of D-053.
4. **Compile flags for K2 kernels:** `-O2 -Minline -gpu=ccXX,nofma`, never `-fast`. `-Minline` is kept, but it does not rescue the real kernels with callees: use the `distribute parallel do` form of rule 7 for them.
5. **Enforcement.** Each K2 kernel is accepted only when its device result equals the gfortran host result bitwise on every output array after every kernel (harness `prototypes/s4_cuda_mass/k2_repro`), and its `-Minfo=mp` output has no "implicit private" line for a variable the kernel meant to keep shared or invariant. The device comparison needs the owner's test GPU (R-39); on the development machine the host-fallback T0 run and the `-Minfo` review are what can be checked.
6. The nvfortran defects were reported to NVIDIA (`k2_repro/mini`; the wrong-result one is confirmed as a compiler bug, TPR #39018). The rules stay in force regardless of a fix and are re-checked on each nvfortran upgrade.

7. **Kernel shape (v0.7, from S4d and NVIDIA's reply; `docs/amrex/s4-cuda-mass-findings.md` §12).**
   - A kernel whose body calls a routine or uses a local array is written with `!$omp target teams distribute parallel do collapse(3)`. This is the default. With nvfortran 26.9 it is within 1.3% of K1 on all seven S4d nests and bitwise equal.
   - Optional single-macro form: `!$omp target teams loop collapse(3) bind(teams,parallel)`. It is an NVHPC extension, so the `s4_omp.inc` macro uses it only under `#ifdef __NVCOMPILER` (gfortran rejects the clause).
   - Plain `teams loop` on a body with a callee is **never used**: 10 to 19 times slower, and wrong in 2 of the 7 nests.
   - The workaround "initialise the private scalar before the call" is **not a rule**: it fixed one nest but not `H_RHO_D_DZD`, with or without `-Minline`.
   - A kernel without callees or local arrays may keep `teams loop` (`bind` changes it by under 1.5%).
   - With `distribute parallel do`, every loop-local scalar must be on the `private` list (otherwise it is shared between threads, a race). This does not conflict with rule 1: only temporaries assigned inside the loop are private.
   - NVIDIA's classification: the slowness is by OpenMP design (`teams loop` with a callee may be run by one thread per team); the wrong result is a confirmed compiler bug (TPR #39018). The `S4_LOOP` macro in `s4_omp.inc` selects the directive per kernel file; the AMReX Integration Lead applies the macro form and re-times.

**GPU build flags (requirement, both styles).** T1 on the device needs contraction and fast-math off: `nvcc --fmad=false` (K1), `-gpu=nofma` (K2), `AMReX_CUDA_FASTMATH=OFF` (the AMReX option defaults to ON, `AMReXCUDAOptions.cmake:204`), no `-use_fast_math`, and host `-ffp-contract=off`. With contraction on, T0 is lost in every case and the two-level blob case fails T1 at about 1.6e-3 (the SUPERBEE limiter branches and clip thresholds amplify ulp changes). These flags go into the CMake presets and the build docs. **D-068:** `AMReX_CUDA_FASTMATH=OFF` is pinned in the driver CUDA configure with FORCE, and the configure fails if it is ON; `docs/tools/kernel_lint.py` checks the pin.

### Data-movement kernels in the C++ driver layer (D-066)
The K1 bullet under "Kernel implementation style" states the ruling. This subsection gives the definition, the conditions and the launch facts it relies on (plan `docs/amrex/stage1-gpu-spike-plan.md`, sections 4 and 7.12).

**Definition.** A *data-movement kernel* is a device loop in the C++ driver layer that only copies, gathers, scatters, packs, unpacks, fills or checksums values: no FDS formula, no property lookup, no branch on a physical state, and no arithmetic other than index arithmetic and the bit operations of a checksum. Examples: the flux read-out pack and the interface-face gather and scatter of the flux hooks, the wall-state upload staging and the device checksum of the wall seam, the SOLID mask copy, ghost and periodic-partner copies. A loop that contains any physics term, however small, is a physics kernel and stays K2.

**Conditions.**
1. The kernel is bitwise testable against a host reference: values are moved or bit-combined, never re-computed. A checksum is an integer XOR/rotate fold, compared bit for bit with the host fold.
2. The GPU build flags rule (K2 `nofma`, nvcc `--fmad=false`, fast math off) applies to these translation units as well, so that a neighbouring physics expression is never fused across the hand-off.
3. Launch and ordering follow the launch model below.

**Reason.** The driver layer is already C++, and AMReX owns its streams, arenas and `ParallelFor`. A K2 kernel for a pure copy would add an nvfortran translation unit and a blocking launch for no maintenance gain, and the single-source argument for K2 (one FDS loop, one kernel) does not apply because no FDS loop stands behind a copy.

**Launch model (measured on one cc 8.9 GPU, nvfortran 26.9, nvcc, AMReX CUDA build).**
- A K2 target launch (`target teams distribute parallel do`, no `nowait`) is blocking: a kernel of about 285 ms made the call return after about 285 ms, and the following device synchronisation took about 3.5 us. A K1 launch returns in about 1.6 to 3 us.
- A K2 launch does not wait for work on other streams: launched right after a kernel on another stream with no sync, it read stale data in 10 of 10 runs. The host must synchronise that stream before a K2 launch that reads what an AMReX stream wrote.
- A `cudaMemcpy` placed after a K2 launch is safe (0 of 10 stale, the launch had completed); after a kernel on a non-blocking stream it was stale in 10 of 10 runs.
- `nowait` K2 launches cost about 2 ms each to issue and are unordered against blocking K2 launches (10 of 10 runs saw an unfinished writer without a task wait). They are not usable.
- Blocking K2 launches issued from different host threads overlap on the device (20 small launches: 25.3 ms on 1 thread, 26.5 ms on 2 threads, 41 ms on 3 threads).
- **Design constraint (D-068): "one stream per box" means one host thread per box for K2 launches.** K2 launches block and do not wait for AMReX streams, so K2 cannot be placed on an AMReX stream.

**Consequences.** Data-movement kernels written as K1 can be enqueued on the development machine stream and overlap with other boxes; K2 kernels are host-synchronous points. The ordering rule for a mixed sequence: stream sync before each K2 launch that consumes a stream result, no extra sync after a K2 launch, and no sync for a copy-back (`cudaMemcpy`) after a K2 launch. The 8-byte checksum copy-back of the wall-state test is the only per-stage sync it adds, and only in the first-run and CI check mode.

### GPU kernel data layout (v0.6, owner-approved; requirement text owned by the Spec & Program Lead)
Applies to K1 and K2 alike, and to every kernel added later.
1. **Flat arrays and scalars only.** A GPU kernel receives raw pointers (or explicit-shape dummies) plus integer bounds, taken from AMReX `MultiFab`/`Array4` data, and scalar parameters. It never receives an FDS derived type, and never touches an allocatable or pointer component of one. Nothing inside a kernel calls `ALLOCATE`, reassigns a pointer or follows a pointer chain.
2. **A shim converts at the call site.** The per-box wrapper (the `POINT_TO_BOX` shim of Option C, or the K1 lambda capture) unpacks whatever the FDS routine reads from derived types into the flat arguments before the launch. Kernel bodies do not change when FDS data structures change; only the shim does.
3. **Small types are packed once.** Small scalar or fixed-size-array types (species/reaction constants, material properties, tables) are packed once at setup into flat device arrays (or passed by value as small structs of scalars) and mapped at setup, not per call.
4. **Nested allocatable components are flattened.** Surface, wall, material and similar records whose components are allocatable arrays (for example `WALL(IW)%...` and per-layer arrays) are flattened to structure-of-arrays with an index table (offset and count per record), built at setup and rebuilt after regrid where the records move. The kernel indexes the flat arrays through the table.
5. **Why this is a rule.** The GPU call-graph survey (`docs/inventory/gpu_callgraph_survey.md`, FireX `36975d7`) finds 829 device-eligible loops; **542 of them have a layout-class blocker** in the body (allocatable/pointer or derived-type component access, derived or indirect write indices, module-global writes and similar), and 178 have no body blocker. **About 294 routines (46,268 lines) need GPU-safe rewrites** (upper bound, lexical, not an effort estimate); 31 reachable routines take a derived-type dummy and 37 a POINTER or ALLOCATABLE dummy. The rule is what keeps that work mechanical.
6. **Not a K1/K2 difference.** CUDA Fortran has the same restriction (no automatic deep copy of a derived type with allocatable components), and OpenMP `target` needs the same explicit `map` of every component. K2 does not remove this work; it keeps the Fortran arithmetic in place once the data is flat.
7. **Unified-memory shortcut (first port only).** `nvfortran -gpu=managed` (CUDA managed memory) could let a first port of a routine run without the flattening, at the cost of page-migration time and of hiding missing-device-data bugs. It is allowed only as a temporary bring-up aid for a routine, never in production and never as the acceptance run; this parallels the existing rejection of AMReX managed memory (`amrex.the_arena_is_managed`) for production. It is [VERIFY] whether the shortcut works at all for FDS's nested types on nvfortran 26.9.

### Upstream FireX merge protocol for ported kernels (v0.6 draft; requirement text owned by the Spec & Program Lead)
**Problem.** FireX is merged into the branch daily (D-034). Upstream loops may name mesh data through derived-type components (`M%U(I,J,K)`, `WC%...`, `SF%...`), while a ported kernel (K1 or K2) receives flat arrays (`U(I,J,K)`) under the data-layout rule. A textual merge of such a hunk into a ported file therefore does not apply cleanly, or applies and silently reads the wrong thing. This is true for K1 and K2 alike. The survey (`gpu_callgraph_survey.md` §8) puts the load at about 5 to 6.6 commits a month touching a candidate routine (median, with and without bulk commits), and 64% of the changed loop hunks (443 of 689) are body-only arithmetic, which is the case this protocol makes cheap.
1. **Naming rule.** A ported kernel keeps the upstream variable names inside its body: the dummy argument is called `U`, not `M_U`, wherever upstream writes `M%U` or reaches `U` through a pointer set by `POINT_TO_MESH`; the same for `RHOP`, `ZZP`, `DX`, and so on. Only the call-site shim knows where the data comes from (it extracts `M%U` and passes it as `U`). Loop bounds and index names also stay as upstream. So an arithmetic-only upstream edit pastes across without change. The style guide shared by K1 and K2 already fixes this (FDS variable names kept); this makes it a merge rule. K1 bodies keep the same names for the same reason.
2. **Mechanical rename step.** A script (`tools/port_rename_hunks.py`, to be written by the implementer who ports the first kernel) rewrites an upstream patch hunk so that it applies to the ported file: it drops the derived-type prefixes (`M%`, `SF%`, `WC%`, `B1%`, and so on) listed in a mapping table (`docs/inventory/port_rename_map.csv`: prefix or component, flat name, kernel or "any"), and refuses (exit code non-zero, hunk printed) on any hunk that names a component not in the table. It never guesses.
3. **Kernel-to-upstream mapping table.** `docs/inventory/port_kernel_map.csv`, one row per ported kernel: kernel file, style (K1 or K2), Fortran routine, file, line range of the loop nest at the last synced upstream commit, shim function, and the commit at which the port was last verified. The daily FireX merge runs a check (`tools/port_merge_check.py`) that diffs the incoming commits against the ranges in the table and lists every ported routine touched upstream, split into body-only hunks (heuristic of `gpu_callgraph_survey.md` §8: no `DO`, `!$OMP`, `CALL`, `USE`, `ALLOCATE`, pointer assignment, declaration or routine name in the hunk) and surface-changing hunks. Body-only hunks go through step 2 and are pasted into the kernel by the kernel's owner. Surface-changing hunks are reviewed by hand. The Legacy Mapper's churn data (`gpu_routine_churn.csv`) and `gpu_callee_routines.csv` seed the first version of the table; every newly ported kernel adds its row as part of its acceptance. Untabled ported routines are a merge-check failure.
4. **Edits that add or change a field a kernel reads.** A new mesh field, wall or surface component, or a new scalar read inside a ported nest is a two-part change: a new argument on the kernel, and the matching extraction in the shim (and, for a nested allocatable record, the structure-of-arrays table, data-layout rule point 4). The merge is not complete until both are made. **Gate:** the kernel's bitwise tests (V&V harness: device or shimmed kernel against the FDS reference at T0/T1, K2 also against the gfortran host result) pass on the merged tree before the merge is committed; a merge that touches a ported routine and skips its bitwise test is not accepted. New arguments get the same header-note treatment as the rest of the kernel (passive scalars, dropped cylindrical terms).
5. **Host-only routines keep the derived types.** A routine that is not ported (host-side setup, output, HVAC, particle bookkeeping, routines with I/O or allocation) keeps its FDS derived-type code and merges as ordinary Fortran. If a flat-array copy is ever costly to keep in sync, a host-only routine may instead stay on derived-type components that are pointer-backed views into the same storage as the flat arrays (a pointer component set at setup to the AMReX array). That fallback is for host code only, never in a device kernel, and adds the alias trap that the W2 alias trick already carries (Role 1 plan, P1); use it only where measurements show the shim copy matters.
6. **Cost accounting.** The protocol turns the K1 re-port count (about 8 to 12.5 routine re-ports a month) into a review of the flagged hunks. It does not remove the review: the per-merge check, the rename script and the mapping table are work items (A-tracked in `signoff-gaps.md`) that must exist before the first kernel is accepted into the branch.

### Wall tables and wall-state ownership (v0.7, rulings W1 and W2)
**W1: wall lists per box or per mesh.** One set of wall tables **per mesh** (per level and per rank), in upstream wall order and with the upstream global `IW`, plus per-box integer index lists `WLIST_EXT` and `WLIST_INT`. Wall kernels loop `DO IWI=1,NWL; IW=WLIST(IWI)` with the body unchanged; for one box per mesh the list is the identity.
- Against one table set per box: it needs a local-to-global `IW` map, host routines would see two numberings, and a regrid would rebuild every table.
- Against offset-and-count per box (data-layout rule point 4): it requires sorting walls by box, which breaks the upstream order that the host routines and the merge protocol rely on.
- Cost: one extra integer gather per wall access ([VERIFY] on nvfortran). The kernel generator's wall-loop form gets `WLIST` and `NWL` arguments (open item for its owner, `signoff-gaps.md`).

**W2: ownership of `UVW_SAVE`, `U_GHOST`, `V_GHOST`, `W_GHOST`.** Their producers (`velo.f90` `MATCH_VELOCITY` and `ccib.f90`) stay on the host. Until the geometry refactor (ADR-003) lands, the AMR Data Layout Implementer (Role 1) owns them. The shim uploads the four arrays once per stage (predictor and corrector each run `MATCH_VELOCITY`) before the wall kernels, with a checksum assertion against the host arrays. *(Wording corrected 2026-10-02: "once per step" was imprecise; D-061.)* Role 1 also owns the table refresh after obstruction events (`REASSIGN_WALL_CELLS`) and after `WALL` reallocation. The AMReX Integration Lead reviews the kernel side. Exit condition: the geometry refactor puts these producers on the device and the upload is deleted.

**Launch model for wall kernels (D-068).** Wall kernels written as K2 are launched with one host thread per box, because K2 launches block and do not wait for AMReX streams; the stream sync before a K2 launch that reads staged wall state is described under "Data-movement kernels in the C++ driver layer".

**Host-written tables.** A table written on the host (for example `B1_RHO_D_DZDN_F` after a host `WALL_BC`) is scattered to the device by the shim before the kernels that read it, and gathered back if a host routine reads a device-written table. Staged plan: `drafts/gpu-staging-plan.md`.

### Upstream merge log
- 2026-10-01: FireX `835588bb54` (11 commits) merged as `bee11f0329` on top of `9a176c9df0`. The source delta is `radi.f90` only; no ported kernel or driver file is touched. Not pushed. Implementers base new commits on `bee11f0329` or a descendant.
- 2026-10-02: checked, no new upstream commits.

### Radiation sweep (D-039, FR-062 ruling)
- One per-box sweep kernel per angle subset. It reads box-face intensities only from the previous exchange (double-buffered faces) and writes its own outgoing faces to the other buffer, so the result does not depend on the order or concurrency in which boxes are swept.
- The face exchange is a host-side step between passes, using one mechanism for faces between boxes on the same rank and on different ranks. With `RADIATION_ITERATIONS=K` there are K sweep passes and K exchanges per radiation step.
- Consequence: at a fixed box layout, radiation is independent of rank count, ranks per GPU, box ownership and processing order; only the box split changes results (D-039). The physics-feeding sums (`RAD_Q_SUM`, `KFST4_SUM`) use the FR-005 (ii) exact accumulation, replacing the `!$OMP CRITICAL` sum at `radi.f90:4119-4122`.
- Optional escalation only if the lag-error check fails at K ≤ 3: a radiation BoxArray with larger boxes (built from geometry and a fixed size, never from rank count), then a globally ordered sweep. Neither is built by default.

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
- Pass order per stage and level (v0.3.3; p1-findings §13.3), with two host-side OR reductions (in AMReX 26.09 `ParallelAllReduce::Or` takes only a single value, so the packed flag vectors use `ParallelAllReduce::Max` over 0/1 integers):
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
**Decided 2026-10-02 (owner): K2 is the default, K1 is a per-kernel fallback only.** The text below is the pre-decision reasoning, kept for the record. Before the decision: K1 and K2 were co-equal candidates and the P1 readability review decided (NFR-044, D-027 (c), D-029). The trade-off to weigh is readability for Fortran-only developers and continuity with FDS's existing OpenMP (favour K2) against upstream support and stream simplicity (favour K1). Portability is not a criterion: NVIDIA is the only GPU target and AMD is out of scope per owner decision 2026-09-25.

**Evidence summary for the K1/K2 choice (v0.6; the owner chose K2 on 2026-10-02).**
- *Performance* (`docs/amrex/s4-cuda-mass-findings.md` §11, one cc 8.9 GPU, one process, unlocked clocks, median of 100 full-stage repetitions): the full 12-kernel mass stage takes 14.20 ms (K1) and 14.17 ms (K2) at 128^3, and 46.90 ms (K1) and 46.68 ms (K2) at 192^3, so K2 is within 0.5% of K1. At 64^3 the stage is 1.65 ms and K2 is 3% slower because of the per-region host cost (about 4-5 us per blocking region); removing the host syncs saves 2-4 us per stage. One kernel, one GPU, one stage: not a whole-code result.
- *Merge burden* (`gpu_callgraph_survey.md` §8, churn of the candidate routines): a K1 port must be re-translated for every upstream commit that touches a ported routine, about 8 (excluding bulk commits) to 12.5 (all commits) routine re-ports a month at the median, if FireX is merged daily as decided (D-034). K2 keeps the Fortran source in place, so those commits merge as ordinary Fortran changes.
- *Single source* (`docs/amrex/s4-k2-single-source.md`): one K2 Fortran file builds for CPU (gfortran, nvfortran host, ifx) and for the GPU with a 3-macro include, bitwise equal to the serial reference on the runs listed there. K1 needs a separate C++ copy of every kernel.
- *Against K2:* two failure modes found only on the GPU (findings §5), now covered by the K2 coding rules but invisible to CPU tests; a compiler dependency on nvfortran with OpenMP offload; the readability reviewers split.
- *Was open before the owner decided:* one closing GPU T1 run; the owner's sign-off on S5. The coordinator's recommendation is K2 as the default with K1 as a per-kernel fallback, conditional on both.

### P1 readability review record (NFR-044)
The mass kernel (`MASS_FINITE_DIFFERENCES`, `mass.f90:20`, and `DENSITY`, `mass.f90:365`) is written in K1 and K2 and reviewed by FDS Fortran developers. The K2 variant is OpenMP-offload Fortran (`nvfortran`, `target teams loop`); an OpenACC variant is not built (v0.3). Per D-029 the review happens only after both variants compile for a real GPU backend (nvfortran/nvc++ + CUDA after A-31); the K1 variant goes first. The K2 variant may also get the optional, non-gating `ifx` compile. Both variants also run on the host (K1 CPU build, K2 host fallback) against the shimmed Fortran kernel.

| Item | Value |
|---|---|
| Review status | Closed 2026-10-02. Package `drafts/s5-readability-review-package.md` (2026-09-29); the two reviewer answers were split (Species & Combustion Lead: K1; FDS Legacy Mapper: K2 with CUDA Fortran for profiled hot loops); the owner decided |
| Reviewers | Project owner; AMR Species & Combustion Lead; FDS Legacy Mapper (owner decision D-048, 2026-09-26) |
| Variants and commits | K1 `s4_mass_k1.H`; K2 `s4_mass_k2.F90` + `s4_mass_k2.H`; branch `s4-cuda-mass`, `b51f4361b3` (`0bbb0c7cb5` plus the S4b fixes) |
| Device compile (backend, compiler versions) | K1: nvcc 13.3 (NVHPC 26.9); K2: nvfortran 26.9, sm_80. nvfortran rejects `has_device_addr`; `is_device_ptr` on the explicit-shape dummies works (v0.5) |
| K2 optional `ifx` compile (non-gating) | TBD |
| Host result vs shimmed kernel | K1 and K2 byte-identical to P1/FDS: 60 comparisons, 69 checks at 1, 2, 4 ranks. On an NVIDIA GPU (cc 8.9) both 29/29 at T1 and T0 bitwise, K2 after the two fixes |
| Outcome (K1 / K2) and reasons | **K2 default, K1 per-kernel fallback (owner decision, 2026-10-02): tests show similar GPU performance at a much lower maintenance cost.** Reviewer positions as received (reviewers split): *Species & Combustion Lead: K1.* Reasons: K1 is faster to read, and K2's two hidden rules are invisible to CPU tests (one gave +306 % density on the GPU) so a Fortran team adding kernels for years will eventually break one; toolchain risk (one vendor compiler, two defects in one kernel, untested pairing with AMReX) is close to decisive; results equal and K1 was faster in the one timing sample. Caveat: chemistry (`DERIVATIVE`, `JACOBIAN`, falloff routines in `chem.f90`) is where K2's Fortran reuse has value, but batched SUNDIALS on the device from an offloaded Fortran region is unproven, so it does not flip the pick; the chemistry port sizing and a CVODE-on-device spike follow the style decision. *FDS Legacy Mapper: K2 with CUDA Fortran for profiled hot loops (OpenMP version kept).* Reasons: kernel authors know the FDS source in Fortran and can check a Fortran kernel line by line; K1 adds a third index shift (`return` in a lambda for `CYCLE`). Condition: the hidden rules must be enforced, not only written (a GPU CI check, a lint rejecting `private` items that the loop only reads, a script flagging unparenthesised multi-term sums), and the bounds plumbing must be generated or templated. K1 stays the fallback if the measured sync cost is unacceptable. *Owner:* decided K2 (see above). The reviewer conditions become requirements for K2: the hidden rules are enforced by a GPU CI check, a lint rejecting `private` items that the loop only reads and a script flagging unparenthesised multi-term sums; bounds plumbing is generated or templated; each kernel header states the passive-scalar and dropped cylindrical-term limits. K1 is used for a kernel only when a K2 version cannot be made correct and acceptably fast, recorded per kernel with the reason. Points both reviewers raised, valid for either style: (a) FDS loops over `N_TOTAL_SCALARS`, the ports over `N_TRACKED_SPECIES`; passive scalars must be stated as handled elsewhere; (b) the ports use scalar `RDX/RDY/RDZ` and drop `R(I)`/`RRN(I)` (uniform Cartesian only), so they are not yet line-by-line ports of the cylindrical terms; each kernel header must say so. |
| Owner sign-off (if K2: any CUDA Fortran kernels) | **Signed off 2026-10-02: K2 default, K1 per-kernel fallback.** CUDA Fortran stays as decided on 2026-09-25: only for profiled hot loops, each with an OpenMP fallback; none is approved yet |

## Layering for maintainability
- **(a) Driver, regrid, pressure: C++.** `AmrCore` subclass, time loop, FillPatch/flux registers, regrid and side-data rebuild orchestration, `FFT::Poisson`/MLMG, particles container, checkpoint. Small and stable once written; needs real C++ skill, owned by the AMReX-side team.
- **(b) Physics kernels: restricted style (K1 or K2).** Where FDS developers work day to day. One kernel = one `ParallelFor` body or one Fortran loop nest under `!$omp target teams loop`, with explicit arguments.
- **(c) Input parsing, setup and output: may stay Fortran on the host.** `read.f90`, setup in `init.f90`, `dump.f90`/`vtkf.f90`. This matches NFR-043's I/O exception and keeps the largest FDS-specific code unchanged. Regrid-time side-data rebuild (`WALL`, `CELL_INDEX`, `EXTERNAL_WALL`, wall records): host-allowed until Phase 11 and device-capable by Phase 11 (owner decision D-047, NFR-043, v0.4).
- **(d) Scope of C++ (D-043):** layer (a) only, plus any K1 fallback kernels (K2 was chosen 2026-10-02, so K1 is the exception, one kernel at a time, with the reason recorded). No other FDS physics moves to C++.

## Pressure path and global reductions (room decisions, 2026-09-25; details in ADR-002 v0.2)
- **Solver choice, per-step selection, shared gauge and eps_H agreement:** see ADR-002 v0.2, "Accepted decisions / Pressure solver" (`amrex::FFT::Poisson` replaces porting `pois.f90`, D-021; FR-037, FR-039, D-012).
- **MLMG order:** `setMaxOrder(2)`, confirmed by P2; the FR-039 fallback to order 3 is not needed. Evidence and the other P2 rulings (mean removal, `average_down_faces`): ADR-002 v0.2. **v0.4:** HYPRE is no longer bottom-solver only; the pressure solve sits behind a solver-agnostic interface with MLMG and an assembled-matrix HYPRE PCG + BoomerAMG backend, default set by comparison A-56 (ADR-002 v1.1, D-046).
- **Refinement ratio ≤ 4, supported {2,4}:** see ADR-002 v0.2, "Refinement ratios" (FR-010, D-030).
- **Global scalar reductions:** exact fixed-point sums from per-box sums, computed domain-wide; see ADR-002 v0.2, "Decomposition requirements" (FR-005 (ii), (v); D-028; R-36). On GPU builds the per-box accumulation runs on the device.
- **On the GPU (NFR-043):** `FFT::Poisson` via cuFFT for single-level uniform runs; composite MLMG on the device otherwise. `PoissonHybrid` (z-stretched only; device branch at `AMReX_FFT_Poisson.H:715-780`) is no longer needed: stretched-grid cases stay FDS-only (owner decision 2026-09-25, D-030), so R-29/A-29 are moot. HYPRE is CPU-only in our builds, including the FireX-pinned v2.32.0-24 `63331f19c` (`HYPRE_USING_CUDA`/`HYPRE_USING_GPU` undefined, `(local GNU third-party library tree)/libs/hypre/63331f19/include/HYPRE_config.h:93, 144`; D-026), so MLMG's native bottom solver is the device default (a GPU HYPRE would need its own build per backend; R-30). v0.4: with the assembled-matrix backend of ADR-002 v1.1, a device HYPRE build becomes a GPU-phase work item. Stretched meshes need no GPU pressure path: AMR mode uses uniform grids on each level and stretched cases stay FDS-only, which closes charter Q11 (c) (ADR-002 v0.2.1).

## Recommendation
**Adopt Option A (C++ AmrCore driver), owner-confirmed via D-027. Migrate through the Option C `POINT_TO_BOX` shim as a CPU-only stage with the three scope decisions above, then extract every kernel onto the device in style K2 (owner decision 2026-10-02), with K1 only as a per-kernel fallback. Reject B definitively.**

The driver decision rests on four points, the first now decisive on its own:
1. GPU: the full time step must run on the device (D-027, NFR-043), and F_Interfaces has no GPU support (`GPU.rst:105-106`).
2. Particles: the fixed `amrex_particle` struct cannot carry FDS particle state.
3. Maintenance: the F_Interfaces layer is maintained but not developed, so every new AMReX feature would need our bindings.
4. ADR-003's EB candidate for complex geometry needs C++ (v0.4: EB is one of two candidates, ADR-003 v1.1).

The pressure path is not among them: F_Interfaces can drive a cell-centred variable-β projection, though `FFT::Poisson` would have needed a C++ wrapper (FR-030).

Confidence: **high** for rejecting B (owner requirement plus upstream documentation); **medium** for the shim's edit count and regrid cost, pending S1/S3; kernel style **decided 2026-10-02: K2 default, K1 per-kernel fallback**. What would change it: nothing short of withdrawing D-027 reopens B; the K1/K2 choice turns on the review, on the first nvfortran compile (device-data mechanism, R-39).

## Rejected alternatives and why
- **K1 as the default kernel style** (owner decision, 2026-10-02): restricted C++ `ParallelFor` kernels for all physics. Rejected because GPU performance is similar to K2 (within 0.5% on the mass stage) while a separate C++ copy of every kernel must be re-translated for each upstream change and is harder for the Fortran developers to maintain. K1 remains the allowed per-kernel fallback.
- **B:** rejected definitively (D-027): no GPU support in F_Interfaces. It also loses on particles, maintenance and EB, and wins only on single language. The kernel refactor is unchanged.
- **D:** static mesh topology; no AMR time interpolation or reflux.
- **E:** V&V and schedule risk; D-027's rewrite happens incrementally instead.
- **Shim as end state:** incompatible with the device by construction (D-027 (b)).
- **OpenACC for K2** (v0.2 candidate; replaced by owner decision, 2026-09-25): the most mature Fortran offload model on nvfortran, and for simple loops its performance is expected to be similar to OpenMP `loop`. Rejected because it is NVIDIA-only in practice, adds a second directive dialect alongside FDS's existing OpenMP, and is less portable. No performance numbers have been measured for FDS kernels in either model.
- **Full OpenMP offload feature set** (complex/nested constructs, `declare target` on module data): harder to read for FDS developers and more compiler variance; K2 is limited to one `target` loop directive per kernel, `teams loop` or (v0.6) `teams distribute parallel do`, plus `collapse`, `private` and the device-address clause.
- **Second `FillBoundary` in D-031 clipping** (between density apply and species clip): one extra exchange plus a device sync per call for a bitwise-identical result (A-37); replaced by the redundant density clip over valid+1; kept only as a prototype cross-check.
- **Single OR over all clip flags in D-031:** wrong order; the density apply must be gated by the reduced density flags before the species gathers read the clipped `RHOP`, so two reductions are needed.
- **`map` of AMReX field data / AMReX managed memory in production:** `map` copies data already in the device arena (breaks IR-005's no-per-step-copy rule); managed memory costs page migration and hides missing-device-data bugs. Managed memory stays a debugging aid only.

## Consequences and risks
- **+** AMReX features without binding debt; ~860 lines of exchange and most MPI calls retire; one route to the device for driver, pressure and kernels.
- **−** Two languages; transitional shim code; every time-step kernel is rewritten (Phase 11 / M11; schedule TBD until S4).
- **R-26 (shim becomes permanent / too costly):** now a timing rule only; extraction itself is required.
- **R-38 (no GPU on the development machine):** NFR-043 is verified here by a real GPU-backend compile and host-fallback runs only (D-029); device correctness and the no-transfer rule wait for hardware.
- **R-39 (offload toolchain shares AMReX's device runtime, unproven):** applies to K2 only (nvfortran OpenMP offload, and any CUDA Fortran kernels); first tested in S4b (v0.5): `is_device_ptr` shares AMReX's device memory and the two K2 defects are recorded in "K2 coding rules"; stream ordering, sync cost and a profile remain open.
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
  - ~~Install the NVIDIA HPC SDK (A-31) so both P1 variants compile for CUDA.~~ Done (A1, S4a, S4b; v0.5).
  - For K2: document how an `nvfortran` OpenMP `target` region is ordered against AMReX's CUDA stream, and measure the per-boundary cost (R-39). In the first compile (A-31), confirm `has_device_addr` on explicit-shape Fortran dummies, or switch to the `is_device_ptr` + `c_f_pointer` fallback.
  - Confirm pointer re-binding rules after `FillBoundary`/regrid.
  - D-031: build the §13.5 two-phase clip (terms once, then gather; density terms on valid+2), rerun the 47 `check_clip.sh` checks, and add the covered-coarse-cell check (only a covered coarse cell clips; uncovered cells not renormalised). A-37 is closed for the unoptimised form.
  - Confirm that the per-box fixed-point accumulation (FR-005 (ii)) can run in device code.
  - F_Interfaces maintenance status is moot for the driver (D-027).
- **AMR Spec & Program Lead:**
  - Update charter Q11, NFR-043, D-029, R-39 and A-30 for the NVIDIA-only target; (v0.3.7) AMD is out of scope, not deferred (A-32 moot); Q11 (c) closed (stretched cases FDS-only); NFR-045/A-36 become optional and non-gating, `ifx` only; R-29/A-29 moot. Q11 (b) stays open.
  - (Asked, v0.3) Update IR-007, NFR-044, R-38 and R-39 from OpenACC to OpenMP-offload Fortran: reworded threading rule (no host-threading `parallel do` inside device kernels; `target` directives only in kernel files), the no-copy device-data rule (`has_device_addr`; `map` only for host scalars), the second-compiler portability check, and CUDA Fortran only for profiled hot loops with OpenMP fallback.
  - Regrid-time side-data rebuild under NFR-043: *answered* (owner, v0.4): host-allowed until Phase 11, device-capable by Phase 11; logged as D-047 in NFR-043; S5 reviewers logged as D-048 (spec v0.4.28).
  - Minimum feature set for the Phase 2 demo. Smokeview: *answered* (ADR-004, FR-072). Stretched grids: *answered* (owner, 2026-09-25): not in AMR mode; stretched cases stay FDS-only.
- **AMR Pressure Solver Lead:** confirm the pressure path sits wholly outside the shim. *Done* (P2, ADR-002 v0.2): maxorder 2 vs 3 study and FFT/MLMG eps_H check (FR-039).
- **AMR V&V Lead:** S1 tolerance class T1 vs FireX baseline, with `OMP_NUM_THREADS=1`; T1 check of each K1/K2 port against its shimmed kernel.
- **GNU/Intel build chiefs:**
  - Mixed C++/Fortran/MPI link on gfortran 14.2 + OpenMPI 5 and oneAPI (NFR-020); add an NVIDIA HPC SDK + CUDA configuration (compile-only here).
  - Optional, non-gating `ifx` compile of K2 kernel files (amdflang dropped).
  - `C_F_POINTER` + bounds-remap idioms on all compilers.
  - AMReX OpenMP vs FDS `-fopenmp`.
  - HYPRE per D-026: FireX pins v2.32.0-24 `63331f19c` (FireX CMake mislabels it "3.0.0"), and AMReX is built against it at `(local AMReX install built against HYPRE 2.32)` (R-30). HYPRE is CPU-only today; v0.4: it is also the assembled-matrix pressure backend (ADR-002 v1.1), so a CUDA HYPRE build is needed for the GPU phase.

## Spike plan
| Spike | Scope | Pass | Overturns if |
|---|---|---|---|
| **S1 Shim feasibility** (2 wk) | Minimal C++ AmrCore (from `Advection_AmrCore`), single level, periodic, no walls. Shim binds `RHO, ZZ, U, V, W, TMP`; runs `DENSITY` and `COMPUTE_VISCOSITY` built from the FireX sources (outside the read-only tree). | T1 vs FireX baseline; edit count per kernel recorded; overhead ≤ 10% | pervasive kernel edits needed, or overhead > 30% (then extract kernels directly, skipping the shim) |
| ~~S2 F_Interfaces counter-spike~~ | Cancelled by D-027: its reopen condition required GPU out of scope. | — | — |
| **S3 Regrid wall-state rebuild** (2 wk, with ADR-003) | Rebuild `WALL`/`BOUNDARY_ONE_D` for a box with one OBST after a synthetic regrid | conserved exactly; rebuild < 10% of step time at a 10-step interval (R-26) | full re-init needed per regrid |
| **S4 GPU cost probe = P1 readability variants** (1–2 wk; v0.4: S4a compile-only gates acceptance, S4b GPU run does not, see "v0.4 changes") | Port the S1 mass kernel to K1 (first) and K2 (OpenMP-offload Fortran, `target teams loop`; no OpenACC variant); compile both for CUDA (D-029, after A-31); K2 optionally also with `ifx` (non-gating); run K1 CPU build and K2 host fallback | both match the shimmed kernel at T1; effort per variant recorded (Phase 11 estimate); sync points per step counted for K2; device-data mechanism (`has_device_addr` or fallback) recorded | K2 cannot share AMReX's device memory (neither `has_device_addr` nor `is_device_ptr`) or stream with nvfortran (drops K2; R-39) |
| **S5 P1 readability review** (NFR-044) | FDS Fortran developers review both variants (K2 = OpenMP-offload Fortran); v0.5: unblocked, package written | outcome recorded in "P1 readability review record" above | — (decides K1 vs K2) |

## Decision needed from the project owner
1. ~~Q8: are C++ components acceptable in the code base?~~ **Decided (D-043):** C++ for the driver and AMReX glue only; physics stays Fortran unless the NFR-044 review picks K1.
2. **Q2** is answered by D-027. **Q11 (a)** answered: NVIDIA is the only GPU target; AMD out of scope (owner, 2026-09-25). **Q11 (c)** closed: AMR mode uses uniform grids on each level and stretched cases stay FDS-only (owner, 2026-09-25). Still open: **Q11 (b)** GPU test hardware.
3. ~~AMD: deferred or dropped?~~ **Decided** (owner, 2026-09-25): AMD is out of scope; the second-compiler check is optional and non-gating (`ifx` only).
4. ~~Q5: refined-level output format?~~ **Decided:** Smokeview format and VTK (FR-072), design in ADR-004 (accepted).
5. ~~Regrid-time side-data rebuild under NFR-043?~~ **Decided** (D-047, owner, 2026-09-26): host-allowed until Phase 11, device-capable by Phase 11 (v0.4).
6. ~~P1 readability review (S5)~~ **Decided 2026-10-02 (owner): K2 default, K1 per-kernel fallback.** Earlier text: reviewers **named** (D-048, owner, 2026-09-26): the project owner, the AMR Species & Combustion Lead and the FDS Legacy Mapper. v0.5: S4a is done and the package is written. Still open: the three answers and the owner's sign-off of the K1/K2 outcome.

## References
FireX `Source/{mesh,init,mass,main,read,dump,vtkf,pres,velo,wall,turb,vege,type,ccib,pois}.f90`, `CMakeLists.txt`, `Build/makefile` @ `36975d765f`. AMReX `Src/F_Interfaces/**`, `Src/FFT/AMReX_FFT_Poisson.H`, `Src/LinearSolvers/MLMG/{AMReX_MLLinOp,AMReX_MLCellLinOp}.H`, `Src/Boundary/AMReX_InterpBndryData.H`, `Docs/sphinx_documentation/source/{Fortran_Chapter,GPU}.rst`, `Tools/CMake/AMReXOptions.cmake`, `Tests/CMakeLists.txt`, `.github/workflows/{cuda,hip}.yml`, `CHANGES.md` @ `99ddfda`. HYPRE v2.32.0-24 `63331f19c` `HYPRE_config.h` (`(local GNU third-party library tree)/libs/hypre/63331f19`; D-026). ERF `Source/Microphysics/Morrison/*`; PeleLMeX, incflo, ERF `Source/`. Tutorials `FortranInterface/Advection_F`, `Amr/Advection_AmrCore`. Teammate docs: `docs/{risks,requirements,roadmap,charter,README,spec-responses}.md`, `docs/inventory/*`, `docs/amrex/{mapping,driver-options}.md`, `docs/pressure/0{0,1}-*.md`.
