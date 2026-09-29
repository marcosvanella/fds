# Wave-2 hiring proposal (draft for owner approval)

Status: draft, 2026-09-29. Prepared by the Spec & Program Lead against spec v0.4.28. Nobody is hired by this document. Effort figures are rough (low confidence): they are derived from the roadmap phase estimates and the prototype record, not from measured implementer throughput. Cost is given as effort and usage load; I have no dollar or token model.

## Principles
- Release in small batches, each gated on a checkable milestone, so usage stays bounded and a failed batch is cheap.
- Every implementer works from the spec (requirements, decisions, ADRs) and reports to its lead: the Chief Architect for data layout, the Pressure Solver Lead for the pressure interface. Implementers write code in the reference tree only under the owner-approved commit rules (D-005, D-034, D-037); no pushes.
- Each role has an acceptance test that already exists in the spec, so "done" is not a judgment call.

## Batch 1 (release now): core data layout and pressure backend interface

### Role 1: Core data layout implementer (Phase 2, milestone M2a then M2)
- **Scope:** `USE_AMREX` CMake option; AMReX-owned level-0 `BoxArray`, `DistributionMapping` and field `MultiFab`s; the per-box `POINT_TO_BOX` shim over unmodified kernels; ghost fill replacing the same-resolution part of `MESH_EXCHANGE`; face-index round trip; non-box level 0 with the static gap mask (IR-002).
- **First deliverables:** (1) level-0 driver and shim reproducing `shunn3_32` and `csmag_32` (kernels bitwise on frozen input, full steps within T2); (2) `shunn3_4mesh_32` against its single-mesh copy (M2a); (3) round-trip and race tests (IR-005, IR-007); (4) `stairwell` and `hallways` setup with the gap mask.
- **Dependencies:** ADR-001 (driver: C++ driver, physics in Fortran, NFR-049); P1 CPU prototype (`amrex/p1-findings.md`); baselines and the A-46/A-47 V&V runs; D-031 ghost depth. Pressure calls can be stubbed with the single-level FFT until Role 2 delivers.
- **Docs and decisions needed:** FR-001, FR-002, FR-005, FR-038, IR-001..007, NFR-010..022, NFR-031, NFR-049; D-021, D-022, D-028, D-031, D-034, D-037, D-043; `inventory/mesh_fields.csv`, `amrex/mapping.md`, `amrex/driver-options.md`; ADR-001 items 5 (D-047) and 6 (D-048).
- **Effort:** about 6-10 implementer-weeks to M2a, another 4-6 to the M2 exit (Phase 2 estimate 4-8 weeks elapsed). Heaviest usage of the batch.

### Role 2: Pressure backend interface implementer (Phase 2 deliverables, then Phase 4)
- **Scope:** the solver-agnostic interface of D-046 (a). The interface defines the discrete operator (7-point plus the P2/B5-validated composite coarse/fine discretisation) and backend fluxes are never used. Backends: AMReX MLMG, and the assembled-matrix HYPRE branch (PCG + BoomerAMG). Common layer: mean removal, pinning and gauge per zone and connected component, the true-residual check, the fine-flux overwrite, per-step selection between FFT, masked MLMG and composite MLMG (FR-039).
- **First deliverables:** (1) interface plus FFT and single-level MLMG backend on a standalone harness, checked against frozen solves; (2) the assembled HYPRE backend for single level and masked domains; (3) CI check that both backends agree within eps_H on frozen solves (A-56 [VERIFY]); (4) masked branch on `hallways` and `stairwell` (A-54, A-55); then in Phase 4 the composite multi-level path.
- **Default backend:** the head-to-head (`pressure/05-backend-head-to-head.md`) recommends MLMG with the FDS BoomerAMG bottom as the default and the assembled HYPRE backend as a selectable alternate, because MLMG's setup at each regrid is 8 to 90 times cheaper. That result is CPU and single level, and the `stairwell` case with dropped all-gap boxes has not run. The default stays an owner decision (D-046).
- **Dependencies:** none for the standalone part; Role 1's `MultiFab` layout for integration; ADR-002 v1.1; the Pressure Lead's A-55 composite masked check before the Phase 4 composite path.
- **Docs and decisions needed:** FR-030..034, FR-036..039; D-021, D-024, D-026, D-028, D-032 (as amended), D-046; NFR-030, NFR-043; `pressure/01`, `04`, `05`; ADR-002 v1.1; R-47, R-62.
- **Effort:** about 5-8 implementer-weeks for the single-level and masked parts, then 4-8 for composite (Phase 4 estimate). Can run alongside Role 1.

**Batch 1 gate for batch 2:** M2a passes (Role 1) and the interface agreement check is green in CI on single-level and masked solves (Role 2).

## Later batches (do not release yet)

- **Batch 2, after M2a: multi-level infrastructure and species (Phase 3).** Role 3, regrid and multi-level transport implementer (tagging, conservative regrid transfer, refluxing; FR-010..013, FR-024, IR-003, IR-008; needs Role 1). Role 4, species and combustion implementer (species transport on levels, then combustion; with the Species & Combustion Lead). Effort about 6-10 implementer-weeks each. The composite pressure path (Role 2, Phase 4) follows.
- **Batch 3, after M4 (composite pressure) and the FR-041b spike G2a: solids (Phase 5).** Role 5, walls and solid-phase implementer (per-level wall records, FR-040 snap rule, FR-041a/b, FR-045..047; with the Solid Phase Lead). About 6-10 implementer-weeks.
- **Batch 4, after M5: radiation (Phase 8).** Role 6, radiation implementer (QR on all levels, FR-060..063 with the S-B coupling of D-039; with the Radiation Lead). About 4-8 implementer-weeks. Particles, output and restart (Phases 7 and 9) and GPU porting (Phase 11) are separate roles, proposed when M5 is reached.

## Suggested order and which held implementers to release
1. Release two implementers now: one to Role 1 and one to Role 2. Roles 3-6 stay held.
2. Release the batch-2 pair when the batch gate above is met, then the solid implementer at M4, then radiation at M5.
3. **Which of the 8 held implementers:** I do not have their list or their skills, so I cannot name them. Please have the CEO map the held implementers against these roles. The selection rule I suggest: Role 1 needs strong C++ and AMReX (`MultiFab`, `AmrCore`) plus Fortran `bind(C)` interop; Role 2 needs numerical linear algebra and HYPRE/MLMG experience. Hold implementers whose strength is radiation, solid or species for batches 2-4.

## Risks to watch
- Two implementers editing the same reference tree: assign disjoint directories, and have the Chief Architect own the first commit (D-005, D-037).
- Role 2's HYPRE backend needs its own composite assembly, which is the least-known effort (R-62, A-56); the estimate is widest there.
- Daily FireX merges start with development (D-034, R-44); the V&V Lead re-pins baselines at each milestone (A-42).
- GPU hardware is a single test-machine GPU (R-38, R-45), so GPU acceptance stays with Phase 11 and is not a batch-1 gate.
