# Architecture sign-off gap list

Owner: Chief Architect. Scope: what still stands between ADR-001 (driver architecture) and sign-off, plus the items that block the rest of the architecture baseline (milestone M1: ADR-001/002/003 accepted). Needs column: **Owner** = needs a project-owner answer; **Spend** = needs builds, runs or installs; **Doc** = writing only.

ADR-001's driver choice (C++ AmrCore driver) is already owner-confirmed (D-027). The only open decision inside ADR-001 is the kernel style, K1 (restricted C++) or K2 (Fortran with OpenMP `target` offload).

Status as of 2026-09-26: A4, A5, B1, B2, B3, B5 (analysis), B6, C1 and C2 are done.

## A. Blocks ADR-001 sign-off

| # | Item | Owner | Needs |
|---|---|---|---|
| A1 | Install the NVIDIA HPC SDK (A-31) so both kernel variants compile for CUDA. Compile-only; no GPU needed. | AMReX Integration Lead, with a build chief | Spend |
| A2 | Spike S4a (compile-only; S4b, the GPU run, does not gate): write the mass kernel as K1 and K2, compile both for CUDA, match the shimmed kernel at T1, record effort and the device-data mechanism (`has_device_addr` or fallback). Depends on A1. | AMReX Integration Lead | Spend |
| A3 | Spike S5, the NFR-044 readability review of both variants by FDS Fortran developers; decides K1 or K2 (D-043: physics stays Fortran unless the review picks K1). Reviewers named by the owner (D-048): the project owner, the AMR Species & Combustion Lead and the FDS Legacy Mapper. The owner signs off the outcome. Depends on A2. | Owner, with the Chief Architect recording | Owner |
| A4 | **Done (owner decision D-047, ADR-001 v0.4, NFR-043): host-allowed until Phase 11, device-capable by Phase 11.** Rule whether regrid-time side-data rebuild (wall records, particle bookkeeping) counts as part of the time step that must run on the GPU (NFR-043, D-027). Recommendation to be drafted: host-allowed through the CPU phases, device-capable by Phase 11. | Chief Architect drafts; owner confirms | Owner (yes/no), Doc |
| A5 | **Done (ADR-001 v0.3.8).** Update ADR-001 text: Q8 item aligned with D-043; D-039 and the FR-062 ruling (independent per-box sweep, host-side exchange) in the kernel-interface section; stale "Decision needed" list cleared (Q5, Q8 answered); ADR-004 cross-reference. | Chief Architect | Doc |

Not required for ADR-001 sign-off: GPU test hardware (Q11 (b)) and device timing (R-39). The kernel style is chosen on readability; device performance is measured during implementation.

## B. Blocks the architecture baseline (ADR-002, ADR-003), not ADR-001

| # | Item | Owner | Needs |
|---|---|---|---|
| B1 | **Done (ADR-003 v0.3).** Fold the FR-040 R3 / FR-034 rulings and the FR-041b ruling into ADR-003; clear the P-a/P-b pending marks; fold in `obst_wall_cface_indexing.csv` (now delivered). | Chief Architect | Doc |
| B2 | **Done (`drafts/ruling-nonbox-level0.md`; ADR-003 v0.3, ADR-002 v0.2.2).** Rule on non-box level 0 (gaps between level-0 meshes become solid level-0 cells). `HVAC_leak_exponent` and about 40 FR-006 cases depend on it. | Chief Architect | Doc |
| B3 | **Done: passes with a solver-configuration change (`drafts/ruling-nonbox-level0.md` §5.1, §5.5).** One pinned cell and one D-032 mean per disconnected part is FDS's own rule (`pres.f90:3369-3406, 3447`). On single-level split cases the true residual including pinned rows is ≤ 9e-9·max\|RHS\| (measured).<br>Required configuration:<br>- a > 0 in masked cells (otherwise the default BiCGStab bottom returns NaN);<br>- HYPRE bottom on level 0 (one pin stops geometric coarsening);<br>- no `setNSolve`;<br>- mean removal per component, not per zone.<br>N5 (zones spanning gap components) is checked in §5.3: none in FDS; setup assertion added. Masked-branch convergence against FDS is in §5.4: one solve, no inter-mesh error; FDS `hallways` averages 6.0 iterations per half-step. The composite multi-level case stays on the E-2 test list. Original item: [VERIFY] one pinned cell per disconnected fluid part under E-2 masking. | AMR Pressure Solver Lead | Doc (code check) |
| B4 | **Done: G2a passes (`solid/03-g2a-cost-check.md`).** R=2 bound on `couch` is 0.06 (1-D solve isolated) to 0.16 (whole WALL timer) against the 25 % gate. Memory at finest-ever depth is 0.1–0.3 GB at R=2 and 0.4–1.3 GB at R=4 against 16 GB. No record-depth cap needed. Caveat: the small 2D charring case reaches 0.40. FR-041b wording deltas sent to the Spec & Program Lead. Original item: G2a cost check for FR-041b (WALL share of step time from stock runs; gate 25 % at R=2; memory bound against 16 GB). If it fails, the Chief Architect rules on a record-depth cap. Solid Lead's FR-041b wording to the Spec & Program Lead. | AMR Solid Phase Lead | Spend (small) |
| B5 | **V&V analysis done (`vv/b5_uglmat_p2_closure.md`, runs in `vv-runs/B5/`); P2 criterion met.** UGLMAT reference A-38 (`ns2d_16_int_1to2_refinement`, N=16/32/64). Like-for-like single solve (FDS UGLMAT-HYPRE step-1 right-hand side, rebuilt from its step-1 `H`, solved by the P2 prototype with the same 2:1 patch; A-38 norm, `H` mean removed, whole domain): at maxorder 2 the prototype-to-UGLMAT `H` error ratio is 0.92 / 0.83 / 0.68, order 2.01 / 2.03 (UGLMAT 1.87 / 1.75), so REC-I2 criterion 2 (≤ 1.1×) passes; criterion 3 difference from FDS 5.4e-3 / 2.5e-3 / 1.2e-3 (≤ 3e-2). Caveat: the step-1 right-hand side is rebuilt from FDS's UGLMAT stencil, not dumped. D-040 is not a P2 check; it needs the Phase-4 time-evolved run. One-cell-thick direction: `ref_ratio_vect = 2 1 2`, blocking factor 1 in y, `LPInfo::setHiddenDirection(1)` gives the true 2-D answer (≤ 7e-13 relative), but MLMG needs 122 / 311 / 546 iterations against 9 in 2-D and fails at `maxiter` 200 for N ≥ 32; see B8. | AMR V&V Lead | Done (analysis) |
| B8 | One-cell-thick direction convergence under `setHiddenDirection` (from B5): find the cause or a mitigation (8 smoothing passes give 99 / 128 iterations; Jacobi and HYPRE bottom do not help; 2 coarse y-cells without a hidden direction converge in 8 iterations at twice the y cells). Until then the one-cell-thick setup stays provisional for production. Does not block ADR-002 if ADR-002 records it as an open implementation item. | AMReX Integration Lead, with the AMR Pressure Solver Lead | Doc (code check), Spend (small) |
| B6 | **Done: no constraint needed (`drafts/ruling-nonbox-level0.md` §5.2).** The spacing is span/5000 (`read.f90:10503, 10597`), sub-mm and independent of every dz. The lookup takes the nearest entry (`func.f90:852-854`), and FDS never takes differences of `P_0`. AMR builds the table once from level-0 extents and evaluates it at each level's cell centres. Gap components share it, and ΔP_zone is one scalar. Original item: [VERIFY] `P_0` ramp spacing. | AMR Pressure Solver Lead | Doc (code check) |
| B7 | **Done: ADR-003 v1.0 accepted, then ADR-002 v1.0 accepted (2026-09-26).** Option C (single global dt, subcycling-ready data) recorded as decided. B8 and the GPU masked branch are recorded in ADR-002 as open implementation items. | Chief Architect; owner informed | Done |

## C. Spec items and owner questions

| # | Item | Owner | Needs |
|---|---|---|---|
| C1 | **Done (spec v0.4.24).** Apply the IR-008 ruling (refinement confined to the region) and the FR-062 wording (`drafts/rulings-IR008-FR062.md`), after the Radiation Lead confirms the wording. | Spec & Program Lead | Doc |
| C2 | **Done (spec v0.4.21/v0.4.22).** Record the D-044 amendment (documentation-only snapshots allowed during the spec phase) in the decision log. | Spec & Program Lead | Doc |
| C3 | Q1: confirm the first gate as kernel-level T0/T1 on explicit kernels with frozen input plus T2 for whole runs (full-run bitwise is infeasible, D-022). Blocks M1's V&V baseline, not ADR-001. | Owner | Owner |
| C4 | Q11 (b): date for the NVIDIA test machine. Blocks device verification (NFR-043), not sign-off. | Owner | Owner |
| C5 | Q7: effort or calendar budget. Needed for the roadmap, not for sign-off. | Owner | Owner |

## D. Runs during implementation (explicitly not blocking sign-off)

- S1 shim feasibility. Its fallback (extract kernels directly, skipping the shim) needs no ADR change.
- S3 regrid wall-state rebuild cost (R-26), with ADR-003's accepted freeze rule (D-010).
- ADR-004 spikes S-A (byte identity), S-B (Fortran I/O on a background thread) and S-C (Smokeview load), under A-53.
- D-031 two-phase clip optimisation and its 47-check rerun.
- FR-063 radiation state across regrids (Phase 8).
- R-39 device cost per kernel boundary, on the NVIDIA machine.
- Non-box level 0 (`drafts/ruling-nonbox-level0.md` §6): the padding [VERIFY] (V&V script check) and the mesh-union fallback [VERIFY] items (i)-(iv) in the Phase 2 prototype (Integration Lead, with the Pressure Solver Lead).
- Masked level-0 branch on the GPU path (HYPRE builds are CPU-only), before the GPU phase.

## Critical path

A1, then A2, then A3 is the only chain that needs spend and owner time; everything else in section A is writing. After A3 the Chief Architect sets ADR-001 to Accepted (writing only). Sections B and C can proceed in parallel.
