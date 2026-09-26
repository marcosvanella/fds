# docs/pressure — AMR pressure/projection design (AMR Pressure Solver Lead)

**Base:** FDS FireX `AMReX` branch @ `36975d765f`, read-only at this repository. All FDS `file:line` citations in this folder
refer to that tree. The earlier draft citations against master `ce1f659` were re-mapped and text-checked (see 00 header).
**Status:** Draft v1, 2026-09-25. Nothing here has been executed with FDS or AMReX (no MPI runtime on the development machine). Measured items are
explicitly labelled.

## Assumptions (for this index)
* The files below are a first spec for review. The programme documents (`docs/requirements.md`, `risks.md`, `roadmap.md`,
  `README.md`, `adr/`, `inventory/`) were read and cross-referenced but **not edited**.
* Phase numbering follows `docs/roadmap.md`: Phase 4 / M4 = composite pressure; Phase 6 = subcycling and sync projections.

## Files
| File | Content |
|---|---|
| `00-fds-pressure-baseline.md` | How FDS does it today: H formulation, D and zones, time step, FFT/ULMAT/GLMAT/UGLMAT, IBM and pressure iteration, HYPRE hooks, tunnel preconditioner, mesh exchanges, data flow, properties to preserve. §11 summarises **what FireX changed vs `ce1f659`** (HYPRE AMG settings, GPU offload, resource sets, `HYPRE_DEVICE_RUN`, moved call sites). |
| `01-amr-mapping-spec.md` | The mapping onto AMReX. (A) exact MAC-projection equivalence and operator choice, including C++ vs Fortran availability; (B) time stepping and subcycling, including what deferring sync projections (R-05) really defers; (C) C/F handling, zone integrals, compatibility, gauge; (D) fate of the pressure iteration; (E) solids; (F) BCs; (G) HYPRE role/versions and the answer to ADR-001; (H) verification plan; (I) **hard acceptance numbers for H vs GLMAT/UGLMAT** (FR-016(c)); (R) verified references. |
| `02-performance-expectations.md` | Cost and memory models for FFT vs MLMG vs HYPRE; measured FFT/stencil/bandwidth proxies on the development machine; NFR-030 impact; AMR-benefit and subcycling trade-off; box-specific expectations; measurement plan. §7 evaluates **AMReX's distributed FFT Poisson solvers** (`FFT::Poisson`/`PoissonHybrid`/`PoissonOpenBC`) as a replacement for `pois.f90` on a uniform level 0. |
| `03-fr034-fr040r3-answers.md` | Answers to the Chief Architect's draft rulings (`adr/drafts/rulings-FR040R3-FR034.md`): FR-034 leak/HVAC owned-face mask, the composite PBAR profile, per-stage consistency with D-032 mean removal, one leak node across levels; FR-040 R3 E-2 β rule and MLMG C/F interpolation across thin walls. Verifies the draft's citations. |
| `04-masked-domain-pressure.md` | **DRAFT, incomplete (stopped by the scope change).** R-47 / FR-037 open point: pressure solve on a masked (non-box) level 0 (40 FR-006 cases). Recommends a single-level MLMG "gap mask"; FFT-based routes are rejected; acceptance criteria and the fit with D-021/D-032 (5)/FR-039. |

## Key recommendations (short form; details and citations in 01)
1. The FDS step is algebraically a MAC projection: `Ũ = uⁿ − δtF`, `∇·(δt∇H) = ∇·Ũ − D*` (corrector: β = δt/2, S = Dⁿ⁺¹). Use a
   cell-centred MAC projection with φ = H. No nodal projector.
2. Roadmap Phase 4: constant-coefficient operator with the lagged baroclinic term and the FDS baroclinic iteration (FR-036).
   Variable-β `∇·(δt/ρ̄ ∇p̃)` is an option (buildable from Fortran via `amrex_abeclaplacian`).
3. One global δt, no subcycling, composite MLMG at both solves (FR-030). Under this scheme no sync projection is required.
   R-05 is essentially a Phase 6 (subcycling) risk.
4. C/F: `average_down_faces` for face velocities (FR-032 exact), `average_down` for D/ρ/H, zone sums over uncovered cells
   (FR-034), per-zone compatibility, and the GLMAT mass-weighted gauge.
5. Solids: E-1 (IBM forcing, FDS default) first, E-2 (remove solids from the operator, UGLMAT-like) next, EB only with ADR-003.
6. Acceptance for H: ε_H = 1e-8 (mean-removed relative L2, single frozen solve, N ≤ 64) where the discretisations are identical.
   On C/F cases vs UGLMAT: FR-032 at machine zero, accuracy against the exact solution not worse than UGLMAT, and a 3e-2
   sanity bound (estimate).
7. MLMG uses `setMaxOrder(2)` (`AMReX_MLLinOp.H:310`; the default is 3 at `:886`), so its Dirichlet/OPEN ghost value is the same linear one
   that `FFT::Poisson` and Crayfishpak use. This is part of the FFT→MLMG acceptance check (FR-039), and the solver is chosen at every step from the hierarchy.
   C/F tangential interpolation stays at order 3 (`AMReX_MLCellLinOp.H:862`). Roadmap P2 tests orders 2 and 3 for observed order ≥ 1.8 at C/F faces;
   if order 2 misses, multi-level runs use order 3 (01 §C.1, REC-I3).

## Open questions / decisions needed
| # | Question | Options / my recommendation | Owner (suggested) |
|---|---|---|---|
| Q1 | FDS stretched meshes (`TRNX/Y/Z`) have no AMReX counterpart | drop in AMR mode (error at input), or replace by refinement. Rec: reject in AMR mode, keep in FFT uniform mode | Chief Architect |
| Q2 | Subcycling (Phase 6) at all? | Estimated gain ≤ 1.33× for N_c ≈ N_f (02 §4). Rec: decide after M4 measurements | Chief Architect + Pressure Lead |
| Q3 | Operator: constant β + baroclinic iteration (A2-a) vs variable β (A2-b) | Rec: A2-a for Phase 4. Revisit if the baroclinic iteration hits `MAX_PRESSURE_ITERATIONS` often | Pressure Lead |
| Q4 | Solids: IBM forcing (E-1) vs masked operator (E-2) vs EB (E-3) | Rec: E-1 then E-2. E-3 only with ADR-003 = EB, which forces C++. E-2 β rule settled (03 Q5): β = 0 on each level's own R1 wall faces; covered and C/F coarse faces always get MLMG's arithmetic fine average (`AMReX_MLABecLaplacian.H:497-510`, `799-815`). The gap-mask subset of E-2 is needed by Phase 2 for masked level 0 (Q19) | Chief Architect (ADR-003) |
| Q5 | Gauge in singular zones | Rec: the GLMAT mass-weighted p̃ gauge per zone (it affects F_B, so it is not cosmetic) | Pressure Lead |
| Q6 | HYPRE version for AMReX: local 3.0.0 (untested by AMReX CI) vs 3.1.0 (tested) | Rec: 3.1.0, or run an AMReX HYPRE smoke test on 3.0.0. Avoid `mac_proj.use_mlhypre` (SStruct; no `setLevelBC`) | Integration Lead |
| Q7 | Tunnel preconditioner | Obsolete under MLMG. Keep `tunnel_demo` as a performance/regression case. Use the hypre bottom solver if coarsening stalls | Pressure Lead |
| Q8 | Driver language (ADR-001) | Pressure is feasible in Fortran for A2-a/A2-b/E-1 (hand-written MAC glue). The overset mask, `MacProjector` and EB need C++ | Chief Architect (ADR-001) |
| Q9 | FFT retention in uniform mode (FR-037) | Estimated MLMG ≈ 6–18× an FFT solve, so NFR-030 is at risk unless f_pres is small. Rec: keep FFT for uniform single-box. Measure f_pres first. See Q16 for doing this with AMReX's FFT instead of `pois.f90` | Integration Lead |
| Q10 | Approve ε_H = 1e-8 and the C/F criteria for FR-016(c) | Needs the requirements owner to copy them into FR-016 (I did not edit requirements.md) | Spec & Program Lead |
| Q11 | Reference builds | References must come from the **FireX `36975d765f`** build: its AMG settings (PMIS/ℓ1-Jacobi) differ from master (HMIS/ℓ1-GS). An MPI runtime must be installed on the development machine to run any FDS reference | Integration / V&V Lead |
| Q12 | Refinement ratios ≠ 2/4, non-nested FDS meshes | **Closed.** `InterpBndryData` asserts ratio ≤ 4 whenever order-3 C/F interpolation is used (`AMReX_InterpBndryData.H:186-188`), and MLMG always uses order 3 there (`IBD_max_order_DEF`, `AMReX_MLCellLinOp.H:862`). So MLMG supports ratios 2 and 4 and aborts on 8. Ratio 3 passes this assert (3 ≤ 4), but its correctness elsewhere in MLMG and in the driver was not tested. Decision: support 2 and 4, reject everything else (3 included) in AMR mode | Inventory / Chief Architect |
| Q13 | Phase numbering: the dispatch said sync projections are "Phase 2"; `risks.md` R-05 and `roadmap.md` say Phase 6 | Rec: roadmap numbering (used in this folder) | Spec & Program Lead |
| Q14 | ADR-002 leaning text says the pressure iteration survives only at solids with masking | Correction: it survives at solids only with IBM forcing (E-1), not with masking (E-2), and always for the baroclinic term (FR-036) unless A2-b is chosen | Chief Architect |
| Q15 | FDS "patch-averaged" interface velocities (`MATCH_VELOCITY`, `velo.f90:2723-2733`) | Replaced by `average_down_faces` at C/F boundaries. Confirm that no other FDS routine depends on the averaged-then-overwritten values | Inventory (Legacy Mapper) |
| Q16 | Replace `pois.f90` by AMReX `FFT::Poisson` for uniform level 0 (02 §7) | **Supported for:** a uniform, unstretched level 0 with one BC type per domain face; BC values lifted into the RHS as Crayfishpak does. Removes the interface error and the interface iteration (equivalent to GLMAT); IBM, mixed-face and baroclinic iterations remain. **Gaps:** stretching (Hybrid: z only), mixed open/closed faces (iteration or a capacitance method), C++ only, multi-node all-to-all cost (estimated ≈ MLMG). Rec (REC-P7): adopt it as the single uniform-mode solver in the AMReX driver and do not port `pois.f90`. Decide after measuring k and f_pres | Chief Architect (with ADR-001) |
| Q17 | FR-034 leak area across levels (Architect draft Ruling 2) | Rec (03 Q1–Q4): one owned-face mask (valid, uncovered gas side) for `NODE_*`, the `MFT` application and `USUM`, with exact sums. FDS already masks `USUM` (`divg.f90:760`) but not `HVAC_BC_IN`. One composite `PBAR = P_0(z) + ΔP_zone(t)` sampled per level. Log the removed mean per zone (it hides face-set mismatches). Multi-zone leakage gates Phase 4, not Phase 6. Amendment applied in 01 §F (REC-F1) | Chief Architect; Spec Lead (FR-034 text) |
| Q18 | MLMG C/F tangential interpolation (order 3, hard-coded) crosses thin walls normal to the interface (`AMReX_InterpBndryData_3D_K.H:23-51`; the mask ignores walls; β does not enter) | Rec (03 Q6): R3-T plus the standard MLMG C/F; no masked interpolant now. Error O(H jump), which is dynamic-pressure scale because zone pressure lives in PBAR. Add a C/F-face flux check to the duct_flow_uglmat_refine acceptance. Fallback: AMReX patch for order-1 tangential interpolation on wall-adjacent faces | Chief Architect; Integration Lead (fallback patch) |
| Q19 | Masked (non-box) level 0, R-47 / FR-037 (40 FR-006 cases, 19 with several zones, OPEN vents on gap faces in the restart anchors and `hallways`) | **Draft only (04 incomplete):** single-level MLMG gap mask (overset 0 plus β = 0 on gap faces, β = 2 plus a known value behind gap-face OPEN vents, one pin per sealed component, D-032 mean removal and gauge). Amend D-032 (5)/FR-039: "FFT on a single uniform level with a box domain". Acceptance: FR-006 T2 vs the default FDS; eps_H only vs FDS UGLMAT on frozen input. Pull the gap mask into Phase 2. Open: `stairwell` cost (fill 0.091) | Chief Architect; Integration Lead; Spec Lead |

## Not verified (carried from 01/02)
* The exact sign convention of `MacProjector`'s velocity update. Native MLMG handling of an incompatible RHS in singular problems
  is now answered by P2 (D-032: reports convergence while the true residual stays at 2.6e-5·max|b|). Per-cell Robin coefficients via `setLevelBC`. Coarse-level β averaging is now verified by reading (03 Q5). Its convergence effect with disconnected components is not. MLMG on a non-covering level 0 with a Neumann C/F BC (04 (a1)). The Fortran E-2 emulation with
  an a-coefficient in solid cells. A divergence-preserving face interpolater in the target AMReX version. MLMG support for
  refinement ratio 3.
* AMReX with the local HYPRE 3.0.0.
* AMReX `FFT::Poisson`: read from the source only, never compiled or run. Transpose and all-to-all costs are estimates from assumed bandwidths (02 §7.0 F3). How `SAVE(IA)` is built in Crayfishpak was not traced.
* All FDS and AMReX timings. Only NumPy proxies were measured. Also unverified: the FDS memory per cell and f_pres.
* Literature listed on the AMReX-Hydro page but not checked: Bell, Colella & Glaz (1989); Almgren, Bell & Szymczak (1996);
  Bell & Marcus (1992).
