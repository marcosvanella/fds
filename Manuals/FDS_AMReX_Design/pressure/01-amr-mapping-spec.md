# 01 — Mapping the FDS pressure / divergence / velocity-correction cycle onto AMReX AMR

**Status:** Draft v1 (first spec). Author role: numerical analyst (pressure/projection).
**Depends on:** `00-fds-pressure-baseline.md` (cited below as "00 §n"). FDS citations use `file:line` in the
**FireX `AMReX` branch at commit `36975d765f`**, checked out read-only at this repository (`Source/*.f90`). They were re-based
from the earlier `ce1f659` draft by diff-hunk mapping with a text check (00 header, 00 §11). `TechGuide/…` is `Manuals/FDS_Technical_Reference_Guide/…`.
**External sources:** AMReX and AMReX-Hydro documentation and source as fetched on 2026-09-25 (URLs given inline). The
literature is listed in §R. Only references whose bibliographic data I checked are listed there.

Legend: **[REC]** = recommendation. **[ALT]** = alternative considered. **[OPEN]** = decision needed (collected in `README.md`).
**[VERIFY]** = a claim about AMReX/HYPRE behaviour that I have not checked in source or by a run, and that must be confirmed in the prototype.

**Programme cross-references (read, not edited; 2026-09-25):**
* `docs/requirements.md`: FR-016 (static two-level equivalence; UGLMAT-HYPRE reference; H tolerance "TBD(Pressure Lead)").
  This spec supplies that number in §I. Also FR-002, FR-023, FR-030…FR-037.
* `docs/risks.md` R-05 (deferring sync projections) and R-20 (zones).
* `docs/roadmap.md`: Phase 4 / M4 "Composite pressure" with the R-05 review gate; Phase 6 "Subcycling and synchronization
  projections". Note: the dispatch message calls the sync-projection phase "Phase 2", but `risks.md`/`roadmap.md` say
  Phase 6. **This spec follows the roadmap's numbering.**
* `docs/README.md`: ADR-002 (global dt, leaning) and ADR-003 (EB vs OBST masking).
* `docs/adr/ADR-001-driver-architecture.md`: it asks the Pressure Lead a question, answered in §G.1.
* `docs/inventory/README.md` and `pressure_fields_access.csv` (Legacy Mapper). Their CODE 5 / `COPY_H_OMESH_TO_MESH` / PRHS
  claims were re-verified and folded into 00 §8.

---

## 0. Assumptions (specific to this spec)

A1. **Scope.** The gas-phase, Cartesian, 3-D path of FDS, with the defaults `SOLVER='FFT'`, `BAROCLINIC=.TRUE.`
    (`cons.f90:229`) and pressure iterations on (`ITERATE_PRESSURE=.TRUE.` unless `VELOCITY_TOLERANCE>100`,
    `read.f90:10157-10163`). Cylindrical, level-set-only and cut-cell (`CC_IBM`) paths are out of scope (00 §0).
A2. **Grids.** Each AMR level is a **uniform** Cartesian grid, as AMReX requires. FDS stretched meshes (`TRNX/TRNY/TRNZ`,
    which Crayfishpak supports in up to two directions, 00 §4) have **no direct counterpart**. This spec assumes they are
    either dropped or replaced by refinement **[OPEN Q1]**.
A3. **Refinement ratio.** Refinement ratio 2 between levels (4 allowed). Arbitrary FDS mesh-to-mesh ratios (e.g. 3:1) and
    non-nested, non-aligned FDS meshes cannot be represented exactly. I have not checked whether AMReX MLMG supports ratio 3.
A4. **Refinement is static in the first multi-level phase (roadmap Phase 4 / milestone M4, "Composite pressure").** The user's `&MESH` blocks define fixed refined regions, like FDS embedded meshes
    (`Verification/Pressure_Solver/dancing_eddies_embed.fds`). Dynamic tagging comes later. The design must not preclude
    it (§B.4).
A5. **Target framework.** AMReX core (MLMG linear solvers) plus the AMReX-Hydro `MacProjector`
    (https://github.com/AMReX-Fluids/AMReX-Hydro/blob/development/Projections/hydro_MacProjector.H). I assume a recent
    AMReX (docs version string "26.10-dev" on https://amrex-codes.github.io/amrex/docs_html/LinearSolvers.html).
A6. **Parallelism.** MPI across boxes, with OpenMP/GPU optional. Nothing in this spec depends on GPU.
A7. **Other physics unchanged.** The divergence `D` (`DIVERGENCE_PART_1`, `divg.f90:22-785`), the flux term `F`
    (`VELOCITY_FLUX`, `velo.f90:563`) and the zone model are assumed to be ported by other teams. This spec fixes only
    **what the pressure step needs from them**: where D lives, how covered cells are treated, and which integrals are formed.
A8. **No results from running AMReX are claimed.** Statements marked [VERIFY] come from reading docs/source, not from execution.

---

## A. Projection formulation and discretisation

### A.1 Exact algebraic equivalence between the FDS step and a MAC projection

FDS predictor (00 §1.2, `velo.f90:1603-1630`; RHS at `pres.f90:250-261`; DDDT at `divg.f90:1610-1619`):

```
u* = uⁿ − δt (F + ∇H),      ∇²H = −∇·F − DDDT,      DDDT = (D* − ∇·uⁿ)/δt
```

Define the face-centred trial velocity `Ũ = uⁿ − δt F`. Then `δt ∇²H = ∇·Ũ − D*`, and `u* = Ũ − δt ∇H` gives
`∇·u* = D*` exactly, up to solver tolerance. This is precisely a **MAC projection** with constant coefficient β = δt,
source `S = D*` and potential `φ = H`:

```
∇·(β ∇φ) = ∇·Ũ − S ,   U = Ũ − β ∇φ .
```

FDS corrector (`velo.f90:1723-1751`, DDDT at `divg.f90:1620-1631`):
`uⁿ⁺¹ = ½(uⁿ + u* − δt(F* + ∇H*))`, with `DDDT = (2Dⁿ⁺¹ − ∇·uⁿ − ∇·u*)/δt`.
Set `Ũ = ½(uⁿ + u*) − (δt/2) F*`, β = δt/2, `S = Dⁿ⁺¹`, `φ = H*` (FDS array `HS`). The same algebra gives
`∇·uⁿ⁺¹ = Dⁿ⁺¹`.

So the FDS step maps onto AMReX-Hydro's `MacProjector` exactly. The API as fetched has `setUMAC`, `setDivU`
(the cell-centred source S), `updateBeta(Real const_beta)` (constant β → `MLPoisson` internally),
`project(phi_in, reltol, atol)` with an initial guess, `setDomainBC`, `setLevelBC`, `setCoarseFineBC` and `getFluxes`
(`hydro_MacProjector.H:26-160`). The RHS it forms is `(∇·Ũ − S)` scaled by 1/β (MLPoisson) or −1 (MLABecLaplacian)
(`hydro_MacProjector.cpp:381-394`). **[VERIFY]** I did not trace the final sign convention of the velocity update in the
`.cpp`. A unit test must confirm `∇·U = S` after `project()`.

Properties this preserves:
* **FDS's self-correcting source.** FDS uses the *actual* discrete `∇·uⁿ` (00 §2.3). The MAC form uses `∇·Ũ` computed from
  the actual `uⁿ`, so any divergence error left in `uⁿ` (from regridding, tolerance or interpolation) is removed at the next
  projection. This matters for AMR (§B.4).
* **Potential equals H, not an increment.** `Ũ` contains no `∇H`, so φ is H itself. The previous H (or HS) is an
  excellent initial guess, and MLMG usually needs few V-cycles from it (estimate in `02` §3).
* **KRES.** H still means `KRES + p̃/ρ` (`velo.f90:283-292`, `velo.f90:3257`). KRES uses face-averaged velocities, so it
  needs one ghost layer of face velocities at C/F boundaries (filled by FillPatch).

**[REC-A1]** Use `Hydro::MacProjector` with face-centred `Ũ`, cell-centred φ = H, constant β (δt or δt/2), `setDivU(D)`
on every level, and a composite multi-level solve (§C). Keep the FDS names H/HS for φ in predictor/corrector.

**What the operator needs from AMReX, and C++ vs Fortran availability.** The projection needs four things: a cell-centred,
multi-level operator (constant coefficient for A2-a; face-centred variable β for A2-b and for β=0 solid faces in E-2); the
face gradient/flux of the solution (velocity update); domain, level and C/F BCs; and optionally an overset mask.
* **AMReX Fortran interface** (checked in `(local AMReX checkout)` @ `99ddfda`). It exposes only `amrex_poisson` and
  `amrex_abeclaplacian`, both cell-centred and multi-level. Bottom solvers are smoother, bicgstab, cg, **hypre** and **petsc**
  (`Src/F_Interfaces/LinearSolvers/AMReX_multigrid_mod.F90:9-14`).
* `amrex_abeclaplacian` has `set_scalars`, `set_acoeffs` and `set_bcoeffs` (`AMReX_abeclaplacian_mod.F90:14-16`). So a variable-
  coefficient `∇·(β∇p)` with face-centred β (e.g. δt/ρ̄_f for A2-b, or β = 0 on solid faces for E-2) **can be built from Fortran**.
* The Fortran linop also has `set_maxorder`, `set_domain_bc`, `set_coarse_fine_bc` and `set_level_bc`
  (`AMReX_linop_mod.F90:12-15`). The multigrid object has `solve`, `get_grad_solution`, `get_fluxes` and `comp_residual`
  (`AMReX_multigrid_mod.F90:24-27`), which is enough to hand-assemble the MAC projection: compute `∇·Ũ − S`, solve, then
  `U = Ũ − β∇φ` from `get_fluxes`.
* **Not available from Fortran:** the `MacProjector` convenience class (AMReX-Hydro, C++ only), the **overset mask** (grep of
  `F_Interfaces/LinearSolvers` finds no mask argument), nodal operators (`MLNodeLaplacian`), EB operators (`MLEBABecLap`) and
  tensor operators. All of these need C++.

**[REC-A1b]** The operator choice does not force the driver language for the roadmap Phase 4 scope:
* A2-a/E-1 needs only `amrex_poisson`. A2-b needs only `amrex_abeclaplacian`.
* E-2 without an overset mask can be emulated in Fortran with `amrex_abeclaplacian`: set β = 0 on faces touching solids, and
  give solid cells a unit `a`-coefficient with zero RHS (α·a·φ decouples them). This emulation is **[VERIFY]**: with a ≠ 0
  somewhere, MLMG will not treat the problem as singular, so per-zone compatibility (§C.4) must be enforced by the application.
* E-3 (EB) and `MacProjector` reuse require C++.

See §G.1 for the ADR-001 implications.

### A.2 Constant coefficient with lagged baroclinic term (FDS form) vs. variable-coefficient projection

FDS keeps the operator constant-coefficient and moves the variable-density effect into the lagged baroclinic term
`F_B = −p̃ ∇(1/ρ)` (`BAROCLINIC_CORRECTION`, `velo.f90:3216-3315`; rationale in `TechGuide/Momentum_Chapter.tex:354`).
It then iterates until the variable-coefficient residual `PRESSURE_ERROR_MAX` (`pres.f90:747-795`) is below
`PRESSURE_TOLERANCE` (exit test `main.f90:1724`; default `20/min(1,Δ)²`, `read.f90:10162`).

| | **A2-a: FDS form (constant β, lagged F_B, outer iteration)** | **A2-b: variable-coefficient `∇·((δt/ρ̄_f)∇p̃)`** |
|---|---|---|
| AMReX operator | `MLPoisson` via `MacProjector::updateBeta(Real)` | `MLABecLaplacian` with face β = δt·2/(ρᵢ+ρᵢ₊₁) (the face density that FDS already uses in the residual check, `pres.f90:747-795`) |
| Fidelity to FDS | identical discrete equations. Results comparable with FDS to round-off/solver tolerance on single-level cases | different momentum split: the solution variable is p̃, and F must exclude F_B and ∇K. Results not directly comparable with FDS |
| Iterations | 1 solve if `PRESSURE_ERROR_MAX < tol`, otherwise repeated (FDS allows up to `MAX_PRESSURE_ITERATIONS=10`, `cons.f90:559`) | 1 solve per half-step (no baroclinic iteration) |
| MLMG convergence | best case (constant coefficient) | degrades with density contrast. Fire contrast ρ_amb/ρ_flame ≈ 1.2/0.2 ≈ 6 is moderate. The literature method accommodates large contrasts (Day & Bell 2000; Almgren et al. 1998) |
| Gauge sensitivity | **F_B depends on the absolute level of p̃**, so the gauge matters (§C.5) | gauge only affects the p̃ output |

**[REC-A2]** Roadmap Phase 4: A2-a, which gives bit-for-bit comparable equations and the easiest V&V against FDS. Implement the
`PRESSURE_ERROR_MAX` diagnostic over the composite grid, since it is needed anyway. Keep A2-b as a flagged experiment
for cases where the baroclinic iteration hits its iteration cap **[OPEN Q3]**.

### A.3 Why not a nodal (approximate) projection or a cell-centred approximate projection

* FDS velocities are MAC-staggered (00 §1.2) and D is imposed exactly on the MAC divergence. AMReX-Hydro's
  `NodalProjector` is an approximate projection for **cell-centred** velocities with nodal φ
  (https://amrex-fluids.github.io/amrex-hydro/docs_html/Projections.html). Using it would force a change of the momentum
  discretisation.
* The Martin & Colella (2000) cell-centred approximate projection has the same objection.
* **[REC-A3]** Use MAC projection only. Nodal or approximate projections are rejected for this port.

### A.4 Single-level FFT fast path (optional)

AMReX-Hydro also contains `Hydro::FFTMacProjector` (`hydro_FFTMacProjector.H`, built on `amrex::FFT::Poisson`). As fetched,
its API takes `Geometry`, low/high BC types and `setUMAC`, with **no `setDivU`**, and it is single-level. It is the
closest analogue of the FDS Crayfishpak solve (`pois.f90:187`), but it would need a source term to be usable here.
**[ALT]** Keep it as a possible uniform-level-0 accelerator or reference. Not needed for M4. Relevant to FR-037 (FFT retention in uniform mode), which can keep the existing Crayfishpak path per box.

**Update (distributed FFT evaluation, 02 §7).** `amrex::FFT::Poisson` (`AMReX_FFT_Poisson.H`, AMReX `99ddfda`) can replace `pois.f90` for uniform level-0 runs inside the AMReX driver. It covers periodic, Neumann (`even`) and Dirichlet (`odd`) BCs, one type per domain face (`AMReX_FFT_Helper.H:59`; `AMReX_FFT_Poisson.H:41-42`, `313-342`). BC values must be lifted into the RHS the same way Crayfishpak already does (`pois.f90:240-270`). One global FFT gives GLMAT's solution (no interface error, no interface iteration). IBM, mixed-face and baroclinic iterations remain. `PoissonHybrid` allows stretching in z only; `PoissonOpenBC` does not apply. **[REC-P7]** One projection path, with the linear solver chosen by configuration (FFT when `finest_level = 0`, the grid is unstretched and every domain face has a single BC type; MLMG otherwise). See 02 §7 and README Q16.

---

## B. Time stepping and subcycling

### B.1 FDS constraints that bind the choice
* FDS uses a **global** time step `DT = MINVAL(DT_NEW)` (`main.f90:715`). Stability is checked per mesh in
  `CHECK_STABILITY` (`velo.f90:3028-3209`).
* FDS uses a two-stage predictor–corrector with **two pressure solves per step** plus iterations (`main.f90:855`, `1090`).
* If CFL/VN is violated after the predictor, the predictor is **re-run** with a reduced DT
  (`CHANGE_TIME_STEP_LOOP`, `main.f90:774-897`; redo at `890-895`).

### B.2 Options

| | **B-1: no subcycling (all levels share δt)** | **B-2: subcycling in time (Berger–Oliger style)** |
|---|---|---|
| Precedent | PeleLMeX (Esclapez et al. 2023: "without level sub-cycling") | IAMR / Almgren et al. 1998. Low-Mach combustion: Pember et al. 1998; Day & Bell 2000 |
| Pressure step | one **composite** MAC projection over all levels per half-step | level projections with coarse-supplied C/F Dirichlet data, then a MAC sync / **sync projection** on the composite grid after fine levels catch up (Almgren et al. 1998) |
| FDS predictor–corrector and DT redo | unchanged: redo the whole hierarchy | must be redesigned per level. A DT reduction on a fine level interacts with the coarse step |
| Zone dP̄/dt (`divg.f90:1517-1521`) | one global integral per step, as in FDS | integrals at different times on different levels need time interpolation and a sync correction for P̄. No precedent in FDS |
| DDDT self-correction | exact per step | covered-coarse and fine divergences are reconciled only at sync |
| Cost | coarse levels take the fine δt. Waste ∝ fraction of cells on coarse levels (`02` §4) | optimal work per level |
| Implementation risk | low | high |

**[REC-B1]** Roadmap Phase 4: **B-1 (no subcycling)**. This matches the recorded ADR-002 leaning and FR-030. Keep the FDS predictor–corrector, the global δt and
`CHANGE_TIME_STEP_LOOP` semantics unchanged, with two composite projections (plus baroclinic iterations) per step.
Revisit B-2 only if profiling shows the coarse-level waste dominates **[OPEN Q2]**.

### B.2a What "deferring sync projections" (R-05) actually defers
In IAMR/Almgren et al. (1998), the sync projection corrects the composite-divergence mismatch that arises **because
levels are advanced with different time steps**. The coarse face velocities used during the coarse step differ from the
time-averaged fine face velocities. With **B-1 there is no such mismatch**:
* one composite solve per half-step enforces the divergence constraint on the composite grid;
* `average_down_faces` (§C.2) makes the coarse C/F face velocities equal to the fine averages at the same time level;
* same-time flux replacement (FR-024; a "reflux" with no time integration) keeps the scalars conservative;
* the zone integrals are formed once per step over uncovered cells (§C.3).

So under B-1, no sync projection is needed for correctness. PeleLMeX (no subcycling) has none. **R-05 is therefore
mainly a subcycling risk** (Phase 6). Under B-1 the remaining C/F risk is the ordinary C/F discretisation error of the
composite operator and of the scalar fluxes. It is measured by the R-05 trigger ("interface-adjacent divergence
residual > 10× interior"), which should read ≈ solver tolerance everywhere under B-1. **[REC-B3]** Keep the R-05 gate
at M4 as a check on this argument, using diagnostics (i) and (iii) of §H.

### B.3 Placement in the step (B-1)
Mirrors `main.f90` (00 §3, §9). Each "for all levels" is a single MultiFab-vector operation:
1. `DENSITY` on all levels, then `average_down` of ρ, Z to covered coarse cells, then FillPatch ghosts.
2. `VELOCITY_FLUX` → F on all faces of all levels.
3. `DIVERGENCE_PART_1` → D on valid cells. Zone sums DSUM/PSUM/USUM over **uncovered** cells only (§C.3), then
   `ParallelAllReduce`, then `DIVERGENCE_PART_2` → dP̄/dt, D*, then `average_down(D)` (§C.2).
4. Pressure iteration loop (§D): `BAROCLINIC_CORRECTION` → (solids forcing, §E) → form Ũ → `MacProjector::project`
   → residual diagnostics → exit tests copied from `main.f90:1724-1741`.
5. Velocity update `U = Ũ − β∇H` (done by the projector). Then `average_down_faces` (§C.2), then `CHECK_STABILITY` on
   each level with a global MIN, then possibly a DT redo.

### B.4 Regridding (needed later; must not be precluded now)
* On regrid, H (or HS) is interpolated as an initial guess only. Face velocities on new fine patches must be interpolated
  from coarse faces. **[VERIFY]** AMReX offers face interpolaters; I did not verify whether a divergence-preserving one is
  available in the target version.
* Thanks to the DDDT self-correction (§A.1), a small divergence error introduced by the interpolation is removed at the
  next predictor projection. It pollutes only the one scalar-advection sub-step in between. **[REC-B2]** Accept this in
  the first regridding phase. Optionally run an extra projection with `S = D` right after regrid for strict V&V runs.

---

## C. Coarse/fine (C/F) handling

### C.1 Composite solve
With B-1, one MLMG composite solve over levels 0…L (`MLMG::solve` on vectors of MultiFabs) gives H on the valid region of
every level. The C/F flux matching (the coarse face flux equals the sum of the fine face fluxes) is built into the composite
residual. This **replaces** both FDS mechanisms:
* the Dirichlet interface iteration of FFT/ULMAT (`pres.f90:115-143`; `TechGuide/Momentum_Chapter.tex:477-505`);
* the GLMAT/UGLMAT two-point C/F flux `A=min(A_int,A_ext)`, `1/(δx_int+δx_ext)` (`pres.f90:5223`, `5229`).

The FDS two-point flux ignores the tangential offset of cell centres (00 §6.2). AMReX interpolates the coarse solution
into fine ghost cells with a polynomial whose order is set by `setMaxOrder` (LinearSolvers docs). **Expect differences
from FDS GLMAT at C/F boundaries.** `Verification/Pressure_Effects/obst_coarse_fine_interface` (with `_exact.csv`) is the
discriminating test (§H).

**Decision: `setMaxOrder(2)` on the MLMG path** (AMReX `99ddfda`; roadmap P2 / FR-039).
* **Setting.** `setMaxOrder` is at `AMReX_MLLinOp.H:310`; the default is `maxorder = 3` (`AMReX_MLLinOp.H:886`).
* **Effect at Dirichlet and OPEN faces.** The Dirichlet ghost value uses `NX = min(blen+1, maxorder)` points (`AMReX_MLLinOp_K.H:124-138`). With order 2 this is the linear
  ghost value `φ_g = 2g − φ₁`, which is the value `FFT::Poisson` (02 §7.2) and Crayfishpak (`pres.f90:454`) use. This makes FFT and MLMG
  agree within ε_H on identical single-level input, and that agreement is part of the acceptance check for switching from FFT to MLMG.
  The driver chooses the solver at every step from the hierarchy (02 §7.9).
* **Effect at C/F faces.** `maxorder` also sets the extrapolation normal to C/F faces, because `applyBC` passes it to the
  boundary kernels (`AMReX_MLCellLinOp.H:999-1092`). The tangential interpolation of coarse values onto the face stays at
  `IBD_max_order_DEF = 3` (`AMReX_MLCellLinOp.H:862`; `AMReX_InterpBndryData.H:126`).
* **P2 test.** P2 must run the refinement study at maxorder 2 and 3 and require an observed order ≥ 1.8 at the C/F interface (§I.3).
  If order 2 misses, keep order 3 for multi-level hierarchies (the FR-039 fallback) and keep order 2 for single-level runs.

### C.2 Face velocities and covered cells
* **Face velocities.** After the projection, coarse faces **on** the C/F boundary and coarse faces **covered** by fine
  cells are overwritten by the area average of the fine faces (`amrex::average_down_faces`). This is the AMR equivalent
  of the FDS averaging of interface normal velocities (`MATCH_VELOCITY`, `velo.f90:2723-2733`; flux averaging in
  `MATCH_VELOCITY_FLUX`, `velo.f90:2951`). The composite divergence is then exactly `S` on uncovered coarse cells
  adjacent to C/F boundaries.
* **Cell data.** D, ρ, H: `average_down` (volume-weighted) into covered coarse cells after each update. D on covered cells
  must equal the average of the fine D. Otherwise the coarse-level smoother sees an inconsistent RHS, and the
  convergence of the composite solve degrades (the residual on covered cells is not used, but the coarse-level
  correction is).
* **Reflux of scalar fluxes.** Needed for conservation of ρ, ρY and energy across C/F boundaries (`FluxRegister`). This is
  outside the pressure step, but D and the zone integrals assume conservative transport.

### C.3 Zone integrals over uncovered cells
FDS forms DSUM and PSUM only over cells **not covered by a lower-numbered (finer, embedded) mesh**, using the
`INTERPOLATED_MESH` mask (`divg.f90:727-753`, test at `divg.f90:734`, `760`; set at `init.f90:1569-1584`). USUM is
summed over zone boundary faces (`divg.f90:757-778`) and reduced globally in `EXCHANGE_DIVERGENCE_INFO`
(`main.f90:2028-2054`).
**[REC-C1]** Compute DSUM/PSUM as Σ over levels of Σ over **valid, uncovered** cells (fine mask from `makeFineMask`,
weighted by cell volume). Compute USUM over boundary faces of the finest level that covers each portion of the boundary.
One `ParallelAllReduce` per zone. This is the exact AMR analogue of the FDS mask.

### C.4 Compatibility (solvability) of singular problems
For a zone with no Dirichlet face, solvability requires `∫_zone (∇·Ũ − S) dV = 0`, i.e.
`∮ Ũ·n = ∫ D`. FDS makes this hold by construction: `dP̄/dt = (DSUM−USUM)/PSUM` and
`D ← D − (R_PBAR − RTRM)·dP̄/dt` (`divg.f90:1517-1521`, `1540`). It then still subtracts the RHS mean before the solve, for
round-off and discretisation mismatch (ULMAT `pres.f90:1683-1737`; GLMAT `pres.f90:3369-3406`).
**[REC-C2]** Keep both. (i) Form S exactly as FDS does, from the composite sums of §C.3. (ii) Before `project()`, subtract
from the RHS the composite, volume-weighted mean **per connected zone**. Do this in application code rather than relying
on solver internals. **[VERIFY]** The AMReX docs document `hypre.adjust_singular_matrix` for the HYPRE path. I did not
verify how native MLMG handles a non-compatible RHS in the singular case.

### C.5 Gauge
F_B uses the absolute level of p̃ (`velo.f90:3254-3272`), so the additive constant in H changes the solution when the
problem is singular. FDS GLMAT/UGLMAT shift H so that the **mass-weighted mean of p̃ is zero**:
`SHIFT_H = Σ Vρ(KRES+H)/Σ Vρ` (`pres.f90:3495-3548`), with a zero-mean fallback (`pres.f90:3549-3557`). The FFT path
instead uses a singular perturbation `POIS_PTB` (`pres.f90:346`).
**[REC-C3]** After every singular solve, apply the GLMAT gauge per zone over uncovered cells. This keeps AMR results
comparable to FDS GLMAT/UGLMAT runs **[OPEN Q5]**.

### C.6 Proper nesting and buffers
AMReX requires proper nesting (with `blocking_factor` and `n_error_buf`). FDS has no such requirement. User meshes that
violate nesting must be adjusted when a `&MESH` layout is converted into an AMR hierarchy (inventory/ADR teams). The
pressure solver itself needs nothing beyond proper nesting.

---

## D. Fate of the FDS pressure-iteration loop

FDS iterates for three distinct reasons (00 §5.2, `main.f90:1601-1745`):

| Reason | FDS mechanism | AMR (REC) |
|---|---|---|
| (1) Mesh-interface mismatch | Dirichlet interface + `VELOCITY_ERROR` at `INTERPOLATED` faces (`pres.f90:802-1076`) | **Eliminated** by the composite solve (§C.1): same-level box boundaries are internal to MLMG, and C/F boundaries are handled by the composite operator |
| (2) Solid no-flux | IBM forcing in `NO_FLUX` (`velo.f90:1348-1563`) with the `WALL_WORK1` correction (`pres.f90:828-832`, `1056`) | Kept only under option E-1. Eliminated under E-2 (§E) |
| (3) Baroclinic lag | `ITERATE_BAROCLINIC_TERM` until `PRESSURE_ERROR_MAX < PRESSURE_TOLERANCE` (`main.f90:1609-1613`, `1724`) | **Kept** under A2-a. Eliminated under A2-b |

**Correction to the ADR-002 leaning text** (`docs/README.md`), which says the pressure iteration "survives only at solid
boundaries, and only if obstruction masking is kept". Two corrections:
* Under E-1, reason (2) keeps the loop alive. Under E-2 (masking), reason (2) is **removed**, because no-flux is exact.
* Independently of solids, reason (3), the baroclinic iteration, survives under A2-a. FR-036 requires its semantics to
  be kept.

So the loop survives for solids only with IBM forcing (not with masking), and always for the baroclinic term unless A2-b
is chosen.

**[REC-D1]** Keep the `PRESSURE_ITERATION_SCHEME` skeleton, its inputs and its exit logic (`MAX_PRESSURE_ITERATIONS`,
`MAX_PREDICTOR_PRESSURE_ITERATIONS` (`read.f90:10199`), `SUSPEND_PRESSURE_ITERATIONS`, `main.f90:1724-1741`) so
that the input files keep their meaning. `VELOCITY_ERROR_MAX` then measures only the solid normal-velocity error (E-1)
and is identically ≈0 under E-2.

**[REC-D2] Inner solver tolerance.** FDS HYPRE uses PCG relative 2-norm tolerance 1e-12 (`imkl.f90:236-263`). MLMG uses
a max-norm residual (LinearSolvers docs). A residual `r` in `∇·U − S` produces a velocity error of order `r·h`
(one cell). Estimated requirement: `r_max·h_finest ≤ 0.01·VELOCITY_TOLERANCE`, so that the solver error is two orders
below the FDS tolerance. With the default `VELOCITY_TOLERANCE = 0.5·Δ` (`read.f90:10161`), this is
`r_max ≤ 0.005` s⁻¹ (numerically independent of Δ). Use `tol_abs` from that bound together with `tol_rel = 1e-10`,
whichever is reached first. *Estimate, to be calibrated on `pressure_iteration3d_default`.*

**[REC-D3] Tunnel preconditioner.** `TUNNEL_POISSON_SOLVER` (`pres.f90:505-697`) is a one-directional global coarse
correction (00 §7). Multigrid provides global coupling on its coarse levels, so the preconditioner has no role in the AMR
code. **Risk:** very long thin domains coarsen poorly in MLMG, because the coarsening stops when a box dimension reaches
its minimum, leaving a large bottom problem. Mitigations: LPInfo agglomeration/consolidation, or the `hypre` bottom solver
(§G). Keep `tunnel_demo` as a regression test **[OPEN Q7]**.

---

## E. Solid obstructions

| | **E-1: IBM forcing on all-cell operator (FDS FFT/GLMAT default)** | **E-2: remove solids from the operator (FDS ULMAT/UGLMAT analogue)** | **E-3: embedded boundaries (MLEBABecLap)** |
|---|---|---|---|
| FDS reference | `NO_FLUX` forcing + iteration (`velo.f90:1348-1563`) | UGLMAT gas-only numbering (`pres.f90:5750-5759`), exact Neumann | cut-cell `CC_IBM` (out of scope) |
| AMReX realisation | `MLPoisson` (constant β) over all cells. F on solid faces forced as in `NO_FLUX` | `MLABecLaplacian`, β_face = 0 on every face touching a solid cell or thin obstruction. Solid cells marked "known" with the **overset mask** (1 = unknown, 0 = known; supported for cell-centred MLABecLaplacian/MLPoisson, LinearSolvers docs; `MacProjector` accepts `a_overset_mask`, `hydro_MacProjector.H:53`) | `MLEBABecLap`, homogeneous Neumann on EB by default (`setEBDirichlet` available) |
| No-flux accuracy | approximate, ≤ `VELOCITY_TOLERANCE` after iterations | exact (to solver tolerance) in one solve | exact on the EB |
| Zones separated by solids | the operator stays connected through solid cells, so there is one null space | **one null space per disconnected zone**: needs per-zone compatibility (§C.4) and a per-zone pin (one mask = 0 cell per zone with a known value) or the per-zone gauge | as E-2 |
| MG robustness | best | variable β with zeros; thin walls can "leak" on coarse MG levels, which slows convergence but does not affect the fine-level answer. **Verified (AMReX 99ddfda):** `prepareForSolve` always calls `averageDownCoeffs` (`AMReX_MLABecLaplacian.H:497-510`). Face β is coarsened by the arithmetic mean of the r² fine subfaces, both between MG levels of one AMR level (`:711-737`) and onto the coarse AMR level under the fine level (`:799-815`; kernel `AMReX_MultiFabUtil_3D_C.H:95-120`). So covered and C/F coarse faces always get the fine average, whatever the driver sets. Under FR-040 R3-T (thin-wall band rule) that average is exactly 0 or 1 in the interface band. Fractional β appears only on M1–M4 covered faces and on coarse MG levels, where walls interior to a coarsened cell vanish. Both affect convergence only. **Rule (FR-040 R3, E-2):** β = 0 on the wall faces of each level's own R1 mask; covered/C/F coarse faces take MLMG's fine average (`03-fr034-fr040r3-answers.md` Q5). **[VERIFY]** convergence when thin walls split a level into disconnected components (per-component pin) | known to work in AMReX codes, higher implementation cost |
| Obstruction create/remove (`GLOBAL_MATRIX_REASSIGN`, `main.f90:1806-1835`) | nothing to rebuild | `updateBeta` + mask update + MLMG re-setup (O(N)) | EB regeneration (expensive) |

**[REC-E1]** Roadmap Phase 4: implement **E-1** first, for fidelity and direct comparability with the FDS default. Implement
**E-2** immediately after as the "tight" option, mirroring how FDS positions UGLMAT. This makes the pressure iteration
unnecessary for solids. **E-3** only together with a future GEOM/cut-cell port **[OPEN Q4]**.
Note for E-2: FDS sets `D = 0` in solid cells (`divg.f90:1547-1555`), which becomes irrelevant for masked cells but must
remain consistent in `average_down` (average over gas cells only, or mask-aware).

---

## F. Boundary conditions

Each FDS pressure BC maps as follows. The key point: in MAC form, a prescribed normal velocity on a boundary face (wall,
supply vent, mirror) is put into `Ũ` directly, and the Poisson BC is then **homogeneous Neumann**. FDS reaches the same
thing through `BXS…` values from `F` (`pres.f90:16-315`).

| FDS boundary | FDS pressure BC | AMReX `LinOpBCType` | Notes |
|---|---|---|---|
| Solid wall, supply/exhaust vent with prescribed velocity, `MIRROR` | Neumann (`init.f90:2478-2495`) | `Neumann` (homogeneous, with the boundary face velocity already in Ũ) | the projector does not modify boundary-face velocity |
| `OPEN` vent / open boundary | Dirichlet on H (`init.f90:2497`) | `Dirichlet`, with ghost values set by `setLevelBC` | H value as in `PRESSURE_SOLVER_COMPUTE_RHS` (`pres.f90:16-315`), incl. `DYNAMIC_PRESSURE` (`pressure_boundary` case) |
| Domain face that is **partly** OPEN and partly wall | the FDS FFT uses one type per face, and an OPEN vent makes the face Dirichlet (`init.f90:2497`). ULMAT/UGLMAT handle it per cell | `Robin` with per-ghost-cell coefficients (a=0 → Neumann, b=0 → Dirichlet). **[VERIFY]** per-cell Robin coefficient arrays via `setLevelBC` | improves on FFT fidelity. Matches UGLMAT |
| `PERIODIC` | periodic (FFT) | `Periodic` (Geometry periodicity) | |
| `INTERPOLATED` (mesh–mesh) | Dirichlet + iteration (`init.f90:2528`) | **none**: internal box boundary or C/F boundary | §C.1 |

**HVAC and leakage** (amended 2026-09-25 for FR-034; the former text said they were "unchanged"). Leak/HVAC faces are prescribed-normal-velocity faces: `CALC_HVAC_BC` sets `U_NORMAL(_S)` at every stage (`wall.f90:1745-1796`, called from `WALL_BC`, `wall.f90:178`). So in MAC form they enter `Ũ` as boundary-face velocities (homogeneous Neumann, first row of the table). They also enter the zone sums (`USUM`, `divg.f90:757-767`, giving dP̄/dt, `divg.f90:1519`) and D through `Q_LEAK` (`divg.f90:543-544`). The zone-pair leak areas (`LEAK_AREA`, `hvac.f90:3374-3375`; localized `AREA`, `hvac.f90:684-685`; pressure scaling `hvac.f90:3585-3609`) are hierarchy-independent scalars and are never summed from faces.

**[REC-F1] Owned-face mask.** Every face-distributed leak/HVAC quantity uses **owned faces only**: faces whose gas-side cell is valid, uncovered (`makeFineMask`) and not solid. This covers `NODE_AREA`, the area-weighted node state, the per-face `MFT` application and the leak/HVAC part of `USUM`. All use `B1%AREA` of the owning level and the FR-005 (ii) exact sum. FDS already masks `USUM` (`divg.f90:760`), but `HVAC_BC_IN` has no covered-cell test (`hvac.f90:2299-2539`), so the mask is new for the node sums and the flux application. Covered coarse leak/HVAC faces get no `CALC_HVAC_BC` flux; their velocity is `average_down_faces` of the fine faces (D-032 (3), §C.2). `PBAR` for `NODE_P`/ρ_F is the owning level's sample of the single composite profile `P_0(z) + ΔP_zone(t)` (`init.f90:452-484`, `mass.f90:548`, `730`). Check: Σ_owned MFT·A = DUCT_MF per step to round-off (FR-020). Details and the D-032 mean-removal diagnostic: `03-fr034-fr040r3-answers.md` Q1–Q4.

---

## G. Linear solver strategy and the role of HYPRE

### G.1 Answer to ADR-001's question to the Pressure Lead
ADR-001 asks whether the pressure design assumes a C++ projection (`MLPoisson`/`MLNodeLaplacian`). **Answer:** the design needs
a **cell-centred** MAC projection. `MLNodeLaplacian` is **not** used (§A.3).

| Pressure capability | C++ (option A) | Fortran interface (option B) |
|---|---|---|
| A2-a constant-coefficient composite solve | `MacProjector` → `MLPoisson` | `amrex_poisson` + hand-coded div/flux update (`AMReX_multigrid_mod.F90:24-27`) |
| A2-b variable β = δt/ρ̄_f | `MacProjector` + `MLABecLaplacian` | `amrex_abeclaplacian` `set_bcoeffs` (`AMReX_abeclaplacian_mod.F90:14-16`) |
| E-2 solids removed | overset mask (`hydro_MacProjector.H:53`) or β = 0 | β = 0 + a-coefficient emulation only (no mask), [VERIFY] |
| E-3 EB (ADR-003) | `MLEBABecLap` | **not available** |
| HYPRE / PETSc bottom | yes | yes (`AMReX_multigrid_mod.F90:9-14`) |
| Whole-hierarchy HYPRE (`mac_proj.use_mlhypre`) | yes (not recommended, §G.2) | no |

**Implication for ADR-001.**
* For the Phase 4 scope (A2-a, E-1, then E-2), **option B is sufficient** for pressure, at the cost of ~100–200 lines of
  hand-written projection glue (estimate). This is the same algebra as `MacProjector::project`.
* **Choosing EB (ADR-003) for solids makes C++ mandatory** for the pressure solve.
* Keeping FFT per box for uniform mode (FR-037) in the legacy Fortran code does not affect ADR-001. **Using AMReX's `FFT::Poisson` as the uniform level-0 solver instead (REC-P7, 02 §7) does:** it has no Fortran interface (`Src/F_Interfaces` has no FFT module) and its classes are C++20 templates (`AMReX_FFT_Poisson.H:43`). Option B would need a ~100–200-line `extern "C"` shim (estimate). That is one more argument for Option A, but it does not block Option B. It also has a GPU path (cuFFT/rocFFT/oneMKL, `AMReX_FFT_Helper.H:23-42`).

### G.2 Solver choices
**[REC-G1] Default: native MLMG** (geometric multigrid, O(N) per V-cycle), with the default bottom solver (bicgstab) or
`cg`, and an initial guess from the previous H (§A.1).

**[REC-G2] HYPRE uses (all optional, via AMReX's interface):**
1. **Bottom solver** (`bottom_solver = hypre`: BoomerAMG/PCG/GMRES; singular fix `hypre.adjust_singular_matrix`): for
   E-2 masking, elongated (tunnel) domains, and later EB. This is the lowest-risk use.
2. **HYPRE as the full solver** (`setMaxCoarseningLevel(0)` + hypre bottom): a *reference* configuration that approximates
   the FDS UGLMAT setup (PCG + BoomerAMG with PMIS coarsening, ℓ1-Jacobi relaxation (type 18), 1 sweep, tol 1e-12;
   `imkl.f90:236-263`; set at `pres.f90:4828-4838`). Use it
   for A/B comparisons with FDS, not for production.
3. **`mac_proj.use_mlhypre`** (`hydro_MacProjector.cpp:80-81`): solves the whole AMR hierarchy with
   `HypreMLABecLap`, which uses the HYPRE **SStruct** interface (`AMReX_HypreMLABecLap.H:175-182`). It asserts no EB
   (`hydro_MacProjector.cpp:92`) and **no `setLevelBC`** (`hydro_MacProjector.cpp:294`), so inhomogeneous Dirichlet
   (OPEN) values cannot be passed. **Not recommended.**

**HYPRE version.** The local build is **v3.0.0** (`(local HYPRE checkout)/CHANGELOG:10`; 32-bit ints, no OpenMP, no GPU; 00 §6.3).
HYPRE 3.0.0 **rewrote the Struct/SStruct code with API changes** and removed FAC/Maxwell/MSG (`(local HYPRE checkout)/CHANGELOG`,
3.0.0 entry). AMReX's own CI (`.github/workflows/hypre.yml`, development branch, fetched 2026-09-25) tests HYPRE
**2.21.0** (GCC 3D), **2.32.0** (GCC EB 2D, nodal EB 2D) and **3.1.0** (CUDA EB 2D), but **not 3.0.0**. The LinearSolvers
docs still say "tested with HYPRE 2.32.0" (as of March 2025).
**FireX coupling.** FireX's own HYPRE use is IJ/ParCSR only, but it now also calls the device APIs (`HYPRE_SetExecutionPolicy`,
`*_Migrate`, `*Initialize_v2`; 00 §11.1), which the local 3.0.0 library exports. If one HYPRE build is shared between the legacy
FDS solvers (reference runs) and AMReX, it must provide both sets of APIs. Any 3.x build does.
**[REC-G3]** For AMReX, build against **3.1.0** (the version AMReX's CI tests), or run an AMReX HYPRE smoke test against
the local 3.0.0 before relying on it. The IJ/ParCSR path used for bottom solves is lower-risk than SStruct
(`use_mlhypre`) **[OPEN Q6]**. For large meshes (>2³¹ rows) a `--enable-bigint`/mixedint build would be needed; not
relevant on the 8-core development machine.

---

## H. Verification plan

All of the following are proposed. **None has been run** (memory budget, A8). Directories are relative to
`Verification` (at `36975d765f`).

### H.1 Tier 0 — operator/solver unit tests (AMReX-only, small)
| Test | Setup | Metric / acceptance |
|---|---|---|
| MMS Poisson, composite | `φ = sin·sin·sin` on 2-level and 3-level grids (ratio 2), Dirichlet and Neumann | ‖φ−φ_exact‖∞ order ≥ 2 away from C/F, ≥ 1.5 in C/F cells (proposed) |
| Projection identity | random Ũ, random S | ‖∇·U − S‖∞ ≤ tol on every uncovered cell of every level after `project()`. Confirms the sign convention [VERIFY §A.1] |
| Free-stream preservation | uniform Ũ, S = 0, through a C/F boundary | U unchanged to 1e-12 |
| Singular / zones | closed box with 2 zones separated by a solid wall (E-2) | converges; the per-zone gauge reproduces the `pres.f90:3495-3548` shift |

### H.2 Tier 1 — FDS equivalence on single level (A2-a, E-1)
For the same input on one uniform level, AMR H and u should agree with FDS `SOLVER='GLMAT'` to solver tolerance.
Cases: `Pressure_Solver/dancing_eddies_1mesh`, `Pressure_Solver/pressure_iteration3d_default`,
`Pressure_Effects/isentropic`, `Pressure_Effects/pressure_rise`.
Metrics: H per §I.2 (ε_H = 1e-8, single frozen-state solve). Max |Δu|/U_ref ≤ 1e-8 after one step. The same number of pressure
iterations as FDS.

### H.3 Tier 2 — multilevel accuracy against analytic solutions
* `Adaptive_Mesh_Refinement/ns2d_16_emb_1to1_refinement`, `ns2d_16_emb_1to2_refinement`, `ns2d_16_int_1to2_refinement`
  (analytic 2-D Navier–Stokes, `VELOCITY_TOLERANCE=1e-6`). Refine by 2× twice and require second-order convergence
  of u in L2 and L∞. The error must not exceed that of the FDS run of the same case.
* `NS_Analytical_Solution/` (same analytic solution, single mesh), as the no-AMR baseline for the above.
* `Pressure_Effects/obst_coarse_fine_interface` versus its `_exact.csv`, to discriminate the AMReX C/F interpolation from
  the FDS two-point flux (§C.1).

### H.4 Tier 3 — FDS multi-mesh cases (compare with FDS results and with the 1-mesh reference)
| Case group | What it tests | Metric |
|---|---|---|
| `Pressure_Solver/dancing_eddies_{default,embed,tight,uglmat_refine,…}` | interface consistency vs 1-mesh | L2(u − u_1mesh) over time ≤ FDS `uglmat` value; max normal-velocity jump across C/F = 0 to tolerance |
| `Pressure_Solver/duct_flow*` | solids + multi-mesh + solver variants | mean flow rate vs FDS; wall normal velocity (E-2: 0; E-1: ≤ tol) |
| `Pressure_Solver/pressure_iteration{2d,3d}_{default,uglmat}` | iteration counts and residuals | iterations ≤ FDS; `PRESSURE_ERROR_MAX < PRESSURE_TOLERANCE` |
| `Pressure_Solver/tunnel_demo`, `tunnel_demo_glmat` | elongated domains (§D REC-D3) | MLMG V-cycles and wall time vs FDS with the preconditioner |
| `Pressure_Solver/stairwell` (11 meshes), `hallways`, `random_obstructions_fft` | many boxes, obstructions | global quantities (`_devc.csv`) within FDS run-to-run tolerance |
| `Pressure_Solver/obst_activation_*`, `Pressure_Effects/zone_break_*`, `zone_shape*`, `ulmat_2zone` | obstruction create/remove, zone merging, per-zone gauge | zone pressure histories vs FDS; no spikes at the matrix-reassign time |
| `Pressure_Effects/isentropic2`, `pressure_boundary`, `thick_orifice_5cm` | dP̄/dt, Dirichlet values, orifice | P̄(t) vs analytic/FDS |

Diagnostics to implement once and use everywhere: (i) max |∇·u − D| per level (uncovered), (ii) max solid normal velocity,
(iii) max C/F flux mismatch (coarse face − fine average), (iv) `PRESSURE_ERROR_MAX` composite, (v) MLMG iterations and
time per solve, (vi) per-zone Σ(∇·u − D)V (should be ≈0).

---

## I. Acceptance tolerance for H against the FDS reference solvers (answers FR-016(c) "TBD(Pressure Lead)")

### I.1 What UGLMAT/GLMAT actually solve to (FireX `36975d765f`)
* **Solver.** HYPRE ParCSR PCG, preconditioned by one BoomerAMG V-cycle, set up per zone in `GET_H_MATRIX_LUDCMP`
  (`SETMAXITER`/`SETTOL`/`SETTWONORM` at `pres.f90:4828-4830`; AMG coarsening/relaxation at `pres.f90:4837-4838`), solved at
  `pres.f90:3461`.
* **Constants** (`imkl.f90:238-240`, `244`, `251`, `261-263`): `HYPRE_SOLVER_MAXIT = 1000`, **`HYPRE_SOLVER_TOL = 1e-12`**,
  `HYPRE_SOLVER_SETTWONORM = 1`, PMIS coarsening, ℓ1-Jacobi relaxation, 1 sweep, AMG tolerance 0, 1 AMG iteration.
  ULMAT-HYPRE uses the same constants (`pres.f90:3021-3031`).
* **Stopping test** (HYPRE `src/krylov/pcg.c:416-419`, `480-484` in `(local HYPRE checkout)`, v3.0.0). With two-norm on and the default
  `a_tol = 0`, PCG stops when **‖r‖₂ ≤ 10⁻¹² ‖b‖₂**, i.e. relative residual of the assembled system, or after 1000 iterations.
* **No convergence check in FDS.** The iteration count and final residual are queried only if `CHECK_POISSON` and
  `HYPRE_SOLVER_SETPRINTLEVEL>0` (`pres.f90:3463-3466`). The print level is the constant 0 (`imkl.f90:241`), so this never happens,
  and a non-converged solve would go unnoticed.
* **Singular zones.** The RHS mean is subtracted (`pres.f90:3369-3406`) and the last unknown is pinned (`pres.f90:3447`). The
  solution is then shifted by the mass-weighted gauge (`pres.f90:3495-3548`).
* **FireX vs master.** FireX uses PMIS/ℓ1-Jacobi where master `ce1f659` uses HMIS/ℓ1-GS (00 §11.1). The two masters therefore
  converge along different paths and agree only to the tolerance. **Reference runs must be produced with the FireX
  `36975d765f` build.** They must run on CPU: the local HYPRE has no GPU backend, so `HYPRE_DEVICE_RUN` is moot.

### I.2 Hard number for discretisation-identical comparisons
Applies when the AMR discrete operator is **the same** as the reference's:
* single-level, multi-box AMR (A2-a, E-1) vs **GLMAT** (FR-002);
* single-level AMR with E-2 masking vs **UGLMAT-HYPRE**.

The comparison is on a **single frozen-state solve**: the same RHS, ideally taken from the FDS run's first step. It is not a
time-marched run, because D and F differ after the first step and chaotic growth then dominates.

Metric (the provisional rule, made precise):
`e_H = ‖(H_AMR − ⟨H_AMR⟩_V) − (H_ref − ⟨H_ref⟩_V)‖_{2,V} / ‖H_ref − ⟨H_ref⟩_V‖_{2,V}`. Here ⟨·⟩_V is the volume-weighted mean
over gas cells of the zone, and ‖·‖_{2,V} is the volume-weighted L2 norm over uncovered valid cells.

**Derivation.** "10 × the looser solver tolerance" (10 × 1e-12 = 1e-11) compares a *residual* tolerance with an *error*. For an
SPD system, rel. error ≤ κ · rel. residual. For the 7-point Neumann Laplacian (mean fixed) on N cells per direction,
κ ≈ λ_max/λ_min ≈ (12/h²)/(π²/L²) = 12N²/π² ≈ 1.2 N² (estimate for a uniform box). Each solver's worst-case error is
therefore 1.2N²·1e-12, and the difference of two solvers is at most twice that: `e_H ≤ 2.4·10⁻¹²·N_max²`.

**[REC-I1] Hard number: ε_H = 1×10⁻⁸** for finest-level grids with **N_max ≤ 64** cells per direction
(2.4e-12 × 64² = 9.8e-9). For larger grids use `ε_H = 2.4·10⁻¹²·N_max²`. This covers all pressure verification inputs listed in
§H: the `ns2d_16_*` cases have N ≤ 32; `pressure_iteration3d` has 16³ per mesh, 32 per direction overall.
Conditions:
* MLMG must reach **‖r‖₂/‖b‖₂ ≤ 1e-12**, the same as HYPRE. MLMG's own test is a max-norm, and
  ‖r‖₂/‖b‖₂ ≤ √N_cells·‖r‖∞/‖b‖∞. So in equivalence tests set `tol_rel = 1e-12` with `set_always_use_bnorm`, then *verify* the
  2-norm residual with `comp_residual` (`AMReX_multigrid_mod.F90:27`).
* If round-off prevents reaching 1e-12, use the measured residual r_m in the bound: ε_H = 10·κ·max(r_m, 1e-12) with κ = 1.2N_max².
* Expected in practice: e_H ≈ 1e-11 to 1e-10, since MG and AMG errors on smooth RHS are far below the κ bound. A result above
  1e-8 indicates an operator or BC mismatch, not solver noise.

### I.3 Coarse–fine cases vs UGLMAT-HYPRE (FR-016(c)): a solver-tolerance rule does not apply
* **The discretisations differ at the C/F faces.** UGLMAT couples coarse and fine cells with a two-point flux
  `A = min(A_int, A_ext)`, `1/(δx_int+δx_ext)` (`pres.f90:5223`, `5229`). The AMReX composite operator uses polynomial C/F
  interpolation (`set_maxorder`, fixed at 2 with a fallback to 3, §C.1). Their H differ by an O(h^p) truncation-error term that is independent of the solver tolerance.
  No tolerance-based threshold can be met, and none should be expected.
* **[REC-I2]** Use three criteria instead:
  1. **Interface velocity (FR-032).** `max |u_c − ⟨u_f⟩_A| / U_ref ≤ 1e-12` on every C/F face. `average_down_faces` makes this
     exact by construction (§C.2), so the target is machine zero.
  2. **Accuracy against the exact solution** (`Adaptive_Mesh_Refinement/ns2d_16_int_1to2_refinement`, `ns2d_16_emb_1to2_refinement`):
     ‖H_AMR − H_exact‖ ≤ 1.1 × ‖H_UGLMAT − H_exact‖, both mean-removed, and observed order ≥ 1.8 under 2× refinement (FR-014).
  3. **Direct AMR vs UGLMAT difference** as a sanity bound: `e_H ≤ 3×10⁻²` on the N_coarse = 16 cases. **Estimate**: second-order
     truncation for the analytic Fourier mode k = 2π/L gives a relative error of about (kh)²/12 = (2π/16)²/12 ≈ 1.3% per scheme,
     so the difference between two schemes is ≤ 2.6%. In general use ε_CF = (2π/N_coarse)²/6. To be calibrated on the first runs.
* **[REC-I3] maxorder (roadmap P2, FR-039).** Run criterion 2 at `setMaxOrder(2)` and `setMaxOrder(3)` (§C.1). Order 2 is
  accepted if its observed order at the C/F interface is ≥ 1.8. Otherwise multi-level runs use order 3, and single-level
  FFT/MLMG equivalence (ε_H, §I.2) keeps using order 2.

---

## R. References (verified bibliographic data)

* A. S. Almgren, J. B. Bell, P. Colella, L. H. Howell, M. L. Welcome, "A conservative adaptive projection method for the
  variable density incompressible Navier–Stokes equations," *J. Comput. Phys.* **142**(1):1–46, 1998.
  doi:10.1006/jcph.1998.5890.
* R. B. Pember, L. H. Howell, J. B. Bell, P. Colella, W. Y. Crutchfield, W. A. Fiveland, J. P. Jessee, "An adaptive projection
  method for unsteady, low-Mach number combustion," *Combust. Sci. Technol.* **140**(1–6):123–168, 1998.
  doi:10.1080/00102209808915770.
* M. S. Day, J. B. Bell, "Numerical simulation of laminar reacting flows with complex chemistry," *Combust. Theory Model.*
  **4**(4):535–556, 2000. doi:10.1088/1364-7830/4/4/309.
* D. F. Martin, P. Colella, "A cell-centered adaptive projection method for the incompressible Euler equations,"
  *J. Comput. Phys.* **163**(2):271–312, 2000. doi:10.1006/jcph.2000.6575.
* A. Nonaka, A. S. Almgren, J. B. Bell, M. J. Lijewski, C. M. Malone, M. Zingale, "MAESTRO: An adaptive low Mach number
  hydrodynamics algorithm for stellar flows," *Astrophys. J. Suppl.* **188**:358–383, 2010.
  doi:10.1088/0067-0049/188/2/358. (Precedent for a time-dependent base-state pressure coupled to a projection, analogous
  to FDS P̄(t).)
* L. Esclapez et al., "PeleLMeX: an AMR low Mach number reactive flow simulation code without level sub-cycling,"
  *J. Open Source Softw.* **8**(90):5450, 2023. doi:10.21105/joss.05450. (Full author list not transcribed.)
* W. Zhang et al., "AMReX: a framework for block-structured adaptive mesh refinement," *J. Open Source Softw.*
  **4**(37):1370, 2019. doi:10.21105/joss.01370.
* R. A. Sweet, "Direct methods for the solution of Poisson's equation on a staggered grid," *J. Comput. Phys.* **12**:422–428,
  1973 (as given in `Manuals/Bibliography/FDS_general.bib:7169-7176` and cited at `TechGuide/Momentum_Chapter.tex:366`).
* M. Vanella, R. McDermott, G. Forney, K. McGrattan, "Fire Dynamics Simulator: Advances in simulation capability for complex
  geometry," Fire and Evacuation Modeling Technical Conference (FEMTC) 2016, paper
  https://media.thunderheadeng.net/femtc/2016_d1-07-vanella-paper.pdf (title/authors checked from the PDF).
* S. Kilian, "The FDS pressure equation: Intuitive understanding and solution strategies," FEMTC 2020,
  https://media.thunderheadeng.net/femtc/2020_d3-10-kilian-paper.pdf (title/author checked from the PDF).
* FDS Technical Reference Guide (NIST SP 1018-1), LaTeX source at `Manuals/FDS_Technical_Reference_Guide/` @ `36975d765f` (identical to `ce1f659`).
* AMReX Linear Solvers docs: https://amrex-codes.github.io/amrex/docs_html/LinearSolvers.html. AMReX-Hydro Projections
  docs: https://amrex-fluids.github.io/amrex-hydro/docs_html/Projections.html. AMReX-Hydro source
  (`Projections/hydro_MacProjector.{H,cpp}`, `hydro_FFTMacProjector.H`) and AMReX source
  (`Src/Extern/HYPRE/AMReX_HypreMLABecLap.H`, `.github/workflows/hypre.yml`), development branches fetched 2026-09-25.
  Line numbers of fetched files refer to that snapshot and may drift.

Not cited because not verified in this session: Bell, Colella & Glaz (1989); Almgren, Bell & Szymczak (1996); Bell & Marcus
(1992). These are listed on the AMReX-Hydro docs page but I did not check their bibliographic data.
