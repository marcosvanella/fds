# ADR-002: Time stepping across AMR levels — subcycling vs single global dt

| Field | Value |
|---|---|
| Status | Proposed (framing memo; leaning recorded, not decided). v0.2 records accepted sub-decisions (see "Accepted decisions") |
| Version | v0.2.1 (2026-09-25): owner decision uniform grids per level in AMR mode (stretched cases FDS-only); D-032 multi-level mean removal; standing true-residual check (FR-039); UGLMAT reference A-38 defined; provisional one-cell-direction setup; see "v0.2.1 changes". v0.2 (2026-09-25): accepted items D-006, D-007, pressure-solver selection and P2 rulings, refinement ratios and D-030, reflux and bit-parity scope, decomposition requirements; see "v0.2 changes" |
| Type | Framing memo |
| Date | 2026-09-25 |
| Owner | Chief Architect, with AMR Pressure Solver Lead |
| Deciders | Project owner; AMR Pressure Solver Lead; AMR Spec & Program Lead |
| Depends on | **ADR-003** (solid treatment in the pressure operator decides whether the pressure iteration survives); ADR-001 (C++ AmrCore driver assumed) |
| **Code base** | **FireX `36975d765f` (2026-09-24), local branch `AMReX`, this repository (read-only)** — the primary code base and acceptance baseline. All `Source/` line numbers are FireX. Manuals/Verification files cited are unchanged from merge-base `ce1f659`. |
| Other evidence | AMReX `99ddfda`; IAMR `0a63b79`, incflo `878827d`, PeleLMeX `5c21556`, ERF `78707d3` |

Cross-references: risks R-01, R-04, R-05, R-14, R-15, R-16, R-17, R-20, R-24, R-30, R-36; requirements FR-001, FR-002, FR-005, FR-010, FR-014, FR-016, FR-020..024, FR-030..039; decisions D-006, D-007, D-012, D-015, D-021, D-022, D-023, D-026, D-028, D-030; P2 `docs/amrex/p2-findings.md`, `prototypes/p2_pressure_mlmg/logs/check_results.md`; P1 `docs/amrex/p1-findings.md` §9; roadmap Phases 4 and 6, gate M4 (`docs/{risks,requirements,roadmap}.md`, v0.1). Pressure baseline memo: `docs/pressure/00-fds-pressure-baseline.md`; pressure mapping spec `docs/pressure/01-amr-mapping-spec.md` (§B, §B.2a, §D, §E, §G.1). Teammate docs cite merge-base `ce1f659` line numbers; FireX changed 22 files (+8,621/−847), so `main.f90`/`pres.f90`/`read.f90` citations there are shifted.

## v0.2.1 changes (2026-09-25)
- **Owner decision (2026-09-25):** AMR mode uses only uniform grids on each level; the 59 stretched-grid cases stay FDS-only (confirms D-030). Closes charter Q11 (c) and the "stretched grids as hard requirement" question (ADR-001 v0.3.7).
- **D-032 (Spec Lead to log):** mean removal on multi-level runs uses the composite volume-weighted mean over uncovered cells only, per pressure zone, with the D-028 fixed-point sum under the `makeFineMask` mask.
- "P3 must check the true residual at least once per run" is replaced by a standing runtime check in the production driver (FR-039).
- UGLMAT reference defined as A-38 (V&V Lead): `ns2d_16_int_1to2_refinement` at N=16/32/64 on the existing FireX binary; no new build; `emb_*` excluded. The GNU Build Chief action is withdrawn.
- Provisional (pending the P2 case): setup for a one-cell-thick direction.

## v0.2 changes (2026-09-25)
- New section "Accepted decisions (v0.2)": D-006 (M2a demo gate) and D-007 (Phase 6 scope rule); `amrex::FFT::Poisson` replaces porting `pois.f90` with per-step solver selection and eps_H agreement (D-021, FR-037, FR-039, D-012); MLMG `setMaxOrder(2)` confirmed by P2; refinement ratios {2,4} (FR-010) and D-030; mass/species reflux (D-023, FR-016) and bit-parity scope (D-022, confirmed here); decomposition requirements (FR-005, D-028); P2 rulings (mean removal, `average_down_faces`, HYPRE as bottom solver only, level-0 coarsening warning, indicative timings, UGLMAT comparison open).
- Placement: D-006 (a demo gate) and the decomposition requirements (FR-005) have no owning ADR; they are recorded here by default. Refinement ratios are a pressure-solver constraint and sit here too. D-009 and D-010 are recorded in ADR-003.
- Existing text brought in line: the composite solve is MLMG with HYPRE only as bottom solver (was "MLMG, or HYPRE"); FFT is used for any single uniform level, not only single-box runs; S1 thresholds follow D-007; S2 criteria follow D-006/D-012/D-022; S3 gating case follows D-015; HYPRE pin follows D-026; the Pressure Lead's HYPRE question is answered.

## Context

**FDS time step (verified in `main.f90`).**
- Explicit second-order predictor–corrector, with **two pressure solves per step**:
  - predictor: `main.f90:745-923`, `PRESSURE_ITERATION_SCHEME` at `:855`;
  - corrector: `main.f90:925-1127`, `PRESSURE_ITERATION_SCHEME` at `:1090`.
- Predictor sequence: density/species (`DENSITY`, `:780`) → exchange (`MESH_EXCHANGE(1)`, `:790`) → momentum fluxes (`VELOCITY_FLUX`, `:819`) → HVAC (`:822-830`) → wall BCs + divergence (`WALL_BC`, `DIVERGENCE_PART_1/2`, `:838-850`) → pressure → `VELOCITY_PREDICTOR` (`:863`).
- The corrector repeats this sequence and adds combustion (`:973`), particles (`MOVE_PARTICLES`, `:987`), radiation (`COMPUTE_RADIATION`, `:1034`) and OBST mass exchange (`:1074-1079`).
- **Pressure iteration** (`main.f90:1601-1745`):
  - Solve H per mesh with FFT (`PRESSURE_SOLVER_FFT`, `pres.f90:318`, Fishpak `H3CZSS` `pois.f90:187`), or ULMAT per mesh (`pres.f90:1423`), or GLMAT/UGLMAT globally (`GLMAT_SOLVER`, `pres.f90:3299`). FireX adds a GPU (`WITH_HYPRE_DEVICE`) path for the HYPRE-based solvers (`pres.f90:1177, 3454`).
  - Exchange `H/HS` + `FVX..` (CODE 5; `main.f90:1617,1636,1663,1687`), compute the velocity error, `MPI_ALLGATHERV` the maxima, and loop until `VELOCITY_TOLERANCE`/`PRESSURE_TOLERANCE` or the max iteration count (`:1726-1731`).
  - The iteration exists to reduce two errors: normal-velocity mismatch at **mesh interfaces**, and penetration at **solid (IBM-masked) boundaries** (User Guide `FDS_User_Guide.tex:9510`).
  - GLMAT removes the interface error. It "produces the exact same pressure solution as the 'FFT' solver would if the domain were one single mesh" (`:9452`), but it still iterates for solids (`:9651`). UGLMAT removes both: unknowns only in gas cells, exact normal-velocity match at mesh boundaries, and it "allows stretching and refinement" (`:9455`).
- **Divergence constraint** (`divg.f90`).
  - `DIVERGENCE_PART_1` (`divg.f90:22`) builds D from thermodynamic sources.
  - The zone integrals `DSUM/PSUM/USUM` are summed over **all meshes** with `MPI_ALLREDUCE` (`main.f90:2038-2040`), giving the background-pressure rate dP̄/dt per zone.
  - `DIVERGENCE_PART_2` (`divg.f90:1431`) then applies it. This is a global, same-time coupling of every mesh, done twice per step.
- **dt control is global and synchronous.**
  - `CHECK_STABILITY` (`velo.f90:3028-3209`) gives `DT_NEW(NM)` per mesh from CFL/VN/particle CFL. `DT = MINVAL(DT_NEW)` after `MPI_ALLGATHERV` (`main.f90:715,738-741,885-891`).
  - A cut on any mesh repeats the whole predictor, including the pressure solve (`CHANGE_TIME_STEP_LOOP`, `main.f90:774-897`).
  - FDS has no local time stepping today (pressure memo §3.1).
- **Sub-models with their own time logic.**
  - Particles already sub-step inside a fluid step so that no particle crosses more than one cell (`DT_P = DT/N_ITER`, `part.f90:1882-1887`).
  - Radiation updates intensities every `TIME_STEP_INCREMENT` steps and `ANGLE_INCREMENT` angle sets (`radi.f90:3855-3865, 4277`), so it counts steps.
  - The HVAC network is solved on rank 0 from vent states on all meshes, in both predictor and corrector (`main.f90:829, 980`).
- **Coarse-fine interfaces already exist statically.** Abutting meshes of different resolution exchange through `EXTERNAL_WALL` `IIO..KKO_MIN:MAX` ranges. Ghosts are copied from a coarser neighbour or averaged from a finer one (`velo.f90:514-547`, `wall.f90:321-330`). There is **no refluxing**: `EXTERNAL_WALL_TYPE%FVN/FVNS` ("flux-limited ∫ρY u_n", `type.f90:478-479`) are declared but have no other reference in `Source/` (verified with `rg`). All meshes share one dt.

**Prior FDS-developer intent (verified; weaker than reported).**
- `Manuals/FDS_Technical_Reference_Guide/samr_notes.tex` (327 lines) is tracked in git (unchanged in FireX). No other `.tex` file references it, and most of its sections are empty headings. **It contains nothing on time stepping or subcycling.** Its substantive content concerns static meshes of different resolution:
  - wall data structures for an embedded mesh (`:50-110`);
  - the H boundary condition at a coarse–fine interface. The coarse→fine direction takes the coarse cell value directly, with "no interpolation between other cells", then interpolates with the previous step's H (`:127-225`). The fine→coarse direction averages (`:227-305`).
- History (from `Utilities/Misc/GoogleCode_Commit_Log.txt`):
  - `samr.f90` was added "for development of structured adaptive mesh refinement" on 2013-12-20 (r17785, line 4493);
  - notes were edited in 2014 (r19785–r20163, lines 2424–2735);
  - `samr.f90` is **not present** in `Source/` (FireX or master).
- Reading: NIST explored SAMR in 2013–14 and stopped at static, same-dt interface coupling. No subcycling design exists to inherit.
- `Verification/Adaptive_Mesh_Refinement/` has 4 inputs:
  - `random_meshes.fds` (overlapping meshes, plume vent) is active (`FDS_Cases.sh:3`).
  - `ns2d_16_emb_1to1`, `ns2d_16_emb_1to2` (8×8 or 16×16 mesh embedded in a 16×16 mesh) and `ns2d_16_int_1to2` (a 16×16 fine mesh surrounded by 12 4×4 coarse meshes) are commented out (`FDS_Cases.sh:4-6`). They test the 2-D analytic NS solution (DNS, periodic) across embedded and abutting refinement.
  - The User Guide says an embedded mesh is **one-way** coupled: "the larger mesh receives no information from the mesh embedded within" (`FDS_User_Guide.tex:1112`; precedence by input order, `:1008`).
  - So there is no existing AMR or two-way refinement mechanism, only static multi-resolution with a global dt (consistent with R-24).

**AMReX precedent (verified in clones).**
- **IAMR subcycles**, with MAC sync (`NavierStokes::mac_sync`, `IAMR/Source/NavierStokes.cpp:1440-1731`) and a level sync projection (`NavierStokesBase::level_sync`, `NavierStokesBase.cpp:1927`). It is built on `AmrLevel` (`NavierStokesBase.H:66-68`), not `AmrCore`.
- **incflo does not subcycle.** It keeps one `m_dt` (`incflo_compute_dt.cpp:304-308`), each stage loops over all levels (`incflo_update_velocity.cpp:26`), and fine→coarse consistency comes from `average_down`/`EB_average_down` of the pressure gradient after a composite projection (`incflo_apply_cc_projection.cpp:513-516`). A code comment says a per-level value "will never change unless we use subcycling" (`incflo/src/utilities/io.cpp:749-750`).
- **PeleLMeX does not subcycle** ("the non-subcycling version of PeleLM", `Docs/sphinx/manual/Model.rst:17, 387-397`). Coarse-fine flux consistency is by averaging down face fluxes, and the state is averaged down each SDC iteration. It has a **closed-chamber** background-pressure algorithm (`m_closed_chamber`, `Source/PeleLMeX.H:2380`; `Model.rst:39`), the closest analogue to FDS pressure zones and dP̄/dt. Its developers moved *away* from PeleLM's subcycling partly because "Synchronization across levels ... is also far simpler" (`Model.rst:397`).
- **ERF subcycles** (`nsubsteps`, `ERF_TimeStep.cpp:246,274,286`). It is mainly compressible, so it is less relevant to a projection method.
- `amrex-tutorials/ExampleCodes/Amr/Advection_AmrCore` has a `do_subcycle` switch (`AmrCoreAdv.cpp:34-37,112`), so AmrCore supports both.

## Decision drivers (ranked)
1. **Correctness of the low-Mach constraint across levels.** ∇·u = D everywhere (FR-023), zone dP̄/dt from composite integrals (FR-034), and conservation (FR-020..024).
2. **Preserve FDS algorithm structure** (predictor/corrector, two pressure solves, dt rejection loop). That keeps V&V applicable (FR-001/002/003).
3. **Implementation risk and time to first multi-level result** (roadmap M3/M4).
4. **Cost:** AMR benefit (NFR-032) vs the fine dt imposed everywhere (R-14).
5. **Sub-model compatibility:** particles, radiation, HVAC, OBST create/remove, devices.

## Options

### Option A — Single global dt, all levels advanced together (no subcycling)
Every level advances with the same dt through FDS's predictor/corrector. Each of the two pressure solves is **one composite multi-level solve** (MLMG, HYPRE only as its bottom solver; `FFT::Poisson` while the hierarchy is a single uniform level, FR-039). Fine→coarse consistency comes from `average_down` of state and face fluxes (flux register in "same-dt" mode).
- **Pros**
  - Identical to FDS's current global-dt model (`main.f90:715,891`), so the step structure, `CHANGE_TIME_STEP_LOOP` and the dt rejection carry over unchanged.
  - Zone integrals `DSUM/PSUM/USUM` stay same-time. They only need composite volumes (covered coarse cells excluded; FR-034).
  - HVAC and radiation step counting are unchanged.
  - Particles need no level-dependent dt; FDS already sub-steps them for the fine CFL (`part.f90:1882-1887`).
  - A composite solve removes the interface velocity error that the iteration exists to reduce (`FDS_User_Guide.tex:9452,9651`). In AMR mode the loop survives for solids **only under IBM forcing on an all-cell operator** (pressure spec §E, E-1), not under masking (E-2) or EB (E-3), and survives for the **baroclinic lag** whenever the FDS constant-coefficient form is kept (A2-a, FR-036) (`docs/pressure/01-amr-mapping-spec.md` §D).
  - Precedent in incflo and PeleLMeX.
  - No sync projection is required for correctness, because all levels are at the same time (pressure spec §B.2a). What remains is coarse-fine flux consistency (reflux/average_down), which Option C builds in.
- **Cons**
  - Coarse levels take fine-level steps. Wasted work grows with level count and coarse-cell fraction (R-14). The size of that penalty depends on *where* the CFL-limiting cell sits. In fire plumes it is often in the refined region, which reduces the penalty, but that is **not yet measured** (Spike S1).
  - Two composite solves per step, times the outer iterations that survive (E-1 solids and/or A2-a baroclinic), may cost more than per-mesh FFT (R-01, R-22).

### Option B — Berger–Oliger subcycling with IAMR-style synchronization
Level l+1 takes r steps per level-l step. Afterwards come refluxing, a MAC sync, a level sync projection, and (for FDS) a dP̄/dt correction.
- **Pros**
  - Best asymptotic cost with deep hierarchies.
  - IAMR and PeleLM prove it for variable-density low-Mach flows.
- **Cons**
  - Every same-time global coupling in FDS breaks mid-step: zone integrals (`main.f90:2038-2040`), the HVAC network (`:829,980`), the global dt rejection (`:774-897`, which would become per-level with parent re-runs), and radiation's step-count cadence (`radi.f90:3855`).
  - The sync machinery is large (`mac_sync` alone is ~290 lines, `NavierStokes.cpp:1440-1731`) and lives on `AmrLevel`, which ADR-001's AmrCore driver does not provide.
  - Particles must be advanced per level and migrate between levels mid-step, which complicates two-way coupling (FR-051).
  - FDS's D ≠ 0 constraint with zone dP̄/dt needs PeleLM-style divergence reflux plus dP0/dt sync. PeleLMeX's authors left this approach for simplicity (`Model.rst:397`).
  - No FDS prior art: `samr_notes.tex` has no time-stepping content.

### Option C — Option A now, "subcycling-ready" data model (the leaning, refined)
Implement A, but from Phase 3 keep per-level old/new state, `t_old/t_new`, and flux registers with time weighting. Keep the level loop factored so that a later subcycled advance (roadmap Phase 6) changes the driver, not the kernels.
- **Pros:** A's low risk now; keeps B open without rework; flux registers are needed anyway for FR-020/021/024.
- **Cons:** some up-front complexity (two time levels per level, ~2× state memory for those fields; NFR-031).

### Option D — Global dt for flow; local sub-stepping only inside stiff sub-models
Chemistry (`MAX_CHEMISTRY_SUBSTEPS`, `read.f90:4644`) and particles (`part.f90:1882-1887`) already sub-step internally. Extend that idea to per-level sub-models only if profiling shows a bottleneck.
- Pros: zero impact on the projection.
- Cons: it does not address the flow CFL penalty. It is complementary to A/C, not an alternative.

## Leaning
**Option C:** a single global dt, no subcycling in the first multi-level phases. A composite MLMG solve for H (HYPRE only as bottom solver) at **both** the predictor and corrector pressure solves, with the per-mesh FFT + interface iteration replaced in AMR mode (FR-030); `amrex::FFT::Poisson` whenever the hierarchy is a single uniform level (FR-037, FR-039, D-021). Flux registers and per-level old/new state are designed in from Phase 3.

Synchronization corrections (divergence reflux and background-pressure dP0/dt sync) are **deferred to roadmap Phase 6**, with the **R-05 review gate at the end of Phase 4 (milestone M4)**. If FR-014/022/023 fail or pass only marginally there, Phase 6 is pulled ahead of Phase 5. This matches the AMR Pressure Solver Lead's written recommendation (**REC-B1**, B-1 no subcycling; `01-amr-mapping-spec.md` §B) and the Spec Lead's roadmap. Per spec §B.2a, under B-1 the deferred work is not a sync *projection* but coarse-fine flux/divergence consistency; **R-05 is mainly a subcycling risk**, and the M4 gate is kept as a check rather than an expected trigger. Note: one relay of that leaning said "Phase 2", but `roadmap.md` places sync projections in Phase 6 and the gate at M4. This memo follows `roadmap.md`.

**Dependency on ADR-003 (explicit):**
The pressure loop has three reasons to exist (spec §D): (1) mesh interfaces — **eliminated** by the composite solve in every case; (2) solids; (3) baroclinic lag.
- ADR-003 keeps FDS OBSTs on Cartesian cells (Option B) **with IBM forcing on an all-cell operator** (spec E-1, the Phase 4 first step per REC-E1): reason (2) survives, and `VELOCITY_TOLERANCE` now measures only solid-boundary penetration (R-15, FR-035).
- ADR-003 Option B **with solids removed from the operator** (spec E-2, `MLABecLaplacian` + β=0 faces/overset mask, the UGLMAT analogue): reason (2) is eliminated; wall normal velocity ≈0 by construction.
- ADR-003 Option A (EB, spec E-3): reason (2) is eliminated.
- Independently of ADR-003, reason (3) keeps the loop alive under A2-a (Phase 4 default, REC-A2; FR-036) and disappears only under A2-b. So "the iteration disappears in AMR mode" is true only for E-2/E-3 **and** A2-b.
- ADR-002 therefore cannot be accepted before ADR-003 (roadmap Phase 1 exit rule; README decision log).

Strongest evidence:
- FDS is already a global-dt code with same-time global couplings (`main.f90:715,891,2038-2040,829`).
- The two AMReX low-Mach codes written most recently (incflo, PeleLMeX) chose no subcycling.
- NIST's own SAMR attempt never reached time-stepping design.

Confidence: medium. The unknown is the cost penalty (R-14), which Spike S1 measures cheaply.

## Accepted decisions (v0.2, 2026-09-25)
These are accepted sub-decisions. The ADR as a whole stays Proposed until ADR-003 is decided (see above).

### Scope and gates
- **D-006, Phase 2 demo gate M2a (accepted).** Full FDS steps through the ADR-001 shim on `shunn3_32` (FR-001, 2-D), `csmag_32` (FR-001, 3-D) and `shunn3_4mesh_32` (FR-002, single-level composite solve vs GLMAT). Gate wording per D-022: bitwise kernel checks on frozen input against single-mesh FDS, full steps within tolerance. Periodic only; no walls, combustion, radiation or particles (in/out list: `spec-responses.md` (a)). No ADR owns demo gates; recorded here because it exercises S2's single-level composite solve.
- **D-007, Phase 6 scope rule (accepted; final scope needs the project owner).** Set by the S1 global-dt penalty on the NFR-032 case: ≤ 1.5× → sync corrections only (subcycling out of the roadmap); 1.5-3× → sync first, subcycling optional with go/no-go at M5; > 3× → subcycling pulled forward to start at M4. Rationale (Pressure Lead): with a global dt and a composite solve the FDS step is a MAC projection (`velo.f90:1603-1630, 1723-1751`), so no sync projection is needed.

### Pressure solver
- **`amrex::FFT::Poisson` replaces porting `pois.f90`** (D-021, FR-037). The solver is chosen per step from the current hierarchy (FR-039): `FFT::Poisson` whenever the hierarchy is a single uniform level, composite MLMG otherwise. D-021's validity limits apply: unstretched, one BC type per domain face; mixed open/closed faces keep parity with FDS FFT (whole-face Dirichlet plus iteration), not with GLMAT.
- **Shared gauge and agreement:** both paths use the same pressure gauge and agree within eps_H = max(1e-8, 2.4e-12·N²) (D-012), for single solves on frozen input only. Multi-step run outputs use the V&V T2/T3 classes.
- **MLMG uses `setMaxOrder(2)`** (FR-039), confirmed by P2 (`docs/amrex/p2-findings.md` §1-2; `prototypes/p2_pressure_mlmg/logs/check_results.md`):
  - Worst FFT vs MLMG mean-removed relative L2 difference: 7.0e-14 (Neumann, periodic and Dirichlet at 64³; Neumann at 128³ and at 34×18×32).
  - maxorder 3 differs from the FFT by 1.8e-3 on Dirichlet (Neumann and periodic are unaffected), so maxorder 2 is confirmed.
  - Coarse-fine observed order (max norm, C/F-adjacent cells): ≥ 1.92 with maxorder 2 (ratio 2 periodic and Dirichlet, ratio 4 periodic); ≥ 1.97 with maxorder 3. Maxorder 3 cuts the interface error by about 40 % but is rejected because of the FFT disagreement. FR-039's fallback (maxorder 3 on multi-level hierarchies) is not needed.
  - Composite-solve fluxes: the coarse C/F flux inside the composite operator equals the fine-face sum to 1.8e-14 relative.
- **P2 rulings:**
  - **Mean removal** before every singular solve is mandatory and uses the D-028 exact fixed-point sum. MLMG reported convergence on a nonzero-mean RHS while the true composite residual stayed at 2.6e-5·max|b| (p2 §2.3), so MLMG's converged flag is not trusted on singular problems.
  - **Multi-level mean (D-032, v0.2.1):** on multi-level runs the mean removed is the composite, volume-weighted mean over **uncovered** cells only, taken per pressure zone (pressure spec §C.3-C.4, REC-C1/REC-C2). It is never a per-level mean and never includes covered cells. The D-028 fixed-point sum is computed with the `makeFineMask` mask applied (`AMReX_MultiFabUtil.H:236-244`).
  - **Standing true-residual check (FR-039, v0.2.1):** on singular solves, after MLMG reports convergence, the production driver computes the true composite residual of the mean-removed problem (uncovered cells only, per D-032). It does so on the first solve of every run and after every regrid. If the residual exceeds the solver tolerance, the run warns and logs it. A debug option runs the check on every solve. P2 already demonstrated the failure mode (p2 §2.3).
  - **Coarse interface fluxes** from `MLMG::getFluxes` are first-order wrong (p2 §2.2) and are always overwritten with `average_down_faces` before the projection (pressure spec §C.2).
  - **HYPRE is used only as MLMG's bottom solver.** P2 matched the default bottom solver to 1.3e-15 with HYPRE v2.32.0-24-g63331f19c (the FireX pin, D-026). If HYPRE is ever used as the whole solver (diagnostics), `hypre.adjust_singular_matrix=1` is mandatory: without it BoomerAMG diverges on closed (singular) domains (p2 §8 item 1). On the GPU (NVIDIA) path MLMG's own bottom solver is the default, since our HYPRE builds are CPU-only.
  - **Level-0 coarsening:** a 34×18×32 level 0 coarsens only to 17×9×16 (p2 §7). FR-010 warns when level 0 cannot be halved at least 3 times in each direction. Preferred fix: remesh level 0; second: HYPRE bottom (CPU only).
  - **Timings are indicative only** (shared machine): at 64³ on 4 ranks, FFT about 0.004 s and MLMG about 0.03 s per solve (p2 §5). NFR-030 needs idle-machine measurements.
  - **Open, the only item keeping P2 open: UGLMAT reference (A-38, V&V Lead; v0.2.1).** Case `Adaptive_Mesh_Refinement/ns2d_16_int_1to2_refinement.fds` (twelve 4×4 coarse meshes around a 16×16 fine mesh, 2:1 interface, fully periodic) at N=16/32/64, with `SOLVER='UGLMAT HYPRE'` added to the existing `&PRES` line (line 22). It runs on the existing FireX binary `(local build directory)/firex-36975d7/ompi_gnu_rel/fds`, linked against the pinned HYPRE, so no new build is needed. Output: error norms of H and u against the analytic solution at each resolution. The P2 comparison case uses the same geometry (REC-I2 criterion 2). The overlapping-mesh `emb_*` cases are excluded for now.

### Grids on each level (owner decision, 2026-09-25; v0.2.1)
- AMR mode uses only uniform grids on each level. The 59 stretched-grid (TRN) cases stay FDS-only, which confirms D-030 (and IR-002's rejection of stretched meshes in AMR mode). Consequences: no stretched-grid pressure path is needed (`PoissonHybrid` unused; R-29/A-29 moot), and charter Q11 (c) (x/y-stretched meshes on the GPU) is closed.
- **One-cell-thick direction (provisional, pending the P2 case):** such a direction (2-D inputs such as `ns2d_16_int_1to2_refinement`, IJK(2)=1) needs `amr.ref_ratio_vect` with 1 in that direction, a `blocking_factor` of 1 there, and `LPInfo::setHiddenDirection` for MLMG (`AMReX_MLLinOp.H:89` setter; used at `:951` and `:1218-1219`). This answers FR-010's open blocking-factor item provisionally.

### Refinement ratios
- **Capped at 4 because of MLMG:** the coarse-fine tangential interpolation is order 3 and asserts ratio ≤ 4 (`AMReX_InterpBndryData.H:186-188`). Supported ratios are {2,4} (FR-010); AMR mode rejects other ratios at input time. No ADR owns refinement constraints; this one is pressure-driven, so it sits here.
- **D-030:** the 31 FireX V&V cases rejected by the {2,4} limit and the 59 stretched-grid (TRN) cases stay FDS-only references. Any needed in AMR mode are remeshed deliberately, never auto-converted.

### Coarse-fine reflux and bit-parity scope
- **Reflux (D-023, FR-016):** mass and species are refluxed at coarse-fine faces; FR-020/021/024 require conservation to round-off. FDS's interface treatment (ghost averaging, first-fine-cell sites, velocity restore) is not copied. P1 evidence: without reflux, tracer mass drifted by up to −4.5e-4 relative over 128 steps, fully explained by the coarse-fine flux mismatch (predicted vs measured to 1e-10 relative; p1 §9, §9.1).
- **Parity scope (D-022, confirmed here):** the comparison against multi-mesh FDS is a tolerance check only (T2/T3). Bit-parity is against single-mesh FDS, for explicit kernels on frozen inputs only, and requires replaying FDS's exact time-step sequence (last-step adjustment at FireX `main.f90:719`).

### Decomposition requirements (FR-005)
Not a time-stepping topic and not owned by ADR-001 or ADR-003; recorded here by default.
- Explicit stages are byte-identical across box split and rank count on frozen inputs (FR-005 (i), which also covers thread count for reduction-free stages).
- Zone sums (`DSUM/PSUM/USUM`), the mean-removal gauge, `FDS_AREA`/`AREA_ADJUST`, and setup-time areas and volumes use exact fixed-point sums, starting from per-box sums and computed domain-wide (FR-005 (ii), (v); D-028; the upstream `FDS_AREA` defect is confirmed by runs).
- The pressure solution agrees within eps_H across box split and rank count (FR-005 (iii)). P2 measured at most 2.4e-16 relative L2 for ranks 1/2/4 × `max_grid_size` 32/16, on every path; not bitwise (p2 §3).
- Runs are bitwise reproducible at fixed rank count, thread count and box layout with `OMP_DYNAMIC=false` (FR-005 (iv)). P2: 3 repeats were bitwise identical on every path (p2 §4).

## Rejected alternatives and why
- **B now:** it breaks four same-time global couplings, and its `AmrLevel`-based precedent does not fit ADR-001's AmrCore driver. It is re-evaluated at Phase 6 with measured data.
- **Keep per-mesh FFT + interface iteration across levels:** FFT needs one tensor-product block per solve and couples blocks only through iterated Dirichlet H (R-01). There is no AMR prior art, and it keeps the interface error that a composite solve removes.
- **Port `pois.f90` (Crayfishpak) per box:** replaced by `amrex::FFT::Poisson` (D-021).
- **MLMG maxorder 3:** better coarse-fine error (about 40 %) but differs from the FFT by 1.8e-3 on Dirichlet faces (P2), breaking the eps_H switch check.
- **HYPRE as the whole solver in production:** diverges on singular matrices without `hypre.adjust_singular_matrix=1` and was much slower in P2 (indicative 0.93-2.95 s vs about 0.03 s for MLMG at 64³, p2 §5); diagnostics only.

## Consequences and risks
- **+** Minimal change to the FDS step, V&V comparability, simpler particles/HVAC/radiation, and a clear Phase 4 target.
- **−** The fine dt applies globally (R-14). NFR-032 may be missed until Phase 6.
- **Risk: deferred sync corrections** (R-05). Under global dt this reduces to coarse-fine flux/divergence consistency and composite zone integrals (spec §B.2a); it is mainly a subcycling risk. Mitigation: flux registers from Phase 3, track the interface divergence residual from Phase 4, M4 gate as a check.
- **Risk: the meaning of pressure-iteration inputs changes** (R-15). Mitigation: document it, and warn when tolerances that only affect interfaces are set.
- **Risk: two composite solves per step are slower than FFT in uniform mode** (R-01/R-22). Mitigation: `FFT::Poisson` whenever the hierarchy is a single uniform level, whatever the box count (FR-037, FR-039). P2 indicative: FFT about 0.004 s vs MLMG about 0.03 s per solve at 64³ on 4 ranks.
- **Risk: no conservation at coarse-fine faces today** (FVN/FVNS unused; R-04). Mitigation: flux registers from Phase 3 (FR-024), enforced by FR-020/021 budgets; mass and species reflux accepted (D-023; P1 drift without reflux up to −4.5e-4). FR-016 intentionally departs from baseline on the coarse side.

## Open questions for teammates
- **AMR Pressure Solver Lead:**
  - *Answered* (spec §G.1): cell-centred MAC projection, not nodal; B-1 no subcycling (REC-B1); A2-a then E-1 → E-2 (REC-A2, REC-E1).
  - *Answered* (v0.2, P2): MLMG is the default, HYPRE only as its bottom solver; HYPRE as the whole solver is diagnostics-only. On the GPU path MLMG's native bottom solver is the default because our HYPRE builds are CPU-only; FireX's `WITH_HYPRE_DEVICE` path (`pres.f90:1177`) stays a uniform-mode FireX option (FR-038).
  - Please confirm "Phase 6 / gate M4" rather than "Phase 2", and re-anchor spec citations to FireX (e.g. spec cites `main.f90:1478-1482`, which is `ce1f659` numbering for the `ITERATE_BAROCLINIC_TERM` set-up now at FireX `main.f90:1609-1613`).
- **AMR V&V Lead:**
  - References for the spikes: `SOLVER='GLMAT'` baseline for single-level multi-box runs (FR-002, same-resolution only, `:9452`); `SOLVER='UGLMAT HYPRE'` for coarse-fine cases, gating on `ns2d_16_int_1to2_refinement` only (`:9455`, FR-016, D-015; `emb_1to2` overlaps, so it is informational only).
  - A-38: run and own the UGLMAT reference (`ns2d_16_int_1to2_refinement`, N=16/32/64, existing FireX binary; see "Pressure solver"). `emb_*` excluded for now.
  - Do the three commented-out AMR cases run on baseline at all (A-09)?
- **FDS Legacy Mapper:** list every piece of step-count or `ICYC`-based logic (radiation increments, dump intervals, device averaging, `DT_RESTART`) that would need per-level semantics under subcycling, so Phase 6 can size it. Deliver `interface_averaging_sites.csv` (pending) for FR-016.
- **AMReX Integration Lead:**
  - Flux-register usage in same-dt mode (`FluxRegister` `CrseInit`/`FineAdd` with unit time weights), and EB variants if ADR-003 goes EB.
  - Can AmrCore host a later subcycled advance without an `AmrLevel`-style sync framework?
- **AMR Spec & Program Lead:**
  - Is a scope cut of Phase 6 to "sync corrections only, no subcycling" acceptable (roadmap Phase 6 note)?
  - Please pick the NFR-032 benefit case so that it is realistic under a global dt.
- **GNU/Intel build chiefs:** the FireX-pinned HYPRE `63331f19c` (v2.32.0-24; FireX CMake mislabels it "3.0.0") linked into AMReX as the MLMG bottom solver, on both toolchains (D-026; P2 linked it statically on GNU). (The earlier GNU Build Chief action for UGLMAT runs is withdrawn: A-38 uses the existing FireX binary.)

## Spike plan
| Spike | Scope | Pass / fail signal |
|---|---|---|
| **S1 Global-dt penalty estimate** (3 days, no new code) | Run baseline anchor cases (`random_meshes`, a plume, `dancing_eddies`) with `&DUMP CFL_FILE=T` (`main.f90:901`, CFL location logging `main.f90:1878-1884`). Record where the limiting cell sits relative to plausible refinement regions. Compute the estimated work of A vs B for 1–2 levels (ratio 2). | Phase 6 scope per D-007 (penalty on the NFR-032 case): ≤ 1.5× → sync corrections only, subcycling out of the roadmap; 1.5-3× → sync first, subcycling optional with go/no-go at M5; > 3× → subcycling pulled forward to start at M4. |
| **S2 Composite solve, single level, multi-box** (1 wk) | MLMG `MLPoisson` for H with FDS RHS on a same-resolution multi-box layout. Compare against baseline `SOLVER='GLMAT'` (FR-002); `shunn3_4mesh_32` is the M2a case (D-006). | `H` within eps_H of GLMAT on a single solve on frozen input (D-012); explicit kernels bitwise vs single-mesh FDS, full steps within tolerance (D-022). The velocity error at box interfaces is at machine zero with no iteration. |
| **S3 Two-level static, same dt** (2 wk) | `ns2d_16_int_1to2_refinement` equivalent (gating, D-015; `emb_1to2` informational only) with composite solve + average_down + flux register (reflux, D-023). References: `UGLMAT HYPRE` baseline and the analytic solution. | FR-016 (fine-side ghost fill of transported/diffused fields bitwise at step 1; interface normal-velocity mismatch at machine zero; agreement with UGLMAT), FR-014 order ≥ 1.8, FR-020/021 budgets. |
| **S4 Zone integral** (3 days) | `zone_break_fast` with a refined patch and composite `DSUM/PSUM/USUM`. | dP̄/dt matches uniform-fine to T2. If it drifts, R-05 fires early. |

A failure of S3's conservation or S4's drift under global dt would argue for pulling sync corrections (Phase 6) forward, not for subcycling itself.

## Decision needed from the project owner
None blocking for the leaning itself. The project owner should, however, confirm the **Phase 6 scope** that D-007's rule produces once S1 reports (sync corrections only vs full subcycling), because it sets whether NFR-032 is achievable at all.

## References
FireX `36975d765f` `Source/{main,divg,pres,pois,velo,part,radi,hvac,type,read}.f90`; `Manuals/FDS_User_Guide/FDS_User_Guide.tex` (1008, 1112, 9452-9457, 9510, 9651); `Manuals/FDS_Technical_Reference_Guide/samr_notes.tex`; `Verification/FDS_Cases.sh:3-6`; `Verification/Adaptive_Mesh_Refinement/*.fds`; `Utilities/Misc/GoogleCode_Commit_Log.txt`. IAMR `Source/NavierStokes{,Base}.cpp`; incflo `src/incflo_compute_dt.cpp`, `src/projection/*`; PeleLMeX `Docs/sphinx/manual/Model.rst`, `Source/PeleLMeX.H`; ERF `Source/TimeIntegration/ERF_TimeStep.cpp`.
