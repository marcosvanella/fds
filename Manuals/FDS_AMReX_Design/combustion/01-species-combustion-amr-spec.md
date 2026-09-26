# 01 — Species & Combustion on AMReX: spec outline

Owner: AMR Species & Combustion Lead · Status: **DRAFT (lean outline, low-spend mode)**, 2026-09-25 · Nothing here is approved.
Reference tree: this repository (FireX `36975d765f`), read-only. Line numbers are from that tree. Read first: `docs/README.md` (D-022, D-023), `requirements.md` (FR-016, FR-020/021/024), ADR-002 v0.2.
Scope for this pass: key citations, design direction, open questions, and candidate verification cases. There is no exhaustive survey yet. Items marked *(unread)* are pointers only.

## 1. Current FDS data flow (FireX)

| Stage | Where | Notes |
|---|---|---|
| Species advection/diffusion, predictor | `main.f90:767` `MASS_FINITE_DIFFERENCES`, `main.f90:780` `DENSITY` | Face fluxes, then update. |
| Corrector | `main.f90:948-949` | Same pair. |
| Combustion, once per step after the corrector | `main.f90:973` `COMBUSTION_LOAD_BALANCED` (`fire.f90:35`) | Skipped unless `N_REACTIONS>0 .OR. INIT_HRRPUV`. |
| Combustion BCs, then divergence | `main.f90:1055-1056` (`COMBUSTION_BC`, `fire.f90:1886`) | |
| Conservative update of `rho*Y` | `mass.f90:451` `ZZS = RHO*ZZ - DT*RHS` | Plus `M_DOT_PPP` at `mass.f90:477-478`, and an explicit `Q_Z` exchange between components 1 and 2 at `mass.f90:491-492` (context not yet read). |
| Density from the species sum | `mass.f90:509` `RHOS = SUM(ZZS(...,1:N_TRACKED_SPECIES))`, then `Y = rhoY/rho` at `mass.f90:531` | Density is not transported separately. It is the sum of the species densities. |
| Realizability | `CHECK_MASS_DENSITY` `mass.f90:775` (clip plus neighbour redistribution, `mass.f90:839`, `:915`); `CLIP_PASSIVE_SCALARS` `mass.f90:966-989` (clips `ZETA` to 0..1, `:984`; called `:540`, `:722`); `GET_REALIZABLE_MF` `func.f90:1556-1571` (clip each value to 0..1, then the most abundant species absorbs the sum error) | The redistribution is stencil-based and conservative within a mesh. The per-cell clip is not conservative. |
| Chemistry driver | `COMBUSTION_GENERAL_LOAD_BALANCED` `fire.f90:89` | Scans for active cells (`fire.f90:117-162`, `CHECK_CHEMICALLY_ACTIVE_STATE` `:365`). With load balancing on, it runs `MPI_ALLGATHER` of active counts (`:189`), ships cells to other ranks (`DISTRIBUTE_CELLS_ACCROSS_MPI_PROCESSES` `:196`/`:560`), calls `COMBUSTION_MODEL` per cell (`:214`), and gathers the results back (`:234`/`:692`). Otherwise it loops over meshes serially (`:246-295`), with OpenMP `SCHEDULE(DYNAMIC)` at `:274`. |
| Per-cell model | `COMBUSTION_MODEL` `fire.f90:771` | Integrators: `FIRE_FORWARD_EULER` `:1498`, `FIRE_RK2` `:1554`, `REACTION_RATE` `:1604`, `CVODE` wrapper `:1058`. Extinction: `EXTINCT_1/2` `:1389/:1422`, auto-ignition `:1349`. Source-term write: `SET_SPECIES_SOURCE_TERM_CELL` `:464`. |
| CVODE/SUNDIALS | `chem.f90` (guarded by `#ifdef WITH_SUNDIALS` at `:7`; `fire.f90:10, 332, 948, 1057`; `read.f90:5598`) | RHS `DERIVATIVE` `:92`, analytic `JACOBIAN` `:333`, per-cell `CVODE_SERIAL` `:753`, error handler `:1015`. |
| Energy and species sources in the divergence | `divg.f90:567, 579` (`DP += ... Q + QR`), `divg.f90:671` (`DP += D_SOURCE`) | Where `D_SOURCE` is formed *(unread)*. |

The key property is that everything in `COMBUSTION_MODEL` is **cell-local**. Its only non-local parts are the active-cell MPI shuffle and the stencil clip in `CHECK_MASS_DENSITY`.

## 2. What changes on an AMReX hierarchy

**Data layout.** `rho*Y` (the conserved variable) and `Y` are stored as `MultiFab` components with `N_TRACKED_SPECIES` components each. Passive scalars, including `ZETA`, are stored as extra components. The lumped-to-primitive matrix (`Z2Y`) and the property tables are read-only device data.

**Refluxing (D-023).** Register the **species face fluxes** in an `amrex::FluxRegister`. There is no separate density flux, because density is the species sum (`mass.f90:509`). Refluxing all tracked species therefore makes density conservative automatically and keeps the sum of `Y` equal to 1 without a separate step. After reflux and `average_down`, recompute `rho = SUM(rhoY)` and then `Y = rhoY/rho`. The FDS velocity restore (`mass.f90:424-436`) and ghost averaging (`wall.f90:319-339`) are **not** ported (D-023). The `M_DOT_PPP` and `Q_Z` terms are volume sources, so they need no reflux. Acceptance follows FR-020/021/024: round-off conservation.

**Realizability under interpolation and averaging.**
- `average_down` of `rhoY` with volume weights, followed by division by the averaged `rho`, gives a convex combination of fine-cell `Y`. It therefore stays inside 0..1 and sums to 1. No clipping is needed there.
- For coarse-to-fine fill (FillPatch and regrid), interpolate `rhoY` and `rho` conservatively, using `cell_cons_interp` with the min/max limiter, which is local extrema preserving per component. Then set `rho` to the sum of the interpolated `rhoY`. Open issue: per-component limiting can leave a sum-of-`Y` error at round-off level. The plan is to apply a `GET_REALIZABLE_MF`-style correction only in newly created fine cells and log how often it fires.
- `CHECK_MASS_DENSITY` redistributes mass across box boundaries through its stencil. On AMReX that needs ghost cells and a rule at coarse-fine faces. Options are (a) restrict redistribution to within one level, with a local clip plus a reflux-consistent correction, or (b) replace it with a flux limiter in the advection stencil, so that clipping becomes a diagnostic only. *Decision needed.*

**Source terms per box and per level.** Chemistry is a `ParallelFor` over the valid cells of each box on each level. Cells covered by a finer level are masked out and then overwritten by `average_down`. The FDS rank-to-rank cell shuffle (`fire.f90:186-236`) is replaced by AMReX load balancing, with a **chemistry cost `MultiFab`** (for example the CVODE step count per cell, from the previous step) feeding `DistributionMapping::makeKnapSack` or SFC. The owner allows this.

**Time stepping (ADR-002).**
- FDS calls combustion once per step, after the corrector (`main.f90:973`).
- With subcycling, each level integrates its chemistry over its own `dt_l`, at the same point in its step. The level's `Q` and `D_SOURCE` then feed that level's divergence.
- With a single global `dt`, there is one chemistry call per level per step.
- Either way, refluxing happens after the fine-level substeps, before `average_down`, and before the next coarse divergence. Also open is whether coarse chemistry in the region covered by the fine level is skipped or recomputed after averaging (the proposal is to skip it).

## 3. Refinement tagging (proposal)
Tag a cell if any of these holds (thresholds are inputs):
- heat release rate per unit volume `Q > q_tag`;
- the relative gradient of fuel or O2 mass fraction, or of mixture fraction, exceeds `g_tag`;
- `|grad T|*dx/T > t_tag`;
- the cell is chemically active per `CHECK_CHEMICALLY_ACTIVE_STATE` (`fire.f90:365`), which is optional.

Add a buffer of 2 or more cells. This must work with ADR-003's static wall refinement. Tagging must not use rank-dependent reductions (FR-005).

## 4. GPU execution of chemistry
- **Explicit paths** (`FIRE_FORWARD_EULER`/`FIRE_RK2`/`REACTION_RATE`) port straightforwardly as per-cell device functions. The hazards are derived-type pointer access (`REACTION%...`), the table lookups, and temporary allocations inside the per-cell path (for example `ALLOCATE(DZ_F0)` at `fire.f90:193, 252`). All of these must be hoisted or turned into flat device arrays.
- **CVODE.** The current path is a per-cell `CVODE_SERIAL` (`chem.f90:753`) called from Fortran. On NVIDIA GPUs there are three options:
  - (a) SUNDIALS batched solve: one `N_Vector` of size `ncell*nspec` per box, with block-diagonal Jacobians through `SUNLinSol_cuSolverSp_batchQR` or MAGMA dense batched solvers. This needs the SUNDIALS CUDA build.
  - (b) A custom per-cell implicit integrator in device code, for example BDF or Rosenbrock, as in the PelePhysics approach.
  - (c) Keep CVODE on the host for small mechanisms. This violates the "full step on GPU" owner decision, so it could only be a fallback.
  - The proposal is to prototype (a) and to keep the analytic `JACOBIAN` (`chem.f90:333`) as the device function.
- **Kernel style (pending the P1 readability review).** With C++ `ParallelFor`, the RHS and Jacobian must be ported to C++ or callable device code. With Fortran OpenMP target offload, `DERIVATIVE`/`JACOBIAN` could stay in Fortran, but a SUNDIALS batched solver called from an offloaded region is unproven. *Spike needed before the choice.*
- **Load imbalance.** The cost varies by orders of magnitude between cells. Two mitigations: the chemistry cost weight in §2, plus compacting active cells into a packed array per box before launching the kernel (on the GPU, the equivalent of the FDS active-cell scan).

## 5. Verification cases (proposal, to confirm with the AMR V&V Lead)
From `Verification/` (directories confirmed present):
- **Chemistry, 0-D reaction rates, and CVODE:** `Species/reactionrate_*` (EDC, Arrhenius with `_cvode` twins, series, fast_slow, lumped_two_air); `Chemistry/ign_delay_*`; `Chemistry/EDC_OneCFDStep_*`, `EDC_MultiCFDStep_*`; `Chemistry/EDC_load_bal_methane_smooke_{serial,parallel}` (this pair exercises load balancing directly).
- **Conservation and realizability:** `Species/mass_balance_reac*.fds`, `mass_balance_gas_volume.fds`, `bound_test_1/2.fds`, `lumped_stoich_*.fds`, `favre_test.fds`.
- **Flames and heat release:** `Species/burke_schumann.fds`, `methane_flame_{simple,lumped,primitive}*.fds`, `hrrpuv_reac_*.fds`, `Fires/simple_test.fds`, `circular_burner.fds`, `HoC_Ideal/NonIdeal.fds`, `tmp_lower_limit_*.fds`.
- **Extinction:** `Extinction/extinction_1/2.fds`.
- **New AMR variants (to create):** two-level versions of `methane_flame_simple`, `burke_schumann` and `mass_balance_reac`, with a fixed refined patch over the flame, then with dynamic tagging. Checks: FR-020/021/024 round-off, T2 against single-level runs, and 0..1 bounds every step.
- Note: `test-plan.md` §FR-020/021/024 cites `species_conservation_1..4`. These did not match under `Verification/Species` (`ls | rg species_conservation`), so the V&V Lead needs to confirm where they live.

## 6. Open questions and risks
1. **(Chief Architect)** `CHECK_MASS_DENSITY` redistribution across box and level boundaries: keep it level-local, or replace it with a limiter (§2)?
2. **(Chief Architect / Integration)** Skip chemistry under covered cells, and where combustion sits in the subcycled step (§2)?
3. **(Integration / Build)** SUNDIALS CUDA build and batched linear solver availability; interop with OpenMP-offload Fortran (§4).
4. **(V&V)** Location of `species_conservation_*`; acceptance criteria for AMR flame cases.
5. **(Chief Architect)** Scope couplings not yet read: HVAC species transport, species/mass bookkeeping in pressure zones (`divg.f90:1496-1499` uses zone pressures), level-set wildfire versus combustion, and `STORE_SPECIES_FLUX` (A-34).
- **Proposed risk (ID TBD):** device-side stiff chemistry (CVODE batched or custom) may be unavailable or slow for the chosen kernel style, which would block "full step on GPU" for finite-rate cases. Mitigation: an early spike, and the explicit paths first.
- **Proposed risk (ID TBD):** per-cell clipping after interpolation or regrid breaks FR-020/021/024 round-off conservation. Mitigation: convex averaging, limited conservative interpolation, and counting any clip that fires.
- **Proposed requirement (ID TBD):** mass fractions stay in 0..1 and sum to 1 within round-off after every FillPatch, regrid, `average_down` and reflux. Verification: a check each step in debug builds.
- **Proposed requirement (ID TBD):** the chemistry cost estimate feeds AMReX load balancing, and results do not depend on the distribution (FR-005).
