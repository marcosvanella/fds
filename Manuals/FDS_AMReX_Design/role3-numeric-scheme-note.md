# Role 3 note: numeric scheme for multi-level transport and regrid

Status: draft for the project owner and the Chief Architect. No code. Spec v0.4.28, FR-010..013, FR-016, FR-020/021, FR-024. Source read at the FDS-AMReX branch (mass.f90, divg.f90, velo.f90, pres.f90, wall.f90, func.f90, main.f90). Anything not read in the source is marked **(assumption)**.

**Words used.** *Level*: one grid resolution; level 1 is finer than level 0 by a ratio of 2 or 4. *Coarse-fine interface*: the faces where a fine patch meets coarse cells. *Covered* coarse cell: one lying under a fine patch. *Face-centred*: stored on cell faces (FDS velocity U, V, W). *Ghost cells*: extra layers outside a block holding copies or interpolated values of neighbours, so a stencil can be evaluated at the block edge. *Restriction* (average-down): fine data to coarse, here a volume average. *Prolongation*: coarse data to fine. *Conservative interpolation*: children are filled so their volume-average equals the parent exactly. *Slope limiter*: caps the interpolation slope so new values stay inside the range of the coarse neighbours (no new extremes, no negatives). *Reflux*: after a step, correct each coarse cell next to the interface by (coarse flux minus sum of fine fluxes) because the two sides moved different amounts of mass.

## 0. What FDS actually does (verified)

- Time step: predictor then corrector, each followed by a pressure solve (main.f90:745-1127; ADR-002 context). Not a pure explicit scheme: each stage ends with an elliptic constraint.
- Mass and species: face value of rho*Z from a flux limiter (`GET_SCALAR_FACE_VALUE`, func.f90:1330, reads one to two cells either side of the face), flux = face value x face-normal velocity, update `ZZS = rho*Z - DT*RHS` (mass.f90:~447-456); rho is the sum of species (mass.f90:497). The corrector averages the two stages (mass.f90:615-625).
- Diffusion is separate: `DEL_RHO_D_DEL_Z` is a cell array built in divg.f90 from two-point face fluxes (`RHO_D_DZDX`, divg.f90:175) and read by the next mass update; it is one stage old (mass.f90:~395-401, `DEL_RHO_D_DEL_Z__0`).
- Energy has no conserved variable. Heat conduction, enthalpy diffusion, radiation and combustion enter the divergence field D (divg.f90 part 1); temperature is derived from rho, Z and the background pressure (e.g. wall.f90:~355). H solves for the velocity whose divergence equals D: right-hand side uses `DDDT = (D_target - div U_current)/DT` (divg.f90:1619, 1631), then `VS = V - DT*(FV + dH/dx)` (velo.f90:1625-1640). Momentum is in advective form, so velocity is not a conserved quantity; only its divergence matters.
- At a mesh interface FDS copies area-averaged neighbour data into the first ghost layer; the second layer is second-order only for equal cell sizes, otherwise a zero-gradient copy (wall.f90:351-386; inventory README correction). Ghost mass fractions are clipped to [0,1] and temperature is rebuilt from the ghost density (wall.f90:~348).

## 1. Time stepping across levels (global DT, no subcycling)

**Where the owner is right.** With one DT from the finest grid on every level there is no time interpolation of coarse boundary data, no storage of fluxes accumulated over substeps, and no fine-vs-coarse time mismatch. The interface problem shrinks to: both sides must use the *same flux* through each shared face in each stage.

**Where it is more than trivial.** (a) Each stage ends in a pressure solve that couples all levels, so "explicit" holds for transport only. (b) The coarse side does not automatically see the same flux. Its flux uses its own stencil and its own face velocity; the fine side uses ghost values interpolated from coarse. The mismatch is a real mass loss: the earlier prototype drifted tracer mass by up to -4.5e-4 relative in 128 steps, fully explained by this (ADR-002, P1 evidence).

**Is reflux needed?** Case by case:

| Case | Verdict |
|---|---|
| Coarse face flux = area-sum of fine fluxes (product face value x velocity, per fine face), same DT, applied in each stage | Mass and species conserved to round-off. A time-accumulating reflux register is **not needed**. A post-step reflux would equal this overwrite algebraically. |
| Coarse side keeps its own flux | Mass drifts. Needs reflux (here a one-step correction, no accumulation). |
| Advection stencil crossing the interface | Affects accuracy (fine faces next to the interface use interpolated coarse ghosts, as FDS does), not conservation, once the coarse flux is overwritten. |
| Diffusion (two-point face flux) | Same overwrite on the diffusive face flux, **before** the cell divergence `DEL_RHO_D_DEL_Z` is formed. |
| Heat conduction / enthalpy diffusion (in D) | Not a mass issue. Same overwrite recommended so the integral of D matches the boundary flux (energy bookkeeping). |
| Radiation, combustion, particle sources | Cell sources, no interface flux. Restriction (average-down) handles them. |
| Pressure projection | Moves no scalar. The coarse face velocity must equal the area-average of the fine ones, which the composite solve provides. In Phase 3 (prescribed velocity) I set it that way by construction. |

**Decision for the plan.** The flux register / refluxing item (R3 as written) is **dropped**. It is replaced by an **interface flux overwrite**: per stage, (1) compute fluxes on all levels, (2) at each covered coarse face set the coarse flux to the area-sum of the fine face fluxes (product form, not product of averages), (3) then update cells; then (4) average-down covered cells. A negative control (overwrite off) remains as a test. Costs: coarse levels step r times more often than their own CFL needs (cheap, since fine cells dominate); coarse CFL is lower, which helps positivity. Positivity risk remains where the overwritten flux removes more than the coarse cell holds; counted and logged, handled by the existing clip policy (D-031, coarse side treats covered faces as walls).

**Spec impact (for the Architect).** FR-024 already allows "refluxing or equivalent". FR-016, D-023 and ADR-002 S3 say "flux register"; they should read "interface flux overwrite". Implementation hook needed from Role 1: an override input for the advective flux product between `MASS_FINITE_DIFFERENCES` and `DENSITY`, and for diffusive face fluxes inside divg. Fallback if kernels cannot change: correct the covered-adjacent coarse cells after the kernel and recompute rho and Y there.

## 2. Regrid: restriction and prolongation by variable class

Regrid happens only between full steps, so stage arrays (`RHOS`, `ZZS`, `US`, `HS`, `FV`) are not transferred; they are rebuilt at the next predictor.

| Class | Prolongation (new fine block) | Restriction (dropped or covered) | Conservation / positivity | Notes |
|---|---|---|---|---|
| rho, species rho*Z (cell) | Where fine data exists: copy. Else limited conservative linear interpolation of each rho*Z_i; rho = sum | Volume average (uniform cells, Cartesian) | Exact. Limiter keeps each rho*Z_i >= 0, hence Y in [0,1] and sum = 1 automatically | clips counted (FR-025) |
| Derived cell fields (TMP, MU, KRES, D-property, Cp) | Recomputed by kernels | Recomputed | n/a | not interpolated |
| Divergence field D, `DEL_RHO_D_DEL_Z` (history) | Recompute with one extra divg pass on changed levels | same | Not conserved quantities | Interpolating would carry sub-cell error into DDDT |
| Face velocity U,V,W | Coincident faces: coarse value; interior fine faces: divergence-preserving face interpolation (AMReX `FaceDivFree`), so each fine cell keeps its coarse cell's divergence | Area average of fine faces on the coarse face (`average_down_faces`) | Net face flux and cell divergence exactly preserved; momentum not conserved (not required in FDS); kinetic energy lost to averaging | Phase 3: analytic velocity re-evaluated on new faces, no transfer |
| Pressure H | Not needed (solved each stage); linear interpolation only as solver initial guess | Average, same use | n/a | background pressure: one table per zone, nothing to move |
| Ghost cells | Never transferred; refilled after regrid by same-level copy, then coarse-to-fine interpolation, then physical boundary fill | n/a | Conservative only through valid cells | Second layer: AMReX gives true interpolated data; FDS uses a zero-gradient copy at level jumps. Architect to confirm which is wanted **(assumption: FR-016 compares layer 1)** |
| Wall / boundary data (B1: RHO_F, ZZ_F, surface mass, solid profiles) | Per-area quantities copied from the parent wall face; extensive totals split by area | Area-weighted mean | Totals conserved by area weighting | Not in Phase 3 (periodic, no walls); Phase 5 (solid-phase owner) |

**Divergence on new fine faces.** The interpolation keeps the coarse divergence, but D varies inside a coarse cell. The next pressure solve removes the leftover because `DDDT` uses the actual current divergence (divg.f90:1619). Cost: a pressure impulse of size ~ dx^2 x error / DT. Recommendation: one extra projection right after each regrid (cheap at regrid intervals), default on once Role 2's composite solver exists; in Phase 3 not exercised.

**New blocks.** Fill order: copy overlapping old fine data (bit-exact), then interpolate the rest from the level below, then ghost fill, then recompute derived fields and D. **Dropped blocks.** Covered coarse cells already hold the fine average from the last average-down, so handing back is a no-op; a debug check compares coarse to the fine average just before the drop. Composite mass and species are unchanged at regrid to round-off (FR-012: <= 1e-12 relative).

**History.** Only end-of-step U,V,W, rho, rho*Z, D and the lagged diffusion cell array carry over. Pressure state is dropped (it is re-solved). Radiation intensities, particles and solid state are outside this note **(assumption)**.

## 3. Recommendation and open risks

Run one global DT from the finest level; replace refluxing with a per-stage interface flux overwrite for mass, species and (for D consistency) heat flux; spend the effort on the transfer table above, which is where the real risk sits, and on the divergence of velocity on new fine faces. Conservation tests: per-step composite budget, per-regrid budget, overwrite-off negative control.

Open risks:
1. Kernel hooks (flux override in `DENSITY` and divg) need Role 1; fallback is a post-kernel correction.
2. Overwritten fluxes can still drive a coarse cell negative; clip policy applies.
3. Leftover divergence after prolongation may excite H noise; measured in R5, extra projection is the remedy.
4. FDS's diffusion lag means regrid needs one extra divg pass.
5. Energy is not a conserved state in FDS; "energy conservation" is only via D consistency, not required by FR-020/021.
6. Not verified: stencil widths of CHARM/SUPERBEE limiters, exact `FaceDivFree` behaviour, radiation exchange at level jumps, cylindrical geometry, wall state.

## Addendum (after reading FDS `MATCH_VELOCITY`, `UVW_SAVE`, `COARSE_MESH_IF`): see `docs/role3-fds-mesh-interface-note.md`
- FDS matches interface velocity by averaging (`velo.f90:2630`) and keeps the unmatched value in `UVW_SAVE` for scalar transport and the DIVG boundary terms (`mass.f90:421-434`, `divg.f90:933,1260`). Our composite path has one face velocity, so both coincide; the K2 kernels keep working if `UVW_SAVE` is set to it.
- FDS already overwrites the coarse **diffusive species** flux with the fine area-weighted sum (`wall.f90:891-958`). It does **not** overwrite the advective flux, nor conduction. Section 1's "no reflux" decision stands; the advective overwrite is new, the diffusive one reproduces FDS.
- Section 2, face velocity: new fine faces on the interface take the coarse value (injection), interior faces `FaceDivFree`; the post-regrid projection is optional because `DDDT` (`divg.f90:1612-1632`) already absorbs divergence error at the next solve (default off in Phase 3).
