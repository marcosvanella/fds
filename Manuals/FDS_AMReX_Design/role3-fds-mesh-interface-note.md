# What FDS already does at mesh interfaces, and what it means for our design (Role 3)

Status: analysis from reading the source and the Technical Reference Guide; nothing was run. Code citations are in `Source/`; "TRG" is `Manuals/FDS_Technical_Reference_Guide` (tracked in git and present in the worktree). **Verified** = read in code; **inferred** = follows from the code but not executed.

**Words.** *Interface*: faces where two FDS meshes touch (here: fine mesh against coarse mesh). *INTERPOLATED boundary*: FDS name for such a face. *Matched velocity*: the normal velocity on the shared face replaced by an average of the two sides. *Unmatched velocity*: the value each side had before that averaging. *Ghost cell*: copy of the neighbour's data just outside a mesh. *Divergence D*: the source term the pressure solve must reproduce as div(u).

## 1. Order of one FDS step (`main.f90`, two-stage predictor/corrector, SSPRK2-like)
Predictor (from line 821): `MASS_FINITE_DIFFERENCES` (face values of rho*Z, 837-841) -> `DENSITY` (rho, rho*Z update, 853) -> exchange 1 (RHOS, ZZS, DS, MU, KRES, Q two cells deep at interfaces, `mesh_exchange.csv`) -> `VELOCITY_FLUX` (momentum terms F, 892) -> `WALL_BC` (911) -> `DIVERGENCE_PART_1` (D from conduction, diffusion, sources, 913) -> `DIVERGENCE_PART_2` (DDDT, 923) -> `PRESSURE_ITERATION_SCHEME` (928) -> `VELOCITY_PREDICTOR` (936) -> exchange 3 (HS, US, VS, WS) -> **`MATCH_VELOCITY`** (987) -> `VELOCITY_BC`. Corrector (from 1002) repeats it with the new values: `MASS_FINITE_DIFFERENCES`+`DENSITY` (1019-1022), parts 1 and 2 (1129, 1158), pressure (1163), `VELOCITY_CORRECTOR` (1168), exchange 6, **`MATCH_VELOCITY`** (1193).

## 2. Matched and unmatched velocity (verified)
- `MATCH_VELOCITY` (`velo.f90:2630-2873`): for each interface face, own face value <- 0.5*(own + other), where "other" is the **area-weighted average of the abutting faces** of the other mesh (`UU_OTHER`, `DA_OTHER`). The coarse face averages the fine faces; the fine faces each average with the one coarse face. The value before the averaging is saved in **`UVW_SAVE(IW)`** (unmatched). Same-size faces also write the average into the local copy of the other mesh (`AREA_RATIO>0.9`).
- Matched velocity (U, V, W / US, VS, WS as stored): momentum terms (`VELOCITY_FLUX`), the divergence of the current velocity in `DDDT` (`divg.f90:1612-1632`), and the pressure update. `MATCH_VELOCITY_FLUX` (`velo.f90:2879`, called `main.f90:1711`) does the same averaging for the momentum flux FV before the pressure solve.
- Unmatched velocity (`UVW_SAVE`): scalar transport, because `DENSITY` copies U or US and **puts `UVW_SAVE` back on interface faces** (`mass.f90:421-434` predictor, `593-607` corrector); and the boundary terms of `DIVERGENCE_PART_1` (`divg.f90:933`, `1260`). So the owner's statement is confirmed. Interior terms of D use the stored (matched) face values (inferred from the pointers at `divg.f90:64-72`).
- Pressure iteration (`main.f90:1674-1818`): H is Dirichlet at interfaces; the solve is repeated (default `ITERATE_PRESSURE`, up to 10 times, tolerance 0.5 cell size) until the two sides' normal velocities agree.

## 3. Scalar fluxes at an interface (verified)
- Ghost data: the coarse mesh gets the area-weighted mean of the fine cells, the fine mesh gets the single coarse value injected (`wall.f90` `ASSIGN_GHOST_VALUE`, 282-388); second layer = copy (zero gradient) unless the faces have equal size. TRG `Mass_Chapter.tex:224-232`: first-order upwinding at refined boundaries.
- Advective: FDS "matches" the mass flux only through those ghost rules plus matched-by-iteration velocity. The coarse side uses (averaged rho) x (its own U), the fine side sum of (rho x U) per fine face: they differ by the covariance of rho and U over the face and by the leftover velocity mismatch (inferred). So FDS is not exactly conservative at a refined interface; our flux overwrite is **not** already there.
- Diffusive species flux: **already overwritten.** `wall.f90:891-958` (`COARSE_MESH_IF`, `EWC%NIC>1`) builds the coarse flux as the area-weighted sum of the fine two-point fluxes (`ARO`), fixes the species sum to zero (DNS/LES), and `divg.f90:204-216` writes it into `RHO_D_DZDX` before the divergence. TRG `Mass_Chapter.tex:234-240` describes it.
- Enthalpy diffusion at the coarse face = face enthalpy x that overwritten species flux (`divg.f90:~303`, inferred from order). Conduction `KDTDX`: no match (INTERPOLATED faces are skipped in `divg.f90:~540`; TRG says FDS "does not explicitly match" the heat flux).

## 4. Our design pieces against FDS
| Our piece | FDS equivalent | Change |
|---|---|---|
| Advective flux overwrite | none (ghost rules only) | keep; overwrite the **product** rho*U; use the unmatched face value `UVW_SAVE` as `DENSITY` does; in the composite path set `UVW_SAVE` = the single face velocity |
| Diffusive species overwrite | yes (`wall.f90:891-958`) | keep our override but reproduce this formula; FDS code reads `OMESH` data that the AMReX port does not have, so it must be bypassed at level interfaces (check) |
| Conduction/enthalpy flux | none for conduction | Phase 3: reproduce FDS, no override; optional later |
| Face-velocity transfer at regrid | interface face average (`MATCH_VELOCITY`) only for existing meshes; no regrid in FDS | new fine faces on the interface take the coarse value (area average then equals coarse); interior faces `FaceDivFree`; FDS 0.5-blend not used |
| Post-regrid projection | `DDDT=(D-div u)/dt` (`divg.f90:1612-1632`) corrects any divergence error at the next solve | optional, default off in Phase 3; run `DIVERGENCE_PART_1` on new levels |
