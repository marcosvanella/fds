# Fine-level OBST/VENT wall cells that conform to the fine grid (design note)

Status: **design note for review, no code.** Author: Role 1 (AMR Data Layout). Co-reviewer: AMR Solid Phase Lead (wall-state questions are marked **[SP?]**).
Base: HEAD `c27dcd7ed9` with patch 0026 (`MESH_OBJ`), hidden fine walls (`27ba12cdbe`, `d5575b2c70`), and Role 2's burner work
(`ffdcfa8203`, `4e30615ef7`, `c27dcd7ed9`, plus `ce1eb49619` on `species-stage2`).

## 1. Problem and where we are today

- **Fine solids are already snapped at fine resolution.** `SideData` comp 0 comes from `pb::obst_list()` (ObstMask.H, FDS `NINT` at the
  level's dx, input XB from `fds_obst_xb` = `OB%X1..Z2`), and `FDS_FINE_SET_SOLID` writes it into the fine `CELL_INDEX`. So the fine
  pressure/transport masks already see fine-snapped OBSTs.
- **What is still coarse:** the fine *wall cells*. `FDS_FINE_WALLS_BUILD` (today) puts a hidden wall on a fine face only where the face lies on
  a **level-0 parent wall face**, and copies the parent's SURF, VENT and surface state. A fine OBST face that does not coincide with a coarse
  wall face (sub-coarse OBST, a step at fine resolution, a VENT edge that snaps differently) gets no wall, or the wrong SURF. Role 2's
  `fds_fine_walls_state` / `fds_fine_walls_un` inherit the same limitation (they select faces with an L0 SOLID parent).
- So the gap is narrower than "run FDS's OBST setup on the fine grid": the cells are right, the **wall inventory and its surface identity** are not.

## 2. Critique of the proposed direction

**(1) Run FDS's OBST/VENT setup on the fine box, with its own CELL/WALL.** Agree with the goal, not with the means. FDS's setup is
spread over `READ_OBST`/`PROC_OBST`/`READ_VENT` (read.f90, per-mesh snapping into `MESHES(NM)%OBSTRUCTION/VENTS`),
`INITIALIZE_MESH_VARIABLES_1/2` (init.f90: wall counting, `INIT_WALL_CELL`, edges, interpolation, Poisson set-up) and
`INIT_WALL_CELL` itself (`M=>MESHES(NM)`, `M%VENTS`, `M%OBSTRUCTION`, `EXTERNAL_WALL`, `SEARCH_OTHER_MESHES`). Re-running it means
`MESH_OBJ` in a dozen read/init routines, a per-box OBSTRUCTION/VENTS table, and fighting the code paths that assume an FDS mesh
(interpolated boundaries, edges, Poisson). That breaks D-084 and the "no init.f90" line the owner held for 0026.
**Recommendation:** a driver builder that does the same arithmetic on the global input lists:
fine-snap each OBST/VENT XB with FDS `NINT` (the ObstMask rule, already shared with the masks), and create a wall on every fine gas face whose
neighbour is a fine solid cell or the domain boundary, with SURF from `OB%SURF_INDEX(IOR)` and VENT overrides from the fine-snapped VENT
planes (the `INIT_WALL_CELL` VENT_SEARCH_LOOP logic: plane match + ALLOW_VENT + RADIUS). Storage stays `ALLOCATE_STORAGE` on `MESH_OBJ`. The
fine box keeps its own `CELL`; its `WALL` stays hidden (see §3, step C for making it live).

**(2) The coarse level keeps a snapped OBST only under covered cells, with flux matching.** Disagree with changing level 0.
Level 0 must stay FDS's own OBST treatment: `aeb2ebf6de` showed that masking L0 drops the MASS_FLUX burner velocity (div 0.075 → 4679), and
L0 bit-identity is a hard requirement. Covered L0 cells are overwritten by average-down anyway, so their OBST content does not reach the answer
except through the coarse pressure solve. Keep L0 untouched. The real risk is a **coarse–fine interface that crosses an OBST whose coarse and
fine snaps differ**: the two levels disagree on what is solid at the interface. Handle that with tagging (point 3), not with a coarse rebuild.
Flux matching at the interface is already there (ADV/DIF override lists); fine solid faces enter them with ADV = DIF = 0 (`ce1eb49619`).

**(3) Tagging keeps sub-coarse OBSTs at the finest level.** Agree, with a precise rule: tag every OBST whose coarse snap differs from its fine
snap (any extent, not just "thinner than a coarse cell"), grown by `NBUF` coarse cells, so that the interface never cuts it. That is the
existing `--tag-obst NBUF` (RegridAmrCore) narrowed to "snap-different" OBSTs. VENTs likewise (a VENT plane whose fine snap differs).
Note the conflict with **D-010**: stateful VENTs (`fds_stateful_vents`) carry a no-refine buffer that keeps them on level 0, while Role 2's burner
runs put the burner on L1 (`--tag-obst 1 --nofine -1`). This note assumes the burner case runs with the D-010 buffer off and the burner
statically refined (allowed by the D-010 static-hierarchy clarification). **[SP?]** confirm.

**(4) New fine wall state at regrid by injection from neighbour or parent.** Partly. For the surfaces in scope (INERT, fixed TMP_FRONT,
HRRPUA/MASS_FLUX burners, OPEN) almost nothing is *state*: `TMP_F`, `U_NORMAL`, `M_DOT_G_PP_ADJUST`, `ZZ_F`, `RHO_F`, `Q_CON_F` are recomputed
every stage from input + gas (FR-041a: "stateless faces are rebuilt from input at each regrid"). What is genuinely carried is small:
`T_IGN` (ramps), `HEAT_TRANS_COEF`/`TMP_F` as an iteration seed (NET_FLUX ADLOOP), `Q_RAD_IN` until fine radiation exists, and B2 wall-function
values (`U_TAU`, `Y_PLUS`). Inject those from the **parent wall** when one exists, else from the **nearest wall of the same OBST face and SURF on the
same level** (deterministic order: face key then index). Thick walls must **not** use injection: FR-041b is decided as Option A
(finest-ever records under a face key, split r² on first refinement). Until those records exist, keep the FR-041a freeze and the error stop
already in `FDS_FINE_WALLS_BUILD` for thermally thick SURFs.

## 3. Data flow (proposed)

```
FDS set-up (unchanged)                               level 0: FDS OBST/VENT/WALL as today (bit-identical)
fds_obst_list (OB%X1..Z2, OB%SURF_INDEX(-3:3), ALLOW_VENT, THIN)   <- extend fds_mesh_query export
fds_vent_list (VT%X1..Z2, IOR, SURF_INDEX, RADIUS/X0..)             <- new export
   | (C++ ObstList: once)                             tagging: snap-different OBST/VENT -> finest level (+NBUF)
bind_level / regrid:
   fine CELL solid  <- SideData (fine NINT snap)      (today)
   FDS_FINE_WALLS_BUILD v2 (fine snap, own inventory):
      faces: fine gas | fine solid  or  fine gas | domain boundary   (coarse-fine faces: none)
      SURF: OBST face SURF, overridden by a matching fine-snapped VENT; domain faces: DEFAULT/VENT
      parent PIW kept when a coincident L0 wall exists, else 0
      ALLOCATE_STORAGE (MESH_OBJ), BC/RDN/AREA from fine metrics
      state seed: parent B1 (PIW>0) or nearest same-SURF fine wall (PIW=0) for T_IGN, seeds, B2   [SP?]
per stage (after L0 s_wall_bc):
   step A/B: fine wall B1 evaluated on fine gas (GAS_VARIABLES + SURFACE_HEAT_TRANSFER; burner M_DOT from SURF)
   fds_fine_walls_state / _un: read the FINE wall's B1 (not the parent's)   -- same C interface (Role 2)
   HRR file: FDS_FINE_WALLS_SUMS unchanged (parent share only where PIW>0; PIW=0 faces have no L0 counterpart)
step C (later): walls live in the fine kernels (WALL_BC on fine boxes), replacing the override-list injection
```

Steps: **A** conforming inventory (hidden), **B** fine-evaluated state feeding Role 2's accessors, **C** live walls (separate decision).

## 4. Files, routines and estimated lines

| file | change | lines |
|---|---|---|
| `driver/fds_mesh_query.f90` | `fds_obst_surf` (SURF_INDEX(-3:3), ALLOW_VENT per OBST), `fds_vent_list` | +60 |
| `driver/fds_fine_level.f90` | `FDS_FINE_WALLS_BUILD` v2: solid-neighbour face scan, fine NINT of OBST/VENT lists, VENT match, PIW=0 path, state seed | ~+150 / −40 |
| `driver/fds_fine_level.f90` | `fds_fine_walls_state`/`_un`: read fine B1; burner `M_DOT_G_PP_ADJUST` from SURF (HRRPUA/ΔH, ramp, T_IGN) for PIW=0 | ~+60 / −20 |
| `driver/fds_fine_level.f90` | `FDS_FINE_WALLS_SUMS`: skip parent share when PIW=0 | +5 |
| `regrid_transport/RegridAmrCore.*`, `pressure_backend/ObstMask.H` | `--tag-obst-snap NBUF`: tag snap-different OBST/VENT; `snap_differs(xb, dx0, dxf)` | +50 |
| `driver/TimeLoop.cpp` | none beyond today's build call (state accessors keep their signature) | 0–10 |
| upstream (patch) | none planned. Option: one line in `INIT_WALL_CELL` (`M=>MESH_OBJ`) if we ever reuse it (not recommended now) | 0 (1) |
| tests | `run_fine_obst_walls.sh`, two cases (§7), notes | +80 |

Total ≈ 400 driver lines, no new upstream lines. Step C (live walls) is not in this count; it is its own note (fine `WALL_BC` with
`MESH_OBJ`-style redirection of the wall loops, ≥ 300 lines and several upstream touch points).

## 5. Rebasing on the upstream WALL/ONE_D layout change

- Upstream touch points stay at 0026's four lines (`M=>MESHES(NM)` selection in `ALLOCATE_STORAGE`, `INITIALIZE_BOUNDARY_PROP1`,
  `HEAT_TRANSFER_COEFFICIENT`; `SURFACE_HEAT_TRANSFER` public). They rebase by re-applying the same pointer selection wherever the new layout
  allocates or reads a wall's mesh.
- All driver access to wall data goes through three routines in `fds_fine_level.f90` (build, state, sums), via `WC%BC_INDEX/B1_INDEX/B2_INDEX`
  and named `B1` fields. If the layout becomes flat SoA (NFR-048 a), only those three bodies change; the C interfaces (`fds_fine_walls_state`
  frozen in `pressure_backend/frozen/fine-wall-state.md`) do not.
- No use of `BOUNDARY_ONE_D`/`THR_D` (thick walls refused), so the ONE_D part of the change does not touch this work.

## 6. L0 bit-identity and D-084

- Level 0 is not modified: no L0 OBST/VENT/WALL change, `MESH_OBJ` stays null outside the fine builder, the new exports are read-only.
  `--levels 1` and the level-0 meshes of a two-level run stay byte-identical (same proof as for 0026: single-level `--run` binaries and
  OMP=1 two-level TWO-LEVEL lines/devc/steps against the previous HEAD).
- Two-level results change only where they should: the fine wall inventory (HRR-file wall terms, Role 2's burner face fluxes) and, with
  `--tag-obst-snap`, the hierarchy. The tag rule is opt-in (default off) so the existing baselines stay valid.
- D-084: one builder, one face rule ("gas next to solid or domain"), one snap rule (FDS NINT, shared with the masks), no second copy of
  READ_OBST/INIT_WALL_CELL, no new upstream code. The thing we deliberately do not do is run FDS's mesh set-up on boxes.

## 7. Test plan

1. **simple_test burner fine** (`--tag-obst 1 --nofine -1 --patch 0 7 0 7`, D-010 buffer off): wall inventory equals the L0 inventory
   refined (each L0 face ↔ 4 fine faces, same SURF, PIW>0 everywhere, since this OBST is coarse-aligned); `fds_fine_walls_state` output
   with fine-evaluated B1 equals today's parent-copied output for `U`, `M_DOT_G_PP_ADJUST`, `ZZ_F`, `TMP_F` to round-off (burner mass flux is
   input-driven); Role 2's `run_fine_burner_check.sh` mass/species gate numbers unchanged; L0 bit-identity.
2. **simple_test with a sub-coarse burner** (burner XB shifted by half a fine cell, e.g. 0.05 m): coarse and fine snaps differ; check PIW=0
   faces exist, burner total mass flux on the fine faces = HRRPUA·A_fine/ΔH within 1e-12 relative, and that `--tag-obst-snap` keeps the
   interface off the OBST.
3. **pyramid with fine-conforming steps** (`pyramid_48_4m.fds` variant with steps of half the coarse cell): counts of fine wall faces per IOR
   equal an FDS single-mesh run at fine dx over the refined region (exact integer match); Q_CON_F sum on the fine walls within resolution of
   the FDS fine run; closed-domain mass PASS; L0 bit-identity.
4. Regrid test: same case with `--regrid-every 5`: inventory rebuilt each time, seeds (T_IGN) carried, no change in L0.

## 8. Relation to Role 2 (burner wall state)

Role 2 owns the consumers (`fds_fine_walls_un`, `fds_fine_walls_state`, the fine pressure RHS and ADV/DIF override lists). This work owns the
inventory and the B1 values behind them. To avoid conflicts: (a) keep both C signatures frozen; (b) change only the *source* of the values
(parent B1 → fine B1) behind a single switch inside `fds_fine_level.f90`, landing after `ce1eb49619` is merged; (c) the selection rule
"parent is an internal SOLID wall with nonzero U_NORMAL/M_DOT" becomes "fine wall with SURF mass flux", which lists the same faces on
coarse-aligned burners (test 1) and adds the PIW=0 faces. `s[2]` (area weight to the parent) is 0 when PIW=0; Role 2's code must not divide by it.

## 9. Questions for the Solid Phase Lead [SP?] and the owner

1. **[SP?]** Is the list of carried quantities for stateless surfaces complete (T_IGN, HTC/TMP_F seed, Q_RAD_IN, B2 U_TAU/Y_PLUS)? Anything
   for HRRPUA with `TAU_Q`/ramps keyed on `T_IGN` beyond copying it?
2. **[SP?]** For a PIW=0 face (no coarse counterpart), is "nearest same-SURF wall of the same OBST face" acceptable as the seed, or should it
   come from the FR-046 face key once it exists?
3. **[SP?]** D-010 vs the burner on L1: confirm that a statically refined burner (hierarchy fixed around it) is within the D-010 clarification,
   and that the no-refine buffer may be off for these cases.
4. **[SP?]** Thick walls stay frozen (error stop) until Option A records exist. Agree that this note must not introduce any injection path
   that a thick wall could later fall into?
5. Owner: confirm level 0 is never re-snapped (point 2 rejected), and that `--tag-obst-snap` is the sanctioned way to keep the interface off
   snap-different OBSTs.
6. Owner: step C (live fine walls in the kernels) is a separate decision; this note stops at hidden, fine-evaluated walls feeding the
   existing override/RHS paths.

## 10. Effort

| step | content | estimate |
|---|---|---|
| A | exports, builder v2, VENT match, snap tagging, inventory tests 1–3 (counts) | 1.5 days |
| B | fine-evaluated B1, burner M_DOT for PIW=0, seeds, switch behind Role 2's accessors, tests 1–2 physics, regrid test 4 | 1.5 days |
| review | Solid Lead questions, Role 2 merge coordination | 0.5 day |
| **A+B** |  | **≈ 3.5 days** |
| C (not proposed now) | live fine walls in WALL_BC/velocity BC/radiation | ≥ 2 weeks, separate note |
