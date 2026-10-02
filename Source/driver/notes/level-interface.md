# Per-level interface of the driver (S9)

Answers Role 3 plan `Source/regrid_transport/notes/plan.md` section 4, items 1-5, `RegridInterface.H` and `flux-override-interface.md` (read-only inputs).
Header status: **stable** for the names listed in section 1 (`Level`, `LevelRegistry`, the `TimeLoop` per-level entry points, `BcStep::exchange` / `CfGhostHook`,
`exact_sum_uncovered` / `exact_sum_hierarchy`). Anything marked "not yet" is declared behaviour that aborts with a message instead of being faked.

## 1. What exists now (all checked, `tests/run_driver_tests.sh` test `levels`, and the unchanged single-level runs)

| Role 3 item | What the driver provides | Files |
|---|---|---|
| 1 hierarchy-aware registry | `Level` (= `Level0` alias, plus `level`, `ref_ratio_from_parent`, `fds_mesh_offset`, `fds_bound()`), `make_layout_level()`, `LevelRegistry` (implements `fdsrt::LevelListener`: one `Level`, `Fields`, `SideData`, covered mask per level; `make_level` / `remake_level` / `clear_level`; level 0 adopted from `TimeLoop`), `TimeLoop::registry()`, `TimeLoop::num_levels()` | `FdsAmr.H`, `LevelRegistry.H/.cpp`, `TimeLoop.H` |
| 2 per-level loop, global dt, exact sums | `advance()` runs every stage through `for_levels`; the stage bodies are public per-level entry points (`stage_viscosity`, `stage_density`, `stage_exchange`, `stage_boundary`, `stage_velocity_flux`, `stage_wall_bc`, `stage_divergence1/2`, `stage_velocity_update`, `stage_velocity_correct`); the per-box dt vectors hold every level (slot = level offset + box), `dt_at_step_start` / `dt_after_pass` take the MIN over all of them, and `global_dt()` reduces per-level values over ranks: one dt, the same in both stages (D-050). `exact_sum_uncovered`, `exact_sum_product_uncovered`, `exact_sum_hierarchy` (D-028 over uncovered cells, one fixed-point scale for the hierarchy, bitwise independent of box split and ranks) | `TimeLoop.cpp`, `ExactSum.H/.cpp` |
| 4 coarse-fine ghost hook | `CfGhostRequest` / `CfGhostHook`, `BcStep::cf_ghost_hook`, `BcStep::exchange(code, predictor)` (same-level fill, then the hook), `TimeLoop::set_cf_ghost_hook(level, hook)`. The OMESH route never sees a coarse-fine face (see 4 below) | `GhostExchange.H/.cpp` |
| 5 per-level SideData rebuild | `LevelRegistry::rebuild_side_data(level)` (called by make/remake; callable again by Role 3), `layout_cell_walls()` provider (domain edge = wall, every other face open, no solid), `set_side_provider()` for a fine-level FDS provider later, `covered_mask(level)` | `LevelRegistry.H/.cpp` |
| 3 face fluxes | not in this change (task (3): `notes/flux-hooks-design.md`) | |

Single-level behaviour is unchanged: `Level0` is an alias, the stage bodies are the former inline statements moved verbatim, level 0 is the only bound level, and the
byte-for-byte comparison of the final fields of the decomposition, shunn3 and csmag cases against the pre-change binary is part of the check list (README "S9").

## 2. What is per level and what stays global (level 0)

Per level (callable with a level index, same calls and order as `advance()`): viscosity + mass finite differences, density (D-031 clip over the level), exchange +
cf hook, boundary step (`after_exchange`), velocity flux, wall BC, divergence part 1/2, velocity predictor / corrector, dt contribution.

Global, level-0 only today (not faked for level > 0): the pressure scheme and the Poisson solve (`pressure_scheme`, composite or per-level solves are Role 2/3's
decision, the driver calls through `pressure_backend/PressureIface.H`), the pressure-zone sums (level-0 mesh data; D-053 keeps FDS order), the outputs and the
mass diagnostic (`mass_row` is a level-0 exact sum; it becomes `exact_sum_hierarchy` once a level is bound), set-up.

Order across levels inside one stage: all levels do the position (for example the exchange) before any level does the next one, coarse to fine. For the exchange
codes this means a fine level's cf hook runs after the coarse level's same-level fill of the same code and before any boundary routine of any level. If Role 3 needs a
different order (for example the coarse level's boundary step before the fine hook), say so: it is a one-line change in `advance()`.

## 3. What a level > 0 needs before a stage can run on it ("fine-level FDS mesh objects")

The kernels address a box by its FDS mesh number (`MESHES(NM)`, `fds_k_*(nm,...)`); `Level::fds_mesh_offset` maps box i of a level to mesh `fds_mesh_offset + i + 1`.
A bound level needs, per box, a `MESH_TYPE` object holding:
1. cell metrics `X Y Z XC YC ZC DX DY DZ DXN DYN DZN RDX RDY RDZ RDXN RDYN RDZN` (uniform: `dx_level`), `R = RRN = 1`, `IBAR/JBAR/KBAR`, `NM`;
2. cell table `CELL`, `CELL_INDEX` (all gas, no solid; obstructions at a fine level are not represented, plan item 5), `WALL`, `EXTERNAL_WALL`, `WALL_INDEX`: the domain-edge walls with
   FDS boundary data (`SURF_ID`, `BC_INDEX`, `ONE_D`), and the faces towards the coarser level or a neighbour box as interface walls (INTERPOLATED_BOUNDARY with no `OMESH`
   neighbour, so that `EXTERNAL_GHOSTS_FILLED` (patches 0003/0004, D-055) skips their `NOM > 0` branches); the driver's `iface()` machinery already turns every non-domain
   face of a box into a no-wall face around the kernels, whatever its neighbour (box or coarse level), so no new kernel logic is needed for the faces themselves;
3. the state arrays are the `Fields` of the level, bound by `fds_shim_bind` (S2 alias) or `POINT_TO_BOX` (patch 0005 DRAFT: not yet oneAPI-validated) per mesh number;
4. the `MESHES` array grows: two options for the Architect: (A) allocate `NMESHES_BASE + N_SPARE` mesh objects at set-up and activate them at regrid (a guarded patch in
   `main.f90`/`init.f90`, because mesh counts are used for `PROCESS`, `OMESH` and the MPI maps), or (B) keep the fine meshes in a second array and switch `POINT_TO_BOX` to it
   by level (patch 0005 grows by one module array; the kernels still receive the mesh number). Role 1 prefers (B): it does not touch the level-0 mesh count loops in
   `main.f90` (these assume NMESHES = level-0 meshes, MPI maps, `OMESH` allocation);
5. `fds_p_mesh_info`, `fds_p_zone_*` per level (the zone tables stay level-0 under D-053 until pressure zones under AMR are decided).
This is the single blocking item for physics on level > 0; nothing else in the driver is fixed to level 0 any more, except the items listed in section 2.

## 4. Coarse-fine faces and the OMESH route

`BcStep::fill_omesh` copies box data between boxes of the **same level** only; with `EXTERNAL_GHOSTS_FILLED` set, `VISCOSITY_BC` (L1399), `VELOCITY_BC` (L1367) and `NO_FLUX` (L1364) skip the `NOM>0` branches (D-055), so at a coarse-fine face FDS reads the ghost cells that Role 3's hook has filled (`UVW_SAVE` is saved before the match as at level 0). A level > 0 requires `ext_ghost` (`BcStep::exchange` aborts otherwise). `FDSTL_EXTGHOST=0` (the OMESH average) is a level-0-only legacy route.

## 5. Date estimates (calendar days from Monday 5 October, working days; estimates, not commitments)

| Item | State | Estimate |
|---|---|---|
| 1 hierarchy-aware registry | done, header stable | done |
| 4 coarse-fine ghost hook | done as a hook; exercised by a unit-test-level fake on level 0 and by Role 3's FillPatch when it arrives | done |
| 5 per-level SideData rebuild | done (layout provider); FDS provider for fine meshes follows item 2 | done |
| 2 per-level loop, exact sums, global dt | structure, dt and masked sums done; execution on level > 0 needs section 3 | 5-8 working days after the Architect picks (A) or (B) (about 9-14 October if decided on 5 October) [estimate] |
| 3 face fluxes out of the kernels + overrides | task (3) of this order; design note first, patch set 0007+ | design 6 October, hooks and the empty-set test 9 October [estimate] |

## 6. Checks added

`tests/test_units.cpp` `test_levels` (1 and 4 ranks): registry make/remake/clear, retired objects readable after remake, covered mask count, layout SideData at box-box, coarse-fine, periodic and domain faces, `exact_sum_uncovered` against a serial 128-bit reference, no-mask equality with `exact_sum`, `exact_sum_hierarchy` bitwise equal for 2 and 4 fine boxes.
