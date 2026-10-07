# S6: OPEN / inhomogeneous Poisson (Role 2)

Architect order after D-092/S6 ACCEPT. Item 5 (4-rank/3-level composite-masked) resumed separately.

## Already built (single-level)

| Piece | Status | Evidence |
|---|---|---|
| `fold_boundary_data(p, BoundaryData, rhs)` | Built | `CommonLayer.cpp`; signs in `PressureIface.H` and `m2-notes.md` §4 (FDS FISHPAK / `pres.f90`) |
| Mixed open/closed faces (N/D/P per face) | Built | FFT / MLMG / HYPRE; `pb_mixed_faces_sweep`, `pb_hypre_mixed_faces_sweep` |
| Fold + mixed exact / MMS | Built | `pb_bcdata_exact`, `pb_bcdata_mms`, `pb_hypre_bcdata_exact` |

## Composite fold (this tip)

`fold_boundary_data(p, bd_levels, rhs_levels)` folds into each level using that level's `h` and `Domain`.
`bd_levels.size()` is `1` (broadcast) or `levels.size()`. Covered cells may receive the fold; backends ignore covered rhs.

Cold-flow: floor INERT (Neumann z-lo) + five OPEN (Dirichlet), constant `H_b`, full-fine two-level → `phi ≈ H_b`.
Harness `mode=comp_bcdata` / ctest `pb_comp_bcdata` (MLMG and HYPRE).

## FDS sign convention (do not invent)

From FDS `pres.f90` (header ~57–58, application ~452–454, ULMAT ~1620–1679), confirmed A-56 fold sign check:

- Neumann `g = dH/dx_d` along the **increasing** coordinate at both low and high faces.
- Dirichlet `H_b` = wall value of `H` (not the ghost).
- Fold into unit-operator rhs `b` (no volume factor):

| face | rhs change |
|---|---|
| low Neumann | `+ g / h` |
| high Neumann | `- g / h` |
| Dirichlet (either side) | `- 2 H_b / h^2` |

Periodic faces take no data (`InvalidInput`). One-cell y drops its term (D-057).

Role 1 maps `BXS..BZF` → `BoundaryData` (`value[f]` MultFab or `constant[f]`) with this convention, then calls fold **before** `solve_pressure`.

## Role 1 coordination (do not edit TimeLoop here)

`TimeLoop.cpp` still `die`s on `BXS..BZF /= 0` (S5 guard). **Role 2 does not remove that die.**

To clear floor INERT + five OPEN on two-level vortex without abort, Role 1 must:

1. Set `p.bc[]` to the mixed open/closed map (OPEN → Dirichlet, INERT/WALL impermeable → Neumann, PERIODIC → Periodic).
2. Build `BoundaryData` from `BXS..BZF` (per level for composite MultFab data, or constants when uniform).
3. Call `fold_boundary_data` on each level's rhs (composite overload), then `solve_pressure`.
4. Only then drop the inhomogeneous `die`.

Partial-face OPEN / Known pad remains the masked path (Q2/Q6); domain-face OPEN is this fold path.


## Role 1: fold call site (pattern)

Copy the cold-flow harness (`mode=comp_bcdata` / `pb_comp_bcdata` in `harness/composite_modes.cpp`):

1. Fill `p.bc[face_index(d,side)]` (OPEN→`Dirichlet`, INERT impermeable→`Neumann`, PERIODIC→`Periodic`). Floor INERT + five OPEN is `bcfaces=DD,DD,ND`.
2. Put FDS `BXS..BZF` into `pb::BoundaryData` (`constant[f]` or `value[f]` MultFab). Signs: header comment on `fold_boundary_data` / table above.
3. Single-level: `fold_boundary_data(p, bd, *rhs_mf, &msg)` then `solve_pressure(p, o)`.
4. Composite: `std::vector<BoundaryData> bds{bd};` (or one entry per level), `std::vector<MultiFab*> rhs_levels` pointing at each level rhs, then `fold_boundary_data(p, bds, rhs_levels, &msg)` then `solve_pressure(p, o, ws)`.
5. Drop the TimeLoop inhomogeneous `die` only after (3)/(4) are wired.

Do not edit TimeLoop from Role 2; this note is the hand-off.

## Out of scope here

- `face_gradient_composite`+mask (later).
- Changing Role 1 TimeLoop `die` (report only).
- Item 5 resumed under Architect order after this ACCEPT (separate tip).
