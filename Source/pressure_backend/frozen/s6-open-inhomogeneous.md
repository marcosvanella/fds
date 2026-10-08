# S6: OPEN / inhomogeneous Poisson (Role 2)

Signatures match tip `bf0ba2d523` and current HEAD (`fold_boundary_data` single-level + composite overloads in `PressureIface.H` / `CommonLayer.cpp`).

## Role 1 wire recipe (URGENT — die lift)

**Rule: call `fold_boundary_data` whenever any of `BXS..BZF` is nonzero on a solved physical face — not only for OPEN.**

| Face kind (FDS → `pb::BC`) | Data in `BXS..BZF` | Fold?
|---|---|---|
| OPEN → `Dirichlet` | wall value `H_b` | **yes** if ≠ 0 (and still fold zeros are no-ops) |
| MIRROR / INERT impermeable / WALL Neumann → `Neumann` | `g = dH/dx_d` (increasing coord) | **yes** if ≠ 0 — **MIRROR top with nonzero BZS/BZF dies without fold** |
| PERIODIC → `Periodic` | none | never (data on periodic → `InvalidInput`) |

Homogeneous Neumann (`g = 0` on all Neumann faces) and homogeneous Dirichlet (`H_b = 0`) need no fold terms (constants stay 0 / null MultFabs). The TimeLoop die fires on `max|BXS..BZF| > 0` on solved faces — that is exactly when fold is mandatory.

### Call order (composite — vortex / two-level)

```cpp
// After PRHS is in each level rhs; BEFORE solve_pressure. Do not lift die until this runs.
pb::BoundaryData bd{};  // value[f]=nullptr, constant[f]=0 by default
for (int f = 0; f < 6; ++f) {
    if (p.bc[f] == pb::BC::Periodic) continue;
    // Map FDS face array → bd.constant[f] or bd.value[f] (FISHPAK signs below).
    // Neumann (MIRROR/INERT): g from BXS/BXF/BYS/BYF/BZS/BZF
    // Dirichlet (OPEN):       H_b from the same arrays
}
std::string msg;
std::vector<pb::BoundaryData> bds{bd};              // or one BoundaryData per level
std::vector<amrex::MultiFab*> rhs_levels = {/* writable MultiFab* per level, same BA/DM as levels[l] */};
pb::Status fs = pb::fold_boundary_data(p, bds, rhs_levels, &msg);
if (fs != pb::Status::Ok) { /* abort with msg */ }

pb::PressureResult r = pb::solve_pressure(p, o, ws);
```

Single-level:

```cpp
pb::Status fs = pb::fold_boundary_data(p, bd, *rhs_mf, &msg);  // rhs BA/DM == p.ba/p.dm
pb::PressureResult r = pb::solve_pressure(p, o);
```

**Signatures (unchanged since `bf0ba2d523`):**

```cpp
Status fold_boundary_data(PressureProblem const& p, BoundaryData const& bd,
                          amrex::MultiFab& rhs, std::string* message = nullptr);
Status fold_boundary_data(PressureProblem const& p,
                          std::vector<BoundaryData> const& bd_levels,
                          std::vector<amrex::MultiFab*>& rhs_levels,
                          std::string* message = nullptr);
```

`BoundaryData`: `std::array<MultiFab const*,6> value{}`, `std::array<Real,6> constant{}` — index `face_index(dir, side)`.

### Vortex_ring MIRROR top (smoke)

`PBZ=L MIRROR` → Neumann on z-hi. If FDS fills `BZF ≠ 0` (or any other face array ≠ 0), the S5 die fires until Role 1 folds that `g` into the rhs. Mapping MIRROR→`pb::BC::Neumann` alone is not enough.

### Fail without fold

Skip fold with nonzero Dirichlet/Neumann data → wrong φ (harness neg control: constant Hb, max|phi−Hb| ≈ |Hb|). Do not use `status==Ok` alone.

### Acceptance template

Fill-in checklist: `frozen/s6-open-coldflow-RESULT.md`. Backend golden: `ctest -R 'pb_comp_bcdata$'` (floor INERT + five OPEN, φ≈Hb).

---

## Already built (single-level)

| Piece | Status | Evidence |
|---|---|---|
| `fold_boundary_data(p, BoundaryData, rhs)` | Built | `CommonLayer.cpp`; signs in `PressureIface.H` and `m2-notes.md` §4 |
| Mixed open/closed faces (N/D/P per face) | Built | FFT / MLMG / HYPRE; `pb_mixed_faces_sweep`, `pb_hypre_mixed_faces_sweep` |
| Fold + mixed exact / MMS | Built | `pb_bcdata_exact`, `pb_bcdata_mms`, `pb_hypre_bcdata_exact` |

## Composite fold

`fold_boundary_data(p, bd_levels, rhs_levels)` folds into each level using that level's `h` and `Domain`.
`bd_levels.size()` is `1` (broadcast) or `levels.size()`. Covered cells may receive the fold; backends ignore covered rhs.

Cold-flow golden: floor INERT (Neumann z-lo) + five OPEN (Dirichlet), constant `H_b`, full-fine two-level → `phi ≈ H_b`.
Harness `mode=comp_bcdata` / ctest `pb_comp_bcdata` (MLMG and HYPRE).

L1 golden: periodic x,y, Neumann z with `g = a` on both faces given as per-level data, fine patch touching z-hi only
(layout 3) → `phi = a z` exact, removed mean at round-off. Negative controls (no fold on L1, on L0) must fail.
Harness `mode=comp_bcdata_l1` / ctest `pb_comp_bcdata_l1` (MLMG and HYPRE).

vortex_ring_re100_wall (2026-10-08): fold is correct on L1 (fine `h`, domain layer only, periodic skipped, signs as below);
all BXS..BZF were ~1e-14 there, so the fold is a no-op in that run. Drift came from the case IC (nonzero `W` on the
MIRROR top face → uniform deficit removed as the rhs mean every solve) and from L1 `DDDT` never being bound (driver).

## FDS sign convention (do not invent)

From FDS `pres.f90` (header ~57–58, application ~452–454, ULMAT ~1620–1679), confirmed A-56 fold sign check:

- Neumann data `g` is `dH/dx_d` along the **increasing** coordinate at both the low and the high face.
- Dirichlet data is `H` at the wall (not the ghost).
- Fold into unit-operator rhs `b` (no volume factor):

| face | rhs change |
|---|---|
| low Neumann | `+ g / h` |
| high Neumann | `- g / h` |
| Dirichlet (either side) | `- 2 H_b / h^2` |

Periodic faces take no data (`InvalidInput`). One-cell y drops its term (D-057).

## Role 1 coordination (die)

`TimeLoop.cpp` historically `die`s on `BXS..BZF /= 0`. **Role 2 does not edit TimeLoop.** Lift the die only after fold is wired for **all** nonzero face data (OPEN and MIRROR/Neumann).

Partial-face OPEN / Known pad remains the masked path (Q2/Q6); domain-face OPEN and domain-face inhomogeneous Neumann are this fold path.

## Out of scope here

- `face_gradient_composite`+mask (later).
- Editing TimeLoop from Role 2.
- 3-level partial-fine masked (PL DEFER; masked-branch-design.md).
- §11.4 true FDS H with nonzero walls: `frozen/s6-fold-fds-h-dump-gap.md` (do not block die lift).
