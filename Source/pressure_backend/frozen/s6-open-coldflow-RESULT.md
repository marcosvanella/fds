# S6 inhomogeneous BC — Role 1 acceptance template (after die lift)

**Role 2 does not lift the TimeLoop die.** Normative wire recipe: `frozen/s6-open-inhomogeneous.md` (signatures = tip `bf0ba2d523` / current HEAD).

## Hard rule

**Fold whenever `max|BXS..BZF| > 0` on a solved face.** That includes:

- OPEN → Dirichlet (`H_b`)
- **MIRROR / INERT / WALL Neumann → Neumann (`g`)** — vortex_ring MIRROR top dies if `BZF≠0` without fold
- Not only the five-OPEN cold-flow

## A. Backend golden (no driver)

```bash
ctest --test-dir /path/to/pb-build -R 'pb_comp_bcdata$' --output-on-failure
```

| Check | Expect | Evidence |
|---|---|---|
| MLMG / HYPRE | Ok, φ≈Hb=1.25, neg control without fold ≈1.25 | |
| Exit | 0 | |

## B. Wire recipe (paste into driver pressure path)

```cpp
// BEFORE solve_pressure; AFTER rhs filled from PRHS
pb::BoundaryData bd{};
for (int f = 0; f < 6; ++f) {
    if (p.bc[f] == pb::BC::Periodic) continue;
    // Neumann (MIRROR/INERT): bd.constant[f] or value[f] = g from BXS..BZF
    // Dirichlet (OPEN):       bd.constant[f] or value[f] = H_b from BXS..BZF
}
std::vector<pb::BoundaryData> bds{bd};
std::vector<amrex::MultiFab*> rhs_levels = { /* per-level writable rhs* */ };
std::string msg;
auto fs = pb::fold_boundary_data(p, bds, rhs_levels, &msg);
// require fs == Ok, then:
auto r = pb::solve_pressure(p, o, ws);
```

Single-level: `fold_boundary_data(p, bd, *rhs_mf, &msg)`.

## C. Checklist (Role 1 fills)

| # | Item | Result | Notes |
|---|---|---|---|
| 1 | Fold wired for **any** nonzero BXS..BZF (OPEN **and** MIRROR/Neumann) | | |
| 2 | Die lift tip hash | | Role 1 only |
| 3 | vortex_ring MIRROR top: no die; solve Ok | | smoke was dying on B≠0 |
| 4 | Optional: floor INERT + five OPEN cold-flow φ≈Hb | | |
| 5 | Neg control / smoke: without fold still wrong if B≠0 | | |
| 6 | `pb_comp_bcdata` green on same backend tip | | |

## Pointers

- `frozen/s6-open-inhomogeneous.md` — full signs + signatures
- `harness/composite_modes.cpp` `run_comp_bcdata`
- §11.4 gap: `frozen/s6-fold-fds-h-dump-gap.md`
