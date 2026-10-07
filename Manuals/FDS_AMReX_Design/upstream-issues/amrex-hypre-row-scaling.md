# Draft upstream issue for AMReX: HYPRE interface row-scales `A` to `D⁻¹A`, breaking PCG when diagonals vary

**State: DRAFT ONLY — ready for the project owner to paste into a GitHub issue. Not sent, not filed, nothing posted anywhere.** Do not open an AMReX-Codes/amrex issue from this tree. No patch file is attached yet; per project decision **D-075 (4)**, a candidate patch is delivered to the owner only after AMReX tests, and the AMReX tree is not edited from this project.

Evidence measured against AMReX 26.09 at commit 99ddfda (`SpMatrix: 32-bit row offsets in the split blocks (#5968)`), HYPRE v2.32.0-24, CPU and GPU harnesses on an NVIDIA test machine (see `docs/pressure/07-a56-backend-comparison-plan.md` §8, especially §8.2–§8.4). Numbers below are distilled from that section; do not treat them as a new measurement.

**Project context (not for the AMReX issue body):** FDS-AMReX Role 2’s assembled HYPRE backend builds its own IJ matrix (does not call `habec_ijmat`); the product workaround for Mf + HYPRE bottom on stretched grids is `hypre_solver=BiCGSTAB` with bottom tolerance 1e-11. This draft targets **AMReX’s** HypreABecLap / Habec path so PCG remains valid when diagonals vary.

---

## Title

`HypreABecLap3` / `habec_ijmat`: row scaling hands HYPRE `D⁻¹A` (unit diagonal), which is nonsymmetric when diagonals vary; PCG + BoomerAMG then fails or diverges on stretched grids

## Summary

AMReX's assembled HYPRE interface for the cell-centred ABecLaplacian (`AMReX_HypreABecLap3.cpp`, Habec kernels, `habec_ijmat`) scales every matrix row by the inverse of its own diagonal (`sten[m] *= diaginv`, centre entry set to 1) and multiplies the right-hand side by the same `diaginv`. For a symmetric discrete operator `A` this hands HYPRE `D⁻¹A`. That matrix has a unit diagonal, but it is **symmetric only when all diagonal entries are equal**. On a stretched (or otherwise non-uniform) grid the diagonals vary, so `D⁻¹A` is nonsymmetric.

The default Krylov method used with this interface in many applications (including FDS's BoomerAMG option set) is PCG, which requires a symmetric positive definite operator and preconditioner. With `D⁻¹A` and a varying diagonal, PCG fails or diverges at moderate stretch ratios, while the same unscaled matrix `A` and a symmetrically scaled matrix `D⁻½ A D⁻½` both converge with the same BoomerAMG options.

Affected: any call path that builds the HYPRE IJ matrix through `habec_ijmat` (and the matching RHS scaling) and then runs PCG, when the stencil diagonal is not constant (stretched cells, and in general non-uniform or masked operators). Uniform-spacing boxes with a constant diagonal are not affected: there `D⁻¹A` remains symmetric and PCG converges as today.

## Where in AMReX

- Driver / interface: `Src/LinearSolvers/MLMG/AMReX_HypreABecLap3.cpp` (and the related HypreABecLap path that uses Habec).
- Kernels: `Src/LinearSolvers/MLMG/AMReX_Habec_*_K.H` (3-D: `AMReX_Habec_3D_K.H`), function `habec_ijmat`.
- Mechanism (read): each row is scaled by `diaginv = 1/sten[centre]`; off-diagonals `sten[m] *= diaginv`; centre set to 1; RHS multiplied by the same `diaginv`. Result passed to HYPRE is `D⁻¹A` with unit diagonal.

## Minimal reproducer (conceptual)

No standalone one-file reproducer is attached to this draft. The measured isolation is:

1. Build the same stretched Neumann Poisson matrix (geometric stretch in z, last/first cell width = R; pure Neumann, one pinned cell) at 64³.
2. Solve with PCG + BoomerAMG (FDS-style options: PMIS coarsening, l1-Jacobi, strong threshold 0.25, interpolation 6, one AMG cycle per PCG iteration) on three matrices:
   - `A` (unscaled, as a direct HYPRE assembly would do),
   - `D⁻¹A` (what AMReX's `habec_ijmat` produces),
   - `D⁻½ A D⁻½` (symmetric scaling).
3. Separately, force `diaginv = 1` in a local copy of `habec_ijmat` and run MLMG with HYPRE as the bottom solver (PCG) on the same problem.

A short AMReX-side unit test suggestion is under **Test suggestion** below. Full tables, scripts and logs are recorded in this project's `docs/pressure/07-a56-backend-comparison-plan.md` §8 (especially §8.2–§8.4).

## Evidence (from pressure/07 §8; cite that doc)

Case: pure Neumann, pinned last cell, volume-weighted mean removal, 64³ unless stated, 4 CPU ranks, solve tolerance 1e-10. Diagonal max/min on the stretched operator: 2.0 at R=1 (boundary rows), 2.6 at R=8, **7.3 at R=64**, 38 at R=512.

### Isolation: same matrix, same AMG, PCG (assembled harness)

| Matrix handed to PCG + BoomerAMG | R=8 | R=64 | R=512 |
|---|---|---|---|
| `A` (unscaled) | 22 it | 27 it, 0.128 s | 39 it, 0.215 s |
| `D⁻¹A` (AMReX `habec_ijmat`) | 107 it, 0.52 s | **1000 it, not converged, 5.2 s** | **1000 it, not converged** |
| `D⁻½ A D⁻½` (symmetric scaling) | 28 it, 0.13 s | 38 it, 0.20 s | 60 it, 0.32 s |

So: the same matrix and options converge in ~27 iterations when the matrix is left unscaled (or symmetrically scaled) and fail when only the one-sided row scaling is applied. Symmetric scaling costs ~40% more iterations at R=64 but converges.

### MLMG + HYPRE bottom (Mf) behaviour with the shipped row scaling

At R=64 the first bottom PCG solve hits the 200-iteration cap: preconditioned residual starts at 0.42, falls to a minimum of 0.285, then grows to 26 by iteration 200. Later outer iterations recover; the first solve accounts for most of the cost (re-solve ~2.0–2.5 s vs ~0.17–0.21 s at R=1/8). At 100³ / 8 ranks, R=64 with PCG does not converge (200 V-cycles, residual ~2.8e20).

A **local** AMReX build forcing `diaginv = 1` in `habec_ijmat` (compile flag; not upstream, not committed) makes the Mf route converge in **3 V-cycles at every R** tested (R = 1, 8, 64, 512: 0.26, 0.28, 0.38, 0.55 s at 64³, 4 ranks).

Ruled out as the primary cause (same doc §8.2): BoomerAMG strength/coarsening on the unscaled matrix (27–39 it through R=512); raising the strength threshold; alternate coarsening/relaxation; `max_coarsening_level` / agglomeration (pinned cell already forces `mg_levels = 1`).

## Workarounds that work today (no AMReX change)

For any level whose diagonal varies (stretched cells; generally any non-uniform or masked operator):

1. Set `hypre.hypre_solver=BiCGSTAB` (GMRES / FlexGMRES also converge).
2. Tighten MLMG's bottom tolerance to `1e-11` (`MLMG::setBottomTolerance`).

With BiCGSTAB + bottom tolerance 1e-11, MLMG converges in **one V-cycle** at every R tested on 64³ (R=1, 8, 64, 512: ~0.19, 0.19, 0.28, 0.35 s). At 100³ the same setting restores convergence at R=64 on CPU and GPU where PCG diverges or is much slower. Keep existing BoomerAMG options.

On a **uniform**-spacing box the diagonal does not vary in this way; PCG remains fine and BiCGSTAB is only extra cost (do not force BiCGSTAB globally without need).

## Proposed fix

Make the HYPRE interface able to hand HYPRE a matrix that stays SPD-compatible with PCG when diagonals vary. Two options; both are small and local to `habec_ijmat` (order of ~5 lines for the skip path; slightly more for symmetric scaling plus matching vector scaling).

### (A) Optional skip of row scaling (preferred for minimal change)

- Leave `A` as-is (no `diaginv` multiply on the stencil; centre stays the true diagonal).
- Do not multiply the RHS by `diaginv` (it becomes multiply-by-1).
- Expose as a runtime / build option (e.g. ParmParse under `hypre.` or a method on the HypreABecLap object) so existing unit-diagonal users are unchanged by default.
- Measured with a local `diaginv=1` build only: Mf + PCG converges in 3 V-cycles at every R (see above). **Not** yet run through the AMReX regression suite.

### (B) Symmetric scaling `D⁻½ A D⁻½` with matching vector scaling

- Scale rows and columns by `1/√d_i` so the matrix HYPRE sees is symmetric with unit diagonal when `A` is SPD.
- Scale the RHS (and unscramble the solution) consistently.
- Measured on the assembled harness: converges (28 / 38 / 60 it at R=8 / 64 / 512) with the same PCG + BoomerAMG options; ~40% more iterations than unscaled `A` at R=64, still far better than row-scaled `D⁻¹A`.

**Preference:** (A) if the goal is “hand HYPRE the physical matrix and let the user’s Krylov choice be valid”; (B) if maintainers want to keep a unit diagonal for numerical reasons but still support PCG. Either should be **opt-in or at least not silently change the default for uniform grids** until tests say otherwise.

### Untested candidate sketch (pseudocode only — not a patch)

```text
// UNTESTED CANDIDATE — illustrative only; not applied to any tree; not run through AMReX tests.
// Location: habec_ijmat (AMReX_Habec_*_K.H) and the matching RHS scaling in HypreABecLap3.
//
// Option A (skip row scaling when requested):
//   Real diaginv = skip_row_scaling ? Real(1.0) : Real(1.0)/sten[centre];
//   for each off-diagonal m: sten[m] *= diaginv;
//   sten[centre] = skip_row_scaling ? sten[centre] : Real(1.0);
//   rhs *= diaginv;   // becomes no-op when skip_row_scaling
//
// Option B (symmetric scaling when requested):
//   Real s_i = 1/sqrt(d_i);  Real s_j = 1/sqrt(d_j) for column neighbour j;
//   a_ij := s_i * a_ij * s_j;  centre := 1;
//   rhs_i := s_i * rhs_i;  and x_physical_i := s_i * x_hypre_i after solve.
```

Caveats: (A) and (B) have **not** been exercised as a runtime option in upstream AMReX tests; the `diaginv=1` result is from a one-flag local build at 99ddfda. Singular / pinned components, overset masks, and GPU Habec paths need the same treatment if the option is added. No FDS- or project-specific code belongs in AMReX.

## Test suggestion

1. Add a stretched-z Neumann (or Dirichlet) ABecLaplacian case with geometric stretch R ≥ 64, HYPRE bottom or HypreABecLap with **PCG**, and assert convergence within a modest iteration cap (e.g. < 100 PCG it for the unscaled or symmetrically scaled path). With today’s row scaling this case fails or hits the cap; with (A) or (B) it should pass.
2. Keep an existing uniform-grid HYPRE + PCG test green with the **default** (current) behaviour, so the option does not change uniform-grid defaults without an explicit switch.
3. If (B) is chosen, assert symmetry of the assembled HYPRE matrix (max |a_ij − a_ji| at round-off) when the option is on.

## What we are NOT asking for

- Do **not** change the default behaviour for uniform grids without an explicit option (or without a clear deprecation path). Uniform-diagonal `D⁻¹A` + PCG works today and should keep working.
- Do **not** add FDS-specific or application-specific code to AMReX.
- Do **not** require dropping PCG as a supported Krylov method for this interface; the request is to make the matrix handed to HYPRE compatible with PCG when diagonals vary.
- This draft does **not** ask maintainers to accept an untested patch; a patch file for the owner follows only after AMReX tests (project D-075 (4)).

## Draft status

- **DRAFT ONLY.** Not sent, not filed, nothing posted to GitHub or elsewhere. Ready for owner paste.
- Written for the project owner to review and, if desired, paste into an AMReX-Codes/amrex issue (title + body from **Title** / sections above; omit the “Project context” paragraph if filing upstream).
- No patch attached yet. Candidate change (~5 lines in `habec_ijmat`, or symmetric scaling with vector matching) to be delivered as a **patch file for the owner** after AMReX tests; this project does not edit the AMReX tree for this item.
- Evidence source: `docs/pressure/07-a56-backend-comparison-plan.md` §8 (D-075 (4)). Refresh date: 2026-10-06 (Role 2 tip after R4/R5 ACCEPT).
- Checklist before filing (owner): (1) confirm AMReX tip still uses the same `habec_ijmat` row-scale; (2) optionally attach a minimal unit-test sketch; (3) still **no** untested patch in the first issue unless AMReX tests already exist.
