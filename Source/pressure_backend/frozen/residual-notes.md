# Residual limit follow-ups (A-67 / A-72, docs/pressure/08 §4.2)

E1 (MLMG half-tol polish) and E2 (mean-removed residual / floor) landed in commits `3252dc32c0` and `a2a4906647`. R6 harness in `a399b4fe9a`. This note records **E3**: composite floor measurement and the level-wise `||A||` change.

## E3: composite floor

### Issue

`residual_floor = 10 * 2^-53 * ||A|| * ||H||_w / ||b||_w`. On a hierarchy the norms use cell-volume weights over uncovered cells, but `||A||` was the **finest** level's `4 sum_d 1/dx_d^2`. Coarse mass then raises the floor by up to the cell-area ratio (4 for ratio 2, 16 for ratio 4).

### Measurement (before the formula change)

`pb_harness mode=hypre_resid`, n=32 coarse (and one n=16), mgs=16, tol_rel=1e-12, 2 ranks, tip before the E3 patch. `check` is the measured residual (round-off level when well converged); `floor/check` is the inflation of interest.

| bc | nlev | ratio | backend | check | floor | limit | floor/check |
|---|---|---|---|---|---|---|---|
| neumann | 2 | 2 | hypre | 8.42e-14 | 6.16e-13 | 1e-12 | 7.3 |
| neumann | 2 | 2 | mlmg | 1.68e-13 | 6.16e-13 | 1e-12 | 3.7 |
| neumann | 2 | 4 | hypre | 2.76e-13 | 2.46e-12 | 2.46e-12 | 8.9 |
| neumann | 2 | 4 | mlmg | 3.48e-14 | 2.46e-12 | 2.46e-12 | 71 |
| neumann | 3 | 2 | hypre | 1.31e-13 | 2.46e-12 | 2.46e-12 | 19 |
| neumann | 3 | 2 | mlmg | 1.69e-13 | 2.46e-12 | 2.46e-12 | 15 |
| neumann | 3 | 4 | hypre | 1.57e-12 | 3.94e-11 | 3.94e-11 | 25 |
| neumann | 3 | 4 | mlmg | 6.36e-14 | 3.94e-11 | 3.94e-11 | 620 |
| periodic | 2 | 2 | hypre | 3.87e-14 | 1.55e-13 | 1e-12 | 4.0 |
| periodic | 2 | 2 | mlmg | 8.89e-14 | 1.55e-13 | 1e-12 | 1.7 |
| periodic | 2 | 4 | hypre | 3.27e-13 | 6.21e-13 | 1e-12 | 1.9 |
| periodic | 2 | 4 | mlmg | 3.12e-14 | 6.21e-13 | 1e-12 | 20 |
| periodic | 3 | 2 | hypre | 1.63e-13 | 6.21e-13 | 1e-12 | 3.8 |
| periodic | 3 | 2 | mlmg | 1.11e-13 | 6.21e-13 | 1e-12 | 5.6 |
| periodic | 3 | 4 | hypre | 1.00e-12 | 9.93e-12 | 9.93e-12 | 9.9 |
| periodic | 3 | 4 | mlmg | 3.69e-14 | 9.93e-12 | 9.93e-12 | 269 |
| neumann | 2 | 2 (n=16) | hypre | 4.41e-14 | 1.55e-13 | 1e-12 | 3.5 |
| neumann | 2 | 2 (n=16) | mlmg | 2.69e-13 | 1.55e-13 | 1e-12 | 0.58 |

Verdict: floor sits well above ~3× measured round-off on several composite cases (worst on ratio 4 and three-level Neumann), so the formula was changed.

### Formula change

Composite only (`CompositeSolve.cpp`): replace finest-only `||A||` by an effective norm with the same volume weights as `||H||_w`:

    anorm_eff = sqrt( sum_l A_l^2 * phi2_l  /  phi2 ),   A_l = 4 sum_d 1/dx_{l,d}^2
    =>  anorm_eff * ||H||_w = sqrt( sum_l A_l^2 * ||H||_{w,l}^2 )

Implemented as one `exact_sum_multi` over the existing H² MultiFabs with weight `vol[l]*A_l^2` (no extra field temps). Single-level path unchanged.

### After (same cases, level-wise ||A||)

| bc | nlev | ratio | backend | check | floor | limit | floor/check | floor reduction |
|---|---|---|---|---|---|---|---|---|
| neumann | 2 | 2 | hypre | 8.42e-14 | 2.00e-13 | 1e-12 | 2.4 | 3.1× |
| neumann | 2 | 2 | mlmg | 1.68e-13 | 2.00e-13 | 1e-12 | 1.2 | 3.1× |
| neumann | 2 | 4 | hypre | 2.76e-13 | 5.45e-13 | 1e-12 | 2.0 | 4.5× |
| neumann | 2 | 4 | mlmg | 3.48e-14 | 5.45e-13 | 1e-12 | 16 | 4.5× |
| neumann | 3 | 2 | hypre | 1.31e-13 | 2.52e-13 | 1e-12 | 1.9 | 9.8× |
| neumann | 3 | 2 | mlmg | 1.69e-13 | 2.52e-13 | 1e-12 | 1.5 | 9.8× |
| neumann | 3 | 4 | hypre | 1.57e-12 | 2.59e-12 | 2.59e-12 | 1.7 | 15× |
| neumann | 3 | 4 | mlmg | 6.36e-14 | 2.59e-12 | 2.59e-12 | 41 | 15× |
| periodic | 2 | 2 | hypre | 3.87e-14 | 6.54e-14 | 1e-12 | 1.7 | 2.4× |
| periodic | 2 | 2 | mlmg | 8.89e-14 | 6.54e-14 | 1e-12 | 0.74 | 2.4× |
| periodic | 2 | 4 | hypre | 3.27e-13 | 2.20e-13 | 1e-12 | 0.67 | 2.8× |
| periodic | 2 | 4 | mlmg | 3.12e-14 | 2.20e-13 | 1e-12 | 7.0 | 2.8× |
| periodic | 3 | 2 | hypre | 1.63e-13 | 9.59e-14 | 1e-12 | 0.59 | 6.5× |
| periodic | 3 | 2 | mlmg | 1.11e-13 | 9.59e-14 | 1e-12 | 0.86 | 6.5× |
| periodic | 3 | 4 | hypre | 1.00e-12 | 1.18e-12 | 1.18e-12 | 1.2 | 8.4× |
| periodic | 3 | 4 | mlmg | 3.69e-14 | 1.18e-12 | 1.18e-12 | 32 | 8.4× |
| neumann | 2 | 2 (n=16) | hypre | 4.41e-14 | 5.02e-14 | 1e-12 | 1.1 | 3.1× |
| neumann | 2 | 2 (n=16) | mlmg | 2.69e-13 | 5.02e-14 | 1e-12 | 0.19 | 3.1× |

Floor/check is now typically O(1)–few (in line with single-level Neumann and the factor `c=10`), except where MLMG's true residual is far below the HYPRE one (ratio-4 three-level / MLMG residual ~1e-14). The limit is no longer inflated by finest-only `||A||` on the ratio-2 and most ratio-4 cases.

Note: MLMG reports `NotConverged` on neumann 3-level ratio 4 (pre-existing, true residual still ~6e-14; same as docs/pressure/07 §14.2); residual_ok stays true.

## D-087 (A-72 E3 CHANGE, docs/pressure/08 §8.1)

Architect ACCEPT WITH CHANGE: keep the level-wise form but restore margin. Composite only:

    anorm = min(2 * anorm_eff, A_finest)

where `anorm_eff = sqrt(sum_l A_l^2 phi2_l / phi2)` is the previous E3 value and `A_finest` is the last (finest) `A_l` of that loop. Same reduction (weight `vol[l]*A_l^2`); no new reduction. Do **not** use `max(old, new)` or bare factor-10.

Shared helper `operator_anorm(dx, dom)` in `CommonLayer.H` serves the single-level path and this loop (D-084).

Expected where the floor binds (n=32, 2 ranks): Neumann 3-level ratio 4 floor/check ~3.3 (floor ~5.2e-12); periodic ~2.4. Tests: `pb_resid_e3_worst` (HYPRE Neumann and periodic, 3-level, ratio 4, n=32) and `pb_resid_e3_near` (near-limit negative control).

