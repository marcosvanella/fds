# 03 · G2a cost check for FR-041b option A (finest-ever records)

Owner: AMR Solid Phase Lead · Status: **result, for the AMR Chief Architect** · 2026-09-26 · Per `adr/drafts/ruling-FR041b-G2.md` §3.

## Verdict
**R=2 passes** the 25% step-time gate for the representative case (`couch`): 0.16 on the loose wall-timer bound, 0.06 on the isolated 1-D solve bound. The memory check passes with wide margin. No record-depth cap (§3.5) is needed on this evidence. One caveat: a small, solid-dominated charring case (`box_burn_away_2D_residue`) exceeds the gate (0.40). Cases like it are where a cap would bite, if it is ever needed.

## Method
- Unmodified FireX release build (GNU/OpenMPI), unchanged Verification inputs. Timing comes from `CHID_cpu.csv` (header main.f90:4146). The `WALL` column is the whole `WALL_BC` timer (wall.f90:278).
- s_wall = WALL / Total, summed over ranks. This is an upper bound because `WALL` also includes gas-side boundary work (near-wall gas variables, heat transfer coefficient, `CALCULATE_ZZ_F`), which does not scale with record count.
- s_solve isolates the 1-D solve. The 1-D solve runs on the corrector every `WALL_INCREMENT` steps (default 2, cons.f90:588; counter at main.f90:1007-1011, gate at wall.f90:92). Each case was rerun with `WALL_INCREMENT=1`. With per-step times w, s_solve = (w_inc1 − w_inc2) / t_inc2, which is the solve's share of the default run.
- Overhead bound = (R²−1) · s, with f=1 (every coarse face refined to the full depth), compounded across levels: a factor of 3 at total ratio R=2 and 15 at R=4. Fine records advance with the owner's dt and increment (spec §4), so each record costs the same as a coarse one.

## Time results

| Case | Ranks | s_wall (max rank) | Bound R=2 / R=4, wall timer | s_solve | Bound R=2 / R=4, solve |
|---|---|---|---|---|---|
| `couch` (flame spread, 8 meshes, 600 s) | 8 | 0.052 (0.081) | **0.16** / 0.79 | 0.021 | **0.06** / 0.32 |
| `box_burn_away1` | 4 | 0.102 (0.138) | 0.31 / 1.53 | ≈0 (noise) | ≈0 |
| `energy_budget_solid` | 1 | 0.103 | 0.31 / 1.54 | 0.012 | 0.04 / 0.18 |
| `box_burn_away_2D_residue` (charring) | 1 | 0.234 | 0.70 / 3.51 | 0.133 | 0.40 / 1.99 |

Only R=2 gates. The worst rank of `couch` on the loose bound is 3 × 0.081 = 0.24, still within 25%.

Caveats:
- The small cases have a high surface-to-volume ratio and short runs (the 2D case totals 4 s of CPU time), so their shares carry large relative noise. `box_burn_away1` shows the wall timer is almost entirely gas-side work.
- The `couch` pair ran on 8 ranks (one per mesh), before the 4-core limit for this shared machine was set. The other cases used 4 ranks or fewer.
- The machine was shared. The last part of the `couch` `WALL_INCREMENT=1` run overlapped other jobs. That inflates its time and biases s_solve upward, which is the conservative direction.
- f=1 is the worst case. Under A2, records exist only inside the refinable region, so the real overhead scales with the refined fraction of faces.

## Memory bound
There are no per-record node counts in FDS output, so the bound is parametric. Fields come from type.f90:176-179 and 217-376, with allocation sizes at func.f90:4791-4820. Here N = `N_CELLS_MAX`, M = materials, L = layers, S = tracked species, P = particle classes.
- Reals per record ≈ (6+M)N + 7 + 2M + ML + M + 11L + 1 + 9S + 32 + 5P + 12. Integers ≈ 18 + 4L + M. Bytes = 8·reals + 4·ints, in a flat layout (NFR-048).
- **Assumed** examples: N=100, M=3, L=3, S=6 gives 8.6 KB per record, so 16 GB holds 2.0 M records. N=30, M=2, L=2, S=3 gives 2.9 KB, so 16 GB holds 5.9 M records.
- The current Fortran derived-type layout adds about (43+2M) allocatable descriptors of about 64 B each, roughly 3 KB per record. This supports the flat SoA layout in NFR-048.
- `couch` (geometry estimate): about 9,850 records (about 9,040 gypsum boundary faces after the open vent and obstruction footprints, plus about 810 obstruction faces). Finest-ever everywhere: 39 k records at R=2 (0.1–0.3 GB) and 158 k at R=4 (0.4–1.3 GB). **Passes**.
- General rule: level-0 record count × (R²)^levels × bytes per record ≤ 16 GB. At R=4 that allows about 125 k (heavy) to 370 k (light) level-0 faces per node or GPU, if every face is refined to full depth.

## Pass signals carried to G2 (A only)
Solid mass and enthalpy change across a regrid ≤ 1e-12 relative. Every record bitwise unchanged across each regrid after its first split (A3).

Note: the bulk plot and restart output of the runs in `vv-runs/g2a` was deleted to free disk. The CSV timing files, logs and inputs used by `share.py` are kept, and the runs can be regenerated from the inputs.
