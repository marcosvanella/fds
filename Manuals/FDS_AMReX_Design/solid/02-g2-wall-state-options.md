# 02 — Wall state on refined levels: option A vs option B (input to ADR-003 Spike G2)

Owner: AMR Solid Phase Lead · Status: **for ruling by the AMR Chief Architect**, 2026-09-26 · Reference tree: FireX `36975d765f`.
Context: FR-041b (cross-level wall-state transfer, after Phase 6), ADR-003 Option B sub-choice, D-010/FR-041a (Phase 5 freeze), D-038 (`VARIABLE_THICKNESS` frozen at setup), R-59, R-60. Detail in `01-solid-phase-amr-spec.md` §2-§5.

## The two options
- **A, finest-ever records.** Each 1-D record (`BOUNDARY_ONE_D` plus the solid parts of `BOUNDARY_PROP1/2`) is keyed to the OBST face patch at the finest resolution that has ever covered it and stored outside the box layout. A coarse owner drives its r² records with its own gas-side inputs and applies the area-sum of their mass, species and heat fluxes; its surface temperature for the gas BC is the area mean of the records (σT⁴-weighted for emission). First-time refinement splits a coarse record into r² identical copies. Records never change on regrid, only their owner.
- **B, remap.** Records live at the owning level's resolution. Coarsening area-averages r² profiles into one; refinement copies one into r².

## Tradeoffs
| | A | B |
|---|---|---|
| Ownership at level changes | Owner changes, records do not. Coarse-to-fine: split once (exact per area); fine-to-coarse: aggregate fluxes each step. | Every change rewrites records. Fine-to-coarse needs a conservative depth remap: profiles have different node counts and layer thicknesses after renoding (wall.f90:2361-2804), per-material densities, and T recovered from enthalpy through a nonlinear `RHO_C_S(T)`. Event state (`T_IGN`, `BURNAWAY`, layer counts) and accumulators (`PART_MASS`, `A_LP_MPUA`, R-60) need one rule each. |
| `VARIABLE_THICKNESS` (D-038) | Consistent by construction: each record keeps the thickness frozen at setup. A split copies the coarse thickness; a coarse owner keeps its fine records' own thicknesses. | Averaging records with different frozen thicknesses gives a thickness no level ever had; the remap must also move the back boundary. Breaks the "frozen at setup" reading of D-038 unless thickness is excluded from averaging, which then breaks mass conservation. |
| Conservation | A regrid changes no record, so solid mass and enthalpy are unchanged by it exactly, not just to round-off. Burn-away heterogeneity is kept through coarse periods. | Mass and enthalpy conservable to round-off with care. Re-refinement is lossy: r² copies of the mean lose spatial detail (burn front, partial burn-away). Surface temperature jumps at regrid. |
| Parity | A face owned by one level from t=0 has one record, identical to FDS on that grid (FR-047). | Same for never-regridded faces. |
| Cost | Solid work and memory stay at the finest-ever record count (R-59). No remap kernel. Records must move between ranks when box ownership changes (small: `NWP`×(1+materials) values each). | Solid work follows current resolution. Regrid-time remap kernel, irregular and branchy (bad for GPU, NFR-048). |
| Code | New keyed store plus aggregation; solid kernel unchanged. | Remap for every solid field; kernel unchanged. |

## Recommendation
**Option A.** It is exact, keeps D-038 without special cases, leaves the kernel untouched, and avoids a remap that is hard to get right for renoded, multi-material, charring profiles. Its only real downside is cost after coarsening (R-59), which is bounded: extra solid work is at most (r² − 1) × the wall share of step time on faces that were refined and then coarsened.

## Smallest spike that settles it (G2a, replaces the B half of G2)
The choice turns on cost alone, so measure cost before building anything:
1. **Measure the wall share.** Run existing FDS on `Fires/box_burn_away1`, `Energy_Budget/energy_budget_solid` and one charring flame-spread input, and read the `WALL` column of `CHID_cpu.csv` (main.f90:4146; `WALL_BC` time, wall.f90:278). This includes gas-side BC work, so it is an upper bound on the 1-D solve share.
2. **Bound A's overhead.** Overhead ≤ (r² − 1) × wall share × the fraction of wall faces that were refined then coarsened. Evaluate for r = 2 and 4 with the fraction at 1 (worst case).
3. **Decide.** If the worst case at r = 2 is below a threshold the Architect sets (suggested: 25 % of step time), choose A and scope G2 to A only: 2-level static box, heated charring OBST, refine, coarsen and refine across its face; pass signals as in ADR-003 (solid energy unchanged across regrid, surface T continuous) plus the A-specific check that every record is bitwise unchanged across each regrid operation. Only if the bound fails does B get prototyped, with the same case and the remap error as the added measurement.

Step 1 needs no new code and fits in the Phase 2 prototype timing work.


**Outcome.** Option A accepted (`adr/drafts/ruling-FR041b-G2.md`). G2a passed at R=2; see `03-g2a-cost-check.md`.
