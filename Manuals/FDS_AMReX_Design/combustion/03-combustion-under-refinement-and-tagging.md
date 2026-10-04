# 03 — Combustion under refinement (D-050, D-058) and species tagging advice for R4

Owner: AMR Species & Combustion Lead · Status: DRAFT for the Architect and Role 3. Inputs: `combustion/01` (spec), `combustion/02` (sign-offs), D-031, D-050, D-058, D-062, FR-011, FR-012, FR-025, `role3-regrid-transport-plan.md` (R3, R4). Source citations are FireX `36975d765f`.

## 1. Species clip accounting

Three different things can change a species amount without a flux. They must not be mixed in one counter, because they have different owners and different budget effects.

| Source | Where | What moves | Budget effect | Counter |
|---|---|---|---|---|
| In-level clip | `CHECK_MASS_DENSITY` (`mass.f90:799-939`), D-031 gather form | A cell outside its allowed range is set to the limit; the excess is taken from, or added to, the six neighbours that are interior targets | Zero inside a level only if every target is a valid, uncovered cell of the same level. A target that is a covered coarse cell is overwritten by `average_down`; a target outside the interior is skipped. Both are a composite budget change | `CLIP_RHO_ZZ(N)` max-reduction (existing), plus new integer counts per level: cells clipped, cells whose scatter target was covered or outside |
| Regrid clip | Limited conservative interpolation in R4 (D-058) | Children are rescaled so their sum equals the parent | Zero by construction (children sum equals parent) | One count per parent, per level, per regrid (FR-025) |
| Level-operation clip | After `FillPatch`, `average_down`, interface flux overwrite | Any residual realizability fix (a species outside [0,1] or sum not 1) | Nonzero, and must be reported | Count and level, warning when nonzero (FR-025) |

Rules proposed:
1. **Flags come from uncovered valid cells only.** A covered coarse cell never starts a clip, because its value is replaced at the next `average_down`. This is the D-031 condition and it keeps the flags independent of the box layout.
2. **The scatter does not cross a level boundary.** The gather writes targets inside the same level. A target that is a covered coarse cell is skipped and the skipped amount is counted in the "outside" counter, so the loss is visible and not hidden in the budget residual. Same for a fine-box target outside the valid region at a coarse-fine interface.
3. **Per-step budget line.** The composite budget (FR-020/021) gets one extra line per species and level: "mass removed or added by clips", the sum over clipped cells of the clip delta, taken with the D-028 exact sum. On a positive-field case (no clipping expected) the line must be exactly 0 and the count 0. On a case that clips, budget error minus the clip line must be round-off. This is the test that tells a clip loss from a flux error.
4. **Regrid.** Per D-062 the transfer works on rho and rho*Z, so the clip on a child is on rho*Z. Rescaling the children to the parent sum keeps the 1e-12 per-regrid limit of FR-012 only if the rescale is applied after the limiting and before `average_down` is repeated. Count it, and test that the count is 0 on the smooth blob.
5. **Open question, unchanged** (spec §6 item 1): keep the redistribution level-local, or replace it by a flux limiter so that clipping is only a diagnostic. I still recommend level-local for Phase 3 and the diagnostic form later, because FDS results at T2 do not depend on it for the anchor cases and the gather already exists.

## 2. The `Z_TEMP` pad of patches 0001 and 0002 under refinement

- The pad matters at **physical-wall faces only** (`mass.f90:147-179` and the `divg.f90` twins): the one-face calls fill three of the four `Z_TEMP` elements, and MP5 reads the fourth (`func.f90:1443-1448`; `Z_TEMP(3)` for A>0, `Z_TEMP(0)` for A<0, and the fifth stencil point is extrapolated from it). The patches pad with 0 so the result is defined.
- **Do not reuse the wall path at coarse-fine interfaces.** An interface face of a fine box must use real ghost data from the `FillPatch` (coarse-to-fine interpolation), not a padded wall call. With the 4-element stencil `(I-1, I, I+1, I+2)` the box needs **2 ghost layers** for scalars where a limiter other than the default is used, and the same two layers are needed for the default limiter on the high side. Role 3 should confirm that the scalar `FillPatch` uses `ngrow ≥ 2`, and fail at input time if the chosen limiter needs more than the allocated ghost width.
- Under the D-050 flux overwrite the coarse face flux is replaced by the area-sum of the fine face fluxes, so a pad error at a wall does not propagate across levels. It stays a local, bounded error at the wall, the same as in single-mesh FDS.
- Until patches 0001 and 0002 are applied upstream, bitwise tests that use MP5 with wall faces compare against undefined memory. Tests for the other limiters (`SUPERBEE`, `CHARM`, `GODUNOV`, `MINMOD`) are not affected. MP5 stays out of the Phase 3 gate set.

## 3. Species tagging advice for R4 (D-058, FR-011)

The working rule is undivided difference with `TAG_KEEP` hysteresis, threshold about 0.05, buffer `N_ERROR_BUF = 2`. For species and combustion I advise the following.

1. **Tag on mass fraction Y (and mixture-fraction-like species), not on rho*Z.** Y is bounded in [0,1], so one absolute threshold has a meaning. rho*Z follows density and would tag the hot plume core everywhere. Derive Y from the mass-weighted state (D-062), then take the undivided difference `max over the six face neighbours of |Y_nb - Y|` per tagged species.
2. **Per-species scale for trace species.** For major species (fuel, O2, CO2, H2O) the 0.05 default is a sensible start. For trace species (CO, soot, smoke, aerosols) 0.05 never triggers. Allow a per-species scale `Y_REF(N)` (default 1, user-settable, no hidden default for soot) and tag on `|dY|/Y_REF`. Do not tag on the background species, because it carries no signal in a mixture with several species.
3. **Combustion zone: tag on HRRPUV, not only on species.** The flame sheet is where fuel and oxidizer both exist, and a single species difference can miss it. Use `HRRPUV > HRRPUV_MIN` (absolute, user value; a relative form needs a global maximum each regrid, which is a reduction and must use the D-028 exact sum). A cell with reaction is tagged for refinement whenever the flame, as a whole, is inside the refined region. This matters because FDS's combustion model depends on the cell size (the mixing time scale in `fire.f90` uses the local `DX` and the `LES` filter width), so a flame that is only partly refined mixes at two different resolutions.
4. **Keep the hysteresis stricter on species than on temperature.** Use `TAG_KEEP ≈ 0.5 × threshold` and a regrid interval of at least the buffer-crossing time (`N_ERROR_BUF ≥ ceil(R · CFL per axis) + 1`, already in `test-plan.md`). A clipped species (section 1) can produce a one-cell spike, so do not tag on a cell that was clipped in the same step; use the value after the realizability pass.
5. **Optional feedback tag.** If a level-operation clip count (section 1, third row) is nonzero in a cell neighbourhood, tag it for refinement at the next regrid. Off by default; gives an honest signal that the resolution is too low there.
6. **Unit cases for R4** (cell by cell against an independent evaluation, as the plan requires): (a) a smooth blob of one species, positive everywhere; (b) a step in a trace species with and without `Y_REF`; (c) a flame-sheet case where fuel and oxidizer overlap in one cell and the species difference alone is below the threshold, so that only the HRRPUV rule tags it; (d) the clipped-cell case, where the tag is taken after the realizability pass.

## 4. What I need from others
- Architect: accept or change rule 2 of section 1 (skipped scatter targets are counted, not redistributed) and the extra budget line.
- Role 3: confirm scalar ghost width `≥ 2`, and the new counters in the regrid report.
- V&V Lead: add the clip-line check to the budget test, and the two combustion tagging cases to the R4 unit set.
