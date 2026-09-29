# Ruling: FR-041b wall state on refined levels (ADR-003 Spike G2)

Owner: AMR Chief Architect. Status: **ACCEPTED (Option A), with G2a as a cost check**. Input: `solid/02-g2-wall-state-options.md`, `solid/01-solid-phase-amr-spec.md` §2-§5. Amends ADR-003 (Option B sub-choice) when folded in.

## 1. Decision

FR-041b uses **Option A, finest-ever records**. Each 1-D wall record is keyed to the OBST face patch at the finest resolution that has ever covered it. The record is stored outside the box layout, under the layout-independent face key of FR-046. A coarse owner drives its r^2 records (r^2 per level jump, compounded across levels) with its own gas-side inputs, and it applies their fluxes summed in a fixed order within the record group. The first refinement of a face splits a record into r^2 identical copies.

Option B (remap) is **rejected** on correctness, not cost:
1. It needs a conservative remap over renoded, multi-material profiles (wall.f90:2361-2804), with T recovered through the nonlinear `RHO_C_S(T)` and rules for event state and accumulators (R-60). There is no FDS reference to verify it against.
2. It averages `VARIABLE_THICKNESS` values into a thickness no level ever had, which contradicts D-038.
3. Re-refinement loses burn-front heterogeneity and makes the surface T jump.
4. It adds an irregular regrid-time kernel on the GPU path (D-027).

Cost is the only open question, and G2a bounds it. If G2a fails, the fallback is **not** Option B. It is a cap on record depth (section 3).

## 2. Conditions on Option A

- A1. Records follow D-038. Split copies inherit the thickness frozen at setup; nothing is recomputed at the finer level.
- A2. Records exist only for faces inside the refinable region at the levels it permits. Outside it, a face keeps one record at its level-0 or static-box resolution. The finest-ever depth is bounded by the local maximum level.
- A3. Determinism (FR-005 (i)): the aggregation order within a record group is fixed by the face key. Record migration between ranks at regrid or load balancing copies records bit for bit. Every record is bitwise unchanged across a regrid, except for the first split.
- A4. Outputs (FR-045, FR-070): budgets and boundary files use the owner's aggregated face value. A wall device at a point reads the finest-ever record containing that point.
- A5. Restart writes records with their face keys and finest-ever depth (replacing dump.f90:3981), independent of layout.
- A6. The solid kernel (`SOLID_HEAT_TRANSFER`) is unchanged (FR-047 (a)). Only the driver loop over records changes.

## 3. Spike G2a (cost check, no new code)

1. From stock FireX runs of `box_burn_away1`, `energy_budget_solid` and one charring spread case, read the WALL column of `CHID_cpu.csv` (main.f90:4146) as a share s of step time. Use the existing baseline runs where they exist.
2. Take the worst-case overhead to be (R^2 - 1) × s × f with f = 1. Here R is the total refinement ratio between the coarsest owner and the finest-ever resolution. Report it for R = 2 and R = 4, since FR-010 allows ratios of 2 and 4.
3. Record memory as well: the record count at the finest-ever depth × the bytes per record (`NWP`×(1+materials) values plus scalars), against 16 GB per node and per GPU.
4. **Threshold:** a worst-case overhead of at most 25% of step time at R = 2, and memory within budget, means A stands and G2 is scoped to A only. At R = 4 the result is reported, not gating.
5. **If it fails:** record depth is capped at max_level-1 relative to the owner (a user parameter, default off). Faces under the cap coarsen through the Option A split in reverse only where the r^2 records are bitwise identical, and are otherwise held. Option B is not built. The Chief Architect rules on the cap design if needed.

G2 pass signals (ADR-003, unchanged, plus A3): solid energy and mass conserved to round-off across refine, coarsen and refine; surface T continuous; every record bitwise unchanged across each regrid after the first split.

## 4. Rejected alternatives

- Option B remap (section 1).
- Building both A and B in G2. The cost question is answered by G2a from existing data, and B fails on correctness.
- Freezing ownership forever (extending D-010 past Phase 6). It blocks refinement that follows a moving fire front.
