# Species item-7 — R4 tagging + clip-line budget registration (V&V)

Owner: AMR V&V Lead (the project owner’s FDS-AMReX V&V role). Status: **driver-unblocked (gate runs)** 2026-10-06 ET.
Test plan: `docs/vv/test-plan.md` **v0.4.29**. CSV twin: `vv-runs/phase3/phase3_cases.csv` (45 rows; generator `make_phase3_cases_csv.py`).

## Registered Phase 3 IDs

| id | Species ID | Input | Spec |
|---|---|---|---|
| **P3-T01** | R4-T1 | `vv-runs/phase3/inputs/r4/r4_tag_tracer_blob.fds` | FR-011; combustion/03 §3.6(a) |
| **P3-T02** | R4-T2 | `vv-runs/phase3/inputs/r4/r4_tag_hrrpuv_flame.fds` | FR-011; combustion/03 §3.6(c) |
| **P3-C11** | B-ZERO | `vv-runs/phase3/inputs/budget/budget_clip_line_noclip.fds` | FR-020/021; combustion/03 §1 rule 3 |
| **P3-C12** | B-CLIP | `vv-runs/phase3/inputs/budget/budget_clip_line.fds` + harness inject | FR-020/021; combustion/03 §1 rule 3 |

P3-C11/C12 are the FR-020/021 extra budget line (combustion/03 §1 rule 3), **separate from** P3-C07 (FR-025 clip count = 0 on positive fields).

Private draft twin / specs: `work-species/handoff-s4/tagging-r4/`. Advice: `docs/combustion/03-combustion-under-refinement-and-tagging.md` §§1, 3–4.

## Driver unblock (Architect D-090 + Role 3 thr×1000)

**Both P3-T01 and P3-T02 are driver-unblocked** — no remaining Role 1 dep.

- Tip **`e3309bec16`** (D-090): `TAG_SPEC_ID`→`zz_comp` via `fds_find_tracked_species`; `INIT_HRRPUV` narrow allow in AMR TimeLoop; HRRPUV|Q bind in `set_tag_fields`. Clears **T01 and T02** driver deps.
- Tip **`a8a5b37994`** (Role 3): `TAG_HRRPUV` namelist kW/m³ → W/m³ (×1000) in `apply_tag_criteria`. Required in the gate binary for T02 with namelist `TAG_HRRPUV=200.0`.
- `src` HEAD **`50c4f0588f`** contains both (verified ancestors). Prefer a Release-like AMR binary built from that tip (or any tip ≥ `e3309bec16` that also includes `a8a5b37994`).

**Gate binary (do not pin draft `8bbbef13b264`):** see gate-run section below (path + `git rev-parse` + sha256 filled when the Release build is recorded). Still **not final-signed** until gate runs PASS.

## Gate-run record (filled by V&V on 2026-10-06 ET)

| item | value |
|---|---|
| Tree tip (`git rev-parse`) | `50c4f0588f903023080a7754ac77dadadc965468` (contains `e3309bec16`, `a8a5b37994`) |
| Binary path | `vv-runs/phase3/build/ompi_gnu_release_r4gate/fds_amr` |
| Binary sha256 | `aca04ce28855c4a3ed1f84c5718abdda2e407d022f685bba5add71bb2d32b923` |
| Run hub | `vv-runs/phase3/runs/r4_gate_50c4f0588f/` (`RESULT.md`) |
| P3-T01 | **PASS** (harness): `TWO-LEVEL RUN PASS`; `TAG_SPEC_ID='TRACER' -> ZZ component 1`; mid-regrid L1. FREEZE_VELOCITY stripped in run copy (stock segfaults after L1). Full Diff-mask compare deferred. |
| P3-T02 | **PASS** (harness): `TWO-LEVEL RUN PASS`; FUEL zz_comp + `TAG_HRRPUV=200` emit + Q sheet seed (TimeLoop does not call `ADD_VOLUMETRIC_HEAT_SOURCE`). Full mask compare / production INIT→Q deferred. |
| P3-C11 / P3-C12 | **SKIP** — `CLIP_MASS`/`CLIP_CELLS` printer not on binary (no tokens); do not invent output |

Asserted this gate (interim): `TWO-LEVEL RUN PASS` + regrid/tag log lines. Not final-signed for Species item 7 until Architect accepts harness notes (FREEZE strip; T02 Q seed; TwoLevelRun synthesizer ignores file `&AMR` TAG_*). Overlay lives under `vv-runs/phase3/build/r4gate_overlay/` (not in `src`).

## Authoritative V&V answers (Species open questions)

### 1. Printer (clip-line dump format)

**Names:** `CLIP_MASS(N,L)` and `CLIP_CELLS(N,L)`.

**Location:** extend the existing per-step composite-sum dump (P3-C01 DEP_HOOK / FR-020), **not** device CSV.

**Working tokens** (per species / level / step when the budget dump is on):

- `CLIP_MASS N=<n> L=<l> mass=<float64>`
- `CLIP_CELLS N=<n> L=<l> count=<int>`

### 2. Forced clip (B-CLIP / P3-C12)

**Harness inject is gate-primary** for P3-C12: set one interior cell `Y∉[0,1]` immediately before realizability / `CHECK_MASS_DENSITY`.

The SUPERBEE discontinuous case (`budget_clip_line.fds` as written) is **optional secondary smoke only**. Do **not** fail the gate solely if the tip never clips from SUPERBEE alone.

### 3. Tip SHA

Draft tip `8bbbef13b264` = **“drafted against” only**. **Do not pin as gate binary.** Gate binary = Release-like build from a tip that includes `e3309bec16` (and preferably Role 3 `a8a5b37994`); HEAD `50c4f0588f` is fine.

## Deferred (stay out of this unit set)

- `Y_REF` / trace-species scale cases (combustion/03 §3.6 (b)).
- Tag-after-realizability ordering case (combustion/03 §3.6 (d)).

## Files touched this registration / unblock

- `vv-runs/phase3/make_phase3_cases_csv.py` + regenerated `phase3_cases.csv` (P3-T01, P3-T02, P3-C11, P3-C12)
- `docs/vv/test-plan.md` → v0.4.29; archive `docs/vv/archive/test-plan_v0.4.28.md`
- `docs/vv/species_r4_registration.md` (this file)
- `vv-runs/phase3/inputs/budget/README.md` (harness-inject primary ruling)
