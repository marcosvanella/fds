# V&V review: draft E_class table (D-077 / S1 host forced-flip survey)

Reviewer: AMR V&V Lead (executor). Subject: Solid Phase Lead draft in
`src-s5sp` / branch `s5-wall-sp` commit `92136ff50a`
(`amrex/s4_mass/s5_gen/README_solve_capture.md` §8c).
Solid Phase tree was not edited. This note is docs-only.

## Verdict

**NEED RE-RUN** before co-signing `E_class` as the Level B bound under D-077 (1) V3.

The draft table numbers match the committed JSON and regenerate from
`solve_flip_survey.py table`. The under-estimation caveat is already in the
README. That is enough to treat the draft as verified *survey evidence*, not
enough to sign it as `E_class`.

## Verification performed

| Check | Result |
|---|---|
| Commit `92136ff50a` on `s5-wall-sp` | Present; HEAD message matches (5 cases, 492 flips, DRAFT E_class) |
| Files present | `README_solve_capture.md`, `solve_flip_survey.py`, `test/test_solve_flip.py`, five `test/data/flip_survey_*.json` |
| Regenerated table from JSON | Byte-identical to README §8c summary and per-case rows |
| Totals | 492 flipped records; 129,555 applicable (call,sign) pairs (~129.5k) |
| Class maxima | D1 7.39 / 5.1 / 33.9 / 605 (shrink_swell +1); D2 9.22 / 5.93 / 748 (shrink_swell −1); D3 node 0.037 K (TMP_F 0.0576); D5a 2.71e-4 K; D5b 8.15e-4 K |
| Q_RAD_IN | Every effect field 0 across all 1,476 imm/max/end samples |
| Cross-wall leak | 0 on every flip |
| `test_solve_flip.py --no-live` | **63/63 PASS**. Claimed 68/68 includes the 5 live tests that need `vv-runs/solve-capture/fds_flip` (absent here; SKIP). Offline suite is green. |

Draft-summary rounding: "2.7e-4" / "8e-4" match 0.000271 / 0.000815. D3 "max dTMP 0.037 K" is the node max; TMP_F max is 0.0576 K (in the table, omitted from the short summary).

## What D-077 / solid/09 requires for a signed E_class

From `solid/09-flip-budget-gate.md` §7 and D-077 (1) V3:

- Forced-flip survey on the host over **P1 and P2**.
- Prefer calls with **near-threshold** decisions; **at least 1000 calls per class where reachable**.
- `E_class` = **twice the largest observed** `Delta_cf`, signed by Solid Phase Lead and V&V Lead.
- Until signed, flips of that class are recorded and the class has **no PASS**.
- D3 / D4 / layer-removal / renoding-check use **derived** Level B bounds where the code supplies them; empirical `E_class` is for D1, remesh-ratio / cell-count (D2), and D5.

V1–V4 (Clopper-Pearson gate, exposed distinct-input denominator, Level A/B, P2 carries gate if P1 short with exposure-tag coverage) are already ratified; this review is only about signing `E_class`.

## Gaps blocking co-sign

1. **No predicate-form `E_class`.** The draft reports maxima. It does not propose
   `E_class = 2 × max_observed` per (class, quantity), nor name which quantities
   enter Level B (node TMP, TMP_F, Q_CON_F, Q_RAD_IN, Q_RAD_OUT, layer mass).

2. **Sample density below the §7 bar.** Fired flips: D1 112, D2 70, D3 40,
   D5a 190, D5b 80. No class reaches 1000 where the pool allows it (D1 48,884
   applicable; D5a 77,800). Evenly spaced picks over-weight quiet late-time
   conduction calls; the README already says maxima may be under-estimated.
   (Complete censuses are fine where the pool is tiny: D2 −1 = 30/30 tagged
   calls; D1 −1 on both shrink_swell cases = 1/1 each.)

3. **P2 missing.** Survey is five host cases only (three pinned, two derived,
   one survey-short of kc). No random corpus. solid/09 requires P1 and P2
   before `E_class` is signed for the device gate.

4. **Q_RAD_IN never exercised.** All cases are `SOLID_PHASE_ONLY`; Q_RAD_IN
   effect is identically 0. No Level B statement on Q_RAD_IN is possible.
   Need at least one case with radiation in (e.g. `back_wall_test`, or another
   pinned wall case with non-zero `B1%Q_RAD_IN`).

5. **D4 not in the empirical table (acceptable if stated).** D4 has a derived
   bound (1e-4 K) and no flip hooks in S1. That is fine under §7 *if* the signed
   package says D4 Level B = derived, not empirical. Delamination / burn-away
   remain rate-only (D-077 A3).

6. **D3 only on a derived, non-pinned input** (`shrink_swell_o2`). Acceptable as
   supplemental empirical evidence; Level B for D3 should lean on the derived
   M_DOT / Y_O2 bound, with the survey as supporting.

7. **SOLID_PHASE_ONLY / no gas feedback.** Survey measures solid-propagated
   effect on one wall. That matches a host Delta_cf survey of the wall solve,
   but must be stated on the signed bound (no claim about gas-coupled
   amplification).

## Exact ask for the Solid Phase Lead

Please deliver a revised package (same worktree / follow-up commit on
`s5-wall-sp` is fine; V&V will not edit that tree) that:

1. **States `E_class` in predicate form** for each empirical class (D1, D2
   remesh/cell-count, D5a, D5b): for each tracked quantity,
   `E_class(qty) = 2 × max |Delta_cf(qty)|` over the survey that will be signed,
   with the under-estimation caveat retained until sample density meets item 3.
2. **Adds ≥1 case with non-zero Q_RAD_IN** and reports dQ_RAD_IN (and leak) for
   it; if still zero after a radiation case, say so with evidence.
3. **Raises sample density** toward ≥1000 fired flips per class where the
   applicable pool allows (priority: D1, D5a). Prefer near-threshold picks
   (TMP_RATIO near k+0.5, REMESH_RATIO near 1, TMP near integer for D5) over
   uniform late-time spacing. Keep full census where the pool is smaller than
   1000.
4. **Adds a P2 (random) forced-flip slice**, or a written schedule that P2 is
   completed and folded into `E_class` before the first device-tier gate run
   (V&V will not co-sign a final `E_class` on P1-only).
5. **Labels Level B sources:** D3/D4 derived; D1/D2/D5 empirical; delamination
   and burn-away rate-only / no Level B.
6. **Keeps** the existing acknowledgements: sample ≠ census, Q_RAD_IN gap until
   item 2, SOLID_PHASE_ONLY, D3 derived-input, flat-table D5a zeros.

After that package lands, V&V will re-verify numbers from JSON, re-run
`test_solve_flip.py` (live if `fds_flip` is present), and co-sign or return
with a short delta list.

## Interim status of this draft

- **Accepted as:** verified S1 host survey evidence at `92136ff50a` (492 flips,
  maxima and caveats as in README §8c).
- **Not accepted as:** signed Level B `E_class` under D-077 (1) V3.
- Device-tier flip-budget gate: classes without a signed `E_class` (or a
  derived Level B bound) have **no PASS** for Level B until the re-run lands.
