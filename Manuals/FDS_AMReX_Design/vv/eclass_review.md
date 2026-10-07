# V&V review: E_class package (D-077 / S1 host forced-flip survey) — re-verify

Reviewer: AMR V&V Lead (executor). Subject: Solid Phase Lead revised package in
`src-s5sp` / branch `s5-wall-sp` commit `73f1bdd6db`
(`amrex/s4_mass/s5_gen/README_solve_capture.md` §8c — predicate form).
Solid Phase tree was not edited. This note is docs-only.

Prior review (NEED RE-RUN) was against tip `92136ff50a` / docs commit `6169b3c`.

## Verdict

**INTERIM CO-SIGN** of the P1-densified predicate `E_class` table in §8c at
`73f1bdd6db`, **with the limits stated below**.

This is **not** a final co-sign of `E_class` for the device-tier flip-budget
gate. Final co-sign waits until full P2 (`N≥50` per (class, sign) on every case
of S) is folded into S and the predicate table is regenerated.

### What is signed (interim)

- Predicate form: `E_class(c, q) = 2 × max |Delta_cf(q)|` over survey S
  (`test/data/flip_survey_*.json`), quantities
  {node TMP, TMP_F, Q_CON_F, Q_RAD_IN, Q_RAD_OUT, layer mass}.
- Empirical Level B numbers for D1 / D2 / D5a / D5b as printed by
  `solve_flip_survey.py table` at this tip (byte-match README §8c).
- Derived Level B labeling: D3 = M_DOT / Y_O2 bound (survey supporting only);
  D4 = 1e-4 K extraction tolerance (no S1 flip hook); delamination / burn-away
  = rate-only / no Level B (D-077 A3).
- Q_RAD_IN exercised on derived `back_wall_test_short` (T_END=5 s): max
  |dQ_RAD_IN| = 8.2×10⁻⁹ W/m² on the flipped wall; cross-wall leak ≤1.5×10⁻³
  reported and **excluded** from Delta_cf (radiation coupling).
- Caveats retained: sample ≠ census for D1/D5a; SOLID_PHASE_ONLY on all but
  back_wall; under-estimation possible outside surveyed calls; E_class(Q_RAD_IN)
  rests on one derived case and three flips.

### What is not signed

- Final device-tier `E_class` until full P2 lands (Solid’s own ask; V&V agrees).
- Any claim that P2 slices already move (or bound) the maxima — delivered P2
  seed-7 slices do not move maxima and are incomplete.
- Gas-coupled amplification beyond the documented cross-wall leak on back_wall.

## Verification performed

| Check | Result |
|---|---|
| Commit `73f1bdd6db` on `s5-wall-sp` | Present; HEAD = `73f1bdd6dbd4194a45d7df45a96d45db0ed4ac79`; message matches (“predicate form, densification, README §8c”) |
| JSON inventory | 17 files `test/data/flip_survey_*.json` (incl. densification, census, back_wall t5/r2/p2, a_fo p2/r2/d1n500/d5a*, kc d5an100) |
| Regenerated table | `/usr/bin/python3 solve_flip_survey.py table test/data/flip_survey_*.json` — three output tables **byte-identical** to README §8c predicate / class-maxima / per-case rows |
| Fired counts | D1 **1471** (≥1000); D2 **900** (full census); D5a **1011** (≥1000); D5b **1462** (full census); D3 **40** supporting only |
| Applicable pools | D1 66884; D2 900; D5a 92200; D5b 1462; D3 509 |
| Predicate E_class (spot) | D1: node TMP 17.5, TMP_F 12.7, Q_CON_F 67.8, Q_RAD_IN 1.64e-08, Q_RAD_OUT 1.51e+03, mass 1.71e-05 (shrink_swell +1). Matches 2× class maxima. |
| Q_RAD_IN case | `flip_survey_back_wall_test_short_t5.json`: `max_q_rad_in` ≈ 1.20×10⁴ W/m² (quantity exercised); three D1 +1 flips at seq 3601 / 10801 / 18000; max \|dQ_RAD_IN\| = **8.20×10⁻⁹**; leaks 1.45×10⁻³ / 9.56×10⁻⁴ / 0 — matches §8c; leak not in Delta_cf |
| P2 schedule | §8c: slices seed 7 delivered (a_fo D1/D5a; back_wall D1); **explicit written schedule** — P2 N≥50/(class,sign) on every case of S before first device-tier gate, then fold and regenerate. States V&V will not final-sign on P1-only. |
| Level B labels | D1/D2/D5 empirical; D3/D4 derived; delam/burn-away rate-only — present in predicate table and caveats |
| SOLID_PHASE_ONLY + leak caveat | Stated; radiation leak excluded from Delta_cf |
| Under-estimation notes | Kept (sample ≠ census for D1/D5a; Q_RAD_IN thin; D5a flat-table zeros) |
| `test/test_solve_flip.py --no-live` | **68/68 PASS** (offline suite green; live binary not required for this check) |

## Prior gaps (from tip `92136ff50a`) — disposition

| # | Prior gap | Status |
|---|---|---|
| 1 | No predicate-form E_class (= 2 × max \|Δ_cf\|) | **CLOSED** |
| 2 | Sample density <1000 where pool allows (D1, D5a) | **CLOSED** (D1 1471, D5a 1011; D2/D5b census) |
| 3 | P2 missing — no final co-sign on P1-only | **OPEN for final**; **satisfied for interim** via delivered slices + written schedule before device gate |
| 4 | Q_RAD_IN never exercised | **CLOSED** (back_wall_test_short T_END=5) |
| 5 | D4 derived OK if stated | **CLOSED** (stated) |
| 6 | D3 derived-input OK as supporting; Level B lean derived | **CLOSED** (stated) |
| 7 | SOLID_PHASE_ONLY must be stated | **CLOSED** (+ radiation-leak caveat) |

## Remaining open before FINAL co-sign

1. **Full P2:** N≥50 per (class, sign) on every case of S, fold into S, regenerate
   predicate table; V&V re-verifies numbers and offline tests, then final co-sign
   (or delta list).
2. Optional hardening (not blocking interim): thicken Q_RAD_IN evidence beyond
   three D1 flips on one derived short case; keep under-estimation caveat until
   P2 is in.

## Exact ask remaining for Solid Phase Lead

Complete the scheduled P2 before the first device-tier flip-budget gate run:
`run CASE --p2 --seed S` with N≥50 per (class, sign) on every case of S, merge
into the survey JSON set, regenerate §8c, and ping V&V for final re-verify.
Until then, treat this interim co-sign as authorizing use of the P1 predicate
numbers as the **working** Level B bound for host/planning work only — **not**
as a PASS for the device-tier gate.

## Interim status of this package

- **Accepted as (interim):** predicate-form Level B `E_class` for D1/D2/D5a/D5b
  (empirical) and labeled derived bounds for D3/D4 / rate-only delam/burn-away,
  at Solid tip `73f1bdd6db`, with Q_RAD_IN exercised and densification meeting
  the ≥1000 / census bar.
- **Not accepted as:** final signed `E_class` for the device-tier flip-budget
  gate under D-077 (1) V3 / solid/09 §7 until full P2 is folded in.
