# Device tolerance for `s5rad_plane3d_ild` scheme 3 (EXP / LIBM-CLASS)

Owner: V&V Lead. Status: **v0.1 PROPOSED pending Architect ratification**. No Architect number exists for this class; this document does not invent one.
Subject: the per-value ulp limit that a device-tier comparison may accept for `s5rad_plane3d_ild` when `RAD_DIFF_SCHEME = 3` (EXPONENTIAL), under `-DS5_LIBM_CLASS`, against the host reference. Schemes 1 and 2 stay bitwise.

Source measurement (Radiation Lead, RTX 4070, nvfortran 26.9, `s5_rad_sweep` at s5-gen `d054f283f0`, 1246 cases): 18 of 1246 cases differ on the device; all 18 have scheme 3; max **5632 ulp** (~1e-12 relative) over `IL` / `UIID` / `UII` / `QR` / `QR_W` / `ILW` / `INRAD_W` / `IL_S` (212 differing fields); schemes 1–2 bitwise; host builds bitwise.

## 1. Proposal (**PROPOSED pending Architect ratification**)

**Device-tier limit for scheme 3 of `s5rad_plane3d_ild`: `maxulp ≤ 8192`** (ordered-integer ulp distance, same rule as `bitcmp.py --ulp` / D-070), judged over the fields the Radiation Lead already reports under `LIBM-CLASS` (`IL`, `UIID`, `UII`, `QR`, `QR_W`, `ILW`, `INRAD_W`, `IL_S`). Schemes 1 and 2 remain **bitwise** (+0/−0 rule) on the device. The host gate stays bitwise for all schemes.

Not set in `vv-runs/gpu_gate/kernel_categories.json` until the Architect ratifies. Until then the gate has no scheme-3 exception: a device log that differs by more than the default `libm` limit of 2 ulp (D-070) would FAIL if judged as ordinary `libm`, and a `LIBM-CLASS` verdict without a ratified limit stays OPEN / not PASS.

Basis:
- Measured max on the reported sweep is **5632 ulp**. 8192 is the next power of two above that (~1.45× margin), about 1.8e-12 relative. Same style of margin used in `docs/vv/dt_class_tolerance.md` (measured max vs proposed class number), not a new policy.
- D-070’s **2 ulp** per-value rule is the default for a single libm call that does not amplify. Here a 1-ulp host/device `EXP` difference is amplified by the plane-to-plane recurrence and the negative-correction clip (Radiation Lead §9.4), so 2 ulp cannot pass the measured cases. Widening the global 2-ulp rule is not proposed.
- Relative form ~1e-11 (about 10× the reported ~1e-12) would also clear 5632 ulp on double values of order 1, but the existing device harness keys on integer `maxulp` (`ULP` line / `bitcmp.py --ulp --limit`). The proposal therefore states an integer ulp limit. If the Architect prefers a relative class bound instead, 1e-11 relative is an equivalent working number on the reported data; it is not adopted here without ratification.
- Flip-budget / Level A (D-077, `solid/09`): Level A still asks device == forced-host within 2 ulp when decisions agree. Scheme-3 sweep differences are continuous field differences from `EXP`, not decision flips; they are a **libm amplification class**, not a flip-budget class. No flip-budget change is proposed.

## 2. How this sits next to existing rulings

| Ruling / pattern | What it says | Relation to this proposal |
|---|---|---|
| D-070 / D-072 (a)(b) | Device `libm` kernels: 2 ulp per value; registry records category; gate owns the comparison | Default stays 2 ulp. Scheme 3 needs a **named exception** (this number) or it fails the default. |
| D-072 (c) / ported | Ported = device run on record within class tolerance | Scheme-3 cases cannot be signed ported under 2 ulp until a class limit exists and is met. |
| `dt_class_tolerance.md` / D-080 | Measured proposal → Architect ratification → gate field | Same process: this file is the measured proposal; gate field waits on ratification. |
| Flip-budget gate (D-077) | Decision flips, Level A/B | Not applicable to EXP amplification along the sweep. |
| Radiation Lead draft `[[libm]]` for `s5rad_plane3d_ild` (D-079, proposed in `s5rad_sweep.toml`) | Category `libm`, status proposed | Category registration is Legacy Mapper / Architect; tolerance number is this V&V proposal. |

## 3. What V&V will do after ratification

1. Record the Architect number (or a rejection) in a change-log line of this file and archive v0.1.
2. Ask Integration to encode the ratified limit in the device-tier judge (scheme-3 / `LIBM-CLASS` path of the sweep device run, or a `ULP` / class field in `kernel_categories.json` if that is the chosen home). V&V does not invent that schema here.
3. Schemes 1–2 stay bitwise; padding and counters stay exact (already reported).

## 4. Out of scope here

- Registering `s5_rad_interp` / `s5_rad_hrr` gate steps (see `vv-runs/gpu_gate/README.md`, pending Integration).
- Changing the host pin counts (done separately: `SWEEP_PIN` → `d054f283f0`).
- Inventing an Architect decision ID or treating this number as accepted.
