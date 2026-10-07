# Named RUNCMP exception for `visc_wale` pure-shear (L1356)

Owner: V&V Lead. Status: **RATIFIED (D-088)**. Measurement date: 2026-10-06.
Subject: a named, measured exception to D-080 `run_class_tolerance = 1e-11` for the pure-shear scenario of kernel `visc_wale` (`COMPUTE_VISCOSITY` CASE WALE, velo.f90:256-279, L1356; quantity `MU`). All other scenarios stay at the class gate. Does not widen D-080 globally.

Source measurement (Wall Loops, test machine, nvfortran 26.9, host ref `-Mnofma -Mnovect -Kieee`; measurement logs archived): device-versus-host `RUNCMP` on MU over coupled steps. Host bitwise already on record (19/19 PASS across 7 flag/thread × 36 cases; 5 mutants caught). L1351 / L1354 already host+device bitwise. Registry entry D-086 (`visc_wale`, `run_compare = class`) notes the pure-shear growth and asks not to silently widen.

**Change log.**
- v0.2 (2026-10-06): Architect **D-088** ratified; status **RATIFIED**; gate commit `3e2b4e0` encodes `run_case_tolerances` for `_s5_` MU = 1e-10. Device RUNCMP on record (`device_logs/stage1_vw/`, peak `_s5_` maxrel 9.3e-12). v0.1 (PROPOSED) archived at `docs/vv/archive/visc_wale_tolerance_v0.1.md`.
- v0.1 (2026-10-06): V&V proposal **(B)** — pure-shear named relative exception `1e-10` on `RUNCMP` maxrel for MU; class `1e-11` unchanged for all other scenarios. Status PROPOSED pending Architect.

## 1. Ratified limit (D-088)

**Named exception for pure-shear only:** `RUNCMP` maxrel on quantity `MU` may be at most **`1e-10` (relative)** for cases whose scenario is pure-shear (`_s5_` in the Wall driver naming: scenario 5 of `s5_visc_wale_run.F90`). Every other scenario keeps D-080 **`1e-11`**.

Match D-085 style: named, measured, reopen if exceeded; do not silently widen the global class number.

Basis:
- Measured peak on pure-shear is **4.274086e-11** (`vw_64x48x32_s5_c2`, 200 steps). Proposed `1e-10` is ≈2.3× that peak (same order of margin as D-085’s 8192 over 5632).
- Second FAIL is **2.814223e-11** (`vw_24x20x16_s5_c2`, 1000 steps), also under `1e-10`.
- Non-pure-shear max across the same logs is **5.946530e-12** (large / 1000 steps) — still under `1e-11`. The class gate is not the problem; only pure-shear is.
- Cause is physical/numerical, not a wrong kernel: for pure shear the WALE invariant `SD2` is analytically 0 (Nicoud & Ducros); residual cancellation noise under host/device coupling grows with steps and with `C_WALE` (c2 = 2.5 fails; c1 = 0.6 stays under). Per-step pure-shear delta claimed ~2.6e-16; other scenarios stay under 1e-11.

## 2. Evidence check (logs verified)

Driver scenarios (`s5_visc_wale_run.F90` `init`): 1 random, 2 zero-vel (DENOM=0), 3 tiny, 4 large, **5 pure shear**, 6 sign changes; `c1`/`c2` = `C_WALE` 0.6 / 2.5. Gate check in `run_device.sh`: `maxrel > 1e-11`.

| Log | Grid × steps | FAIL cases (only) | maxrel | All other cases |
|---|---|---|---|---|
| `dev_small.log` | 1³ / 4×3×2 / 7×5×6 × 5 | none | max 2.21e-15 | 36/36 PASS |
| `dev_200.log` | 24×20×16 × 200 | none | pure-shear peak 9.30e-12 | 12/12 PASS |
| `dev_big.log` | 64×48×32 × 200 | `vw_64x48x32_s5_c2` | **4.274086e-11** | 11/12 PASS |
| `dev_1000.log` | 24×20×16 × 1000 | `vw_24x20x16_s5_c2` | **2.814223e-11** | 11/12 PASS |

Only `_s5_` (pure-shear) fails the class gate. Pure-shear `c1` stays under (`6.80e-12` and `2.72e-12` on the two FAIL runs). Claimed headline numbers match the logs.

## 3. Is maxrel ill-conditioned because `|MU_host|≈0`?

**No.** On pure shear, `NU_EDDY` is theoretically 0, so `MU = MU_DNS + RHOP*NU_EDDY ≈ MU_DNS`. In the Wall driver, `MU_DNS` (`mud`) is drawn as `1e-5*(0.5+U[0,1])` ∈ **[5e-6, 1.5e-5]** — small, but not ≈0. Finite maxrel values (no Inf/NaN) confirm the denominator is nonzero.

Estimated absolute `|ΔMU|` at the FAIL peaks ≈ maxrel × O(1e-5) → **O(1e-16 to 6e-16)**, consistent with cancellation residue in theoretically-zero `SD2` under coupling, not with a vanishing-denominator artifact. An absolute MU floor is not required for the exception; keep the gate metric as relative maxrel (same as D-080 / `RUNCMP`).

## 4. Options considered

| Option | Verdict | Why |
|---|---|---|
| **(A)** Exempt pure-shear from the 1e-11 class gate | Reject | Drops the cancellation scenario from the gate entirely; loses the check that growth stays bounded. |
| **(B)** Pure-shear own looser tolerance (**this proposal**: `1e-10` relative) | **Recommend** | Named, measured exception (D-085 pattern). Class `1e-11` stays for all other scenarios. Reopen if a later run exceeds `1e-10`. |
| **(C)** Cap size/steps for this kernel’s gate | Reject | Hides the known growth mode instead of naming it; 1000-step runs are already inside the D-080 envelope. Silent narrowing of coverage. |

Prefer named, measured exceptions over silent widen (project preference; matches D-085 vs raising global 2-ulp).

## 5. How this sits next to existing rulings

| Ruling / pattern | Relation |
|---|---|
| D-080 `run_class_tolerance = 1e-11` | Unchanged as the class default. This is a **named exception** for one scenario of one kernel. |
| D-085 scheme-3 `maxulp ≤ 8192` | Same process: measured peak → named limit with margin → reopen if exceeded. |
| D-086 `visc_wale` libm entry | Already notes pure-shear can exceed 1e-11 and forbids silent widen; this file is the V&V number for that reopen. |
| Per-value 2-ulp / `ULP` line | Unchanged. `RUNCMP` is the step/run-level line only. |

## 6. What V&V and Integration do next (after Architect)

1. Architect ratified the pure-shear number as `1e-10` relative on MU (D-088).
2. Integration encodes the named exception in the device-tier judge / `kernel_categories.json` (case match on pure-shear / `_s5_`, or equivalent). V&V does not invent that schema here.
3. If any later pure-shear `RUNCMP` exceeds the ratified number, reopen; do not silently raise it.

## 7. Out of scope here

- Changing D-080’s global `1e-11`.
- Exempting pure-shear from host bitwise or from the `ULP` line.
- Caps on grid size or step count for the `visc_wale` gate (option C).
