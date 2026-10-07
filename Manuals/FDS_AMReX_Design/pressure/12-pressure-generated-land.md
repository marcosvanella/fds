# 12. Generated forms of the pressure-lead loops (post s5-gen FF)

Status: four nests landed as generated K2 kernels with sidecar entries and bitwise tests; three groups stay host code with the generator gap named. Pins: survey FireX 36975d7 (sidecar line numbers); s5-gen tip `6dca8b2209` (Wall L1356 `visc_wale` `8dcb39998c` included; Mesh Data particle tables `0185448a9c` included). Host code for every loop remains on the private `pressure-lead-loops` branch (docs 10, 11).

## 1. Summary

| Loop | velo.f90 (survey) | Disposition | Sidecar name | Tests |
|---|---|---|---|---|
| L1351 | 115-121 CONSMAG MU | GENERATED KERNEL | `mu_consmag` | bitwise PASS |
| L1354 | 187-196 Deardorff MU | GENERATED KERNEL (`C_DEARDORFF` by value) | `mu_deardorff` | bitwise PASS |
| L1360 | 413-435 COMPUTE_STRAIN_RATE | GENERATED KERNEL (explicit `private` list closes the PRIVATE-list gap) | `strain_rate` | bitwise PASS |
| L1379 | 970-978 CORIOLIS UP/VP/WP | GENERATED KERNEL | `coriolis_up` | bitwise PASS |
| L1365 | NO_FLUX FVX/FVY/FVZ | HOST CODE; generated form needs an obstruction-box builder (as `s5_pvf.py`) plus `idempotent` | none | host bitwise PASS (doc 10) |
| L1385-L1387 | DIRECT_FORCE FVX/FVY/FVZ | HOST CODE; generator refuses constant subscript `FVEC(n)` | none | host bitwise PASS (doc 10) |
| L1381-L1383 | CORIOLIS FVX/FVY/FVZ | HOST CODE; generator refuses constant subscript `OVEC(n)` | none | host bitwise PASS (doc 11) |

## 2. Landed kernels

- Sidecar (`s5_markers.toml`): four `[[kernel]]` entries, `layout = "exact"` (FDS-index dummies on the field arrays as allocated: no gather, no repack, no copy). `C_DEARDORFF` added to `[policy.scalars_by_value]`. `mu_deardorff` and `strain_rate` carry an explicit `private` list of their scalar temporaries.
- Generated signatures (all `bind(c)`, cell nest `K,J,I` over `1..KBAR,1..JBAR,1..IBAR`):
  - `s5gen_mu_consmag(IBAR,JBAR,KBAR,CSD2,MU,MU_DNS,RHOP,STRAIN_RATE)`
  - `s5gen_mu_deardorff(IBAR,JBAR,KBAR,C_DEARDORFF,LES_FILTER_WIDTH,MU,MU_DNS,RHOP,UP,UP_HAT,VP,VP_HAT,WP,WP_HAT)`
  - `s5gen_strain_rate(IBAR,JBAR,KBAR,RDX,RDY,RDZ,STRAIN_RATE,UU,VV,WW)`
  - `s5gen_coriolis_up(IBAR,JBAR,KBAR,UP,UU,VP,VV,WP,WW)`
- Golden signatures: four new keys, no drift in existing keys.
- Simplicity / data layout (D-084): each kernel is the upstream statement list and nothing more (no extra FLOP, branch or temporary). One pass per nest, unit stride in `I`, read-only inputs and a single written field. `exact` layout keeps the FDS arrays in place, so there is no host or device copy. `coriolis_up` writes UP/VP/WP in one sweep over UU/VV/WW. Its FV updates (L1381-L1383) stay host code for now, so those arrays are read in a second sweep there. Fusing the two would remove that sweep once `OVEC(n)` is accepted.

## 3. Tests (`test/test_pressure_lead.py`, box CPU)

- Reference: the upstream loop text, copied by line range from 36975d7 (SHA-256 pinned per slice), compiled in a routine with the generated dummies. Inputs are shared; outputs are compared bit for bit over the whole array, ghosts included. Grid 13x11x9.
- Matrix: six flag sets (O0, O2, O0omp, O2omp, O2omp_off, O2omp_dpd) x callee switch dpd/bind x 1/4/8 threads (OpenMP sets), `-ffp-contract=off`.
- Mutants of the generated file, all caught: CONSMAG product reassociated; Deardorff `VP_HAT` read as `UP_HAT`; strain divergence reassociated; strain DUDY stencil; Coriolis `UU(I-1)` and `WW(K-1)` stencils.
- Result: 36/36 PASS (2 text checks, 28 bitwise runs, 6 mutants) on the regenerated file at the s5-gen tip.

## 4. Gaps kept as host code

- L1385-L1387 (DIRECT_FORCE): s5gen refuses a rank-1 array with a constant subscript (`FVEC(1)`): "rank-1 array FVEC subscript 1 is not a loop variable +/- constant". Feature needed: constant-subscript rank-1 read passed as a by-value scalar (`FVEC(i)` hoisted to an argument).
- L1381-L1383 (CORIOLIS FV): same refusal for `OVEC(2)` / `OVEC(3)`. Same feature.
- L1365 (NO_FLUX): outer loop over obstructions with per-obstruction box bounds and CELL_INDEX alias; needs an obstruction-box builder in the style of `s5_pvf.py` and an `idempotent` entry. Not attempted in this land. The Wall owner's `wall_no_flux` (L1366, NO_FLUX wall loop) landed with `unique` + `unique_check`; the obstruction loop could follow the same marker pattern once a box builder exists.

## 5. Register notes for the Mapper

- L1351, L1354, L1360, L1379: GENERATED (sidecar names above). Host versions on the private branch stay as cross-checks.
- L1385-L1387, L1381-L1383: HOST; blocker = constant-subscript rank-1 scalar (`FVEC(n)`, `OVEC(n)`).
- L1365: HOST; blocker = obstruction-box builder + idempotent.
- Shared-file commit: `6dca8b2209` on s5-gen (rebase of `b7ebb72338` onto `8dcb39998c`), under the s5-gen lock: sidecar, regenerated K2/H/args/report, golden (four new keys only; no drift including Wall `visc_wale`), and `test/test_pressure_lead.py`. Private branch `pressure-lead-s5gen` matches that tip.
- If s5-gen moves on generator inputs again: re-apply the four sidecar entries and the `C_DEARDORFF` policy name, regenerate with `--update-golden`, check that the only new golden keys are the four above, then rerun `test/test_pressure_lead.py`.
