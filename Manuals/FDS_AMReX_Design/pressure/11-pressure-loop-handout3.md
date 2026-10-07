# 11. Dispositions of the third pressure hand-out (higher-share open list)

Status: complete for this hand-out. L1360, L1354, L1351, L1379 and L1381-L1383 are host code on the private `pressure-lead-loops` branch with bitwise tests. COMPRESSION_WAVE L1276-L1287 left to Mesh Data / Mapper. No shared generator sidecar until `s5-wall-next` / FF.

Scope: higher-share open velocity-side loops claimed by the AMR Pressure Solver Lead. Line numbers: survey (FireX 36975d7) and head `8ad685f435` where they moved (`velo.f90` +9 for viscosity/strain, +12 for Coriolis). Nest bodies are identical at both pins. Dispositions: HOST CODE now; GENERATED KERNEL deferred pending generator features and `s5-wall-next`.

## 1. Summary

| Loop | file:lines (survey → head) | What | Disposition | Tests |
|---|---|---|---|---|
| L1360 | velo.f90:413-435 → 422-444 | COMPUTE_STRAIN_RATE default model | HOST CODE now (PRIVATE-list gap, same as L1355) | bitwise PASS |
| L1354 | velo.f90:187-196 → 196-205 | COMPUTE_VISCOSITY Deardorff MU nest | HOST CODE now (translatable-now; generated later) | bitwise PASS |
| L1351 | velo.f90:115-121 → 124-130 | COMPUTE_VISCOSITY CONSMAG MU nest | HOST CODE now (translatable-now; generated later) | bitwise PASS |
| L1379 | velo.f90:970-978 → 982-990 | CORIOLIS_FORCE cell-centred UP,VP,WP | HOST CODE now (translatable-now; generated later) | bitwise PASS |
| L1381 | velo.f90:992-1000 → 1004-1012 | CORIOLIS FVX += 2(OVEC(2)WBAR-OVEC(3)VBAR) | HOST CODE now; GENERATED KERNEL needs constant-subscript `OVEC(n)` | bitwise PASS |
| L1382 | velo.f90:1006-1014 → 1018-1026 | CORIOLIS FVY update | HOST CODE now; same feature note | bitwise PASS |
| L1383 | velo.f90:1020-1028 → 1032-1040 | CORIOLIS FVZ update | HOST CODE now; same feature note | bitwise PASS |

## 2. Evidence per loop

### L1360 (velo.f90:422-444 at head): HOST CODE (PRIVATE-list gap)

- Default-model strain rate `STRAIN_RATE(I,J,K)` from face velocities `UU/VV/WW` and metrics `RDX/RDY/RDZ`, with `ONTH = 1/3` (prec.f90). Locals `DUDX`..`DWDY` are assigned every pass; the generator front end rejects them as possibly live after the loop (same PRIVATE-list gap as the Mapper's L1355 Vreman nest).
- Implemented: `compute_strain_rate` in `StrainRateLoops.{H,cpp}`.

### L1351 (velo.f90:124-130): HOST CODE (translatable-now)

- CONSMAG/DYNSMAG: `MU = MU_DNS + RHOP*CSD2*STRAIN_RATE`. Plain cell nest; front end accepts. Host now; sidecar after `s5-wall-next`.
- Implemented: `viscosity_consmag`.

### L1354 (velo.f90:196-205): HOST CODE (translatable-now)

- Deardorff MU nest: `DELTA = LES_FILTER_WIDTH`; `KSGS` from `(UP-UP_HAT)**2` etc.; `NU_EDDY = C_DEARDORFF*DELTA*SQRT(KSGS)`; `MU = MU_DNS + RHOP*NU_EDDY`. Inputs `UP`/`UP_HAT` are formed by earlier nests and `TEST_FILTER` (not this loop).
- Implemented: `viscosity_deardorff`.

### L1379 (velo.f90:982-990): HOST CODE (translatable-now)

- Cell-centred `UP/VP/WP = 0.5*(face + neighbour face)` into WORK7-9. Front end accepts. Same formula as the Deardorff pre-filter nest (not claimed here).
- Implemented: `coriolis_up_vp_wp`.

### L1381, L1382, L1383: HOST CODE now; generator feature note

- Face nests: `FVX += 2*(OVEC(2)*WBAR - OVEC(3)*VBAR)` and the Y/Z analogues. `OVEC` defaults to zero (cons.f90), so the nests are skipped unless `&WIND CORIOLIS_VECTOR` is set.
- Generator gap: constant-subscript rank-1 `OVEC(n)` as by-value scalars — the same feature that unlocks L1385-L1387 (`FVEC`). Documented in `pressure-velocity-h-loops.md`. Host translations pass `ovec1/2/3` as scalars.
- Implemented: `coriolis_fvx` / `coriolis_fvy` / `coriolis_fvz`.

## 3. Implementation and tests

Private branch `pressure-lead-loops`, `Source/pressure_lead/`.

| Artifact | Loops |
|---|---|
| `StrainRateLoops.{H,cpp}` | L1360 |
| `ViscosityLoops.{H,cpp}` | L1351, L1354 |
| `CoriolisLoops.{H,cpp}` | L1379, L1381-L1383 |
| `tests/run_strain_tests.sh` | L1360 |
| `tests/run_visc_tests.sh` | L1351, L1354 |
| `tests/run_coriolis_tests.sh` | L1379, L1381-L1383 |

Reference = verbatim upstream nests, pinned by first/last statement and SHA-256. Compile with `-ffp-contract=off`.

### Test results

- **L1360:** 28 bitwise runs, 0 differences; 5 synthetic cases; mutants 1-2 caught.
- **L1351/L1354:** 28 bitwise runs, 0 differences; 6 cases (CONSMAG and Deardorff, JBAR=1, C_DEARDORFF=0); mutants 1-2 caught.
- **L1379/L1381-L1383:** 28 bitwise runs, 0 differences; 5 cases (including zero OVEC); mutants 1-4 caught.
- Real-FDS dumps: not in this round (Coriolis off by default; strain/viscosity need a dedicated dump hook).

## 4. Notes for other owners

- **Legacy Mapper:** L1360/L1354/L1351/L1379/L1381-L1383 → in progress (host on private branch). Generated form waits on `s5-wall-next` and, for L1360, the PRIVATE-list feature; for L1381-L1383, the constant-subscript `OVEC(n)` feature (shared with L1385-L1387).
- **COMPRESSION_WAVE L1276-L1287:** not touched (Mesh Data / Mapper).
- **No blockers other than the wait on `s5-wall-next` / FF for shared/generator commits.**

## 5. Commits

- Code (private branch pressure-lead-loops, not pushed): e36141be1e.
- Docs (this file and the README row): 3621d15 (content); hash-stamp follow-up.
