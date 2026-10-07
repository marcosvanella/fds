# 10. Dispositions of the second pressure hand-out (Legacy Mapper round-next)

Status: complete for this hand-out. L1365 and L1385-L1387 are host code on the private `pressure-lead-loops` branch with bitwise tests; L1391 retired; L1130 and L1121 host set-up. No shared generator sidecar entry is committed until `s5-wall-next` lands (Legacy Mapper instruction).

Scope: the loops of `docs/amrex/loop-handout-round-next.md` accepted for the AMR Pressure Solver Lead, plus L1130 and L1121 named in the same hand-off. Line numbers: survey (FireX 36975d7) and, where they moved, the FDS-AMReX head `8ad685f435` (FireX 1fb5c300fa merged). `velo.f90` shifted by +15 between those commits for the NO_FLUX and DIRECT_FORCE nests; `pres.f90` is unchanged (`git diff 36975d7 HEAD -- Source/pres.f90` is empty). Dispositions: RETIRED, HOST CODE, GENERATED KERNEL (deferred until the shared wall merge).

## 1. Summary

| Loop | file:lines (survey → head) | What | Disposition | Tests |
|---|---|---|---|---|
| L1365 | velo.f90:1404-1459 → 1419-1474 | NO_FLUX: FVX/FVY/FVZ inside obstructions | HOST CODE now (private branch); GENERATED KERNEL later (IDEMPOTENT, near `pvf_fvx`) after `s5-wall-next` | bitwise PASS (section 3) |
| L1385 | velo.f90 nest survey 896-903 / head 908-915 | DIRECT_FORCE FVX | HOST CODE now (private); GENERATED KERNEL later (scalar `FVEC(1)` gap) | bitwise PASS (section 3) |
| L1386 | survey 917-924 / head 929-936 | DIRECT_FORCE FVY | HOST CODE now (private); same | bitwise PASS (section 3) |
| L1387 | survey 938-945 / head 950-957 | DIRECT_FORCE FVZ | HOST CODE now (private); same | bitwise PASS (section 3) |
| L1391 | velo.f90:1243-1253 → 1255-1265 | VELOCITY_FLUX_CYLINDRICAL vorticity/stress | RETIRED (cylindrical deferred, D-042) | none needed |
| L1130 | pres.f90:5580-5601 | GET_H_REGFACES WALL_LOOP_1 (P4 tables) | HOST CODE (set-up only) | none needed |
| L1121 | pres.f90:4012-4026 | CHECK_UNSUPPORTED_MESH SETS_LOOP | HOST CODE (set-up only; not reached on the AMR FFT path) | none needed |

No shared file, no generator sidecar, no file of another owner and no FDS source is edited in this round. Waiting on `s5-wall-next` before any `shared:` / sidecar commit for L1365 or L1385-L1387.

## 2. Evidence per loop

### L1365 (velo.f90:1419-1474 at head): HOST CODE now; IDEMPOTENT generated form later

- What: per obstruction `N`, three face nests write `FVX`/`FVY`/`FVZ` so that the velocity inside solids is driven toward zero against the `HP` gradient (velo.f90:1423-1472). `PREDICTOR` switches `DUUDT` between `-RFODT*U` and `-RFODT*(U+US)` (and the V/W analogues). Both neighbouring cells must be `CELL%SOLID`.
- Outer loop is over obstructions, not a K,J,I nest: the generator front end rejects it as a single kernel (work-list, generator gap). The accepted shape is the same host-launch-per-box rule as `pvf_fvx` (`markers/pvf_kernels.toml`), with `idempotent = true`: overlapping obstructions that share a face write the same bits (the stored value depends only on the face and the fields, not on which obstruction wrote it). Justified in `pressure-velocity-h-loops.md` (L1365).
- AMR path today: the driver already calls the unmodified Fortran `NO_FLUX` on level 0 (`fds_p_no_flux` / `fds_step.f90:366`, used in the pressure iteration). Fine-level boxes cannot run that routine.
- This round: C++ host translations `no_flux_fvx` / `no_flux_fvy` / `no_flux_fvz` in `Source/pressure_lead` (private branch), one call per obstruction box, OpenMP over the face nest. Sidecar entry and generated kernel wait for `s5-wall-next`.

### L1385, L1386, L1387 (velo.f90:908-957): HOST CODE now

- What: `FVX -= RRHO*FVEC(1)*TIME_RAMP_FACTOR*SIN_THETA` (and FVY with `FVEC(2)*COS_THETA`, FVZ with `FVEC(3)`), `RRHO = 2/(RHOP+RHOP_neighbour)`. Perfect K,J,I nests. Ramp factor and angle are evaluated on the host before the nest (velo.f90:899-905).
- Generator gap: constant-subscript rank-1 `FVEC(1)` as a by-value scalar (hand-out). Front end does not yet accept that form.
- Called from `VELOCITY_FLUX` when `ANY(ABS(FVEC)>TWENTY_EPSILON_EB) .OR. CTRL_DIRECT_FORCE` (velo.f90:851). `FVEC` defaults to zero (cons.f90:345), so the nests are skipped unless `&WIND FORCE_VECTOR` / geostrophic / pressure-gradient force is set.
- This round: C++ host translations `direct_force_fvx` / `direct_force_fvy` / `direct_force_fvz` in `Source/pressure_lead`. Generated kernels wait on the scalar-gap feature and on `s5-wall-next` for any shared sidecar commit.

### L1391 (velo.f90:1255-1265 at head): RETIRED

- Cylindrical vorticity `OMY` and stress `TXZ` inside `VELOCITY_FLUX_CYLINDRICAL`. Cylindrical meshes are deferred in AMR mode (D-042; FR-004 rejects them with refinement). The driver aborts on them in every mode it runs (`FdsAmr.cpp:20`, `TimeLoop.cpp:290` `chk(IP[14], "CYLINDRICAL")`). `fds_kernels.f90` records that `VELOCITY_FLUX_CYLINDRICAL` is never called under that rejection. Nothing replaces L1391 in the AMR path; the FDS-only build keeps it.
- Uniform-mode parity with L1210 is out of scope for the AMR path. If cylindrical returns: the nest is a plain K,J,I cell loop (translatable-now for the generator; no feature gap of its own), but that is a later decision.

### L1130 (pres.f90:5580-5601): HOST CODE (set-up, P4)

- `WALL_LOOP_1` of `GET_H_REGFACES`: marks `LOG_INTWC` for internal `SOLID_BOUNDARY` faces, then the routine allocates and fills `REGFACE_*_H` (the K,J,I counting nests that follow). Blocker family P4 (table needed).
- Called only from GLMAT/UGLMAT set-up (`pres.f90:3871`, stage 3 of the GLMAT path). The AMR driver refuses every pressure method other than FFT (`TimeLoop.cpp:290` `chk(IP[13] != 1, "pressure solver other than FFT")`), so this loop is never reached on the AMR path today.
- Disposition: stays host (set-up / table build). Not a device candidate. No C++ translation in this round.

### L1121 (pres.f90:4012-4026): HOST CODE (set-up only)

- `SETS_LOOP` of `CHECK_UNSUPPORTED_MESH`: for each disjoint mesh set, a wall loop over `N_EXTERNAL_WALL_CELLS` looking for a Dirichlet/OPEN face, then `MPI_ALLREDUCE` on `DIRI_SET`. Set-up only; stops the run if a set has no Dirichlet face.
- Called only from GLMAT set-up (`pres.f90:3709`). Same FFT-only gate as L1130: not reached on the AMR path.
- Disposition: host (set-up). Could be retired if GLMAT/UGLMAT are permanently out of the AMR path; until that is recorded, host is the safe class. No C++ translation.

## 3. Implementation and tests (this round)

Private branch `pressure-lead-loops` under `Source/pressure_lead/` (owner-only directory; product CMake does not build it yet).

| Artifact | Loops |
|---|---|
| `PlViews.H` | shared Fortran-ordered views (`F1`, `F3`, `I3`, `Box3`) |
| `NoFluxLoops.{H,cpp}` | L1365 (`no_flux_fvx` / `fvy` / `fvz`), host launch per obstruction box |
| `DirectForceLoops.{H,cpp}` | L1385-L1387; Fortran left-associative `RRHO*FVEC*RAMP*(SIN|COS)` |
| `tests/run_noflux_tests.sh` | L1365 |
| `tests/run_dforce_tests.sh` | L1385-L1387 |

Reference = verbatim upstream line ranges (pinned by first/last statement and SHA-256). Compile with `-ffp-contract=off`.

### Test results

- **L1365:** 28 bitwise runs (6 flag sets × dpd/bind × 1/4/8 threads), 0 differences; reference files byte-identical across builds; 8 synthetic cases (predictor/corrector, overlapping obst = idempotent, JBAR=1, partial solid, empty solid, non-unit metrics); mutants 1-5 caught (5/5).
- **L1385-L1387:** 28 bitwise runs, 0 differences; reference files byte-identical; 6 synthetic cases (each axis alone, all three, zero force, JBAR=1, varying density); mutants 1-4 caught (4/4).
- Real-FDS dumps for these loops: not in this round (DIRECT_FORCE is off unless `&WIND FORCE_VECTOR` is set; NO_FLUX needs an obstruction dump hook). The residual dumps of doc 09 remain the real-FDS precedent for the directory.

## 4. Notes for other owners

- **Legacy Mapper / register:** L1365 and L1385-L1387 → in progress (host on private branch; generated form waiting on `s5-wall-next` and, for DIRECT_FORCE, the scalar `FVEC(i)` feature). L1391 → not planned (retired). L1130, L1121 → not planned (host, set-up).
- **Generator owners:** do not commit L1365 / L1385-L1387 sidecar entries until `s5-wall-next` lands. L1365 wants `idempotent = true` and a host-per-box launch like `pvf_*`. L1385-L1387 need by-value scalar treatment of `FVEC(1..3)`.
- **Driver owner:** level-0 already calls Fortran `NO_FLUX`; fine-level boxes will need the C++ host (or later generated) form when multi-level velocity flux is wired.
- **No blockers other than the wait on `s5-wall-next` for shared/generator commits.**

## 5. Commits

- Code (private branch `pressure-lead-loops`, not pushed): `4ca156050e` (L1365, L1385-L1387 host code and tests; `PlViews.H`). Prior residual commits `2a31d4562d`, `8bb1b207d5` unchanged.
- Docs (this file and the README row): `be27948`.
