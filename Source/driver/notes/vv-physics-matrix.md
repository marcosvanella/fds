# Multi-level physics V&V matrix (2 and 3 levels)

Owner: V&V Lead. Written 2026-10-14 after the owner ruling of the same day: physics verification is top priority on 2 and 3 levels.
This note changes no code. Anything not confirmed in the repo, docs, or run logs is marked **TBD**. Paths are relative to `/workspace/fds-amr/` unless they start with `Source/`.
Paths prefixed `Source/` are in `src/`.

Status key: **OPEN** = known open failure; **NOT RUN** = never run on 2 or 3 levels at the tip; **PASS** = passed on the stated tip; **SKIP** = blocked on a missing printer or input.

## 1. Conservation (FR-012 / FR-020 / FR-021)

| id | case | levels / regrid | metric | pass limit | owner | script | status |
|---|---|---|---|---|---|---|---|
| P3-C01 | `blob2d_static_mr_amr` (test-plan §5.12.2) | 2, static | M_tot relative change, per step and cumulative | per step <=1e-12; cumulative <=1e-10 | Role 3 | existing Phase 3 harness (path TBD in `vv-runs/phase3/`) | Phase 3 status TBD; NOT RUN on 3 levels |
| P3-C02 | same | 2, static | M_s per species | same limits | Species Lead | same | as C01 |
| P3-C03 | isothermal frozen blob | 2 | E_h closure | 1e-12 per step / 1e-10 cumulative | Species Lead | same | TBD |
| P3-C05 | FR-024 face-flux identity, ADV+DIF | 2 | in-code sum; permuted sum | bitwise; <=1e-14 | Role 3 | `Source/driver/tests/check_shared_face_flux.py` | TBD |
| P3-C06 | composite sums across 1/2/4 ranks, 2 `max_grid_size` | 2 | composite sums | bitwise (`EXACT_SUMS=T`) | Role 1 | same harness | TBD |
| P3-C07 | realizability | 2 | clip count; max abs(ΣY−1) | 0; <=1e-14 | Species Lead | same | TBD |
| P3-C08 | negative control: overwrite OFF | 2 | cumulative species imbalance | >1e-6 (must fail) | V&V | same | TBD |
| P3-C11 | `vv-runs/phase3/inputs/budget/budget_clip_line_noclip.fds` | 2 | CLIP_MASS, CLIP_CELLS | exactly 0 | Species Lead | harness | **PASS** on `9b93ac753c` (CLIP_MASS) |
| P3-C12 | `vv-runs/phase3/inputs/budget/budget_clip_line.fds` + harness inject | 2 | abs(ΔM−flux−CLIP_MASS)/M_ref per step | <=~1e-12 | Species Lead | harness | **PASS** on `9b93ac753c` (CLIP_MASS) |
| MC-R2 | `vv-runs/owner-show/inputs/blob_16_l0.fds --two-level-run --regrid-every 5` | 2, regrid 5 | composite mass relative change per regrid and cumulative | per regrid <=1e-12 (FR-012) | Role 3 | `vv-runs/owner-show/run_show.sh` | **PASS** on `bf0ba2d523` and on `791e3865c1`: worst 3.37e-15, final −1.86e-15 (`vv-runs/physics-matrix/RESULTS.md`) |
| MC-R3 | same blob, 3 levels | 3, regrid 5 | M_tot and each M_s per regrid / step / cumulative | 1e-12 / 1e-12 / 1e-10 | Role 3 | `--levels 3 --patch 2 13 2 13 --patch2 4 11 4 11 --regrid-every 5` (`vv-runs/physics-matrix/RESULTS.md`) | **PASS** on `791e3865c1`: worst per step 1.34e-14, final −9.4e-15 (`vv-runs/physics-matrix/RESULTS.md`); also PASS on `650cc56f5a`/`591f85a589` |
| MC-S3 | 2 and 3 levels, each species, regrid 5 | 2, 3 | M_s per species through regrids | same as FR-021 | Species Lead | driver RESULT rho*Z lines (`vv-runs/physics-matrix/RESULTS.md`) | **PASS** on `791e3865c1`: 2 levels Z1 3.1e-5, Z2 4.4e-7; 3 levels Z1 1.7e-4, Z2 1.9e-6 (residual/limit) (`vv-runs/physics-matrix/RESULTS.md`) |
| ISO-DR | isothermal derefinement (T uniform, OBST re-snap) | 2, 3; refine then derefine | T, Z, P0 and div u unchanged; pressure zones stay separate on every level; every OBST present on every level (an unresolved OBST becomes zero-thickness and never disappears); mass change recorded only | relative 1e-12 on T, Z and P0; 1e-9 on div u; total mass diagnostic only (CEO 2026-10-14) | Role 3 + Role 1 (OBST) | NEW: `vv-runs/physics-matrix/inputs/isodr_16_l0.fds`, **BLOCKED**: with uniform Z there is no tag, so the driver aborts with "expected finest level 1, got 0". Needs a force-then-derefine harness. Zone/OBST case template: `Source/driver/notes/fine-obst-walls.md` (thin XB=1.80,1.83 wall, 2 and 3 levels, `--regrid-every 5`) | **PASS** on `791e3865c1`, 2 and 3 levels, with `--init-patch` (commit `211641acd1`) and the `DEREFINE after step 5` print (`805030f653`): T, Z, P0 relative change 0 (<=1e-12); div 2.2e-16 / 6.5e-16 (<=1e-9); mass diagnostic only (`vv-runs/physics-matrix/RESULTS.md`). Scope: uniform at-rest state; OBST re-snap and pressure-zone separation not yet exercised |
| ZONE-W | thin wall, burner on one side (`fine-obst-walls.md` §check) | 2, 3; regrid 5 | D_PBAR_DT in zone without burner | 0 to round-off | Role 2 | NEW | NOT RUN |
| FAIL-Z1 | 3-level case, Z1 | 3 | Z1 ratio vs reference | 1.0 at round-off | Role 3 + Species Lead | TBD | **OPEN**: Z1 off by 1.77x. Role 3 log `scratch/role3-3lev/f3/run.log` (simple_test, 3 levels, 1760 steps). RESULT: species balance/limit 1.771, mass residual/limit 0.226 (worst over steps), fine div 7.7e-10, gate FAIL (`role3-3lev/NOTES.md`). Role 3 says 591f85a589 explains the old 0.226 mass ratio; a 3-level rerun on 591f85a589 is pending |
| FAIL-PYR | pyramid owner show (`vv-runs/owner-show/pyramid_p1_64ecfec176_full/`, tip `64ecfec176`) | 2 (nlev=2 per `movies/oct14/pyramid_render.log`) | L0 mass change | FR-020 | Role 3 | `run_show.sh` | **FIXED by `591f85a589`** (Role 3): the culprit `4d8e5f6377` primed the zeroed OBST-face velocities (U zeroed, US not) when the projection was OFF, and step-1 wall velocity was read as inflow. Pyramid L0 10 steps: residual/limit 6.2e-9 PASS. Bisect: `scratch/role3-bisect/` (p4d8, p40b, `PROPOSED-pyramid-l0-ic-zero.patch`); notes `scratch/role3-gate/NOTES.md`. The fix also explains the old 0.226 mass and ~0.6 Z2 ratios. V&V regate on 2 levels: TBD |
| FAIL-VR | `vv-runs/owner-show/vortex_ring_re100_wall/` (overnight 2026-10-06, binary `build/ompi_gnu_release_s6fold/fds_amr`) | 2, adaptive | composite mass drift | 1e-10 cumulative | Role 3 (investigating) | `run_show.sh` | **OPEN**: about 5.7% at step 1157, when the run aborted because the disk was full (`RESULT.md` §Disk rescue). Rerun in a new outdir `vortex_ring_re100_wall_pe20/` |

| ST-TT2 | simple_test `--tag-temp`, 2 levels, after the `591f85a589` fix (Role 3, `scratch/role3-gate/` r2l) | 2, regrid 5 | fine div; mass and species residual/limit | 1e-9; <=1 | Role 3 | driver RESULT | **PASS** as reported by Role 3: div 1.6e-10, mass 1.1e-6, species 5e-7. Not regated by V&V |

## 2. Enthalpy / divergence

| id | case | levels / regrid | metric | pass limit | owner | script | status |
|---|---|---|---|---|---|---|---|
| DIV-ISO2 | isothermal closed box | 2 | fine max abs(div u − D) | <=1e-9, with Poisson L2 residual <=1e-12 | Role 2 | NEW (input NEW) | NOT RUN |
| DIV-ISO3 | same | 3 | same | same | Role 2 | NEW | NOT RUN |
| DIV-ADI2 | adiabatic closed box, heat source | 2 | div u − D per level; D_PBAR_DT consistent | 1e-9 / 1e-12 | Role 2 + Species Lead | NEW | NOT RUN |
| DIV-ADI3 | same | 3 | same | same | Role 2 | NEW | NOT RUN |
| DIV-E2L | Role 3 energy 2-level (`scratch/role3-energy2l/amr2/run.log`) | 2 | L1 max abs(div u − D) | 1e-9 | Role 3 | existing run dir | **PASS**, informally: step 182 L1 4.31e-11, L0 9.76e-12. The same log shows composite mass change −2.27e-7, constant from step 35 (check against FR-020) |
| DIV-BLOB | blob adaptive (`RESULT.md` item 11) | 2, regrid 5 | worst max abs(div u − D) | 1e-9 | Role 3 | `run_show.sh` | **PASS** on `791e3865c1` (2 levels 6.7e-13 at step 2; 3 levels 2.0e-12 at step 1, neither a regrid step; driver gate PASS). It was 5.8e-4 / 7.9e-4 on `591f85a589` and 1.436 on `bf0ba2d523` (`vv-runs/physics-matrix/RESULTS.md`) |
| FAIL-DIVG | fine level, all cases since commit `161671e402` | 2 | growth of fine div u − D | 1e-9 | Role 2 + Role 3 | bisect from `161671e402` (NEW) | **CLOSED** (CEO 2026-10-14): the growth came from Role 2's uncommitted WIP, so it was not a regression |
| FAIL-VR-DIV | `vortex_ring_re100_wall` | 2 | L1 max abs(div u − D) | 1e-9 | Role 3 | `run_show.sh` | **OPEN**: about 12–19.5, spikes to about 23 (`RESULT.md`) |

## 3. Radiation

| id | case | levels / regrid | metric | pass limit | owner | script | status |
|---|---|---|---|---|---|---|---|
| RAD-Q2 | `ild_3d_gray`-style box, fine single-mesh FDS reference | 2 | composite Q_RADI vs fine single mesh (L2/max rel) | TBD (proposal: same as the D-073b gate; limit TBD) | Radiation Lead | `Source/driver/tests/radiation/compare_runs.py` | NOT RUN |
| RAD-Q3 | same | 3 | same | TBD | Radiation Lead | same | NOT RUN |
| RAD-CF | intensity continuity at C/F faces (spec `docs/radiation/01-radiation-amr-spec.md`: fine ghost IL prolonged from coarse, same angle set; the C/F behaviour is marked [VERIFY] there) | 2, 3 | jump in IL across C/F per angle | TBD | Radiation Lead | NEW | NOT RUN |
| RAD-CFO | thin OBST at a C/F face (spec §candidates) | 2 | Q_RADI, IL | TBD | Radiation Lead | NEW | NOT RUN |
| RAD-S3 | device scheme-3 `s5rad_plane3d_ild` (D-085, `vv/rad_scheme3_tolerance.md`) | single | maxulp | <=8192 ulp | Radiation Lead | `vv-runs/gpu_gate/kernel_categories.json` judge | **PASS**: landed (`docs/radiation/03` :189). Not a multi-level test |
| RAD-D073b | D-073b gate | TBD | TBD | TBD | Radiation Lead | TBD | TBD: no D-073b text found under `docs/` |
| RAD-FIRE | `race_test_1_r4`, `race_test_4_r4` (A-35 4:1), fire + radiation + C/F | 2 | Q_RADI and HRR vs FDS | TBD | Radiation Lead | `compare_runs.py` | NOT RUN |

## 4. Accuracy

| id | case | levels / regrid | metric | pass limit | owner | script | status |
|---|---|---|---|---|---|---|---|
| ACC-ST-HRR | simple_test (`work-species/simple_test/*/simple_test_s2*.fds`) vs FDS coarse and fine | 2 | HRR(t), Q_* budget columns | between coarse and fine, or within TBD % of fine | Species Lead | `Source/driver/tests/radiation/compare_runs.py` (the brief's `scratch/r1/compare/` does not exist) | NOT RUN |
| ACC-ST-PROF | same | 2, 3 | T and velocity profiles (devc) | TBD | V&V | same + `Source/driver/tests/compare_run.py` | NOT RUN |
| FAIL-BURN | simple_test with the burner on the fine level (`--tag-obst 1 --nofine -1 --patch 0 7 0 7`), `8de7213eff` + patch 0026. Data on an NVIDIA test machine (cited as given, not readable from the box): `bf5` (two-level), `tall` (patch 0 7 0 15), `ff` (fine FDS, default FFT), `ugf` / `ugc` (clean reference, plain FDS `&PRES SOLVER='UGLMAT'`), `uglcmp.py`, `pdt`, `ffprofT`. Coarse references on the box: `work-species/simple_test/` | 2 | E 0-2 s; mean HRR for t>=2 s; E 0-1 s | fixed tolerance against the clean reference (UGLMAT = FFT with VELOCITY_TOLERANCE=1E-7, MAX_PRESSURE_ITERATIONS=200). Tolerance value TBD (proposed: mean HRR t>=2 s within 5%; E 0-2 s within a fixed percent to be set by the owner). The earlier scatter-band criterion is withdrawn | Species Lead (+ Role 2) | the test machine `uglcmp.py`; `Source/driver/tests/run_fine_burner_check.sh` | **OPEN**. The two-level run FAILS E 0-2 s: bf5 155.0 kJ, 16% low; tall patch 0 7 0 15 132.4 kJ, 29% low. It PASSES mean HRR for t>=2 s: +1.6% (tall +4.4%). The patch top is ruled out: the tall run matches bf5 to 0.1 kW up to 1.5 s. Against the clean reference, E 0-1 s is 35.6 kJ two-level (41.5 with a forced first dt) versus 50.6 kJ (default FFT gave 54.5). The deficit is about 9-15 kJ, all from 0.4 s on (departs after ~0.3 s). The earlier default-FFT fine-FDS scatter band (E 0-2 s 185.35-185.62 kJ) is an artifact: default FFT has a burner-face pressure artifact before 0.35 s. Earlier numbers vs default-FFT fine FDS (E 0-2 s: coarse 197.7, L0 port 197.7, fine 185.4, two-level 155.0 kJ; 0-5 s: 671.0 / 670.9 / 665.0 / 641.3 kJ; mean HRR t>=2 s 157.8 / 157.7 / 159.9 / 162.1 kW) are superseded as the reference. Gates at 5 s PASS: mass 0.519x, species 0.558x, fine div 7.8e-12 |
| ACC-VR | vortex ring `Source/driver/tests/cases/vortex_ring_{16_l0,32_l0,32_4m,32_8m,48_4m,48_8m}.fds` | 2, 3 | vorticity, core trajectory, KE vs FDS 32/48 | TBD | Role 3 | `compare_runs.py` + `vv-runs/owner-show/movies` core log | PRIMARY PASS on `50c4f0588f` (smoke only; no accuracy metric) |
| ACC-GC | vortex ring 16/32/48 and simple_test, AMR vs uniform | 2, 3 | observed order | about 2 (TBD) | V&V | NEW | NOT RUN |
| ACC-QB | Q_* budget closure (P3-C04/C10, FR-022 form) | 2 | closure | max(1.05·abs(base), 0.005·max abs(Q_TOTAL)) | Species Lead | existing (`Energy_Budget_*` cases) | NOT RUN, Phase 4 |

**Row count:** Conservation 18, Enthalpy/div 8, Radiation 7, Accuracy 6 (total 39; 3 are OPEN failures (FAIL-Z1, FAIL-BURN, FAIL-VR/FAIL-VR-DIV not rerun), FAIL-PYR is FIXED, DIV-BLOB and ISO-DR now PASS, 1 is CLOSED).

## 5. Execution order and assignments

1. **Triage the open failures (day 0, cheap).** FAIL-DIVG is closed (Role 2 WIP). Role 3 owns the vortex_ring_re100_wall mass/div results (FAIL-VR, FAIL-VR-DIV). Role 3 and the Species Lead take FAIL-Z1. Role 3 takes FAIL-PYR. The V&V Lead first recovers the log and input for each of these (all TBD now).
2. **Conservation gates.** MC-R2 is already passing and becomes the regression floor. Then MC-R3, MC-S3, P3-C01..C08 on 2 and 3 levels, ISO-DR and ZONE-W. Role 3 and Role 1 own OBST/zone representation; Role 2 owns ZONE-W. The V&V Lead runs the gates.
3. **Divergence.** Role 2 runs DIV-ISO/ADI at 2 and then 3 levels. This needs the Poisson L2 residual 1e-12 printed alongside, (FAIL-DIVG closed). DIV-BLOB (1.436) is regated once that fix lands.
4. **Radiation.** The Radiation Lead runs RAD-Q2/Q3, RAD-CF, RAD-CFO and RAD-FIRE, and must first define the tolerances and locate D-073b.
5. **Accuracy.** The Species Lead runs ACC-ST, Role 2 runs FAIL-BURN, Role 3 runs ACC-VR, and the V&V Lead runs ACC-GC and the composite report through `compare_runs.py`.

## 6. Gaps (NEW / TBD)

- New inputs or scripts are needed for: 3-level runner (MC-R3/MC-S3), ISO-DR, ZONE-W, isothermal and adiabatic closed boxes, RAD-CF/RAD-CFO, ACC-GC.
- `scratch/r1/compare/compare_runs.py` does not exist. The real script is `Source/driver/tests/radiation/compare_runs.py`.
- Logs for Z1 1.77x and pyramid 51,562x (Role 3, `4d8e5f6377`) and the burner HRR table (Species Lead, `8de7213eff`, the test machine) are not on the box yet. Build fixed at `650cc56f5a`.
- D-073b is not found under `docs/`. Tolerances for radiation and accuracy are TBD.
- The Phase 3 harness path for P3-C01..C08 is not pinned here (TBD).
