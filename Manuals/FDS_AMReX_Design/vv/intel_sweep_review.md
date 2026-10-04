# Review of the Intel setup-only sweep (941 inputs): the 16 failures against GNU

Owner: V&V Lead. Status: v0.1. Date of the check: 2026-10-04.
Subject: the Intel host-gate result of the setup-only sweep (T_END=0 copies of all 941 inputs of `scope_case_list.csv`, FireX 36975d765f, Intel Release `impi_intel_rel`, sha256 `cfaa5d06...ad70`): 875 PASS, 16 FAIL, 50 SKIP.

## 1. Result

**All 16 failures are expected. None is Intel-specific.** Each of the 16 inputs was run again with the GNU Release reference (same source commit, same method). All 16 stop at set-up on GNU with the same message as on Intel (16 of 16 identical after the normalisation the sweep applies to digits). A sample of 20 inputs that pass on Intel also pass on GNU (20 of 20). So the sweep gives no finding against the Intel build; the 16 messages are valid entries for the G1 (IR-001) reference, with the reasons in section 3.

| Group | Inputs | Count | Reason (same on GNU) | Class |
|---|---|---|---|---|
| Restart "b" inputs | `clocks_restart_b`, `device_restart_b`, `geom_restart_b`, `geom_ls_restart_b`, `restart_test1b`, `restart_ulmat_b` | 6 | `ERROR(1050)`: the `.smv` file of the "a" run does not exist. A setup-only copy of the "a" input stops before it writes restart files, so no "b" input can pass in this method, on any compiler. | expected (needs files from a previous run) |
| `GEOM` keywords not in this build | `geom_azim`, `geom_elev`, `geom_scale` | 3 | `ERROR(101)` on the `&GEOM` line: the inputs use `AZIM`, `ELEV`, `SCALE`, `XYZ0`, which are not in the `&GEOM` namelist of this source (`Source/geom.f90`, `NAMELIST /GEOM/`). | expected (feature not in this build) |
| `MISC` keywords not in this build | `cube_cc_compute`, `sphere_cc_compute` (`COMPUTE_CUTCELLS_ONLY`), `zero_thick_roof` (`GEOM_DEFAULT_THICKNESS`) | 3 | `ERROR(101)` on the `&MISC` line: neither keyword exists in the source (no hit in `Source/`). | expected (feature not in this build) |
| `PART` keyword not in this build | `cloud_drag` | 1 | `ERROR(101)` on the `&PART` line: the input sets `EVAPORATE`, which is not in the `&PART` namelist of `Source/read.f90`. | expected (feature not in this build) |
| Turbulence model not in this build | `rng_32`, `rng_64` | 2 | `ERROR(129)`: `TURBULENCE_MODEL='RNG'` is not recognised (no `RNG` in the source). | expected (feature not in this build) |
| Reaction input does not balance | `lumped_stoich_soot` | 1 | `Problem with REAC 1. Unbalanced stoichiometry`: the `NU` values do not balance H, C and O (errors 0.45, 0.009, 0.25 atoms); the build refuses the reaction. | expected (the input is rejected by this build's own check) |
| | | **16** | | **16 expected, 0 Intel-specific** |

## 2. How it was checked

1. **Intel result.** Read from the test machine scratch folder of the Intel work (`intel_validation/hostgate/`): the sweep script `setup_sweep/setup_sweep.py`, the result table `setup_sweep/setup_sweep_results.tsv` (941 rows: path, class, status, reason, ranks, wall seconds; 875 PASS, 16 FAIL, 50 SKIP, no duplicate path) and the failing logs `setup_sweep/fail_logs/`, plus the "Task 4" part of the folder README. Nothing there was changed. The README states that no GNU sweep exists, and I confirmed it: `vv-runs/baseline`, `vv-runs/refbin` and the build folders hold only full-length baselines, no setup-only sweep. So a GNU run was needed.
2. **Method on the Intel side.** Each input is copied into a mirror of the repository layout (symlinks to the other files of its directory and to `Utilities/`, so `../../Utilities/...` paths resolve), `&TIME` is patched to `T_END=0.`, ranks are the `firebot_ranks` of the list (empty = 1), one core per rank, 180 s timeout. PASS = exit 0 and `STOP: Set-up only`. SKIP = more ranks than the 6 cores.
3. **GNU run, same method.** The Intel script was reused with three changes only: Open MPI instead of Intel MPI (`mpirun --bind-to none`); 3 cores (`taskset`, `nice -n 10`) on the development machine; the two 4-rank inputs (`device_restart_b`, `restart_ulmat_b`) ran oversubscribed on the 3 cores (`--oversubscribe`) instead of being skipped. Binary: `vv-runs/refbin/gnu_ompi_firex-36975d7/fds`, GNU Release, sha256 `f4214a0a...3238` (verified), source commit 36975d765f. Inputs: `git archive 36975d7 Verification Utilities` from the repository, which holds the same 941 `.fds` files as the list. The script, the file lists, both result tables and the GNU failing logs are in `vv-runs/analysis/setup_sweep_gnu_check/`.
4. **Inputs run on GNU:** the 16 failing inputs, and 20 passing ones chosen to be near the failing ones and to vary: 3 `Complex_Geometry`, 3 `Restart` (two "a" inputs and one base case), 2 `Turbulence`, 2 `Species`, 2 `Sprinklers_and_Sprays`, 3 multi-rank (`shunn3_4mesh_128`, `ht3d_sphere_48`, `dancing_eddies_tight_no_precon`, 4 ranks), and others from the DEFERRED, OUT and UNCLEAR classes (`geom_extruded_poly`, `leak_geom`, `part_drag_prof_vy`, `anca-couce-fig1_5K`, `plate_view_factor_cyl_30`). The list is in `sample20.txt`.
5. **Cause of the "feature not in this build" stops.** For each keyword I searched the namelist declarations of the source commit (`git grep` at 36975d7): `AZIM`, `ELEV`, `SCALE` appear in the geometry type but not in `NAMELIST /GEOM/`; `COMPUTE_CUTCELLS_ONLY`, `GEOM_DEFAULT_THICKNESS` and `RNG` do not appear in `Source/` at all; `EVAPORATE` is not in `NAMELIST /PART/`. A control: `cloud_drag` with a comma added after `DRAG_COEFFICIENT=1.27` (the only other suspect) fails in the same way on GNU, which leaves `EVAPORATE` as the unknown keyword. A keyword the namelist does not know gives the read error that becomes `ERROR(101)`; this is the same in both Fortran runtimes.

## 3. Results

GNU, the 16 inputs (all `FAIL` at set-up, no crash, no timeout, GNU exit status 0 for all of them):

| Input | Ranks | Intel message = GNU message |
|---|---|---|
| `Complex_Geometry/cube_cc_compute` | 1 | yes: `ERROR(101) Problem with MISC line` |
| `Complex_Geometry/geom_azim` | 1 | yes: `ERROR(101) GEOM ID=geom1 Check &GEOM input line` |
| `Complex_Geometry/geom_elev` | 1 | yes: same |
| `Complex_Geometry/geom_scale` | 1 | yes: same |
| `Complex_Geometry/sphere_cc_compute` | 1 | yes: `ERROR(101) Problem with MISC line` |
| `Complex_Geometry/zero_thick_roof` | 1 | yes: same |
| `Restart/clocks_restart_b` | 1 | yes: `ERROR(1050)` file `clocks_restart_a.smv` does not exist |
| `Restart/geom_ls_restart_b` | 1 | yes: `ERROR(1050)` `geom_ls_restart.smv` |
| `Restart/geom_restart_b` | 1 | yes: `ERROR(1050)` `geom_restart.smv` |
| `Restart/device_restart_b` | 4 | yes: `ERROR(1050)` `device_restart_a.smv` |
| `Restart/restart_test1b` | 1 | yes: `ERROR(1050)` `restart_test1a.smv` |
| `Restart/restart_ulmat_b` | 4 | yes: `ERROR(1050)` `restart_ulmat.smv` |
| `Species/lumped_stoich_soot` | 1 | yes: `Problem with REAC 1. Unbalanced stoichiometry` |
| `Sprinklers_and_Sprays/cloud_drag` | 1 | yes: `ERROR(101) Problem with PART line` |
| `Turbulence/rng_32` | 1 | yes: `ERROR(129) TURBULENCE_MODEL, RNG, is not recognized` |
| `Turbulence/rng_64` | 1 | yes: same |

GNU, the 20 passing inputs: 20 of 20 PASS (`STOP: Set-up only`), Intel PASS for the same 20. The two "a" restart inputs in the sample pass on both, which confirms the reading of the six "b" stops.

## 4. Points for the owners of the lists

1. **Eight of the 16 are class IN in `scope_case_list.csv`** (`clocks_restart_b`, `device_restart_b`, `restart_test1b`, `restart_ulmat_b`, `rng_32`, `rng_64`, `cloud_drag`, `lumped_stoich_soot`), the other eight DEFERRED. The four IN restart "b" inputs are covered by the full-length restart chains (`restart_test1a/1b` in the GNU and Intel baselines), not by a setup-only copy. `rng_32`, `rng_64`, `cloud_drag` and `lumped_stoich_soot` are in scope but cannot pass set-up on the reference build on either compiler. They should be recorded as known set-up stops in the G1 reference, and the owner of the list should decide whether they are reclassified or whether the build is expected to accept them. This is not an Intel question.
2. **Exit status.** The GNU runs exit with status 0 on all 16 set-up errors, so the exit status alone cannot detect these stops; the sweep's PASS rule (exit 0 and the `STOP: Set-up only` line) does. The Intel table does not record the exit status of the failing runs, so I did not compare it. For the G1 message diff, compare the message text, as done here.
3. **Restart "b" inputs in a setup-only sweep** will always stop with `ERROR(1050)`. If the G1 reference should be free of them, the sweep can skip inputs whose `RESTART=.TRUE.` file is the product of another input, or list them as expected stops (this review lists them).

## 5. What this does not cover

- **The 50 SKIP inputs** (35 need 8 ranks, 5 need 16, 5 need 24, one each 10, 13, 32, 36, 64) were not run on either compiler, so the sweep says nothing about them. They need a machine with enough cores (or an oversubscribed set-up-only run, which the Intel README says was not attempted).
- **GNU was not run on all 941.** The check is the 16 failures plus 20 passes. An input that fails on GNU and passes on Intel could exist among the other 855 passing inputs; the full GNU sweep is cheap on an idle machine (the Intel one took 4.5 minutes on 6 cores) if the G1 reference is to be a two-sided message diff.
- **Load.** The development machine was heavily loaded by other jobs during the GNU run (load average 48 to 67 on 8 cores), so the GNU wall times are not comparable with the Intel ones and are not reported. Only the status and the message matter here.
- **The Intel flat-layout run** (36 inputs failed with `CATF file ../../Utilities/... not found`) is superseded by the layout used in the table above and was not used.
