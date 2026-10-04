# R4 part 2 and R2b: registry integration, moving-blob gate rows, GPU tagging path, domain-edge status

Inventory and status. Numbers are from the box CPU build (gfortran, host) unless stated.

## 1. Registry integration: what is done, what runs against Role 1's entry points, what is blocked

Done and tested (ctest, 1 and 4 ranks):
- `RegistryTransfer` (Role 1) as the `LevelDataTransfer` of `RegridAmrCore`: `fill_initial_level` at t = 0 (initial conditions evaluated directly on the new level), `fill_new_level` / `fill_remade_level`
  (conservative limited prolongation of rho*Z, rho = species sum, faces by the D-060 path, bitwise copy of the previous fine data over the overlap), `hierarchy_done` (mass-weighted average-down, D-062).
  Agreement with `CellTransfer` / `FaceTransfer` / `LevelOps`: `test_species_avgdown`, `test_blob_registry` (real registry, mock transport).
- D-060 face path (`FaceTransfer`, normal-linear prolongation, divergence of every all-new fine cell equals its parent's to round-off), D-062 mass-weighted restriction, D-063 post-regrid projection with the
  pressure_backend composite MLMG (`PressureBackendSolver`, setting `POST_REGRID_PROJECTION`, AUTO / ON / OFF).
- `LevelRegistry::begin_regrid` / `end_regrid`: `RegridAmrCore::regrid_dynamic` opens and closes one bracket per regrid. `test_blob_registry` now asserts it after every run: brackets = regrids, bracket closed,
  nothing left retired, registry levels = hierarchy levels, `make_level` count covers the finest level.
- `LevelRegistry::fill_initial_level`: used by `init_from_tags` for every new level at t = 0 (`RegistryTransfer::initial` set in the test).
- Blob regrid test on the real registry with the mock conservative transport: 3-D 32^3, 2-D 48x1x48 three levels, 2-D ratio 4, ns2d_16-style 16x1x16.

Role 1 side as found on the branch (git log) and in the shared working tree:
- Committed: `bind_level` / `unbind_level` (levels bound in order; only the top bound level can be rebound), `FluxStages` D-061 order, `RegistryTransfer`, registry bracket, flux read-out and override hooks.
- Not committed (working tree of Role 1, seen read only): `TwoLevelRun` (`--two-level-run`: a static patch made by `RegridAmrCore`, listener = `LevelRegistry`, data = `RegistryTransfer` with a `derive` hook, bound by `bind_level`,
  `install_post_regrid_projection` for the initial projection), `CompositeOps`, a mirror fill for the domain-edge ghost layers of fine boxes (`mirror_domain_edges`), `HypreStub`. This is the first user of the regrid hook.

Blocked on Role 1 (nothing to wire on our side until it is committed):
1. Fine-box domain-edge ghost layers (A-60): working-tree change only, not on the branch (item 4 below).
2. A dynamic regrid during a time loop run: `bind_level` rebinds only the top level, so a regrid that remakes a level below an existing finer one needs unbind of the finer levels first; `TwoLevelRun` uses
   `REGRID_INTERVAL=0` (one static patch). No time-loop call of `RegridAmrCore::regrid_dynamic` between steps, and no call of `install_post_regrid_projection` from the time loop itself.
3. Derived fields of a new level in the registry transfer: the stage arrays D, DS, RSUM, MU, KRES, H, HS are not transferred; `RegistryTransfer::derive` (EOS, `run_divergence_part1`) is only in the working tree.
4. Known gaps listed in `driver/notes/level-binding.md`: predictor DS differs from level 0 by up to 0.4 % for non-uniform species on a bound level, `INTERPOLATED_MESH` mask, one pressure zone, composite D_PBAR_DT.
5. Real FDS stages in a moving-blob regrid run: needs 2 and 3. Until then the moving blob uses the mock transport (item 2).

## 2. Moving-blob gate rows (V&V plan section 5.12; thresholds as written there, none assumed)

Run: `ctest -R moving_blob_p3` (`test_blob_registry --gates`, 1 and 4 ranks; the whole `test_blob_registry` includes the same checks). Transport: the mock conservative upwind scheme on the real registry,
because the real FDS stages need a driver build with fine-level binding and dynamic regrids (item 1, blocked); the real-stage static-patch rows stay in `tests/run_e2e_driver.sh` (skipped without a driver build).

| Row | Check | Threshold | Result |
|---|---|---|---|
| P3-B01 | composite mass and species per step; cumulative; clips; containment | 1e-12; 1e-10; 0; 0 | per step <= 2e-16, run <= 3e-16, clips 0, containment 0 (ratio 2 cases) |
| P3-C07 | max abs(sum Z - 1) | 1e-14 | <= 5.6e-16 |
| P3-B02 | relative L2 of rho*Z_1 against the uniform-fine run, AMR <= 0.5 x uniform coarse | 0.5 | 0.34 (3-D), 0.35 (2-D three levels); ratio 4 reported only |
| P3-B03 | circular centroid against the uniform-fine run | 0.1 fine cell | 0.039, 0.007 |
| P3-B04 | repeat determinism (hierarchy and data bitwise); ranks by `run_regrid_rank_check.sh` | identical | identical |
| P3-B05 | `REGRID_INTERVAL=0` control: containment > 0 and discrimination exceeded | inverted | containment 2947 cells, L2 ratio 0.79 > 0.5 |
| P3-F03 | no tags, AMR mode: no level-1 boxes, data bitwise equal to the single-level run | bitwise | pass |
| P3-F04 | everything tagged: level 1 bitwise equal to the uniform fine run; level 0 equal to the average-down of level 1 | bitwise; 1e-14 | 0; 1.4e-16 |
Not covered here (not possible with this transport): P3-C03 enthalpy (isothermal blob, no energy equation in the mock), P3-C05/C06 flux identity and EXACT_SUMS decomposition (FluxTests and the rank-hash test),
P3-B09 corner front distance (`corner_check`: bitwise outside n+1 cells for n steps of the upwind scheme; the 4n of the driver stencil is for the real stages), P3-D/P3-X rows.

## 3. GPU path of the tagging kernels (K2)

- `rt_tag_kernels.F90` now follows the pattern of the generated K2 files (`s4_omp.inc` / `s5gen_k2.F90`): macro tails `S4_LOOP` (`target teams loop collapse(3)`), `S5_LOOP_CALLEE`
  (`target teams distribute parallel do collapse(3)`, the region that calls the declare-target `cell_value` / `ndiff`), `S4_DEV` (`is_device_ptr` on nvfortran, `has_device_addr` otherwise); host build `parallel do collapse(2)`.
  Switches `S4_OFFLOAD` (or `RT_OFFLOAD`), `S4K2_COLL`, `S4K2_THREADS`, `S4K2_TEAMS`, `S5_FORCE_DPD`. Bounds and direction flags are copied to scalars before each region (array dummies in a region would be mapped as copies).
  `rt_tag_count` is host only (the approved clause list has no sum reduction).
- Build-time choice for the library: `-DRT_TAG_OFFLOAD=ON [-DRT_TAG_OFFLOAD_FLAGS="-mp=gpu;-gpu=mem:managed,nofma"]` compiles `rt_tag_kernels.F90` with `RT_OFFLOAD` (same entry points and arguments); `tag_kernel_check` builds both variants
  and compares them bitwise (ctest `regrid_transport_tagkernel_host_vs_offload`). Device use needs an AMReX with a GPU backend; the installed one is CPU only (the configure step warns).
- `docs/tools/k2_ci_check.py` on this file: K2-01 to K2-06 and the clause lint pass; the three K2-07 findings remain (header note naming an upstream routine and lines): FDS has no tagging routine, so none exists to name. The
  header states this; a registry rule for kernels without an upstream routine is a decision for the tool owner.
- Laptop run (RTX 4070, nvfortran 26.9, `-O2 -mp=gpu -gpu=mem:managed,nofma`, no fast math): 38 cases, output lines identical to the gfortran host build, 38 kernel launches seen, `OMP_TARGET_OFFLOAD=MANDATORY` passes.

## 4. Domain-edge ghost layers of fine boxes (E3b / E1)

Role 1's fix (A-60) is not on the branch (no commit since the E3b analysis; only a working-tree change). `MASS_TOL` (5e-6) and `E1_ZZ_TOL` (1e-4) in `tests/run_e2e_driver.sh` are unchanged. When it lands: rebuild the driver,
rerun `run_e2e_driver.sh` on 1 and 4 ranks, set `MASS_TOL` to `E3A_TOL` (1e-12) and `E1_ZZ_TOL` to 1e-12.
