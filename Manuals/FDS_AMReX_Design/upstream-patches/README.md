# Upstream patches (decision D-051)

Proposed changes to upstream FDS source (FireX or master) that come out of the AMR/GPU work, for the owner to review and commit upstream himself. Nobody else commits or pushes them upstream.

## Rules for a patch
1. One short-scope change per patch (one marker set, one loop cleanup, one routine). No mixed patches.
2. Delivered as a file in this folder: `git format-patch` output or a unified diff against FireX, applicable to the upstream source file as it is, without our `WITH_AMREX` code.
3. A short plain-language rationale in the patch description: what changes, why the GPU or AMR work needs it, what it does not change.
4. A behavior-unchanged check: the exact command or test and its result (for example the V&V baseline comparison against the reference binary, bitwise or within the stated tolerance), run on the patched upstream file.
5. A target line: **FireX**, **master**, or **both**. Master is merged into FireX periodically, so write the patch to apply cleanly to both and say which one it was tested against.
6. Name: `NNNN-short-name.patch` (four-digit sequence), and add a row to the index below in the same change.

## Standard behavior-unchanged check (V&V Lead proposal, accepted as the default)
1. Control build: the reference GNU Release binary `vv-runs/refbin/gnu_ompi_firex-36975d7/fds`. Patched build: the same FireX source plus the patch, same flags, HYPRE and toolchain. Compare patched against control, never against older recorded baselines. A master-targeted patch is also built once on master against a master control.
2. Cases (about 5 minutes): ns2d_16 on 1 rank, obst_activation_default on 4 ranks, shunn3_4mesh_32 on 4 ranks, csmag_32 on 1 rank, plus one case that exercises the changed lines, with a coverage line in the report. A patch that touches shared code also runs the Tier 1 set (39 runs, about 17 minutes).
3. Pass criterion: bitwise (`cmp`-identical outputs apart from wall-clock and CPU columns and timing lines), run with `SIG_FIGS=17`. A patch declared to change results states which outputs change and why, and gets a stated tolerance metric instead.
4. Loop or bounds patches: one extra run of the exercising case on a `-fcheck=all` build, which must be clean.
The V&V Lead is writing `vv-runs/tools/patch_check.sh` so each patch can cite the command and its output.

## Index
| No. | File | Target | Source file(s) | Summary | Behavior-unchanged check | Author | State (proposed / reviewed / committed upstream by owner) |
|---|---|---|---|---|---|---|---|
| 0001 | `0001-divg-species-ztemp-pad.patch` | both (tested on master) | `Source/divg.f90` (`SPECIES_ADVECTION_PART_1_NEW`, 12 statements) | Assign all four `Z_TEMP` elements in the two wall loops (pad with `0._EB`, as mass.f90 does); removes a stale-scratch read under `FLUX_LIMITER='MP5'` | SUPERBEE and CHARM: output files byte-identical to the unpatched binary; MP5: changes (expected). Full V&V baseline comparison not run | GPU Legacy Mapper stream | proposed |
| 0002 | `0002-divg-enthalpy-ztemp-pad.patch` | both (tested on master) | `Source/divg.f90` (`ENTHALPY_ADVECTION_NEW`, 6 statements) | Same fix in the enthalpy wall loop | SUPERBEE, CHARM, MP5: byte-identical on the test case (MP5 effect not demonstrated). Full V&V baseline comparison not run | GPU Legacy Mapper stream | proposed |
| 0003 | `0003-mass-check-density-two-pass.patch` | both (tested on FireX 36975d765f; applies to master ce1f659cd4) | `Source/mass.f90` (`CHECK_MASS_DENSITY`, species loop, 35 lines added and 8 removed) | Split the species clipping loop into a per-cell pass (stores the seven amounts and a flag per cell) and a serial scatter pass in the original cell order; makes the first pass device-eligible | FDS run (propane, four OBSTs, OPEN vents): restart, smoke3d, slice, devc and hrr files byte-identical to the unpatched build, clipping exercised; routine harness, 60 meshes, byte-identical at -O0 -fcheck=all, -O2, -O3, -O2 -fopenmp; 8 mutants caught. Full V&V baseline comparison not run | GPU Legacy Mapper stream | proposed (sign-off request in `0003-mass-check-density-two-pass.signoff.md`, not yet sent) |

## Blocked-loop family reviews (decision D-051, item 1)
Each blocked-loop family is signed off by a domain lead before its rewrite becomes a patch.

| Family | Reviewer | Status |
|---|---|---|
| Zone sums with a pressure-zone index (`USUM(IPZ)`, `DSUM`, `PSUM`) (P1), `CONNECTED_ZONES` (P2), solid-cell DP correction L0394 (P3), `LOG_INTWC` (P4) | Pressure Solver Lead | signed off 2026-10-02: P1 accept under D-053 (measure the single-thread device add pass before making it the GPU default; privatise `IPZ`; only uncovered cells and owned faces under AMR; HVAC `U_NORMAL` final first), P2 accept, P3 accept with change (two-pass CSR gather in ascending wall index; `BOUNDARY_PROP1(WC%BC_INDEX)` vs `B1_INDEX` is a candidate upstream patch, debug assert until then), P4 accept |
| `DELTA_RHO_ZZ`/`DELTA_RHO` scatter (S1), wall nests with pointer scratch (S2), wall scatter into `U_DOT_DEL_RHO_*` (S3), `SETTLING_VELOCITY` (S4) | Species & Combustion Lead | signed off 2026-10-02: S1 accept (two-pass gather, same operation order, no FMA, interior targets only), S2 accept with change (uniqueness assertion only for writes that pass their flow guards; pad the fourth `Z_TEMP` element with 0 until patches 0001/0002 land, exclude MP5 divg from bit tests meanwhile), S3 accept (fixed per-cell wall key independent of decomposition; zero the result first), S4 accept (no other kernel uses `WORK7-9` between the three kernels). Details in `combustion/` notes of the lead |
| `CELL_COUNTER` (SP1), ghost mirror (SP2), DP accumulate (SP3), `WALL_BC` (SP4) | Solid Phase Lead | signed off 2026-10-02 (`solid/04-blocked-loop-signoff.md`): SP1 accept, SP2 change (the target-distinct assertion would fail at corners and thin obstructions; gather per target cell and take the highest wall index that passes the SOLID/EXTERIOR test), SP3 accept (exclude NULL, INTERPOLATED, OPEN walls), SP4 change (many-to-one sums are guarded upstream by `OMP CRITICAL`; per-wall scratch and ascending-wall-index sums per obstruction and per gas cell, no device atomics; `MASS_FLUX_VAR` statistically equal, not bitwise) |
| `CHECK_MASS_DENSITY` species clipping loop, two-pass split (patch 0003) | Species & Combustion Lead | sign-off request written (`0003-mass-check-density-two-pass.signoff.md`), not yet sent |
| Radiation loops | Radiation Lead | requested 2026-10-02 |
| `CHECK_STABILITY` (V1) and `CHECK_DIVERGENCE` (V2) reductions | V&V Lead | signed off 2026-10-02: V1 accept (merge per-box results by value, then larger global K,J,I; device `**ONTH` last-bit difference accepted if the location matches for distinct maxima; tie test without `pow`), V2 accept (run still stops on NaN through the existing check; minimum keeps first-wins, smaller global index; test for it), OpenMP merge not translated, accept |
