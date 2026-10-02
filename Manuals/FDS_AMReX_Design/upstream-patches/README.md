# Upstream patches (decision D-051)

Proposed changes to upstream FDS source (FireX or master) that come out of the AMR/GPU work, for the owner to review and commit upstream himself. Nobody else commits or pushes them upstream.

## Rules for a patch
1. One short-scope change per patch (one marker set, one loop cleanup, one routine). No mixed patches.
2. Delivered as a file in this folder: `git format-patch` output or a unified diff against FireX, applicable to the upstream source file as it is, without our `WITH_AMREX` code.
3. A short plain-language rationale in the patch description: what changes, why the GPU or AMR work needs it, what it does not change.
4. A behavior-unchanged check: the exact command or test and its result (for example the V&V baseline comparison against the reference binary, bitwise or within the stated tolerance), run on the patched upstream file.
5. A target line: **FireX**, **master**, or **both**. Master is merged into FireX periodically, so write the patch to apply cleanly to both and say which one it was tested against.
6. Name: `NNNN-short-name.patch` (four-digit sequence), and add a row to the index below in the same change.

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
| Zone sums with a pressure-zone index (`USUM(IPZ)`, `DSUM`, `PSUM`) | Pressure Solver Lead | requested 2026-10-02 |
| `DELTA_RHO_ZZ` scatter, species and combustion loops | Species & Combustion Lead | requested 2026-10-02 |
| Solid-phase counters and wall loops (`CELL_COUNTER` and similar) | Solid Phase Lead | requested 2026-10-02 |
| `CHECK_MASS_DENSITY` species clipping loop, two-pass split (patch 0003) | Species & Combustion Lead | sign-off request written (`0003-mass-check-density-two-pass.signoff.md`), not yet sent |
| Radiation loops | Radiation Lead | requested 2026-10-02 |
| `CHECK_STABILITY` reductions and other cross-cutting reductions | V&V Lead | requested 2026-10-02 |
