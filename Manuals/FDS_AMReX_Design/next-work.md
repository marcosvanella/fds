# Next work per role (kept by the Chief Architect)

Owner rule: nobody sits idle. Each row says what the role is on, what is ready next, and what blocks it. The Chief Architect updates this file at each milestone and reports to the project coordinator when a role has no ready work. A role that finishes its list tells the Chief Architect.

Status as of the weekend work restart, from the latest reports.

Decisions since the last list: D-056 fine-level mesh objects by option B (condition met for 0005/0006; 0007-0009 still need oneAPI validation), D-057 2-D/singular pressure BC mapping owner, D-058 R4 tagging and regrid design, D-059 shared corner ghost cells accepted, D-060 normal-linear face prolongation, D-061 flux override ordering and W2 once per stage.

| Role | Ready work, in order | Blocked or waiting on |
|---|---|---|
| Role 1, Data Layout | (1) two-level end-to-end run of the periodic 2-D case (`ns2d_16` 1-to-2 refinement) with Role 3's override lists wired to the flux hooks; check composite mass and the max abs(div u - D); (2) hook the registry entry points to Role 3's regrid transfer; (3) FR-072 output plan with the Legacy Mapper's list of mesh-looping output routines (ADR-004); (4) first level-0 + level-1 timing run for the Integration Lead | oneAPI validation of 0007-0009 before physics use |
| Role 2, Pressure Backend | (1) M2: assembled-matrix HYPRE PCG+BoomerAMG backend and the MLMG backend behind the solver-agnostic interface; (2) rebuild entry point after regrid (D-058); (3) 2-D and singular-case mapping check (D-057) using Role 1's `notes/fft-thin-direction-check.md` | A-56 comparison from the Pressure Lead for the default |
| Role 3, Regrid and Transport | (1) R2b: flux-override adapter against Role 1's hooks (patch 0009) and a two-level conservation run; (2) R4 part 2: registry integration of the transfer operators (begin/end_regrid, fill_initial_level); (3) moving-blob test at the driver level; (4) GPU path for tagging kernels (K2, host loops proven) | Role 1 two-level run |
| AMReX Integration Lead | (1) review Role 1's flux-hooks design (`notes/flux-hooks-design.md`) and report; (2) stage-1 GPU spike on the test machine (standing approval), including the wall-seam GPU needs in `notes/wall-seam-design.md`; (3) apply the rule 7 macro form in `s4_omp.inc` and re-time | none |
| Pressure Solver Lead | (1) A-56 backend comparison (MLMG vs HYPRE, eps_H) plan and first run; (2) review the 2-D mapping (D-057); (3) single-thread device add measurement for zone sums (P1 note) | none |
| Legacy Mapper | (1) list of every output routine that loops over meshes for FR-072/ADR-004; (2) send the patch 0003 sign-off request to the Species Lead; (3) O3 (L1402) uniqueness-assertion review routing; (4) check that patch 0007-0009 touch no other code path that uses `MESHES(NM)` with a fine mesh number | none |
| Species and Combustion Lead | (1) review patch 0003; (2) combustion under refinement note for D-050/D-058 (species clip accounting, `Z_TEMP` pad of 0001/0002); (3) species tagging threshold advice for R4 | none |
| Solid Phase Lead | (1) fine-level solid phase plan (1-D wall conduction under refinement; OBST deferred); (2) review the generated SP2/SP4 kernels once the Wall Loops Engineer lands them | generator output for (2) |
| Radiation Lead | (1) blocked-loop sign-off for the radiation family (still open); (2) FR-062 sweep implementation notes for Phase 4; (3) radiation under regrid (what must be rebuilt) | none |
| V&V Lead | (1) Phase 3 acceptance case list for two-level runs (conservation, FR-016 baseline, blob) and thresholds; (2) FDS-only marking for the 23 pressure-code-0 inputs; (3) 2-D (`ns2d_16`) test for D-057; (4) Debug build for the patch-check tool; (5) multi-box merge tests for `CHECK_STABILITY`/`CHECK_DIVERGENCE` | Debug build |
| Spec and Program Lead | (1) carry D-056 to D-061 into the spec and IR/NFR text; (2) chase the Radiation and Species sign-offs; (3) update the roadmap for dates given by Roles 1 and 3 | none |
| GPU Generator Engineer | (1) coverage items for VELOCITY_FLUX, mass flux, predictor, corrector, DIV_PART_2; (2) P1/V1 summation-order statements | none |
| GPU Wall Loops Engineer | (1) wall loops with the Solid and Species conditions (SP2 per-target gather, SP4 scratch sums, S2 guarded assertion); (2) O3 (L1402) | none |
| GPU Mesh Data Loops Engineer | (1) `PATCH_VELOCITY_FLUX` (L1390); (2) remaining O2 edge and exchange loops (D-055 excludes `NOM>0` branches) | none |
| GNU Build Chief | (1) Release and Debug reference binaries current; (2) validate patches 0005-0009 with gfortran `-fcheck=all` and the Debug build; (3) runtime checks of the OFF build | none |
| Intel Build Chief | (1) validate patches 0007-0009 with 0005/0006 under oneAPI; (2) kernelcheck, decomposition check, csmag check under oneAPI; (3) same-compiler baseline option for the output gate | none |
