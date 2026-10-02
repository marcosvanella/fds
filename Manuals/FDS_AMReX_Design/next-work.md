# Next work per role (kept by the Chief Architect)

Owner request (2026-10-02): idle team members always get relevant work. Each row says whether the role is busy or free, what is ready for it, and what blocks it. The Chief Architect updates this file at each milestone and reports to the project coordinator when a role has no ready work. A role that finishes its list tells the Chief Architect, who either adds ready work or reports it.

Status as of 2026-10-02, from the latest reports (re-check with the role before relying on it).

| Role | State | Ready work, in order | Blocked or waiting on |
|---|---|---|---|
| Role 1, Data Layout | busy | (1) decision on fine-level mesh objects: option B ruled (D-056); (2) flux read-out and override hooks (design due 6 Oct); (3) W1/W2 wall-state seam; (4) kernel loops on level 1 once the fine-level mesh objects exist (about 9-14 Oct) | oneAPI validation of patch 0005 and 0006 |
| Role 2, Pressure Backend | busy | (1) M2 plan; (2) 2-D case boundary mapping check for the FFT backend (D-057); (3) HYPRE backend per ADR-002 v1.1 | A-56 comparison from the Pressure Lead |
| Role 3, Regrid and Transport | partly free | (1) R2 AmrCore subclass (waits on Role 1 items 1-5, now delivered); (2) fine-face and ghost-fill unit tests against the interface header; (3) overwrite averaging code (area-sum of fine fluxes) as a stand-alone tested function | flux hook design from Role 1 (6 Oct) for end-to-end use |
| Integration Lead | busy | WP1 remaining non-wall kernels on real fields; review of the flux hook; CUDA AMReX install with FFT and MPI (standing approval for test-machine runs) | WP3 generator coverage (Legacy Mapper and generator engineers) |
| Pressure Solver Lead | free after sign-offs | (1) A-56 backend comparison plan and run (MLMG vs HYPRE, eps_H); (2) review of the 2-D boundary mapping (D-057); (3) measure the single-thread device add for zone sums (P1 note) | none |
| Legacy Mapper | busy | amend `amrex/blocked-loop-families.md` with the sign-off conditions; patch 0003 sign-off request to the Species Lead; O3 (L1402) uniqueness assertion review; list of output routines that loop over meshes for FR-072 | none |
| Species and Combustion Lead | free after sign-off | (1) review patch 0003; (2) combustion-on-AMR spec follow-up for D-050 wording; (3) check MP5 divg `Z_TEMP` pad (patches 0001/0002) | none |
| Solid Phase Lead | free after sign-off | (1) fine-level solid phase plan (1-D wall conduction under refinement, OBST deferred); (2) SP4 `WALL_BC` scratch-sum review of the generated kernel | generator output |
| Radiation Lead | free | (1) blocked-loop sign-off for the radiation family (open); (2) FR-062 sweep implementation notes for Phase 4 | none |
| V&V Lead | busy | (1) `vv-runs/tools/patch_check.sh`; (2) Debug build refresh; (3) multi-box merge tests for `CHECK_STABILITY` and `CHECK_DIVERGENCE` (V1/V2 conditions); (4) 2-D case test (ns2d_16) for D-057 | Debug build |
| Spec and Program Lead | busy | spec v0.4.31 follow-through; chase the Radiation and Species sign-offs; carry D-056 and D-057 into IR/NFR text | none |
| GPU Generator Engineer | busy | P1/V1 summation-order statements; coverage items for VELOCITY_FLUX, mass flux, predictor, corrector, DIV_PART_2 | none |
| GPU Wall Loops Engineer | busy | wall loops with the Solid and Species conditions (SP2 per-target gather, SP4 scratch sums, S2 guarded assertion); O3 (L1402) | none |
| GPU Mesh Data Loops Engineer | busy | `PATCH_VELOCITY_FLUX` (L1390) first, then O2 edge and exchange loops that remain (D-055) | exchange-buffer layout for later items |
| GNU Build Chief | check | keep the Release and Debug reference binaries current for V&V (Debug rebuild is blocking `patch_check.sh`) | none |
| Intel Build Chief | ready | oneAPI validation of patches 0005 and 0006 on the test machine (standing approval for runs) | none; this unblocks D-056 |
