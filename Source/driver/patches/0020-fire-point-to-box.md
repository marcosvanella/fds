# Patch 0020: `fire.f90` combustion kernels find their box through `POINT_TO_BOX` (D-056)

## What it does
- The same macro as patch 0008, placed before `MODULE FIRE`. It applies only under `WITH_AMREX`:
  `#define POINT_TO_MESH POINT_TO_BOX`.
  - Every `CALL POINT_TO_MESH(NM)` of the file therefore reaches a fine-level box: COMBUSTION_LOAD_BALANCED, COMBUSTION_GENERAL_LOAD_BALANCED, the load-balance routines, COMBUSTION_BC and CONDENSATION_EVAPORATION.
  - Before this patch, a fine-mesh COMBUSTION call stopped with `POINT_TO_MESH: mesh N is not a level-0 mesh`.
- A new module function `FIRE_MESH(NM)` returns the mesh object: `MESHES(NM)` on level 0, or `FINE_LEVEL(L)%BOX(IB)` for a fine box.
  - The 12 direct `MESHES(NM)%` reads/writes (N_CUTCELL_MESH, NCHEM_ACTIVE_CELLS, NCHEM_ACTIVE_CC) now go through a local pointer `MB=>FIRE_MESH(NM)`, set right after the `POINT_TO_MESH(NM)` call.
  - `MB` is declared in 4 routines and assigned after 6 `POINT_TO_MESH` calls.
- On level 0, `FIRE_MESH(NM)` is `MESHES(NM)`, so results are unchanged.
- The companion driver change (same commit, `driver/fds_fine_level.f90`) has BUILD_FINE_BOX allocate `CHEM_ACTIVE_CELLS` the way init.f90 does. Without it, COMBUSTION_GENERAL_LOAD_BALANCED segfaults on a fine box at the first chemically active cell.

## Checks (GNU Release, 1 rank)
- **Level 0 is bit-identical.** simple_test (reacting, RADIATION off), `--levels 1`, 20 steps, binaries built before and after: the plotfiles (step 20 and final), `st_multi_level.csv`, `simple_test_steps.csv` and the s3d/smv files are byte-identical.
- **Fine combustion runs.** With a per-mesh FDS_K_COMB (Species' stage2b WIP, adapted to HEAD in a scratch copy and not committed), a two-level simple_test run (`--tag-vort auto --blocking 2`, 20 steps) runs combustion on the L1 boxes and completes all 20 steps.
  - There is no POINT_TO_MESH stop.
  - The known two-level health issue remains: fine div 0.87.
- **Pyramid** two-level, 20 steps, 4 ranks: PASS (L1 div 9.55e-12, mass 5.3e-13).
