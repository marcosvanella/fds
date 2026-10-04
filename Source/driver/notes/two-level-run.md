# Two-level run (S14)

## What runs
`fds_amr ns2d_16_l0.fds --two-level-run --patch 4 11 4 11 --steps N` (ns2d_16, level 0 = 16 x 1 x 16, level 1 = coarse cells 4..11, ratio (2,1,2), one fine box 16 x 1 x 16 cells).
Per step: viscosity; predictor density (`stage_density_all`: ADV read-out of every level, Role 3's override lists, finest-first update); exchange 1/4; velocity flux; wall BC;
DIVERGENCE_PART_1 with the DIF overwrite; pressure scheme (BAROCLINIC, MATCH_VELOCITY_FLUX on level 0, NO_FLUX, FV average-down, level-0 PRHS by FDS, fine PRHS in C++, one composite solve,
H gradient at the fine C/F faces); velocity predictor; the same for the corrector; face-velocity average-down after each velocity update.

## Findings that fixed the first runs
* NaN in D/DS (the "D goes NaN after the first prime pass"): not the D_PBAR_DT zone term (ns2d has no pressure zone). A new fine level carries only the transferred state; the stage arrays
  RSUM, MU, KRES, D, DS, H, HS were zero, DIVERGENCE_PART_1 divides by RSUM, and the average-down of the cell scalars wrote the zero into the covered coarse cells too (64 cells of level 0).
  FDS skips the density/EOS update in cycle 1 (ICYC<=1), so nothing rebuilt them. Fix: `RegistryTransfer::derive` in `TwoLevelRun.cpp` injects the parent value into the new level for
  D, DS, RSUM, MU, KRES, H, HS. (Role 3's direct `bind_level` tests must do the same; an EOS rebuild of RSUM/TMP from the fine ZZ is the better derive for a non-uniform state.)
* FVX/FVY/FVZ are cell-shaped arrays holding the flux at the UPPER face of the cell (FDS index). The fine PRHS used the lower-face convention: the composite right-hand side had a mean of 1.6
  and the predictor divergence error was 1.06; with (FVX(I-1)-FVX(I))/dx the error is 1.6e-2. The FV average-down has to use the same convention (`composite_average_down_flux`: nodal face n = FV(n-1)).
* TWO_D: the PRHS has no y term when the level has one cell in y.

## Design constraint for GPU work
K2 kernel launches are blocking and are not issued on AMReX streams. Per-box concurrency therefore means one host thread per box; the driver has to synchronise the AMReX stream
(`amrex::Gpu::streamSynchronize()`) before any K2 launch that consumes a result written on a stream (for example the PRHS, the average-down or the ghost fills), and again before AMReX
reads what a K2 launch wrote.

## Numbers (ns2d_16_l0, patch coarse cells 4..11, ratio (2,1,2), 1 rank, flux overwrite on)
* Level 0 periodic ghost corruption (the last open cause of the composite PRHS mean): `after_exchange(6)` run for a fine box (fill_omesh, then VELOCITY_BC) overwrote the periodic z ghost rows of the
  LEVEL-0 U (and the x ghost columns of W). ROOT CAUSE (S14.4, `notes/fine-velocity-bc.md`): the `BcStep` of a fine level was built from the registry `Level`, whose `fds_mesh_offset` is 0, so it called
  the FDS routines with mesh number box index + 1 instead of NM0 + box index + 1: they ran on level-0 mesh objects. Fixed by `BcStep::set_mesh_offset(nm0)` in `bind_level`; save_uvw and VELOCITY_BC (fill_omesh is level 0 only)
  run on a fine level again (default `FDSTL_SKIPAFT=0`), results identical to the skip. Reproducer: `tests/run_stage_boundary_fine_check.sh`.
* Result after the fixes: 40 steps (t = 0.42): composite mass 3.947843523 -> same, relative change 1.3e-14; rho*Z the same; max|div u - D| over the uncovered cells 1.8e-12 (level 0 6e-14,
  level 1 1.8e-12), removed mean of the Poisson right-hand side 4e-13 (relative to rms 4e-17), 100 composite solves with MLMG.
* Debug switches: `FDSTL_GHOSTDBG=1` (RANGE/DIVERR/RHSSUM/WRAP lines; 2: field ranges; 3: worst cell and Poisson residual), `FDSTL_TLDIAG=1` (PROBE lines of the level-0 ghost rows).

## Input converter in the driver (D-076)
`main.cpp` calls `prepare_amr_input` (Role 3's `InputConverter`, patch `regrid_transport/notes/driver-patch-main-input-converter.patch`): a multi-level `&MESH` input with an `&AMR` line is converted to a level-0 input (finer meshes removed, one cover mesh per hole, written as `<stem>_amr_level0.fds` in the working directory by rank 0); a level-0-only input is passed through unchanged and no file is written. `assemble_level0`, `print_pressure_bc` and the end-of-run `fds_setup` use the converted name; `build_level0` sees level-0 meshes only. The converted hierarchy and parameters (`conv.hierarchy`, `conv.params`) are not yet consumed by the step loop. Checked with `tests/cases/ns2d_13mesh_amr.fds` (13 meshes: 12 level-0 meshes, 1 finer mesh removed, 1 cover mesh added; 13 boxes, one cell size): `--run` of 10 steps gives the same T and DT as the one-mesh `ns2d_16_l0` and the same UVEL device record; the single-level runs of `run_level_bind_check.sh` stay bitwise equal to the reference executable.

## Controls, multi-box fine level, 4 ranks (S14.4; 40 steps, t = 0.42, fine patch coarse cells 4..11, ratio (2,1,2), composite mass 3.947843523)
| Run | fine boxes | ranks | mass change (relative) | max abs(div u - D) |
|---|---|---|---|---|
| overwrite ON | 1 | 1 | 1.3e-14 | 1.84e-12 |
| `--no-overwrite` | 1 | 1 | 1.3e-14 | 1.84e-12 |
| `--maxsize 8` | 4 | 1 | 1.3e-14 | 1.84e-12 |
| `--maxsize 4` | 16 | 1 | 1.3e-14 | 1.84e-12 |
| `ns2d_16_4m` (4 level-0 meshes, one per rank), default boxes | 4 | 4 | -1.2e-15 | 1.83e-12 |
| `ns2d_16_4m`, `--maxsize 4` | 16 | 4 | -1.2e-15 | 1.84e-12 |
* The overwrite-off control is not discriminating on `ns2d_16`: constant density and one species give a constant advective flux, so the interface flux overwrite changes nothing. The discriminating controls: `driver_unit_tests` TWOLEVEL (overwrite ON composite rho change 2.2e-15, OFF 7.2e-4 over 45 steps) and Role 3's `run_e2e_driver.sh` E3b (ON 6e-16, OFF 4.35e-5).
* The multi-box and 4-rank runs need no same-level seam FV averaging on this case (fine boxes of one level share the faces through the AMReX exchange); a case with a flux through a seam of two fine boxes that are not exchanged is not covered.
* Before S14.4 the 4-rank runs with 16 fine boxes crashed (segfault in `fds_p_save_uvw`, then in `fds_g_fill_om`): the BcStep mesh-number offset, `notes/fine-velocity-bc.md`.
* Role 3's regrid tests against the committed `RegistryTransfer` (standalone `regrid_transport` build): all 14 ctest entries pass (celltransfer, regrid_core 1/4 ranks and rank check, facetransfer, species_avgdown, blob_registry 1/4 ranks and rank check, moving_blob_p3 1/4 ranks); `test_regrid_core` worst composite change per regrid 1.87e-16, clips 0, 6 of 6 regrids changed the grids, hierarchy hash 9678b9b0b9139a7b, data hash 1cbb18482dc6ee2d.

## Variable-density two-level case (discriminating overwrite control)

The constant-density runs above cannot tell the interface flux overwrite (D-050, vv/test-plan.md 5.12.6) from its absence: with uniform density the
coarse and fine mass fluxes agree and the drift stays at round-off. `tests/cases/blob_16_l0.fds` (one mesh) and `blob_16_4m.fds` (four 8x1x8 meshes,
one per rank) are the 16x1x16 doubly periodic unit square with a uniform wind (U0=1, W0=0.5) and a 28/44 species blob (`&INIT`, TRACER, 0.2 to 0.4 in x,
0.4 to 0.6 in z), so density and the species mass fractions vary across the coarse/fine interface. The fine patch is coarse cells 4..11 in x and z (at np 4 it spans all four meshes). `tests/run_two_level_density_check.sh <build> [work] [steps1=100] [steps4=40]` runs it three ways at
1 and 4 ranks (ON; OFF = `--no-overwrite`; REF = a fine level that covers the whole domain) and checks ON drift <= 1e-12, OFF drift >= 1e-5, and the
last-step max|div u - D| of the two-level run within 3x of REF. `run_driver_tests.sh` runs it with `DRIVER_DENSITY_CHECK=1` and always runs
`run_stage_boundary_fine_check.sh`.

| run | steps | mass change | rho*Z1 change | rho*Z2 change |
|---|---|---|---|---|
| np 1 ON | 100 | -9.1e-15 | -1.1e-14 | -2.2e-15 |
| np 1 OFF | 100 | 2.5e-4 (worst 1.4e-3) | 5.7e-4 | -3.8e-3 |
| np 1 REF | 100 | -1.7e-14 | -1.8e-14 | -2.8e-15 |
| np 4 ON | 40 | 1.6e-15 | 2.0e-15 | -6.1e-16 |
| np 4 OFF | 40 | -1.4e-3 | 2.8e-3 | -5.6e-2 |
| np 4 REF | 40 | 3.5e-15 | 3.8e-15 | -4.0e-16 |

max|div u - D| over uncovered cells: the worst-over-steps value is 1.436 in all six runs. It is the first step (the baroclinic term is lagged when
`PRESSURE_ITERATIONS` is not iterated; the single-level driver run of `blob_32` shows the same large error, so it is not a two-level effect) and does not depend on the
overwrite. The last step is the comparison that matters: np 1 0.633 (two-level) vs 0.623 (REF), np 4 0.664 vs 0.709. The overwrite changes conservation
(round-off versus 1e-3 to 6e-2) and not the divergence error, so the two numbers are reported together.
