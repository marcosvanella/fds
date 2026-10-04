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
* Level 0 periodic ghost corruption (the last open cause of the composite PRHS mean): `after_exchange(6)` run for a fine box (fill_omesh, then VELOCITY_BC of the fine mesh) overwrote the
  periodic z ghost rows of the LEVEL-0 U (and the x ghost columns of W); the first predictor FVX/FVZ at the periodic faces then differed from the single-level run (sums 91.8 vs 10.6 at the two faces of a
  periodic pair) and the sum of PRHS did not telescope (6.38 = removed mean 1.6 x domain volume). A level > 0 has no external wall cell, so `BcStep::after_exchange` no longer runs fill_omesh or
  VELOCITY_BC on it (`FDSTL_SKIPAFT=0` restores them). The cause inside the fine-mesh VELOCITY_BC is not understood; to be checked when fine walls are built.
* Result after the fixes: 40 steps (t = 0.42): composite mass 3.947843523 -> same, relative change 1.3e-14; rho*Z the same; max|div u - D| over the uncovered cells 1.8e-12 (level 0 6e-14,
  level 1 1.8e-12), removed mean of the Poisson right-hand side 4e-13 (relative to rms 4e-17), 100 composite solves with MLMG.
* Debug switches: `FDSTL_GHOSTDBG=1` (RANGE/DIVERR/RHSSUM/WRAP lines; 2: field ranges; 3: worst cell and Poisson residual), `FDSTL_TLDIAG=1` (PROBE lines of the level-0 ghost rows).
