# fds_csmag32_periodic: FDS-derived frozen pressure case (M4)

One real FDS right-hand side and FDS's own solved H, for checking the backends against FDS itself. Single mesh, 32 x 32 x 32 cells,
uniform cells (dx = dy = dz = 0.56549/32), all six faces periodic, no obstructions.

## Provenance
| Item | Value |
|---|---|
| FDS | FireX 36975d765f (`FDS-6.11.1-1244-g36975d765f`), GNU 14.2 / Open MPI 5.0.7, Release, OpenMP on, HYPRE and SUNDIALS off (not used by the case), 1 rank, 1 thread |
| Patch | upstream patch 0009 (`docs/upstream-patches/0009-pressure-rhs-dump.patch`, rationale in `0009-pressure-rhs-dump.signoff.md`), enabled with `FDS_PDUMP_STEPS=1,3,10` |
| Input | `fds_csmag32.fds` (here): the V&V case `csmag_32` with `FISHPAK_BC=0,0,0` (FFT Poisson solver, periodic in x, y, z), `CHID='fds_csmag32'` and `SIG_FIGS=17` added to `&DUMP`. It also needs `cbc32_uvw.csv` (870916 bytes, not committed): `vv-runs/baseline/gnu_ompi_firex-36975d7/csmag_32__fishpak_bc000/cbc32_uvw.csv` |
| Step, stage | step (ICYC) 3, t = 4.8932e-2 s, dt = 1.7363e-2 s; stage P (predictor, H) and stage C (corrector, HS) |
| Solver in FDS | `PRES_FLAG` 0 (FFT, FISHPAK `H3CZSS`), `IPS` 0, `LBC=MBC=NBC=0` (periodic) |
| Control | the dump-enabled run is bitwise identical to the unpatched build in all 40 comparable output files (see the patch sign-off) |

## Files (`fds-pdump-1` format, written by the FDS patch)
`fds_csmag32_pdump_n000003_<P|C>_m001_` followed by `meta.txt` (key = values: extents, cells, face types, PRES_FLAG, step, stage, T, DT, face coordinates),
`rhs.bin` (FDS `PRHS` before the solve), `phi.bin` (FDS `H`/`HS` after the solve), `rho.bin` (`RHO`/`RHOS`), `kres.bin` (`KRES`).
All `.bin`: little-endian float64, valid cells only, 32768 values, I fastest. About 2.1 MB in total. (The dump also writes `zone.bin`
and `bc.bin`; they are trivial for this case, one pressure zone and no boundary data, and are not kept.)

## Test
`pb_fds_frozen` (`harness/m4_modes.cpp`, registered in `tests/m4_fds_frozen.cmake` as `pb_fds_frozen_csmag32_periodic_{P,C}_np{1,2}`):
rebuilds the problem from `meta.txt`, solves the stored `rhs.bin` with FFT, MLMG and HYPRE, removes the volume-weighted (here arithmetic) mean of
every H (the additive constant is a gauge), and compares with `phi.bin` and between backends. Pass: relative L2 difference <= eps_H = 1e-8.
FDS's own H is also put into the interface's 7-point operator (`RESIDUAL` line) to show that FDS and the backends solve the same discrete operator.

Result (step 3, corrector, 1 rank; the predictor and 2 ranks agree at the same level):

| comparison | relative L2 | relative max |
|---|---|---|
| FFT vs FDS | 3.5e-15 | 1.2e-15 |
| MLMG vs FDS | 2.3e-12 | 7.8e-13 |
| HYPRE vs FDS | 1.6e-13 | 2.0e-13 |
| MLMG vs FFT | 2.3e-12 | 7.8e-13 |
| HYPRE vs FFT | 1.6e-13 | 2.0e-13 |
| FDS H in the interface operator: max abs residual / max abs (rhs - mean) | 6.0e-16 | |

## Regenerate
Apply patch 0009 to a scratch copy of FireX 36975d765f, build, then in a directory with `fds_csmag32.fds` and `cbc32_uvw.csv`:
`FDS_PDUMP_STEPS=3 mpirun -np 1 fds fds_csmag32.fds`.

## Multi-mesh
Not committed. The patch writes one file set per mesh (`_m<NM>_`), and a multi-mesh periodic case such as `shunn3_4mesh_32` with
`FDS_PDUMP_STEPS=n` gives four sets (checked: the run completes and its outputs are unchanged), but `pb_fds_frozen` stops on `nmeshes > 1`.
FDS's multi-mesh FFT (default) and GLMAT/ULMAT differ from one global solve by the accepted coupling difference, so a multi-mesh comparison
needs the gathered global right-hand side and an agreed tolerance first.
