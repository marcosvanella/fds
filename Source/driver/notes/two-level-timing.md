# Two-level driver timing (r6 item 4)

For the Integration Lead and for `docs/amrex/stage1-gpu-spike-plan.md` section 10 ("Measured against estimated"). **run** = measured here; **ESTIMATE** = a model, basis given.
CPU only (the box: gfortran 14.2 `-O2`, Open MPI 5.0.7, AMReX 26.09). The box has **one shared core** (load average 15 to 43 during the runs), so the numbers are wall
seconds under load: best of 3 runs per configuration, 40 steps each. The 1-rank rows can be compared with each other. The **4-rank rows put 4 ranks on one core and are only
indicative**; their spread is large (for example ns2d16 np4 REF 235 ms per step against ns2d32 np4 REF 77 ms).

## What was timed
Cases: `ns2d_16_l0` (periodic analytic 2-D NS, constant density) and `blob_16_l0` (variable-density tracer blob, two species), generated at N x 1 x N coarse cells, N = 16, 32, 64,
as one mesh (1 rank) or 2 x 2 meshes (4 ranks). Per case and rank count:
- **single**: the single-level CPU driver on the coarse case (`fds_amr <case> --run`); pressure = level-0 `pressure_backend` FFT.
- **two**: `fds_amr <case> --two-level-run --patch N/4 3N/4-1 N/4 3N/4-1` (ratio 2: a fine patch over the middle half, N x N fine cells = 25 % of the domain at fine
  resolution); composite pressure = `pressure_backend` composite MLMG, flux overwrite ON, cf ghost hooks.
- **ref**: the same two-level run with the patch = the whole domain (`--patch 0 N-1 0 N-1`).
- **fine**: the single-level CPU driver on the uniform fine case (2N x 1 x 2N), the true uniform fine run.
All two-level and REF runs printed `TWO-LEVEL RUN PASS` (composite mass change at round-off).

## Instrumentation (the cheapest possible)
AMReX is built without TinyProfiler (`AMReX_TINY_PROFILE` OFF), so no `BL_PROFILE`. The driver's existing accumulators were reused: `t_loop`, `t_pressure`,
`t_poisson_solve`, `t_fds_outputs` and the `g_prof` slots 0 to 3 (FillBoundary, `fill_omesh`, `match_periodic_faces`, FDS bc routines). Added: three `g_prof` slots
(`Fields.H`; one `ProfScope` = two `amrex::second()` calls per scope, a few scopes per step) and `TimeLoop::write_perf()` so `--two-level-run` writes the same
`<chid>_driver_perf.csv` and, with `FDSTL_PROFILE=1`, `<chid>_driver_profile.txt`:
- slot 4 **avg-down**: `composite_average_down_faces` (end of predictor and corrector) and `composite_average_down_flux` (in the pressure scheme);
- slot 5 **cf ghost fill**: the cf scalar ghost hook (`BcStep::exchange`) and `composite_fill_fine_face_ghosts` (codes 3 and 6);
- slot 6 **flux overwrite**: the ADV read-out (`compute_stage_fluxes`) and the DIF read-out (`readout_dif`) of every level plus the override lists. `readout_dif` *is* the
  DIVERGENCE_PART_1 pass of the level, so this slot includes one regular DIVERGENCE_PART_1 per level; the re-runs with the lists (`apply_flux_divergence`,
  `run_divergence_part1`) are in "host stages".
Not measured: host/device copy counts and bytes (CPU build; no device copies) and launches (no kernels in this build). Hand-offs are therefore zero by construction here.

## Per-stage breakdown (run, ms per step, best of 3)
"pressure scheme" = `t_pressure` (the whole PRESSURE_ITERATION_SCHEME incl. the solve); "of which solve" = `t_poisson_solve` (FFT for single and fine, composite MLMG for two
and ref). "host stages" = step - pressure scheme - outputs - avg-down - cf ghost fill - flux overwrite (the FDS stages outside the pressure scheme, same-level ghost exchange
included). FillBoundary, `fill_omesh` and the FDS bc routines run inside and outside the pressure scheme and are listed, not subtracted.

| run | np | step | pressure scheme | of which solve | avg-down | cf ghost fill | flux overwrite | FillBoundary | fill_omesh | FDS bc routines | outputs | host stages |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ns2d16_np1_single | 1 | 3.71 | 0.99 | 0.40 | 0.000 | 0.000 | 0.000 | 0.027 | 0.728 | 0.151 | 0.209 | 2.51 |
| ns2d16_np1_two | 1 | 14.65 | 9.61 | 9.34 | 0.072 | 0.852 | 1.331 | 0.048 | 0.436 | 0.199 | 0.490 | 2.29 |
| ns2d16_np1_ref | 1 | 29.64 | 19.79 | 19.37 | 0.172 | 2.416 | 2.515 | 0.122 | 0.561 | 0.236 | 0.096 | 4.65 |
| ns2d16_np1_fine | 1 | 7.89 | 1.99 | 0.67 | 0.000 | 0.000 | 0.000 | 0.038 | 1.861 | 0.605 | 0.132 | 5.77 |
| ns2d16_np4_single | 4 | 22.40 | 8.47 | 1.43 | 0.000 | 0.000 | 0.000 | 1.313 | 12.403 | 0.165 | 0.764 | 13.17 |
| ns2d16_np4_two | 4 | 57.55 | 41.45 | 38.84 | 0.107 | 2.601 | 2.104 | 1.669 | 7.452 | 0.126 | 0.370 | 10.92 |
| ns2d16_np4_ref | 4 | 235.11 | 153.77 | 140.88 | 1.586 | 12.659 | 5.644 | 16.230 | 44.015 | 0.375 | 1.250 | 60.20 |
| ns2d16_np4_fine | 4 | 109.35 | 52.10 | 3.08 | 0.000 | 0.000 | 0.000 | 7.914 | 109.178 | 0.477 | 1.204 | 56.05 |
| ns2d32_np1_single | 1 | 14.85 | 3.99 | 1.32 | 0.000 | 0.000 | 0.000 | 0.136 | 3.335 | 0.993 | 0.411 | 10.44 |
| ns2d32_np1_two | 1 | 74.33 | 44.18 | 42.72 | 0.190 | 4.705 | 8.590 | 0.261 | 2.626 | 1.601 | 0.263 | 16.40 |
| ns2d32_np1_ref | 1 | 128.76 | 75.70 | 73.88 | 0.632 | 10.764 | 14.589 | 0.851 | 2.103 | 1.472 | 0.333 | 26.75 |
| ns2d32_np1_fine | 1 | 34.91 | 7.40 | 1.56 | 0.000 | 0.000 | 0.000 | 0.063 | 8.087 | 3.091 | 0.319 | 27.20 |
| ns2d32_np4_single | 4 | 26.49 | 14.64 | 3.02 | 0.000 | 0.000 | 0.000 | 0.943 | 17.840 | 0.423 | 0.408 | 11.44 |
| ns2d32_np4_two | 4 | 53.86 | 36.84 | 33.05 | 0.463 | 2.383 | 4.066 | 0.881 | 9.639 | 0.453 | 0.327 | 9.78 |
| ns2d32_np4_ref | 4 | 76.99 | 50.35 | 46.08 | 0.646 | 5.519 | 6.266 | 2.452 | 9.068 | 0.386 | 0.376 | 13.84 |
| ns2d32_np4_fine | 4 | 67.87 | 32.51 | 3.26 | 0.000 | 0.000 | 0.000 | 1.369 | 48.222 | 1.551 | 0.293 | 35.07 |
| ns2d64_np1_single | 1 | 29.94 | 6.66 | 1.56 | 0.000 | 0.000 | 0.000 | 0.066 | 7.430 | 2.424 | 0.189 | 23.09 |
| ns2d64_np1_two | 1 | 189.53 | 101.60 | 97.24 | 0.653 | 14.399 | 23.836 | 0.674 | 6.609 | 4.592 | 0.332 | 48.71 |
| ns2d64_np1_ref | 1 | 451.71 | 249.80 | 241.93 | 2.515 | 40.855 | 49.597 | 2.856 | 7.381 | 7.179 | 0.452 | 108.49 |
| ns2d64_np1_fine | 1 | 185.36 | 33.24 | 6.81 | 0.000 | 0.000 | 0.000 | 0.161 | 33.835 | 16.573 | 1.029 | 151.09 |
| ns2d64_np4_single | 4 | 57.57 | 25.83 | 2.57 | 0.000 | 0.000 | 0.000 | 1.365 | 39.477 | 1.239 | 0.689 | 31.06 |
| ns2d64_np4_two | 4 | 315.34 | 199.32 | 167.98 | 1.170 | 16.795 | 18.444 | 6.692 | 72.368 | 1.807 | 0.454 | 79.15 |
| ns2d64_np4_ref | 4 | 259.38 | 158.66 | 147.37 | 4.214 | 18.628 | 23.378 | 7.995 | 28.306 | 1.993 | 0.653 | 53.85 |
| ns2d64_np4_fine | 4 | 179.96 | 75.88 | 8.59 | 0.000 | 0.000 | 0.000 | 2.029 | 114.800 | 5.043 | 1.269 | 102.80 |
| blob16_np1_single | 1 | 3.22 | 0.69 | 0.29 | 0.000 | 0.000 | 0.000 | 0.028 | 0.557 | 0.155 | 0.284 | 2.24 |
| blob16_np1_two | 1 | 18.63 | 11.44 | 11.13 | 0.079 | 1.089 | 2.127 | 0.159 | 0.500 | 0.207 | 0.300 | 3.60 |
| blob16_np1_ref | 1 | 31.53 | 18.18 | 17.74 | 0.181 | 2.322 | 4.250 | 0.132 | 0.512 | 0.240 | 0.058 | 6.54 |
| blob16_np1_fine | 1 | 10.06 | 2.18 | 0.65 | 0.000 | 0.000 | 0.000 | 0.036 | 2.069 | 0.587 | 0.069 | 7.81 |
| blob16_np4_single | 4 | 5.50 | 2.66 | 0.54 | 0.000 | 0.000 | 0.000 | 0.145 | 3.829 | 0.054 | 0.045 | 2.80 |
| blob16_np4_two | 4 | 22.26 | 14.98 | 13.51 | 0.306 | 1.284 | 1.863 | 0.676 | 3.080 | 0.179 | 0.064 | 3.76 |
| blob16_np4_ref | 4 | 22.32 | 15.07 | 13.88 | 0.254 | 1.467 | 1.808 | 0.667 | 2.713 | 0.080 | 0.060 | 3.66 |
| blob16_np4_fine | 4 | 13.67 | 6.90 | 0.71 | 0.000 | 0.000 | 0.000 | 0.225 | 9.811 | 0.277 | 0.158 | 6.61 |
| blob32_np1_single | 1 | 10.83 | 2.14 | 0.67 | 0.000 | 0.000 | 0.000 | 0.037 | 2.304 | 0.590 | 0.062 | 8.63 |
| blob32_np1_two | 1 | 55.50 | 29.35 | 28.26 | 0.191 | 4.156 | 8.726 | 0.139 | 1.795 | 0.876 | 0.104 | 12.98 |
| blob32_np1_ref | 1 | 145.14 | 77.43 | 75.24 | 0.784 | 14.866 | 19.447 | 0.647 | 2.439 | 1.305 | 0.203 | 32.41 |
| blob32_np1_fine | 1 | 40.13 | 7.00 | 1.47 | 0.000 | 0.000 | 0.000 | 0.063 | 8.058 | 2.564 | 0.225 | 32.90 |
| blob32_np4_single | 4 | 11.42 | 4.98 | 0.71 | 0.000 | 0.000 | 0.000 | 0.314 | 7.568 | 0.175 | 0.060 | 6.38 |
| blob32_np4_two | 4 | 37.73 | 23.51 | 21.20 | 0.184 | 2.133 | 3.378 | 0.728 | 6.404 | 0.308 | 0.096 | 8.43 |
| blob32_np4_ref | 4 | 54.59 | 32.73 | 30.09 | 0.257 | 4.060 | 6.404 | 0.741 | 6.124 | 0.284 | 0.080 | 11.07 |
| blob32_np4_fine | 4 | 27.97 | 10.89 | 1.01 | 0.000 | 0.000 | 0.000 | 0.200 | 17.636 | 0.669 | 0.100 | 16.98 |
| blob64_np1_single | 1 | 37.72 | 6.81 | 1.48 | 0.000 | 0.000 | 0.000 | 0.058 | 7.860 | 2.367 | 0.201 | 30.71 |
| blob64_np1_two | 1 | 207.99 | 98.35 | 94.08 | 0.640 | 13.646 | 34.236 | 0.643 | 6.526 | 4.139 | 0.322 | 60.80 |
| blob64_np1_ref | 1 | 557.44 | 294.01 | 285.51 | 2.541 | 39.865 | 75.868 | 3.137 | 8.369 | 7.025 | 0.482 | 144.68 |
| blob64_np1_fine | 1 | 215.43 | 33.29 | 5.66 | 0.000 | 0.000 | 0.000 | 0.159 | 34.863 | 14.760 | 1.001 | 181.14 |
| blob64_np4_single | 4 | 36.62 | 14.33 | 1.73 | 0.000 | 0.000 | 0.000 | 0.562 | 22.295 | 0.762 | 0.118 | 22.17 |
| blob64_np4_two | 4 | 137.75 | 79.44 | 69.81 | 1.065 | 7.375 | 16.424 | 1.916 | 20.765 | 1.435 | 0.177 | 33.27 |
| blob64_np4_ref | 4 | 208.42 | 120.52 | 114.24 | 1.237 | 13.532 | 24.928 | 4.280 | 15.284 | 1.604 | 0.364 | 47.84 |
| blob64_np4_fine | 4 | 99.16 | 33.64 | 2.89 | 0.000 | 0.000 | 0.000 | 0.724 | 55.599 | 3.333 | 0.345 | 65.18 |

## Two-level against single-level and REF (ms per step)
Model (**ESTIMATE**): two-level = single coarse + 0.25 x uniform fine (the fine patch is 25 % of the uniform fine cells; one dt for both levels; no cost for the AMR machinery).
"AMR overhead per fine-patch cell" = (two - single) / fine-patch cells - uniform fine cost per cell, per step. Note: two, ref and fine all step with the fine dt (one dt over
levels, no subcycling), single with the coarse dt, so per simulated second the single-level coarse run is a further 2x cheaper.

| case | np | coarse cells | fine patch cells | single coarse | two-level run | model ESTIMATE | run / model | REF (patch = domain) | uniform fine single-level | AMR overhead per fine-patch cell | uniform fine cost per cell | two-level / uniform fine |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ns2d16 | 1 | 256 | 256 | 3.71 | 14.65 | 5.68 | 2.6x | 29.64 | 7.89 | 35.01 us | 7.71 us | 1.86x |
| ns2d16 | 4 | 256 | 256 | 22.40 | 57.55 | 49.74 | 1.2x | 235.11 | 109.35 | 30.49 us | 106.79 us | 0.53x |
| ns2d32 | 1 | 1024 | 1024 | 14.85 | 74.33 | 23.58 | 3.2x | 128.76 | 34.91 | 49.56 us | 8.52 us | 2.13x |
| ns2d32 | 4 | 1024 | 1024 | 26.49 | 53.86 | 43.46 | 1.2x | 76.99 | 67.87 | 10.15 us | 16.57 us | 0.79x |
| ns2d64 | 1 | 4096 | 4096 | 29.94 | 189.53 | 76.29 | 2.5x | 451.71 | 185.36 | 27.65 us | 11.31 us | 1.02x |
| ns2d64 | 4 | 4096 | 4096 | 57.57 | 315.34 | 102.56 | 3.1x | 259.38 | 179.96 | 51.95 us | 10.98 us | 1.75x |
| blob16 | 1 | 256 | 256 | 3.22 | 18.63 | 5.74 | 3.2x | 31.53 | 10.06 | 50.38 us | 9.82 us | 1.85x |
| blob16 | 4 | 256 | 256 | 5.50 | 22.26 | 8.92 | 2.5x | 22.32 | 13.67 | 52.09 us | 13.35 us | 1.63x |
| blob32 | 1 | 1024 | 1024 | 10.83 | 55.50 | 20.86 | 2.7x | 145.14 | 40.13 | 33.82 us | 9.80 us | 1.38x |
| blob32 | 4 | 1024 | 1024 | 11.42 | 37.73 | 18.41 | 2.0x | 54.59 | 27.97 | 18.87 us | 6.83 us | 1.35x |
| blob64 | 1 | 4096 | 4096 | 37.72 | 207.99 | 91.58 | 2.3x | 557.44 | 215.43 | 28.42 us | 13.15 us | 0.97x |
| blob64 | 4 | 4096 | 4096 | 36.62 | 137.75 | 61.41 | 2.2x | 208.42 | 99.16 | 18.64 us | 6.05 us | 1.39x |

## Reading (1 rank)
- **The composite MLMG solve is the largest item**: 45 to 64 % of the two-level step (ns2d64 97 of 190 ms, blob64 94 of 208 ms) against 0.3 to 7 ms for the FFT of the
  single-level runs (single coarse and uniform fine). This is where the two-level time goes; the remedies (an FFT-type composite solve, a cheaper MLMG set-up) are Role 2's.
- **Flux overwrite slot** 9 to 16 % (blob64 34 ms; includes the regular DIVERGENCE_PART_1 pass), **cf ghost fill** 6 to 8 % (blob64 13.6 ms), **avg-down** at most 0.5 %.
- **FDS host stages** of the two-level step: 61 ms (blob64) against 31 ms single coarse and 181 ms uniform fine: the fine patch's host work scales with its cells, as expected.
- **Two-level vs model**: run 2.3x to 3.2x the ESTIMATE; at N = 64 the two-level step (25 % of the domain refined) costs as much as the uniform fine run (0.97x to 1.02x).
- **REF** (patch = domain) is 2.4x to 3.8x the uniform fine single-level run: the composite path (MLMG + overwrite + cf fill) costs 2.4 to 3.8 times a plain fine step.
- **AMR overhead per fine-patch cell**: 28 us (blob64) and 28 us (ns2d64) per step, against 11 to 13 us per cell for a uniform fine step.

## Row for section 10 of the spike plan (Integration Lead to fold in; Role 1 does not edit docs/amrex)
| two-level CPU step, 1 rank | blob 64^2 coarse + 64^2 fine patch: 208 ms (composite MLMG 94, flux overwrite 34, cf ghost fill 14, avg-down 0.6, host stages 61) vs single coarse 38 ms, uniform fine 215 ms; model 92 ms | run (model ESTIMATE) | box, 1 shared core, best of 3, `notes/two-level-timing.md` |

## Commands
```
source Source/driver/tests/env.sh
Source/driver/tests/run_two_level_timing.sh <build-dir> <work-dir> "16 32 64" 40 3     # FDSTL_PROFILE=1 is set by the script; writes <work-dir>/timing.csv
python3 Source/driver/tests/two_level_timing_table.py <work-dir>/timing.csv             # the two tables above
# one run by hand, e.g. the two-level blob at N = 64 (case generated by the script as blob_64_1m.fds):
FDSTL_PROFILE=1 mpirun -np 1 <build-dir>/fds_amr blob_64_1m.fds --two-level-run --patch 16 47 16 47 --steps 40 --outdir o --chid t   # -> o/t_driver_profile.txt
```
Next (not done): a quiet machine (the laptop CPU, 24 cores) for the 4-rank rows; device numbers wait for the device-AMR hooks (WP9, WP9b, WP12).
