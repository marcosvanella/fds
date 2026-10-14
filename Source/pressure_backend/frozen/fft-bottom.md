# FFT bottom for the composite MLMG (opt-in, off by default)

`PressureOptions::mlmg_fft_bottom` (`pressure.mlmg_fft_bottom`, env `PB_MLMG_FFT_BOTTOM=1`): AMReX `BottomSolver::custom` on an MLPoisson subclass
(`FftBottomPoisson`, CompositeSolve.cpp) that solves the coarsest MG level of AMR level 0 exactly with `amrex::FFT::Poisson`. Conditions: unmasked,
no face_neumann, no one-cell direction, bottom level <= 1 box per rank, mlmg_krylov BiCGSTAB or GMRES. Otherwise the default bottom.
Test: pb_comp_fft_bottom (same solution as the default bottom, used only on request).

Benchmark, the test machine, 4 ranks x 1 OMP, vortex_ring_48_4m two levels (--tag-vort auto --maxsize 8 --blocking 4), 30 steps (61 solves), medians per
solve from <chid>_pressure.csv, ms/step = t_loop / steps (FDSTL_PROFILE). Load average 1.9 -> 4.8 during the runs (24 cores).

| case | arm | iters (med/max) | ms/solve (med) | ms/step | true rel2 (med/max) | max div |
|---|---|---|---|---|---|---|
| static | V-cycle | 9/9 | 20.4 | 260 | 1.1e-10/7.1e-10 | 6.5e-13 |
| static | BiCGStab | 4/5 | 29.3 | 283 | 3.3e-12/9.0e-12 | 9.3e-14 |
| static | V-cycle + FFT bottom | 200/200 (stagnates) | 492 | 1240 | 9.2e-08/2.6e-07 | 4.2e-10 |
| static | BiCGStab + FFT bottom | 6/6 | 44.0 | 336 | 1.7e-12/2.4e-12 | 2.5e-14 |
| adaptive (regrid 5) | V-cycle | 9/10 | 100.1 | 766 | 8.0e-11/8.7e-11 | 2.1e-12 |
| adaptive | BiCGStab | 4/5 | 148.6 | 867 | 4.5e-12/2.0e-11 | 6.7e-12 |
| adaptive | V-cycle + FFT bottom | 200/200 (stagnates) | 2172 | 4879 | 7.3e-08/1.4e-07 | 3.5e-10 |
| adaptive | BiCGStab + FFT bottom | 5/6 | 185.9 | 948 | 6.9e-12/1.0e-11 | 6.3e-13 |

Result: no win. The coarsest MG level is a few cells per direction, so the default bottom is already cheap and near exact; the FFT plan and
redistribution cost more. Plain V-cycles with the FFT bottom stagnate near 1e-7 on this case (not on the harness cases: not diagnosed), so
the option is restricted to Krylov. Variant (2) (FFT on the whole of level 0, MLMG for the fine corrections) was not built.
