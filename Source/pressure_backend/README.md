# pressure_backend (Role 2): solver-agnostic pressure interface, milestone M1

Single level, CPU. FFT and single-level MLMG backends behind `PressureIface.H`; the common layer wraps them.

| file | content |
|---|---|
| `PressureIface.H/.cpp` | problem/options/result types, selector, `solve_pressure()` |
| `PressureBackend.H`, `FFTBackend.cpp`, `MLMGBackend.cpp` | backend contract; `FFT::Poisson`; `MLPoisson` at `setMaxOrder(2)` |
| `CommonLayer.H/.cpp` | components and pins, mean removal, gauge, reference operator, true residual |
| `ExactSum.H` | private mask-aware volume-weighted exact fixed-point sum (decomposition independent) |
| `harness/` | standalone CMake project and key=value driver `pb_harness` (no FDS sources) |
| `tests/` | CTest registrations and `pb_test.py` |
| `frozen/` | see `frozen/README.md`; `mean-removal-vs-fds.md` (note), `stretched_study.py`, `fds_dump_hook.py`, `fds_cases/` (FDS study) |

Mean removal and gauge (`CommonLayer.H`): exact sums, per singular component. On uniform cells `remove_mean` is the
ULMAT arithmetic mean removal of the volume-scaled RHS; `apply_gauge` takes optional `gauge_weight` (rho) and
`gauge_offset` (KRES) for the FDS gauge. Per-cell volumes (`MeanKind`) are implemented for the future volume-scaled
backends. See `frozen/mean-removal-vs-fds.md`.

Selector: FFT only on one level with no masked or covered cells and uniformly closed (Neumann/periodic) or
uniformly open faces. MLMG on request for the same problems. Masked, composite, variable-coefficient, non-uniform cell widths, cylindrical, and mixed
open/closed single-level requests return `Status::NotBuilt` with a message; `phi` is left untouched.

Build and test (AMReX with MPI, OMP, FFT, LSOLVERS; the install config also needs a Fortran compiler in CMake):
```
cmake -S Source/pressure_backend/harness -B /workspace/pb-build \
  -DCMAKE_PREFIX_PATH=/workspace/amrex-install-hypre232 \
  -DHYPRE_ROOT=/workspace/firemodels-gnu/libs/hypre/63331f19 \
  -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_BUILD_TYPE=Release
cmake --build /workspace/pb-build -j6
(cd /workspace/pb-build && ctest --output-on-failure)
```
Single run: `mpirun -np 4 pb_harness n_cell="64 64 64" bc=neumann max_grid_size=32 backends="fft mlmg"`.
