# pressure_backend (Role 2): solver-agnostic pressure interface, milestones M1 and M2 (first composite solve)

CPU. FFT and single-level MLMG backends behind `PressureIface.H`; the common layer wraps them. A first composite
(multi-level) solve is built on `amrex::MLMG` (see "Composite solve" below and `frozen/composite-notes.md`).

| file | content |
|---|---|
| `PressureIface.H/.cpp` | problem/options/result types, selector, `solve_pressure()` |
| `PressureBackend.H`, `FFTBackend.cpp`, `MLMGBackend.cpp` | backend contract; `FFT::Poisson`; `MLPoisson` at `setMaxOrder(2)` |
| `Composite.H`, `CompositeSolve.cpp` | composite (multi-level) solve: selector and layout checks, `PressureWorkspace`, `solve_pressure` composite path, `face_gradient_composite` |
| `CommonLayer.H/.cpp` | components and pins, mean removal, gauge, reference operator, true residual |
| `ExactSum.H` | private mask-aware volume-weighted exact fixed-point sum (decomposition independent) |
| `harness/` | standalone CMake project and key=value driver `pb_harness` (no FDS sources); `main.cpp` single-level modes, `composite_modes.cpp` composite modes (`comp`, `comp_ns2d`, `comp_sel`, `comp_ws`) |
| `tests/` | CTest registrations and `pb_test.py` |
| `frozen/` | see `frozen/README.md`; `mean-removal-vs-fds.md` (note), `stretched_study.py`, `fds_dump_hook.py`, `fds_cases/` (FDS study) |

Mean removal and gauge (`CommonLayer.H`): exact sums, per singular component. On uniform cells `remove_mean` is the
ULMAT arithmetic mean removal of the volume-scaled RHS; `apply_gauge` takes optional `gauge_weight` (rho) and
`gauge_offset` (KRES) for the FDS gauge. Per-cell volumes (`MeanKind`) are implemented for the future volume-scaled
backends. See `frozen/mean-removal-vs-fds.md`.

Composite solve (milestone M2, first version): fill `PressureProblem::levels` (one `PressureLevel` per level: `ba`, `dm`,
`geom`, `ref_ratio` to the next coarser level, `rhs`, `phi`; coarse to fine). When `levels` is non-empty it is
authoritative and the single-level `ba/dm/geom/rhs/phi/nlevels` are ignored. No covered-cell mask is passed: covered cells are
the next finer level's BoxArray coarsened by its `ref_ratio`. `solve_pressure(p, opts, &workspace)` returns `phi` on every level
(valid cells; covered coarse cells = average-down of the fine solution) and the composite true residual on the uncovered cells.
The workspace (`PressureWorkspace::rebuild`, D-058) holds the hierarchy-dependent setup so repeated solves on one hierarchy do not
rebuild it. `face_gradient_composite` returns face-centred dphi/dn per level with the coarse faces at C/F faces replaced by the
average of the fine gradients. Only MLMG does composite; unsupported requests are `NotBuilt` (see the header and
`frozen/composite-notes.md`).

Selector: FFT only on one level with no masked or covered cells and uniformly closed (Neumann/periodic) or
uniformly open faces. MLMG on request for the same problems. Masked, variable-coefficient, non-uniform cell widths, cylindrical, mixed
open/closed, and multi-level requests through the single-level fields return `Status::NotBuilt` with a message; `phi` is left untouched.

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
