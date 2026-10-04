# pressure_backend (Role 2): solver-agnostic pressure interface, milestones M1 and M2 (first composite solve)

CPU. FFT and single-level MLMG backends behind `PressureIface.H`; the common layer wraps them. A first composite
(multi-level) solve is built on `amrex::MLMG` (see "Composite solve" below and `frozen/composite-notes.md`).

| file | content |
|---|---|
| `PressureIface.H/.cpp` | problem/options/result types, selector, `solve_pressure()` |
| `PressureBackend.H`, `FFTBackend.cpp`, `MLMGBackend.cpp` | backend contract; `FFT::Poisson`; `MLPoisson` at `setMaxOrder(2)` |
| `Composite.H`, `CompositeSolve.cpp` | composite (multi-level) solve: selector and layout checks, `PressureWorkspace`, `solve_pressure` composite path, `face_gradient_composite` |
| `HypreBackend.H/.cpp` | HYPRE assembled-matrix backend (IJ matrix of the MLMG composite operator, PCG/GMRES/BiCGSTAB + BoomerAMG, identity-row pin); explicit request only; see `frozen/hypre-notes.md` |
| `CommonLayer.H/.cpp` | components and pins, mean removal, gauge, reference operator, true residual |
| `ExactSum.H` | private mask-aware volume-weighted exact fixed-point sum (decomposition independent) |
| `harness/` | standalone CMake project and key=value driver `pb_harness` (no FDS sources); `main.cpp` single-level modes, `composite_modes.cpp` composite modes (`comp`, `comp_ns2d`, `comp_sel`, `comp_ws`) |
| `tests/` | CTest registrations and `pb_test.py` |
| `frozen/` | see `frozen/README.md`; `hypre-notes.md` (HYPRE backend), `mixed-nd-hierarchy-note.md` (where the mixed N/D excess error sits), `mean-removal-vs-fds.md` (note), `stretched_study.py`, `fds_dump_hook.py`, `fds_cases/` (FDS study) |

Mean removal and gauge (D-067; `PressureIface.H` header comment, `CommonLayer.H`): exact sums, per singular component.
Default mean removal is the composite volume-weighted mean over the uncovered cells (`PressureProblem::mean_kind =
MeanKind::Volume`); the FDS arithmetic removal of the volume-scaled right-hand side is the runtime parity switch
`MeanKind::ScaledArithmetic` (same operation on uniform cells up to rounding). The gauge `sum(rho*V*(KRES - H)) = 0` is always
applied to a singular component, after any backend; `gauge_weight` (rho, default 1) and `gauge_offset` (KRES, default 0) are
`PressureProblem` fields on a single level and `PressureLevel` fields in a composite hierarchy. See `frozen/mean-removal-vs-fds.md`
(section 8) and `frozen/composite-notes.md`.

Composite solve (milestone M2, first version): fill `PressureProblem::levels` (one `PressureLevel` per level: `ba`, `dm`,
`geom`, `ref_ratio` to the next coarser level, `rhs`, `phi`; coarse to fine). When `levels` is non-empty it is
authoritative and the single-level `ba/dm/geom/rhs/phi/nlevels` are ignored. No covered-cell mask is passed: covered cells are
the next finer level's BoxArray coarsened by its `ref_ratio`. `solve_pressure(p, opts, &workspace)` returns `phi` on every level
(valid cells; covered coarse cells = average-down of the fine solution) and the composite true residual on the uncovered cells.
The workspace (`PressureWorkspace::rebuild`, D-058) holds the hierarchy-dependent setup so repeated solves on one hierarchy do not
rebuild it. `face_gradient_composite` returns face-centred dphi/dn per level with the coarse faces at C/F faces replaced by the
average of the fine gradients. Only MLMG does composite; unsupported requests are `NotBuilt` (see the header and
`frozen/composite-notes.md`).

Selector: FFT only on one level with no masked or covered cells; any per-face mix of Neumann, Dirichlet and periodic
faces (each direction: PP, NN, DD, ND, DN), see frozen/m2-notes.md. MLMG on request for the same problems. A Dirichlet face in a one-cell x or z
direction, masked, variable-coefficient, non-uniform cell widths, cylindrical, and multi-level requests through the single-level fields return `Status::NotBuilt` with a message; `phi` is left untouched.

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

## FDS pressure-area loops (see frozen/fds-loops-notes.md)

`FdsPressureLoops.H/.cpp` (CPU host code, namespace `fdsloops`): FDS `pres.f90` loops L1211 (PRHS divergence, IPS 1/4/7), L1207
(`P = RHOP*(HP-KRES)` over the full box), L1220-L1222 (H boundary fill) and L1209 (Poisson boundary arrays), same operation order as the Fortran.
Harness `pb_fds_loops`; ctests `pb_fdsloops_*` (bitwise against the verbatim upstream loops, 18 mutants, drift check, real-FDS dump archive).

## Build options (PB_WITH_HYPRE)

The HYPRE backend is optional at compile time. Everything that names HYPRE is behind the compile definition `PB_WITH_HYPRE`:
`HypreBackend.H` and `HypreBackend.cpp` (both empty without it), `make_hypre_backend` (`PressureBackend.H`), the `hypre` / `hsys` members of
`PbWorkspaceImpl.H`, and the HYPRE branches of `PressureIface.cpp` and `CompositeSolve.cpp`. `PressureIface.H` (`BackendKind::HYPRE`,
`HypreOptions`, the `hypre_*` result fields, `PressureWorkspace::hypre_built()`) is unchanged and needs no HYPRE header.

- **Without `PB_WITH_HYPRE`** (the driver build): no HYPRE type or symbol is referenced, so a link against AMReX alone (no HYPRE) needs no stub.
  `BackendKind::HYPRE` is rejected by the selector, for a single level and for a composite request, before anything else is looked at:
  `Status::NotBuilt` with `message == pb::kHypreNotBuilt` (`"built without HYPRE"`); `PressureWorkspace::hypre_built()` is always false.
  Compiling `HypreBackend.cpp` without the definition is harmless (empty object), so a source list may keep it.
- **With `PB_WITH_HYPRE`**: add `HypreBackend.cpp` to the sources, link an AMReX built with HYPRE, and define it for **every**
  `pressure_backend` translation unit (the layout of `PressureWorkspace::Impl` depends on it; mixing the two settings in one program is an
  ODR violation).
- **Harness** (`harness/CMakeLists.txt`): option `PB_WITH_HYPRE`, default ON when the installed AMReX has HYPRE (`AMReX_HYPRE_FOUND`); the
  library target `pressure_backend` carries the definition publicly. The harness tests the HYPRE backend, so configuring it with the option
  OFF is an error. The build without the definition is covered by `pb_nohypre_obj` (the driver's pressure_backend source list plus
  `HypreBackend.cpp`, compiled without the definition) and two ctests: `pb_nohypre_symbols` (no HYPRE symbol in any of those objects, by `nm`)
  and `pb_nohypre_check` (HYPRE request is `NotBuilt` / `kHypreNotBuilt` single level, composite and with a workspace; FFT and MLMG still solve).
- **Driver** (`Source/driver/CMakeLists.txt`, owned by Role 1): drop `HypreStub.cpp` and the `if (EXISTS .../HypreBackend.cpp)` block. Do not
  define `PB_WITH_HYPRE` unless `HypreBackend.cpp` is also in the sources and the AMReX is built with HYPRE.

## M2 additions (see frozen/m2-notes.md)

- `PressureOptions::trigger` / `full_checks_on` (FR-039): full checks versus the cheap path.
- `PressureWorkspace` caches the FFT plan (keyed by ba/dm/geom/bc); `PbWorkspaceImpl.H` is private. The class layout changed: recompile dependants.
- Harness modes `trigger1`, `comp_trigger`, `fftcache`; ctests `pb_trigger_single_*`, `pb_comp_trigger`, `pb_fftcache_*`.
- `fold_boundary_data` / `BoundaryData`: inhomogeneous Dirichlet/Neumann wall data folded into the rhs (sign convention in the header and frozen/m2-notes.md).
- Mixed open/closed faces (any per-direction NN/DD/ND/DN/PP) on FFT, MLMG and composite; `effective_bc()` for a one-cell y.

## HYPRE backend (see frozen/hypre-notes.md)

- `PressureOptions::backend = BackendKind::HYPRE` (never chosen by Auto; the backend option lives on `PressureOptions`, not on `PressureProblem`), settings in `PressureOptions::hypre` (`HypreOptions`: Krylov, BoomerAMG parameters), tolerance and iteration cap are `tol_rel` / `max_iter`.
- Same discrete operator as MLMG (7-point, second order, MLMG's C/F ghost and reflux), assembled from host buffers; single level and composite (ratios 2 and 4, n levels, one-cell y). Mean removal and gauge are the common layer's. Same `NotBuilt` rules as MLMG.
- The set-up (matrix and AMG hierarchy) is cached in `PressureWorkspace` and dropped by `rebuild`; `PressureWorkspace::hypre_built()`.
- Harness: `backends="fft mlmg hypre"` (single level), `mode=hypre_cmp` (composite HYPRE vs MLMG), `mode=hypre_op` (operator oracle against MLMG), `mode=hypcache`, `backend=hypre` on `comp`, `comp_ns2d`, `comp_sel`, `comp_ws`, `mixed1`, `bcdata`; `mode=err_map` (error by distance to C/F and domain faces). ctests `pb_hypre_*` (21, plus the residual-check tests below).

## Residual check (see frozen/hypre-notes.md, "Residual check")

`PressureResult::residual_rel2` / `residual_relmax` are the full ||b - L H|| / ||b|| (2-norm, max), always reported. For a non-singular
component, or when the backend applied no pin (MLMG, FFT), that is the checked value (`residual_check`). When HYPRE pinned
an unknown (singular problem) the pin cell is excluded: the residual of a conservative operator with a mean-free right-hand side sums to
zero, so the pin row holds the sum of all other residuals (up to sqrt(N) times their 2-norm) and is not an independent equation. The warning is on
`residual_rel2_nopin = ||r_nopin||_2 / ||b||_2` (scaled like MLMG's); also reported are
`residual_relmax_nopin`, the pin row (`residual_pin_rel`, `residual_pin_abs`), `residual_check`, `residual_limit` and the extras
`residual_rel2_mr` (mean-removed; equals the full value) and `residual_backward`. `verbose >= 2` logs full and non-pin values.

**Size-scaled limit (FR-031, A-67; docs/pressure/07 sections 14-16).** The checked value is compared with
`residual_limit = max(residual_tol, residual_floor)`, `residual_floor = 10 * 2^-53 * ||A|| ||H||_2 / ||b||_2`, `||A|| = 4 sum_d 1/dx_d^2`
(finest level, directions of more than one cell), `||H||_2` and `||b||_2` the weighted 2-norms of the returned solution and of the
mean-removed right-hand side. The floor is computed for every component, singular or not. `residual_tol` stays 1e-12, so the limit equals
`residual_tol` while the floor is below it (smooth unit-cube data: about 65 to 100 cells per direction, depending on the problem) and grows like N^2 above that. Warn iff
`residual_check > residual_limit`; the warning text names the limit, `residual_tol` and the floor. `residual_ok` means `residual_check <= residual_limit`.
Harness: `mode=rescheck`, `mode=hypre_resid`; ctests `pb_hypre_resid_*`, `pb_resid_check_*`, `pb_hypre_resid_comp`.
