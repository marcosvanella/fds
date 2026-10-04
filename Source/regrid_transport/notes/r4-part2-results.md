# R4 part 2: results of the face path, mass-weighted restriction, dynamic blob and GPU tagging checks

## Restriction of species (steering item)
- Average-down of species uses rho and rho*Z over covered cells and derives Z = (rho*Z)/rho. `LevelOps::average_down_species` is used by
  `average_down_registry` and by the coarse-fine ghost hook (`DriverAdapter.cpp`); no other restriction path touches mass fractions.
- `tests/test_species_avgdown.cpp`: bare level, registry and hook; negative control with a linear Z average changes species mass (1e-1 relative on the
  test data) and the test fails against it; RHO/ZZ/ZZS agree bitwise on covered cells with the driver's `RegistryTransfer::hierarchy_done`.

## Face path (D-060)
- `prolong_faces_level`, `average_down_faces` (all components), `max_divergence_error_local` in `FaceTransfer`; `tests/test_facetransfer.cpp`.
- All-new fine cells have the parent's divergence to about 1e-14. Next to retained fine faces (seams) max|div u - D| is of the size of the velocity
  gradient times the cell width; reported by `test_blob_registry`, not asserted. The projection of the pressure solve removes it each step.

## Dynamic blob end to end (`tests/test_blob_registry.cpp`)
- Real RegridAmrCore, LevelRegistry and RegistryTransfer; mock conservative upwind transport of rho*Z with a prescribed constant-in-time velocity;
  2-D and 3-D, ratios 2 and 4, up to 3 levels, regrid every 4 steps.
- Composite mass and species change per regrid below double resolution; no negative rho*Z, no clips; same hierarchy and data hash at 1, 2 and 4 ranks
  (`run_regrid_rank_check.sh`); AMR closer to the uniform-fine run than uniform-coarse; cells beyond the dependence cone equal the uniform-fine run bitwise.
- Negative controls: no interface flux overwrite drifts mass; linear Z restriction breaks mass; no tag buffer leaves feature cells off the finest level.
- `&AMR_REGION` must be given: the default region is empty and all tags are discarded.

## GPU tagging kernels
- `tests/tag_kernel_check.F90` built twice (host default, `RT_OFFLOAD`), each against an independent cell-by-cell reference, outputs compared bitwise
  (`regrid_transport_tagkernel_host_vs_offload`). Without an accelerator the offload build runs its target regions on the host.
- Device run: `nvfortran -fast -Mpreprocess -mp=gpu -gpu=mem:managed -DRT_OFFLOAD rt_tag_kernels.F90 tag_kernel_check.F90`; 38 cases identical to the host
  build (kernels launched on the GPU). CMake: `-DRT_OFFLOAD_FLAGS="-mp=gpu;-gpu=mem:managed"` with nvfortran.
- Using the offload kernels inside the C++ library needs an AMReX build with GPU support (device-resident TagBox/FArrayBox); the installed AMReX is CPU only.

## Post-regrid composite projection (rulings D-062, D-063)
- D-062: restriction, average-down and prolongation of species act on rho and rho*Z only, Z derived (see above; unchanged).
- D-063: a regrid that makes new fine faces next to retained ones leaves max|div u - D| of the size of the velocity gradient at those seams (about 1.0 in the 3-D case,
  3.5 for ratio 4, measured with a velocity that is exactly divergence free on every level). The driver then projects over the whole composite hierarchy:
  Lap phi = div u - D on the uncovered cells, u -= grad phi, fine faces overwrite coarse faces. No band re-prolongation.
- Hook: `PostRegridProjection.{H,cpp}`: `project_after_regrid(levels, solver, options)`, `ProjectionLevel` (non-owning view of one level: geometry, ratio to the coarser level,
  three face MultiFabs, optional D, covered mask), abstract `CompositePoissonSolver` with `rebuild(levels)` (hierarchy changed, no pointers kept) and `solve(rhs, grad, tol_rel,
  tol_abs, max_iter)`. The solver returns the gradient of its own solution on the faces of every level (its own coarse-fine treatment), which is what makes div (u - grad) = D hold
  to the solver tolerance. Role 2's `PressureWorkspace::rebuild` and `PressureLevel` vector do not exist yet in `Source/pressure_backend` (only the single-level `PressureIface.H`);
  the names here are placeholders to be aligned or wrapped when those headers land. A solve is skipped when max|div u - D| is already within the bound.
- Acceptance: `ProjectionReport::accepted` = max|div u - D| after <= accept_abs + accept_rel * (value before). Switch: `ProjectionOptions::enabled` (false only measures; the
  Phase 3 default stays off per the ruling) and, in the tests, `--report-only` (numbers printed, bound not asserted) against the default (bound asserted).
- Mock solver `tests/MlmgCompositeSolver.H` (amrex::MLPoisson/MLMG, periodic). It needs two workarounds worth passing to Role 2: the face gradient along a hidden (single-cell)
  direction is not zero and must be zeroed (it broke the telescoping of the composite divergence), and MLMG with a hidden direction and three AMR levels did not converge in the
  2-D 48x1x48 case (stagnated at about 5 %); 2-D with two levels (ratios 2 and 4) and 3-D with two levels converge in 5 to 28 iterations.
- Results (mock solver, 1 and 4 ranks identical): 3-D 32^3 two levels: max|div u - D| 1.0 -> 5e-9 (bound 1e-9 u/dx = 1.9e-7); 2-D ratio 4: 3.5 -> 4e-8 (bound 3.8e-7). Largest
  velocity change 0.1 % of u (3-D) and 0.4 % (ratio 4) at the seam. Profile of the change against distance to the nearest source cell (3-D, in cells of the level): 3.8e-3, 1.0e-3, 4.9e-4,
  2.8e-4, 2.0e-4, 1.6e-4, 1.3e-4, 1.0e-4 at 0..7 and beyond: it falls by a factor 4 in the first cells and then slowly (algebraic tail, a dipole-like source), not exponentially. Ratio 4: 1.3e-2 down to 5e-4.
  Composite mass and species per regrid and over the run are unchanged by the projection (the velocity does not enter the flux overwrite balance): per regrid 0, run 3e-16 with and without.
  The AMR result stays closer to the uniform-fine run than the uniform-coarse one (L1 of rho*Z_1 in 3-D 1.83e-3 against 6.35e-3).
- Negative controls: projection off fails the bound in every regrid; a solver that corrects nothing fails it; the tests print both. (`regrid_transport_blob_registry*`.)
- NOT GATE (no projection support yet, reported only): (a) the 2-D three-level case with the mock solver (does not converge, above); (b) every case the mock/real composite solver
  does not take: masked cells (OBST), non-box level 0, domain walls other than periodic/Neumann, stretched or cylindrical meshes. The ctests carry the label `projection_mock_not_gate`;
  the transport checks in the same tests stay gates.

## Phase 4 note: GPU run of the offload tagging test
The offload build of the tagging kernels was run on one NVIDIA GPU with nvfortran in a standalone program. For Phase 4, use the Integration Lead's stage-1 spike GPU AMReX build to run
the offload tagging test inside the library (device-resident TagBox and FArrayBox data) and compare with the host path bitwise. GPU-enabled AMReX is not needed in Phase 3 (host only, D-047).
