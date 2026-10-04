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
  to the solver tolerance. Bound to Role 2's composite MLMG, see the next section. A solve is skipped when max|div u - D| is already within the bound.
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
  does not take: masked cells (OBST), non-box level 0, domain walls other than periodic/Neumann, stretched or cylindrical meshes. The ctests carry the label `projection_real_solver`;
  the transport checks in the same tests stay gates.

### Bound to the real composite solver (`PressureBackendSolver`)
- `PressureBackendSolver.{H,cpp}` implements `CompositePoissonSolver` with `pb::PressureProblem::levels` (`pb::PressureLevel`), `pb::PressureWorkspace::rebuild`, `pb::solve_pressure(..., workspace)`
  and `pb::face_gradient_composite` (single-cell directions zeroed; fine-flux overwrite of the D-032(3) kind included). Periodic geometry directions get periodic BC, others Neumann; the
  pressure result `status` other than Ok becomes `SolveReport.ok = false` with the message. `tol_abs` is not used by the backend. CMake builds the `pressure_backend` sources it needs
  as a private static library (no file there is edited) and defines `RT_HAVE_PRESSURE_BACKEND`; without FFT/LSOLVERS in AMReX the real solver is absent and the tests fall back to the mock.
- The AMReX-MLMG mock stays test-only (`--mock-solver`); with it the 2-D three-level case is still NOT GATE.
- Hook: `RegridAmrCore::regrid_dynamic` runs the registered hook after the new hierarchy is complete and before it returns (so before the next dt), when the policy asks for it.
  `DriverAdapter.H`: `install_post_regrid_projection(core, registry, state)` registers a hook that builds the projection levels from the registry, runs `project_after_regrid` and
  copies the projected face velocities to the stage copies (US, VS, WS) if the registry has them. Driver wiring of the install call belongs to the time loop (not edited here).
- Setting `POST_REGRID_PROJECTION` in `&AMR`: `'AUTO'` (default), `'ON'`, `'OFF'` (also .TRUE./.FALSE.). AUTO projects when some level >= 1 keeps old cells and gains new cells in the
  same regrid (the D-063 condition: new fine faces next to retained ones); with no retained fine cells it is not needed. ON: after every regrid that changed the grids. OFF: never.
  The ruling text states the condition but not the key name or its values; the name and the three values are an implementation choice to be confirmed.
  `OFF` makes the parser add a warning (`Report::warnings`) saying that the setting is a diagnostic one and not valid for a gate run; the caller prints it into the .out file with the other
  input warnings, and the gate scripts must refuse a case that carries it. `&MISC EXACT_SUMS` (logical, default F) is read by the same parser so that the line is accepted without error
  (`AmrParams::exact_sums`); the AMR code does not use it.
- Results with the real solver, tolerance 1e-12 relative, acceptance bound 1e-9 u/dx_fine, 1 and 4 ranks the same bound values:
  3-D 32^3 two levels ratio 2: max|div u - D| 1.0 -> 7e-14, at most 9 iterations; 2-D 48x1x48 three levels ratio 2: 1.0 -> 2.5e-13 (8 iterations; converges where the mock did not);
  2-D 16x1x16 (ns2d_16 style, one-cell direction y, two levels, ratio 2): 1.4 -> 6e-13; 2-D 32x1x32 ratio 4 (not a gate): 3.5 -> 1e-12.
  The check also runs after every regrid whether the hook ran or AUTO skipped it: worst value right after any regrid equals the post-projection value in all cases, so AUTO leaves no seam.
- Seam velocity change (largest by level, u = 3): 3-D 1.6e-3 and 3.6e-3 (0.1 %); 2-D three levels 5e-4, 1.4e-3, 2.1e-3; 16-cell case 7e-3, 1.7e-2 (the cells are 4 times coarser);
  ratio 4: 4e-3, 1.3e-2. Profile against distance to the source (cells of the level, 3-D): 3.6e-3, 9.6e-4, 3.8e-4, 2.4e-4, 1.8e-4, 1.5e-4, 1.2e-4, 1.0e-4 at 0..7 and beyond: a factor 4 in the first
  cells, then a slow tail; the 48-cell three-level case has a bump at distance 3 (4.9e-4 after 2.0e-4), where a second level boundary lies.
- Conservation is unaffected: composite mass and species change per regrid 0 with and without the projection, whole run 3e-16 (without 2e-16); AMR remains closer to the uniform-fine run than
  uniform-coarse. Controls: projection off fails the bound; a solver that corrects nothing fails it.
- Known limits: ratio 2 and 4 only; no masks or OBST; the gradient on the fine side of coarse-fine faces is first order at ratio-4 corner patches (up to 12 %), so the ratio-4 numbers are
  reported, not asserted, and gate numbers come from ratio 2; the one-cell direction is handled by the backend with an extruded 4-cell periodic copy (4 times the cells); level 0 must be a box with periodic or Neumann walls.
- The 16-cell case checks the tag-buffer property by report only (blob about 3 coarse cells wide: two tail cells sit one regrid behind); the larger cases assert it.

## Phase 4 note: GPU run of the offload tagging test
The offload build of the tagging kernels was run on one NVIDIA GPU with nvfortran in a standalone program. For Phase 4, use the Integration Lead's stage-1 spike GPU AMReX build to run
the offload tagging test inside the library (device-resident TagBox and FArrayBox data) and compare with the host path bitwise. GPU-enabled AMReX is not needed in Phase 3 (host only, D-047).
