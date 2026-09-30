# Review of `Source/pressure_backend/PressureIface.H` from the driver side (S5)

For the Architect and Role 2. Review only: nothing in `Source/pressure_backend/` was edited. Kernel-facing rules for the driver: passive scalars are handled by
`Fields.cpp` (they do not enter the pressure problem); only uniform Cartesian metrics are used (no R/RRN, no CYLINDRICAL/TRN*, IR-002).

## Verdict
The signature (`pb::solve_pressure(PressureProblem&, PressureOptions&)`, the `PressureProblem` fields `ba/dm/geom/bc/rhs/phi`) is **usable as is** for the S5 case
(`shunn3_32`: FFT backend, x and z periodic, y Neumann with one cell, homogeneous boundary data). It is called from `TimeLoop.cpp` (`solve_poisson`) with
`verbose = 0`. `FFT::Poisson` was not needed as a stub. Results on the 80 steps of `shunn3_32`: residual `||b - L phi||` about 1e-15, `removed_mean` 0 or about
1e-17 relative, pressure iteration counts equal to FDS (10/10, then 7/7).

## What the driver had to do around the interface
1. `Geometry` periodicity must match `bc`. FDS uses MBC = 0 (periodic) in y for this case while the level has y non-periodic and one cell; the driver maps a
   1-cell non-periodic direction with FDS code 0 to Neumann (the y term is zero either way). Dies on any other combination.
2. `phi` needs one ghost layer, `rhs` none, same BoxArray/DistributionMapping. The driver copies the result into H/HS and then builds the ghost layer of H
   itself (level fill plus the physical-face assignment block that ends `PRESSURE_SOLVER_FFT`, kept in the driver as `fds_p_h_ghost`).

## What is missing for S6 (request)
1. **Inhomogeneous boundary data.** FDS passes BXS..BZF (non-zero Neumann or Dirichlet data from open vents and moving walls) into the solve. `PressureProblem` has
   no input for them ("Dirichlet: homogeneous in M1, values folded into the RHS" but no helper does the fold). The driver aborts when any boundary datum on a
   solved physical face is non-zero. Request: an RHS-fold helper `fold_boundary(problem, bxs..bzf)` or explicit boundary-value arrays in `PressureProblem`,
   with the FDS sign convention stated (FDS FISHPAK_BC: which faces carry the derivative and which the value, and the 1/2 dx factors).
2. **Mixed open/closed faces** (for example an open vent on one side and walls elsewhere, Dirichlet in one direction and Neumann in another) return `NotBuilt`
   for an unmasked single level. Most FDS cases with vents need exactly this. Request: FFT support for Dirichlet/Neumann/periodic per direction (the FFT
   Poisson of AMReX supports them direction by direction) or the MLMG path for the mixed case.
3. **Obstructions** (`cell_class` masks) are `NotBuilt`; the S5 driver refuses cases with OBSTs.
4. **Performance.** The FFT plan is rebuilt on every call and FDS solves up to 10 times per step. Request: plan caching keyed by `ba/dm/geom/bc`.
5. **Decomposition independence (D-028).** The distributed FFT returns results whose last bits depend on the box layout and the rank count (the driver's own
   exact sums are bitwise independent, see the driver README). If whole-run bit identity across decompositions is wanted, the solve would have to be
   made layout independent (a fixed reduction order in the transposes) or the requirement stated as T2 only.
6. **Multi-mesh FDS cases.** With several FDS meshes FDS solves one Poisson problem per mesh with interpolated boundary data and iterates; the driver replaces
   this by one level-wide periodic/Neumann solve. That is a different discretisation of the coupling (it is the intended AMR behaviour); the comparison
   with the multi-mesh baseline is therefore a T2 statement, not a bitwise one.

No signature change is required for S6 except items 1 and 2.
