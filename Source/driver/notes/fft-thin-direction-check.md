# FFT backend: one-cell ("thin") directions and the driver-side boundary mapping (S9.5, task A(ii), D-057)

Driver-side notes for Role 2, who owns the 2-D/singular boundary mapping check of the FFT backend (D-057); `pressure_backend` is not touched. What the driver passes, what it assumes of
the backend, which verification inputs are affected, and what is open on Role 2's side. Code: `PressureBcMap.H` (pure, unit tested), `TimeLoop::Impl::pressure_bc_map`, mode `fds_amr case.fds --pressure-bc`.

## 1. Facts the mapping relies on
1. **FDS** (`read.f90` 703): `JBAR == 1` sets `TWO_D`; the FFT/FISHPAK 2-D operator (`H2CZSS`) has no y term. A one-cell y direction is therefore NOT part of FDS's pressure operator, whatever the y boundary code is.
   A one-cell x or z direction is solved by FDS with all its terms, including the Dirichlet (open) coupling of a one-cell direction.
2. **amrex::FFT::Poisson** (AMReX 26.09, `AMReX_FFT_Poisson.H`, eigenvalue set-up, 3-D class, about line 320): `dxfac[idim] = 0` for every direction with `Domain().length(idim) == 1`,
   whatever the boundary type of that direction (Neumann/even, Dirichlet/odd or periodic). A one-cell direction contributes nothing to the operator in this backend. (`FFTBackend.cpp` always uses `FFT::Poisson`, never `PoissonHybrid`; the driver requests `BackendKind::FFT`.)
3. Consequence: for the y direction of a TWO_D case the FFT operator equals FDS's, for any boundary type of y. For a one-cell x or z direction with a Dirichlet face the FFT operator lacks the term FDS has: wrong pressure. With Neumann or periodic faces in a one-cell direction the term is zero in both: exact.
4. `pb::select_backend` (`PressureIface.cpp`, about line 60) accepts only face sets with zero or six Dirichlet faces ("single level with mixed open/closed domain faces is not built"), counted over all six faces, the ignored y faces included.

## 2. Driver mapping (`map_pressure_bc_direction`, per direction; FDS code of the mesh at the low face / at the high face)
| Case | Result |
|---|---|
| no box at a face | error |
| level periodic in this direction (single mesh code 0, or several meshes) | `PP` |
| code 0 on a non-periodic direction, more than one cell | error (as before) |
| code 0, one cell | `NN` |
| code > 4 (cylindrical axis) | error (as before) |
| more than one cell | low face D if code 1 or 2, high face D if code 1 or 4, else N (as before) |
| one cell, y of a TWO_D case (FDS ignores it) | `NN`, with a reported note if a Dirichlet face was mapped (never needed in the verification inputs, see 3) |
| one cell, x or z, any Dirichlet face (code 1, 2, 4) | **error** with a message (S5 policy: refuse rather than solve a different operator) |
| one cell, x or z, N or periodic | `NN` / `PP`, exact |
| TWO_D y ignored and both solved directions fully Dirichlet (`DD`,`DD`) | y takes `DD` (`ignored_direction_follows_open`), so that item 4 sees six open faces instead of a "mixed" set; no effect on the operator by item 2 |

Run-time report: `fds_amr case.fds --pressure-bc` prints `PBC <name> twod=<0/1> n=<NXxNYxNZ> codes=<LBC>,<MBC>,<NBC> x=<..> y=<..> z=<..>` (two letters per direction: low and high face, `P`/`D`/`N`, `ERROR [reason]`) plus the notes; the same strings are what `pb::PressureProblem::bc` receives.
Unit tests: `driver_unit_tests`, test `pressure_bc_map` (21 checks, 1 and 4 ranks).

## 3. Sweep over the verification inputs (reproduce: `tests/run_pressure_bc_sweep.sh <fds_amr> notes/thin-direction-case-list.txt <work>`; table `notes/fft-thin-direction-cases.csv`)
Selection: the 286 `Verification/` inputs (not `archive/`) whose `&MESH IJK` has a direction of one cell. **All 286 have a one-cell y; none has a one-cell x or z**, so the refusal rule of the table never fires on the verification set. In every case that runs through the set-up the y code is 3 (Neumann both faces), so the Dirichlet-to-Neumann note never fires either.

Result of the set-up-only sweep (strings x,y,z; counts of inputs):
- **Mapped and accepted by the face rule of item 4** (the case may still stop later for other S5-scope reasons: radiation, inhomogeneous boundary data, chemistry, ...): `PP,NN,PP` 91 (Complex_Geometry 35, NS_Analytical_Solution 12, Scalar_Analytical_Solution 44); `NN,NN,NN` 86 (Chemistry 33, Pyrolysis 12, Radiation 16, Species 13, WUI 7, Fires 2, Flowfields, HVAC, Pressure_Effects 1 each); `PP,NN,NN` 25 (Complex_Geometry 15, Scalar_Analytical_Solution 5, Turbulence 3, Miscellaneous 2); `DD,DD,DD` 1 (`Fires/tmp_lower_limit_default`: y takes D, last row of the table).
- **Mixed open/closed faces, refused by `select_backend` for Role 2 (read from `PressureIface.cpp`, not run end to end because these inputs stop earlier for other reasons; 31)**: `ND,NN,NN` 17 (Pressure_Solver 11, Energy_Budget 2, Fires 2, Complex_Geometry 1, Flowfields 1); `DD,NN,NN` 4 (Flowfields); `NN,NN,ND` 4 (Flowfields 2, Pyrolysis 2); `DD,NN,ND` 3 (Atmospheric_Effects `lee_waves`, Fires 1, Species 1); `ND,NN,ND` 2 (Pressure_Solver); `DN,NN,DD` 1 (`stack_effect`). The mixed D/N strings themselves are what FDS has; the FFT strings `even`/`odd` per face would express them, so this is a question of the selector's scope, not of the thin direction.
- **FDS code 0 on non-periodic directions (23)**: `soborot_*` (21, Scalar_Analytical_Solution, PERIODIC_TEST=13 with prescribed velocity) and `Species/bound_test_1`, `_2`: FDS reports pressure code 0 for all three directions although the level is not periodic (no pressure solve is set up). The driver aborts with "FDS pressure code 0 (periodic) on a non-periodic domain direction" exactly as before this change; to run these cases the driver needs a "no pressure solve" mode (FREEZE_VELOCITY analogue), not a boundary mapping. Open item for the Architect.
- **No `PBC` line (29)**: refused earlier at level-0 assembly (11 CYLINDRICAL, IR-002; 10 refined/overlapping/non-lattice mesh sets: `ns2d_16_*`, `random_meshes`, `geom_channel2`, `geom_channel_tmp*`, `check_kappa`, `dancing_eddies_embed`, `dancing_eddies_uglmat_refine`), 7 Pyrolysis TGA-only runs (FDS stops before the time loop), 1 `rotated_cube_0deg_256_stm` whose set-up did not finish in 90 s on the box (the 128 and 256 `_obs` siblings give `PP,NN,PP`).

## 4. Statement on an open boundary in a one-cell direction
- In the y direction of a TWO_D case: FDS ignores it and so does the FFT backend (item 2): an open (Dirichlet) y boundary has no effect on the pressure solve in either; the driver maps it to N and reports it (no verification input has it). The open boundary still acts on the velocity and scalar boundary conditions through the FDS boundary routines as always.
- In x or z with one cell and a Dirichlet face: the driver refuses (error). This needs either a backend that keeps the one-cell Dirichlet term (e.g. the MLMG backend with the 7-point operator, or the FFT with `dxfac` kept for a Dirichlet one-cell direction) or a different decomposition; none is built.

## 5. Caveats for Role 2 (stage-1 plan section 7.6, repeated; the check itself is Role 2's)
- `PoissonHybrid` is wrong for singular problems and the driver never selects it; `FFT::Poisson` is used for all closed (all `NN`, singular, mean-zero) and all-periodic cases. Whether the singular all-Neumann problem (86 cases above) gets the same mean-zero convention as FDS is the item to check on the backend side (not checked by the driver work); the driver supplies `NN` per face and the right-hand side as FDS forms it, without a solvability correction of its own.
- A TWO_D case is a 3-D `FFT::Poisson` problem with `ny = 1`: item 2 gives the 2-D operator exactly; no extra 2-D code path is needed on the driver side.
