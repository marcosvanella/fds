# D-057: FDS boundary types to FFT / MLMG pressure BCs, one-cell directions, singular problems

Check of the mapping that the driver's `PressureBcMap` hands to `pb::PressureProblem::bc`, against the backends in this directory.
Test: ctest `pressure_bc_map` (harness modes `bcmap`, `mixed1`, `fftraw`; dense numpy references; 1 and 2 ranks).
Inputs: the driver note `Source/driver/notes/fft-thin-direction-check.md` (read only). Related: m2-notes.md (mixed faces,
`effective_bc`), composite tests `pb_comp_faces_plane2d_*`.

## Answers

1. **What the driver passes for one-cell directions.**
   - One-cell y of a TWO_D case: `NN` (or `PP` when the level is periodic in y), except when both solved directions are fully
     open (`DD`,`DD` in x and z): then y takes `DD` so that the face set is "six open". A Dirichlet y face is never produced for
     the verification set otherwise (all 286 thin inputs have a one-cell y with code 3).
   - One-cell x or z: Neumann or periodic faces pass as `NN` / `PP`; any Dirichlet face is refused by the driver before the solve.
2. **Can D / DD / ND occur in a one-cell direction?** Yes, in exactly one place: the y faces of the fully open 2-D box
   (`DD,DD,DD`, `Fires/tmp_lower_limit_default`). A Dirichlet face in a one-cell x or z cannot reach the interface from the
   driver (refused), and if it does the interface refuses it too (NotBuilt, message names the direction).
3. **How the interface treats it.**
   - Raw `amrex::FFT::Poisson` sets the factor of every length-1 direction to 0 whatever its boundary type. This equals the
     FDS operator (no y term in TWO_D) for any y type, so the FFT result with y = PP, NN, DD, ND or DN is the 2-D solution
     (tested: relative difference to the dense 2-D operator 7e-14 for all five). It disagrees with the *full* 3-D operator
     (which keeps the term -2 phi/dy^2 of a Dirichlet face) by a relative 4e2 for DD in y (tested, negative control). So raw FFT
     is right for y only in the sense of FDS; it must not be compared with a full operator.
   - The interface therefore defines the problem as the 2-D one: `effective_bc()` turns the Dirichlet faces of a one-cell y into
     Neumann for the backend (FFT and MLMG), the composite solve, the reference operator of the true-residual check, and
     the singular-component test. Before this change the all-open 2-D box was reported with a false true-residual warning because
     the reference operator contained the y term; now the residual is at round-off (6e-15 FFT, 3e-13 MLMG).
   - A one-cell y with Dirichlet faces and closed x, z faces is **singular** (the y term does not exist, nothing pins the level):
     mean removal and the zero-mean gauge apply. With x or z open it is non-singular.
   - A Dirichlet face in a one-cell x or z direction: **NotBuilt** (single level, FFT or MLMG, and composite): raw FFT drops
     a term that FDS keeps (tested: relative error 4e2 against the full dense operator). With N or P faces in a one-cell x or z the result is exact
     (tested against the full dense operator, 1e-13).
   - The interface never changes the driver's strings: `DD,DD,DD` in, effective `DD,NN,DD` in the solve.

## Sweep strings (inputs counted in the driver note), one-cell y, grid 16 x 1 x 16

| string (x,y,z) | inputs | backend (Auto) | singular | result |
|---|---|---|---|---|
| PP,NN,PP | 91 | FFT | yes | Ok, equals the 2-D dense operator |
| NN,NN,NN | 86 | FFT | yes | Ok |
| PP,NN,NN | 25 | FFT | yes | Ok |
| DD,DD,DD | 1 | FFT | no | Ok; y acts as N (effective_bc); y term dropped |
| ND,NN,NN | 17 | FFT | no | Ok (was NotBuilt: mixed faces) |
| DD,NN,NN | 4 | FFT | no | Ok |
| NN,NN,ND | 4 | FFT | no | Ok |
| DD,NN,ND | 3 | FFT | no | Ok |
| ND,NN,ND | 2 | FFT | no | Ok |
| DN,NN,DD | 1 | FFT | no | Ok |
| FDS pressure code 0 on non-periodic directions | 23 | not reached | | refused by the driver (Architect decision); the interface itself refuses a Periodic BC on a non-periodic Geometry (InvalidInput) |
| no PBC line | 29 | not reached | | refused in level-0 assembly by the driver |

Every Ok row was also solved by an explicit MLMG request (same answers, 1e-13 to the dense reference) and on 1 and 2 ranks.
Extra strings with a Dirichlet y face (never produced by the driver): `NN,DD,NN` and `NN,ND,PP` are singular, `DD,DN,NN` is not.

## Other checks

- Singular problems (all-N, all-P, N/P mixes, and a one-cell y with D faces and closed x, z) use `FFT::Poisson`; the code
  never uses `PoissonHybrid` (the test greps `FFTBackend.cpp`). Mean removal is the D-067 default (volume-weighted).
- Periodicity guards: a Periodic BC on a non-periodic Geometry, a Periodic BC on one face only, and a Neumann BC on a periodic
  Geometry are `InvalidInput`.

## Open on the driver / Pressure Lead side

- FDS solves a one-cell x or z with Dirichlet faces with all its terms; this is not built (MLMG's 7-point operator would
  keep the term and could be enabled for that case after a dense check; not requested).
