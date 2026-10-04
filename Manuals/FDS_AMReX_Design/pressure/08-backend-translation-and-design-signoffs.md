# 08. Sign-off reviews: pres.f90 host translations (L1211, L1207, L1220-L1222, L1209), the A-68 masked-branch design and the A-67 residual limit

**Status: reviews complete for the items below (code read line by line, tests re-run, three probes added); the masked product-level comparison with HYPRE is blocked, see section 6.** Everything marked "measured" was run by me on the development machine; everything marked "read" is from source only.
Reviewed code: the backend implementer's tree at the pressure-backend tip `f0755644fb` (loop translations unchanged since `99c55244e1`, reader fix `a1b17a25f5`; residual limit `f30ec6b471`; masked design note `55cbbf9432`; masked single-level code `919d75aa41` to `f0755644fb`). Reference: FireX `36975d765f` `Source/pres.f90` (the file in the tree is byte-identical to that commit; the loop tests pin the line ranges by hash).
Accepted spec used as the yardstick: docs 04 (masked-domain pressure), 06 (gauge per zone and connected component, composite volume-weighted mean, D-067) and 07 section 15 (`residual_limit`, measured constant, solution error).

## 1. Answers in short

| Item | Verdict | Required change or test (details in the section named) |
|---|---|---|
| L1211 `pres_compute_rhs_div` (pres.f90:250-260) | **ACCEPT** as a host reference for the 3-D non-cylindrical branch | tests R1 (a real-FDS 2-D and a real-FDS periodic dump) and R2 (face-array index mapping) before it is wired into the AMR path; the cylindrical branch and the transposed IPS layouts are not covered by design (section 2.2) |
| L1207 `pres_p_from_h` (pres.f90:758-764) | **ACCEPT** the function; **D-067 ordering NOT VERIFIED** (no caller exists) | test R3: gauge, then H ghost fill, then P, checked on a singular component; note that this loop is the residual-check P and does not produce `1/RHOP` (section 2.3) |
| L1220, L1221, L1222 `pres_h_bc_x/y/z` (pres.f90:450-462, 466-477, 481-492) | **ACCEPT WITH CHANGE** | periodic wrap needs a whole-domain guard (change C1); the tunnel-preconditioner add-back between the solve and the fill is not represented and must be refused (C2); real-FDS coverage is two code sets only (test R4) (section 2.4) |
| L1209 `pres_poisson_boundary_arrays` (pres.f90:65-228) | **ACCEPT WITH CHANGE** | the wind-array coverage check is weaker than the reads: a high-z open wall reads index KBP1 and the check passes a view that ends at KBAR (change C3, reproduced with an address sanitizer); real-FDS coverage gaps (test R5) (section 2.5) |
| A-68 masked-branch design (note and single-level MLMG code) | **ACCEPT WITH CHANGE** for the single-level MLMG path; **composite masks need a design addition before implementation** | D1 to D8 in section 3 (cache the labels, per-(zone, component) gauge, solid cells out of the C/F interpolation, split-component rule, solid means OBST only under E-2, pins per component for HYPRE, scale tests, vent values from L1209) |
| A-67 size-scaled residual limit | **ACCEPT WITH CHANGE** | the MLMG overshoot of the default tolerance is not handled (change E1); the floor should use the gauge-free H (E2); composite floor and Known-face norm are untested (E3, E4); one real solve above 1e-12 is needed (test R6) (section 4) |

Evidence in one line each: the 93 synthetic bitwise cases and the 23 real-FDS files pass on a fresh build (0 differing elements); all 18 mutants are caught; the three review probes (wind guard, real-FDS coverage census, gauge sensitivity of the floor) are in section 5; the gauge probe is a numpy study of the rounding floor, not a run of the product code.

## 2. Loop-by-loop review against FireX 36975d765f `pres.f90`

### 2.1 What was compared, and how the tests were judged

1. Code: `Source/pressure_backend/FdsPressureLoops.H/.cpp`, `harness/loops_modes.cpp`, `tests/fds_loops_test.py`, `tests/fds_loops/ref_loops.f90.in` and `fds_loops_gen.py`, `tests/m5_fds_loops.cmake`, `frozen/fds-loops-notes.md`, `frozen/fds_loops_cases/fds_loops_cases.tar`, `frozen/fds_loops_dump_hook.py`. Claim register: `docs/amrex/loop_claims.csv` and `loop_work_list.csv` rows L1207, L1209, L1211, L1220 to L1222; the legacy loop map `docs/amrex/pressure-velocity-h-loops.md` and the implementer's status note `docs/amrex/pressure-loops-r2-status.md`.
2. Method: every statement of the Fortran loop was matched to the C++ statement (operand order, sign, index, branch order, which array), the bounds of every array were taken from `init.f90` (not from the notes), and the test data was read for what it does and does not exercise.
3. Re-run (measured, box, g++ and gfortran, `-ffp-contract=off` on both sides as in the test):

| Check | Result |
|---|---|
| 93 synthetic case files against the verbatim Fortran loops (L1211 10, L1207 5, L1220-L1222 70, L1209 8) | all elements bitwise equal (L1211 2,077; L1207 2,027; H fill 19,320; L1209 1,460 elements) |
| 23 files written by a real FDS run (L1211 6, L1207 6, H fill 6, L1209 5) | all elements bitwise equal (L1211 7,326; L1207 9,600; H fill 9,600; L1209 3,430 elements) |
| 18 mutants (`PB_FDSLOOPS_MUTANT` 1 to 18) | 18 of 18 caught |
| Reference text | `pres.f90` in the tree equals `36975d765f` (no diff); the pinned ranges 250-260, 450-462, 466-477, 481-492, 758-764 and 65-228 start and end on the statements the generator expects |

4. Limits of "bitwise": both sides are compiled without floating-point contraction and without fast-math. A production FDS build with contraction (for example an Intel build with `-xHost`, or gfortran with `-march=native`) can differ in the last bit; the real-FDS files came from a gfortran `-O1` scratch build. So the bitwise result is a statement about the arithmetic order and the logic, which is what these tests are for, not a promise about a differently contracted FDS binary. Acceptance of the AMR route against FDS must stay on eps_H or kernel-level tolerance, not on bitwise.
5. None of the five functions has a call site in the product route (read: only `harness/loops_modes.cpp` calls them; the status note says the driver calls the Fortran `fds_p_h_ghost` for physical-face H ghosts and that the C++ version "serves the fine levels", but no fine-level code calls it). Every verdict below is therefore a verdict on the translation, and each "required test" is a test the call site must pass when it is added.

### 2.2 L1211, PRHS divergence (pres.f90:250-260), `pres_compute_rhs_div`

| Point | Finding |
|---|---|
| Statements | `TRM1=(FVX(I-1,J,K)-FVX(I,J,K))*RDX(I)`, `TRM2` in y with `RDY(J)`, `TRM3` in z with `RDZ(K)`, `TRM4=-DDDT(I,J,K)`, `PRHS=TRM1+TRM2+TRM3+TRM4`: identical, same association; the function contains no contraction-prone form that differs from the Fortran |
| Cells touched | `DO K=1,KBAR; J=1,JBAR; I=1,IBAR`: the function loops over the box it is given (FDS: `fds_interior`), reads one face on the low side (`FVX(I-1)`), and writes only `PRHS(I,J,K)` of the box. Checked: the view guards require `FVX(il-1:ih)`, `FVY(jl-1:jh)`, `FVZ(kl-1:kh)`, `DDDT`, `RDX/RDY/RDZ` and `PRHS` to cover what is touched |
| Array bounds | FDS: `FVX/FVY/FVZ/DDDT(0:IBP1,0:JBP1,0:KBP1)`, `RDX(0:IBP1)` etc., `PRHS(ITRN,JTRN,KTRN)` with lower bound 1 and `ITRN=IBAR` if x is periodic, `JTRN=1` if `JBAR=1` (init.f90:2376-2384). The tests use exactly these shapes, including both `ITRN` variants and `JBAR=1` and a 1x1x1 mesh |
| IPS cases | `CASE(:1,4,7)` is the natural layout (IPS 0 and 1 included, which the notes call "1/4/7"); IPS 2, 3, 5, 6 store `PRHS` transposed (L1212-L1214) and are not provided. That is correct for the AMR route (natural layout) and is stated in the header |
| Cylindrical | the `IF (CYLINDRICAL)` loop (pres.f90:236-247, `R(I-1)*FVX..., RRN(I)`) is a different loop and is not translated. The function has no cylindrical argument, so a caller that forgot the refusal would silently get a Cartesian right-hand side. The product interface refuses `cylindrical` with `NotBuilt` (`PressureIface.cpp`, `CompositeSolve.cpp`), so today this cannot happen through the API |
| 2-D | for `JBAR=1` FDS still evaluates `TRM2` with `RDY(1)`; the function does too (tested: five sizes include `JBAR=1`). In a `TWO_D` FDS case `FVY` is zero so the term vanishes; the driver note says the one-cell y direction has "no y term". Both are the same number only if `FVY` is exactly zero (or the term is skipped by design); this has to hold at the call site |
| RHOP versus RHO | not used in this loop |
| Staggered data | FDS `FVX(I)` is the flux on the high-x face of cell `I`. An AMReX x-face array is indexed by the low face of cell `i` (face `i` = FDS `FVX(i-1)`). `F3::from_array4` wraps an array "as is", so wrapping a nodal Array4 and calling the function with FDS indices would be wrong by one cell and silent. No test covers this mapping |
| Tests | 10 synthetic (five sizes times two `ITRN` layouts) and 6 real-FDS files, all 3-D with `ITRN=IBAR+1`, `FVY` nonzero. Not real-FDS: 2-D (`JBAR=1`), periodic x (`ITRN=IBAR`), stretched `RDX` |

**Verdict: ACCEPT.** Tests to ask for: **R1** one real-FDS dump of a `TWO_D` input and one of an input periodic in x (same hook, two small inputs; the hook already works) comparing L1211, L1207 and the H fill; **R2** at the first call site, an index-mapping test: `PRHS` from AMReX face arrays equals the existing driver `PRHS` (cell for cell, bitwise) on a random field, including a box that does not start at 0.

### 2.3 L1207, P = RHOP*(HP-KRES) (pres.f90:758-764), `pres_p_from_h`, and the D-067 ordering

| Point | Finding |
|---|---|
| Statement | `P(I,J,K)=RHOP(I,J,K)*(HP(I,J,K)-KRES(I,J,K))`: identical (the mutants that distribute the product or flip the sign are caught) |
| Cells touched | `DO K=0,KBP1; J=0,JBP1; I=0,IBP1`: the full box with ghost cells in all three directions (`fds_with_ghosts`). Tested: the output array starts at a sentinel (`-7.77e77`) and every element, ghosts and edges included, is compared |
| Array bounds | `RHO/RHOS(-1:IBP1+1,...)` (init.f90:525-526), `H/HS/KRES(0:IBP1,...)`, `WORK7(0:IBP1,...)`; the function needs only 0:IBP1 and the tests give `RHO` its true lower bound -1 |
| RHOP versus RHO | `RHOP=>RHO, HP=>H` in the predictor, `RHOS, HS` in the corrector (pres.f90:716-722). The function takes whichever arrays the caller passes. The pairing is the caller's responsibility and is not tested (no caller) |
| What this loop is | it sits inside `IF (ITERATE_BAROCLINIC_TERM)` of `PRESSURE_SOLVER_CHECK_RESIDUALS`: it forms `P` (`WORK7`) for the residual check of the baroclinic iteration. The baroclinic term itself uses its own loop in `velo.f90` (`BAROCLINIC_CORRECTION`, 3254-3261) with the same `P` formula **and** `RRHO=1/RHOP` in the same loop. The work list calls the two "same formula"; that is true for `P` only, so this function cannot replace the `velo.f90` loop without a second output |
| Ghost data it consumes | `HP` ghost cells must already be filled by L1220-L1222 (the notes' trap 4), and `KRES` and `RHOP` ghosts must be set. In FDS the interior-ghost corner and edge cells of `HP` are never filled by L1220-L1222 (they loop over interior tangential indices only), so the corner and edge values of `P` carry whatever `HP` held; the residual check never reads them (7-point stencil), and the test compares them bitwise only because both sides start from the same random array |

**D-067 ordering (mean removal, solve, gauge, `P=rho_s(H-KRES)`, baroclinic pass 2).** The function is a pure map of its inputs, so the order can only be right or wrong at the caller. What the order requires, in terms of this loop:

1. The gauge shift (`sum rho V (KRES - H) = 0` per zone and component, doc 06 section 5) must be applied to the interior `H` **before** the ghost fill. For a singular component the ghost fill (`H(0) = H(1) - DXI*BXS`, wrap for periodic) is shift-invariant, so a fill done after the gauge gives ghosts that carry the shift; a fill done before it gives ghosts that are off by the gauge constant, and `P` at the ghosts is then off by `rho*c` while interior `P` is right. The baroclinic pass 2 differences `P` across the wall face, so the error enters there.
2. `P` must be formed from the gauged `H` (otherwise the constant that doc 06 section 5 item 2 shows to matter is lost).
3. `RHOP`/`HP` must be the pair of the same stage.

Nothing in the repository tests this order for L1207 (the existing gauge tests check the gauge sum on `H`, not `P` or its ghosts).

**Verdict: ACCEPT** the translation. **D-067 ordering: NOT VERIFIED.** Test **R3** (at the first call site, singular all-Neumann component with nonuniform `rho`): (a) after gauge, fill, then L1207, `sum_V P = 0` over the component to round-off; (b) the same with fill-before-gauge fails by `c * sum_V rho` (negative control); (c) ghost `P` equals `rho_ghost*(H_ghost - KRES_ghost)` with the gauged interior. Also ask for a second output (`RRHO`) or a separate function if the baroclinic loop is to use the same code.

### 2.4 L1220, L1221, L1222, H boundary fill (pres.f90:450-462, 466-477, 481-492)

| Point | Finding |
|---|---|
| Statements and order | the six `IF` lines of each loop are in the Fortran order, so later lines overwrite earlier ones where several apply (for example x code 1 writes the low ghost twice, as in FDS). Compared line by line: low/high Neumann (`HP(0)=HP(1)-DXI*BXS`, `HP(IBP1)=HP(IBAR)+DXI*BXF`), Dirichlet (`-HP+2*B`), codes 5 and 6 (x only, `HP(0)=HP(1)`), periodic. All identical |
| Cells touched | x: `K=1..KBAR, J=1..JBAR`, writes `HP(0,J,K)` and `HP(IBP1,J,K)`; y: `K, I`; z: `J, I`. Edge and corner ghosts are not written, in FDS or in the C++ |
| Array bounds | `BXS/BXF(JDIM,KBP1)` with `JDIM=1` if `JBAR=1`; `BYS/BYF(IBP1,KBP1)`; `BZS/BZF(IBP1,JDIM)` (init.f90:2448-2453); the guards require the interior ranges |
| Boundary types | the function takes the FISHPAK codes (0 periodic, 1 DD, 2 DN, 3 NN, 4 ND, 5/6 axis codes in x). Which wall gets Neumann or Dirichlet is decided upstream (init.f90:2605-2629, one code per face) and by L1209; all faces of a mesh therefore share one pressure type per face, and mixed faces become Dirichlet over the whole face (init.f90:2497-2554). That is FDS behaviour and is not relevant to the loop itself |
| Periodic | `LBC==0` wraps `HP(0)=HP(IBAR)` inside the box it is given. It is correct only if the development machine is the whole extent in that direction (the header says so). A per-box call on a decomposed level silently wraps inside the box |
| Tunnel preconditioner | in FDS, between the solve and these loops, `IF (TUNNEL_PRECONDITIONER)` adds `H_BAR` to `HP` and `BXS_BAR/BXF_BAR` to `BXS/BXF` (pres.f90:436-445). Not represented anywhere; the driver does not refuse it that I could find (read: no `TUNNEL_PRECONDITIONER` handling in the pressure or regrid code; the frozen FDS case records `tunnel_preconditioner = F`) |
| Tests | 70 synthetic cases: every x code 0 to 6 against y/z codes 0 to 4, one 3-D mesh and one `JBAR=1` mesh, whole `H` array compared. Real FDS: 6 files, but only the code sets (LBC,MBC,NBC) = (3,3,4) and (1,3,4); no periodic, no code 2 or 5/6 |

**Verdict: ACCEPT WITH CHANGE.** **C1** add an explicit "whole extent" argument (or an assertion that `b` equals the domain extent in a periodic direction) so a per-box call cannot wrap inside the box. **C2** refuse `TUNNEL_PRECONDITIONER` in the input converter or the driver (one line, with a test), or translate the add-back. Test **R4**: the same two real-FDS dumps as R1 plus one with a Dirichlet low face (a vent on the low-x face) so that codes 0 and 2 meet real data.

### 2.5 L1209, Poisson boundary arrays BXS..BZF (pres.f90:65-228), `pres_poisson_boundary_arrays`

| Point | Finding |
|---|---|
| Wall loop | one pass over the external wall cells in list order; a later wall overwrites an earlier one that names the same entry (tested: 42 to 69 duplicate face entries per synthetic file). The function takes the list in the order the caller gives, so the caller must give FDS's `WALL` order |
| Neumann branch | `BXS(J,K)=HX(0)*(-FVX(0,J,K)+DUNDT)`, `BXF=HX(IBP1)*(-FVX(IBAR,J,K)-DUNDT)`, same in y and z: identical; the sign of the z-high `DUNDT` is mutant 18 |
| Dirichlet, not open (solid) | `0.5*(HP(0,J,K)+HP(1,J,K))+WALL_WORK1(IW)` and the five siblings: identical |
| Dirichlet, interpolated | `(DX_OTHER*HP(1,J,K)+DX(1)*HP(0,J,K))/(DX(1)+DX_OTHER)+WALL_WORK1`: identical; `DX_OTHER` is read from the neighbour mesh at `IIO_MIN/JJO_MIN/KKO_MIN` in FDS and arrives as `PoissonWall::other_d` from the caller. `NIC>1` is not reached (D-072) |
| Dirichlet, open | `TSI` rule (`ABS(T_IGN-T_BEGIN)<=20 eps` and `PRESSURE_RAMP_INDEX>=1`), `P_EXTERNAL=ramp*DYNAMIC_PRESSURE`, `VEL_EDDY` chosen by the **vent's** `ABS(IOR)` (not the wall's), `H0=0.5*(U0**2+V0**2+W0**2)`, wind `H0` with the sign pair per IOR, then the six `SELECT` branches comparing `UU(0,J,K)<0`, `UU(IBAR,J,K)>0`, ...: identical, including `<` against `>` per side (mutants 14 to 17) |
| Array bounds | `U/US(-1:IBP1,..)`, `V/VS(..,-1:JBP1,..)`, `W/WS(..,..,-1:KBP1)` (the cause of the former OPEN-wall difference, fixed in the reader; the function itself never read index -1 and needs `UU(0:IBAR,...)` etc.), `U_WIND/V_WIND/W_WIND(0:KBP1)` (init.f90:718-720), `HX(0:IBP1)`, `DX(1:IBAR)`, `RDXN(0:IBAR)` |
| Wall indices | `BC%II/JJ/KK` is the wall (ghost) cell, not the gas cell (init.f90:3330-3340, `IIG` is the gas cell). So an `IOR=-3` wall has `KK=KBP1` and an `IOR=3` wall has `KK=0`; the tangential indices are 1..N. The real-FDS files confirm it: `KK` of the open top walls is 11 for `KBAR=10` and 9 for `KBAR=8` |
| Wind reads | with `OPEN_WIND_BOUNDARY`, FDS reads `W_WIND(K)` with `K=KK` for a z wall, that is index 0 or `KBP1`. The function requires `u_wind/v_wind/w_wind` to cover **1..KBAR** only and then reads at `KK`. A caller whose wind view starts at 1 passes the check and reads one element past the end. Reproduced (probe, section 5.1): an address sanitizer reports a heap-buffer-overflow at the `W_WIND(K)` read for an `IOR=-3` open wall with `KK=KBP1`. The synthetic tests do not see it because they draw `KK` in 1..KBAR for every wall and give the wind arrays the lower bound 0 |
| Periodic and null walls | a wall whose `PRESSURE_BC_TYPE` is neither `NEUMANN` nor `DIRICHLET` (periodic faces, unset) writes nothing, as in FDS; the 8 synthetic files include type 0 walls. Real FDS: 26 and 18 NULL walls (Neumann type) only |
| Solid obstruction walls | an obstruction at a domain face is an external wall with a solid `BOUNDARY_TYPE`: Neumann branch (covered by real FDS: 414 solid Neumann walls in the closed run, `DUNDT /= 0` on 13 Neumann walls of that run) or, on a Dirichlet face, the pseudo-Dirichlet average plus `WALL_WORK1` (synthetic only: no real-FDS wall of that kind; `WALL_WORK1 /= 0` for 1 wall) |

Real-FDS coverage census (measured from the archive, five L1209 files): wall kinds present are NULL-Neumann, SOLID-Neumann, OPEN-Dirichlet (IOR -3 in the closed run; IOR 1, -1 and -3 in the two-mesh run, wind on, pressure ramp index >= 1 for 120 walls) and INTERPOLATED-Dirichlet. Not present in any real-FDS file: SOLID-Dirichlet (the `WALL_WORK1` average), OPEN walls with IOR +3, +-2, synthetic eddy (`N_EDDY>0`; the dump hook stops on eddies), a `T_IGN` different from `T_BEGIN`, periodic faces. The synthetic files cover all of them but with randomly placed cells (so KK is never a ghost index there).

**Verdict: ACCEPT WITH CHANGE.** **C3** make the wind guard match the reads: require `u_wind`, `v_wind`, `w_wind` to cover `0..KBP1` (or validate `KK` per wall and require exactly the indices read); add a synthetic case with `KK=0` and `KK=KBP1` for z walls with wind on and a wind view whose lower bound is 1 (must throw). Test **R5**: a real-FDS dump with (a) a vent plus an obstruction on one Dirichlet face (SOLID-Dirichlet with `WALL_WORK1`), (b) open vents on y faces and on the +z face, (c) a synthetic-eddy inflow, (d) a ramped open vent with `T_IGN /= T_BEGIN`.

## 3. A-68 masked-branch design (note `frozen/masked-notes.md`, commit `55cbbf9432`; code `919d75aa41` to `f0755644fb`)

The design note was written first; the single-level MLMG code and the `pb_masked` test (`ctest`, 12x10x8 cells, 1 to 3 ranks, six cases against an independent dense numpy reference) landed afterwards. The review covers the note and, where a claim can be checked, the code.

### 3.1 What matches the accepted spec

| Spec point (docs 04, 06, 01 section E) | Design | Verdict |
|---|---|---|
| Gap cells and obstructions removed from the operator, exact no-flux walls | `Solid` class: no unknown, no equation, every face coefficient 0 | matches (the OBST reading needs D5) |
| OPEN vent on a gas-gap face: known value behind the face, face coefficient 2 | `Known` class: flux `2(g - phi)/dx^2`, `g` is the value on the shared face; in the MLMG form `alpha = sum 2/dx^2` on the gas cell and `+2g/dx^2` in the source | matches doc 04 section 2(a2) and FDS's linear ghost `-H + 2*B` (pres.f90:450-492) |
| One null space per connected component, per-component compatibility | flood fill over faces with coefficient /= 0, periodic wrap included, ids ranked by the lowest global cell index (decomposition independent) | matches doc 01 section E and the FDS rule "one mean and one pin per indefinite connected zone" |
| Per-component pin: lowest global index | recorded as `ComponentInfo::pin`; MLMG does not apply it (semi-definite system, mean removed before, gauge after) | matches; HYPRE needs it applied (D6) |
| Singular unless connected to Dirichlet | open if a gas cell touches a Dirichlet domain face after `effective_bc`, or a `Known` cell | matches, including D-057 (a Dirichlet face of a one-cell y direction does not make a component open) |
| D-067 per component: volume-weighted mean removal, then `rho V` gauge, residual over gas cells only | `remove_mean` and `apply_gauge` unchanged and driven by the labels; residual, `||b||`, `||H||` and the floor over gas cells; `Solid` returns 0, `Known` returns `g` | matches doc 06 section 5 items 1, 2, 4 and the draft ADR wording |
| FFT legality | an all-zero or null `cell_class` takes the unchanged unmasked path (bitwise); any mask makes FFT, HYPRE and composite `NotBuilt`; `Auto` picks MLMG | matches doc 04 section 5 (FFT only on a box domain) |
| Covered cells on the single-level API | `NotBuilt` | safe |

Measured by the implementer (`pb_masked`): solution error against the dense reference 5e-13 to 5.8e-12, true residual 1.2e-13 to 7.7e-13 against limits of 1e-12 to 1.9e-12, labels, singular flags, cell counts and removed means exact, results identical on 1, 2 and 3 ranks. Read, not re-run by me (the AMReX build is not available on the loaded box).

### 3.2 Required changes (D1 to D8)

| # | Item | Finding | Required change or test |
|---|---|---|---|
| D1 | Component identification cost | `label_components_masked` runs min-label relaxation with sweeps inside each box and a ghost exchange plus a global reduction per round, on every solve (note: "no cache yet"). The number of rounds grows with the number of box-to-box hops along the longest component path; the stairwell union is 18M bounding-box cells with long passages. Not measured at that size. The global index is an `int` (asserted below 2^31 cells) | cache the labels, the singular flags and the pins in the workspace keyed on the mask (static for a run, new after a regrid); report rounds and time on the stairwell geometry (three components with the `cutk` cut) next to the solve time; keep the 2^31 assertion as a refusal with a message |
| D2 | Gauge key | components are computed from the mask alone; the driver's zone ids are `NotBuilt`. Doc 06 section 5.1 fixes the gauge per zone **and** component | when zone ids exist, the label is the pair (zone, component); a component that spans two zones must be an error (they are one Poisson problem with two background pressures); test: two zones separated by a one-cell solid slab, and a component cut by a zone boundary |
| D3 | Solid cells in the C/F interpolation (composite extension) | `Solid` and `Known` cells carry `phi = 0` (identity row). At a coarse-fine boundary the fine ghost values are interpolated from coarse cells; a stencil that includes a `Solid` coarse cell would read that dummy 0 and corrupt the flux of a gas fine cell next to the corner | rule for the composite masked branch: a coarse `Solid` or `Known` cell must not enter an interpolation stencil (use the nearest gas value, or drop to constant order for that stencil); the face coefficient on a C/F face is 0 when the fine-side neighbour is `Solid`; test: a patch whose boundary touches a solid corner and edge |
| D4 | Component split across refinement levels | the note says components "would have to be found on the composite graph" but gives no rule. Required rule: a component is a set of uncovered gas cells of all levels connected by same-level faces with coefficient /= 0 and by C/F faces between an uncovered coarse gas cell and an uncovered fine gas cell; covered coarse cells are not part of it. One mean removal, one gauge and one pin per component with the composite exact sum (level volumes as weights, doc 06 section 5.1). Per-level labelling would split a corridor whose middle is refined into two or three components and apply different gauges | state the rule in the note; test: a corridor whose middle third is covered by a fine patch (one composite component, one gauge, `sum rho V (KRES - H) = 0` over the composite set) and a sealed pocket inside a patch |
| D5 | What `Solid` means | the note says "OBST or gap cell". The accepted spec (doc 04 assumption A4, doc 01 REC-E1) keeps obstructions on the FDS FFT path with forcing (E-1) for fidelity and offers removing them from the operator (E-2) as the tight option; gap cells are masked in both | say in the note that OBST-as-`Solid` is E-2: the `NO_FLUX`/`WALL_WORK1` iteration for those faces is then switched off (doc 01 section D), `D = 0` in solid cells is irrelevant, and the average-down over a masked fine patch uses gas cells only (doc 01 section E). Gap-only masks are the case the 40 FR-006 inputs need |
| D6 | Pins for the assembled HYPRE branch | the next step in the note is "identity rows for `Solid` and `Known`, `Known` neighbours to the right-hand side". An assembled singular system also needs one pin row per sealed component; `residual_sums` and `evaluate_residual` take a single pin cell | extend the pin handling to a list of pins (the residual excludes all pin rows) and test two sealed components with two pins |
| D7 | Scale, convergence margin, decomposition | tests are 12x10x8 cells on 1 to 3 ranks with one `max_grid_size`; the slowest case needs 199 of the 200 allowed iterations (Dirichlet face plus a thin slab). The study found 3 MLMG iterations on the masked hallways and stairwell with a HYPRE bottom solver (doc 05 and doc 06 section 8); the product path uses a BiCGStab+CG bottom solver with 8 smoothing sweeps instead | run the product masked path on the hallways and the 3-component stairwell (same geometry and RHS as doc 06 section 8) and report iterations and time against the study's Mf rows; vary `max_grid_size` (FR-005); either raise `max_iter` for masked solves or add the Krylov outer loop the note mentions, and add a test that fails when the iteration count is within 5 % of the cap |
| D8 | Where the known values come from | FDS's OPEN value depends on the flow direction (`KRES` of the gas cell when the normal velocity points out, `H0` otherwise, pres.f90:176-210), so `known_value` must be rebuilt at every pressure iteration from the OPEN branch of L1209 (`BXS..BZF` of the vent walls). The other Dirichlet branches of L1209 (solid pseudo-Dirichlet average plus `WALL_WORK1`, interpolated) do not exist under E-2 (exact Neumann on solid walls, composite instead of mesh interfaces) | document the link; test: `known_value` taken from `pres_poisson_boundary_arrays` on a flipping-sign normal velocity changes `g` and the solution accordingly. The doc 04 open point (does UGLMAT put the OPEN value on the face?) is still `[VERIFY]` |

Other observations (no change required): the tolerance is scaled by `|rhs|_inf/|s|_inf` so MLMG's norm and the common-layer residual agree when `Known` terms dominate (reasonable, and the negative control warns); isolated sealed cells get `H = 0` before the gauge (the gauge then sets them to the weighted `KRES`, as the mean removal removes their whole right-hand side); the MLMG coarse-level leak through averaged coefficients changes the rate, not the answer.

**Verdict: ACCEPT WITH CHANGE** for the single-level MLMG masked path (D1, D2, D5, D7, D8 before it carries real cases); **the composite and HYPRE masked branches need D3, D4 and D6 written down before implementation.**

## 4. A-67 size-scaled residual limit (commit `f30ec6b471`; `CommonLayer.cpp` `evaluate_residual`, `CompositeSolve.cpp`)

### 4.1 Compared with the accepted rule (doc 07 section 15; FR-031 text in the requirements)

| Point | Accepted spec | Code | Result |
|---|---|---|---|
| Formula | `max(residual_tol, 10 * 2^-53 * ||A|| * ||H||_2 / ||b||_2)` | `residual_limit = max(o.residual_tol, residual_floor)`, `residual_floor = kResidualRoundoff * 2^-53 * anorm * sqrt(phi2) / sqrt(b2)`, `kResidualRoundoff = 10` | identical |
| Norms | plain Euclidean norms over cells; `||A|| = 4 sum_d 1/dx_d^2` over directions of more than one cell | `phi.norm2`, `rhs.norm2` (unweighted) on a single level; `anorm = 4 sum 1/dx_d^2` over directions longer than one cell; composite: finest-level `dx` | identical; composite uses the volume-weighted sums of the composite residual (`W.vol`) in the ratio, see E3 |
| Compared quantity | residual relative to the right-hand side with the pin row excluded | `residual_check = rel2_nopin` for a singular component, `rel2` otherwise | identical |
| Which components | singular and non-singular alike | floor computed for both; test checks `floor > 0` for Neumann and Dirichlet | identical |
| Small N | base applies up to about N = 77 | limit equals `residual_tol` while the floor is below it; a zero right-hand side falls back to the absolute floor (`bn` replaced by 1) and a zero field gives floor 0 | correct; at N = 48 and 96 the tests see the base and the floor-capable regime respectively |
| Reporting | `residual_floor`, `residual_backward`, `residual_limit` in every case; warning names limit, base and floor | all three reported; warning text "exceeds L (residual_tol T, round-off floor F)" | done |
| Masked | norms over gas cells | `phi` zeroed on non-gas cells, `sumw` = gas count | done |
| Negative controls | loose tolerances still warn | `tol_rel` 1e-6 on MLMG and HYPRE (Dirichlet, N = 48) and `max_iter` 4 warn; synthetic formula cases for singular and non-singular (1.5e-12 passes under a 2.2e-12 limit, 3e-12 warns, a base above the floor wins) | done |
| A-56 eps_H criterion | unchanged | untouched | correct |

### 4.2 Required changes (E1 to E4) and one test

| # | Finding | Required change |
|---|---|---|
| E1 | The spec says the MLMG overshoot of the default tolerance needs one of two answers: run MLMG at `tol_rel` 5e-13, or accept a base of 1.5e-12 (section 15 finding 3: 3 of 24 default MLMG solves reached 1.18e-12 to 1.32e-12, at N = 48 Dirichlet, N = 64 Neumann and one more). The code keeps `residual_tol = 1e-12` **and** `tol_rel = 1e-12`, so at N <= 77, where the floor is below the base, those solves exceed the limit and warn although the solution is fine | do one of the two (the cheapest consistent form: MLMG runs to `0.5 * tol_rel`, as the HYPRE path already runs its Krylov solver to `0.1 * tol_rel`); this is derived from my measured table, not re-run on the new tip |
| E2 | The code evaluates the residual and the floor on the **gauge-fixed** `H`. Rounding in `L H` grows with an additive constant, so the floor follows it. Probe (section 5.3): periodic 96^3, smooth field, an offset of 100 times the rms of `H` gives a residual of 2.05e-12 by rounding alone (3.8e-14 without offset); the gauge-fixed floor is 4.8e-11, the mean-removed floor 4.8e-13. In FDS the `KRES` part of the gauge is of order 1 to 50 against pressure variations of order 0.01 to 1, so ratios of 10 to 1000 are realistic. The implementer's own masked table shows it already: the two `rho` gauge cases at 12x10x8 report limits of 1.5e-12 and 1.9e-12 (above the 1e-12 base at a tiny N) | recommended, not blocking: evaluate the residual and the floor on the mean-removed `H` of each singular component (the system the solver solved), so the verdict does not depend on the gauge constant; if kept as is, state that the limit grows with the gauge offset |
| E3 | Composite floor is unmeasured (also stated in the requirement text). The ratio `||H||_w/||b||_w` uses the level volumes as weights while `||A||` is the finest level's, so mass in coarse regions raises the floor by up to the ratio of the cell areas (4 for ratio 2, 16 for ratio 4) | measure it on the composite cases of doc 07 section 14.2 (record `check`, `floor`, `limit` per solve) before the limit is trusted for composites; if the floor sits more than about 3 times above the measured round-off level, use the level-wise `||A||` weighted the same way |
| E4 | `||A|| = 4 sum 1/dx^2` bounds the stencil of an unmasked cell. A masked cell with a `Known` face has row sum `(5 + 2)/dx^2` per direction instead of `6/dx^2` (up to 8 sum 1/dx^2 with `Known` cells on both sides); the margin `c = 10` against the measured constant (at most 5.6) covers the usual 1.17 factor | no change for now; if masked cases with many vent cells report a floor close to the check, replace `||A||` by the largest row sum of the masked operator |
| R6 | Test | `pb_resid_limit` covers N = 48 and 96 on real solves (the floor never exceeds the base there) and the over-base regime only with synthetic sums. Add one real solve in the regime that motivated the change: HYPRE, Neumann, smooth data, N = 160 (residual 1.5e-12 in my measurement, limit 3.7e-12 by the closed form); mark it as a long test |

**Verdict: ACCEPT WITH CHANGE** (E1 required; E2 recommended; E3 before composite use; R6 for the real-solve evidence).

## 5. Review probes (scratch, reproducible)

All files are under `scratch/pressure-signoff/r2review/` (sources, logs; binaries and case files are regenerable and not kept).

1. **Wind guard (L1209).** `wind_probe.cpp` builds a 3x3x3 context, wind on, one open `IOR=-3` wall with `KK=KBP1`, and `u_wind/v_wind/w_wind` views of lower bound 1 and extent `KBAR` (they satisfy the function's own `covers(1, KBAR)` check). Built with an address sanitizer against the tip's `FdsPressureLoops.cpp`: heap-buffer-overflow, a read of 8 bytes immediately after the 24-byte wind array, at the `W_WIND(K)` read in the open branch. No exception is thrown.
2. **Real-FDS coverage census.** `cov.py` reads the L1209, H-fill and L1211 records of `frozen/fds_loops_cases/fds_loops_cases.tar` and prints wall kinds, orientations, ghost indices, wind, ramp, eddy, `T_IGN` windows, duplicates and codes (numbers in sections 2.4 and 2.5).
3. **Gauge sensitivity of the floor.** `gauge_probe.py` (numpy, periodic 7-point operator, FFT solution): residual by rounding and both floors for offsets of 0, 1, 1e2, 1e4, 1e6 times the rms of `H` at N = 32, 64 and 96 (N = 96: residual 3.8e-14, 4.4e-14, 2.1e-12, 1.8e-10, 2.1e-8; gauge-fixed floor 4.8e-13 to 4.8e-7; mean-removed floor 4.8e-13 throughout).
4. **Re-run of the implementer's loop tests and mutants.** `mut.sh` (18 injected mutants, all caught); the build used `g++ -O2 -ffp-contract=off` and `gfortran -O2 -ffp-contract=off`.

## 6. Parts 3 and 4 of the request: status

**Part 3, masked product-level comparison with HYPRE: not run, blocked on the A-68 implementation.** Checked in the backend implementer's tree: the masked branch has landed only for the single-level MLMG backend (code in three commits ending at the masked-test commit; the design note still has uncommitted status edits). The note's own status table lists the assembled HYPRE backend with masks, composite solves with masks and driver-supplied component ids as `NotBuilt`, and the selector refuses masks for FFT. A product-level HYPRE answer on a masked case therefore does not exist yet, and the comparison was skipped without waiting. What would be possible now, but was not asked and was not run: product masked MLMG against the study harness's assembled masked HYPRE on the hallways and stairwell geometry (this is also the D7 request). The action table in the README keeps the masked product-level comparison open.

**Part 4, cause of the two unstable `strD8` CPU timing rows: cause not proven; ruled-out list and open hypothesis recorded in doc 07 section 7.4.** Per-run logs show identical iteration counts (54 and 21) and identical error norms in all repeats, stable first solves, equal memory, held clocks and cool package in the slowest run, an idle GPU, and no sibling threads; the slow runs are explained by none of these. The open hypothesis is foreign load on a pinned core during the run (the quiet gate samples only before the start). The deciding rerun (per-CPU load sampled during the run, eight Mb repeats, per-solve times) was not done because the quiet gate failed at the single check made for it: a foreign process held one pinned performance core.

## 7. Left open by this review

- Tests R1 to R6 and changes C1 to C3, D1 to D8, E1 to E4 are requests to the backend implementer; none of them was implemented by me (the backend tree is read-only for this work).
- The masked product-level comparison with HYPRE is blocked (section 6).
- Not re-run by me: the AMReX-based tests (`pb_masked`, `pb_resid_limit`); the development machine has one core and a load far above one, so the review relies on the implementer's reported numbers for those and on my own runs for the loop tests and the probes.
