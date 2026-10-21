# 17 — FFT as the default single-level pressure solver (delta design memo)

Status: **DESIGN** (no implementation), 2026-10-21. Written as a delta against the src tip `0c47fd77c6`, AMReX 26.09 (`99ddfda`), FDS reference `ce1f659`. Line numbers are at those commits. Nothing pushed.

Owner decision: a single-level run (`--levels 1`, no regrid to finer levels) on a cuboid domain with no OBST must use `amrex::FFT::Poisson` by default; MLMG stays the default otherwise; an explicit override exists.

## 0. What exists, and why CSIRO runs on MLMG today

| Piece | Where | What it does now |
|---|---|---|
| FFT backend | `Source/pressure_backend/FFTBackend.cpp` 11-19 (`fft_bc`), 31-40 (`prepare`), 42-53 (`solve`) | Single-level `FFT::Poisson<MultiFab>`; **already per face**: each of the six faces maps independently Neumann→`even`, Periodic→`periodic`, Dirichlet→`odd`; plan cached on BA/DM/Geometry/BC (24-28). |
| Single-level selector | `PressureIface.cpp` 94-144 (`select_backend`) | Refuses stretched (107-115), variable coefficients (116-118), masked+FFT (120), covered cells (124), Dirichlet in a one-cell x/z (126-133). `Auto` → FFT when unmasked (142). **The code does not test for mixed open/closed faces**; the header comment `PressureIface.H` 367-368 and the `FFTBackend.cpp` line 1 comment ("uniform closed or uniform open faces") are stale. |
| Composite selector | `CompositeSolve.cpp` 289-370 (`select_composite`) | `requested == FFT` → NotBuilt (292-294); `Auto` → HYPRE if masked (367), else **MLMG** (368). Used whenever `PressureProblem::levels` is non-empty (`PressureIface.cpp` 100), **including one level**. |
| Driver, `--run` (S5 loop) | `driver/TimeLoop.cpp` 1354-1370 (`solve_level0`), 1362 | Hard-codes `BackendKind::FFT`; BC from `bc_type[]` (PressureBcMap.H 31-57, already mixed per face). |
| Driver, `--two-level-run --levels N` | `TimeLoop.cpp` 1524-1600 (`solve_composite`), 1575-1590 | Always fills `p.levels` (also for N=1) and requests `Auto`. |

**Correction to the premise.** The CSIRO `--levels 1` case (e.g. `csiro_uniform_0p625`) is not on MLMG because of a mixed-face refusal: it is on MLMG because `--two-level-run` always builds a composite problem (`p.levels.size() == 1`) and `select_composite` never returns FFT (CompositeSolve.cpp 292-294, 366-368). The FFT per-face machinery already handles a Neumann floor with Dirichlet sides and top (`fft_bc`, FFTBackend.cpp 11-19, used per face at 35-38). The delta is therefore mainly **routing** plus three gaps: partly open faces (`face_neumann`), the override flag, and the log line.

## 1. BC mapping and inhomogeneous data

### 1.1 Mapping (unchanged, confirmed)
| FDS | FDS code (pres.f90 55-60, PressureBcMap.H 4) | pb::BC | FFT::Boundary | Transform |
|---|---|---|---|---|
| OPEN (and `P_EXTERNAL`, mesh-interpolated faces at a domain face) | 1/2/4 Dirichlet side | Dirichlet | `odd` | DST |
| MIRROR, INERT / solid wall, other SURF walls on the domain face | 3, or the Neumann side of 2/4 | Neumann | `even` | DCT |
| PERIODIC | 0 | Periodic | `periodic` | DFT |

AMReX 26.09 supports this per face for cell-centred data, Dirichlet/Neumann **at the face** (half-cell from the first cell centre): `AMReX_FFT_Poisson.H` 41 (constructor takes one `pair<Boundary,Boundary>` per direction), 313 (periodic), 330-341 (odd-odd, odd-even, even-odd → wave-number offset 0.5, i.e. the DST-II/DCT-II mixed cases), 343-357 (eigenvalue `k2`; `k2 == 0` is left undivided at 353). Constraints: a periodic direction must be periodic on both sides (51-52; `PoissonHybrid` asserts it at 189-192). Nodal data is not offered for Poisson; the project uses cell-centred H, so this matches FDS's staggering (H at cell centres, BXS etc. at the face: pres.f90 452-455 sets the ghost as `2 BXS - H(1)` for Dirichlet and `H(1) - dx BXS` for Neumann).

### 1.2 Inhomogeneous OPEN H and wall Neumann terms
- FDS: `BXS..BZF` hold the face value of H (Dirichlet; pres.f90 105-106, 125-128, 188-196) or dH/dn (Neumann, `-F_n - d(u_n)/dt`; pres.f90 81-83) and Crayfishpak `H3CZSS` (pres.f90 346; pois.f90 187-) folds them into the right-hand side before the transform.
- FDS-AMReX: identical idea, already shared by every backend: `pb::fold_boundary_data` (CommonLayer.cpp 587-595, `fold_one_level` 532-585) adds `-2 g / h^2` for a Dirichlet face and `±g / h` for a Neumann face to the boundary-cell RHS, then solves the homogeneous problem. Called for the single-level path by `fold_then_solve_level0` (TimeLoop.cpp 1256-1265) and for the composite one by `fold_composite` (1267-1280). This is exactly a homogeneous-BC RHS correction, so `FFT::Poisson` (homogeneous `odd`/`even`) needs **no change**: the FFT route must keep calling the fold before `solve_pressure`.
- Can it match FDS's own FFT? On a **single mesh** yes, to round-off: same 7-point operator, same face-centred Dirichlet/Neumann closure, same folded data (Crayfishpak's tri-diagonal in the last direction vs a full 3-D transform only changes round-off). Measured precedent: MLMG vs FFT 7e-14 (doc 07 §7.1). Against **multi-mesh** stock FDS it cannot match: FDS solves per mesh and leaves interface velocity errors controlled by `PRESSURE_ITERATIONS`; this is an INFO comparison only (section 4).
- Partly open faces (`PressureProblem::face_neumann`, set by `set_face_neumann`, TimeLoop.cpp 1215-1225: an OPEN face with `SOLID_BOUNDARY` patches) cannot be expressed as one FFT boundary type per face: **excluded** from the FFT route (stay on MLMG). CSIRO Case_C064_LS has fully OPEN sides/top and a fully wall floor (ROAD/IGNITION/LINE vents are wall SURFs), so `face_neumann` is empty.

### 1.3 Stretched grids
`FFT::Poisson` needs uniform spacing; `PoissonHybrid` allows a non-uniform z only (`AMReX_FFT_Poisson.H` 163-, 195-196 z non-periodic, 227-237 `dz`, 458-459). FDS-AMReX does **not** build stretched grids on any backend (PressureIface.cpp 107-115; CompositeSolve.cpp 303-309). The FFT route therefore requires uniform dx in all three directions; `PoissonHybrid` is out of scope until stretched grids exist (note: FDS's Crayfishpak does support a stretched direction via `HX`, pres.f90 346).

## 2. Pure Neumann / periodic: zero mode and gauge
- FDS: `H3CZSS` detects the singular case (all directions Neumann/periodic, `ELMBDA == 0`; pois.f90 340-375, `ISING = 1` at 375) and removes the weighted mean `PERT = S1/S3` of F (pois.f90 376-407, PERT at 396); no post-solve normalisation of H beyond that.
- AMReX: the `k2 == 0` mode is not divided (AMReX_FFT_Poisson.H 353-355), so the solution carries `scale × (mean of rhs)` as a constant; with a compatible RHS that is round-off.
- FDS-AMReX: the common layer removes the RHS mean before every backend and fixes the constant after (D-067; PressureIface.H 27-35, 66-69; PressureIface.cpp 293-294 `remove_mean`; gauge `sum(rho V (KRES - H)) = 0`). The FFT route inherits this unchanged; no backend-specific gauge. CSIRO (Dirichlet sides) is non-singular, so this matters only for the MMS all-Neumann/periodic cases.

## 3. Auto-detection rule
### 3.1 Predicate (driver, `solve_composite`, before building `p.levels`)
```
fft_single_level =
       nl == 1                                   // TimeLoop: lc.size()==1, i.e. --levels 1 (to.max_level == 0) -> no regrid possible
    && ro.pressure_backend in {"", "auto"}        // no explicit backend
    && pb::obst_list().empty()                    // no OBST of any thickness (ObstList.H), incl. zero-thickness
    && (ccls.empty() || !ccls[0])                 // no masked level-0 cells
    && cthin.empty()                              // no thin (zero-thickness) faces
    && all p.face_neumann[f].empty()              // every domain face uniform per bc
    && geometry uniform (cell_width all equal), Cartesian, not cylindrical
    && PressureBcMap gave no error for any direction
```
Implementation: when true, the driver solves through the single-level `PressureProblem` (empty `levels`, `p.ba/p.dm/p.geom` from level 0, `bc_type`, `fold_then_solve_level0`) exactly as `solve_level0` (TimeLoop.cpp 1354-1370), with `o.backend = Auto`. `select_backend` then returns FFT (PressureIface.cpp 142) and its existing refusals act as the final check; on refusal fall back to the composite path (MLMG) and log why. Keep the check in the driver (it owns OBST/thin/face data), and add a pure helper `pb::fft_single_level_ok(PressureProblem const&, std::string* why)` in `PressureIface.{H,cpp}` for the backend-side part (uniform, unmasked, no face_neumann, no covered cells, D-057) so tests can call it.

### 3.2 Override
Existing flag `--pressure-backend auto|mlmg|hypre` (main.cpp 62; TimeLoop.cpp 1324-1328): add `fft`. `--pressure-backend mlmg` forces MLMG on a single level; `fft` forces FFT and dies with the selector message if not eligible. Note `study_options` turns any of these flags into study mode (true residual every solve, TimeLoop.cpp 1322, 1333); the override must set `o.backend` without forcing study mode (split the backend choice out of `study`).

### 3.3 Log line (rank 0, once per run and after any change)
```
PRESSURE SOLVER: FFT (amrex::FFT::Poisson, single level, BC x=DD y=DD z=ND)
PRESSURE SOLVER: MLMG (composite, <reason>)          e.g. reason "OBST present", "--levels 3", "face_neumann on zmin", "--pressure-backend mlmg"
```
BC letters from `bc_string` (PressureBcMap.H 25).

## 4. Verification plan
References: (R1) **single-mesh stock FDS** at `ce1f659` (Crayfishpak), primary; (R2) MLMG at `tol_rel` 1e-13 (`--pressure-backend mlmg --pressure-tol 1e-13`); INFO only: multi-mesh stock FDS (reports the size of the interface error with default `PRESSURE_ITERATIONS`).

1. **MMS** (harness `pb_harness`, single level): `phi = prod cos/sin` chosen per face type so the exact solution satisfies inhomogeneous face data; BC sets: NNN, DDD, PPP, ND/DN per direction, CSIRO-type (x DD, y DD, z ND), PPN, PDN. Grids 16³,32³,64³,128³: expect L2/Linf order 2.00±0.05; FFT vs MLMG(1e-13) ≤ 1e-11 relative; FFT vs stock FDS single mesh (same RHS dumped) ≤ 1e-12 relative. 1, 2, 4 ranks.
2. **CSIRO no-OBST, uniform 1.25 m, single mesh**: Case_C064_LS-derived, full domain 120 × 120 × 40 m → 96 × 96 × 32 = 295k cells (one MESH; under 2 GB, single-mesh stock FDS ≈ 1-2 h to 120 s on one P-core: fits). OPEN x/y/zmax, wall floor, no OBST. Runs: FFT (Auto), MLMG (R2), stock FDS single mesh (R1), stock FDS 4×4 meshes (INFO). Compare per step: max|div u − D| (gate 1e-9), DEVC time series (max rel. diff), fire-front x (HRR-weighted leading edge at 30/60/90/120 s; criterion |Δx| ≤ 1 cell vs R1 and R2). 1/2/4 ranks.
3. **Timing per solve**: setup_s and solve_s from the pressure CSV (`--pressure-csv`, but timed in a non-study run too via the log), median of ≥ 100 solves, FFT vs MLMG at 1/2/4 ranks.
4. **Machine**: owner-provided NVIDIA test machine (i9-13900HX, 8 P-cores = logical CPUs 0,2,…,14; E-cores 16-31). Pin with `mpirun --bind-to core --cpu-set 0,2,4,6` (4 ranks) and `OMP_NUM_THREADS=1`. Another engineer's 16-rank run may be on the machine: check `top`/`ps` first, use only free P-cores, never exceed 8 P-core ranks in total, and do not time while it overlaps (record load in the log).
5. **Gates / identity**: divergence gate as today (1e-9). Bit-identity: FFT results across 1/2/4 ranks are **not** expected to be bit-identical (the distributed transform changes summation order); expect ≤ 1e-13 relative in H and identical DEVC to 1e-10; the D-067 mean/gauge sums are exact-sum and decomposition independent. Within one rank count, runs must be bit-identical (deterministic FFTW plans: require `FFTW_ESTIMATE`/fixed plan, check `FFT::Info`).

## 5. AMReX version and build
- `amrex::FFT` landed in 24.11 (#4193); `FFT::Poisson` Neumann/Dirichlet in 24.12 (#4202); `PoissonHybrid` user dz in 24.12 (#4229); 26.01 fixed `PoissonHybrid` nx/ny = 1 (#4829). 26.09 has all of it (`Src/FFT/AMReX_FFT_Poisson.H`, CHANGES.md).
- Build: `AMReX_FFT=ON`; CPU needs FFTW3 (found via `AMReX_FFTW_SEARCH`, OpenMP variant optional); CUDA uses cuFFT, HIP rocFFT, SYCL oneMKL. **heFFTe is not needed** (AMReX's own distributed transposes).
- Current builds: `amrex-install-hypre232` (used by the src build) has `AMReX_FFT ON`, `AMReX_GIT_VERSION 26.09` (AMReXConfig.cmake), FFTW3 found (`/usr/lib/x86_64-linux-gnu/libfftw3*`). The CUDA install must be checked the same way before GPU runs (not checked here).

## 6. Implementation breakdown (backend implementer)
| # | File / function | Change | Test | Acceptance |
|---|---|---|---|---|
| 1 | `pressure_backend/PressureIface.{H,cpp}`: new `fft_single_level_ok(p, why)`; fix stale comments H 367-368 and `FFTBackend.cpp` 1 | Reuse the checks of `select_backend` 103-133; add `face_neumann` all empty | `tests/test_units.cpp`: eligible CSIRO BC set, each refusal reason | all reasons distinct strings |
| 2 | `driver/TimeLoop.cpp` `solve_composite` (1524) | If predicate 3.1 holds, call a single-level path (factor out of `solve_level0` 1354-1370: build `p`, `fold_then_solve_level0`, `solve_pressure` with workspace `pws0`) | CSIRO `--levels 1` | log says FFT; div gate passes |
| 3 | `driver/main.cpp` 62, `TimeLoop.cpp` 1324-1328 | `--pressure-backend fft`; backend choice applied without enabling study mode | CLI test | `mlmg` forces MLMG, `fft` on OBST case dies with selector message |
| 4 | log | Lines of 3.3 on rank 0 at first solve and on change | grep in `pb_test.py` | exact text |
| 5 | `harness` + `tests/pb_test.py` | MMS mixed-BC cases of 4.1, rank 1/2/4 | ctest | order ≥ 1.95; FFT–MLMG ≤ 1e-11 |
| 6 | V&V script | CSIRO 96×96×32 FFT/MLMG/stock-FDS + timing table | manual run on test machine | 4.2 criteria; FFT solve faster than MLMG |
| 7 | docs | Update `Source/pressure_backend/README.md`, this doc's status | — | — |

Risks / blockers: (a) the composite path keeps per-level state (`l0wall`, regrid hooks, `dirichlet_fine_h_ghosts`) that the single-level route must not skip for level 0 — review `solve_composite` 1524-1600 for level-0 side effects (e.g. `l0_wall_fv`) before routing; (b) cross-rank bit-identity is not available with FFT (rule change for V&V gates that assume it); (c) `--run` already uses FFT on all inputs including OBST cases without a mask (TimeLoop.cpp 1362), worth a separate check; (d) partly open faces stay on MLMG.
