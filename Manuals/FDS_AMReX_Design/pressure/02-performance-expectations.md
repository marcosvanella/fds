# 02 — Performance expectations: FDS FFT/HYPRE pressure solve vs AMReX MLMG composite solve

**Status:** Draft v1. Author role: numerical analyst (pressure/projection).
**Base:** FDS FireX `AMReX` branch @ `36975d765f` (this repository, read-only); citations as in 00/01.
**Labels used throughout:** **[MEASURED]** = run on the development machine during this session, with the command and conditions given.
**[ESTIMATE]** = derived; the derivation is shown next to it. **[UNVERIFIED]** = plausible but neither measured nor derived here.

---

## 0. Assumptions

P1. **Hardware.** The development machine: 8 cores ("Intel(R) Xeon(R) Processor", 1 thread/core, 320 MiB L3 reported by `lscpu`), 16 GB RAM with
    ~2.5 GB free (`free -m`: 2458–2551 MB free, ~6.4 GB "available" during this session). No GPU (programme decision D-004).
P2. **No FDS or AMReX runs were possible.** The only prebuilt FDS binaries are the old master (`ce1f659`) builds under `(local FDS master checkout)/Build/*`
    (none under `Build`). They link against
    `libmpi.so.40`, which is not installed (`ldd` shows "not found"). There is no `mpirun` on the development machine. Installing an MPI runtime
    was out of scope for this task. So every FDS/AMReX number below is an **estimate**.
P3. **Micro-benchmarks are proxies.** They use NumPy 2.2.4, single thread (`OMP_NUM_THREADS=1`), on the development machine. NumPy's FFT
    (pocketfft) is a stand-in for Crayfishpak. The NumPy stencil is a stand-in for a compiled AMReX kernel, and is slower
    than compiled code by a factor I did not measure.
P4. **Scope.** Uniform-grid and static two-level (ratio 2) configurations, no subcycling (01 REC-B1), A2-a operator (constant
    coefficient) unless noted. Pressure step only. D, F and transport costs are out of scope.
P5. **Solver settings as recommended in 01:** MLMG tolerance per 01 REC-D2, initial guess from the previous H. FDS reference:
    FFT with iterations (default) or UGLMAT-HYPRE (01 §I.1).

---

## 1. What FDS spends today

| Solver | Work per solve | Solves per step | Notes |
|---|---|---|---|
| FFT (default) | O(N log N) per mesh, direct (Tech Guide `Appendices.tex:2976-2981`; Crayfishpak `H3CZSS`, `pois.f90:187`) | 2 × (pressure iterations): at least 1, at most `MAX_PRESSURE_ITERATIONS=10` (`cons.f90:559`) per half-step. The loop runs while the interface/solid velocity error > `VELOCITY_TOLERANCE` or `PRESSURE_ERROR_MAX > PRESSURE_TOLERANCE` (`main.f90:1724-1741`) | Serial per mesh (00 §8). Parallel only across meshes. Tunnel preconditioner adds a 1-D global solve on rank 0 (`pres.f90:505-697`) |
| ULMAT | sparse direct (PARDISO) or PCG-AMG per mesh | as FFT | exact at solids, still iterates at interfaces |
| GLMAT/UGLMAT (HYPRE) | PCG + 1 BoomerAMG V-cycle per iteration, O(nnz) each; stop at ‖r‖₂ ≤ 1e-12‖b‖₂ or 1000 its (`imkl.f90:238-240`) | 2 × (iterations for baroclinic/solids only) | AMG setup only when matrices are rebuilt (`GLOBAL_MATRIX_REASSIGN`, `main.f90:1806-1835`). FireX adds optional GPU offload of the HYPRE solve (00 §11.2), irrelevant here (no GPU) |

The FireX build writes a `PRES` column (timer 5) and a new `PACK` column (timer 12, MPI pack/unpack) to `CHID_cpu.csv`
(`main.f90:4146`). That is the instrument for measuring the pressure fraction `f_pres` of the step. **`f_pres` was not measured (P2).**

### 1.1 [MEASURED] FFT proxy on the development machine
Command: `OMP_NUM_THREADS=1 python3 bench.py` (script: `numpy.fft.rfftn` + `irfftn` of a float64 cube, best of 5).

| Grid | fwd+inv FFT time | per cell |
|---|---|---|
| 64³ (262,144 cells) | 3.93 ms | 15.0 ns |
| 128³ (2,097,152 cells) | 28.7 ms | 13.7 ns |

A Poisson solve by FFT costs one forward and one inverse transform plus a diagonal scaling, so **~14–15 ns/cell per solve,
single core**, is the proxy for one FDS FFT pressure solve. Crayfishpak uses real fast trigonometric transforms plus a
tridiagonal solve in the stretched direction, so its constant differs. [UNVERIFIED] It is plausibly within a factor of 2.

### 1.2 [MEASURED] Stencil and bandwidth proxies
Same script and conditions:
* 7-point Laplacian apply (NumPy slicing, allocates the result): 64³ → 1.51 ms (**5.8 ns/cell**); 128³ → 14.5 ms (**6.9 ns/cell**).
* Triad-like stream (`x = 3c; x += b`, 20 M doubles, counting 40 B/element): **24.1 GB/s** effective, single thread.
  Peak RSS of this test was ≈ 0.5 GB (three 160 MB arrays); nothing larger was run, to respect the memory budget.

---

## 2. MLMG cost model

### 2.1 Work per V-cycle [ESTIMATE]
Per fine level, one V-cycle does ν₁+ν₂ smoothing sweeps (AMReX default 2+2 red-black GS
[UNVERIFIED default for the target version]), one residual, one restriction and one prolongation. That is ≈ 7–8 stencil-like
passes over the level. Coarser MG levels add a geometric factor Σ 8⁻ᵏ = 8/7.
→ **≈ 8–9 fine-grid passes per V-cycle.**

Cost per pass, compiled and memory-bound: 3 streams × 8 B = 24 B/cell (φ read/write, rhs read). At the measured 24 GB/s
single-thread, that is ≈ 1 ns/cell ideal. Ghost-cell exchange, box overhead (32³ boxes, (34/32)³ ≈ 1.2) and the fact that
GS is not a pure stream give an assumed practical 2–3 ns/cell per pass.
→ **≈ 15–25 ns/cell per V-cycle, single core** (range from 8 passes × 2 ns to 9 passes × 3 ns, rounded).

### 2.2 Iterations [ESTIMATE]
* For the constant-coefficient cell-centred Poisson problem, geometric MG typically reduces the residual by about 10× per
  V-cycle [UNVERIFIED for AMReX defaults; to be measured].
* With the previous H as initial guess, the initial residual is already small. Reaching the 01 REC-D2 tolerance
  (r_max ≤ 0.005 s⁻¹, or relative 1e-10) is estimated at **6–10 V-cycles** for a production solve.
* Equivalence tests at 1e-12 (01 §I.2) need **10–14 V-cycles**.
* The composite multi-level solve adds C/F interpolation and reflux inside the residual. I assume ~10–20% overhead per V-cycle
  versus single-level [UNVERIFIED].

### 2.3 Per-solve comparison [ESTIMATE]
| | ns/cell per solve (single core) | derivation |
|---|---|---|
| FDS FFT (proxy) | ~14 | §1.1 measured proxy |
| MLMG, production tolerance | ~90–250 | 6–10 V-cycles × 15–25 ns |
| MLMG, 1e-12 equivalence tests | ~150–350 | 10–14 V-cycles × 15–25 ns |
| HYPRE PCG-AMG (UGLMAT today) | ~300–1500 [UNVERIFIED] | 10–30 PCG its × (AMG V-cycle ≈ 2–4× a GMG V-cycle, from operator complexity) |

**Consequence.** On a *uniform single box*, one MLMG solve is expected to cost **~6–18× one FFT solve**. FDS multi-mesh runs
often need several FFT solves per half-step because of the interface iteration. The composite solve needs none of those.
So the multi-mesh ratio is smaller: r ≈ (6–18)/(mean pressure iterations).

### 2.4 Impact on the whole step (NFR-030, R-22)
Uniform-mode slowdown = `1 + f_pres·(r − 1)`, where r is the per-step pressure cost ratio (MLMG/FFT).
NFR-030 (≤ 1.25× baseline) therefore requires `f_pres·(r−1) ≤ 0.25`:

| f_pres (unmeasured) | max allowed r |
|---|---|
| 0.05 | 6.0 |
| 0.10 | 3.5 |
| 0.20 | 2.25 |

With r ≈ 6–18 for single-box uniform runs, **NFR-030 is unlikely to be met by MLMG alone** unless f_pres ≤ ~0.03. This supports
keeping the FFT path in uniform single-box mode (FR-037) and measuring `f_pres` early from the FireX `_cpu.csv` `PRES` column
(**first action once an MPI runtime is available**). For multi-mesh uniform runs with many FDS pressure iterations, MLMG may
break even.

---

## 3. Memory per cell for the pressure step [ESTIMATE]

Counted in doubles per cell per level. Face arrays count ≈ 1 double per cell per direction. The factor 1.2 is ghost/box
overhead for 32³ boxes with one ghost layer; 8/7 is the MG hierarchy factor.

| Item | A2-a (`MLPoisson`) | A2-b / E-2 (`MLABecLaplacian`) |
|---|---|---|
| Ũ / U face velocities (owned by momentum, already exist in FDS as U,V,W / US,VS,WS) | 0 (not counted) | 0 |
| φ = H, S = D (exist in FDS; D is `DS`/`D`) | 0 | 0 |
| Projector RHS | 1 | 1 |
| MLMG internal (sol, rhs, res, cor, cor-hold) × 8/7 | ≈ 5.7 | ≈ 5.7 |
| Face β (3) × 8/7 and a-coefficient (1) × 8/7 | — | ≈ 4.6 |
| Subtotal × 1.2 overhead | **≈ 8 doubles ≈ 64 B/cell** | **≈ 14 doubles ≈ 110 B/cell** |
| Overset mask (int, E-2 via C++) | — | +4–5 B/cell |

For comparison: FDS's own pressure work arrays (`PRHS`, `BXS…BZF`, Crayfishpak `SAVE` arrays; 00 §4/§8) are O(1) doubles per
cell. MLMG therefore adds roughly **60–120 B/cell** on top of the FDS state.

HYPRE BoomerAMG (reference solver only): a 7-point matrix is ≈ 7 × 12 B = 84 B/row. With operator complexity 1.5–3
[UNVERIFIED typical PMIS range] plus PCG vectors, that is **≈ 250–500 B/cell**. This is a second reason not to use full-HYPRE
for production on the development machine.

**On the development machine (≈ 2.5 GB free):** the pressure-solver add-on alone would allow ~20–40 M cells, so it is **not** the constraint.
The FDS state per cell is much larger. [UNVERIFIED] It is of order 1 KB/cell or more (`docs/inventory/mesh_fields.csv` lists
1,314 mesh fields, most of them not per-cell). **Keep V&V runs on the development machine to ≤ ~1 M cells in total** until a real peak-RSS
measurement exists.

---

## 4. AMR benefit and subcycling trade-off [ESTIMATE]

Example plume case (for FR-014/NFR-032 sizing only):
* Uniform fine: 128³ = 2.10 M cells.
* Two-level AMR: coarse 64³ (0.26 M) plus a fine patch covering 1/8 of the domain at ratio 2 (0.26 M). Total 0.52 M cells,
  **4.0× fewer cells**.

Without subcycling (01 REC-B1), δt is set by the fine level in both runs. Work per step scales with cells times an AMR
overhead of 1.2–1.5× [UNVERIFIED; FillPatch, reflux, composite MLMG]. The speed-up is therefore **≈ 2.7–3.3×**, which meets
the NFR-032 target (≤ 50% of uniform-fine time) if accuracy is maintained.

Waste of B-1 relative to subcycling: coarse cells advance with the fine δt. The wasted fraction of cell-updates is
`N_c(1 − 1/r)/(N_c + N_f)`. For the example (N_c = N_f, r = 2) that is **25%**, so subcycling would give ≤ 1.33× more.
Subcycling pays off only when coarse cells dominate (N_c ≫ N_f) or with more levels. This quantifies [OPEN Q2] and
roadmap Phase 6's "measured speed-up vs global dt" (`docs/roadmap.md`).

---

## 5. Parallel scaling notes [ESTIMATE / UNVERIFIED]

* **FFT (FDS):** each mesh is solved independently (perfectly parallel across meshes), but global coupling comes only from
  the interface iteration, and its iteration count grows with the number of meshes in a chain. That is the motivation for the
  tunnel preconditioner (00 §7).
* **MLMG:** global coupling on coarse MG levels. Agglomeration/consolidation (LPInfo) keeps the coarse levels efficient. On
  one 8-core node, communication is intra-node and cheap. The limiter is **memory bandwidth**: the measured single-thread
  stream is 24 GB/s. Node bandwidth was not measured, so 8-core speed-up for a bandwidth-bound V-cycle is expected to be
  4–6×, not 8× [UNVERIFIED].
* **UGLMAT-HYPRE:** FireX gathers rows to "resource-set" masters (00 §11.3). With `FDS_RANKS_PER_GPU` unset each rank is its
  own master, so CPU scaling is as on master.

### 5.1 Concrete expectation on the development machine (8 cores, CPU) [ESTIMATE]
| Case | Cells | MLMG per solve, 1 core | 8 cores (÷4–6) | Pressure per step (2 solves × ~1.5 baroclinic iterations) |
|---|---|---|---|---|
| `ns2d_16_*` (2-D, ≤ 32²) | ≤ 1 k | < 1 ms | — | negligible |
| `pressure_iteration3d_default` (8 × 16³) | 32 k | 3–8 ms | ~1–2 ms | ~5–25 ms |
| 64³ uniform | 262 k | 25–65 ms | 5–15 ms | 15–45 ms |
| 128³ uniform | 2.1 M | 190–520 ms | 35–130 ms | 100–400 ms. Memory ≈ 2.1 M × 110 B ≈ 0.23 GB for the solver, fits; whole-FDS state likely does not (§3) |

---

## 6. Measurement plan (to replace the estimates above)

1. Install an MPI runtime (or build FDS FireX against one) and run `Timing_Benchmarks/openmp_test64a.fds` (NFR-030 anchor) and
   `Pressure_Solver/pressure_iteration3d_default.fds` with `SOLVER='FFT'`, `'GLMAT'` and `'UGLMAT'`. Record `CHID_cpu.csv`
   (`PRES`, `PACK`, Total) and the pressure-iteration count (`VELOCITY_ERROR_FILE=T`, format at `main.f90:1714-1719`). This gives
   `f_pres` and the per-solve FFT/HYPRE baselines.
2. AMReX MLMG stand-alone: a 64³ and 128³ Poisson problem with 1 and 8 ranks. Record V-cycles and time for relative tolerances
   1e-10 and 1e-12, and ns/cell per V-cycle, to calibrate §2.1–2.3.
3. Two-level static case (`ns2d_16_int_1to2_refinement` extended to 3-D) to measure the composite overhead (§2.2) and memory
   (peak RSS, `/usr/bin/time -v` if installed).
4. Keep every run ≤ ~1 GB RSS on the development machine. Watch `free -m` before each run.

---

## 7. AMReX distributed FFT Poisson solvers as the uniform level-0 solver (replacing `pois.f90`)

**Question (AMR Chief Architect).** Can `amrex::FFT::Poisson`, `PoissonHybrid` or `PoissonOpenBC`
(`(local AMReX checkout)/Src/FFT/AMReX_FFT_Poisson.H`, AMReX `99ddfda`) replace Crayfishpak (`pois.f90`) for uniform level-0 runs
inside the AMReX driver? The goal is a single pressure path (one projection wrapper) with the linear solver chosen per run:
FFT on a uniform level 0, composite MLMG when AMR is active.

**Short answer.**
* `FFT::Poisson` can do it when level 0 is uniform (no stretching) and each of the six domain faces has a single BC type.
* `PoissonHybrid` adds stretching in z only.
* `PoissonOpenBC` does not apply: it is a free-space solver, not an FDS OPEN boundary.
* All three are C++ only.

Citations in this section are `AMReX_FFT_Poisson.H:<line>` unless another file is named. FDS citations are against FireX `36975d765f`.

### 7.0 Assumptions for this section
* **F1. Code reading only.** I read the headers and did not compile or run anything: the development machine is busy with an AMReX build and has no MPI runtime (P2). Every AMReX behaviour below comes from the source at `99ddfda`.
* **F2. "Uniform level 0"** means one cell size per direction over the whole domain box. AMReX `Geometry` imposes this anyway.
* **F3. Network and copy bandwidths are assumed, not measured.**
  * Inter-node: a 100 Gb/s NIC per node shared by 32–64 ranks, i.e. 0.2–0.4 GB/s per rank.
  * Intra-node: 5–10 GB/s effective copy bandwidth per rank, taken as ¼–½ of the measured 24 GB/s single-thread triad (§1.2) to allow for pack, MPI shared-memory copy and unpack.
* **F4. FDS behaviour** is as documented in 00; I re-read the cited lines for this section.

### 7.1 What the three classes provide (verified in source)
* **`Boundary` enum.** `AMREX_ENUM( Boundary, periodic, even, odd );` (`AMReX_FFT_Helper.H:59`).
  * Each direction takes a (low, high) pair (`:41-42`). Periodic must be set on both sides (`AMReX_FFT_R2X.H:239-244`).
  * The ghost fill `fill_physbc` (declared `:19`, defined `:863-919`) sets ghost = +interior for `even` and ghost = −interior for `odd` (`:910-914`). So `even` is zero normal gradient (Neumann) and `odd` is zero value at the face (Dirichlet).
  * The spectral offsets 0, ½ and 1 (`:328-342`) are those of cell-centred DCT/DST transforms with the boundary on the cell face, which matches FDS.
* **`Poisson`.**
  * Uses R2C when all directions are periodic and R2X otherwise (`:48-58`).
  * Solves ∇·∇φ = rhs (`:81-91`) using the eigenvalues of the 7-point stencil, `dxfac·(cos a − 1)` per direction (`:317-357`).
  * One cell size per direction; δx ≠ δy ≠ δz is allowed (`:317-320`). Cell-centred data only (`:304`).
* **`PoissonHybrid`.**
  * FFT in x and y, tridiagonal solve in z (`:156-161`; `solve_z` at `:636-847`).
  * x and y must each be periodic on both sides or on neither (`:185-194`). z must not be periodic (`:195-196`).
  * Non-uniform dz is supported through `solve(soln, rhs, dz)` (`:549-579`), with coefficients `2/(dz_k(dz_k+dz_{k±1}))` (`TriA`/`TriC`, `:453-474`). Generic tridiagonal functors are also accepted (`:255-256`).
* **`PoissonOpenBC`** (3-D only, `:101-154`, `:373-431`).
  * Convolves with the free-space Green's function (a 1/r kernel with 2×2×2 Gauss quadrature, `:396-422`) on a doubled domain (`AMReX_FFT_OpenBCSolver.H:201-231`).
  * This models an isolated system with potential → 0 at infinity. The FDS OPEN condition is a Dirichlet value on the boundary face (`pres.f90:145-224`). **Not applicable to FDS**, and it would cost about 8× the transform volume.

### 7.2 Q1: are the BCs homogeneous only, and how are FDS values lifted?
**Yes, homogeneous only.**
* No constructor or `solve` takes boundary values (`:41-42`, `:91`, `:221-237`).
* After every solve, `fill_physbc` overwrites the ghost cells with homogeneous images (`:363`, `:630`, `:910-914`).

**FDS already lifts its boundary values into the RHS inside Crayfishpak.** `H3CZSS` adds them to the RHS of the first and last cell (`pois.f90:240-270`):
* Dirichlet: `F(1,J,K) −= 2·BDXS·SAVE(IA)` (`pois.f90:243`).
* Neumann: `F(1,J,K) += SAVE(IA)·DX·BDXS` (`pois.f90:251`); high side at `259` and `267`.

Here `SAVE(IA)` is taken to be the 1/δx² coefficient (an assumption from the form; not traced through `S3CFIS`). The AMReX driver would do the same, cell by cell, before calling `solve`:

| FDS BC (from `BXS…BZF`) | RHS lift (low x face, cell 1) | Ghost value after the solve (replacing AMReX's homogeneous fill) |
|---|---|---|
| Dirichlet value g (OPEN; solid cell on a Dirichlet face) | rhs₁ −= 2g/δx² | H₀ = 2g − H₁ (as `pres.f90:454`) |
| Neumann gradient g | rhs₁ += g/δx, with the sign flipped on the high face as in `pois.f90:267` | H₀ = H₁ − δx·g (form of `pres.f90:452`) |

* **Neumann values need no lift in the projection form.** In 01 §F the prescribed normal velocity enters the trial velocity, so the Neumann lift is zero. Only Dirichlet (OPEN) values need lifting.
* **Both steps are cell-local.** Estimated at ~30 lines of kernel code, from counting 6 faces × 2 statements plus loops.
* **Per-cell values are fully supported; only the BC *type* is per face.**

### 7.3 Q2: coverage of every FDS boundary case
**How FDS assigns BCs.**
* **FFT:** the BC type is per mesh and per face.
  * The type is held in `LBC/MBC/NBC`: 0 periodic, 1–4 the Dirichlet/Neumann pairs (`init.f90:2476-2491`, `cons.f90:97-101`). The default is Neumann–Neumann (`init.f90:2493-2495`).
  * **Any OPEN vent on a face makes the whole mesh face Dirichlet** (`init.f90:2497-2524`). Every mesh–mesh (interpolated) face is also Dirichlet (`init.f90:2528-2556`).
  * Each wall cell's `PRESSURE_BC_TYPE` is then copied from its face type (`init.f90:2598-2629`).
* **Per-cell values carry the physics.**
  * Neumann cells: `∂H/∂n = −F_n ∓ DUNDT` (`pres.f90:75-93`).
  * Solid or mirror cells on a Dirichlet face: `BXS = ½(H_ghost + H₁) + WALL_WORK1` (`pres.f90:97-113`).
  * `WALL_WORK1` comes from the velocity-error iteration (`pres.f90:802-1076`), whose target is `U_NORMAL` (`1044-1050`).
* **GLMAT** instead applies homogeneous Neumann per cell at external solid or mirror cells, even on a Dirichlet face (`pres.f90:3612-3632`, comment "Set Homogeneous Neumann in external SOLID_BOUNDARY").
* **ULMAT** gives up the FFT when a face's BC type is not uniform (`pres.f90:1221-1265`; the test is at `1258`).

| FDS case | FDS FFT today | `FFT::Poisson` over level 0 | Verdict |
|---|---|---|---|
| Periodic | One mesh: `FISHPAK_BC=0` (`init.f90:2558-2560`), allowed only for a single mesh (`read.f90:10167-10168`). Several meshes: PERIODIC vents become INTERPOLATED boundaries (`init.f90:3151`, `3251-3254`), i.e. Dirichlet plus iteration | `periodic` pair, exact | Exact; **better** than FDS with several meshes |
| Solid wall / MIRROR over the whole face | Neumann (`LBC`=3) | `even` plus lift (zero lift in projection form) | Exact |
| OPEN over the whole face | Dirichlet (`init.f90:2497-2524`), values from `pres.f90:145-224` | `odd` plus per-cell lift (§7.2) | Exact |
| Prescribed velocity (supply/exhaust) on a Neumann face | Per-cell Neumann values (`pres.f90:75-93`) | `even`, with u_n in the trial velocity | Exact |
| Face partly OPEN, partly closed | Whole *mesh* face Dirichlet; closed cells emulated by Dirichlet plus iteration (`pres.f90:97-113`, `802-1076`) | Whole *domain* face `odd`, with the same emulation and iteration | See note below |
| Mesh–mesh interface | Dirichlet plus iteration (`pres.f90:115-143`) | Not a boundary: one global solve | Error and iteration removed |
| Domain that is not a box (union of meshes) | Allowed | The level-0 domain is a box, so uncovered parts must be blocked (OBST) | Costs IBM iterations |

**Mixed faces.** `FFT::Poisson` gives parity with FDS FFT, but the emulated region grows from one mesh face to the whole domain face, so more iterations are possible (a qualitative estimate). Exact alternatives:
* GLMAT-style per-cell Neumann cannot be expressed in `FFT::Poisson`.
* A capacitance-matrix correction: M extra unknowns (M = closed cells on Dirichlet faces), with CG on the Schur complement at one FFT solve per CG iteration. This is an estimate; AMReX does not provide it.
* Use MLMG for such cases. Per-cell mixed BCs in MLMG are themselves unverified (01 §F).

### 7.4 Q3: stretched grids
* **FDS** allows `TRNX/TRNY/TRNZ` stretching in at most two directions per mesh (`init.f90:2356-2366`, `ERROR(425)` at `2366`), using Crayfishpak `H3CZIS/H3CSIS`.
* **AMReX:** `Poisson` is uniform only. `PoissonHybrid` handles stretched z only, and z must be non-periodic (`:195-196`, `:549-579`).
  * Its z coefficients `2/(dz_k(dz_k+dz_{k±1}))` are the standard cell-centred finite-volume form. This equals `RDZ·RDZN` with `DZN = ½(dz_k+dz_{k+1})`; that is the FDS form I assume, as I did not re-read FDS's stretched Laplacian for this section.
* **Hybrid only helps if the whole driver carries a stretched z,** but `Geometry`, and therefore MLMG, FillPatch and `MacProjector`, is uniform.
* **Consequence:** stretched runs stay outside the AMReX driver (01/README Q1). Hybrid is worth revisiting only if z-stretched uniform runs become a requirement.

### 7.5 Q4: is one global FFT equivalent to GLMAT, and which iterations remain?
**Confirmed, with qualifications.**
* **What GLMAT is.**
  * GLMAT numbers every cell, solids included (`pres.f90:5739-5748`, "Classic IBM"), in one global system (`pres.f90:3698`, `N_ZONE_GLOBMAT = 0`).
  * Mesh interfaces are ordinary internal rows, so no interface BC is needed (`pres.f90:3613-3616`: "that's the whole point of a global solve").
  * `ITERATE_PRESSURE` is left true only for diagnostics (`pres.f90:3712-3714`).
* **Why a global FFT matches it.** One FFT over a uniform level 0 has the same unknowns (all cells), the same constant-coefficient 7-point operator and no interfaces.
  * It therefore reproduces GLMAT's solution **with or without obstructions**, because both treat solids by IBM forcing on all cells. The condition is that each external face has a single BC type.
* **Remaining differences:**
  1. On mixed faces, GLMAT uses per-cell Neumann (`pres.f90:3618-3632`); the FFT uses Dirichlet emulation plus iteration.
  2. `NO_FLUX` uses `DHFCT = 0` at external walls for GLMAT and 1 for FFT (`velo.f90:1483-1487`), so the forced F differs in those cells.
  3. GLMAT runs PCG to 1e-12 (01 §I.1); the FFT is a direct solve.
* **Removed:**
  * the mesh-interface discretisation error (the Dirichlet interface, `pres.f90:115-143`);
  * the interface part of the velocity-error iteration (INTERPOLATED cells in `COMPUTE_VELOCITY_ERROR`, `pres.f90:802-1076`);
  * the reason for the tunnel preconditioner (`pres.f90:505-697`).
* **Still present:**
  * the solid (IBM) iteration at obstructions (`velo.f90:1404-1540`, `pres.f90:802-1076`);
  * the mixed-face emulation iteration (§7.3);
  * the baroclinic iteration (`main.f90:1609-1613`; 01 §D).

  This is the same iteration logic as the MLMG path with E-1 (01 §D), so FFT and MLMG can share one iteration driver.

### 7.6 Q5: cost [ESTIMATE unless marked]
**Communication read from the code.**
* R2X starts from x-pencils (`AMReX_FFT_R2X.H:263`), redistributes x→y and y→z (`AMReX_FFT_R2X.H:423-457`; `ParallelCopy` at `AMReX_FFT_R2X.H:789-811`), and reverses both on the way back (`AMReX_FFT_R2X.H:960-962`).
* It also copies in from the user's MultiFab (`AMReX_FFT_R2X.H:780`) and copies out (`AMReX_FFT_R2X.H:759`, `866-868`, `923-925`).
* So one solve does **4 transposes plus 2 layout copies = up to 6 full-field redistributions**. The two layout copies become local if the level-0 boxes are chosen to match the x-pencils.
* Each redistribution moves ≈ 8 B/cell: real data in the r2r path, or complex data over half the cells.

| Solver | Compute (ns/cell/solve, per core) | Communication (ns/cell/solve) | Total per solve | Derivation |
|---|---|---|---|---|
| FDS per-mesh FFT | 14–15 **[MEASURED proxy, §1.1]** | none in the solve; `MESH_EXCHANGE(5)` plus velocity error each iteration | ≈ 15 × k, with k = 1–10 pressure iterations | §1 |
| Global FFT, one node | ≈ 15 (same transforms; FFTW r2r on CPU) | 6 × (8 B × ~3 copy passes at 5–10 GB/s) ≈ 6 × (2.4–4.8) ≈ 15–30 | **≈ 30–45** | F3 |
| Global FFT, many nodes | ≈ 15 | 6 × 8 B / (0.2–0.4 GB/s) ≈ 6 × (20–40) ≈ 120–240 | **≈ 135–255** (latency excluded: ~√P partners per transpose, because the pencils are decomposed in 2 dimensions, `AMReX_FFT_R2X.H:263`) | F3 |
| Composite MLMG | 90–250 (§2.3) | halos plus coarse levels (small) | ≥ 90–250 | §2.3 |

**Consequences.**
* **One node.**
  * A global FFT solve costs ≈ 2–3× one per-mesh FFT solve, but it replaces the k solves of the FDS interface iteration. It breaks even at k ≈ 2–3, and it removes the interface error.
  * It is ≈ 2–8× cheaper than MLMG per solve.
  * In the §2.4 formula, r ≈ (2–3)/k, so **NFR-030 is met for uniform runs whenever k ≥ 2–3 today**. For k = 1 (a single mesh, no iterations) the slowdown is 1 + f_pres(r − 1) with r ≈ 2–3, so f_pres ≤ 0.12–0.25 would be needed.
* **Many nodes.** The all-to-all volume makes the global FFT about as expensive as MLMG, and it scales worse (all-to-all per solve versus halo exchange). For large multi-node uniform runs, FDS's per-mesh FFT may remain faster per step if k is small.
* **Memory.** R2X aliases two work MultiFabs (`AMReX_FFT_R2X.H:405-417`), about 2–3 doubles per cell ≈ 16–24 B/cell (plus the spectral MultiFab for Hybrid, `:489-537`). That is below MLMG's ≈ 64 B/cell (§3).
* **Accuracy.** A direct solve leaves a residual at round-off, so it meets the 01 §I.2 ε_H test with margin (estimate: relative residual ~1e-14–1e-13 for N ≤ 128).
* **Gauge and compatibility (inference from code).**
  * For all-Neumann or all-periodic problems, `Poisson` does not divide the k² = 0 mode; it only rescales it (`:353-356`). The incompatible part of the RHS is dropped and reappears as the solution's constant. Apply the 01 §C.5 gauge afterwards.
  * `PoissonHybrid` instead doubles one diagonal entry of the singular mode (`:699-702`, `:748-750`, `:811-813`). The caller must remove the RHS mean first (01 §C.4).

### 7.7 Q6: C++ only (ADR-001)
* **No Fortran binding.**
  * `Src/F_Interfaces` contains only `AmrCore`, `Base`, `LinearSolvers`, `Octree` and `Particle`, and a case-insensitive search for "fft" there returns nothing.
  * The classes are C++ templates with C++20 `requires` clauses (`:43`, `:119-120`).
* **A Fortran driver would need a C++ shim.** It would be an `extern "C"` function that takes the `MultiFab` pointer held by `amrex_multifab`, builds the `Geometry` and BC pairs, and calls `solve`. Estimated at ~100–200 lines, including the lift from §7.2.
* **This is one more C++ item** alongside `MacProjector`, the overset mask, EB and nodal operators (01 §G.1). It is an argument for ADR-001 **Option A**, but not a blocker for Option B.

### 7.8 Q7: GPU path
* **Yes.** `FFT::Poisson` uses cuFFT, rocFFT or oneMKL DFT backends (`AMReX_FFT_Helper.H:23-42`; FFTW on CPU).
  * Kernels run as device lambdas through `ParallelFor` (`:344-357`, `:686-696`). `fill_physbc` has a GPU path (`:894-896`), and there is a device-vector dz overload (`:237`, `:549-559`).
* **Caveats.**
  * On GPU, the r2r transforms (Neumann/Dirichlet directions) are emulated with complex FFTs of length 2n–4n (`AMReX_FFT_Helper.H:804-822`). Those directions therefore cost ~2–4× a periodic transform (derived from the lengths).
  * Hybrid's GPU tridiagonal path is marked "TODO … optimize" (`:715-718`) and runs one thread per (i, j) column with strided z access.
* **FireX today:** its GPU pressure path is HYPRE only (00 §11.2).
* **Not relevant on the development machine,** which has no GPU (D-004).

### 7.9 Recommendation
**[REC-P7]** Use `FFT::Poisson` as the uniform level-0 linear solver inside the AMReX driver instead of porting `pois.f90`. It can be called through AMReX-Hydro `FFTMacProjector` (01 §A.4) or through a thin wrapper that adds the source term.
* **One projection path.** The solver is selected by configuration:
  * FFT when `finest_level = 0`, the grid is unstretched, and each domain face has a single BC type;
  * composite MLMG otherwise.
* **Keep FDS's per-mesh FFT** only in the legacy Fortran code, as the NFR-030 baseline, until f_pres and k are measured.
* **Do not use `PoissonOpenBC`.** Consider `PoissonHybrid` only if z-stretched runs become a requirement.
* **Work needed:**
  * the Dirichlet lift and ghost re-fill (§7.2);
  * zero-mode and gauge handling (§7.6);
  * the C++ shim if ADR-001 picks Option B (§7.7);
  * mixed faces left on the existing iteration (§7.3) until a capacitance or MLMG alternative is chosen.
* **Measurement to add once MPI is available:**
  * the k distribution from `VELOCITY_ERROR_FILE` for the §6 cases;
  * `FFT::Poisson` 64³/128³ on 1 and 8 ranks, recording the time split between transforms and `ParallelCopy`, to calibrate the table above.

---

## Appendix A. Benchmark script used for §1.1–1.2 (verbatim)

```python
import numpy as np, time, os
def t(f, n=5):
    f(); best=1e9
    for _ in range(n):
        s=time.perf_counter(); f(); best=min(best,time.perf_counter()-s)
    return best
print("numpy", np.__version__, "threads env OMP", os.environ.get("OMP_NUM_THREADS"))
for N in (64,128):
    a=np.random.rand(N,N,N)
    tf=t(lambda: np.fft.irfftn(np.fft.rfftn(a), s=a.shape))
    print(f"FFT fwd+inv {N}^3: {tf*1e3:.2f} ms  -> {tf/N**3*1e9:.1f} ns/cell")
    def lap():
        r=np.empty_like(a)
        r[1:-1,1:-1,1:-1]=(a[2:,1:-1,1:-1]+a[:-2,1:-1,1:-1]+a[1:-1,2:,1:-1]+a[1:-1,:-2,1:-1]+a[1:-1,1:-1,2:]+a[1:-1,1:-1,:-2]-6*a[1:-1,1:-1,1:-1])
        return r
    tl=t(lap)
    print(f"7-pt stencil apply {N}^3 (numpy): {tl*1e3:.2f} ms -> {tl/N**3*1e9:.1f} ns/cell")
    del a
M=20_000_000
b=np.random.rand(M); c=np.random.rand(M); x=np.empty(M)
def triad(): np.multiply(c,3.0,out=x); np.add(x,b,out=x)
tt=t(triad)
# bytes: mult reads c writes x (16B), add reads x,b writes x (24B) =40B/elem
print(f"triad-like: {tt*1e3:.1f} ms -> {40*M/tt/1e9:.1f} GB/s effective (single thread numpy)")
```
