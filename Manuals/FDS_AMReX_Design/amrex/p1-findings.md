# P1 findings: unmodified FDS mass kernels on a static two-level AMReX hierarchy

Scope: the goal was to run FDS `MASS_FINITE_DIFFERENCES` and `DENSITY` (mass.f90, FireX `36975d765f`) box by box on a static two-level AMReX `AmrCore` hierarchy (AMReX `99ddfda`), without editing any FDS file. The code is in `prototypes/p1_mass_shim/` (see its README for build and run commands). All runs used GNU 14.2 and OpenMPI 5.0.7 on the shared machine, with at most 2 ranks and 64^3 cells. File:line citations refer to `src/Source` at `36975d765f` unless a prototype file is named.

## 1. Summary

**What works:**
- Both kernels run unmodified on every box of a two-level hierarchy, on 1 or 2 ranks.
- Two shim variants give bitwise-identical results: aliasing FAB memory, and copying in and out.
- Mass is conserved to round-off on one level: relative change about 1e-15 after 8–128 steps.
- The -O0 debug build and the -O3 release build give bitwise-identical results.
- With FP traps on (invalid, zero, overflow), no trap fires in any single-level or two-level run.

**Decomposition check (requested):** the single-level result is bitwise identical for max_grid_size 32 vs 8 and for 1 vs 2 ranks, as long as `CHECK_MASS_DENSITY` does not clip species. That covers the smooth blob with uniform or sheared velocity, and the FDS initial field with the GODUNOV limiter. It is **not** bitwise when species clipping happens, which is the case for the SUPERBEE limiter on the sharp FDS initial field. The cause is a box-wide pass inside FDS (§6.2). This is a property of FDS itself, not of the shim. **Resolved by D-031 (§13):** with the gather-form clip (`p1.clip=gather`) this case is bitwise identical across layouts and rank counts, and bitwise identical to single-mesh FDS, also with 64 boxes on 2 ranks.

**Match with FDS: bitwise.**
- P1 on a single 32^3 box, started from an FDS restart, reproduces FDS restarted from the same file bitwise over 8 steps (RHO, ZZ and TMP), with both SUPERBEE and GODUNOV.
- With GODUNOV, P1 on 64 boxes over 2 ranks also matches the FDS single-mesh run bitwise.
- The ulp-level differences I saw at first came from FDS's end-of-run time-step rule (main.f90:719), not from the shim (§5).

**Two levels:**
- Without reflux, the tracer mass drifts by up to 4.5e-4 (relative) over 128 steps while the blob crosses the coarse/fine interface.
- The drift is entirely coarse/fine flux mismatch: the mismatch summed from the kernels' own face fluxes (`p1_box_flux`) equals the measured mass change to 10 significant digits (§9.1). A FluxRegister reflux would therefore remove it.
- There is no spike at the interface: errors in the coarse cells next to the fine patch lie between the single-level 32^3 and 64^3 errors.

**Cost (64^3, 1 rank, 8 steps):**
- The kernels take 168 ns per cell per stage with a single box. They take 65 % longer with 8^3 boxes, because FDS works on padded arrays.
- Shim overhead per call is 0.5–5.6 % of kernel time with alias and 15–31 % with copy.
- Building the side data once takes 0.06–0.17 s.

## 2. Which FDS files are compiled, and why

The FDS files are compiled in place and unmodified into `libfds_closure.a` (prototype `CMakeLists.txt:15-18`). They use the FDS release flags `-O3 -cpp -std=f2018 -frecursive -ffpe-summary=none -fall-intrinsics -fopenmp`, which are the same as the FDS CMake build's `flags.make`, and no defines.

The set is the transitive USE-closure of `mass.f90`: 15 files, 82,312 of about 181k FDS lines. The only undefined externals left are libc, libm, MPI, libgfortran and OpenMP.

| file | lines | needed because |
|---|---|---|
| prec.f90 | 60 | `PRECISION_PARAMETERS` (mass.f90:5) |
| cons.f90 | 1011 | `GLOBAL_CONSTANTS` (mass.f90:6): flags, limiter, species counts, RHOMIN/MAX, T_USED |
| type.f90, mesh.f90 | 2214, 918 | `MESH_POINTERS` / `MESHES` / `POINT_TO_MESH` (mass.f90:7) |
| imkl.f90, devc.f90, prop.f90 | 617, 209, 8935 | modules USEd by type/mesh/func (MKL interface stubs, device and property types) |
| func.f90 | 7221 | `COMP_FUNCTIONS` (CURRENT_TIME, :8), `MATH_FUNCTIONS` (GET_SCALAR_FACE_VALUE, :22), `PHYSICAL_FUNCTIONS` (:23, :367) |
| gsmv.f90, data.f90 | 657, 3262 | USEd by func and ccib |
| turb.f90 | 3835 | `MANUFACTURED_SOLUTIONS` for DENSITY (:368), used only when PERIODIC_TEST=7 |
| soot.f90 | 599 | `SOOT_ROUTINES` SETTLING_VELOCITY (:369), used only with depositing species |
| ccib.f90, geom.f90 | 24046, 27735 | `CC_SCALARS` (:370, SET_EXIMADVFLX_3D and ROTATED_CUBE_RHS_ZZ, used only with CC_IBM or PERIODIC_TEST 21-23). ccib pulls in all of geom. These two files are 63 % of the closure and most of the 2-minute build. |

**Implication:** mass.f90 cannot be compiled alone. Half of the closure is there only to satisfy USE statements for branches the prototype never takes. Nothing needs stubbing, though.

## 3. The POINT_TO_MESH blocker and the two workarounds

The original plan was to re-point the `MESH_POINTERS` module pointers at FAB data. **That cannot work.**
- Both kernels call `POINT_TO_MESH(NM)` themselves (mass.f90:42, :397).
- `POINT_TO_MESH` re-associates every pointer to the ALLOCATABLE components of `MESHES(NM)` (mesh.f90:505-916; components are declared ALLOCATABLE at mesh.f90:16-354, e.g. U at :18 and RHO at :39).
- mass.f90:477 and :639 read `ALLOCATED(MESHES(NM)%M_DOT_PPP)` directly.

So box data must arrive through `MESHES(NM)%X`. There is one MESHES entry per AMReX box: NM = level offset + box index + 1, and `NMESHES` is set to the total. I found two ways to do this without editing FDS, selected with `p1.shim_mode`:

- **W1, copy** (standard Fortran):
  - Make bounds-remapped views of the FAB data (`RHO_V(-1:,-1:,-1:) => F3` after `C_F_POINTER`, prototype `p1_shim.f90:263-270`).
  - Copy them into FDS-shaped allocatables in `MESHES(NM)` before the kernels, and copy the outputs back afterwards (`p1_shim.f90:272-275, 312-318`).
- **W2, alias** (the default; not standard Fortran):
  - `p1_alias.c` writes the descriptor of the unallocated ALLOCATABLE component so that it points at FAB memory with FDS lower bounds, and nulls it afterwards (`p1_shim.f90:277-284, 320-323`).
  - This works with gfortran 14 and gives bitwise-identical results to W1. Other compilers are not guaranteed to allow it.

Both keep the FDS lower bounds. The kernels depend on them: `FX_P(LBOUND(FX,1):,...)` remaps at mass.f90:81-83 and :212-214, and the fixed loop ranges run from -1 to IBP1+1.

**Mapping between FAB layout and FDS arrays** (MultiFab ghost counts chosen so that the FDS bounds come out exactly):
- **Cell-centred:** rho, rhos, tmp (1 comp) and zz, zzs (NS comps) use ng=2, which gives FDS `-1:IBP1+1`. FDS cell I is global cell lo+I-1.
- **Face-centred:** U, V, W are nodal in their own direction with ng=1. That gives exactly `U(-1:IBP1,0:JBP1,0:KBP1)`. FDS face I is global face lo+I.
- **FX/FY/FZ cannot be AMReX face MultiFabs.** Their bounds are `0:IBP1` in every direction (init.f90:611-613), so the normal direction has an asymmetric ghost. They stay as per-box scratch in MESHES.

## 4. Per-box pointers and side data

Everything below is set by `p1_box_setup` (once per box, static grid) or `p1_box_stage` (every call).

**Bound to FAB memory every call** (alias or copy):
- RHO, RHOS, TMP: `-1:IBP1+1` in each direction.
- ZZ, ZZS: `-1:IBP1+1` ×3, `1:NS`.
- U, V, W and US, VS, WS: US/VS/WS are bound to the same data as U/V/W, which assumes a frozen velocity (see §8).

**Per-box scratch in MESHES(NM)**, allocated once with FDS bounds:
- `WORK_PAD(-1:IBP1+1)^3`.
- `WORK_U/V/W` with U/V/W bounds.
- `FX/FY/FZ(0:IBP1,0:JBP1,0:KBP1, 0:NS)`. These must be zeroed like init.f90:611-613 does. The MW-correction loop (mass.f90:327-350) reads rows J=0 and K=0 of FX, which GET_SCALAR_FACE_VALUE never writes (range 1:JBAR at mass.f90:84). The values computed there are never used, but they must be finite for trap-clean runs.
- `DEL_RHO_D_DEL_Z` and `SWORK4 (0:IBP1)^3 × NS`. DEL_RHO_D_DEL_Z = 0 because P1 has no diffusion.
- `WORK4`, `WORK5`, `RSUM`, `PRESSURE_ZONE` = 0.
- `PBAR`, `PBAR_S(0:KBP1,0:N_ZONE)` = P_INF, and `D_PBAR_DT(_S)` = 0.

**Grid metrics:** `DX, RDX, RRN(0:IBP1)`, `R(0:IBAR)`, `DY, RDY`, `DZ, RDZ`.

**Box sizes:** IBAR/JBAR/KBAR and IBP1/JBP1/KBP1 (plus IBM1, …). Boxes must be at least 2 cells in every direction.

**Cells:** `CELL_INDEX(0:IBP1)^3` and `CELL(:)` for every cell including the ghost ring, with SOLID = F and EXTERIOR = T on the ring.

**Walls (steering point 1):** there is one EXTERNAL_WALL entry per box-face cell, NEXT = 2(JK+IK+IJ), in IOR order +1, -1, +2, -2, +3, -3 (`p1_shim.f90:184-232`). Each entry has:
- `WALL%BOUNDARY_TYPE = INTERPOLATED_BOUNDARY`, `BC_INDEX = IW`, `B1_INDEX = 1`;
- `BOUNDARY_COORD` II/JJ/KK, IIG/JJG/KKG, II2/JJ2/KK2 and IOR, following FDS conventions;
- `EXTERNAL_WALL%BOUNDARY_TYPE_PREVIOUS = INTERPOLATED` and `NOM = 0`;
- **`CELL(gas cell)%WALL_INDEX(-IOR) = IW`**, the CHECK_MASS_DENSITY barrier (mass.f90:831-836, 907-912).

There is a single dummy `BOUNDARY_PROP1(1)` with ZZ_F allocated. It is pointer-assigned at mass.f90:97 but never read for interpolated walls.

**UVW_SAVE (steering point 2):** `UVW_SAVE(IW)` = the box's own face-normal velocity at that face, refreshed on every call (`p1_shim.f90:287-300`). That makes DENSITY's restore (mass.f90:429-434, 601-606) a no-op. MATCH_VELOCITY (velo.f90:2669, 2732-2774) is not emulated.

**Left unallocated:** `M_DOT_PPP` and `D_SOURCE`. The `ALLOCATED` test at mass.f90:477 and :639 then skips the source term.

**Global state set once** (`p1_init_globals`):
- N_TRACKED_SPECIES = N_TOTAL_SCALARS = 2, N_PASSIVE_SCALARS = 0.
- N_ZONE = 0. It has **no default** (cons.f90:695), and the PBAR loops use it.
- GRAVITATIONAL_SETTLING and THERMOPHORETIC_SETTLING forced to F. **Both default to .TRUE.** (cons.f90:231, 233), which would call SETTLING_VELOCITY if any species were depositing.
- PERIODIC_TEST = 0; SOLID_PHASE_ONLY, CC_IBM and STORE_SPECIES_FLUX = F.
- FLUX_LIMITER_MW_CORRECTION = T with N_LOWER_SCALARS = 0 (init.f90:605-609).
- I_FLUX_LIMITER.
- RHOMIN/RHOMAX = 0.01 / 100, the FDS `&CLIP` defaults (read.f90:10381-10382). No RHO clipping happened in any run.
- ICYC = 2, because DENSITY returns immediately for ICYC ≤ 1 (mass.f90:390).
- T_USED allocated.
- SPECIES_MIXTURE(:)%MW (32, 64), DEPOSITING = F, and MWR_Z = 1/MW (read.f90:3535).
- NMESHES, and `MESHES(1:nbox_total)` allocated on every rank. Only the boxes a rank owns have components.

**Per-call flags:** PREDICTOR, CORRECTOR, FIRST_PASS = T.

**Memory:**
- CELL_TYPE is 384 bytes per cell, and each wall-cell entry (WALL + EXTERNAL_WALL + BOUNDARY_COORD) is 540 bytes. Both were measured with `storage_size`.
- Side data therefore comes to about **1000 B per cell for 16^3 boxes and 760 B per cell for 32^3 boxes**, roughly 4 MB per 16^3 box. The FAB state (7 cell arrays with ng=2 plus 3 face arrays) is about 150 B per cell.
- The CELL array alone is more than half of the side data.

## 5. Match with FDS (single 32^3 box vs FDS single mesh)

**Reference setup** (`fds_ref/make_inputs.sh`):
- 32^3 fully periodic unit cube.
- `FREEZE_VELOCITY`, DNS mode, `DT=1/64` with `LOCK_TIME_STEP`, DIFFUSIVITY = 0.
- GASA (MW 32, background) and GASB (MW 64), with three nested boxes of GASB (0.2 / 0.4 / 0.6).
- Constant ramps set the velocity. After FDS's initial projection it is not uniform (U from 0.981 to 1.019), but it is divergence-free to 2e-14 and exactly 1 / 0.5 / 0.25 on the boundary faces.

**Procedure** (`fds_ref/check_fds.sh`):
1. FDS runs 2 steps and writes a restart. DENSITY is skipped at ICYC = 1.
2. FDS is **restarted** from that file and runs 8 more steps.
3. P1 starts from the same restart (FDS RHO, ZZ and full U/V/W) and runs 8 steps with FDS's T/DT sequence.

**Results after 8 steps:**

| case | result vs FDS |
|---|---|
| SUPERBEE, 1 box, alias | **bitwise** (RHO, ZZ, TMP) |
| SUPERBEE, 1 box, copy | **bitwise** |
| GODUNOV, 1 box | **bitwise** |
| GODUNOV, 64 boxes (max_grid_size 8), 2 ranks | **bitwise**, even though FDS ran a single mesh |
| SUPERBEE, 64 boxes, 2 ranks | differs by up to 5e-5, from box-local species clipping (§6.2); **bitwise with the D-031 gather clip (§13)** |

**What the "T/DT sequence" means.** At the start of every step FDS applies `IF ((T+DT+DT_END_FILL)>T_END) DT = MAX(T_END-T+TWENTY_EPSILON_EB,DT_END_MINIMUM)` (main.f90:719). DT_END_FILL is 1e-6 (cons.f90:429), so this fires on the **last step of every run**:
- The last step of the 2-step run used DT = 1/64 + 20ε = 0.01562500000000444 and ended at T = 1/32 + 20ε.
- The restart stores T and DT, and LOCK_TIME_STEP keeps that DT. A restarted run therefore steps with 0.01562500000000444, and its own last step is shortened again (to 0.015624999999995559 for a 2-step continuation).
- P1 reproduces this with `p1.t0`, `p1.dt` and `p1.t_end` (`p1_driver.cpp`, same formula).

**How I found it:**
1. An independent numpy re-implementation of the two kernels for the GODUNOV case, following mass.f90's order of operations (`fds_ref/emulate_godunov.py`), matched P1 bitwise.
2. That emulation matched FDS for one step after a restart, where the DT rule happens to give exactly 1/64 again, but not for two steps.
3. With the DT sequence above, it matched FDS for two steps.

**Earlier misleading comparisons:**
- Comparing against a *from-scratch* FDS run (whose step 2 used exactly 1/64, unlike the restart file's) gave differences of up to 2.5e-14 (SUPERBEE) and 1.3e-15 (GODUNOV). That was the whole of the earlier "unexplained" mismatch.
- The first reference used DIFFUSIVITY = 1e-30. That gives DEL_RHO_D_DEL_Z ≈ 1e-28, which also changes results at the ulp level wherever the flux divergence is small, so the reference now uses DIFFUSIVITY = 0. FDS honours `D_USER >= 0` (data.f90:2595, 2757). The old runs are in `fds_ref/old_d1e-30/`.
- FDS's same-resolution ghost value is `MAX(0,MIN(1,(ρZ)/ρ))` (wall.f90:333-334, 359). It never mattered here, because the periodic faces lie in pure background. P1 can emulate it (`p1.fds_ghost_rounding=1`) but does not by default (§6).

## 6. Decomposition independence (steering point 3)

`check_decomp.sh` runs a single level on 32^3 for 8 steps with max_grid_size 32 on 1 rank as the reference. It compares that against max_grid_size 8 on 1 rank (64 boxes), max_grid_size 32 on 2 ranks (2 boxes) and max_grid_size 8 on 2 ranks.

| case | 32/1 vs 8/1 | vs 32/2 | vs 8/2 |
|---|---|---|---|
| blob, uniform velocity, SUPERBEE, alias | bitwise | bitwise | bitwise |
| blob, sheared velocity (`vel_type=1`), SUPERBEE | bitwise | bitwise | bitwise |
| blob, SUPERBEE, copy mode | bitwise | bitwise | bitwise |
| FDS initial field + FDS velocity, GODUNOV | bitwise | bitwise | bitwise |
| FDS initial field + FDS velocity, SUPERBEE | **differs** (5e-5 after 8 steps) | differs (4e-8) | differs |

In addition:
- The GODUNOV case on 64 boxes over 2 ranks is also bitwise identical to FDS's own single-mesh result (§5).
- The two-level runs are bitwise identical between max_grid_size 16 and 8, on both levels.
- The alias-mode result at max_grid_size 32 on 1 rank matches copy mode at 8/1 and 32/2 bitwise.

### 6.1 Why the own-face UVW_SAVE policy works

At a box face, both neighbouring boxes hold the same face velocity in the nodal MultiFab, and each sets UVW_SAVE to that value. So the flux through the face is computed from identical inputs on both sides, and the face-value stencil (GET_SCALAR_FACE_VALUE, func.f90:1330, 2 cells upwind and 1 downwind) only sees exact ghost copies.

### 6.2 Why SUPERBEE on sharp data breaks it: CHECK_MASS_DENSITY is box-global

The first difference shows up at step 3, and it covers every blob cell of one box, not just cells near the box face. The cause is in `CHECK_MASS_DENSITY` (mass.f90:775-963):
- If **any** cell in the box has ρZ < 0 or ρZ > ρ, then `CLIP_RHO_ZZ(N)` becomes true (mass.f90:888-889).
- In that case the early return at mass.f90:943 is skipped, and the absorb/renormalise loop at mass.f90:947-960 runs over **every cell of the box**. `RHO_ZZ(N) = RHO_ZZ(N) + RHOP - SUM(RHO_ZZ)` changes values by an ulp even in cells that were never clipped.
- Which cells get that treatment therefore depends on which box the clipped cell sits in. SUPERBEE then amplifies those ulps at the sharp steps, reaching 5e-5 after 8 steps.
- On top of this, the clipping redistribution itself stops at box faces (the WALL_INDEX barrier, mass.f90:831-836, 907-912).

The evidence:
- GODUNOV (monotone at this CFL, so no clipping) on the same data is bitwise.
- CENTRAL (heavy clipping) differs, and also changes tracer mass by 1e-3, because clipping is not conservative.

FDS multi-mesh runs have exactly the same property, since this code is per mesh. **Bitwise decomposition independence is therefore only achievable when no species clipping happens**, unless the driver changes how CHECK_MASS_DENSITY is scoped. That would need an FDS change. (D-031, §13, does exactly that in a prototype copy, without editing FDS, and is bitwise equal to single-mesh FDS.)

Emulating FDS's ghost rounding (wall.f90:333-334) would also break decomposition independence, at the ulp level (tested with the blob: max_grid_size 32 vs 8 differ by 6.7e-16). The driver should use plain copies, as P1 does by default.

## 7. Floating-point traps and ghost coverage

**Debug build:** FDS debug flags (`-O0 -finit-real=snan -fcheck=all -ffpe-trap=... -fbounds-check`), with traps switched on through AMReX (`amrex.fpe_trap_invalid/zero/overflow=1`). The main program is C++, so gfortran's `-ffpe-trap` alone never enables them.

**Clean runs** (no trap, no bounds error):
- single level, alias and copy;
- FDS initial field on 2 ranks;
- two-level blob on 2 ranks, alias and copy.

**Negative control `p1.poison_corners=1`** (sNaN in edge and corner ghosts):
- With traps on, it traps at **mass.f90:73**, the loop `RHO_Z_P = RHOP*ZZP` over the full `-1:IBP1+1` range, called from mass.f90:81.
- With traps off, the result is **bitwise identical** to the unpoisoned run.
- So edge and corner ghosts are read, and must hold finite values, but they do not affect the result. Only the 2 face-adjacent ghost layers matter.

**Negative control `p1.zero_fx=0`:** no trap. The fresh heap pages happened to be zero, so this does not show that zeroing FX is unnecessary. FDS zeroes FX itself (init.f90:611-613).

## 8. Timings (64^3, 1 rank, 8 steps = 16 stages, -O3)

| max_grid_size (boxes) | mode | side-data setup (once) | bind | kernels | unbind | (bind+unbind)/kernels | ghost fill |
|---|---|---|---|---|---|---|---|
| 64 (1) | alias | 0.064 s | 0.003 | 0.702 | 0.00003 | 0.5 % | 0.003 |
| 32 (8) | alias | 0.073 | 0.012 | 0.903 | 0.0001 | 1.4 % | 0.009 |
| 16 (64) | alias | 0.099 | 0.028 | 0.815 | 0.0003 | 3.5 % | 0.024 |
| 8 (512) | alias | 0.169 | 0.065 | 1.161 | 0.0007 | 5.6 % | 0.078 |
| 64 (1) | copy | 0.068 | 0.094 | 0.747 | 0.026 | 16 % | 0.004 |
| 32 (8) | copy | 0.078 | 0.109 | 0.902 | 0.026 | 15 % | 0.011 |
| 16 (64) | copy | 0.096 | 0.160 | 0.790 | 0.032 | 24 % | 0.025 |
| 8 (512) | copy | 0.187 | 0.302 | 1.100 | 0.038 | 31 % | 0.080 |

Notes:
- **Kernels:** 168 ns per cell per stage with one box. With 8^3 boxes it is about 65 % more, because the FDS loops cover padded ranges (`-1:IBP1+1`), and a 12^3 padded region is 3.4× an 8^3 box.
- **Bind time in alias mode** is mostly the UVW_SAVE loop over box faces.
- **Setup** is the one-off cost of CELL, walls and scratch. It would be paid again on every regrid.
- Timings are single samples on a shared machine, so treat them as ±10 %.

## 9. Two-level results (static patch, 2 ranks)

**Setup:**
- Coarse grid 32^3. The fine patch covers coarse cells 8..23, i.e. fine cells 16..47.
- Blob Z = 0.05 + 0.5·exp(-r²/0.1²), starting at the centre, uniform velocity (1, 0.5, 0.25), dt = 1/256 on both levels (no subcycling), 128 steps.
- By T = 0.25 the blob has moved to the +x interface; by T = 0.5 it has left the patch.

**Ghost fill:**
- Level 0 uses FillBoundary.
- Level 1 uses FillPatchTwoLevels with piecewise-constant interpolation. That is the same "injection from the one coarse cell" that FDS uses on the fine side.

**Averaging down:** after each step, rho and ρZ are averaged down conservatively, and Z = ρZ/ρ is reset only on covered coarse cells.

**Mass drift without reflux:**
- Rho: -1.1e-5 at T = 0.5.
- Tracer (relative to total tracer mass): -1.1e-4 at step 16, -4.4e-4 at step 32, -2.0e-4 at step 64, -4.5e-4 at step 80, -4.2e-4 at step 128. It rises and falls as the blob flank crosses interface faces, which is the behaviour of a coarse/fine flux mismatch.
- The single-level runs conserve both to 1e-15.
- The drift is about 0.8 % of the blob's excess mass. **A FluxRegister is needed.**

### 9.1 Direct flux-mismatch measurement (`p1.cf_diag=1`)

- **Method:** after each box's DENSITY call, the driver reads the face fluxes FX·WORK_U·area that FDS just used (`p1_box_flux`). It does this for every box face on the patch surface: fine boxes on level 1, and covered coarse boxes on level 0. This needs `amr.max_grid_size=8`, so coarse boxes do not straddle the patch boundary. The driver checks this and warns if they do.
- **What is summed:** Σ 0.5·dt·(F_coarse,out − F_fine,out) over both stages and all steps. Each stage's flux gets weight 0.5·dt because FDS's predictor/corrector (Heun) update applies 0.5·dt·(F^n + F^*) per step.
- **Expectation:** without reflux, the composite (coarse uncovered + fine) mass should change by exactly this amount.
- **Same case as §9, 2 ranks, mgs 8.** The relative drift matches the mgs 16 run: rho −1.11e-5 and tracer −4.19e-4 at step 128.

| steps | quantity | predicted Σ0.5·dt·ΔF | measured Δmass |
|---|---|---|---|
| 16 | rho | −3.856587781e-06 | −3.856587781e-06 |
| 16 | tracer ρZ | −7.713175563e-06 | −7.713175563e-06 |
| 128 | rho | −1.520881505e-05 | −1.520881506e-05 |
| 128 | tracer ρZ | −3.041763011e-05 | −3.041763011e-05 |

- **Conclusion:** the mass error is entirely the coarse/fine flux mismatch, to 1e-10 relative. Nothing else loses mass: not the ghost fill, not average-down, not the kernels. Correcting the coarse cells next to the patch by these fluxes (reflux) restores conservation.
- **The rho mismatch is half the tracer mismatch** at both sampled steps (to 10 digits). So the background species' mismatch is −0.5× the tracer's. I have not analysed why the ratio is exactly this; it probably comes from how the limited face density and mass fractions combine in this blob setup.

**Smoothness** (`analysis/two_level_check.py`; mean / max |Z − exact| on coarse cells):

| T | region | two-level | single 32^3 | single 64^3 (coarsened) |
|---|---|---|---|---|
| 0.25 | covered | 4.1e-4 / 1.4e-2 | 1.2e-3 / 7.2e-2 | 3.9e-4 / 1.4e-2 |
| 0.25 | 1-cell ring outside patch | 4.0e-4 / 2.6e-2 | 1.0e-3 / 7.7e-2 | 2.9e-4 / 1.8e-2 |
| 0.25 | 2nd ring | 7.7e-4 / 2.7e-2 | 9.8e-4 / 3.4e-2 | 3.6e-4 / 1.4e-2 |
| 0.50 | rest of coarse grid | 6.1e-4 / 8.1e-2 | 6.9e-4 / 1.0e-1 | 3.4e-4 / 2.8e-2 |

- The fine level on the patch at T = 0.25 has mean error 4.8e-4, against 4.6e-4 for the 64^3 run on the same cells.
- The error jump between neighbouring coarse cells at the interface is 3.6e-2 for two-level and 5.5e-2 for 32^3.
- **So there is no interface artefact:** two-level errors next to the interface lie between the 32^3 and 64^3 errors.

## 10. Problems hit

1. **Pointer remapping alone is impossible** because of POINT_TO_MESH (§3). This needed W1/W2.
2. **Missing or unsafe defaults:** N_ZONE has no default (cons.f90:695). The settling flags default to T (cons.f90:231, 233). DENSITY does nothing at ICYC ≤ 1 (mass.f90:390). T_USED is an unallocated allocatable.
3. **FX shape:** FX has asymmetric bounds, so it cannot be a MultiFab, and it must be zeroed.
4. **Direct component access:** mass.f90:477 and :639 read `ALLOCATED(MESHES(NM)%M_DOT_PPP)` directly.
5. **Diffusivity in the reference:** DIFFUSIVITY = 1e-30 in the first FDS reference was not negligible at the ulp level, so the reference now uses 0. FDS honours `D_USER >= 0` (data.f90:2595, 2757).
6. **FP traps:** they are not active unless AMReX enables them (C++ main). hwloc/libxml2 hits a NaN inside MPI_Init when traps are on, avoided with `HWLOC_LIBXML=0`.
7. **Build environment:** `pkg-config` is missing, so AMReX's FFTW search must be pointed at a CMake shim (README).
8. **FDS time-step rule:** FDS shortens or lengthens the last step of every run by about 20ε (main.f90:719), and a restart carries the modified DT. Until I found this, P1 appeared to differ from FDS by up to 2.5e-14. Any bitwise regression test against FDS has to reproduce the T/DT sequence (§5).

## 11. Implications for the driver design

1. **One MESHES entry per AMReX box** is the unit FDS kernels understand. Two options:
   - Keep W2 (alias; lowest overhead at 0.5–6 %, but non-standard and gfortran-specific).
   - Make a small FDS-side change so a box can be bound without going through ALLOCATABLE components: e.g. a POINT_TO_BOX hook, or POINTER components. That removes the non-standard trick and the copies (W1 costs 15–31 %).
2. **Side data is large and rebuilt on every regrid.**
   - About 0.8–1 kB per cell, dominated by the 384-byte CELL_TYPE, plus 540 B per box-face cell for walls. Setup is 0.06–0.17 s per 64^3 at 1–512 boxes.
   - The driver should build it per box in parallel, and consider a slimmer CELL, since only SOLID/EXTERIOR/WALL_INDEX are used by mass.f90.
3. **Box size:** avoid small boxes. The kernel cost follows the padded box volume (+65 % at 8^3). max_grid_size ≥ 16 looks reasonable.
4. **Ghost cells:**
   - Fill 2 face-adjacent layers with plain copies at box faces. Keep edges and corners finite; they are read but do not affect results.
   - Do not emulate FDS's (ρZ)/ρ ghost rounding: it costs decomposition independence.
   - On the coarse side of a c/f interface, P1 uses averaged-down covered coarse cells (volume average of 8 fine cells). FDS uses an area average of the 4 face-adjacent fine cells (wall.f90:319-339). The difference is small but real.
5. **UVW_SAVE own-face policy (no MATCH_VELOCITY)** gives decomposition-independent single-level results whenever no species clipping happens (§6).
6. **Species clipping is per box and partly box-global** (§6.2). Results will depend on the box layout whenever clipping happens, exactly as FDS results depend on the mesh layout. Decide whether that is acceptable, or plan an FDS change to CHECK_MASS_DENSITY. → Decided as D-031; the gather form in §13 removes the layout dependence within a level.
7. **Reflux is required:**
   - Measured drift is up to 4.5e-4 relative tracer mass in 128 steps.
   - The hook is `p1_box_flux` (FX·WORK_U·area per box face and species, available after DENSITY). A FluxRegister would take 0.5·(predictor + corrector) face fluxes at c/f faces and reflux level 0 before average_down (marked `// FLUXREG:` in `p1_driver.cpp`).
   - Averaging rho and ρZ down conservatively already works.
8. **Parallelism:**
   - No OpenMP threading over boxes. The kernels work through global MESH_POINTERS module variables, and DENSITY updates `T_USED(3)` (mass.f90:355, 768) without protection.
   - MPI over boxes works as-is; MESHES must be sized to the global box count on every rank.
9. **Regression testing against FDS:** bitwise comparison is possible, but it has to reproduce FDS's T/DT sequence, including the end-of-run DT rule (main.f90:719) and the fact that a restart carries the modified DT. `fds_ref/check_fds.sh` is a working template.
10. **Diffusion is out of P1's scope:** DEL_RHO_D_DEL_Z comes from DIVERGENCE_PART_1 (divg.f90). Even ~1e-28 values change results at the ulp level, so a later prototype adding divg must produce it box by box consistently.

## 12. Not done / open

- No FluxRegister/reflux was implemented. The flux mismatch that it would correct is measured and accounts for all the drift (§9.1).
- No regridding (static grids by design) and no subcycling.
- `zero_fx=0` as a negative control is inconclusive (§7).
- No FDS multi-mesh reference: the FDS comparison is FDS single-mesh against P1 with 1 or 64 boxes.

## 13. Gather-form clipping, D-031

**Result: pass.** With `p1.clip=gather`, SUPERBEE on the sharp FDS initial field (species clipping in 14 of 16 stages) is bitwise identical for max_grid_size 32 vs 8 and 1 vs 2 ranks. It is also bitwise identical to single-mesh FDS, both with 1 box and with 64 boxes on 2 ranks. That second case differed by 5.2e-5 in RHO (6.9e-5 in ZZ) with FDS's own clip. All numbers come from `check_clip.sh` (exit 0; output in `prototypes/p1_mass_shim/runs/clip/check_clip_output.txt`). Line numbers below were re-verified on FireX `36975d765f`.

### 13.1 What FDS does and why it depends on the layout

`CHECK_MASS_DENSITY` (mass.f90:775-963) is called from inside `DENSITY` (predictor :518, corrector :700), once per mesh. It runs in this order:

1. **Density pass.**
   - Loop over interior cells. A cell outside [RHOMIN, RHOMAX] sets CLIP_RHOMIN/CLIP_RHOMAX (:815/:819).
   - Weights from its 6 neighbours: MASS_N(d) = 0 wherever `CELL%WALL_INDEX(d)/=0` (:831-836).
   - Scatter into DELTA_RHO (:840-846).
   - Apply to the interior iff a flag is set (:853-854).
2. **Species pass, per species.**
   - The same clip against [0, RHOP], where RHOP is the density-clipped value. It sets CLIP_RHO_ZZ(N) (:889).
   - Weights at :907-912, CONST at :915, scatter into DELTA_RHO_ZZ at :916-922.
   - Interior apply iff CLIP_RHO_ZZ(N) (:927, :931-937).
3. **Early return** if nothing clipped (:943).
4. **Renormalisation** of every interior cell so that ΣρZ = ρ (:947-961).

The layout dependence has two causes:
- **(a)** Every mesh-boundary face carries a wall cell, so MASS_N = 0 across box faces. A clipped cell next to a box face spreads its mass over fewer neighbours, which gives a different CONST.
- **(b)** Steps 1-apply, 2-apply and 4 are gated by per-mesh flags. Unclipped cells therefore change, or not, depending on which box they share with a clipped cell.

### 13.2 What was built (algorithm unchanged, loop structure changed)

**Files**, all in `prototypes/p1_mass_shim/`; nothing under this repository was edited:

- **`mass_split.f90`** (module `P1_DENSITY_SPLIT`), generated by `gen_mass_split.py`.
  - `DENSITY` cannot be run without its internal `CALL CHECK_MASS_DENSITY`, because it is a CONTAINed routine called unconditionally. So a prototype copy of `DENSITY` is split at that call:
    - `DENSITY_PRE_CLIP` = mass.f90:365-513 (predictor, up to `RHOS=SUM(ZZS)`) + :584-695 (corrector, up to `RHO=SUM(ZZ)`).
    - `DENSITY_POST_CLIP` = :365-397 (declarations, early returns, POINT_TO_MESH) + :522-582 + :704-768 + `CLIP_PASSIVE_SCALARS` :966-989.
  - Changes against the source:
    - The line ranges are copied verbatim.
    - The two `CALL CHECK_MASS_DENSITY` lines (:518, :700) are left out.
    - All `!$OMP` directive lines are dropped (P1 runs 1 thread).
    - The routines are renamed.
  - Nothing else changed. `CHECK_MASS_DENSITY` itself is **not** copied. It is re-expressed in `clip_gather.f90`.
- **`clip_gather.f90`** (module `P1_CLIP_GATHER`).
  - Explicit-shape arguments only, with no module pointers or globals. Arrays are indexed by **global** AMReX cell index. Array bounds (FAB incl. ghosts) and loop bounds are passed separately, and RHOP, RHO_ZZ, the mask and DELTA each get their own bounds.
  - `CLIP_TERMS` recomputes a cell's clip exactly as FDS writes it: the same expressions and operation order for VC1/VC/MASS_C/MASS_N/`SUM(MASS_N)`/CONST (:801-839, :875-915), and the same SOLID/range test order (density: range first, :809-811; species: SOLID first, :885-888).
  - `GATHER` visits, for each target cell, the sources in the order the K,J,I scatter delivers to that cell: (k-1), (j-1), (i-1), self, (i+1), (j+1), (k+1). It accumulates `DEL = DEL + CONST*SUM_MASS_N/VC(0)` (self) or `DEL = DEL - CONST*MASS_N(d)/VC(d)` (neighbour, d = direction of the target seen from the source). DEL starts from 0 as in :785/:871, so every floating-point operation is the one FDS performs, in the same order.
  - Only single-mesh-FDS source cells contribute (mask component 7, below).
- **Build.** Both files are compiled with the FDS flags (CMake `FDS_COMMON_FLAGS` + `-O3`, or the debug set).
- **`p1_shim.f90`.** `p1_box_stage` gets a `PART` argument: 0 = unmodified `MASS_FINITE_DIFFERENCES + DENSITY` (the FDS clip, still the default); 1 = `MASS_FINITE_DIFFERENCES + DENSITY_PRE_CLIP`; 2 = `DENSITY_POST_CLIP`. There is also a new `p1_level_wall_mask`.
- **`p1_driver.cpp`.**
  - Selects the clip with `p1.clip=fds|gather`.
  - `clip_level()` runs the passes below.
  - Also added: per-stage diagnostics, the `p1.clip_ghost=fill2|redundant` ghost variants, and the `p1.clip_poison` negative control.

**Face mask** (Legacy Mapper's correction; `p1_level_wall_mask`, `P1Amr::build_clip_mask`). It is an iMultiFab with 8 components and ng=2, so the face mask and SOLID are defined on valid+2:
- **Comp 0:** SOLID.
- **Comps 1-6:** WALL_INDEX(-1,+1,-2,+2,-3,+3).
- **Comp 7:** 1 iff the cell is a clip source in the single-mesh loop 1:IBAR.
- **Where it comes from.** The mask is **not** built from boundary types. It is taken from the WALL_INDEX that single-mesh FDS has for the whole domain:
  - The driver allocates one extra `MESHES` entry per level (NM = nbox_total+lev+1) and runs the unchanged `P1_BOX_SETUP` on the level's full domain box as a single mesh.
  - That gives every domain-boundary face a wall cell, periodic faces included, the analogue of init.f90:76-107. This matters for parity: single-mesh FDS never redistributes across the domain boundary, even in a periodic direction.
  - The mask is then copied from the periodic-image cell.
- **Faces between boxes of the same level** are interior faces of that single mesh, so their WALL_INDEX is already 0. Nothing has to be zeroed, and the shim's own per-box WALL_INDEX (an EXTERNAL_WALL on every box face) is not used.
- **Obstructions:** P1 has none. With OBSTs, the same route would pick up exactly the exposed faces, both sides of thin obstructions, and skip permanently covered faces (init.f90:199-201), because it is FDS's own wall-cell set for the domain.
- **Memory:** the domain mesh costs about 16 MB at 32^3 (the CELL array). A production driver should build the 8-int (or bit-packed) mask directly from FDS's global wall data rather than a full domain mesh.
- **Levels > 0** (Chief Architect's rule):
  - Faces towards cells that the level's BoxArray does not cover (coarse/fine faces) stay nonzero, like an FDS mesh interface.
  - Coarse-interpolated ghosts are not sources (comp 7 = 0).
  - Only faces between boxes of the **same** level are open.
  - Flags are OR-reduced per level.
  - The reason: gathering from interpolated ghosts would move mass across levels without a matching reflux. D-031's layout independence is therefore a property of how boxes are split *within* a level at a fixed level structure. Results are not expected to stay the same when the level structure changes.

### 13.3 Pass and reduction sequence (per stage, per level)

This follows the Chief Architect's sequencing correction.

| # | step | FDS lines | notes |
|---|---|---|---|
| 0 | `DENSITY_PRE_CLIP` on every box | 365-513 / 584-695 | ρZ and ρ (RHOS/ZZS in the predictor, RHO/ZZ in the corrector) are now current in valid cells only |
| 1 | **FillBoundary #1**: RHOP (1 comp) and RHO_ZZ (NS comps), ng=2 | – | the extra exchange |
| 2 | density gather on the valid box → DELTA_RHO; local CLIP_RHOMIN/MAX from **valid cells only** | 799-849 | |
| 3 | **OR-reduction #1** (`amrex::ParallelAllReduce::Or`) of CLIP_RHOMIN, CLIP_RHOMAX over all boxes and ranks | – | host-side |
| 4 | iff reduced density flag: density apply on the interior (valid) only | 853-854 | |
| 5 | iff reduced density flag: **FillBoundary #2** of RHOP | – | needed because step 4 updates only valid cells, and the species pass reads the clipped RHOP of ghost cells (as QMAX of source ghosts). RHO_ZZ is not changed by the density step when NS>1, so it is not refilled. NS=1 case (:858-861): RHO_ZZ(:,1)=RHOP iff the reduced flag is set, then return |
| 6 | species gather for every N (valid box) → DELTA_RHO_ZZ(:,N); local CLIP_RHO_ZZ(N) from valid cells only | 870-925 | all N can be gathered before the reduction, because species N reads only RHO_ZZ(:,N) and RHOP |
| 7 | **OR-reduction #2** of CLIP_RHO_ZZ(1:NS) | – | host-side |
| 8 | species apply for each N iff **reduced** CLIP_RHO_ZZ(N) | 927, 931-937 | |
| 9 | renormalisation iff reduced (CLIP_RHOMIN ∨ CLIP_RHOMAX ∨ any CLIP_RHO_ZZ) | 947-961 (replaces the per-box return at 943) | |
| 10 | `DENSITY_POST_CLIP` on every box | 522-582 / 704-768 | |

Gating by the reduced flags matters for bitwise parity, not just for tidiness. Applying unconditionally would, for example, turn a −0.0 ρZ into +0.0, and renormalising when single-mesh FDS would not changes ulps (see the mutation test in 13.4).

**Ghost variants** (`p1.clip_ghost`):
- **`fill2`** (default) is the sequence above.
- **`redundant`** (the target form):
  - Step 1 copies RHOP into a temporary MultiFab with **ng=3** and fills it; RHO_ZZ is filled with ng=2.
  - Steps 2/4 gather and apply the density clip on **valid+1**, so step 5 is not needed. The clipped valid values are copied back.
  - The stencil (mass.f90:775-922) needs ng=3 for RHOP/RHOS, ng=2 for RHO_ZZ, and the face mask/SOLID on valid+2.
  - The temporary is needed because the shim's FDS-bounds arrays (RHO/RHOS/ZZ/ZZS bound as `-1:IBP1+1`, i.e. ng=2, required by POINT_TO_MESH and the kernels' fixed loop ranges) cannot hold a third ghost layer without changing the FDS array bounds.
  - Flags still come from valid cells only: a redundant clip in any ghost (same-level, coarse/fine or periodic image) never sets a flag.
- **Both variants are bitwise identical** in every case tested (section 6 of `check_clip.sh`: SUPERBEE and rhomax=1.85, g32/np1 and g8/np2, with and without poisoned ghosts, vs FDS; 64^3 blob g8/g16/g32; two levels). **The bitwise-equality check was done.** All other results in this section use `fill2` unless stated.

### 13.4 Acceptance results (`check_clip.sh`, 8 steps, 32^3 periodic, ≤2 ranks)

| check | variant | result |
|---|---|---|
| SUPERBEE, FDS initial field (check_decomp setup), gather: g32/np1 vs g8/np1, g32/np2, g8/np2 | fill2 | **bitwise** (RHO, ZZ1, ZZ2) |
| same, gather g32/np1 vs `p1.clip=fds` g32/np1 (1 box = single-mesh FDS behaviour) | fill2 | **bitwise** |
| control: `p1.clip=fds` g32/np1 vs g8/np2 | – | still differs (as in §6) |
| FDS restart case (check_fds setup, FDS T/DT sequence), gather, 1 box, 1 rank vs FDS | fill2 | **bitwise** (RHO, ZZ, TMP) |
| same, gather, **64 boxes, 2 ranks vs FDS** | fill2 and redundant | **bitwise** (FDS clip: RHO 5.2e-5, ZZ 6.9e-5 on 7745-8033 cells) |
| density clip active, `p1.rhomax=1.85` (CLIP_RHOMAX, 2520 cells over 16 stages) and `p1.rhomin=1.333` (CLIP_RHOMIN, 154354 cells): gather g32/np1 vs g8/np1, g32/np2, g8/np2, and vs fds-clip 1 box | fill2 (+ redundant for rhomax) | **bitwise** |
| two levels, static patch, blob σ=0.03, dt=1/64, 12 steps, 7970 species-clipped fine cells: gather mgs 16/np1 vs 8/np1, 16/np2, 8/np2 (and redundant 8/np2), levels 0 and 1 | both | **bitwise**. FDS clip: level 1 differs, max \|ΔZ\| 0.99 in this CFL-1 stress case (not investigated further) |
| `check_decomp.sh` (P1 regression, default `p1.clip=fds`) | – | exit 0, unchanged (all bitwise, SUPERBEE/FDS-init still EXPECTED_DIFF) |
| `fds_ref/check_fds.sh` (P1 regression) | – | 4 BITWISE, EXPECTED_DIFF case unchanged (RHO 5.222e-05) |

- **No "first differing cell" analysis was needed:** nothing differs.
- **The clip fires.** In the FDS restart case, species clipping happens in 14 of 16 stages (1 to 168 cells per species per stage; 489 per species, 978 in total over the 8 steps). The counts are identical for 1 box and 64 boxes/2 ranks. No density clipping happens with the default &CLIP limits (as in §4), which is why the rhomax/rhomin cases were added.
- **Mass.**
  - Per stage, the gather clip changes each species' total mass by at most 8.7e-16 relative, i.e. round-off, the same as FDS's clip on this data. FDS's clip is not exactly conservative in general (the `MIN(1,…)` in CONST, the final clamps, the `SUM_MASS_N<=TWO_EPSILON_EB` skip). The gather form inherits exactly FDS's behaviour, because the fields are bitwise equal to FDS's.
  - End-of-run relative change: rho −2.4e-15 / tracer −1.5e-15 (1 box) vs −3.3e-16 / −2.2e-16 (64 boxes). The fields are bitwise equal; only the summation order of `MultiFab::sum` differs.
  - With density clipping, species mass does change (the density clip moves mass, and renormalisation follows it), exactly as FDS.
- **Only valid cells may set clip flags** (Chief Architect's condition (a)). Test: `p1.clip_poison=1` sets RHOP=1e3 and tracer ρZ=−1 in every non-source ghost after the pre-clip fill.
  - GODUNOV (no real clipping), 64 boxes/2 ranks: no flag is ever set, and the result is bitwise equal to FDS, in both variants.
  - SUPERBEE and rhomax=1.85 with poison: bitwise equal to the unpoisoned runs and FDS.
  - Mutation check (done by hand, not in the script):
    - Removing *only* the valid-cell guard on the flag count is **not** caught, because a second guard also holds: non-source ghosts are never evaluated as sources. Same-level ghosts that are sources duplicate valid cells of other boxes, whose flags the OR already includes.
    - Removing both guards **is** caught in the redundant variant: flags are set in all 16 stages, renormalisation runs when FDS would not, and RHO differs from FDS by up to 1.1e-15.
    - The fill2 variant never clips ghosts, so it is immune.
- **Debug build.** `-O0 -fcheck=all -fbounds-check`, FP traps on: g8/np2 gather in both variants, rhomax 100 and 1.85, poisoned. No trap or bounds error, and results are bitwise equal to the -O3 build.

### 13.5 Cost (64^3, blob σ=0.03 with clipping, 8 steps = 16 stages, -O3, shared machine, single samples, ±10 %)

Times in s; "ghost fill" is P1's existing per-stage fill; the "fds" box loop includes FDS's own clip.

| mgs (boxes), 1 rank | fds: loop / box loop | gather fill2: loop / kernels / extra FillBoundary / clip passes | gather redundant: loop / extra fill / clip | existing ghost fill |
|---|---|---|---|---|
| 64 (1) | 0.81 / 0.81 | 1.01 / 0.68 / 0.0036 / 0.32 | 1.04 / 0.0083 / 0.33 | 0.0031 |
| 32 (8) | 0.96 / 0.95 | 1.13 / 0.79 / 0.0090 / 0.33 | 1.18 / 0.016 / 0.35 | 0.0067 |
| 16 (64) | 0.87 / 0.84 | 1.31 / 0.91 / 0.028 / 0.35 | 1.36 / 0.037 / 0.38 | 0.017 |
| 8 (512) | 1.30 / 1.22 | 1.81 / 1.32 / 0.081 / 0.36 | 1.97 / 0.10 / 0.45 | 0.058 |

- **Extra FillBoundary cost** (fill2, no density clip): 1.2-1.6× P1's existing ghost fill, i.e. 0.4-4.5 % of the loop time. The data volume is about the same: RHOP + NS ρZ components with ng=2.
- **With density clipping** (`rhomax=1.7`, 15 of 16 stages clip):
  - The second exchange adds a little on 1 rank: fill2 extra fill 0.010/0.031/0.095 s at mgs 32/16/8, against 0.017/0.036/0.097 s for redundant.
  - On 2 ranks at mgs 16: fill2 0.033 s vs redundant 0.030 s.
  - So on this CPU development machine the redundant variant only pays off once the second exchange costs latency (many ranks, GPU).
- **The clip passes themselves** cost about 20 ms per stage at 64^3, about 80-110 ns per cell per stage. That is about 2.6× FDS's own scatter clip (≈0.12 s over 16 stages, from the box-loop difference at mgs 64), and it raises the loop time by 25-40 %.
  - The reason: this prototype recomputes `CLIP_TERMS` for 7 sources per target (×(1 + NS) passes) through a non-inlined call with explicit-shape arguments, even though almost no cell clips.
  - An obvious and still bitwise-safe optimisation: first compute a per-cell "clipped" byte and the 7 terms CONST·MASS_N(d)/VC(d) on valid+1, then gather only where a neighbour clipped. It was not done tonight.

### 13.6 What a GPU version would need

- **A C++ port of `CLIP_TERMS`/`GATHER`** in `amrex::ParallelFor`, since the Fortran cannot run on the device.
  - The gather form is already race-free (each thread writes only its own cell), deterministic, and has no atomics. That is the main point of D-031 for GPUs: the FDS scatter would need atomics and would lose bitwise reproducibility.
- **Keep the literal operation order** of every expression (`DY*DZ` then `DX*`, `CONST*MASS_N/VC`, `SUM` over MASS_N(-3:3) left to right, `MAXLOC` first-index tie rule in the renormalisation). Compile with FMA contraction off (`-fmad=false` for nvcc, `-ffp-contract=off` for hipcc/clang) if bitwise parity with CPU FDS is required. The CPU build here uses no FMA (plain `-O3`, no `-march`).
- **The mask as device data.** It is 8 ints per cell now; 1 byte would do (6 face bits + SOLID + source). It has to be rebuilt on regrid from the level's BoxArray and FDS's global wall/OBST data (the domain-mesh trick in 13.2 is only for the prototype).
- **Two host synchronisations per stage per level** for the two OR reductions. Each is a device reduction (`ReduceOps`) followed by `ParallelAllReduce::Or`. The per-species flags should be packed into one reduction. The density flag also decides whether FillBoundary #2 runs; the redundant variant removes that exchange, but the density-flag reduction is still needed before the apply.
- **Early exits:** most cells return after one comparison, so warps diverge only near clipped cells. The two-phase form from 13.5 (flag + term arrays, then gather) suits GPUs better than recomputing 7× per target.
- **Scratch:** DELTA_RHO (ng=1 for the redundant variant) and DELTA_RHO_ZZ (NS comps), plus the ng=3 RHOP temporary for the redundant variant. These should come from The_Arena, not per-box allocation.
