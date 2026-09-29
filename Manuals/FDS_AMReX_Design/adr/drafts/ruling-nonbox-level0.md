# Ruling: non-box level 0 (meshes whose union is not a box)

Owner: Chief Architect, with the AMReX Integration Lead. Status: ruling final; amended 2026-09-26 by §6 (level-0 padding N1a, N6 fallback) and by the A-46 vent count (N2, N3 wording). Inputs: `pressure/04-masked-domain-pressure.md` (Pressure Solver Lead proposal, draft), requirements IR-002, FR-006, FR-037, FR-039, D-021, D-032 (5), risk R-47. Affects 40 IN cases of the FR-006 set (`vv/scope_case_list.csv`, `nonbox_domain = 1`), including the restart anchors, `hallways`, `simple_duct` and `HVAC_leak_exponent`.

## 1. Rulings

- **N1. Level 0 covers the bounding box; gaps are static solid cells.** In AMR mode level 0 is the bounding box of the level-0 `&MESH` blocks. Cells in no `&MESH` block are gap cells: solid for the whole run, on every level, never refinable (they are excluded from the refinable region, IR-008), never removable by DEVC/CTRL. Gap cells are a separate mask class from OBST cells, so OBST create/remove and the FR-040 checks never touch them. Gap boundaries lie on level-0 cell faces, so they refine exactly on every finer level. The AMReX level-0 domain is this bounding box padded to a multiple of the level-0 blocking factor (§6.1, N1a); padding cells are gap cells.
- **N2. Gas-gap faces reproduce FDS exterior walls.** In FDS every mesh face with no neighbour is an exterior wall with the default surface (`init.f90:76-107`). A gas-gap face therefore gets a wall record with `DEFAULT_SURF_INDEX`, and any `&VENT` that FDS snaps onto that mesh face applies there exactly as in FDS, including `OPEN` vents. The V&V count (A-46, `vv/gap_face_vents.md`) finds gap-face vents in 21 of the 40 cases: `OPEN` vents in 11 (`bi_dir`, `velocity_bc_test`, `qfan_multi`, `back_wall_test`, `hallways`, `shrink_swell`, `hot_spheres`, `device_restart_a/b/base_case`, `LS4_ember_ignition`), only non-`OPEN` vents (MIRROR, velocity and other surfaces) in 10, and no PERIODIC vent. The part of a gap-face vent backed by an OBST stays a solid wall, as in FDS. Where `&SURF DEFAULT=T` is set (15 cases), gap walls get that surface, not INERT. These wall records are ordinary wall state (FR-041); they never move, because the gap never refines.
- **N3. Pressure: single-level masked MLMG ("gap mask", proposal option (a2)).** Overset mask 0 on gap cells; beta = 0 on gas-gap faces; the Dirichlet rule with the known value behind `OPEN` vents on gap faces; one deterministic pinned cell per sealed component; D-032 mean removal and gauge per zone and connected component. Obstructions inside the meshes keep the FFT-path E-1 treatment. **Solver configuration amended by §5.5** (a > 0 in masked cells, HYPRE bottom on level 0, no N-solve); checked in §5.1 and §5.4.
  - FR-039 selection gets a third branch: `FFT::Poisson` for a single level on a box domain; single-level masked MLMG for a single level with gap cells; composite MLMG once a refined level exists, with the same gap mask on level 0. The mask is static, so a masked run never switches solver.
  - D-021 is unchanged in substance (a masked level 0 is already outside its validity domain). D-032 (5) reads "FFT whenever the hierarchy is a single uniform level on a box domain (no gap cells)".
  - FR-039's "FFT vs MLMG on identical input" check is replaced for masked domains by a single-solve comparison against FDS `SOLVER='UGLMAT'` on derived copies of `hallways`, `device_restart_a` and `simple_duct`, within eps_H (proposal §3 item 2), including its [VERIFY] on where UGLMAT places the `OPEN` value. `hallways` and `device_restart_a` test `OPEN` placement on a gap face. `simple_duct` has no gap-face vent (its HVAC and leak vents sit on interior planes and attach to obstructions, A-46), so it tests sealed gap walls only. No replacement case is needed.
- **N4. Phase.** The gap mask moves to Phase 2, because the FR-006 uniform subset gates the Phase 2 exit. Its scope there is small: one level, static mask, no coarse/fine interface, no OBST masking.
- **N5. Zones.** Mean removal and the pinned cell are per connected component, inside each zone. If a `&ZONE` can span several gap-separated components, each component gets its own pin and mean (Pressure Solver Lead verifies `&ZONE XYZ` handling in `read.f90`; proposal O5). **Checked (§5.3):** in FDS a sealed zone cannot span gap-separated components, so one ΔP_zone and one `D_PBAR_DT` per zone stay correct. Compatibility, pin and gauge are applied per singular operator component whatever the zone ID, and a setup assertion checks that each sealed zone is exactly one component.
- **N6. Large bounding boxes (`stairwell`, fill 0.091, about 18 M bounding-box cells).** **Estimated** (`amrex/n6-stairwell-memory-estimate.md`): a covering level 0 needs 15-20 GB against NFR-031's 12 GB, so the fallback is needed. **The fallback is ruled in §6.2**: a level 0 built from the mesh union, with no overset mask, under setup-checked conditions. This replaces the earlier wording "omit all-gap boxes", which would keep gap cells, and therefore a mask, inside a non-covering level 0.
- **N7. Performance (NFR-030).** Masked cases are measured first. If they miss NFR-030 only because of the bounding-box overhead, they get a recorded exception rather than a new solver.

## 2. Rationale
- It keeps FDS's discretisation on gap walls exactly: exterior walls with no-flux pressure faces, the default surface for heat and momentum, and vents on those faces.
- It is the only evaluated route that handles `OPEN` vents on faces inside the bounding box and keeps disjoint zones disjoint.
- It reuses the E-2 masking machinery and MLMG, already validated in P2; no new numerical method.
- It gives a same-discretisation reference (UGLMAT) for an eps_H check.

## 3. Rejected alternatives
- **Non-covering level 0 as the default** (solver BoxArray over the meshes only, Neumann coarse/fine condition on the gaps). It cannot express `OPEN` vents on gap faces. AmrCore supports a non-covering level 0 through `PostProcessBaseGrids`, but MLMG supports it only without an overset mask (§6.2). Kept only as the N6 fallback, in the mesh-union form of §6.2.
- **Bounding-box FFT with FDS obstruction iteration.** Approximate no-flux where FDS is exact, extra iterations on the largest walls, merged zones, and no `OPEN` on interior faces.
- **Capacitance or immersed correction to FFT.** New, unvalidated code at or above MLMG cost with no accuracy gain.
- **Gap cells as gas with zero velocity forced.** Merges disjoint zones and is not FDS's discretisation.
- **Rejecting non-box inputs in AMR mode (FR-004).** Loses 40 IN cases, including the restart anchors.

## 4. Follow-ups
- Spec & Program Lead: IR-002 (gap cells, N1-N2), FR-037/FR-039/D-032 (5) wording (N3), FR-006 phase note (N4), R-47 mitigation chosen, `HVAC_leak_exponent` IN.
- V&V Lead (A-46): list every vent on a gap-facing mesh face across the 40 cases, split into `OPEN` and others. **Done** (`vv/gap_face_vents.md`), folded into N2/N3.
- Pressure Solver Lead: MLMG convergence with large masked regions (proposal O3); N5 check. **Done, §5.** Still open for the Pressure Solver Lead: the composite (multi-level) masked solve (§5.1, last item) and the O6 timing.
- AMReX Integration Lead: N6 estimate. **Done** (`amrex/n6-stairwell-memory-estimate.md`), ruled in §6. Still open: the §6.1 and §6.2 [VERIFY] items in the Phase 2 prototype.

## 5. Pressure checks for sign-off (results)

Four checks: (a) one pinned cell per disconnected part (gap list B3; R3-T ruling §1.3 (i)); (b) the `P_0` ramp (B6; FR-034 ruling §2.2 (h)); (c) N5; (d) convergence of the masked branch against FDS (O3).

Evidence:
- Source reads: FireX `36975d765f` and AMReX `99ddfda`.
- The owner-provided FDS pressure-solver lecture slides (04_Pressure).
- Small runs under `scratch/pressure-signoff/`:
  - `mlmg_mask/` is a single-level masked MLMG test (AMReX `99ddfda` with the HYPRE 2.32 bottom solver, np ≤ 2). It uses the gap layouts of `hallways` and `device_restart_a` on their bounding boxes, plus two split test boxes. The right-hand side is synthetic, with its mean removed per component. The residual reported is an independent 7-point true residual over all gas cells, pinned rows included.
  - `fds_restart/` and `fds_hallways/` are FDS reference runs on derived copies of the two inputs (np = 2).

How numbers are labelled: FDS numbers are measured, MLMG iteration counts and residual ratios are measured, and conversions are derived.

### 5.1 (a) Pinned cell per disconnected component: PASS, given the §5.5 configuration

**This is FDS's own rule.**
- FDS builds one matrix per connected-zone group. UGLMAT marks the matrix indefinite when the group has no Dirichlet face (`pres.f90:4995-5007`); ULMAT does the same per mesh (`pres.f90:2784-2785`).
- Each indefinite matrix gets:
  - the arithmetic mean of the RHS removed over its unknowns (`pres.f90:3369-3406`, `1683-1700`);
  - its last unknown pinned (matrix at `pres.f90:2965-2975`, `4749-4760`; `F_H(last) = 0` at `pres.f90:1748`, `3447`);
  - a mass-weighted gauge shift after the solve (`pres.f90:3495-3540`).
- Zones join a group only through gas faces or open/interpolated boundaries (`MERGE_PRESSURE_ZONES`, `divg.f90:1295-1321`), so a group is one gas-connected component.
- The slides agree that an all-Neumann matrix is singular (slide 19) and give one matrix per pressure zone per mesh (slides 34, 39). They do not describe the pin.
- Difference: the source uses one matrix per connected-zone group, not per zone (`pres.f90:1205-1213`). The source wins, and the rule is unaffected.

**It is consistent with D-032.**
- With the mean removed per component, the pinned equation that is dropped holds automatically. Its residual is minus the sum of the component's other residuals (derived).
- Measured on every converged case:
  - true residual, pinned rows included: ≤ 9e-9 · max|RHS|;
  - per-component residual sums: ≤ 4e-9.
- Overset-mask pins and a-coefficient pins give the same solution within 1e-13 to 1.6e-11 relative L2. This held on the `split`, `wall`, `hallways` (sealed and open) and `device_restart_a` layouts (sealed and open); eps_H ≥ 1e-8.

**It works when a level splits into parts.** Three split cases converge with the §5.5 configuration:
- `split`: two sealed boxes separated by a 6-cell gap.
- `wall`: one box cut by a β = 0 plane.
- The `device_restart_a` walls (y = 6.8–7.2 m) masked with the door closed. This gives one sealed part of 4,354 cells plus one part holding the OPEN vent. FDS reports the same 4,354 cells for its auto zone 1.

**How parts are identified.** Use a flood fill over the faces with β ≠ 0 of the operator actually solved, with OPEN-known cells counted as Dirichlet. Do not use the zone index. At level 0 of an FDS input the two agree (§5.3), and the setup asserts it.

**Failures found; §5.5 exists because of them:**
1. **NaN from the default bottom solver.** MLMG's default BiCGStab returns NaN everywhere when masked cells have a = 0 (measured on `split`).
   - `MLABecLaplacian::normalize` divides by the diagonal α·a + β·Σb with no overset guard (`AMReX_MLABecLaplacian.H:1340-1400`; kernel `AMReX_MLABecLap_3D_K.H:61-76`). BiCGStab calls it at `AMReX_MLCGSolver.H:218, 270, 302`.
   - A gap cell has a = 0 and all its b = 0, so the diagonal is zero.
   - With a > 0 in the masked cells the same case converges.
2. **One pinned cell stops level-0 MG coarsening.** A cell with mask 0 inside gas makes MLMG stop at the first coarse level that mixes known and unknown cells (`AMReX_MLCellABecLap.H:192-241`), so the solve becomes bottom-solver only. Measured on sealed `hallways` (524,288 bounding-box cells, one pin):

   | Bottom solver | Result |
   |---|---|
   | BiCGStab | 50 MLMG iterations, 38.2 s |
   | CG | no convergence (resid/bnorm 0.59 after 100) |
   | HYPRE (BoomerAMG through AMReX's overset-aware IJ interface, `AMReX_HypreABecLap3.cpp`) | 3 iterations, 0.80 s |

3. **N-solve does not converge.** `MLMG::setNSolve(1)` (`AMReX_MLABecLaplacian.H:1402-1515`) failed on `split`, `wall` and `hallways`: resid/bnorm 0.01–0.31 after 100 iterations.
4. **a-coefficient pins need aligned walls under geometric coarsening.** They converge only when every β = 0 face lies on a 2^L boundary of the coarsening:
   - wall on x-face 32: 24 V-cycles to 1e-10;
   - wall on x-face 29: stalls (resid/bnorm 0.14 after 100 V-cycles).

   This confirms R3-T §1.3 (i)'s warning that coarse MG levels reconnect the parts. Face-β averaging (`AMReX_MLABecLaplacian.H:711-737`) drops a wall that falls inside a coarse cell.

**Not run.** A composite (multi-level) masked solve. On fine AMR levels the overset mask is coarsened without the mixed-cell check (`AMReX_MLCellABecLap.H:256-290`). Components then have to be found on the composite graph. **[VERIFY]** in the Phase 4 prototype.

### 5.2 (b) `P_0` ramp: PASS, no new constraint

**What FDS does.**
- The table is one global reserved ramp over [ZSW, ZFW]: `ZS_MIN`/`ZF_MAX` without HVAC (`read.f90:2230-2231`); with HVAC, also `DZS_MAX`/`DZF_MAX` and the duct-node heights (`read.f90:2227-2228`). The ramp itself is created at `read.f90:2301-2312`.
- `ZS_MIN`/`ZF_MAX` are the min/max over all meshes (`read.f90:775-776`).
- The table has `NUMBER_INTERPOLATION_POINTS` = 5000 entries (`read.f90:10503`), with `RDT` = 5000/span (`read.f90:10597`). It is integrated at `init.f90:465-469`.
- `EVALUATE_RAMP` returns the nearest entry, `NINT(position·RDT)`, with no interpolation (`func.f90:852-854`).

**What this means for spacing.**
- The spacing is span/5000: 0.56 mm for `device_restart_a` (2.8 m) and 0.8 mm for `hallways` (4 m). It does not depend on any mesh dz.
- The nearest-entry error is ≤ ρ·g·span/10000, about 5e-3 Pa for a 4 m span (derived).
- FDS never takes differences of `P_0`:
  - the hydrostatic divergence term uses ρ_0·g analytically (`divg.f90:684`);
  - the momentum terms average ρ_0 (`velo.f90:827, 1327`).
- `PBAR` starts as the same lookup at the same z (`init.f90:1051-1054`), so ambient cancellation is exact.
- A finer dz only means that, below span/5000, two vertically adjacent cells can share an entry. That is harmless.

**AMR rule** (refines FR-034 §2.2 (h)):
- Build the table once, from the level-0 global extents (the union of the `&MESH` z-extents, which is the unpadded bounding box; never the §6.1 padded domain). With HVAC, use the level-0 `DZS_MAX`/`DZF_MAX`, not the finest dz.
- On every level, evaluate `P_0(z_k^ℓ)` with the same `EVALUATE_RAMP` at that level's cell centres, and set `PBAR_ℓ = P_0 + ΔP_zone`. ΔP_zone is one scalar (`mass.f90:548, 730`), the same on all levels and all components.
- Gap cells do not change the extents, so all gap components share the one function.
- Do not interpolate `P_0` between levels and do not rebuild the table per level. Either would break the exact cancellation and parity with FDS on level 0.

### 5.3 (c) N5, zones spanning gap components: PASS, with a documented rule and a setup assertion

**In FDS, a sealed zone cannot span gap-separated components.**
- Only the first `XYZ` of a `&ZONE` is used (`read.f90:13690`).
- `ASSIGN_PRESSURE_ZONE` floods within one mesh. It passes through gas and removable solids and stops at non-removable solids and thin OBST faces (`func.f90:5558-5729`).
- `ZONE_BOUNDARY_EXCHANGE` crosses only INTERPOLATED boundaries and sets zone 0 at OPEN boundaries (`main.f90:2836-2900`; `2879-2881`, `2889-2894`).
- Overlapping explicit zones stop with ERROR(872) (`main.f90:2636-2640`).
- Every leftover sealed region gets its own auto zone (`main.f90:2652-2680`).
- Gap faces are exterior walls, not interpolated boundaries, so no flood crosses a gap.

**Consequences.**
- Each sealed zone is one component. One ΔP_zone with one `D_PBAR_DT` (`divg.f90:1519`, added to the divergence at `1523-1540`) is exactly that component's compatibility condition.
- Zones merged at run time merge only through gas (`divg.f90:1295-1321`), so a merged group is still one component.
- Zone 0 can span several components, but each of them has an OPEN (Dirichlet) face, so none is singular. Zone 0 has no `D_PBAR_DT` (the loops start at `IPZ = 1`, `divg.f90:1505, 1518`).
- Exception: level-set mode sets `NO_PRESSURE_ZONES` (`read.f90:1981-1984`; `main.f90:2620-2626`). All cells are then zone 0 with no flood fill. A sealed gap component with no OPEN face would be singular and have no `D_PBAR_DT`. The per-component rule below still pins it and removes its mean (not run).

**What goes wrong if a zone did span two components** (measured demonstration):
- One mean over both components leaves each component incompatible.
- MLMG still reports convergence (resid/bnorm ≈ 1e-11). The pinned rows absorb the whole imbalance:

  | Case | True max residual | Per-component residual sums |
  |---|---|---|
  | `split` / `wall` | 3.3e3 · max\|RHS\| | ±6.6e3 (synthetic units) |
  | masked `device_restart_a` walls | 209 · max\|RHS\| | ±358 (synthetic units) |

- The pin hides the error, so it has to be excluded by construction; MLMG will not detect it.

**Rule.**
- Mean removal, pin and gauge apply per singular operator component, whatever the zone ID. `D_PBAR_DT` stays per zone.
- Setup assertion: every sealed zone (ID ≥ 1) is exactly one operator component, and every component without an OPEN face lies inside one sealed zone. It cannot fail for FDS input; it guards the AMR zone and mask code.
- No extra constraint is needed. The mask is static, so components are recomputed only when zones merge, as in FDS.

### 5.4 (d) Masked single-level MLMG against FDS: PASS on the tested cases

**FDS reference** (measured; FFT with default `&PRES`; `VELOCITY_ERROR_FILE`; np = 2):

| Case | Velocity tolerance | Pressure iterations per half-step | Final max velocity error | Half-steps ending above tolerance | Max pressure residual (tolerance) |
|---|---|---|---|---|---|
| `device_restart_a`, 0–40 s, 463 half-steps | 0.2 m/s | 1 in every half-step | 0.083 m/s | 0 | 1.35 s⁻² (125) |
| `hallways`, 0–10 s, 2,341 half-steps | 0.03125 m/s | mean 6.0; 932 hit the cap of 10 | 0.104 m/s | 911 | 1.36e4 s⁻² (5.12e3) |

In every recorded iteration (463 and 14,095 of them) the cell with the largest velocity error lies on a mesh-boundary face. It is inter-mesh error.

**Masked MLMG** (same gap geometry on the bounding box; OPEN-known cells; synthetic RHS; np = 2; relative tolerance 1e-10):

| Case | Configuration | MG levels | Iterations | True max residual / max\|RHS\| | Solve time |
|---|---|---|---|---|---|
| `device_restart_a`, OPEN | overset + HYPRE bottom | 1 (nz = 7 cannot coarsen) | 3 | 4.8e-13 | 0.08 s |
| `hallways`, OPEN | overset, geometric MG, BiCGStab bottom | 5 | 20 (14 to 1e-6) | 4.4e-11 | 0.15 s |
| `hallways`, OPEN | overset + HYPRE, max coarsening level 0 | 1 | 3 | 8.6e-12 | 1.2 s |
| `hallways`, sealed (one pin) | overset + HYPRE bottom | 1 | 3 | 1.3e-9 | 0.80 s |

**Why this passes.**
- A single-level solve has one velocity per face, so there is no inter-mesh velocity error at all. Gap walls are exact (β = 0).
- The error FDS iterates on in these cases is therefore absent, and one solve replaces FDS's 1–10 iterations.
- The residual meets FDS's pressure tolerance for any RHS magnitude below about 1e14 s⁻² (derived: 5.12e3 / 4.4e-11).
- OBSTs inside meshes keep the E-1 iteration (N3). The masked solve does not change that iteration count.

**Not covered.** Composite masked solves, run time against NFR-030 (O6), and iteration counts on an FDS-generated RHS. The solver is linear, so the residual ratios carry over.

### 5.5 Required configuration for the masked branch (amends N3)
1. Overset mask 0 on gap cells, OPEN-known cells and pinned cells. Set a > 0 in every mask-0 cell (for example a·α = 1/dx²); the value has no effect under the mask but avoids the zero-diagonal NaN.
2. Level-0 bottom solver: HYPRE BoomerAMG through the overset-aware IJ path. Geometric coarsening then runs as far as the mask allows. Where a pin or unaligned gap walls stop it, HYPRE solves the level directly. That is equivalent to the `setMaxCoarseningLevel(0)` + HYPRE reference of 01 §G.2, which is the closest to FDS UGLMAT.
3. Do not use `setNSolve`. Use a-coefficient pins (the Fortran fallback of 01 REC-A1b) only with the maximum coarsening level set so that every β = 0 face and gap boundary lies on a 2^L boundary.
4. Components:
   - find them with a flood fill over faces with β ≠ 0; pin deterministically (lowest global cell index);
   - remove the D-032 volume-weighted mean per component, using the exact sum (FDS uses the arithmetic mean over unknowns, `pres.f90:3385-3404`; the two are identical on a uniform level 0);
   - after the solve, apply a per-component gauge as in `pres.f90:3495-3540`;
   - run the FR-039 true-residual check with the pinned rows included.
5. Add the setup assertion of §5.3.

## 6. Amendment (2026-09-26): level-0 padding and the N6 fallback

Inputs: `amrex/n6-stairwell-memory-estimate.md` (Integration Lead) and `vv/gap_face_vents.md` (A-46). The pressure results of §5 are unchanged.

### 6.1 N1a. Level-0 domain padded for divisibility

**Problem.** `stairwell`'s bounding box is 189×173×549 cells. No direction divides by 2, so AMReX accepts only blocking factor 1 (`AMReX_AmrMesh.cpp:1252-1259`), MLMG cannot coarsen level 0, and refined grids have no alignment.

**Ruling.**
1. For a non-box level 0, the AMReX level-0 domain is the bounding box of the level-0 `&MESH` blocks, extended in each direction to the smallest multiple of the level-0 blocking factor bf0. For `stairwell` with bf0 = 8 this gives 192×176×552.
2. Padding goes on the high side. It goes on the low side of a direction whose high bounding-box face carries an `OPEN` or PERIODIC vent. If both faces of a direction carry one, that direction is not padded; it runs with blocking factor 1 in that direction and FR-010's coarsening warning.
3. Padding cells are gap cells (N1): static solid, never refined. A face between gas and padding is a gas-gap face under N2. FDS treats the same face as a mesh face with no neighbour and gives it the same default surface, so the discretisation does not change. Vents on these faces apply as in FDS. The normal flux of a non-`OPEN` vent enters the pressure right-hand side as on any wall.
4. Padding is internal and is by whole level-0 cells, so every grid plane is unchanged. Anything FDS derives from mesh extents uses the unpadded `&MESH` union: the `P_0` table (§5.2), `XB` snapping of OBSTs, vents and devices (D-009), and the output meshes (ADR-004, written per `&MESH`). IR-008 blocking-factor alignment is measured from the padded domain's low corner.
5. Box domains are never padded. They keep the FFT path and the existing FR-010 rules.

**Rationale.** It restores MLMG coarsening and bf0-aligned refinement at no measurable cost: 19.5 GB padded against 19.7 GB unpadded for covering, and 3.15 against 3.19 GB for the fallback (the padded layout aligns better with 16³ boxes).

**Rejected.**
- Unpadded level 0 with blocking factor 1: no MG coarsening (the bottom solver runs on 18 M cells) and no alignment for refined grids.
- Padding both sides symmetrically: it moves two faces instead of one, for no gain.
- Asking users to remesh: FR-006 runs the 40 cases unmodified.

**[VERIFY] (V&V Lead, script check).** List any of the 40 cases whose padded direction carries `OPEN` or PERIODIC vents on both bounding-box faces. None is expected.

### 6.2 N6 fallback: level 0 built from the mesh union

**Problem.** Dropping only the all-gap boxes (the earlier N6 wording) keeps gap cells inside the kept boxes, and those cells need the overset mask. MLMG asserts that a non-covering level 0 has no overset mask (`AMReX_MLABecLaplacian.H:847-850`, `AMReX_MLPoisson.H:233-236`). The assert is compiled out of release builds. In release, MLMG then marks the problem singular even with a pinned cell and applies one global singular fix. It needed 13-18 iterations instead of 3 and matched the covering solution only after mean removal. With several components that one global fix is wrong in principle.

**Ruling.**
1. An overset mask on a non-covering level 0 is not allowed in any build. We do not ship on behaviour that our own debug builds reject.
2. The fallback level-0 BoxArray is the union of the level-0 `&MESH` blocks, chopped by `max_grid_size` and installed through `AmrMesh::PostProcessBaseGrids` (`AMReX_AmrMesh.H:412-421`). The Geometry domain stays the §6.1 padded bounding box. Level 0 then holds no gap cells, so there is no gap mask and no gap pin.
   - MLMG applies homogeneous Neumann on the uncovered boundary (`AMReX_MLLinOp.H:251-266`, `1238-1240`, `1514`). That is exactly the gas-gap no-flux of N2/N3.
   - The normal flux of a non-`OPEN` gap-face vent enters the right-hand side as on any wall.
   - `OPEN` vents on bounding-box faces use the domain-face treatment of the pressure spec (per-ghost-cell Robin coefficients, `pressure/01-amr-mapping-spec.md` §D table), not the overset mask. Its [VERIFY] applies.
   - FR-039's third branch covers both forms (covering masked, and mesh-union). The form is fixed at setup and never switches during a run.
3. The fallback is used only when all of these setup checks pass. Otherwise the case runs covering or gets a recorded exception (N7):
   - (a) No `OPEN` vent on a gas-gap face. `stairwell` passes: both its vents are on bounding-box faces (A-46). §6.1 padding moves its `Extract` vent (not `OPEN`) onto a gap face, which is allowed.
   - (b) No overset mask anywhere on level 0. OBSTs inside the meshes use E-1 (IBM forcing on the all-cell operator), not E-2 masking, while the fallback is in use. E-1 is the first step of REC-E1 anyway. E-2 on fallback cases waits for MLMG support of overset masks with an uncovered level 0.
   - (c) At most one singular operator component (a component with no Dirichlet face). MLMG's single global singular treatment is then correct, together with D-032 mean removal and the per-component gauge (§5.5 item 4). The flood fill of §5.1 counts the components.
   - (d) The Integration Lead's trigger: covered cells are below 30 % of the padded bounding box, and the covering estimate exceeds 75 % of the usable memory U. U is 12 GB on the workstation (NFR-031); on a GPU it is 3/4 of device memory divided by the ranks that share it. On the workstation only `stairwell` meets this among the 40 cases.
4. Driver obligations.
   - Ghost cells of level-0 boxes that fall outside the union lie inside the domain, so neither `FillBoundary` nor `PhysBCFunct` fills them. After every FillPatch the driver sets them to the gap-wall ghost values defined by the N2 wall records, which are the same values the covering design holds in its gap cells.
   - Nesting for level 1 and above comes from the level-0 grids (`AMReX_AmrMesh.cpp:631-644`). IR-008 already keeps refinement inside the meshes.
5. Box size. Boxes follow the mesh extents and can fall below the 16-cell guidance (`requirements.md:443`; `stairwell` gets 13-14-cell boxes). This is allowed for fallback cases as a recorded NFR exception. Level-0 MG coarsening is limited by the mesh extents. On the CPU path the HYPRE bottom solver (§5.5 item 2) covers the rest; the GPU path is item 7 (iv).
6. Memory. `stairwell` level 0 is about 1.7-2.0 GB (Integration Lead estimate), well inside U.
7. **[VERIFY] in the Phase 2 prototype** (Integration Lead, with the Pressure Solver Lead):
   - (i) On a low-fill case small enough for both forms, the mesh-union solve matches covering masked MLMG within eps_H after mean removal.
   - (ii) A debug build runs clean, with no assert.
   - (iii) Level-1 grids tagged over level-0 boxes that are not bf0-aligned pass proper nesting and IR-008 alignment.
   - (iv) Bottom-solver convergence on the GPU path, where our HYPRE builds are CPU-only (ADR-002).

**Rejected.**
- Dropping only the all-gap boxes and keeping the mask (the earlier N6 wording). It relies on an assert being compiled out, and its singular treatment is wrong with several components.
- Patching or relaxing the AMReX assert. That is an upstream change outside our control, and the release behaviour it would bless is itself wrong with several components. It can be revisited if upstream adds support for an overset mask with Neumann on the uncovered boundary. The Integration Lead may ask upstream, but nothing depends on the answer.
- `a` = 1 in gap cells without a mask. MLMG aborted in the tiny run ("MLMG failed").
- A covering level 0 for `stairwell` on the workstation: 15-20 GB against 12 GB.
- A run-time exception as the first choice. The fallback is cheap and keeps the case in scope. The exception remains for cases that fail (a)-(c).
