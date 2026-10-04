# Interface flux overwrite: how a stage is wired (Role 3 side, R2b)

Plain words: a *face flux* is the amount of a quantity crossing one cell face per unit area and time. At the edge of a fine patch the coarse faces take the
area-weighted mean of the fine face fluxes that tile them (the "overwrite"); with one time step on all levels this is all that is needed, no flux register.

## Code
- `FluxOverrideOps` builds the lists for one coarse level from the fine level's face fluxes (periodic wrap, shared faces listed in every box that holds them).
- `FluxStageRunner::set_overrides(FluxAccess&, kind)` does it for every level, finest first, and ALWAYS calls `set_flux_override` on every level (the finest gets an
  empty set, which clears the lists of the previous stage). `overwrite=false` clears everything: the negative control and the FDS interface behaviour of the FR-016 baseline.
- `GhostShare` counts the covered coarse cells that serve two or more ghost faces (D-059); the number must be printed at each regrid (the regrid wiring owns that call).
- `tests/mock_transport.H` is a small conservative transport code with the stage structure of one FDS SSPRK2 step, used to test all of the above without the driver's
  fine-level binding (test `test_flux_stage`, 1 and 4 ranks): composite mass conserved to ~1e-17 with the overwrite ON, 1e-3..1e-4 drift with it OFF; uniform state stays
  uniform; three levels; 3-D patch with a lug; patch on periodic edges; full-domain level 1 bitwise equal to the single-level run; override equals the independent mean and
  keeps the species sum zero.

## Stage order the driver loop must follow (per stage, predictor and corrector)
1. `stage_viscosity` (COMPUTE_VISCOSITY + MASS_FINITE_DIFFERENCES) on all levels
2. `flux_readout_adv` on all levels            -> `runner.set_overrides(fs, Adv)`
3. `stage_density` on all levels (consumes the ADV lists)
4. `stage_exchange(code 1 or 4)` and `stage_boundary` on all levels, coarse first (the coarse-fine ghost hook runs here), `stage_velocity_flux`, `stage_wall_bc`
5. `flux_readout_dif` on all levels (= DIVERGENCE_PART_1)   -> `runner.set_overrides(fs, Dif)`
6. `flux_apply_dif` on all levels (re-run where a list exists)
7. average-down of the cell fields (`average_down_hierarchy`)

## Finding for Role 1 (FluxStages as committed)
`FluxStages::apply_flux_divergence(level)` runs `stage_density` and then immediately `flux_readout_dif(level)` (DIVERGENCE_PART_1). In FDS order DIVERGENCE_PART_1 comes
AFTER the exchange of RHOS/ZZS (steps 4), and with several levels, finest-first, a fine level's DIF read-out would run before the coarse level has its updated density and
before the coarse-fine ghost hook has filled the fine ghost cells of that stage. The DIF read-out of a level must therefore not be part of `apply_flux_divergence`; it
belongs to step 5 above. Either split `apply_flux_divergence` (density only) and move the DIF read-out into `run_divergence_part1` / a new phase, or the stage loop calls the
primitive `TimeLoop` entry points directly as listed. Role 3's runner only needs `stage_flux` and `set_flux_override` from `FluxAccess` and works with either.

## Requirements on the velocities
The advective override value is face value x the single interface-face velocity. The coarse interface-face velocity must equal the mean of the fine faces that tile it
(prescribed fields: build the coarse field by average-down of the fine face velocities, or from a stream function sampled at the common nodes, as the mock test does).

## R2b: wired against the committed hooks (S12 FluxStages, patches 0007-0009) — status: run against the validated patch set (0005-0009, commit 1f3f4de4c4) with gfortran 14 and oneAPI ifx/Intel MPI; same numbers
Role 1 split the stage as asked (`apply_flux_divergence` = density only, `readout_dif` / `run_divergence_part1` separate), so the runner is used as designed; no change to `FluxStageRunner`.
`DriverModes.cpp` (`fds_amr <case> --rt-e2e transport|ghost|list`, entered through the hook of `docs/upstream-patches/0006-r2b-driver-e2e-hook.patch`) runs one transport-only step as
predictor then corrector, each: `stage_state`; `stage_viscosity` (all levels, coarse first); `compute_stage_fluxes` (ADV read-out); runner ADV overrides; `apply_flux_divergence` (finest first);
`stage_exchange` (code 1 / 4: the coarse-fine ghost hook of `install_cf_ghost_hooks` runs here) and `stage_boundary`; `stage_velocity_flux`; `stage_init_divergence`; `stage_wall_bc`;
`readout_dif` (all levels); runner DIF overrides; `run_divergence_part1` (finest first). The velocity is prescribed and constant (no pressure solve); a level 1 is made with the registry transfer
and `bind_level(1)`, and `install_cf_ghost_hooks` runs after `bind_level` (the other order loses the hook). One corrector-form pass without the density update (`prime`) fills the lagged
DEL_RHO_D_DEL_Z / D before the first step. There is no average-down inside the step (the FDS ghost-rule values of the covered cells stay intact); the composite sums skip covered cells.

What the runs showed (gfortran and oneAPI agree to the digits printed; 1 and 4 ranks agree):
- FDS skips the density update while ICYC <= 1 (DENSITY_PRE_CLIP), so a driver that starts its own step counter must start the first step at cycle 2. A first version of this mode did not and
  the ADV read-out was identically zero; the symptom was a transport test that "conserved" because nothing moved. `Env::icyc` starts at 1 now.
- Advection only (species diffusivity set to 1e-12): overwrite ON conserves the composite species to 2e-14 over 6 steps, OFF drifts 4.3e-5 (E3a).
- With the default diffusivity: ON 4.4e-7, OFF 4.3e-5 (E3b, 6 steps, 100 times better but NOT round-off). **Cause found, not ours (Role 1, level binding):** a fine box of the 2D cases has
  one cell in y, so its y faces are domain-edge faces; the fine level has no wall cells there (listed gap "fine-box domain-edge wall cells" in `Source/driver/notes/level-binding.md`) and its y
  ghost layer is never filled (ZZ, RHO, RHO_D read 0 at J=0 and JBAR+1 in every stage). `DIVERGENCE_PART_1` then differences against it: RHO_D_DZDY(:,0,:) = -5e-5 on the smooth profile, a spurious
  species flux through the y faces (level 0 has MIRROR walls there and gives exactly 0). The interface overwrites are exact: per stage the composite sum of the advective divergence is 0 (to 1e-17
  relative), the coarse DIF override equals the fine read-out face by face (W/E/B/T sums agree to all digits), and the composite sum of DEL_RHO_D_DEL_Z equals the sum of the fine y-edge
  flux terms (-2.2e-7 per stage on the smooth profile). Why it looked like a corrector effect: the prime pass runs on the symmetric initial profile, where the y-edge terms cancel over the
  species-sum correction; once the profile has advected they do not, so the stage sums are 0 at the first predictor and grow after that (about u * D * dt^2 per step, quadratic in the step count).
  Experiments (smooth profile, u = 1, 6 steps): zeroing the fine y-edge diffusive flux in the divergence (scratch overlay, not committed) gives 1.3e-15 instead of 4.4e-7; the same case in 3D
  (triply periodic, fine box y range 4..11, so no edge) gives 5e-15, and a 3D fine box covering all of y (periodic images fill the ghosts) 4e-14. Tests: E3c (the 3D interior box, with overwrite OFF as negative
  control, 4.6e-6) and E1c below. What Role 1 needs: wall cells (or an equivalent mirror/periodic/open treatment, as level 0 has) for fine-box faces on the domain edge, or the fill of those
  ghost layers. No change in `regrid_transport` can fix it (the fine DIV1 reads the ghost itself); E3b keeps its loose bound until then and must be tightened to E3A_TOL.
- E1 (full-domain level 1 against the uniform fine run): RHO equal to 1e-15, ZZ differs by 5e-6 per step in every cell (4e-5 after 6 steps), also for pure diffusion with zero velocity. This is
  the same defect (the full-domain fine box in 2D has the unbound y edge, the uniform fine reference run has MIRROR walls): the 3D triply periodic version (E1c, 32^3 fine against 16^3 + full level 1)
  agrees to 3e-16 for RHO and both ZZ. The "KNOWN GAP" of `level-binding.md` (predictor DS differs by up to 0.4 % on a level with non-uniform species) is very probably the same cause (not rerun here).
  The 2D E1 keeps RHO to 1e-12 and ZZ to 1e-4 until Role 1 binds the edge; then tighten E1_ZZ_TOL.
- Level-1 D / DS fields read back from the registry hold a uniform -0.24 1/s (level 0: 1e-15) after DIVERGENCE_PART_1 on a uniform-temperature, uniform-pressure state: the composite zone
  term (D_PBAR_DT over the levels of a zone) is not bound yet (listed in `level-binding.md`); the first prime pass also leaves D as NaN on level 1. Not used by the species update.
- E2: a uniform state stays uniform (exact) across the interface for 6 steps. E5: corner blobs run, D-059 count printed (32 covered coarse cells serve one ghost face, 16 serve two faces; patch of 8 x 8 coarse cells).
- E4 (`ghost` mode, `ns2d_16_l0` with the dumps of the instrumented FDS run of `ns2d_16_int_1to2_refinement`, step 2 and 3): the ghost layers 1 and 2 written by a real `stage_exchange` agree
  with the FDS values to 5e-15 (RHO, ZZ, TMP, RSUM; corrector and predictor forms), except the hole-corner cells, reported separately: KRES at 4 hole-corner cells differs (42 %, the D-059 shared cells).
- Findings for Role 1: (1) the `IBAR_MAX`/`JBAR_MAX`/`KBAR_MAX` raise in `fds_fine_level.f90` (heap corruption in VELOCITY_FLUX when a fine box is larger than every level-0 mesh) is in the working tree now;
  (2) `FDS_G_FILL_OM` crashes for several fine boxes on one rank (`docs/upstream-patches/0007-r2b-fill-om-bounds.patch`); (3) `stage_boundary(1, 3|6)` aborts in `fds_p_save_uvw`
  ("not a level-0 FDS mesh"), so the velocity-matching codes are run on level 0 only here (the prescribed velocity is constant); (4) the level-1 D / DS and the NaN above.
- `tests/run_e2e_driver.sh <build>` runs E1 to E5 on 1 and 4 ranks (E4 with `RT_FR016_DUMP=<prefix>`), then Role 1's two-level prescribed-velocity test (`driver_unit_tests`, which uses our override lists).
