# Flux read-out and override hooks: design (task B, S9.4)

Status: design note for the Architect and Role 3, written before any code. Basis: `Source/regrid_transport/notes/flux-override-interface.md` (Role 3 request, sections
numbered below as "RO §n"), `RegridInterface.H` `FluxAccess`, D-050 (single global dt, no reflux, the coarse face flux is overwritten by the fine flux at every stage),
the ruling "no heat override in phase 3". Line numbers are those of `Source/mass.f90` and `Source/divg.f90` at the current HEAD.

## 1. What is read and what is overridden

| Kind | Quantity (unit kg/(m2 s) of rho*Z_n) | Where it exists in FDS | Hook |
|---|---|---|---|
| ADV | stage product `FX(I,J,K,n)*UU(I,J,K)` (and `FY*VV`, `FZ*WW`), faces I = 0..IBAR | not stored (inline in the mass update, `mass.f90` 448-450 predictor and 625-627 corrector); `ADV_FX` exists only with `STORE_SPECIES_FLUX` and after the corrector holds the stage average (`mass.f90` 654) | driver side only: the density kernel is the driver's generated copy `fds_density_split.f90` (`tools/gen_density_split.py`), so the hook is a change of the generator, **no patch of `mass.f90`** |
| DIF | `RHO_D_DZDX/Y/Z(I,J,K,n)` (`SWORK1..3`) after the species-sum fix (`divg.f90` ~249) and the wall corrections (~196-225), i.e. exactly what the divergence reads (`DEL_RHO_D_DEL_Z`, ~411) | `DIF_FX` exists only with `STORE_SPECIES_FLUX` (`divg.f90` 269-278, 336-349) and holds the NEGATIVE of `RHO_D_DZDX` | **guarded patch 0008** (`divg.f90`, one call under `#ifdef WITH_AMREX`) |
| heat | enthalpy diffusion flux `H_S*RHO_D_DZDX` (`divg.f90` ~303) and conduction `KDTDX` (~516) | formed inside the species loop after the species flux | no override now (RO §2, ruling): see section 7 |

Sign and unit convention of the read-out and of the override value: the value the divergence loop reads (`RHO_D_DZDX` itself, not `DIF_FX`). ADV is the product with the
UNMATCHED face velocity, as FDS has it in `DENSITY` (the driver copy restores `UVW_SAVE` on interface walls, `fds_density_split.f90` 86-91, 181-186). A fine box that has
a coarse-fine face has no wall there (the driver's `iface()` turns every non-domain face into a no-wall face), and a coarse box has no wall under or next to a fine patch,
so the FDS mesh-interface branches (`wall.f90` 891 `COARSE_MESH_IF`, `divg.f90` 204-216 `EWC%NIC>1`) are never reached at a level interface (answer to RO §6.3).
Fine boxes get `NIC = 1` interface walls (`level-interface.md` §3).

## 2. Registered arrays (answer to RO §6.1: yes, no `STORE_SPECIES_FLUX` side effects)

`StageFlux` per level, owned by `LevelRegistry`: for kind in {ADV, DIF} and dir in {x,y,z} one nodal `MultiFab` (`convert(valid_box, nodal dir)`, ncomp = N_TOTAL_SCALARS, 0 ghost,
same BoxArray and DistributionMapping as the level), allocated on first use of `compute_stage_fluxes` and remade with the level (retired objects as for `Fields`). AMReX face
`a = lo + I` is FDS face I (the HIGH face of cell I): the Fortran view of a FAB uses lower bounds `(0, 1, 1)` for x (resp. `(1,0,1)`, `(1,1,0)`) and goes through the same
`FDS_HOOK_SET_VIEW` route as the box views (rank-4 variant: scalar index 1..N last). `STORE_SPECIES_FLUX` stays `.FALSE.`: no `ADV_F*`/`DIF_F*` arrays are allocated by `init.f90`.
`stage_flux(level, kind, dir)` returns these. Memory: 6 nodal arrays x N_TOTAL_SCALARS doubles per cell, about 6/1 of one scalar field per scalar (tens of MB for 128^3, 1 scalar).

## 3. Phases and where each is executed (three phases of RO §4, finest first)

`FluxAccess` is implemented by a `TimeLoop` adapter (`FluxStages`, new, over the existing per-level stage entry points). For stage `predictor`:
1. `compute_stage_fluxes(level, predictor)` for every level:
   - ADV: after `stage_mfd` (FX computed by the existing `MASS_FINITE_DIFFERENCES` call) a new generated entry `DENSITY_ADV_READOUT(NM)` builds the same `UU,VV,WW` as `DENSITY_PRE_CLIP` (including the `UVW_SAVE` restore) and writes `FX*UU` etc. into the registered ADV arrays through the hook module. It reads only; the cell update does not run.
   - DIF: `stage_divergence1` as today, with the hook (patch 0008) in READ mode: at the hook point the hook copies `RHO_D_DZDX..Z(0:IBAR,..,1:N_TOTAL_SCALARS)` of the box into the registered DIF arrays and returns. The rest of `DIVERGENCE_PART_1` runs as at level 0.
2. Role 3 reads `stage_flux` and calls `set_flux_override(level, kind, per_local_box)` (coarse faces only, `FluxOverride` of RO §3), finest first. The driver stores the lists in a hook-module table indexed by box (Fortran view `n_ovr`, `idx(3,n_ovr)`, `val(N_TOTAL_SCALARS, n_ovr)`; index conversion `I = a - lo_x` for the normal direction, `J = b - lo_y + 1` tangential, done once at `set_flux_override`). Every entry is range-checked (inside `0..IBAR`, `1..JBAR`, `1..KBAR`; a ghost face is an error).
3. `apply_flux_divergence(level, predictor)` for every level, finest first:
   - ADV: the density update of `DENSITY_PRE_CLIP`/`DENSITY_POST_CLIP`. With `n_ovr(box) = 0` the generated code is the present loop text, unchanged (the loop is selected by `IF (N_ADV_OVR>0)`); with overrides a second copy of the loop reads a per-face product array that equals `FX*UU` except at the listed faces (a separate loop because a materialised product need not round like the inline expression if the compiler contracts to FMA; the empty-set path is therefore the original loop and bitwise unchanged by construction).
   - DIF: for a level with a non-empty DIF override set, `DIVERGENCE_PART_1` is run again for the level with the hook in OVERRIDE mode (replace the listed `RHO_D_DZDX` values for all species at the hook point, before the heat loop and before `DEL_RHO_D_DEL_Z` is formed). A level with an empty set (finest level, or the whole run with no override) is not touched: its phase 1 result is final. Cost: the divergence kernel runs twice on levels that have a finer level. This relies on `DIVERGENCE_PART_1` being re-runnable with identical inputs (test T3). Optimisation if needed later: a derived split copy of `DIVERGENCE_PART_1` at the hook point (as for `DENSITY`), which would avoid the repetition.

Because the species flux override happens before the heat loop at `divg.f90` ~303, the enthalpy diffusion flux `H_S*RHO_D_DZDX` of an overridden face is automatically `h_face(coarse T) x overridden flux`, which is what FDS does at its own mesh interfaces (RO §2); conduction `KDTDX` is not overridden.

Level order in the driver: phase 1 runs coarse to fine (any order, no dependency); phases 2 and 3 run fine to coarse. The cf-ghost hook has filled the fine ghost cells before phase 1, so a fine face next to the interface is computed with the interpolated ghost like any other face (RO §2 "coarse value injected as ghost").

## 4. Answers to the open questions of RO §6
1. Registered arrays without side effects: yes, section 2.
2. Shared faces between two boxes of one level: both boxes compute the same value bitwise today for DIF (`.5*(a+b)*dz*rdxn` is symmetric and the ghost values are bitwise copies) and for ADV (the stage kernels are decomposition independent bitwise, `run_decomp_check.sh`; `tests/check_shared_face_flux.py` compares the two copies of a shared face directly). An override entry is applied by both boxes with the same value if Role 3 lists the face in both boxes; a face listed in one box only would make the two copies differ, so the driver (check at `set_flux_override`) requires the list of a shared coarse face to appear in every box that holds it, and aborts with a message otherwise.
3. The `COARSE_MESH_IF`/`NIC>1` branches are skipped at level interfaces (section 1); `UVW_SAVE` on a fine box's interface wall is saved before the match as at level 0 (`level-interface.md` section 4); it must hold the interface-face velocity (coarse = area average of the fine faces), which is a requirement on Role 3's ghost hook (`CfGhostRequest` result) that the driver will check in the T2 setting, so that matched and unmatched velocity coincide.
4. Ghost faces never feed the divergence: the loops read faces `I-1` and `I` for `I = 1..IBAR`, i.e. faces `0..IBAR` (`mass.f90` 448-450, `divg.f90` 411-413). Overrides on faces outside that range are an error.

## 5. Bitwise-unchanged guarantee and tests
- T1 (empty set): a run of the driver with the adapter's three phases and no overrides is bitwise equal to the same run without the adapter (final fields and step logs, `bitcmp.sh` style), 1 and 4 ranks, `dec2`, `dec4`, a case with `STORE_SPECIES_FLUX` off.
- T1b: `USE_AMREX=OFF` with patch 0008: `check_off_bitwise.sh` (the macro is undefined, so the preprocessed source is unchanged apart from blank lines).
- T2 (no-op override): an override list that carries the kernel's own values (read from phase 1) at the box-box faces of one level must reproduce the empty-set result bitwise (ADV and DIF), 1 and 4 ranks. Tests the index conversion and the second loop copy of the density update.
- T3 (re-run): `DIVERGENCE_PART_1` run twice per stage gives bitwise the same fields as once (needed for the DIF phase 3).
- T4: read-out content: `stage_flux(ADV)` equals `FX*UU` recomputed on the C++ side from `FX` (read through `fds_k_xfer` mode 3) and `U`; `stage_flux(DIF)` equals `-DIF_FX` of a run with `STORE_SPECIES_FLUX` on (reference only; that run is not the production configuration).
- T5: sum rule: for the DIF read-out the sum over the tracked species of each face is zero to rounding (species-sum fix), so a conservative override keeps it.

## 6. Patches and files
- 0007: DRAFT grown `POINT_TO_BOX` (D-056 option B), `patches/0007-*.md`; not part of this task, listed because the fine-level hook needs boxes of level > 0.
- 0008: `divg.f90`, one `#ifdef WITH_AMREX` block in `DIVERGENCE_PART_1` after the species-sum correction (`CALL FDS_HOOK_DIF_FLUX(NM, PREDICTOR, ...)`, plus a `USE FDS_AMREX_HOOKS` line). The hook returns at once when the box has no registered DIF array and no override: with the macro undefined the patch changes nothing.
- Driver files: `fds_flux_hooks.f90` (hook module: registered views, override tables, READ/OVERRIDE state), `tools/gen_density_split.py` (read-out entry and the override copy of the update loop), `FluxStages.{H,cpp}` (adapter over `TimeLoop` stage entries, implements `FluxAccess`), `LevelRegistry` (StageFlux arrays), tests T1-T5.

## 7. Heat flux, second step
Only if the V&V lead later asks for matching of conduction or of `h x F` per fine face: a second hook at `divg.f90` ~303 (the `H_RHO_D_DZDX` arrays, `WORK5..7`) with a `FluxKind::HEAT` and the same list structure, and at ~516 for `KDTDX` (`WORK1`). Not designed further: the ruling is that heat follows FDS.

## 8. Dates (calendar days, estimates, not commitments)
Design (this note) 6 October. ADV read-out and override with T1, T2, T4 (driver side only): 12 October, behind the level > 0 kernel work of D-056 which has priority. Patch 0008, DIF read-out and override with T1b, T2, T3, T5: 14 October. The real 2-level exercise waits for Role 3's `set_flux_override` calls and for oneAPI validation of 0005/0007 for fine boxes; at level 0 the tests above run with fake overrides.
