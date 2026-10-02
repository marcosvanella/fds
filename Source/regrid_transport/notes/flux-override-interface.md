# Interface flux overwrite: flux read-out and override input (request to the driver)

Status: proposal from Role 3 for Role 1, reviewed by the Integration Lead. Basis: D-050 / ADR-002 v1.2 (single global dt, no
reflux; the coarse face flux is overwritten by the fine fluxes at every stage). Plain-language terms: a *face flux* is the amount of
a quantity crossing one cell face per unit area and time; the *coarse-fine interface* is the set of coarse faces on the outline of
a fine patch; "overwrite" means the coarse side uses, for those faces, the fine-side value instead of its own.

## 1. What we compute (Role 3)
For every coarse face on the interface of a level l+1 patch, in every stage (predictor and corrector), for each scalar n:
`F_coarse = (1/A_c) * sum over the fine faces f tiling that coarse face of ( A_f * F_f )`.
With uniform Cartesian cells every fine face has the same area, so this is the plain mean of the `r1*r2` fine face values
(`r1,r2` = refinement ratios in the two directions tangential to the face; 1 in a hidden single-cell direction). F is a flux per area,
the *product* face value times velocity, not a product of averages. Levels are handled finest first. Role 3 also handles the
exchange between ranks; Role 1 only sees local boxes.

## 2. Read-out: stage face fluxes (Role 1 provides)
Per box and per direction d in {x,y,z}, a face-centred array, scalar index n = 1..N_TOTAL_SCALARS as the component, no ghost layer:
- **ADV**: the stage product used in the mass update, `FX(I,J,K,n)*UU(I,J,K)` (and FY*VV, FZ*WW). Not the stored `ADV_FX` after the
  corrector, which is the average of the two stages (mass.f90:654). Unit kg/(m2 s) of rho*Z_n. On interface faces `UU` is the
  *unmatched* face velocity: `DENSITY` restores `UVW_SAVE` there (mass.f90:421-434, 593-607; driver copy fds_density_split.f90:86-91,
  181-186), while the momentum terms use the matched value (`MATCH_VELOCITY`, velo.f90:2630). In the multi-level path an interface
  face has one velocity (coarse = area average of the fine faces), so `UVW_SAVE` must equal that velocity; matched and unmatched then coincide.
- **DIF**: the species diffusive face flux `RHO_D_DZDX(I,J,K,n)` (and Y, Z) after the divg.f90 species-sum correction and the wall
  corrections, i.e. the values the divergence loop reads. Same unit. FDS already overwrites the coarse diffusive species flux with
  the area-weighted sum of the fine two-point fluxes at its own mesh interfaces (wall.f90:891-958, `EWC%NIC>1`; divg.f90:204-216; TRG
  Mass_Chapter.tex:234). Our override must give the same numbers (fine face flux with the coarse value injected as ghost, species-sum
  fix as wall.f90:955-959). That FDS branch reads `OMESH` data the AMReX port does not fill, so Role 1 bypasses it at level interfaces. (`DIF_FX` in the registry already exists; confirm it holds
  these values and is filled in AMR mode without `STORE_SPECIES_FLUX` side effects.)
- **HEAT (optional, later)**: FDS does not match conduction `KDTDX` at refined interfaces (skipped for INTERPOLATED faces, divg.f90:~540)
  and forms the enthalpy diffusion flux as face enthalpy x the overwritten species flux (divg.f90:~303). Phase 3 reproduces FDS: no
  HEAT override. If wanted later: `KDTDX` and `H_RHO_D_DZDX` per species n at the same positions. Not needed for mass and species.

Indexing. FDS face I in direction x is the HIGH face of cell I, valid I = 0..IBAR, with J = 1..JBAR, K = 1..KBAR (analogous for y, z).
In the registry's AMReX index (Fields.H index map) the same face is face `a = lo_x + I`, the LOW face of cell `a`; the array is nodal
in direction d, box `convert(valid_box, nodal d)`, so FDS I = a - lo_x. The array is per box (one FAB per box of the level), same
BoxArray and DistributionMapping as the level; shared faces between two boxes of one level are held by both boxes.

## 3. Override input (Role 3 provides, Role 1 applies)
A sparse face list per box (the interface is a surface, so a mask array would be mostly empty):
```
struct FluxOverride {            // one set per box, per stage kind (ADV, DIF, later HEAT)
  int dir;                       // 0,1,2: direction of the face normal
  std::vector<std::array<int,3>> face;   // AMReX face index (low face of cell a), box-local list, sorted by (k,j,i)
  std::vector<double> value;     // size face.size()*nscal, face-major: value[f*nscal + n-1]
};
```
Fortran view per box and direction: `n_ovr`, `idx(3,n_ovr)` in FDS box indices (I = a - lo_x for the normal direction, J = b - lo_y + 1 for
the tangential ones), `val(N_TOTAL_SCALARS, n_ovr)`. Each entry replaces the stage product for all scalars at that face. Every face
in the list lies in the box's valid face range and on a coarse-fine interface (the list holds coarse faces only).

## 4. When it is applied
- **ADV**: in the mass update after `MASS_FINITE_DIFFERENCES` and before the divergence of the flux inside `DENSITY` (predictor and
  corrector). The kernel replaces `FX*UU` by `val` at listed faces; the `R(I)`, `RDX`, `RRN` factors are unchanged.
- **DIF**: in `divg.f90`, after the species-sum correction (line ~249) and the wall corrections (~200-225), before `DEL_RHO_D_DEL_Z`
  (~411) is formed. Because every fine face has species sum zero there, the overwritten coarse face also has it.
- **HEAT**: before `DP` collects `DIV_DIFF_HEAT_FLUX` (~394) and `DELKDELT` (~564).
- The override must be set before the stage and cleared (or replaced) at the next stage. The covered coarse cells are replaced by the
  average-down afterwards; the overwrite only has to be right for the uncovered coarse cells next to the interface.

## 5. Bitwise-unchanged guarantee
With `n_ovr = 0` for every box and direction, the kernels must give results bitwise identical to the current ones (same operation
order, no extra arithmetic). Suggested form: the override loop sits behind `IF (N_OVR>0)` and writes into the product only at listed
faces; for ADV the product is computed exactly as `FX(I,J,K,N)*UU(I,J,K)` in both cases. Test: `check_off_bitwise.sh` with an empty set,
and decomposition independence with a non-empty set on a single-level case where the override equals the kernel's own value (must be
a no-op bitwise).

## 6. Questions for Role 1 and the Integration Lead
1. Is the stage product available per stage as a registered array (no `STORE_SPECIES_FLUX` side effects), or should it be added as
   an extra output (cell arrays unchanged)?
2. Shared faces between level-0 boxes: do both boxes compute bitwise the same value today? The override relies on both using the same
   replacement value, which is guaranteed only for listed faces.
3. Confirm the `COARSE_MESH_IF` branch (wall.f90:891) is skipped, and that `UVW_SAVE` is set to the face velocity, at level interfaces.
4. Do `ghost` faces (index -1 and IBP1) ever feed the divergence? We assume not, so overrides never target them.
