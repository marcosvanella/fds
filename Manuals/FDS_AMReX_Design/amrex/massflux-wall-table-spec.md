# Wall-table spec for the mass-flux wall nests (L0880, L0882, L0401, L0403)

Scope: the four wall loops that overwrite species face fluxes after `GET_SCALAR_FACE_VALUE`:

| Loop | Routine | Nest | Writes | Uses wall-face value? |
|---|---|---|---|---|
| L0880 | `MASS_FINITE_DIFFERENCES` | mass.f90:93-188 (inside `SPECIES_LOOP`, mass.f90:65-192) | `FX/FY/FZ(:,:,:,N)`, N=1..N_TOTAL_SCALARS | yes (`B1%RHO_F*B1%ZZ_F(N)` or 0) |
| L0882 | same | mass.f90:224-320 | `FX/FY/FZ(:,:,:,0)` | yes (`B1%RHO_F/MW_F` or 0) |
| L0401 | `SPECIES_ADVECTION_PART_1_NEW` | divg.f90:1021-1084 | `FX_ZZ/FY_ZZ/FZ_ZZ(:,:,:,N)` | no (off-wall overwrite only) |
| L0403 | same | divg.f90:1120-1183 | `FX_ZZ/FY_ZZ/FZ_ZZ(:,:,:,0)` | no (off-wall overwrite only) |

All line numbers are at the pinned baseline (36975d7). `NWE = N_EXTERNAL_WALL_CELLS`, `NWI = N_INTERNAL_WALL_CELLS`,
`NS = N_TOTAL_SCALARS` (policy rename, s5_markers.toml:15-18). Conventions follow docs/adr/drafts/gpu-generator-design.md,
section "WALL flat tables": table name = prefix + upstream component name; integers `integer(c_int)`, reals `real(eb)`,
masks `integer(c_int)` (nonzero = true); the driver fills them.

## 1. What each nest reads

Per wall `IW = 1..NWE+NWI`, skipping `WC%BOUNDARY_TYPE==NULL_BOUNDARY` (mass.f90:95, 226; divg.f90:1023, 1122):

1. Wall-face phase (L0880 mass.f90:108-129, L0882 mass.f90:239-261):
   - `WC%BOUNDARY_TYPE` against `SOLID_BOUNDARY`, `INTERPOLATED_BOUNDARY`;
   - ghost cell `(II,JJ,KK)`: `CELL(IC)%SOLID` and `CELL(IC)%EXTERIOR` (the "thin obstruction" test);
   - `BC%IIG,JJG,KKG,IOR` (the face written is `F*(IIG-1,..)` for IOR>0, `F*(IIG,..)` for IOR<0);
   - `B1%RHO_F`; L0880 also `B1%ZZ_F(N)`; L0882 `B1%ZZ_F(1:N_TRACKED_SPECIES)` for `MW_F`.
2. Off-wall phase (all four nests; mass.f90:133-186, 271-318; divg.f90:1029-1084, 1129-1180), only when
   `BOUNDARY_TYPE` is neither `INTERPOLATED_BOUNDARY` nor `OPEN_BOUNDARY`:
   - `BC%II,JJ,KK,IOR`;
   - guard `.NOT.(CELL(CELL_INDEX(cell))%WALL_INDEX(+-n)>0)` on the neighbouring gas cell:
     IOR=+1: cell (II+1,JJ,KK) flag +1 (mass.f90:140, divg.f90:1036); IOR=-1: cell (II-1,JJ,KK) flag -1 (mass.f90:150);
     IOR=+-2, +-3 the same on y and z. The six cases are the same pattern with flag index equal to IOR.

## 2. Fields required, with verification

| # | Name (shape, type) | Meaning | Status in the worktree | Upstream source |
|---|---|---|---|---|
| 1 | `W_BOUNDARY_TYPE(NWE+NWI)` int | `WALL(IW)%BOUNDARY_TYPE` | exists | type.f90:494 (`WALL_TYPE`), `BOUNDARY_TYPE=0` default |
| 2 | `BC_II, BC_JJ, BC_KK, BC_IIG, BC_JJG, BC_KKG, BC_IOR (NWE+NWI)` int | ghost cell, gas cell, orientation | exist | mass.f90:99-105; init.f90:3298 uses `-IOR` |
| 3 | `B1_RHO_F(NWE+NWI)` real | `BOUNDARY_PROP1(WC%B1_INDEX)%RHO_F` | exists (golden entry `B1_RHO_F:wall(NWE+NWI):real:in`) | type.f90:325 |
| 4 | `B1_ZZ_F(NWE+NWI,NS)` real | `B1%ZZ_F(1:N_TRACKED_SPECIES)` | **new** | type.f90:296 (allocatable, `1:N_TRACKED_SPECIES`); filled init.f90:3553 |
| 5 | `WALL_INDEX(0:IBAR+1,0:JBAR+1,0:KBAR+1,-3:3)` int | `CELL(CELL_INDEX(i,j,k))%WALL_INDEX(n)` | **new** | type.f90:2182 (`DIMENSION(-3:3)`, default 0); set only at init.f90:3298 |
| 6 | `EXTERIOR(0:IBAR+1,0:JBAR+1,0:KBAR+1)` int mask | `CELL(...)%EXTERIOR` | **new** | type.f90:2177 |
| 7 | `SOLID(0:IBAR+1,0:JBAR+1,0:KBAR+1)` int mask | `CELL(...)%SOLID` | exists (s5_markers.toml:38 flatten rule) | type.f90:2176 |
| 8 | face-write check (section 4) | no two walls write one face | **new** | see section 4 |

### 2.1 `B1_ZZ_F(NWE+NWI,NS)` (row 4)

- Column `N` is `B1%ZZ_F(N)`. The driver fills columns `1..N_TRACKED_SPECIES`.
- `NS = N_TOTAL_SCALARS = N_TRACKED_SPECIES + N_PASSIVE_SCALARS` (read.f90:3032). Upstream `ZZ_F` has only
  `N_TRACKED_SPECIES` entries, but L0880 indexes it with `N` up to `N_TOTAL_SCALARS` (mass.f90:65, 119-124). That is an
  out-of-bounds read upstream whenever `N_PASSIVE_SCALARS>0`. `N_PASSIVE_SCALARS` is declared with default 0 (cons.f90:469)
  and no assignment exists in Source/ at this commit, so today `NS = N_TRACKED_SPECIES`.
  Requirement: the driver must abort (or zero-fill columns above `N_TRACKED_SPECIES`, to be stated) if `NS>N_TRACKED_SPECIES`;
  the kernel keeps the upstream indexing `B1_ZZ_F(IW,N)` verbatim.
- L0882 needs no extra field for `MW_F`: it is `1/DOT_PRODUCT(MWR_Z,ZZ_GET)` (func.f90:1700-1707) evaluated in the kernel from
  row `IW` of `B1_ZZ_F` and the existing table `MWR_Z(NS)` (s5_markers.toml:31). A per-wall `B1_MW_F` precomputed by the driver would NOT be bitwise
  safe unless it used the same `DOT_PRODUCT` order; do not offer it.
- Dimension order must match the existing `B1_RHO_D_DZDN_F(NWE+NWI,NS)`: wall index first.

### 2.2 `WALL_INDEX(0:IBAR+1,0:JBAR+1,0:KBAR+1,-3:3)` (row 5)

- Preferred form: the raw upstream integer (value = wall index `IW>=1`, or 0). Then the nest text `WALL_INDEX(n)>0` and the
  `WALL_INDEX(n)==0` forms (mass.f90:907-912, in the DELTA_RHO_ZZ nest) stay verbatim and the generator needs one flatten rule.
- Accepted alternative (the offered "six flags"): the same shape with value `(WALL_INDEX(n)>0)` as 0/1. Equivalent for every use in
  my nests because `WALL_INDEX` is never negative: its only writer is init.f90:3298 (`= IW`, `IW>=1`), default 0
  (type.f90:2182); the other uses in ccib.f90 test `/=0` or `>0`.
- Index 0 of the last dimension is unused (always 0). Slots -3..-1 and 1..3 are the six neighbours; the nests use `n = IOR`.
- Range of the first three indices: all `CELL_INDEX` arguments of the guards lie in `0..IBAR+1` etc. Checked per case:
  IOR=+1 reads cell `II+1`; external walls have `II=0` (so `II+1=1`), internal walls have gas cell `IIG=II+1<=IBAR` so
  `II+1<=IBAR`. IOR=-1 reads `II-1`, internal `II>=2` (gas cell `IIG=II-1>=1`), external `II=IBAR+1`. Same for y, z.
  The masks therefore never need an index outside `0..IBAR+1` (resp. J, K).
- Refresh rule: upstream only writes it in the wall-cell initialisation (init.f90:3298), which also runs when obstructions are
  created or removed. The driver must rebuild the table whenever the wall cells are (re)initialised.
- The table is NOT redundant with `BC_IIG/BC_IOR`: the guard looks at a cell that may be owned by another wall; it has to be a cell-indexed table.

### 2.3 `EXTERIOR` mask (row 6)

- Same shape as `SOLID`. `CELL%EXTERIOR` defaults `.FALSE.` (type.f90:2177). It is set `.TRUE.` only inside `BLOCK_CELL` when
  `OBST_INDEX==0` (func.f90:5521), and the only caller with `OBST_INDEX==0` is read.f90:11697-11702, which blocks the six outer ghost
  slabs (`I=0`, `I=IBP1`, `J=0`, `J=JBP1`, `K=0`, `K=KBP1`) with `IVAL=1`. It is never reset to `.FALSE.`, so it is static after
  read: the driver can build it once per mesh.
- `EXTERIOR` stays `.TRUE.` on ghost cells whose `SOLID` is later cleared (interpolated and open external boundaries,
  init.f90:3281, 5087, 5092). The thin-obstruction test is `SOLID_BOUNDARY .AND. .NOT.SOLID .AND. .NOT.EXTERIOR` and must use the
  live `SOLID` mask.
- `SOLID` (row 7) is dynamic: init.f90:3128, 3281, 5077-5092 and `BLOCK_CELL` (func.f90:5515-5518) change it. The existing mask
  already has this contract; the requirement is only that it is refreshed after obstruction creation/removal (before the nest).
- `W_THIN` (golden entry `W_THIN:wall(NWE+NWI):logical:in`; upstream `WC%THIN`, init.f90:3315) is NOT a substitute. It is evaluated once
  at wall initialisation from `SOLID(ICG)`, `SOLID(IC)`, `EXTERIOR(IC)`, whereas mass.f90:108 and :239 test the live ghost `SOLID` and
  also require `BOUNDARY_TYPE==SOLID_BOUNDARY`. Using `W_THIN` would change results after obstruction removal.

## 3. Index ranges used by the kernels

| Quantity | Range | Reason |
|---|---|---|
| `IW` | `1..NWE+NWI` | loop bounds (mass.f90:92) |
| `BC_IOR` | `+-1,+-2,+-3` | `SELECT CASE`; any other value is a no-op in both phases (no `CASE DEFAULT`) |
| face stores `F*(IIG-1,..)`/`F*(IIG,..)` | `FX(0:IBAR,..)`, `FY(.,0:JBAR,.)`, `FZ(..,0:KBAR)` | `IIG` is the gas cell, `1<=IIG<=IBAR` |
| off-wall stores `FX(II+1)`, `FX(II-2)` | `1..IBAR-1` resp. `0..IBAR-1`, same for y,z | from the ranges in 2.2 |
| velocity reads `UU(II+1)`, `UU(II-2)` | `0..IBAR+1` | `UU` is `(0:IBAR+1,..)` |
| scalar reads `RHO_Z_P(II+1:II+2)`, `RHO_Z_P(II-2:II-1)` | `-1..IBAR+2` | `WORK_PAD` is `(-1:IBP1+1,..)`; for `IBAR=1`, `II+2=3=IBP1+1` still inside |
| `N` | `1..NS` (L0880, L0401), `0` (L0882, L0403) | `FX/FY/FZ` and `FX_ZZ/..` lower species bound is 0 only with `FLUX_LIMITER_MW_CORRECTION` (init.f90:606-608) |

The driver should check the first two rows' implied bounds once (`1<=BC_IIG<=IBAR` etc.) when it builds the tables.

## 4. Face-write conflict check (what "no two walls write the same off-wall face" must mean)

Each loop is `!$OMP DO` upstream and writes face arrays from many walls, so correctness of a one-thread-per-wall kernel depends on
the set of faces written. A purely static test "no two walls target one off-wall face" is both too strong and incomplete:

- Too strong: in a two-cell gas gap, wall A (IOR=+1, ghost `II`) and wall C (IOR=-1, ghost `II+3`) both target `FX(II+1)`. Their
  guards `UU(II+1)>0` and `UU(II+1)<0` are mutually exclusive at run time, so at most one store happens. This pairing is legal and
  must pass.
- Incomplete: wall-face stores and off-wall stores can collide too (e.g. a one-cell gap: off-wall target of wall A `FX(II+1)` equals
  the wall-face of wall B with IOR=-1 and ghost `II+2`); upstream avoids it through the `WALL_INDEX` guard, not through disjoint targets.

Required check, as a predicate over the built tables (per mesh, per axis, keyed by the face `(axis,i,j,k)`):

- Define W-writers: walls with `BOUNDARY_TYPE!=NULL` and, in the `ELSE` branch, `!=INTERPOLATED`, plus all thin-path walls
  (`SOLID_BOUNDARY` with live `SOLID=0,EXTERIOR=0` at the ghost cell). Their face is `IIG-1`/`IIG` as in 1. Store value: 0 (thin path) or
  `RHO_F*ZZ_F(N)` / `RHO_F/MW_F`.
- Define O-writers: walls with type not in {NULL, INTERPOLATED, OPEN}, IOR in `+-1..+-3`, whose neighbour guard flag is 0. Their face is
  `II+1`/`II-2` as in 1. They store only when the sign condition on the velocity holds.
- Pass when, for every face: (a) at most one W-writer, or every W-writer is thin-path (all store 0: duplicate stores are idempotent;
  this is the thin-obstruction case where both sides write the same face); (b) at most two O-writers, and if two, their IOR on that
  axis are opposite; (c) a face with a W-writer and an O-writer is reported (the guard should have removed the O-writer; if it did not,
  the serial order would matter and the kernel is not exact).
- Fail with a listing (IW pairs, face, axis) instead of silently serialising. Because the guard flags (row 5) are static between
  wall initialisations, the check can be evaluated once per rebuild of the tables, not per step.
- Not required to be checked per step: the velocity signs (data-dependent) and the species index (every species and the `0` slot
  write distinct array planes, so the conflict structure is identical for all `N`).

For (a) with mixed boundary types, two W-writers on one face where one is thin-path (stores 0) and the other is not (stores
`RHO_F*ZZ_F`) would make the result order-dependent. The check must flag that mixed case; upstream relies on it not occurring.

## 5. Additional findings that affect the nests (not wall-table items)

1. `Z_TEMP` is not fully assigned in the divergence nests. In L0401/L0403 the constructors are 3 elements wide:
   IOR=+1: `Z_TEMP(0:2)=(/RHO_Z_P(II+1),RHO_Z_P(II+1:II+2)/)` (divg.f90:1037, 1136), IOR=-1: `Z_TEMP(1:3)=...` (divg.f90:1047, 1146),
   and the same for y and z. The mass.f90 versions (L0880/L0882) are 4 elements wide with `DUMMY=0` (mass.f90:141, 151, 158, 165,
   172, 179; 273-311). `GET_SCALAR_FACE_VALUE` with `LIMITER=MP5_LIMITER` reads `U(I+IP2)` (A>0, IOR=+1: `Z_TEMP(3)`) resp.
   `U(I+IM1)` (A<0, IOR=-1: `Z_TEMP(0)`) (func.f90:1443-1448). With MP5 the divg nests therefore read an element never assigned in that
   iteration (per-thread scratch from an earlier wall, or uninitialised); verified unchanged in the current upstream worktree and in master. Other limiters never read it (func.f90:1346-1436). The
   plan for the generated kernels is to initialise the scratch to 0 and to compare bitwise only under the non-MP5 limiters for
   L0401/L0403, recording MP5 as upstream-undefined there.
2. L0403 builds `Z_TEMP` from `RHO_Z_P` (divg.f90:1136-1177) while its global-face call uses `RHO_RMW` (divg.f90:1111-1113), unlike L0882,
   which uses `RHO_RMW` in both (mass.f90:215, 273-311). This is **not** a numerical difference: `RHO_Z_P` is aimed at `WORK_PAD`
   (divg.f90:991) and `RHO_RMW` is aimed at the same array (divg.f90:1094), which holds `RHO/MW` at that point. The kernel keeps the
   upstream text and reads the `WORK_PAD` table once under either name. (An earlier reading of this as a possible defect is withdrawn.)
3. The off-wall call has the argument order `GET_SCALAR_FACE_VALUE(A=U_TEMP, U=Z_TEMP, F=F_TEMP, 1,1,1,1,1,1, IOR, LIMITER)`
   (func.f90:1330), with scratch `F_WORK(0:3,0:3,0:3)` and `U_WORK,Z_WORK(-1:3,-1:3,-1:3)` (mass.f90:27-28). The kernels need the
   callee interface from the Mesh Data stream; the wall tables above do not depend on it.

## 6. What I need back

1. Confirmation of rows 4-6 names and shapes (or the exact alternative names), and whether row 5 is raw or 0/1.
2. A statement of the driver's refresh points for rows 5 and 7 (after wall-cell initialisation and obstruction creation/removal).
3. The section 4 check, with the three outcomes (a)-(c), reachable as a table-only function so my kernels can assert it in their tests.
4. Whether `NS>N_TRACKED_SPECIES` aborts or zero-fills (2.1).
