# Upstream issue candidates (firemodels/fds): drafts, not filed

Prepared 2026-09-25 from the FDS mesh-data inventory.

**Bases.**
- FireX: `Source`, branch `AMReX`, `36975d765fcead401e14b094a04f910ac42eab8a`.
- Upstream master: `(local FDS master checkout)/Source` at `ce1f659cd419aa4f7bb6f60c77e3b009ebbfe824` (2026-08-31).
- Line numbers are given for both.

**Status.**

| # | title | status |
|---|---|---|
| 1 | Shape-OBST area adjustment depends on mesh decomposition and MPI process count | **confirmed by runs** on both FireX and master binaries |
| 2 | Thin-OBST collapse in y compares against the x origin (XS instead of YS) | found by reading, **not run-confirmed** |
| 3 | Thin-OBST collapse in z uses the x indices (I1/I2 instead of K1/K2) | found by reading, **not run-confirmed** |
| 4 | Thin-OBST collapse mixes metres and cell indices (x, y, z) | found by reading, **not run-confirmed** |

The brief's "candidates 2 and 3" are the A-20 slips. The unit mix is listed separately as #4 so that each issue has a single fix. File #2–#4 together if preferred.

---

## Candidate 1: `SHAPE` OBST area adjustment depends on the mesh decomposition and MPI process count (MULTIPLIER%FDS_AREA never reset or reduced)

### Summary

A `SHAPE='CYLINDER'` OBST built from a `MULT` voxel array straddles a mesh boundary. The total mass flux (and hence HRR) from its surfaces then changes when the same domain is:
- split into two meshes, or
- run on one versus two MPI processes.

In the reproducer below, the burning side area is 1.5× the true value with 2 meshes on 1 process and 2× with 2 meshes on 2 processes. The top area is 0.75× with 2 meshes on 1 process. The result also depends on where the split lies.

### Affected code (identical in FireX and master)

| what | FireX 36975d7 | master ce1f659 |
|---|---|---|
| `ADJUST_OBST_SHAPE_AREA`: accumulate `MR%FDS_AREA` (pass 1), then `B1%AREA_ADJUST = OB%SHAPE_AREA/MR%FDS_AREA` (pass 2) | init.f90:339-432 (accumulate 355-390, adjust 394-430) | init.f90:339-432 (byte-identical routine) |
| called once per mesh, inside the per-mesh `INITIALIZE_MESH_VARIABLES_1` | init.f90:917, main.f90:246-248 | init.f90:917 |
| `FDS_AREA(6)=0` default; only initialisation, never reset | type.f90:1856 | type.f90:1852 |
| global flag `OBST_SHAPE_AREA_ADJUST`: multi-mesh allowed only for CYLINDER; overwritten by each shaped OBST line | read.f90:10947-10954; cons.f90:267 | read.f90:10814-10821 |
| cylinder top/bottom `SHAPE_AREA(1)` = circle ∩ **current mesh** | read.f90:10962, 10967 | read.f90:10829 |
| cylinder side `SHAPE_AREA(2)` = 2πRH (**whole** cylinder) | read.f90:10966 | same block |
| `OB%MULT_INDEX = MULT_INDEX`: all voxels of one MULT share one `MULTIPLIER` entry | read.f90:11251 | read.f90:11118 |

### Mechanism

- `MULTIPLIER(:)` is a global array on each process, and `FDS_AREA` is summed into it mesh by mesh (init.f90:366-387).
- The sum is never zeroed between meshes and never summed across MPI processes.
- It is used right away, in the same call, as the denominator of `AREA_ADJUST` (init.f90:407-426).

The consequences are as follows.
- **Several meshes on one process.** Mesh k divides by the faces of meshes 1..k. The first mesh gets too large an adjustment; later meshes are diluted by earlier meshes' faces.
- **One mesh per process.** Each mesh divides by only its own faces.
- **What the target needs:**
  - For the side (whole-shape `SHAPE_AREA(2)`), the target needs the global face sum.
  - For the top and bottom (per-mesh `SHAPE_AREA(1)`/`(3)`), it needs the per-mesh face sum.
  - The code gives neither consistently.
- **Shared MULT.** The accumulator is keyed by MULT rather than by OBST line, so two shaped OBSTs that share one `MULT_ID` dilute each other even on a single mesh (variant e).

### Minimal input

The base case is variant a. Variants b–e change only the `MESH` lines, noted after the listing.

```
&HEAD CHID='a_1mesh', TITLE='MULTIPLIER FDS_AREA shape-area adjust decomposition test' /
&MESH IJK=20,20,15, XB=0.0,0.4,0.0,0.4,0.0,0.3 /
&TIME T_END=0.1 /
&DUMP DT_HRR=0.02, DT_DEVC=0.02 /
&SPEC ID='METHANE' /
&SPEC ID='PROPANE' /
&RAMP ID='one', T=0.0, F=1.0 /
&RAMP ID='one', T=100.0, F=1.0 /
&SURF ID='top',  MASS_FLUX(1)=0.01, SPEC_ID(1)='METHANE', RAMP_MF(1)='one', COLOR='RED' /
&SURF ID='side', MASS_FLUX(1)=0.01, SPEC_ID(1)='PROPANE', RAMP_MF(1)='one', COLOR='BLUE' /
&SURF ID='wall', COLOR='GRAY' /
&MULT ID='voxels', DX=0.02, DY=0.02, I_UPPER=9, J_UPPER=9 /
&OBST XB=0.10,0.12,0.10,0.12,0.00,0.10, MULT_ID='voxels', SURF_IDS='top','side','wall',
      SHAPE='CYLINDER', RADIUS=0.1, HEIGHT=0.1, XYZ=0.2,0.2,0.0 /
&VENT MB='XMIN', SURF_ID='OPEN' /
&VENT MB='XMAX', SURF_ID='OPEN' /
&VENT MB='YMIN', SURF_ID='OPEN' /
&VENT MB='YMAX', SURF_ID='OPEN' /
&VENT MB='ZMAX', SURF_ID='OPEN' /
&TAIL /
```

- **(b) 2 meshes, 1 process:** replace the MESH line with `&MESH IJK=10,20,15, XB=0.0,0.2,0.0,0.4,0.0,0.3 /` and `&MESH IJK=10,20,15, XB=0.2,0.4,0.0,0.4,0.0,0.3 /`, and run `mpiexec -n 1`.
- **(c) 2 meshes, 2 processes:** as (b), with `MPI_PROCESS=0` / `MPI_PROCESS=1`, run `mpiexec -n 2`.
- **(d) off-centre split, 1 process:** `IJK=7,...,XB=0.0,0.14,...` and `IJK=13,...,XB=0.14,0.4,...`.
- **(e) 1 mesh, two cylinders sharing one MULT:** `&MESH IJK=40,20,15, XB=0.0,0.8,0.0,0.4,0.0,0.3 /`, plus a second OBST with `XB=0.50,0.52,0.10,0.12,0.00,0.10`, `XYZ=0.6,0.2,0.0`, same MULT and SURFs.

The top faces emit METHANE and the sides emit PROPANE, so `MLR_METHANE` and `MLR_PROPANE` in `CHID_hrr.csv` measure the effective top and side areas separately (× 0.01 kg/m²/s). The bottom sits on the floor and is not exposed. Each run takes about 1–2 s.

### Expected vs actual

**Expected**, independent of decomposition:
- MLR_METHANE = 0.01·πR² = 3.1416e-4 kg/s
- MLR_PROPANE = 0.01·2πRH = 6.2832e-4 kg/s
- Variant e: twice both values.

Values are from `CHID_hrr.csv` at t = 0.1 s. The FireX and master binaries give the same numbers to all printed digits.

| variant | meshes / processes | MLR_METHANE (top) | ratio to exact | MLR_PROPANE (side) | ratio to exact |
|---|---|---|---|---|---|
| a | 1 / 1 | 3.1415927E-04 | 1.000 | 6.2831853E-04 | 1.000 |
| b | 2 / 1 (split x=0.2) | 2.3561945E-04 | **0.750** | 9.4247780E-04 | **1.500** |
| c | 2 / 2 (split x=0.2) | 3.1415927E-04 | 1.000 | 1.2566371E-03 | **2.000** |
| d | 2 / 1 (split x=0.14) | 2.7374480E-04 | **0.871** | 1.0681415E-03 | **1.700** |
| e | 1 / 1, 2 cylinders, shared MULT (exact 6.2832e-4 / 1.2566e-3) | 3.1415927E-04 | **0.500** | 6.2831853E-04 | **0.500** |

**The mechanism predicts b–d exactly.**
- A face-count model reproduces all six b–d values to 8 digits. It uses the voxel faces per mesh (top 0.016/0.016 m² and side 0.04/0.04 m² for the x=0.2 split; top 0.0048/0.0272 and side 0.024/0.056 for x=0.14), `FDS_AREA` accumulated over meshes 1..k, and the per-mesh circle ∩ mesh `SHAPE_AREA(1)`.
- For example, (b) top = πR²/2 + (πR²/2)·½ = 0.75πR², and side = 2πRH + 2πRH·½ = 3πRH. (c) side = 2πRH per mesh.

### Impact

- **Who is affected:** any shaped OBST whose MULT voxels span more than one mesh. With multiple meshes that means CYLINDER only (read.f90:10950); with one mesh, any shape that shares a MULT.
- **What changes:** the burning rate, HRR, and every other flux scaled by `B1%AREA_ADJUST` (e.g. wall.f90:1316; dump.f90:10239-10412) depend on the mesh layout and the rank count. Combined with the global `OBST_SHAPE_AREA_ADJUST` flag, even the choice of which shaped OBST line comes last in the input matters.
- **Validation cases (*uncertain*).** The `Validation/FM_Burner` cases (e.g. `FM_15cm_Burner_CH4_5mm.fds`) use a MULT-voxel cylinder over a 208-mesh layout that splits it at x=0, y=0 and in z.
  - Only the top burns (`SURF_IDS='burner','wall','wall'`). With one mesh per process the top is per-mesh consistent, so those cases are probably unaffected as usually run.
  - They would change if meshes were grouped onto fewer processes.

### Suggested fix

1. **Key the accumulator by shape, not by MULT.** Store it on the OBST line: give each shaped OBST line a shape index, and give `OBSTRUCTION` a `SHAPE_INDEX` alongside `MULT_INDEX`. This also fixes variant e.
2. **Split the routine into two phases.**
   - Phase 1: loop over all meshes of the process, zeroing the accumulator once before the mesh loop.
   - Then `MPI_ALLREDUCE(MPI_SUM)` of the shape-area accumulators.
   - Phase 2: loop again to set `AREA_ADJUST`.
   - Move the call out of the per-mesh `INITIALIZE_MESH_VARIABLES_1`, or split it into two routines called from main.f90 on either side of the reduction.
3. **Use whole-shape target areas throughout.** Replace the per-mesh `CIRCLE_CELL_INTERSECTION_AREA(...,M%XS,M%XF,M%YS,M%YF)` at read.f90:10962 with the whole-shape value (πR², or the intersection with the domain/XB). The global denominator of fix 2 then gives exact totals for top, side and bottom on any decomposition.
4. **Make the enable flag per OBST.** Replace the global `OBST_SHAPE_AREA_ADJUST` with a per-OBST flag, e.g. `OB%SHAPE_AREA_ADJUST`, so it is no longer decided by the last shaped OBST line read.

Once this is done, the NMESHES==1 restriction for SPHERE/CONE/BOX at read.f90:10950 could probably be dropped (*uncertain*).

**Note for FireX/AMR.** The global shape-area sum must be accumulated exactly (fixed-point) to be bit-identical across box layouts; see `global_reductions.csv`.

### Related observation (reading only, *uncertain*)

The BOX branch of pass 2 pairs faces and target areas inconsistently (init.f90:421-426 in both bases).
- `IOR=-1` uses `SHAPE_AREA(1)` (=LENGTH·WIDTH), but `IOR=+1` uses `SHAPE_AREA(2)`.
- `IOR=+2` uses `SHAPE_AREA(1)`, and `IOR=-2` uses `SHAPE_AREA(3)`.
- So opposite faces of the box get different target areas.

This path is active only on single-mesh runs with a default ORIENTATION. It was not run-tested; report it in the same issue as a question.

---

## Candidate 2: thin-OBST collapse in y uses the x origin `XS` (READ_OBST)

**Status: found by reading, not run-confirmed.**

| | FireX 36975d7 | master ce1f659 |
|---|---|---|
| routine | read.f90 `READ_OBST` (10742-11726), block 11185-11198 | same block, 11052-11065 |
| guard (correct) | read.f90:11192 `GINV(XB4-YS,2,NM)-GINV(XB3-YS,2,NM)<0.25_EB/RDETA` | read.f90:11059 |
| **slip** | read.f90:11193 `IF(GINV(XB3-XS,2,NM)-REAL(OB%J1,EB) < REAL(OB%J2,EB) - GINV(XB4-YS,2,NM))` | read.f90:11060 |

**What happens.**
- A thin OBST in y (thinner than ¼ cell, not THICKENed) is collapsed to one face, J1 or J2, chosen by which side it is closer to.
- The left-hand distance uses `XB3-XS` (the y coordinate minus the mesh **x** origin) where `XB3-YS` is intended.
- So the chosen face depends on the mesh's x origin. The same thin OBST can collapse to different y faces in meshes with different XS.

**Expected fix:** `GINV(XB3-YS,2,NM)`, together with the unit fix of #4.

**Suggested test (not run):** a thin plate at y = j·Δy + 0.1Δy, in two single-mesh runs whose domains differ only in XS (e.g. XB x-range 0..1 vs 5..6). The OBST face y-location in the .smv/.out output should agree.

---

## Candidate 3: thin-OBST collapse in z uses the x indices `I1/I2` (READ_OBST)

**Status: found by reading, not run-confirmed.**

| | FireX 36975d7 | master ce1f659 |
|---|---|---|
| guard (correct) | read.f90:11209 (`...<0.25_EB/RDZETA .AND. OB%K1 /= OB%K2`) | read.f90:11076 |
| **slip** | read.f90:11210 `IF(GINV(XB5-ZS,3,NM)-REAL(OB%I1,EB) < REAL(OB%I2,EB) - GINV(XB6-ZS,3,NM))` | read.f90:11077 |
| assignments (correct, K) | read.f90:11211-11213 `OB%K2=OB%K1` / `OB%K1=OB%K2` | read.f90:11078-11080 |

**What happens.**
- The z collapse compares the z position against the **x** indices OB%I1/OB%I2 instead of OB%K1/OB%K2.
- The choice between the bottom and top face therefore depends on the OBST's x-position.

**Expected fix:** `REAL(OB%K1,EB)` / `REAL(OB%K2,EB)`, together with the unit fix of #4.

**Suggested test (not run):** a horizontal thin plate at z = k·Δz + 0.1Δz, placed at two x positions far apart. The collapsed z face should be identical.

---

## Candidate 4: thin-OBST collapse mixes metres and cell indices

**Status: found by reading, not run-confirmed.**

| | FireX 36975d7 | master ce1f659 |
|---|---|---|
| x | read.f90:11176 `GINV(XB1-XS,1,NM)-REAL(OB%I1,EB) < REAL(OB%I2,EB) - GINV(XB2-XS,1,NM)` | read.f90:11043 |
| y | read.f90:11193 | read.f90:11060 |
| z | read.f90:11210 | read.f90:11077 |
| index definition for comparison | read.f90:11157-11158 `OB%I1 = NINT(GINV(XB1-XS,1,NM)*RDXI)` | same block |
| GINV returns a length (m) | func.f90:6722-6750 | func.f90 |

**What happens.**
- `GINV(...)` returns a position in metres; cell indices are `GINV(...)*RDXI` (read.f90:11157).
- The collapse test subtracts an integer index from a length in metres. Unless Δx = 1 m, it does not measure "distance to face I1 vs I2", so the nearer face is not reliably chosen.
- The ¼-cell guard on the line before (read.f90:11175, `<0.25_EB/RDXI`) is dimensionally correct.

**Expected fix:** multiply by `RDXI`/`RDETA`/`RDZETA`, e.g. `GINV(XB1-XS,1,NM)*RDXI-REAL(OB%I1,EB) < REAL(OB%I2,EB)-GINV(XB2-XS,1,NM)*RDXI`. Do the same for y (with `YS`, #2) and for z (with `K1/K2`, #3).

**Note for FireX.** Level-0 bit-parity with master requires reproducing #2–#4 as is. Fix them only behind a flag.

---

### Reproduction artefacts (candidate 1)

- **Directory:** `(local project tools directory)/inventory/runs/fds_area/`. It holds `a_1mesh.fds`, `b_2mesh_1proc.fds`, `c_2mesh_2proc.fds`, `d_2mesh_1proc_split014.fds` and `e_1mesh_2cyl_sharedmult.fds`, with their `*_hrr.csv`, `*.out` and `*.log`. The same runs on the master binary are in `master/`.
- **Binaries:**
  - FireX: `(local build directory)/firex-36975d7/ompi_gnu_rel/fds`, revision FDS-6.11.1-1244-g36975d765f-AMReX, GCC 14.2.0, Open MPI 5.0.7.
  - master: `(local FDS master checkout)/Build/ompi_gnu_linux/fds_ompi_gnu_linux`, revision date matching ce1f659.
- **Run settings:** `OMP_NUM_THREADS=1`, `mpirun -np 1|2`; each run takes 1–2 s wall time.
