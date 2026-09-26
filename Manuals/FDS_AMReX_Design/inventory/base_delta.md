# Base delta: ce1f659 (`(local FDS master checkout)/Source`) → FireX (`Source`)

**New reference base:**
- Path: this repository
- Branch: `AMReX` (git worktree, read-only)
- HEAD: `36975d765fcead401e14b094a04f910ac42eab8a`, "Merge pull request #16596 from cxp484/FireX", 2026-09-24

**Old base:** `(local FDS master checkout)` at ce1f659.

**Size:**
- Source/ holds **34** `.f90` files (33 + new `vtkf.f90`) and 180,840 lines. (The 35 quoted earlier counted one file too many; `find Source -name '*.f90' | wc -l` = 34.)
- ce1f659 had 33 files and 173,066 lines.

**How it was produced:**
- `tools/inventory/remap.py` does line-level `difflib` over each file pair.
- `fortran_scope.collect_globals` was run on both trees (all 100 derived types and every module-level declaration).
- `routine_index.py` was run on both trees, for routine-level body comparison ignoring indentation.
- Format: "old → new" is `file:line` on ce1f659 → `file:line` on FireX.

## 1. Derived types (MESH_TYPE and all nested types)

No type was added or removed; all 100 types are present in both trees. No existing component was removed or re-declared (same type, kind, dims and attributes).

**MESH_TYPE:** 431 → 444 components (+13). All the new ones are output bookkeeping for VTK/ParaView and unique slices:

| added component | FireX decl | notes |
|---|---|---|
| `NP`, `NC` (INTEGER) | mesh.f90:159 | Points/cells per mesh for VTK. Read at vtkf.f90:845-846. No MESH_POINTERS alias. |
| `N_SLCF_O`, `N_SLCF_VTK` | mesh.f90:323 | `N_SLCF_VTK` is set only on `MESHES(1)` (read.f90:16166) and read from `MESHES(1)` (vtkf.f90:1434, 3718). |
| `N_UNIQUE_SLCF` | mesh.f90:325 | |
| `UNIQUE_SLICE_NAMES`, `ALL_SLICE_NAMES`, `ALL_SLICE_QUANTITIES` | mesh.f90:326 | |
| `ALL_SLICE_TOPOLOGIES` | mesh.f90:327 | |
| `EMPTY_UNIQUE_SLICE`, `UNIQUE_SLICE_IS_SL3D`, `UNIQUE_SLICE_IN_DOMAIN` | mesh.f90:328 | `UNIQUE_SLICE_IN_DOMAIN` is allocated on `MESHES(1)` only (read.f90:16003). |
| `UNIQUE_SLCF_AGL` | mesh.f90:329 | |

**Other types:**
- `SLICE_TYPE` gains `SLCF_NAME` (type.f90:1603).
- `ZONE_SOLVE_TYPE` (type.f90:1791 → 1792) gains `NUNKH_LOCAL_RS`, `TOT_NNZ_H_RS`, `NUNKH_LOC_RS` and `UNKH_IND_RS` (type.f90:1816-1818). These are the "resource set" (GPU) matrix-gather counters used by the pressure path (§4).

**Line shift:** the MESH_TYPE block moves from mesh.f90:16-346 to mesh.f90:16-354. `MESHES` moves from mesh.f90:348 to 356.

**Effect on mesh_fields.csv:** 1,314 → 1,328 rows. The 14 new rows are categorised `other` under the sub-group "output / VTK and unique-slice bookkeeping (new in FireX)".

**Module-level additions** (these feed `module_globals.csv`):
- GLOBAL_CONSTANTS: +40 variables.
  - Resource-set communicators and ranks: `FDS_RANKS_PER_GPU` (cons.f90:740), `MPI_COMM_RS`, `MPI_COMM_RS_MASTERS` (cons.f90:743/746), `MY_RANK_RS`, `N_MPI_RS`, `MASTERS_RS` (cons.f90:741-750).
  - `HYPRE_DEVICE_RUN` (cons.f90:569).
  - VTK/STL/ParaView file names and switches.
- OUTPUT_CLOCKS: +22 VTK clocks and counters.
- FDS program: `MESHES_PER_PROCESS`, `N_WRITTEN`, `ERROR` (main.f90:92-95).
- New module `VTK_FDS_INTERFACE` (vtkf.f90:3-4443).

## 2. POINT_TO_MESH / MESH_POINTERS

**MESH_POINTERS:**
- Moves from mesh.f90:353 to mesh.f90:361.
- Goes from 396 to 405 declared pointers. The +9 are `N_SLCF_O` (mesh.f90:480) and the eight unique-slice arrays (mesh.f90:497-500).

**POINT_TO_MESH:**
- Moves from mesh.f90:494-896 to mesh.f90:505-916.
- Goes from 393 to 402 remaps (mesh.f90:513-914), the +9 matching the new pointers: mesh.f90:864 and 907-914.
- No other lines changed.
- `MTR`, `MSR` and `WEM` are still declared and never remapped.

**Call sites:** 200 → 214 live `CALL POINT_TO_MESH` sites, in 163 routines instead of 152.
- The +14 are in dump.f90 (+5) and vtkf.f90 (+9).
- **New:** two calls pass a literal mesh number, `CALL POINT_TO_MESH(1)` at vtkf.f90:1875 and vtkf.f90:2482. Together with the `MESHES(1)%...` fields above, mesh 1 now acts as a holder for global VTK slice metadata. This is relevant to box renumbering on regrid.

## 3. MESH_EXCHANGE codes / POST_RECEIVES

- **MESH_EXCHANGE** moves from main.f90:2975-3793 to main.f90:3117-3975.
  - The set of CODEs, the package contents and the pack/unpack logic are **unchanged**. A whitespace-insensitive diff shows only three things:
    1. a routine-local `USE MEMORY_FUNCTIONS, ONLY: PACK_PARTICLE,PACK_WALL,...` (main.f90:3119);
    2. a new local `T_NOW_SUB`;
    3. 21 timing lines accumulating into `T_USED(12)` around each pack/unpack block.
  - The PKG1 unpack block was re-indented (e.g. CODE 1/4 unpack main.f90:3542-3572 → 3704-3736).
- **POST_RECEIVES** moves from main.f90:2803-2969 to main.f90:2945-3111. No non-whitespace change.
- **CODE 5 pack** moves from main.f90:3206-3247 to 3352-3395; **CODE 5 unpack** from main.f90:3576-3604 to 3740-3770. The extra lines are the timing instrumentation.
- **Renamed routines** (a typo fix): `ALLOCATE_RADIAITON_RECV_PKG`/`_SEND_PKG` (main.f90:3796/3837) became `ALLOCATE_RADIATION_RECV_PKG`/`_SEND_PKG` (main.f90:3980/4024).
- **New collective-exchange helpers** for VTK: `EXCHANGE_NOBST_INFO`, `EXCHANGE_NSLICE_INFO` and `EXCHANGE_NPATCH_INFO` (main.f90:4721-4965), plus program-level functions `NP(NM)` and `NC(NM)` (main.f90:4778-4794). These function names shadow nothing in MESH_TYPE, which is only reached as `MESHES(NM)%NP`.

## 4. Pressure path

- **`PRESSURE_ITERATION_SCHEME`** moves from main.f90:1470-1614 to main.f90:1601-1745 and is **byte-identical**. Old → new call sites:

  | call | ce1f659 | FireX |
  |---|---|---|
  | CODE 5 (CC_IBM only) | main.f90:1486 | main.f90:1617 |
  | CODE 5 after BAROCLINIC_CORRECTION | 1505 | 1636 |
  | CODE 5 after GLMAT_SOLVER | 1532 | 1663 |
  | `COPY_H_OMESH_TO_MESH` | 1533 | 1664 |
  | CODE 5 before COMPUTE_VELOCITY_ERROR | 1556 | 1687 |
  | MPI_ALLGATHERV of errors | 1571 | 1702 |

- **GLMAT_SOLVER** (pres.f90:3241-3551 → 3299-3658): the global-matrix solve now runs on **resource-set (RS) masters only**.
  - Every rank in an RS `MPI_GATHERV`s its RHS `F_H` to the RS master (pres.f90:3414-3419).
  - Only `MY_RANK_RS==0` solves, inside `MASTER_IF` (pres.f90:3422-3482). The PARDISO/cluster communicator changes from `MPI_COMM_WORLD` to `MPI_COMM_RS_MASTERS` (pres.f90:3436, 3440).
  - HYPRE vectors migrate to and from the device (pres.f90:3455-3472).
  - The solution is sent back to the RS ranks with `MPI_SCATTERV` (pres.f90:3486-3491).
- **GET_H_MATRIX_LUDCMP** (pres.f90:4254-4592 → 4372-4876):
  - the matrix graph is gathered to RS masters (MPI_GATHERV, pres.f90:4511, 4570, 4572; `MASTER_IF` pres.f90:4613);
  - the HYPRE sub-communicator is split from `MPI_COMM_RS_MASTERS` (pres.f90:4733).
- **ULMAT** (local-matrix solver):
  - ULMAT_SOLVER_SETUP sets the HYPRE execution policy and memory location for device or host (pres.f90:1176-1187).
  - The HYPRE vectors migrate device↔host in ULMAT_SOLVE_ZONE (pres.f90:1753-1770) and in ULMAT_H_MATRIX_SOLVER_SETUP (pres.f90:3011-3044).
  - PRHS and the BXS..BZF re-allocation move from pres.f90:1260-1283 to 1280-1303; the shapes are unchanged.
- **RS definition:**
  - `DEFINE_RS_COMM_INFO` (main.f90:5096-5146) is called at main.f90:153.
  - It reads the environment variable `FDS_RANKS_PER_GPU` (main.f90:5108) and splits `MPI_COMM_WORLD` by `MY_RANK/FDS_RANKS_PER_GPU` (main.f90:5123-5124).
  - `HYPRE_DEVICE_RUN` defaults to .TRUE. (cons.f90:569) and is a PRES namelist input (read.f90:10069).
- **FFT path unchanged:** PRESSURE_SOLVER_FFT and PRESSURE_SOLVER_COMPUTE_RHS show no body change. TUNNEL_POISSON_SOLVER has 4 changed lines.

## 5. Other notable routine-level changes

- 72 routines were added: VTK/HDF5 output in dump.f90 and vtkf.f90 (59 in vtkf.f90), `FIND_WALL_INDEX` (func.f90:5449), and `SORT_COLUMNS_BY_FIRST_ROW` (init.f90:2155).
- 3 were removed or renamed:
  - `MASS_FINITE_DIFFERENCES_NEW` became `MASS_FINITE_DIFFERENCES` (mass.f90:20 in both; called at main.f90:670/851 → 767/948);
  - the two `ALLOCATE_RADIAITON_*` routines (§3).
- 76 routines changed, ignoring indentation. The largest changes are:
  - WRITE_SMOKEVIEW_FILE, READ_SLCF;
  - GET_H_MATRIX_LUDCMP;
  - COMPUTE_RADIATION/RADIATION_FVM (radi.f90:3740-5002 → 3744-5144);
  - READ_DUMP;
  - INITIALIZE_HT3D_WALL_CELLS;
  - DUMP_MESH_OUTPUTS;
  - GLMAT_SOLVER.

  The full list can be reproduced with `routine_index.py` on both trees.

## 6. Line-number remap for citations already given to the team

Status key:
- `same`: same line number and identical text.
- `moved`: identical text at a new line number.
- Ranges map both endpoints.
- Nothing cited was changed or removed.

| ce1f659 citation | FireX | status | what |
|---|---|---|---|
| mesh.f90:16-346 | mesh.f90:16-354 | start same, end moved; +13 components inside | MESH_TYPE |
| mesh.f90:348 | mesh.f90:356 | moved | `MESHES` declaration |
| mesh.f90:353 | mesh.f90:361 | moved | MODULE MESH_POINTERS |
| mesh.f90:494-896 | mesh.f90:505-916 | moved; +9 remaps inside | POINT_TO_MESH |
| mesh.f90:502-894 | mesh.f90:513-906 (last remap now mesh.f90:914) | moved | remap lines |
| main.f90:637 | main.f90:734 | moved | the only main.f90 POINT_TO_MESH call |
| main.f90:1939 | main.f90:2077 | moved | `ALLOCATE(MESHES(NM)%OMESH(NMESHES))` |
| main.f90:3206-3247 | main.f90:3352-3395 | moved; +timing lines | CODE 5 pack |
| main.f90:3576-3604 | main.f90:3740-3770 | moved; +timing lines | CODE 5 unpack |
| main.f90:1486 / 1505 / 1532 / 1533 / 1556 / 1571 | main.f90:1617 / 1636 / 1663 / 1664 / 1687 / 1702 | moved | pressure-iteration calls |
| main.f90:1470-1614 | main.f90:1601-1745 | moved (identical body) | PRESSURE_ITERATION_SCHEME |
| main.f90:2975-3793 | main.f90:3117-3975 | moved | MESH_EXCHANGE |
| main.f90:2803-2969 | main.f90:2945-3111 | moved | POST_RECEIVES |
| main.f90:812 / 1009 | main.f90:909 / 1105 | moved | CODE 3 (after predictor) / CODE 6 (after corrector) |
| init.f90:522-529 | init.f90:522-529 | same | 2-ghost-layer arrays |
| init.f90:533-535 | init.f90:533-535 | same | U/V/W allocation |
| init.f90:543-545 | init.f90:543-545 | same | FVX/FVY/FVZ |
| init.f90:2379 | init.f90:2445 | moved | PRHS |
| init.f90:2382-2387 | init.f90:2448-2453 | moved | BXS..BZF |
| init.f90:449 | init.f90:449 | same | P_0 |
| pres.f90:1260-1283 | pres.f90:1280-1303 | moved | ULMAT re-allocation of BXS..BZF/PRHS |
| pres.f90:3936 | pres.f90:4054 | moved | SUBROUTINE COPY_H_OMESH_TO_MESH |
| wall.f90:282-388 | wall.f90:282-388 | same | ASSIGN_GHOST_VALUE |
| velo.f90:516-546 | velo.f90:516-546 | same | OMESH MU/KRES/D/DS averaging |
| velo.f90:1394 | velo.f90:1394 | same | NO_FLUX reads OM_HP |
| type.f90:462-482 | type.f90:462-482 | same | EXTERNAL_WALL_TYPE |
| type.f90:478-479 | type.f90:478-479 | same | `FVN`, `FVNS` (EXTERNAL_WALL_TYPE) |
| ccib.f90:370 | ccib.f90:370 | same | "Assumes POINT_TO_MESH(NM) has been called" |

**Every other citation:** the checkpoint-1 README citations were re-mapped automatically with `remap.py`. All resolved to identical source text (status same or moved). README.md now cites FireX line numbers throughout.
