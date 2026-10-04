# Patch 0010 (DRAFTED - needs Build Chief validation): `main.f90`, AMR level-0 abort guard for coarse/fine faces (D-065 Q1, D-076)

**Status: drafted (extended with the `NIC_CHECK` hook, D-077), not built.** Cut against `Source/main.f90` of the FDS-AMReX branch tree (blob `2e1327bbe8`, identical at tips 991a759f79 and 0353957a93). `patch -p1 --dry-run` is clean on that file (see "Apply check"). Nobody has compiled it yet; the GNU Build Chief and the Intel Build Chief validate it (section "Validation request"). It continues the driver patch series (`Source/driver/patches/0001`-`0009`); the file is kept here because the source trees are read-only for the author. Apply from the top of the source tree with `patch -p1` or `git apply`.

## What it does
Four additions in `Source/main.f90`, all inside `#ifdef WITH_AMREX` (no existing statement is edited):
1. `INITIALIZE_MESH_EXCHANGE_1`, declaration part (patched main.f90:2151-2159): `USE COMP_FUNCTIONS, ONLY: SHUTDOWN` and the locals `AMR_N_BAD`, `AMR_IW1`, `AMR_MESSAGE`.
2. `INITIALIZE_MESH_EXCHANGE_1`, end of the routine, after the `MPI_PARTICLE_EXCHANGE` block and before `END SUBROUTINE` (patched main.f90:2397-2424; inserted before original line 2374, the `END SUBROUTINE`): a loop over `M%EXTERNAL_WALL(1:M%N_EXTERNAL_WALL_CELLS)`. It skips walls with `EWC%NOM<1` (NIC is set only for walls that abut another mesh, main.f90:2184) and counts the walls with `EWC%NIC>1`. If any, it writes the message below and calls `SHUTDOWN(...,PROCESS_0_ONLY=.FALSE.)` once for the mesh (first offending wall is named, the count is given).
3. The call loop in the main set-up sequence (patched main.f90:333-344; the loop itself is original 326-328): `CALL STOP_CHECK(1)` right after the `DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX / CALL INITIALIZE_MESH_EXCHANGE_1(NM)` loop, then the `NIC_CHECK` reduction and print (item 4).
4. The `NIC_CHECK` hook (D-077 (2), V&V test P-3). Declaration `INTEGER :: AMR_NIC_CHECKED=0,AMR_NIC_TOTAL=0` in the declaration part of the main routine (patched main.f90:115-117; the routine is `FDS_SETUP` under `WITH_AMREX`, which has `SAVE`, so the host variable is visible in the internal routine). In the guard loop of `INITIALIZE_MESH_EXCHANGE_1`, `AMR_NIC_CHECKED` is increased by one for every `EXTERNAL_WALL` entry with `NOM>0` (the same set the guard tests, right after the `NOM<1` skip, patched main.f90:2409). After the `STOP_CHECK(1)` the call loop does `AMR_NIC_TOTAL = AMR_NIC_CHECKED`, `MPI_ALLREDUCE(MPI_IN_PLACE,...,MPI_SUM)` over `MPI_COMM_WORLD` when `N_MPI_PROCESSES>1` (the same call form as `STOP_CHECK`, original main.f90:1976, executed once by every rank), and rank 0 writes the line below. The line is reached only when no rank aborted, because `STOP_CHECK(1)` ends all ranks first.

## Expected behaviour
| Run | Result |
|---|---|
| `WITH_AMREX` undefined | no change at all (patch compiled out) |
| `WITH_AMREX`, every level-0 face NIC=1 | exactly one line on stderr from rank 0, `NIC_CHECK level0: <n> walls checked`, with `<n>` = number of `EXTERNAL_WALL` entries with `NOM>0` summed over all level-0 meshes and ranks; then the driver's `FDS-AMReX level 0: ... box(es)` line on stdout |
| `WITH_AMREX`, any NIC>1 level-0 wall | `ERROR(9001)` line(s), `ERROR: FDS was improperly set-up - FDS stopped`, no `NIC_CHECK` line, no `level 0:` line, exit status 0 |

`NIC_CHECK level0: <n> walls checked` has no leading blank and no rank suffix. It goes to stderr (`LU_ERR`, the unit of every other FDS console message, for example `Starting FDS ...` at main.f90:180 and the stop message at main.f90:2047), not to stdout: unit 6 (`LU_OUTPUT`, cons.f90:602) is reconnected to `CHID.out` by `INITIALIZE_DIAGNOSTIC_FILE` (dump.f90:3068-3070, called at main.f90:612), so a Fortran write to unit 6 depends on timing and on buffering, while stderr is unbuffered. Tests read it from `stderr.txt`. The driver's own lines (`FDS-AMReX level 0: ...`, FdsAmr.cpp:87) are C++ and go to stdout.

The count includes periodic faces (a periodic domain face finds its partner mesh at the opposite end, init.f90:3151-3160, so `NOM>0`) and an internal mesh face seen from both sides (once per mesh). Domain-edge walls with `NOM=0` are not counted.

Hand counts (by reading the inputs and `read.f90:700`, `N_EXTERNAL_WALL_CELLS = 2*(IBAR*JBAR+IBAR*KBAR+JBAR*KBAR)`):
- `shunn3_4mesh_32` (four meshes `16,1,16`, 2 x 2 in x and z; `PBX=-1/1` and `PBZ=-1/1` are `PERIODIC`, no vent on y): per mesh the x faces give 2 x (JBAR x KBAR = 16) = 32 and the z faces 2 x (IBAR x JBAR = 16) = 32, all with `NOM>0` (neighbour or periodic partner); the y faces (2 x 256) have `NOM=0`. 4 x 64 = **256**.
- `dec4_np4.fds` (16 meshes `8,1,8`, 4 x 4 in x and z, `PERIODIC` vents on x and z, nothing on y): per mesh 2 x 8 + 2 x 8 = 32 with `NOM>0`. 16 x 32 = **512**.
Both agree with the expected values 256 and 512 (V&V P-3). Not run; the Build Chiefs confirm.

## Why
D-065 Q1 retires the `EWC%NIC>1` branches (species diffusive flux match, `wall.f90:897`; coarse flux overwrite, `divg.f90:219`) in the AMR route. D-072 (f) and D-076 make that a design requirement: the AMR route never has an interpolated boundary between level-0 meshes of different resolution (a coarse/fine jump is a refinement level, handled by ghost fill and the interface flux overwrite). Reading the code cannot prove the whole route, so the requirement is enforced at run time: any level-0 mesh with a NIC>1 external wall aborts the set-up with a clear message, before any exchange or time step. Without the guard such an input would reach the retired branches silently. Until the input pre-pass (converter) is wired, inputs with finer meshes stop here instead of at `FdsAmr.cpp:30` (which checks the cell sizes only after `FDS_SETUP` returns).

## Message
Written to the error unit by every rank that has an offending mesh (`SHUTDOWN` with `PROCESS_0_ONLY=.FALSE.` prefixes a blank line and appends ` (MPI Process:<rank>, CHID: <chid>)`):

`ERROR(9001): AMR mode needs equal-resolution level-0 meshes. MESH <NM>: external wall cell <IW> (IOR=<IOR>) abuts MESH <NOM> with NIC=<NIC> (<count> such wall cells in this mesh). Finer meshes are converted into AMR levels by the input pre-pass; run this input without AMR.`

It names the mesh (`NM`), the wall (`IW`, `IOR`, other mesh), NIC, and states the requirement and the conversion route. Code 9001 is unused in the tree (checked: no `ERROR(90xx)`).

## Multi-rank behaviour (cannot hang)
- `INITIALIZE_MESH_EXCHANGE_1` contains no MPI call (checked by reading original main.f90:2132-2374; the guard adds none). The only collectives the patch adds are the unconditional `STOP_CHECK(1)` and, after it, the `MPI_ALLREDUCE` of the `NIC_CHECK` count (reached only if `STOP_CHECK` returned, which is the same on every rank); every rank reaches each exactly once, whatever number of meshes it holds, so no rank can wait for a message from a rank that has stopped.
- `SHUTDOWN` (func.f90:120-139) only prints and sets `STOP_STATUS = SETUP_STOP`; it does not stop. `EXTERNAL_WALL` exists on the owning rank only, so the offending rank is not necessarily rank 0.
- The new `STOP_CHECK(1)` (original main.f90:1969-1993) is executed unconditionally by every rank, once, after the loop. With `N_MPI_PROCESSES>1` it does `MPI_ALLREDUCE(MPI_MAX)` of `STOP_STATUS`, so every rank sees `SETUP_STOP`; then `END_FDS` runs on all ranks (rank 0 prints `ERROR: FDS was improperly set-up - FDS stopped`, then `MPI_FINALIZE` and `STOP`). When nothing is wrong `STOP_STATUS` is 0 on all ranks and `STOP_CHECK` returns at once (two extra `MPI_ALLREDUCE` at set-up, with more than one rank, only in the `WITH_AMREX` build: `STOP_CHECK` and the `NIC_CHECK` sum). This is the pattern of the existing set-up-only stop at original main.f90:303-308.
- Exit status: `END_FDS` ends with a plain `STOP` (original main.f90:2073), so the process exit status is 0, as for every other FDS set-up error (for example `ERROR(431)`, init.f90:3185-3187). A test must key on the message text and on the absence of the driver line `level 0: <n> box(es)` and of the `NIC_CHECK` line, not on the exit status. In `FDS_SETUP` the `STOP` happens before the driver finalizes AMReX, the same as the existing set-up-only stop; this is accepted for an abort path.

## Where it is a no-op
- `WITH_AMREX` undefined (`USE_AMREX=OFF`): all four hunks are removed by the preprocessor. Checked: `gfortran -cpp -E -P` of main.f90 with the macro undefined, before and after the patch, blank lines dropped: identical (4055 lines, same method after the `NIC_CHECK` hook was added). With `-DWITH_AMREX` the preprocessed file (`-E -P`, blank lines kept) grows by 46 lines (5169 to 5215; 36 lines with blank lines dropped); 35 of the 46 were the guard alone, 11 are the `NIC_CHECK` hook.
- `WITH_AMREX` defined and no NIC>1 wall (all equal-resolution meshes, one mesh, periodic single mesh, any input with only NIC=1 faces): no `ERROR(9001)` message, one `NIC_CHECK` line, `STOP_CHECK` returns at once, results unchanged. Fine-level boxes are not in `MESHES` and have `N_EXTERNAL_WALL_CELLS=0` (`fds_fine_level.f90:162-165`), and the routine runs only for `MESHES(NM)`, so it checks level-0 meshes only.
- The guard is independent of the species count. The NIC>1 branches need 2 or more tracked species to matter, but the guard refuses every NIC>1 level-0 face (D-076: any), which is stricter and also covers the velocity and pressure-side interpolation.

## Risk
Low. An input with a coarse/fine face that ran before in the AMR build (it could only get to `FdsAmr.cpp:30` and abort there) now aborts earlier with this message. `STOP_CHECK` is placed before `MPI_INITIALIZATION_CHORES(4)`; `END_FDS` frees `N_REQ1`... requests, which are 0 at that point, as at the set-up stop at original main.f90:303-308 (to be confirmed by the negative run on 2 and 6 ranks).

## How to revert
Delete the four hunks (`patch -R -p1`), or build without `-DWITH_AMREX`. No other file is touched.

## Controls
Files: `inputs/0010-neg-control-converter-pass.fds` (negative, runs through stock `fds_amr`), `inputs/0010-neg-control-race_test_1.fds` (negative, needs the bypass build); positive controls below.

**Negative control: `inputs/0010-neg-control-race_test_1.fds`** (copy of `Verification/Thread_Check/race_test_1.fds`; only `CHID` and a `TITLE` differ). The stock `fds_amr` refuses this input in the pre-pass before FDS set-up (ratio 5); it exercises the guard only on a scratch build with the pre-pass bypassed, see "Converter-passing negative control" (fallback) below. For the stock build use `inputs/0010-neg-control-converter-pass.fds`. Read from the input: six meshes, `&REAC FUEL='N-HEXANE'` (more than one tracked species). Mesh 3, `IJK=30,30,20` on `-0.15,0.15,-0.15,0.15,0,0.2`, has cell size 0.01; meshes 1, 2, 4, 5, 6 have cell size 0.05, so the ratio is 5 and the coarse wall cell has NIC = 5 x 5 = 25. Derived by hand from the mesh geometry (not yet run):

| Mesh | Face on mesh 3 | NOM | Offending wall cells | NIC | IOR of the first |
|---|---|---|---|---|---|
| 1 | high x | 3 | 24 | 25 | -1 |
| 2 | low x | 3 | 24 | 25 | 1 |
| 4 | low z | 3 | 36 | 25 | 3 |
| 5 | high y | 3 | 24 | 25 | -2 |
| 6 | low y | 3 | 24 | 25 | 2 |

Mesh 3 (the fine side) has NIC=1 on all faces and prints nothing. Expected result of `mpirun -np N fds_amr amr_nic_guard_neg_race_test_1.fds` (N = 1, 2, 6) on the patched `USE_AMREX=ON` build:
- stderr: one `ERROR(9001): ...` line per offending mesh (five with N=6, one per rank that holds a mesh from the table; with N=1 only the first mesh processed, then the rest are silent because `SHUTDOWN` prints only while `STOP_STATUS` is not yet `SETUP_STOP`), each with the numbers of the table row (the `IW` and `IOR` are those of the lowest-numbered offending wall of that mesh; the counts 24/24/36/24/24 and NIC=25 are fixed by the geometry), followed once (rank 0) by `ERROR: FDS was improperly set-up - FDS stopped (CHID: amr_nic_guard_neg_race_test_1)`;
- no `level 0: ... box(es)` line on stdout (the driver is never reached);
- all ranks end within the 300 s test timeout (no hang), exit status 0;
- on an `USE_AMREX=OFF` build the same input runs normally (no guard), which confirms the macro guard.

## Converter-passing negative control

**Why a second negative control.** The first control (`inputs/0010-neg-control-race_test_1.fds`) never reaches the guard through the stock `fds_amr`. `main.cpp:150` calls `prepare_amr_input` before `fds_setup(0,...)` (main.cpp:155). The converter rejects mesh 3 (cell size 0.01 against 0.05, ratio 5): `refinement ratio 5 between mesh 1 and mesh 3 is not supported; AMR mode supports ratios 2 and 4 (FR-010)` (`Hierarchy.cpp:170-171`, called from `group_meshes` at `InputConverter.cpp:296`). It is still a valid control for the guard itself on a build where the pre-pass is bypassed (see "Fallback" below).

**Can any plain mesh layout pass the converter and keep a level-0 NIC>1 face? No (by design).** Read from the converter source:
- Every mesh is classified only by `(XB range)/IJK` (`Hierarchy.cpp:24-34`, `InputConverter.cpp:165-166`). A mesh whose cell size is not an integer fraction (ratio 2 or 4) of the coarsest one is an error (`Hierarchy.cpp:203-217`, tolerance `kTol = 1e-3` cells, `Hierarchy.cpp:13`); with unequal cell sizes and no `&AMR` line it is an error too (`Hierarchy.cpp:173-177`); a touching pair with a ratio other than 2 or 4 is an error (`Hierarchy.cpp:165-171`). A mesh that is accepted at a finer level is removed from the level-0 text and replaced by cover meshes of the level-0 cell size (`InputConverter.cpp:355-372`).
- What is left on level 0 has the level-0 cell size, lower corners on the level-0 lattice (`Hierarchy.cpp:244-249` and `261-266`), no overlaps (`Hierarchy.cpp:223-227`) and no gaps (`Hierarchy.cpp:285-296`). The working tree also re-reads the converted text and checks that equality in `verify_equal_level0` (`InputConverter.cpp:374-390`, `395-439`; an uncommitted addition at the time of writing). With equal, lattice-aligned meshes, the two sample points that FDS puts at 0.025 and 0.975 of a wall cell width (`init.f90:3147-3170`; `(ITER*0.95-0.475)*(X(I)-X(I-1))` at `init.f90:3161`) always land in one cell of the neighbour, so `EWC%NIC = 1` (`main.f90:2184`). To split the two samples, the neighbour's cell edges would have to be offset by about 0.025 of a cell or more; the lattice test allows 1e-3 of a cell (`Hierarchy.cpp:13`, `263`). A shift in a single-cell direction (where the lattice is not tested, `Hierarchy.cpp:259`) would make the two samples find different meshes and stop FDS with `ERROR(431)` (`init.f90:3185`) before the guard, so it does not reach it either.
- So the guard cannot fire through the converter for any layout described by `IJK` and `XB` alone. It is a defence in depth for what the converter cannot see.

**What the converter cannot see: `&TRNX`/`&TRNY`/`&TRNZ` grid stretching.** The converter reads only `IJK`, `XB`, `MPI_PROCESS`, `ID`, `MULT_ID` of a `&MESH` line (`InputConverter.cpp:165-174`) and no `TRN*` key or group anywhere (no `TRN` in any file of `Source/regrid_transport`, tests included). FDS builds the physical cell edges from the transformation (`read.f90:1184-1216` linear form, `func.f90:6647-6676` `G`), so a mesh with six equal computational cells can have physical cells of different widths. Layout below: the race_test_1 six-mesh domain, mesh 3 given the coarse cell size (`IJK=6,6,4`, 0.05, so no ratio-5 pair), and mesh 5 stretched in x so that a part of its y-face meets mesh 3/4 across cells of a different width.

File: `inputs/0010-neg-control-converter-pass.fds` (copy of `inputs/0010-neg-control-race_test_1.fds`; `CHID` and `TITLE` changed, `T_END=0.02`, the six `&MESH` lines kept except mesh 3 `IJK=6,6,4` and `TRNX_ID='XSTRETCH'` on mesh 5, three `&TRNX` lines, and the line `&AMR MAX_LEVEL=0, BLOCKING_FACTOR=2 /`). The `&AMR` line is needed because the 34 x 18 level-0 domain is not divisible by the default blocking factor 8 (`Hierarchy.cpp:377-385`: error "not divisible by BLOCKING_FACTOR 8 on level 0"); `MAX_LEVEL=0` adds no refinement, and `BLOCKING_FACTOR` is read at `AmrInput.cpp:209`. The rest of the input (`&REAC` N-HEXANE, `&SURF`, `&OBST`, `&VENT`, `&DEVC`) is unchanged.

| Mesh | `IJK` | Extent (x; y; z) | Cell size seen by the converter | Physical cells in FDS |
|---|---|---|---|---|
| 1 | 14,18,32 | -0.85..-0.15; -0.45..0.45; 0..1.6 | 0.05 | 0.05 |
| 2 | 14,18,32 | 0.15..0.85; -0.45..0.45; 0..1.6 | 0.05 | 0.05 |
| 3 | 6,6,4 | -0.15..0.15; -0.15..0.15; 0..0.2 | 0.05 | 0.05 |
| 4 | 6,6,28 | -0.15..0.15; -0.15..0.15; 0.2..1.6 | 0.05 | 0.05 |
| 5 | 6,6,32, `TRNX_ID` | -0.15..0.15; -0.45..-0.15; 0..1.6 | 0.05 | y, z 0.05; x edges -0.15, -0.10, -0.075, -0.05, 0.05, 0.10, 0.15 (widths 0.05, 0.025, 0.025, 0.10, 0.05, 0.05) |
| 6 | 6,6,32 | -0.15..0.15; 0.15..0.45; 0..1.6 | 0.05 | 0.05 |

The `&TRNX` points (linear form, interior break points only, absolute coordinates; `read.f90:1184-1198` subtracts the mesh start): `(CC,PC) = (-0.10,-0.10), (0.00,-0.05), (0.05,0.05)`, slopes 1, 0.5, 2, 1. The end cells keep 0.05, so the faces of mesh 5 to meshes 1 and 2 are not changed.

**Why the converter accepts it (file:line).** Checked with the converter library of the committed tip (`adbf1225e2`, no `verify_equal_level0`) and with the working tree (with it), both built as the `fds_amr_convert_input` tool (`ConvertInputTool.cpp`; it calls `convert_input`, the function main.cpp:97 calls): exit 0, `meshes 6, level 0 meshes 6, finer meshes removed 0 (0 lines), cover meshes added 0, levels 1`, and one warning (FR-010 coarsening: the domain 34x18x32 can be halved once, `Hierarchy.cpp:301-317`). Steps:
1. `find_group_spans` (`InputConverter.cpp:104-148`) finds all groups; `&TRNX` groups are found as spans but no key of them is parsed (only MESH, MULT, VENT, AMR, AMR_REGION, MISC and the `MESH_ID` test of other groups are read, `InputConverter.cpp:157`, `205-228`, `277-291`, `318-334`). `parse_meshes` (`InputConverter.cpp:150-203`) gives six meshes from `IJK`/`XB`. `TRNX_ID` is never read.
2. `group_meshes` (`Hierarchy.cpp:81-299`): cell size 0.05 in x, y, z for every mesh (`Hierarchy.cpp:31`), so `dx0 = 0.05` (`105-107`); no pair differs (`128-136` `continue`), `any_finer` stays false, so the missing-`&AMR` rule `173` does not apply either way; every mesh has ratio 1 and level 0 (`191-218`); no overlap (`223-227`); level-0 domain 34 x 18 x 32 on the lattice (`244-249`); index boxes aligned (`252-281`); no gaps (`285-296`).
3. `build_static_hierarchy` accepts blocking factor 2 (34 and 18 are divisible by 2; the error text is at `Hierarchy.cpp:377-385`). Without the `&AMR` line it fails there with the default 8; the committed-tip tool and the working-tree tool give the same result.
4. In the working tree, `verify_equal_level0` (`InputConverter.cpp:389`) sees six meshes of cell size 0.05 on the lattice, disjoint, tiling 34 x 18 x 32: pass.
5. `level0_text` equals the input except that the `&AMR` line becomes a `!` comment line (`emit_level0_text`, `InputConverter.cpp:230-268`), so `prepare_amr_input` writes `<stem>_amr_level0.fds` into the working directory (main.cpp:105-113) and FDS reads that file (main.cpp:155). The `&TRNX` lines and `TRNX_ID` go through unchanged. The driver prints `FDS-AMReX: input converter: 6 meshes, 6 level-0 meshes, 0 finer meshes removed, 0 cover meshes added, hierarchy levels 1` (main.cpp:152) before the set-up.

**FDS accepts the layout, and what an unpatched build does.** The converted file was run (`mpirun -np 1`) on an existing unpatched driver build (revision string `FDS-6.11.1-1244-g36975d765f-FDS-AMReX`; it has no converter, so the converted file was given to it, which is what the stock build passes to `fds_setup`). FDS set-up completed through `INITIALIZE_MESH_EXCHANGE_1`, the Poisson initialisation and the first mesh outputs without `ERROR(431)` or any other error; then the driver stopped with `amrex::Abort::0::M2a: TRNX/TRNY/TRNZ (nonuniform) meshes are not supported (IR-002)` (`FdsSetup.cpp:31`, after `fds_setup` returned). On a patched build the guard fires inside `fds_setup`, before that line. Two side findings: (1) the driver's TRN test (`fds_mesh_query.f90:38`) looks at `TRNX_ID` of the `&MESH` line only, so a `&TRNX` line that selects its mesh by `MESH_NUMBER` (`read.f90:1034-1038`) escapes it; (2) the input converter does not read `TRN*` at all, so it can neither refuse nor convert such inputs. A converter rule "refuse any input with a `&TRNX`, `&TRNY` or `&TRNZ` group in AMR mode" (owner Role 3) would close this route and would make the guard unreachable from the pre-pass for good.

**NIC by hand (face by face).** Wall cells are numbered per mesh in the order of `init.f90:76-107`: `IOR=1` (low x, `I=0`), `IOR=-1`, `IOR=2` (low y, `J=0`), `IOR=-2`, `IOR=3`, `IOR=-3`; within a face K is the outer loop and the first tangential index the inner one. For a wall cell the two sample points are the wall-cell centre moved by -0.475 and +0.475 of the wall cell width in each tangential direction and by `MESH_SEPARATION_DISTANCE` (1e-3, `read.f90:916`) into the neighbour (`init.f90:3147-3170`); `NIC` = (number of neighbour cells spanned in I) x (in J) x (in K) (`main.f90:2184`).
- Mesh 3 and 4 low-y face (`IOR=2`, y = -0.15) against mesh 5 (cells in x: mesh 3/4 cell `i` is [-0.15+0.05(i-1), -0.15+0.05i]). `i=1` [-0.15,-0.10]: samples -0.14875 and -0.10125, both in mesh 5 cell 1 [-0.15,-0.10]: NIC=1. `i=2` [-0.10,-0.05]: samples -0.09875 and -0.05125, in mesh 5 cells 2 [-0.10,-0.075] and 3 [-0.075,-0.05]: spans 2 cells in I, NIC = 2 x 1 x 1 = **2**. `i=3` [-0.05,0] and `i=4` [0,0.05]: both samples in mesh 5 cell 4 [-0.05,0.05]: NIC=1. `i=5,6`: mesh 5 cells 5, 6: NIC=1. K: all cells have dz = 0.05 aligned, span 1. So the offenders of mesh 3 are `i=2`, `k=1..4` (4 cells) and of mesh 4 `i=2`, `k=1..28` (28 cells).
- Mesh 5 high-y face (`IOR=-2`, y = -0.15) against meshes 3 (z < 0.2) and 4 (z > 0.2). Cells `i=1,2,3` (widths 0.05, 0.025, 0.025): both samples of each lie in one neighbour cell: NIC=1 (cell 1 [-0.15,-0.10] is one neighbour cell; cells 2 and 3 lie inside neighbour cell [-0.10,-0.05]). Cell `i=4` [-0.05,0.05]: samples -0.0475 and +0.0475, in neighbour cells [-0.05,0] and [0,0.05]: NIC = **2**. Cells 5, 6: NIC=1. The offenders of mesh 5 are `i=4`, `k=1..4` (neighbour mesh 3, 4 cells) and `k=5..32` (neighbour mesh 4, 28 cells): 32 cells.
- Alignment tests do not trip: the checks at `init.f90:3194-3221` compare the spanned neighbour width with the wall cell width (0.05 against 0.05 for mesh 3/4 `i=2`; 0.10 against 0.10 for mesh 5 `i=4`) and the mixed-ratio tests (`3214-3221`) need the other direction to differ too, which it does not (dz equal).
- All other faces have NIC=1 (meshes 1, 2 and 6 have none; the x-faces of mesh 5 meet mesh 1 and 2 through the y and z cells, which are 0.05 on both sides; mesh 3 to mesh 4 is an equal-cell z-face).
The same table was reproduced by a small script that applies the sample-point rule to every external wall cell of the six meshes in FDS order (neighbour found by the first mesh that contains the point, bounds inclusive, as `func.f90:5424`): NIC>1 on exactly 4 + 28 + 32 = 64 external wall cells, 0 elsewhere, NOM>0 walls per mesh 576, 576, 132, 708, 576, 576 (3144 in all; this total is not printed on this control because the guard stops the run first).

**Expected offender table (the lines `ERROR(9001)` prints; first offending wall of each mesh, count over the mesh):**

| Mesh | First offending wall `IW` | `IOR` | Abuts (`NOM`) | `NIC` | Offending wall cells in the mesh |
|---|---|---|---|---|---|
| 3 | 50 (= 48 + 2: face `IOR=2` starts at 49, `k=1`, `i=2`) | 2 | 5 | 2 | 4 |
| 4 | 338 (= 336 + 2: face `IOR=2` starts at 337) | 2 | 5 | 2 | 28 |
| 5 | 580 (= 576 + 4: face `IOR=-2` starts at 577, `k=1`, `i=4`) | -2 | 3 | 2 | 32 |

The face start indices follow from the face sizes: mesh 3 (6,6,4): 24 per x-face, 24 per y-face, so `IOR=2` starts at 49; mesh 4 (6,6,28): 168 per x-face, so `IOR=2` starts at 337; mesh 5 (6,6,32): 192 per x-face and 192 per y-face, so `IOR=-2` starts at 577. Meshes 1, 2, 6 print nothing. Expected `err.txt` lines (rank suffix and `CHID` are appended by `SHUTDOWN`; with one rank only the first mesh in rank order, mesh 3, prints, because `SHUTDOWN` prints only while `STOP_STATUS` is not yet `SETUP_STOP`; with six ranks all three lines appear):

```
ERROR(9001): AMR mode needs equal-resolution level-0 meshes. MESH 3: external wall cell 50 (IOR=2) abuts MESH 5 with NIC=2 (4 such wall cells in this mesh). Finer meshes are converted into AMR levels by the input pre-pass; run this input without AMR.
ERROR(9001): AMR mode needs equal-resolution level-0 meshes. MESH 4: external wall cell 338 (IOR=2) abuts MESH 5 with NIC=2 (28 such wall cells in this mesh). Finer meshes are converted ...
ERROR(9001): AMR mode needs equal-resolution level-0 meshes. MESH 5: external wall cell 580 (IOR=-2) abuts MESH 3 with NIC=2 (32 such wall cells in this mesh). Finer meshes are converted ...
```

followed once by `ERROR: FDS was improperly set-up - FDS stopped (CHID: amr_nic_guard_neg_conv_pass)`; no `NIC_CHECK` line; no `level 0:` line on stdout; the converter line above is on stdout; exit status 0. Not run on a patched build; numbers are derived by hand and by the script.

**Run (same form as the first control, on 1, 2 and 6 ranks):**
`mkdir -p <w> && cp <docs>/upstream-patches/inputs/0010-neg-control-converter-pass.fds <w>/amr_nic_guard_neg_conv_pass.fds && cd <w> && timeout 300 mpirun --bind-to none --oversubscribe -np <N> <bld>/fds_amr amr_nic_guard_neg_conv_pass.fds > out.txt 2> err.txt; echo $?`
Pass: the lines of the table above in `err.txt` (at least one; one per rank that holds an offending mesh), `grep -c "input converter: 6 meshes, 6 level-0 meshes" out.txt` is 1, `grep -c "level 0:" out.txt` is 0, no `NIC_CHECK`. Fail: `M2a: TRNX/TRNY/TRNZ` in `err.txt` (guard missing from the build), or `input converter ERROR` (the converter learned to refuse `TRN*`; then use the fallback).

**Fallback: scratch-build bypass of the pre-pass (for `inputs/0010-neg-control-race_test_1.fds` and any input the converter refuses).** In a scratch copy of the tree only, in `Source/driver/main.cpp` replace line 150

`if (!prepare_amr_input(argv[1], conv, fds_input_s)) { MPI_Abort(MPI_COMM_WORLD, 1); return 1; }`

by the single line

`(void)conv;`

(`fds_input_s` then stays `argv[1]`, the original file; main.cpp:151 and 155 are unchanged, the line at 152 prints `0 meshes`). `conv` is then a default `ConvertResult` and is unused. The converter is skipped completely and `fds_setup(0,...)` reads the original input, so the ratio-5 control reaches `INITIALIZE_MESH_EXCHANGE_1` and the table of the first control applies. This build is a test instrument only; the reference tree and the shipped driver keep the pre-pass. If a later change makes the converter refuse `TRN*` inputs (see above), this bypass is the only way to exercise the guard from a stock-input pair, and the guard stays a defence in depth: the pre-pass makes it unreachable by design for inputs it can read.

**Positive controls (guard silent; equal-resolution multi-mesh with 2 species):**
- `shunn3_4mesh_32` (`Verification/Scalar_Analytical_Solution`, 4 meshes of `IJK=16,1,16`, cell size 0.125 in every mesh, `&SPEC` BACKGROUND and SCALAR), 4 ranks: `Source/driver/tests/check_setup_amr.sh <build-dir> shunn3_4mesh_32 4 4` must print PASS (`level 0: 4 box(es)` and `STOP: Set-up only`).
- `Source/driver/tests/cases/dec4_np4.fds` (16 meshes of `IJK=8,1,8`, all cell size 0.0625, two `&SPEC`), 4 ranks: must print `level 0: 16 box(es)`, no `ERROR(9001)`.
- Not a positive control: `race_test_1_r4.fds` (regrid_transport cases) still has a ratio-4 face (mesh 3 is `24,24,16` on 0.3 x 0.3 x 0.2, cell size 0.0125 against 0.05); the guard also aborts it.

## Apply check
`git -C <src> show HEAD:Source/main.f90 > main.f90` in a scratch directory with the layout `Source/main.f90`, then `patch -p1 --dry-run < 0010-main-amr-level0-nic-guard.patch`: clean, no offset, no fuzz, four hunks (re-checked after the `NIC_CHECK` hook was added; the patched file is byte-identical to the file the patch was cut from, 5278 lines). Line numbers in this note are those of the original file unless marked "patched". Fortran lines are at most 132 characters (the new lines were checked); the file compiles under the default free-form limit.

## Validation request (GNU Build Chief and Intel Build Chief)
Do not edit the reference tree. Work on a copy of the FDS-AMReX tree outside it.
1. Copy the tree, apply: `cd <copy> && patch -p1 < <docs>/upstream-patches/0010-main-amr-level0-nic-guard.patch`. Expect no rejects.
2. Configure and build the AMR driver with your toolchain exactly as for patches 0003-0009 (`cmake -S <copy> -B <bld> -DUSE_AMREX=ON ...`, which compiles the Fortran sources with `-DWITH_AMREX`; for GNU use the options in `Source/driver/tests/env.sh`, for oneAPI the ones used for 0009). Report any warning from `main.f90` (expect none; with `-check all -traceback -fpe0 -init=snan` on ifx and `-fcheck=all` with FP traps on Debug gfortran). Also build `USE_AMREX=OFF` and run `Source/driver/tests/check_off_bitwise.sh <copy> <bld-off> shunn3_4mesh_32 4` (expect PASS: output bitwise identical to the baseline).
3. Negative control, on 1, 2 and 6 ranks:
   `mkdir -p <w> && cp <docs>/upstream-patches/inputs/0010-neg-control-race_test_1.fds <w>/amr_nic_guard_neg_race_test_1.fds && cd <w> && timeout 300 mpirun --bind-to none --oversubscribe -np <N> <bld>/fds_amr amr_nic_guard_neg_race_test_1.fds > out.txt 2> err.txt; echo $?`
   Expect: `grep -c "ERROR(9001)" err.txt` is at least 1 and each line matches the table above; `grep "improperly set-up" err.txt` once; `grep -c "level 0:" out.txt` is 0; the command returns well inside the timeout (a timeout, exit 124, is a FAIL: hang); exit status 0.
3b. Converter-passing negative control, stock `fds_amr`, 1, 2 and 6 ranks: the command and pass criteria in "Converter-passing negative control" (Run). Expect the three `ERROR(9001)` lines of its table, the converter line on stdout, no `level 0:` line. For the ratio-5 input of step 3 use a scratch build with the one-line bypass of the same section (the stock pre-pass stops it first).
4. Same input on the `USE_AMREX=OFF` `fds` build, `-np 6`: no `ERROR(9001)`; it may be stopped after a few steps, only the start is of interest.
5. Positive controls: `check_setup_amr.sh <bld> shunn3_4mesh_32 4 4` (PASS) and `dec4_np4.fds` on 4 ranks (`level 0: 16 box(es)`, no `ERROR(9001)`).
6. `NIC_CHECK` line (D-077 (2), V&V P-3): on the `USE_AMREX=ON` build, `grep -c "NIC_CHECK level0:" stderr.txt` must be exactly 1 on each positive control and the line must read `NIC_CHECK level0: 256 walls checked` for `shunn3_4mesh_32` on 4 ranks and `NIC_CHECK level0: 512 walls checked` for `dec4_np4.fds` on 4 ranks (hand counts above); the number must not change with the rank count (also run `shunn3_4mesh_32` on 1 and 2 ranks: 256 each time). On the negative control (1, 2 and 6 ranks) there must be no `NIC_CHECK` line, in `stderr.txt` or `stdout.txt`. On the `USE_AMREX=OFF` build there is no such line. If a count differs, report the case, the rank count and the number printed.
7. Report: toolchain and flags, the output of steps 3 to 5 (the `err.txt` lines), and whether the message text and numbers match the table. If a rank hangs on 2 or 6 ranks, report the rank that has not reached `STOP_CHECK`.
