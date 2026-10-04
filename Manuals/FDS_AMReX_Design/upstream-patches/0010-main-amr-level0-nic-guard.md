# Patch 0010 (DRAFTED - needs Build Chief validation): `main.f90`, AMR level-0 abort guard for coarse/fine faces (D-065 Q1, D-076)

**Status: drafted, not built.** Cut against `Source/main.f90` of the FDS-AMReX branch tree (blob `2e1327bbe8`, identical at tips 991a759f79 and 0353957a93). `patch -p1 --dry-run` is clean on that file (see "Apply check"). Nobody has compiled it yet; the GNU Build Chief and the Intel Build Chief validate it (section "Validation request"). It continues the driver patch series (`Source/driver/patches/0001`-`0009`); the file is kept here because the source trees are read-only for the author. Apply from the top of the source tree with `patch -p1` or `git apply`.

## What it does
Three additions in `Source/main.f90`, all inside `#ifdef WITH_AMREX` (no existing statement is edited):
1. `INITIALIZE_MESH_EXCHANGE_1`, declaration part (patched main.f90:2141-2149): `USE COMP_FUNCTIONS, ONLY: SHUTDOWN` and the locals `AMR_N_BAD`, `AMR_IW1`, `AMR_MESSAGE`.
2. `INITIALIZE_MESH_EXCHANGE_1`, end of the routine, after the `MPI_PARTICLE_EXCHANGE` block and before `END SUBROUTINE` (patched main.f90:2387-2413; inserted before original line 2374, the `END SUBROUTINE`): a loop over `M%EXTERNAL_WALL(1:M%N_EXTERNAL_WALL_CELLS)`. It skips walls with `EWC%NOM<1` (NIC is set only for walls that abut another mesh, main.f90:2184) and counts the walls with `EWC%NIC>1`. If any, it writes the message below and calls `SHUTDOWN(...,PROCESS_0_ONLY=.FALSE.)` once for the mesh (first offending wall is named, the count is given).
3. The call loop in the main set-up sequence (patched main.f90:330-334; the loop itself is original 326-328): `CALL STOP_CHECK(1)` right after the `DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX / CALL INITIALIZE_MESH_EXCHANGE_1(NM)` loop.

## Why
D-065 Q1 retires the `EWC%NIC>1` branches (species diffusive flux match, `wall.f90:897`; coarse flux overwrite, `divg.f90:219`) in the AMR route. D-072 (f) and D-076 make that a design requirement: the AMR route never has an interpolated boundary between level-0 meshes of different resolution (a coarse/fine jump is a refinement level, handled by ghost fill and the interface flux overwrite). Reading the code cannot prove the whole route, so the requirement is enforced at run time: any level-0 mesh with a NIC>1 external wall aborts the set-up with a clear message, before any exchange or time step. Without the guard such an input would reach the retired branches silently. Until the input pre-pass (converter) is wired, inputs with finer meshes stop here instead of at `FdsAmr.cpp:30` (which checks the cell sizes only after `FDS_SETUP` returns).

## Message
Written to the error unit by every rank that has an offending mesh (`SHUTDOWN` with `PROCESS_0_ONLY=.FALSE.` prefixes a blank line and appends ` (MPI Process:<rank>, CHID: <chid>)`):

`ERROR(9001): AMR mode needs equal-resolution level-0 meshes. MESH <NM>: external wall cell <IW> (IOR=<IOR>) abuts MESH <NOM> with NIC=<NIC> (<count> such wall cells in this mesh). Finer meshes are converted into AMR levels by the input pre-pass; run this input without AMR.`

It names the mesh (`NM`), the wall (`IW`, `IOR`, other mesh), NIC, and states the requirement and the conversion route. Code 9001 is unused in the tree (checked: no `ERROR(90xx)`).

## Multi-rank behaviour (cannot hang)
- `INITIALIZE_MESH_EXCHANGE_1` contains no MPI call (checked by reading original main.f90:2132-2374; the guard adds none). The only collective the patch adds is the unconditional `STOP_CHECK(1)`, which every rank reaches exactly once, whatever number of meshes it holds, so no rank can wait for a message from a rank that has stopped.
- `SHUTDOWN` (func.f90:120-139) only prints and sets `STOP_STATUS = SETUP_STOP`; it does not stop. `EXTERNAL_WALL` exists on the owning rank only, so the offending rank is not necessarily rank 0.
- The new `STOP_CHECK(1)` (original main.f90:1969-1993) is executed unconditionally by every rank, once, after the loop. With `N_MPI_PROCESSES>1` it does `MPI_ALLREDUCE(MPI_MAX)` of `STOP_STATUS`, so every rank sees `SETUP_STOP`; then `END_FDS` runs on all ranks (rank 0 prints `ERROR: FDS was improperly set-up - FDS stopped`, then `MPI_FINALIZE` and `STOP`). When nothing is wrong `STOP_STATUS` is 0 on all ranks and `STOP_CHECK` returns at once (one extra `MPI_ALLREDUCE` at set-up, only in the `WITH_AMREX` build). This is the pattern of the existing set-up-only stop at original main.f90:303-308.
- Exit status: `END_FDS` ends with a plain `STOP` (original main.f90:2073), so the process exit status is 0, as for every other FDS set-up error (for example `ERROR(431)`, init.f90:3185-3187). A test must key on the message text and on the absence of the driver line `level 0: <n> box(es)`, not on the exit status. In `FDS_SETUP` the `STOP` happens before the driver finalizes AMReX, the same as the existing set-up-only stop; this is accepted for an abort path.

## Where it is a no-op
- `WITH_AMREX` undefined (`USE_AMREX=OFF`): all three hunks are removed by the preprocessor. Checked: `gfortran -cpp -E -P` of main.f90 with the macro undefined, before and after the patch, blank lines dropped: identical (4055 lines). With `-DWITH_AMREX` the preprocessed file grows by 80 lines.
- `WITH_AMREX` defined and no NIC>1 wall (all equal-resolution meshes, one mesh, periodic single mesh, any input with only NIC=1 faces): no message, `STOP_CHECK` returns at once, results unchanged. Fine-level boxes are not in `MESHES` and have `N_EXTERNAL_WALL_CELLS=0` (`fds_fine_level.f90:162-165`), and the routine runs only for `MESHES(NM)`, so it checks level-0 meshes only.
- The guard is independent of the species count. The NIC>1 branches need 2 or more tracked species to matter, but the guard refuses every NIC>1 level-0 face (D-076: any), which is stricter and also covers the velocity and pressure-side interpolation.

## Risk
Low. An input with a coarse/fine face that ran before in the AMR build (it could only get to `FdsAmr.cpp:30` and abort there) now aborts earlier with this message. `STOP_CHECK` is placed before `MPI_INITIALIZATION_CHORES(4)`; `END_FDS` frees `N_REQ1`... requests, which are 0 at that point, as at the set-up stop at original main.f90:303-308 (to be confirmed by the negative run on 2 and 6 ranks).

## How to revert
Delete the three hunks (`patch -R -p1`), or build without `-DWITH_AMREX`. No other file is touched.

## Controls
Files: `inputs/0010-neg-control-race_test_1.fds` (negative); positive controls below.

**Negative control: `inputs/0010-neg-control-race_test_1.fds`** (copy of `Verification/Thread_Check/race_test_1.fds`; only `CHID` and a `TITLE` differ). Read from the input: six meshes, `&REAC FUEL='N-HEXANE'` (more than one tracked species). Mesh 3, `IJK=30,30,20` on `-0.15,0.15,-0.15,0.15,0,0.2`, has cell size 0.01; meshes 1, 2, 4, 5, 6 have cell size 0.05, so the ratio is 5 and the coarse wall cell has NIC = 5 x 5 = 25. Derived by hand from the mesh geometry (not yet run):

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

**Positive controls (guard silent; equal-resolution multi-mesh with 2 species):**
- `shunn3_4mesh_32` (`Verification/Scalar_Analytical_Solution`, 4 meshes of `IJK=16,1,16`, cell size 0.125 in every mesh, `&SPEC` BACKGROUND and SCALAR), 4 ranks: `Source/driver/tests/check_setup_amr.sh <build-dir> shunn3_4mesh_32 4 4` must print PASS (`level 0: 4 box(es)` and `STOP: Set-up only`).
- `Source/driver/tests/cases/dec4_np4.fds` (16 meshes of `IJK=8,1,8`, all cell size 0.0625, two `&SPEC`), 4 ranks: must print `level 0: 16 box(es)`, no `ERROR(9001)`.
- Not a positive control: `race_test_1_r4.fds` (regrid_transport cases) still has a ratio-4 face (mesh 3 is `24,24,16` on 0.3 x 0.3 x 0.2, cell size 0.0125 against 0.05); the guard also aborts it.

## Apply check
`git -C <src> show HEAD:Source/main.f90 > main.f90` in a scratch directory with the layout `Source/main.f90`, then `patch -p1 --dry-run < 0010-main-amr-level0-nic-guard.patch`: clean, no offset, no fuzz. The patched file is byte-identical to the file the patch was cut from. Line numbers in this note are those of the original file unless marked "patched". Fortran lines are at most 132 characters (the new lines were checked); the file compiles under the default free-form limit.

## Validation request (GNU Build Chief and Intel Build Chief)
Do not edit the reference tree. Work on a copy of the FDS-AMReX tree outside it.
1. Copy the tree, apply: `cd <copy> && patch -p1 < <docs>/upstream-patches/0010-main-amr-level0-nic-guard.patch`. Expect no rejects.
2. Configure and build the AMR driver with your toolchain exactly as for patches 0003-0009 (`cmake -S <copy> -B <bld> -DUSE_AMREX=ON ...`, which compiles the Fortran sources with `-DWITH_AMREX`; for GNU use the options in `Source/driver/tests/env.sh`, for oneAPI the ones used for 0009). Report any warning from `main.f90` (expect none; with `-check all -traceback -fpe0 -init=snan` on ifx and `-fcheck=all` with FP traps on Debug gfortran). Also build `USE_AMREX=OFF` and run `Source/driver/tests/check_off_bitwise.sh <copy> <bld-off> shunn3_4mesh_32 4` (expect PASS: output bitwise identical to the baseline).
3. Negative control, on 1, 2 and 6 ranks:
   `mkdir -p <w> && cp <docs>/upstream-patches/inputs/0010-neg-control-race_test_1.fds <w>/amr_nic_guard_neg_race_test_1.fds && cd <w> && timeout 300 mpirun --bind-to none --oversubscribe -np <N> <bld>/fds_amr amr_nic_guard_neg_race_test_1.fds > out.txt 2> err.txt; echo $?`
   Expect: `grep -c "ERROR(9001)" err.txt` is at least 1 and each line matches the table above; `grep "improperly set-up" err.txt` once; `grep -c "level 0:" out.txt` is 0; the command returns well inside the timeout (a timeout, exit 124, is a FAIL: hang); exit status 0.
4. Same input on the `USE_AMREX=OFF` `fds` build, `-np 6`: no `ERROR(9001)`; it may be stopped after a few steps, only the start is of interest.
5. Positive controls: `check_setup_amr.sh <bld> shunn3_4mesh_32 4 4` (PASS) and `dec4_np4.fds` on 4 ranks (`level 0: 16 box(es)`, no `ERROR(9001)`).
6. Report: toolchain and flags, the output of steps 3 to 5 (the `err.txt` lines), and whether the message text and numbers match the table. If a rank hangs on 2 or 6 ranks, report the rank that has not reached `STOP_CHECK`.
