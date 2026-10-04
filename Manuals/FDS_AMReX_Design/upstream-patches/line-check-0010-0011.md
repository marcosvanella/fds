# Line-number check of output patches 0010 and 0011

Scope: the two patches that Role 1 item 7 (`next-work.md`) waits for. "0010" is `0010-main-amr-level0-nic-guard.patch` (note `0010-main-amr-level0-nic-guard.md`). "0011" is `UP-0011-pressure-rhs-dump.patch` (note `UP-0011-pressure-rhs-dump.signoff.md`, the pressure RHS dump). Reference tree: the FDS-AMReX branch, read with `git show HEAD:<path>`; `Source/main.f90` has blob `2e1327bbe8` there (the same blob as at tips 991a759f79 and 0353957a93). "Original" line numbers are those of the unpatched file, "patched" those after applying the patch. Scratch copies were used for every apply test and removed.

## Verdict
| Patch | Applies | Line numbers in the patch hunks | Line numbers in the note | Corrections needed |
|---|---|---|---|---|
| 0010 | clean on FDS-AMReX HEAD, no offset, no fuzz | all hunk headers correct | all citations correct except one number | 1 (preprocessor growth, below) |
| 0011 (UP-0011) | clean on FireX 36975d765f; offsets +73 on the FDS-AMReX HEAD and on bee11f0329; offsets -131 (main.f90) and -4 (pres.f90 hunk 2) on master ce1f659cd4 | all hunk headers correct | routine name and patch number need correcting | 3 (below) |

## Patch 0010
Apply: `patch -p1 --dry-run` on `Source/main.f90` at FDS-AMReX HEAD: "checking file Source/main.f90", three hunks, no offset message. The patched file has 5267 lines (original 5226, 41 added).

Lines checked (original unless marked patched), all match the note:
| Citation in the note | Checked content | Result |
|---|---|---|
| main.f90:326-328 | `DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX` / `CALL INITIALIZE_MESH_EXCHANGE_1(NM)` / `ENDDO` | match |
| patched main.f90:330-334 | `#ifdef WITH_AMREX` (330), `CALL STOP_CHECK(1)` (333), `#endif` (334) | match |
| main.f90:334 (implied: STOP_CHECK before `MPI_INITIALIZATION_CHORES(4)`) | `CALL MPI_INITIALIZATION_CHORES(4)` is original line 334 | match |
| main.f90:303-308 | `IF (SETUP_ONLY .OR. CHECK_MESH_ALIGNMENT .OR. STOP_STATUS/=0)` (303) ... `CALL STOP_CHECK(1)` (306) `ENDIF` (307) | match (the range includes one trailing blank line) |
| main.f90:2132-2374 | routine `INITIALIZE_MESH_EXCHANGE_1` starts at 2132, `END SUBROUTINE` at 2374; the only `MPI_` token in the range is the logical `MPI_PARTICLE_EXCHANGE` (2360), no MPI call | match |
| patched main.f90:2141-2149 | `USE COMP_FUNCTIONS, ONLY: SHUTDOWN` block (2141-2143) and the `AMR_N_BAD,AMR_IW1` / `AMR_MESSAGE` block (2146-2149) | match |
| patched main.f90:2387-2413 | `#ifdef WITH_AMREX` (2387) ... `#endif` (2413); `END SUBROUTINE INITIALIZE_MESH_EXCHANGE_1` follows at 2415 | match |
| "inserted before original line 2374" | original 2374 is the `END SUBROUTINE`; the insertion follows the `MPI_PARTICLE_EXCHANGE` block (2360-2371) | match |
| main.f90:2184 | `EWC%NIC = (EWC%IIO_MAX-EWC%IIO_MIN+1)*...` (2184); set inside the `NOM` search over `EXTERNAL_WALL(IW)` | match |
| main.f90:1969-1993 | `STOP_CHECK`: `MPI_ALLREDUCE(...STOP_STATUS...MPI_MAX)` (1976), `IF (END_CODE==1) CALL END_FDS` (1991) | match |
| main.f90:2073 | the plain `STOP` of `END_FDS` (routine 1998-2075) | match |
| `N_REQ1`... are 0 at that point | main.f90:108 `N_REQ1=0,...`; freed at 2055-2064 | match |
| func.f90:120-139 | `SUBROUTINE SHUTDOWN` (120) to `END SUBROUTINE SHUTDOWN` (139); only prints and sets `STOP_STATUS = SETUP_STOP` | match |
| init.f90:3185-3187 | `ERROR(431)` write (3185), `STOP_STATUS = SETUP_STOP` (3186), `IERR = 1` (3187) | match |
| fds_fine_level.f90:162-165 | `M%EXTERNAL_WALL(0)` (163), `M%N_EXTERNAL_WALL_CELLS=0` (165) | match |
| wall.f90:897 | `COARSE_MESH_IF: IF (EWC%NIC>1) THEN` | match |
| divg.f90:219 | `IF (EWC%NIC>1) THEN` (comment "overwrite coarse mesh diffusive flux") | match |
| FdsAmr.cpp:30 | `amrex::Abort("M2a: all meshes must have the same cell size ...")` | match |
| "code 9001 is unused" | no `ERROR(9001)` or any four-digit `ERROR(90xx)` in `Source`; the nearest are `ERROR(901)` to `ERROR(910)` | match |
| Fortran lines at most 132 characters (new lines) | patched lines 2384-2415 and 330-334: none over 132; the two lines over 132 in the file (patched 121 and 122, original 121 and 122) are older | match |
| "`gfortran -cpp -E -P` without `WITH_AMREX`: identical, 4055 lines" | blank lines dropped: both files give 4055 and `diff` is empty | match |

Mismatch (1):
- The note states that with `-DWITH_AMREX` the preprocessed file "grows by 80 lines". Measured with `gfortran -cpp -DWITH_AMREX -E -P`: original 5169 lines, patched 5204 lines (all lines, +35; non-blank lines +29). Correct value: +35 lines (+29 without blank lines). No behaviour consequence.

Not checked (needs a build, left to the Build Chiefs as the note says): the compile, the message text on 1, 2 and 6 ranks, the hand-derived counts 24/24/36/24/24 and NIC=25 of the negative control.

## Patch 0011 (UP-0011-pressure-rhs-dump.patch)
Apply, hunk by hunk (`patch -p1 --dry-run`, files `Source/main.f90` and `Source/pres.f90` exported with `git show`):
| Base | main.f90 hunks (-1650,6 and -1669,6) | pres.f90 hunks (-9,7 and -1076,6) |
|---|---|---|
| FireX 36975d765f | clean, no offset | clean, no offset |
| FDS-AMReX HEAD | both succeed at +73 (1723, 1748) | clean, no offset |
| bee11f0329 | both succeed at +73 (1723, 1748) | clean, no offset |
| master ce1f659cd4 | both succeed at -131 (1519, 1544) | hunk 2 at -4 (1072); hunk 1 clean |

The offsets +73 and -131 agree with the sign-off note; the -4 offset of pres.f90 hunk 2 on master is not in the note.

Lines checked (FDS-AMReX HEAD, patched file):
| Citation | Checked content | Result |
|---|---|---|
| hunk 1, after the `PRESSURE_SOLVER_COMPUTE_RHS` loop and before `SELECT CASE(PRES_FLAG)` | original main.f90:1723 `CALL PRESSURE_SOLVER_COMPUTE_RHS(T,DT,NM)`; new call block patched 1729; `SELECT CASE(PRES_FLAG)` follows (patched 1734, original 1728) | match |
| hunk 2, after the `END SELECT` and before the residual check | new call patched 1754; the comment "Check the residuals of the Poisson solution" is original 1745 | match |
| "`PRESSURE_ITERATION_LOOP`" (note, item 1) | `PRESSURE_ITERATION_LOOP` (main.f90:1696) is the label of the `DO` construct inside `SUBROUTINE PRESSURE_ITERATION_SCHEME` (main.f90:1674); both blocks are inside it | match, but see correction 1 |
| 1 line in the PUBLIC list | pres.f90:11-12, `PRESSURE_SOLVER_DUMP` added after `COMPUTE_VELOCITY_ERROR` | match |
| routine added at the end of module PRES, "about 190 lines" | hunk -1076,6 +1076,194: 188 lines added before `END MODULE PRES` (original pres.f90:1079); the routine is patched 1090-1244, the helper follows | match |
| added lines at most 132 characters | none longer | match |
| driver files `pdump_<icyc>_<P|C>_{rhs,phi}.bin` of Role 1 | `FDSTL_PDUMP` at Source/driver/TimeLoop.cpp:834, file names at 841 | match |
| reader `Source/pressure_backend/harness/m4_modes.cpp` (`pb_fds_frozen`) and `frozen/fds_csmag32_periodic/README.md` | both files exist; `pb_fds_frozen` named at m4_modes.cpp:2 | match |

Mismatches and corrections:
1. Number: the patch header says `[PATCH 0009]`, the sign-off note is titled "patch 0009", and `Source/pressure_backend/frozen/fds_csmag32_periodic/README.md:10` says "upstream patch 0009", while the file and the index (`README.md` row UP-0011) use UP-0011. Patch 0009 also exists in the driver series (`Source/driver/patches/0001`-`0009`), so "0009" is ambiguous. Correct all three to "UP-0011". The first two are in this folder; the third is in the source tree (hand to its owner, Role 2).
2. Wording: the note names `PRESSURE_ITERATION_LOOP` as if it were a routine. It is the loop label inside `PRESSURE_ITERATION_SCHEME` (main.f90:1674, loop at 1696). Suggested text: "`PRESSURE_ITERATION_SCHEME`, inside the `PRESSURE_ITERATION_LOOP` construct".
3. Offset: add "pres.f90 hunk 2 at offset -4 on master ce1f659cd4" to the applicability statement.

Not checked: the behaviour-unchanged runs (csmag_32, shunn3_4mesh_32, ns2d_16), which the note reports and which need builds; the compile of the added routine against the HEAD `pres.f90` module scope.

## Resolution
- Patch 0010: the preprocessor growth statement in the note is corrected; after the `NIC_CHECK` hook was added the measured growth with `-DWITH_AMREX` is 46 lines (5169 to 5215, blank lines kept), and the no-`WITH_AMREX` result is unchanged (4055 lines, identical).
- UP-0011: the patch header and the sign-off note title now say UP-0011; the wording on `PRESSURE_ITERATION_LOOP` (a loop label in `PRESSURE_ITERATION_SCHEME`) and the pres.f90 hunk 2 offset of -4 on master are in the note. Open outside this folder: `Source/pressure_backend/frozen/fds_csmag32_periodic/README.md:10` (source tree, owner Role 2) still says "upstream patch 0009".
