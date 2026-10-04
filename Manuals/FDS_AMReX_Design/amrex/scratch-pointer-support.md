# Scratch pointers, plane aliases and the wall-index table in the generator (legacy_ptr.py)

Scope: the four wall loops that fill the first off-wall advective flux through a scratch pointer: `SPECIES_ADVECTION_PART_1_NEW` (divg.f90 at 36975d7:
1021-1084 `WALL_LOOP_2`, species N; 1120-1183 `WALL_LOOP_3`, density slot 0) and `MASS_FINITE_DIFFERENCES` (mass.f90:93-188 `WALL_LOOP_2`, species N, which also writes the wall-face
flux with the thin-obstruction test; 224-320 `WALL_LOOP_3`, slot 0, wall-face value `B1%RHO_F/MW_F` through `GET_MOLECULAR_WEIGHT`). The mass.f90 loops also need the wall tables of
`massflux-wall-table-spec.md` (`B1_RHO_F`, `B1_ZZ_F`, `WALL_INDEX`, `SOLID`, `EXTERIOR`). Code: `amrex/s4_mass/s5_gen/legacy_ptr.py`, enabled per kernel by `scratch = true` in the sidecar.

## 1. What upstream does

Per thread, three pointers are aimed at fixed-size private arrays before the loop (divg.f90:1016-1019: `U_TEMP=>U_WORK`, `F_TEMP=>F_WORK`, `Z_TEMP=>Z_WORK`, with
`U_WORK,Z_WORK(-1:3,-1:3,-1:3)` and `F_WORK(0:3,0:3,0:3)`, divg.f90:980-982). Inside the loop a one-dimensional strip of the scalar is copied into
`Z_TEMP` by an array-constructor assignment (divg.f90:1037, 1047, ...), the face velocity into `U_TEMP(1,1,1)`, the routine
`GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,IOR,LIMITER)` (func.f90:1330-1487) is called for one face, and `F_TEMP(1,1,1)` is stored. The
guard reads `CELL(CELL_INDEX(i,j,k))%WALL_INDEX(n)`.

## 2. Rewrites (done on the parsed loop before the normal analysis; the upstream text is not edited)

| upstream | generated |
|---|---|
| `U_TEMP(1,1,1)` | `US` (private REAL scalar) |
| `F_TEMP(1,1,1)` | `FS` (private REAL scalar) |
| `Z_TEMP(lo:hi,1,1)` / `(1,lo:hi,1)` / `(1,1,lo:hi) = (/.../)` | `ZS(lo:hi) = (/.../)`, expanded to scalar assignments by `Emitter.assign_constructor`; `ZS` is private `ZS(0:3)` |
| 3-wide constructor (`0:2` or `1:3`) | padded with `0._EB` to `0:3` (trailing for `+IOR`, leading for `-IOR`) |
| `CALL GET_SCALAR_FACE_VALUE(U_TEMP,Z_TEMP,F_TEMP,1,1,1,1,1,1,IOR,LIM)` | `CALL GET_SCALAR_FACE_VALUE_PT(US,ZS,FS,LIM)` |
| local `PARAMETER` with a numeric literal (`DUMMY=0._EB`, mass.f90:30) | the literal (the generic path has no rule for a PARAMETER declared inside the routine and prints the bare name, which does not compile) |
| `CELL(CELL_INDEX(a,b,c))%WALL_INDEX(n)` | `WALL_INDEX(a,b,c,n)`, read-only INTEGER table `(0:IBAR+1,0:JBAR+1,0:KBAR+1,-3:3)` (policy.arrays entry; type.f90:2182) |

Checks (each is a generator error with its own message, tested in `test/test_legacy_scratch.py`): the scratch pointers must be aimed at local TARGET arrays with
literal bounds by an association above the range (an association inside the range, a mesh array as target, a rank mismatch are refused); the three call
arguments must be three distinct scratch pointers in the roles A, Z, F; arguments 4-9 must be the literal 1 and IOR a literal 1, 2 or 3; the constructor before
a call must fill the axis `|IOR|`; only the element (1,1,1) of `U_TEMP`/`F_TEMP` may be used; every other use of a scratch pointer is refused; `n` of
`WALL_INDEX(n)` must be a literal in -3..-1, 1..3 and the table is never written.

Why the pad: the caller fills three of the four elements and the fourth keeps the value an earlier wall left. Only MP5 reads it. The pad makes the generated code
deterministic and equals the `DUMMY=0` of mass.f90 and upstream patches UP-0001/UP-0002. For the limiters 0-4 it changes nothing (tested bitwise against the
upstream routine).

## 3. Callee interface assumed

`SUBROUTINE GET_SCALAR_FACE_VALUE_PT(A,Z,F,LIMITER)`: `A` real in (face velocity), `Z(0:3)` real in (four cells along the normal, face between `Z(1)` and
`Z(2)`), `F` real out, `LIMITER` integer by value, `!$omp declare target`. The sidecar keys `scratch_callee` and `scratch_call` rename the callee and the upstream
routine. The tests use a stand-in with this interface (limiters 0-4) and also the real callee derived by `s5_gsfv.py` (limiters 0-5; section 5). Two facts from the tests:

* the generic call analysis (`leaf_events`) treats a dummy without INTENT as in/out, so a `VALUE` dummy (`LIMITER`) made `I_FLUX_LIMITER` look like a loop-carried
  scalar. `legacy_ptr` sets intent IN on VALUE dummies of this callee only; the general fix is one condition in `leaf_events` (a VALUE dummy is an input) and
  belongs to the owner of that function.
* wall regions use the directive `S5_LOOP_WALL` (always `distribute parallel do`), so `S5_CALLEE_DPD` and `S5_CALLEE_BIND` give identical code for these kernels.

## 4. Plane aliases (the FX_P remap)

`FX_P(LBOUND(FX_ZZ,1):,LBOUND(FX_ZZ,2):,LBOUND(FX_ZZ,3):) => FX_ZZ(:,:,:,N)` (divg.f90:1009-1011, 1108-1110; mass.f90:81-83, 212-214) keeps the index values:
`FX_P(i,j,k)` is `FX_ZZ(i,j,k,N)`. `leg_plane_alias` resolves a pointer at a CALL to `(array, index)` from the nearest remap before the call in the same
statement list, and refuses anything but this exact form (same array in the three bounds and in the target, three bare `:`, a scalar literal or name as the last
index). `leg_rewrite_planes` replaces each use of `FX_P` as a CALL argument by the section `FX_ZZ(:, :, :, N)`, removes the remap statements, and refuses any other use
(element reference, call without remap, unused remap). Tested on the real units (divg.f90 `SPECIES_ADVECTION_PART_1_NEW`, mass.f90 `MASS_FINITE_DIFFERENCES`) and on
refusal cases in `test/test_legacy_plane.py`. The result serves any whole-field callee: an explicit-shape dummy `F(lb1:ub1,lb2:ub2,lb3:ub3)` receives the contiguous
plane by sequence association; a per-face form writes `FX_ZZ(I,J,K,N) = FS` itself.

## 5. Interaction with the Mesh Data hooks (`s5_gsfv.py`)

`s5_gsfv.install(s5gen)` applies four hooks (`Gen.callee`, `Emitter.call`, `Emitter.part_ref`, `Emitter.pointer_stmt`) that do the same caller-side rewrite at emission time. Both
mechanisms were run together (`test/test_legacy_scratch.py`, second half: generator started through `test/leg_gen_hooks.py`, real source, no stand-in):

* For a kernel with `scratch = true` the scratch pointers and the call are rewritten before emission, so `part_ref`, `call` and `pointer_stmt` find nothing to do. Only the
  `Gen.callee` hook is used: the call is named `GET_SCALAR_FACE_VALUE_PT`, which does not exist in the source tree; `legacy_ptr` asks `Gen.callee` for the upstream name
  `GET_SCALAR_FACE_VALUE`, takes the record whose `emit_name` is the one-face name, and registers it under that name with a unit and scope parsed from the derived text (the
  upstream unit has 11 dummies in another order, so the intents of a 4-argument call must come from the derived routine). If the hook patch changes the record's keys (`emit_name`,
  `text`, `file`) the registration must follow.
* With the real callee (limiters 0-5, MP5 included) the four nests match the upstream loop bitwise. For the divg.f90 nests the reference for MP5 is the loop with the pad of upstream
  patch UP-0001 (the verbatim loop reads a stale element there); for the mass.f90 nests the verbatim loop.
* Without `scratch = true` the hook-only path leaves `US`, `FS`, `ZS` without a declaration: `Gen.gsfv_privates` is set by `emit_call` but nothing in `s5gen.py` reads it. `legacy_ptr`
  declares them (synthetic private locals), so the combination works; the hook patch needs an equivalent when the hooks are used without `legacy_ptr`.
* Review points for the hook patch: (a) one mechanism must own the pad of the unassigned element (both pad to `0._EB`; a second pad is a no-op); (b) `ZS`/`US`/`FS` must not be renamed;
  (c) the intent of the callee's `LIMITER` (the derived routine declares it `INTENT(IN)`; the announced interface says by value, which needs intent IN in the generic call analysis);
  (d) the six per-limiter whole-field kernels take the plane from `leg_rewrite_planes` or an equivalent.

## 6. Status

* The four wall nests generate and match the upstream loop text bitwise (reference = verbatim upstream loop text calling the upstream pointer-based routine; limiters 0-4 with the stand-in, 0-5 with the derived callee,
  species 1-3, three data sets, six flag sets, both callee switches, 4 and 8 threads, mutants, negative cases): `test/test_legacy_scratch.py`. Test data keep the faces written by
  different walls disjoint (the upstream `!$OMP DO` over walls has no ordering either; two walls writing one face would make the reference itself order dependent).
* The kernel entries are NOT in the shared `s5_markers.toml`: with the real source the callee `GET_SCALAR_FACE_VALUE_PT` does not exist yet, so adding them would stop the whole
  generation. They are added with the callee (the entries are in `test/test_legacy_scratch.py`, function `markers`, plus `policy.arrays.WALL_INDEX` and `policy.assoc` for UU, VV, WW).
* MP5 is covered only in the run with the derived callee.
* The two nests of divg.f90 and mass.f90 whose range holds a loop over species (`SPECIES_LOOP`) are kernels per species N (N is an argument), as the other species kernels.
