!> \brief Hook module of the guarded FDS edits (USE_AMREX=ON only; patches 0003, 0004, 0005 in Source/driver/patches).
!>
!> EXTERNAL_GHOSTS_FILLED: when .TRUE. the ghost cells across a box boundary (EXTERNAL_WALL(IW)%NOM>0) have been filled from AMReX data before the boundary
!> routines run, so the NOM>0 branches of VISCOSITY_BC, NO_FLUX, VELOCITY_BC, MATCH_VELOCITY, MATCH_VELOCITY_FLUX (velo.f90) and ASSIGN_GHOST_VALUE (wall.f90)
!> must not overwrite them with OMESH values ("filled externally", ADR-001 Option C). Default .FALSE.: FDS behaves as before. S4 keeps it .FALSE. (the driver
!> fills OMESH and runs the unmodified routines, GhostExchange.cpp); the TRUE path is exercised from S5 on, once the box views below are in use.
!>
!> BOX_VIEW(NM): pointer views of the AMReX data of box NM, with the FDS index bounds (lower bounds as in MESHES(NM)%X), so that POINT_TO_BOX (patch 0005,
!> mesh.f90) can point the MESH_POINTERS module pointers at box data without going through the MESHES(NM)%X allocatables. FDS_HOOK_SET_VIEW builds a view from
!> the C address of a FAB with the standard C_F_POINTER plus pointer bounds remapping (Fortran 2008, no compiler extension).
!>
!> Kernel-facing file rules (M2a):
!> (a) Passive scalars (N_TOTAL_SCALARS beyond the tracked species) are handled by Fields.cpp: ZZ/ZZS views carry the last extent N_TOTAL_SCALARS as given by the caller.
!> (b) Only uniform Cartesian metrics are used; R(I)/RRN(I) are not involved and CYLINDRICAL/TRN* meshes are rejected at level-0 assembly (IR-002).
MODULE FDS_AMREX_HOOKS

USE ISO_C_BINDING
USE PRECISION_PARAMETERS, ONLY: EB

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE

LOGICAL, PUBLIC, SAVE :: EXTERNAL_GHOSTS_FILLED = .FALSE.

!> Time, time step and cycle number that the driver hands to FDS_SETUP(MODE=3) (end-of-step outputs, patch 0006).
REAL(EB), PUBLIC, SAVE :: OUT_T = 0._EB, OUT_DT = 0._EB
INTEGER,  PUBLIC, SAVE :: OUT_ICYC = 0
!> Set .TRUE. by the patched FDS_SETUP (patch 0006) when set-up returns; the driver calls MODE=3 only if it is .TRUE. (an unpatched main.f90 would run END_FDS).
LOGICAL,  PUBLIC, SAVE :: STEP_OUTPUTS_PATCHED = .FALSE.
!> Patch 0011: when .TRUE. the patched FDS_SETUP(MODE=3) also runs the per-mesh dump loop of MAIN_LOOP (DUMP_MESH_OUTPUTS: slices, boundary files, Plot3D, ...). Default .FALSE.
!> (set by fds_hook_set_mesh_dumps), so the behaviour of patch 0006 is unchanged until the driver asks.
LOGICAL,  PUBLIC, SAVE :: OUT_MESH_DUMPS = .FALSE.
!> Patch 0012: number of output meshes (the meshes of the FDS output files, which may exceed NMESHES once boxes of finer levels are output meshes). 0 (default) = NMESHES.
!> Must be set (fds_hook_set_out_meshes) before fds_setup(0), because ASSIGN_FILE_NAMES sizes the output tables.
INTEGER,  PUBLIC, SAVE :: N_OUT_MESHES = 0

TYPE, PUBLIC :: BOX_VIEW_TYPE
   REAL(EB), POINTER, DIMENSION(:,:,:)   :: U=>NULL(),V=>NULL(),W=>NULL(),US=>NULL(),VS=>NULL(),WS=>NULL(),D=>NULL(),DS=>NULL(),H=>NULL(),HS=>NULL(), &
                                            KRES=>NULL(),FVX=>NULL(),FVY=>NULL(),FVZ=>NULL(),RHO=>NULL(),RHOS=>NULL(),MU=>NULL(),TMP=>NULL(),Q=>NULL(),RSUM=>NULL()
   REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZ=>NULL(),ZZS=>NULL()
END TYPE BOX_VIEW_TYPE

TYPE(BOX_VIEW_TYPE), ALLOCATABLE, TARGET, PUBLIC, SAVE :: BOX_VIEW(:)

PUBLIC :: FDS_HOOK_SET_FLAG,FDS_HOOK_SET_VIEW,FDS_HOOK_SET_STEP,FDS_HOOK_SET_MESH_DUMPS,FDS_HOOK_SET_OUT_MESHES,FDS_HOOK_STEP_OUTPUTS,FDS_HOOK_FINE_ABORT,FDS_HOOK_FINE_GUARD,FDS_HOOK_SET_FINE_READY,FDS_HOOK_SHADOW,SHADOW_PROC,SHADOW_IF,FDS_HOOK_BIND_VIEW,FDS_HOOK_L0_ONLY

!> D-056 (option B): the set of kernel wrappers (fds_kernels.f90 entry names) that may run on a fine-level mesh number (NM > NMESHES). Empty until patch 0007 is validated and the
!> kernels are switched to POINT_TO_BOX: every wrapper aborts on a fine number now. Set by FDS_HOOK_FINE_READY (also callable from the draft driver file fds_fine_mesh_b.f90).
LOGICAL, SAVE :: FINE_READY = .FALSE.

!> Shadow hook (fine-level box test, D-056 option B): when a procedure is registered (draft file fds_fine_box_b.f90, option FDS_AMR_FINE_B_DRAFT), every kernel wrapper calls it before (PHASE=0)
!> and after (PHASE=1) the kernel on a level-0 box, so that the same kernel can be run on a fine-level copy through POINT_TO_BOX and compared. Not registered in a normal run.
ABSTRACT INTERFACE
   SUBROUTINE SHADOW_IF(PHASE,KIND,NM,T,DT,EST,DTNEW,ICHG)
      IMPORT :: EB
      INTEGER, INTENT(IN) :: PHASE,KIND,NM,EST,ICHG
      REAL(EB), INTENT(IN) :: T,DT,DTNEW
   END SUBROUTINE SHADOW_IF
END INTERFACE
PROCEDURE(SHADOW_IF), POINTER, SAVE :: SHADOW_PROC => NULL()

CONTAINS

!> Set or clear EXTERNAL_GHOSTS_FILLED.
SUBROUTINE FDS_HOOK_SET_FLAG(FLAG) BIND(C,NAME='fds_hook_set_flag')
INTEGER(C_INT), VALUE :: FLAG
EXTERNAL_GHOSTS_FILLED = (FLAG/=0)
END SUBROUTINE FDS_HOOK_SET_FLAG

!> 1 when main.f90 carries patch 0006 (FDS_SETUP(MODE=3) end-of-step outputs), else 0.
FUNCTION FDS_HOOK_STEP_OUTPUTS() BIND(C,NAME='fds_hook_step_outputs') RESULT(F)
INTEGER(C_INT) :: F
F = MERGE(1,0,STEP_OUTPUTS_PATCHED)
END FUNCTION FDS_HOOK_STEP_OUTPUTS

!> Set T, DT and ICYC for the next FDS_SETUP(MODE=3) call.
SUBROUTINE FDS_HOOK_SET_STEP(T,DT,ICYC) BIND(C,NAME='fds_hook_set_step')
REAL(C_DOUBLE), VALUE :: T,DT
INTEGER(C_INT), VALUE :: ICYC
OUT_T = T ; OUT_DT = DT ; OUT_ICYC = ICYC
END SUBROUTINE FDS_HOOK_SET_STEP

!> Switch the per-mesh dump loop of FDS_SETUP(MODE=3) on (FLAG/=0) or off (patch 0011).
SUBROUTINE FDS_HOOK_SET_MESH_DUMPS(FLAG) BIND(C,NAME='fds_hook_set_mesh_dumps')
INTEGER(C_INT), VALUE :: FLAG
OUT_MESH_DUMPS = (FLAG/=0)
END SUBROUTINE FDS_HOOK_SET_MESH_DUMPS

!> Set the number of output meshes N_OUT_MESHES (patch 0012); effective only if called before fds_setup(0).
SUBROUTINE FDS_HOOK_SET_OUT_MESHES(N) BIND(C,NAME='fds_hook_set_out_meshes')
INTEGER(C_INT), VALUE :: N
N_OUT_MESHES = N
END SUBROUTINE FDS_HOOK_SET_OUT_MESHES

!> Make BOX_VIEW(NM)%<WHICH> a view of the C array P: rank 3 (LB(1:3), EXT(1:3)) or, for ZZ (WHICH=21) and ZZS (WHICH=22), rank 4 (LB(1:4), EXT(1:4)).
!> WHICH: 1 U 2 V 3 W 4 US 5 VS 6 WS 7 D 8 DS 9 H 10 HS 11 KRES 12 FVX 13 FVY 14 FVZ 15 RHO 16 RHOS 17 MU 18 TMP 19 Q 20 RSUM. NM: 1-based box number;
!> BOX_VIEW is allocated on first use with NMAX entries. Returns 0, or 1 for an unknown WHICH or NM out of range.
FUNCTION FDS_HOOK_SET_VIEW(NMAX,NM,WHICH,LB,EXT,P) BIND(C,NAME='fds_hook_set_view') RESULT(IERR)
INTEGER(C_INT), VALUE :: NMAX,NM,WHICH
INTEGER(C_INT), INTENT(IN) :: LB(4),EXT(4)
TYPE(C_PTR), VALUE :: P
INTEGER(C_INT) :: IERR
IERR = 1
IF (.NOT.ALLOCATED(BOX_VIEW)) ALLOCATE(BOX_VIEW(NMAX))
IF (NM<1 .OR. NM>SIZE(BOX_VIEW)) RETURN
IERR = FDS_HOOK_BIND_VIEW(BOX_VIEW(NM),WHICH,LB,EXT,P)
END FUNCTION FDS_HOOK_SET_VIEW

!> The body of FDS_HOOK_SET_VIEW for any view object V: BOX_VIEW(NM) of a level-0 box, FINE_LEVEL(L)%VIEW(IB) of a fine-level box (draft/fds_fine_box_b.f90). 0, or 1 for an unknown WHICH.
FUNCTION FDS_HOOK_BIND_VIEW(V,WHICH,LB,EXT,P) RESULT(IERR)
TYPE(BOX_VIEW_TYPE), INTENT(INOUT) :: V
INTEGER(C_INT), VALUE :: WHICH
INTEGER(C_INT), INTENT(IN) :: LB(4),EXT(4)
TYPE(C_PTR), VALUE :: P
INTEGER(C_INT) :: IERR
REAL(EB), POINTER, DIMENSION(:) :: FLAT
IERR = 1
IF (WHICH<1 .OR. WHICH>22) RETURN
IF (WHICH>=21) THEN
   CALL C_F_POINTER(P,FLAT,[INT(EXT(1))*EXT(2)*EXT(3)*EXT(4)])
   IF (WHICH==21) V%ZZ (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1,LB(4):LB(4)+EXT(4)-1) => FLAT
   IF (WHICH==22) V%ZZS(LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1,LB(4):LB(4)+EXT(4)-1) => FLAT
   IERR = 0
   RETURN
ENDIF
CALL C_F_POINTER(P,FLAT,[INT(EXT(1))*EXT(2)*EXT(3)])
SELECT CASE(WHICH)
   CASE( 1) ; V%U   (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 2) ; V%V   (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 3) ; V%W   (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 4) ; V%US  (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 5) ; V%VS  (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 6) ; V%WS  (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 7) ; V%D   (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 8) ; V%DS  (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 9) ; V%H   (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(10) ; V%HS  (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(11) ; V%KRES(LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(12) ; V%FVX (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(13) ; V%FVY (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(14) ; V%FVZ (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(15) ; V%RHO (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(16) ; V%RHOS(LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(17) ; V%MU  (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(18) ; V%TMP (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(19) ; V%Q   (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(20) ; V%RSUM(LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
END SELECT
IERR = 0
END FUNCTION FDS_HOOK_BIND_VIEW

!> Declare (FLAG/=0) that the kernel wrappers may run on fine-level mesh numbers (only after patch 0007 is validated and the kernels call POINT_TO_BOX; off by default).
SUBROUTINE FDS_HOOK_SET_FINE_READY(FLAG) BIND(C,NAME='fds_hook_set_fine_ready')

INTEGER(C_INT), VALUE :: FLAG

FINE_READY = (FLAG/=0)

END SUBROUTINE FDS_HOOK_SET_FINE_READY

!> Stop the run with a clear message: mesh number NM is not a level-0 FDS mesh and WHERE would read MESHES(NM) (D-056 option B).
SUBROUTINE FDS_HOOK_FINE_ABORT(WHERE,NM)

CHARACTER(*), INTENT(IN) :: WHERE
INTEGER, INTENT(IN) :: NM

WRITE(0,'(A,A,A,I0,A)') 'fds_amr ERROR in ',TRIM(WHERE),': mesh number ',NM,' is not a level-0 FDS mesh. Boxes of refinement level > 0 are fine-level mesh objects (D-056 option B, ' // &
   'FINE_LEVEL(:), patch 0007) and must be reached through POINT_TO_BOX; this routine reads MESHES(NM) and has not been made fine-ready.'
ERROR STOP 1

END SUBROUTINE FDS_HOOK_FINE_ABORT

!> Guard at the entry of every kernel wrapper: a fine-level mesh number (NM > NMESHES) aborts with a clear message unless fine boxes have been declared ready (FINE_READY).
SUBROUTINE FDS_HOOK_FINE_GUARD(WHERE,NM,NMESHES_L0)

CHARACTER(*), INTENT(IN) :: WHERE
INTEGER(C_INT), INTENT(IN) :: NM
INTEGER, INTENT(IN) :: NMESHES_L0

IF (NM>NMESHES_L0 .AND. .NOT.FINE_READY) CALL FDS_HOOK_FINE_ABORT(WHERE,INT(NM))

END SUBROUTINE FDS_HOOK_FINE_GUARD

!> Guard at the entry of a driver routine that is level-0 only (it reads MESHES(NM), FDS pressure/vent/wall structures of a level-0 mesh): a number above NMESHES always aborts, fine-ready or not.
SUBROUTINE FDS_HOOK_L0_ONLY(WHERE,NM,NMESHES_L0)

CHARACTER(*), INTENT(IN) :: WHERE
INTEGER(C_INT), INTENT(IN) :: NM
INTEGER, INTENT(IN) :: NMESHES_L0

IF (NM>NMESHES_L0) CALL FDS_HOOK_FINE_ABORT(WHERE,INT(NM))

END SUBROUTINE FDS_HOOK_L0_ONLY

!> Call of the registered shadow procedure (no-op when none is registered).
SUBROUTINE FDS_HOOK_SHADOW(PHASE,KIND,NM,T,DT,EST,DTNEW,ICHG)

INTEGER, INTENT(IN) :: PHASE,KIND,NM,EST,ICHG
REAL(EB), INTENT(IN) :: T,DT,DTNEW

IF (ASSOCIATED(SHADOW_PROC)) CALL SHADOW_PROC(PHASE,KIND,NM,T,DT,EST,DTNEW,ICHG)

END SUBROUTINE FDS_HOOK_SHADOW

END MODULE FDS_AMREX_HOOKS
