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

TYPE, PUBLIC :: BOX_VIEW_TYPE
   REAL(EB), POINTER, DIMENSION(:,:,:)   :: U=>NULL(),V=>NULL(),W=>NULL(),US=>NULL(),VS=>NULL(),WS=>NULL(),D=>NULL(),DS=>NULL(),H=>NULL(),HS=>NULL(), &
                                            KRES=>NULL(),FVX=>NULL(),FVY=>NULL(),FVZ=>NULL(),RHO=>NULL(),RHOS=>NULL(),MU=>NULL(),TMP=>NULL(),Q=>NULL(),RSUM=>NULL()
   REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZ=>NULL(),ZZS=>NULL()
END TYPE BOX_VIEW_TYPE

TYPE(BOX_VIEW_TYPE), ALLOCATABLE, TARGET, PUBLIC, SAVE :: BOX_VIEW(:)

PUBLIC :: FDS_HOOK_SET_FLAG,FDS_HOOK_SET_VIEW,FDS_HOOK_SET_STEP,FDS_HOOK_STEP_OUTPUTS

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

!> Make BOX_VIEW(NM)%<WHICH> a view of the C array P: rank 3 (LB(1:3), EXT(1:3)) or, for ZZ (WHICH=21) and ZZS (WHICH=22), rank 4 (LB(1:4), EXT(1:4)).
!> WHICH: 1 U 2 V 3 W 4 US 5 VS 6 WS 7 D 8 DS 9 H 10 HS 11 KRES 12 FVX 13 FVY 14 FVZ 15 RHO 16 RHOS 17 MU 18 TMP 19 Q 20 RSUM. NM: 1-based box number;
!> BOX_VIEW is allocated on first use with NMAX entries. Returns 0, or 1 for an unknown WHICH or NM out of range.
FUNCTION FDS_HOOK_SET_VIEW(NMAX,NM,WHICH,LB,EXT,P) BIND(C,NAME='fds_hook_set_view') RESULT(IERR)
INTEGER(C_INT), VALUE :: NMAX,NM,WHICH
INTEGER(C_INT), INTENT(IN) :: LB(4),EXT(4)
TYPE(C_PTR), VALUE :: P
INTEGER(C_INT) :: IERR
REAL(EB), POINTER, DIMENSION(:) :: FLAT
IERR = 1
IF (.NOT.ALLOCATED(BOX_VIEW)) ALLOCATE(BOX_VIEW(NMAX))
IF (NM<1 .OR. NM>SIZE(BOX_VIEW)) RETURN
IF (WHICH<1 .OR. WHICH>22) RETURN
IF (WHICH>=21) THEN
   CALL C_F_POINTER(P,FLAT,[INT(EXT(1))*EXT(2)*EXT(3)*EXT(4)])
   IF (WHICH==21) BOX_VIEW(NM)%ZZ (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1,LB(4):LB(4)+EXT(4)-1) => FLAT
   IF (WHICH==22) BOX_VIEW(NM)%ZZS(LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1,LB(4):LB(4)+EXT(4)-1) => FLAT
   IERR = 0
   RETURN
ENDIF
CALL C_F_POINTER(P,FLAT,[INT(EXT(1))*EXT(2)*EXT(3)])
SELECT CASE(WHICH)
   CASE( 1) ; BOX_VIEW(NM)%U   (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 2) ; BOX_VIEW(NM)%V   (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 3) ; BOX_VIEW(NM)%W   (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 4) ; BOX_VIEW(NM)%US  (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 5) ; BOX_VIEW(NM)%VS  (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 6) ; BOX_VIEW(NM)%WS  (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 7) ; BOX_VIEW(NM)%D   (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 8) ; BOX_VIEW(NM)%DS  (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 9) ; BOX_VIEW(NM)%H   (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(10) ; BOX_VIEW(NM)%HS  (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(11) ; BOX_VIEW(NM)%KRES(LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(12) ; BOX_VIEW(NM)%FVX (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(13) ; BOX_VIEW(NM)%FVY (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(14) ; BOX_VIEW(NM)%FVZ (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(15) ; BOX_VIEW(NM)%RHO (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(16) ; BOX_VIEW(NM)%RHOS(LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(17) ; BOX_VIEW(NM)%MU  (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(18) ; BOX_VIEW(NM)%TMP (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(19) ; BOX_VIEW(NM)%Q   (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(20) ; BOX_VIEW(NM)%RSUM(LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
END SELECT
IERR = 0
END FUNCTION FDS_HOOK_SET_VIEW

END MODULE FDS_AMREX_HOOKS
