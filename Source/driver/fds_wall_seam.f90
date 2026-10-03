!> \brief Host side of the wall-table seam of ADR-001 (rulings W1 and W2): per-box wall index lists, the staged copy of the four host-produced wall arrays that a
!> device kernel would read, and the checksum assertion against the host arrays. notes/wall-seam-design.md has the design and what still needs the GPU.
!>
!> W1: one set of wall tables per mesh object in upstream order with the global IW (the driver has one mesh object per box, so the lists are the identity ranges):
!> WLIST_EXT = 1..N_EXTERNAL_WALL_CELLS and WLIST_INT = N_EXTERNAL_WALL_CELLS+1..N_WALL_CELLS. The tables are rebuilt by FDS_WSEAM_REFRESH when the wall counts change
!> (REASSIGN_WALL_CELLS, reallocation of WALL); the function reports whether a rebuild happened, so that a device shim knows to re-gather its tables.
!> W2: UVW_SAVE, U_GHOST, V_GHOST, W_GHOST (producers MATCH_VELOCITY and ccib.f90) stay on the host. FDS_WSEAM_UPLOAD gathers them over WLIST_EXT into STAGE(4,:) (the
!> stand-in for the host-to-device copy of the shim) and asserts that the staged copy and the host arrays have the same checksum; FDS_WSEAM_CHECK asserts that the host arrays
!> still match the staged copy (nothing wrote them behind the upload). The "device" side of the copy is host memory here: the real transfer is GPU work.
MODULE FDS_WALL_SEAM

USE ISO_C_BINDING
USE PRECISION_PARAMETERS, ONLY: EB
USE FDS_BOX_OBJ, ONLY: BOX_OBJ
USE MESH_VARIABLES, ONLY: MESH_TYPE

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE

TYPE WSEAM_T
   INTEGER :: N_WALL=-1,NWL_EXT=0,NWL_INT=0
   INTEGER, ALLOCATABLE :: WLIST_EXT(:),WLIST_INT(:)
   REAL(EB), ALLOCATABLE :: STAGE(:,:)       !< (4,NWL_EXT): UVW_SAVE, U_GHOST, V_GHOST, W_GHOST of the external walls
   INTEGER(C_INT64_T) :: CHK=0               !< checksum of the staged copy
   LOGICAL :: STAGED=.FALSE.
END TYPE WSEAM_T

TYPE(WSEAM_T), ALLOCATABLE, TARGET, SAVE :: WS(:)
INTEGER, SAVE :: NWS=0

CONTAINS

SUBROUTINE ENSURE(NM)
INTEGER, INTENT(IN) :: NM
TYPE(WSEAM_T), ALLOCATABLE :: TMP(:)
INTEGER :: I
IF (NM>NWS) THEN
   ALLOCATE(TMP(NM))
   DO I=1,NWS
      CALL MOVE_ALLOC(WS(I)%WLIST_EXT,TMP(I)%WLIST_EXT)
      CALL MOVE_ALLOC(WS(I)%WLIST_INT,TMP(I)%WLIST_INT)
      CALL MOVE_ALLOC(WS(I)%STAGE,TMP(I)%STAGE)
      TMP(I)%N_WALL=WS(I)%N_WALL ; TMP(I)%NWL_EXT=WS(I)%NWL_EXT ; TMP(I)%NWL_INT=WS(I)%NWL_INT
      TMP(I)%CHK=WS(I)%CHK ; TMP(I)%STAGED=WS(I)%STAGED
   ENDDO
   CALL MOVE_ALLOC(TMP,WS)
   NWS=NM
ENDIF
END SUBROUTINE ENSURE

!> Order-sensitive bitwise checksum of the four arrays over the list.
FUNCTION CHECKSUM(M,LIST,N) RESULT(C)
TYPE(MESH_TYPE), INTENT(IN) :: M
INTEGER, INTENT(IN) :: N,LIST(:)
INTEGER(C_INT64_T) :: C
INTEGER :: I,IW
C=0_C_INT64_T
DO I=1,N
   IW=LIST(I)
   C=IEOR(ISHFTC(C,7),TRANSFER(M%UVW_SAVE(IW),C))
   C=IEOR(ISHFTC(C,7),TRANSFER(M%U_GHOST(IW),C))
   C=IEOR(ISHFTC(C,7),TRANSFER(M%V_GHOST(IW),C))
   C=IEOR(ISHFTC(C,7),TRANSFER(M%W_GHOST(IW),C))
ENDDO
END FUNCTION CHECKSUM

FUNCTION CHECKSUM_STAGE(S) RESULT(C)
TYPE(WSEAM_T), INTENT(IN) :: S
INTEGER(C_INT64_T) :: C
INTEGER :: I,K
C=0_C_INT64_T
DO I=1,S%NWL_EXT
   DO K=1,4
      C=IEOR(ISHFTC(C,7),TRANSFER(S%STAGE(K,I),C))
   ENDDO
ENDDO
END FUNCTION CHECKSUM_STAGE

!> (Re)build the index lists of box NM. Returns 1 when the lists were (re)built (first call or changed wall counts), 0 when they were current.
FUNCTION FDS_WSEAM_REFRESH(NM) BIND(C,NAME='fds_wseam_refresh') RESULT(CHANGED)
INTEGER(C_INT), VALUE :: NM
INTEGER(C_INT) :: CHANGED
TYPE(MESH_TYPE), POINTER :: M
INTEGER :: I
CALL ENSURE(NM)
M => BOX_OBJ(NM)
CHANGED=0
IF (WS(NM)%N_WALL==M%N_WALL_CELLS .AND. WS(NM)%NWL_EXT==M%N_EXTERNAL_WALL_CELLS) RETURN
CHANGED=1
WS(NM)%N_WALL=M%N_WALL_CELLS
WS(NM)%NWL_EXT=M%N_EXTERNAL_WALL_CELLS
WS(NM)%NWL_INT=M%N_WALL_CELLS-M%N_EXTERNAL_WALL_CELLS
IF (ALLOCATED(WS(NM)%WLIST_EXT)) DEALLOCATE(WS(NM)%WLIST_EXT)
IF (ALLOCATED(WS(NM)%WLIST_INT)) DEALLOCATE(WS(NM)%WLIST_INT)
IF (ALLOCATED(WS(NM)%STAGE))     DEALLOCATE(WS(NM)%STAGE)
ALLOCATE(WS(NM)%WLIST_EXT(MAX(1,WS(NM)%NWL_EXT)),WS(NM)%WLIST_INT(MAX(1,WS(NM)%NWL_INT)),WS(NM)%STAGE(4,MAX(1,WS(NM)%NWL_EXT)))
DO I=1,WS(NM)%NWL_EXT ; WS(NM)%WLIST_EXT(I)=I ; ENDDO
DO I=1,WS(NM)%NWL_INT ; WS(NM)%WLIST_INT(I)=WS(NM)%NWL_EXT+I ; ENDDO
WS(NM)%STAGE=0._EB ; WS(NM)%STAGED=.FALSE.
END FUNCTION FDS_WSEAM_REFRESH

!> Counts of box NM (builds the lists when needed): NEXT = NWL_EXT, NINT = NWL_INT.
SUBROUTINE FDS_WSEAM_COUNTS(NM,NEXT,NINT) BIND(C,NAME='fds_wseam_counts')
INTEGER(C_INT), VALUE :: NM
INTEGER(C_INT), INTENT(OUT) :: NEXT,NINT
INTEGER(C_INT) :: DUMMY
DUMMY=FDS_WSEAM_REFRESH(NM)
NEXT=WS(NM)%NWL_EXT ; NINT=WS(NM)%NWL_INT
END SUBROUTINE FDS_WSEAM_COUNTS

!> W2 upload of box NM: gather, then assert staged checksum = host checksum. Returns 0 on success, 1 when the assertion fails.
FUNCTION FDS_WSEAM_UPLOAD(NM) BIND(C,NAME='fds_wseam_upload') RESULT(RC)
INTEGER(C_INT), VALUE :: NM
INTEGER(C_INT) :: RC
TYPE(MESH_TYPE), POINTER :: M
INTEGER :: I,IW,DUMMY
INTEGER(C_INT64_T) :: CH
DUMMY=FDS_WSEAM_REFRESH(NM)
M => BOX_OBJ(NM)
RC=0
DO I=1,WS(NM)%NWL_EXT
   IW=WS(NM)%WLIST_EXT(I)
   WS(NM)%STAGE(1,I)=M%UVW_SAVE(IW) ; WS(NM)%STAGE(2,I)=M%U_GHOST(IW)
   WS(NM)%STAGE(3,I)=M%V_GHOST(IW)  ; WS(NM)%STAGE(4,I)=M%W_GHOST(IW)
ENDDO
CH=CHECKSUM(M,WS(NM)%WLIST_EXT,WS(NM)%NWL_EXT)
WS(NM)%CHK=CHECKSUM_STAGE(WS(NM))
WS(NM)%STAGED=.TRUE.
IF (CH/=WS(NM)%CHK) RC=1
END FUNCTION FDS_WSEAM_UPLOAD

!> Assert that the host arrays of box NM still match the staged copy. Returns 0 when they do (or nothing is staged), 1 otherwise.
FUNCTION FDS_WSEAM_CHECK(NM) BIND(C,NAME='fds_wseam_check') RESULT(RC)
INTEGER(C_INT), VALUE :: NM
INTEGER(C_INT) :: RC
TYPE(MESH_TYPE), POINTER :: M
RC=0
IF (NM>NWS) RETURN
IF (.NOT. WS(NM)%STAGED) RETURN
M => BOX_OBJ(NM)
IF (CHECKSUM(M,WS(NM)%WLIST_EXT,WS(NM)%NWL_EXT)/=WS(NM)%CHK) RC=1
END FUNCTION FDS_WSEAM_CHECK

!> Test hook (negative control): perturb the first staged-host value of box NM after an upload, as a host producer that wrote late would. Returns 0, or 2 when the box has no external wall.
FUNCTION FDS_WSEAM_POKE(NM) BIND(C,NAME='fds_wseam_poke') RESULT(RC)
INTEGER(C_INT), VALUE :: NM
INTEGER(C_INT) :: RC
TYPE(MESH_TYPE), POINTER :: M
CALL ENSURE(NM)
RC=2
IF (WS(NM)%NWL_EXT<1) RETURN
M => BOX_OBJ(NM)
M%UVW_SAVE(WS(NM)%WLIST_EXT(1))=M%UVW_SAVE(WS(NM)%WLIST_EXT(1))+1.E-9_EB
RC=0
END FUNCTION FDS_WSEAM_POKE

END MODULE FDS_WALL_SEAM
