!> \brief Flux read-out and override hooks of the AMReX driver (notes/flux-hooks-design.md): the stage face fluxes of the species equations, ADV (the product face value x
!> face velocity that the density update differences) and DIF (RHO_D_DZDX/Y/Z as DIVERGENCE_PART_1 hands them to the divergence), are copied into arrays registered by the driver,
!> and listed face values can be replaced before they are used (the coarse face flux is overwritten by the fine flux, D-050).
!>
!> Indexing: box number NM (level-0 mesh number or fine-level number, as the kernels). Face I of a flux array is FDS face I (the HIGH face of cell I): x faces I=0..IBAR,
!> J=1..JBAR, K=1..KBAR; y faces J=0..JBAR; z faces K=0..KBAR; scalar index n=1..N_TOTAL_SCALARS. The registered views have exactly these index ranges (lower bounds (0,1,1,1), (1,0,1,1), (1,1,0,1)).
!> Nothing is read or written for a box that has no mode set: the hooks return at once, the kernels then behave as without them (the DIF call is a guarded FDS patch, 0009; the ADV call is in
!> the generated density copy, fds_density_split.f90).
!>
!> Modes (per box and kind): 0 off, 1 read-out only (ADV: the density update is skipped for the box, nothing else changes), 2 override (the list of the box is applied; an empty list changes nothing).
MODULE FDS_FLUX_HOOKS

USE ISO_C_BINDING
USE PRECISION_PARAMETERS, ONLY: EB
USE GLOBAL_CONSTANTS, ONLY: N_TOTAL_SCALARS
USE MESH_POINTERS, ONLY: IBAR,JBAR,KBAR
USE FDS_BOX_OBJ, ONLY: BOX_OBJ
USE MESH_VARIABLES, ONLY: MESH_TYPE

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE
PUBLIC :: FDS_HOOK_DIF_FLUX,FDS_HOOK_ADV,FDS_FLUX_ACTIVE,PRX,PRY,PRZ

INTEGER, PARAMETER :: KADV=0,KDIF=1

TYPE VIEW4
   REAL(EB), POINTER :: A(:,:,:,:) => NULL()
END TYPE VIEW4

TYPE OVR_T
   INTEGER :: N=0                           !< number of listed faces
   INTEGER, ALLOCATABLE :: DIR(:),IDX(:,:)   !< direction 1..3, FDS face index (I,J,K)
   REAL(EB), ALLOCATABLE :: VAL(:,:)         !< VAL(n,entry)
END TYPE OVR_T

TYPE FLUX_BOX_T
   TYPE(VIEW4) :: V(0:1,3)                  !< registered read-out arrays (kind, direction)
   INTEGER :: MODE(0:1)=0
   TYPE(OVR_T) :: OVR(0:1)
   REAL(EB), ALLOCATABLE :: PROD(:,:,:,:,:)  !< PROD(:,:,:,n,dir): ADV face products of the override loop (0:IBAR+1 ... per direction)
END TYPE FLUX_BOX_T

TYPE(FLUX_BOX_T), ALLOCATABLE, TARGET, SAVE :: FB(:)
REAL(EB), POINTER, SAVE :: PRX(:,:,:,:) => NULL(),PRY(:,:,:,:) => NULL(),PRZ(:,:,:,:) => NULL()   !< ADV face products for the override copy of the density loop (set by FDS_HOOK_ADV)

CONTAINS

!> Allocate the per-box tables for box numbers 1..NTOT (idempotent; entries are kept).
SUBROUTINE FDS_FLUX_RESERVE(NTOT) BIND(C,NAME='fds_flux_reserve')
INTEGER(C_INT), VALUE :: NTOT
TYPE(FLUX_BOX_T), ALLOCATABLE :: TMP(:)
INTEGER :: I
IF (.NOT.ALLOCATED(FB)) THEN ; ALLOCATE(FB(NTOT)) ; RETURN ; ENDIF
IF (SIZE(FB)>=NTOT) RETURN
ALLOCATE(TMP(NTOT))
DO I=1,SIZE(FB) ; TMP(I)%V = FB(I)%V ; TMP(I)%MODE = FB(I)%MODE ; TMP(I)%OVR = FB(I)%OVR ; ENDDO
CALL MOVE_ALLOC(TMP,FB)
END SUBROUTINE FDS_FLUX_RESERVE

!> Register the read-out array of (box NM, kind 0 ADV / 1 DIF, direction 0..2): LB(4) lower bounds, EXT(4) extents, P the C data (all components, Fortran order). Returns 0, or 1 for a bad argument.
FUNCTION FDS_FLUX_REGISTER(NM,KIND,DIR,LB,EXT,P) BIND(C,NAME='fds_flux_register') RESULT(IERR)
INTEGER(C_INT), VALUE :: NM,KIND,DIR
INTEGER(C_INT), INTENT(IN) :: LB(4),EXT(4)
TYPE(C_PTR), VALUE :: P
INTEGER(C_INT) :: IERR
REAL(EB), POINTER :: FLAT(:)
IERR = 1
IF (.NOT.ALLOCATED(FB)) RETURN
IF (NM<1 .OR. NM>SIZE(FB) .OR. KIND<0 .OR. KIND>1 .OR. DIR<0 .OR. DIR>2) RETURN
CALL C_F_POINTER(P,FLAT,[INT(EXT(1))*EXT(2)*EXT(3)*EXT(4)])
FB(NM)%V(KIND,DIR+1)%A(LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1,LB(4):LB(4)+EXT(4)-1) => FLAT
IERR = 0
END FUNCTION FDS_FLUX_REGISTER

!> Set the mode of (box NM, kind). Returns 0, or 1 for a bad argument.
FUNCTION FDS_FLUX_SET_MODE(NM,KIND,MODE) BIND(C,NAME='fds_flux_set_mode') RESULT(IERR)
INTEGER(C_INT), VALUE :: NM,KIND,MODE
INTEGER(C_INT) :: IERR
IERR = 1
IF (.NOT.ALLOCATED(FB)) RETURN
IF (NM<1 .OR. NM>SIZE(FB) .OR. KIND<0 .OR. KIND>1 .OR. MODE<0 .OR. MODE>2) RETURN
FB(NM)%MODE(KIND) = MODE
IERR = 0
END FUNCTION FDS_FLUX_SET_MODE

!> Replace the override list of (box NM, kind) by N entries: DIR(i) 0..2, IDX(3,i) the FDS face index, VAL(NSC,i) the values (NSC must be N_TOTAL_SCALARS). N = 0 clears the list.
!> Every entry is range-checked against the box (x: I 0..IBAR, J 1..JBAR, K 1..KBAR; y, z likewise): returns 0, 1 bad argument, 2 face outside the valid range (list not stored).
FUNCTION FDS_FLUX_SET_OVERRIDE(NM,KIND,N,DIR,IDX,NSC,VAL) BIND(C,NAME='fds_flux_set_override') RESULT(IERR)
INTEGER(C_INT), VALUE :: NM,KIND,N,NSC
INTEGER(C_INT), INTENT(IN) :: DIR(*),IDX(3,*)
REAL(C_DOUBLE), INTENT(IN) :: VAL(NSC,*)
INTEGER(C_INT) :: IERR
TYPE(MESH_TYPE), POINTER :: M
INTEGER :: E,D,LO(3),HI(3)
IERR = 1
IF (.NOT.ALLOCATED(FB)) RETURN
IF (NM<1 .OR. NM>SIZE(FB) .OR. KIND<0 .OR. KIND>1 .OR. N<0) RETURN
IF (N>0 .AND. NSC/=N_TOTAL_SCALARS) RETURN
M => BOX_OBJ(NM)
DO E=1,N
   D = DIR(E)+1
   IF (D<1 .OR. D>3) RETURN
   LO = [1,1,1] ; HI = [M%IBAR,M%JBAR,M%KBAR]
   LO(D) = 0
   IF (ANY(IDX(:,E)<LO) .OR. ANY(IDX(:,E)>HI)) THEN ; IERR = 2 ; RETURN ; ENDIF
ENDDO
IF (ALLOCATED(FB(NM)%OVR(KIND)%DIR)) DEALLOCATE(FB(NM)%OVR(KIND)%DIR,FB(NM)%OVR(KIND)%IDX,FB(NM)%OVR(KIND)%VAL)
FB(NM)%OVR(KIND)%N = N
IF (N>0) THEN
   ALLOCATE(FB(NM)%OVR(KIND)%DIR(N),FB(NM)%OVR(KIND)%IDX(3,N),FB(NM)%OVR(KIND)%VAL(NSC,N))
   FB(NM)%OVR(KIND)%DIR = DIR(1:N)+1 ; FB(NM)%OVR(KIND)%IDX = IDX(1:3,1:N) ; FB(NM)%OVR(KIND)%VAL = VAL(1:NSC,1:N)
ENDIF
IERR = 0
END FUNCTION FDS_FLUX_SET_OVERRIDE

!> True when the box has a mode set for the kind (the generated density copy asks this before it calls FDS_HOOK_ADV).
LOGICAL FUNCTION FDS_FLUX_ACTIVE(NM,KIND)
INTEGER, INTENT(IN) :: NM,KIND
FDS_FLUX_ACTIVE = .FALSE.
IF (.NOT.ALLOCATED(FB)) RETURN
IF (NM<1 .OR. NM>SIZE(FB)) RETURN
FDS_FLUX_ACTIVE = (FB(NM)%MODE(KIND)/=0)
END FUNCTION FDS_FLUX_ACTIVE

!> DIF hook (called from DIVERGENCE_PART_1 after the species-sum correction, patch 0009): mode 1 copies the faces into the registered arrays, mode 2 replaces the listed values.
!> N0: lower bound of the species index of the arrays (0 or 1; species 1..N_TOTAL_SCALARS are used).
SUBROUTINE FDS_HOOK_DIF_FLUX(NM,N0,RX,RY,RZ)
INTEGER, INTENT(IN) :: NM,N0
REAL(EB), INTENT(INOUT) :: RX(0:,0:,0:,N0:),RY(0:,0:,0:,N0:),RZ(0:,0:,0:,N0:)
INTEGER :: E,N,D,I,J,K
IF (.NOT.ALLOCATED(FB)) RETURN
IF (NM<1 .OR. NM>SIZE(FB)) RETURN
SELECT CASE(FB(NM)%MODE(KDIF))
   CASE(1)
      IF (ASSOCIATED(FB(NM)%V(KDIF,1)%A)) THEN
         DO N=1,N_TOTAL_SCALARS ; DO K=1,KBAR ; DO J=1,JBAR ; DO I=0,IBAR ; FB(NM)%V(KDIF,1)%A(I,J,K,N) = RX(I,J,K,N) ; ENDDO ; ENDDO ; ENDDO ; ENDDO
      ENDIF
      IF (ASSOCIATED(FB(NM)%V(KDIF,2)%A)) THEN
         DO N=1,N_TOTAL_SCALARS ; DO K=1,KBAR ; DO J=0,JBAR ; DO I=1,IBAR ; FB(NM)%V(KDIF,2)%A(I,J,K,N) = RY(I,J,K,N) ; ENDDO ; ENDDO ; ENDDO ; ENDDO
      ENDIF
      IF (ASSOCIATED(FB(NM)%V(KDIF,3)%A)) THEN
         DO N=1,N_TOTAL_SCALARS ; DO K=0,KBAR ; DO J=1,JBAR ; DO I=1,IBAR ; FB(NM)%V(KDIF,3)%A(I,J,K,N) = RZ(I,J,K,N) ; ENDDO ; ENDDO ; ENDDO ; ENDDO
      ENDIF
   CASE(2)
      DO E=1,FB(NM)%OVR(KDIF)%N
         D = FB(NM)%OVR(KDIF)%DIR(E) ; I = FB(NM)%OVR(KDIF)%IDX(1,E) ; J = FB(NM)%OVR(KDIF)%IDX(2,E) ; K = FB(NM)%OVR(KDIF)%IDX(3,E)
         SELECT CASE(D)
            CASE(1) ; RX(I,J,K,1:N_TOTAL_SCALARS) = FB(NM)%OVR(KDIF)%VAL(1:N_TOTAL_SCALARS,E)
            CASE(2) ; RY(I,J,K,1:N_TOTAL_SCALARS) = FB(NM)%OVR(KDIF)%VAL(1:N_TOTAL_SCALARS,E)
            CASE(3) ; RZ(I,J,K,1:N_TOTAL_SCALARS) = FB(NM)%OVR(KDIF)%VAL(1:N_TOTAL_SCALARS,E)
         END SELECT
      ENDDO
END SELECT
END SUBROUTINE FDS_HOOK_DIF_FLUX

!> ADV hook (called from the generated density copy once UU,VV,WW hold the velocities the update uses, interface-wall values restored): the stage product FX*UU (FY*VV, FZ*WW).
!> Mode 1: copies it into the registered arrays and sets READ_ONLY (the caller skips the update). Mode 2 with a non-empty list: builds the product arrays PRX,PRY,PRZ, replaces
!> the listed faces and sets OVR (the caller runs its second copy of the update loop on PRX,PRY,PRZ); with an empty list OVR stays false and the caller's own loop runs unchanged.
SUBROUTINE FDS_HOOK_ADV(NM,N0,LU,LV,LW,UU,VV,WW,FX,FY,FZ,READ_ONLY,OVR)
INTEGER, INTENT(IN) :: NM,N0,LU(3),LV(3),LW(3)   !< lower bounds of UU, VV, WW (the work arrays WORK_U.. start at -1 in their own direction)
REAL(EB), INTENT(IN) :: UU(LU(1):,LU(2):,LU(3):),VV(LV(1):,LV(2):,LV(3):),WW(LW(1):,LW(2):,LW(3):),FX(0:,0:,0:,N0:),FY(0:,0:,0:,N0:),FZ(0:,0:,0:,N0:)
LOGICAL, INTENT(OUT) :: READ_ONLY,OVR
INTEGER :: E,N,D,I,J,K,NS
READ_ONLY = .FALSE. ; OVR = .FALSE.
IF (.NOT.ALLOCATED(FB)) RETURN
IF (NM<1 .OR. NM>SIZE(FB)) RETURN
SELECT CASE(FB(NM)%MODE(KADV))
   CASE(1)
      READ_ONLY = .TRUE.
      IF (ASSOCIATED(FB(NM)%V(KADV,1)%A)) THEN
         DO N=1,N_TOTAL_SCALARS ; DO K=1,KBAR ; DO J=1,JBAR ; DO I=0,IBAR ; FB(NM)%V(KADV,1)%A(I,J,K,N) = FX(I,J,K,N)*UU(I,J,K) ; ENDDO ; ENDDO ; ENDDO ; ENDDO
      ENDIF
      IF (ASSOCIATED(FB(NM)%V(KADV,2)%A)) THEN
         DO N=1,N_TOTAL_SCALARS ; DO K=1,KBAR ; DO J=0,JBAR ; DO I=1,IBAR ; FB(NM)%V(KADV,2)%A(I,J,K,N) = FY(I,J,K,N)*VV(I,J,K) ; ENDDO ; ENDDO ; ENDDO ; ENDDO
      ENDIF
      IF (ASSOCIATED(FB(NM)%V(KADV,3)%A)) THEN
         DO N=1,N_TOTAL_SCALARS ; DO K=0,KBAR ; DO J=1,JBAR ; DO I=1,IBAR ; FB(NM)%V(KADV,3)%A(I,J,K,N) = FZ(I,J,K,N)*WW(I,J,K) ; ENDDO ; ENDDO ; ENDDO ; ENDDO
      ENDIF
   CASE(2)
      IF (FB(NM)%OVR(KADV)%N==0) RETURN
      NS = N_TOTAL_SCALARS
      IF (ALLOCATED(FB(NM)%PROD)) THEN
         IF (ANY(SHAPE(FB(NM)%PROD)/=[IBAR+1,JBAR+1,KBAR+1,NS,3])) DEALLOCATE(FB(NM)%PROD)
      ENDIF
      IF (.NOT.ALLOCATED(FB(NM)%PROD)) ALLOCATE(FB(NM)%PROD(0:IBAR,0:JBAR,0:KBAR,NS,3))
      FB(NM)%PROD = 0._EB
      DO N=1,NS
         DO K=1,KBAR ; DO J=1,JBAR ; DO I=0,IBAR ; FB(NM)%PROD(I,J,K,N,1) = FX(I,J,K,N)*UU(I,J,K) ; ENDDO ; ENDDO ; ENDDO
         DO K=1,KBAR ; DO J=0,JBAR ; DO I=1,IBAR ; FB(NM)%PROD(I,J,K,N,2) = FY(I,J,K,N)*VV(I,J,K) ; ENDDO ; ENDDO ; ENDDO
         DO K=0,KBAR ; DO J=1,JBAR ; DO I=1,IBAR ; FB(NM)%PROD(I,J,K,N,3) = FZ(I,J,K,N)*WW(I,J,K) ; ENDDO ; ENDDO ; ENDDO
      ENDDO
      DO E=1,FB(NM)%OVR(KADV)%N
         D = FB(NM)%OVR(KADV)%DIR(E) ; I = FB(NM)%OVR(KADV)%IDX(1,E) ; J = FB(NM)%OVR(KADV)%IDX(2,E) ; K = FB(NM)%OVR(KADV)%IDX(3,E)
         FB(NM)%PROD(I,J,K,1:NS,D) = FB(NM)%OVR(KADV)%VAL(1:NS,E)
      ENDDO
      CALL POINT_PROD(FB(NM))
      OVR = .TRUE.
END SELECT
END SUBROUTINE FDS_HOOK_ADV

SUBROUTINE POINT_PROD(B)
TYPE(FLUX_BOX_T), INTENT(IN), TARGET :: B
PRX(0:,0:,0:,1:) => B%PROD(:,:,:,:,1)
PRY(0:,0:,0:,1:) => B%PROD(:,:,:,:,2)
PRZ(0:,0:,0:,1:) => B%PROD(:,:,:,:,3)
END SUBROUTINE POINT_PROD

END MODULE FDS_FLUX_HOOKS
