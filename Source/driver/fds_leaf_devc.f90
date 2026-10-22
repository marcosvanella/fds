!> \brief Leaf-cell queries for the gas SPATIAL_STATISTIC / point devices of dump.f90 (UPDATE_DEVICES_1): which level-0 or fine cells are covered by a finer level, the fine boxes to add,
!> and the finest cell at a point. Registered by fds_leaf_register (TimeLoop::bind_level); without a fine level the hooks are never set and the devices are the level-0 ones unchanged.
!> A fine box belongs to the level-0 mesh whose extent the device range was clipped to (SDV%X1..Z2); a box of another rank than that mesh is added by its own rank (dump.f90 LEAF_FOREIGN_DEVICES).
MODULE FDS_LEAF_DEVC

USE ISO_C_BINDING
USE PRECISION_PARAMETERS, ONLY: EB
USE MESH_VARIABLES, ONLY: MESHES
USE MESH_POINTERS, ONLY: FINE_LEVEL
USE FDS_AMREX_HOOKS, ONLY: LEAF_COVERED_PROC,LEAF_BOX_PROC,LEAF_FINEST_PROC

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE
PUBLIC :: FDS_LEAF_REGISTER

CONTAINS

SUBROUTINE FDS_LEAF_REGISTER() BIND(C,NAME='fds_leaf_register')
LEAF_COVERED_PROC => LEAF_COVERED
LEAF_BOX_PROC     => LEAF_BOX
LEAF_FINEST_PROC  => LEAF_FINEST
END SUBROUTINE FDS_LEAF_REGISTER

!> Level of box number NM (0: a level-0 mesh) and its fine level index.
INTEGER FUNCTION BOX_LEVEL(NM)
INTEGER, INTENT(IN) :: NM
INTEGER :: LEV
BOX_LEVEL = 0
IF (NM<=SIZE(MESHES) .OR. .NOT.ALLOCATED(FINE_LEVEL)) RETURN
DO LEV=1,SIZE(FINE_LEVEL)
   IF (NM>FINE_LEVEL(LEV)%NM0 .AND. NM<=FINE_LEVEL(LEV)%NM0+FINE_LEVEL(LEV)%N_BOXES) THEN ; BOX_LEVEL = LEV ; RETURN ; ENDIF
ENDDO
END FUNCTION BOX_LEVEL

!> True when the point (X,Y,Z) of box NM lies inside a box of the next finer level.
LOGICAL FUNCTION LEAF_COVERED(NM,X,Y,Z)
INTEGER, INTENT(IN) :: NM
REAL(EB), INTENT(IN) :: X,Y,Z
INTEGER :: LEV,IB
LEAF_COVERED = .FALSE.
IF (.NOT.ALLOCATED(FINE_LEVEL)) RETURN
LEV = BOX_LEVEL(NM)+1
IF (LEV>SIZE(FINE_LEVEL)) RETURN
DO IB=1,FINE_LEVEL(LEV)%N_BOXES
   ASSOCIATE(R=>FINE_LEVEL(LEV)%BOX(IB))
   IF (X>R%XS .AND. X<R%XF .AND. Y>R%YS .AND. Y<R%YF .AND. Z>R%ZS .AND. Z<R%ZF) THEN ; LEAF_COVERED = .TRUE. ; RETURN ; ENDIF
   END ASSOCIATE
ENDDO
END FUNCTION LEAF_COVERED

!> The IB-th fine box (levels in order) as box number NB (0 past the last) and the cell range IR = (I1,I2,J1,J2,K1,K2) of its cells inside the extent X1..Z2, rounded as the level-0
!> device range (read.f90); an empty intersection gives I1>I2.
INTEGER FUNCTION LEAF_BOX(IB,X1,X2,Y1,Y2,Z1,Z2,IR)
INTEGER, INTENT(IN) :: IB
REAL(EB), INTENT(IN) :: X1,X2,Y1,Y2,Z1,Z2
INTEGER, INTENT(OUT) :: IR(6)
INTEGER :: LEV,N
LEAF_BOX = 0 ; IR = (/1,0,1,0,1,0/)
IF (.NOT.ALLOCATED(FINE_LEVEL)) RETURN
N = IB
DO LEV=1,SIZE(FINE_LEVEL)
   IF (N>FINE_LEVEL(LEV)%N_BOXES) THEN ; N = N-FINE_LEVEL(LEV)%N_BOXES ; CYCLE ; ENDIF
   ASSOCIATE(R=>FINE_LEVEL(LEV)%BOX(N))
   LEAF_BOX = FINE_LEVEL(LEV)%NM0+N
   IF (X2<R%XS .OR. X1>R%XF .OR. Y2<R%YS .OR. Y1>R%YF .OR. Z2<R%ZS .OR. Z1>R%ZF) RETURN
   CALL RANGE(X1,X2,R%XS,R%IBAR,R%RDXI,IR(1),IR(2))
   CALL RANGE(Y1,Y2,R%YS,R%JBAR,R%RDETA,IR(3),IR(4))
   CALL RANGE(Z1,Z2,R%ZS,R%KBAR,R%RDZETA,IR(5),IR(6))
   END ASSOCIATE
   RETURN
ENDDO
CONTAINS
   SUBROUTINE RANGE(A1,A2,AS,NB,RD,I1,I2)
   REAL(EB), INTENT(IN) :: A1,A2,AS,RD
   INTEGER, INTENT(IN) :: NB
   INTEGER, INTENT(OUT) :: I1,I2
   I1 = NINT(MAX(A1-AS,0._EB)*RD)+1 ; I2 = MIN(NB,NINT((A2-AS)*RD))
   IF (I1>I2 .AND. I2>=1 .AND. A1>=A2) I1 = I2   ! a zero-thickness device takes the layer below its plane, as at level 0 (a box that only touches the extent adds nothing)
   END SUBROUTINE RANGE
END FUNCTION LEAF_BOX

!> Finest fine box containing (X,Y,Z): box number (0 if none) and the cell indices of the point in it.
INTEGER FUNCTION LEAF_FINEST(X,Y,Z,I,J,K)
REAL(EB), INTENT(IN) :: X,Y,Z
INTEGER, INTENT(OUT) :: I,J,K
INTEGER :: LEV,IB
LEAF_FINEST = 0 ; I = 0 ; J = 0 ; K = 0
IF (.NOT.ALLOCATED(FINE_LEVEL)) RETURN
DO LEV=SIZE(FINE_LEVEL),1,-1
   DO IB=1,FINE_LEVEL(LEV)%N_BOXES
      ASSOCIATE(R=>FINE_LEVEL(LEV)%BOX(IB))
      IF (X>R%XS .AND. X<R%XF .AND. Y>R%YS .AND. Y<R%YF .AND. Z>R%ZS .AND. Z<R%ZF) THEN
         LEAF_FINEST = FINE_LEVEL(LEV)%NM0+IB
         I = MIN(R%IBAR,INT((X-R%XS)*R%RDXI)+1) ; J = MIN(R%JBAR,INT((Y-R%YS)*R%RDETA)+1) ; K = MIN(R%KBAR,INT((Z-R%ZS)*R%RDZETA)+1)
         RETURN
      ENDIF
      END ASSOCIATE
   ENDDO
ENDDO
END FUNCTION LEAF_FINEST

END MODULE FDS_LEAF_DEVC
