!> \brief The FDS mesh object of a box by mesh number: MESHES(NM) for a level-0 box and, once patch 0007 is applied (FDS_FINE_B defined by the CMake auto-detection),
!> FINE_LEVEL(L)%BOX(IB) for a fine-level box number (D-056 option B). Driver-owned Fortran uses BOX_OBJ(NM) wherever it used to write MESHES(NM) and may be called for a fine box.
MODULE FDS_BOX_OBJ

USE ISO_C_BINDING
USE MESH_VARIABLES, ONLY: MESHES,MESH_TYPE
#ifdef FDS_FINE_B
USE FDS_AMREX_HOOKS, ONLY: FDS_HOOK_BIND_VIEW
USE MESH_POINTERS, ONLY: FINE_LEVEL
#endif

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE
PUBLIC :: BOX_OBJ,HAS_M_DOT_PPP,FDS_FINE_B_SET_VIEW

CONTAINS

!> Pointer to the mesh object of box NM; for a number that belongs to no level (or a fine number without patch 0007) the run stops with a message.
FUNCTION BOX_OBJ(NM) RESULT(M)

INTEGER, INTENT(IN) :: NM
TYPE(MESH_TYPE), POINTER :: M
INTEGER :: LEV

M => NULL()
IF (NM>=1 .AND. NM<=SIZE(MESHES)) THEN
   M => MESHES(NM)
   RETURN
ENDIF
#ifdef FDS_FINE_B
IF (ALLOCATED(FINE_LEVEL)) THEN
   DO LEV=1,SIZE(FINE_LEVEL)
      IF (NM>FINE_LEVEL(LEV)%NM0 .AND. NM<=FINE_LEVEL(LEV)%NM0+FINE_LEVEL(LEV)%N_BOXES) THEN
         M => FINE_LEVEL(LEV)%BOX(NM-FINE_LEVEL(LEV)%NM0)
         RETURN
      ENDIF
   ENDDO
ENDIF
#endif
WRITE(0,'(A,I0,A)') 'fds_amr ERROR in BOX_OBJ: mesh number ',NM,' belongs to no level (level-0 FDS meshes, or fine-level boxes with patch 0007).'
ERROR STOP 1

END FUNCTION BOX_OBJ

!> ALLOCATED(MESHES(NM)%M_DOT_PPP) for any box (the generated density copy asks this).
FUNCTION HAS_M_DOT_PPP(NM) RESULT(L)

INTEGER, INTENT(IN) :: NM
LOGICAL :: L
TYPE(MESH_TYPE), POINTER :: M

M => BOX_OBJ(NM)
L = ALLOCATED(M%M_DOT_PPP)

END FUNCTION HAS_M_DOT_PPP

!> Bind a field of fine-level box IB of level L (the second array, FINE_LEVEL(L)%VIEW(IB)) to AMReX data: the counterpart of fds_hook_set_view for level > 0. WHICH, LB, EXT, P as there.
!> POINT_TO_BOX(NM) with NM = FINE_LEVEL(L)%NM0+IB then points the stage functions at that data. Returns 0, 1 for an unknown WHICH, 2 for an unknown level or box,
!> 3 when the build has no fine-level support (patch 0007 not applied).
FUNCTION FDS_FINE_B_SET_VIEW(L,IB,WHICH,LB,EXT,P) BIND(C,NAME='fds_fine_b_set_view') RESULT(IERR)

INTEGER(C_INT), VALUE :: L,IB,WHICH
INTEGER(C_INT), INTENT(IN) :: LB(4),EXT(4)
TYPE(C_PTR), VALUE :: P
INTEGER(C_INT) :: IERR

#ifdef FDS_FINE_B
IERR = 2
IF (.NOT.ALLOCATED(FINE_LEVEL)) RETURN
IF (L<1 .OR. L>SIZE(FINE_LEVEL)) RETURN
IF (IB<1 .OR. IB>FINE_LEVEL(L)%N_BOXES) RETURN
IF (.NOT.ALLOCATED(FINE_LEVEL(L)%VIEW)) ALLOCATE(FINE_LEVEL(L)%VIEW(FINE_LEVEL(L)%N_BOXES))
IERR = FDS_HOOK_BIND_VIEW(FINE_LEVEL(L)%VIEW(IB),WHICH,LB,EXT,P)
#else
IERR = 3
#endif

END FUNCTION FDS_FINE_B_SET_VIEW

END MODULE FDS_BOX_OBJ
