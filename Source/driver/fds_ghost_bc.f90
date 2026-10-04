!> \brief Boundary-condition step of the AMR driver for box-to-box and periodic boundaries of the same level (M2a, S4).
!>
!> FDS gets the values of the cells across a mesh boundary (box-to-box interface or periodic image) from OMESH(NOM), the exchanged copy of the neighbouring mesh, and
!> the boundary routines VISCOSITY_BC, VELOCITY_BC, MATCH_VELOCITY and WALL_BC (ASSIGN_GHOST_VALUE) read it. The AMR driver keeps that mechanism for the M2a
!> same-level case: FDS_G_FILL_OM is the block copy of MESH_EXCHANGE (data of box NOM, ghost cells included, sent by the driver to every rank), the boundary
!> routines below are the unmodified FDS routines, so the ghost values they write are FDS's own, bitwise, by construction. The copy is the only place that knows
!> where the neighbour data comes from; the guarded patches of Source/driver/patches (0003, 0004) are the hooks for replacing it by AMReX ghost data.
!>
!> Kernel-facing file rules (M2a):
!> (a) Passive scalars (N_TOTAL_SCALARS beyond the tracked species) are handled by Fields.cpp: ZZ/ZZS carry ncomp=N_TOTAL_SCALARS and the copy here
!>     loops over the last extent, whatever its size.
!> (b) Only uniform Cartesian metrics are used; R(I)/RRN(I) are not involved and CYLINDRICAL/TRN* meshes are rejected at level-0 assembly (IR-002).
!> Boundary table: README.md, "Fortran/C++ boundary".
MODULE FDS_GHOST_BC

USE ISO_C_BINDING
USE PRECISION_PARAMETERS, ONLY: EB
USE GLOBAL_CONSTANTS, ONLY: PREDICTOR,CORRECTOR
USE MESH_VARIABLES, ONLY: MESHES,MESH_TYPE
USE TYPES, ONLY: OMESH_TYPE
USE VELO, ONLY: MATCH_VELOCITY,MATCH_VELOCITY_FLUX,VELOCITY_BC,VISCOSITY_BC
USE WALL_ROUTINES, ONLY: WALL_BC

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE

PUBLIC :: FDS_G_FILL_OM,FDS_G_PHASE,FDS_G_MATCH,FDS_G_MATCH_FLUX,FDS_G_VELOCITY_BC,FDS_G_VISCOSITY_BC,FDS_G_WALL_BC,FDS_G_MU_EDGES,FDS_G_MU_EDGES_DOM

CONTAINS

!> Copy of one array of box NOM into OMESH(NOM) of box NM: the block copy that the same-rank branch of MESH_EXCHANGE does (OM%X(IMIN:IMAX,..) = MESHES(NOM)%X(IMIN:IMAX,..)),
!> whole allocated OMESH extent, ghost cells of NOM included. DATA is the native FAB of box NOM (FDS lower bounds LB(3), extents EXT(3), NC components; RHO/RHOS
!> have three ghost layers there), sent to every rank by the driver. WHICH: 1 MU 2 RHO 3 RHOS 4 U 5 V 6 W 7 US 8 VS 9 WS 10 H 11 HS 12 FVX 13 FVY 14 FVZ 15 D 16 DS
!> 17 KRES 18 Q 19 ZZ 20 ZZS. Returns 0, or 1 if NM has no OMESH(NOM) array of that kind (not a neighbour).
FUNCTION FDS_G_FILL_OM(NM,NOM,WHICH,LB,EXT,NC,DATA) BIND(C,NAME='fds_g_fill_om') RESULT(IERR)
INTEGER(C_INT), VALUE :: NM,NOM,WHICH,NC
INTEGER(C_INT), INTENT(IN) :: LB(3),EXT(3)
REAL(C_DOUBLE), INTENT(IN) :: DATA(*)
INTEGER(C_INT) :: IERR
TYPE(OMESH_TYPE), POINTER :: OM
INTEGER :: I,J,K,N,IS,JS,KS
IERR = 1
IF (.NOT.ALLOCATED(MESHES(NM)%OMESH)) RETURN
IF (NOM<LBOUND(MESHES(NM)%OMESH,1) .OR. NOM>UBOUND(MESHES(NM)%OMESH,1)) RETURN   ! NOM is not a mesh known to NM (for example a fine box)
OM => MESHES(NM)%OMESH(NOM)
SELECT CASE(WHICH)
   CASE( 1) ; IF (ALLOCATED(OM%MU))   CALL C3(OM%MU)
   CASE( 2) ; IF (ALLOCATED(OM%RHO))  CALL C3(OM%RHO)
   CASE( 3) ; IF (ALLOCATED(OM%RHOS)) CALL C3(OM%RHOS)
   CASE( 4) ; IF (ALLOCATED(OM%U))    CALL C3(OM%U)
   CASE( 5) ; IF (ALLOCATED(OM%V))    CALL C3(OM%V)
   CASE( 6) ; IF (ALLOCATED(OM%W))    CALL C3(OM%W)
   CASE( 7) ; IF (ALLOCATED(OM%US))   CALL C3(OM%US)
   CASE( 8) ; IF (ALLOCATED(OM%VS))   CALL C3(OM%VS)
   CASE( 9) ; IF (ALLOCATED(OM%WS))   CALL C3(OM%WS)
   CASE(10) ; IF (ALLOCATED(OM%H))    CALL C3(OM%H)
   CASE(11) ; IF (ALLOCATED(OM%HS))   CALL C3(OM%HS)
   CASE(12) ; IF (ALLOCATED(OM%FVX))  CALL C3(OM%FVX)
   CASE(13) ; IF (ALLOCATED(OM%FVY))  CALL C3(OM%FVY)
   CASE(14) ; IF (ALLOCATED(OM%FVZ))  CALL C3(OM%FVZ)
   CASE(15) ; IF (ALLOCATED(OM%D))    CALL C3(OM%D)
   CASE(16) ; IF (ALLOCATED(OM%DS))   CALL C3(OM%DS)
   CASE(17) ; IF (ALLOCATED(OM%KRES)) CALL C3(OM%KRES)
   CASE(18) ; IF (ALLOCATED(OM%Q))    CALL C3(OM%Q)
   CASE(19) ; IF (ALLOCATED(OM%ZZ))   CALL C4(OM%ZZ)
   CASE(20) ; IF (ALLOCATED(OM%ZZS))  CALL C4(OM%ZZS)
END SELECT

CONTAINS

INTEGER FUNCTION IX(I,J,K,N)
INTEGER, INTENT(IN) :: I,J,K,N
IX = 1 + (I-LB(1)) + EXT(1)*((J-LB(2)) + EXT(2)*((K-LB(3)) + EXT(3)*(N-1)))
END FUNCTION IX

LOGICAL FUNCTION INSIDE(I,J,K)
INTEGER, INTENT(IN) :: I,J,K
INSIDE = I>=LB(1) .AND. I<LB(1)+EXT(1) .AND. J>=LB(2) .AND. J<LB(2)+EXT(2) .AND. K>=LB(3) .AND. K<LB(3)+EXT(3)
END FUNCTION INSIDE

SUBROUTINE C3(A)
REAL(EB), ALLOCATABLE, INTENT(INOUT) :: A(:,:,:)
DO K=LBOUND(A,3),UBOUND(A,3)
   DO J=LBOUND(A,2),UBOUND(A,2)
      DO I=LBOUND(A,1),UBOUND(A,1)
         IF (INSIDE(I,J,K)) A(I,J,K) = DATA(IX(I,J,K,1))
      ENDDO
   ENDDO
ENDDO
IERR = 0
END SUBROUTINE C3

SUBROUTINE C4(A)
REAL(EB), ALLOCATABLE, INTENT(INOUT) :: A(:,:,:,:)
DO N=1,MIN(UBOUND(A,4),NC)
   DO K=LBOUND(A,3),UBOUND(A,3)
      DO J=LBOUND(A,2),UBOUND(A,2)
         DO I=LBOUND(A,1),UBOUND(A,1)
            IF (INSIDE(I,J,K)) A(I,J,K,N) = DATA(IX(I,J,K,N))
         ENDDO
      ENDDO
   ENDDO
ENDDO
IERR = 0
END SUBROUTINE C4

END FUNCTION FDS_G_FILL_OM

!> Set the stage flag that FDS keeps in module variables (PREDICTOR/CORRECTOR) before a boundary-condition call.
SUBROUTINE FDS_G_PHASE(PRED) BIND(C,NAME='fds_g_phase')
INTEGER(C_INT), VALUE :: PRED
PREDICTOR = (PRED/=0) ; CORRECTOR = .NOT.PREDICTOR
END SUBROUTINE FDS_G_PHASE

!> MATCH_VELOCITY(NM), unmodified (PREDICTOR/CORRECTOR as set by FDS_K_STATE); OMESH must have been filled by FDS_G_FILL_OMESH.
SUBROUTINE FDS_G_MATCH(NM) BIND(C,NAME='fds_g_match')
INTEGER(C_INT), VALUE :: NM
CALL MATCH_VELOCITY(NM)
END SUBROUTINE FDS_G_MATCH

!> MATCH_VELOCITY_FLUX(NM), unmodified (returns at once for a single FDS mesh).
SUBROUTINE FDS_G_MATCH_FLUX(NM) BIND(C,NAME='fds_g_match_flux')
INTEGER(C_INT), VALUE :: NM
CALL MATCH_VELOCITY_FLUX(NM)
END SUBROUTINE FDS_G_MATCH_FLUX

!> VELOCITY_BC(T,NM,EST), unmodified.
SUBROUTINE FDS_G_VELOCITY_BC(T,NM,EST) BIND(C,NAME='fds_g_velocity_bc')
REAL(C_DOUBLE), VALUE :: T
INTEGER(C_INT), VALUE :: NM,EST
CALL VELOCITY_BC(T,NM,APPLY_TO_ESTIMATED_VARIABLES=(EST/=0))
END SUBROUTINE FDS_G_VELOCITY_BC

!> VISCOSITY_BC(NM,EST), unmodified.
SUBROUTINE FDS_G_VISCOSITY_BC(NM,EST) BIND(C,NAME='fds_g_viscosity_bc')
INTEGER(C_INT), VALUE :: NM,EST
CALL VISCOSITY_BC(NM,APPLY_TO_ESTIMATED_VARIABLES=(EST/=0))
END SUBROUTINE FDS_G_VISCOSITY_BC

!> WALL_BC(T,DT,NM), unmodified.
SUBROUTINE FDS_G_WALL_BC(T,DT,NM) BIND(C,NAME='fds_g_wall_bc')
REAL(C_DOUBLE), VALUE :: T,DT
INTEGER(C_INT), VALUE :: NM
CALL WALL_BC(T,DT,NM)
END SUBROUTINE FDS_G_WALL_BC

!> Copy of MU and KRES into the edge cells of the domain (the closing lines of COMPUTE_VISCOSITY, velo.f90): the clamped copies that VELOCITY_FLUX reads
!> through the edge averages of MU. In a time step they are written by the viscosity kernel; the kernel check calls this after its boundary replay because
!> the replay restores no kernel output.
SUBROUTINE FDS_G_MU_EDGES(NM) BIND(C,NAME='fds_g_mu_edges')
INTEGER(C_INT), VALUE :: NM
TYPE(MESH_TYPE), POINTER :: M
INTEGER :: IBAR,JBAR,KBAR,IBP1,JBP1,KBP1
M => MESHES(NM)
IBAR=M%IBAR ; JBAR=M%JBAR ; KBAR=M%KBAR ; IBP1=IBAR+1 ; JBP1=JBAR+1 ; KBP1=KBAR+1
CALL EDG(M%MU)
CALL EDG(M%KRES)
CONTAINS
SUBROUTINE EDG(A)
REAL(EB), INTENT(INOUT) :: A(0:,0:,0:)
A(   0,0:JBP1,   0) = A(   1,0:JBP1,1)
A(IBP1,0:JBP1,   0) = A(IBAR,0:JBP1,1)
A(IBP1,0:JBP1,KBP1) = A(IBAR,0:JBP1,KBAR)
A(   0,0:JBP1,KBP1) = A(   1,0:JBP1,KBAR)
A(0:IBP1,   0,   0) = A(0:IBP1,   1,1)
A(0:IBP1,JBP1,0)    = A(0:IBP1,JBAR,1)
A(0:IBP1,JBP1,KBP1) = A(0:IBP1,JBAR,KBAR)
A(0:IBP1,0,KBP1)    = A(0:IBP1,   1,KBAR)
A(0,   0,0:KBP1)    = A(   1,   1,0:KBP1)
A(IBP1,0,0:KBP1)    = A(IBAR,   1,0:KBP1)
A(IBP1,JBP1,0:KBP1) = A(IBAR,JBAR,0:KBP1)
A(0,JBP1,0:KBP1)    = A(   1,JBAR,0:KBP1)
END SUBROUTINE EDG
END SUBROUTINE FDS_G_MU_EDGES

!> Time-loop form of FDS_G_MU_EDGES: only the edge cells that lie on an edge of the DOMAIN are written. FDS does the clamped copy at the end of COMPUTE_VISCOSITY;
!> the driver's full ghost fill (periodic images, box neighbours) runs after it and would replace those edge cells by the periodic/neighbour values, which FDS
!> does not use (the edge averages of MU and KRES in VELOCITY_FLUX read them). The interface edges of a decomposed mesh are left to the ghost fill so that the
!> result does not depend on the box layout. MASK bits: 1 low-x, 2 high-x, 4 low-y, 8 high-y, 16 low-z, 32 high-z are set when the box face lies on the domain boundary
!> (periodic or not); a statement is applied when both sides it names are domain sides. WHICH bit 0 = MU, bit 1 = KRES (both normally; one is dropped only by the fault-injection
!> switch FDSTL_SKIP_FIX of the regression test).
SUBROUTINE FDS_G_MU_EDGES_DOM(NM,MASK,WHICH) BIND(C,NAME='fds_g_mu_edges_dom')
INTEGER(C_INT), VALUE :: NM,MASK,WHICH
TYPE(MESH_TYPE), POINTER :: M
INTEGER :: IBAR,JBAR,KBAR,IBP1,JBP1,KBP1
LOGICAL :: XL,XH,YL,YH,ZL,ZH
M => MESHES(NM)
IBAR=M%IBAR ; JBAR=M%JBAR ; KBAR=M%KBAR ; IBP1=IBAR+1 ; JBP1=JBAR+1 ; KBP1=KBAR+1
XL=BTEST(MASK,0) ; XH=BTEST(MASK,1) ; YL=BTEST(MASK,2) ; YH=BTEST(MASK,3) ; ZL=BTEST(MASK,4) ; ZH=BTEST(MASK,5)
IF (BTEST(WHICH,0)) CALL EDG(M%MU)
IF (BTEST(WHICH,1)) CALL EDG(M%KRES)
CONTAINS
SUBROUTINE EDG(A)
REAL(EB), INTENT(INOUT) :: A(0:,0:,0:)
IF (XL.AND.ZL) A(   0,0:JBP1,   0) = A(   1,0:JBP1,1)
IF (XH.AND.ZL) A(IBP1,0:JBP1,   0) = A(IBAR,0:JBP1,1)
IF (XH.AND.ZH) A(IBP1,0:JBP1,KBP1) = A(IBAR,0:JBP1,KBAR)
IF (XL.AND.ZH) A(   0,0:JBP1,KBP1) = A(   1,0:JBP1,KBAR)
IF (YL.AND.ZL) A(0:IBP1,   0,   0) = A(0:IBP1,   1,1)
IF (YH.AND.ZL) A(0:IBP1,JBP1,0)    = A(0:IBP1,JBAR,1)
IF (YH.AND.ZH) A(0:IBP1,JBP1,KBP1) = A(0:IBP1,JBAR,KBAR)
IF (YL.AND.ZH) A(0:IBP1,0,KBP1)    = A(0:IBP1,   1,KBAR)
IF (XL.AND.YL) A(0,   0,0:KBP1)    = A(   1,   1,0:KBP1)
IF (XH.AND.YL) A(IBP1,0,0:KBP1)    = A(IBAR,   1,0:KBP1)
IF (XH.AND.YH) A(IBP1,JBP1,0:KBP1) = A(IBAR,JBAR,0:KBP1)
IF (XL.AND.YH) A(0,JBP1,0:KBP1)    = A(   1,JBAR,0:KBP1)
! The last four statements of FDS write the eight corner cells from the z-ghost cells A(.,.,0) and A(.,.,KBP1), which at that point of COMPUTE_VISCOSITY hold the
! mirror of the adjacent gas cell (WALL_LOOP_2), not the periodic image that the driver has there. So a corner is the adjacent interior corner cell.
IF (XL.AND.YL.AND.ZL) A(   0,   0,   0) = A(   1,   1,   1)
IF (XL.AND.YL.AND.ZH) A(   0,   0,KBP1) = A(   1,   1,KBAR)
IF (XH.AND.YL.AND.ZL) A(IBP1,   0,   0) = A(IBAR,   1,   1)
IF (XH.AND.YL.AND.ZH) A(IBP1,   0,KBP1) = A(IBAR,   1,KBAR)
IF (XH.AND.YH.AND.ZL) A(IBP1,JBP1,   0) = A(IBAR,JBAR,   1)
IF (XH.AND.YH.AND.ZH) A(IBP1,JBP1,KBP1) = A(IBAR,JBAR,KBAR)
IF (XL.AND.YH.AND.ZL) A(   0,JBP1,   0) = A(   1,JBAR,   1)
IF (XL.AND.YH.AND.ZH) A(   0,JBP1,KBP1) = A(   1,JBAR,KBAR)
END SUBROUTINE EDG
END SUBROUTINE FDS_G_MU_EDGES_DOM

END MODULE FDS_GHOST_BC
