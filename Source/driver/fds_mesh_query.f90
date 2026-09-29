!> \brief Read-only queries of the FDS mesh set-up data for the C++ driver (bind(C)); no kernel is called from here.
!>
!> Kernel-facing file rules (M2a):
!> (a) Passive scalars (N_TOTAL_SCALARS beyond the tracked species) are not touched here. They are handled by Fields.cpp (S2), which
!>     sizes the RHO_ZZ component list from N_TOTAL_SCALARS; this file only reports the two counts.
!> (b) Only uniform Cartesian metrics are used. R(I)/RRN(I) are dropped; CYLINDRICAL and TRNX/TRNY/TRNZ meshes are reported
!>     to the driver, which rejects them (IR-002).
!>
!> Boundary table entry: see README.md, "Fortran/C++ boundary".

MODULE FDS_MESH_QUERY

USE ISO_C_BINDING, ONLY: C_INT,C_DOUBLE
USE MESH_VARIABLES, ONLY: MESHES
USE GLOBAL_CONSTANTS, ONLY: INTERPOLATED_BOUNDARY,NMESHES,PROCESS,PERIODIC_DOMAIN_X,PERIODIC_DOMAIN_Y,PERIODIC_DOMAIN_Z,CYLINDRICAL,N_TRACKED_SPECIES, &
                            N_TOTAL_SCALARS,MY_RANK,N_MPI_PROCESSES

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE

CONTAINS

!> \brief Number of FDS meshes (the same on every MPI rank)
FUNCTION FDS_GET_NMESHES() BIND(C,NAME='fds_get_nmeshes') RESULT(N)
INTEGER(C_INT) :: N
N = NMESHES
END FUNCTION FDS_GET_NMESHES

!> \brief Cell counts, physical extent, owning rank and a nonuniform-grid flag of mesh NM (1-based FDS index)
SUBROUTINE FDS_GET_MESH(NM,IJK,XB,RANK,NONUNIFORM) BIND(C,NAME='fds_get_mesh')
INTEGER(C_INT), VALUE :: NM
INTEGER(C_INT), INTENT(OUT) :: IJK(3),RANK,NONUNIFORM
REAL(C_DOUBLE), INTENT(OUT) :: XB(6)
IJK = (/MESHES(NM)%IBAR,MESHES(NM)%JBAR,MESHES(NM)%KBAR/)
XB  = (/MESHES(NM)%XS,MESHES(NM)%XF,MESHES(NM)%YS,MESHES(NM)%YF,MESHES(NM)%ZS,MESHES(NM)%ZF/)
RANK = PROCESS(NM)
NONUNIFORM = 0
IF (MESHES(NM)%TRNX_ID/='null' .OR. MESHES(NM)%TRNY_ID/='null' .OR. MESHES(NM)%TRNZ_ID/='null') NONUNIFORM = 1
END SUBROUTINE FDS_GET_MESH

!> \brief Domain-level flags: periodicity per direction, cylindrical option, species counts, MPI size
SUBROUTINE FDS_GET_DOMAIN(PERIODIC,CYL,N_TRACKED,N_TOTAL,NRANKS) BIND(C,NAME='fds_get_domain')
INTEGER(C_INT), INTENT(OUT) :: PERIODIC(3),CYL,N_TRACKED,N_TOTAL,NRANKS
PERIODIC = 0
IF (PERIODIC_DOMAIN_X) PERIODIC(1) = 1
IF (PERIODIC_DOMAIN_Y) PERIODIC(2) = 1
IF (PERIODIC_DOMAIN_Z) PERIODIC(3) = 1
CYL = 0
IF (CYLINDRICAL) CYL = 1
N_TRACKED = N_TRACKED_SPECIES
N_TOTAL   = N_TOTAL_SCALARS
NRANKS    = N_MPI_PROCESSES
END SUBROUTINE FDS_GET_DOMAIN

!> \brief Cell wall data of the local mesh NM for the C++ side-data rebuild (SideData). FLAGS(0:6,I,J,K) for the valid cells:
!> 0 = CELL%SOLID (0/1); 1..6 = wall code of the faces (-x,+x,-y,+y,-z,+z) = WALL_INDEX(-1),(1),(-2),(2),(-3),(3):
!> 0 no wall cell, 1 wall cell, 2 mesh-to-mesh interface (INTERPOLATED_BOUNDARY with a neighbouring mesh NOM>0; the driver opens
!> these in the face mask). Periodic faces keep code 1 (single-mesh FDS has a wall cell there too, init.f90:76-107).
SUBROUTINE FDS_GET_CELL_WALLS(NM,N,FLAGS) BIND(C,NAME='fds_get_cell_walls')
INTEGER(C_INT), VALUE :: NM
INTEGER(C_INT), INTENT(IN) :: N(3)
INTEGER(C_INT), INTENT(OUT) :: FLAGS(0:6,N(1),N(2),N(3))
INTEGER :: I,J,K,IC,IW,F
INTEGER, PARAMETER :: CODES(6) = (/-1,1,-2,2,-3,3/)
DO K=1,N(3) ; DO J=1,N(2) ; DO I=1,N(1)
   IC = MESHES(NM)%CELL_INDEX(I,J,K)
   FLAGS(0,I,J,K) = 0
   IF (MESHES(NM)%CELL(IC)%SOLID) FLAGS(0,I,J,K) = 1
   DO F=1,6
      IW = MESHES(NM)%CELL(IC)%WALL_INDEX(CODES(F))
      IF (IW==0) THEN
         FLAGS(F,I,J,K) = 0
      ELSEIF (IW<=MESHES(NM)%N_EXTERNAL_WALL_CELLS .AND. MESHES(NM)%WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY &
              .AND. MESHES(NM)%EXTERNAL_WALL(IW)%NOM>0) THEN
         FLAGS(F,I,J,K) = 2
      ELSE
         FLAGS(F,I,J,K) = 1
      ENDIF
   ENDDO
ENDDO ; ENDDO ; ENDDO
END SUBROUTINE FDS_GET_CELL_WALLS

END MODULE FDS_MESH_QUERY
