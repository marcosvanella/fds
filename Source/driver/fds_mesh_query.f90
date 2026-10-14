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

USE ISO_C_BINDING, ONLY: C_INT,C_DOUBLE,C_CHAR,C_NULL_CHAR
USE MESH_VARIABLES, ONLY: MESHES
USE TYPES, ONLY: SPECIES_MIXTURE
USE GLOBAL_CONSTANTS, ONLY: INTERPOLATED_BOUNDARY,NULL_BOUNDARY,NMESHES,PROCESS,PERIODIC_DOMAIN_X,PERIODIC_DOMAIN_Y,PERIODIC_DOMAIN_Z,CYLINDRICAL,N_TRACKED_SPECIES, &
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

!> \brief Wall type of the six faces (-x,+x,-y,+y,-z,+z) of local mesh NM for the fine-level velocity BC (GhostExchange mirror_domain_edges, patch 0015): the
!> WALL BOUNDARY_TYPE (cons.f90: 1 SOLID/INERT, 3 MIRROR) shared by every external wall cell of that face, 0 when they differ, -1 when there is none. NULL walls are skipped.
SUBROUTINE FDS_GET_FACE_BC(NM,CODE) BIND(C,NAME='fds_get_face_bc')
INTEGER(C_INT), VALUE :: NM
INTEGER(C_INT), INTENT(OUT) :: CODE(6)
INTEGER :: IW,IO,F,BT
CODE = -1
DO IW=1,MESHES(NM)%N_EXTERNAL_WALL_CELLS
   BT = MESHES(NM)%WALL(IW)%BOUNDARY_TYPE
   IF (BT==NULL_BOUNDARY) CYCLE
   IO = MESHES(NM)%BOUNDARY_COORD(MESHES(NM)%WALL(IW)%BC_INDEX)%IOR
   F = 2*ABS(IO) - MERGE(1,0,IO>0)   ! IOR=+1 is the wall of the -x face (its normal points into the mesh)
   IF (CODE(F)==-1) CODE(F) = BT
   IF (CODE(F)/=BT) CODE(F) = 0
ENDDO
END SUBROUTINE FDS_GET_FACE_BC

!> \brief 0-based ZZ component for tracked species / tracer ID (SPECIES_MIXTURE%ID), or -1 if unknown.
!> NAME is a NUL-terminated C string. Passive search is exact on TRIM(ID) (case-sensitive, as FDS namelist IDs).
FUNCTION FDS_FIND_TRACKED_SPECIES(NAME) BIND(C,NAME='fds_find_tracked_species') RESULT(IDX)
CHARACTER(KIND=C_CHAR), INTENT(IN) :: NAME(*)
INTEGER(C_INT) :: IDX
CHARACTER(LEN=60) :: ID
INTEGER :: I,N
ID = ' '
DO I=1,60
   IF (NAME(I)==C_NULL_CHAR) EXIT
   ID(I:I) = NAME(I)
ENDDO
IDX = -1
IF (.NOT.ALLOCATED(SPECIES_MIXTURE)) RETURN
DO N=1,N_TRACKED_SPECIES
   IF (TRIM(SPECIES_MIXTURE(N)%ID)==TRIM(ID)) THEN
      IDX = N - 1   ! C++ ZZ component (BACKGROUND first => 0)
      RETURN
   ENDIF
ENDDO
END FUNCTION FDS_FIND_TRACKED_SPECIES

!> \brief Number of OBSTs of the meshes of this rank (OBSTRUCTION(1:N_OBST); a multi-mesh OBST counts once per mesh it was clipped to)
FUNCTION FDS_OBST_COUNT() BIND(C,NAME='fds_obst_count') RESULT(N)
INTEGER(C_INT) :: N
INTEGER :: NM
N = 0
DO NM=1,NMESHES
   IF (PROCESS(NM)/=MY_RANK) CYCLE
   N = N + MESHES(NM)%N_OBST
ENDDO
END FUNCTION FDS_OBST_COUNT

!> \brief XB of those OBSTs as READ_OBST stored them (OB%X1..Z2: input XB clipped to the mesh, input nudges kept), 6 per OBST, same order
SUBROUTINE FDS_OBST_XB(N,XB) BIND(C,NAME='fds_obst_xb')
INTEGER(C_INT), VALUE :: N
REAL(C_DOUBLE), INTENT(OUT) :: XB(6,N)
INTEGER :: NM,I,M
M = 0
DO NM=1,NMESHES
   IF (PROCESS(NM)/=MY_RANK) CYCLE
   DO I=1,MESHES(NM)%N_OBST
      M = M + 1
      IF (M>N) RETURN
      ASSOCIATE (OB=>MESHES(NM)%OBSTRUCTION(I))
         XB(:,M) = [OB%X1,OB%X2,OB%Y1,OB%Y2,OB%Z1,OB%Z2]
      END ASSOCIATE
   ENDDO
ENDDO
END SUBROUTINE FDS_OBST_XB

!> \brief D-010: XB of the stateful VENTs of the local meshes (SURF with a specified mass flux (MASS_FLUX, HRRPUA, MLRPUA), a specified heat flux or a pyrolysis model),
!> 6 per VENT. With XB = NULL (N = 0) only the count is returned. These walls stay on level 0 (no-refine buffer, RegridAmrCore::ManualTagsPlacement).
FUNCTION FDS_STATEFUL_VENTS(N,XB) BIND(C,NAME='fds_stateful_vents') RESULT(M)
USE TYPES, ONLY: SURFACE
USE GLOBAL_CONSTANTS, ONLY: SPECIFIED_MASS_FLUX,NET_FLUX_BC,CONVECTIVE_FLUX_BC,PYROLYSIS_NONE
INTEGER(C_INT), VALUE :: N
REAL(C_DOUBLE), INTENT(OUT) :: XB(6,*)
INTEGER(C_INT) :: M
INTEGER :: NM,I
M = 0
DO NM=1,NMESHES
   IF (PROCESS(NM)/=MY_RANK) CYCLE
   DO I=1,MESHES(NM)%N_VENT
      ASSOCIATE (VT=>MESHES(NM)%VENTS(I))
      ASSOCIATE (SF=>SURFACE(VT%SURF_INDEX))
      IF (SF%SPECIES_BC_INDEX/=SPECIFIED_MASS_FLUX .AND. SF%THERMAL_BC_INDEX/=NET_FLUX_BC .AND. SF%THERMAL_BC_INDEX/=CONVECTIVE_FLUX_BC &
          .AND. SF%PYROLYSIS_MODEL==PYROLYSIS_NONE) CYCLE
      M = M + 1
      IF (M<=N) XB(:,M) = [VT%X1,VT%X2,VT%Y1,VT%Y2,VT%Z1,VT%Z2]
      END ASSOCIATE
      END ASSOCIATE
   ENDDO
ENDDO
END FUNCTION FDS_STATEFUL_VENTS

END MODULE FDS_MESH_QUERY
