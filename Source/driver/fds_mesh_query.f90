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
USE GLOBAL_CONSTANTS, ONLY: NMESHES,PROCESS,PERIODIC_DOMAIN_X,PERIODIC_DOMAIN_Y,PERIODIC_DOMAIN_Z,CYLINDRICAL,N_TRACKED_SPECIES, &
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

END MODULE FDS_MESH_QUERY
