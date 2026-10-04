#!/usr/bin/env python3
"""Scratch-only hook: dump per-cell pressure-solve data from ULMAT_SOLVE_ZONE and GLMAT_SOLVER.
Usage: fds_dump_hook.py path/to/scratch/Source/pres.f90   (edits the given scratch copy in place)
Env at run time: FDS_DBG_RHS_OFFSET (real, default 0) adds offset*cell_volume to F_H before the mean removal;
FDS_DBG_FORCE_ULMAT (any value) makes ULMAT build its matrix instead of falling back to FFT on a regular mesh;
FDS_DBG_FIXMEAN (any value) makes the single-rank whole-domain GLMAT path remove the mean like the multi-rank path does."""
import sys
p = sys.argv[1]
s = open(p).read()

def rep(old, new, count=1):
    global s
    assert s.count(old) == count, (old, s.count(old))
    s = s.replace(old, new)

# ---------------- ULMAT ----------------
decl_u = """REAL(EB), POINTER, DIMENSION(:,:,:) :: HP,RHOP
INTEGER, SAVE :: DBG_NCALL=0
INTEGER :: DBG_LU,DBG_STAT
CHARACTER(80) :: DBG_FN,DBG_ENV
REAL(EB) :: DBG_OFF,DBG_VOL
REAL(EB), ALLOCATABLE, DIMENSION(:) :: DBG_FH0,DBG_FH1,DBG_XH0,DBG_XH1
"""
# first occurrence is in ULMAT_SOLVE_ZONE (line ~1482); the decl at 706 and 5887 are other routines
idx = s.index("REAL(EB), POINTER, DIMENSION(:,:,:) :: HP,RHOP\n", s.index("SUBROUTINE ULMAT_SOLVE_ZONE(NM,IPZ)"))
s = s[:idx] + decl_u + s[idx+len("REAL(EB), POINTER, DIMENSION(:,:,:) :: HP,RHOP\n"):]

rep("""! For indefinite matrices substract mean of source F_H:
H_INDEFINITE_IF_1 : IF (ZM%MTYPE==SYMM_INDEFINITE ) THEN
""", """DBG_NCALL = DBG_NCALL + 1
DBG_OFF = 0._EB
CALL GET_ENVIRONMENT_VARIABLE('FDS_DBG_RHS_OFFSET',DBG_ENV,STATUS=DBG_STAT)
IF (DBG_STAT==0) READ(DBG_ENV,*) DBG_OFF
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (MUNKH(I,J,K)<=0) CYCLE
         DBG_VOL = ((1._EB-CYL_FCT)*DY(J) + CYL_FCT*RC(I))*DX(I)*DZ(K)
         ZM%F_H(MUNKH(I,J,K)) = ZM%F_H(MUNKH(I,J,K)) + DBG_OFF*DBG_VOL
      ENDDO
   ENDDO
ENDDO
ALLOCATE(DBG_FH0(ZM%NUNKH),DBG_FH1(ZM%NUNKH),DBG_XH0(ZM%NUNKH),DBG_XH1(ZM%NUNKH))
DBG_FH0 = ZM%F_H
! For indefinite matrices substract mean of source F_H:
H_INDEFINITE_IF_1 : IF (ZM%MTYPE==SYMM_INDEFINITE ) THEN
""")
rep("""! Solve the system...

LIBRARY_SELECT: SELECT CASE(ULMAT_SOLVER_LIBRARY)""", """DBG_FH1 = ZM%F_H
! Solve the system...

LIBRARY_SELECT: SELECT CASE(ULMAT_SOLVER_LIBRARY)""")
rep("""! For indefinite matrices, substract mean of solution X_H:
H_INDEFINITE_IF_2 : IF (ZM%MTYPE==SYMM_INDEFINITE ) THEN
SUM_XH(1:2) = 0._EB; MEAN_XH = 0._EB""".replace("SUM_XH(1:2)", "   SUM_XH(1:2)"),
"""DBG_XH0 = ZM%X_H
! For indefinite matrices, substract mean of solution X_H:
H_INDEFINITE_IF_2 : IF (ZM%MTYPE==SYMM_INDEFINITE ) THEN
   SUM_XH(1:2) = 0._EB; MEAN_XH = 0._EB""")
rep("""IF (ZM%MTYPE==SYMM_INDEFINITE) THEN
   SUM_GAUGE = 0._EB
   IF (PREDICTOR) THEN
      RHOP => RHO""", """DBG_XH1 = ZM%X_H
IF (ZM%MTYPE==SYMM_INDEFINITE) THEN
   SUM_GAUGE = 0._EB
   IF (PREDICTOR) THEN
      RHOP => RHO""")
rep("""IF(CC_IBM) CALL GET_H_GUARD_CUTCELL(IPZ,HP)

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW
""", """IF(CC_IBM) CALL GET_H_GUARD_CUTCELL(IPZ,HP)

IF (DBG_NCALL<=4) THEN
   IF (PREDICTOR) THEN; RHOP => RHO; ELSE; RHOP => RHOS; ENDIF
   WRITE(DBG_FN,'(A,I0,A)') 'dbg_ULMAT_',DBG_NCALL,'.txt'
   OPEN(NEWUNIT=DBG_LU,FILE=TRIM(DBG_FN),STATUS='REPLACE',ACTION='WRITE')
   WRITE(DBG_LU,'(A,3I6,I3,I8)') 'DIMS',IBAR,JBAR,KBAR,MERGE(1,0,PREDICTOR),ZM%NUNKH
   WRITE(DBG_LU,'(A,I3,ES25.16)') 'MTYPE',ZM%MTYPE,DBG_OFF
   WRITE(DBG_LU,'(A,*(ES25.16))') 'DX ',DX(1:IBAR)
   WRITE(DBG_LU,'(A,*(ES25.16))') 'DY ',DY(1:JBAR)
   WRITE(DBG_LU,'(A,*(ES25.16))') 'DZ ',DZ(1:KBAR)
   WRITE(DBG_LU,'(A,*(ES25.16))') 'DXN',DXN(0:IBAR)
   WRITE(DBG_LU,'(A,*(ES25.16))') 'DYN',DYN(0:JBAR)
   WRITE(DBG_LU,'(A,*(ES25.16))') 'DZN',DZN(0:KBAR)
   WRITE(DBG_LU,'(A)') 'COLS I J K PRHS FH0 FH1 XH0 XH1 XHFIN HP RHOP KRES'
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (MUNKH(I,J,K)<=0) CYCLE
            IROW = MUNKH(I,J,K)
            WRITE(DBG_LU,'(3I6,9ES25.16)') I,J,K,PRHS(I,J,K),DBG_FH0(IROW),DBG_FH1(IROW),DBG_XH0(IROW),DBG_XH1(IROW),&
                                         ZM%X_H(IROW),HP(I,J,K),RHOP(I,J,K),KRES(I,J,K)
         ENDDO
      ENDDO
   ENDDO
   CLOSE(DBG_LU)
ENDIF
DEALLOCATE(DBG_FH0,DBG_FH1,DBG_XH0,DBG_XH1)

T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW
""")

# ---------------- GLMAT ----------------
rep("""INTEGER :: IERR

! INTEGER  :: JCOL
! REAL(EB) :: LHS""", """INTEGER :: IERR
INTEGER, SAVE :: DBG_NCALL=0
INTEGER :: DBG_LU,DBG_STAT
CHARACTER(80) :: DBG_FN,DBG_ENV
REAL(EB) :: DBG_OFF,DBG_VOL
REAL(EB), ALLOCATABLE, DIMENSION(:) :: DBG_FH0,DBG_FH1,DBG_XH0,DBG_XH1

! INTEGER  :: JCOL
! REAL(EB) :: LHS""")
rep("""IPZ_LOOP : DO IPZ=0,N_ZONE_GLOBMAT

   ZSL => ZONE_SOLVE(IPZ)

   IF (ZSL%NUNKH_TOTAL==0) CYCLE
""", """DBG_NCALL = DBG_NCALL + 1
DBG_OFF = 0._EB
CALL GET_ENVIRONMENT_VARIABLE('FDS_DBG_RHS_OFFSET',DBG_ENV,STATUS=DBG_STAT)
IF (DBG_STAT==0) READ(DBG_ENV,*) DBG_OFF

IPZ_LOOP : DO IPZ=0,N_ZONE_GLOBMAT

   ZSL => ZONE_SOLVE(IPZ)

   IF (ZSL%NUNKH_TOTAL==0) CYCLE
""")
rep("""      CALL GET_FH_FROM_PRHS_AND_BCS(NM,DT,CYL_FCT,UNKH,ZSL%NUNKH_LOCAL,IPZ,ZSL%F_H)
   ENDDO MESH_LOOP_1
""", """      CALL GET_FH_FROM_PRHS_AND_BCS(NM,DT,CYL_FCT,UNKH,ZSL%NUNKH_LOCAL,IPZ,ZSL%F_H)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (CCVAR(I,J,K,UNKH)<=0 .OR. ZONE_SOLVE(PRESSURE_ZONE(I,J,K))%CONNECTED_ZONE_PARENT/=IPZ) CYCLE
               IROW = CCVAR(I,J,K,UNKH) - ZSL%UNKH_IND(NM_START)
               DBG_VOL = ((1._EB-CYL_FCT)*DY(J) + CYL_FCT*RC(I))*DX(I)*DZ(K)
               ZSL%F_H(IROW) = ZSL%F_H(IROW) + DBG_OFF*DBG_VOL
            ENDDO
         ENDDO
      ENDDO
   ENDDO MESH_LOOP_1
   IF (ALLOCATED(DBG_FH0)) DEALLOCATE(DBG_FH0,DBG_FH1,DBG_XH0,DBG_XH1)
   ALLOCATE(DBG_FH0(ZSL%NUNKH_LOCAL),DBG_FH1(ZSL%NUNKH_LOCAL),DBG_XH0(ZSL%NUNKH_LOCAL),DBG_XH1(ZSL%NUNKH_LOCAL))
   DBG_FH0 = ZSL%F_H(1:ZSL%NUNKH_LOCAL)
""")
rep("""   ! WRITE(LU_ERR,*) 'SUM_FH=',SUM(F_H),H_MATRIX_INDEFINITE
""", """   DBG_FH1 = ZSL%F_H(1:ZSL%NUNKH_LOCAL)
   ! WRITE(LU_ERR,*) 'SUM_FH=',SUM(F_H),H_MATRIX_INDEFINITE
""")
rep("""   DEALLOCATE(DISPL)

   IF (ZSL%MTYPE==SYMM_INDEFINITE) THEN
      IF ((.NOT.PRES_ON_WHOLE_DOMAIN""", """   DEALLOCATE(DISPL)
   DBG_XH0 = ZSL%X_H(1:ZSL%NUNKH_LOCAL)

   IF (ZSL%MTYPE==SYMM_INDEFINITE) THEN
      IF ((.NOT.PRES_ON_WHOLE_DOMAIN""")
rep("""   ! WRITE(LU_ERR,*) 'SUM_XH=',SUM(X_H),SUM(A_H(1:IA_H(NUNKH_LOCAL+1)))

   ! Dump result back to mesh containers:""", """   DBG_XH1 = ZSL%X_H(1:ZSL%NUNKH_LOCAL)
   ! WRITE(LU_ERR,*) 'SUM_XH=',SUM(X_H),SUM(A_H(1:IA_H(NUNKH_LOCAL+1)))

   ! Dump result back to mesh containers:""")
# write after mesh loop 2
rep("""      ENDDO WALL_CELL_LOOP_2
   ENDDO MESH_LOOP_2
ENDDO IPZ_LOOP
""", """      ENDDO WALL_CELL_LOOP_2
   ENDDO MESH_LOOP_2
   IF (DBG_NCALL<=4 .AND. IPZ==0) THEN
      NM = LOWER_MESH_INDEX
      CALL POINT_TO_MESH(NM)
      IF (PREDICTOR) THEN; HP => H; RHOP => RHO; ELSE; HP => HS; RHOP => RHOS; ENDIF
      WRITE(DBG_FN,'(A,I0,A)') 'dbg_GLMAT_',DBG_NCALL,'.txt'
      OPEN(NEWUNIT=DBG_LU,FILE=TRIM(DBG_FN),STATUS='REPLACE',ACTION='WRITE')
      WRITE(DBG_LU,'(A,3I6,I3,I8)') 'DIMS',IBAR,JBAR,KBAR,MERGE(1,0,PREDICTOR),ZSL%NUNKH_LOCAL
      WRITE(DBG_LU,'(A,I3,ES25.16)') 'MTYPE',ZSL%MTYPE,DBG_OFF
      WRITE(DBG_LU,'(A,*(ES25.16))') 'DX ',DX(1:IBAR)
      WRITE(DBG_LU,'(A,*(ES25.16))') 'DY ',DY(1:JBAR)
      WRITE(DBG_LU,'(A,*(ES25.16))') 'DZ ',DZ(1:KBAR)
      WRITE(DBG_LU,'(A,*(ES25.16))') 'DXN',DXN(0:IBAR)
      WRITE(DBG_LU,'(A,*(ES25.16))') 'DYN',DYN(0:JBAR)
      WRITE(DBG_LU,'(A,*(ES25.16))') 'DZN',DZN(0:KBAR)
      WRITE(DBG_LU,'(A)') 'COLS I J K PRHS FH0 FH1 XH0 XH1 XHFIN HP RHOP KRES'
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (CCVAR(I,J,K,UNKH)<=0) CYCLE
               IROW = CCVAR(I,J,K,UNKH) - ZSL%UNKH_IND(NM_START)
               WRITE(DBG_LU,'(3I6,9ES25.16)') I,J,K,PRHS(I,J,K),DBG_FH0(IROW),DBG_FH1(IROW),DBG_XH0(IROW),DBG_XH1(IROW),&
                                            ZSL%X_H(IROW),HP(I,J,K),RHOP(I,J,K),KRES(I,J,K)
            ENDDO
         ENDDO
      ENDDO
      CLOSE(DBG_LU)
   ENDIF
ENDDO IPZ_LOOP
""")

# ---------------- force ULMAT instead of FFT on a single regular mesh ----------------
rep("""INTEGER :: I,J,K,IPZ,ICC,JCC,IW,IOR,ZBTYPE_LAST(-3:3),WALL_BTYPE,NZIM,IPZIM,IZERO,JDIM,IERR
""", """INTEGER :: I,J,K,IPZ,ICC,JCC,IW,IOR,ZBTYPE_LAST(-3:3),WALL_BTYPE,NZIM,IPZIM,IZERO,JDIM,IERR
CHARACTER(10) :: DBG_FENV
INTEGER :: DBG_FSTAT
""")
rep("""   IF (M%N_INTERNAL_CFACE_CELLS>0) ZM%USE_FFT=.FALSE.
""", """   IF (M%N_INTERNAL_CFACE_CELLS>0) ZM%USE_FFT=.FALSE.
   CALL GET_ENVIRONMENT_VARIABLE('FDS_DBG_FORCE_ULMAT',DBG_FENV,STATUS=DBG_FSTAT)
   IF (DBG_FSTAT==0) ZM%USE_FFT=.FALSE.
""")

# ---------------- optional scratch fix of the single-rank whole-domain GLMAT mean (SUM_FH(2) never set) ----------------
rep("""         IF (N_MPI_PROCESSES>1) CALL MPI_ALLREDUCE(SUM_FH(1),SUM_FH(2),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
""", """         IF (N_MPI_PROCESSES>1) CALL MPI_ALLREDUCE(SUM_FH(1),SUM_FH(2),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
         IF (N_MPI_PROCESSES==1) THEN   ! scratch: FDS_DBG_FIXMEAN set -> behave as with several ranks
            CALL GET_ENVIRONMENT_VARIABLE('FDS_DBG_FIXMEAN',DBG_ENV,STATUS=DBG_STAT)
            IF (DBG_STAT==0) SUM_FH(2)=SUM_FH(1)
         ENDIF
""")

open(p, "w").write(s)
print("hook applied")
