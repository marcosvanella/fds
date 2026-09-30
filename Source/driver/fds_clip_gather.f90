!> Layout-independent "gather" form of the FDS routine CHECK_MASS_DENSITY (D-031; mass.f90 lines 775-963 of the reference tree), used between
!> DENSITY_PRE_CLIP and DENSITY_POST_CLIP (fds_density_split.f90). Derived from the P1 prototype (clip_gather.f90) with the symbols renamed
!> FDS_CLIP_*; the FDS source is not edited. The algorithm and every expression of FDS are kept; only the loop structure changes.
!> Kernel-facing rules (M2a): (a) passive scalars: the clip acts on the N_TRACKED_SPECIES first components of ZZ/ZZS only, as FDS does; the
!>     remaining components (passive scalars) are carried by Fields.cpp (ncomp = N_TOTAL_SCALARS) and are clipped by CLIP_PASSIVE_SCALARS in
!>     DENSITY_POST_CLIP. (b) only uniform Cartesian metrics are used: DX, DY, DZ are constant per direction, R(I)/RRN(I) are not used.
!>
!> FDS's algorithm is kept; only the loop structure changes.
!>  * FDS scatters: every clipped cell adds +CONST*SUM_MASS_N/VC(0) to itself and -CONST*MASS_N(d)/VC(d) to each
!>    neighbour (mass.f90:840-846, 916-922), K,J,I ascending. Here every cell GATHERS its own term and the terms of its
!>    6 face neighbours (ghost cells included), in exactly the order the K,J,I scatter delivers them to that cell:
!>    (k-1), (j-1), (i-1), self, (i+1), (j+1), (k+1). Each neighbour's CONST/MASS_N/VC is recomputed with the same
!>    expressions (CLIP_TERMS below), so the result is bitwise the single-mesh FDS value.
!>  * The WALL_INDEX barrier (mass.f90:831-836, 907-912) is replaced by an explicit face mask MASK(:,:,:,1:6)
!>    (faces -1,+1,-2,+2,-3,+3): the WALL_INDEX single-mesh FDS computes for the whole level, with faces between boxes
!>    of the same level zero (they are not walls in single-mesh FDS). MASK(:,:,:,0) = SOLID, MASK(:,:,:,7) = 1 if the
!>    cell is a clip SOURCE in single-mesh FDS, i.e. a real cell of this level (not a periodic image outside the domain,
!>    not a coarse-interpolated fine ghost). Only sources contribute to their neighbours.
!>  * The box-global early return / renormalisation (mass.f90:927, 943, 947-961) and the density apply (853-854) are
!>    split into separate passes gated by flags OR-reduced over all boxes and ranks by the driver.
!> RHOP may have a different ghost width than RHO_ZZ and MASK (RLO:RHI): the native RHO/RHOS FAB has 3 ghost layers, ZZ/ZZS and the mask 2.
!> Explicit-shape arguments only; array bounds ALO:AHI are GLOBAL AMReX cell indices of the FAB (incl. 2 ghosts),
!> loop bounds LLO:LHI the valid box. No module pointers, no globals.
MODULE FDS_CLIP_GATHER

USE ISO_C_BINDING
USE PRECISION_PARAMETERS, ONLY: EB, TWO_EPSILON_EB

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE
PUBLIC :: FDS_CLIP_DENSITY, FDS_CLIP_DENSITY_APPLY, FDS_CLIP_SPECIES, FDS_CLIP_SPECIES_APPLY, FDS_CLIP_RENORM, &
          FDS_CLIP_SPECIES_ONE, FDS_CLIP_SUMS

CONTAINS

!> CONST, MASS_N, VC of cell (I,J,K) as FDS computes them for the clip of Q (Q=RHOP for density with limits
!> QMIN/QMAX; Q=RHO_ZZ(:,:,:,N) for species with QMIN=RHO_ZZ_MIN=0, QMAX=RHOP(I,J,K)). Same expressions and
!> evaluation order as mass.f90:801-807,822-839 (density) and 875-881,887-915 (species). ACTIVE=.FALSE. where FDS
!> CYCLEs (in range, SOLID, SUM_MASS_N<=TWO_EPSILON_EB). CLIPPED=.TRUE. where FDS sets the CLIP flag.
!> Q and DX/DY/DZ have bounds QLO:QHI, MASK has MLO:MHI (they may differ: the redundant variant uses ng=3 for RHOP).
PURE SUBROUTINE CLIP_TERMS(QLO,QHI,MLO,MHI,I,J,K,Q,QMIN,QMAX,MASK,DX,DY,DZ,SPECIES,ACTIVE,CLIPPED,CONST,SUM_MASS_N,MASS_N,VC)
INTEGER, INTENT(IN) :: QLO(3),QHI(3),MLO(3),MHI(3),I,J,K
REAL(EB), INTENT(IN) :: Q(QLO(1):QHI(1),QLO(2):QHI(2),QLO(3):QHI(3)),QMIN,QMAX
INTEGER(C_INT), INTENT(IN) :: MASK(MLO(1):MHI(1),MLO(2):MHI(2),MLO(3):MHI(3),0:7)
REAL(EB), INTENT(IN) :: DX(QLO(1):QHI(1)),DY(QLO(2):QHI(2)),DZ(QLO(3):QHI(3))
LOGICAL, INTENT(IN) :: SPECIES
LOGICAL, INTENT(OUT) :: ACTIVE,CLIPPED
REAL(EB), INTENT(OUT) :: CONST,SUM_MASS_N,MASS_N(-3:3),VC(-3:3)
REAL(EB) :: VC1(-3:3),Q_CUT,SIGN_FACTOR,MASS_C
ACTIVE = .FALSE. ; CLIPPED = .FALSE. ; CONST = 0._EB ; SUM_MASS_N = 0._EB ; MASS_N = 0._EB ; VC = 0._EB
IF (SPECIES) THEN                                                   ! mass.f90:885,888 order: SOLID first
   IF (MASK(I,J,K,0)/=0) RETURN
   IF (Q(I,J,K)>=QMIN .AND. Q(I,J,K)<=QMAX) RETURN
ELSE                                                                ! mass.f90:809,811 order: range first
   IF (Q(I,J,K)>=QMIN .AND. Q(I,J,K)<=QMAX) RETURN
   IF (MASK(I,J,K,0)/=0) RETURN
ENDIF
CLIPPED = .TRUE.
IF (Q(I,J,K)<QMIN) THEN
   Q_CUT = QMIN ; SIGN_FACTOR = 1._EB
ELSE
   Q_CUT = QMAX ; SIGN_FACTOR = -1._EB
ENDIF
VC1( 0)  = DY(J)  *DZ(K)
VC1(-1)  = VC1( 0)
VC1( 1)  = VC1( 0)
VC1(-2)  = DY(J-1)*DZ(K)
VC1( 2)  = DY(J+1)*DZ(K)
VC1(-3)  = DY(J)  *DZ(K-1)
VC1( 3)  = DY(J)  *DZ(K+1)
VC( 0)  = DX(I)  * VC1( 0)
VC(-1)  = DX(I-1)* VC1(-1)
VC( 1)  = DX(I+1)* VC1( 1)
VC(-2)  = DX(I)  * VC1(-2)
VC( 2)  = DX(I)  * VC1( 2)
VC(-3)  = DX(I)  * VC1(-3)
VC( 3)  = DX(I)  * VC1( 3)
MASS_C = ABS(Q_CUT-Q(I,J,K))*VC(0)
IF (MASK(I,J,K,1)==0) MASS_N(-1) = ABS(MIN(QMAX,MAX(QMIN,Q(I-1,J,K)))-Q_CUT)*VC(-1)
IF (MASK(I,J,K,2)==0) MASS_N( 1) = ABS(MIN(QMAX,MAX(QMIN,Q(I+1,J,K)))-Q_CUT)*VC( 1)
IF (MASK(I,J,K,3)==0) MASS_N(-2) = ABS(MIN(QMAX,MAX(QMIN,Q(I,J-1,K)))-Q_CUT)*VC(-2)
IF (MASK(I,J,K,4)==0) MASS_N( 2) = ABS(MIN(QMAX,MAX(QMIN,Q(I,J+1,K)))-Q_CUT)*VC( 2)
IF (MASK(I,J,K,5)==0) MASS_N(-3) = ABS(MIN(QMAX,MAX(QMIN,Q(I,J,K-1)))-Q_CUT)*VC(-3)
IF (MASK(I,J,K,6)==0) MASS_N( 3) = ABS(MIN(QMAX,MAX(QMIN,Q(I,J,K+1)))-Q_CUT)*VC( 3)
SUM_MASS_N = SUM(MASS_N)
IF (SUM_MASS_N<=TWO_EPSILON_EB) RETURN
CONST = SIGN_FACTOR*MIN(1._EB,MASS_C/SUM_MASS_N)
ACTIVE = .TRUE.
END SUBROUTINE CLIP_TERMS


!> Gather the correction DELTA(I,J,K) for every cell of the loop region LLO:LHI (the valid box, or valid+1 for the
!> redundant density variant). SPECIES: Q=RHO_ZZ(:,N), QMIN=0, QMAX=RHOP(cell); otherwise Q=RHOP, QMIN=RHOMIN,
!> QMAX=RHOMAX. NCLIP/NLO/NHI count only cells of the VALID box VLO:VHI (only valid cells may set clip flags).
SUBROUTINE GATHER(QLO,QHI,RLO,RHI,MLO,MHI,LLO,LHI,DLO,DHI,VLO,VHI,Q,RHOP,QMIN_IN,QMAX_IN,MASK,DX,DY,DZ,SPECIES,DELTA,NCLIP,NLO,NHI)
INTEGER, INTENT(IN) :: QLO(3),QHI(3),RLO(3),RHI(3),MLO(3),MHI(3),LLO(3),LHI(3),DLO(3),DHI(3),VLO(3),VHI(3)
REAL(EB), INTENT(IN) :: Q(QLO(1):QHI(1),QLO(2):QHI(2),QLO(3):QHI(3)),RHOP(RLO(1):RHI(1),RLO(2):RHI(2),RLO(3):RHI(3))
REAL(EB), INTENT(IN) :: QMIN_IN,QMAX_IN
INTEGER(C_INT), INTENT(IN) :: MASK(MLO(1):MHI(1),MLO(2):MHI(2),MLO(3):MHI(3),0:7)
REAL(EB), INTENT(IN) :: DX(QLO(1):QHI(1)),DY(QLO(2):QHI(2)),DZ(QLO(3):QHI(3))
LOGICAL, INTENT(IN) :: SPECIES
REAL(EB), INTENT(INOUT) :: DELTA(DLO(1):DHI(1),DLO(2):DHI(2),DLO(3):DHI(3))
INTEGER, INTENT(OUT) :: NCLIP,NLO,NHI
INTEGER :: I,J,K,S,II,JJ,KK,D
REAL(EB) :: CONST,SUM_MASS_N,MASS_N(-3:3),VC(-3:3),QMAX,DEL
LOGICAL :: ACTIVE,CLIPPED,VALID
! neighbour offsets in FDS K,J,I scatter-delivery order; D = direction index of the TARGET cell as seen from the source
INTEGER, PARAMETER :: OI(7)=[0,0,-1,0,1,0,0], OJ(7)=[0,-1,0,0,0,1,0], OK(7)=[-1,0,0,0,0,0,1], OD(7)=[3,2,1,0,-1,-2,-3]
NCLIP = 0 ; NLO = 0 ; NHI = 0
DO K=LLO(3),LHI(3)
   DO J=LLO(2),LHI(2)
      DO I=LLO(1),LHI(1)
         VALID = I>=VLO(1) .AND. I<=VHI(1) .AND. J>=VLO(2) .AND. J<=VHI(2) .AND. K>=VLO(3) .AND. K<=VHI(3)
         DEL = 0._EB                                              ! DELTA_RHO(_ZZ) = 0 (mass.f90:785, 871)
         DO S=1,7
            II = I+OI(S) ; JJ = J+OJ(S) ; KK = K+OK(S) ; D = OD(S)
            IF (MASK(II,JJ,KK,7)==0) CYCLE                        ! not a cell of the single-mesh loop 1:IBAR
            QMAX = QMAX_IN ; IF (SPECIES) QMAX = RHOP(II,JJ,KK)
            CALL CLIP_TERMS(QLO,QHI,MLO,MHI,II,JJ,KK,Q,QMIN_IN,QMAX,MASK,DX,DY,DZ,SPECIES,ACTIVE,CLIPPED,CONST,SUM_MASS_N,MASS_N,VC)
            IF (S==4 .AND. CLIPPED .AND. VALID) THEN
               NCLIP = NCLIP+1
               IF (Q(I,J,K)<QMIN_IN) THEN ; NLO = NLO+1 ; ELSE ; NHI = NHI+1 ; ENDIF
            ENDIF
            IF (.NOT.ACTIVE) CYCLE
            IF (D==0) THEN
               DEL = DEL + CONST*SUM_MASS_N/VC( 0)                 ! mass.f90:840 / 916
            ELSE
               DEL = DEL - CONST*MASS_N(D)/VC(D)                   ! mass.f90:841-846 / 917-922
            ENDIF
         ENDDO
         DELTA(I,J,K) = DEL
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE GATHER


!> Pass 1 (density, mass.f90:799-849): DRHO = gathered DELTA_RHO on the loop region LLO:LHI (valid box for the
!> 2nd-FillBoundary variant, valid+1 for the redundant variant); FLAGS(1:2) = local CLIP_RHOMIN, CLIP_RHOMAX
!> (mass.f90:815,819) from VALID cells VLO:VHI only; NCL = clipped valid cells. RHOP/DX: bounds QLO:QHI, MASK: MLO:MHI,
!> DRHO: DLO:DHI.
SUBROUTINE FDS_CLIP_DENSITY(QLO,QHI,MLO,MHI,LLO,LHI,DLO,DHI,VLO,VHI,RHOP,MASK,DX,DY,DZ,RHO_MIN,RHO_MAX,DRHO,FLAGS,NCL) &
           BIND(C,NAME='fds_clip_density')
INTEGER(C_INT), INTENT(IN) :: QLO(3),QHI(3),MLO(3),MHI(3),LLO(3),LHI(3),DLO(3),DHI(3),VLO(3),VHI(3)
REAL(C_DOUBLE), INTENT(IN) :: RHOP(QLO(1):QHI(1),QLO(2):QHI(2),QLO(3):QHI(3))
INTEGER(C_INT), INTENT(IN) :: MASK(MLO(1):MHI(1),MLO(2):MHI(2),MLO(3):MHI(3),0:7)
REAL(C_DOUBLE), INTENT(IN) :: DX(QLO(1):QHI(1)),DY(QLO(2):QHI(2)),DZ(QLO(3):QHI(3))
REAL(C_DOUBLE), VALUE :: RHO_MIN,RHO_MAX
REAL(C_DOUBLE), INTENT(INOUT) :: DRHO(DLO(1):DHI(1),DLO(2):DHI(2),DLO(3):DHI(3))
INTEGER(C_INT), INTENT(OUT) :: FLAGS(2),NCL
INTEGER :: NC,NLO,NHI
CALL GATHER(QLO,QHI,QLO,QHI,MLO,MHI,LLO,LHI,DLO,DHI,VLO,VHI,RHOP,RHOP,RHO_MIN,RHO_MAX,MASK,DX,DY,DZ,.FALSE.,DRHO,NC,NLO,NHI)
FLAGS(1) = MERGE(1,0,NLO>0) ; FLAGS(2) = MERGE(1,0,NHI>0) ; NCL = NC
END SUBROUTINE FDS_CLIP_DENSITY


!> Density apply (mass.f90:853-854) over LLO:LHI (valid, or valid+1 in the redundant variant); the driver calls it only
!> if the REDUCED CLIP_RHOMIN .OR. CLIP_RHOMAX is set.
SUBROUTINE FDS_CLIP_DENSITY_APPLY(QLO,QHI,LLO,LHI,DLO,DHI,RHOP,DRHO,RHO_MIN,RHO_MAX) BIND(C,NAME='fds_clip_density_apply')
INTEGER(C_INT), INTENT(IN) :: QLO(3),QHI(3),LLO(3),LHI(3),DLO(3),DHI(3)
REAL(C_DOUBLE), INTENT(INOUT) :: RHOP(QLO(1):QHI(1),QLO(2):QHI(2),QLO(3):QHI(3))
REAL(C_DOUBLE), INTENT(IN) :: DRHO(DLO(1):DHI(1),DLO(2):DHI(2),DLO(3):DHI(3))
REAL(C_DOUBLE), VALUE :: RHO_MIN,RHO_MAX
RHOP(LLO(1):LHI(1),LLO(2):LHI(2),LLO(3):LHI(3)) = &
   MIN(RHO_MAX,MAX(RHO_MIN,RHOP(LLO(1):LHI(1),LLO(2):LHI(2),LLO(3):LHI(3))+DRHO(LLO(1):LHI(1),LLO(2):LHI(2),LLO(3):LHI(3))))
END SUBROUTINE FDS_CLIP_DENSITY_APPLY


!> N_TRACKED_SPECIES==1 (mass.f90:858-861), gated by the reduced density flag.
SUBROUTINE FDS_CLIP_SPECIES_ONE(RLO,RHI,ALO,AHI,LLO,LHI,RHOP,RHO_ZZ) BIND(C,NAME='fds_clip_species_one')
INTEGER(C_INT), INTENT(IN) :: RLO(3),RHI(3),ALO(3),AHI(3),LLO(3),LHI(3)
REAL(C_DOUBLE), INTENT(IN) :: RHOP(RLO(1):RHI(1),RLO(2):RHI(2),RLO(3):RHI(3))
REAL(C_DOUBLE), INTENT(INOUT) :: RHO_ZZ(ALO(1):AHI(1),ALO(2):AHI(2),ALO(3):AHI(3),1)
RHO_ZZ(LLO(1):LHI(1),LLO(2):LHI(2),LLO(3):LHI(3),1) = RHOP(LLO(1):LHI(1),LLO(2):LHI(2),LLO(3):LHI(3))
END SUBROUTINE FDS_CLIP_SPECIES_ONE


!> Pass 2 (species N, mass.f90:870-925) on the valid box LLO:LHI: DZZ = gathered DELTA_RHO_ZZ(:,N); FLAG = local
!> CLIP_RHO_ZZ(N) (mass.f90:889) from valid cells. RHOP (bounds RLO:RHI) must be the density-clipped RHOP, current in
!> the first ghost layer. RHO_ZZ/DX: bounds QLO:QHI (2 ghosts); MASK: MLO:MHI.
SUBROUTINE FDS_CLIP_SPECIES(QLO,QHI,RLO,RHI,MLO,MHI,LLO,LHI,NS,N,RHOP,RHO_ZZ,MASK,DX,DY,DZ,DZZ,FLAG,NCL) &
           BIND(C,NAME='fds_clip_species')
INTEGER(C_INT), INTENT(IN) :: QLO(3),QHI(3),RLO(3),RHI(3),MLO(3),MHI(3),LLO(3),LHI(3)
INTEGER(C_INT), VALUE :: NS,N
REAL(C_DOUBLE), INTENT(IN) :: RHOP(RLO(1):RHI(1),RLO(2):RHI(2),RLO(3):RHI(3))
REAL(C_DOUBLE), INTENT(IN) :: RHO_ZZ(QLO(1):QHI(1),QLO(2):QHI(2),QLO(3):QHI(3),NS)
INTEGER(C_INT), INTENT(IN) :: MASK(MLO(1):MHI(1),MLO(2):MHI(2),MLO(3):MHI(3),0:7)
REAL(C_DOUBLE), INTENT(IN) :: DX(QLO(1):QHI(1)),DY(QLO(2):QHI(2)),DZ(QLO(3):QHI(3))
REAL(C_DOUBLE), INTENT(INOUT) :: DZZ(LLO(1):LHI(1),LLO(2):LHI(2),LLO(3):LHI(3))
INTEGER(C_INT), INTENT(OUT) :: FLAG,NCL
INTEGER :: NC,NLO,NHI
CALL GATHER(QLO,QHI,RLO,RHI,MLO,MHI,LLO,LHI,LLO,LHI,LLO,LHI,RHO_ZZ(:,:,:,N),RHOP,0._EB,0._EB,MASK,DX,DY,DZ,.TRUE.,DZZ,NC,NLO,NHI)
FLAG = MERGE(1,0,NC>0) ; NCL = NC
END SUBROUTINE FDS_CLIP_SPECIES


!> Species apply (mass.f90:927-937); the driver calls it only if the REDUCED CLIP_RHO_ZZ(N) is set.
SUBROUTINE FDS_CLIP_SPECIES_APPLY(RLO,RHI,ALO,AHI,LLO,LHI,NS,N,RHOP,RHO_ZZ,DZZ) BIND(C,NAME='fds_clip_species_apply')
INTEGER(C_INT), INTENT(IN) :: RLO(3),RHI(3),ALO(3),AHI(3),LLO(3),LHI(3)
INTEGER(C_INT), VALUE :: NS,N
REAL(C_DOUBLE), INTENT(IN) :: RHOP(RLO(1):RHI(1),RLO(2):RHI(2),RLO(3):RHI(3))
REAL(C_DOUBLE), INTENT(INOUT) :: RHO_ZZ(ALO(1):AHI(1),ALO(2):AHI(2),ALO(3):AHI(3),NS)
REAL(C_DOUBLE), INTENT(IN) :: DZZ(LLO(1):LHI(1),LLO(2):LHI(2),LLO(3):LHI(3))
REAL(EB), PARAMETER :: RHO_ZZ_MIN = 0._EB
INTEGER :: I,J,K
DO K=LLO(3),LHI(3)
   DO J=LLO(2),LHI(2)
      DO I=LLO(1),LHI(1)
         RHO_ZZ(I,J,K,N) = MIN(RHOP(I,J,K),MAX(RHO_ZZ_MIN,RHO_ZZ(I,J,K,N)+DZZ(I,J,K)))
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE FDS_CLIP_SPECIES_APPLY


!> Renormalisation (mass.f90:947-961); the driver calls it only if the REDUCED
!> CLIP_RHOMIN .OR. CLIP_RHOMAX .OR. ANY(CLIP_RHO_ZZ) is set (replaces the per-box early return at mass.f90:943).
SUBROUTINE FDS_CLIP_RENORM(RLO,RHI,ALO,AHI,LLO,LHI,NS,NT,RHOP,RHO_ZZ,MASK) BIND(C,NAME='fds_clip_renorm')
INTEGER(C_INT), INTENT(IN) :: RLO(3),RHI(3),ALO(3),AHI(3),LLO(3),LHI(3)
INTEGER(C_INT), VALUE :: NS,NT   ! NS = extent of the species dimension of RHO_ZZ, NT = N_TRACKED_SPECIES (the passive scalars are not renormalised)
REAL(C_DOUBLE), INTENT(IN) :: RHOP(RLO(1):RHI(1),RLO(2):RHI(2),RLO(3):RHI(3))
REAL(C_DOUBLE), INTENT(INOUT) :: RHO_ZZ(ALO(1):AHI(1),ALO(2):AHI(2),ALO(3):AHI(3),NS)
INTEGER(C_INT), INTENT(IN) :: MASK(ALO(1):AHI(1),ALO(2):AHI(2),ALO(3):AHI(3),0:7)
INTEGER :: I,J,K,N
REAL(EB) :: SUM_RHO_ZZ,RHO_ZZ_TEST
DO K=LLO(3),LHI(3)
   DO J=LLO(2),LHI(2)
      DO I=LLO(1),LHI(1)
         IF (MASK(I,J,K,0)/=0) CYCLE
         SUM_RHO_ZZ = SUM(RHO_ZZ(I,J,K,1:NT))
         N = MAXLOC(RHO_ZZ(I,J,K,1:NT),1)
         RHO_ZZ_TEST = RHO_ZZ(I,J,K,N) + RHOP(I,J,K) - SUM_RHO_ZZ
         IF (RHO_ZZ_TEST<0._EB .OR. RHO_ZZ_TEST>RHOP(I,J,K)) THEN  ! Renormalize the original set of RHO_ZZ
            RHO_ZZ(I,J,K,1:NT) = RHOP(I,J,K) * RHO_ZZ(I,J,K,1:NT)/SUM_RHO_ZZ
         ELSE  ! Absorb mass deficit/excess into largest RHO_ZZ
            RHO_ZZ(I,J,K,N) = RHO_ZZ_TEST
         ENDIF
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE FDS_CLIP_RENORM


!> Diagnostic: per-box sum over the valid box of RHO_ZZ(:,N)*DX*DY*DZ for each species (mass per species).
SUBROUTINE FDS_CLIP_SUMS(ALO,AHI,LLO,LHI,NS,RHO_ZZ,DX,DY,DZ,MS) BIND(C,NAME='fds_clip_sums')
INTEGER(C_INT), INTENT(IN) :: ALO(3),AHI(3),LLO(3),LHI(3)
INTEGER(C_INT), VALUE :: NS
REAL(C_DOUBLE), INTENT(IN) :: RHO_ZZ(ALO(1):AHI(1),ALO(2):AHI(2),ALO(3):AHI(3),NS)
REAL(C_DOUBLE), INTENT(IN) :: DX(ALO(1):AHI(1)),DY(ALO(2):AHI(2)),DZ(ALO(3):AHI(3))
REAL(C_DOUBLE), INTENT(INOUT) :: MS(NS)
INTEGER :: I,J,K,N
DO N=1,NS ; DO K=LLO(3),LHI(3) ; DO J=LLO(2),LHI(2) ; DO I=LLO(1),LHI(1)
   MS(N) = MS(N) + RHO_ZZ(I,J,K,N)*DX(I)*DY(J)*DZ(K)
ENDDO ; ENDDO ; ENDDO ; ENDDO
END SUBROUTINE FDS_CLIP_SUMS

END MODULE FDS_CLIP_GATHER
