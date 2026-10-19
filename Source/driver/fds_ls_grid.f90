!> \brief Level set on one global 2-D grid (the finest-level resolution), advanced independently of the gas-phase meshes.
!>
!> The numerics are a copy of LEVEL_SET_FIRESPREAD / LEVEL_SET_SPREAD_RATE / LEVEL_SET_ADVECT_FLUX / IGNITE_GRID_CELL of Source/vege.f90 (flat terrain, one
!> grid, no CC_IBM) with the mesh pointers replaced by the arrays of this module and the wall state (T_IGN, TAU_LS, BURN_DURATION, AREA_ADJUST) kept per LS cell.
!> The arithmetic order is the one of vege.f90, so a grid equal to the level-0 mesh reproduces the original bit for bit (shadow test, FDS_LS_SHADOW_*).
!> Every rank holds the whole grid and advances it with the same inputs (the wind at the reference height, LS_REF_WIND), so no halo exchange between ranks is needed.
MODULE FDS_LS_GRID

USE ISO_C_BINDING
USE PRECISION_PARAMETERS, ONLY: EB,TWENTY_EPSILON_EB,PI,PIO2
USE GLOBAL_CONSTANTS, ONLY: PREDICTOR,LEVEL_SET_MODE,LEVEL_SET_COUPLED_WIND,LEVEL_SET_COUPLED_FIRE,I_FLUX_LIMITER,CHID,NMESHES,SOLID_BOUNDARY,STORE_FIRE_ARRIVAL,STORE_FIRE_RESIDENCE,N_TRACKED_SPECIES
USE TYPES, ONLY: SURFACE,SURFACE_TYPE
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE PHYSICAL_FUNCTIONS, ONLY: GET_WIND_AT_HEIGHT
USE MESH_VARIABLES, ONLY: MESHES,MESH_TYPE

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE
PUBLIC :: LS_DRIVE,LS_APPLY_L0,LS_FLUX_ACTIVE,LS_BLOCK_FLUX,LS_ACTIVE,LS_INIT,LS_STEP,LS_REF_WIND,LS_REF_HEIGHT,LS_NX,LS_NY,LS_PHI,LS_PHI1,LS_T_IGN,LS_TAU,LS_BURN_DUR,LS_AREA_ADJ,LS_SURF,LS_DX,LS_DY, &
   LS_UREF,LS_VREF,LS_U,LS_V,LS_SHADOW_MODE,FDS_LS_SET_MODE,FDS_LS_SHADOW_RUN,FDS_LS_SHADOW_CMP,FDS_LS_PREPARE,FDS_LS_WIND_FROM_PAIRS,FDS_LS_REF_HEIGHTS,FDS_LS_DRIVE_STEP,LS_RESTRICT_PHI

REAL(EB), PARAMETER :: PHI_LS_MIN=-1._EB, PHI_LS_MAX=1._EB

LOGICAL :: LS_ACTIVE=.FALSE.          !< the LS grid exists
LOGICAL :: LS_WIND_SET=.FALSE.        !< LS_UREF/LS_VREF were set by FDS_LS_WIND_FROM_PAIRS for the coming stage
INTEGER :: LS_RATIO=1                 !< LS cell = level-0 cell / LS_RATIO
LOGICAL :: LS_DRIVE=.FALSE.           !< --ls-fine: the LS grid drives the fire (the original LEVEL_SET_FIRESPREAD is not called; level-0 walls get the LS cell means)
REAL(EB) :: LS_X0=0._EB,LS_Y0=0._EB      !< lower corner of the LS grid
INTEGER, ALLOCATABLE :: LS_MI0(:),LS_MJ0(:)   !< level-0 cell offset of mesh NM in the LS grid (level-0 cells)
LOGICAL, ALLOCATABLE :: LS_MGROUND(:)        !< mesh NM sits on the ground (lowest ZS) and owns LS cells
REAL(EB), ALLOCATABLE :: LS_FLUX(:,:)    !< fuel mass flux of every LS cell at the current stage, kg/m2/s (LS_APPLY_L0)
REAL(EB) :: LS_FLUX_T = -HUGE(1._EB)      !< stage time (T+DT predictor, T corrector) of LS_FLUX
LOGICAL :: LS_FLUX_PRED = .FALSE.
INTEGER :: LS_SHADOW_MODE=0           !< --ls-fine-shadow: run beside the original LEVEL_SET_FIRESPREAD and compare (level-0 equal grid)
INTEGER :: LS_NX=0,LS_NY=0,LIMITER_LS=1
REAL(EB), ALLOCATABLE :: LS_DX(:),LS_DY(:),RDXL(:),RDYL(:)
REAL(EB), ALLOCATABLE, TARGET :: LS_PHI(:,:),LS_PHI1(:,:),FLUX0(:,:),FLUX1(:,:),SR_X(:,:),SR_Y(:,:),WORK1(:,:),WORK2(:,:)
REAL(EB), ALLOCATABLE :: PHI_WS(:,:),THETA_ELPS(:,:),UMF(:,:),ROS_HEAD(:,:)
REAL(EB), ALLOCATABLE :: LS_UREF(:,:),LS_VREF(:,:)    !< wind at the reference height of the cell (input, raw)
REAL(EB), ALLOCATABLE :: LS_U(:,:),LS_V(:,:)          !< midflame wind (output, as U_LS, V_LS of vege.f90)
REAL(EB), ALLOCATABLE :: LS_T_IGN(:,:),LS_TAU(:,:),LS_BURN_DUR(:,:),LS_AREA_ADJ(:,:)
INTEGER, ALLOCATABLE :: LS_SURF(:,:)

!> shadow-test counters
INTEGER :: SH_STEPS=0,SH_BAD_PHI=0,SH_BAD_WALL=0,SH_CELLS_IGN=0
REAL(EB) :: SH_MAXDIFF=0._EB

CONTAINS


!> \brief Allocate the grid; DX,DY per cell; SURF index per cell; T_IGN etc. initial values (as the wall cells have them at the start)
SUBROUTINE LS_INIT(NX,NY,DX,DY,SURF,T_IGN0,BURN_DUR0,AREA_ADJ0,TAU0)

INTEGER, INTENT(IN) :: NX,NY,SURF(NX,NY)
REAL(EB), INTENT(IN) :: DX(NX),DY(NY),T_IGN0(NX,NY),BURN_DUR0(NX,NY),AREA_ADJ0(NX,NY),TAU0(NX,NY)

LS_NX=NX ; LS_NY=NY
ALLOCATE(LS_DX(NX),LS_DY(NY),RDXL(NX),RDYL(NY))
LS_DX=DX ; LS_DY=DY ; RDXL=1._EB/DX ; RDYL=1._EB/DY
ALLOCATE(LS_PHI(-1:NX+2,-1:NY+2),LS_PHI1(-1:NX+2,-1:NY+2)) ; LS_PHI=PHI_LS_MIN ; LS_PHI1=PHI_LS_MIN
ALLOCATE(FLUX0(NX,NY),FLUX1(NX,NY),SR_X(NX,NY),SR_Y(NX,NY)) ; FLUX0=0._EB ; FLUX1=0._EB ; SR_X=0._EB ; SR_Y=0._EB
ALLOCATE(WORK1(0:NX,0:NY),WORK2(0:NX,0:NY)) ; WORK1=0._EB ; WORK2=0._EB
ALLOCATE(PHI_WS(NX,NY),THETA_ELPS(NX,NY),UMF(NX,NY),ROS_HEAD(NX,NY)) ; PHI_WS=0._EB ; THETA_ELPS=0._EB ; UMF=0._EB ; ROS_HEAD=0._EB
ALLOCATE(LS_UREF(NX,NY),LS_VREF(NX,NY)) ; LS_UREF=0._EB ; LS_VREF=0._EB
ALLOCATE(LS_U(0:NX+1,0:NY+1),LS_V(0:NX+1,0:NY+1)) ; LS_U=0._EB ; LS_V=0._EB
ALLOCATE(LS_T_IGN(NX,NY),LS_TAU(NX,NY),LS_BURN_DUR(NX,NY),LS_AREA_ADJ(NX,NY))
LS_T_IGN=T_IGN0 ; LS_BURN_DUR=BURN_DUR0 ; LS_AREA_ADJ=AREA_ADJ0 ; LS_TAU=TAU0
ALLOCATE(LS_SURF(NX,NY)) ; LS_SURF=SURF
LIMITER_LS=I_FLUX_LIMITER
IF (LIMITER_LS>3) LIMITER_LS=1
LS_ACTIVE=.TRUE.

END SUBROUTINE LS_INIT


!> \brief Reference height of the wind sample of an LS cell (m above the ground), as in LEVEL_SET_FIRESPREAD
REAL(EB) FUNCTION LS_REF_HEIGHT(I,J)
INTEGER, INTENT(IN) :: I,J
LS_REF_HEIGHT = SURFACE(LS_SURF(I,J))%VEG_LSET_WIND_HEIGHT
IF (SURFACE(LS_SURF(I,J))%VEG_LSET_WIND_HEIGHT<0._EB) LS_REF_HEIGHT = 6.1_EB
END FUNCTION LS_REF_HEIGHT


!> \brief Wind at the height HREF above the ground from a column of cell-centre samples (z above ground ascending, u, v already averaged to the cell centre).
!> Same selection and interpolation as the coupled-wind block of LEVEL_SET_FIRESPREAD; the samples may come from cells of different size (leaf cells of the column).
SUBROUTINE LS_REF_WIND(NS,ZS,US,VS,HREF,UO,VO)

INTEGER, INTENT(IN) :: NS
REAL(EB), INTENT(IN) :: ZS(NS),US(NS),VS(NS),HREF
REAL(EB), INTENT(OUT) :: UO,VO
INTEGER :: KDUM,KWIND
REAL(EB) :: ZWIND(2),U_Z(2),V_Z(2)

IF (ZS(NS)<HREF) THEN      ! sample height above the top cell: the top cell
   UO = US(NS) ; VO = VS(NS)
ELSE
   KWIND = 0
   DO KDUM=1,NS
      IF (ZS(KDUM)>=HREF) THEN
         KWIND = KDUM
         ZWIND(2) = ZS(KWIND)
         IF (KWIND>1) ZWIND(1) = ZS(KWIND-1)
         EXIT
      ENDIF
   ENDDO
   U_Z(2) = US(KWIND) ; V_Z(2) = VS(KWIND)
   IF (KWIND==1) THEN   ! wind from the first cell: plug flow near the ground
      U_Z(1) = U_Z(2) ; V_Z(1) = V_Z(2) ; ZWIND(1) = 0._EB
   ELSE
      U_Z(1) = US(KWIND-1) ; V_Z(1) = VS(KWIND-1)
   ENDIF
   UO = U_Z(1) + (HREF-ZWIND(1))/(ZWIND(2)-ZWIND(1))*(U_Z(2)-U_Z(1))
   VO = V_Z(1) + (HREF-ZWIND(1))/(ZWIND(2)-ZWIND(1))*(V_Z(2)-V_Z(1))
ENDIF

END SUBROUTINE LS_REF_WIND


!> \brief One RK stage of the level set (LEVEL_SET_FIRESPREAD of vege.f90): input LS_UREF, LS_VREF (wind at the reference height).
SUBROUTINE LS_STEP(T,DT)

REAL(EB), INTENT(IN) :: T,DT
INTEGER :: IIG,JJG
REAL(EB) :: UMF_TMP,PHI_W_X,PHI_W_Y,UMF_X,UMF_Y,UMF_MAG,WIND_FACTOR,REF_WIND_HEIGHT,WW
TYPE(SURFACE_TYPE), POINTER :: SF

DO JJG=1,LS_NY
   DO IIG=1,LS_NX

      SF => SURFACE(LS_SURF(IIG,JJG))

      IF (.NOT. SF%VEG_LSET_SPREAD) CYCLE

      ! Ignite landscape at user specified location and time

      IF (LS_PHI(IIG,JJG)<0._EB .AND. T>SF%VEG_LSET_IGNITE_T) LS_PHI(IIG,JJG) = PHI_LS_MAX

      REF_WIND_HEIGHT = LS_REF_HEIGHT(IIG,JJG)

      IF (LEVEL_SET_COUPLED_WIND) THEN
         LS_U(IIG,JJG) = LS_UREF(IIG,JJG)
         LS_V(IIG,JJG) = LS_VREF(IIG,JJG)
      ELSE
         CALL GET_WIND_AT_HEIGHT(T,REF_WIND_HEIGHT,LS_U(IIG,JJG),LS_V(IIG,JJG),WW)
      ENDIF

      UMF_TMP = 1._EB

      ! Wind at midflame height (UMF). From Andrews 2012, USDA FS Gen Tech Rep. RMRS-GTR-266 (with added SI conversion)
      IF (SF%VEG_LSET_WIND_HEIGHT<0._EB) &
         UMF_TMP = 1.83_EB / LOG((20.0_EB + 1.18_EB * SF%VEG_LSET_HT) /(0.43_EB * SF%VEG_LSET_HT))  ! Bova et al., Eq. A2

      ! Factor 60 converts U from m/s to m/min which is used in elliptical model.

      UMF_X = UMF_TMP * LS_U(IIG,JJG) * 60.0_EB
      UMF_Y = UMF_TMP * LS_V(IIG,JJG) * 60.0_EB
      UMF_MAG = SQRT(UMF_X**2 + UMF_Y**2)

      ! Adjust wind for output slices
      LS_U(IIG,JJG) = UMF_TMP * LS_U(IIG,JJG)
      LS_V(IIG,JJG) = UMF_TMP * LS_V(IIG,JJG)

      ! Compute wind factor affecting spread rate R(U) = R_0*(1+WIND_FACTOR)

      IF (SF%I_RAMP_LS_WIND>0) THEN
         WIND_FACTOR = EVALUATE_RAMP(UMF_MAG/60._EB,SF%I_RAMP_LS_WIND)
      ELSE
         WIND_FACTOR = SF%C_ROTH * ((3.281_EB * UMF_MAG)**SF%B_ROTH) * SF%BETA_ROTH ! Bova et al., Eq. A1
      ENDIF

      IF (UMF_MAG>TWENTY_EPSILON_EB) THEN
         PHI_W_X = WIND_FACTOR*UMF_X/UMF_MAG
         PHI_W_Y = WIND_FACTOR*UMF_Y/UMF_MAG
      ELSE
         PHI_W_X = 0.0_EB
         PHI_W_Y = 0.0_EB
      ENDIF

      ! Flat terrain: no slope term

      PHI_WS(IIG,JJG) = SQRT(PHI_W_X**2 + PHI_W_Y**2)
      THETA_ELPS(IIG,JJG) = ATAN2(PHI_W_Y,PHI_W_X)

      UMF(IIG,JJG) = SQRT(UMF_X**2 + UMF_Y**2)  ! U in Eq. A6, Bova et al.

      ! The following two lines convert ATAN2 output to compass system (0 to 2 pi CW from +Y-axis)

      THETA_ELPS(IIG,JJG) = PIO2 - THETA_ELPS(IIG,JJG)
      IF (THETA_ELPS(IIG,JJG) < 0.0_EB) THETA_ELPS(IIG,JJG) = 2.0_EB*PI + THETA_ELPS(IIG,JJG)

      IF (SF%VEG_LSET_ROS_00 > 0._EB) ROS_HEAD(IIG,JJG) = SF%VEG_LSET_ROS_00*(1._EB + PHI_WS(IIG,JJG))  ! Bova et al., Eq. A3

   ENDDO
ENDDO

! Runge-Kutta Scheme

CALL SPREAD_RATE
CALL ADVECT_FLUX

IF (PREDICTOR) THEN
   LS_PHI1(1:LS_NX,1:LS_NY) = LS_PHI(1:LS_NX,1:LS_NY) - DT*FLUX0(1:LS_NX,1:LS_NY)
   LS_PHI1 = MAX(PHI_LS_MIN,MIN(PHI_LS_MAX,LS_PHI1))
ELSE
   LS_PHI(1:LS_NX,1:LS_NY) = LS_PHI(1:LS_NX,1:LS_NY) - 0.5_EB*DT*(FLUX0(1:LS_NX,1:LS_NY) + FLUX1(1:LS_NX,1:LS_NY))
   LS_PHI = MAX(PHI_LS_MIN,MIN(PHI_LS_MAX,LS_PHI))
ENDIF

! Ignition time and burn history of newly ignited cells

IF (.NOT.PREDICTOR) THEN
   DO JJG=1,LS_NY
      DO IIG=1,LS_NX
         SF => SURFACE(LS_SURF(IIG,JJG))
         IF (.NOT. SF%VEG_LSET_SPREAD) CYCLE
         IF (LS_PHI(IIG,JJG)>=0._EB .AND. LS_T_IGN(IIG,JJG)>9.E5_EB) CALL IGNITE_GRID_CELL(IIG,JJG,SF)
      ENDDO
   ENDDO
ENDIF

CONTAINS

SUBROUTINE IGNITE_GRID_CELL(I,J,SFC)

INTEGER, INTENT(IN) :: I,J
TYPE(SURFACE_TYPE), POINTER :: SFC
REAL(EB) :: CROSSING_DISTANCE,FIRE_DEPTH,TAU_2,ROS_MAG

LS_T_IGN(I,J) = T
ROS_MAG = MAX(0.0001_EB,SQRT(SR_X(I,J)**2 + SR_Y(I,J)**2))  ! Rate Of Spread magnitude
IF (ABS(SR_X(I,J))<TWENTY_EPSILON_EB) THEN
   CROSSING_DISTANCE = LS_DY(J)
ELSEIF (ABS(SR_Y(I,J))<TWENTY_EPSILON_EB) THEN
   CROSSING_DISTANCE = LS_DX(I)
ELSE
   CROSSING_DISTANCE = SQRT( (MIN(LS_DX(I),ABS(SR_X(I,J)/SR_Y(I,J))*LS_DY(J)))**2 + &
                             (MIN(LS_DY(J),ABS(SR_Y(I,J)/SR_X(I,J))*LS_DX(I)))**2 )
ENDIF
FIRE_DEPTH = SFC%BURN_DURATION*ROS_MAG
IF (FIRE_DEPTH >= CROSSING_DISTANCE) THEN
   LS_TAU(I,J) = CROSSING_DISTANCE/ROS_MAG
   TAU_2 = (FIRE_DEPTH - CROSSING_DISTANCE)/ROS_MAG
ELSE
   LS_TAU(I,J) = FIRE_DEPTH/ROS_MAG
   TAU_2 = (CROSSING_DISTANCE - FIRE_DEPTH)/ROS_MAG
ENDIF

LS_BURN_DUR(I,J) = 2._EB*LS_TAU(I,J) + TAU_2
! Adjust mass flux to conserve total energy output
IF (LEVEL_SET_COUPLED_FIRE) LS_AREA_ADJ(I,J) = LS_AREA_ADJ(I,J) * SFC%BURN_DURATION/(LS_TAU(I,J)+TAU_2)

END SUBROUTINE IGNITE_GRID_CELL

END SUBROUTINE LS_STEP


!> \brief Components of the spread rate vector (LEVEL_SET_SPREAD_RATE of vege.f90, flat terrain)
SUBROUTINE SPREAD_RATE

INTEGER :: I,J,IM1,IP1,JM1,JP1
REAL(EB) :: DPHIDX,DPHIDY,F_EAST,F_WEST,F_NORTH,F_SOUTH,MAG_F
REAL(EB) :: COS_THETA,SIN_THETA,XSF,YSF,UMF_DUM
REAL(EB) :: AROS,A_ELPS,A_ELPS2,BROS,B_ELPS2,B_ELPS,C_ELPS,DENOM,ROS_TMP,LB,LBD,HB
REAL(EB), DIMENSION(:) :: NORMAL_FIRELINE(2)
REAL(EB), POINTER, DIMENSION(:,:) :: PHI_LS_P

IF (PREDICTOR) THEN
   PHI_LS_P => LS_PHI
ELSE
   PHI_LS_P => LS_PHI1
ENDIF

SR_X = 0.0_EB ; SR_Y = 0.0_EB

FLUX_JLOOP: DO J=1,LS_NY

   JM1 = J-1
   JP1 = J+1

   DO I=1,LS_NX

      IM1 = I-1
      IP1 = I+1

      F_EAST  = 0.5_EB*( PHI_LS_P(I,J) + PHI_LS_P(IP1,J) )
      F_WEST  = 0.5_EB*( PHI_LS_P(I,J) + PHI_LS_P(IM1,J) )
      F_NORTH = 0.5_EB*( PHI_LS_P(I,J) + PHI_LS_P(I,JP1) )
      F_SOUTH = 0.5_EB*( PHI_LS_P(I,J) + PHI_LS_P(I,JM1) )

      DPHIDX = (F_EAST-F_WEST)   * RDXL(I)
      DPHIDY = (F_NORTH-F_SOUTH) * RDYL(J)

      MAG_F = SQRT(DPHIDX**2 + DPHIDY**2)

      IF (MAG_F > 0._EB) THEN   !components of unit vector normal to PHI contours
         NORMAL_FIRELINE(1) = -DPHIDX/MAG_F
         NORMAL_FIRELINE(2) = -DPHIDY/MAG_F
         XSF =  DPHIDY
         YSF = -DPHIDX
      ELSE
         NORMAL_FIRELINE = 0._EB
         XSF=0._EB
         YSF=0._EB
      ENDIF

      ! ROS does not change with wind or slope
      IF (SURFACE(LS_SURF(I,J))%VEG_LSET_ROS_FIXED) THEN
         SR_X(I,J) = SURFACE(LS_SURF(I,J))%VEG_LSET_ROS_00*NORMAL_FIRELINE(1) !spread rate components
         SR_Y(I,J) = SURFACE(LS_SURF(I,J))%VEG_LSET_ROS_00*NORMAL_FIRELINE(2)
      ELSE
         ! Effective wind direction (theta) is clockwise from y-axis (Richards 1990)
         COS_THETA = COS(THETA_ELPS(I,J))
         SIN_THETA = SIN(THETA_ELPS(I,J))

         ROS_TMP = ROS_HEAD(I,J)

         ! Magnitude of wind speed at midflame height must be in units of m/s here

         UMF_DUM = UMF(I,J)/60.0_EB

         ! Length to breadth ratio of ellipse based on effective UMF (Bova et al., Eq. A6)

         LB = 0.936_EB * EXP(0.2566_EB * UMF_DUM) + 0.461_EB * EXP(-0.1548_EB * UMF_DUM) - 0.397_EB

         ! Constraint LB max = 8 from Finney 2004

         LB = MAX(1.0_EB,MIN(LB,8.0_EB))  ! (Bova et al., Eq. A7)

         ! Head to back ratio based on LB

         LBD = SQRT(LB**2 - 1.0_EB)
         HB = (LB + LBD) / (LB - LBD)

         ! A_ELPS and B_ELPS notation is consistent with Farsite and Richards

         B_ELPS =  0.5_EB * (ROS_TMP + ROS_TMP/HB)
         B_ELPS2 = B_ELPS**2
         A_ELPS =  B_ELPS / LB
         A_ELPS2=  A_ELPS**2
         C_ELPS =  B_ELPS - (ROS_TMP/HB)

         AROS  = XSF*COS_THETA - YSF*SIN_THETA
         BROS  = XSF*SIN_THETA + YSF*COS_THETA
         DENOM = A_ELPS2*BROS**2 + B_ELPS2*AROS**2

         IF (DENOM > 0._EB) THEN
            DENOM = 1._EB / SQRT(DENOM)
         ELSE
            DENOM = 0._EB
         ENDIF

         SR_X(I,J) = DENOM * ( A_ELPS2*COS_THETA*BROS - B_ELPS2*SIN_THETA*AROS) + C_ELPS*SIN_THETA  ! Bova et al., Eq. A8
         SR_Y(I,J) = DENOM * (-A_ELPS2*SIN_THETA*BROS - B_ELPS2*COS_THETA*AROS) + C_ELPS*COS_THETA  ! Bova et al., Eq. A9

      ENDIF

   ENDDO

ENDDO FLUX_JLOOP

END SUBROUTINE SPREAD_RATE


!> \brief Advective flux of the level set equation (LEVEL_SET_ADVECT_FLUX of vege.f90)
SUBROUTINE ADVECT_FLUX

INTEGER :: I,J
REAL(EB), DIMENSION(:) :: Z(4)
REAL(EB), POINTER, DIMENSION(:,:) :: PHI_LS_P,FLUX_LS_P,F_X,F_Y
REAL(EB) :: DPHIDX,DPHIDY,SR_X_AVG,SR_Y_AVG

F_X => WORK1
F_Y => WORK2

IF (PREDICTOR) THEN
   PHI_LS_P => LS_PHI
   FLUX_LS_P => FLUX0
ELSE
   PHI_LS_P => LS_PHI1
   FLUX_LS_P => FLUX1
ENDIF

DO J=1,LS_NY
   DO I=0,LS_NX
      Z(1) = PHI_LS_P(I-1,J)
      Z(2) = PHI_LS_P(I  ,J)
      Z(3) = PHI_LS_P(I+1,J)
      Z(4) = PHI_LS_P(I+2,J)
      SR_X_AVG = 0.5_EB*(SR_X(MIN(I+1,LS_NX),J)+SR_X(MAX(1,I),J))
      F_X(I,J) = SCALAR_FACE_VALUE_LS(SR_X_AVG,Z,LIMITER_LS)
   ENDDO
ENDDO

DO J=0,LS_NY
   DO I=1,LS_NX
      Z(1) = PHI_LS_P(I,J-1)
      Z(2) = PHI_LS_P(I,J  )
      Z(3) = PHI_LS_P(I,J+1)
      Z(4) = PHI_LS_P(I,J+2)
      SR_Y_AVG = 0.5_EB*(SR_Y(I,MIN(J+1,LS_NY))+SR_Y(I,MAX(1,J)))
      F_Y(I,J) = SCALAR_FACE_VALUE_LS(SR_Y_AVG,Z,LIMITER_LS)
   ENDDO
ENDDO

DO J=1,LS_NY
   DO I=1,LS_NX
      IF (.NOT. SURFACE(LS_SURF(I,J))%VEG_LSET_SPREAD) CYCLE
      DPHIDX = (F_X(I,J)-F_X(I-1,J))*RDXL(I)
      DPHIDY = (F_Y(I,J)-F_Y(I,J-1))*RDYL(J)
      FLUX_LS_P(I,J) = SR_X(I,J)*DPHIDX + SR_Y(I,J)*DPHIDY
   ENDDO
ENDDO

END SUBROUTINE ADVECT_FLUX


REAL(EB) FUNCTION SCALAR_FACE_VALUE_LS(SR_XY,Z,LIMITER)

INTEGER, INTENT(IN) :: LIMITER
REAL(EB) :: SR_XY
REAL(EB), INTENT(IN), DIMENSION(4) :: Z
REAL(EB) :: B,DZLOC,DZUP,R,ZUP,ZDWN

IF (SR_XY > 0._EB) THEN  ! flow is left to right
   DZLOC = Z(3)-Z(2)
   DZUP  = Z(2)-Z(1)
   IF (ABS(DZLOC) > 0._EB) THEN
      R = DZUP/DZLOC
   ELSE
      R = 0._EB
   ENDIF
   ZUP  = Z(2)
   ZDWN = Z(3)
ELSE  ! flow is right to left
   DZLOC = Z(3)-Z(2)
   DZUP  = Z(4)-Z(3)
   IF (ABS(DZLOC) > 0._EB) THEN
      R = DZUP/DZLOC
   ELSE
      R = 0._EB
   ENDIF
   ZUP  = Z(3)
   ZDWN = Z(2)
ENDIF

B = 0._EB
IF (LIMITER==2) THEN
   B = MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
ELSEIF (LIMITER==3) THEN
   B = MAX(0._EB,MIN(1._EB,R))
ENDIF

SCALAR_FACE_VALUE_LS = ZUP + 0.5_EB * B * ( ZDWN - ZUP )

END FUNCTION SCALAR_FACE_VALUE_LS


!> \brief Build the LS grid from level-0 mesh 1 (one mesh, flat terrain): LS cell = level-0 cell / R, surface index and wall state injected from the level-0 ground wall of the cell.
SUBROUTINE LS_INIT_FROM_MESH(R)

USE MPI_F08
USE GLOBAL_CONSTANTS, ONLY: MY_RANK,PROCESS,N_MPI_PROCESSES,N_SURF
INTEGER, INTENT(IN) :: R
TYPE(MESH_TYPE), POINTER :: M
INTEGER :: I,J,IC,JC,IG,JG,IW,NX,NY,NM,NXC,NYC,IERR
REAL(EB) :: ZB,X0,X1,Y0,Y1,DX0,DY0
REAL(EB), ALLOCATABLE :: DXY(:),G(:,:,:)   ! DXY: LS spacings in x then y; G(:,:,1:6): T_IGN, BURN_DURATION, AREA_ADJUST, TAU_LS, surface index, coverage count
INTEGER, ALLOCATABLE :: SF0(:,:)
LOGICAL, ALLOCATABLE :: SEEN(:)
REAL(EB) :: E_ROTH,BETA_OP_ROTH
TYPE(SURFACE_TYPE), POINTER :: SFP

IF (LS_ACTIVE) RETURN
IF (LEVEL_SET_MODE==2 .OR. LEVEL_SET_MODE==5) CALL ABORT_INIT('the LS grid supports LEVEL_SET_MODE 1, 3 and 4 only')
IF (NMESHES>1 .AND. .NOT.LS_DRIVE) CALL ABORT_INIT('the LS shadow test supports one level-0 mesh only')
! the ground meshes (lowest ZS) tile the LS grid; uniform level-0 cells of one size in all of them
ZB = MINVAL(MESHES(1:NMESHES)%ZS)
X0 = HUGE(1._EB) ; X1 = -HUGE(1._EB) ; Y0 = X0 ; Y1 = X1 ; DX0 = -1._EB ; DY0 = -1._EB
ALLOCATE(LS_MI0(NMESHES),LS_MJ0(NMESHES),LS_MGROUND(NMESHES))
DO NM=1,NMESHES
   M => MESHES(NM)
   LS_MGROUND(NM) = (M%ZS<=ZB)
   IF (.NOT.LS_MGROUND(NM)) CYCLE
   X0 = MIN(X0,M%XS) ; X1 = MAX(X1,M%XF) ; Y0 = MIN(Y0,M%YS) ; Y1 = MAX(Y1,M%YF)
   IF (DX0<0._EB) THEN ; DX0 = (M%XF-M%XS)/REAL(M%IBAR,EB) ; DY0 = (M%YF-M%YS)/REAL(M%JBAR,EB) ; ENDIF
   IF (ABS((M%XF-M%XS)/REAL(M%IBAR,EB)-DX0)>1.E-6_EB*DX0 .OR. ABS((M%YF-M%YS)/REAL(M%JBAR,EB)-DY0)>1.E-6_EB*DY0) &
      CALL ABORT_INIT('the LS grid needs level-0 ground meshes of one cell size')
ENDDO
NXC = NINT((X1-X0)/DX0) ; NYC = NINT((Y1-Y0)/DY0)
DO NM=1,NMESHES
   LS_MI0(NM) = NINT((MESHES(NM)%XS-X0)/DX0) ; LS_MJ0(NM) = NINT((MESHES(NM)%YS-Y0)/DY0)
ENDDO
NX = NXC*R ; NY = NYC*R
ALLOCATE(DXY(NX+NY),G(NX,NY,6),SF0(NX,NY))
DXY = -HUGE(1._EB) ; G = 0._EB
! every rank fills the LS cells of its ground meshes; one exact reduction makes the map global (each cell has one owner: SUM of one value and zeros, MAX for the spacings)
DO NM=1,NMESHES
   IF (.NOT.LS_MGROUND(NM) .OR. PROCESS(NM)/=MY_RANK) CYCLE
   M => MESHES(NM)
   IF (.NOT.ALLOCATED(M%CELL_INDEX)) CYCLE
   IF (MAXVAL(M%DX(1:M%IBAR))-MINVAL(M%DX(1:M%IBAR))>1.E-6_EB*DX0 .OR. MAXVAL(M%DY(1:M%JBAR))-MINVAL(M%DY(1:M%JBAR))>1.E-6_EB*DY0) &
      CALL ABORT_INIT('the LS grid needs uniform level-0 cells (no stretched meshes)')
   DO I=1,M%IBAR*R
      IG = LS_MI0(NM)*R+I
      DXY(IG) = M%DX((I-1)/R+1)/REAL(R,EB)
      IF (R==1) DXY(IG) = M%DX(I)
   ENDDO
   DO J=1,M%JBAR*R
      JG = LS_MJ0(NM)*R+J
      DXY(NX+JG) = M%DY((J-1)/R+1)/REAL(R,EB)
      IF (R==1) DXY(NX+JG) = M%DY(J)
   ENDDO
   DO J=1,M%JBAR*R
      DO I=1,M%IBAR*R
         IC = (I-1)/R+1 ; JC = (J-1)/R+1 ; IG = LS_MI0(NM)*R+I ; JG = LS_MJ0(NM)*R+J
         IF (ABS(M%Z_LS(IC,JC))>0._EB .OR. M%K_LS(IC,JC)/=1) CALL ABORT_INIT('the LS grid supports flat terrain only')
         G(IG,JG,5) = REAL(M%LS_SURF_INDEX(IC,JC),EB)
         IW = M%CELL(M%CELL_INDEX(IC,JC,M%K_LS(IC,JC)))%WALL_INDEX(-3)
         G(IG,JG,1) = M%BOUNDARY_PROP1(M%WALL(IW)%B1_INDEX)%T_IGN
         G(IG,JG,2) = M%BOUNDARY_PROP1(M%WALL(IW)%B1_INDEX)%BURN_DURATION
         G(IG,JG,3) = M%BOUNDARY_PROP1(M%WALL(IW)%B1_INDEX)%AREA_ADJUST
         G(IG,JG,4) = M%BOUNDARY_PROP2(M%WALL(IW)%B2_INDEX)%TAU_LS
         G(IG,JG,6) = G(IG,JG,6) + 1._EB
      ENDDO
   ENDDO
ENDDO
IF (N_MPI_PROCESSES>1) THEN
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,DXY,NX+NY,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERR)
   CALL MPI_ALLREDUCE(MPI_IN_PLACE,G,6*NX*NY,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
ENDIF
IF (ANY(G(:,:,6)/=1._EB)) CALL ABORT_INIT('the level-0 ground meshes do not tile the LS grid (gap or overlap)')
SF0 = NINT(G(:,:,5))
! surface constants that INITIALIZE_LEVEL_SET_FIRESPREAD_2 (vege.f90) sets only for the surfaces of the rank's own meshes: set them for every spreading surface of the
! global map, same formulas (identical values where the rank had them already)
ALLOCATE(SEEN(0:N_SURF)) ; SEEN = .FALSE.
DO J=1,NY
   DO I=1,NX
      SEEN(SF0(I,J)) = .TRUE.
   ENDDO
ENDDO
DO I=0,N_SURF
   IF (.NOT.SEEN(I)) CYCLE
   SFP => SURFACE(I)
   IF (.NOT.SFP%VEG_LSET_SPREAD) CYCLE
   SFP%VEG_LSET_HT = MAX(0.001_EB,SFP%VEG_LSET_HT)
   SFP%B_ROTH = 0.15988_EB * (SFP%VEG_LSET_SIGMA**0.54_EB)
   SFP%C_ROTH = 7.47_EB * EXP(-0.8711_EB * (SFP%VEG_LSET_SIGMA**0.55_EB))
   E_ROTH = 0.715_EB * EXP(-0.01094_EB * SFP%VEG_LSET_SIGMA)
   BETA_OP_ROTH = 0.20395_EB * (SFP%VEG_LSET_SIGMA**(-0.8189_EB))
   SFP%BETA_ROTH = (SFP%VEG_LSET_BETA/BETA_OP_ROTH)**(-E_ROTH)
ENDDO
CALL LS_INIT(NX,NY,DXY(1:NX),DXY(NX+1:NX+NY),SF0,G(:,:,1),G(:,:,2),G(:,:,3),G(:,:,4))
LS_RATIO = R
LS_X0 = X0 ; LS_Y0 = Y0
ALLOCATE(LS_FLUX(NX,NY)) ; LS_FLUX = 0._EB
IF (LS_DRIVE .AND. (STORE_FIRE_ARRIVAL .OR. STORE_FIRE_RESIDENCE)) &
   CALL ABORT_INIT('--ls-fine does not support fire arrival / residence outputs (the level-0 walls carry LS cell means)')

CONTAINS
SUBROUTINE ABORT_INIT(MSG)
CHARACTER(*), INTENT(IN) :: MSG
WRITE(0,'(A)') 'ERROR: '//MSG
STOP 1
END SUBROUTINE ABORT_INIT

END SUBROUTINE LS_INIT_FROM_MESH



!> \brief Level-0 PHI_LS of mesh NM (cell means of the LS cells under it) for the PHI plotfile / tagging and the "LEVEL SET VALUE" output; the ghost cells that lie in
!> another ground mesh get that mesh's mean too (the global grid has it), ghost cells outside the domain are left alone
SUBROUTINE LS_RESTRICT_PHI(NM)

INTEGER, INTENT(IN) :: NM
TYPE(MESH_TYPE), POINTER :: M
INTEGER :: I,J,II,JJ,R,IG,JG,ILO,IHI,JLO,JHI
REAL(EB) :: S

IF (.NOT.LS_ACTIVE .OR. NM>SIZE(LS_MGROUND)) RETURN
IF (.NOT.LS_MGROUND(NM)) RETURN
M => MESHES(NM) ; R = LS_RATIO
ILO = MERGE(0,1,LS_MI0(NM)>0) ; IHI = MERGE(M%IBP1,M%IBAR,(LS_MI0(NM)+M%IBAR)*R<LS_NX)
JLO = MERGE(0,1,LS_MJ0(NM)>0) ; JHI = MERGE(M%JBP1,M%JBAR,(LS_MJ0(NM)+M%JBAR)*R<LS_NY)
DO J=JLO,JHI
   DO I=ILO,IHI
      IG = LS_MI0(NM)+I ; JG = LS_MJ0(NM)+J
      S = 0._EB
      DO JJ=(JG-1)*R+1,JG*R
         DO II=(IG-1)*R+1,IG*R
            S = S + LS_PHI(II,JJ)
         ENDDO
      ENDDO
      M%PHI_LS(I,J) = S/REAL(R*R,EB)
   ENDDO
ENDDO

END SUBROUTINE LS_RESTRICT_PHI


!> \brief --ls-fine: one stage of the LS grid with the wind gathered for this stage (C++ density(), after ls_fine_gather_wind; once per rank and stage)
SUBROUTINE FDS_LS_DRIVE_STEP(T,DT) BIND(C,NAME='fds_ls_drive_step')
REAL(C_DOUBLE), VALUE :: T,DT
IF (.NOT.LS_DRIVE) RETURN
IF (.NOT.LS_WIND_SET) THEN
   WRITE(0,'(A)') 'ERROR: ls-fine: the wind of the LS grid was not gathered'
   STOP 1
ENDIF
LS_WIND_SET = .FALSE.
CALL LS_STEP(T,DT)
LS_FLUX_T = -HUGE(1._EB)   ! the LS state changed: LS_APPLY_L0 recomputes LS_FLUX
END SUBROUTINE FDS_LS_DRIVE_STEP


!> \brief Fraction of the plateau fuel flux of an LS cell at this stage: the trapezoid of WALL_BC (wall.f90) from the LS cell's T_IGN, TAU_LS and BURN_DURATION
REAL(EB) FUNCTION LS_FRAC(I,J,T,DT)

INTEGER, INTENT(IN) :: I,J
REAL(EB), INTENT(IN) :: T,DT
REAL(EB) :: TSI

LS_FRAC = 0._EB
IF (T<LS_T_IGN(I,J) .OR. LS_T_IGN(I,J)+LS_BURN_DUR(I,J)<T) RETURN
IF (PREDICTOR) THEN
   TSI = T+DT-LS_T_IGN(I,J)
ELSE
   TSI = T-LS_T_IGN(I,J)
ENDIF
IF (TSI<=LS_TAU(I,J)) THEN
   LS_FRAC = TSI/LS_TAU(I,J)
ELSEIF (TSI<=LS_BURN_DUR(I,J)-LS_TAU(I,J)) THEN
   LS_FRAC = 1._EB
ELSE
   LS_FRAC = 1._EB-(TSI-(LS_BURN_DUR(I,J)-LS_TAU(I,J)))/LS_TAU(I,J)
ENDIF

END FUNCTION LS_FRAC


!> \brief Before WALL_BC of level-0 mesh NM (FDS_G_WALL_BC): the fuel flux of every LS cell at this stage (LS_FLUX), and on every level-0 spreading wall the LS cell means:
!> B2%PHI_LS = mean PHI; the wall trapezoid of WALL_BC is made to return mean(flux)/(MASS_FLUX*ADJUST_FAC) of the wall: T_IGN = T_stage - F*LAM, TAU_LS = LAM, BURN_DURATION huge
!> (so TSI/TAU = F, any value of F in 0..1) and AREA_ADJUST = max(1,total) with F = total/AREA_ADJUST. M_DOT_G_PP_ADJUST (and the burning rate output) is then the mean LS flux.
!> No flux is touched when LEVEL_SET_COUPLED_FIRE is off (only PHI_LS).
SUBROUTINE LS_APPLY_L0(T,DT,NM)

REAL(EB), INTENT(IN) :: T,DT
INTEGER, INTENT(IN) :: NM
REAL(EB), PARAMETER :: LAM=1.E6_EB
TYPE(MESH_TYPE), POINTER :: M
TYPE(SURFACE_TYPE), POINTER :: SF
INTEGER :: I,J,II,JJ,IW,R
REAL(EB) :: SM,SP,MF0,TOT,AA,FR,TEFF,RR

IF (.NOT.LS_DRIVE .OR. .NOT.LS_ACTIVE .OR. NM>SIZE(LS_MGROUND)) RETURN
TEFF = MERGE(T+DT,T,PREDICTOR)
IF (LEVEL_SET_COUPLED_FIRE .AND. (TEFF/=LS_FLUX_T .OR. (PREDICTOR.NEQV.LS_FLUX_PRED))) THEN   ! once per stage on every rank (the replicated grid; fine walls read it)
   LS_FLUX_T = TEFF ; LS_FLUX_PRED = PREDICTOR
   DO J=1,LS_NY
      DO I=1,LS_NX
         SF => SURFACE(LS_SURF(I,J))
         LS_FLUX(I,J) = 0._EB
         IF (SF%VEG_LSET_SPREAD .AND. ABS(SUM(SF%MASS_FLUX(1:N_TRACKED_SPECIES)))>TWENTY_EPSILON_EB) &
            LS_FLUX(I,J) = LS_FRAC(I,J,T,DT)*SUM(SF%MASS_FLUX(1:N_TRACKED_SPECIES))*SF%M_DOT_G_PP_ADJUST_FAC*LS_AREA_ADJ(I,J)
      ENDDO
   ENDDO
ENDIF
IF (.NOT.LS_MGROUND(NM)) RETURN
M => MESHES(NM) ; R = LS_RATIO ; RR = REAL(R*R,EB)
DO J=1,M%JBAR
   DO I=1,M%IBAR
      IW = M%CELL(M%CELL_INDEX(I,J,1))%WALL_INDEX(-3)
      IF (IW<1) CYCLE
      IF (M%WALL(IW)%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE
      SF => SURFACE(M%WALL(IW)%SURF_INDEX)
      IF (.NOT.SF%VEG_LSET_SPREAD) CYCLE
      SM = 0._EB ; SP = 0._EB
      DO JJ=(LS_MJ0(NM)+J-1)*R+1,(LS_MJ0(NM)+J)*R
         DO II=(LS_MI0(NM)+I-1)*R+1,(LS_MI0(NM)+I)*R
            SP = SP + LS_PHI(II,JJ)
            IF (LEVEL_SET_COUPLED_FIRE) SM = SM + LS_FLUX(II,JJ)
         ENDDO
      ENDDO
      ASSOCIATE(B1=>M%BOUNDARY_PROP1(M%WALL(IW)%B1_INDEX),B2=>M%BOUNDARY_PROP2(M%WALL(IW)%B2_INDEX))
      B2%PHI_LS = SP/RR
      IF (LEVEL_SET_COUPLED_FIRE) THEN
         MF0 = SUM(SF%MASS_FLUX(1:N_TRACKED_SPECIES))*SF%M_DOT_G_PP_ADJUST_FAC
         IF (ABS(MF0)>TWENTY_EPSILON_EB) THEN
            TOT = SM/RR/MF0
            B2%TAU_LS = LAM
            IF (TOT>0._EB) THEN
               AA = MAX(1._EB,TOT) ; FR = TOT/AA
               B1%T_IGN = TEFF-FR*LAM ; B1%BURN_DURATION = 1.E30_EB ; B1%AREA_ADJUST = AA
            ELSE
               B1%T_IGN = -1._EB ; B1%BURN_DURATION = -1._EB ; B1%AREA_ADJUST = 1._EB
            ENDIF
         ENDIF
      ENDIF
      END ASSOCIATE
   ENDDO
ENDDO

END SUBROUTINE LS_APPLY_L0


LOGICAL FUNCTION LS_FLUX_ACTIVE()
LS_FLUX_ACTIVE = LS_DRIVE .AND. LS_ACTIVE .AND. LEVEL_SET_COUPLED_FIRE
END FUNCTION LS_FLUX_ACTIVE


!> \brief Mean LS flux of the LS cells under a fine wall centred at (X,Y) with sides DXL, DYL (the wall of a cell of a finer level)
REAL(EB) FUNCTION LS_BLOCK_FLUX(X,Y,DXL,DYL)

REAL(EB), INTENT(IN) :: X,Y,DXL,DYL
INTEGER :: ILO,IHI,JLO,JHI,I,J
REAL(EB) :: S

ILO = MAX(1,NINT((X-0.5_EB*DXL-LS_X0)/LS_DX(1))+1) ; IHI = MIN(LS_NX,NINT((X+0.5_EB*DXL-LS_X0)/LS_DX(1)))
JLO = MAX(1,NINT((Y-0.5_EB*DYL-LS_Y0)/LS_DY(1))+1) ; JHI = MIN(LS_NY,NINT((Y+0.5_EB*DYL-LS_Y0)/LS_DY(1)))
S = 0._EB
DO J=JLO,JHI
   DO I=ILO,IHI
      S = S + LS_FLUX(I,J)
   ENDDO
ENDDO
LS_BLOCK_FLUX = S/REAL(MAX(1,(IHI-ILO+1)*(JHI-JLO+1)),EB)

END FUNCTION LS_BLOCK_FLUX


!> \brief C entry: the LS PHI of the global grid (x fastest, NX*NY values; 0 and returns 1 if the grid does not exist yet)
INTEGER(C_INT) FUNCTION FDS_LS_GET_PHI(NX,NY,BUF) BIND(C,NAME='fds_ls_get_phi')
INTEGER(C_INT), VALUE :: NX,NY
REAL(C_DOUBLE), INTENT(OUT) :: BUF(NX*NY)
FDS_LS_GET_PHI = 1
IF (.NOT.LS_ACTIVE) RETURN
IF (NX/=LS_NX .OR. NY/=LS_NY) RETURN
BUF = RESHAPE(LS_PHI(1:NX,1:NY),(/NX*NY/))
FDS_LS_GET_PHI = 0
END FUNCTION FDS_LS_GET_PHI


!> \brief C entry: make sure the LS grid exists (ratio R to the level-0 mesh); returns its size and the largest wind reference height of the spreading cells
SUBROUTINE FDS_LS_PREPARE(R,NX,NY,HMAX) BIND(C,NAME='fds_ls_prepare')

INTEGER(C_INT), VALUE :: R
INTEGER(C_INT), INTENT(OUT) :: NX,NY
REAL(C_DOUBLE), INTENT(OUT) :: HMAX
INTEGER :: I,J

CALL LS_INIT_FROM_MESH(R)
NX = LS_NX ; NY = LS_NY ; HMAX = 0._EB
DO J=1,LS_NY
   DO I=1,LS_NX
      IF (SURFACE(LS_SURF(I,J))%VEG_LSET_SPREAD) HMAX = MAX(HMAX,LS_REF_HEIGHT(I,J))
   ENDDO
ENDDO

END SUBROUTINE FDS_LS_PREPARE


!> \brief C entry: reference height of every LS column (x fastest), -1 where the cell does not spread
SUBROUTINE FDS_LS_REF_HEIGHTS(NX,NY,H) BIND(C,NAME='fds_ls_ref_heights')
INTEGER(C_INT), VALUE :: NX,NY
REAL(C_DOUBLE), INTENT(OUT) :: H(NX,NY)
INTEGER :: I,J
H = -1._EB
DO J=1,MIN(NY,LS_NY)
   DO I=1,MIN(NX,LS_NX)
      IF (SURFACE(LS_SURF(I,J))%VEG_LSET_SPREAD) H(I,J) = LS_REF_HEIGHT(I,J)
   ENDDO
ENDDO
END SUBROUTINE FDS_LS_REF_HEIGHTS


!> \brief C entry: wind at the reference height of every spreading LS cell from the two leaf samples of its column that bracket the height (pre-reduced by LsFine.cpp:
!> ZLO/ULO/VLO the highest sample below it, -1E300 if none; ZHI/UHI/VHI the lowest at or above it, 1E300 if none; heights above the ground, cell-centre wind).
!> LS_REF_WIND gives the same result as with the whole column. NBAD counts spreading columns without any sample.
SUBROUTINE FDS_LS_WIND_FROM_PAIRS(ZLO,ULO,VLO,ZHI,UHI,VHI,NBAD) BIND(C,NAME='fds_ls_wind_from_pairs')

REAL(C_DOUBLE), INTENT(IN), DIMENSION(LS_NX,LS_NY) :: ZLO,ULO,VLO,ZHI,UHI,VHI
INTEGER(C_INT), INTENT(OUT) :: NBAD
INTEGER :: I,J,NS
REAL(EB) :: ZC(2),UC(2),VC(2)

NBAD = 0
DO J=1,LS_NY
   DO I=1,LS_NX
      IF (.NOT.SURFACE(LS_SURF(I,J))%VEG_LSET_SPREAD) CYCLE
      NS = 0
      IF (ZLO(I,J)>-1.E299_EB) THEN ; NS = NS+1 ; ZC(NS) = ZLO(I,J) ; UC(NS) = ULO(I,J) ; VC(NS) = VLO(I,J) ; ENDIF
      IF (ZHI(I,J)< 1.E299_EB) THEN ; NS = NS+1 ; ZC(NS) = ZHI(I,J) ; UC(NS) = UHI(I,J) ; VC(NS) = VHI(I,J) ; ENDIF
      IF (NS==0) THEN
         NBAD = NBAD+1 ; LS_UREF(I,J) = 0._EB ; LS_VREF(I,J) = 0._EB
      ELSE
         CALL LS_REF_WIND(NS,ZC(1:NS),UC(1:NS),VC(1:NS),LS_REF_HEIGHT(I,J),LS_UREF(I,J),LS_VREF(I,J))
      ENDIF
   ENDDO
ENDDO
LS_WIND_SET = .TRUE.

END SUBROUTINE FDS_LS_WIND_FROM_PAIRS


!> \brief --ls-fine-shadow: C entry that switches the shadow mode on
SUBROUTINE FDS_LS_SET_MODE(MODE) BIND(C,NAME='fds_ls_set_mode')
INTEGER(C_INT), VALUE :: MODE
LS_SHADOW_MODE = MODE
IF (MODE==2) THEN   ! --ls-fine: drive
   LS_DRIVE = .TRUE. ; LS_SHADOW_MODE = 0
ENDIF
END SUBROUTINE FDS_LS_SET_MODE


!> \brief Shadow test, before the original LEVEL_SET_FIRESPREAD of level-0 mesh NM: build the LS grid from the mesh on the first call (grid = the mesh),
!> sample the wind column of every cell the way the original does, and advance the LS grid one stage. Level-0 mesh 1 only (one mesh, flat terrain).
SUBROUTINE FDS_LS_SHADOW_RUN(T,DT,NM)

REAL(EB), INTENT(IN) :: T,DT
INTEGER, INTENT(IN) :: NM
TYPE(MESH_TYPE), POINTER :: M
INTEGER :: I,J,K,NS
REAL(EB), ALLOCATABLE :: ZS(:),US(:),VS(:)

IF (NM/=1 .OR. NMESHES/=1) RETURN
M => MESHES(NM)
IF (.NOT.LS_ACTIVE) CALL LS_INIT_FROM_MESH(1)
IF (LS_RATIO/=1) CALL ABORT_LS('ls-fine-shadow: the LS grid must equal the level-0 mesh (ratio 1)')
IF (LEVEL_SET_MODE==2 .OR. LEVEL_SET_MODE==5) CALL ABORT_LS('ls-fine-shadow: LEVEL_SET_MODE 2 and 5 are not supported')
IF (.NOT.LS_WIND_SET) THEN    ! no gathered wind for this stage: sample the level-0 mesh the way the original does
   ALLOCATE(ZS(M%KBAR),US(M%KBAR),VS(M%KBAR))
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         IF (.NOT.SURFACE(LS_SURF(I,J))%VEG_LSET_SPREAD) CYCLE
         NS = 0
         DO K=M%K_LS(I,J),M%KBAR
            NS = NS+1
            ZS(NS) = M%ZC(K)-M%Z_LS(I,J)
            US(NS) = 0.5_EB*(M%U(I-1,J,K)+M%U(I,J,K))
            VS(NS) = 0.5_EB*(M%V(I,J-1,K)+M%V(I,J,K))
         ENDDO
         CALL LS_REF_WIND(NS,ZS,US,VS,LS_REF_HEIGHT(I,J),LS_UREF(I,J),LS_VREF(I,J))
      ENDDO
   ENDDO
ENDIF
LS_WIND_SET = .FALSE.
CALL LS_STEP(T,DT)

CONTAINS
SUBROUTINE ABORT_LS(MSG)
CHARACTER(*), INTENT(IN) :: MSG
WRITE(0,'(A)') 'ERROR: '//MSG
STOP 1
END SUBROUTINE ABORT_LS

END SUBROUTINE FDS_LS_SHADOW_RUN


!> \brief Shadow test, after the original LEVEL_SET_FIRESPREAD: compare PHI (both stages), T_IGN, TAU_LS, BURN_DURATION, AREA_ADJUST of every cell bit for bit
!> and write the running totals to <CHID>_ls_shadow.txt.
SUBROUTINE FDS_LS_SHADOW_CMP(NM)

INTEGER, INTENT(IN) :: NM
TYPE(MESH_TYPE), POINTER :: M
INTEGER :: I,J,IW,IU
REAL(EB) :: A,B

IF (NM/=1 .OR. NMESHES/=1 .OR. .NOT.LS_ACTIVE) RETURN
M => MESHES(NM)
DO J=-1,M%JBAR+2
   DO I=-1,M%IBAR+2
      IF (PREDICTOR) THEN
         A = LS_PHI1(I,J) ; B = M%PHI1_LS(I,J)
      ELSE
         A = LS_PHI(I,J) ; B = M%PHI_LS(I,J)
      ENDIF
      IF (A/=B) THEN
         SH_BAD_PHI = SH_BAD_PHI+1 ; SH_MAXDIFF = MAX(SH_MAXDIFF,ABS(A-B))
      ENDIF
   ENDDO
ENDDO
IF (.NOT.PREDICTOR) THEN
   SH_CELLS_IGN = 0
   DO J=1,M%JBAR
      DO I=1,M%IBAR
         IW = M%CELL(M%CELL_INDEX(I,J,M%K_LS(I,J)))%WALL_INDEX(-3)
         IF (LS_T_IGN(I,J)<9.E5_EB) SH_CELLS_IGN = SH_CELLS_IGN+1
         IF (LS_T_IGN(I,J)/=M%BOUNDARY_PROP1(M%WALL(IW)%B1_INDEX)%T_IGN .OR. LS_BURN_DUR(I,J)/=M%BOUNDARY_PROP1(M%WALL(IW)%B1_INDEX)%BURN_DURATION .OR. &
             LS_AREA_ADJ(I,J)/=M%BOUNDARY_PROP1(M%WALL(IW)%B1_INDEX)%AREA_ADJUST .OR. LS_TAU(I,J)/=M%BOUNDARY_PROP2(M%WALL(IW)%B2_INDEX)%TAU_LS) SH_BAD_WALL = SH_BAD_WALL+1
      ENDDO
   ENDDO
   SH_STEPS = SH_STEPS+1
   OPEN(NEWUNIT=IU,FILE=TRIM(CHID)//'_ls_shadow.txt',STATUS='REPLACE')
   WRITE(IU,'(A,I0,A,I0,A,I0,A,ES10.3,A,I0)') 'LS shadow: stages compared with the original, steps ',SH_STEPS,', PHI values differing ',SH_BAD_PHI, &
      ', cells with differing T_IGN/TAU/BURN_DURATION/AREA_ADJUST ',SH_BAD_WALL,', max |dPHI| ',SH_MAXDIFF,', ignited cells ',SH_CELLS_IGN
   CLOSE(IU)
ENDIF

END SUBROUTINE FDS_LS_SHADOW_CMP

END MODULE FDS_LS_GRID
