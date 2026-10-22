!> \brief Fine-level mesh objects of the AMR driver (D-056 option B, needs patch 0007 in Source/mesh.f90; compiled when FDS_FINE_B is defined, see CMakeLists.txt).
!> BUILD_FINE_BOX makes the mesh object of one axis-aligned uniform Cartesian fine box (moved here unchanged from draft/fds_fine_box_b.f90, which now uses it).
!> FDS_FINE_LEVEL_CREATE / FDS_FINE_LEVEL_DESTROY are the C entry points of TimeLoop::bind_level: they make and drop FINE_LEVEL(LEV) for the boxes of this rank.
!> The data views of the boxes are bound by FDS_FINE_B_SET_VIEW (fds_box_obj.f90). Not validated with the oneAPI compilers; DRAFT with patches 0007 to 0009.
MODULE FDS_FINE_LEVEL

USE ISO_C_BINDING
USE PRECISION_PARAMETERS, ONLY: EB
USE GLOBAL_CONSTANTS, ONLY: N_TOTAL_SCALARS,FLUX_LIMITER_MW_CORRECTION,N_TRACKED_SPECIES,TMPA,RHOA,RSUM0,U0,V0,W0,PARTICLE_DRAG,CHECK_VN,REAC_SOURCE_CHECK, &
   STORE_SPECIES_FLUX,N_ZONE,STRATIFICATION,I_RAMP_TMP0_Z,I_RAMP_P0_Z,P_INF,STORE_PRESSURE_POISSON_RESIDUAL,STORE_CUTCELL_DIVERGENCE,STORE_CARTESIAN_DIVERGENCE,N_LP_ARRAY_INDICES,N_REACTIONS,RADIATION,N_LAGRANGIAN_CLASSES,PARTICLE_FILE,TMPM,NULL_BOUNDARY,DEFAULT_SURF_INDEX
USE MESH_VARIABLES, ONLY: MESHES,MESH_TYPE
USE FDS_LS_GRID, ONLY: LS_FLUX_ACTIVE,LS_BLOCK_FLUX
USE GLOBAL_CONSTANTS, ONLY: LOWER_MESH_INDEX,TURB_MODEL,C_SMAGORINSKY,CONSMAG,DYNSMAG,T_BEGIN
USE GLOBAL_CONSTANTS, ONLY: IBAR_MAX,JBAR_MAX,KBAR_MAX,TMPA4,SIGMA,NMESHES,PI
USE RADCONS, ONLY: UIIDIM
USE PHYSICAL_FUNCTIONS, ONLY: LES_FILTER_WIDTH_FUNCTION,GET_WIND_AT_HEIGHT
USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
USE MESH_POINTERS, ONLY: FINE_LEVEL,FINE_LEVEL_TYPE
USE FDS_AMREX_HOOKS, ONLY: HRR_DET,HRR_DET_ALLOC,HRR_QW,HRR_MW,HRR_MADD,HRR_QENTH,HRR_QENTH_ON
USE FDS_AMREX_HOOKS, ONLY: FDS_HOOK_SET_FINE_READY,FINE_DEVC_PROC,LEAF_WALL_COUNT_PROC,LEAF_WALL_PROC,LEAF_WALL_COVERED_PROC,LEAF_BOXWALL_COUNT_PROC,LEAF_BOXWALL_PROC,RAD_WALLS_PROC,FINE_WALLS_PTR_PROC
USE FDS_KERNELS, ONLY: FDS_K_RESERVE_BOXES

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE
PUBLIC :: FDS_FINE_WALLS_PARENT_XCHG,FDS_FINE_SET_LISTS,BUILD_FINE_BOX,FDS_FINE_LEVEL_CREATE,FDS_FINE_LEVEL_DESTROY,FDS_FINE_SET_SOLID,FDS_FINE_WALLS_BUILD,FDS_FINE_WALLS_RAD_XCHG,FDS_FINE_WALLS_SUMS,FDS_FINE_FUEL_MDOT

!> Hidden fine-level wall cells (patch 0026): per local box a separate mesh object W holding only WALL/BOUNDARY_COORD/BOUNDARY_PROP1/2 (the box itself keeps WALL(0:0), so the
!> kernels never see them), and per wall the level-0 parent wall (mesh PNM, wall PIW) whose face contains it. Used for the composite HRR-file wall terms only.
TYPE FWALL_BOX_TYPE
   TYPE(MESH_TYPE), POINTER :: W=>NULL()
   INTEGER :: NW=0
   INTEGER, ALLOCATABLE :: PNM(:),PIW(:)
   INTEGER, ALLOCATABLE :: RNM(:),RIW(:)   !< ruling 2 interim remap of a SOLID wall without parent (PNM = 0): the nearest level-0 wall of the same surface, boundary type and IOR in the same face plane (0: none)
   INTEGER, ALLOCATABLE :: OWN(:)          !< level-0 mesh whose extent contains the gas cell of the wall (devices are visited per level-0 mesh)
   INTEGER, ALLOCATABLE :: LIDX(:)         !< indices of the leaf walls of the box (built by REFRESH_LEAF)
   LOGICAL, ALLOCATABLE :: LEAF(:)         !< the gas cell of the wall is not covered by the next level (leaf sums); refreshed by REFRESH_LEAF
   REAL(EB), ALLOCATABLE :: PV(:,:)   !< radiation state pulled from the level L-1 wall of the same face (TMP_F, EMISSIVITY, Q_RAD_OUT, T_IGN), FDS_FINE_WALLS_RAD_XCHG
   LOGICAL, ALLOCATABLE :: HAVEP(:)
   INTEGER(8), ALLOCATABLE :: RKEY(:)      !< key of the level-0 parent face when the parent mesh is on another rank (global cell, IOR; 0 if none), FDS_FINE_WALLS_PARENT_XCHG
   REAL(EB), ALLOCATABLE :: PS(:,:)        !< gathered level-0 state (PARENT_STATE layout, then T_IGN, Q_RAD_IN, RHO_D_F(1:NZ)), valid where HAVES: the parent on another rank (RKEY/=0),
                                           !< or, for a wall without parent, the nearest gathered level-0 wall of the same surface and IOR (FDS_FINE_WALLS_PARENT_XCHG)
   LOGICAL, ALLOCATABLE :: HAVES(:)
   REAL(EB), ALLOCATABLE :: PQ(:)          !< leaf wall over a covered thermally thick level-0 wall: Q_CON_F*AREA*DT of this step's predictor (FDS_FINE_L0Q_COMMIT)
END TYPE
TYPE FWALL_LEVEL_TYPE
   TYPE(FWALL_BOX_TYPE), ALLOCATABLE :: B(:)
END TYPE
TYPE(FWALL_LEVEL_TYPE), ALLOCATABLE, SAVE :: FWALL(:)
INTEGER, ALLOCATABLE, SAVE :: LW_OFF(:),LW_LEV(:),LW_BOX(:),LW_W(:)   !< leaf walls grouped by level-0 owner mesh (CSR, built by REFRESH_LEAF): LEAF_WALL_COUNT/LEAF_WALL are O(1)
REAL(EB), SAVE :: ORPHAN_E = 0._EB   !< ruling 2: Q_CON_F*AREA*DT of the leaf SOLID hidden walls without a level-0 wall (PNM = 0; no solid takes it), predictor stages
INTEGER, SAVE :: ORPHAN_NMAX = 0     !< largest number of such walls of one box in one predictor stage (local)
REAL(EB), SAVE :: REMAP_E = 0._EB    !< the part of them remapped (RNM/RIW) into a covered thick level-0 solid through L0Q (WARNING)
INTEGER, SAVE :: REMAP_NMAX = 0
INTEGER, SAVE :: BUILD_VER=0,LEAF_VER=-1   !< bumped at every level create/destroy and wall build; LEAF is recomputed when they differ

!> Global OBST and VENT lists (set once by the driver, fds_fine_set_lists): the per-level wall builder snaps them at the dx of each box (FDS rounding, as ObstMask.H).
INTEGER, SAVE :: NOB=0,NVT=0,FVT=0   ! FVT: VENT of the last FACE_SURFACE call (0 none)
LOGICAL, SAVE :: BED_RUN=.FALSE.   ! FDS_FINE_SET_BED
REAL(EB), SAVE :: QENTH_V = 0._EB   !< FDS_SET_QENTH
LOGICAL, SAVE :: QENTH_SET = .FALSE.
REAL(EB), ALLOCATABLE, SAVE :: OBX(:,:),VTX(:,:)   ! XB(6,*); VENT XB, XYZ, SPREAD_RATE, INPUT_AREA (-1: no AREA_ADJUST) (11,*)
INTEGER, SAVE :: G0(3)=0,NG(3)=0   ! global cell index (0-based) of cell (1,1,1) of the current box, number of cells of the domain at its level
INTEGER, ALLOCATABLE, SAVE :: THD(:),THP(:),THT(:,:)   ! thin OBST (owner rule): direction (0 none), plane and transverse cell ranges at the current level
INTEGER, SAVE :: NTHIN = 0   ! thin-face walls found by the last builder pass (WALLDIAG)
REAL(EB), ALLOCATABLE, SAVE :: GROUGH(:,:)   !< ROUGHNESS of the level-0 ground wall of every level-0 column of the domain (GROUND_ROUGH_INIT)
REAL(EB), SAVE :: DXLEV(3,20) = 0._EB   ! cell size of fine level 1.. (set by FDS_FINE_LEVEL_CREATE)
INTEGER, ALLOCATABLE, SAVE :: OBI(:,:),VTS(:,:)   ! OBST and VENT XB snapped at the current box's dx (FDS NINT rounding)
INTEGER, ALLOCATABLE, SAVE :: OBS(:,:),VTI(:,:)    ! OBST SURF_INDEX(-3:3) (7,*); VENT IOR,SURF_INDEX,BOUNDARY_TYPE (3,*)

LOGICAL :: FWALL_T0 = .FALSE.   ! FDS_FINE_WALLS_QCON STAGE 0: FWALL_STATE keeps a thick wall's initial (parent) HTC and Q_CON_F
INTEGER :: FWALL_STAGE = -1   ! FDS_FINE_WALLS_QCON STAGE (0 t=0, 1 predictor, 2 corrector) while it runs: FWALL_STATE writes the OPEN ghost cells; -1 otherwise (diagnostic callers)
CONTAINS

SUBROUTINE FDS_FINE_SET_LISTS(NO,OXB,OSF,NV,VXB,VI) BIND(C,NAME='fds_fine_set_lists')
INTEGER(C_INT), VALUE :: NO,NV
REAL(C_DOUBLE), INTENT(IN) :: OXB(6,*),VXB(11,*)
INTEGER(C_INT), INTENT(IN) :: OSF(7,*),VI(3,*)
IF (ALLOCATED(OBX)) DEALLOCATE(OBX,OBS,VTX,VTI)
ALLOCATE(OBX(6,NO),OBS(7,NO),VTX(11,NV),VTI(3,NV))
NOB = NO ; NVT = NV
IF (NO>0) THEN ; OBX = OXB(:,1:NO) ; OBS = OSF(:,1:NO) ; ENDIF
IF (NV>0) THEN ; VTX = VXB(:,1:NV) ; VTI = VI(:,1:NV) ; ENDIF
END SUBROUTINE FDS_FINE_SET_LISTS

!> Make FINE_LEVEL(LEV) with NBOX boxes (LEV = the next level above the existing ones, or the top level again after FDS_FINE_LEVEL_DESTROY). Box IB has IJK(:,IB) cells
!> on [XB(1,IB),XB(2,IB)] x [XB(3,IB),XB(4,IB)] x [XB(5,IB),XB(6,IB)]; its mesh object is built only when LOCAL(IB)/=0 (the boxes of this rank). NM0 returns the offset of the
!> mesh numbers (box IB is mesh NM0+IB). Per-mesh-number kernel arrays are made large enough and the kernel wrappers are declared fine-ready.
!> IERR: 0 ok, 1 bad level (not the top level), 2 no fine-level support in this tree, 3 the level-0 meshes of this rank have more than one pressure zone (a fine box gets the single zone of level 0).
FUNCTION FDS_FINE_LEVEL_CREATE(LEV,NBOX,IJK,XB,LOCAL,NM0) BIND(C,NAME='fds_fine_level_create') RESULT(IERR)

INTEGER(C_INT), VALUE :: LEV,NBOX
INTEGER(C_INT), INTENT(IN) :: IJK(3,NBOX),LOCAL(NBOX)
REAL(C_DOUBLE), INTENT(IN) :: XB(6,NBOX)
INTEGER(C_INT), INTENT(OUT) :: NM0
INTEGER(C_INT) :: IERR
TYPE(FINE_LEVEL_TYPE), ALLOCATABLE :: TMP(:)
INTEGER :: L,IB,N0,NM,Z0,NMREF,I,J,K,KP
LOGICAL :: Z0_SET

IERR = 0 ; NM0 = 0
CALL GROUND_ROUGH_INIT()   ! collective, once
IF (LEV<1) THEN ; IERR = 1 ; RETURN ; ENDIF
IF (.NOT.ALLOCATED(FINE_LEVEL)) ALLOCATE(FINE_LEVEL(0))
IF (LEV<SIZE(FINE_LEVEL)) THEN ; IF (ANY(FINE_LEVEL(LEV+1:)%N_BOXES/=0)) THEN ; IERR = 1 ; RETURN ; ENDIF ; ENDIF           ! only the top level may be (re)made: the mesh numbers of the levels above would move
IF (LEV>SIZE(FINE_LEVEL)+1) THEN ; IERR = 1 ; RETURN ; ENDIF
IF (LEV==SIZE(FINE_LEVEL)+1) THEN
   ALLOCATE(TMP(LEV))
   DO L=1,SIZE(FINE_LEVEL)
      CALL MOVE_ALLOC(FINE_LEVEL(L)%BOX,TMP(L)%BOX) ; CALL MOVE_ALLOC(FINE_LEVEL(L)%VIEW,TMP(L)%VIEW)
      TMP(L)%NM0 = FINE_LEVEL(L)%NM0 ; TMP(L)%N_BOXES = FINE_LEVEL(L)%N_BOXES
   ENDDO
   CALL MOVE_ALLOC(TMP,FINE_LEVEL)
ENDIF
BUILD_VER = BUILD_VER+1
N0 = SIZE(MESHES)
DO L=1,LEV-1
   N0 = N0 + FINE_LEVEL(L)%N_BOXES
ENDDO
! pressure zone of the fine boxes: the one zone that all cells of the level-0 meshes of this rank have (a case with several zones is refused)
Z0 = 0 ; Z0_SET = .FALSE. ; NMREF = 0
DO NM=1,SIZE(MESHES)
   IF (.NOT.ALLOCATED(MESHES(NM)%PRESSURE_ZONE)) CYCLE
   IF (NMREF==0 .AND. ALLOCATED(MESHES(NM)%PBAR)) NMREF = NM
   ASSOCIATE(R=>MESHES(NM))
   DO K=1,R%KBAR ; DO J=1,R%JBAR ; DO I=1,R%IBAR   ! gas cells only: an OBST cell carries no zone of its own
      IF (R%CELL(R%CELL_INDEX(I,J,K))%SOLID) CYCLE
      IF (.NOT.Z0_SET) THEN ; Z0 = R%PRESSURE_ZONE(I,J,K) ; Z0_SET = .TRUE. ; ENDIF
   ENDDO ; ENDDO ; ENDDO
   END ASSOCIATE
ENDDO
FINE_LEVEL(LEV)%NM0 = N0 ; FINE_LEVEL(LEV)%N_BOXES = NBOX
IF (NBOX>0 .AND. LEV<=20) DXLEV(:,LEV) = (/XB(2,1)-XB(1,1),XB(4,1)-XB(3,1),XB(6,1)-XB(5,1)/)/REAL(IJK(:,1),EB)
IF (ALLOCATED(FINE_LEVEL(LEV)%BOX)) DEALLOCATE(FINE_LEVEL(LEV)%BOX)
IF (ALLOCATED(FINE_LEVEL(LEV)%VIEW)) DEALLOCATE(FINE_LEVEL(LEV)%VIEW)
ALLOCATE(FINE_LEVEL(LEV)%BOX(NBOX),FINE_LEVEL(LEV)%VIEW(NBOX))
DO IB=1,NBOX
   ASSOCIATE(B=>FINE_LEVEL(LEV)%BOX(IB))   ! extents of every box, also of the other ranks' (radi.f90 RAD_COVERED)
   B%XS=XB(1,IB) ; B%XF=XB(2,IB) ; B%YS=XB(3,IB) ; B%YF=XB(4,IB) ; B%ZS=XB(5,IB) ; B%ZF=XB(6,IB)
   END ASSOCIATE
   IF (LOCAL(IB)/=0) CALL BUILD_FINE_BOX(FINE_LEVEL(LEV)%BOX(IB),IJK(1,IB),IJK(2,IB),IJK(3,IB),XB(1,IB),XB(2,IB),XB(3,IB),XB(4,IB),XB(5,IB),XB(6,IB))
   IF (LOCAL(IB)/=0) THEN
      FINE_LEVEL(LEV)%BOX(IB)%PARTICLE_TAG = N0+IB   ! as init.f90 (PARTICLE_TAG = NM): unique tags, stepped by NMESHES in part.f90
      FINE_LEVEL(LEV)%BOX(IB)%PRESSURE_ZONE = Z0
      CALL FINE_BOX_ZONES(FINE_LEVEL(LEV)%BOX(IB))   ! D7: every cell inherits the zone of its level-0 cell (zones are never renumbered, so the chain l-1,...,0 ends there)
      ! zone-wide background pressure state (PBAR, PBAR_S, D_PBAR_DT, D_PBAR_DT_S, U_LEAK, R_PBAR) of a new box = the one of the level-0 mesh (without STRATIFICATION PBAR is uniform in K;
      ! with STRATIFICATION the box keeps its own P_0 profile and takes only the zone rates)
      IF (NMREF>0) THEN
         ASSOCIATE(M=>FINE_LEVEL(LEV)%BOX(IB),R=>MESHES(NMREF))
            M%D_PBAR_DT = R%D_PBAR_DT ; M%D_PBAR_DT_S = R%D_PBAR_DT_S ; M%U_LEAK = R%U_LEAK
            IF (.NOT.STRATIFICATION) THEN
               DO K=0,M%KBP1
                  M%PBAR(K,:) = R%PBAR(1,:) ; M%PBAR_S(K,:) = R%PBAR_S(1,:) ; M%R_PBAR(K,:) = R%R_PBAR(1,:)
               ENDDO
            ELSE   ! injected children carry the parent's rho and T: the background profile must be the parent's, piecewise constant per parent cell
               DO K=0,M%KBP1
                  KP = 1 ; DO WHILE (KP<R%KBAR .AND. M%ZC(K)>R%Z(KP)) ; KP = KP+1 ; ENDDO
                  M%P_0(K) = R%P_0(KP) ; M%TMP_0(K) = R%TMP_0(KP) ; M%RHO_0(K) = R%RHO_0(KP)
                  M%PBAR(K,:) = R%PBAR(KP,:) ; M%PBAR_S(K,:) = R%PBAR_S(KP,:) ; M%R_PBAR(K,:) = R%R_PBAR(KP,:)
               ENDDO
            ENDIF
         END ASSOCIATE
      ENDIF
   ENDIF
ENDDO
CALL FDS_K_RESERVE_BOXES(INT(N0+NBOX,C_INT))
CALL FDS_HOOK_SET_FINE_READY(1_C_INT)
NM0 = INT(N0,C_INT)

END FUNCTION FDS_FINE_LEVEL_CREATE

!> D7/D8: pressure zone of every cell of fine box M = the zone of the level-0 cell that contains its centre (a cell over an OBST of level 0 keeps the box's default zone).
!> Ghost cells copy the nearest interior cell.
SUBROUTINE FINE_BOX_ZONES(M)
USE TRAN, ONLY: GET_IJK
TYPE(MESH_TYPE), INTENT(INOUT) :: M
INTEGER :: I,J,K,L,NM0,I0,J0,K0
REAL(EB) :: XI,YJ,ZK
DO K=1,M%KBAR ; DO J=1,M%JBAR ; DO I=1,M%IBAR
   NM0 = 0
   DO L=1,SIZE(MESHES)
      IF (.NOT.ALLOCATED(MESHES(L)%PRESSURE_ZONE)) CYCLE
      IF (M%XC(I)>MESHES(L)%XS .AND. M%XC(I)<MESHES(L)%XF .AND. M%YC(J)>MESHES(L)%YS .AND. M%YC(J)<MESHES(L)%YF .AND. &
          M%ZC(K)>MESHES(L)%ZS .AND. M%ZC(K)<MESHES(L)%ZF) THEN ; NM0 = L ; EXIT ; ENDIF
   ENDDO
   IF (NM0==0) CYCLE
   CALL GET_IJK(M%XC(I),M%YC(J),M%ZC(K),NM0,XI,YJ,ZK,I0,J0,K0)
   IF (MESHES(NM0)%CELL(MESHES(NM0)%CELL_INDEX(I0,J0,K0))%SOLID) CYCLE
   M%PRESSURE_ZONE(I,J,K) = MESHES(NM0)%PRESSURE_ZONE(I0,J0,K0)
ENDDO ; ENDDO ; ENDDO
M%PRESSURE_ZONE(0,1:M%JBAR,1:M%KBAR) = M%PRESSURE_ZONE(1,1:M%JBAR,1:M%KBAR) ; M%PRESSURE_ZONE(M%IBP1,1:M%JBAR,1:M%KBAR) = M%PRESSURE_ZONE(M%IBAR,1:M%JBAR,1:M%KBAR)
M%PRESSURE_ZONE(:,0,:) = M%PRESSURE_ZONE(:,1,:) ; M%PRESSURE_ZONE(:,M%JBP1,:) = M%PRESSURE_ZONE(:,M%JBAR,:)
M%PRESSURE_ZONE(:,:,0) = M%PRESSURE_ZONE(:,:,1) ; M%PRESSURE_ZONE(:,:,M%KBP1) = M%PRESSURE_ZONE(:,:,M%KBAR)
END SUBROUTINE FINE_BOX_ZONES

!> Drop the boxes and views of the top level LEV (before it is made again after a regrid, or when it is cleared). IERR 1 when LEV is not the top level.
FUNCTION FDS_FINE_LEVEL_DESTROY(LEV) BIND(C,NAME='fds_fine_level_destroy') RESULT(IERR)

INTEGER(C_INT), VALUE :: LEV
INTEGER(C_INT) :: IERR

IERR = 0
IF (.NOT.ALLOCATED(FINE_LEVEL)) THEN ; IERR = 1 ; RETURN ; ENDIF
IF (LEV>SIZE(FINE_LEVEL) .OR. LEV<1) THEN ; IERR = 1 ; RETURN ; ENDIF
IF (ANY(FINE_LEVEL(LEV+1:)%N_BOXES/=0)) THEN ; IERR = 1 ; RETURN ; ENDIF   ! top bound level: the levels above were destroyed first (FINE_LEVEL keeps its size)
IF (ALLOCATED(FINE_LEVEL(LEV)%BOX)) DEALLOCATE(FINE_LEVEL(LEV)%BOX)
IF (ALLOCATED(FINE_LEVEL(LEV)%VIEW)) DEALLOCATE(FINE_LEVEL(LEV)%VIEW)
FINE_LEVEL(LEV)%N_BOXES = 0
CALL FWALL_DROP(LEV)

END FUNCTION FDS_FINE_LEVEL_DESTROY

!> CELL%SOLID of the interior cells of local box IB of level LEV from S (SideData comp 0 on valid+2): a solid cell points at the shared solid cell CELL(UBOUND).
SUBROUTINE FDS_FINE_SET_SOLID(LEV,IB,S) BIND(C,NAME='fds_fine_set_solid')
INTEGER(C_INT), VALUE :: LEV,IB
INTEGER(C_INT), INTENT(IN) :: S(-1:FINE_LEVEL(LEV)%BOX(IB)%IBAR+2,-1:FINE_LEVEL(LEV)%BOX(IB)%JBAR+2,-1:FINE_LEVEL(LEV)%BOX(IB)%KBAR+2)
ASSOCIATE(M=>FINE_LEVEL(LEV)%BOX(IB))
WHERE (S(1:M%IBAR,1:M%JBAR,1:M%KBAR)/=0) M%CELL_INDEX(1:M%IBAR,1:M%JBAR,1:M%KBAR) = UBOUND(M%CELL,1)
END ASSOCIATE
END SUBROUTINE FDS_FINE_SET_SOLID

!> Allocate the arrays of a fine-level mesh object of NI x NJ x NK cells on [XS,XF]x[YS,YF]x[ZS,ZF] (axis-aligned, uniform, Cartesian, no obstructions) with the contents the kernels read
!> (level-interface.md section 3): metrics (the arithmetic of READ_MESH for an untransformed grid), cell tables (ghost ring SOLID/EXTERIOR as READ_OBST sets it), empty edge and wall tables
!> (a fine box has no WALL cells: its non-domain faces are interface faces whose ghost values the driver fills, and there is no OMESH neighbour), and the state, field and work arrays of
!> INITIALIZE_MESH_VARIABLES_1 (initial values as there), including the Lagrangian particle storage of init.f90 (W7, PV-1: stock LAGRANGIAN_PARTICLE records; BOUNDARY_* storage grows
!> in ALLOCATE_STORAGE from the empty arrays). PARTICLE_TAG is set by the caller (mesh number of the box).
SUBROUTINE BUILD_FINE_BOX(M,NI,NJ,NK,XS,XF,YS,YF,ZS,ZF)

TYPE(MESH_TYPE), INTENT(INOUT) :: M
INTEGER, INTENT(IN) :: NI,NJ,NK
REAL(EB), INTENT(IN) :: XS,XF,YS,YF,ZS,ZF
INTEGER :: I,J,K,IC,NC,IBP1,JBP1,KBP1,N_LOWER_SCALARS

N_LOWER_SCALARS = MERGE(0,1,FLUX_LIMITER_MW_CORRECTION)   ! as in INITIALIZE_MESH_VARIABLES_1
! the largest box dimensions are sizes of automatic arrays in the kernels (VELOCITY_FLUX: GX(0:IBAR_MAX), GY, GZ): a fine box larger than every level-0 mesh must raise them
IBAR_MAX=MAX(IBAR_MAX,NI,NJ,NK) ; JBAR_MAX=MAX(JBAR_MAX,NJ) ; KBAR_MAX=MAX(KBAR_MAX,NK)
M%IBAR=NI ; M%JBAR=NJ ; M%KBAR=NK ; M%IBM1=NI-1 ; M%JBM1=NJ-1 ; M%KBM1=NK-1 ; M%IBP1=NI+1 ; M%JBP1=NJ+1 ; M%KBP1=NK+1
IBP1=M%IBP1 ; JBP1=M%JBP1 ; KBP1=M%KBP1
M%XS=XS ; M%XF=XF ; M%YS=YS ; M%YF=YF ; M%ZS=ZS ; M%ZF=ZF
M%DXI=(XF-XS)/REAL(NI,EB) ; M%DETA=(YF-YS)/REAL(NJ,EB) ; M%DZETA=(ZF-ZS)/REAL(NK,EB)
M%RDXI=1._EB/M%DXI ; M%RDETA=1._EB/M%DETA ; M%RDZETA=1._EB/M%DZETA
ALLOCATE(M%X(0:NI),M%XC(0:NI+1),M%HX(0:NI+1),M%DX(0:NI+1),M%RDX(0:NI+1),M%DXN(0:NI),M%RDXN(0:NI),M%R(0:NI),M%RC(0:NI+1),M%RRN(0:NI+1))
ALLOCATE(M%Y(0:NJ),M%YC(0:NJ+1),M%HY(0:NJ+1),M%DY(0:NJ+1),M%RDY(0:NJ+1),M%DYN(0:NJ),M%RDYN(0:NJ))
ALLOCATE(M%Z(0:NK),M%ZC(0:NK+1),M%HZ(0:NK+1),M%DZ(0:NK+1),M%RDZ(0:NK+1),M%DZN(0:NK),M%RDZN(0:NK))
M%R=1._EB ; M%RC=1._EB ; M%RRN=1._EB
CALL AXIS(NI,XS,XF,M%DXI,M%X,M%XC,M%HX,M%DX,M%RDX,M%DXN,M%RDXN,M%DXMIN)
CALL AXIS(NJ,YS,YF,M%DETA,M%Y,M%YC,M%HY,M%DY,M%RDY,M%DYN,M%RDYN,M%DYMIN)
CALL AXIS(NK,ZS,ZF,M%DZETA,M%Z,M%ZC,M%HZ,M%DZ,M%RDZ,M%DZN,M%RDZN,M%DZMIN)

! cell tables: the ghost ring gets cell numbers (READ_OBST order), is SOLID and EXTERIOR; interior cells use index 0 (CELL(0), non-solid)
ALLOCATE(M%CELL_INDEX(0:IBP1,0:JBP1,0:KBP1)) ; M%CELL_INDEX = 0
NC = 0
DO K=0,KBP1
   DO J=0,JBP1
      DO I=0,IBP1
         IF (I<=1 .OR. I>=NI .OR. J<=1 .OR. J>=NJ .OR. K<=1 .OR. K>=NK) THEN ; NC = NC+1 ; ENDIF
      ENDDO
   ENDDO
ENDDO
CALL NUMBER_CELLS(M,NI,NJ,NK,NC)
ALLOCATE(M%CELL(0:NC+1)) ; M%CELL(NC+1)%SOLID=.TRUE.   ! CELL(NC+1): the shared interior solid cell (FDS_FINE_SET_SOLID)
DO K=0,KBP1
   DO J=0,JBP1
      DO I=0,IBP1
         IC = M%CELL_INDEX(I,J,K)
         IF (IC>0) THEN
            M%CELL(IC)%I=I ; M%CELL(IC)%J=J ; M%CELL(IC)%K=K
            IF (I==0 .OR. I==IBP1 .OR. J==0 .OR. J==JBP1 .OR. K==0 .OR. K==KBP1) THEN ; M%CELL(IC)%SOLID=.TRUE. ; M%CELL(IC)%EXTERIOR=.TRUE. ; ENDIF
         ENDIF
      ENDDO
   ENDDO
ENDDO
ALLOCATE(M%INTERPOLATED_MESH(1:NI,1:NJ,1:NK)) ; M%INTERPOLATED_MESH = 0     ! no cell of a fine box is covered by a finer mesh (yet)
ALLOCATE(M%EDGE(0:0)) ; ALLOCATE(M%WALL(0:0),M%EXTERNAL_WALL(0),M%THIN_WALL(0:10))
M%WALL(0)%BOUNDARY_TYPE = NULL_BOUNDARY ; M%WALL(0)%SURF_INDEX = DEFAULT_SURF_INDEX
M%N_WALL_CELLS=0 ; M%N_WALL_CELLS_DIM=0 ; M%N_INTERNAL_WALL_CELLS=0 ; M%N_EXTERNAL_WALL_CELLS=0
ALLOCATE(M%UVW_SAVE(0),M%U_GHOST(0),M%V_GHOST(0),M%W_GHOST(0),M%WALL_WORK1(0),M%WALL_WORK2(0))

! state, field and work arrays (init.f90, INITIALIZE_MESH_VARIABLES_1)
ALLOCATE(M%TMP(-1:IBP1+1,-1:JBP1+1,-1:KBP1+1)) ; M%TMP = TMPA
ALLOCATE(M%RHO(-1:IBP1+1,-1:JBP1+1,-1:KBP1+1)) ; M%RHO = RHOA
ALLOCATE(M%RHOS(-1:IBP1+1,-1:JBP1+1,-1:KBP1+1)) ; M%RHOS = RHOA
ALLOCATE(M%ZZ(-1:IBP1+1,-1:JBP1+1,-1:KBP1+1,N_TOTAL_SCALARS)) ; M%ZZ = 0._EB
ALLOCATE(M%ZZS(-1:IBP1+1,-1:JBP1+1,-1:KBP1+1,N_TOTAL_SCALARS)) ; M%ZZS = 0._EB
ALLOCATE(M%WORK_PAD(-1:IBP1+1,-1:JBP1+1,-1:KBP1+1)) ; M%WORK_PAD = 0._EB
ALLOCATE(M%U(-1:IBP1,0:JBP1,0:KBP1),M%US(-1:IBP1,0:JBP1,0:KBP1),M%WORK_U(-1:IBP1,0:JBP1,0:KBP1)) ; M%U=0._EB ; M%US=0._EB ; M%WORK_U=0._EB
ALLOCATE(M%V(0:IBP1,-1:JBP1,0:KBP1),M%VS(0:IBP1,-1:JBP1,0:KBP1),M%WORK_V(0:IBP1,-1:JBP1,0:KBP1)) ; M%V=0._EB ; M%VS=0._EB ; M%WORK_V=0._EB
ALLOCATE(M%W(0:IBP1,0:JBP1,-1:KBP1),M%WS(0:IBP1,0:JBP1,-1:KBP1),M%WORK_W(0:IBP1,0:JBP1,-1:KBP1)) ; M%W=0._EB ; M%WS=0._EB ; M%WORK_W=0._EB
ALLOCATE(M%FVX(0:IBP1,0:JBP1,0:KBP1),M%FVY(0:IBP1,0:JBP1,0:KBP1),M%FVZ(0:IBP1,0:JBP1,0:KBP1)) ; M%FVX=0._EB ; M%FVY=0._EB ; M%FVZ=0._EB
ALLOCATE(M%FVX_B(0:IBP1,0:JBP1,0:KBP1),M%FVY_B(0:IBP1,0:JBP1,0:KBP1),M%FVZ_B(0:IBP1,0:JBP1,0:KBP1)) ; M%FVX_B=0._EB ; M%FVY_B=0._EB ; M%FVZ_B=0._EB
IF (PARTICLE_DRAG) THEN
   ALLOCATE(M%FVX_D(0:IBP1,0:JBP1,0:KBP1),M%FVY_D(0:IBP1,0:JBP1,0:KBP1),M%FVZ_D(0:IBP1,0:JBP1,0:KBP1)) ; M%FVX_D=0._EB ; M%FVY_D=0._EB ; M%FVZ_D=0._EB
ENDIF
IF (N_LP_ARRAY_INDICES>0) THEN   ! as init.f90 (particle radiation and averaged droplet arrays)
   ALLOCATE(M%QR_W(0:IBP1,0:JBP1,0:KBP1)) ; M%QR_W=0._EB
   ALLOCATE(M%AVG_DROP_DEN(0:IBP1,0:JBP1,0:KBP1,N_LP_ARRAY_INDICES),M%AVG_DROP_AREA(0:IBP1,0:JBP1,0:KBP1,N_LP_ARRAY_INDICES))
   ALLOCATE(M%AVG_DROP_TMP(0:IBP1,0:JBP1,0:KBP1,N_LP_ARRAY_INDICES),M%AVG_DROP_RAD(0:IBP1,0:JBP1,0:KBP1,N_LP_ARRAY_INDICES))
   M%AVG_DROP_DEN=0._EB ; M%AVG_DROP_AREA=0._EB ; M%AVG_DROP_TMP=TMPM ; M%AVG_DROP_RAD=0._EB
ENDIF
ALLOCATE(M%H(0:IBP1,0:JBP1,0:KBP1),M%HS(0:IBP1,0:JBP1,0:KBP1)) ; M%H=0.5_EB*(U0**2+V0**2+W0**2) ; M%HS=0.5_EB*(U0**2+V0**2+W0**2)
ALLOCATE(M%KRES(0:IBP1,0:JBP1,0:KBP1),M%DDDT(0:IBP1,0:JBP1,0:KBP1),M%D(0:IBP1,0:JBP1,0:KBP1),M%DS(0:IBP1,0:JBP1,0:KBP1)) ; M%KRES=0._EB ; M%DDDT=0._EB ; M%D=0._EB ; M%DS=0._EB
ALLOCATE(M%MU(0:IBP1,0:JBP1,0:KBP1),M%MU_DNS(0:IBP1,0:JBP1,0:KBP1)) ; M%MU=0._EB ; M%MU_DNS=0._EB
IF (CHECK_VN) THEN ; ALLOCATE(M%D_Z_MAX(0:IBP1,0:JBP1,0:KBP1)) ; M%D_Z_MAX=0._EB ; ENDIF
ALLOCATE(M%STRAIN_RATE(0:IBP1,0:JBP1,0:KBP1)) ; M%STRAIN_RATE=0._EB
ALLOCATE(M%CSD2(0:IBP1,0:JBP1,0:KBP1)) ; M%CSD2=0._EB
ALLOCATE(M%Q(0:IBP1,0:JBP1,0:KBP1)) ; M%Q=0._EB
ALLOCATE(M%MIX_TIME(0:IBP1,0:JBP1,0:KBP1)) ; M%MIX_TIME=0._EB
ALLOCATE(M%PRESSURE_ZONE(0:IBP1,0:JBP1,0:KBP1)) ; M%PRESSURE_ZONE=0   ! single (ambient) pressure zone; FDS_FINE_LEVEL_CREATE refuses cases with N_ZONE>0
ALLOCATE(M%DEL_RHO_D_DEL_Z(0:IBP1,0:JBP1,0:KBP1,N_TOTAL_SCALARS)) ; M%DEL_RHO_D_DEL_Z=0._EB
ALLOCATE(M%RSUM(0:IBP1,0:JBP1,0:KBP1)) ; M%RSUM=RSUM0
ALLOCATE(M%FX(0:IBP1,0:JBP1,0:KBP1,N_LOWER_SCALARS:N_TOTAL_SCALARS),M%FY(0:IBP1,0:JBP1,0:KBP1,N_LOWER_SCALARS:N_TOTAL_SCALARS),M%FZ(0:IBP1,0:JBP1,0:KBP1,N_LOWER_SCALARS:N_TOTAL_SCALARS))
M%FX=0._EB ; M%FY=0._EB ; M%FZ=0._EB
IF (STORE_SPECIES_FLUX) THEN
   ALLOCATE(M%ADV_FX(0:IBP1,0:JBP1,0:KBP1,1:N_TOTAL_SCALARS),M%ADV_FY(0:IBP1,0:JBP1,0:KBP1,1:N_TOTAL_SCALARS),M%ADV_FZ(0:IBP1,0:JBP1,0:KBP1,1:N_TOTAL_SCALARS))
   ALLOCATE(M%DIF_FX(0:IBP1,0:JBP1,0:KBP1,1:N_TOTAL_SCALARS),M%DIF_FY(0:IBP1,0:JBP1,0:KBP1,1:N_TOTAL_SCALARS),M%DIF_FZ(0:IBP1,0:JBP1,0:KBP1,1:N_TOTAL_SCALARS))
   ALLOCATE(M%DIF_FXS(0:IBP1,0:JBP1,0:KBP1,1:N_TOTAL_SCALARS),M%DIF_FYS(0:IBP1,0:JBP1,0:KBP1,1:N_TOTAL_SCALARS),M%DIF_FZS(0:IBP1,0:JBP1,0:KBP1,1:N_TOTAL_SCALARS))
   M%ADV_FX=0._EB ; M%ADV_FY=0._EB ; M%ADV_FZ=0._EB ; M%DIF_FX=0._EB ; M%DIF_FY=0._EB ; M%DIF_FZ=0._EB ; M%DIF_FXS=0._EB ; M%DIF_FYS=0._EB ; M%DIF_FZS=0._EB
ENDIF
IF (STORE_PRESSURE_POISSON_RESIDUAL) THEN ; ALLOCATE(M%PP_RESIDUAL(0:IBP1,0:JBP1,0:KBP1)) ; M%PP_RESIDUAL=0._EB ; ENDIF
IF (STORE_CUTCELL_DIVERGENCE) THEN ; ALLOCATE(M%CCVELDIV(0:IBP1,0:JBP1,0:KBP1)) ; M%CCVELDIV=0._EB ; ENDIF
IF (STORE_CARTESIAN_DIVERGENCE) THEN ; ALLOCATE(M%CARTVELDIV(0:IBP1,0:JBP1,0:KBP1)) ; M%CARTVELDIV=0._EB ; ENDIF
IF (N_REACTIONS>0 .OR. N_LP_ARRAY_INDICES>0) THEN
   ALLOCATE(M%D_SOURCE(0:IBP1,0:JBP1,0:KBP1)) ; M%D_SOURCE=0._EB
   ALLOCATE(M%M_DOT_PPP(0:IBP1,0:JBP1,0:KBP1,1:N_TRACKED_SPECIES)) ; M%M_DOT_PPP=0._EB
ENDIF
ALLOCATE(M%CHEM_ACTIVE_CELLS(IBP1*JBP1*KBP1,4)) ; M%CHEM_ACTIVE_CELLS=-1   ! as init.f90: COMBUSTION_LOAD_BALANCED on a fine box (fire.f90, patch 0020)
IF (REAC_SOURCE_CHECK) THEN
   ALLOCATE(M%REAC_SOURCE_TERM(0:IBP1,0:JBP1,0:KBP1,N_TRACKED_SPECIES),M%Q_REAC(0:IBP1,0:JBP1,0:KBP1,N_REACTIONS)) ; M%REAC_SOURCE_TERM=0._EB ; M%Q_REAC=0._EB
ENDIF
ALLOCATE(M%CHI_R(0:IBP1,0:JBP1,0:KBP1),M%QR(0:IBP1,0:JBP1,0:KBP1),M%KAPPA_GAS(0:IBP1,0:JBP1,0:KBP1),M%UII(0:IBP1,0:JBP1,0:KBP1)) ; M%CHI_R=0._EB ; M%QR=0._EB ; M%KAPPA_GAS=0._EB ; M%UII=4._EB*SIGMA*TMPA4
M%NLP = 0 ; M%RAD_CALL_COUNTER = MESHES(LOWER_MESH_INDEX)%RAD_CALL_COUNTER ; M%ANGLE_INC_COUNTER = MESHES(LOWER_MESH_INDEX)%ANGLE_INC_COUNTER ; M%RAD_FRESH = .TRUE.   ! in step with a local level-0 mesh
IF (RADIATION) THEN ; ALLOCATE(M%UIID(0:IBP1,0:JBP1,0:KBP1,1:UIIDIM),M%CELL_ILW(1:UBOUND(M%CELL,1),1:3),M%OMESH(NMESHES)) ; M%UIID=4._EB*SIGMA*TMPA4/REAL(UIIDIM,EB) ; ENDIF   ! as init.f90 (radiation on a fine box, patch 0022)
ALLOCATE(M%WORK1(0:IBP1,0:JBP1,0:KBP1),M%WORK2(0:IBP1,0:JBP1,0:KBP1),M%WORK3(0:IBP1,0:JBP1,0:KBP1),M%WORK4(0:IBP1,0:JBP1,0:KBP1),M%WORK5(0:IBP1,0:JBP1,0:KBP1))
ALLOCATE(M%WORK6(0:IBP1,0:JBP1,0:KBP1),M%WORK7(0:IBP1,0:JBP1,0:KBP1),M%WORK8(0:IBP1,0:JBP1,0:KBP1),M%WORK9(0:IBP1,0:JBP1,0:KBP1),M%IWORK1(0:IBP1,0:JBP1,0:KBP1))
ALLOCATE(M%SWORK1(0:IBP1,0:JBP1,0:KBP1,N_LOWER_SCALARS:N_TOTAL_SCALARS),M%SWORK2(0:IBP1,0:JBP1,0:KBP1,N_LOWER_SCALARS:N_TOTAL_SCALARS))
ALLOCATE(M%SWORK3(0:IBP1,0:JBP1,0:KBP1,N_LOWER_SCALARS:N_TOTAL_SCALARS),M%SWORK4(0:IBP1,0:JBP1,0:KBP1,1:N_TOTAL_SCALARS))
M%WORK1=0._EB ; M%WORK2=0._EB ; M%WORK3=0._EB ; M%WORK4=0._EB ; M%WORK5=0._EB ; M%WORK6=0._EB ; M%WORK7=0._EB ; M%WORK8=0._EB ; M%WORK9=0._EB ; M%IWORK1=0
M%SWORK1=0._EB ; M%SWORK2=0._EB ; M%SWORK3=0._EB ; M%SWORK4=0._EB
ALLOCATE(M%U_WIND(0:KBP1),M%V_WIND(0:KBP1),M%W_WIND(0:KBP1))
DO K=0,KBP1   ! as COMPUTE_WIND_COMPONENTS: the prevailing-wind profile at the box's own cell heights
   CALL GET_WIND_AT_HEIGHT(T_BEGIN,M%ZC(K),M%U_WIND(K),M%V_WIND(K),M%W_WIND(K))
ENDDO
ALLOCATE(M%LES_FILTER_WIDTH(0:IBP1,0:JBP1,0:KBP1)) ; M%LES_FILTER_WIDTH=0._EB
DO K=0,KBP1 ; DO J=0,JBP1 ; DO I=0,IBP1   ! as INITIALIZE_MESH_VARIABLES_1: filter width and the constant-Smagorinsky coefficient
   M%LES_FILTER_WIDTH(I,J,K) = LES_FILTER_WIDTH_FUNCTION(M%DX(I),M%DY(J),M%DZ(K))
   IF (TURB_MODEL==CONSMAG .OR. TURB_MODEL==DYNSMAG) M%CSD2(I,J,K) = (C_SMAGORINSKY*M%LES_FILTER_WIDTH(I,J,K))**2
ENDDO ; ENDDO ; ENDDO

! background atmosphere (INITIALIZE_ATMOSPHERE, cell-centred part) and background pressure (PBAR, PBAR_S, R_PBAR, D_PBAR_DT, D_PBAR_DT_S, U_LEAK)
ALLOCATE(M%P_0(0:KBP1),M%TMP_0(0:KBP1),M%RHO_0(0:KBP1))
IF (STRATIFICATION) THEN
   DO K=0,KBP1
      M%TMP_0(K) = TMPA*EVALUATE_RAMP(M%ZC(K),I_RAMP_TMP0_Z)
      M%P_0(K)   = EVALUATE_RAMP(M%ZC(K),I_RAMP_P0_Z)
      M%RHO_0(K) = M%P_0(K)/(M%TMP_0(K)*RSUM0)
   ENDDO
ELSE
   M%TMP_0(:) = TMPA ; M%P_0(:) = P_INF ; M%RHO_0(:) = M%P_0(:)/(M%TMP_0(:)*RSUM0)
ENDIF
ALLOCATE(M%PBAR(0:KBP1,0:N_ZONE),M%PBAR_S(0:KBP1,0:N_ZONE),M%R_PBAR(0:KBP1,0:N_ZONE),M%D_PBAR_DT(0:N_ZONE),M%D_PBAR_DT_S(0:N_ZONE),M%U_LEAK(0:N_ZONE))
M%D_PBAR_DT=0._EB ; M%D_PBAR_DT_S=0._EB ; M%U_LEAK=0._EB ; M%R_PBAR=0._EB
DO K=0,KBP1 ; M%PBAR(K,:)=M%P_0(K) ; M%PBAR_S(K,:)=M%P_0(K) ; ENDDO
ALLOCATE(M%BOUNDARY_COORD(0),M%BOUNDARY_PROP1(0),M%BOUNDARY_PROP2(0),M%BOUNDARY_RADIA(0))
M%NLP=0 ; M%NLPDIM=50   ! particle storage as init.f90 (INITIALIZE_MESH_VARIABLES_1)
IF (N_LAGRANGIAN_CLASSES>0) THEN ; ALLOCATE(M%PARTICLE_LAST(N_LAGRANGIAN_CLASSES)) ; M%PARTICLE_LAST=0 ; ENDIF
IF (PARTICLE_FILE) ALLOCATE(M%LAGRANGIAN_PARTICLE(M%NLPDIM))
M%CLIP_RHOMIN=.FALSE. ; M%CLIP_RHOMAX=.FALSE. ; M%DT_RESTRICT_COUNT=0

END SUBROUTINE BUILD_FINE_BOX

!> One axis of the grid exactly as READ_MESH makes it for an untransformed uniform mesh (G(X)=X).
SUBROUTINE AXIS(N,S,E,DXI,X,XC,HX,DX,RDX,DXN,RDXN,DXMIN)

INTEGER, INTENT(IN) :: N
REAL(EB), INTENT(IN) :: S,E,DXI
REAL(EB), INTENT(OUT) :: X(0:N),XC(0:N+1),HX(0:N+1),DX(0:N+1),RDX(0:N+1),DXN(0:N),RDXN(0:N),DXMIN
INTEGER :: I

DXMIN = 1000000._EB ; DX = 0._EB ; RDX = 0._EB ; HX = 1._EB ; XC = 0._EB
DO I=0,N
   X(I) = S + REAL(I,EB)*DXI
   IF (I>0) THEN
      DX(I) = X(I) - X(I-1)
      DXMIN = MIN(DXMIN,DX(I))
      RDX(I) = 1._EB/DX(I)
   ENDIF
ENDDO
X(0) = S ; X(N) = E
DX(0) = DX(1) ; DX(N+1) = DX(N) ; RDX(0) = 1._EB/DX(1) ; RDX(N+1) = 1._EB/DX(N)
DO I=0,N
   DXN(I) = 0.5_EB*(DX(I)+DX(I+1))
   RDXN(I) = 1._EB/DXN(I)
   IF (I>0) XC(I) = 0.5_EB*(X(I)+X(I-1))
ENDDO
XC(0) = S - 0.5_EB*DX(0) ; XC(N+1) = E + 0.5_EB*DX(N+1)

END SUBROUTINE AXIS

!> The cell numbering of the ghost ring in the order of READ_OBST (low/high x layers, then y, then z).
SUBROUTINE NUMBER_CELLS(M,NI,NJ,NK,NC)

TYPE(MESH_TYPE), INTENT(INOUT) :: M
INTEGER, INTENT(IN) :: NI,NJ,NK
INTEGER, INTENT(INOUT) :: NC
INTEGER :: I,J,K,CNT

CNT = 0
DO K=0,NK+1
   DO J=0,NJ+1
      DO I=0,1 ; CALL NUM(I,J,K) ; ENDDO
      DO I=NI,NI+1 ; CALL NUM(I,J,K) ; ENDDO
   ENDDO
ENDDO
DO K=0,NK+1
   DO I=0,NI+1
      DO J=0,1 ; CALL NUM(I,J,K) ; ENDDO
      DO J=NJ,NJ+1 ; CALL NUM(I,J,K) ; ENDDO
   ENDDO
ENDDO
DO J=0,NJ+1
   DO I=0,NI+1
      DO K=0,1 ; CALL NUM(I,J,K) ; ENDDO
      DO K=NK,NK+1 ; CALL NUM(I,J,K) ; ENDDO
   ENDDO
ENDDO
NC = CNT

CONTAINS
   SUBROUTINE NUM(I,J,K)
   INTEGER, INTENT(IN) :: I,J,K
   IF (M%CELL_INDEX(I,J,K)==0) THEN ; CNT = CNT+1 ; M%CELL_INDEX(I,J,K) = CNT ; ENDIF
   END SUBROUTINE NUM

END SUBROUTINE NUMBER_CELLS
!> Drop the hidden walls of level LEV (with the boxes, FDS_FINE_LEVEL_DESTROY).
SUBROUTINE FWALL_DROP(LEV)
INTEGER, INTENT(IN) :: LEV
INTEGER :: IB
BUILD_VER = BUILD_VER+1
IF (.NOT.ALLOCATED(FWALL)) RETURN
IF (LEV>SIZE(FWALL)) RETURN
IF (.NOT.ALLOCATED(FWALL(LEV)%B)) RETURN
DO IB=1,SIZE(FWALL(LEV)%B)
   IF (ASSOCIATED(FWALL(LEV)%B(IB)%W)) DEALLOCATE(FWALL(LEV)%B(IB)%W)
ENDDO
DEALLOCATE(FWALL(LEV)%B)
END SUBROUTINE FWALL_DROP

!> Build the hidden wall cells of local box IB of level LEV (after FDS_FINE_SET_SOLID; again at every regrid). A fine gas cell gets a wall on its face in direction D (+-1,+-2,+-3)
!> when that face lies on the D-face of its level-0 parent gas cell and the parent has a SOLID or OPEN wall there (domain boundary or OBST face). Coarse-fine faces have no
!> parent wall, so they get none. The wall takes the parent's surface; a thermally thick (stateful) surface on a fine box stops the run. Storage via ALLOCATE_STORAGE (MESH_OBJ).
SUBROUTINE FDS_FINE_WALLS_BUILD(LEV,IB) BIND(C,NAME='fds_fine_walls_build')
USE MESH_VARIABLES, ONLY: MESH_OBJ
USE TYPES, ONLY: SURFACE,WALL_TYPE,BOUNDARY_COORD_TYPE
USE GLOBAL_CONSTANTS, ONLY: SOLID_BOUNDARY,OPEN_BOUNDARY,MIRROR_BOUNDARY,THERMALLY_THICK,LOGLAW_HTC_MODEL,INERT_SURF_INDEX,T_BEGIN
USE MEMORY_FUNCTIONS, ONLY: ALLOCATE_STORAGE
USE TRAN, ONLY: GET_IJK
USE GLOBAL_CONSTANTS, ONLY: MY_RANK,PROCESS
INTEGER(C_INT), VALUE :: LEV,IB
INTEGER :: I,J,K,D,NM0,I0,J0,K0,IW0,IW,N,NM,L,N_SKIP,SI,BT,NMP,NN,GK(3),RNM0,RIW0
LOGICAL :: LOCAL0,WALLQ
INTEGER, ALLOCATABLE :: T(:,:),RT(:,:)
REAL(EB) :: XI,YJ,ZK,FF,F0,TOL,FA,DXY(3)
TYPE(MESH_TYPE), POINTER :: R,WM,M0
TYPE(WALL_TYPE), POINTER :: WC,WC0
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC

BUILD_VER = BUILD_VER+1
IF (.NOT.ALLOCATED(FWALL)) ALLOCATE(FWALL(0))
IF (SIZE(FWALL)<LEV) THEN
   BLOCK
   TYPE(FWALL_LEVEL_TYPE), ALLOCATABLE :: TMP(:)
   ALLOCATE(TMP(LEV))
   DO L=1,SIZE(FWALL) ; IF (ALLOCATED(FWALL(L)%B)) CALL MOVE_ALLOC(FWALL(L)%B,TMP(L)%B) ; ENDDO
   CALL MOVE_ALLOC(TMP,FWALL)
   END BLOCK
ENDIF
IF (.NOT.ALLOCATED(FWALL(LEV)%B)) ALLOCATE(FWALL(LEV)%B(FINE_LEVEL(LEV)%N_BOXES))
FINE_DEVC_PROC => FINE_DEVC   ! fine wall-flux devices (patch 0027 hook in FDS_SETUP(3))
LEAF_WALL_COUNT_PROC => LEAF_WALL_COUNT ; LEAF_WALL_PROC => LEAF_WALL ; LEAF_WALL_COVERED_PROC => LEAF_WALL_COVERED ; LEAF_BOXWALL_COUNT_PROC => LEAF_BOXWALL_COUNT ; LEAF_BOXWALL_PROC => LEAF_BOXWALL ; RAD_WALLS_PROC => FRAD_WALLS ; FINE_WALLS_PTR_PROC => FWALL_POINTERS   ! fine walls in the wall-device statistics
R => FINE_LEVEL(LEV)%BOX(IB)
NM = FINE_LEVEL(LEV)%NM0 + IB
! Per-level wall inventory (same code for every level >= 1): the OBST and VENT lists snapped at this box's own dx (FDS rounding NINT((X-X0)/DX), as ObstMask.H), no parent copy.
! Gas cell (I,J,K) gets a wall on face D when the cell across it is an OBST cell or outside the domain; coarse-fine / box-box faces get none. Surface: the last VENT of the list that covers
! the face (plane and extents snapped at this dx, IOR = -D), else the OBST face surface (last OBST containing the cell across), else SURF_DEFAULT. PNM/PIW = the level-0 wall with the same
! surface whose face contains this one (0 if none: new fine faces and faces whose surface differs; they have no parent state).
CALL SNAP_SETUP(R,LEV)
NTHIN = 0
ALLOCATE(T(12,6*R%IBAR*R%JBAR*R%KBAR),RT(2,6*R%IBAR*R%JBAR*R%KBAR)) ; N = 0 ; N_SKIP = 0 ; RT = 0   ! I,J,K,D,NM0,IW0,SURF_INDEX,BOUNDARY_TYPE of each wall, global level-0 cell of a parent on another rank
DO K=1,R%KBAR ; DO J=1,R%JBAR ; DO I=1,R%IBAR
   IF (R%CELL(R%CELL_INDEX(I,J,K))%SOLID) CYCLE
   NM0 = 0
   DO L=1,MIN(NMESHES,SIZE(MESHES))   ! level-0 meshes of every rank (extents are known everywhere; the WALL data are not)
      IF (R%XC(I)>MESHES(L)%XS .AND. R%XC(I)<MESHES(L)%XF .AND. R%YC(J)>MESHES(L)%YS .AND. R%YC(J)<MESHES(L)%YF .AND. &
          R%ZC(K)>MESHES(L)%ZS .AND. R%ZC(K)<MESHES(L)%ZF) THEN ; NM0 = L ; EXIT ; ENDIF
   ENDDO
   IF (NM0==0) THEN ; WRITE(0,'(A)') 'fds_amr ERROR in FDS_FINE_WALLS_BUILD: no level-0 mesh under a fine cell' ; ERROR STOP 1 ; ENDIF
   M0 => MESHES(NM0)
   ! The parent mesh may belong to another rank (the fine-level boxes follow the AMReX distribution, not the level-0 owner): its WALL array is not stored here, so such a wall has no parent state here.
   LOCAL0 = PROCESS(NM0)==MY_RANK
   IF (LOCAL0) CALL GET_IJK(R%XC(I),R%YC(J),R%ZC(K),NM0,XI,YJ,ZK,I0,J0,K0)
   DO D=-3,3
      IF (D==0) CYCLE
      CALL FACE_SURFACE(R,I,J,K,D,SI,BT,WALLQ)
      IF (.NOT.WALLQ) CYCLE
      IF (BT/=SOLID_BOUNDARY .AND. BT/=OPEN_BOUNDARY .AND. BT/=MIRROR_BOUNDARY) CYCLE   ! MIRROR: the radiation reflects there (no wall = ambient intensity)
      IF (SURFACE(SI)%TMP_GAS_FRONT>0._EB .OR. &
          ANY(SURFACE(SI)%VEL_T/=0._EB) .OR. SURFACE(SI)%HEAT_TRANSFER_MODEL==LOGLAW_HTC_MODEL) THEN
         WRITE(0,'(3A,I0,A)') 'fds_amr ERROR in FDS_FINE_WALLS_BUILD: unsupported surface ',TRIM(SURFACE(SI)%ID),' on level ',LEV, &
                         ' (levels above 0 support inert / fixed-temperature surfaces without VEL_T, TMP_GAS_FRONT or the log-law HTC model only)'
         ERROR STOP 1
      ENDIF
      IW0 = 0 ; NMP = 0 ; RNM0 = 0 ; RIW0 = 0
      IF (LOCAL0) THEN
         SELECT CASE(D)
            CASE(-1) ; FF=R%X(I-1) ; F0=M0%X(I0-1) ; TOL=R%DXI
            CASE( 1) ; FF=R%X(I)   ; F0=M0%X(I0)   ; TOL=R%DXI
            CASE(-2) ; FF=R%Y(J-1) ; F0=M0%Y(J0-1) ; TOL=R%DETA
            CASE( 2) ; FF=R%Y(J)   ; F0=M0%Y(J0)   ; TOL=R%DETA
            CASE(-3) ; FF=R%Z(K-1) ; F0=M0%Z(K0-1) ; TOL=R%DZETA
            CASE( 3) ; FF=R%Z(K)   ; F0=M0%Z(K0)   ; TOL=R%DZETA
         END SELECT
         IF (ABS(FF-F0)<=1.E-3_EB*TOL .AND. .NOT.M0%CELL(M0%CELL_INDEX(I0,J0,K0))%SOLID) THEN
            IW0 = M0%CELL(M0%CELL_INDEX(I0,J0,K0))%WALL_INDEX(D)
            IF (IW0>0) THEN
               IF (M0%WALL(IW0)%SURF_INDEX==SI .AND. M0%WALL(IW0)%BOUNDARY_TYPE==BT) NMP = NM0
            ENDIF
            IF (NMP==0) IW0 = 0
            IF (NMP==0 .AND. BT==SOLID_BOUNDARY) CALL REMAP_TARGET(M0,I0,J0,K0,D,SI,BT,R%XC(I),R%YC(J),R%ZC(K),RNM0,RIW0)   ! ruling 2 interim
         ENDIF
      ENDIF
      GK = -1
      IF (.NOT.LOCAL0) CALL REMOTE_PARENT_CELL(R,I,J,K,D,NM0,GK)
      IF (.NOT.LOCAL0) N_SKIP = N_SKIP+1   ! a hidden wall over another rank's mesh (cells without walls there are allowed: test FDS_AMR_DEBUG_FOREIGN_BOXES)
      N = N+1 ; T(:,N) = (/I,J,K,D,NMP,IW0,SI,BT,GK,FVT/)
      IF (RNM0>0) THEN ; RT(1,N) = NM0 ; RT(2,N) = RIW0 ; ENDIF
   ENDDO
ENDDO ; ENDDO ; ENDDO

IF (N_SKIP>0) THEN   ! hidden walls: cannot happen, fine boxes are cut at the level-0 mesh boundaries and owned by the rank of their mesh (RegridAmrCore::chop_at_level0, TwoLevelRun owner_dm)
   WRITE(0,'(A,I0,A,I0,A,I0,A)') 'fds_amr ERROR in FDS_FINE_WALLS_BUILD (rank ',MY_RANK,'): level ',LEV,' box ',IB,' has hidden walls over a level-0 mesh of another rank (fine boxes must lie in one level-0 mesh of their own rank)'
   ERROR STOP 1
ENDIF
ASSOCIATE(FB=>FWALL(LEV)%B(IB))
IF (ASSOCIATED(FB%W)) DEALLOCATE(FB%W)
ALLOCATE(FB%W) ; WM => FB%W
FB%NW = N
IF (ALLOCATED(FB%PNM)) DEALLOCATE(FB%PNM,FB%PIW)
ALLOCATE(FB%PNM(N),FB%PIW(N)) ; FB%PNM = T(5,1:N) ; FB%PIW = T(6,1:N)
IF (ALLOCATED(FB%RNM)) DEALLOCATE(FB%RNM,FB%RIW)
ALLOCATE(FB%RNM(N),FB%RIW(N)) ; FB%RNM = RT(1,1:N) ; FB%RIW = RT(2,1:N)
IF (ALLOCATED(FB%OWN)) DEALLOCATE(FB%OWN,FB%LEAF)
ALLOCATE(FB%OWN(N),FB%LEAF(N)) ; FB%OWN = 0 ; FB%LEAF = .FALSE.
IF (ALLOCATED(FB%PV)) DEALLOCATE(FB%PV,FB%HAVEP)
ALLOCATE(FB%PV(4,N),FB%HAVEP(N)) ; FB%PV = 0._EB ; FB%HAVEP = .FALSE.
IF (ALLOCATED(FB%RKEY)) DEALLOCATE(FB%RKEY,FB%PS,FB%HAVES)
IF (ALLOCATED(FB%PQ)) DEALLOCATE(FB%PQ)
ALLOCATE(FB%PQ(N)) ; FB%PQ = 0._EB
ALLOCATE(FB%RKEY(N),FB%PS(8+3*N_TRACKED_SPECIES,N),FB%HAVES(N)) ; FB%RKEY = 0 ; FB%PS = 0._EB ; FB%HAVES = .FALSE.   ! filled by FDS_FINE_WALLS_PARENT_XCHG every stage
DO IW=1,N ; IF (T(9,IW)>=0) FB%RKEY(IW) = KEYOF(T(9:11,IW),-T(4,IW)) ; ENDDO   ! filled by FDS_FINE_WALLS_RAD_XCHG before each radiation update
ALLOCATE(WM%WALL(0:N),WM%WALL_WORK1(N),WM%WALL_WORK2(N)) ; WM%WALL_WORK1 = 0._EB ; WM%WALL_WORK2 = 0._EB ; WM%N_WALL_CELLS_DIM = N ; WM%N_WALL_CELLS = N ; WM%N_EXTERNAL_WALL_CELLS = 0 ; WM%N_INTERNAL_WALL_CELLS = N
MESH_OBJ => WM
DO IW=1,N
   CALL ALLOCATE_STORAGE(NM,SURF_INDEX=MERGE(INERT_SURF_INDEX,T(7,IW),SURFACE(T(7,IW))%THERMAL_BC_INDEX==THERMALLY_THICK),WALL_INDEX=IW)   ! thick: no ONE_D above level 0 (D-N5); the face takes the level-0 TMP_F
ENDDO
MESH_OBJ => NULL()
DO IW=1,N
   I=T(1,IW) ; J=T(2,IW) ; K=T(3,IW) ; D=T(4,IW)
   WC => WM%WALL(IW)
   WC%SURF_INDEX = T(7,IW) ; WC%BOUNDARY_TYPE = T(8,IW) ; WC%VENT_INDEX = 0
   IF (WC%BR_INDEX>0) THEN   ! wall intensities: the level-0 parent wall's (it has converged its initialization sweeps; ALLOCATE_STORAGE starts at 0, so the wall would reflect nothing and the gas next to it would see a deficit), else ambient equilibrium
      DO NN=1,SIZE(WM%BOUNDARY_RADIA(WC%BR_INDEX)%BAND)
         WM%BOUNDARY_RADIA(WC%BR_INDEX)%BAND(NN)%ILW = SIGMA*TMPA4/PI
         IF (T(5,IW)>0) THEN ; IF (MESHES(T(5,IW))%WALL(T(6,IW))%BR_INDEX>0) WM%BOUNDARY_RADIA(WC%BR_INDEX)%BAND(NN)%ILW = &
            MESHES(T(5,IW))%BOUNDARY_RADIA(MESHES(T(5,IW))%WALL(T(6,IW))%BR_INDEX)%BAND(NN)%ILW ; ENDIF
      ENDDO
   ENDIF
   IF (T(5,IW)>0) WC%VENT_INDEX = MESHES(T(5,IW))%WALL(T(6,IW))%VENT_INDEX
   BC => WM%BOUNDARY_COORD(WC%BC_INDEX)
   BC%IOR = -D ; BC%IIG = I ; BC%JJG = J ; BC%KKG = K ; BC%II = I ; BC%JJ = J ; BC%KK = K
   SELECT CASE(ABS(D))
      CASE(1) ; BC%II = I+SIGN(1,D) ; BC%II2 = I+2*SIGN(1,D) ; BC%JJ2 = J ; BC%KK2 = K
      CASE(2) ; BC%JJ = J+SIGN(1,D) ; BC%JJ2 = J+2*SIGN(1,D) ; BC%II2 = I ; BC%KK2 = K
      CASE(3) ; BC%KK = K+SIGN(1,D) ; BC%KK2 = K+2*SIGN(1,D) ; BC%II2 = I ; BC%JJ2 = J
   END SELECT
   BC%NVEC = 0._EB ; BC%NVEC(ABS(D)) = -REAL(SIGN(1,D),EB)
   BC%X = R%XC(I) ; BC%Y = R%YC(J) ; BC%Z = R%ZC(K)
   SELECT CASE(D)
      CASE(-1) ; BC%X = R%X(I-1)
      CASE( 1) ; BC%X = R%X(I)
      CASE(-2) ; BC%Y = R%Y(J-1)
      CASE( 2) ; BC%Y = R%Y(J)
      CASE(-3) ; BC%Z = R%Z(K-1)
      CASE( 3) ; BC%Z = R%Z(K)
   END SELECT
   IF (OWN_BURNER(WC%SURF_INDEX)) THEN   ! D-010-exempt burner: own area and stock ignition time (init.f90: T_BEGIN + distance from XYZ / SPREAD_RATE) at this face
      CALL FWALL_GEOM(R,BC,FF,WM%BOUNDARY_PROP1(WC%B1_INDEX)%AREA)
      IF (T(12,IW)>0) THEN ; IF (VTX(10,T(12,IW))>0._EB) WM%BOUNDARY_PROP1(WC%B1_INDEX)%T_IGN = T_BEGIN + NORM2((/BC%X,BC%Y,BC%Z/)-VTX(7:9,T(12,IW)))/VTX(10,T(12,IW)) ; ENDIF
      IF (T(12,IW)>0) THEN ; IF (VTX(11,T(12,IW))>0._EB) THEN   ! init.f90 AREA_ADJUST = AREA_MULTIPLIER * INPUT_AREA / FDS_AREA, the vent area snapped at this dx
         FA = 1._EB ; DXY = (/R%DX(1),R%DY(1),R%DZ(1)/)
         DO NN=1,3 ; IF (NN/=ABS(D)) FA = FA*(VTS(2*NN,T(12,IW))-VTS(2*NN-1,T(12,IW)))*DXY(NN) ; ENDDO
         WM%BOUNDARY_PROP1(WC%B1_INDEX)%AREA_ADJUST = SURFACE(WC%SURF_INDEX)%AREA_MULTIPLIER*VTX(11,T(12,IW))/FA
      ENDIF ; ENDIF
   ENDIF
   DO L=1,MIN(NMESHES,SIZE(MESHES))   ! level-0 mesh holding the gas cell
      IF (R%XC(I)>MESHES(L)%XS .AND. R%XC(I)<MESHES(L)%XF .AND. R%YC(J)>MESHES(L)%YS .AND. R%YC(J)<MESHES(L)%YF .AND. R%ZC(K)>MESHES(L)%ZS .AND. R%ZC(K)<MESHES(L)%ZF) THEN
         FB%OWN(IW) = L ; EXIT
      ENDIF
   ENDDO
ENDDO
END ASSOCIATE
IF (WALL_DIAG()) WRITE(0,'(A,I0,A,I0,A,3I6,A,3I5,A,I0,A,6I7,A,I0,A,I0)') 'WALLDIAG level ',LEV,' box ',IB,' g0',G0,' n',R%IBAR,R%JBAR,R%KBAR,' nwall ',N,' perD(-1,1,-2,2,-3,3)', &
   COUNT(T(4,1:N)==-1),COUNT(T(4,1:N)==1),COUNT(T(4,1:N)==-2),COUNT(T(4,1:N)==2),COUNT(T(4,1:N)==-3),COUNT(T(4,1:N)==3),' withparent ',COUNT(T(5,1:N)>0),' thin ',NTHIN
DEALLOCATE(T,RT)
END SUBROUTINE FDS_FINE_WALLS_BUILD

!> True if the point lies inside a box of fine level LEV (extents of the boxes of all ranks; half-open, see FINE_OWNER_BOX).
LOGICAL FUNCTION POINT_COVERED(LEV,X,Y,Z)
INTEGER, INTENT(IN) :: LEV
REAL(EB), INTENT(IN) :: X,Y,Z
POINT_COVERED = FINE_OWNER_BOX(LEV,X,Y,Z)>0
END FUNCTION POINT_COVERED

!> Box of fine level LEV that owns the point (0 if none): half-open intervals [XS,XF) x [YS,YF) x [ZS,ZF), so a point on a face shared by two boxes (a particle on a
!> box edge) has exactly one owner, the box on its high side; boxes of a level do not overlap. Cell centres never lie on a box face, so leaf flags of cells and walls are as before.
INTEGER FUNCTION FINE_OWNER_BOX(LEV,X,Y,Z)
INTEGER, INTENT(IN) :: LEV
REAL(EB), INTENT(IN) :: X,Y,Z
INTEGER :: IB
FINE_OWNER_BOX = 0
IF (.NOT.ALLOCATED(FINE_LEVEL)) RETURN
IF (LEV<1 .OR. LEV>SIZE(FINE_LEVEL)) RETURN
IF (.NOT.ALLOCATED(FINE_LEVEL(LEV)%BOX)) RETURN
DO IB=1,FINE_LEVEL(LEV)%N_BOXES
   ASSOCIATE(B=>FINE_LEVEL(LEV)%BOX(IB))
   IF (X>=B%XS .AND. X<B%XF .AND. Y>=B%YS .AND. Y<B%YF .AND. Z>=B%ZS .AND. Z<B%ZF) THEN ; FINE_OWNER_BOX = IB ; RETURN ; ENDIF
   END ASSOCIATE
ENDDO
END FUNCTION FINE_OWNER_BOX

!> Leaf flag of every hidden wall: its gas cell is not inside a box of the next level (leaf-only sums: each location counted at the finest level covering it).
SUBROUTINE REFRESH_LEAF()
INTEGER :: LEV,IB,IW
IF (LEAF_VER==BUILD_VER .OR. .NOT.ALLOCATED(FWALL)) RETURN
DO LEV=1,SIZE(FWALL)
   IF (.NOT.ALLOCATED(FWALL(LEV)%B)) CYCLE
   DO IB=1,SIZE(FWALL(LEV)%B)
      IF (.NOT.ASSOCIATED(FWALL(LEV)%B(IB)%W)) CYCLE
      ASSOCIATE(R=>FINE_LEVEL(LEV)%BOX(IB),FB=>FWALL(LEV)%B(IB))
      DO IW=1,FB%NW
         ASSOCIATE(BC=>FB%W%BOUNDARY_COORD(FB%W%WALL(IW)%BC_INDEX))
         FB%LEAF(IW) = .NOT.POINT_COVERED(LEV+1,R%XC(BC%IIG),R%YC(BC%JJG),R%ZC(BC%KKG))
         END ASSOCIATE
      ENDDO
      END ASSOCIATE
   ENDDO
ENDDO
! indexes: per box the leaf walls, and all leaf walls grouped by the level-0 mesh that owns them (counting sort, two passes)
BLOCK
INTEGER :: NMX,NM,P,N
INTEGER, ALLOCATABLE :: FP(:)
NMX = MAX(NMESHES,SIZE(MESHES))
IF (ALLOCATED(LW_OFF)) DEALLOCATE(LW_OFF,LW_LEV,LW_BOX,LW_W)
ALLOCATE(LW_OFF(0:NMX)) ; LW_OFF = 0
DO P=1,2
   IF (P==2) THEN
      DO NM=1,NMX ; LW_OFF(NM) = LW_OFF(NM)+LW_OFF(NM-1) ; ENDDO
      ALLOCATE(LW_LEV(MAX(1,LW_OFF(NMX))),LW_BOX(MAX(1,LW_OFF(NMX))),LW_W(MAX(1,LW_OFF(NMX))))
      ALLOCATE(FP(NMX)) ; FP = LW_OFF(0:NMX-1)   ! running fill pointer of each mesh (mesh NM occupies LW_OFF(NM-1)+1..LW_OFF(NM))
   ENDIF
   DO LEV=1,SIZE(FWALL)
      IF (.NOT.ALLOCATED(FWALL(LEV)%B)) CYCLE
      DO IB=1,SIZE(FWALL(LEV)%B)
         IF (.NOT.ASSOCIATED(FWALL(LEV)%B(IB)%W)) CYCLE
         IF (P==1) THEN
            IF (ALLOCATED(FWALL(LEV)%B(IB)%LIDX)) DEALLOCATE(FWALL(LEV)%B(IB)%LIDX)
            N = 0
            DO IW=1,FWALL(LEV)%B(IB)%NW ; IF (FWALL(LEV)%B(IB)%LEAF(IW)) N = N+1 ; ENDDO
            ALLOCATE(FWALL(LEV)%B(IB)%LIDX(N)) ; N = 0
            DO IW=1,FWALL(LEV)%B(IB)%NW ; IF (FWALL(LEV)%B(IB)%LEAF(IW)) THEN ; N = N+1 ; FWALL(LEV)%B(IB)%LIDX(N) = IW ; ENDIF ; ENDDO
         ENDIF
         DO N=1,SIZE(FWALL(LEV)%B(IB)%LIDX)
            IW = FWALL(LEV)%B(IB)%LIDX(N) ; NM = FWALL(LEV)%B(IB)%OWN(IW)
            IF (NM<1 .OR. NM>NMX) CYCLE
            IF (P==1) THEN ; LW_OFF(NM) = LW_OFF(NM)+1
            ELSE ; FP(NM) = FP(NM)+1 ; LW_LEV(FP(NM)) = LEV ; LW_BOX(FP(NM)) = IB ; LW_W(FP(NM)) = IW ; ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO
END BLOCK
! covered level-0 SOLID walls (L0Q, conditions 1 and 2): gas cell under level 1; E is kept (an uncovered wall's next solid update takes it)
BLOCK
USE FDS_AMREX_HOOKS, ONLY: L0Q
USE GLOBAL_CONSTANTS, ONLY: SOLID_BOUNDARY,MY_RANK,PROCESS
INTEGER :: NM,NW
TYPE(MESH_TYPE), POINTER :: M
IF (.NOT.ALLOCATED(L0Q)) ALLOCATE(L0Q(NMESHES))
DO NM=1,NMESHES
   IF (PROCESS(NM)/=MY_RANK .OR. .NOT.ALLOCATED(MESHES(NM)%WALL)) CYCLE
   M => MESHES(NM) ; NW = M%N_EXTERNAL_WALL_CELLS+M%N_INTERNAL_WALL_CELLS
   IF (.NOT.ALLOCATED(L0Q(NM)%COV)) THEN ; ALLOCATE(L0Q(NM)%COV(NW),L0Q(NM)%E(NW)) ; L0Q(NM)%E = 0._EB ; ENDIF
   DO IW=1,NW
      ASSOCIATE(BC=>M%BOUNDARY_COORD(M%WALL(IW)%BC_INDEX))
      L0Q(NM)%COV(IW) = M%WALL(IW)%BOUNDARY_TYPE==SOLID_BOUNDARY .AND. POINT_COVERED(1,M%XC(BC%IIG),M%YC(BC%JJG),M%ZC(BC%KKG))
      END ASSOCIATE
   ENDDO
ENDDO
END BLOCK
LEAF_VER = BUILD_VER
END SUBROUTINE REFRESH_LEAF

!> Start of each corrector (TimeLoop): the predictor fluxes of the leaf hidden walls over covered thermally thick level-0 walls go to their solids (L0Q E, condition 1),
!> each fine face once; L0Q_GAS_CUM counts them (gas side of the closure check).
SUBROUTINE FDS_FINE_L0Q_COMMIT() BIND(C,NAME='fds_fine_l0q_commit')
USE FDS_AMREX_HOOKS, ONLY: L0Q,L0Q_COVERED,L0Q_GAS_CUM
INTEGER :: LEV,IB,IW,NM0,IW0
IF (.NOT.ALLOCATED(FWALL)) RETURN
CALL REFRESH_LEAF()
DO LEV=1,SIZE(FWALL)
   IF (.NOT.ALLOCATED(FWALL(LEV)%B)) CYCLE
   DO IB=1,SIZE(FWALL(LEV)%B)
      IF (.NOT.ASSOCIATED(FWALL(LEV)%B(IB)%W)) CYCLE
      ASSOCIATE(FB=>FWALL(LEV)%B(IB))
      DO IW=1,FB%NW
         IF (FB%PQ(IW)==0._EB) CYCLE
         NM0 = FB%PNM(IW) ; IW0 = FB%PIW(IW)
         IF (NM0==0) THEN ; NM0 = FB%RNM(IW) ; IW0 = FB%RIW(IW) ; ENDIF   ! ruling 2 interim remap
         IF (FB%LEAF(IW) .AND. L0Q_COVERED(NM0,IW0)) THEN
            L0Q(NM0)%E(IW0) = L0Q(NM0)%E(IW0) + FB%PQ(IW) ; L0Q_GAS_CUM = L0Q_GAS_CUM + FB%PQ(IW)
         ENDIF
         FB%PQ(IW) = 0._EB
      ENDDO
      END ASSOCIATE
   ENDDO
ENDDO
END SUBROUTINE FDS_FINE_L0Q_COMMIT

!> Closure check (FDSTL_QCON): S(1) fine energy committed to the covered level-0 solids, S(2) energy the level-0 solids took, S(3) the part of S(2) taken by walls a regrid uncovered (local; the caller sums over ranks)
SUBROUTINE FDS_L0Q_SUMS(S) BIND(C,NAME='fds_l0q_sums')
USE FDS_AMREX_HOOKS, ONLY: L0Q_GAS_CUM,L0Q_SOL_CUM,L0Q_FLUSH_CUM
REAL(C_DOUBLE), INTENT(OUT) :: S(3)
S = (/L0Q_GAS_CUM,L0Q_SOL_CUM,L0Q_FLUSH_CUM/)
END SUBROUTINE FDS_L0Q_SUMS

!> Ruling 2: fine SOLID faces without a level-0 wall: local largest per-box count in one predictor stage, and their energy (J). FDS_FINE_ORPHAN_OUT writes the reduced values to the .out.
!> S(3:4) the remapped ones (interim, WARNING), S(5) the convective energy into the level-0 thick ground solids (L0Q_GND_CUM), S(6) the same of |Q_CON_F|.
SUBROUTINE FDS_FINE_ORPHAN_SUMS(S) BIND(C,NAME='fds_fine_orphan_sums')
USE FDS_AMREX_HOOKS, ONLY: L0Q_GND_CUM,L0Q_GND_ABS
REAL(C_DOUBLE), INTENT(OUT) :: S(6)
S = (/REAL(ORPHAN_NMAX,EB),ORPHAN_E,REAL(REMAP_NMAX,EB),REMAP_E,L0Q_GND_CUM,L0Q_GND_ABS/)
END SUBROUTINE FDS_FINE_ORPHAN_SUMS
SUBROUTINE FDS_FINE_ORPHAN_OUT(S) BIND(C,NAME='fds_fine_orphan_out')
USE GLOBAL_CONSTANTS, ONLY: LU_OUTPUT,MY_RANK
REAL(C_DOUBLE), INTENT(IN) :: S(6)
IF (MY_RANK/=0) RETURN
WRITE(LU_OUTPUT,'(/A,I0,A,ES24.16,A)') ' AMR fine walls without a level-0 wall: ',NINT(S(1)),' (largest per-box count in a predictor stage, summed over ranks), energy ',S(2),' J'
IF (S(3)>0._EB) WRITE(LU_OUTPUT,'(A,I0,A,ES24.16,A,ES24.16,A)') ' WARNING: AMR fine walls remapped to the nearest level-0 wall of the same surface: ',NINT(S(3)), &
   ' (largest per-box count in a predictor stage, summed over ranks), energy ',S(4),' J; level-0 ground convective energy ',S(5),' J (of |flux| ',S(6),' J)'
END SUBROUTINE FDS_FINE_ORPHAN_OUT

!> Ruling 2 interim: the nearest level-0 wall of mesh M0 with surface SI, boundary type BT and IOR -D (WALL_INDEX(D)) in the face plane of cell (I0,J0,K0), within 3 cells tangentially,
!> nearest to the fine cell centre (X,Y,Z); RIW = 0 if none.
SUBROUTINE REMAP_TARGET(M0,I0,J0,K0,D,SI,BT,X,Y,Z,RNM,RIW)
TYPE(MESH_TYPE), POINTER :: M0
INTEGER, INTENT(IN) :: I0,J0,K0,D,SI,BT
REAL(EB), INTENT(IN) :: X,Y,Z
INTEGER, INTENT(OUT) :: RNM,RIW
INTEGER :: A,B,I,J,K,IW,R
REAL(EB) :: DIST,BEST
RNM = 0 ; RIW = 0 ; BEST = HUGE(1._EB)
DO R=1,3
   DO A=-R,R ; DO B=-R,R
      IF (MAX(ABS(A),ABS(B))/=R) CYCLE
      SELECT CASE(ABS(D))
         CASE(1) ; I=I0 ; J=J0+A ; K=K0+B
         CASE(2) ; I=I0+A ; J=J0 ; K=K0+B
         CASE DEFAULT ; I=I0+A ; J=J0+B ; K=K0
      END SELECT
      IF (I<1 .OR. I>M0%IBAR .OR. J<1 .OR. J>M0%JBAR .OR. K<1 .OR. K>M0%KBAR) CYCLE
      IF (M0%CELL(M0%CELL_INDEX(I,J,K))%SOLID) CYCLE
      IW = M0%CELL(M0%CELL_INDEX(I,J,K))%WALL_INDEX(D)
      IF (IW<=0) CYCLE
      IF (M0%WALL(IW)%SURF_INDEX/=SI .OR. M0%WALL(IW)%BOUNDARY_TYPE/=BT) CYCLE
      DIST = (M0%XC(I)-X)**2+(M0%YC(J)-Y)**2+(M0%ZC(K)-Z)**2
      IF (DIST<BEST) THEN ; BEST = DIST ; RIW = IW ; RNM = 1 ; ENDIF
   ENDDO ; ENDDO
   IF (RIW>0) RETURN
ENDDO
END SUBROUTINE REMAP_TARGET

!> Env FDS_AMR_WALL_DIAG set: one WALLDIAG line per box at every wall build (wall counts per direction, for the tests).
LOGICAL FUNCTION WALL_DIAG()
INTEGER :: L
LOGICAL, SAVE :: DONE=.FALSE.,ON=.FALSE.
IF (.NOT.DONE) THEN ; CALL GET_ENVIRONMENT_VARIABLE('FDS_AMR_WALL_DIAG',LENGTH=L) ; ON = L>0 ; DONE = .TRUE. ; ENDIF
WALL_DIAG = ON
END FUNCTION WALL_DIAG

!> Snap the OBST and VENT lists at the dx of box R (FDS: NINT((X-X0)/DX) from the domain corner), and the box offset G0 and domain size NG at this level.
SUBROUTINE SNAP_SETUP(R,LEV)
TYPE(MESH_TYPE), INTENT(IN) :: R
INTEGER, INTENT(IN) :: LEV
REAL(EB) :: X0(3),X1(3),DX(3),DXK(3)
INTEGER :: O,D,KL,LS,NZ,DZ,LO(3),HI(3),RR(3),T1,T2
LOGICAL :: OK
X0 = (/MINVAL(MESHES%XS),MINVAL(MESHES%YS),MINVAL(MESHES%ZS)/) ; X1 = (/MAXVAL(MESHES%XF),MAXVAL(MESHES%YF),MAXVAL(MESHES%ZF)/)
DX = (/R%DXI,R%DETA,R%DZETA/)
G0 = NINT((/R%XS-X0(1),R%YS-X0(2),R%ZS-X0(3)/)/DX) ; NG = NINT((X1-X0)/DX)
IF (ALLOCATED(OBI)) DEALLOCATE(OBI,VTS)
ALLOCATE(OBI(6,NOB),VTS(6,NVT))
DO O=1,NOB ; DO D=1,3 ; OBI(2*D-1:2*D,O) = NINT((OBX(2*D-1:2*D,O)-X0(D))/DX(D)) ; ENDDO ; ENDDO
DO O=1,NVT ; DO D=1,3 ; VTS(2*D-1:2*D,O) = NINT((VTX(2*D-1:2*D,O)-X0(D))/DX(D)) ; ENDDO ; ENDDO
! Thin OBSTs (same rule as ObstMask.H thin_faces): coarsest level KL<=LEV at which the snap has zero thickness in exactly one direction and is thick in the other two. From there on the OBST has
! no solid cells; it is the face set at plane lo*r over the transverse extent snapped at KL, refined by r (nested). Rods and points keep the per-level snap.
IF (ALLOCATED(THD)) DEALLOCATE(THD,THP,THT)
ALLOCATE(THD(NOB),THP(NOB),THT(4,NOB)) ; THD = 0 ; THP = 0 ; THT = 0
DO O=1,NOB
   DO KL=0,LEV
      IF (KL==0) THEN ; DXK = (/MESHES(1)%DXI,MESHES(1)%DETA,MESHES(1)%DZETA/) ; ELSE ; DXK = DXLEV(:,KL) ; ENDIF
      IF (ANY(DXK<=0._EB)) CYCLE
      NZ = 0 ; DZ = 0 ; OK = .TRUE.
      DO D=1,3
         LO(D) = NINT((OBX(2*D-1,O)-X0(D))/DXK(D)) ; HI(D) = NINT((OBX(2*D,O)-X0(D))/DXK(D))
         IF (HI(D)==LO(D)) THEN ; NZ = NZ+1 ; DZ = D ; ELSEIF (HI(D)<LO(D)) THEN ; OK = .FALSE. ; ENDIF
      ENDDO
      IF (OK .AND. NZ==1) THEN
         LS = KL ; EXIT
      ENDIF
   ENDDO
   IF (KL>LEV) CYCLE
   RR = NINT(DXK/DX)   ! DXK is the level-LS cell size at the exit
   IF (LEV>=1) THEN
      T1 = MOD(DZ,3)+1 ; T2 = MOD(DZ+1,3)+1
      THD(O) = DZ ; THP(O) = LO(DZ)*RR(DZ)
      THT(1,O) = LO(T1)*RR(T1) ; THT(2,O) = HI(T1)*RR(T1) ; THT(3,O) = LO(T2)*RR(T2) ; THT(4,O) = HI(T2)*RR(T2)
      OBI(2*DZ,O) = OBI(2*DZ-1,O)   ! no solid cells from here on
      IF (THP(O)<=0 .OR. THP(O)>=NG(DZ)) THD(O) = 0   ! domain boundary plane: the domain wall exists anyway
   ENDIF
ENDDO
END SUBROUTINE SNAP_SETUP

!> Wall of gas cell (I,J,K) of box R on its face D (+-1..3): WALLQ true when the cell across is an OBST cell (snapped at this level) or outside the domain; then SI, BT = its surface
!> and boundary type (last covering VENT, else the face surface of the last OBST containing the cell across, else SURF_DEFAULT).
SUBROUTINE FACE_SURFACE(R,I,J,K,D,SI,BT,WALLQ)
USE GLOBAL_CONSTANTS, ONLY: DEFAULT_SURF_INDEX,SOLID_BOUNDARY
TYPE(MESH_TYPE), INTENT(IN) :: R
INTEGER, INTENT(IN) :: I,J,K,D
INTEGER, INTENT(OUT) :: SI,BT
LOGICAL, INTENT(OUT) :: WALLQ
INTEGER :: A,SG,G(3),GN(3),OWN,O,P,TT,THIN,T1,T2
LOGICAL :: IN
A = ABS(D) ; SG = SIGN(1,D)
G = G0 + (/I-1,J-1,K-1/) ; GN = G ; GN(A) = G(A)+SG
OWN = 0
IF (GN(A)>=0 .AND. GN(A)<NG(A)) THEN
   DO O=NOB,1,-1
      IF (ALL(GN>=OBI(1:5:2,O)) .AND. ALL(GN<=OBI(2:6:2,O)-1)) THEN ; OWN = O ; EXIT ; ENDIF
   ENDDO
   WALLQ = OWN>0
ELSE
   WALLQ = .TRUE.
ENDIF
P = G(A) + MERGE(1,0,SG>0)
THIN = 0   ! zero-thickness OBST face on this face (last OBST wins over the solid cell across)
DO O=NOB,1,-1
   IF (THD(O)/=A) CYCLE
   IF (THP(O)/=P) CYCLE
   T1 = MOD(A,3)+1 ; T2 = MOD(A+1,3)+1
   IF (G(T1)<THT(1,O) .OR. G(T1)>THT(2,O)-1 .OR. G(T2)<THT(3,O) .OR. G(T2)>THT(4,O)-1) CYCLE
   THIN = O ; EXIT
ENDDO
IF (THIN>0) THEN ; IF (THIN>OWN) THEN ; WALLQ = .TRUE. ; OWN = THIN ; NTHIN = NTHIN+1 ; ENDIF ; ENDIF
SI = DEFAULT_SURF_INDEX ; BT = SOLID_BOUNDARY ; FVT = 0
IF (.NOT.WALLQ) RETURN
IF (OWN>0) SI = OBS(-D+4,OWN)
DO O=NVT,1,-1
   IF (VTI(1,O)/=-D) CYCLE
   IF (VTS(2*A-1,O)/=P .OR. VTS(2*A,O)/=P) CYCLE
   IN = .TRUE.
   DO TT=1,3
      IF (TT==A) CYCLE
      IF (G(TT)<VTS(2*TT-1,O) .OR. G(TT)>VTS(2*TT,O)-1) IN = .FALSE.
   ENDDO
   IF (IN) THEN ; SI = VTI(2,O) ; BT = VTI(3,O) ; FVT = O ; EXIT ; ENDIF
ENDDO
END SUBROUTINE FACE_SURFACE

!> Radiation state of the hidden walls of level LEV (>= 1), pulled from the walls of level LEV-1 only (one-way, read-only toward the fine box; collective over all ranks, called by the driver
!> before the radiation update of level LEV, so after every regrid and every step). Every rank packs the solid/open walls of its level LEV-1 data (level 0: the WALL array of its meshes; level >= 1:
!> the hidden walls of its boxes) as records (global cell index at that level, IOR, TMP_F, EMISSIVITY, Q_RAD_OUT, T_IGN, SURF_INDEX, BOUNDARY_TYPE), all ranks gather them, and each local box
!> wall finds its parent record by (parent cell, IOR) when its face lies on the parent cell face and the surface and boundary type agree; otherwise it has no parent (HAVEP false).
SUBROUTINE FDS_FINE_WALLS_RAD_XCHG(LEV) BIND(C,NAME='fds_fine_walls_rad_xchg')
USE MPI_F08
USE GLOBAL_CONSTANTS, ONLY: MY_RANK,PROCESS,N_MPI_PROCESSES
INTEGER(C_INT), VALUE :: LEV
INTEGER, PARAMETER :: NF=10
INTEGER :: LP,NM0,IB,IW,N,NT,IERR,I,J,K,A,RATIO,P,G(3),GP(3),D,NTOT
INTEGER, ALLOCATABLE :: CNT(:),DSP(:),IX(:)
INTEGER(8), ALLOCATABLE :: KEY(:)
INTEGER(8) :: KQ
REAL(EB), ALLOCATABLE :: SB(:,:),RB(:,:)
REAL(EB) :: X0(3),DXP,DXF
TYPE(MESH_TYPE), POINTER :: R,WM
IF (.NOT.ALLOCATED(FINE_LEVEL)) RETURN   ! the same on every rank (the driver creates every level on every rank): no collective is skipped
IF (.NOT.ALLOCATED(FWALL)) ALLOCATE(FWALL(0))
LP = LEV-1
X0 = (/MINVAL(MESHES(1:NMESHES)%XS),MINVAL(MESHES(1:NMESHES)%YS),MINVAL(MESHES(1:NMESHES)%ZS)/)
! pack the walls of level LP owned by this rank
N = 0
ALLOCATE(SB(NF,0:0))
DO P=1,2   ! 1: count, 2: fill
   IF (P==2) THEN ; DEALLOCATE(SB) ; ALLOCATE(SB(NF,MAX(1,N))) ; ENDIF
   N = 0
   IF (LP==0) THEN
      DO NM0=1,NMESHES
         IF (PROCESS(NM0)/=MY_RANK .OR. .NOT.ALLOCATED(MESHES(NM0)%WALL)) CYCLE
         DO IW=1,MESHES(NM0)%N_EXTERNAL_WALL_CELLS+MESHES(NM0)%N_INTERNAL_WALL_CELLS
            ASSOCIATE(WC=>MESHES(NM0)%WALL(IW))
            IF (WC%B1_INDEX<1) CYCLE
            N = N+1
            IF (P==2) CALL PACK(MESHES(NM0),WC,MESHES(NM0)%BOUNDARY_COORD(WC%BC_INDEX),MESHES(NM0)%BOUNDARY_PROP1(WC%B1_INDEX),SB(:,N))
            END ASSOCIATE
         ENDDO
      ENDDO
   ELSEIF (SIZE(FWALL)>=LP) THEN
     IF (ALLOCATED(FWALL(LP)%B)) THEN
      DO IB=1,SIZE(FWALL(LP)%B)
         WM => FWALL(LP)%B(IB)%W
         IF (.NOT.ASSOCIATED(WM)) CYCLE
         DO IW=1,FWALL(LP)%B(IB)%NW
            N = N+1
            IF (P==2) CALL PACK(FINE_LEVEL(LP)%BOX(IB),WM%WALL(IW),WM%BOUNDARY_COORD(WM%WALL(IW)%BC_INDEX),WM%BOUNDARY_PROP1(WM%WALL(IW)%B1_INDEX),SB(:,N))
         ENDDO
      ENDDO
     ENDIF
   ENDIF
ENDDO
! gather the records of all ranks
ALLOCATE(CNT(0:N_MPI_PROCESSES-1),DSP(0:N_MPI_PROCESSES-1))
CNT(MY_RANK) = N*NF
IF (N_MPI_PROCESSES>1) CALL MPI_ALLGATHER(MPI_IN_PLACE,1,MPI_INTEGER,CNT,1,MPI_INTEGER,MPI_COMM_WORLD,IERR)
DSP(0) = 0 ; DO I=1,N_MPI_PROCESSES-1 ; DSP(I) = DSP(I-1)+CNT(I-1) ; ENDDO
NTOT = SUM(CNT)/NF
ALLOCATE(RB(NF,MAX(1,NTOT)))
IF (N_MPI_PROCESSES>1) THEN
   CALL MPI_ALLGATHERV(SB,N*NF,MPI_DOUBLE_PRECISION,RB,CNT,DSP,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
ELSE
   RB(:,1:N) = SB(:,1:N)
ENDIF
! sort the records by key
ALLOCATE(KEY(NTOT),IX(NTOT))
DO I=1,NTOT ; KEY(I) = KEYOF(NINT(RB(1:3,I)),NINT(RB(4,I))) ; IX(I) = I ; ENDDO
CALL SORT_KEYS(NTOT,KEY,IX)
! parent lookup of every wall of the local boxes of level LEV
DXP = 0._EB   ! cell size of level LP: from a box of this rank (a box of another rank has no cells here), the largest over the ranks
IF (LP==0) THEN
   DXP = (MESHES(1)%XF-MESHES(1)%XS)/MESHES(1)%IBAR
ELSE
   DO IB=1,FINE_LEVEL(LP)%N_BOXES
      IF (FINE_LEVEL(LP)%BOX(IB)%IBAR>0) DXP = MAX(DXP,(FINE_LEVEL(LP)%BOX(IB)%XF-FINE_LEVEL(LP)%BOX(IB)%XS)/FINE_LEVEL(LP)%BOX(IB)%IBAR)
   ENDDO
   IF (N_MPI_PROCESSES>1) CALL MPI_ALLREDUCE(MPI_IN_PLACE,DXP,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERR)
ENDIF
IF (SIZE(FWALL)<LEV) RETURN
IF (.NOT.ALLOCATED(FWALL(LEV)%B)) RETURN
DO IB=1,SIZE(FWALL(LEV)%B)
   WM => FWALL(LEV)%B(IB)%W
   IF (.NOT.ASSOCIATED(WM)) CYCLE
   R => FINE_LEVEL(LEV)%BOX(IB)
   DXF = (R%XF-R%XS)/R%IBAR ; RATIO = NINT(DXP/DXF)
   FWALL(LEV)%B(IB)%HAVEP = .FALSE.
   DO IW=1,FWALL(LEV)%B(IB)%NW
      ASSOCIATE(WC=>WM%WALL(IW),BC=>WM%BOUNDARY_COORD(WM%WALL(IW)%BC_INDEX))
      G = (/NINT((R%XC(BC%IIG)-X0(1))/DXF-0.5_EB),NINT((R%YC(BC%JJG)-X0(2))/DXF-0.5_EB),NINT((R%ZC(BC%KKG)-X0(3))/DXF-0.5_EB)/)
      A = ABS(BC%IOR) ; D = -BC%IOR
      IF (D<0 .AND. MOD(G(A),RATIO)/=0) CYCLE   ! the fine face is inside the parent cell, not on its face
      IF (D>0 .AND. MOD(G(A),RATIO)/=RATIO-1) CYCLE
      GP = G/RATIO
      KQ = KEYOF(GP,BC%IOR)
      I = FIND_KEY(NTOT,KEY,KQ)
      IF (I==0) CYCLE
      K = IX(I)
      IF (NINT(RB(9,K))/=WC%SURF_INDEX .OR. NINT(RB(10,K))/=WC%BOUNDARY_TYPE) CYCLE
      FWALL(LEV)%B(IB)%PV(:,IW) = RB(5:8,K) ; FWALL(LEV)%B(IB)%HAVEP(IW) = .TRUE.
      END ASSOCIATE
   ENDDO
ENDDO

CONTAINS

SUBROUTINE PACK(M,WC,BC,B1,REC)
USE TYPES, ONLY: WALL_TYPE,BOUNDARY_COORD_TYPE,BOUNDARY_PROP1_TYPE
TYPE(MESH_TYPE), INTENT(IN) :: M
TYPE(WALL_TYPE), INTENT(IN) :: WC
TYPE(BOUNDARY_COORD_TYPE), INTENT(IN) :: BC
TYPE(BOUNDARY_PROP1_TYPE), INTENT(IN) :: B1
REAL(EB), INTENT(OUT) :: REC(NF)
REAL(EB) :: DXM
DXM = (M%XF-M%XS)/M%IBAR
REC(1) = NINT((M%XC(BC%IIG)-X0(1))/DXM-0.5_EB) ; REC(2) = NINT((M%YC(BC%JJG)-X0(2))/DXM-0.5_EB) ; REC(3) = NINT((M%ZC(BC%KKG)-X0(3))/DXM-0.5_EB)
REC(4) = BC%IOR ; REC(5) = B1%TMP_F ; REC(6) = B1%EMISSIVITY ; REC(7) = B1%Q_RAD_OUT ; REC(8) = B1%T_IGN ; REC(9) = WC%SURF_INDEX ; REC(10) = WC%BOUNDARY_TYPE
END SUBROUTINE PACK

END SUBROUTINE FDS_FINE_WALLS_RAD_XCHG

INTEGER(8) FUNCTION KEYOF(GG,IOR)
INTEGER, INTENT(IN) :: GG(3),IOR
KEYOF = ((INT(GG(3),8)*100000_8+INT(GG(2),8))*100000_8+INT(GG(1),8))*8_8+INT(IOR+3,8)
END FUNCTION KEYOF

SUBROUTINE SORT_KEYS(NN,KY,ID)   ! shell sort of KY with the index list ID
INTEGER, INTENT(IN) :: NN
INTEGER(8), INTENT(INOUT) :: KY(NN)
INTEGER, INTENT(INOUT) :: ID(NN)
INTEGER :: GAP,II,JJ,IT
INTEGER(8) :: KT
GAP = NN/2
DO WHILE (GAP>0)
   DO II=GAP+1,NN
      KT = KY(II) ; IT = ID(II) ; JJ = II
      DO WHILE (JJ>GAP)
         IF (KY(JJ-GAP)<=KT) EXIT
         KY(JJ) = KY(JJ-GAP) ; ID(JJ) = ID(JJ-GAP) ; JJ = JJ-GAP
      ENDDO
      KY(JJ) = KT ; ID(JJ) = IT
   ENDDO
   GAP = GAP/2
ENDDO
END SUBROUTINE SORT_KEYS

INTEGER FUNCTION FIND_KEY(NN,KY,KK)   ! index of KK in the sorted KY, 0 if absent
INTEGER, INTENT(IN) :: NN
INTEGER(8), INTENT(IN) :: KY(NN),KK
INTEGER :: LO,HI,MD
FIND_KEY = 0 ; LO = 1 ; HI = NN
DO WHILE (LO<=HI)
   MD = (LO+HI)/2
   IF (KY(MD)==KK) THEN ; FIND_KEY = MD ; RETURN ; ENDIF
   IF (KY(MD)<KK) THEN ; LO = MD+1 ; ELSE ; HI = MD-1 ; ENDIF
ENDDO
END FUNCTION FIND_KEY


!> Global level-0 cell (0-based, from the domain corner) of the parent gas cell of hidden fine wall (I,J,K,D) of box R when its level-0 mesh NM0 is on another rank
!> (uniform level-0 cells; mesh extents are known on every rank). GK = -1 when the fine face does not lie on a face of that parent cell.
SUBROUTINE REMOTE_PARENT_CELL(R,I,J,K,D,NM0,GK)
TYPE(MESH_TYPE), INTENT(IN) :: R
INTEGER, INTENT(IN) :: I,J,K,D,NM0
INTEGER, INTENT(OUT) :: GK(3)
REAL(EB) :: D0(3),XS(3),GX0(3),FF,TOL,F0
INTEGER :: IJ0(3),A
GK = -1
XS = (/MESHES(NM0)%XS,MESHES(NM0)%YS,MESHES(NM0)%ZS/)
D0 = (/(MESHES(NM0)%XF-XS(1))/MESHES(NM0)%IBAR,(MESHES(NM0)%YF-XS(2))/MESHES(NM0)%JBAR,(MESHES(NM0)%ZF-XS(3))/MESHES(NM0)%KBAR/)
GX0 = (/MINVAL(MESHES(1:NMESHES)%XS),MINVAL(MESHES(1:NMESHES)%YS),MINVAL(MESHES(1:NMESHES)%ZS)/)
IJ0 = FLOOR(((/R%XC(I),R%YC(J),R%ZC(K)/)-XS)/D0)+1
A = ABS(D)
SELECT CASE(A)
   CASE(1) ; FF = MERGE(R%X(I),R%X(I-1),D>0) ; TOL = R%DXI
   CASE(2) ; FF = MERGE(R%Y(J),R%Y(J-1),D>0) ; TOL = R%DETA
   CASE(3) ; FF = MERGE(R%Z(K),R%Z(K-1),D>0) ; TOL = R%DZETA
END SELECT
F0 = XS(A)+REAL(MERGE(IJ0(A),IJ0(A)-1,D>0),EB)*D0(A)
IF (ABS(FF-F0)>1.E-3_EB*TOL) RETURN
GK = NINT((XS+(REAL(IJ0,EB)-0.5_EB)*D0-GX0)/D0-0.5_EB)
END SUBROUTINE REMOTE_PARENT_CELL

!> Level-0 wall state of the hidden fine walls without a local parent (collective; the driver calls it on every rank after the level-0 WALL_BC of each stage and after a bind).
!> Every rank packs all its level-0 SOLID walls under level 1, keyed by (global gas cell, IOR), with the face centre and the PS layout; all ranks gather them. A local fine wall with
!> a parent on another rank (RKEY) copies its parent's record into PS; a fine wall without parent copies the nearest record of the same surface and IOR (ties: lowest key), the
!> same on every rank count. A burner face (PYROLYSIS_SPECIFIED) without a parent state is an error: its flux would be missing from the gas and the leaf sums.
SUBROUTINE FDS_FINE_WALLS_PARENT_XCHG() BIND(C,NAME='fds_fine_walls_parent_xchg')
USE MPI_F08
USE GLOBAL_CONSTANTS, ONLY: MY_RANK,PROCESS,N_MPI_PROCESSES,SOLID_BOUNDARY,PYROLYSIS_SPECIFIED
USE TYPES, ONLY: SURFACE
INTEGER :: NZ,NPS,NF,NM0,IW,N,P,I,K,IERR,NTOT,LEV,IB
INTEGER, ALLOCATABLE :: CNT(:),DSP(:),IX(:)
INTEGER(8), ALLOCATABLE :: KEY(:)
REAL(EB), ALLOCATABLE :: SB(:,:),RB(:,:)
REAL(EB) :: D0(3),GX0(3),DB,DD
TYPE(MESH_TYPE), POINTER :: M
IF (.NOT.ALLOCATED(FINE_LEVEL)) RETURN   ! the same on every rank
NZ = N_TRACKED_SPECIES ; NPS = 8+3*NZ ; NF = 7+NPS   ! record: global gas cell (3), IOR, face centre (3), PS layout
GX0 = (/MINVAL(MESHES(1:NMESHES)%XS),MINVAL(MESHES(1:NMESHES)%YS),MINVAL(MESHES(1:NMESHES)%ZS)/)
N = 0
DO P=1,2   ! 1: count, 2: fill
   IF (P==2) ALLOCATE(SB(NF,MAX(1,N)))
   N = 0
   DO NM0=1,NMESHES
      IF (PROCESS(NM0)/=MY_RANK .OR. .NOT.ALLOCATED(MESHES(NM0)%WALL)) CYCLE
      M => MESHES(NM0)
      D0 = (/(M%XF-M%XS)/M%IBAR,(M%YF-M%YS)/M%JBAR,(M%ZF-M%ZS)/M%KBAR/)
      DO IW=1,M%N_EXTERNAL_WALL_CELLS+M%N_INTERNAL_WALL_CELLS
         IF (M%WALL(IW)%BOUNDARY_TYPE/=SOLID_BOUNDARY .OR. M%WALL(IW)%B1_INDEX<1) CYCLE
         ASSOCIATE(BC=>M%BOUNDARY_COORD(M%WALL(IW)%BC_INDEX),B1=>M%BOUNDARY_PROP1(M%WALL(IW)%B1_INDEX))
         IF (.NOT.POINT_COVERED(1,M%XC(BC%IIG),M%YC(BC%JJG),M%ZC(BC%KKG))) CYCLE
         N = N+1
         IF (P==2) THEN
            SB(1:3,N) = REAL(NINT(((/M%XC(BC%IIG),M%YC(BC%JJG),M%ZC(BC%KKG)/)-GX0)/D0-0.5_EB),EB) ; SB(4,N) = BC%IOR
            SB(5:7,N) = (/BC%X,BC%Y,BC%Z/)
            SB(8:7+NPS,N) = (/B1%U_NORMAL_S,B1%U_NORMAL,B1%AREA,B1%RHO_F,B1%TMP_F,REAL(M%WALL(IW)%SURF_INDEX,EB),B1%ZZ_F(1:NZ),B1%M_DOT_G_PP_ADJUST(1:NZ), &
                              B1%T_IGN,B1%Q_RAD_IN,B1%RHO_D_F(1:NZ)/)
         ENDIF
         END ASSOCIATE
      ENDDO
   ENDDO
ENDDO
ALLOCATE(CNT(0:N_MPI_PROCESSES-1),DSP(0:N_MPI_PROCESSES-1))
CNT = 0 ; CNT(MY_RANK) = N*NF
CALL MPI_ALLGATHER(MPI_IN_PLACE,1,MPI_INTEGER,CNT,1,MPI_INTEGER,MPI_COMM_WORLD,IERR)
DSP(0) = 0 ; DO I=1,N_MPI_PROCESSES-1 ; DSP(I) = DSP(I-1)+CNT(I-1) ; ENDDO
NTOT = SUM(CNT)/NF
ALLOCATE(RB(NF,MAX(1,NTOT)))
CALL MPI_ALLGATHERV(SB,N*NF,MPI_DOUBLE_PRECISION,RB,CNT,DSP,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IERR)
ALLOCATE(KEY(NTOT),IX(NTOT))
DO I=1,NTOT ; KEY(I) = KEYOF(NINT(RB(1:3,I)),NINT(RB(4,I))) ; IX(I) = I ; ENDDO
CALL SORT_KEYS(NTOT,KEY,IX)
IF (.NOT.ALLOCATED(FWALL)) RETURN
DO LEV=1,SIZE(FWALL)
   IF (.NOT.ALLOCATED(FWALL(LEV)%B)) CYCLE
   DO IB=1,SIZE(FWALL(LEV)%B)
      IF (.NOT.ASSOCIATED(FWALL(LEV)%B(IB)%W)) CYCLE
      ASSOCIATE(FB=>FWALL(LEV)%B(IB))
      FB%HAVES = .FALSE.
      DO IW=1,FB%NW
         IF (FB%W%WALL(IW)%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE
         K = 0
         IF (FB%RKEY(IW)/=0) THEN   ! parent on another rank: unreachable (FDS_FINE_WALLS_BUILD stops on fine cells over another rank's mesh)
            WRITE(0,'(A,I0,A,I0,A,I0,A)') 'fds_amr ERROR in FDS_FINE_WALLS_PARENT_XCHG (rank ',MY_RANK,'): level ',LEV,' box ',IB, &
               ' has a hidden wall whose level-0 parent wall is not on this rank (fine boxes must lie in one level-0 mesh of their own rank)'
            ERROR STOP 1
         ELSEIF (FB%PNM(IW)==0) THEN   ! no parent: nearest record of the same surface and IOR, in key order (rank independent)
            DB = HUGE(1._EB)
            ASSOCIATE(BC=>FB%W%BOUNDARY_COORD(FB%W%WALL(IW)%BC_INDEX))
            DO I=1,NTOT
               IF (NINT(RB(13,IX(I)))/=FB%W%WALL(IW)%SURF_INDEX .OR. NINT(RB(4,IX(I)))/=BC%IOR) CYCLE
               DD = (RB(5,IX(I))-BC%X)**2+(RB(6,IX(I))-BC%Y)**2+(RB(7,IX(I))-BC%Z)**2
               IF (DD<DB) THEN ; DB = DD ; K = IX(I) ; ENDIF
            ENDDO
            END ASSOCIATE
         ENDIF
         IF (K>0) THEN ; IF (NINT(RB(13,K))/=FB%W%WALL(IW)%SURF_INDEX) K = 0 ; ENDIF
         IF (K>0) THEN ; FB%PS(:,IW) = RB(8:7+NPS,K) ; FB%HAVES(IW) = .TRUE. ; ENDIF
         IF (SURFACE(FB%W%WALL(IW)%SURF_INDEX)%PYROLYSIS_MODEL==PYROLYSIS_SPECIFIED .AND. FB%PNM(IW)==0 .AND. .NOT.(FB%RKEY(IW)/=0 .AND. FB%HAVES(IW))) THEN
            WRITE(0,'(A,I0,A,I0,A,I0,3A)') 'fds_amr ERROR in FDS_FINE_WALLS_PARENT_XCHG (rank ',MY_RANK,'): level ',LEV,' box ',IB,' has a face of burner surface ', &
               TRIM(SURFACE(FB%W%WALL(IW)%SURF_INDEX)%ID),' without a level-0 parent wall state (no level-0 face of that surface under it); refine the burner whole or keep it on level 0 (--nofine)'
            ERROR STOP 1
         ENDIF
      ENDDO
      END ASSOCIATE
   ENDDO
ENDDO
END SUBROUTINE FDS_FINE_WALLS_PARENT_XCHG

!> Level-0 parent state of hidden fine wall IW of box FB: V = U_NORMAL_S, U_NORMAL, AREA, RHO_F, TMP_F, SURF_INDEX, ZZ_F(1:NZ), M_DOT_G_PP_ADJUST(1:NZ) of the parent SOLID
!> wall, from the local level-0 mesh or (parent on another rank) from PS. False if the wall has none.
LOGICAL FUNCTION PARENT_STATE(FB,IW,V)
USE GLOBAL_CONSTANTS, ONLY: SOLID_BOUNDARY
TYPE(FWALL_BOX_TYPE), INTENT(IN) :: FB
INTEGER, INTENT(IN) :: IW
REAL(EB), INTENT(OUT) :: V(6+2*N_TRACKED_SPECIES)
INTEGER :: NM0,IW0,NZ
NZ = N_TRACKED_SPECIES
PARENT_STATE = .FALSE. ; V = 0._EB
IF (OWN_BURNER(FB%W%WALL(IW)%SURF_INDEX)) THEN   ! D-010-exempt burner (PV-13): the face's own stock state, not the parent's
   ASSOCIATE(B1=>FB%W%BOUNDARY_PROP1(FB%W%WALL(IW)%B1_INDEX))
   V = (/B1%U_NORMAL_S,B1%U_NORMAL,B1%AREA,B1%RHO_F,B1%TMP_F,REAL(FB%W%WALL(IW)%SURF_INDEX,EB),B1%ZZ_F(1:NZ),B1%M_DOT_G_PP_ADJUST(1:NZ)/)
   END ASSOCIATE
   PARENT_STATE = .TRUE. ; RETURN
ENDIF
NM0 = FB%PNM(IW) ; IW0 = FB%PIW(IW)
IF (NM0>0) THEN
   IF (MESHES(NM0)%WALL(IW0)%BOUNDARY_TYPE/=SOLID_BOUNDARY) RETURN
   ASSOCIATE(B1=>MESHES(NM0)%BOUNDARY_PROP1(MESHES(NM0)%WALL(IW0)%B1_INDEX))
   V = (/B1%U_NORMAL_S,B1%U_NORMAL,B1%AREA,B1%RHO_F,B1%TMP_F,REAL(MESHES(NM0)%WALL(IW0)%SURF_INDEX,EB),B1%ZZ_F(1:NZ),B1%M_DOT_G_PP_ADJUST(1:NZ)/)
   END ASSOCIATE
   PARENT_STATE = .TRUE.
ELSEIF (FB%RKEY(IW)/=0) THEN   ! parent on another rank (a wall without parent has no parent state: its PS is for the sums only)
   IF (.NOT.FB%HAVES(IW)) RETURN
   V = FB%PS(1:6+2*N_TRACKED_SPECIES,IW) ; PARENT_STATE = .TRUE.
ENDIF
END FUNCTION PARENT_STATE

!> Stock wall BC (wall.f90 WALL_CELL_BC) on the exempt-burner hidden walls of fine box NM (NM > NMESHES), after the box's WALL_BC: their own ignition time, flux and surface state.
SUBROUTINE FDS_FINE_BURNER_BC(T,DT,NM) BIND(C,NAME='fds_fine_burner_bc')
USE MESH_VARIABLES, ONLY: MESH_OBJ
USE MESH_POINTERS, ONLY: WALL,BOUNDARY_COORD,BOUNDARY_PROP1,BOUNDARY_PROP2,BOUNDARY_RADIA,POINT_TO_BOX
USE WALL_ROUTINES, ONLY: WALL_CELL_BC
REAL(C_DOUBLE), VALUE :: T,DT
INTEGER(C_INT), VALUE :: NM
INTEGER :: LEV,IB,IW
TYPE(MESH_TYPE), POINTER :: WM
IF (.NOT.ALLOCATED(FWALL) .OR. .NOT.ALLOCATED(FINE_LEVEL)) RETURN
DO LEV=1,MIN(SIZE(FWALL),SIZE(FINE_LEVEL))
   IB = NM-FINE_LEVEL(LEV)%NM0
   IF (IB<1 .OR. IB>FINE_LEVEL(LEV)%N_BOXES .OR. .NOT.ALLOCATED(FWALL(LEV)%B)) CYCLE
   WM => FWALL(LEV)%B(IB)%W
   IF (.NOT.ASSOCIATED(WM)) RETURN
   IF (.NOT.ANY((/(OWN_BURNER(WM%WALL(IW)%SURF_INDEX),IW=1,FWALL(LEV)%B(IB)%NW)/))) RETURN
   CALL POINT_TO_BOX(NM) ; MESH_OBJ => WM
   WALL=>WM%WALL ; BOUNDARY_COORD=>WM%BOUNDARY_COORD ; BOUNDARY_PROP1=>WM%BOUNDARY_PROP1 ; BOUNDARY_PROP2=>WM%BOUNDARY_PROP2 ; BOUNDARY_RADIA=>WM%BOUNDARY_RADIA
   DO IW=1,FWALL(LEV)%B(IB)%NW
      IF (OWN_BURNER(WM%WALL(IW)%SURF_INDEX)) CALL WALL_CELL_BC(T,DT,NM,IW)
   ENDDO
   MESH_OBJ => NULL() ; CALL POINT_TO_BOX(NM) ; RETURN
ENDDO
END SUBROUTINE FDS_FINE_BURNER_BC

!> D-010-exempt burner surface (fds_mesh_query FDS_STATEFUL_VENTS): HRRPUA/MLRPUA with a specified temperature. Its fine faces carry their own stock state.
LOGICAL FUNCTION OWN_BURNER(SI)
USE TYPES, ONLY: SURFACE
USE GLOBAL_CONSTANTS, ONLY: SPECIFIED_TEMPERATURE
INTEGER, INTENT(IN) :: SI
OWN_BURNER = BED_RUN .AND. (SURFACE(SI)%HRRPUA>0._EB .OR. SURFACE(SI)%MLRPUA>0._EB) .AND. SURFACE(SI)%THERMAL_BC_INDEX==SPECIFIED_TEMPERATURE &
             .AND. .NOT.SURFACE(SI)%VEG_LSET_SPREAD
END FUNCTION OWN_BURNER

!> The run has a particle bed (W2a): the D-010 exemption and the own state of exempt burner faces apply (runs without a bed are unchanged).
SUBROUTINE FDS_FINE_SET_BED() BIND(C,NAME='fds_fine_set_bed')
BED_RUN = .TRUE.
END SUBROUTINE FDS_FINE_SET_BED

!> ROUGHNESS of the level-0 ground wall of every level-0 column of the domain (collective, once: called by FDS_FINE_LEVEL_CREATE on every rank). Each rank fills the
!> columns of its ground meshes, a MAX reduction makes the map global (roughness >= 0, 0 where there is no ground wall). Uniform level-0 cells assumed.
SUBROUTINE GROUND_ROUGH_INIT()
USE MPI_F08
USE GLOBAL_CONSTANTS, ONLY: MY_RANK,PROCESS,N_MPI_PROCESSES
USE TYPES, ONLY: SURFACE
INTEGER :: NM,I,J,IW,NX,NY,IERR,I0,J0
REAL(EB) :: D0(2),GX0(3)
IF (ALLOCATED(GROUGH)) RETURN
GX0 = (/MINVAL(MESHES(1:NMESHES)%XS),MINVAL(MESHES(1:NMESHES)%YS),MINVAL(MESHES(1:NMESHES)%ZS)/)
D0 = (/(MESHES(1)%XF-MESHES(1)%XS)/MESHES(1)%IBAR,(MESHES(1)%YF-MESHES(1)%YS)/MESHES(1)%JBAR/)
NX = NINT((MAXVAL(MESHES(1:NMESHES)%XF)-GX0(1))/D0(1)) ; NY = NINT((MAXVAL(MESHES(1:NMESHES)%YF)-GX0(2))/D0(2))
ALLOCATE(GROUGH(NX,NY)) ; GROUGH = 0._EB
DO NM=1,NMESHES
   IF (PROCESS(NM)/=MY_RANK .OR. MESHES(NM)%ZS>GX0(3) .OR. .NOT.ALLOCATED(MESHES(NM)%CELL_INDEX)) CYCLE
   I0 = NINT((MESHES(NM)%XS-GX0(1))/D0(1)) ; J0 = NINT((MESHES(NM)%YS-GX0(2))/D0(2))
   DO J=1,MESHES(NM)%JBAR
      DO I=1,MESHES(NM)%IBAR
         IW = MESHES(NM)%CELL(MESHES(NM)%CELL_INDEX(I,J,1))%WALL_INDEX(-3)
         IF (IW>0) GROUGH(I0+I,J0+J) = SURFACE(MESHES(NM)%BOUNDARY_PROP1(MESHES(NM)%WALL(IW)%B1_INDEX)%SURF_INDEX)%ROUGHNESS
      ENDDO
   ENDDO
ENDDO
IF (N_MPI_PROCESSES>1) CALL MPI_ALLREDUCE(MPI_IN_PLACE,GROUGH,NX*NY,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERR)
END SUBROUTINE GROUND_ROUGH_INIT

!> Radiation on fine box NM (radi.f90 COMPUTE_RADIATION, same code for every level): ON=1 points the wall variables of MESH_POINTERS to the hidden walls of the box and copies the surface
!> state radiation reads (TMP_F, EMISSIVITY, Q_RAD_OUT, T_IGN of the level-0 parent wall; AREA of the fine face), so the face ghost intensity is the wall's, as at level 0; ON=0 points back to the box.
SUBROUTINE FRAD_WALLS(NM,ON)
USE MESH_POINTERS, ONLY: WALL,BOUNDARY_COORD,BOUNDARY_PROP1,BOUNDARY_RADIA,WALL_WORK1,WALL_WORK2,N_WALL_CELLS,N_INTERNAL_WALL_CELLS,N_EXTERNAL_WALL_CELLS,POINT_TO_BOX
USE TYPES, ONLY: SURFACE
INTEGER, INTENT(IN) :: NM,ON
INTEGER :: LEV,IB,IW
REAL(EB) :: RDN,AREA
TYPE(MESH_TYPE), POINTER :: WM
IF (ON==0) THEN ; CALL POINT_TO_BOX(NM) ; RETURN ; ENDIF
IF (.NOT.ALLOCATED(FWALL) .OR. .NOT.ALLOCATED(FINE_LEVEL)) RETURN
DO LEV=1,MIN(SIZE(FWALL),SIZE(FINE_LEVEL))
   IB = NM-FINE_LEVEL(LEV)%NM0
   IF (IB<1 .OR. IB>FINE_LEVEL(LEV)%N_BOXES) CYCLE
   IF (.NOT.ALLOCATED(FWALL(LEV)%B)) RETURN
   WM => FWALL(LEV)%B(IB)%W
   IF (.NOT.ASSOCIATED(WM)) RETURN
   DO IW=1,FWALL(LEV)%B(IB)%NW
      ASSOCIATE(B1=>WM%BOUNDARY_PROP1(WM%WALL(IW)%B1_INDEX),BC=>WM%BOUNDARY_COORD(WM%WALL(IW)%BC_INDEX),SF=>SURFACE(WM%WALL(IW)%SURF_INDEX))
      CALL FWALL_GEOM(FINE_LEVEL(LEV)%BOX(IB),BC,RDN,AREA)
      B1%AREA = AREA
      IF (OWN_BURNER(WM%WALL(IW)%SURF_INDEX)) THEN   ! exempt burner: its own stock state
      ELSEIF (FWALL(LEV)%B(IB)%HAVEP(IW)) THEN   ! the wall of level LEV-1 on the same face (one-way pull, FDS_FINE_WALLS_RAD_XCHG)
         B1%TMP_F = FWALL(LEV)%B(IB)%PV(1,IW) ; B1%EMISSIVITY = FWALL(LEV)%B(IB)%PV(2,IW) ; B1%Q_RAD_OUT = FWALL(LEV)%B(IB)%PV(3,IW) ; B1%T_IGN = FWALL(LEV)%B(IB)%PV(4,IW)
      ELSE   ! no parent wall (new fine face, OBST unresolved above): surface values
         B1%TMP_F = SF%TMP_FRONT ; B1%EMISSIVITY = SF%EMISSIVITY ; B1%Q_RAD_OUT = 0._EB
      ENDIF
      END ASSOCIATE
   ENDDO
   WALL=>WM%WALL ; BOUNDARY_COORD=>WM%BOUNDARY_COORD ; BOUNDARY_PROP1=>WM%BOUNDARY_PROP1 ; BOUNDARY_RADIA=>WM%BOUNDARY_RADIA
   WALL_WORK1=>WM%WALL_WORK1 ; WALL_WORK2=>WM%WALL_WORK2
   N_WALL_CELLS=>WM%N_WALL_CELLS ; N_INTERNAL_WALL_CELLS=>WM%N_INTERNAL_WALL_CELLS ; N_EXTERNAL_WALL_CELLS=>WM%N_EXTERNAL_WALL_CELLS
   RETURN
ENDDO
END SUBROUTINE FRAD_WALLS

!> FINE_WALLS_PTR_PROC: ON=1 points WALL, BOUNDARY_COORD, BOUNDARY_PROP1/2 and the wall counts at the hidden walls of fine box NM (no state change); ON=0 restores the box.
SUBROUTINE FWALL_POINTERS(NM,ON)
USE MESH_POINTERS, ONLY: WALL,BOUNDARY_COORD,BOUNDARY_PROP1,BOUNDARY_PROP2,N_WALL_CELLS,N_INTERNAL_WALL_CELLS,N_EXTERNAL_WALL_CELLS,POINT_TO_BOX
INTEGER, INTENT(IN) :: NM,ON
INTEGER :: LEV,IB
TYPE(MESH_TYPE), POINTER :: WM
IF (ON==0) THEN ; CALL POINT_TO_BOX(NM) ; RETURN ; ENDIF
IF (.NOT.ALLOCATED(FWALL) .OR. .NOT.ALLOCATED(FINE_LEVEL)) RETURN
DO LEV=1,MIN(SIZE(FWALL),SIZE(FINE_LEVEL))
   IB = NM-FINE_LEVEL(LEV)%NM0
   IF (IB<1 .OR. IB>FINE_LEVEL(LEV)%N_BOXES) CYCLE
   IF (.NOT.ALLOCATED(FWALL(LEV)%B)) RETURN
   WM => FWALL(LEV)%B(IB)%W
   IF (.NOT.ASSOCIATED(WM)) RETURN
   WALL=>WM%WALL ; BOUNDARY_COORD=>WM%BOUNDARY_COORD ; BOUNDARY_PROP1=>WM%BOUNDARY_PROP1 ; BOUNDARY_PROP2=>WM%BOUNDARY_PROP2
   N_WALL_CELLS=>WM%N_WALL_CELLS ; N_INTERNAL_WALL_CELLS=>WM%N_INTERNAL_WALL_CELLS ; N_EXTERNAL_WALL_CELLS=>WM%N_EXTERNAL_WALL_CELLS
   RETURN
ENDDO
END SUBROUTINE FWALL_POINTERS

!> Face velocity of the hidden OBST walls of local box IB of level LEV (e.g. a MASS_FLUX burner on the fine level): for each wall with a nonzero parent U_NORMAL_S (PRED=1)
!> or U_NORMAL, the fine cell and wall direction IJKD(:,N) = I,J,K,D and UF(N) = -SIGN(IOR)*U_NORMAL of the level-0 parent wall (the normal velocity in +axis direction, as FDS sets
!> it on a wall face). A domain-boundary SOLID wall (a vent) counts like an OBST wall; the driver re-applies its velocity after the vel_bc mirror. N <= NMAX.
SUBROUTINE FDS_FINE_WALLS_UN(LEV,IB,PRED,NMAX,N,IJKD,UF) BIND(C,NAME='fds_fine_walls_un')
USE GLOBAL_CONSTANTS, ONLY: SOLID_BOUNDARY
INTEGER(C_INT), VALUE :: LEV,IB,PRED,NMAX
INTEGER(C_INT), INTENT(OUT) :: N,IJKD(4,NMAX)
REAL(C_DOUBLE), INTENT(OUT) :: UF(NMAX)
INTEGER :: IW
REAL(EB) :: UN,V(6+2*N_TRACKED_SPECIES)
TYPE(MESH_TYPE), POINTER :: WM
N = 0
IF (.NOT.ALLOCATED(FWALL)) RETURN
IF (SIZE(FWALL)<LEV) RETURN
IF (.NOT.ALLOCATED(FWALL(LEV)%B)) RETURN
WM => FWALL(LEV)%B(IB)%W
IF (.NOT.ASSOCIATED(WM)) RETURN
DO IW=1,FWALL(LEV)%B(IB)%NW
   IF (.NOT.PARENT_STATE(FWALL(LEV)%B(IB),IW,V)) CYCLE   ! no level-0 parent state
   ASSOCIATE(BC=>WM%BOUNDARY_COORD(WM%WALL(IW)%BC_INDEX))
   UN = MERGE(V(1),V(2),PRED==1)*LS_FLUX_SCALE(FINE_LEVEL(LEV)%BOX(IB),BC,SUM(V(7+N_TRACKED_SPECIES:6+2*N_TRACKED_SPECIES)),NINT(V(6)))   ! as FDS_FINE_WALLS_STATE (the FV the projection used)
   IF (UN==0._EB .OR. N>=NMAX) CYCLE
   N = N+1
   IJKD(:,N) = (/BC%IIG,BC%JJG,BC%KKG,-BC%IOR/)
   UF(N) = -SIGN(1,BC%IOR)*UN
   END ASSOCIATE
ENDDO
END SUBROUTINE FDS_FINE_WALLS_UN

!> Parent wall state of the fine burner faces of local box IB of level LEV (frozen/fine-wall-state.md): every hidden fine wall whose level-0 parent is an internal (OBST)
!> SOLID wall with a nonzero U_NORMAL(_S) or M_DOT_G_PP_ADJUST. IJKD(:,N) = fine box-local I,J,K (1-based) of the gas cell and D (+-1..3, the face of that cell; IOR=-D).
!> S(:,N), NS = 5+2*NZ: 1 UN face velocity in +axis direction (m/s; -SIGN(IOR)*U_NORMAL_S if PRED=1, else U_NORMAL), 2 fine face area (m2), 3 fine area / parent AREA,
!> 4 RHO_F (kg/m3), 5 TMP_F (K), 6:5+NZ ZZ_F, 6+NZ:5+2*NZ M_DOT_G_PP_ADJUST (kg/m2/s, >0 into the gas). Values are the parent's at call time. Returns 0, 1 if NZ/=N_TRACKED_SPECIES or N>NMAX.
INTEGER(C_INT) FUNCTION FDS_FINE_WALLS_STATE(LEV,IB,PRED,NZ,NMAX,N,IJKD,S) BIND(C,NAME='fds_fine_walls_state')
USE GLOBAL_CONSTANTS, ONLY: N_TRACKED_SPECIES,SOLID_BOUNDARY
INTEGER(C_INT), VALUE :: LEV,IB,PRED,NZ,NMAX
INTEGER(C_INT), INTENT(OUT) :: N,IJKD(4,NMAX)
REAL(C_DOUBLE), INTENT(OUT) :: S(5+2*NZ,NMAX)
INTEGER :: IW
REAL(EB) :: UN,RDN,AREA,RAT,V(6+2*N_TRACKED_SPECIES)
TYPE(MESH_TYPE), POINTER :: WM
N = 0 ; FDS_FINE_WALLS_STATE = 0
IF (NZ/=N_TRACKED_SPECIES) THEN ; FDS_FINE_WALLS_STATE = 1 ; RETURN ; ENDIF
IF (.NOT.ALLOCATED(FWALL)) RETURN
IF (SIZE(FWALL)<LEV) RETURN
IF (.NOT.ALLOCATED(FWALL(LEV)%B)) RETURN
WM => FWALL(LEV)%B(IB)%W
IF (.NOT.ASSOCIATED(WM)) RETURN
DO IW=1,FWALL(LEV)%B(IB)%NW
   IF (.NOT.PARENT_STATE(FWALL(LEV)%B(IB),IW,V)) CYCLE   ! no level-0 parent state
   ASSOCIATE(BC=>WM%BOUNDARY_COORD(WM%WALL(IW)%BC_INDEX))
   UN = MERGE(V(1),V(2),PRED==1)
   IF (UN==0._EB .AND. ALL(V(7+NZ:6+2*NZ)==0._EB)) CYCLE
   RAT = LS_FLUX_SCALE(FINE_LEVEL(LEV)%BOX(IB),BC,SUM(V(7+NZ:6+2*NZ)),NINT(V(6)))
   IF (RAT==0._EB) CYCLE
   IF (N>=NMAX) THEN ; FDS_FINE_WALLS_STATE = 1 ; RETURN ; ENDIF
   N = N+1
   CALL FWALL_GEOM(FINE_LEVEL(LEV)%BOX(IB),BC,RDN,AREA)
   IJKD(:,N) = (/BC%IIG,BC%JJG,BC%KKG,-BC%IOR/)
   S(1:5,N) = (/-SIGN(1,BC%IOR)*UN*RAT,AREA,AREA/V(3),V(4),V(5)/)
   S(6:5+NZ,N) = V(7:6+NZ) ; S(6+NZ:5+2*NZ,N) = V(7+NZ:6+2*NZ)*RAT
   END ASSOCIATE
ENDDO
END FUNCTION FDS_FINE_WALLS_STATE

!> N_TRACKED_SPECIES for the C side (fds_fine_walls_state's NZ).
INTEGER(C_INT) FUNCTION FDS_N_TRACKED_SPECIES() BIND(C,NAME='fds_n_tracked_species')
FDS_N_TRACKED_SPECIES = N_TRACKED_SPECIES
END FUNCTION FDS_N_TRACKED_SPECIES

!> Ambient gas state of level-0 cell layer K (clamped to the mesh): RHO, TMP and the NZ tracked mass fractions, as an inflow face of FDS sets them (TMP_0, PBAR, ZZ0).
SUBROUTINE FDS_AMBIENT_STATE(K,NZ,RHO,TMP,ZZ) BIND(C,NAME='fds_ambient_state')
USE PHYSICAL_FUNCTIONS, ONLY: GET_SPECIFIC_GAS_CONSTANT
USE TYPES, ONLY: SPECIES_MIXTURE
INTEGER(C_INT), VALUE :: K,NZ
REAL(C_DOUBLE), INTENT(OUT) :: RHO,TMP,ZZ(*)
REAL(EB) :: RSUM,Z0(1:N_TRACKED_SPECIES)
INTEGER :: KK
KK = MAX(1,MIN(K,MESHES(LOWER_MESH_INDEX)%KBAR))   ! a level-0 mesh of this rank (one level-0 mesh in z: TMP_0, PBAR depend on z only)
Z0 = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0
CALL GET_SPECIFIC_GAS_CONSTANT(Z0,RSUM)
TMP = MESHES(LOWER_MESH_INDEX)%TMP_0(KK) ; RHO = MESHES(LOWER_MESH_INDEX)%PBAR(KK,0)/(RSUM*TMP) ; ZZ(1:NZ) = Z0(1:NZ)
END SUBROUTINE FDS_AMBIENT_STATE

!> Tangential ghost value G of a fine cell next to the ground, as VELOCITY_BC sets it for a WALL_MODEL_BC edge: G = SLIP_COEF*U with SLIP_COEF of WALL_MODEL at the cell centre height
!> DZ/2 (roughness of the level-0 wall of column IC,JC). The ground stress is FINE_GROUND_STRESS's; the ghost feeds the gradients (WALE near-wall viscosity, advection). MUC is unused
!> (a stress-matched ghost MUC*(U-G)/DZ = RHO*U_TAU**2 clipped to -|U| when MU was small, so the near-wall WALE MU stayed ~4x below stock over the burner).
SUBROUTINE FDS_WALL_GHOST(U,RHO,TMP,MUC,DZ,IC,JC,G) BIND(C,NAME='fds_wall_ghost')
USE TYPES, ONLY: SURFACE
USE TURBULENCE, ONLY: WALL_MODEL
USE GLOBAL_CONSTANTS, ONLY: I_MAX_TEMP,MU_RSQMW_Z,RSQ_MW_Z
REAL(C_DOUBLE), VALUE :: U,RHO,TMP,MUC,DZ
INTEGER(C_INT), VALUE :: IC,JC
REAL(C_DOUBLE), INTENT(OUT) :: G
REAL(EB) :: SLIP,UTAU,YPLUS,NU,ROUGH
INTEGER :: IW,I,J,ITMP
ROUGH = 0._EB   ! roughness of the level-0 ground wall of global level-0 column IC,JC (any mesh, any rank: GROUND_ROUGH_INIT)
IF (ALLOCATED(GROUGH)) THEN ; I = MAX(1,MIN(IC,SIZE(GROUGH,1))) ; J = MAX(1,MIN(JC,SIZE(GROUGH,2))) ; ROUGH = GROUGH(I,J) ; ENDIF
ITMP = MIN(I_MAX_TEMP,NINT(TMP))
NU = MU_RSQMW_Z(ITMP,1)/RSQ_MW_Z(1)/RHO
CALL WALL_MODEL(SLIP,UTAU,YPLUS,NU,ROUGH,0.5_EB*DZ,ABS(U))
G = SLIP*U   ! stock VELOCITY_BC (WALL_MODEL_BC, VEL_T = 0): VEL_GHOST = SLIP_COEF*VEL_GAS; the ground stress itself is FINE_GROUND_STRESS's
END SUBROUTINE FDS_WALL_GHOST

!> 1 if the case has a prevailing wind at the open boundaries (OPEN_WIND_BOUNDARY).
INTEGER(C_INT) FUNCTION FDS_OPEN_WIND() BIND(C,NAME='fds_open_wind')
USE GLOBAL_CONSTANTS, ONLY: OPEN_WIND_BOUNDARY
FDS_OPEN_WIND = MERGE(1,0,OPEN_WIND_BOUNDARY)
END FUNCTION FDS_OPEN_WIND

!> The default starting time step of read.f90 (CFL_MAX*CELL_SIZE/VEL_CHAR) for the smallest cell size CS of all levels, also as DT_INITIAL (init.f90, the
!> RESTRICT_TIME_STEP cap); -1 and no change if the input gave DT or a time ramp.
REAL(C_DOUBLE) FUNCTION FDS_DT_DEFAULT(CS) BIND(C,NAME='fds_dt_default')
USE GLOBAL_CONSTANTS, ONLY: CFL_MAX,DT_DEFAULT_VEL_CHAR,DT_INITIAL
REAL(C_DOUBLE), VALUE :: CS
FDS_DT_DEFAULT = -1._EB
IF (DT_DEFAULT_VEL_CHAR>0._EB) THEN
   FDS_DT_DEFAULT = CFL_MAX*CS/DT_DEFAULT_VEL_CHAR
   DT_INITIAL = FDS_DT_DEFAULT
ENDIF
END FUNCTION FDS_DT_DEFAULT

!> t=0 start state of fine box NM (init.f90 INITIALIZE_ATMOSPHERE): the background profile at the box's own heights (a box made by a regrid during the run keeps its
!> parent's piecewise-constant one); returns P_0 in P0(0:KBP1) for RHO = P_0(K)/(TMP*RSUM) (init.f90).
SUBROUTINE FDS_FINE_T0_PROFILE(NM,P0) BIND(C,NAME='fds_fine_t0_profile')
USE GLOBAL_CONSTANTS, ONLY: STRATIFICATION,I_RAMP_TMP0_Z,I_RAMP_P0_Z,P_INF
USE FDS_BOX_OBJ, ONLY: BOX_OBJ
USE TYPES, ONLY: INITIALIZATION,INITIALIZATION_TYPE,SPECIES_MIXTURE
USE PHYSICAL_FUNCTIONS, ONLY: GET_SPECIFIC_GAS_CONSTANT
USE GLOBAL_CONSTANTS, ONLY: N_INIT
USE MESH_POINTERS, ONLY: POINT_TO_BOX,TMPV=>TMP,ZZV=>ZZ,ZZSV=>ZZS,RSUMV=>RSUM
INTEGER(C_INT), VALUE :: NM
REAL(C_DOUBLE), INTENT(OUT) :: P0(*)
TYPE(MESH_TYPE), POINTER :: M
TYPE(INITIALIZATION_TYPE), POINTER :: IN
INTEGER :: K,I,J,N,NS
REAL(EB) :: VF(N_TRACKED_SPECIES),ZZ_GET(N_TRACKED_SPECIES)
M => BOX_OBJ(NM)
DO K=0,M%KBP1
   IF (STRATIFICATION) THEN
      M%TMP_0(K) = TMPA*EVALUATE_RAMP(M%ZC(K),I_RAMP_TMP0_Z)
      M%P_0(K)   = EVALUATE_RAMP(M%ZC(K),I_RAMP_P0_Z)
   ELSE
      M%TMP_0(K) = TMPA ; M%P_0(K) = P_INF
   ENDIF
   M%RHO_0(K) = M%P_0(K)/(M%TMP_0(K)*RSUM0)
   M%PBAR(K,:) = M%P_0(K) ; M%PBAR_S(K,:) = M%P_0(K)
   P0(K+1) = M%P_0(K)
ENDDO
! init.f90: TMP(:,:,K) = TMP_0(K), then the &INIT blocks (TMP, ZZ, RSUM) at the box's own cell centres; RHO is set by the caller
CALL POINT_TO_BOX(NM)   ! the level data (the box object's own arrays are not the registry's)
DO K=1,M%KBAR
   TMPV(1:M%IBAR,1:M%JBAR,K) = M%TMP_0(K)
ENDDO
DO N=1,N_INIT
   IN => INITIALIZATION(N)
   IF ((IN%NODE_ID/='null')) CYCLE
   IF (.NOT. (IN%ADJUST_INITIAL_CONDITIONS)) CYCLE
   DO K=1,M%KBAR
      DO J=1,M%JBAR
         DO I=1,M%IBAR
            IF (M%XC(I)<IN%X1.OR.M%XC(I)>IN%X2.OR.M%YC(J)<IN%Y1.OR.M%YC(J)>IN%Y2.OR.M%ZC(K)<IN%Z1.OR.M%ZC(K)>IN%Z2) CYCLE
            IF (IN%VOLUME_FRACTIONS_SPECIFIED) THEN
               VF = 0._EB
               DO NS=2,N_TRACKED_SPECIES
                  IF (IN%RAMP_VF_Z_INDEX(NS)>0) THEN
                     VF(NS) = EVALUATE_RAMP(M%ZC(K),IN%RAMP_VF_Z_INDEX(NS))
                  ELSE
                     VF(NS) = IN%VOLUME_FRACTION(NS)
                  ENDIF
               ENDDO
               VF(1) = 1._EB - SUM(VF)
               ZZV(I,J,K,1:N_TRACKED_SPECIES) = VF(1:N_TRACKED_SPECIES)*SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%MW / &
                                             SUM(VF(1:N_TRACKED_SPECIES)*SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%MW)
            ELSEIF (IN%MASS_FRACTIONS_SPECIFIED) THEN
               DO NS=2,N_TRACKED_SPECIES
                  IF (IN%RAMP_MF_Z_INDEX(NS)>0) THEN
                     ZZV(I,J,K,NS) = EVALUATE_RAMP(M%ZC(K),IN%RAMP_MF_Z_INDEX(NS))
                  ELSE
                     ZZV(I,J,K,NS) = IN%MASS_FRACTION(NS)
                  ENDIF
               ENDDO
               ZZV(I,J,K,1) = 1._EB - SUM(ZZV(I,J,K,2:N_TRACKED_SPECIES))
            ENDIF
            ZZSV(I,J,K,1:N_TRACKED_SPECIES) = ZZV(I,J,K,1:N_TRACKED_SPECIES)
            IF (IN%RAMP_TMP_Z_INDEX>0) THEN
               TMPV(I,J,K) = TMPM + EVALUATE_RAMP(M%ZC(K),IN%RAMP_TMP_Z_INDEX)
            ELSEIF (IN%TEMPERATURE>0._EB) THEN
               TMPV(I,J,K) = IN%TEMPERATURE
            ENDIF
            ZZ_GET(1:N_TRACKED_SPECIES) = ZZV(I,J,K,1:N_TRACKED_SPECIES)
            CALL GET_SPECIFIC_GAS_CONSTANT(ZZ_GET,RSUMV(I,J,K))
         ENDDO
      ENDDO
   ENDDO
ENDDO
END SUBROUTINE FDS_FINE_T0_PROFILE

!> INITIAL_NOISE (init.f90) random numbers of mesh number NSEED, cell size CS: the seed and the wasted calls as stock, then N amplitudes
!> NOISE_VELOCITY*(-1+2*RN)*CS in stock's loop order; returns 0 if NOISE is off.
INTEGER(C_INT) FUNCTION FDS_NOISE_RN(NSEED,CS,N,BUF) BIND(C,NAME='fds_noise_rn')
USE GLOBAL_CONSTANTS, ONLY: NOISE,NOISE_VELOCITY,RND_SEED
INTEGER(C_INT), VALUE :: NSEED,N
REAL(C_DOUBLE), VALUE :: CS
REAL(C_DOUBLE), INTENT(OUT) :: BUF(*)
INTEGER :: SIZE_RND,I
INTEGER, ALLOCATABLE :: SEED_RND(:)
REAL :: RN2
CALL RANDOM_SEED(SIZE=SIZE_RND)
ALLOCATE(SEED_RND(SIZE_RND)) ; SEED_RND = 2819 * 13*NSEED + RND_SEED
CALL RANDOM_SEED(PUT=SEED_RND) ; DEALLOCATE(SEED_RND)
DO I=1,NSEED ; CALL RANDOM_NUMBER(RN2) ; ENDDO
FDS_NOISE_RN = 0
IF (.NOT.NOISE) RETURN
DO I=1,N ; CALL RANDOM_NUMBER(RN2) ; BUF(I) = NOISE_VELOCITY*(-1._EB+2._EB*REAL(RN2,EB))*CS ; ENDDO
FDS_NOISE_RN = 1
END FUNCTION FDS_NOISE_RN

!> The starting H of init.f90, 0.5*(U0**2+V0**2+W0**2).
REAL(C_DOUBLE) FUNCTION FDS_H_INIT() BIND(C,NAME='fds_h_init')
FDS_H_INIT = 0.5_EB*(U0**2+V0**2+W0**2)
END FUNCTION FDS_H_INIT

!> The prevailing wind (U_WIND, V_WIND, W_WIND of FDS) at height Z and time T.
SUBROUTINE FDS_WIND_AT_HEIGHT(T,Z,UW,VW,WW) BIND(C,NAME='fds_wind_at_height')
USE PHYSICAL_FUNCTIONS, ONLY: GET_WIND_AT_HEIGHT
REAL(C_DOUBLE), VALUE :: T,Z
REAL(C_DOUBLE), INTENT(OUT) :: UW,VW,WW
REAL(EB) :: A,B,C
CALL GET_WIND_AT_HEIGHT(T,Z,A,B,C)
UW = A ; VW = B ; WW = C
END SUBROUTINE FDS_WIND_AT_HEIGHT

!> RDN and AREA of hidden wall BC of fine box R (as INIT_WALL_CELL).
!> --ls-fine: factor of the parent's fuel flux on a hidden fine wall of a spreading surface: (mean flux of the LS cells under the wall)/(flux of the level-0 parent). 1 when the LS grid
!> is not driving; 0 where the LS cells under the wall do not burn (the face is then no burner face).
REAL(EB) FUNCTION LS_FLUX_SCALE(R,BC,MP,SURF_INDEX)   ! MP: the parent's SUM(M_DOT_G_PP_ADJUST)
USE TYPES, ONLY: BOUNDARY_COORD_TYPE,SURFACE
TYPE(MESH_TYPE), INTENT(IN) :: R
TYPE(BOUNDARY_COORD_TYPE), INTENT(IN) :: BC
REAL(EB), INTENT(IN) :: MP
INTEGER, INTENT(IN) :: SURF_INDEX
LS_FLUX_SCALE = 1._EB
IF (.NOT.LS_FLUX_ACTIVE()) RETURN
IF (.NOT.SURFACE(SURF_INDEX)%VEG_LSET_SPREAD .OR. ABS(BC%IOR)/=3) RETURN
IF (MP<=0._EB) THEN ; LS_FLUX_SCALE = 0._EB ; RETURN ; ENDIF
LS_FLUX_SCALE = LS_BLOCK_FLUX(BC%X,BC%Y,R%DX(BC%IIG),R%DY(BC%JJG))/MP
END FUNCTION LS_FLUX_SCALE

SUBROUTINE FWALL_GEOM(R,BC,RDN,AREA)
USE TYPES, ONLY: BOUNDARY_COORD_TYPE
TYPE(MESH_TYPE), INTENT(IN) :: R
TYPE(BOUNDARY_COORD_TYPE), INTENT(IN) :: BC
REAL(EB), INTENT(OUT) :: RDN,AREA
SELECT CASE(ABS(BC%IOR))
   CASE(1) ; RDN = R%RDXN(MIN(BC%II,BC%IIG)) ; AREA = R%DY(BC%JJG)*R%DZ(BC%KKG)
   CASE(2) ; RDN = R%RDYN(MIN(BC%JJ,BC%JJG)) ; AREA = R%DX(BC%IIG)*R%DZ(BC%KKG)
   CASE(3) ; RDN = R%RDZN(MIN(BC%KK,BC%KKG)) ; AREA = R%DX(BC%IIG)*R%DY(BC%JJG)
END SELECT
END SUBROUTINE FWALL_GEOM

!> Fuel mass loss rate (kg/s, rank-local) over the LEAF burner walls of a two-level run: R(1) level-0 SOLID walls (dump.f90 WALL_LOOP2 sum, SUM_N M_DOT_G_PP_ADJUST*AREA),
!> R(2) the covered share of those (fine area/parent area), R(3) the hidden fine walls of level LEV (parent's ADJUST x fine AREA). Leaf = R(1)-R(2)+R(3). Diagnostic only.
!> R(4:6): the same three sums of the burning area (walls with a nonzero fuel rate), leaf burner area = R(4)-R(5)+R(6) (REGRID diagnostic).
SUBROUTINE FDS_FINE_FUEL_MDOT(LEV,R) BIND(C,NAME='fds_fine_fuel_mdot')
USE GLOBAL_CONSTANTS, ONLY: SOLID_BOUNDARY,MY_RANK,PROCESS
INTEGER(C_INT), VALUE :: LEV   ! unused: leaf over all levels
REAL(C_DOUBLE), INTENT(OUT) :: R(6)
INTEGER :: IB,IW,NM0,IW0,L
REAL(EB) :: RDN,AREA,MD
TYPE(MESH_TYPE), POINTER :: WM
R = 0._EB
DO NM0=1,SIZE(MESHES)
   IF (.NOT.ALLOCATED(MESHES(NM0)%CELL_INDEX)) CYCLE
   IF (PROCESS(NM0)/=MY_RANK) CYCLE   ! rank-local sum: the WALL array of another rank's mesh is not stored here
   DO IW=1,MESHES(NM0)%N_EXTERNAL_WALL_CELLS+MESHES(NM0)%N_INTERNAL_WALL_CELLS
      IF (MESHES(NM0)%WALL(IW)%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE
      ASSOCIATE(B1=>MESHES(NM0)%BOUNDARY_PROP1(MESHES(NM0)%WALL(IW)%B1_INDEX))
      MD = SUM(B1%M_DOT_G_PP_ADJUST(1:N_TRACKED_SPECIES))*B1%AREA
      R(1) = R(1) + MD
      IF (MD/=0._EB) R(4) = R(4) + B1%AREA
      IF (LEAF_WALL_COVERED(NM0,IW)) THEN ; R(2) = R(2) + MD ; IF (MD/=0._EB) R(5) = R(5) + B1%AREA ; ENDIF
      END ASSOCIATE
   ENDDO
ENDDO
IF (.NOT.ALLOCATED(FWALL)) RETURN
CALL REFRESH_LEAF()
DO L=1,SIZE(FWALL)
IF (.NOT.ALLOCATED(FWALL(L)%B)) CYCLE
DO IB=1,SIZE(FWALL(L)%B)
   WM => FWALL(L)%B(IB)%W
   IF (.NOT.ASSOCIATED(WM)) CYCLE
   DO IW=1,FWALL(L)%B(IB)%NW
      IF (.NOT.FWALL(L)%B(IB)%LEAF(IW)) CYCLE
      IF (OWN_BURNER(WM%WALL(IW)%SURF_INDEX)) THEN
         R(3) = R(3) + SUM(WM%BOUNDARY_PROP1(WM%WALL(IW)%B1_INDEX)%M_DOT_G_PP_ADJUST(1:N_TRACKED_SPECIES))*WM%BOUNDARY_PROP1(WM%WALL(IW)%B1_INDEX)%AREA
         IF (SUM(WM%BOUNDARY_PROP1(WM%WALL(IW)%B1_INDEX)%M_DOT_G_PP_ADJUST(1:N_TRACKED_SPECIES))/=0._EB) R(6) = R(6) + WM%BOUNDARY_PROP1(WM%WALL(IW)%B1_INDEX)%AREA
         CYCLE
      ENDIF
      IF (FWALL(L)%B(IB)%RKEY(IW)/=0 .AND. FWALL(L)%B(IB)%HAVES(IW)) THEN   ! parent on another rank: its state from PS
         CALL FWALL_GEOM(FINE_LEVEL(L)%BOX(IB),WM%BOUNDARY_COORD(WM%WALL(IW)%BC_INDEX),RDN,AREA)
         MD = SUM(FWALL(L)%B(IB)%PS(7+N_TRACKED_SPECIES:6+2*N_TRACKED_SPECIES,IW))
         MD = MD*LS_FLUX_SCALE(FINE_LEVEL(L)%BOX(IB),WM%BOUNDARY_COORD(WM%WALL(IW)%BC_INDEX),MD,NINT(FWALL(L)%B(IB)%PS(6,IW)))
         R(3) = R(3) + MD*AREA ; IF (MD/=0._EB) R(6) = R(6) + AREA
         CYCLE
      ENDIF
      NM0 = FWALL(L)%B(IB)%PNM(IW) ; IW0 = FWALL(L)%B(IB)%PIW(IW)   ! a wall without parent has no fuel flux (FDS_FINE_WALLS_STATE)
      IF (NM0==0) CYCLE
      IF (MESHES(NM0)%WALL(IW0)%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE
      ASSOCIATE(B1=>MESHES(NM0)%BOUNDARY_PROP1(MESHES(NM0)%WALL(IW0)%B1_INDEX),BC=>WM%BOUNDARY_COORD(WM%WALL(IW)%BC_INDEX))
      MD = SUM(B1%M_DOT_G_PP_ADJUST(1:N_TRACKED_SPECIES))
      CALL FWALL_GEOM(FINE_LEVEL(L)%BOX(IB),BC,RDN,AREA)
      MD = MD*LS_FLUX_SCALE(FINE_LEVEL(L)%BOX(IB),BC,MD,MESHES(NM0)%WALL(IW0)%SURF_INDEX)
      R(3) = R(3) + MD*AREA ; IF (MD/=0._EB) R(6) = R(6) + AREA
      END ASSOCIATE
   ENDDO
ENDDO
ENDDO
END SUBROUTINE FDS_FINE_FUEL_MDOT

!> Gas-side wall state of leaf hidden wall IW of box IB of level LEV at time T (as WALL_BC does for a wall cell of a mesh): the level-0 parent's surface state
!> (TMP_F, U_NORMAL, ZZ_F, ...; the solid stays on level 0), the fine geometry (RDN, AREA), the near-wall gas of the fine cell, and Q_CON_F
!> (SURFACE_HEAT_TRANSFER; OPEN: the INFLOW_OUTFLOW arithmetic, its ghost writes only at t=0, FWALL_T0). Pointers are left on box NM.
SUBROUTINE FWALL_STATE(LEV,IB,IW,T)
USE MESH_VARIABLES, ONLY: MESH_OBJ
USE MESH_POINTERS, ONLY: POINT_TO_BOX,UP=>U,VP=>V,WP=>W,TMP_0P=>TMP_0,PBAR_P=>PBAR,TMPG=>TMP,RHOG=>RHO,ZZG=>ZZ,RSUMG=>RSUM,RHOSG=>RHOS,ZZSG=>ZZS,USG=>US,VSG=>VS,WSG=>WS
USE TYPES, ONLY: SURFACE,SURFACE_TYPE,WALL_TYPE,BOUNDARY_COORD_TYPE,BOUNDARY_PROP1_TYPE,SPECIES_MIXTURE
USE GLOBAL_CONSTANTS, ONLY: SPECIFIED_TEMPERATURE,NET_FLUX_BC,CONVECTIVE_FLUX_BC,INFLOW_OUTFLOW,NO_CONVECTION,THERMALLY_THICK
USE PHYSICAL_FUNCTIONS, ONLY: GET_SPECIFIC_GAS_CONSTANT,HEAT_TRANSFER_COEFFICIENT
USE WALL_ROUTINES, ONLY: SURFACE_HEAT_TRANSFER
USE PRECISION_PARAMETERS, ONLY: TWENTY_EPSILON_EB
INTEGER, INTENT(IN) :: LEV,IB,IW
REAL(EB), INTENT(IN) :: T
INTEGER :: NM,NM0
REAL(EB) :: RDN,AREA,UN,RSUM_F,RAT
TYPE(MESH_TYPE), POINTER :: R,WM,M0
TYPE(WALL_TYPE), POINTER :: WC,WC0
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
TYPE(SURFACE_TYPE), POINTER :: SF
WM => FWALL(LEV)%B(IB)%W ; R => FINE_LEVEL(LEV)%BOX(IB) ; NM = FINE_LEVEL(LEV)%NM0 + IB
WC => WM%WALL(IW) ; BC => WM%BOUNDARY_COORD(WC%BC_INDEX) ; B1 => WM%BOUNDARY_PROP1(WC%B1_INDEX)
NM0 = FWALL(LEV)%B(IB)%PNM(IW) ; SF => SURFACE(WC%SURF_INDEX)
CALL FWALL_GEOM(R,BC,RDN,AREA)
CALL POINT_TO_BOX(NM)
IF (NM0>0) THEN   ! local parent: its wall data
   M0 => MESHES(NM0) ; WC0 => M0%WALL(FWALL(LEV)%B(IB)%PIW(IW))
   IF (WC%B2_INDEX>0 .AND. WC0%B2_INDEX>0) WM%BOUNDARY_PROP2(WC%B2_INDEX) = M0%BOUNDARY_PROP2(WC0%B2_INDEX)
   IF (.NOT.OWN_BURNER(WC%SURF_INDEX)) B1 = M0%BOUNDARY_PROP1(WC0%B1_INDEX)
ELSE   ! gathered state (parent on another rank, or the nearest level-0 wall of the same surface and IOR), else the surface input the fine wall was made with
   IF (.NOT.OWN_BURNER(WC%SURF_INDEX) .AND. FWALL(LEV)%B(IB)%HAVES(IW)) CALL PS_TO_B1(FWALL(LEV)%B(IB)%PS(:,IW),B1)
ENDIF
MESH_OBJ => WM
RAT = 1._EB
IF (.NOT.OWN_BURNER(WC%SURF_INDEX)) RAT = LS_FLUX_SCALE(R,BC,SUM(B1%M_DOT_G_PP_ADJUST(1:N_TRACKED_SPECIES)),WC%SURF_INDEX)
B1%RDN = RDN ; B1%AREA = AREA
IF (RAT/=1._EB) THEN ; B1%M_DOT_G_PP_ADJUST = B1%M_DOT_G_PP_ADJUST*RAT ; B1%U_NORMAL = B1%U_NORMAL*RAT ; B1%U_NORMAL_S = B1%U_NORMAL_S*RAT ; ENDIF
CALL GAS_VARIABLES(BC,B1)
SELECT CASE(SF%THERMAL_BC_INDEX)
   CASE(SPECIFIED_TEMPERATURE,NET_FLUX_BC,CONVECTIVE_FLUX_BC)
      CALL SURFACE_HEAT_TRANSFER(NM,T,SF,BC,B1,WALL_INDEX=IW)
   CASE(THERMALLY_THICK)   ! the level-0 solid's TMP_F, this face's gas and HTC (wall.f90 WALL_CELL_LOOP_0); the solid takes the leaf mean (L0Q, FDS_FINE_L0Q_COMMIT)
      IF (FWALL_T0) THEN ; MESH_OBJ => NULL() ; RETURN ; ENDIF   ! t=0 (main.f90 WALL_BC in INITIALIZATION_PHASE): no solid update yet, the parent's initial HTC and Q_CON_F (0)
      B1%HEAT_TRANS_COEF = HEAT_TRANSFER_COEFFICIENT(NM,T,B1%TMP_G-B1%TMP_F,SF,WALL_INDEX_IN=IW)
      B1%Q_CON_F = B1%HEAT_TRANS_COEF*(B1%TMP_G-B1%TMP_F)
   CASE(NO_CONVECTION)
      B1%TMP_F = B1%TMP_G
   CASE(INFLOW_OUTFLOW)
      SELECT CASE(BC%IOR)
         CASE( 1) ; UN =  UP(BC%II  ,BC%JJ,BC%KK)
         CASE(-1) ; UN = -UP(BC%II-1,BC%JJ,BC%KK)
         CASE( 2) ; UN =  VP(BC%II,BC%JJ  ,BC%KK)
         CASE(-2) ; UN = -VP(BC%II,BC%JJ-1,BC%KK)
         CASE( 3) ; UN =  WP(BC%II,BC%JJ,BC%KK  )
         CASE(-3) ; UN = -WP(BC%II,BC%JJ,BC%KK-1)
      END SELECT
      IF (UN>TWENTY_EPSILON_EB) THEN
         B1%TMP_F = TMP_0P(BC%KK) ; B1%ZZ_F = SPECIES_MIXTURE(1:N_TRACKED_SPECIES)%ZZ0
      ELSE
         B1%TMP_F = B1%TMP_G ; B1%ZZ_F = B1%ZZ_G
      ENDIF
      CALL GET_SPECIFIC_GAS_CONSTANT(B1%ZZ_F,RSUM_F)
      B1%RHO_F = PBAR_P(BC%KK,B1%PRESSURE_ZONE)/(RSUM_F*B1%TMP_F)
      B1%Q_CON_F = 2._EB*B1%K_G*(B1%TMP_G-B1%TMP_F)*B1%RDN
      IF (FWALL_STAGE>=0) THEN   ! wall.f90 INFLOW_OUTFLOW ghost cell values (predictor: RHOS, ZZS; corrector and t=0: RHO, ZZ)
         TMPG(BC%II,BC%JJ,BC%KK) = B1%TMP_F ; RSUMG(BC%II,BC%JJ,BC%KK) = RSUM_F
         IF (FWALL_STAGE==1) THEN
            ZZSG(BC%II,BC%JJ,BC%KK,1:N_TRACKED_SPECIES) = B1%ZZ_F(1:N_TRACKED_SPECIES) ; RHOSG(BC%II,BC%JJ,BC%KK) = B1%RHO_F
         ELSE
            ZZG(BC%II,BC%JJ,BC%KK,1:N_TRACKED_SPECIES) = B1%ZZ_F(1:N_TRACKED_SPECIES) ; RHOG(BC%II,BC%JJ,BC%KK) = B1%RHO_F
         ENDIF
      ENDIF
END SELECT
MESH_OBJ => NULL()
END SUBROUTINE FWALL_STATE

!> Hidden walls of fine box NM at time T, every stage after WALL_BC (patch: the fine DIVERGENCE_PART_1 reads them): leaf walls take FWALL_STATE (gas-side Q_CON_F),
!> the others (gas cell covered by the next level, counted there) Q_CON_F = 0, so a face is counted once, at its leaf level.
!> STAGE 1 (predictor of a time step): a leaf wall over a covered thermally thick level-0 wall keeps Q_CON_F*AREA*DT in PQ (each pass overwrites it) for that solid.
SUBROUTINE FDS_FINE_WALLS_QCON(NM,T,DT,STAGE) BIND(C,NAME='fds_fine_walls_qcon')
USE FDS_AMREX_HOOKS, ONLY: L0Q_COVERED
USE GLOBAL_CONSTANTS, ONLY: INITIALIZATION_PHASE,THERMALLY_THICK,SOLID_BOUNDARY
USE TYPES, ONLY: SURFACE
INTEGER(C_INT), VALUE :: NM,STAGE
REAL(C_DOUBLE), VALUE :: T,DT
INTEGER :: LEV,IB,IW,NORPH,NREMAP
REAL(EB) :: RDN,AREA
TYPE(MESH_TYPE), POINTER :: WM
IF (.NOT.ALLOCATED(FWALL)) RETURN
CALL REFRESH_LEAF()
FWALL_T0 = STAGE==0   ! the t=0 call (TimeLoop t0_divergence_levels)
FWALL_STAGE = STAGE
DO LEV=1,SIZE(FWALL)
   IF (.NOT.ALLOCATED(FWALL(LEV)%B)) CYCLE
   IB = NM - FINE_LEVEL(LEV)%NM0
   IF (IB<1 .OR. IB>SIZE(FWALL(LEV)%B)) CYCLE
   IF (.NOT.ASSOCIATED(FWALL(LEV)%B(IB)%W)) RETURN
   WM => FWALL(LEV)%B(IB)%W ; NORPH = 0 ; NREMAP = 0
   DO IW=1,FWALL(LEV)%B(IB)%NW
      IF (FWALL(LEV)%B(IB)%LEAF(IW)) THEN
         CALL FWALL_STATE(LEV,IB,IW,T)
         IF (STAGE==1 .AND. .NOT.INITIALIZATION_PHASE) THEN
            FWALL(LEV)%B(IB)%PQ(IW) = 0._EB
            IF (FWALL(LEV)%B(IB)%PNM(IW)==0 .AND. WM%WALL(IW)%BOUNDARY_TYPE==SOLID_BOUNDARY) THEN   ! fine face without a level-0 wall (ruling 2)
               CALL FWALL_GEOM(FINE_LEVEL(LEV)%BOX(IB),WM%BOUNDARY_COORD(WM%WALL(IW)%BC_INDEX),RDN,AREA)
               IF (SURFACE(WM%WALL(IW)%SURF_INDEX)%THERMAL_BC_INDEX==THERMALLY_THICK .AND. L0Q_COVERED(FWALL(LEV)%B(IB)%RNM(IW),FWALL(LEV)%B(IB)%RIW(IW))) THEN   ! interim remap
                  FWALL(LEV)%B(IB)%PQ(IW) = WM%BOUNDARY_PROP1(WM%WALL(IW)%B1_INDEX)%Q_CON_F*AREA*DT
                  REMAP_E = REMAP_E + FWALL(LEV)%B(IB)%PQ(IW) ; NREMAP = NREMAP+1
               ELSE
                  ORPHAN_E = ORPHAN_E + WM%BOUNDARY_PROP1(WM%WALL(IW)%B1_INDEX)%Q_CON_F*AREA*DT ; NORPH = NORPH+1
               ENDIF
            ENDIF
            IF (SURFACE(WM%WALL(IW)%SURF_INDEX)%THERMAL_BC_INDEX==THERMALLY_THICK .AND. L0Q_COVERED(FWALL(LEV)%B(IB)%PNM(IW),FWALL(LEV)%B(IB)%PIW(IW))) THEN
               CALL FWALL_GEOM(FINE_LEVEL(LEV)%BOX(IB),WM%BOUNDARY_COORD(WM%WALL(IW)%BC_INDEX),RDN,AREA)
               FWALL(LEV)%B(IB)%PQ(IW) = WM%BOUNDARY_PROP1(WM%WALL(IW)%B1_INDEX)%Q_CON_F*AREA*DT
            ENDIF
         ENDIF
      ELSE
         WM%BOUNDARY_PROP1(WM%WALL(IW)%B1_INDEX)%Q_CON_F = 0._EB
      ENDIF
   ENDDO
   ORPHAN_NMAX = MAX(ORPHAN_NMAX,NORPH) ; REMAP_NMAX = MAX(REMAP_NMAX,NREMAP)
   EXIT
ENDDO
FWALL_T0 = .FALSE. ; FWALL_STAGE = -1
END SUBROUTINE FDS_FINE_WALLS_QCON

!> Level-0 wall state of a PS record (FDS_FINE_WALLS_PARENT_XCHG layout) into the fine wall's B1
SUBROUTINE PS_TO_B1(P,B1X)
USE TYPES, ONLY: BOUNDARY_PROP1_TYPE
REAL(EB), INTENT(IN) :: P(:)
TYPE(BOUNDARY_PROP1_TYPE), INTENT(INOUT) :: B1X
INTEGER :: NZ
NZ = N_TRACKED_SPECIES
B1X%U_NORMAL_S = P(1) ; B1X%U_NORMAL = P(2) ; B1X%RHO_F = P(4) ; B1X%TMP_F = P(5) ; B1X%ZZ_F(1:NZ) = P(7:6+NZ) ; B1X%M_DOT_G_PP_ADJUST(1:NZ) = P(7+NZ:6+2*NZ)
B1X%T_IGN = P(7+2*NZ) ; B1X%Q_RAD_IN = P(8+2*NZ) ; B1X%RHO_D_F(1:NZ) = P(9+2*NZ:8+3*NZ)
END SUBROUTINE PS_TO_B1

!> Near-wall gas state of a wall (NEAR_SURFACE_GAS_VARIABLES for a wall without VEL_T or TMP_GAS_FRONT), from the end-of-step fields of the box the pointers point to.
SUBROUTINE GAS_VARIABLES(BCX,B1X)
USE TYPES, ONLY: BOUNDARY_COORD_TYPE,BOUNDARY_PROP1_TYPE
USE MESH_POINTERS, ONLY: UP=>U,VP=>V,WP=>W,ZZM=>ZZ,TMPM=>TMP,RHOM=>RHO
TYPE(BOUNDARY_COORD_TYPE), INTENT(IN) :: BCX
TYPE(BOUNDARY_PROP1_TYPE), INTENT(INOUT) :: B1X
REAL(EB) :: UB,VB,WB
INTEGER :: I,J,K
I = BCX%IIG ; J = BCX%JJG ; K = BCX%KKG
UB = 0.5_EB*(UP(I,J,K)+UP(I-1,J,K)) ; VB = 0.5_EB*(VP(I,J,K)+VP(I,J-1,K)) ; WB = 0.5_EB*(WP(I,J,K)+WP(I,J,K-1))
SELECT CASE(ABS(BCX%IOR))
   CASE(1) ; B1X%U_TANG = SQRT(VB**2+WB**2)
   CASE(2) ; B1X%U_TANG = SQRT(UB**2+WB**2)
   CASE(3) ; B1X%U_TANG = SQRT(UB**2+VB**2)
END SELECT
B1X%TMP_G = TMPM(I,J,K) ; B1X%RHO_G = RHOM(I,J,K) ; B1X%ZZ_G(1:N_TRACKED_SPECIES) = ZZM(I,J,K,1:N_TRACKED_SPECIES)
END SUBROUTINE GAS_VARIABLES

!> Energy check of the gas-side wall flux (run_wall_qcon_check.sh). LEV=0: the level-0 SOLID walls whose gas cell is not covered by level 1; LEV>0: the hidden walls
!> of level LEV. S(1) solid side sum(Q_CON_F*AREA) over leaf walls, S(2) gas side: what DIVERGENCE_PART_1 adds, sum(AREA_ADJUST*Q_CON_F*RDN*cell volume),
!> S(3) sum|Q_CON_F*AREA| of the non-leaf walls (must be 0 above level 0: a face counts once). Local walls only (the caller sums over ranks).
SUBROUTINE FDS_WALLS_QCON_SUMS(LEV,S) BIND(C,NAME='fds_walls_qcon_sums')
USE TYPES, ONLY: WALL_TYPE,BOUNDARY_COORD_TYPE,BOUNDARY_PROP1_TYPE
USE GLOBAL_CONSTANTS, ONLY: SOLID_BOUNDARY,MY_RANK,PROCESS
INTEGER(C_INT), VALUE :: LEV
REAL(C_DOUBLE), INTENT(OUT) :: S(3)
INTEGER :: NM0,IW,IB
REAL(EB) :: RDN,AREA,V
LOGICAL :: LF
TYPE(MESH_TYPE), POINTER :: M,R
TYPE(WALL_TYPE), POINTER :: WC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
S = 0._EB
IF (LEV==0) THEN
   DO NM0=1,NMESHES
      IF (PROCESS(NM0)/=MY_RANK) CYCLE
      M => MESHES(NM0)
      DO IW=1,M%N_EXTERNAL_WALL_CELLS+M%N_INTERNAL_WALL_CELLS
         WC => M%WALL(IW) ; IF (WC%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE
         BC => M%BOUNDARY_COORD(WC%BC_INDEX) ; B1 => M%BOUNDARY_PROP1(WC%B1_INDEX)
         LF = .TRUE. ; IF (ALLOCATED(FWALL)) LF = .NOT.POINT_COVERED(1,M%XC(BC%IIG),M%YC(BC%JJG),M%ZC(BC%KKG))
         V = M%DX(BC%IIG)*M%DY(BC%JJG)*M%DZ(BC%KKG)
         IF (LF) THEN ; S(1) = S(1) + B1%Q_CON_F*B1%AREA ; S(2) = S(2) + B1%AREA_ADJUST*B1%Q_CON_F*B1%RDN*V
         ELSE ; S(3) = S(3) + ABS(B1%Q_CON_F*B1%AREA) ; ENDIF
      ENDDO
   ENDDO
   RETURN
ENDIF
IF (.NOT.ALLOCATED(FWALL)) RETURN
IF (LEV>SIZE(FWALL)) RETURN
IF (.NOT.ALLOCATED(FWALL(LEV)%B)) RETURN
CALL REFRESH_LEAF()
DO IB=1,SIZE(FWALL(LEV)%B)
   IF (.NOT.ASSOCIATED(FWALL(LEV)%B(IB)%W)) CYCLE
   M => FWALL(LEV)%B(IB)%W ; R => FINE_LEVEL(LEV)%BOX(IB)
   DO IW=1,FWALL(LEV)%B(IB)%NW
      WC => M%WALL(IW) ; IF (WC%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE
      BC => M%BOUNDARY_COORD(WC%BC_INDEX) ; B1 => M%BOUNDARY_PROP1(WC%B1_INDEX)
      CALL FWALL_GEOM(R,BC,RDN,AREA) ; V = R%DX(BC%IIG)*R%DY(BC%JJG)*R%DZ(BC%KKG)
      IF (FWALL(LEV)%B(IB)%LEAF(IW)) THEN ; S(1) = S(1) + B1%Q_CON_F*AREA ; S(2) = S(2) + B1%AREA_ADJUST*B1%Q_CON_F*B1%RDN*V
      ELSE ; S(3) = S(3) + ABS(B1%Q_CON_F*AREA) ; ENDIF
   ENDDO
ENDDO
END SUBROUTINE FDS_WALLS_QCON_SUMS

!> Wall terms of the HRR file (UPDATE_HRR's Q_CONV, Q_COND, Q_DIFF) on the hidden walls of the local boxes of level LEV at time T, and the same terms of their level-0 parent
!> walls, weighted by the covered fraction of the parent face. The fine walls get the parent's surface state (B1 copied: TMP_F, U_NORMAL, ZZ_F, RHO_F, RHO_D_F, K_G, ...) and
!> their own near-wall gas (GAS_VARIABLES) and Q_CON_F (SURFACE_HEAT_TRANSFER; OPEN walls: the INFLOW_OUTFLOW arithmetic without its ghost writes).
!> Q(1:3) fine Q_CONV,Q_COND,Q_DIFF; Q(4:6) covered parent share; Q(7) fine sum Q_CON_F*A; Q(8) parent share of it; Q(9) max |Q_CON_F - parent Q_CON_F| when the fine wall is
!> evaluated on the parent's own data (identity check of the fine-wall path, only when CHK/=0; HTC still reads the hidden wall's BC slot, which holds the
!> parent's BC during the check); Q(10) number of fine walls; Q(11) fine non-OPEN wall area; Q(12) level-0 SOLID wall area (CHK/=0 only); Q(13) Q_CON_F*A of the leaf OPEN
!> fine walls, Q(14) of the covered level-0 OPEN walls (W).
SUBROUTINE FDS_FINE_WALLS_SUMS(LEV,TT,CHK,Q) BIND(C,NAME='fds_fine_walls_sums')
USE MESH_VARIABLES, ONLY: MESH_OBJ
USE OUTPUT_DATA, ONLY: M_DOT,N_Q_DOT
USE GLOBAL_CONSTANTS, ONLY: MY_RANK,PROCESS
USE MESH_POINTERS, ONLY: POINT_TO_BOX,UP=>U,VP=>V,WP=>W,ZZM=>ZZ,TMP_0P=>TMP_0,PBAR_P=>PBAR,TMPM=>TMP,RHOM=>RHO
USE TYPES, ONLY: SURFACE,SURFACE_TYPE,WALL_TYPE,BOUNDARY_COORD_TYPE,BOUNDARY_PROP1_TYPE,SPECIES_MIXTURE
USE GLOBAL_CONSTANTS, ONLY: SPECIFIED_TEMPERATURE,NET_FLUX_BC,CONVECTIVE_FLUX_BC,INFLOW_OUTFLOW,NO_CONVECTION,SOLID_BOUNDARY,OPEN_BOUNDARY
USE PHYSICAL_FUNCTIONS, ONLY: GET_SENSIBLE_ENTHALPY,GET_SENSIBLE_ENTHALPY_Z,GET_SPECIFIC_GAS_CONSTANT
USE WALL_ROUTINES, ONLY: SURFACE_HEAT_TRANSFER
USE PRECISION_PARAMETERS, ONLY: TWENTY_EPSILON_EB
INTEGER(C_INT), VALUE :: LEV,CHK
REAL(C_DOUBLE), VALUE :: TT
REAL(C_DOUBLE), INTENT(OUT) :: Q(14)
INTEGER :: IB,IW,NM,NM0,NMB
REAL(EB) :: T,RDN,AREA,FR,UN,QV(3),Q0(3),RSUM_F,RAT
TYPE(MESH_TYPE), POINTER :: R,WM,M0
TYPE(WALL_TYPE), POINTER :: WC,WC0
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC,BC0
TYPE(BOUNDARY_COORD_TYPE) :: BCS
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1,B10
TYPE(SURFACE_TYPE), POINTER :: SF

Q = 0._EB ; T = TT
IF (HRR_DET) CALL HRR_DET_ALLOC(N_Q_DOT,N_TRACKED_SPECIES,NMESHES)
IF (.NOT.ALLOCATED(FWALL)) RETURN
IF (LEV>SIZE(FWALL)) RETURN
IF (.NOT.ALLOCATED(FWALL(LEV)%B)) RETURN
IF (CHK/=0) THEN   ! Q(12): level-0 SOLID wall area of the local meshes (diagnostic: fine-covered fraction = Q(11)/Q(12))
   DO NM0=1,SIZE(MESHES)
      IF (.NOT.ALLOCATED(MESHES(NM0)%CELL_INDEX)) CYCLE
      IF (PROCESS(NM0)/=MY_RANK) CYCLE
      DO IW=1,MESHES(NM0)%N_EXTERNAL_WALL_CELLS+MESHES(NM0)%N_INTERNAL_WALL_CELLS
         IF (MESHES(NM0)%WALL(IW)%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE
         Q(12) = Q(12) + MESHES(NM0)%BOUNDARY_PROP1(MESHES(NM0)%WALL(IW)%B1_INDEX)%AREA
      ENDDO
   ENDDO
ENDIF
CALL REFRESH_LEAF()
IF (LEV==1) THEN   ! the level-0 walls under level 1 are not leaf: their terms (Q(4:6)) are removed by the caller; every leaf wall of every level is added (Q(1:3))
   DO NM0=1,SIZE(MESHES)
      IF (.NOT.ALLOCATED(MESHES(NM0)%CELL_INDEX) .OR. PROCESS(NM0)/=MY_RANK .OR. .NOT.ALLOCATED(MESHES(NM0)%WALL)) CYCLE
      M0 => MESHES(NM0)
      DO IW=1,M0%N_EXTERNAL_WALL_CELLS+M0%N_INTERNAL_WALL_CELLS
         WC0 => M0%WALL(IW)
         IF (WC0%BOUNDARY_TYPE/=SOLID_BOUNDARY .AND. WC0%BOUNDARY_TYPE/=OPEN_BOUNDARY) CYCLE
         BC0 => M0%BOUNDARY_COORD(WC0%BC_INDEX) ; B10 => M0%BOUNDARY_PROP1(WC0%B1_INDEX)
         IF (.NOT.POINT_COVERED(1,M0%XC(BC0%IIG),M0%YC(BC0%JJG),M0%ZC(BC0%KKG))) CYCLE
         CALL POINT_TO_BOX(NM0)
         CALL WALL_TERMS(BC0,B10,B10%AREA,Q0)
         IF (OWN_BURNER(WC0%SURF_INDEX)) THEN
            M_DOT(1:N_TRACKED_SPECIES) = M_DOT(1:N_TRACKED_SPECIES) - B10%M_DOT_G_PP_ADJUST(1:N_TRACKED_SPECIES)*B10%AREA
            IF (HRR_DET) THEN
               HRR_MADD = HRR_MADD - B10%M_DOT_G_PP_ADJUST(1:N_TRACKED_SPECIES)*B10%AREA
               HRR_MW(:,NM0) = HRR_MW(:,NM0) - B10%M_DOT_G_PP_ADJUST(1:N_TRACKED_SPECIES)*B10%AREA
            ENDIF
         ENDIF
         IF (HRR_DET) HRR_QW(:,NM0) = HRR_QW(:,NM0) - Q0
         Q(4:6) = Q(4:6) + Q0
         Q(8) = Q(8) + B10%Q_CON_F*B10%AREA
         IF (WC0%BOUNDARY_TYPE==OPEN_BOUNDARY) Q(14) = Q(14) + B10%Q_CON_F*B10%AREA
      ENDDO
   ENDDO
ENDIF
DO IB=1,SIZE(FWALL(LEV)%B)
   IF (.NOT.ASSOCIATED(FWALL(LEV)%B(IB)%W)) CYCLE
   WM => FWALL(LEV)%B(IB)%W ; R => FINE_LEVEL(LEV)%BOX(IB) ; NM = FINE_LEVEL(LEV)%NM0 + IB
   MESH_OBJ => WM
   IF (HRR_DET) NMB = L0_MESH_OF(R)
   DO IW=1,FWALL(LEV)%B(IB)%NW
      IF (.NOT.FWALL(LEV)%B(IB)%LEAF(IW)) CYCLE
      WC => WM%WALL(IW) ; BC => WM%BOUNDARY_COORD(WC%BC_INDEX) ; B1 => WM%BOUNDARY_PROP1(WC%B1_INDEX)
      NM0 = FWALL(LEV)%B(IB)%PNM(IW)
      SF => SURFACE(WC%SURF_INDEX)
      CALL FWALL_GEOM(R,BC,RDN,AREA)
      IF (NM0>0) THEN   ! local parent: its wall data
         M0 => MESHES(NM0)
         WC0 => M0%WALL(FWALL(LEV)%B(IB)%PIW(IW))
         BC0 => M0%BOUNDARY_COORD(WC0%BC_INDEX) ; B10 => M0%BOUNDARY_PROP1(WC0%B1_INDEX)
         ! identity check: the fine-wall path run on the parent's own wall data must give the parent's Q_CON_F
         IF (CHK/=0 .AND. .NOT.OWN_BURNER(WC%SURF_INDEX) .AND. (SF%THERMAL_BC_INDEX==SPECIFIED_TEMPERATURE .OR. SF%THERMAL_BC_INDEX==NET_FLUX_BC .OR. &
             SF%THERMAL_BC_INDEX==CONVECTIVE_FLUX_BC)) THEN
            CALL POINT_TO_BOX(NM0)
            B1 = B10 ; BCS = BC ; BC = BC0
            CALL SURFACE_HEAT_TRANSFER(NM0,T,SF,BC,B1,WALL_INDEX=IW)
            Q(9) = MAX(Q(9),ABS(B1%Q_CON_F-B10%Q_CON_F))
            BC = BCS
         ENDIF
      ENDIF
      CALL FWALL_STATE(LEV,IB,IW,T)   ! parent surface state, own near-wall gas and Q_CON_F
      MESH_OBJ => WM
      IF (OWN_BURNER(WC%SURF_INDEX)) THEN   ! exempt burner: own state; the HRR-file MLR takes the leaf face instead of the covered parent
         M_DOT(1:N_TRACKED_SPECIES) = M_DOT(1:N_TRACKED_SPECIES) + B1%M_DOT_G_PP_ADJUST(1:N_TRACKED_SPECIES)*AREA
         IF (HRR_DET) THEN
            HRR_MADD = HRR_MADD + B1%M_DOT_G_PP_ADJUST(1:N_TRACKED_SPECIES)*AREA
            HRR_MW(:,NMB) = HRR_MW(:,NMB) + B1%M_DOT_G_PP_ADJUST(1:N_TRACKED_SPECIES)*AREA
         ENDIF
      ENDIF
      CALL WALL_TERMS(BC,B1,AREA,QV)
      Q(1:3) = Q(1:3) + QV
      IF (HRR_DET) HRR_QW(:,NMB) = HRR_QW(:,NMB) + QV
      Q(7) = Q(7) + B1%Q_CON_F*AREA
      IF (WC%BOUNDARY_TYPE==OPEN_BOUNDARY) Q(13) = Q(13) + B1%Q_CON_F*AREA
      Q(10) = Q(10) + 1._EB
      IF (SF%THERMAL_BC_INDEX/=INFLOW_OUTFLOW) Q(11) = Q(11) + AREA
   ENDDO
   MESH_OBJ => NULL()
ENDDO

CONTAINS

!> Level-0 mesh that contains the centre of box RB (fine boxes lie in one level-0 mesh since f2253885b6).
INTEGER FUNCTION L0_MESH_OF(RB)
TYPE(MESH_TYPE), INTENT(IN) :: RB
REAL(EB) :: XC,YC,ZC
INTEGER :: K
XC = 0.5_EB*(RB%XS+RB%XF) ; YC = 0.5_EB*(RB%YS+RB%YF) ; ZC = 0.5_EB*(RB%ZS+RB%ZF)
L0_MESH_OF = 1
DO K=1,NMESHES
   IF (XC>=MESHES(K)%XS .AND. XC<MESHES(K)%XF .AND. YC>=MESHES(K)%YS .AND. YC<MESHES(K)%YF .AND. ZC>=MESHES(K)%ZS .AND. ZC<MESHES(K)%ZF) THEN
      L0_MESH_OF = K ; RETURN
   ENDIF
ENDDO
END FUNCTION L0_MESH_OF



!> UPDATE_HRR's wall-loop terms of one wall on the mesh the pointers point to: -u_n rho_f h_s A, -q_con A, -sum h_s,a rho D dZ/dn A.
SUBROUTINE WALL_TERMS(BCX,B1X,A,QX)
TYPE(BOUNDARY_COORD_TYPE), INTENT(IN) :: BCX
TYPE(BOUNDARY_PROP1_TYPE), INTENT(IN) :: B1X
REAL(EB), INTENT(IN) :: A
REAL(EB), INTENT(OUT) :: QX(3)
REAL(EB) :: U_N,H_S,H_A,HJ,ZG(1:N_TRACKED_SPECIES)
INTEGER :: N
SELECT CASE(BCX%IOR)
   CASE( 1) ; U_N = -UP(BCX%IIG-1,BCX%JJG,BCX%KKG)
   CASE(-1) ; U_N =  UP(BCX%IIG,BCX%JJG,BCX%KKG)
   CASE( 2) ; U_N = -VP(BCX%IIG,BCX%JJG-1,BCX%KKG)
   CASE(-2) ; U_N =  VP(BCX%IIG,BCX%JJG,BCX%KKG)
   CASE( 3) ; U_N = -WP(BCX%IIG,BCX%JJG,BCX%KKG-1)
   CASE DEFAULT ; U_N =  WP(BCX%IIG,BCX%JJG,BCX%KKG)
END SELECT
ZG = B1X%ZZ_F(1:N_TRACKED_SPECIES)
CALL GET_SENSIBLE_ENTHALPY(ZG,H_S,B1X%TMP_F)
HJ = 0._EB
IF (N_TRACKED_SPECIES>1) THEN
   DO N=1,N_TRACKED_SPECIES
      CALL GET_SENSIBLE_ENTHALPY_Z(N,B1X%TMP_F,H_A)
      HJ = HJ + 2._EB*H_A*B1X%RHO_D_F(N)*(ZZM(BCX%IIG,BCX%JJG,BCX%KKG,N)-B1X%ZZ_F(N))*B1X%RDN
   ENDDO
ENDIF
QX(1) = -U_N*B1X%RHO_F*H_S*A
QX(2) = -B1X%Q_CON_F*A
QX(3) = -HJ*A
END SUBROUTINE WALL_TERMS

END SUBROUTINE FDS_FINE_WALLS_SUMS

!> Fine wall-flux devices (milestone 4, registered as FINE_DEVC_PROC, called by FDS_SETUP(3) after UPDATE_DEVICES_1): a point wall DEVC
!> (SPATIAL_STATISTIC null, DV%WALL_INDEX>0) whose face point lies on a hidden level-1 wall of its level-0 wall gets that fine wall's value.
!> Supported: RADIATIVE (1), CONVECTIVE (2), TOTAL (10) HEAT FLUX, WALL TEMPERATURE (5), HEAT TRANSFER COEFFICIENT (21), as SOLID_PHASE_OUTPUT.
!> The fine B1 is the one FDS_FINE_WALLS_SUMS evaluated this step (the radiative part is the parent's: no fine radiation).
!> Composite Q_ENTH (W) of this rank for the HRR file, applied by FINE_DEVC after UPDATE_HRR: the driver sums it in box order on rank 0 (0 elsewhere), so it
!> replaces the level-0 per-rank sum of UPDATE_HRR, whose rank grouping carried round-off (rank-invariant HRR file).
SUBROUTINE FDS_SET_QENTH(V) BIND(C,NAME='fds_set_qenth')
REAL(C_DOUBLE), VALUE :: V
QENTH_V = V ; QENTH_SET = .TRUE.
FINE_DEVC_PROC => FINE_DEVC   ! FINE_DEVC applies it: register it on every rank (FDS_FINE_WALLS_BUILD registers it only on a rank with fine boxes, so a rank without
                              ! any kept its level-0 Q_ENTH in the HRR tally on top of rank 0's composite one)
END SUBROUTINE FDS_SET_QENTH

SUBROUTINE FINE_DEVC(T)
USE DEVICE_VARIABLES, ONLY: DEVICE,DEVICE_TYPE,SUBDEVICE_TYPE,N_DEVC
USE TYPES, ONLY: WALL_TYPE,BOUNDARY_COORD_TYPE,BOUNDARY_PROP1_TYPE
USE FDS_AMREX_HOOKS, ONLY: OUT_ICYC
USE GLOBAL_CONSTANTS, ONLY: TMPM
USE OUTPUT_DATA, ONLY: Q_DOT
REAL(EB), INTENT(IN) :: T
INTEGER :: N,IB,IW,Q,NM0,LEV
REAL(EB) :: V,V0,H
TYPE(DEVICE_TYPE), POINTER :: DV
TYPE(SUBDEVICE_TYPE), POINTER :: SDV
TYPE(MESH_TYPE), POINTER :: WM,R
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
IF (QENTH_SET) THEN ; Q_DOT(9) = QENTH_V ; HRR_QENTH = QENTH_V ; HRR_QENTH_ON = .TRUE. ; QENTH_SET = .FALSE. ; ENDIF
IF (.NOT.ALLOCATED(FWALL)) RETURN
IF (SIZE(FWALL)<1) RETURN
CALL REFRESH_LEAF()
DEVC_LOOP: DO N=1,N_DEVC
   DV => DEVICE(N)
   IF (DV%WALL_INDEX<=0 .OR. DV%SPATIAL_STATISTIC/='null') CYCLE
   Q = -DV%QUANTITY_INDEX(1)
   IF (Q/=1 .AND. Q/=2 .AND. Q/=5 .AND. Q/=10 .AND. Q/=21) CYCLE
   NM0 = DV%MESH
   IF (NM0<1 .OR. NM0>SIZE(MESHES)) CYCLE
   IF (DV%SUBDEVICE_INDEX(NM0)==0) CYCLE
   SDV => DV%SUBDEVICE(DV%SUBDEVICE_INDEX(NM0))
   DO LEV=SIZE(FWALL),1,-1   ! the finest level first: the leaf wall at the device point
   IF (.NOT.ALLOCATED(FWALL(LEV)%B)) CYCLE
   DO IB=1,SIZE(FWALL(LEV)%B)
      IF (.NOT.ASSOCIATED(FWALL(LEV)%B(IB)%W)) CYCLE
      WM => FWALL(LEV)%B(IB)%W ; R => FINE_LEVEL(LEV)%BOX(IB)
      DO IW=1,FWALL(LEV)%B(IB)%NW
         IF (.NOT.FWALL(LEV)%B(IB)%LEAF(IW) .OR. FWALL(LEV)%B(IB)%OWN(IW)/=NM0) CYCLE
         BC => WM%BOUNDARY_COORD(WM%WALL(IW)%BC_INDEX)
         IF (DV%IOR/=0 .AND. DV%IOR/=BC%IOR) CYCLE
         H = 0.5_EB*MAX(R%DXI,R%DETA,R%DZETA)*(1._EB+1.E-6_EB)
         IF (ABS(BC%X-DV%X)>H .OR. ABS(BC%Y-DV%Y)>H .OR. ABS(BC%Z-DV%Z)>H) CYCLE
         B1 => WM%BOUNDARY_PROP1(WM%WALL(IW)%B1_INDEX)
         V = FINE_WALL_VALUE(Q,B1)
         V0 = SDV%VALUE_1 ; SDV%VALUE_1 = V
         IF (MOD(OUT_ICYC,50)==0) WRITE(0,'(A,ES13.5,1X,A,I0,2(A,ES14.6))') 'FINE-DEVC t',T,TRIM(DV%ID)//' level ',LEV,' L0',V0,' fine',V
         CYCLE DEVC_LOOP
      ENDDO
   ENDDO
   ENDDO
ENDDO DEVC_LOOP
END SUBROUTINE FINE_DEVC

!> Value of wall quantity Q (positive solid-phase index) of a hidden fine wall with state B1: the five quantities FINE_DEVC supports (radiative, convective, total heat flux, wall
!> temperature, HTC); other quantities are the parent's (copied state), so their wall statistics are leaf already.
REAL(EB) FUNCTION FINE_WALL_VALUE(Q,B1)
USE TYPES, ONLY: BOUNDARY_PROP1_TYPE
USE GLOBAL_CONSTANTS, ONLY: TMPM
INTEGER, INTENT(IN) :: Q
TYPE(BOUNDARY_PROP1_TYPE), INTENT(IN) :: B1
SELECT CASE(Q)
   CASE(1)  ; FINE_WALL_VALUE = (B1%Q_RAD_IN-B1%Q_RAD_OUT)*0.001_EB
   CASE(2)  ; FINE_WALL_VALUE = B1%Q_CON_F*0.001_EB
   CASE(5)  ; FINE_WALL_VALUE = B1%TMP_F-TMPM
   CASE(10) ; FINE_WALL_VALUE = (B1%Q_RAD_IN-B1%Q_RAD_OUT+B1%Q_CON_F)*0.001_EB
   CASE(21) ; FINE_WALL_VALUE = B1%HEAT_TRANS_COEF
END SELECT
END FUNCTION FINE_WALL_VALUE

!> BURNING RATE (Q=7) and HRRPUA (Q=9) of fine wall IW as dump.f90 SOLID_PHASE_OUTPUT, on its PARENT_STATE (own state on an exempt burner, else the parent's); 0 without one.
REAL(EB) FUNCTION FUEL_VALUE(Q,FB,IW)
USE TYPES, ONLY: REACTION
USE GLOBAL_CONSTANTS, ONLY: N_REACTIONS
INTEGER, INTENT(IN) :: Q,IW
TYPE(FWALL_BOX_TYPE), INTENT(IN) :: FB
REAL(EB) :: PV(6+2*N_TRACKED_SPECIES)
INTEGER :: NR
FUEL_VALUE = 0._EB
IF (.NOT.PARENT_STATE(FB,IW,PV)) RETURN
DO NR=1,N_REACTIONS
   FUEL_VALUE = FUEL_VALUE + PV(6+N_TRACKED_SPECIES+REACTION(NR)%FUEL_SMIX_INDEX)*MERGE(1._EB,REACTION(NR)%HOC_COMPLETE*0.001_EB,Q==7)
ENDDO
END FUNCTION FUEL_VALUE

!> Hooks of the wall-device statistics (see FDS_AMREX_HOOKS): QI is DV%QUANTITY_INDEX(1) (negative for a solid-phase quantity).
!> Leaf wall devices (any number of levels): the hidden walls whose gas cell is not covered by the next level and lies in level-0 mesh NM are counted in the statistics of NM
!> instead of the covered level-0 walls (LEAF_WALL_COVERED). A fine wall over a level-0 mesh of another rank is not visited here.
INTEGER FUNCTION LEAF_WALL_COUNT(NM,QI)
INTEGER, INTENT(IN) :: NM,QI
LEAF_WALL_COUNT = 0
IF (QI/=-1 .AND. QI/=-2 .AND. QI/=-5 .AND. QI/=-7 .AND. QI/=-9 .AND. QI/=-10 .AND. QI/=-21) RETURN
IF ((QI==-7 .OR. QI==-9) .AND. .NOT.BED_RUN) RETURN   ! leaf burner sums: bed runs (exempt burners under a fine level)
IF (.NOT.ALLOCATED(FWALL)) RETURN
CALL REFRESH_LEAF()
IF (NM>=1 .AND. NM<=UBOUND(LW_OFF,1)) LEAF_WALL_COUNT = LW_OFF(NM)-LW_OFF(NM-1)
END FUNCTION LEAF_WALL_COUNT

INTEGER FUNCTION LEAF_WALL(NM,K,QI,WC,BC,B1,V,IWP)
USE TYPES, ONLY: WALL_TYPE,BOUNDARY_COORD_TYPE,BOUNDARY_PROP1_TYPE
INTEGER, INTENT(IN) :: NM,K,QI
TYPE(WALL_TYPE), POINTER :: WC
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
REAL(EB), INTENT(OUT) :: V
INTEGER, INTENT(OUT) :: IWP
INTEGER :: LEV,IB,IW
CALL REFRESH_LEAF()
LEAF_WALL = 0 ; V = 0._EB ; IWP = 0
IF (NM<1 .OR. NM>UBOUND(LW_OFF,1)) RETURN
IF (K<1 .OR. K>LW_OFF(NM)-LW_OFF(NM-1)) RETURN
LEV = LW_LEV(LW_OFF(NM-1)+K) ; IB = LW_BOX(LW_OFF(NM-1)+K) ; IW = LW_W(LW_OFF(NM-1)+K)
ASSOCIATE(WM=>FWALL(LEV)%B(IB)%W)
WC => WM%WALL(IW) ; BC => WM%BOUNDARY_COORD(WC%BC_INDEX) ; B1 => WM%BOUNDARY_PROP1(WC%B1_INDEX)
END ASSOCIATE
V = FINE_WALL_VALUE(-QI,B1) ; IWP = FWALL(LEV)%B(IB)%PIW(IW) ; LEAF_WALL = FINE_LEVEL(LEV)%NM0+IB
IF (QI==-7 .OR. QI==-9) V = FUEL_VALUE(-QI,FWALL(LEV)%B(IB),IW)
END FUNCTION LEAF_WALL

!> Leaf walls of fine box NBX (global box number): count, and the K-th (BOUNDARY_COORD, owner level-0 mesh, parent wall index); O(1) per call (BND_STAGE_RATE).
INTEGER FUNCTION LEAF_BOXWALL_COUNT(NBX)
INTEGER, INTENT(IN) :: NBX
INTEGER :: LEV,IB
LEAF_BOXWALL_COUNT = 0
IF (.NOT.ALLOCATED(FWALL)) RETURN
CALL REFRESH_LEAF()
DO LEV=1,SIZE(FWALL)
   IF (.NOT.ALLOCATED(FWALL(LEV)%B)) CYCLE
   IB = NBX-FINE_LEVEL(LEV)%NM0
   IF (IB>=1 .AND. IB<=SIZE(FWALL(LEV)%B)) THEN
      IF (ASSOCIATED(FWALL(LEV)%B(IB)%W)) LEAF_BOXWALL_COUNT = SIZE(FWALL(LEV)%B(IB)%LIDX)
      RETURN
   ENDIF
ENDDO
END FUNCTION LEAF_BOXWALL_COUNT

SUBROUTINE LEAF_BOXWALL(NBX,K,BC,OWNM,IWP,B1)
USE TYPES, ONLY: BOUNDARY_COORD_TYPE,BOUNDARY_PROP1_TYPE
INTEGER, INTENT(IN) :: NBX,K
TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
INTEGER, INTENT(OUT) :: OWNM,IWP
TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
INTEGER :: LEV,IB,IW
OWNM = 0 ; IWP = 0
DO LEV=1,SIZE(FWALL)
   IF (.NOT.ALLOCATED(FWALL(LEV)%B)) CYCLE
   IB = NBX-FINE_LEVEL(LEV)%NM0
   IF (IB>=1 .AND. IB<=SIZE(FWALL(LEV)%B)) THEN
      IW = FWALL(LEV)%B(IB)%LIDX(K)
      BC => FWALL(LEV)%B(IB)%W%BOUNDARY_COORD(FWALL(LEV)%B(IB)%W%WALL(IW)%BC_INDEX)
      OWNM = FWALL(LEV)%B(IB)%OWN(IW) ; IWP = FWALL(LEV)%B(IB)%PIW(IW)
      B1 => FWALL(LEV)%B(IB)%W%BOUNDARY_PROP1(FWALL(LEV)%B(IB)%W%WALL(IW)%B1_INDEX)
      IF (OWN_BURNER(FWALL(LEV)%B(IB)%W%WALL(IW)%SURF_INDEX)) IWP = -1
      RETURN
   ENDIF
ENDDO
END SUBROUTINE LEAF_BOXWALL

!> Level-0 wall IW of mesh NM is not leaf: its gas cell lies inside a level-1 box (its finer walls, if any, are counted by LEAF_WALL).
LOGICAL FUNCTION LEAF_WALL_COVERED(NM,IW)
INTEGER, INTENT(IN) :: NM,IW
LEAF_WALL_COVERED = .FALSE.
ASSOCIATE(M=>MESHES(NM),BC=>MESHES(NM)%BOUNDARY_COORD(MESHES(NM)%WALL(IW)%BC_INDEX))
LEAF_WALL_COVERED = POINT_COVERED(1,M%XC(BC%IIG),M%YC(BC%JJG),M%ZC(BC%KKG))
END ASSOCIATE
END FUNCTION LEAF_WALL_COVERED

END MODULE FDS_FINE_LEVEL
