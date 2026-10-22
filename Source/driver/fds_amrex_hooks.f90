!> \brief Hook module of the guarded FDS edits (USE_AMREX=ON only; patches 0003, 0004, 0005 in Source/driver/patches).
!>
!> EXTERNAL_GHOSTS_FILLED: when .TRUE. the ghost cells across a box boundary (EXTERNAL_WALL(IW)%NOM>0) have been filled from AMReX data before the boundary
!> routines run, so the NOM>0 branches of VISCOSITY_BC, NO_FLUX, VELOCITY_BC, MATCH_VELOCITY, MATCH_VELOCITY_FLUX (velo.f90) and ASSIGN_GHOST_VALUE (wall.f90)
!> must not overwrite them with OMESH values ("filled externally", ADR-001 Option C). Default .FALSE.: FDS behaves as before. S4 keeps it .FALSE. (the driver
!> fills OMESH and runs the unmodified routines, GhostExchange.cpp); the TRUE path is exercised from S5 on, once the box views below are in use.
!>
!> BOX_VIEW(NM): pointer views of the AMReX data of box NM, with the FDS index bounds (lower bounds as in MESHES(NM)%X), so that POINT_TO_BOX (patch 0005,
!> mesh.f90) can point the MESH_POINTERS module pointers at box data without going through the MESHES(NM)%X allocatables. FDS_HOOK_SET_VIEW builds a view from
!> the C address of a FAB with the standard C_F_POINTER plus pointer bounds remapping (Fortran 2008, no compiler extension).
!>
!> Kernel-facing file rules (M2a):
!> (a) Passive scalars (N_TOTAL_SCALARS beyond the tracked species) are handled by Fields.cpp: ZZ/ZZS views carry the last extent N_TOTAL_SCALARS as given by the caller.
!> (b) Only uniform Cartesian metrics are used; R(I)/RRN(I) are not involved and CYLINDRICAL/TRN* meshes are rejected at level-0 assembly (IR-002).
MODULE FDS_AMREX_HOOKS

USE ISO_C_BINDING
USE PRECISION_PARAMETERS, ONLY: EB
USE TYPES, ONLY: WALL_TYPE,BOUNDARY_COORD_TYPE,BOUNDARY_PROP1_TYPE

IMPLICIT NONE (TYPE,EXTERNAL)
PRIVATE

LOGICAL, PUBLIC, SAVE :: EXTERNAL_GHOSTS_FILLED = .FALSE.

!> Verification only (env FDSTL_VERIFY_STOCK_EDGE_LAG=1, off by default): stock VELOCITY_BC's ground-edge (vegetation) branch reads the tangential velocity of an
!> interface ghost cell before the same edge loop refreshes it, i.e. the value left by the previous call. fds_p_lag_save (fds_step.f90) keeps U,V of box NM (US,VS when EST=1)
!> just before the driver's ghost fill; STOCK_EDGE_LAG_VT gives that edge the saved values in the ghost cells, which reproduces stock bitwise.
!> The switch also keeps the diagonal ghosts (outside the box in two or three directions, across a mesh-corner line) out of the fill (GhostExchange.cpp
!> ghost_exchange): stock does not refresh them (cell fields, US,VS,WS) or refreshes them a step late (U,V,W, VELOCITY_BC corrector branch).
LOGICAL, PUBLIC, SAVE :: VERIFY_STOCK_EDGE_LAG = .FALSE.
!> The OPEN-boundary branch of the same loop has the same lag: its inflow test (U_NORM) reads interface ghosts that a later edge, or a later orientation of the same
!> edge, refreshes (STOCK_EDGE_LAG_UN). Read and write order is the key IE*8+P of the edge loop, P = 3*(I_SGN+1)/2+IS-1 (SIGN_LOOP, ORIENTATION_LOOP).
TYPE :: LAG_T
   REAL(EB), ALLOCATABLE :: U(:,:,:),V(:,:,:),W(:,:,:)
   INTEGER, ALLOCATABLE :: KU(:,:,:),KV(:,:,:),KW(:,:,:)   !< slot 0: first key of stock VELOCITY_BC's 'ghost cell values at the edge' write of that ghost U, V, W (HUGE: none)
END TYPE LAG_T
TYPE(LAG_T), ALLOCATABLE, SAVE :: LAG(:,:)
PUBLIC :: LAG_T,LAG,FDS_HOOK_LAG_SAVE,STOCK_EDGE_LAG_VT,STOCK_EDGE_LAG_UN

!> Time, time step and cycle number that the driver hands to FDS_SETUP(MODE=3) (end-of-step outputs, patch 0006).
REAL(EB), PUBLIC, SAVE :: OUT_T = 0._EB, OUT_DT = 0._EB
INTEGER,  PUBLIC, SAVE :: OUT_ICYC = 0
!> Set .TRUE. by the patched FDS_SETUP (patch 0006) when set-up returns; the driver calls MODE=3 only if it is .TRUE. (an unpatched main.f90 would run END_FDS).
LOGICAL,  PUBLIC, SAVE :: STEP_OUTPUTS_PATCHED = .FALSE.
!> Patch 0011: when .TRUE. the patched FDS_SETUP(MODE=3) also runs the per-mesh dump loop of MAIN_LOOP (DUMP_MESH_OUTPUTS: slices, boundary files, Plot3D, ...). Default .FALSE.
!> (set by fds_hook_set_mesh_dumps), so the behaviour of patch 0006 is unchanged until the driver asks.
LOGICAL,  PUBLIC, SAVE :: OUT_MESH_DUMPS = .FALSE.
!> Patch 0012: number of output meshes (the meshes of the FDS output files, which may exceed NMESHES once boxes of finer levels are output meshes). 0 (default) = NMESHES.
!> Must be set (fds_hook_set_out_meshes) before fds_setup(0), because ASSIGN_FILE_NAMES sizes the output tables.
INTEGER,  PUBLIC, SAVE :: N_OUT_MESHES = 0

TYPE, PUBLIC :: BOX_VIEW_TYPE
   REAL(EB), POINTER, DIMENSION(:,:,:)   :: U=>NULL(),V=>NULL(),W=>NULL(),US=>NULL(),VS=>NULL(),WS=>NULL(),D=>NULL(),DS=>NULL(),H=>NULL(),HS=>NULL(), &
                                            KRES=>NULL(),FVX=>NULL(),FVY=>NULL(),FVZ=>NULL(),RHO=>NULL(),RHOS=>NULL(),MU=>NULL(),TMP=>NULL(),Q=>NULL(),RSUM=>NULL(),DDDT=>NULL()
   REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZ=>NULL(),ZZS=>NULL()
END TYPE BOX_VIEW_TYPE

TYPE(BOX_VIEW_TYPE), ALLOCATABLE, TARGET, PUBLIC, SAVE :: BOX_VIEW(:)

PUBLIC :: FDS_HOOK_SET_FLAG,FDS_HOOK_SET_VIEW,FDS_HOOK_SET_STEP,FDS_HOOK_SET_MESH_DUMPS,FDS_HOOK_SET_OUT_MESHES,FDS_HOOK_STEP_OUTPUTS,FDS_HOOK_FINE_ABORT,FDS_HOOK_FINE_GUARD,FDS_HOOK_SET_FINE_READY,FDS_HOOK_SHADOW,SHADOW_PROC,SHADOW_IF,FDS_HOOK_BIND_VIEW,FDS_HOOK_L0_ONLY

!> D-056 (option B): the set of kernel wrappers (fds_kernels.f90 entry names) that may run on a fine-level mesh number (NM > NMESHES). Empty until patch 0007 is validated and the
!> kernels are switched to POINT_TO_BOX: every wrapper aborts on a fine number now. Set by FDS_HOOK_FINE_READY (also callable from the draft driver file fds_fine_mesh_b.f90).
LOGICAL, SAVE :: FINE_READY = .FALSE.

!> Shadow hook (fine-level box test, D-056 option B): when a procedure is registered (draft file fds_fine_box_b.f90, option FDS_AMR_FINE_B_DRAFT), every kernel wrapper calls it before (PHASE=0)
!> and after (PHASE=1) the kernel on a level-0 box, so that the same kernel can be run on a fine-level copy through POINT_TO_BOX and compared. Not registered in a normal run.
ABSTRACT INTERFACE
   SUBROUTINE SHADOW_IF(PHASE,KIND,NM,T,DT,EST,DTNEW,ICHG)
      IMPORT :: EB
      INTEGER, INTENT(IN) :: PHASE,KIND,NM,EST,ICHG
      REAL(EB), INTENT(IN) :: T,DT,DTNEW
   END SUBROUTINE SHADOW_IF
END INTERFACE
PROCEDURE(SHADOW_IF), POINTER, SAVE :: SHADOW_PROC => NULL()

!> Fine wall-flux devices (patch 0027): when registered (fds_fine_level.f90, with hidden fine walls), FDS_SETUP(3) calls it after the
!> UPDATE_GLOBAL_OUTPUTS loop, so that wall devices on a fine-covered face read the hidden fine wall instead of level 0.
ABSTRACT INTERFACE
   SUBROUTINE FINE_DEVC_IF(T)
      IMPORT :: EB
      REAL(EB), INTENT(IN) :: T
   END SUBROUTINE FINE_DEVC_IF
END INTERFACE
PROCEDURE(FINE_DEVC_IF), POINTER, PUBLIC, SAVE :: FINE_DEVC_PROC => NULL()

!> Leaf-cell hooks of the gas devices (UPDATE_DEVICES_1, registered by fds_leaf_devc.f90 when a fine level exists): LEAF_COVERED(NM,X,Y,Z) is true for a point of box NM inside a finer box;
!> LEAF_BOX(IB,X1,X2,Y1,Y2,Z1,Z2,IR) is the IB-th fine box number (0 after the last) with the cell range IR(6) of the extent in it; LEAF_FINEST(X,Y,Z,I,J,K) the finest fine box at a point (0 if none).
ABSTRACT INTERFACE
   LOGICAL FUNCTION LEAF_COVERED_IF(NM,X,Y,Z)
      IMPORT :: EB
      INTEGER, INTENT(IN) :: NM
      REAL(EB), INTENT(IN) :: X,Y,Z
   END FUNCTION LEAF_COVERED_IF
   INTEGER FUNCTION LEAF_BOX_IF(IB,X1,X2,Y1,Y2,Z1,Z2,IR)
      IMPORT :: EB
      INTEGER, INTENT(IN) :: IB
      REAL(EB), INTENT(IN) :: X1,X2,Y1,Y2,Z1,Z2
      INTEGER, INTENT(OUT) :: IR(6)
   END FUNCTION LEAF_BOX_IF
   INTEGER FUNCTION LEAF_FINEST_IF(X,Y,Z,I,J,K)
      IMPORT :: EB
      REAL(EB), INTENT(IN) :: X,Y,Z
      INTEGER, INTENT(OUT) :: I,J,K
   END FUNCTION LEAF_FINEST_IF
END INTERFACE
PROCEDURE(LEAF_COVERED_IF), POINTER, PUBLIC, SAVE :: LEAF_COVERED_PROC => NULL()
PROCEDURE(LEAF_BOX_IF), POINTER, PUBLIC, SAVE :: LEAF_BOX_PROC => NULL()
PROCEDURE(LEAF_FINEST_IF), POINTER, PUBLIC, SAVE :: LEAF_FINEST_PROC => NULL()

!> Leaf wall devices (UPDATE_DEVICES_1 SPATIAL_STATISTIC over walls; registered with the hidden fine walls, fds_fine_level.f90): for a quantity index QI whose hidden fine wall value
!> differs from its parent's (the FINE_DEVC set) LEAF_WALL_COUNT(NM,QI) is the number of hidden fine walls on level-0 mesh NM, LEAF_WALL(NM,K,QI,WC,BC,B1,V,IWP) the K-th (its WALL, BOUNDARY_COORD and
!> BOUNDARY_PROP1 entries, value V, parent wall index IWP, box number as result) and LEAF_WALL_COVERED(NM,IW) tells that level-0 wall IW has fine children (its share is then the fine walls').
ABSTRACT INTERFACE
   INTEGER FUNCTION LEAF_WALL_COUNT_IF(NM,QI)
      INTEGER, INTENT(IN) :: NM,QI
   END FUNCTION LEAF_WALL_COUNT_IF
   INTEGER FUNCTION LEAF_WALL_IF(NM,K,QI,WC,BC,B1,V,IWP)
      IMPORT :: EB,WALL_TYPE,BOUNDARY_COORD_TYPE,BOUNDARY_PROP1_TYPE
      INTEGER, INTENT(IN) :: NM,K,QI
      TYPE(WALL_TYPE), POINTER :: WC
      TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
      TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
      REAL(EB), INTENT(OUT) :: V
      INTEGER, INTENT(OUT) :: IWP
   END FUNCTION LEAF_WALL_IF
   LOGICAL FUNCTION LEAF_WALL_COVERED_IF(NM,IW)
      INTEGER, INTENT(IN) :: NM,IW
   END FUNCTION LEAF_WALL_COVERED_IF
   INTEGER FUNCTION LEAF_BOXWALL_COUNT_IF(NBX)
      INTEGER, INTENT(IN) :: NBX
   END FUNCTION LEAF_BOXWALL_COUNT_IF
   SUBROUTINE LEAF_BOXWALL_IF(NBX,K,BC,OWNM,IWP,B1)
      IMPORT :: BOUNDARY_COORD_TYPE,BOUNDARY_PROP1_TYPE
      INTEGER, INTENT(IN) :: NBX,K
      TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
      INTEGER, INTENT(OUT) :: OWNM,IWP
      TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
   END SUBROUTINE LEAF_BOXWALL_IF
END INTERFACE
!> Leaf walls of one fine box (O(1) per wall, BND_STAGE_RATE): count and K-th (BOUNDARY_COORD, owner level-0 mesh, parent wall index; -1 and the wall's own B1 on an exempt burner).
PROCEDURE(LEAF_BOXWALL_COUNT_IF), POINTER, PUBLIC, SAVE :: LEAF_BOXWALL_COUNT_PROC => NULL()
PROCEDURE(LEAF_BOXWALL_IF), POINTER, PUBLIC, SAVE :: LEAF_BOXWALL_PROC => NULL()
PROCEDURE(LEAF_WALL_COUNT_IF), POINTER, PUBLIC, SAVE :: LEAF_WALL_COUNT_PROC => NULL()
PROCEDURE(LEAF_WALL_IF), POINTER, PUBLIC, SAVE :: LEAF_WALL_PROC => NULL()
PROCEDURE(LEAF_WALL_COVERED_IF), POINTER, PUBLIC, SAVE :: LEAF_WALL_COVERED_PROC => NULL()

!> Radiation walls of a fine box (registered by fds_fine_level.f90): RAD_WALLS_PROC(NM,1) makes the wall pointers of COMPUTE_RADIATION (WALL, BOUNDARY_*, WALL_WORK, wall counts)
!> point to the hidden fine walls of box NM (their B1 state copied from the level-0 parent wall); RAD_WALLS_PROC(NM,0) points them back to the box.
ABSTRACT INTERFACE
   SUBROUTINE RAD_WALLS_IF(NM,ON)
      INTEGER, INTENT(IN) :: NM,ON
   END SUBROUTINE RAD_WALLS_IF
END INTERFACE
PROCEDURE(RAD_WALLS_IF), POINTER, PUBLIC, SAVE :: RAD_WALLS_PROC => NULL()
!> Fine box NM: ON=1 points the WALL/BOUNDARY_* mesh pointers and wall counts at its hidden walls (velo.f90 near-wall viscosity), ON=0 restores the box.
PROCEDURE(RAD_WALLS_IF), POINTER, PUBLIC, SAVE :: FINE_WALLS_PTR_PROC => NULL()
!> Ground walls of the ghost cells of fine box NM (fds_fine_ground_ghosts, filled from the neighbour fine box): GSI(I,J) = SURF_INDEX of the ground
!> wall of cell (I,J) (0: a fine cell without one; -1: not a fine cell, i.e. coarse-fine or outside the domain), GRDN its RDN. Indices 0:IBP1, 0:JBP1.
TYPE, PUBLIC :: FACE_OWN_TYPE
   INTEGER, ALLOCATABLE :: GSI(:,:)
   REAL(EB), ALLOCATABLE :: GRDN(:,:)
END TYPE FACE_OWN_TYPE
TYPE(FACE_OWN_TYPE), ALLOCATABLE, PUBLIC, SAVE :: FACE_OWN(:)

!> Rank-invariant HRR file (driver runs, HRR_DET): the step terms of each level-0 mesh (UPDATE_HRR: HRR_QM, HRR_MM; leaf wall corrections: HRR_QW, HRR_MW) are
!> tallied per mesh (HRR_QDM, HRR_MDM); rank-level terms that have no mesh (particles; composite Q_ENTH of rank 0) go to HRR_QDR, HRR_MDR. At an HRR dump the
!> tallies are reduced (one source per entry) and summed in mesh order, then rank order, on rank 0 (main.f90). HRR_QADD, HRR_MADD: the per-rank driver
!> additions already in Q_DOT, M_DOT (removed from the rank residual).
LOGICAL,  PUBLIC, SAVE :: HRR_DET = .FALSE., HRR_QENTH_ON = .FALSE.
REAL(EB), PUBLIC, SAVE :: HRR_QENTH = 0._EB
REAL(EB), ALLOCATABLE, PUBLIC, SAVE :: HRR_QM(:,:),HRR_MM(:,:),HRR_QW(:,:),HRR_MW(:,:),HRR_QDM(:,:),HRR_MDM(:,:),HRR_QDR(:),HRR_MDR(:),HRR_QADD(:),HRR_MADD(:)
PUBLIC :: HRR_DET_ALLOC

!> Covered level-0 walls (Solid Phase Lead conditions 1, 2): per level-0 mesh NM and wall IW, COV = a SOLID wall whose gas cell lies under level 1 (REFRESH_LEAF).
!> Its 1-D solid (SOLID_HEAT_TRANSFER) is driven by the leaf fine faces over it: E = sum Q_CON_F*AREA*DT of their predictor stages (last pass), committed at the start
!> of each corrector (FDS_FINE_L0Q_COMMIT), taken and zeroed by the next solid update. The level-0 wall term of the covered gas cell is left out of
!> DIVERGENCE_PART_1 (the fine level carries it). GAS_CUM / SOL_CUM: committed fine energy and energy taken by the level-0 solids (local, run_wall_qcon_check.sh).
TYPE L0Q_TYPE
   LOGICAL, ALLOCATABLE :: COV(:)
   REAL(EB), ALLOCATABLE :: E(:)
END TYPE
TYPE(L0Q_TYPE), ALLOCATABLE, PUBLIC, SAVE :: L0Q(:)
REAL(EB), PUBLIC, SAVE :: L0Q_GAS_CUM = 0._EB, L0Q_SOL_CUM = 0._EB, L0Q_FLUSH_CUM = 0._EB   !< FLUSH: the part taken by uncovered walls (ruling 1)
REAL(EB), PUBLIC, SAVE :: L0Q_GND_CUM = 0._EB, L0Q_GND_ABS = 0._EB   !< convective energy into the level-0 thick ground (IOR 3) solids, Q_CON_F*AREA*DT_BC per update (ruling 2 remap ratio)
PUBLIC :: L0Q_COVERED,L0Q_PENDING

CONTAINS

!> Level-0 wall IW of mesh NM is covered by level 1 (L0Q table; .FALSE. without fine levels)
LOGICAL FUNCTION L0Q_COVERED(NM,IW)
INTEGER, INTENT(IN) :: NM,IW
L0Q_COVERED = .FALSE.
IF (.NOT.ALLOCATED(L0Q)) RETURN
IF (NM<1 .OR. NM>SIZE(L0Q)) RETURN
IF (.NOT.ALLOCATED(L0Q(NM)%COV)) RETURN
IF (IW<1 .OR. IW>SIZE(L0Q(NM)%COV)) RETURN
L0Q_COVERED = L0Q(NM)%COV(IW)
END FUNCTION L0Q_COVERED

!> Fine energy still held for level-0 wall IW of mesh NM (0 if none): after a regrid uncovered the wall, its next solid update takes it (ruling 1)
REAL(EB) FUNCTION L0Q_PENDING(NM,IW)
INTEGER, INTENT(IN) :: NM,IW
L0Q_PENDING = 0._EB
IF (.NOT.ALLOCATED(L0Q)) RETURN
IF (NM<1 .OR. NM>SIZE(L0Q)) RETURN
IF (.NOT.ALLOCATED(L0Q(NM)%E)) RETURN
IF (IW<1 .OR. IW>SIZE(L0Q(NM)%E)) RETURN
L0Q_PENDING = L0Q(NM)%E(IW)
END FUNCTION L0Q_PENDING

!> Allocate (once) the HRR_DET arrays: NQ terms, NS tracked species, NM level-0 meshes.
SUBROUTINE HRR_DET_ALLOC(NQ,NS,NM)
INTEGER, INTENT(IN) :: NQ,NS,NM
IF (ALLOCATED(HRR_QM)) RETURN
ALLOCATE(HRR_QM(NQ,NM),HRR_MM(NS,NM),HRR_QW(3,NM),HRR_MW(NS,NM),HRR_QDM(NQ,NM),HRR_MDM(NS,NM),HRR_QDR(NQ),HRR_MDR(NS),HRR_QADD(NQ),HRR_MADD(NS))
HRR_QM = 0._EB ; HRR_MM = 0._EB ; HRR_QW = 0._EB ; HRR_MW = 0._EB ; HRR_QDM = 0._EB ; HRR_MDM = 0._EB
HRR_QDR = 0._EB ; HRR_MDR = 0._EB ; HRR_QADD = 0._EB ; HRR_MADD = 0._EB
END SUBROUTINE HRR_DET_ALLOC

!> Set or clear EXTERNAL_GHOSTS_FILLED.
SUBROUTINE FDS_HOOK_SET_FLAG(FLAG) BIND(C,NAME='fds_hook_set_flag')
INTEGER(C_INT), VALUE :: FLAG
EXTERNAL_GHOSTS_FILLED = (FLAG/=0)
END SUBROUTINE FDS_HOOK_SET_FLAG

!> 1 when main.f90 carries patch 0006 (FDS_SETUP(MODE=3) end-of-step outputs), else 0.
FUNCTION FDS_HOOK_STEP_OUTPUTS() BIND(C,NAME='fds_hook_step_outputs') RESULT(F)
INTEGER(C_INT) :: F
F = MERGE(1,0,STEP_OUTPUTS_PATCHED)
END FUNCTION FDS_HOOK_STEP_OUTPUTS

!> Set T, DT and ICYC for the next FDS_SETUP(MODE=3) call.
SUBROUTINE FDS_HOOK_SET_STEP(T,DT,ICYC) BIND(C,NAME='fds_hook_set_step')
REAL(C_DOUBLE), VALUE :: T,DT
INTEGER(C_INT), VALUE :: ICYC
OUT_T = T ; OUT_DT = DT ; OUT_ICYC = ICYC
END SUBROUTINE FDS_HOOK_SET_STEP

!> Switch the per-mesh dump loop of FDS_SETUP(MODE=3) on (FLAG/=0) or off (patch 0011).
SUBROUTINE FDS_HOOK_SET_MESH_DUMPS(FLAG) BIND(C,NAME='fds_hook_set_mesh_dumps')
INTEGER(C_INT), VALUE :: FLAG
OUT_MESH_DUMPS = (FLAG/=0)
END SUBROUTINE FDS_HOOK_SET_MESH_DUMPS

!> Set the number of output meshes N_OUT_MESHES (patch 0012); effective only if called before fds_setup(0).
SUBROUTINE FDS_HOOK_SET_OUT_MESHES(N) BIND(C,NAME='fds_hook_set_out_meshes')
INTEGER(C_INT), VALUE :: N
N_OUT_MESHES = N
END SUBROUTINE FDS_HOOK_SET_OUT_MESHES

!> Make BOX_VIEW(NM)%<WHICH> a view of the C array P: rank 3 (LB(1:3), EXT(1:3)) or, for ZZ (WHICH=21) and ZZS (WHICH=22), rank 4 (LB(1:4), EXT(1:4)).
!> WHICH: 1 U 2 V 3 W 4 US 5 VS 6 WS 7 D 8 DS 9 H 10 HS 11 KRES 12 FVX 13 FVY 14 FVZ 15 RHO 16 RHOS 17 MU 18 TMP 19 Q 20 RSUM 23 DDDT (21 ZZ, 22 ZZS). NM: 1-based box number;
!> BOX_VIEW is allocated on first use with NMAX entries. Returns 0, or 1 for an unknown WHICH or NM out of range.
FUNCTION FDS_HOOK_SET_VIEW(NMAX,NM,WHICH,LB,EXT,P) BIND(C,NAME='fds_hook_set_view') RESULT(IERR)
INTEGER(C_INT), VALUE :: NMAX,NM,WHICH
INTEGER(C_INT), INTENT(IN) :: LB(4),EXT(4)
TYPE(C_PTR), VALUE :: P
INTEGER(C_INT) :: IERR
IERR = 1
IF (.NOT.ALLOCATED(BOX_VIEW)) ALLOCATE(BOX_VIEW(NMAX))
IF (NM<1 .OR. NM>SIZE(BOX_VIEW)) RETURN
IERR = FDS_HOOK_BIND_VIEW(BOX_VIEW(NM),WHICH,LB,EXT,P)
END FUNCTION FDS_HOOK_SET_VIEW

!> The body of FDS_HOOK_SET_VIEW for any view object V: BOX_VIEW(NM) of a level-0 box, FINE_LEVEL(L)%VIEW(IB) of a fine-level box (draft/fds_fine_box_b.f90). 0, or 1 for an unknown WHICH.
FUNCTION FDS_HOOK_BIND_VIEW(V,WHICH,LB,EXT,P) RESULT(IERR)
TYPE(BOX_VIEW_TYPE), INTENT(INOUT) :: V
INTEGER(C_INT), VALUE :: WHICH
INTEGER(C_INT), INTENT(IN) :: LB(4),EXT(4)
TYPE(C_PTR), VALUE :: P
INTEGER(C_INT) :: IERR
REAL(EB), POINTER, DIMENSION(:) :: FLAT
IERR = 1
IF (WHICH<1 .OR. WHICH>23) RETURN
IF (WHICH==21 .OR. WHICH==22) THEN
   CALL C_F_POINTER(P,FLAT,[INT(EXT(1))*EXT(2)*EXT(3)*EXT(4)])
   IF (WHICH==21) V%ZZ (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1,LB(4):LB(4)+EXT(4)-1) => FLAT
   IF (WHICH==22) V%ZZS(LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1,LB(4):LB(4)+EXT(4)-1) => FLAT
   IERR = 0
   RETURN
ENDIF
CALL C_F_POINTER(P,FLAT,[INT(EXT(1))*EXT(2)*EXT(3)])
SELECT CASE(WHICH)
   CASE( 1) ; V%U   (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 2) ; V%V   (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 3) ; V%W   (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 4) ; V%US  (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 5) ; V%VS  (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 6) ; V%WS  (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 7) ; V%D   (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 8) ; V%DS  (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE( 9) ; V%H   (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(10) ; V%HS  (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(11) ; V%KRES(LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(12) ; V%FVX (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(13) ; V%FVY (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(14) ; V%FVZ (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(15) ; V%RHO (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(16) ; V%RHOS(LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(17) ; V%MU  (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(18) ; V%TMP (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(19) ; V%Q   (LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(20) ; V%RSUM(LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
   CASE(23) ; V%DDDT(LB(1):LB(1)+EXT(1)-1,LB(2):LB(2)+EXT(2)-1,LB(3):LB(3)+EXT(3)-1) => FLAT
END SELECT
IERR = 0
END FUNCTION FDS_HOOK_BIND_VIEW

!> Declare (FLAG/=0) that the kernel wrappers may run on fine-level mesh numbers (only after patch 0007 is validated and the kernels call POINT_TO_BOX; off by default).
SUBROUTINE FDS_HOOK_SET_FINE_READY(FLAG) BIND(C,NAME='fds_hook_set_fine_ready')

INTEGER(C_INT), VALUE :: FLAG

FINE_READY = (FLAG/=0)

END SUBROUTINE FDS_HOOK_SET_FINE_READY

!> Stop the run with a clear message: mesh number NM is not a level-0 FDS mesh and WHERE would read MESHES(NM) (D-056 option B).
SUBROUTINE FDS_HOOK_FINE_ABORT(WHERE,NM)

CHARACTER(*), INTENT(IN) :: WHERE
INTEGER, INTENT(IN) :: NM

WRITE(0,'(A,A,A,I0,A)') 'fds_amr ERROR in ',TRIM(WHERE),': mesh number ',NM,' is not a level-0 FDS mesh. Boxes of refinement level > 0 are fine-level mesh objects (D-056 option B, ' // &
   'FINE_LEVEL(:), patch 0007) and must be reached through POINT_TO_BOX; this routine reads MESHES(NM) and has not been made fine-ready.'
ERROR STOP 1

END SUBROUTINE FDS_HOOK_FINE_ABORT

!> Guard at the entry of every kernel wrapper: a fine-level mesh number (NM > NMESHES) aborts with a clear message unless fine boxes have been declared ready (FINE_READY).
SUBROUTINE FDS_HOOK_FINE_GUARD(WHERE,NM,NMESHES_L0)

CHARACTER(*), INTENT(IN) :: WHERE
INTEGER(C_INT), INTENT(IN) :: NM
INTEGER, INTENT(IN) :: NMESHES_L0

IF (NM>NMESHES_L0 .AND. .NOT.FINE_READY) CALL FDS_HOOK_FINE_ABORT(WHERE,INT(NM))

END SUBROUTINE FDS_HOOK_FINE_GUARD

!> Guard at the entry of a driver routine that is level-0 only (it reads MESHES(NM), FDS pressure/vent/wall structures of a level-0 mesh): a number above NMESHES always aborts, fine-ready or not.
SUBROUTINE FDS_HOOK_L0_ONLY(WHERE,NM,NMESHES_L0)

CHARACTER(*), INTENT(IN) :: WHERE
INTEGER(C_INT), INTENT(IN) :: NM
INTEGER, INTENT(IN) :: NMESHES_L0

IF (NM>NMESHES_L0) CALL FDS_HOOK_FINE_ABORT(WHERE,INT(NM))

END SUBROUTINE FDS_HOOK_L0_ONLY

!> Call of the registered shadow procedure (no-op when none is registered).
SUBROUTINE FDS_HOOK_SHADOW(PHASE,KIND,NM,T,DT,EST,DTNEW,ICHG)

INTEGER, INTENT(IN) :: PHASE,KIND,NM,EST,ICHG
REAL(EB), INTENT(IN) :: T,DT,DTNEW

IF (ASSOCIATED(SHADOW_PROC)) CALL SHADOW_PROC(PHASE,KIND,NM,T,DT,EST,DTNEW,ICHG)

END SUBROUTINE FDS_HOOK_SHADOW

SUBROUTINE FDS_HOOK_LAG_SAVE(NM,EST,NMAX,U,V,W)
INTEGER, INTENT(IN) :: NM,EST,NMAX
REAL(EB), POINTER, INTENT(IN) :: U(:,:,:),V(:,:,:),W(:,:,:)
VERIFY_STOCK_EDGE_LAG = .TRUE.
IF (.NOT.ALLOCATED(LAG)) ALLOCATE(LAG(0:1,NMAX))
LAG(EST,NM)%U = U ; LAG(EST,NM)%V = V ; LAG(EST,NM)%W = W
END SUBROUTINE FDS_HOOK_LAG_SAVE

!> Value of component IC (1 U, 2 V, 3 W) at (I,J,K) that stock's edge loop sees at read key KEYR: the saved (previous-call) value when that ghost is written later.
REAL(EB) FUNCTION STOCK_LAG_VAL(NM,E,KEYR,IC,I,J,K,CUR) RESULT(R)
INTEGER, INTENT(IN) :: NM,E,KEYR,IC,I,J,K
REAL(EB), INTENT(IN) :: CUR
R = CUR
IF (.NOT.ALLOCATED(LAG)) RETURN
IF (.NOT.ALLOCATED(LAG(E,NM)%U) .OR. .NOT.ALLOCATED(LAG(0,NM)%KU)) RETURN
SELECT CASE(IC)
   CASE(1) ; IF (LAG(0,NM)%KU(I,J,K)>KEYR .AND. LAG(0,NM)%KU(I,J,K)<HUGE(1)) R = LAG(E,NM)%U(I,J,K)
   CASE(2) ; IF (LAG(0,NM)%KV(I,J,K)>KEYR .AND. LAG(0,NM)%KV(I,J,K)<HUGE(1)) R = LAG(E,NM)%V(I,J,K)
   CASE(3) ; IF (LAG(0,NM)%KW(I,J,K)>KEYR .AND. LAG(0,NM)%KW(I,J,K)<HUGE(1)) R = LAG(E,NM)%W(I,J,K)
END SELECT
END FUNCTION STOCK_LAG_VAL

!> VEL_T of the ground-edge branch (read key KEYR) at gas cell (I,J,K).
REAL(EB) FUNCTION STOCK_EDGE_LAG_VT(NM,EST,KEYR,I,J,K,UC,VC) RESULT(R)
INTEGER, INTENT(IN) :: NM,KEYR,I,J,K
LOGICAL, INTENT(IN) :: EST
REAL(EB), INTENT(IN) :: UC,VC
INTEGER :: E
E = MERGE(1,0,EST)
R = SQRT(STOCK_LAG_VAL(NM,E,KEYR,1,I,J,K,UC)**2 + STOCK_LAG_VAL(NM,E,KEYR,2,I,J,K,VC)**2)
END FUNCTION STOCK_EDGE_LAG_VT

!> U_NORM of the OPEN-boundary branch (velo.f90 OPEN_AND_WIND_BC, same expressions) with the reads that stock makes before their refresh replaced by the saved values.
REAL(EB) FUNCTION STOCK_EDGE_LAG_UN(NM,EST,KEYR,IEC,IOR,II,JJ,KK,IBAR,JBAR,KBAR,UU,VV,WW,UNC) RESULT(UN)
INTEGER, INTENT(IN) :: NM,KEYR,IEC,IOR,II,JJ,KK,IBAR,JBAR,KBAR
LOGICAL, INTENT(IN) :: EST
REAL(EB), POINTER, INTENT(IN) :: UU(:,:,:),VV(:,:,:),WW(:,:,:)
REAL(EB), INTENT(IN) :: UNC
INTEGER :: E
E = MERGE(1,0,EST)
UN = UNC
SELECT CASE(IEC)
   CASE(1)
      IF (JJ==0    .AND. IOR== 2) UN = 0.5_EB*(G(2,II,   0,KK) + G(2,II,   0,KK+1))
      IF (JJ==JBAR .AND. IOR==-2) UN = 0.5_EB*(G(2,II,JBAR,KK) + G(2,II,JBAR,KK+1))
      IF (KK==0    .AND. IOR== 3) UN = 0.5_EB*(G(3,II,JJ,0)    + G(3,II,JJ+1,   0))
      IF (KK==KBAR .AND. IOR==-3) UN = 0.5_EB*(G(3,II,JJ,KBAR) + G(3,II,JJ+1,KBAR))
   CASE(2)
      IF (II==0    .AND. IOR== 1) UN = 0.5_EB*(G(1,   0,JJ,KK) + G(1,   0,JJ,KK+1))
      IF (II==IBAR .AND. IOR==-1) UN = 0.5_EB*(G(1,IBAR,JJ,KK) + G(1,IBAR,JJ,KK+1))
      IF (KK==0    .AND. IOR== 3) UN = 0.5_EB*(G(3,II,JJ,   0) + G(3,II+1,JJ,   0))
      IF (KK==KBAR .AND. IOR==-3) UN = 0.5_EB*(G(3,II,JJ,KBAR) + G(3,II+1,JJ,KBAR))
   CASE(3)
      IF (II==0    .AND. IOR== 1) UN = 0.5_EB*(G(1,   0,JJ,KK) + G(1,   0,JJ+1,KK))
      IF (II==IBAR .AND. IOR==-1) UN = 0.5_EB*(G(1,IBAR,JJ,KK) + G(1,IBAR,JJ+1,KK))
      IF (JJ==0    .AND. IOR== 2) UN = 0.5_EB*(G(2,II,   0,KK) + G(2,II+1,   0,KK))
      IF (JJ==JBAR .AND. IOR==-2) UN = 0.5_EB*(G(2,II,JBAR,KK) + G(2,II+1,JBAR,KK))
END SELECT
CONTAINS
REAL(EB) FUNCTION G(IC,I,J,K)
INTEGER, INTENT(IN) :: IC,I,J,K
SELECT CASE(IC)
   CASE(1) ; G = STOCK_LAG_VAL(NM,E,KEYR,1,I,J,K,UU(I,J,K))
   CASE(2) ; G = STOCK_LAG_VAL(NM,E,KEYR,2,I,J,K,VV(I,J,K))
   CASE DEFAULT ; G = STOCK_LAG_VAL(NM,E,KEYR,3,I,J,K,WW(I,J,K))
END SELECT
END FUNCTION G
END FUNCTION STOCK_EDGE_LAG_UN

END MODULE FDS_AMREX_HOOKS
