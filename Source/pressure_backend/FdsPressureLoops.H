// FdsPressureLoops.H: host (CPU) translations of FDS pressure-area loops (Role 2, pressure backend).
//
// Each function is the loop text of FDS pres.f90 (FireX 36975d7; line ranges below) with the same operations in the same order, so the
// results are bitwise equal to FDS for the same inputs. They are NOT generator kernels (the S5 generator sidecars are owned by other
// engineers); they are independent, hand-written CPU references and host code, tagged
//     HAND-WRITTEN <loop id> <file:lines> (FireX 36975d7)
// at the function, and checked bitwise against the verbatim upstream lines compiled with gfortran (tests/fds_loops_test.py,
// frozen/fds-loops-notes.md). Claim register: docs/amrex/loop_claims.csv, owner AMR Pressure Backend Implementer.
//
//   L1211  pres_compute_rhs_div        pres.f90:250-260   PRESSURE_SOLVER_COMPUTE_RHS, PRHS, 3-D branch for IPS 1/4/7
//   L1207  pres_p_from_h               pres.f90:758-764   PRESSURE_SOLVER_CHECK_RESIDUALS, P = RHOP*(HP-KRES) over the full box
//   L1220  pres_h_bc_x                 pres.f90:450-462   PRESSURE_SOLVER_FFT, H boundary fill in x (LBC codes)
//   L1221  pres_h_bc_y                 pres.f90:466-477   same in y (MBC codes)
//   L1222  pres_h_bc_z                 pres.f90:481-492   same in z (NBC codes)
//   L1209  pres_poisson_boundary_arrays pres.f90:65-228   PRESSURE_SOLVER_COMPUTE_RHS, BXS..BZF from the external wall cells
//   L1215  (retired) pres.f90:394-400: the copy of PRHS into HP after the FFT solve has no counterpart, pressure_backend writes H.
//
// Arrays are Fortran-ordered views (I fastest) with explicit lower bounds and extents, so an FDS array (0:IBP1,...) or an AMReX
// Array4 of one box is described without copying. Index conventions are those of FDS: cell (I,J,K) with I=1..IBAR, ghost cells 0 and
// IBP1. A Box3 selects the cell region; FDS itself is Box3{1,1,1, IBAR,JBAR,KBAR}. Nothing here depends on AMReX.
#ifndef PB_FDS_PRESSURE_LOOPS_H
#define PB_FDS_PRESSURE_LOOPS_H

#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>

namespace pressure_backend { namespace fdsloops {

struct F1 {   // Fortran rank-1 view: element i lives at p[i - lo]
    double* p = nullptr; int lo = 0; int n = 0;
    double& operator()(int i) const { return p[i - lo]; }
    bool covers(int a, int b) const { return p && a >= lo && b <= lo + n - 1; }
};
struct F2 {   // rank-2 view, first index fastest
    double* p = nullptr; int lo[2] = {1, 1}; int n[2] = {0, 0};
    double& operator()(int i, int j) const { return p[(i - lo[0]) + static_cast<std::ptrdiff_t>(n[0]) * (j - lo[1])]; }
    bool covers(int a0, int b0, int a1, int b1) const { return p && a0 >= lo[0] && b0 <= lo[0] + n[0] - 1 && a1 >= lo[1] && b1 <= lo[1] + n[1] - 1; }
};
struct F3 {   // rank-3 view, first index fastest
    double* p = nullptr; int lo[3] = {0, 0, 0}; int n[3] = {0, 0, 0};
    double& operator()(int i, int j, int k) const {
        return p[(i - lo[0]) + static_cast<std::ptrdiff_t>(n[0]) * ((j - lo[1]) + static_cast<std::ptrdiff_t>(n[1]) * (k - lo[2]))];
    }
    bool covers(int a0, int b0, int a1, int b1, int a2, int b2) const {
        return p && a0 >= lo[0] && b0 <= lo[0] + n[0] - 1 && a1 >= lo[1] && b1 <= lo[1] + n[1] - 1 && a2 >= lo[2] && b2 <= lo[2] + n[2] - 1;
    }
    // View of an amrex::Array4-like object (members p, jstride, kstride, begin.x/y/z, end.x/y/z; unit i stride, no component offset).
    // Only valid when jstride == nx and kstride == nx*ny, which holds for a FArrayBox component.
    template <class A4> static F3 from_array4(A4 const& a) {
        F3 f; f.p = const_cast<double*>(a.p);
        f.lo[0] = a.begin.x; f.lo[1] = a.begin.y; f.lo[2] = a.begin.z;
        f.n[0] = a.end.x - a.begin.x; f.n[1] = a.end.y - a.begin.y; f.n[2] = a.end.z - a.begin.z;
        return f;
    }
};

struct Box3 { int lo[3]; int hi[3]; };   // inclusive cell range

// Thrown for a feature the translation does not cover (the product route reports it as NotBuilt, never silently as another feature).
struct NotBuilt : std::runtime_error { using std::runtime_error::runtime_error; };

// What the H boundary fill (L1220-L1222) needs to know about the whole mesh besides the box it fills.
struct HFillOptions {
    Box3 domain;                              // the whole cell extent of the mesh in all three directions (FDS: fds_interior)
    bool tunnel_preconditioner = false;       // FDS TUNNEL_PRECONDITIONER: its add-back (pres.f90:436-445) is not translated -> NotBuilt
};

inline Box3 fds_interior (int ibar, int jbar, int kbar) { return Box3{{1, 1, 1}, {ibar, jbar, kbar}}; }
inline Box3 fds_with_ghosts (int ibar, int jbar, int kbar) { return Box3{{0, 0, 0}, {ibar + 1, jbar + 1, kbar + 1}}; }

// ---- L1211 ----------------------------------------------------------------------------------------------------------------------
// HAND-WRITTEN L1211 pres.f90:250-260 (FireX 36975d7), non-cylindrical branch for IPS 1/4/7:
//   TRM1=(FVX(I-1,J,K)-FVX(I,J,K))*RDX(I); TRM2=(FVY(I,J-1,K)-FVY(I,J,K))*RDY(J); TRM3=(FVZ(I,J,K-1)-FVZ(I,J,K))*RDZ(K);
//   TRM4=-DDDT(I,J,K); PRHS(I,J,K)=TRM1+TRM2+TRM3+TRM4   over the cells of b.
// prhs is indexed (I,J,K) with its own extents (FDS: ITRN,JTRN,KTRN, which can be IBAR when x is periodic). Only IPS 1/4/7 (natural
// layout); IPS 2,3,5,6 store PRHS transposed (L1212-L1214) and are not provided: the AMR route keeps the natural layout.
// Throws std::invalid_argument if a view does not cover what the loop touches.
void pres_compute_rhs_div (Box3 const& b, F3 const& fvx, F3 const& fvy, F3 const& fvz, F3 const& dddt,
                           F1 const& rdx, F1 const& rdy, F1 const& rdz, F3 const& prhs);

// ---- L1207 ----------------------------------------------------------------------------------------------------------------------
// HAND-WRITTEN L1207 pres.f90:758-764 (FireX 36975d7): P(I,J,K)=RHOP(I,J,K)*(HP(I,J,K)-KRES(I,J,K)) over the cells of b.
// FDS passes the box with ghost cells, K=0..KBP1, J=0..JBP1, I=0..IBP1 (fds_with_ghosts): the ghost values of HP (after the H
// boundary fill, L1220-L1222), RHOP and KRES must be set. RHOP is allocated -1:IBP1+1 in FDS; any view that covers b works.
void pres_p_from_h (Box3 const& b, F3 const& rhop, F3 const& hp, F3 const& kres, F3 const& p);

// ---- L1220, L1221, L1222 ---------------------------------------------------------------------------------------------------------
// HAND-WRITTEN L1220 pres.f90:450-462, L1221 pres.f90:466-477, L1222 pres.f90:481-492 (FireX 36975d7): ghost layer of H next to
// the faces of b from the Poisson boundary arrays and the boundary codes (FDS: b = fds_interior, ghost cells 0 and IBAR+1, ...).
// Codes (FISHPAK/Crayfishpak): 0 periodic, 1 Dirichlet-Dirichlet, 2 Dirichlet low / Neumann high, 3 Neumann-Neumann,
// 4 Neumann low / Dirichlet high; the x direction (LBC) also has 5 (Dirichlet high, H(0)=H(1)) and 6 (Neumann high, H(0)=H(1)).
// bxs/bxf are indexed (J,K) over the J,K range of b, bys/byf (I,K), bzs/bzf (I,J). dxi, deta, dzeta are the Poisson coordinate
// spacings (FDS DXI, DETA, DZETA). Periodic (code 0) wraps inside b, so b must be the whole domain extent in that direction:
// opt.domain gives that extent and the function throws std::invalid_argument if code 0 is asked for with b.lo/b.hi different from
// opt.domain.lo/hi in the direction (a per-box call on a decomposed level would silently wrap inside the box). With
// opt.tunnel_preconditioner the functions throw NotBuilt (the add-back of H_BAR, BXS_BAR, BXF_BAR between the solve and the fill is
// not translated). The statements are kept in the FDS order (later ones overwrite earlier ones where several apply).
void pres_h_bc_x (Box3 const& b, HFillOptions const& opt, int lbc, double dxi, F2 const& bxs, F2 const& bxf, F3 const& hp);
void pres_h_bc_y (Box3 const& b, HFillOptions const& opt, int mbc, double deta, F2 const& bys, F2 const& byf, F3 const& hp);
void pres_h_bc_z (Box3 const& b, HFillOptions const& opt, int nbc, double dzeta, F2 const& bzs, F2 const& bzf, F3 const& hp);

// ---- L1209 ----------------------------------------------------------------------------------------------------------------------
// Constants of FDS cons.f90 used by the loop.
namespace fdsconst {
constexpr int NULL_BOUNDARY = 0, SOLID_BOUNDARY = 1, OPEN_BOUNDARY = 2, INTERPOLATED_BOUNDARY = 6;
constexpr int DIRICHLET = 1, NEUMANN = 2;
constexpr double TWENTY_EPSILON_EB = 20.0 * 2.220446049250313e-16;   // 20*EPSILON(1._EB)
}

struct PoissonVent {      // the VENTS(WC%VENT_INDEX) fields the loop reads (OPEN boundaries)
    int pressure_ramp_index = 0;
    double dynamic_pressure = 0.0;
    int n_eddy = 0;       // VT%N_EDDY
    int ior = 0;          // VT%IOR (the vent's, not the wall cell's)
    F2 u_eddy, v_eddy, w_eddy;   // VT%U_EDDY(J,K), V_EDDY(I,K), W_EDDY(I,J), read only if n_eddy > 0
};

struct PoissonWall {      // one external wall cell (FDS WALL(IW), EXTERNAL_WALL(IW), BOUNDARY_COORD, BOUNDARY_PROP1)
    int i = 0, j = 0, k = 0;          // BC%II, JJ, KK
    int ior = 0;                      // BC%IOR
    int pressure_bc_type = 0;         // EWC%PRESSURE_BC_TYPE (fdsconst::NEUMANN, DIRICHLET)
    int boundary_type = 0;            // WC%BOUNDARY_TYPE
    double dundt = 0.0;               // EWC%DUNDT
    double wall_work1 = 0.0;          // WALL_WORK1(IW)
    // INTERPOLATED boundary: the neighbour mesh cell width in the normal direction at the overlap, i.e. MESHES(NOM)%DX(EWC%IIO_MIN)
    // for |IOR|=1, %DY(JJO_MIN) for 2, %DZ(KKO_MIN) for 3. The caller takes it from the driver's ghost-exchange data (no neighbour read here).
    double other_d = 0.0;
    // OPEN boundary
    double t_ign = 0.0;               // B1%T_IGN
    double rho_f = 1.0;               // B1%RHO_F
    PoissonVent const* vent = nullptr;   // VT; required for OPEN walls
};

struct PoissonBcContext {  // the mesh data of pres.f90:65-228 (UU,VV,WW,HP are the predictor or corrector set chosen by the caller)
    int ibar = 0, jbar = 0, kbar = 0;
    F3 fvx, fvy, fvz, hp, kres, uu, vv, ww;
    F1 hx, hy, hz;                    // HX(0:IBP1), HY(0:JBP1), HZ(0:KBP1)
    F1 dx, dy, dz, rdxn, rdyn, rdzn;  // DX(1:IBAR).., RDXN(0:IBAR).. (RDXN(0), RDXN(IBAR)), ...
    F1 u_wind, v_wind, w_wind;        // U_WIND(0:KBP1) etc.: read at the wall cell's KK (0 and KBP1 for z walls) if open_wind_boundary; must cover 0..KBP1
    bool tunnel_preconditioner = false;   // FDS TUNNEL_PRECONDITIONER: not translated, NotBuilt
    double u0 = 0.0, v0 = 0.0, w0 = 0.0;
    bool open_wind_boundary = false;
    double t = 0.0, dt = 0.0, t_begin = 0.0;
    // EVALUATE_RAMP(TSI,RAMP_INDEX) of FDS math.f90, supplied by the caller (callee of the loop; not translated here).
    std::function<double(double, int)> evaluate_ramp;
    F2 bxs, bxf, bys, byf, bzs, bzf;  // BXS(J,K),BXF; BYS(I,K),BYF; BZS(I,J),BZF (outputs; entries not set by a wall keep their value)
};

// HAND-WRITTEN L1209 pres.f90:65-228 (FireX 36975d7): Poisson boundary arrays BXS..BZF from the external wall cells, in the order of
// the wall list (a later wall writes over an earlier one that names the same entry, as in FDS). Neumann: dH/dn = -F_n - dun/dt;
// Dirichlet: average of the ghost and gas H plus WALL_WORK1 (solid), the interpolated average with the neighbour width, or the open
// boundary value P_EXTERNAL/RHO_F + KRES or + H0 (with the wind option). The NIC>1 / neighbour branches are not reached in the AMR
// route (D-072): the neighbour width comes in through PoissonWall::other_d.
void pres_poisson_boundary_arrays (PoissonBcContext const& c, PoissonWall const* walls, int nwalls);

}}   // namespace pressure_backend::fdsloops

#endif
