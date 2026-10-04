// FdsPressureLoops.cpp: see FdsPressureLoops.H. Hand-written host translations of FDS pres.f90 loops (FireX 36975d7).
// Keep every expression in the operand order of the Fortran text: the bitwise tests compare with the verbatim loops.
// Build flags: -ffp-contract=off (set in harness/CMakeLists.txt), no fast-math.
//
// PB_FDSLOOPS_MUTANT=n (test builds only, tests/m5_fds_loops.cmake) changes one operand, sign or index so the bitwise test
// can show that it notices (mutation check of generator-howto.md section 4); 0 / undefined is the real code.
#include "FdsPressureLoops.H"

#ifndef PB_FDSLOOPS_MUTANT
#define PB_FDSLOOPS_MUTANT 0
#endif
#define MUT(n) (PB_FDSLOOPS_MUTANT == (n))

namespace pressure_backend { namespace fdsloops {

namespace {
void need (bool ok, const char* what) { if (!ok) throw std::invalid_argument(std::string("fdsloops: view does not cover ") + what); }
void check_hfill (Box3 const& b, HFillOptions const& opt, int dir, int code, const char* name)
{
#if PB_FDSLOOPS_MUTANT != 21
    if (opt.tunnel_preconditioner) throw NotBuilt("fdsloops: TUNNEL_PRECONDITIONER is not built (the H_BAR/BXS_BAR/BXF_BAR add-back before the H boundary fill is not translated)");
#endif
#if PB_FDSLOOPS_MUTANT != 19
    if (code == 0 && (b.lo[dir] != opt.domain.lo[dir] || b.hi[dir] != opt.domain.hi[dir]))
        throw std::invalid_argument(std::string("fdsloops: periodic wrap in ") + name + " needs the box to be the whole domain extent in that direction");
#endif
}
}

// ---- L1211 ----------------------------------------------------------------------------------------------------------------------
// HAND-WRITTEN L1211 pres.f90:250-260 (FireX 36975d7)
void pres_compute_rhs_div (Box3 const& b, F3 const& fvx, F3 const& fvy, F3 const& fvz, F3 const& dddt,
                           F1 const& rdx, F1 const& rdy, F1 const& rdz, F3 const& prhs)
{
    const int il = b.lo[0], ih = b.hi[0], jl = b.lo[1], jh = b.hi[1], kl = b.lo[2], kh = b.hi[2];
    need(fvx.covers(il - 1, ih, jl, jh, kl, kh), "FVX(I-1:I,J,K)");
    need(fvy.covers(il, ih, jl - 1, jh, kl, kh), "FVY(I,J-1:J,K)");
    need(fvz.covers(il, ih, jl, jh, kl - 1, kh), "FVZ(I,J,K-1:K)");
    need(dddt.covers(il, ih, jl, jh, kl, kh), "DDDT");
    need(rdx.covers(il, ih) && rdy.covers(jl, jh) && rdz.covers(kl, kh), "RDX/RDY/RDZ");
    need(prhs.covers(il, ih, jl, jh, kl, kh), "PRHS");
    for (int K = kl; K <= kh; ++K) {
        for (int J = jl; J <= jh; ++J) {
            for (int I = il; I <= ih; ++I) {
                double TRM1, TRM2, TRM3, TRM4;
#if PB_FDSLOOPS_MUTANT == 1
                TRM1 = (fvx(I, J, K) - fvx(I - 1, J, K)) * rdx(I);   // sign of the x difference
#else
                TRM1 = (fvx(I - 1, J, K) - fvx(I, J, K)) * rdx(I);
#endif
                TRM2 = (fvy(I, J - 1, K) - fvy(I, J, K)) * rdy(J);
                TRM3 = (fvz(I, J, K - 1) - fvz(I, J, K)) * rdz(K);
                TRM4 = -dddt(I, J, K);
#if PB_FDSLOOPS_MUTANT == 2
                prhs(I, J, K) = TRM1 + TRM2 + (TRM3 + TRM4);   // different association
#elif PB_FDSLOOPS_MUTANT == 3
                prhs(I, J, K) = TRM1 + TRM2 + TRM3 - TRM4;     // sign of DDDT
#elif PB_FDSLOOPS_MUTANT == 4
                prhs(I, J, K) = TRM1 + TRM3 + TRM2 + TRM4;     // y and z terms swapped in the sum
#else
                prhs(I, J, K) = TRM1 + TRM2 + TRM3 + TRM4;
#endif
            }
        }
    }
}

// ---- L1207 ----------------------------------------------------------------------------------------------------------------------
// HAND-WRITTEN L1207 pres.f90:758-764 (FireX 36975d7)
void pres_p_from_h (Box3 const& b, F3 const& rhop, F3 const& hp, F3 const& kres, F3 const& p)
{
    need(rhop.covers(b.lo[0], b.hi[0], b.lo[1], b.hi[1], b.lo[2], b.hi[2]), "RHOP");
    need(hp.covers(b.lo[0], b.hi[0], b.lo[1], b.hi[1], b.lo[2], b.hi[2]), "HP");
    need(kres.covers(b.lo[0], b.hi[0], b.lo[1], b.hi[1], b.lo[2], b.hi[2]), "KRES");
    need(p.covers(b.lo[0], b.hi[0], b.lo[1], b.hi[1], b.lo[2], b.hi[2]), "P");
#if PB_FDSLOOPS_MUTANT == 5
    const int il = b.lo[0] + 1, ih = b.hi[0] - 1, jl = b.lo[1] + 1, jh = b.hi[1] - 1, kl = b.lo[2] + 1, kh = b.hi[2] - 1;   // interior only (ghost cells skipped)
#else
    const int il = b.lo[0], ih = b.hi[0], jl = b.lo[1], jh = b.hi[1], kl = b.lo[2], kh = b.hi[2];
#endif
    for (int K = kl; K <= kh; ++K)
        for (int J = jl; J <= jh; ++J)
            for (int I = il; I <= ih; ++I) {
#if PB_FDSLOOPS_MUTANT == 6
                p(I, J, K) = rhop(I, J, K) * (hp(I, J, K) + kres(I, J, K));   // sign of KRES
#elif PB_FDSLOOPS_MUTANT == 7
                p(I, J, K) = rhop(I, J, K) * hp(I, J, K) - rhop(I, J, K) * kres(I, J, K);   // distributed
#else
                p(I, J, K) = rhop(I, J, K) * (hp(I, J, K) - kres(I, J, K));
#endif
            }
}

// ---- L1220, L1221, L1222 ---------------------------------------------------------------------------------------------------------
// HAND-WRITTEN L1220 pres.f90:450-462 (FireX 36975d7)
void pres_h_bc_x (Box3 const& b, HFillOptions const& opt, int LBC, double DXI, F2 const& bxs, F2 const& bxf, F3 const& hp)
{
    check_hfill(b, opt, 0, LBC, "x");
    const int i0 = b.lo[0], i1 = b.hi[0], jl = b.lo[1], jh = b.hi[1], kl = b.lo[2], kh = b.hi[2];   // FDS: i0=1, i1=IBAR
    need(hp.covers(i0 - 1, i1 + 1, jl, jh, kl, kh), "HP(ILO-1:IHI+1,J,K)");
    need(bxs.covers(jl, jh, kl, kh) && bxf.covers(jl, jh, kl, kh), "BXS/BXF(J,K)");
    for (int K = kl; K <= kh; ++K) {
        for (int J = jl; J <= jh; ++J) {
            if (LBC == 3 || LBC == 4)             hp(i0 - 1, J, K) = hp(i0, J, K) - DXI * bxs(J, K);
#if PB_FDSLOOPS_MUTANT == 8
            if (LBC == 3 || LBC == 2 || LBC == 6) hp(i1 + 1, J, K) = hp(i1, J, K) - DXI * bxf(J, K);   // sign at the high face
#else
            if (LBC == 3 || LBC == 2 || LBC == 6) hp(i1 + 1, J, K) = hp(i1, J, K) + DXI * bxf(J, K);
#endif
            if (LBC == 1 || LBC == 2)             hp(i0 - 1, J, K) = -hp(i0, J, K) + 2.0 * bxs(J, K);
#if PB_FDSLOOPS_MUTANT == 9
            if (LBC == 1 || LBC == 4)             hp(i1 + 1, J, K) = -hp(i1, J, K) + 2.0 * bxf(J, K);   // code 5 dropped
#else
            if (LBC == 1 || LBC == 4 || LBC == 5) hp(i1 + 1, J, K) = -hp(i1, J, K) + 2.0 * bxf(J, K);
#endif
            if (LBC == 5 || LBC == 6)             hp(i0 - 1, J, K) = hp(i0, J, K);
            if (LBC == 0) {
                hp(i0 - 1, J, K) = hp(i1, J, K);
                hp(i1 + 1, J, K) = hp(i0, J, K);
            }
        }
    }
}

// HAND-WRITTEN L1221 pres.f90:466-477 (FireX 36975d7)
void pres_h_bc_y (Box3 const& b, HFillOptions const& opt, int MBC, double DETA, F2 const& bys, F2 const& byf, F3 const& hp)
{
    check_hfill(b, opt, 1, MBC, "y");
    const int il = b.lo[0], ih = b.hi[0], j0 = b.lo[1], j1 = b.hi[1], kl = b.lo[2], kh = b.hi[2];   // FDS: j0=1, j1=JBAR
    need(hp.covers(il, ih, j0 - 1, j1 + 1, kl, kh), "HP(I,JLO-1:JHI+1,K)");
    need(bys.covers(il, ih, kl, kh) && byf.covers(il, ih, kl, kh), "BYS/BYF(I,K)");
    for (int K = kl; K <= kh; ++K) {
        for (int I = il; I <= ih; ++I) {
            if (MBC == 3 || MBC == 4) hp(I, j0 - 1, K) = hp(I, j0, K) - DETA * bys(I, K);
            if (MBC == 3 || MBC == 2) hp(I, j1 + 1, K) = hp(I, j1, K) + DETA * byf(I, K);
            if (MBC == 1 || MBC == 2) hp(I, j0 - 1, K) = -hp(I, j0, K) + 2.0 * bys(I, K);
#if PB_FDSLOOPS_MUTANT == 10
            if (MBC == 1 || MBC == 4) hp(I, j1 + 1, K) = -hp(I, j1, K) + 2.0 * bys(I, K);   // wrong array (BYS for BYF)
#else
            if (MBC == 1 || MBC == 4) hp(I, j1 + 1, K) = -hp(I, j1, K) + 2.0 * byf(I, K);
#endif
            if (MBC == 0) {
                hp(I, j0 - 1, K) = hp(I, j1, K);
                hp(I, j1 + 1, K) = hp(I, j0, K);
            }
        }
    }
}

// HAND-WRITTEN L1222 pres.f90:481-492 (FireX 36975d7)
void pres_h_bc_z (Box3 const& b, HFillOptions const& opt, int NBC, double DZETA, F2 const& bzs, F2 const& bzf, F3 const& hp)
{
    check_hfill(b, opt, 2, NBC, "z");
    const int il = b.lo[0], ih = b.hi[0], jl = b.lo[1], jh = b.hi[1], k0 = b.lo[2], k1 = b.hi[2];   // FDS: k0=1, k1=KBAR
    need(hp.covers(il, ih, jl, jh, k0 - 1, k1 + 1), "HP(I,J,KLO-1:KHI+1)");
    need(bzs.covers(il, ih, jl, jh) && bzf.covers(il, ih, jl, jh), "BZS/BZF(I,J)");
    for (int J = jl; J <= jh; ++J) {
        for (int I = il; I <= ih; ++I) {
            if (NBC == 3 || NBC == 4) hp(I, J, k0 - 1) = hp(I, J, k0) - DZETA * bzs(I, J);
#if PB_FDSLOOPS_MUTANT == 11
            if (NBC == 3 || NBC == 2) hp(I, J, k1 + 1) = hp(I, J, k1) + DZETA * bzs(I, J);   // wrong array (BZS for BZF)
#else
            if (NBC == 3 || NBC == 2) hp(I, J, k1 + 1) = hp(I, J, k1) + DZETA * bzf(I, J);
#endif
            if (NBC == 1 || NBC == 2) hp(I, J, k0 - 1) = -hp(I, J, k0) + 2.0 * bzs(I, J);
            if (NBC == 1 || NBC == 4) hp(I, J, k1 + 1) = -hp(I, J, k1) + 2.0 * bzf(I, J);
            if (NBC == 0) {
                hp(I, J, k0 - 1) = hp(I, J, k1);
                hp(I, J, k1 + 1) = hp(I, J, k0);
            }
        }
    }
}

// ---- L1209 ----------------------------------------------------------------------------------------------------------------------
// HAND-WRITTEN L1209 pres.f90:65-228 (FireX 36975d7)
void pres_poisson_boundary_arrays (PoissonBcContext const& c, PoissonWall const* walls, int nwalls)
{
    using namespace fdsconst;
#if PB_FDSLOOPS_MUTANT != 21
    if (c.tunnel_preconditioner) throw NotBuilt("fdsloops: TUNNEL_PRECONDITIONER is not built (the H_BAR/BXS_BAR/BXF_BAR add-back is not translated)");
#endif
    const int IBAR = c.ibar, JBAR = c.jbar, KBAR = c.kbar, IBP1 = IBAR + 1, JBP1 = JBAR + 1, KBP1 = KBAR + 1;
    need(c.hp.covers(0, IBP1, 0, JBP1, 0, KBP1) && c.kres.covers(1, IBAR, 1, JBAR, 1, KBAR), "HP(0:IBP1,..) / KRES");
    need(c.fvx.covers(0, IBAR, 0, JBP1, 0, KBP1) && c.fvy.covers(0, IBP1, 0, JBAR, 0, KBP1) && c.fvz.covers(0, IBP1, 0, JBP1, 0, KBAR), "FVX/FVY/FVZ");
    need(c.hx.covers(0, IBP1) && c.hy.covers(0, JBP1) && c.hz.covers(0, KBP1), "HX/HY/HZ");
    need(c.dx.covers(1, IBAR) && c.dy.covers(1, JBAR) && c.dz.covers(1, KBAR), "DX/DY/DZ");
    need(c.rdxn.covers(0, IBAR) && c.rdyn.covers(0, JBAR) && c.rdzn.covers(0, KBAR), "RDXN/RDYN/RDZN");
    need(c.bxs.covers(1, JBAR, 1, KBAR) && c.bxf.covers(1, JBAR, 1, KBAR) && c.bys.covers(1, IBAR, 1, KBAR) && c.byf.covers(1, IBAR, 1, KBAR) &&
         c.bzs.covers(1, IBAR, 1, JBAR) && c.bzf.covers(1, IBAR, 1, JBAR), "BXS..BZF");

    bool any_open = false;
    for (int IW = 0; IW < nwalls; ++IW) {
        const PoissonWall& w = walls[IW];
        const int a = w.ior < 0 ? -w.ior : w.ior;
        if (a < 1 || a > 3) throw std::invalid_argument("fdsloops: wall IOR must be +-1, +-2 or +-3");
        const bool okJK = w.j >= 1 && w.j <= JBAR && w.k >= 1 && w.k <= KBAR, okIK = w.i >= 1 && w.i <= IBAR && w.k >= 1 && w.k <= KBAR,
                   okIJ = w.i >= 1 && w.i <= IBAR && w.j >= 1 && w.j <= JBAR;
        if (!(a == 1 ? okJK : a == 2 ? okIK : okIJ)) throw std::invalid_argument("fdsloops: wall cell index outside the mesh");
        if (w.boundary_type == OPEN_BOUNDARY && w.pressure_bc_type == DIRICHLET) any_open = true;
    }
    if (any_open) {
        need(c.uu.covers(0, IBAR, 1, JBAR, 1, KBAR) && c.vv.covers(1, IBAR, 0, JBAR, 1, KBAR) && c.ww.covers(1, IBAR, 1, JBAR, 0, KBAR), "UU/VV/WW");
        need(static_cast<bool>(c.evaluate_ramp), "evaluate_ramp (a callable)");
#if PB_FDSLOOPS_MUTANT == 20
        if (c.open_wind_boundary) need(c.u_wind.covers(1, KBAR) && c.v_wind.covers(1, KBAR) && c.w_wind.covers(1, KBAR), "U_WIND/V_WIND/W_WIND(K)");   // the former guard
#else
        // FDS: U_WIND, V_WIND, W_WIND(0:KBP1); the loop reads them at the wall cell's KK, which is 0 or KBP1 for a z wall
        if (c.open_wind_boundary) need(c.u_wind.covers(0, KBP1) && c.v_wind.covers(0, KBP1) && c.w_wind.covers(0, KBP1), "U_WIND/V_WIND/W_WIND(0:KBP1)");
#endif
    }

    const F3& HP = c.hp; const F3& KRES = c.kres; const F3& UU = c.uu; const F3& VV = c.vv; const F3& WW = c.ww;
    const F3& FVX = c.fvx; const F3& FVY = c.fvy; const F3& FVZ = c.fvz;
    const F2& BXS = c.bxs; const F2& BXF = c.bxf; const F2& BYS = c.bys; const F2& BYF = c.byf; const F2& BZS = c.bzs; const F2& BZF = c.bzf;

    for (int IW = 0; IW < nwalls; ++IW) {
        const PoissonWall& WC = walls[IW];
        const int I = WC.i, J = WC.j, K = WC.k, IOR = WC.ior;

        // Apply pressure gradients at NEUMANN boundaries: dH/dn = -F_n - d(u_n)/dt
        if (WC.pressure_bc_type == NEUMANN) {
            switch (IOR) {
                case  1: BXS(J, K) = c.hx(0)    * (-FVX(0, J, K)    + WC.dundt); break;
                case -1: BXF(J, K) = c.hx(IBP1) * (-FVX(IBAR, J, K) - WC.dundt); break;
                case  2: BYS(I, K) = c.hy(0)    * (-FVY(I, 0, K)    + WC.dundt); break;
                case -2: BYF(I, K) = c.hy(JBP1) * (-FVY(I, JBAR, K) - WC.dundt); break;
                case  3: BZS(I, J) = c.hz(0)    * (-FVZ(I, J, 0)    + WC.dundt); break;
                #if PB_FDSLOOPS_MUTANT == 18
                case -3: BZF(I, J) = c.hz(KBP1) * (-FVZ(I, J, KBAR) + WC.dundt); break;
#else
                case -3: BZF(I, J) = c.hz(KBP1) * (-FVZ(I, J, KBAR) - WC.dundt); break;
#endif
                default: break;
            }
        }

        // Apply pressures at DIRICHLET boundaries, depending on the specific type
        if (WC.pressure_bc_type == DIRICHLET) {

            if (WC.boundary_type != OPEN_BOUNDARY && WC.boundary_type != INTERPOLATED_BOUNDARY) {
                // Solid boundary that uses a Dirichlet BC: average of the last computed ghost and adjacent gas cell pressures.
                switch (IOR) {
                    case  1: BXS(J, K) = 0.5 * (HP(0, J, K)    + HP(1, J, K))    + WC.wall_work1; break;
                    #if PB_FDSLOOPS_MUTANT == 12
                    case -1: BXF(J, K) = 0.5 * (HP(IBAR, J, K) + HP(IBP1, J, K)); break;   // WALL_WORK1 dropped
#else
                    case -1: BXF(J, K) = 0.5 * (HP(IBAR, J, K) + HP(IBP1, J, K)) + WC.wall_work1; break;
#endif
                    case  2: BYS(I, K) = 0.5 * (HP(I, 0, K)    + HP(I, 1, K))    + WC.wall_work1; break;
                    case -2: BYF(I, K) = 0.5 * (HP(I, JBAR, K) + HP(I, JBP1, K)) + WC.wall_work1; break;
                    case  3: BZS(I, J) = 0.5 * (HP(I, J, 0)    + HP(I, J, 1))    + WC.wall_work1; break;
                    case -3: BZF(I, J) = 0.5 * (HP(I, J, KBAR) + HP(I, J, KBP1)) + WC.wall_work1; break;
                    default: break;
                }
            }

            // Interpolated boundary: average of the neighbouring cells from the previous time step (HP of the neighbour mesh has
            // already been copied to the external cells in NO_FLUX).
            if (WC.boundary_type == INTERPOLATED_BOUNDARY) {
                const double OD = WC.other_d;   // DX_OTHER, DY_OTHER or DZ_OTHER of the case
                switch (IOR) {
                    case  1: BXS(J, K) = (OD * HP(1, J, K)    + c.dx(1)    * HP(0, J, K))    / (c.dx(1)    + OD) + WC.wall_work1; break;
                    case -1: BXF(J, K) = (OD * HP(IBAR, J, K) + c.dx(IBAR) * HP(IBP1, J, K)) / (c.dx(IBAR) + OD) + WC.wall_work1; break;
                    #if PB_FDSLOOPS_MUTANT == 13
                    case  2: BYS(I, K) = (OD * HP(I, 1, K)    + c.dy(1)    * HP(I, 0, K))    / (c.dy(1)) + WC.wall_work1; break;   // OD missing in the denominator
#else
                    case  2: BYS(I, K) = (OD * HP(I, 1, K)    + c.dy(1)    * HP(I, 0, K))    / (c.dy(1)    + OD) + WC.wall_work1; break;
#endif
                    case -2: BYF(I, K) = (OD * HP(I, JBAR, K) + c.dy(JBAR) * HP(I, JBP1, K)) / (c.dy(JBAR) + OD) + WC.wall_work1; break;
                    case  3: BZS(I, J) = (OD * HP(I, J, 1)    + c.dz(1)    * HP(I, J, 0))    / (c.dz(1)    + OD) + WC.wall_work1; break;
                    case -3: BZF(I, J) = (OD * HP(I, J, KBAR) + c.dz(KBAR) * HP(I, J, KBP1)) / (c.dz(KBAR) + OD) + WC.wall_work1; break;
                    default: break;
                }
            }

            // OPEN (passive opening to exterior of domain) boundary. Apply inflow/outflow BC.
            if (WC.boundary_type == OPEN_BOUNDARY) {
                if (!WC.vent) throw std::invalid_argument("fdsloops: OPEN wall without a PoissonVent");
                const PoissonVent& VT = *WC.vent;
                double TSI;
                const double d_ign = WC.t_ign - c.t_begin;
                #if PB_FDSLOOPS_MUTANT == 16
                if ((d_ign < 0.0 ? -d_ign : d_ign) <= TWENTY_EPSILON_EB) {   // ramp-index condition dropped
#else
                if ((d_ign < 0.0 ? -d_ign : d_ign) <= TWENTY_EPSILON_EB && VT.pressure_ramp_index >= 1) {
#endif
                    TSI = c.t;
                } else {
                    TSI = c.t - c.t_begin;
                }
                const double TIME_RAMP_FACTOR = c.evaluate_ramp(TSI, VT.pressure_ramp_index);
                const double P_EXTERNAL = TIME_RAMP_FACTOR * VT.dynamic_pressure;

                // Synthetic eddy method for OPEN inflow boundaries
                double VEL_EDDY = 0.0;
                if (VT.n_eddy > 0) {
                    #if PB_FDSLOOPS_MUTANT == 14
                    switch (IOR < 0 ? -IOR : IOR) {   // the wall's IOR instead of the vent's
#else
                    switch (VT.ior < 0 ? -VT.ior : VT.ior) {
#endif
                        case 1: need(VT.u_eddy.covers(J, J, K, K), "VT%U_EDDY(J,K)"); VEL_EDDY = VT.u_eddy(J, K); break;
                        case 2: need(VT.v_eddy.covers(I, I, K, K), "VT%V_EDDY(I,K)"); VEL_EDDY = VT.v_eddy(I, K); break;
                        case 3: need(VT.w_eddy.covers(I, I, J, J), "VT%W_EDDY(I,J)"); VEL_EDDY = VT.w_eddy(I, J); break;
                        default: break;
                    }
                }

                // Wind inflow boundary conditions
                double H0 = 0.5 * (c.u0 * c.u0 + c.v0 * c.v0 + c.w0 * c.w0);
                if (c.open_wind_boundary) {
                    switch (IOR) {
                        case  1: H0 = HP(1, J, K)    + 0.5 / (c.dt * c.rdxn(0))    * (c.u_wind(K) + VEL_EDDY - UU(0, J, K));    break;
                        #if PB_FDSLOOPS_MUTANT == 15
                        case -1: H0 = HP(IBAR, J, K) + 0.5 / (c.dt * c.rdxn(IBAR)) * (c.u_wind(K) + VEL_EDDY - UU(IBAR, J, K)); break;
#else
                        case -1: H0 = HP(IBAR, J, K) - 0.5 / (c.dt * c.rdxn(IBAR)) * (c.u_wind(K) + VEL_EDDY - UU(IBAR, J, K)); break;
#endif
                        case  2: H0 = HP(I, 1, K)    + 0.5 / (c.dt * c.rdyn(0))    * (c.v_wind(K) + VEL_EDDY - VV(I, 0, K));    break;
                        case -2: H0 = HP(I, JBAR, K) - 0.5 / (c.dt * c.rdyn(JBAR)) * (c.v_wind(K) + VEL_EDDY - VV(I, JBAR, K)); break;
                        case  3: H0 = HP(I, J, 1)    + 0.5 / (c.dt * c.rdzn(0))    * (c.w_wind(K) + VEL_EDDY - WW(I, J, 0));    break;
                        case -3: H0 = HP(I, J, KBAR) - 0.5 / (c.dt * c.rdzn(KBAR)) * (c.w_wind(K) + VEL_EDDY - WW(I, J, KBAR)); break;
                        default: break;
                    }
                }

                switch (IOR) {
                    case 1:
                        if (UU(0, J, K) < 0.0) BXS(J, K) = P_EXTERNAL / WC.rho_f + KRES(1, J, K);
                        else                   BXS(J, K) = P_EXTERNAL / WC.rho_f + H0;
                        break;
                    case -1:
                        if (UU(IBAR, J, K) > 0.0) BXF(J, K) = P_EXTERNAL / WC.rho_f + KRES(IBAR, J, K);
                        else                      BXF(J, K) = P_EXTERNAL / WC.rho_f + H0;
                        break;
                    case 2:
                        if (VV(I, 0, K) < 0.0) BYS(I, K) = P_EXTERNAL / WC.rho_f + KRES(I, 1, K);
                        else                   BYS(I, K) = P_EXTERNAL / WC.rho_f + H0;
                        break;
                    case -2:
                        #if PB_FDSLOOPS_MUTANT == 17
                        if (VV(I, JBAR, K) >= 0.0) BYF(I, K) = P_EXTERNAL / WC.rho_f + KRES(I, JBAR, K);
#else
                        if (VV(I, JBAR, K) > 0.0) BYF(I, K) = P_EXTERNAL / WC.rho_f + KRES(I, JBAR, K);
#endif
                        else                      BYF(I, K) = P_EXTERNAL / WC.rho_f + H0;
                        break;
                    case 3:
                        if (WW(I, J, 0) < 0.0) BZS(I, J) = P_EXTERNAL / WC.rho_f + KRES(I, J, 1);
                        else                   BZS(I, J) = P_EXTERNAL / WC.rho_f + H0;
                        break;
                    case -3:
                        if (WW(I, J, KBAR) > 0.0) BZF(I, J) = P_EXTERNAL / WC.rho_f + KRES(I, J, KBAR);
                        else                      BZF(I, J) = P_EXTERNAL / WC.rho_f + H0;
                        break;
                    default: break;
                }
            }
        }
    }
}

}}   // namespace pressure_backend::fdsloops
