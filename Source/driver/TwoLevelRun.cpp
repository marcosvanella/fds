#include <cstdlib>
#include <memory>
// TwoLevelRun.cpp: see TwoLevelRun.H.
#include "LsFine.H"
#include "ObstList.H"
#include "TwoLevelRun.H"
#include "CompositeOps.H"

#include <AMReX_ParallelContext.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <cstdio>
#include <sstream>

#include "DriverAdapter.H"
#include "Hierarchy.H"
#include "PostRegridProjection.H"
#include "PressureBackendSolver.H"
#include "RegistryTransfer.H"
#include "RegridAmrCore.H"
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include "TimeLoopWiring.H"

extern "C" void fds_p_set_legacy_save_guard(int on);   // fds_step.f90
extern "C" void fds_p_bnd_rate(int ns, double* r);
extern "C" void fds_fine_set_bed();
extern "C" int fds_open_wind();
extern "C" double fds_dt_default(double cs);
extern "C" double fds_h_init();
extern "C" int fds_noise_rn(int nseed, double cs, int n, double* buf);   // fds_fine_level.f90
extern "C" void fds_fine_orphan_sums(double* s);   // fds_fine_level.f90 (ruling 2)
extern "C" void fds_fine_orphan_out(const double* s);   // s[6]: unmapped count, energy; remapped count, energy; level-0 ground convective energy, of |Q_CON_F|
extern "C" void fds_wind_at_height(double t, double z, double* u, double* v, double* w);   // fds_fine_level.f90
extern "C" void fds_fine_fuel_mdot(int lev, double* r);   // fds_fine_level.f90: fuel MLR and burning area over leaf walls, r[6] (diagnostic)
extern "C" void fds_p_pbar(int nm, int nmax, double* out, int* nout);
extern "C" int fds_p_reac_coef(int ns, double* c);   // fds_step.f90: species production per unit heat release (one reaction), else 0
extern "C" void fds_k_consts(double* tend, double* dtfill, double* dtmin, double* dt0, int* nzone);   // T_END etc.   // fds_step.f90: level-0 wall inflow rate of the last step (stage mean): total, per tracked species, through-flow
extern "C" void fds_k_phi(int nm, double* buf);
extern "C" int fds_ls_get_phi(int nx, int ny, double* buf);   // fds_ls_grid.f90: PHI of the global LS grid   // fds_kernels.f90: PHI_LS of a level-0 mesh
extern "C" int fds_find_tracked_species(const char* name);  // fds_mesh_query.f90; 0-based ZZ index or -1

namespace fdsamr {

namespace {

// Composite sums over the uncovered cells of all levels: volume integral of RHO and of RHO*ZZ_n. Each box sums its own cells; the per-box sums of all ranks
// (one owner per entry, so the reduction adds zeros only) are then added in (level, box) order: the same bits at any rank count (as 4048b2c87b, RTE factor).
std::vector<double> composite_mass(LevelRegistry& reg, int nl, int ns)
{
    const int nc = 1 + ns;
    std::vector<int> off(nl + 1, 0);
    for (int l = 0; l < nl; ++l) off[l + 1] = off[l] + static_cast<int>(reg.fields(l)["RHO"].size());
    std::vector<double> pb(static_cast<std::size_t>(off[nl]) * nc, 0.0);
    for (int l = 0; l < nl; ++l) {
        const Level& lv = reg.level(l);
        Fields& F = reg.fields(l);
        const amrex::iMultiFab* cov = reg.covered_mask(l);
        const double vol = lv.dx[0] * lv.dx[1] * lv.dx[2];
        for (amrex::MFIter mfi(F["RHO"]); mfi.isValid(); ++mfi) {
            auto r = F["RHO"].const_array(mfi);
            auto z = F["ZZ"].const_array(mfi);
            amrex::Array4<const int> cv, sd;
            if (cov) cv = cov->const_array(mfi);
            if (l > 0) sd = reg.side_data(l).mask().const_array(mfi);   // fine OBST cells hold the burner wall state (TimeLoop fine_burner_state), not gas
            double* m = pb.data() + static_cast<std::size_t>(off[l] + mfi.index()) * nc;
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                if ((cov && cv(i, j, k) != 0) || (l > 0 && sd(i, j, k, 0) != 0)) return;
                m[0] += r(i, j, k) * vol;
                for (int n = 0; n < ns; ++n) m[1 + n] += r(i, j, k) * z(i, j, k, n) * vol;
            });
        }
    }
    amrex::ParallelAllReduce::Sum(pb.data(), static_cast<int>(pb.size()), amrex::ParallelContext::CommunicatorSub());
    std::vector<double> m(nc, 0.0);
    for (int b = 0; b < off[nl]; ++b) for (int c = 0; c < nc; ++c) m[c] += pb[static_cast<std::size_t>(b) * nc + c];
    return m;
}

// Species check 7 (P0 W12): min and max of every tracked mass fraction over the uncovered gas cells of all levels (fine OBST cells masked, as composite_mass).
std::pair<double, double> leaf_z_range(LevelRegistry& reg, int nl, int ns)
{
    double zmin = 0.0, zmax = 0.0;
    if (ns > 0) { zmin = 1e300; zmax = -1e300; }
    for (int l = 0; l < nl && ns > 0; ++l) {
        Fields& F = reg.fields(l);
        const amrex::iMultiFab* cov = reg.covered_mask(l);
        for (amrex::MFIter mfi(F["ZZ"]); mfi.isValid(); ++mfi) {
            auto z = F["ZZ"].const_array(mfi);
            amrex::Array4<const int> cv, sd;
            if (cov) cv = cov->const_array(mfi);
            if (l > 0) sd = reg.side_data(l).mask().const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                if ((cov && cv(i, j, k) != 0) || (l > 0 && sd(i, j, k, 0) != 0)) return;
                for (int n = 0; n < ns; ++n) { zmin = std::min(zmin, z(i, j, k, n)); zmax = std::max(zmax, z(i, j, k, n)); }
            });
        }
    }
    amrex::ParallelAllReduce::Min(zmin, amrex::ParallelContext::CommunicatorSub());
    amrex::ParallelAllReduce::Max(zmax, amrex::ParallelContext::CommunicatorSub());
    return {zmin, zmax};
}

// Composite heat release rate (W): Q*vol over the uncovered level-0 cells and every finer gas cell (each point of the domain counted once; fine OBST cells masked).
double composite_hrr(LevelRegistry& reg, int nl)
{
    double q = 0.0;
    for (int l = 0; l < nl; ++l) {
        Fields& F = reg.fields(l);
        if (!F.has("Q")) return 0.0;
        const Level& lv = reg.level(l);
        const amrex::iMultiFab* cov = reg.covered_mask(l);
        const double vol = lv.dx[0] * lv.dx[1] * lv.dx[2];
        for (amrex::MFIter mfi(F["Q"]); mfi.isValid(); ++mfi) {
            auto a = F["Q"].const_array(mfi);
            amrex::Array4<const int> cv, sd;
            if (cov) cv = cov->const_array(mfi);
            if (l > 0) sd = reg.side_data(l).mask().const_array(mfi);   // fine OBST cells (as composite_mass)
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { if (!(cov && cv(i, j, k) != 0) && !(l > 0 && sd(i, j, k, 0) != 0)) q += a(i, j, k) * vol; });
        }
    }
    amrex::ParallelAllReduce::Sum(q, amrex::ParallelContext::CommunicatorSub());
    return q;
}


// Diagnostic (FDSTL_LEAFDIAG=N, default off): one line per level over the leaf gas cells (uncovered, not fine OBST): species mass sum(rho Z vol), Z and sum(Z)-1 range, TMP and RHO range, Q*V,
// and the cell of max TMP / min RHO / max Q. old_ba != nullptr: only the cells of levels >= 1 that were not in the old level (new cells), plus their max |TMP, RHO, Z - parent value|.
void leaf_diag(LevelRegistry& reg, const char* what, int s, double t, int ns, const std::vector<amrex::BoxArray>* old_ba = nullptr)
{
    for (int l = 0; l < reg.num_levels(); ++l) {
        if (old_ba && l < 1) continue;
        Fields& F = reg.fields(l);
        const amrex::iMultiFab* cov = reg.covered_mask(l);
        const double vol = reg.level(l).dx[0] * reg.level(l).dx[1] * reg.level(l).dx[2];
        const amrex::IntVect rr = l > 0 ? reg.level(l).ref_ratio_from_parent : amrex::IntVect(1);
        std::unique_ptr<amrex::MultiFab> pT, pR, pZ;   // parent values on the coarsened fine boxes
        if (old_ba) {
            const amrex::BoxArray cba = amrex::BoxArray(F["TMP"].boxArray()).coarsen(rr);
            pT.reset(new amrex::MultiFab(cba, F["TMP"].DistributionMap(), 1, 0)); pT->ParallelCopy(reg.fields(l - 1)["TMP"], 0, 0, 1);
            pR.reset(new amrex::MultiFab(cba, F["TMP"].DistributionMap(), 1, 0)); pR->ParallelCopy(reg.fields(l - 1)["RHO"], 0, 0, 1);
            pZ.reset(new amrex::MultiFab(cba, F["TMP"].DistributionMap(), ns, 0)); pZ->ParallelCopy(reg.fields(l - 1)["ZZ"], 0, 0, ns);
        }
        // v: 0 n, 1 QV, 2 Tmin, 3 Tmax, 4 Rmin, 5 Rmax, 6 Smin, 7 Smax, 8 dT, 9 dR, 10 dZ, then m(ns), Zmin(ns), Zmax(ns)
        std::vector<double> v(11 + 3 * ns, 0.0);
        v[2] = v[4] = v[6] = v[11 + ns] = 1e300; v[3] = v[5] = v[7] = v[11 + 2 * ns] = -1e300;
        for (int n = 0; n < ns; ++n) { v[11 + ns + n] = 1e300; v[11 + 2 * ns + n] = -1e300; }
        amrex::IntVect lt(-1), lr(-1), lq(-1); double qmax = -1e300;
        for (amrex::MFIter mfi(F["RHO"]); mfi.isValid(); ++mfi) {
            auto r = F["RHO"].const_array(mfi); auto z = F["ZZ"].const_array(mfi); auto tm = F["TMP"].const_array(mfi); auto q = F["Q"].const_array(mfi);
            amrex::Array4<const int> cv, sd; amrex::Array4<const amrex::Real> a, b, c;
            if (cov) cv = cov->const_array(mfi);
            if (l > 0) sd = reg.side_data(l).mask().const_array(mfi);
            if (old_ba) { a = pT->const_array(mfi); b = pR->const_array(mfi); c = pZ->const_array(mfi); }
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                if ((cov && cv(i, j, k) != 0) || (l > 0 && sd(i, j, k, 0) != 0)) return;
                if (old_ba && (*old_ba)[l].contains(amrex::IntVect(i, j, k))) return;
                ++v[0]; v[1] += q(i, j, k) * vol;
                if (tm(i, j, k) < v[2]) v[2] = tm(i, j, k);
                if (tm(i, j, k) > v[3]) { v[3] = tm(i, j, k); lt = amrex::IntVect(i, j, k); }
                if (r(i, j, k) < v[4]) { v[4] = r(i, j, k); lr = amrex::IntVect(i, j, k); }
                v[5] = std::max(v[5], r(i, j, k));
                if (q(i, j, k) > qmax) { qmax = q(i, j, k); lq = amrex::IntVect(i, j, k); }
                double sz = 0.0;
                for (int n = 0; n < ns; ++n) {
                    sz += z(i, j, k, n); v[11 + n] += r(i, j, k) * z(i, j, k, n) * vol;
                    v[11 + ns + n] = std::min(v[11 + ns + n], z(i, j, k, n)); v[11 + 2 * ns + n] = std::max(v[11 + 2 * ns + n], z(i, j, k, n));
                }
                v[6] = std::min(v[6], sz - 1.0); v[7] = std::max(v[7], sz - 1.0);
                if (old_ba) {
                    const amrex::IntVect pc = amrex::coarsen(amrex::IntVect(i, j, k), rr);
                    v[8] = std::max(v[8], std::abs(tm(i, j, k) - a(pc)));
                    v[9] = std::max(v[9], std::abs(r(i, j, k) - b(pc)));
                    for (int n = 0; n < ns; ++n) v[10] = std::max(v[10], std::abs(z(i, j, k, n) - c(pc, n)));
                }
            });
        }
        // sums, then mins, then maxes over ranks (locations are those of rank-local extremes: exact for one rank)
        std::vector<double> sm = {v[0], v[1], v[8], v[9], v[10]}, mn = {v[2], v[4], v[6]}, mx = {v[3], v[5], v[7]};
        for (int n = 0; n < ns; ++n) { sm.push_back(v[11 + n]); mn.push_back(v[11 + ns + n]); mx.push_back(v[11 + 2 * ns + n]); }
        const double dT = v[8], dR = v[9], dZ = v[10];
        amrex::ParallelAllReduce::Sum(sm.data(), 2, amrex::ParallelContext::CommunicatorSub());
        amrex::ParallelAllReduce::Sum(sm.data() + 5, ns, amrex::ParallelContext::CommunicatorSub());
        amrex::ParallelAllReduce::Min(mn.data(), static_cast<int>(mn.size()), amrex::ParallelContext::CommunicatorSub());
        amrex::ParallelAllReduce::Max(mx.data(), static_cast<int>(mx.size()), amrex::ParallelContext::CommunicatorSub());
        double dm3[3] = {dT, dR, dZ}; amrex::ParallelAllReduce::Max(dm3, 3, amrex::ParallelContext::CommunicatorSub());
        if (sm[0] == 0) continue;
        std::ostringstream os; os << std::setprecision(10);
        os << "[leafdiag] " << what << " step " << s << " t " << t << " L" << l << " n " << sm[0] << " QV " << sm[1] << " TMP " << mn[0] << " " << mx[0] << " @" << lt << " RHO " << mn[1] << " " << mx[1]
           << " @" << lr << " Qmax @" << lq << " sumZ-1 " << mn[2] << " " << mx[2];
        for (int n = 0; n < ns; ++n) os << " | Z" << n + 1 << " m " << sm[5 + n] << " min " << mn[3 + n] << " max " << mx[3 + n];
        if (old_ba) os << " | vs parent max|dT| " << dm3[0] << " max|dRHO| " << dm3[1] << " max|dZ| " << dm3[2];
        amrex::Print() << os.str() << "\n";
    }
}

void leaf_boxes(LevelRegistry& reg, const char* what, int s)
{
    for (int l = 1; l < reg.num_levels(); ++l) {
        std::ostringstream os; os << "[leafdiag] boxes " << what << " step " << s << " L" << l << " n " << reg.level(l).ba.size() << ":";
        for (int i = 0; i < static_cast<int>(reg.level(l).ba.size()); ++i) os << " " << reg.level(l).ba[i];
        amrex::Print() << os.str() << "\n";
    }
}

std::vector<fdsrt::ProjectionLevel> projection_levels(LevelRegistry& reg, int nl)
{
    std::vector<fdsrt::ProjectionLevel> lv;
    for (int l = 0; l < nl; ++l) {
        Fields& F = reg.fields(l);
        fdsrt::ProjectionLevel pl;
        pl.geom = reg.level(l).geom;
        pl.ref_ratio = l > 0 ? reg.level(l).ref_ratio_from_parent : amrex::IntVect(1);
        pl.vel = {&F["U"], &F["V"], &F["W"]};
        pl.D = &F["D"];
        pl.covered = reg.covered_mask(l);
        lv.push_back(pl);
    }
    return lv;
}

}  // namespace

namespace {

void fill_phi(LevelRegistry& reg, int nl);

// Cell-centred plot Multifab: RHO, U, V, W + OMEGA_Y = ∂u/∂z − ∂w/∂x (FDS edge form, 1/s), the mean of the cell's four y-edges.
// + TMP (K) and, with --tag-phi, the level-set PHI (k-constant, refreshed here). Owner viz: ω_y on mid-y xz + dark-grey mesh/AMR boxes. VORT (|ω|) stays TAG_VORT-only via ensure().
void write_rho_vel_plotfile(const std::string& path, LevelRegistry& reg, TimeLoop& loop, int nl)
{
    const bool phi = reg.fields(0).has("PHI");
    if (phi) fill_phi(reg, nl);
    amrex::Vector<amrex::MultiFab> mf(nl);
    amrex::Vector<const amrex::MultiFab*> pmf(nl);
    amrex::Vector<amrex::Geometry> geom(nl);
    amrex::Vector<int> level_steps(nl, loop.icyc());
    amrex::Vector<amrex::IntVect> ref_ratio(std::max(nl - 1, 0));
    for (int l = 0; l < nl; ++l) {
        const Level& Lv = reg.level(l);
        Fields& F = reg.fields(l);
        // Ghosts incl. coarse-fine and wall: the corrector-end exchange (code 6) on U, V, W, H, restored afterwards so a dump never changes the run.
        std::vector<amrex::MultiFab> keep;
        for (const char* n : {"U", "V", "W", "H"}) { keep.emplace_back(F[n].boxArray(), F[n].DistributionMap(), 1, F[n].nGrow()); amrex::MultiFab::Copy(keep.back(), F[n], 0, 0, 1, F[n].nGrow()); }
        loop.stage_exchange(l, 6, false);
        const int nz = F["ZZ"].nComp(), c0 = phi ? 7 : 6;   // after the base fields: HRRPUV, ZZ components, LEVEL, BOX_ID
        mf[l].define(Lv.ba, Lv.dm, c0 + 3 + nz, 0);
        amrex::MultiFab::Copy(mf[l], F["RHO"], 0, 0, 1, 0);
        amrex::Array<const amrex::MultiFab*, 3> faces{&F["U"], &F["V"], &F["W"]};
        amrex::average_face_to_cellcenter(mf[l], 1, faces);
        amrex::MultiFab::Copy(mf[l], F["TMP"], 0, 5, 1, 0);
        if (phi) amrex::MultiFab::Copy(mf[l], F["PHI"], 0, 6, 1, 0);
        amrex::MultiFab::Copy(mf[l], F["Q"], 0, c0, 1, 0); mf[l].mult(0.001, c0, 1);   // HRRPUV, kW/m3 (FDS Q is W/m3)
        amrex::MultiFab::Copy(mf[l], F["ZZ"], 0, c0 + 1, nz, 0);
        mf[l].setVal(l, c0 + 1 + nz, 1);
        for (amrex::MFIter mfi(mf[l]); mfi.isValid(); ++mfi) mf[l][mfi].setVal<amrex::RunOn::Host>(mfi.index(), mfi.validbox(), c0 + 2 + nz, 1);
        const double fz = 0.25 / Lv.dx[2], fx = 0.25 / Lv.dx[0];
        for (amrex::MFIter mfi(mf[l]); mfi.isValid(); ++mfi) {
            auto a = mf[l].array(mfi);
            auto u = F["U"].const_array(mfi);
            auto w = F["W"].const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                // edge (i,k) = (u(i,k)-u(i,k-1))/dz - (w(i,k)-w(i-1,k))/dx (AMReX face indices); sum over edges i..i+1, k..k+1 telescopes to this
                a(i, j, k, 4) = fz * (u(i, j, k + 1) + u(i + 1, j, k + 1) - u(i, j, k - 1) - u(i + 1, j, k - 1)) - fx * (w(i + 1, j, k) + w(i + 1, j, k + 1) - w(i - 1, j, k) - w(i - 1, j, k + 1));
            });
        }
        int q = 0;
        for (const char* n : {"U", "V", "W", "H"}) { amrex::MultiFab::Copy(F[n], keep[q], 0, 0, 1, F[n].nGrow()); ++q; }
        pmf[l] = &mf[l];
        geom[l] = Lv.geom;
        if (l + 1 < nl) ref_ratio[l] = reg.level(l + 1).ref_ratio_from_parent;
    }
    amrex::Vector<std::string> names{"RHO", "U", "V", "W", "OMEGA_Y", "TMP"};
    if (phi) names.push_back("PHI");
    names.push_back("HRRPUV");
    for (int n = 0; n < reg.fields(0)["ZZ"].nComp(); ++n) names.push_back("ZZ" + std::to_string(n + 1));   // tracked species mass fractions, FDS order
    names.push_back("LEVEL"); names.push_back("BOX_ID");
    amrex::WriteMultiLevelPlotfile(path, nl, pmf, names, geom, loop.time(), level_steps, ref_ratio);
    amrex::Print() << "PLOTFILE wrote " << path << " (RHO,U,V,W,OMEGA_Y,TMP[,PHI],HRRPUV,ZZ1..,LEVEL,BOX_ID cell-centred; OMEGA_Y = du/dz - dw/dx [1/s] from the 4 cell edges, " << nl << " level(s))\n";
}

}  // namespace



namespace {

// Cell-centred undivided |ω| from face U,V,W (FDS dump.f90 vorticity without /dx).
// Level-set PHI (k-constant): level 0 from the FDS PHI_LS of its mesh, finer levels by injection from the level below (ghosts included).
void fill_phi(LevelRegistry& reg, int nl)
{
    reg.fields(0).ensure("PHI");
    amrex::MultiFab& P0 = reg.fields(0)["PHI"];
    for (amrex::MFIter mfi(P0); mfi.isValid(); ++mfi) {
        auto p = P0.array(mfi);
        const amrex::Box vb = mfi.validbox();
        std::vector<double> buf(static_cast<size_t>(vb.length(0) + 2) * (vb.length(1) + 2));
        fds_k_phi(mfi.index() + 1, buf.data());
        const int ilo = vb.smallEnd(0), jlo = vb.smallEnd(1), nbx = vb.length(0) + 2, nby = vb.length(1) + 2;   // buf: cells 0:IBP1 x 0:JBP1 of the mesh (box-local, x fastest)
        amrex::LoopOnCpu(mfi.fabbox(), [&](int i, int j, int k) {
            const int ib = std::min(std::max(i - ilo + 1, 0), nbx - 1), jb = std::min(std::max(j - jlo + 1, 0), nby - 1);
            p(i, j, k) = buf[ib + static_cast<size_t>(jb) * nbx];
        });
    }
    std::vector<double> lsphi;   // --ls-fine: the finer levels take PHI from the global LS grid (mean of the LS cells of the level cell), so tagging follows the fine front
    int lsnx = 0, lsny = 0, lsR = 1;
    if (ls_fine_config().mode >= 3 && nl > 1) {
        lsR = ls_fine_config().ratio; lsnx = reg.level(0).geom.Domain().length(0) * lsR; lsny = reg.level(0).geom.Domain().length(1) * lsR;
        lsphi.resize(static_cast<size_t>(lsnx) * lsny);
        if (fds_ls_get_phi(lsnx, lsny, lsphi.data()) != 0) lsphi.clear();
    }
    int cum = 1;
    for (int l = 1; l < nl; ++l) {
        reg.fields(l).ensure("PHI");
        const amrex::IntVect r = reg.level(l).ref_ratio_from_parent;
        cum *= r[0];
        if (!lsphi.empty()) {
            const int q = lsR / cum;
            amrex::MultiFab& Pl = reg.fields(l)["PHI"];
            for (amrex::MFIter mfi(Pl); mfi.isValid(); ++mfi) {
                auto pf = Pl.array(mfi);
                amrex::LoopOnCpu(mfi.fabbox(), [&](int i, int j, int k) {
                    const int ic = std::min(std::max(i, 0), lsnx / q - 1), jc = std::min(std::max(j, 0), lsny / q - 1);
                    double s = 0.0;
                    for (int jj = 0; jj < q; ++jj) for (int ii = 0; ii < q; ++ii) s += lsphi[static_cast<size_t>(ic * q + ii) + static_cast<size_t>(jc * q + jj) * lsnx];
                    pf(i, j, k) = s / (q * q);
                });
            }
            continue;
        }
        amrex::MultiFab& Pf = reg.fields(l)["PHI"];
        amrex::BoxArray cba = reg.level(l).ba; cba.coarsen(r);
        amrex::MultiFab Pc(cba, reg.level(l).dm, 1, 1);
        Pc.setVal(0.0);
        Pc.ParallelCopy(reg.fields(l - 1)["PHI"], 0, 0, 1, 0, 1, reg.level(l - 1).geom.periodicity());
        for (amrex::MFIter mfi(Pf); mfi.isValid(); ++mfi) {
            auto pf = Pf.array(mfi); auto pc = Pc.const_array(mfi);
            amrex::LoopOnCpu(mfi.fabbox(), [&](int i, int j, int k) { pf(i, j, k) = pc(amrex::coarsen(i, r[0]), amrex::coarsen(j, r[1]), amrex::coarsen(k, r[2])); });
        }
    }
}

// ωx = (W_{j+1}-W_j) - (V_{k+1}-V_k), cyclic. Thin directions contribute 0.
void fill_undivided_vorticity(LevelRegistry& reg, int nl)
{
    for (int l = 0; l < nl; ++l) {
        Fields& F = reg.fields(l);
        F.ensure("VORT");
        const amrex::Periodicity per = reg.level(l).geom.periodicity();
        for (const char* n : {"U", "V", "W"}) F.fill_ghosts(n);   // one source per seam ghost
        const amrex::Box dom = reg.level(l).geom.Domain();
        const int dx_on = (dom.length(0) > 1) ? 1 : 0;
        const int dy_on = (dom.length(1) > 1) ? 1 : 0;
        const int dz_on = (dom.length(2) > 1) ? 1 : 0;
        amrex::MultiFab& Vo = F["VORT"];
        const amrex::MultiFab& U = F["U"];
        const amrex::MultiFab& V = F["V"];
        const amrex::MultiFab& W = F["W"];
        for (amrex::MFIter mfi(Vo); mfi.isValid(); ++mfi) {
            auto vo = Vo.array(mfi);
            auto u = U.const_array(mfi);
            auto v = V.const_array(mfi);
            auto w = W.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                const double ox = dy_on * (w(i, j + 1, k) - w(i, j, k)) - dz_on * (v(i, j, k + 1) - v(i, j, k));
                const double oy = dz_on * (u(i, j, k + 1) - u(i, j, k)) - dx_on * (w(i + 1, j, k) - w(i, j, k));
                const double oz = dx_on * (v(i + 1, j, k) - v(i, j, k)) - dy_on * (u(i, j + 1, k) - u(i, j, k));
                vo(i, j, k) = std::sqrt(ox * ox + oy * oy + oz * oz);
            });
        }
    }
}

}  // namespace

int two_level_run(TimeLoop& loop, const Level0& l0, double dt_setup, const TwoLevelOptions& o)
{
    int fails = 0;
    LevelRegistry& reg = loop.registry();
    const int ns = l0.dom.n_total;
    auto probe = [&](const char* tag) {   // FDSTL_TLDIAG: periodic z ghost of U on level 0 at x face 5 (must equal the opposite valid row)
        if (!std::getenv("FDSTL_TLDIAG")) return;
        const amrex::MultiFab& u = reg.fields(0)["U"];
        for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) {
            auto a = u.const_array(mfi);
            amrex::Print() << "PROBE " << tag << ": U(5,0,-1) " << a(5, 0, -1) << " U(5,0,15) " << a(5, 0, 15) << " U(5,0,16) " << a(5, 0, 16) << " U(5,0,0) " << a(5, 0, 0) << "\n";
        }
    };
    probe("entry");
    const amrex::Box dom = l0.geom.Domain();

    // ---- the hierarchy: &AMR / &AMR_REGION -> AmrParams, level-0 meshes -> Hierarchy -> RegridAmrCore; levels 1..max_level by regrid_dynamic (nested tags)
    if (o.max_level < 0 || o.max_level > 8) { amrex::Print() << "MULTI-LEVEL: max_level must be 0..8, got " << o.max_level << "\n"; return 1; }
    const char* tag = (o.max_level == 1) ? "TWO-LEVEL" : "MULTI-LEVEL";
    amrex::IntVect rr(o.ratio, o.ratio, o.ratio);
    if (ls_fine_config().mode >= 3) { int r = 1; for (int l = 0; l < o.max_level; ++l) r *= o.ratio; ls_fine_config().ratio = r; }   // --ls-fine: LS cell = finest possible cell
    for (int d = 0; d < 3; ++d) if (dom.length(d) == 1) rr[d] = 1;
    const int plo[3] = {o.patch[0], dom.smallEnd(1), o.patch[2]}, phi_[3] = {o.patch[1], dom.bigEnd(1), o.patch[3]};
    amrex::Box patch(amrex::IntVect(plo[0], plo[1], plo[2]), amrex::IntVect(phi_[0], phi_[1], phi_[2]));
    const int p2lo[3] = {o.patch2[0], dom.smallEnd(1), o.patch2[2]}, p2hi[3] = {o.patch2[1], dom.bigEnd(1), o.patch2[3]};
    amrex::Box patch2(amrex::IntVect(p2lo[0], p2lo[1], p2lo[2]), amrex::IntVect(p2hi[0], p2hi[1], p2hi[2]));
    if (o.max_level >= 2 && !patch.contains(patch2)) {
        amrex::Print() << tag << ": patch2 must lie inside patch (level-0 cells)\n";
        return 1;
    }
    load_obst_list();   // the bed list (W2a) is needed for the refinable region below
    if (o.bed.size() == 6) pb::bed_list().push_back({o.bed[0], o.bed[1], o.bed[2], o.bed[3], o.bed[4], o.bed[5]});
    if (!pb::bed_list().empty()) { load_nofine(1); fds_fine_set_bed(); }   // D-010 exemption of analytic burners and their own fine-face state, bed runs only
    std::ostringstream amr;
    amr << "&AMR MAX_LEVEL=" << o.max_level << ", REF_RATIO=" << o.ratio << ", BLOCKING_FACTOR=" << o.blocking << ", MAX_GRID_SIZE=" << o.maxsize
        << ", N_PROPER=1, N_ERROR_BUF=" << (std::isnan(o.tag_phi) ? 1 : o.tag_phi_buf) << ", REGRID_INTERVAL=" << o.regrid_every << ", POST_REGRID_PROJECTION='" << o.projection << "'";
    // Adaptive TAG_*: vortex-ring TAG_VORT path, else blob TAG_SPEC_*. Static (regrid_every==0) uses nested set_tag_function only.
    const bool use_tag_phi = !std::isnan(o.tag_phi);
    const bool use_tag_vort = o.tag_vort_auto || !std::isnan(o.tag_vort);
    const bool use_tag_temp = !std::isnan(o.tag_temp), use_tag_obst = o.tag_obst >= 0;
    const bool adaptive_free = use_tag_vort || use_tag_phi || (o.regrid_every > 0 && (use_tag_temp || use_tag_obst)) || !pb::bed_list().empty();   // criteria that move the fine region by themselves: the refinable region is the whole domain
    double tag_vort_thr = o.tag_vort;
    if (o.regrid_every > 0 && use_tag_vort) {
        fill_undivided_vorticity(reg, 1);   // L0 |ω| before AmrParams (auto thr needs max)
        if (o.tag_vort_auto) {
            // Owner overnight: 0.45*max hugs ring core tighter than 0.35 (less fat L1).
            tag_vort_thr = 0.45 * reg.fields(0)["VORT"].max(0, 0);   // MultiFab::max is global; level 0 by design: a tagging threshold fixed before any fine level exists, not an output statistic
            amrex::Print() << tag << ": TAG_VORT auto thr = 0.45 * max|VORT|_L0 = " << tag_vort_thr << "\n";
        }
        amr << ", TAG_VORT=" << std::setprecision(17) << tag_vort_thr << ", TAG_KEEP=0.85";
    } else if (o.regrid_every > 0 && !use_tag_temp && !use_tag_obst && !use_tag_phi) {   // default blob criterion only when no other adaptive criterion (--tag-temp, --tag-obst) is asked for
        amr << ", TAG_SPEC_ID='TRACER', TAG_SPEC_DIFF=0.05, TAG_KEEP=0.8";
    }
    if (o.regrid_every > 0 && use_tag_phi) amr << ", TAG_PHI_DIFF=" << std::setprecision(17) << o.tag_phi << ", TAG_KEEP=0.85";
    if (o.regrid_every > 0 && !std::isnan(o.tag_temp)) amr << ", TAG_TMP_ABOVE=" << std::setprecision(17) << o.tag_temp - 293.15;   // --tag-temp (TMPA 293.15 = TagFields::tmp_ambient)
    amr << " /\n";
    // Vortex adaptive: whole-domain refinable region so the ring can leave the static patch.
    // --levels 1 (max_level 0): level 0 alone through the same loop and plotfiles, no &AMR_REGION (its LEVEL must be >= 1)
    if (o.max_level > 0) {
        amr << std::setprecision(17) << "&AMR_REGION XB=";
        if (adaptive_free) {
            for (int d = 0; d < 3; ++d) amr << l0.geom.ProbLo(d) << "," << (d == 2 && !std::isnan(o.fine_zmax) ? std::min(l0.geom.ProbHi(d), o.fine_zmax) : l0.geom.ProbHi(d)) << (d < 2 ? "," : "");   // --fine-zmax caps the refinable region in z
        } else {
            for (int d = 0; d < 3; ++d) amr << l0.geom.ProbLo(d) + plo[d] * l0.dx[d] << "," << l0.geom.ProbLo(d) + (phi_[d] + 1) * l0.dx[d] << (d < 2 ? "," : "");
        }
        amr << ", LEVEL=" << o.max_level << " /\n";
    }
    fdsrt::Report rep;
    fdsrt::AmrParams ap = fdsrt::parse_amr_params(amr.str(), rep);
    std::vector<fdsrt::MeshInput> meshes;
    for (const MeshInfo& m : l0.mesh) {
        fdsrt::MeshInput mi;
        for (int d = 0; d < 3; ++d) mi.ijk[d] = m.ijk[d];
        for (int q = 0; q < 6; ++q) mi.xb[q] = m.xb[q];
        mi.rank = m.rank;
        meshes.push_back(mi);
    }
    std::array<bool, 3> per{l0.dom.periodic[0] != 0, l0.dom.periodic[1] != 0, l0.dom.periodic[2] != 0};
    fdsrt::Hierarchy h;
    if (!rep.ok() || !fdsrt::build_hierarchy_from_meshes(meshes, ap, per, h, rep) || !rep.ok()) {
        for (const auto& e : rep.errors) amrex::Print() << tag << " hierarchy error: " << e << "\n";
        return 1;
    }
    for (const auto& w : rep.warnings) amrex::Print() << tag << " hierarchy warning: " << w << "\n";
    fdsrt::RegridAmrCore core(h, ap);
    RegistryTransfer tr(reg, l0.dom.n_tracked);
    // Stage arrays a regrid does not transfer (D, DS, RSUM, MU, KRES, H, HS) of a new fine level: parent cell value (piecewise constant), any level >= 1 (D-053).
    tr.derive = [&](int level) {
        if (level < 1) return;
        const Level& lf = reg.level(level);
        amrex::BoxArray cba = lf.ba; cba.coarsen(lf.ref_ratio_from_parent);
        const amrex::IntVect ratio = lf.ref_ratio_from_parent;
        for (const char* nm : {"D", "DS", "RSUM", "MU", "KRES", "H", "HS"}) {
            if (!reg.fields(level).has(nm) || !reg.fields(level - 1).has(nm)) continue;
            amrex::MultiFab& Df = reg.fields(level)[nm];
            amrex::MultiFab ct(cba, lf.dm, Df.nComp(), 0);
            ct.setVal(0.0);
            ct.ParallelCopy(reg.fields(level - 1)[nm], 0, 0, Df.nComp(), 0, 0, reg.level(level - 1).geom.periodicity());
            for (amrex::MFIter mfi(Df); mfi.isValid(); ++mfi) {
                auto d = Df.array(mfi); auto c = ct.const_array(mfi); const int nc = Df.nComp();
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { for (int n = 0; n < nc; ++n) d(i, j, k, n) = c(amrex::coarsen(i, ratio[0]), amrex::coarsen(j, ratio[1]), amrex::coarsen(k, ratio[2]), n); });
            }
        }
    };
    core.set_data_transfer(&tr);
    // D-090: map TAG_SPEC_ID -> ZZ component before set_tag_fields (no hardcode). Abort if the name is unknown.
    int zz_comp_tag = 0;
    if (!ap.tag_spec_id.empty()) {
        const int idx = fds_find_tracked_species(ap.tag_spec_id.c_str());
        if (idx < 0) {
            amrex::Abort("D-090: TAG_SPEC_ID='" + ap.tag_spec_id + "' is not a tracked species or tracer (SPECIES_MIXTURE%ID)");
        }
        zz_comp_tag = idx;
        amrex::Print() << tag << ": TAG_SPEC_ID='" << ap.tag_spec_id << "' -> ZZ component " << zz_comp_tag << "\n";
    }
    // Registry fields for ErrorEst TAG_* (preferred when AmrParams has a criterion). Fallback nested set_tag_function for static.
    core.set_tag_fields([&](int lev) {
        Fields& F = reg.fields(lev);
        const amrex::Periodicity per = reg.level(lev).geom.periodicity();
        fdsrt::TagFields f;
        if (F.has("TMP")) { F["TMP"].FillBoundary(per); f.tmp = &F["TMP"]; }
        if (F.has("RHO")) { F["RHO"].FillBoundary(per); f.rho = &F["RHO"]; }
        if (F.has("ZZ")) {
            F["ZZ"].FillBoundary(per);
            f.zz = &F["ZZ"];
            f.zz_comp = zz_comp_tag;   // D-090: from TAG_SPEC_ID via fds_find_tracked_species
        }
        // D-090: TAG_HRRPUV field — HRRPUV if registered, else FDS Q (W/m3; thr×1000 in apply_tag_criteria)
        if (F.has("HRRPUV")) f.hrrpuv = &F["HRRPUV"];
        else if (F.has("Q")) f.hrrpuv = &F["Q"];
        if (F.has("VORT")) { F["VORT"].FillBoundary(per); f.vort = &F["VORT"]; }
        if (F.has("PHI")) { F["PHI"].FillBoundary(per); f.phi = &F["PHI"]; }
        for (int d = 0; d < 3; ++d) f.dirs[d] = (dom.length(d) == 1) ? 0 : 1;
        return f;
    });
    // Nested tags: level 0 tags `patch` (creates level 1); levels 1..max_level-1 tag the refined image of `patch2` (creates levels 2..max_level).
    core.set_tag_function([&](int lev, amrex::TagBoxArray& tags, amrex::Real) {
        if (lev >= o.max_level) return;
        amrex::IntVect cum(1, 1, 1);
        for (int l = 0; l < lev; ++l) for (int d = 0; d < 3; ++d) cum[d] *= rr[d];
        const amrex::Box src = (lev == 0) ? patch : patch2;
        amrex::Box tagbox(amrex::IntVect(src.smallEnd(0) * cum[0], src.smallEnd(1) * cum[1], src.smallEnd(2) * cum[2]),
                          amrex::IntVect(src.bigEnd(0) * cum[0] + cum[0] - 1, src.bigEnd(1) * cum[1] + cum[1] - 1, src.bigEnd(2) * cum[2] + cum[2] - 1));
        for (amrex::MFIter mfi(tags); mfi.isValid(); ++mfi) {
            auto a = tags.array(mfi);
            const amrex::Box b = mfi.validbox() & tagbox;
            if (b.ok()) amrex::LoopOnCpu(b, [&](int i, int j, int k) { a(i, j, k) = amrex::TagBox::SET; });
        }
    });
    load_obst_list();
    core.set_tag_obst(o.tag_obst);
    core.set_tag_obst_snap(o.tag_obst_snap != 0 && o.max_level >= 1);
    core.set_nofine(o.nofine);
    if (o.nofine >= 0 && !pb::nofine_list().empty()) amrex::Print() << tag << ": D-010 " << pb::nofine_list().size() << " stateful VENT(s) stay on level 0 (no-refine buffer " << o.nofine << " level-0 cells)\n";
    // fine boxes are cut at the level-0 mesh boundaries (RegridAmrCore::chop_at_level0); each goes to the rank of the mesh under it (FDS PROCESS), so its parent walls,
    // ambient state and devices are local
    const fdsrt::RegridAmrCore::DmFn owner_dm = [&](int lev, const amrex::BoxArray& ba) {
        amrex::IntVect r(1);
        for (int l = 0; l < lev; ++l) r *= core.refRatio(l);
        amrex::Vector<int> pmap(ba.size(), 0);
        for (int i = 0; i < static_cast<int>(ba.size()); ++i) {
            const amrex::Box b = amrex::coarsen(ba[i], r);
            for (int m = 0; m < static_cast<int>(l0.ba.size()); ++m)
                if (l0.ba[m].contains(b.smallEnd())) { pmap[i] = l0.dm[m]; break; }
        }
        // Test only (FDS_AMR_DEBUG_FOREIGN_BOXES set): a fine box that touches no level-0 domain face (no hidden walls in an OBST-free case) goes to the next
        // rank, so the leaf devices over another rank's mesh (dump.f90 FOREIGN pass, main.f90 LEAF_SUM_IN_RANK_ORDER) are exercised.
        static const bool dbg_foreign = std::getenv("FDS_AMR_DEBUG_FOREIGN_BOXES") != nullptr;
        const int np = amrex::ParallelDescriptor::NProcs();
        if (dbg_foreign && np > 1) {
            const amrex::Box dom = core.Geom(0).Domain();
            int nmoved = 0;
            for (int i = 0; i < static_cast<int>(ba.size()); ++i) {
                const amrex::Box b = amrex::coarsen(ba[i], r);
                bool inner = true;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) inner = inner && b.smallEnd(d) > dom.smallEnd(d) && b.bigEnd(d) < dom.bigEnd(d);
                if (!inner) continue;
                pmap[i] = (pmap[i] + 1) % np; ++nmoved;
                amrex::Print() << "DEBUG-FOREIGN level " << lev << " box " << i << " " << ba[i] << " to rank " << pmap[i] << "\n";
            }
            amrex::Print() << "DEBUG-FOREIGN level " << lev << ": " << nmoved << " of " << ba.size() << " boxes on a non-parent rank\n";
        }
        return amrex::DistributionMapping(pmap);
    };
    core.init_static(reg, false, owner_dm, &l0.dm);
    probe("after init_static");

    // ---- level 0 alone: one corrector-form pass for the D of the initial state (projection target and parent value of new levels' D)
    loop.set_state(0.0, dt_setup, 0);
    probe("before first prime");
    loop.prime_levels();
    probe("after first prime");
    fdsrt::PressureBackendSolver solver;
    std::copy(loop.pressure_bc(), loop.pressure_bc() + 6, solver.face_bc.begin());   // OPEN faces Dirichlet as in the time-step solve
    solver.face_neumann = loop.pressure_face_neumann();   // mixed face: Neumann over its SOLID cells, as in the time-step solve
    fdsrt::PostRegridProjectionState pst;
    pst.solver = &solver;
    pst.D = [&](int l) { return &reg.fields(l)["D"]; };
    pst.options.accept_abs = 1.0e-9;
    pst.before = [&](std::vector<fdsrt::ProjectionLevel>& lv) { composite_zero_solid_faces(reg, 0, static_cast<int>(lv.size()), {"U", "V", "W"}, 0); };   // faces in an OBST (prolonged, or the level-0 IC)
    fdsrt::install_post_regrid_projection(core, reg, pst);   // D-063: the hook of RegridAmrCore::regrid_dynamic
    loop.zero_closed_wall_normal();   // IC: Neumann faces carry no net volume flux into the first projection
    // One AmrCore::regrid grows at most one new fine level; repeat until max_level (or stuck).
    core.set_force_tag_function((o.init_patch || o.patch_set) && o.regrid_every > 0);   // --init-patch or an explicit --patch/--patch2: patch/patch2 forced at step 0 only (adaptive runs would otherwise ignore them)
    for (int pass = 0; pass < o.max_level; ++pass) {
        if (use_tag_vort) fill_undivided_vorticity(reg, reg.num_levels());
        if (use_tag_phi) fill_phi(reg, reg.num_levels());
        core.regrid_dynamic(0.0);
        if (core.finestLevel() >= o.max_level) break;
        if (pass > 0 && !core.last_outcome().changed) break;
    }
    core.set_force_tag_function(false);
    probe("after regrid_dynamic");
    const bool stock_ic = true;   // every level's initial velocity is the stock one (set below; init.f90 always takes U_WIND(K), U0 without a wind), none interpolated
    {   // init.f90: U(:,:,K) = U_WIND(K) at the cell height of each mesh; a fine level is set as stock sets its own mesh, not prolonged from level 0
        for (int l = 1; l <= core.finestLevel(); ++l) {
            const amrex::Geometry& g = core.Geom(l);
            const char* vn[6] = {"U", "V", "W", "US", "VS", "WS"};   // init.f90 US = U before INITIAL_NOISE: the estimate keeps the noise-free profile
            for (int n = 0; n < 6; ++n) {
                const int d = n % 3;
                amrex::MultiFab& mf = reg.fields(l)[vn[n]];
                for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
                    auto a = mf.array(mfi);
                    amrex::LoopOnCpu(mfi.fabbox(), [&](int i, int j, int k) {
                        double w[3]; fds_wind_at_height(0.0, (k + 0.5) * g.CellSize(2) + g.ProbLo(2), &w[0], &w[1], &w[2]);
                        a(i, j, k) = w[d];
                    });
                }
            }
            // main.f90 INITIAL_NOISE on every fine box as on a mesh of that number (its index in the level): same seed, wasted calls and loop order as init.f90
            const amrex::Real* h = g.CellSize(); const double cs = std::cbrt(h[0] * h[1] * h[2]);
            amrex::MultiFab& U = reg.fields(l)["U"]; amrex::MultiFab& V = reg.fields(l)["V"]; amrex::MultiFab& W = reg.fields(l)["W"];
            for (amrex::MFIter mfi(reg.fields(l)["RHO"]); mfi.isValid(); ++mfi) {
                const amrex::Box b = mfi.validbox(); const int i0 = b.smallEnd(0), j0 = b.smallEnd(1), k0 = b.smallEnd(2);
                const int ib = b.length(0), jb = b.length(1), kb = b.length(2);
                const int n = ib * (jb - 1) * (kb - 1) + (ib - 1) * jb * (kb - 1) + (ib - 1) * (jb - 1) * kb;
                std::vector<double> rn(std::max(n, 1));
                if (!fds_noise_rn(mfi.index() + 1, cs, n, rn.data())) continue;
                auto u = U.array(mfi); auto v = V.array(mfi); auto w = W.array(mfi);
                const double rdx = 1.0 / h[0], rdy = 1.0 / h[1], rdz = 1.0 / h[2];
                int q = 0;   // FDS face index F -> F + lo, cell index C -> C - 1 + lo
                for (int K = 1; K <= kb - 1; ++K) for (int J = 1; J <= jb - 1; ++J) for (int I = 1; I <= ib; ++I) { const double r = rn[q++];
                    w(I - 1 + i0, J - 1 + j0, K + k0) -= r * rdy; w(I - 1 + i0, J + j0, K + k0) += r * rdy; v(I - 1 + i0, J + j0, K - 1 + k0) += r * rdz; v(I - 1 + i0, J + j0, K + k0) -= r * rdz; }
                for (int K = 1; K <= kb - 1; ++K) for (int J = 1; J <= jb; ++J) for (int I = 1; I <= ib - 1; ++I) { const double r = rn[q++];
                    w(I - 1 + i0, J - 1 + j0, K + k0) -= r * rdx; w(I + i0, J - 1 + j0, K + k0) += r * rdx; u(I + i0, J - 1 + j0, K - 1 + k0) += r * rdz; u(I + i0, J - 1 + j0, K + k0) -= r * rdz; }
                for (int K = 1; K <= kb; ++K) for (int J = 1; J <= jb - 1; ++J) for (int I = 1; I <= ib - 1; ++I) { const double r = rn[q++];
                    v(I - 1 + i0, J + j0, K - 1 + k0) -= r * rdx; v(I + i0, J + j0, K - 1 + k0) += r * rdx; u(I + i0, J - 1 + j0, K - 1 + k0) += r * rdy; u(I + i0, J + j0, K - 1 + k0) -= r * rdy; }
            }
        }
        // faces touching an OBST cell carry no velocity (as pst.before does for the projected IC): the noise is added after the regrid, and a noisy face inside a burner OBST
        // left |div u - D| = 0.118 on its solid cells at step 1 (SEAM_burner_L2_m4r4)
        composite_zero_solid_faces(reg, 1, core.finestLevel() + 1, {"U", "V", "W"}, 0);
        composite_zero_solid_faces(reg, 1, core.finestLevel() + 1, {"US", "VS", "WS"}, 0);
    }
    const int finest = core.finestLevel();
    if (finest != o.max_level) {
        amrex::Print() << tag << ": expected finest level " << o.max_level << ", got " << finest << " (tagging / blocking / region)\n";
        return 1;
    }
    if (!pb::bed_list().empty() && finest >= 1) {   // W2a (PV-10): particles are inserted on L_ins = the finest level, which must cover every bed cell
        const amrex::Geometry& g = core.Geom(finest);
        bool ok = true;
        for (auto const& xb : pb::bed_list()) { const amrex::Box b = pb::obst_tag_box(g, xb, 0) & g.Domain(); ok = ok && (!b.ok() || core.boxArray(finest).contains(b)); }
        amrex::Print() << tag << ": particle bed level L_ins " << finest << ", dx " << g.CellSize(0) << " " << g.CellSize(1) << " " << g.CellSize(2) << (ok ? "" : ": does NOT cover the whole bed") << "\n";
        if (!ok) return 1;
    }
    amrex::Print() << tag << " regrid_dynamic: changed " << core.last_outcome().changed << ", new cells " << core.last_outcome().new_cells << ", retained " << core.last_outcome().retained_cells
                   << ", post-regrid hook calls " << pst.calls;
    if (pst.calls > 0) amrex::Print() << " (max|div u - D| before " << pst.last.div_before << ", after " << pst.last.div_after << ", accepted " << pst.last.accepted << ")";
    amrex::Print() << "\n";
    for (int l = 1; l <= finest; ++l) loop.bind_level(l);   // mesh offset + derive of parent stage arrays already done by RegistryTransfer::derive at make_level
    fdsrt::CfHookStats cfs;
    fdsrt::install_cf_ghost_hooks(loop, finest, fdsrt::ThermoProvider(), &cfs);
    loop.set_flux_overwrite(o.overwrite);
    probe("after bind/hooks");
    for (int l = 1; l <= finest; ++l)
        amrex::Print() << tag << ": level " << l << " = " << reg.level(l).ba.size() << " box(es), " << reg.level(l).ba.numPts() << " cells, ratio from parent ("
                       << reg.level(l).ref_ratio_from_parent[0] << "," << reg.level(l).ref_ratio_from_parent[1] << "," << reg.level(l).ref_ratio_from_parent[2] << ")\n";
    amrex::Print() << tag << ": levels 0.." << finest << ", blocking " << o.blocking << ", flux overwrite " << (o.overwrite ? "on" : "off") << ", projection " << o.projection
                   << ", bind_only " << (o.bind_only ? "yes" : "no") << "\n";
    const int nl = finest + 1;
    if (o.bind_only) {
        amrex::Print() << tag << " BIND-ONLY PASS: registry has " << reg.num_levels() << " level(s), TimeLoop bound " << loop.num_levels() << "\n";
        return (reg.num_levels() == nl && loop.num_levels() == nl) ? 0 : 1;
    }
    if (o.boundary_test != 0) {
        // stage_boundary(1,3) and (1,6): the boundary routines FDS runs after the exchanges of the predictor and the corrector, on the bound fine level (fds_p_save_uvw on a fine box number, and with
        // FDSTL_SKIPAFT=0 also fill_omesh / VELOCITY_BC). Test 2 restores the pre-S14.1 guard (level-0 numbers only): the call must then stop the run. Fine-level calls must leave level 0 alone.
        if (o.boundary_test == 2) fds_p_set_legacy_save_guard(1);
        std::vector<amrex::MultiFab> keep;
        const char* vn[3] = {"U", "V", "W"};
        for (int d = 0; d < 3; ++d) { const amrex::MultiFab& m = reg.fields(0)[vn[d]]; keep.emplace_back(m.boxArray(), m.DistributionMap(), m.nComp(), m.nGrow()); amrex::MultiFab::Copy(keep.back(), m, 0, 0, m.nComp(), m.nGrow()); }
        amrex::Print() << "STAGE-BOUNDARY-FINE: calling stage_boundary(1,3) and stage_boundary(1,6)" << (o.boundary_test == 2 ? " with the pre-S14.1 level-0-only guard" : "") << "\n";
        loop.stage_boundary(1, 3);
        loop.stage_boundary(1, 6);
        double worst = 0.0;
        for (int d = 0; d < 3; ++d) {
            const amrex::MultiFab& m = reg.fields(0)[vn[d]];
            amrex::MultiFab diff(m.boxArray(), m.DistributionMap(), m.nComp(), m.nGrow());
            amrex::MultiFab::Copy(diff, m, 0, 0, m.nComp(), m.nGrow());
            amrex::MultiFab::Subtract(diff, keep[d], 0, 0, m.nComp(), m.nGrow());
            const double w = diff.norminf(0, m.nGrow());
            amrex::Print() << "STAGE-BOUNDARY-FINE: level 0 " << vn[d] << " (ghost layers included) changed by max |diff| " << w << "\n";
            worst = std::max(worst, w);
        }
        amrex::Print() << "STAGE-BOUNDARY-FINE: returned, level 0 U/V/W max change " << worst << (worst == 0.0 ? " (PASS)" : " (FAIL)") << "\n";
        return worst == 0.0 ? 0 : 1;
    }

    // ---- the fine level's own D, KRES, MU, ... by one corrector-form pass, then the post-regrid projection of the two-level velocity with that D (the body of the hook, Role 3 PostRegridProjection)
    double dt_start = dt_setup;   // read.f90: default DT = CFL_MAX * CHARACTERISTIC_CELL_SIZE / VEL_CHAR, the cell size being the smallest of all meshes, so the finest level's
    {                             // (FDS set it from level 0 alone); a DT or time ramp from the input stays
        double cs = 1.0e30;
        for (int l = 0; l <= core.finestLevel(); ++l) { const amrex::Real* h = core.Geom(l).CellSize(); cs = std::min(cs, std::cbrt(h[0] * h[1] * h[2])); }
        const double d = fds_dt_default(cs);
        if (d > 0.0) dt_start = d;
        amrex::Print() << tag << ": starting dt " << dt_start << " (cell size " << cs << ")\n";
    }
    loop.set_state(0.0, dt_start, 0);
    loop.t0_stock_rho_levels();
    loop.prime_levels();
    // main.f90 inserts the particles (:634) after the initial radiation iterations (:597-622): the bed on L_ins after rad_init_levels (W2: the &INIT beds; with --bed too)
    auto insert_bed = [&] { if (!pb::bed_list().empty() && finest >= 1) { loop.insert_bed(finest); loop.particle_report(); } };
    const bool bed_first = std::getenv("FDSTL_BED_BEFORE_RAD_INIT") != nullptr;   // negative control (run_rad_init_csiro_check.sh)
    if (bed_first) insert_bed();
    loop.rad_init_levels();
    if (!bed_first) insert_bed();
    loop.t0_divergence_levels();   // the fine levels' t=0 D the stock way, after their initial radiation (not prolonged from level 0)
    probe("after second prime");
    if (std::getenv("FDSTL_TLDIAG")) {
        for (int l = 0; l < nl; ++l)
            for (const char* nm : {"RHO", "TMP", "ZZ", "D", "DS", "DDDT", "KRES", "MU", "U", "V", "W"}) {
                if (!reg.fields(l).has(nm)) continue;
                const amrex::MultiFab& mf = reg.fields(l)[nm];
                amrex::Print() << "TLDIAG prime [covered cells included] L" << l << " " << nm << " min " << mf.min(0, 0) << " max " << mf.max(0, 0) << (mf.contains_nan() ? " NAN" : "") << "\n";
            }
    }
    {
        std::vector<fdsrt::ProjectionLevel> lv = projection_levels(reg, nl);
        fdsrt::ProjectionOptions po = pst.options;
        po.enabled = o.projection == "ON" || (o.projection == "AUTO" && !stock_ic);   // AUTO: the start state is stock's, so project only when a level's initial velocity is interpolated
        loop.zero_closed_wall_normal();   // fine faces prolonged onto INERT/MIRROR, then the projection
        if (po.enabled) pst.before(lv);   // OFF: leave the level-0 IC as FDS set it (zeroing U alone, without a projection or US, breaks L0 mass)
        const fdsrt::ProjectionReport pr = fdsrt::project_after_regrid(lv, &solver, po);
        amrex::Print() << tag << " initial projection (" << o.projection << "): ran " << pr.ran << " solved " << pr.solved << " max|div u - D| before " << pr.div_before << " after " << pr.div_after
                       << " (level 0 " << (pr.div_after_level.empty() ? 0.0 : pr.div_after_level[0]) << ", level 1 " << (pr.div_after_level.size() > 1 ? pr.div_after_level[1] : 0.0) << "), accepted " << pr.accepted
                       << ", rhs mean/|rhs| " << pr.rhs_sum_rel << ", solver " << pr.solver << "\n";
        if (po.enabled && pr.ran && pr.solved) {
            for (int l = 0; l < nl; ++l) {
                Fields& F = reg.fields(l);
                const char* st[3] = {"US", "VS", "WS"}; const char* vn[3] = {"U", "V", "W"};
                for (int d = 0; d < 3; ++d) { amrex::MultiFab::Copy(F[st[d]], F[vn[d]], 0, 0, 1, 0); F.fill_ghosts(st[d]); F.fill_ghosts(vn[d]); }
            }
            loop.prime_levels();   // ghost faces of the projected velocity, D of the projected state
        }
        if (po.enabled && pr.ran && !pr.accepted) { ++fails; amrex::Print() << tag << ": the initial projection was not accepted\n"; }
    }

    {   // init.f90: H = HS = 0.5 (U0^2+V0^2+W0^2) at the first step (stock solves no pressure before it); the set-up passes above left the solve's H on every level,
        // and pres.f90 takes the start-of-step H (HP) as the wind inflow value of the OPEN faces
        const double h0 = fds_h_init();
        for (int l = 0; l < nl; ++l) for (const char* n : {"H", "HS"}) { reg.fields(l)[n].setVal(h0); }
    }
    if (std::getenv("FDSTL_T0RAD")) loop.rad_t0_report();
    // ---- the run
    const std::vector<double> m0 = composite_mass(reg, nl, ns);
    double zrange_min = 0.0, zrange_max = 0.0;   // species check 7: leaf 0 <= Z <= 1 (within 1e-10) after every step
    // Open-domain mass balance: inflow = the mass DENSITY moved through the level-0 wall faces (OPEN vents, MASS_FLUX/HRRPUA walls), summed over the steps.
    // Gate (CEO ruling): open case (through-flow > 1e-10 of the initial mass) |change - inflow| <= 1e-6 x cumulative through-flow; closed case <= 1e-10 x initial mass.
    // Open case: limit 1e-6 x (initial mass + cumulative through-flow), the mass that has been in the domain (1e-6 x through-flow alone was ~1e5 times tighter than intended early on).
    // Species: the reaction source sum_n C(n) * Q is added to the inflow (COMBUSTION's M_DOT_PPP of step s enters RHO*Z in DENSITY of step s+1 with that step's DT).
    std::vector<double> inflow(1 + ns, 0.0), r(2 + ns), reac(ns, 0.0), rc(std::max(ns, 1), 0.0);
    double through = 0.0;
    // particle mass that became gas: DENSITY of step s takes the particle sources of step s-1, so the gas change of step s matches the particle loss up to step s-1
    const double pm0 = loop.particle_mass();
    double pm_prev = pm0;
    auto mass_limit = [&] { return through > 1e-10 * m0[0] ? 1e-6 * (m0[0] + through) : 1e-10 * m0[0]; };
    const bool have_rc = ns > 0 && fds_p_reac_coef(ns, rc.data()) != 0;
    double hrr_prev = composite_hrr(reg, nl), hrr_sum = 0.0, hrr_last = 0.0;
    double t_end = 0.0, dtfill = 0.0, dtmin = 0.0, dt0 = 0.0; int nzone = 0;
    fds_k_consts(&t_end, &dtfill, &dtmin, &dt0, &nzone);
    std::string logname = o.outdir + "/" + o.chid + "_multi_level.csv";
    std::FILE* lf = amrex::ParallelDescriptor::IOProcessor() ? std::fopen(logname.c_str(), "w") : nullptr;
    if (lf) { std::fprintf(lf, "step,t,dt,mass,mass_rel_change,divmax"); for (int l = 0; l < nl; ++l) std::fprintf(lf, ",div_l%d", l); std::fprintf(lf, ",removed_rel,inflow%s\n", nl >= 2 ? ",mlr_l0,mlr_leaf" : ""); }
    double worst_mass = 0.0, worst_mass_ratio = 0.0, worst_div = 0.0, worst_fine = 0.0;   // worst_fine: max|div u - D| over levels >= 1
    int untracked = 0;                                               // steps with fine div > max(2*L0, 1e-10)
    std::vector<double> bylev;
    std::vector<double> fm_worst[2];
    long stale_regs = 0;
    int nl_now = nl;
    double deref_worst = 0.0, deref_div = 0.0; int deref_n = 0;   // ISO-DR gate: worst relative change of T/Z/p0 over derefining regrids, worst div u - D after them
    const int ldiag = std::getenv("FDSTL_LEAFDIAG") ? std::max(1, std::atoi(std::getenv("FDSTL_LEAFDIAG"))) : 0;
    int s_done = 0; bool t_end_hit = false;
    for (int s = 1; s <= o.steps; ++s) {
        s_done = s;
        { const double pm = loop.particle_mass(); inflow[0] += pm_prev - pm; pm_prev = pm; }   // particle mass lost up to step s-1 (gas source of this step)
        const std::vector<double> psrc = ns > 0 ? loop.particle_src(ns) : std::vector<double>();   // per-species particle source rates DENSITY of this step applies (species balance)
        if (!loop.advance()) { amrex::Print() << tag << ": STOP at step " << s << " (instability or non-finite DT)\n"; ++fails; break; }
        if (const char* fd = std::getenv("FDSTL_LEVEL_CSV")) {   // diagnostic: the finest level's fields after step s, one CSV per rank (name,comp,i,j,k,value), for cell-by-cell comparisons
            Fields& F = reg.fields(nl - 1);
            std::FILE* fp = std::fopen((std::string(fd) + "/level_s" + std::to_string(s) + "_r" + std::to_string(amrex::ParallelDescriptor::MyProc()) + ".csv").c_str(), "w");
            for (const char* n : {"RHO", "TMP", "U", "V", "W", "H", "D", "ZZ"}) {
                if (!fp || !F.has(n)) continue;
                const amrex::MultiFab& mf = F[n];
                for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
                    const auto a = mf.const_array(mfi);
                    for (int c = 0; c < mf.nComp(); ++c)
                        amrex::LoopOnCpu(std::getenv("FDSTL_LEVEL_CSV_GHOST") ? amrex::grow(mfi.validbox(), std::min(1, mf.nGrow())) & mfi.fabbox() : mfi.validbox(), [&](int i, int j, int k) { std::fprintf(fp, "%s,%d,%d,%d,%d,%.17g\n", n, c, i, j, k, a(i, j, k, c)); });
                }
            }
            if (fp) std::fclose(fp);
        }
        fds_p_bnd_rate(ns, r.data());   // global (summed over the boxes in box order on every rank)
        for (int n = 0; n <= ns; ++n) inflow[n] += r[n] * loop.steps().back().dt;
        for (int n = 0; n < ns; ++n) inflow[1 + n] += psrc[n] * loop.steps().back().dt;   // particle off-gas (incl. char-oxidation air uptake) per species
        through += r[1 + ns] * loop.steps().back().dt;
        if (have_rc) for (int n = 0; n < ns; ++n) reac[n] += rc[n] * hrr_prev * loop.steps().back().dt;
        // Mid-run adaptive regrid (Role 3 call shape): gate on AmrParams.regrid_interval; unbind top-down via finest_bound;
        // regrid_dynamic brackets LevelListener — do NOT call reg.begin_regrid/end_regrid here.
        {
            const int ri = core.amr_params().regrid_interval;
            if (ri > 0 && s % ri == 0) {
                // ISO-DR gate: snapshot of the level-0 T, Z and the background pressure; compared after the regrid when it removed a level (derefinement keeps them: owner rule)
                const int fb0 = core.finestLevel();
                const amrex::MultiFab* ft0 = reg.fields(0).has("TMP") ? &reg.fields(0)["TMP"] : nullptr;
                const amrex::MultiFab* fz0 = reg.fields(0).has("ZZ") ? &reg.fields(0)["ZZ"] : nullptr;
                std::unique_ptr<amrex::MultiFab> snap_t, snap_z;
                // Leaf field of T and Z at level-0 resolution: the level-0 cells under a finer box hold the average of the finest data (levels averaged down, finest first), not the coarse solution.
                auto leaf_snap = [&](const char* name) {
                    std::vector<std::unique_ptr<amrex::MultiFab>> S;
                    for (int l = 0; l <= fb0; ++l) { const amrex::MultiFab& f = reg.fields(l)[name]; S.emplace_back(new amrex::MultiFab(f.boxArray(), f.DistributionMap(), f.nComp(), 0)); amrex::MultiFab::Copy(*S[l], f, 0, 0, f.nComp(), 0); }
                    for (int l = fb0; l >= 1; --l) amrex::average_down(*S[l], *S[l - 1], 0, S[l]->nComp(), reg.level(l).ref_ratio_from_parent);
                    return std::move(S[0]);
                };
                if (ft0) snap_t = leaf_snap("TMP");
                std::unique_ptr<amrex::MultiFab> snap_r; if (reg.fields(0).has("RHO")) snap_r = leaf_snap("RHO");
                if (fz0) {   // Z is restricted as the regrid does it: average of rho*Z over average of rho
                    const int nz = fz0->nComp(); std::vector<std::unique_ptr<amrex::MultiFab>> R, Q;
                    for (int l = 0; l <= fb0; ++l) {
                        const amrex::MultiFab& z = reg.fields(l)["ZZ"]; const amrex::MultiFab& r = reg.fields(l)["RHO"];
                        R.emplace_back(new amrex::MultiFab(z.boxArray(), z.DistributionMap(), 1, 0)); amrex::MultiFab::Copy(*R[l], r, 0, 0, 1, 0);
                        Q.emplace_back(new amrex::MultiFab(z.boxArray(), z.DistributionMap(), nz, 0)); amrex::MultiFab::Copy(*Q[l], z, 0, 0, nz, 0);
                        for (int n = 0; n < nz; ++n) amrex::MultiFab::Multiply(*Q[l], *R[l], 0, n, 1, 0);
                    }
                    for (int l = fb0; l >= 1; --l) { amrex::average_down(*Q[l], *Q[l - 1], 0, nz, reg.level(l).ref_ratio_from_parent); amrex::average_down(*R[l], *R[l - 1], 0, 1, reg.level(l).ref_ratio_from_parent); }
                    for (int n = 0; n < nz; ++n) amrex::MultiFab::Divide(*Q[0], *R[0], 0, n, 1, 0);
                    snap_z = std::move(Q[0]);
                }
                std::vector<std::vector<double>> snap_p(reg.level(0).ba.size());
                for (int i = 0; i < static_cast<int>(snap_p.size()); ++i) if (reg.level(0).dm[i] == amrex::ParallelDescriptor::MyProc()) { snap_p[i].assign(4096, 0.0); int n = 0; fds_p_pbar(i + 1, 4096, snap_p[i].data(), &n); snap_p[i].resize(n); }
                std::vector<amrex::BoxArray> old_ba;
                if (ldiag) { leaf_diag(reg, "pre-regrid", s, loop.time(), ns); leaf_boxes(reg, "pre-regrid", s); for (int l = 0; l < reg.num_levels(); ++l) old_ba.push_back(reg.level(l).ba); }
                const int fb = loop.finest_bound();
                // leaf fuel rate and burning area before / after the regrid (re-snapped burner faces change the area, so the fuel rate HRRPUA*A changes with it)
                auto leaf_fuel = [&](double* a) { double r[6]; fds_fine_fuel_mdot(1, r); amrex::ParallelAllReduce::Sum(r, 6, amrex::ParallelContext::CommunicatorSub()); a[0] = r[0] - r[1] + r[2]; a[1] = r[3] - r[4] + r[5]; };
                double fuel0[2], fuel1[2]; leaf_fuel(fuel0);
                loop.part_stash();   // W8: migrate, deposit, prime, project, deposit, prime
                for (int l = fb; l >= 1; --l) loop.unbind_level(l);
                if (use_tag_vort) fill_undivided_vorticity(reg, reg.num_levels());
                if (use_tag_phi) fill_phi(reg, reg.num_levels());
                // auto: thr stays 0.45*max|VORT|_L0 of now (AmrParams thr fixed, so rescale VORT by thr/(0.45*max))
                if (o.tag_vort_auto) { const double sc = tag_vort_thr / (0.45 * reg.fields(0)["VORT"].max(0, 0)); for (int l = 0; l < reg.num_levels(); ++l) reg.fields(l)["VORT"].mult(sc); }
                const bool changed = core.regrid_dynamic(loop.time());
                const int finest = core.finestLevel();
                bool rebound = false;
                for (int l = 1; l <= finest; ++l) { loop.bind_level(l); rebound = true; }
                if (finest >= 1) {
                    fdsrt::install_cf_ghost_hooks(loop, finest, fdsrt::ThermoProvider(), &cfs);
                    loop.set_flux_overwrite(o.overwrite);
                }
                loop.part_place();
                loop.part_deposit();   // before the prime: the projection must see the particle divergence
                if (changed || rebound) loop.prime_levels();
                if (finest > fb0) {   // a new level: the hook projected onto the D it inherited from its parent; its own D comes from the prime, so project again and prime (as the initial setup does)
                    std::vector<fdsrt::ProjectionLevel> lv = projection_levels(reg, finest + 1);
                    pst.before(lv);
                    const fdsrt::ProjectionReport pr = fdsrt::project_after_regrid(lv, &solver, pst.options);
                    if (pr.ran && pr.solved) {
                        for (int l = 0; l <= finest; ++l) {
                            Fields& F = reg.fields(l);
                            const char* st[3] = {"US", "VS", "WS"}; const char* vn[3] = {"U", "V", "W"};
                            for (int d = 0; d < 3; ++d) { amrex::MultiFab::Copy(F[st[d]], F[vn[d]], 0, 0, 1, 0); F.fill_ghosts(st[d]); F.fill_ghosts(vn[d]); }
                        }
                        loop.part_deposit();   // the prime rebuilds D after the projection: the particle sources must be in place for it
                        loop.prime_levels();
                    }
                }
                if (ldiag) { leaf_diag(reg, "post-regrid", s, loop.time(), ns); leaf_boxes(reg, "post-regrid", s); old_ba.resize(std::max<std::size_t>(old_ba.size(), reg.num_levels())); leaf_diag(reg, "new-cells", s, loop.time(), ns, &old_ba); }
                if (finest < fb0) {   // derefinement: relative change of T, each Z and the background pressure at level 0, and div u - D of the remaining levels
                    amrex::BoxArray cov_new; if (finest >= 1) { cov_new = reg.level(1).ba; cov_new.coarsen(reg.level(1).ref_ratio_from_parent); }   // still covered: not a leaf cell at level 0
                    auto relmax = [&](const amrex::MultiFab& a, const amrex::MultiFab& b0, int nc) {
                        std::vector<double> w(nc, 0.0);
                        for (amrex::MFIter mfi(a); mfi.isValid(); ++mfi) {
                            auto x = a.const_array(mfi); auto y = b0.const_array(mfi);
                            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { if (cov_new.contains(amrex::IntVect(i, j, k))) return; for (int n = 0; n < nc; ++n) w[n] = std::max(w[n], std::abs(x(i, j, k, n) - y(i, j, k, n)) / std::max(std::abs(y(i, j, k, n)), 1e-300)); });
                        }
                        amrex::ParallelAllReduce::Max(w.data(), nc, amrex::ParallelContext::CommunicatorSub());
                        return w;
                    };
                    const double dT = snap_t ? relmax(reg.fields(0)["TMP"], *snap_t, 1)[0] : 0.0;
                    const double dR = snap_r ? relmax(reg.fields(0)["RHO"], *snap_r, 1)[0] : 0.0;
                    std::vector<double> dZ = snap_z ? relmax(reg.fields(0)["ZZ"], *snap_z, snap_z->nComp()) : std::vector<double>();
                    double dP = 0.0;
                    for (int i = 0; i < static_cast<int>(snap_p.size()); ++i) if (!snap_p[i].empty()) {
                        std::vector<double> pn(4096, 0.0); int n = 0; fds_p_pbar(i + 1, 4096, pn.data(), &n);
                        for (int q = 0; q < n && q < static_cast<int>(snap_p[i].size()); ++q) dP = std::max(dP, std::abs(pn[q] - snap_p[i][q]) / std::max(std::abs(snap_p[i][q]), 1e-300));
                    }
                    amrex::ParallelAllReduce::Max(dP, amrex::ParallelContext::CommunicatorSub());
                    double dD = 0.0;
                    { std::vector<fdsrt::ProjectionLevel> plv = projection_levels(reg, finest + 1); dD = fdsrt::composite_divergence_error(plv); }
                    double dZmax = 0.0; for (double z : dZ) dZmax = std::max(dZmax, z);
                    deref_worst = std::max({deref_worst, dT, dR, dZmax, dP}); deref_div = std::max(deref_div, dD); ++deref_n;
                    amrex::Print() << tag << " DEREFINE after step " << s << ": relative change T " << dT << ", RHO " << dR << ", Z";
                    for (std::size_t n = 0; n < dZ.size(); ++n) amrex::Print() << " " << dZ[n];
                    amrex::Print() << ", background pressure " << dP << " (each <= 1e-12), max|div u - D| " << dD << " (<= 1e-9)\n";
                }
                leaf_fuel(fuel1);
                nl_now = finest + 1;
                amrex::Print() << tag << " mid-regrid after step " << s << " t " << loop.time()
                               << " changed " << changed << " finest " << finest
                               << " boxes L1 " << (finest >= 1 ? reg.level(1).ba.size() : 0) << " cells L1 " << (finest >= 1 ? reg.level(1).ba.numPts() : 0)
                               << " leaf fuel rate " << std::setprecision(10) << fuel0[0] << " -> " << fuel1[0] << " kg/s, burning area " << fuel0[1] << " -> " << fuel1[1] << " m2\n";
                // Mid-run dump shares the same write helper (Architect: regrid dump path); only without a --plot-every cadence (else it dumped every regrid).
                if (!o.plotfile.empty() && o.plot_every <= 0) {
                    const std::string plt = o.plotfile + "/plt_" + o.chid + "_s" + std::to_string(s);
                    write_rho_vel_plotfile(plt, reg, loop, nl_now);
                }
            }
        }
        // Owner-show cadence (D-089 pack): --plot-every N with --plotfile; dump when s % N == 0 (1-based step). Final dump unchanged; mid-regrid dump off.
        if (o.plot_every > 0 && !o.plotfile.empty() && (s % o.plot_every) == 0) {
            const std::string plt = o.plotfile + "/plt_" + o.chid + "_s" + std::to_string(s);
            write_rho_vel_plotfile(plt, reg, loop, nl_now);
        }
        // coarse-fine flux matching on the first steps after a regrid (prime + 2 steps): per kind and scalar the coarse-minus-fine-sum mismatch the override removed,
        // and read-out arrays / runner grids not of the current hierarchy (stale flux registers; gate: 0)
        {
            const int ri = core.amr_params().regrid_interval;
            const fdsrt::StageStats* fm = nl_now > 1 ? loop.flux_match() : nullptr;
            if (fm && ri > 0 && s % ri == std::min(2, ri - 1)) {
                for (int k = 0; k < 2; ++k) {
                    std::vector<double> p = fm->pre[k];
                    int np_ = static_cast<int>(p.size());   // a rank without fine boxes has an empty vector: agree on the length, pad with 0 (the values are >= 0), then one collective
                    amrex::ParallelAllReduce::Max(np_, amrex::ParallelContext::CommunicatorSub());
                    p.resize(np_, 0.0);
                    if (np_ > 0) amrex::ParallelAllReduce::Max(p.data(), np_, amrex::ParallelContext::CommunicatorSub());
                    long st = fm->stale[k];
                    amrex::ParallelAllReduce::Sum(st, amrex::ParallelContext::CommunicatorSub());
                    stale_regs += st;
                    if (fm_worst[k].size() < p.size()) fm_worst[k].resize(p.size(), 0.0);
                    amrex::Print() << tag << " FLUXMATCH step " << s << (k ? " DIF" : " ADV") << " stale " << st << " coarse-fine before override:";
                    for (std::size_t n = 0; n < p.size(); ++n) { fm_worst[k][n] = std::max(fm_worst[k][n], p[n]); amrex::Print() << " " << (n + 1 == p.size() ? "rho" : "Z" + std::to_string(n + 1)) << " " << p[n]; }
                    amrex::Print() << "\n";
                }
                loop.reset_flux_match();
            }
        }
        const std::vector<double> m = composite_mass(reg, nl_now, ns);
        { const auto zr = leaf_z_range(reg, nl_now, ns); zrange_min = std::min(zrange_min, zr.first); zrange_max = std::max(zrange_max, zr.second); }
        const double dm = (m[0] - m0[0] - inflow[0]) / m0[0];   // mass balance residual: change minus boundary inflow
        std::vector<fdsrt::ProjectionLevel> lv = projection_levels(reg, nl_now);
        // Solid cells are left out of the check: a burner face gives its solid cell div u = U_wall/dz, not a gas cell's D (and a level-0 OBST cell that no fine level covers is not solved: masked).
        std::vector<std::unique_ptr<amrex::iMultiFab>> skip(nl_now);
        for (int l = 0; l < nl_now; ++l) {
            const amrex::iMultiFab& sm = reg.side_data(l).mask();
            skip[l].reset(new amrex::iMultiFab(sm.boxArray(), sm.DistributionMap(), 1, 0));
            for (amrex::MFIter mfi(*skip[l]); mfi.isValid(); ++mfi) {
                auto o = skip[l]->array(mfi); auto s0 = sm.const_array(mfi);
                amrex::Array4<const int> c; if (lv[l].covered) c = lv[l].covered->const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { o(i, j, k) = (s0(i, j, k) != 0 || (lv[l].covered && c(i, j, k) != 0)) ? 1 : 0; });
            }
            lv[l].covered = skip[l].get();
        }
        const double dv = fdsrt::composite_divergence_error(lv, &bylev);
        loop.pressure_csv_div(s, bylev);   // P2 solver-study CSV (no-op without the flags)
        worst_mass = std::max(worst_mass, std::abs(dm));
        worst_mass_ratio = std::max(worst_mass_ratio, std::abs(m[0] - m0[0] - inflow[0]) / mass_limit());
        worst_div = std::max(worst_div, dv);
        if (bylev.size() > 1) {
            const double f = *std::max_element(bylev.begin() + 1, bylev.end());
            worst_fine = std::max(worst_fine, f);
            untracked += f > std::max(2.0 * bylev[0], 1e-10);
        }
        const double dts = loop.steps().empty() ? 0.0 : loop.steps().back().dt;
        // CLIP_MASS / CLIP_CELLS (P3-C11/C12): every dump step when P3-C01 CSV path is active (emit zeros too).
        double fm[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};   // leaf fuel MLR (diagnostic): collective, so every rank takes part (lf is open on rank 0 only)
        if (nl >= 2) { fds_fine_fuel_mdot(1, fm); amrex::ParallelAllReduce::Sum(fm, 6, amrex::ParallelContext::CommunicatorSub()); }
        if (lf) {
            const int ntr = l0.dom.n_tracked;
            for (int l = 0; l < nl_now; ++l)
                for (int n = 1; n <= ntr; ++n) {
                    amrex::Print() << "CLIP_MASS N=" << n << " L=" << l << " mass=" << std::setprecision(17) << loop.clip_mass(n, l) << "\n";
                    amrex::Print() << "CLIP_CELLS N=" << n << " L=" << l << " count=" << loop.clip_cells(n, l) << "\n";
                }
            std::fprintf(lf, "%d,%.9g,%.9g,%.17g,%.6e,%.6e", s, loop.time(), dts, m[0], dm, dv);
            for (int l = 0; l < nl_now; ++l) std::fprintf(lf, ",%.6e", l < static_cast<int>(bylev.size()) ? bylev[l] : 0.0);
            std::fprintf(lf, ",%.3e,%.9e", loop.composite_report().removed_rel, inflow[0]);
            if (nl >= 2) std::fprintf(lf, ",%.9e,%.9e", fm[0], fm[0] - fm[1] + fm[2]);
            std::fprintf(lf, "\n");
        }
        if (o.log_every > 0 && (s % o.log_every == 0 || s == o.steps)) {
            amrex::Print() << tag << " step " << s << " t " << loop.time() << " dt " << dts << " composite mass change " << dm << " max|div u - D| " << dv;
            for (int l = 0; l < nl_now; ++l) amrex::Print() << " L" << l << " " << (l < static_cast<int>(bylev.size()) ? bylev[l] : 0.0);
            amrex::Print() << "\n";
        }
        if (ldiag && s % ldiag == 0) leaf_diag(reg, "step", s, loop.time(), ns);
        hrr_prev = composite_hrr(reg, nl_now);
        hrr_sum += hrr_prev * loop.steps().back().dt; hrr_last = hrr_prev;
        if (o.log_every > 0 && (s % o.log_every == 0 || s == o.steps)) amrex::Print() << tag << " HRR at step " << s << " composite " << hrr_prev / 1000.0 << " kW\n";
        if (loop.time() >= t_end) { amrex::Print() << tag << ": T_END " << t_end << " reached at step " << s << " t " << loop.time() << "\n"; t_end_hit = true; break; }
    }
    if (!t_end_hit && fails == 0 && s_done == o.steps) amrex::Print() << tag << ": STEP CAP " << o.steps << " reached before T_END (t " << loop.time() << " < T_END " << t_end << ")\n";
    if (lf) std::fclose(lf);
    loop.write_perf();
    if (!o.plotfile.empty() && fails == 0) {
        const std::string plt = o.plotfile + "/plt_" + o.chid + "_final";
        write_rho_vel_plotfile(plt, reg, loop, loop.num_levels());
    }
    const std::vector<double> m1 = composite_mass(reg, loop.num_levels(), ns);
    const fdsamr::TimeLoop::CompositeReport& cr = loop.composite_report();
    amrex::Print() << tag << " RESULT steps " << loop.steps().size() << " composite mass " << m0[0] << " -> " << m1[0] << " boundary inflow " << inflow[0] << " (through-flow " << through << ") balance residual / initial mass " << (m1[0] - m0[0] - inflow[0]) / m0[0] << " (worst over steps " << worst_mass << "), residual / limit worst " << worst_mass_ratio << "\n";
    // every species, also those starting at zero: |change - inflow| against the same limit as the mass
    double worst_species = 0.0;
    for (int n = 0; n < ns; ++n) {
        const double res = m1[1 + n] - m0[1 + n] - inflow[1 + n] - reac[n] - loop.clip_species(n);   // the D-031 clip moves rho*Z between species (FDS renormalisation)
        worst_species = std::max(worst_species, std::abs(res) / mass_limit());
        amrex::Print() << tag << " RESULT composite rho*Z" << n + 1 << " " << m0[1 + n] << " -> " << m1[1 + n] << " boundary inflow " << inflow[1 + n] << " reaction source " << reac[n] << " clip " << loop.clip_species(n) << " balance residual " << res << " (/ limit " << std::abs(res) / mass_limit() << ")\n";
    }
    amrex::Print() << tag << " RESULT composite HRR final " << hrr_last / 1000.0 << " kW, time mean " << (loop.time() > 0.0 ? hrr_sum / loop.time() / 1000.0 : 0.0) << " kW, t " << loop.time() << " (T_END " << t_end << ")\n";
    amrex::Print() << tag << " RESULT max|div u - D| over uncovered cells, worst over steps " << worst_div << "\n";
    // PASS gate (owner ruling): worst fine div <= 1e-9 (Poisson residual tol); the L0-tracking step count is info only; |mass rel| <= 1e-10
    for (int k = 0; k < 2; ++k) {
        amrex::Print() << tag << " RESULT flux matching " << (k ? "DIF" : "ADV") << " worst coarse-fine mismatch before override:";
        for (std::size_t n = 0; n < fm_worst[k].size(); ++n) amrex::Print() << " " << (n + 1 == fm_worst[k].size() ? "rho" : "Z" + std::to_string(n + 1)) << " " << fm_worst[k][n];
        amrex::Print() << "; after override 0 by construction (hook writes the list values)\n";
    }
    amrex::Print() << tag << " RESULT flux registers: stale (old-hierarchy) read-out arrays " << stale_regs << " (gate: 0)\n";
    const bool deref_ok = deref_n == 0 || (deref_worst <= 1e-12 && deref_div <= 1e-9);
    if (deref_n > 0) amrex::Print() << tag << " RESULT derefinement invariance (" << deref_n << " derefining regrids): worst relative change of T, Z, background pressure " << deref_worst << " (<= 1e-12), worst max|div u - D| after " << deref_div << " (<= 1e-9): " << (deref_ok ? "PASS" : "FAIL") << "\n";
    bool finite = std::isfinite(hrr_last) && std::isfinite(loop.time());
    for (double x : m1) finite = finite && std::isfinite(x);
    if (!finite) amrex::Print() << tag << " RESULT non-finite composite mass, HRR or time: FAIL\n";
    double worst_pred = 0.0; for (const StepRecord& r : loop.steps()) worst_pred = std::max(worst_pred, r.pred_div_fine);
    {   // ruling 2: fine SOLID faces without a level-0 wall (no solid takes their flux): count and energy, to stdout and the .out
        double orp[6]; fds_fine_orphan_sums(orp); amrex::ParallelDescriptor::ReduceRealSum(orp, 6);
        amrex::Print() << tag << " RESULT fine walls without a level-0 wall: " << orp[0] << " (largest per-box stage count, summed over ranks), energy " << std::setprecision(17) << orp[1] << " J\n";
        if (orp[2] > 0) amrex::Print() << tag << " RESULT WARNING fine walls remapped to the nearest level-0 wall of the same surface (interim): " << orp[2] << " (largest per-box stage count, summed over ranks), energy " << orp[3]
                                       << " J, level-0 ground convective energy " << orp[4] << " J (of |flux| " << orp[5] << " J), |remapped| / |ground| " << (orp[4] != 0 ? std::abs(orp[3] / orp[4]) : 0.0)
                                       << ", |remapped| / ground |flux| " << (orp[5] != 0 ? std::abs(orp[3]) / orp[5] : 0.0) << " (< 0.01)\n";
        fds_fine_orphan_out(orp);
    }
    amrex::Print() << tag << " RESULT worst fine-level predictor max|div US - DS| " << worst_pred << " (<= 1e-9)\n";
    const bool zok = zrange_min >= -1e-10 && zrange_max <= 1.0 + 1e-10;
    amrex::Print() << tag << " RESULT species check 7 (leaf 0 <= Z <= 1 within 1e-10, every step): min " << zrange_min << " max " << zrange_max << (zok ? " PASS" : " FAIL") << "\n";
    const bool gate = finite && worst_fine <= 1e-9 && worst_pred <= 1e-9 && worst_mass_ratio <= 1.0 && worst_species <= 1.0 && stale_regs == 0 && deref_ok && zok;
    amrex::Print() << tag << " RESULT gate " << (gate ? "PASS" : "FAIL") << ": worst fine-level div " << worst_fine << " (<= 1e-9), steps with fine div > max(2*L0, 1e-10) " << untracked << " (info only)"
                   << ", worst |mass balance| / limit " << worst_mass_ratio << ", final |species balance| / limit " << worst_species << " (<= 1; limit " << (through > 1e-10 * m0[0] ? "1e-6 x (initial mass + through-flow), open" : "1e-10 x initial mass, closed") << ")\n";
    amrex::Print() << tag << " RESULT composite pressure: " << cr.solves << " solves, backend " << cr.backend << ", removed mean (largest) " << cr.removed_mean << ", relative to rms " << cr.removed_rel
                   << ", fine H ghost cells set " << cr.h_ghost_cells << "\n";
    if (cr.res_checked > 0) amrex::Print() << tag << " RESULT pressure residual check: " << cr.res_checked << " solves checked, largest residual " << cr.res_check_max << ", largest limit " << cr.res_limit_max << " (round-off floor " << cr.res_floor_max << "), " << cr.res_failed << " above the limit\n";
    amrex::Print() << tag << " RESULT cf scalar ghost hook: calls " << cfs.calls << ", fine ghost cells " << cfs.fine_ghost_cells << ", covered cells " << cfs.covered_cells << ", conflicts " << cfs.conflicts << "\n";
    return fails + !gate;
}

}  // namespace fdsamr
