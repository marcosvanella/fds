// Residual-limit modes of the harness (A-67 follow-ups E1, E2, R6, E3; frozen/residual-notes.md).
//   mode=manuf    single-level manufactured problem as in the A-56 residual study: H* smooth (+ optional grid noise), b := L H* with the
//                 backend's operator and ghost rules, solve through the product API, print one RESID line (check, floor, limit, error).
//                 Options: n, mgs, bcpairs (e.g. NN,NN,NN | DD,DD,DD | PP,PP,PP), backend (fft|mlmg|hypre), noise, tol_rel, max_iter,
//                 gauge_k (constant KRES-like gauge offset = gauge_k * rms(H*)), expect_warn (-1 skip, 0 none, 1 at least one).
// Every mode prints CHECK PASS|FAIL lines and returns the number of failures.
#include "PressureIface.H"
#include "CommonLayer.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <cmath>
#include <iomanip>

using namespace amrex;

namespace {

int g_rfail = 0;
void rcheck (bool ok, std::string const& what)
{
    Print() << "CHECK " << (ok ? "PASS " : "FAIL ") << what << "\n";
    if (!ok) { ++g_rfail; }
}

double hnoise (int i, int j, int k)
{
    unsigned long long x = (unsigned long long)(i+1)*73856093ULL ^ (unsigned long long)(j+1)*19349663ULL ^ (unsigned long long)(k+1)*83492791ULL;
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33;
    return double(x >> 11)/9007199254740992.0 - 0.5;
}

pb::BC bc_of (char c)
{
    return c == 'N' ? pb::BC::Neumann : c == 'D' ? pb::BC::Dirichlet : pb::BC::Periodic;
}

pb::BackendKind backend_of (std::string const& s)
{
    return s == "fft" ? pb::BackendKind::FFT : s == "mlmg" ? pb::BackendKind::MLMG : s == "hypre" ? pb::BackendKind::HYPRE : pb::BackendKind::Auto;
}

void run_manuf (ParmParse& pp)
{
    int n = 32, mgs = 32, maxit = 200, expect_warn = -1;
    std::string bcs = "NN,NN,NN", be = "mlmg";
    double tol = 1e-12, noise = 0.0, gauge_k = 0.0;
    pp.query("n", n); pp.query("mgs", mgs); pp.query("bcpairs", bcs); pp.query("backend", be); pp.query("tol_rel", tol);
    pp.query("noise", noise); pp.query("max_iter", maxit); pp.query("gauge_k", gauge_k); pp.query("expect_warn", expect_warn);
    pb::PressureProblem p;
    bool singular = true; int per[3] = {0, 0, 0};
    for (int d = 0; d < 3; ++d) {
        const std::string t = bcs.substr(3*d, 2);
        for (int s = 0; s < 2; ++s) {
            p.bc[pb::face_index(d, s)] = bc_of(t[s]);
            if (t[s] == 'D') { singular = false; }
            if (t[s] == 'P') { per[d] = 1; }
        }
    }
    Box dom(IntVect(0), IntVect(n-1));
    RealBox rb({0.,0.,0.}, {1.,1.,1.});
    p.geom = Geometry(dom, rb, CoordSys::cartesian, Array<int,3>{per[0], per[1], per[2]});
    p.ba = BoxArray(dom); p.ba.maxSize(mgs); p.dm = DistributionMapping(p.ba);
    MultiFab ps(p.ba, p.dm, 1, 1), rhs(p.ba, p.dm, 1, 0), phi(p.ba, p.dm, 1, 1), kres(p.ba, p.dm, 1, 0);
    ps.setVal(0.0); phi.setVal(0.0);
    const double PI = 3.14159265358979323846, h = 1.0/n, ih2 = 1.0/(h*h);
    for (MFIter mfi(ps); mfi.isValid(); ++mfi) {
        auto a = ps.array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const double x = (i+.5)*h, y = (j+.5)*h, z = (k+.5)*h;
            a(i,j,k) = std::cos(PI*x)*std::cos(2*PI*y)*std::cos(PI*z) + 0.5*std::sin(1.5*PI*x)*std::sin(PI*y)*std::cos(3*PI*z) + 0.3 + noise*hnoise(i,j,k);
        });
    }
    const double sm = ps.sum(0) / double(dom.numPts());
    if (singular) { ps.plus(-sm, 0, 1, 0); }
    const double prms = ps.norm2(0) / std::sqrt(double(dom.numPts()));
    ps.FillBoundary(p.geom.periodicity());
    for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
        auto b = rhs.array(mfi); auto a = ps.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const double c = a(i,j,k); double s = 0;
            for (int d = 0; d < 3; ++d) {
                for (int sd = 0; sd < 2; ++sd) {
                    int nb[3] = {i, j, k}; nb[d] += sd ? 1 : -1;
                    double v;
                    if ((nb[d] < 0 || nb[d] >= n) && !per[d]) { v = (p.bc[pb::face_index(d, sd)] == pb::BC::Neumann) ? c : -c; }
                    else { v = a(nb[0], nb[1], nb[2]); }
                    s += v - c;
                }
            }
            b(i,j,k) = s*ih2;
        });
    }
    kres.setVal(gauge_k * prms);
    p.rhs = &rhs; p.phi = &phi;
    if (gauge_k != 0.0) { p.gauge_offset = &kres; }
    pb::PressureOptions o; o.verbose = 0; o.tol_rel = tol; o.max_iter = maxit; o.removed_mean_warn = 1e9; o.backend = backend_of(be);
    const double t0 = amrex::second();
    pb::PressureResult r = pb::solve_pressure(p, o);
    double wall = amrex::second() - t0; ParallelDescriptor::ReduceRealMax(wall);
    MultiFab e(p.ba, p.dm, 1, 0), pv(p.ba, p.dm, 1, 0);
    MultiFab::Copy(e, phi, 0, 0, 1, 0); MultiFab::Copy(pv, ps, 0, 0, 1, 0); MultiFab::Subtract(e, pv, 0, 0, 1, 0);
    if (singular) { e.plus(-e.sum(0)/double(dom.numPts()), 0, 1, 0); }
    const double err = e.norm2(0) / pv.norm2(0);
    const double epsH = std::max(1e-8, 2.4e-12*double(n)*n);
    // Residual of the RETURNED (gauge-fixed) H, evaluated independently of the product check: shows what the gauge offset does to H itself.
    pb::ResidualNorms tn = pb::true_residual(p, phi, rhs);
    Print() << std::setprecision(5) << "RESID n=" << n << " bc=" << bcs << " backend=" << be << " noise=" << noise << " gauge_k=" << gauge_k
            << " status=" << pb::to_string(r.status) << " iters=" << r.backend_status.iterations << " own=" << r.backend_status.own_residual
            << " check=" << r.residual_check << " floor=" << r.residual_floor << " limit=" << r.residual_limit << " ok=" << (r.residual_ok ? 1 : 0)
            << " nwarn=" << r.warnings.size() << " raw=" << r.residual_rel2 << " nopin=" << r.residual_rel2_nopin
            << " returned_H_rel2=" << tn.rel2 << " err_rel2=" << err << " epsH=" << epsH << " wall_s=" << wall << " msg=\"" << r.message << "\"\n";
    if (expect_warn != 1) {   // a negative control (expect_warn=1) deliberately under-solves
        rcheck(r.status == pb::Status::Ok, "manuf: status Ok");
        // A large gauge_k adds KRES-like representation noise to the returned H; the E2 test checks residual independence, not that
        // the gauged H still matches H* to eps_H. Skip the error check when an offset is applied.
        if (gauge_k == 0.0) { rcheck(err <= epsH, "manuf: solution error within eps_H"); }
    }
    if (expect_warn == 0) { rcheck(r.residual_ok && r.warnings.empty(), "manuf: no residual warning"); }
    if (expect_warn == 1) { rcheck(!r.residual_ok && !r.warnings.empty(), "manuf: residual warning raised (negative control)"); }
    (void)singular;
}

} // namespace

int run_resid_mode (std::string const& mode, ParmParse& pp)
{
    if (mode == "manuf") { run_manuf(pp); }
    else { return -1; }
    return g_rfail;
}
