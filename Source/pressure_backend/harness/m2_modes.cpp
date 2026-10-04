// Harness modes for the M2 items on the single-level path (declared in main.cpp):
//   mode=trigger1  FR-039 trigger points (full checks versus the cheap path), single level, FFT and MLMG.
//   mode=fftcache  FFT plan caching in PressureWorkspace: bitwise identical results, plan reuse counters, timing.
// (Further modes are added below: boundary-condition mixes, boundary data, the D-057 mapping check.)
#include "PressureIface.H"
#include "CommonLayer.H"
#include "ExactSum.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>

using namespace amrex;

namespace {

int g_m2fail = 0;
void mcheck (bool ok, std::string const& what)
{
    Print() << "CHECK " << (ok ? "PASS " : "FAIL ") << what << "\n";
    if (!ok) { ++g_m2fail; }
}

// "NN,DD,ND" -> six face types (x lo, x hi, y lo, ...: stored as face_index(dir, side)). Letters: N Neumann, D Dirichlet, P periodic.
bool parse_pairs (std::string const& s, std::array<pb::BC,6>& bc)
{
    std::vector<std::string> tok;
    std::string cur;
    for (char c : s) { if (c == ',' || c == ' ') { if (!cur.empty()) { tok.push_back(cur); cur.clear(); } } else { cur += c; } }
    if (!cur.empty()) { tok.push_back(cur); }
    if (tok.size() != 3) { return false; }
    for (int d = 0; d < 3; ++d) {
        if (tok[d].size() != 2) { return false; }
        for (int side = 0; side < 2; ++side) {
            const char c = tok[d][side];
            pb::BC b;
            if (c == 'N') { b = pb::BC::Neumann; } else if (c == 'D') { b = pb::BC::Dirichlet; } else if (c == 'P') { b = pb::BC::Periodic; } else { return false; }
            bc[pb::face_index(d, side)] = b;
        }
    }
    return true;
}

struct S1 {
    Box domain;
    Geometry geom;
    BoxArray ba;
    DistributionMapping dm;
    std::array<pb::BC,6> bc;
    Real dx = 0;
};

S1 make_s1 (Vector<int> const& n, std::string const& pairs, int mgs)
{
    S1 s;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(parse_pairs(pairs, s.bc), "bcpairs must look like NN,DD,ND (letters N D P)");
    const int nmax = std::max({n[0], n[1], n[2]});
    s.dx = Real(1.0) / nmax;
    s.domain = Box(IntVect(0), IntVect(n[0]-1, n[1]-1, n[2]-1));
    RealBox rb({0.,0.,0.}, {n[0]*s.dx, n[1]*s.dx, n[2]*s.dx});
    Array<int,3> per{s.bc[0] == pb::BC::Periodic, s.bc[1] == pb::BC::Periodic, s.bc[2] == pb::BC::Periodic};
    s.geom = Geometry(s.domain, rb, CoordSys::cartesian, per);
    s.ba = BoxArray(s.domain);
    s.ba.maxSize(mgs);
    s.dm = DistributionMapping(s.ba);
    return s;
}

Real smooth_rhs (Real x, Real y, Real z)
{
    const Real pi = Real(3.141592653589793238462643383279502884);
    const Real r2 = (x-Real(0.31))*(x-Real(0.31)) + (y-Real(0.57))*(y-Real(0.57)) + (z-Real(0.45))*(z-Real(0.45));
    return std::exp(-r2/Real(0.02)) + Real(0.5)*std::sin(Real(3.)*pi*x)*std::cos(Real(2.)*pi*y+Real(0.4))*(z+Real(0.25));
}

void fill_smooth (S1 const& s, MultiFab& rhs, double offset)
{
    for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
        auto const& a = rhs.array(mfi);
        const Real dx = s.dx;
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { a(i,j,k) = smooth_rhs((i+0.5)*dx, (j+0.5)*dx, (k+0.5)*dx) + offset; });
    }
}

std::vector<double> gather_s (MultiFab const& mf, Box const& domain)
{
    BoxArray ba1(domain);
    DistributionMapping dm1(Vector<int>{ParallelDescriptor::IOProcessorNumber()});
    MultiFab g(ba1, dm1, 1, 0);
    g.ParallelCopy(mf, 0, 0, 1);
    std::vector<double> v;
    if (ParallelDescriptor::IOProcessor()) { v.assign(g[0].dataPtr(), g[0].dataPtr() + domain.numPts()); }
    return v;
}

pb::PressureProblem problem_of (S1 const& s, MultiFab& rhs, MultiFab& phi)
{
    pb::PressureProblem p;
    p.ba = s.ba; p.dm = s.dm; p.geom = s.geom; p.bc = s.bc; p.rhs = &rhs; p.phi = &phi;
    return p;
}

using Clock = std::chrono::steady_clock;
double secs (Clock::time_point a, Clock::time_point b) { return std::chrono::duration<double>(b - a).count(); }

// --------------------------------------------------------------------------------------------------------------
// mode=trigger1
// --------------------------------------------------------------------------------------------------------------
void run_trigger1 (ParmParse& pp)
{
    Vector<int> n{16, 16, 16}; pp.queryarr("n_cell", n);
    int mgs = 8; pp.query("mgs", mgs);
    std::string bcs = "NN,NN,NN"; pp.query("bcpairs", bcs);
    S1 s = make_s1(n, bcs, mgs);
    MultiFab rhs(s.ba, s.dm, 1, 0), phi(s.ba, s.dm, 1, 1);
    fill_smooth(s, rhs, 0.3);                  // incompatible: nonzero mean, so the compatibility diagnostic has something to see
    pb::PressureProblem p = problem_of(s, rhs, phi);
    for (std::string const name : {"fft", "mlmg"}) {
        pb::PressureOptions o; o.verbose = 0; o.removed_mean_warn = 1.0e-3;
        o.backend = (name == "fft") ? pb::BackendKind::FFT : pb::BackendKind::MLMG;
        const std::string tag = name + ": ";
        auto run = [&] (pb::PressureOptions const& oo, pb::PressureWorkspace* ws, std::vector<double>& out) {
            phi.setVal(0.0);
            pb::PressureResult r = pb::solve_pressure(p, oo, ws);
            out = gather_s(phi, s.domain);
            return r;
        };
        std::vector<double> ref, v;
        // defaults: full checks on every solve, as before the option existed
        pb::PressureResult r0 = run(o, nullptr, ref);
        mcheck(r0.status == pb::Status::Ok && r0.full_checks && r0.residual_checked && r0.triggers == pb::SolveDebug, tag + "default options: full checks, residual checked, trigger = Debug");
        mcheck(r0.components[0].removed_rel > 1e-3 && !r0.warnings.empty(), tag + "default options: compatibility diagnostic measured and warning raised");
        // Routine solve without a workspace: cheap path
        pb::PressureOptions oc = o; oc.trigger = pb::SolveRoutine;
        pb::PressureResult r1 = run(oc, nullptr, v);
        mcheck(r1.status == pb::Status::Ok && !r1.full_checks && !r1.residual_checked && r1.warnings.empty() && r1.components[0].removed_rel == 0.0, tag + "Routine, no workspace: cheap path (no residual, no diagnostic, no warning)");
        mcheck(v == ref, tag + "cheap path solution is bitwise the solution of the full path");
        mcheck(r1.components[0].removed_mean == r0.components[0].removed_mean && r1.components[0].gauge_shift == r0.components[0].gauge_shift, tag + "cheap path reports the same removed mean and gauge constant (bitwise)");
        // workspace bookkeeping
        pb::PressureWorkspace ws;
        mcheck(ws.num_solves() == 0 && ws.auto_triggers() == pb::SolveFirst, tag + "fresh workspace: next solve is a FirstSolve");
        pb::PressureResult a = run(oc, &ws, v);
        mcheck(a.full_checks && (a.triggers & pb::SolveFirst) && a.residual_checked, tag + "workspace, first solve: full checks (FirstSolve detected by the workspace)");
        mcheck(v == ref, tag + "workspace first solve bitwise equal");
        pb::PressureResult b2 = run(oc, &ws, v);
        mcheck(!b2.full_checks && b2.triggers == pb::SolveRoutine && !b2.residual_checked, tag + "workspace, second Routine solve: cheap path");
        mcheck(v == ref, tag + "workspace cheap solve bitwise equal to the full-check solve");
        pb::PressureOptions of = oc; of.trigger = pb::SolveFirstAfterRegrid;
        pb::PressureResult c2 = run(of, &ws, v);
        mcheck(c2.full_checks && c2.residual_checked && (c2.triggers & pb::SolveFirstAfterRegrid), tag + "caller flag FirstAfterRegrid: full checks");
        pb::PressureOptions od = oc; od.trigger = pb::SolveDebug;
        mcheck(run(od, &ws, v).full_checks, tag + "caller flag Debug: full checks");
        mcheck(!run(oc, &ws, v).full_checks, tag + "back to Routine: cheap path");
        // an explicit rebuild marks a regrid
        std::string msg;
        pb::Status rs = ws.rebuild(p, &msg);
        mcheck(rs == pb::Status::Ok && ws.regrid_pending() && ws.auto_triggers() == pb::SolveFirstAfterRegrid, tag + "rebuild(): next solve is FirstAfterRegrid (" + msg + ")");
        mcheck(run(oc, &ws, v).full_checks, tag + "first solve after rebuild(): full checks");
        mcheck(!run(oc, &ws, v).full_checks && !ws.regrid_pending(), tag + "the following Routine solve is cheap again");
        // a changed layout is detected by the solve itself and counts as a regrid
        {
            S1 s2 = make_s1(n, bcs, mgs >= 8 ? mgs/2 : mgs*2);
            MultiFab rhs2(s2.ba, s2.dm, 1, 0), phi2(s2.ba, s2.dm, 1, 1);
            fill_smooth(s2, rhs2, 0.3);
            pb::PressureProblem p2 = problem_of(s2, rhs2, phi2);
            phi2.setVal(0.0);
            pb::PressureResult d = pb::solve_pressure(p2, oc, &ws);
            mcheck(name == "mlmg" || (d.full_checks && d.workspace_rebuilt && !d.backend_status.plan_reused), tag + "new box layout with the same workspace: detected, full checks, plan rebuilt");
            std::vector<double> v2 = gather_s(phi2, s2.domain);
            mcheck(v2.size() == ref.size(), tag + "layout change solve returned a field");
        }
        // mask of triggers
        pb::PressureWorkspace ws2;
        pb::PressureOptions om = oc; om.full_checks_on = pb::SolveDebug;          // only Debug solves are checked
        mcheck(!run(om, &ws2, v).full_checks, tag + "full_checks_on = Debug only: FirstSolve alone does not trigger the full checks");
        om.trigger = pb::SolveDebug;
        mcheck(run(om, &ws2, v).full_checks, tag + "full_checks_on = Debug only: Debug does");
        pb::PressureOptions on = oc; on.full_checks_on = 0u;
        mcheck(!run(on, &ws2, v).full_checks, tag + "full_checks_on = 0: never");
        // cheap path keeps the error handling
        {
            pb::PressureProblem bad = p; bad.bc[0] = pb::BC::Dirichlet;     // mixed faces on a backend/problem that may be unsupported
            pb::PressureProblem pc = p; pc.cylindrical = true;
            phi.setVal(7.0);
            pb::PressureResult e = pb::solve_pressure(pc, oc, &ws2);
            mcheck(e.status == pb::Status::NotBuilt && !e.message.empty() && phi.min(0) == 7.0, tag + "cheap path: NotBuilt is still reported, phi untouched");
            pb::PressureOptions ot = oc; ot.max_iter = 1; ot.tol_rel = 1e-14;
            if (name == "mlmg") {
                pb::PressureResult nc = run(ot, nullptr, v);
                mcheck(nc.status == pb::Status::NotConverged, tag + "cheap path: NotConverged is still reported");
            }
        }
        Print() << "TRIGGER " << name << " ok\n";
    }
}

// --------------------------------------------------------------------------------------------------------------
// mode=fftcache: repeated FFT solves on one layout with and without the workspace
// --------------------------------------------------------------------------------------------------------------
void run_fftcache (ParmParse& pp)
{
    Vector<int> n{64, 64, 64}; pp.queryarr("n_cell", n);
    int mgs = 32; pp.query("mgs", mgs);
    int nsolve = 10; pp.query("nsolve", nsolve);
    std::string bcs = "NN,NN,NN"; pp.query("bcpairs", bcs);
    S1 s = make_s1(n, bcs, mgs);
    MultiFab rhs(s.ba, s.dm, 1, 0), phi(s.ba, s.dm, 1, 1);
    fill_smooth(s, rhs, 0.0);
    pb::PressureProblem p = problem_of(s, rhs, phi);
    pb::PressureOptions o; o.backend = pb::BackendKind::FFT; o.verbose = 0; o.trigger = pb::SolveRoutine;
    auto solve_n = [&] (pb::PressureWorkspace* ws, std::vector<double>& last, bool& all_reused_after_first) {
        std::vector<double> first;
        all_reused_after_first = true;
        ParallelDescriptor::Barrier();
        const auto t0 = Clock::now();
        for (int it = 0; it < nsolve; ++it) {
            phi.setVal(0.0);
            pb::PressureResult r = pb::solve_pressure(p, o, ws);
            if (r.status != pb::Status::Ok) { mcheck(false, "solve status Ok"); }
            if (ws && it > 0 && !r.backend_status.plan_reused) { all_reused_after_first = false; }
        }
        ParallelDescriptor::Barrier();
        const double t = secs(t0, Clock::now());
        last = gather_s(phi, s.domain);
        double tm = t; ParallelDescriptor::ReduceRealMax(tm);
        return tm;
    };
    std::vector<double> v_fresh, v_ws, v_ws2;
    bool re1 = true, re2 = true;
    solve_n(nullptr, v_fresh, re1);                           // warm-up (AMReX caches, FFT libraries)
    const double t_fresh = solve_n(nullptr, v_fresh, re1);
    pb::PressureWorkspace ws;
    const double t_ws = solve_n(&ws, v_ws, re1);              // includes the one plan build
    const double t_ws2 = solve_n(&ws, v_ws2, re2);            // all reuse
    mcheck(v_fresh == v_ws && v_ws == v_ws2, "workspace solutions are bitwise identical to the fresh-plan solutions");
    mcheck(re1 && re2, "every solve after the first reused the cached plan");
    mcheck(ws.fft_plan_builds() == 1, "exactly one plan was built for " + std::to_string(2*nsolve) + " solves (builds " + std::to_string(ws.fft_plan_builds()) + ", reuses " + std::to_string(ws.fft_plan_reuses()) + ")");
    // invalidate on rebuild
    mcheck(ws.matches(p) && ws.fft_plan_built(), "workspace matches the layout and holds a plan");
    ParallelDescriptor::Barrier();
    const auto tb0 = Clock::now();
    const pb::Status rbs = ws.rebuild(p);
    ParallelDescriptor::Barrier();
    double t_build = secs(tb0, Clock::now()); ParallelDescriptor::ReduceRealMax(t_build);
    mcheck(rbs == pb::Status::Ok && ws.fft_plan_builds() == 2, "rebuild() discards and rebuilds the plan");
    {   // another layout / another BC invalidate
        S1 s2 = make_s1(n, bcs, std::max(8, mgs/2));
        MultiFab r2(s2.ba, s2.dm, 1, 0), p2(s2.ba, s2.dm, 1, 1); fill_smooth(s2, r2, 0.0); p2.setVal(0.0);
        pb::PressureProblem q = problem_of(s2, r2, p2);
        mcheck(!ws.matches(q), "plan key differs for another box layout");
        pb::PressureResult r = pb::solve_pressure(q, o, &ws);
        mcheck(r.status == pb::Status::Ok && !r.backend_status.plan_reused && ws.fft_plan_builds() == 3, "a solve on another layout rebuilds the plan");
        std::vector<double> v2 = gather_s(p2, s2.domain);
        mcheck(v2 == v_fresh, "result on the other box layout equals the first within bitwise (same FFT, same data)" );
    }
    {
        S1 s3 = make_s1(n, (bcs == "PP,PP,PP") ? "NN,NN,NN" : "PP,PP,PP", mgs);
        MultiFab r3(s3.ba, s3.dm, 1, 0), p3(s3.ba, s3.dm, 1, 1); fill_smooth(s3, r3, 0.0); p3.setVal(0.0);
        pb::PressureProblem q = problem_of(s3, r3, p3);
        mcheck(!ws.matches(q), "plan key differs for other boundary types");
    }
    Print() << std::setprecision(6) << "FFTCACHE n=" << n[0] << "x" << n[1] << "x" << n[2] << " nranks=" << ParallelDescriptor::NProcs()
            << " mgs=" << mgs << " nsolve=" << nsolve << " fresh_s=" << t_fresh << " ws_first_s=" << t_ws << " ws_reuse_s=" << t_ws2
            << " per_solve_fresh_ms=" << 1e3*t_fresh/nsolve << " per_solve_reuse_ms=" << 1e3*t_ws2/nsolve
            << " plan_build_ms=" << 1e3*t_build << " speedup=" << t_fresh/t_ws2 << " plan_builds=" << ws.fft_plan_builds() << " plan_reuses=" << ws.fft_plan_reuses() << "\n";
}

// --------------------------------------------------------------------------------------------------------------
// mode=mixed1: one solve with the given face types, rhs and solution dumped for the dense numpy reference.
//   keys: n_cell, bcpairs (NN,DD,ND ...), mgs, backend=fft|mlmg|auto, out=<prefix>
//   prints MIXED status=<> backend=<> msg=<> and the library's own true-residual rel2.
// --------------------------------------------------------------------------------------------------------------
Real hash_noise (int i, int j, int k)
{
    std::uint64_t h = 1469598103934665603ULL;
    for (int v : {i, j, k}) { h ^= std::uint64_t(v + 1000); h *= 1099511628211ULL; h ^= h >> 29; }
    return Real(double(h % 2000003ULL) / 1000001.5 - 1.0);
}

void write_field (std::string const& fname, std::vector<double> const& v)
{
    if (!ParallelDescriptor::IOProcessor()) { return; }
    std::ofstream f(fname, std::ios::binary);
    f.write(reinterpret_cast<const char*>(v.data()), std::streamsize(v.size()*sizeof(double)));
}

void mixed_one (Vector<int> const& n, std::string const& bcs, int mgs, std::string const& be, std::string const& out)
{
    S1 s = make_s1(n, bcs, mgs);
    MultiFab rhs(s.ba, s.dm, 1, 0), phi(s.ba, s.dm, 1, 1);
    for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
        auto const& a = rhs.array(mfi);
        const Real dx = s.dx;
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { a(i,j,k) = smooth_rhs((i+0.5)*dx, (j+0.5)*dx, (k+0.5)*dx) + 0.5*hash_noise(i,j,k) + 0.3; });
    }
    phi.setVal(0.0);
    pb::PressureProblem p = problem_of(s, rhs, phi);
    pb::PressureOptions o; o.verbose = 0; o.removed_mean_warn = 1.0e9;   // the rhs is deliberately incompatible
    o.backend = (be == "fft") ? pb::BackendKind::FFT : (be == "mlmg") ? pb::BackendKind::MLMG : pb::BackendKind::Auto;
    o.tol_rel = 1.0e-13; o.max_iter = 100;
    pb::PressureResult r = pb::solve_pressure(p, o);
    Print() << std::setprecision(6) << "MIXED bcpairs=" << bcs << " status=" << pb::to_string(r.status) << " backend=" << (r.backend.empty() ? std::string("none") : r.backend)
            << " singular=" << (r.components.empty() ? -1 : int(r.components[0].singular))
            << " true_rel2=" << r.residual_rel2 << " msg=[" << r.message << "]\n";
    if (r.status == pb::Status::Ok && !out.empty()) {
        write_field(out + "_rhs.bin", gather_s(rhs, s.domain));
        write_field(out + "_phi.bin", gather_s(phi, s.domain));
    }
}

void run_mixed1 (ParmParse& pp)
{
    Vector<int> n{6, 5, 4}; pp.queryarr("n_cell", n);
    int mgs = 4; pp.query("mgs", mgs);
    std::string bcs = "NN,NN,NN"; pp.query("bcpairs", bcs);
    std::string be = "auto"; pp.query("backend", be);
    std::string out; pp.query("out", out);
    int sweep = 0; pp.query("sweep", sweep);
    if (!sweep) { mixed_one(n, bcs, mgs, be, out); return; }
    // all 5^3 combinations of {PP, NN, DD, ND, DN} per direction, in the order the python reference expects
    const char* L[5] = {"PP", "NN", "DD", "ND", "DN"};
    int idx = 0;
    for (int a = 0; a < 5; ++a) for (int b = 0; b < 5; ++b) for (int c = 0; c < 5; ++c) {
        mixed_one(n, std::string(L[a]) + "," + L[b] + "," + L[c], mgs, be, out + "_" + std::to_string(idx++));
    }
}

} // namespace

int run_m2_mode (std::string const& mode, ParmParse& pp)
{
    if (mode == "trigger1") { run_trigger1(pp); }
    else if (mode == "fftcache") { run_fftcache(pp); }
    else if (mode == "mixed1") { run_mixed1(pp); }
    else { return -1; }
    return g_m2fail;
}
