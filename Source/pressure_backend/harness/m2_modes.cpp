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
#include <AMReX_FFT_Poisson.H>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <functional>
#include <memory>
#include <sstream>

using namespace amrex;

namespace {

int g_m2fail = 0;
void mcheck (bool ok, std::string const& what)
{
    Print() << "CHECK " << (ok ? "PASS " : "FAIL ") << what << "\n";
    if (!ok) { ++g_m2fail; if (!ParallelDescriptor::IOProcessor()) { amrex::AllPrint() << "CHECK FAIL (rank " << ParallelDescriptor::MyProc() << ") " << what << "\n"; } }
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
// mode=hypcache: HYPRE set-up (matrix + BoomerAMG) cached in the PressureWorkspace on the single-level path
// --------------------------------------------------------------------------------------------------------------
void run_hypcache (ParmParse& pp)
{
    Vector<int> n{32, 32, 32}; pp.queryarr("n_cell", n);
    int mgs = 16; pp.query("mgs", mgs);
    int nsolve = 5; pp.query("nsolve", nsolve);
    std::string bcs = "NN,NN,NN"; pp.query("bcpairs", bcs);
    S1 s = make_s1(n, bcs, mgs);
    MultiFab rhs(s.ba, s.dm, 1, 0), phi(s.ba, s.dm, 1, 1);
    fill_smooth(s, rhs, 0.0);
    pb::PressureProblem p = problem_of(s, rhs, phi);
    pb::PressureOptions o; o.backend = pb::BackendKind::HYPRE; o.verbose = 0; o.trigger = pb::SolveRoutine;
    std::vector<double> fresh, wsv;
    double t_fresh = 0.0, t_first = 0.0, t_reuse = 0.0;
    int iters = 0; bool all_reused = true;
    for (int it = 0; it < nsolve; ++it) {      // no workspace: set-up built and dropped every solve
        phi.setVal(0.0); const double t0 = amrex::second();
        pb::PressureResult r = pb::solve_pressure(p, o, nullptr);
        t_fresh += amrex::second() - t0;
        mcheck(r.status == pb::Status::Ok && r.backend == "HYPRE", "fresh solve status Ok on HYPRE");
        iters = r.backend_status.iterations;
        fresh = gather_s(phi, s.domain);
    }
    pb::PressureWorkspace ws;
    mcheck(!ws.hypre_built(), "empty workspace holds no HYPRE set-up");
    for (int it = 0; it < nsolve; ++it) {
        phi.setVal(0.0); const double t0 = amrex::second();
        pb::PressureResult r = pb::solve_pressure(p, o, &ws);
        const double dt = amrex::second() - t0;
        mcheck(r.status == pb::Status::Ok, "workspace solve status Ok");
        if (it == 0) { t_first = dt; mcheck(!r.backend_status.plan_reused && ws.hypre_built(), "first solve builds the set-up"); }
        else { t_reuse += dt; if (!r.backend_status.plan_reused) { all_reused = false; } }
        wsv = gather_s(phi, s.domain);
        mcheck(wsv == fresh, "workspace solution is bitwise identical to the fresh-set-up solution");
    }
    mcheck(all_reused, "every later solve reused the cached HYPRE set-up");
    mcheck(ws.matches(p) && ws.hypre_built(), "workspace matches the layout and holds a HYPRE set-up");
    mcheck(ws.rebuild(p) == pb::Status::Ok && !ws.hypre_built(), "rebuild() invalidates the cached HYPRE set-up");
    {
        phi.setVal(0.0);
        pb::PressureResult r = pb::solve_pressure(p, o, &ws);
        mcheck(r.status == pb::Status::Ok && !r.backend_status.plan_reused && ws.hypre_built(), "the next solve after rebuild() builds a new set-up");
        mcheck(gather_s(phi, s.domain) == fresh, "solution after rebuild is bitwise identical");
    }
    {   // another box layout: the key differs, set-up rebuilt for the new layout (no pointers into the old MultiFabs)
        S1 s2 = make_s1(n, bcs, std::max(8, mgs/2));
        MultiFab r2(s2.ba, s2.dm, 1, 0), p2(s2.ba, s2.dm, 1, 1); fill_smooth(s2, r2, 0.0); p2.setVal(0.0);
        pb::PressureProblem q = problem_of(s2, r2, p2);
        mcheck(!ws.matches(q), "workspace key differs for another box layout");
        pb::PressureResult r = pb::solve_pressure(q, o, &ws);
        mcheck(r.status == pb::Status::Ok && !r.backend_status.plan_reused && r.workspace_rebuilt, "a solve on another layout rebuilds the set-up and flags it");
        std::vector<double> v2 = gather_s(p2, s2.domain);
        double dn = 0.0, nn = 0.0;
        for (std::size_t q2 = 0; q2 < v2.size() && q2 < fresh.size(); ++q2) { dn += (v2[q2]-fresh[q2])*(v2[q2]-fresh[q2]); nn += fresh[q2]*fresh[q2]; }
        ParallelDescriptor::ReduceRealSum(dn); ParallelDescriptor::ReduceRealSum(nn);    // gathered on the I/O rank only
        mcheck(std::sqrt(dn/nn) <= 1e-8, "result on the other box layout agrees within eps_H (rel L2 " + std::to_string(std::sqrt(dn/nn)) + ")");
    }
    {   // options are part of the key
        pb::PressureOptions o2 = o; o2.hypre.relax_type = 6;
        pb::PressureResult r = pb::solve_pressure(p, o2, &ws);
        mcheck(r.status == pb::Status::Ok && !r.backend_status.plan_reused, "other HYPRE options rebuild the set-up");
    }
    Print() << std::setprecision(6) << "HYPCACHE n=" << n[0] << "x" << n[1] << "x" << n[2] << " nranks=" << ParallelDescriptor::NProcs() << " mgs=" << mgs << " iters=" << iters
            << " per_solve_fresh_ms=" << 1e3*t_fresh/nsolve << " first_ms=" << 1e3*t_first << " per_solve_reuse_ms=" << 1e3*t_reuse/std::max(1, nsolve-1) << "\n";
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
    o.backend = (be == "fft") ? pb::BackendKind::FFT : (be == "mlmg") ? pb::BackendKind::MLMG : (be == "hypre") ? pb::BackendKind::HYPRE : pb::BackendKind::Auto;
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

// --------------------------------------------------------------------------------------------------------------
// mode=bcdata: inhomogeneous boundary data folded into the right-hand side (fold_boundary_data).
//   part=exact: discrete test. rhs := L_full(phi) built with explicit ghost values from data, fold, solve the homogeneous
//               problem, compare with phi (up to a constant if singular).
//   part=mms:   manufactured solution u = cos(1.3x+.2) cosh(.7y) sin(.9z+.4) on [0,1]^3 with exact wall data; error vs u.
// --------------------------------------------------------------------------------------------------------------
Real bdata_fn (int f, int i, int j, int k)
{
    return (Real(0.2) + Real(0.3)*f) * hash_noise(i + 7*f, j + 3*f, k) + Real(0.2) + Real(0.05)*f;
}

// Slab MultiFab of the face-adjacent layer with another distribution than the problem.
std::unique_ptr<MultiFab> make_slab (S1 const& s, int f, std::function<Real(int,int,int)> const& fn)
{
    const int d = f % 3, side = f / 3;
    Box layer = s.domain;
    if (side == 0) { layer.setBig(d, s.domain.smallEnd(d)); } else { layer.setSmall(d, s.domain.bigEnd(d)); }
    BoxArray ba(layer); ba.maxSize(3);
    auto const old = DistributionMapping::strategy();
    DistributionMapping::strategy(DistributionMapping::ROUNDROBIN);
    DistributionMapping dm(ba);
    DistributionMapping::strategy(old);
    auto mf = std::make_unique<MultiFab>(ba, dm, 1, 0);
    for (MFIter mfi(*mf); mfi.isValid(); ++mfi) {
        auto const& a = mf->array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { a(i,j,k) = fn(i,j,k); });
    }
    return mf;
}

void run_bcdata (ParmParse& pp)
{
    std::string part = "exact"; pp.query("part", part);
    Vector<int> n{8, 7, 6}; pp.queryarr("n_cell", n);
    int mgs = 4; pp.query("mgs", mgs);
    std::string bcs = "ND,DN,NN"; pp.query("bcpairs", bcs);
    std::string be = "fft"; pp.query("backend", be);
    S1 s = make_s1(n, bcs, mgs);
    pb::PressureOptions o; o.verbose = 0; o.removed_mean_warn = 1.0e9; o.tol_rel = 1.0e-13; o.max_iter = 100;
    o.backend = (be == "mlmg") ? pb::BackendKind::MLMG : (be == "hypre") ? pb::BackendKind::HYPRE : pb::BackendKind::FFT;
    const Real h = s.dx;
    bool singular = true;
    for (int f = 0; f < 6; ++f) { singular = singular && s.bc[f] != pb::BC::Dirichlet; }
    const Real pi = Real(3.141592653589793238462643383279502884);
    (void)pi;
    if (part == "exact") {
        // unknown field and data (slabs on odd faces, constants on even faces)
        MultiFab phi_ex(s.ba, s.dm, 1, 0), rhs(s.ba, s.dm, 1, 0), phi(s.ba, s.dm, 1, 1);
        for (MFIter mfi(phi_ex); mfi.isValid(); ++mfi) {
            auto const& a = phi_ex.array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { a(i,j,k) = hash_noise(i,j,k) + Real(0.3)*i - Real(0.2)*j; });
        }
        pb::BoundaryData bd;
        std::vector<std::unique_ptr<MultiFab>> slabs;
        for (int f = 0; f < 6; ++f) {
            if (s.bc[f] == pb::BC::Periodic) { continue; }
            if (f % 2 == 1) { slabs.push_back(make_slab(s, f, [f] (int i, int j, int k) { return bdata_fn(f, i, j, k); })); bd.value[f] = slabs.back().get(); }
            else { bd.constant[f] = Real(0.4) + Real(0.1)*f; }
        }
        // data value at a boundary cell (i,j,k) of face f
        auto data_at = [&] (int f, int i, int j, int k) { return bd.value[f] ? bdata_fn(f, i, j, k) : bd.constant[f]; };
        // explicit full operator with ghost values from the data; evaluated from a gathered copy of phi (rank-local lookup via a global array)
        std::vector<double> G = gather_s(phi_ex, s.domain);
        G.resize(s.domain.numPts());
        ParallelDescriptor::Bcast(G.data(), int(s.domain.numPts()), ParallelDescriptor::IOProcessorNumber());
        const int nx = n[0], ny = n[1];
        auto at = [&] (int i, int j, int k) { return Real(G[i + nx*(j + ny*k)]); };
        for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
            auto const& r = rhs.array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const int c[3] = {i, j, k};
                const Real v0 = at(i, j, k);
                Real sum = 0;
                for (int d = 0; d < 3; ++d) {
                    if (d == 1 && n[1] == 1) { continue; }
                    for (int side = 0; side < 2; ++side) {
                        int t[3] = {i, j, k}; t[d] += side == 0 ? -1 : 1;
                        Real nb;
                        const int f = pb::face_index(d, side);
                        if (t[d] >= 0 && t[d] < n[d]) { nb = at(t[0], t[1], t[2]); }
                        else if (s.bc[f] == pb::BC::Periodic) { t[d] = (t[d] + n[d]) % n[d]; nb = at(t[0], t[1], t[2]); }
                        else if (s.bc[f] == pb::BC::Dirichlet) { nb = Real(2)*data_at(f, i, j, k) - v0; }
                        else { nb = v0 + (side == 0 ? Real(-1) : Real(1)) * h * data_at(f, i, j, k); }
                        sum += nb - v0;
                    }
                    (void)c;
                }
                // each direction contributes (nb_lo + nb_hi - 2 v0)/h^2: both sides were summed above
                r(i,j,k) = sum / (h*h);
            });
        }
        // fold + negative control without fold
        MultiFab rhs0(s.ba, s.dm, 1, 0); MultiFab::Copy(rhs0, rhs, 0, 0, 1, 0);
        pb::PressureProblem p = problem_of(s, rhs, phi);
        std::string msg;
        pb::Status fs = pb::fold_boundary_data(p, bd, rhs, &msg);
        mcheck(fs == pb::Status::Ok, "fold_boundary_data Ok " + msg);
        phi.setVal(0.0);
        pb::PressureResult r = pb::solve_pressure(p, o);
        mcheck(r.status == pb::Status::Ok, "solve status Ok (" + r.message + ")");
        auto ph = gather_s(phi, s.domain), pe = gather_s(phi_ex, s.domain);
        double err = 0, mx = 0;
        if (ParallelDescriptor::IOProcessor()) {
            double shift = 0;
            if (singular) { for (std::size_t q = 0; q < ph.size(); ++q) { shift += ph[q] - pe[q]; } shift /= double(ph.size()); }
            for (std::size_t q = 0; q < ph.size(); ++q) { err = std::max(err, std::abs(ph[q] - shift - pe[q])); mx = std::max(mx, std::abs(pe[q])); }
        }
        ParallelDescriptor::ReduceRealMax(err); ParallelDescriptor::ReduceRealMax(mx);
        Print() << std::setprecision(6) << "BCDATA part=exact bcpairs=" << bcs << " backend=" << be << " singular=" << int(singular) << " max_err_rel=" << err/mx << "\n";
        mcheck(err <= 1.0e-10 * mx, "folded solve reproduces the field whose data-ghost operator built the rhs (max error " + std::to_string(err/mx) + ")");
        // negative control: without the fold the solution is wrong
        MultiFab phi2(s.ba, s.dm, 1, 1); phi2.setVal(0.0);
        pb::PressureProblem p2 = problem_of(s, rhs0, phi2);
        pb::solve_pressure(p2, o);
        auto ph2 = gather_s(phi2, s.domain);
        double err2 = 0;
        if (ParallelDescriptor::IOProcessor()) {
            double shift = 0;
            if (singular) { for (std::size_t q = 0; q < ph2.size(); ++q) { shift += ph2[q] - pe[q]; } shift /= double(ph2.size()); }
            for (std::size_t q = 0; q < ph2.size(); ++q) { err2 = std::max(err2, std::abs(ph2[q] - shift - pe[q])); }
        }
        ParallelDescriptor::ReduceRealMax(err2);
        mcheck(err2 > 1.0e-3 * mx, "negative control: without the fold the error is large (" + std::to_string(err2/mx) + ")");
        // error handling
        {
            pb::BoundaryData bp; pb::PressureProblem pq = p;
            for (int f = 0; f < 6; ++f) { if (s.bc[f] == pb::BC::Periodic) { bp.constant[f] = 1.0; std::string m2; MultiFab rr(s.ba, s.dm, 1, 0); rr.setVal(5.0);
                    mcheck(pb::fold_boundary_data(pq, bp, rr, &m2) == pb::Status::InvalidInput && rr.max(0) == 5.0 && rr.min(0) == 5.0, "data on a periodic face: InvalidInput, rhs untouched"); break; } }
        }
    } else {
        // manufactured solution
        auto u = [] (Real x, Real y, Real z) { return std::cos(1.3*x + 0.2)*std::cosh(0.7*y)*std::sin(0.9*z + 0.4); };
        auto du = [] (int d, Real x, Real y, Real z) {
            if (d == 0) { return -1.3*std::sin(1.3*x + 0.2)*std::cosh(0.7*y)*std::sin(0.9*z + 0.4); }
            if (d == 1) { return 0.7*std::cos(1.3*x + 0.2)*std::sinh(0.7*y)*std::sin(0.9*z + 0.4); }
            return 0.9*std::cos(1.3*x + 0.2)*std::cosh(0.7*y)*std::cos(0.9*z + 0.4);
        };
        MultiFab rhs(s.ba, s.dm, 1, 0), phi(s.ba, s.dm, 1, 1);
        for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
            auto const& r = rhs.array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { r(i,j,k) = Real(-2.01)*u((i+0.5)*h, (j+0.5)*h, (k+0.5)*h); });
        }
        pb::BoundaryData bd;
        std::vector<std::unique_ptr<MultiFab>> slabs;
        for (int f = 0; f < 6; ++f) {
            const int d = f % 3, side = f / 3;
            const bool dir = (s.bc[f] == pb::BC::Dirichlet);
            auto fn = [=] (int i, int j, int k) {
                Real x = (i+0.5)*h, y = (j+0.5)*h, z = (k+0.5)*h;
                const Real w = side == 0 ? Real(0) : Real(n[d])*h;
                if (d == 0) { x = w; } else if (d == 1) { y = w; } else { z = w; }
                return dir ? u(x, y, z) : du(d, x, y, z);
            };
            slabs.push_back(make_slab(s, f, fn)); bd.value[f] = slabs.back().get();
        }
        pb::PressureProblem p = problem_of(s, rhs, phi);
        std::string msg;
        mcheck(pb::fold_boundary_data(p, bd, rhs, &msg) == pb::Status::Ok, "fold_boundary_data Ok");
        phi.setVal(0.0);
        pb::PressureResult r = pb::solve_pressure(p, o);
        mcheck(r.status == pb::Status::Ok, "solve status Ok (" + r.message + ")");
        MultiFab ex(s.ba, s.dm, 1, 0);
        for (MFIter mfi(ex); mfi.isValid(); ++mfi) {
            auto const& a = ex.array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { a(i,j,k) = u((i+0.5)*h, (j+0.5)*h, (k+0.5)*h); });
        }
        auto ph = gather_s(phi, s.domain), pe = gather_s(ex, s.domain);
        double e2 = 0, shift = 0;
        if (ParallelDescriptor::IOProcessor()) {
            if (singular) { for (std::size_t q = 0; q < ph.size(); ++q) { shift += ph[q] - pe[q]; } shift /= double(ph.size()); }
            for (std::size_t q = 0; q < ph.size(); ++q) { const double e = ph[q] - shift - pe[q]; e2 += e*e; }
            e2 = std::sqrt(e2 / double(ph.size()));
        }
        ParallelDescriptor::ReduceRealMax(e2);
        Print() << std::setprecision(8) << "BCDATA part=mms n=" << n[0] << " bcpairs=" << bcs << " backend=" << be << " singular=" << int(singular) << " err_l2=" << e2 << "\n";
    }
}

// --------------------------------------------------------------------------------------------------------------
// mode=fftraw: amrex::FFT::Poisson called directly (no common layer, no effective BC) for the D-057 mapping note.
//   keys: n_cell, bcpairs, out. Same right-hand side as mixed1. Prints FFTRAW and dumps rhs / phi.
// mode=bcmap: selector / status checks of the D-057 mapping (one-cell directions, periodicity guards).
// --------------------------------------------------------------------------------------------------------------
void run_fftraw (ParmParse& pp)
{
    Vector<int> n{16, 1, 16}; pp.queryarr("n_cell", n);
    std::string bcs = "NN,NN,NN"; pp.query("bcpairs", bcs);
    std::string out; pp.query("out", out);
    S1 s = make_s1(n, bcs, 8);
    MultiFab rhs(s.ba, s.dm, 1, 0), phi(s.ba, s.dm, 1, 1);
    for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
        auto const& a = rhs.array(mfi);
        const Real dx = s.dx;
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { a(i,j,k) = smooth_rhs((i+0.5)*dx, (j+0.5)*dx, (k+0.5)*dx) + 0.5*hash_noise(i,j,k) + 0.3; });
    }
    phi.setVal(0.0);
    auto fb = [] (pb::BC b) { return b == pb::BC::Neumann ? FFT::Boundary::even : b == pb::BC::Periodic ? FFT::Boundary::periodic : FFT::Boundary::odd; };
    Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM> fbc{
        std::make_pair(fb(s.bc[0]), fb(s.bc[3])), std::make_pair(fb(s.bc[1]), fb(s.bc[4])), std::make_pair(fb(s.bc[2]), fb(s.bc[5]))};
    FFT::Poisson<MultiFab> solver(s.geom, fbc);
    solver.solve(phi, rhs);
    Print() << "FFTRAW bcpairs=" << bcs << " n=" << n[0] << "x" << n[1] << "x" << n[2] << "\n";
    if (!out.empty()) { write_field(out + "_rhs.bin", gather_s(rhs, s.domain)); write_field(out + "_phi.bin", gather_s(phi, s.domain)); }
}

void run_bcmap (ParmParse&)
{
    struct Case { const char* what; int nx, ny, nz; const char* pairs; pb::Status st; const char* frag; int singular; };
    const Case cases[] = {
        // the strings of the driver sweep, with the one-cell y (n = 16 x 1 x 16)
        {"PP,NN,PP (91 inputs)", 16, 1, 16, "PP,NN,PP", pb::Status::Ok, "", 1},
        {"NN,NN,NN (86 inputs)", 16, 1, 16, "NN,NN,NN", pb::Status::Ok, "", 1},
        {"PP,NN,NN (25 inputs)", 16, 1, 16, "PP,NN,NN", pb::Status::Ok, "", 1},
        {"DD,DD,DD (1 input: fully open 2-D box, y follows)", 16, 1, 16, "DD,DD,DD", pb::Status::Ok, "", 0},
        {"ND,NN,NN (17 inputs)", 16, 1, 16, "ND,NN,NN", pb::Status::Ok, "", 0},
        {"DD,NN,NN (4 inputs)", 16, 1, 16, "DD,NN,NN", pb::Status::Ok, "", 0},
        {"NN,NN,ND (4 inputs)", 16, 1, 16, "NN,NN,ND", pb::Status::Ok, "", 0},
        {"DD,NN,ND (3 inputs)", 16, 1, 16, "DD,NN,ND", pb::Status::Ok, "", 0},
        {"ND,NN,ND (2 inputs)", 16, 1, 16, "ND,NN,ND", pb::Status::Ok, "", 0},
        {"DN,NN,DD (1 input)", 16, 1, 16, "DN,NN,DD", pb::Status::Ok, "", 0},
        // Dirichlet faces in the one-cell y (never produced by the driver for the verification set): act as Neumann
        {"one-cell y: D faces only in y -> singular (y term does not exist)", 16, 1, 16, "NN,DD,NN", pb::Status::Ok, "", 1},
        {"one-cell y: ND in y with closed x, z -> singular", 16, 1, 16, "NN,ND,PP", pb::Status::Ok, "", 1},
        {"one-cell y: DN in y, x open -> nonsingular", 16, 1, 16, "DD,DN,NN", pb::Status::Ok, "", 0},
        // one-cell x or z
        {"one-cell x, Neumann (exact)", 1, 16, 16, "NN,NN,NN", pb::Status::Ok, "", 1},
        {"one-cell x, periodic (exact)", 1, 16, 16, "PP,NN,NN", pb::Status::Ok, "", 1},
        {"one-cell z, periodic (exact)", 16, 16, 1, "NN,NN,PP", pb::Status::Ok, "", 1},
        {"one-cell x, DD refused", 1, 16, 16, "DD,NN,NN", pb::Status::NotBuilt, "one-cell x", -1},
        {"one-cell x, ND refused", 1, 16, 16, "ND,NN,NN", pb::Status::NotBuilt, "one-cell x", -1},
        {"one-cell x, DN refused", 1, 16, 16, "DN,DD,DD", pb::Status::NotBuilt, "one-cell x", -1},
        {"one-cell z, DD refused", 16, 16, 1, "NN,NN,DD", pb::Status::NotBuilt, "one-cell z", -1},
        {"one-cell z, DN refused", 16, 16, 1, "DD,DD,DN", pb::Status::NotBuilt, "one-cell z", -1},
    };
    for (auto const& c : cases) {
        for (pb::BackendKind req : {pb::BackendKind::Auto, pb::BackendKind::MLMG}) {
            S1 s = make_s1(Vector<int>{c.nx, c.ny, c.nz}, c.pairs, 8);
            MultiFab rhs(s.ba, s.dm, 1, 0), phi(s.ba, s.dm, 1, 1);
            fill_smooth(s, rhs, 0.3); phi.setVal(7.0);
            pb::PressureProblem p = problem_of(s, rhs, phi);
            pb::PressureOptions o; o.verbose = 0; o.backend = req; o.removed_mean_warn = 1.0e9;
            pb::PressureResult r = pb::solve_pressure(p, o);
            const std::string tag = std::string(c.what) + (req == pb::BackendKind::Auto ? " [auto]" : " [mlmg]");
            bool good = (r.status == c.st);
            if (c.st == pb::Status::Ok) {
                good = good && r.backend == (req == pb::BackendKind::Auto ? "FFT" : "MLMG") && !r.components.empty() && int(r.components[0].singular) == c.singular;
                good = good && r.residual_checked && r.residual_ok && r.warnings.empty();
            } else {
                good = good && r.message.find(c.frag) != std::string::npos && phi.min(0) == 7.0 && phi.max(0) == 7.0;
            }
            Print() << "  " << tag << ": status=" << pb::to_string(r.status) << " backend=" << r.backend << " singular=" << (r.components.empty() ? -1 : int(r.components[0].singular))
                    << " true_rel2=" << r.residual_rel2 << " nwarn=" << r.warnings.size() << (r.message.empty() ? "" : " msg=\"" + r.message + "\"") << "\n";
            mcheck(good, tag);
        }
    }
    // periodicity guards: the interface refuses a BC that does not match the Geometry (a driver "code 0 on a non-periodic direction" never gets here)
    {
        S1 s = make_s1(Vector<int>{16, 1, 16}, "NN,NN,NN", 8);
        MultiFab rhs(s.ba, s.dm, 1, 0), phi(s.ba, s.dm, 1, 1);
        fill_smooth(s, rhs, 0.0); phi.setVal(7.0);
        pb::PressureProblem p = problem_of(s, rhs, phi);
        pb::PressureOptions o; o.verbose = 0;
        p.bc[0] = pb::BC::Periodic; p.bc[3] = pb::BC::Periodic;      // periodic BC, non-periodic Geometry
        pb::PressureResult r = pb::solve_pressure(p, o);
        mcheck(r.status == pb::Status::InvalidInput && r.message.find("eriodic") != std::string::npos && phi.max(0) == 7.0, "periodic BC on a non-periodic Geometry: InvalidInput (" + r.message + ")");
        p.bc[0] = pb::BC::Periodic; p.bc[3] = pb::BC::Neumann;
        r = pb::solve_pressure(p, o);
        mcheck(r.status == pb::Status::InvalidInput, "periodic on one face only: InvalidInput (" + r.message + ")");
        S1 sp = make_s1(Vector<int>{16, 1, 16}, "PP,NN,NN", 8);
        MultiFab rhs2(sp.ba, sp.dm, 1, 0), phi2(sp.ba, sp.dm, 1, 1); fill_smooth(sp, rhs2, 0.0); phi2.setVal(7.0);
        pb::PressureProblem q = problem_of(sp, rhs2, phi2);
        q.bc[0] = pb::BC::Neumann; q.bc[3] = pb::BC::Neumann;           // periodic Geometry, Neumann BC
        r = pb::solve_pressure(q, o);
        mcheck(r.status == pb::Status::InvalidInput, "periodic Geometry with a Neumann BC: InvalidInput (" + r.message + ")");
    }
    // effective_bc() itself
    {
        Box d1(IntVect(0), IntVect(15, 0, 15)), d3(IntVect(0), IntVect(15));
        std::array<pb::BC,6> a{pb::BC::Dirichlet, pb::BC::Dirichlet, pb::BC::Dirichlet, pb::BC::Neumann, pb::BC::Dirichlet, pb::BC::Periodic};
        auto e1 = pb::effective_bc(a, d1), e3 = pb::effective_bc(a, d3);
        mcheck(e1[1] == pb::BC::Neumann && e1[4] == pb::BC::Neumann && e1[0] == pb::BC::Dirichlet && e1[2] == pb::BC::Dirichlet && e1[5] == pb::BC::Periodic, "effective_bc: one-cell y turns only the y Dirichlet faces into Neumann");
        mcheck(e3 == a, "effective_bc: unchanged when y has more than one cell");
    }
}

} // namespace


// --------------------------------------------------------------------------------------------------------------
// mode=masked: masked single-level problem (frozen/masked-notes.md). Random Solid cells (solid_frac), an optional full Solid plane at
// k = slab that splits the box, random Known cells above the slab (known_frac, value 0.3 + smooth), a noisy right-hand side that is also
// non-zero on excluded cells, optional rho/KRES gauge fields. Writes <out>_{class,known,rhs,phi,label,rho,kres}.bin (doubles, x fastest)
// and prints MASKED lines; tests/pb_test.py masked rebuilds the problem densely with numpy and compares.
// --------------------------------------------------------------------------------------------------------------
void run_masked (ParmParse& pp)
{
    Vector<int> n{12, 10, 8}; pp.queryarr("n_cell", n);
    int mgs = 8; pp.query("mgs", mgs);
    std::string bcs = "NN,NN,ND"; pp.query("bcpairs", bcs);
    double solid_frac = 0.10; pp.query("solid_frac", solid_frac);
    double known_frac = 0.0; pp.query("known_frac", known_frac);
    int slab = -1; pp.query("slab", slab);
    int seed = 1; pp.query("seed", seed);
    int gauge_rho = 0; pp.query("gauge_rho", gauge_rho);
    std::string be = "mlmg"; pp.query("backend", be);
    std::string out; pp.query("out", out);
    double tol = 1.0e-12; pp.query("tol_rel", tol);
    int max_iter = 200; pp.query("max_iter", max_iter);
    S1 s = make_s1(n, bcs, mgs);
    iMultiFab cls(s.ba, s.dm, 1, 0);
    MultiFab kv(s.ba, s.dm, 1, 0), rhs(s.ba, s.dm, 1, 0), phi(s.ba, s.dm, 1, 1), rho(s.ba, s.dm, 1, 0), kres(s.ba, s.dm, 1, 0);
    for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
        auto const& c = cls.array(mfi);
        auto const& g = kv.array(mfi);
        auto const& r = rhs.array(mfi);
        auto const& ro = rho.array(mfi);
        auto const& kr = kres.array(mfi);
        const Real dx = s.dx;
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const double h1 = 0.5*(hash_noise(i + 7*seed, j, k) + 1.0);          // in [0,1]
            const double h2 = 0.5*(hash_noise(i, j + 11*seed, k + 3) + 1.0);
            int cc = int(pb::CellGas);
            if (k == slab) { cc = pb::CellSolid; }
            else if (h1 < solid_frac) { cc = pb::CellSolid; }
            else if (known_frac > 0.0 && k > slab && h2 < known_frac) { cc = pb::CellKnown; }
            c(i,j,k) = cc;
            const Real x = (i+0.5)*dx, y = (j+0.5)*dx, z = (k+0.5)*dx;
            g(i,j,k) = Real(0.3) + smooth_rhs(x, y, z);
            r(i,j,k) = smooth_rhs(x, y, z) + Real(0.5)*hash_noise(i, j, k) + Real(0.2);
            ro(i,j,k) = Real(1.0) + Real(0.5)*Real(std::sin(3.0*x + y)) * Real(std::cos(2.0*z));
            kr(i,j,k) = Real(0.5)*(x*x + y*y + z*z);
        });
    }
    phi.setVal(0.0);
    pb::PressureProblem p = problem_of(s, rhs, phi);
    p.cell_class = &cls;
    p.known_value = &kv;
    if (gauge_rho) { p.gauge_weight = &rho; p.gauge_offset = &kres; }
    pb::PressureOptions o; o.verbose = 0; o.removed_mean_warn = 1.0e9;   // the rhs is deliberately incompatible
    o.backend = (be == "fft") ? pb::BackendKind::FFT : (be == "hypre") ? pb::BackendKind::HYPRE : (be == "mlmg") ? pb::BackendKind::MLMG : pb::BackendKind::Auto;
    o.tol_rel = tol; o.max_iter = max_iter;
    pb::ComponentMap cm = pb::label_components(p);
    pb::PressureResult r = pb::solve_pressure(p, o);
    Print() << std::setprecision(8) << "MASKED status=" << pb::to_string(r.status) << " backend=" << (r.backend.empty() ? std::string("none") : r.backend)
            << " iters=" << r.backend_status.iterations << " ncomp=" << r.components.size() << " residual_ok=" << int(r.residual_ok) << " nwarn=" << r.warnings.size()
            << " true_rel2=" << r.residual_rel2 << " limit=" << r.residual_limit << " msg=[" << r.message << "]\n";
    for (auto const& c : r.components) {
        Print() << std::setprecision(17) << "MASKEDCOMP id=" << c.id << " singular=" << int(c.singular) << " ncells=" << c.ncells << " pin=" << c.pin[0] << "," << c.pin[1] << "," << c.pin[2]
                << " removed_mean=" << c.removed_mean << " gauge_shift=" << c.gauge_shift << "\n";
    }
    if (r.status == pb::Status::Ok || r.status == pb::Status::NotConverged) {
        // exact checks done here: excluded cells hold 0 / g
        double e_solid = 0.0, e_known = 0.0;
        for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
            auto const& f = phi.const_array(mfi);
            auto const& c = cls.const_array(mfi);
            auto const& g = kv.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                if (c(i,j,k) == pb::CellSolid) { e_solid = std::max(e_solid, std::abs(double(f(i,j,k)))); }
                if (c(i,j,k) == pb::CellKnown) { e_known = std::max(e_known, std::abs(double(f(i,j,k)) - double(g(i,j,k)))); }
            });
        }
        ParallelDescriptor::ReduceRealMax(e_solid); ParallelDescriptor::ReduceRealMax(e_known);
        mcheck(e_solid == 0.0 && e_known == 0.0, "masked: H is 0 on Solid cells and g on Known cells");
    }
    if (!out.empty()) {
        MultiFab clsd(s.ba, s.dm, 1, 0), lab(s.ba, s.dm, 1, 0);
        for (MFIter mfi(clsd); mfi.isValid(); ++mfi) {
            auto const& d = clsd.array(mfi); auto const& c = cls.const_array(mfi); auto const& l = lab.array(mfi); auto const& cl = cm.label.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { d(i,j,k) = c(i,j,k); l(i,j,k) = cl(i,j,k); });
        }
        write_field(out + "_class.bin", gather_s(clsd, s.domain));
        write_field(out + "_known.bin", gather_s(kv, s.domain));
        write_field(out + "_rhs.bin", gather_s(rhs, s.domain));
        write_field(out + "_phi.bin", gather_s(phi, s.domain));
        write_field(out + "_label.bin", gather_s(lab, s.domain));
        write_field(out + "_rho.bin", gather_s(rho, s.domain));
        write_field(out + "_kres.bin", gather_s(kres, s.domain));
    }
}

int run_m2_mode (std::string const& mode, ParmParse& pp)
{
    if (mode == "trigger1") { run_trigger1(pp); }
    else if (mode == "fftcache") { run_fftcache(pp); }
    else if (mode == "hypcache") { run_hypcache(pp); }
    else if (mode == "mixed1") { run_mixed1(pp); }
    else if (mode == "bcdata") { run_bcdata(pp); }
    else if (mode == "fftraw") { run_fftraw(pp); }
    else if (mode == "bcmap") { run_bcmap(pp); }
    else if (mode == "masked") { run_masked(pp); }
    else { return -1; }
    return g_m2fail;
}
