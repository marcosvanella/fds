// Standalone harness for the pressure backend interface (no FDS sources). Arguments are key=value.
//
// mode=solve    (default) build a single-level problem, solve it with each backend in `backends`, print one
//               RESULT line per backend (iterations, residuals, removed mean, pin, hash of H) and CMP lines
//               (rel. L2 against the first backend, or against `ref_file`) with the eps_H verdict.
// mode=gen      write the synthetic RHS and the FFT reference H (gauge-fixed) as raw float64, x fastest.
// mode=selector check that unsupported requests return an explicit "not built" status.
// mode=exactsum exact-sum and mean-removal decomposition checks on a deterministic wide-range field.
// mode=meankind per-cell volume mean removal (both MeanKind) and the rho*volume, KRES gauge on a synthetic stretched
//               volume field; writes raw files for the independent numpy check (see tests/pb_test.py, meankind).
// mode=comp, comp_ns2d, comp_sel, comp_ws, comp_gauge  composite (multi-level) pressure solve tests, see composite_modes.cpp.
// mode=trigger1, fftcache, ...  M2 single-level modes (m2_modes.cpp).
// mode=diff     compare two raw fields: a=<file> b=<file> n_cell="nx ny nz" (rel. L2, max abs, eps_H verdict).
#include "PressureIface.H"
#include "CommonLayer.H"
#include "ExactSum.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>

using namespace amrex;

namespace {

int g_fail = 0;
void check (bool ok, std::string const& what)
{
    Print() << "CHECK " << (ok ? "PASS " : "FAIL ") << what << "\n";
    if (!ok) { ++g_fail; }
}

pb::BC parse_bc (std::string const& s)
{
    if (s == "neumann") { return pb::BC::Neumann; }
    if (s == "periodic") { return pb::BC::Periodic; }
    if (s == "dirichlet") { return pb::BC::Dirichlet; }
    amrex::Abort("bc must be neumann|periodic|dirichlet");
    return pb::BC::Neumann;
}

// Same smooth non-separable RHS as prototype P2 (cell centres on the unit cube).
Real rhs_func (Real x, Real y, Real z)
{
    const Real pi = Real(3.141592653589793238462643383279502884);
    Real r2 = (x-Real(0.31))*(x-Real(0.31)) + (y-Real(0.57))*(y-Real(0.57)) + (z-Real(0.45))*(z-Real(0.45));
    return std::exp(-r2/Real(0.02)) + Real(0.5)*std::sin(Real(3.)*pi*x)*std::cos(Real(2.)*pi*y+Real(0.4))*(z+Real(0.25));
}

// Variable density and a KRES-like offset for the gauge tests (cell centres on the unit cube).
Real rho_func (Real x, Real y, Real z)
{
    const Real pi = Real(3.141592653589793238462643383279502884);
    return Real(1.2)*(Real(1.0) + Real(0.3)*std::sin(Real(2.)*pi*x)*std::cos(pi*y) + Real(0.2)*z);
}
Real kres_func (Real x, Real y, Real z) { return Real(0.05)*(x*x + Real(2.)*y - z); }

struct Setup {
    Box domain;
    Geometry geom;
    BoxArray ba;
    DistributionMapping dm;
    std::array<pb::BC,6> bc;
    int nmax = 0;
    Real dx = 0;
    double eps_H = 0;
};

Setup make_setup (ParmParse& pp)
{
    Setup s;
    Vector<int> n_cell{64, 64, 64};
    pp.queryarr("n_cell", n_cell);
    int mgs = 32; pp.query("max_grid_size", mgs);
    std::string bcs = "neumann"; pp.query("bc", bcs);
    s.bc.fill(parse_bc(bcs));
    s.nmax = std::max({n_cell[0], n_cell[1], n_cell[2]});
    s.dx = Real(1.0) / s.nmax;
    s.domain = Box(IntVect(0), IntVect(AMREX_D_DECL(n_cell[0]-1, n_cell[1]-1, n_cell[2]-1)));
    RealBox rb({0.,0.,0.}, {n_cell[0]*s.dx, n_cell[1]*s.dx, n_cell[2]*s.dx});
    const int per = (s.bc[0] == pb::BC::Periodic) ? 1 : 0;
    s.geom = Geometry(s.domain, rb, CoordSys::cartesian, Array<int,3>{per, per, per});
    s.ba = BoxArray(s.domain);
    s.ba.maxSize(mgs);
    s.dm = DistributionMapping(s.ba);
    s.eps_H = std::max(1.0e-8, 2.4e-12*double(s.nmax)*double(s.nmax));
    return s;
}

void fill_rhs (Setup const& s, MultiFab& rhs)
{
    for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
        auto const& a = rhs.array(mfi);
        const Real dx = s.dx;
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            a(i,j,k) = rhs_func((i+0.5)*dx, (j+0.5)*dx, (k+0.5)*dx);
        });
    }
}

// The driver supplies a compatible RHS (FR-030: S formed from DSUM/USUM). For the synthetic field the
// harness emulates that by removing the mean once with the common layer; the common layer's own removal
// then finds only round-off (removed-mean diagnostic quiet).
void make_compatible (Setup const& s, MultiFab& rhs)
{
    pb::PressureProblem p;
    p.ba = s.ba; p.dm = s.dm; p.geom = s.geom; p.bc = s.bc;
    pb::ComponentMap cm = pb::label_components(p);
    if (cm.comps[0].singular) { pb::remove_mean(rhs, cm, nullptr, s.dx*s.dx*s.dx); }
}

// Gather a MultiFab to rank 0 as one Fab (x fastest).
std::vector<double> gather (MultiFab const& mf, Box const& domain)
{
    BoxArray ba1(domain);
    DistributionMapping dm1(Vector<int>{ParallelDescriptor::IOProcessorNumber()});
    MultiFab g(ba1, dm1, 1, 0);
    g.ParallelCopy(mf, 0, 0, 1);
    std::vector<double> v;
    if (ParallelDescriptor::IOProcessor()) {
        v.assign(g[0].dataPtr(), g[0].dataPtr() + domain.numPts());
    }
    return v;
}

std::uint64_t fnv1a (std::vector<double> const& v)
{
    std::uint64_t h = 1469598103934665603ULL;
    auto const* p = reinterpret_cast<const unsigned char*>(v.data());
    for (std::size_t i = 0; i < v.size()*sizeof(double); ++i) { h ^= p[i]; h *= 1099511628211ULL; }
    return h;
}

std::string hex64 (std::uint64_t h) { std::ostringstream o; o << std::hex << std::setw(16) << std::setfill('0') << h; return o.str(); }
std::string hexd (double x) { std::uint64_t u; std::memcpy(&u, &x, 8); return hex64(u); }

void write_raw (std::vector<double> const& v, std::string const& fname)
{
    if (!ParallelDescriptor::IOProcessor()) { return; }
    std::ofstream ofs(fname, std::ios::binary);
    ofs.write(reinterpret_cast<const char*>(v.data()), static_cast<std::streamsize>(v.size()*sizeof(double)));
}

void read_raw (std::string const& fname, Box const& domain, MultiFab& mf)
{
    BoxArray ba1(domain);
    DistributionMapping dm1(Vector<int>{ParallelDescriptor::IOProcessorNumber()});
    MultiFab g(ba1, dm1, 1, 0);
    if (ParallelDescriptor::IOProcessor()) {
        std::ifstream ifs(fname, std::ios::binary);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ifs.good(), "cannot open " + fname);
        ifs.read(reinterpret_cast<char*>(g[0].dataPtr()), static_cast<std::streamsize>(domain.numPts()*sizeof(double)));
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ifs.good(), "short read " + fname);
    }
    mf.ParallelCopy(g, 0, 0, 1);
}

struct Diff { double max_abs, rel_max, rel_l2; };
Diff compare (MultiFab const& a, MultiFab const& b)
{
    MultiFab d(a.boxArray(), a.DistributionMap(), 1, 0);
    MultiFab::Copy(d, a, 0, 0, 1, 0);
    MultiFab::Subtract(d, b, 0, 0, 1, 0);
    Diff r;
    r.max_abs = d.norm0(0);
    r.rel_max = r.max_abs / b.norm0(0);
    r.rel_l2 = d.norm2(0) / b.norm2(0);
    return r;
}

void print_diff (std::string const& label, Diff const& d, double eps_H)
{
    Print() << std::setprecision(6) << "CMP " << label << " rel_l2=" << d.rel_l2 << " max_abs=" << d.max_abs
            << " rel_max=" << d.rel_max << " eps_H=" << eps_H << " " << ((d.rel_l2 <= eps_H) ? "PASS" : "FAIL") << "\n";
    if (!(d.rel_l2 <= eps_H)) { ++g_fail; }
}

pb::BackendKind parse_backend (std::string const& s)
{
    if (s == "fft") { return pb::BackendKind::FFT; }
    if (s == "mlmg") { return pb::BackendKind::MLMG; }
    if (s == "hypre") { return pb::BackendKind::HYPRE; }
    if (s == "auto") { return pb::BackendKind::Auto; }
    amrex::Abort("backend must be fft|mlmg|auto");
    return pb::BackendKind::Auto;
}

// ---------------------------------------------------------------------------------------------
void run_solve (ParmParse& pp)
{
    Setup s = make_setup(pp);
    std::vector<std::string> backends;
    if (!pp.queryarr("backends", backends) || backends.empty()) { backends = {"fft", "mlmg"}; }
    std::string out; pp.query("out", out);
    std::string rhs_file; pp.query("rhs_file", rhs_file);
    std::string ref_file; pp.query("ref_file", ref_file);
    pb::PressureOptions o;
    pp.query("tol_rel", o.tol_rel);
    {   // HYPRE settings: krylov=auto|pcg|gmres|bicgstab coarsen= relax= sweeps= strong= interp= agg=
        std::string k;
        if (pp.query("krylov", k)) { o.hypre.krylov = k == "pcg" ? pb::HypreKrylov::PCG : k == "gmres" ? pb::HypreKrylov::GMRES : k == "bicgstab" ? pb::HypreKrylov::BiCGSTAB : pb::HypreKrylov::Auto; }
        pp.query("coarsen", o.hypre.coarsen_type); pp.query("relax", o.hypre.relax_type); pp.query("sweeps", o.hypre.num_sweeps);
        pp.query("strong", o.hypre.strong_threshold); pp.query("interp", o.hypre.interp_type); pp.query("agg", o.hypre.agg_levels);
        pp.query("max_iter", o.max_iter);
    }
    int rm = 1; pp.query("remove_mean", rm); o.remove_mean = (rm != 0);
    pp.query("verbose", o.verbose);
    double inject = 0.0; pp.query("rhs_offset", inject);   // adds a constant to the RHS before the solve
    int gauge_rho = 0; pp.query("gauge_rho", gauge_rho);   // 1: rho-weighted gauge with KRES offset (FDS SYMM_INDEFINITE)
    std::string mean_kind = "volume"; pp.query("mean_kind", mean_kind);   // volume (D-067 default) | scaled (FDS parity switch)
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mean_kind == "volume" || mean_kind == "scaled", "mean_kind must be volume|scaled");

    Print() << "PB solve: n_cell=" << s.domain.length(0) << "x" << s.domain.length(1) << "x" << s.domain.length(2)
            << " bc=" << pb::to_string(s.bc[0]) << " nboxes=" << s.ba.size() << " nranks=" << ParallelDescriptor::NProcs()
            << " eps_H=" << s.eps_H << "\n";

    MultiFab rhs(s.ba, s.dm, 1, 0);
    if (rhs_file.empty()) { fill_rhs(s, rhs); make_compatible(s, rhs); } else { read_raw(rhs_file, s.domain, rhs); }
    if (inject != 0.0) { rhs.plus(inject, 0, 1, 0); }

    MultiFab rho(s.ba, s.dm, 1, 0), kres(s.ba, s.dm, 1, 0);
    if (gauge_rho) {
        for (MFIter mfi(rho); mfi.isValid(); ++mfi) {
            auto const& r = rho.array(mfi); auto const& k_ = kres.array(mfi);
            const Real dx = s.dx;
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                r(i,j,k) = rho_func((i+0.5)*dx, (j+0.5)*dx, (k+0.5)*dx);
                k_(i,j,k) = kres_func((i+0.5)*dx, (j+0.5)*dx, (k+0.5)*dx);
            });
        }
        if (!out.empty()) { write_raw(gather(rho, s.domain), out + "_rho.bin"); write_raw(gather(kres, s.domain), out + "_kres.bin"); }
    }

    MultiFab ref;
    if (!ref_file.empty()) { ref.define(s.ba, s.dm, 1, 0); read_raw(ref_file, s.domain, ref); }

    std::vector<std::pair<std::string, std::unique_ptr<MultiFab>>> sol;
    for (auto const& name : backends) {
        pb::PressureProblem p;
        p.ba = s.ba; p.dm = s.dm; p.geom = s.geom; p.bc = s.bc;
        p.rhs = &rhs;
        p.mean_kind = (mean_kind == "scaled") ? pb::MeanKind::ScaledArithmetic : pb::MeanKind::Volume;
        if (gauge_rho) { p.gauge_weight = &rho; p.gauge_offset = &kres; }
        auto phi = std::make_unique<MultiFab>(s.ba, s.dm, 1, 1);
        phi->setVal(0.0);
        p.phi = phi.get();
        o.backend = parse_backend(name);
        const double t_solve0 = amrex::second();
        pb::PressureResult r = pb::solve_pressure(p, o);
        const double t_solve = amrex::second() - t_solve0;
        Print() << std::setprecision(6) << "RESULT " << name << " status=" << pb::to_string(r.status)
                << " backend=" << r.backend << " iters=" << r.backend_status.iterations
                << " own_res=" << r.backend_status.own_residual
                << " true_rel2=" << r.residual_rel2 << " true_relmax=" << r.residual_relmax
                << " residual_ok=" << (r.residual_ok ? 1 : 0) << " nwarn=" << r.warnings.size() << " seconds=" << t_solve;
        if (r.backend == "HYPRE") { Print() << " hypre_setup=" << r.backend_status.hypre_setup_seconds << " hypre_solve=" << r.backend_status.hypre_solve_seconds << " hypre_method=" << r.backend_status.hypre_method; }
        for (auto const& c : r.components) {
            Print() << " comp" << c.id << "_singular=" << (c.singular ? 1 : 0)
                    << " comp" << c.id << "_pin=" << c.pin[0] << "," << c.pin[1] << "," << c.pin[2]
                    << " comp" << c.id << "_removed_mean=" << hexd(c.removed_mean)
                    << " comp" << c.id << "_removed_rel=" << c.removed_rel
                    << " comp" << c.id << "_gauge_shift=" << hexd(c.gauge_shift);
        }
        std::vector<double> v = gather(*phi, s.domain);
        Print() << " hash=" << hex64(fnv1a(v)) << "\n";
        if (r.status != pb::Status::Ok) { ++g_fail; continue; }
        if (!out.empty()) { write_raw(v, out + "_" + name + ".bin"); }
        MultiFab h(s.ba, s.dm, 1, 0);
        MultiFab::Copy(h, *phi, 0, 0, 1, 0);
        if (!ref_file.empty()) { print_diff(name + "_vs_ref", compare(h, ref), s.eps_H); }
        sol.emplace_back(name, std::make_unique<MultiFab>(std::move(h)));
    }
    if (ref_file.empty()) {
        for (std::size_t i = 1; i < sol.size(); ++i) {
            print_diff(sol[i].first + "_vs_" + sol[0].first, compare(*sol[i].second, *sol[0].second), s.eps_H);
        }
    }
}

void run_gen (ParmParse& pp)
{
    Setup s = make_setup(pp);
    std::string out; pp.get("out", out);
    MultiFab rhs(s.ba, s.dm, 1, 0);
    fill_rhs(s, rhs);
    make_compatible(s, rhs);
    pb::PressureProblem p;
    p.ba = s.ba; p.dm = s.dm; p.geom = s.geom; p.bc = s.bc; p.rhs = &rhs;
    MultiFab phi(s.ba, s.dm, 1, 1); phi.setVal(0.0); p.phi = &phi;
    pb::PressureOptions o; o.backend = pb::BackendKind::FFT;
    pb::PressureResult r = pb::solve_pressure(p, o);
    AMREX_ALWAYS_ASSERT(r.status == pb::Status::Ok);
    // Store the RHS as the driver would supply it (unmodified) and H (gauge-fixed by the common layer).
    write_raw(gather(rhs, s.domain), out + "_rhs.bin");
    write_raw(gather(phi, s.domain), out + "_H.bin");
    Print() << "GEN " << out << " rhs_fnv=" << hex64(fnv1a(gather(rhs, s.domain)))
            << " H_fnv=" << hex64(fnv1a(gather(phi, s.domain))) << "\n";
}

void run_diff (ParmParse& pp)
{
    Setup s = make_setup(pp);
    std::string a, b; pp.get("a", a); pp.get("b", b);
    MultiFab fa(s.ba, s.dm, 1, 0), fb(s.ba, s.dm, 1, 0);
    read_raw(a, s.domain, fa); read_raw(b, s.domain, fb);
    print_diff(a + "_vs_" + b, compare(fa, fb), s.eps_H);
}

// ---------------------------------------------------------------------------------------------
void run_selector ()
{
    Box dom(IntVect(0), IntVect(7));
    RealBox rb({0.,0.,0.}, {1.,1.,1.});
    auto make = [&] (bool per, pb::BC bc) {
        pb::PressureProblem p;
        p.geom = Geometry(dom, rb, CoordSys::cartesian, Array<int,3>{per?1:0, per?1:0, per?1:0});
        p.ba = BoxArray(dom); p.ba.maxSize(4);
        p.dm = DistributionMapping(p.ba);
        p.bc.fill(bc);
        return p;
    };
    MultiFab rhs(BoxArray(dom).maxSize(4), DistributionMapping(BoxArray(dom).maxSize(4)), 1, 0);
    rhs.setVal(0.0);
    auto expect = [&] (std::string const& what, pb::PressureProblem p, pb::Status st, pb::BackendKind want,
                       pb::BackendKind req) {
        p.rhs = &rhs;
        MultiFab phi(p.ba, p.dm, 1, 1); phi.setVal(7.0); p.phi = &phi;
        pb::PressureOptions o; o.verbose = 0; o.backend = req;
        pb::Selection sel = pb::select_backend(p, req);
        pb::PressureResult r = pb::solve_pressure(p, o);
        bool ok = (r.status == st) && (st != pb::Status::NotBuilt || (!sel.ok && !r.message.empty()))
                  && (st == pb::Status::NotBuilt || sel.kind == want);
        if (st == pb::Status::NotBuilt) { ok = ok && phi.min(0) == 7.0 && phi.max(0) == 7.0; }   // untouched
        Print() << "  " << what << ": status=" << pb::to_string(r.status)
                << (r.message.empty() ? "" : " msg=\"" + r.message + "\"") << "\n";
        check(ok, "selector " + what);
    };
    using K = pb::BackendKind;
    expect("single level neumann -> FFT", make(false, pb::BC::Neumann), pb::Status::Ok, K::FFT, K::Auto);
    expect("single level periodic -> FFT", make(true, pb::BC::Periodic), pb::Status::Ok, K::FFT, K::Auto);
    expect("single level dirichlet -> FFT", make(false, pb::BC::Dirichlet), pb::Status::Ok, K::FFT, K::Auto);
    expect("explicit MLMG request", make(false, pb::BC::Neumann), pb::Status::Ok, K::MLMG, K::MLMG);
    {   pb::PressureProblem p = make(false, pb::BC::Neumann); p.nlevels = 2;
        expect("composite (2 levels) not built", p, pb::Status::NotBuilt, K::Auto, K::Auto);
        expect("composite with explicit MLMG not built", p, pb::Status::NotBuilt, K::Auto, K::MLMG); }
    {   pb::PressureProblem p = make(false, pb::BC::Neumann); p.bc[pb::face_index(0,1)] = pb::BC::Dirichlet;
        expect("mixed open/closed faces -> FFT (M2)", p, pb::Status::Ok, K::FFT, K::Auto);
        expect("mixed open/closed faces with explicit MLMG", p, pb::Status::Ok, K::MLMG, K::MLMG); }
    {   pb::PressureProblem p = make(false, pb::BC::Neumann);
        iMultiFab cls(p.ba, p.dm, 1, 0); cls.setVal(0);
        for (MFIter mfi(cls); mfi.isValid(); ++mfi) { if (mfi.index() == 0) { const Box vb = mfi.validbox(); const IntVect sm = vb.smallEnd(); cls[mfi].setVal<RunOn::Host>(1, Box(sm, sm)); } }
        p.cell_class = &cls;
        expect("masked cell not built", p, pb::Status::NotBuilt, K::Auto, K::Auto);
        iMultiFab zero(p.ba, p.dm, 1, 0); zero.setVal(0);
        p.cell_class = &zero;
        expect("all-zero cell_class is unmasked", p, pb::Status::Ok, K::FFT, K::Auto); }
    {   pb::PressureProblem p = make(false, pb::BC::Neumann);
        iMultiFab unc(p.ba, p.dm, 1, 0); unc.setVal(1);
        for (MFIter mfi(unc); mfi.isValid(); ++mfi) { if (mfi.index() == 0) { const Box vb = mfi.validbox(); const IntVect sm = vb.smallEnd(); unc[mfi].setVal<RunOn::Host>(0, Box(sm, sm)); } }
        p.uncovered = &unc;
        expect("covered cell not built", p, pb::Status::NotBuilt, K::Auto, K::Auto); }
    {   pb::PressureProblem p = make(false, pb::BC::Neumann); p.cylindrical = true;
        expect("cylindrical geometry not built", p, pb::Status::NotBuilt, K::Auto, K::Auto); }
    {   pb::PressureProblem p = make(false, pb::BC::Neumann);
        for (int d = 0; d < 3; ++d) { p.cell_width[d].assign(8, Real(0.125)); }
        expect("explicit uniform cell widths -> FFT", p, pb::Status::Ok, K::FFT, K::Auto);
        p.cell_width[2][3] = Real(0.1875);
        expect("non-uniform cell widths (stretched z) not built", p, pb::Status::NotBuilt, K::Auto, K::Auto);
        expect("non-uniform cell widths with explicit MLMG not built", p, pb::Status::NotBuilt, K::Auto, K::MLMG);
        for (int d = 0; d < 3; ++d) { p.cell_width[d].assign(8, Real(0.2)); }
        expect("uniform cell widths that disagree with geometry are invalid", p, pb::Status::InvalidInput, K::FFT, K::Auto); }
    {   pb::PressureProblem p = make(false, pb::BC::Neumann); MultiFab a(p.ba, p.dm, 1, 0); a.setVal(1.0);
        p.cell_coef_a = &a;
        expect("cell coefficient a not built", p, pb::Status::NotBuilt, K::Auto, K::Auto); }
    // HYPRE (assembled-matrix backend): explicit request only, never selected by Auto; same support rules as MLMG.
    {
        expect("explicit HYPRE request (neumann)", make(false, pb::BC::Neumann), pb::Status::Ok, K::HYPRE, K::HYPRE);
        expect("explicit HYPRE request (periodic)", make(true, pb::BC::Periodic), pb::Status::Ok, K::HYPRE, K::HYPRE);
        expect("explicit HYPRE request (dirichlet)", make(false, pb::BC::Dirichlet), pb::Status::Ok, K::HYPRE, K::HYPRE);
        pb::PressureProblem p = make(false, pb::BC::Neumann); p.nlevels = 2;
        expect("composite (2 levels) with explicit HYPRE and no level data not built", p, pb::Status::NotBuilt, K::Auto, K::HYPRE);
        p = make(false, pb::BC::Neumann); p.bc[pb::face_index(0,1)] = pb::BC::Dirichlet;
        expect("mixed open/closed faces with explicit HYPRE", p, pb::Status::Ok, K::HYPRE, K::HYPRE);
        p = make(false, pb::BC::Neumann);
        iMultiFab cls(p.ba, p.dm, 1, 0); cls.setVal(0);
        for (MFIter mfi(cls); mfi.isValid(); ++mfi) { if (mfi.index() == 0) { const Box vb = mfi.validbox(); const IntVect sm = vb.smallEnd(); cls[mfi].setVal<RunOn::Host>(1, Box(sm, sm)); } }
        p.cell_class = &cls;
        expect("masked cell with explicit HYPRE not built", p, pb::Status::NotBuilt, K::Auto, K::HYPRE);
        p = make(false, pb::BC::Neumann);
        iMultiFab unc(p.ba, p.dm, 1, 0); unc.setVal(1);
        for (MFIter mfi(unc); mfi.isValid(); ++mfi) { if (mfi.index() == 0) { const Box vb = mfi.validbox(); const IntVect sm = vb.smallEnd(); unc[mfi].setVal<RunOn::Host>(0, Box(sm, sm)); } }
        p.uncovered = &unc;
        expect("covered cell with explicit HYPRE not built", p, pb::Status::NotBuilt, K::Auto, K::HYPRE);
        p = make(false, pb::BC::Neumann); p.cylindrical = true;
        expect("cylindrical geometry with explicit HYPRE not built", p, pb::Status::NotBuilt, K::Auto, K::HYPRE);
        p = make(false, pb::BC::Neumann);
        for (int d = 0; d < 3; ++d) { p.cell_width[d].assign(8, Real(0.125)); }
        p.cell_width[2][3] = Real(0.1875);
        expect("non-uniform cell widths with explicit HYPRE not built", p, pb::Status::NotBuilt, K::Auto, K::HYPRE);
        p = make(false, pb::BC::Neumann); MultiFab a(p.ba, p.dm, 1, 0); a.setVal(1.0);
        p.cell_coef_a = &a;
        expect("cell coefficient a with explicit HYPRE not built", p, pb::Status::NotBuilt, K::Auto, K::HYPRE);
    }
}

// ---------------------------------------------------------------------------------------------
// Deterministic wide-range field, a function of the global index only.
double wide_value (int i, int j, int k)
{
    std::uint64_t h = 1469598103934665603ULL;
    for (int v : {i, j, k}) { h ^= std::uint64_t(v + 1000); h *= 1099511628211ULL; h ^= h >> 29; }
    const double m = double(h >> 11) / 9007199254740992.0 - 0.5;
    const int e = int((h >> 3) % 40) - 20;
    return std::ldexp(m, e) + 1.0e-3*std::sin(0.1*i + 0.2*j + 0.3*k);
}

void run_exactsum (ParmParse& pp)
{
    Setup s = make_setup(pp);
    MultiFab f(s.ba, s.dm, 1, 0);
    iMultiFab unc(s.ba, s.dm, 1, 0), lab(s.ba, s.dm, 1, 0);
    for (MFIter mfi(f); mfi.isValid(); ++mfi) {
        auto const& a = f.array(mfi); auto const& u = unc.array(mfi); auto const& l = lab.array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            a(i,j,k) = wide_value(i, j, k);
            u(i,j,k) = ((i + 2*j + 3*k) % 7 == 0) ? 0 : 1;      // "covered" pattern
            l(i,j,k) = (i < s.domain.length(0)/2) ? 0 : ((j % 5 == 0) ? -1 : 1);   // two components, some excluded
        });
    }
    pb::ExactSumResult r = pb::exact_sum(f, 0, 0.125, 2, &unc, &lab);
    Print() << "EXACTSUM sum0=" << hexd(r.sum[0]) << " sum1=" << hexd(r.sum[1]) << " n0=" << r.count[0]
            << " n1=" << r.count[1] << "\n";
    // Reference in long double (not the bitwise claim; a sanity bound only).
    long double ref0 = 0, ref1 = 0;
    for (int k = 0; k <= s.domain.bigEnd(2); ++k) for (int j = 0; j <= s.domain.bigEnd(1); ++j)
        for (int i = 0; i <= s.domain.bigEnd(0); ++i) {
            if ((i + 2*j + 3*k) % 7 == 0) { continue; }
            const int c = (i < s.domain.length(0)/2) ? 0 : ((j % 5 == 0) ? -1 : 1);
            if (c == 0) { ref0 += 0.125L*wide_value(i,j,k); } else if (c == 1) { ref1 += 0.125L*wide_value(i,j,k); }
        }
    check(std::abs(double(ref0) - r.sum[0]) <= 1e-13*std::abs(double(ref0)) + 1e-300, "exact sum component 0 matches long double sum");
    check(std::abs(double(ref1) - r.sum[1]) <= 1e-13*std::abs(double(ref1)) + 1e-300, "exact sum component 1 matches long double sum");

    // Mean removal through the common layer: two singular components, mask-aware, idempotent.
    pb::ComponentMap cm;
    cm.label.define(s.ba, s.dm, 1, 0);
    cm.label.ParallelCopy(lab, 0, 0, 1);
    for (int c = 0; c < 2; ++c) { pb::ComponentInfo ci; ci.id = c; ci.singular = true; cm.comps.push_back(ci); }
    pb::remove_mean(f, cm, &unc, Real(0.125));
    std::vector<double> v1 = gather(f, s.domain);
    Print() << "MEANREMOVAL hash1=" << hex64(fnv1a(v1)) << " removed0=" << hexd(cm.comps[0].removed_mean)
            << " removed1=" << hexd(cm.comps[1].removed_mean) << "\n";
    pb::ExactSumResult after = pb::exact_sum(f, 0, 0.125, 2, &unc, &lab);
    const double scale = 0.125 * f.norm0(0) * double(s.domain.numPts());
    check(std::abs(after.sum[0]) <= 1e-15*scale && std::abs(after.sum[1]) <= 1e-15*scale, "mean after removal is at round-off");
    pb::remove_mean(f, cm, &unc, Real(0.125));
    std::vector<double> v2 = gather(f, s.domain);
    Print() << "MEANREMOVAL hash2=" << hex64(fnv1a(v2)) << "\n";
    check(ParallelDescriptor::IOProcessor() ? (v1 == v2) : true, "mean removal is idempotent (bitwise)");
}


// ---------------------------------------------------------------------------------------------
// Per-cell volumes of a stretched mesh (varies with j and k): v = dx^3 * f(j,k).
Real vol_func (Setup const& s, int j, int k)
{
    const Real nz = Real(s.domain.length(2));
    return s.dx*s.dx*s.dx*(Real(0.5) + (k + Real(0.5))/nz)*(Real(1.0) + Real(0.3)*std::sin(Real(j)));
}

void run_meankind (ParmParse& pp)
{
    Setup s = make_setup(pp);
    std::string out; pp.get("out", out);
    double offset = 0.5; pp.query("rhs_offset", offset);
    MultiFab b0(s.ba, s.dm, 1, 0), v(s.ba, s.dm, 1, 0), rho(s.ba, s.dm, 1, 0), kres(s.ba, s.dm, 1, 0), phi0(s.ba, s.dm, 1, 0);
    for (MFIter mfi(b0); mfi.isValid(); ++mfi) {
        auto const& b = b0.array(mfi); auto const& vv = v.array(mfi); auto const& r = rho.array(mfi);
        auto const& kr = kres.array(mfi); auto const& ph = phi0.array(mfi);
        const Real dx = s.dx;
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const Real x = (i+0.5)*dx, y = (j+0.5)*dx, z = (k+0.5)*dx;
            b(i,j,k) = rhs_func(x, y, z) + offset;     // incompatible: nonzero mean
            vv(i,j,k) = vol_func(s, j, k);
            r(i,j,k) = rho_func(x, y, z); kr(i,j,k) = kres_func(x, y, z);
            ph(i,j,k) = std::cos(Real(3.)*x + y) + z*z + Real(0.7);   // arbitrary field to gauge
        });
    }
    pb::PressureProblem p;
    p.ba = s.ba; p.dm = s.dm; p.geom = s.geom; p.bc = s.bc;
    write_raw(gather(b0, s.domain), out + "_b0.bin");
    write_raw(gather(v, s.domain), out + "_v.bin");
    write_raw(gather(rho, s.domain), out + "_rho.bin");
    write_raw(gather(kres, s.domain), out + "_kres.bin");
    write_raw(gather(phi0, s.domain), out + "_phi0.bin");
    struct K { const char* name; pb::MeanKind kind; };
    for (K k : {K{"V", pb::MeanKind::Volume}, K{"S", pb::MeanKind::ScaledArithmetic}}) {
        MultiFab b(s.ba, s.dm, 1, 0); MultiFab::Copy(b, b0, 0, 0, 1, 0);
        pb::ComponentMap cm = pb::label_components(p);
        pb::remove_mean(b, cm, nullptr, 0.0, &v, k.kind);
        std::vector<double> v1 = gather(b, s.domain);
        Print() << "MEANKIND " << k.name << " removed_mean=" << hexd(cm.comps[0].removed_mean)
                << " removed_rel=" << cm.comps[0].removed_rel << " hash=" << hex64(fnv1a(v1)) << "\n";
        write_raw(v1, out + "_b" + k.name + ".bin");
        pb::remove_mean(b, cm, nullptr, 0.0, &v, k.kind);
        std::vector<double> v2 = gather(b, s.domain);
        check(ParallelDescriptor::IOProcessor() ? (v1 == v2) : true, std::string("mean removal idempotent (bitwise), kind ") + k.name);
    }
    {   // gauge with volume field only, and with rho*volume and KRES
        pb::ComponentMap cm = pb::label_components(p);
        MultiFab g1(s.ba, s.dm, 1, 0); MultiFab::Copy(g1, phi0, 0, 0, 1, 0);
        pb::apply_gauge(g1, cm, nullptr, 0.0, &v);
        write_raw(gather(g1, s.domain), out + "_gV.bin");
        MultiFab g2(s.ba, s.dm, 1, 0); MultiFab::Copy(g2, phi0, 0, 0, 1, 0);
        pb::apply_gauge(g2, cm, nullptr, 0.0, &v, &rho, &kres);
        write_raw(gather(g2, s.domain), out + "_gRK.bin");
        Print() << "MEANKIND gauge shifts V=" << hexd(0.0) << " RK=" << hexd(cm.comps[0].gauge_shift) << "\n";
    }
}

} // namespace

int run_composite_mode (std::string const& mode, ParmParse& pp);   // composite_modes.cpp
int run_m2_mode (std::string const& mode, ParmParse& pp);          // m2_modes.cpp

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ParmParse pp;
        std::string mode = "solve"; pp.query("mode", mode);
        if (mode == "solve") { run_solve(pp); }
        else if (mode == "gen") { run_gen(pp); }
        else if (mode == "diff") { run_diff(pp); }
        else if (mode == "selector") { run_selector(); }
        else if (mode == "exactsum") { run_exactsum(pp); }
        else if (mode == "meankind") { run_meankind(pp); }
        else {
            int cf = run_m2_mode(mode, pp);
            if (cf < 0) { cf = run_composite_mode(mode, pp); }
            if (cf < 0) { amrex::Abort("mode must be solve|gen|diff|selector|exactsum|meankind|comp|comp_ns2d|comp_sel|comp_ws|comp_gauge|comp_trigger|trigger1|fftcache|..."); }
            g_fail += cf;
        }
    }
    int fails = g_fail;
    amrex::Finalize();
    return fails == 0 ? 0 : 1;
}
