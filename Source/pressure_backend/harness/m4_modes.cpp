// M4: FDS-derived frozen cases. Separate executable pb_fds_frozen (own main; the shared pb_harness main.cpp is not touched):
//   pb_fds_frozen mode=fds_frozen  prefix=<path prefix of one FDS dump>  [backends="fft mlmg hypre"] [max_grid_size=16] [eps_H=1e-8]
//
// Reads one stage of one mesh written by the FDS-side dump (upstream patch 0006, format fds-pdump-1; see
// frozen/fds_csmag32_periodic/README.md): <prefix>meta.txt, <prefix>rhs.bin (FDS PRHS, double, I fastest, valid cells) and
// <prefix>phi.bin (FDS H or HS after its own Poisson solve). The problem is rebuilt from the metadata (cell counts, extents,
// face types), solved with each requested backend on the dumped right-hand side, and every H is compared with FDS's H after
// removing the volume-weighted mean of both (uniform cells: the arithmetic mean), because the additive constant of H is a
// gauge. The verdict is the relative L2 difference against eps_H (default 1e-8). Single mesh, unmasked, uniform cells only.
// CMP lines: <backend>_vs_fds and <backend>_vs_<first backend>; RESIDUAL line: residual of FDS's own H in the 7-point operator
// of the interface (the discrete operator the backends solve), relative to max|rhs - mean|.
#include "PressureIface.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>

using namespace amrex;

namespace {

int g_m4fail = 0;
void m4check (bool ok, std::string const& what)
{
    Print() << "CHECK " << (ok ? "PASS " : "FAIL ") << what << "\n";
    if (!ok) { ++g_m4fail; }
}

struct Meta {
    std::map<std::string, std::vector<std::string>> kv;
    std::vector<std::string> const& get (std::string const& k) const
    {
        auto it = kv.find(k);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(it != kv.end(), "fds_frozen: key missing in meta file: " + k);
        return it->second;
    }
};

Meta read_meta (std::string const& fname)
{
    Meta m;
    std::ifstream ifs(fname);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ifs.good(), "fds_frozen: cannot open " + fname);
    std::string line;
    while (std::getline(ifs, line)) {
        const auto eq = line.find('=');
        if (eq == std::string::npos) { continue; }
        std::string key = line.substr(0, eq);
        key.erase(std::remove(key.begin(), key.end(), ' '), key.end());
        std::istringstream is(line.substr(eq + 1));
        std::vector<std::string> v; std::string t;
        while (is >> t) { v.push_back(t); }
        m.kv[key] = v;
    }
    return m;
}

std::vector<double> read_doubles (std::string const& fname, std::size_t n)
{
    std::ifstream ifs(fname, std::ios::binary);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ifs.good(), "fds_frozen: cannot open " + fname);
    std::vector<double> v(n);
    ifs.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(n*sizeof(double)));
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ifs.good(), "fds_frozen: short read " + fname);
    return v;
}

void scatter (std::vector<double> const& v, Box const& domain, MultiFab& mf)
{
    BoxArray ba1(domain);
    DistributionMapping dm1(Vector<int>{ParallelDescriptor::IOProcessorNumber()});
    MultiFab g(ba1, dm1, 1, 0);
    if (ParallelDescriptor::IOProcessor()) { std::copy(v.begin(), v.end(), g[0].dataPtr()); }
    mf.ParallelCopy(g, 0, 0, 1);
}

std::vector<double> gather (MultiFab const& mf, Box const& domain)
{
    BoxArray ba1(domain);
    DistributionMapping dm1(Vector<int>{ParallelDescriptor::IOProcessorNumber()});
    MultiFab g(ba1, dm1, 1, 0);
    g.ParallelCopy(mf, 0, 0, 1);
    std::vector<double> v;
    if (ParallelDescriptor::IOProcessor()) { v.assign(g[0].dataPtr(), g[0].dataPtr() + domain.numPts()); }
    return v;
}

double mean_of (std::vector<double> const& v)
{
    long double s = 0.0L;
    for (double x : v) { s += x; }
    return double(s / static_cast<long double>(v.size()));
}

struct Diff { double rel_l2, rel_max, max_abs; };
// a and b with their (arithmetic = volume-weighted, uniform cells) means removed; relative to b.
Diff mean_free_diff (std::vector<double> const& a, std::vector<double> const& b)
{
    const double ma = mean_of(a), mb = mean_of(b);
    long double d2 = 0, b2 = 0; double dm = 0, bm = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const double x = (a[i] - ma) - (b[i] - mb), y = b[i] - mb;
        d2 += static_cast<long double>(x)*x; b2 += static_cast<long double>(y)*y;
        dm = std::max(dm, std::abs(x)); bm = std::max(bm, std::abs(y));
    }
    return {double(std::sqrt(d2/b2)), dm/bm, dm};
}

pb::BC parse_face (std::string const& s)
{
    if (s == "periodic") { return pb::BC::Periodic; }
    if (s == "neumann") { return pb::BC::Neumann; }
    if (s == "dirichlet") { return pb::BC::Dirichlet; }
    amrex::Abort("fds_frozen: face type '" + s + "' is not supported (periodic|neumann|dirichlet)");
    return pb::BC::Neumann;
}

void run_fds_frozen (ParmParse& pp)
{
    std::string prefix; pp.get("prefix", prefix);
    std::vector<std::string> backends;
    if (!pp.queryarr("backends", backends) || backends.empty()) { backends = {"fft", "mlmg", "hypre"}; }
    int mgs = 16; pp.query("max_grid_size", mgs);
    double eps_H = 1.0e-8; pp.query("eps_H", eps_H);

    const Meta meta = read_meta(prefix + "meta.txt");
    AMREX_ALWAYS_ASSERT(meta.get("format")[0] == "fds-pdump-1");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::stoi(meta.get("nmeshes")[0]) == 1, "fds_frozen: single-mesh cases only");
    const int nx = std::stoi(meta.get("ijk")[0]), ny = std::stoi(meta.get("ijk")[1]), nz = std::stoi(meta.get("ijk")[2]);
    const auto& xb = meta.get("xb");
    const double lx = std::stod(xb[1]) - std::stod(xb[0]), ly = std::stod(xb[3]) - std::stod(xb[2]), lz = std::stod(xb[5]) - std::stod(xb[4]);
    const auto& fb = meta.get("face_bc_xmin_xmax_ymin_ymax_zmin_zmax");
    AMREX_ALWAYS_ASSERT(fb.size() == 6);
    std::array<pb::BC,6> bc;
    for (int dir = 0; dir < 3; ++dir) { for (int side = 0; side < 2; ++side) { bc[pb::face_index(dir, side)] = parse_face(fb[2*dir + side]); } }
    for (int dir = 0; dir < 3; ++dir) {   // periodic is a pair
        AMREX_ALWAYS_ASSERT((bc[pb::face_index(dir,0)] == pb::BC::Periodic) == (bc[pb::face_index(dir,1)] == pb::BC::Periodic));
    }
    {   // inhomogeneous Poisson boundary data would need the S6 extension: the case must have none
        double bm = 0.0; for (auto const& s : meta.get("bc_data_maxabs_xs_xf_ys_yf_zs_zf")) { bm = std::max(bm, std::stod(s)); }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(bm == 0.0, "fds_frozen: nonzero Poisson boundary data (BXS..BZF)");
    }
    {   // uniform cell widths
        for (int dir = 0; dir < 3; ++dir) {
            const auto& c = meta.get(dir == 0 ? "x" : dir == 1 ? "y" : "z");
            const int n = (dir == 0 ? nx : dir == 1 ? ny : nz);
            AMREX_ALWAYS_ASSERT(int(c.size()) == n + 1);
            const double h0 = std::stod(c[1]) - std::stod(c[0]);
            for (int i = 1; i <= n; ++i) {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::abs((std::stod(c[i]) - std::stod(c[i-1])) - h0) < 1.0e-9*h0, "fds_frozen: non-uniform cells are not built");
            }
        }
    }
    const std::size_t ncell = std::size_t(nx)*ny*nz;
    const std::vector<double> rhs_v = read_doubles(prefix + "rhs.bin", ncell), fds_v = read_doubles(prefix + "phi.bin", ncell);

    Box domain(IntVect(0), IntVect(nx-1, ny-1, nz-1));
    RealBox rb({std::stod(xb[0]), std::stod(xb[2]), std::stod(xb[4])}, {std::stod(xb[1]), std::stod(xb[3]), std::stod(xb[5])});
    Array<int,3> per{bc[0] == pb::BC::Periodic, bc[1] == pb::BC::Periodic, bc[2] == pb::BC::Periodic};
    Geometry geom(domain, rb, CoordSys::cartesian, per);
    BoxArray ba(domain); ba.maxSize(mgs);
    DistributionMapping dm(ba);
    MultiFab rhs(ba, dm, 1, 0);
    scatter(rhs_v, domain, rhs);

    Print() << "FDS_FROZEN prefix=" << prefix << " ijk=" << nx << "x" << ny << "x" << nz << " L=" << lx << "," << ly << "," << lz
            << " bc=" << fb[0] << "," << fb[1] << "," << fb[2] << "," << fb[3] << "," << fb[4] << "," << fb[5]
            << " step=" << meta.get("step")[0] << " stage=" << meta.get("stage")[0] << " pres_flag=" << meta.get("pres_flag")[0]
            << " nranks=" << ParallelDescriptor::NProcs() << " nboxes=" << ba.size() << " eps_H=" << eps_H << "\n";

    // The residual of FDS's own H in the 7-point interface operator (rank 0, serial). Ghost rule: N p, D -p, P wrap.
    if (ParallelDescriptor::IOProcessor()) {
        const double h[3] = {lx/nx, ly/ny, lz/nz};
        const int n[3] = {nx, ny, nz};
        const double mrhs = mean_of(rhs_v);
        auto at = [&](int i, int j, int k) { return std::size_t(i) + std::size_t(nx)*(j + std::size_t(ny)*k); };
        auto val = [&](int i, int j, int k, int dir, int side) {   // value of FDS H one cell beyond the domain
            int q[3] = {i, j, k};
            const pb::BC b = bc[pb::face_index(dir, side)];
            if (b == pb::BC::Periodic) { q[dir] = side ? 0 : n[dir]-1; return fds_v[at(q[0], q[1], q[2])] * 1.0; }
            const double p = fds_v[at(i, j, k)];
            return b == pb::BC::Neumann ? p : -p;
        };
        double rmax = 0.0, fmax = 0.0;
        for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i) {
            const double p = fds_v[at(i,j,k)];
            double L = 0.0;
            const int ijk[3] = {i, j, k};
            for (int dir = 0; dir < 3; ++dir) {
                int lo[3] = {i, j, k}, hi[3] = {i, j, k};
                lo[dir] -= 1; hi[dir] += 1;
                const double pl = (ijk[dir] > 0)       ? fds_v[at(lo[0], lo[1], lo[2])] : val(i, j, k, dir, 0);
                const double ph = (ijk[dir] < n[dir]-1) ? fds_v[at(hi[0], hi[1], hi[2])] : val(i, j, k, dir, 1);
                L += (pl - 2.0*p + ph)/(h[dir]*h[dir]);
            }
            rmax = std::max(rmax, std::abs(L - (rhs_v[at(i,j,k)] - mrhs)));
            fmax = std::max(fmax, std::abs(rhs_v[at(i,j,k)] - mrhs));
        }
        Print() << std::setprecision(6) << "RESIDUAL fds_H_in_interface_operator max|L H - (rhs - mean)|/max|rhs - mean| = " << rmax/fmax
                << " (rhs mean " << mrhs << ", max|rhs - mean| " << fmax << ")\n";
    }

    std::vector<std::pair<std::string, std::vector<double>>> sol;
    for (auto const& name : backends) {
        pb::PressureProblem p;
        p.ba = ba; p.dm = dm; p.geom = geom; p.bc = bc; p.rhs = &rhs;
        MultiFab phi(ba, dm, 1, 1); phi.setVal(0.0); p.phi = &phi;
        pb::PressureOptions o;
        o.backend = name == "fft" ? pb::BackendKind::FFT : name == "mlmg" ? pb::BackendKind::MLMG : pb::BackendKind::HYPRE;
        const pb::PressureResult r = pb::solve_pressure(p, o);
        Print() << std::setprecision(6) << "RESULT " << name << " status=" << pb::to_string(r.status) << " backend=" << r.backend
                << " iters=" << r.backend_status.iterations << " true_rel2=" << r.residual_rel2 << " true_relmax=" << r.residual_relmax;
        for (auto const& c : r.components) { Print() << " comp" << c.id << "_singular=" << (c.singular ? 1 : 0) << " removed_rel=" << c.removed_rel; }
        Print() << "\n";
        m4check(r.status == pb::Status::Ok, name + " status Ok");
        if (r.status != pb::Status::Ok) { continue; }
        sol.emplace_back(name, gather(phi, domain));
    }
    auto verdict = [&](std::string const& label, Diff const& d) {
        Print() << std::setprecision(6) << "CMP " << label << " rel_l2=" << d.rel_l2 << " rel_max=" << d.rel_max << " max_abs=" << d.max_abs
                << " eps_H=" << eps_H << " " << (d.rel_l2 <= eps_H ? "PASS" : "FAIL") << "\n";
        if (!(d.rel_l2 <= eps_H)) { ++g_m4fail; }
    };
    if (ParallelDescriptor::IOProcessor()) {
        for (auto const& s : sol) { verdict(s.first + "_vs_fds", mean_free_diff(s.second, fds_v)); }
        for (std::size_t i = 1; i < sol.size(); ++i) { verdict(sol[i].first + "_vs_" + sol[0].first, mean_free_diff(sol[i].second, sol[0].second)); }
    }
}

} // namespace

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ParmParse pp;
        std::string mode = "fds_frozen"; pp.query("mode", mode);
        if (mode == "fds_frozen") { run_fds_frozen(pp); }
        else { amrex::Abort("pb_fds_frozen: mode must be fds_frozen"); }
    }
    const int fails = g_m4fail;
    amrex::Finalize();
    return fails == 0 ? 0 : 1;
}
