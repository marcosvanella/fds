// Harness modes for the composite (multi-level) pressure solve. Declared in main.cpp.
//
// mode=comp       manufactured-solution test on a refined-patch hierarchy (levels 1..3, ratio 2 or 4, optional
//                 full-domain fine level, optional 2-D plane with a one-cell y): composite solve, errors against
//                 the exact solution, comparison with the uniform single-level solve at the finest resolution,
//                 decomposition and workspace-reuse checks, face gradients and the divergence consistency check.
// mode=comp_ns2d the Role 1 ns2d_16 two-level hierarchy (16x1x16 periodic, central 8x8 coarse cells refined by
//                 (2,1,2)), geometry replicated from the case files (no driver code).
// mode=comp_sel   composite selector / validation checks (NotBuilt and InvalidInput messages).
// mode=comp_ws    workspace rebuild checks (D-058).
#include "PressureIface.H"
#include "ExactSum.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>

using namespace amrex;

namespace {

int g_cfail = 0;
void ccheck (bool ok, std::string const& what)
{
    Print() << "CHECK " << (ok ? "PASS " : "FAIL ") << what << "\n";
    if (!ok) { ++g_cfail; }
}

pb::BC cparse_bc (std::string const& s)
{
    if (s == "neumann") { return pb::BC::Neumann; }
    if (s == "periodic") { return pb::BC::Periodic; }
    if (s == "dirichlet") { return pb::BC::Dirichlet; }
    amrex::Abort("bc must be neumann|periodic|dirichlet");
    return pb::BC::Neumann;
}

const Real kPi = Real(3.141592653589793238462643383279502884);

std::uint64_t fnv (std::vector<double> const& v, std::uint64_t h = 1469598103934665603ULL)
{
    auto const* p = reinterpret_cast<const unsigned char*>(v.data());
    for (std::size_t i = 0; i < v.size()*sizeof(double); ++i) { h ^= p[i]; h *= 1099511628211ULL; }
    return h;
}
std::string hx (std::uint64_t h) { std::ostringstream o; o << std::hex << std::setw(16) << std::setfill('0') << h; return o.str(); }

// Gather the valid cells of a (possibly partial or face-centred) MultiFab into one Fab on rank 0 (x fastest); cells
// not in the BoxArray are 0.
std::vector<double> gather_box (MultiFab const& mf, Box const& box)
{
    BoxArray ba1(box);
    DistributionMapping dm1(Vector<int>{ParallelDescriptor::IOProcessorNumber()});
    MultiFab g(ba1, dm1, 1, 0);
    g.setVal(0.0);
    g.ParallelCopy(mf, 0, 0, 1);
    std::vector<double> v;
    if (ParallelDescriptor::IOProcessor()) { v.assign(g[0].dataPtr(), g[0].dataPtr() + box.numPts()); }
    return v;
}

// --------------------------------------------------------------------------------------------------------------
// Hierarchy description
// --------------------------------------------------------------------------------------------------------------
struct Cfg {
    int n = 32;                    // coarse cells per refined direction
    int nlev = 2;
    int ratio = 2;
    bool full = false;             // level 1 covers the whole domain
    bool plane2d = false;          // one cell in y (ratio 1 in y)
    pb::BC bc = pb::BC::Neumann;
    std::array<pb::BC,3> bcd = {pb::BC::Neumann, pb::BC::Neumann, pb::BC::Neumann};   // per direction (set from bc, or bcs3 for a neumann/periodic mix)
    int mgs = 16;                  // max grid size on every level
    int kx = 2, ky = 1, kz = 2;    // mode numbers of the manufactured solution
    double ly = 0.0;               // physical y extent (plane2d); 0 = dx*0.7
    void set_bc (pb::BC b) { bc = b; bcd = {b, b, b}; }
    int layout = 0;                // level 1 patch: 0 middle half, 1 corner at the low end of the domain (touches the domain faces, wraps when periodic), 2 two separate patches (nlev = 2)
};

struct Hier {
    std::vector<Geometry> geom;
    std::vector<BoxArray> ba;
    std::vector<DistributionMapping> dm;
    std::vector<IntVect> ratio;                         // to the coarser level
    std::vector<std::unique_ptr<MultiFab>> rhs, phi;
    std::vector<std::unique_ptr<iMultiFab>> unc;        // independent covered mask (harness's own)
};

struct Modes { double lambda = 0.0; };

Real exact_fn (Cfg const& c, Real x, Real y, Real z)
{
    const Real xs[3] = {x, y, z};
    const int ks[3] = {c.kx, c.ky, c.kz};
    Real r = 1.0;
    for (int d = 0; d < 3; ++d) {
        if (ks[d] == 0) { continue; }
        if (c.bcd[d] == pb::BC::Neumann) { r *= std::cos(kPi*ks[d]*xs[d]); }
        else if (c.bcd[d] == pb::BC::Dirichlet) { r *= std::sin(kPi*ks[d]*xs[d]); }
        else { const Real ph[3] = {0.3, 0.7, 0.1}; r *= std::sin(2*kPi*ks[d]*xs[d] + ph[d]); }
    }
    return r;
}
Real lambda_of (Cfg const& c)
{
    const int ks[3] = {c.kx, c.ky, c.kz};
    Real l = 0.0;
    for (int d = 0; d < 3; ++d) {
        const Real w = (c.bcd[d] == pb::BC::Periodic) ? 2*kPi : kPi;
        l += w*w*Real(ks[d]*ks[d]);
    }
    return l;
}
Real rhs_fn (Cfg const& c, Real x, Real y, Real z) { return -lambda_of(c) * exact_fn(c, x, y, z); }

DistributionMapping make_dm (BoxArray const& ba, int kind)
{
    if (kind == 0) { return DistributionMapping(ba); }
    auto const old = DistributionMapping::strategy();
    DistributionMapping::strategy(kind == 1 ? DistributionMapping::ROUNDROBIN : DistributionMapping::KNAPSACK);
    DistributionMapping dm(ba);
    DistributionMapping::strategy(old);
    return dm;
}

Hier build_hier (Cfg const& c, int dmkind)
{
    Hier h;
    const int ny = c.plane2d ? 1 : c.n;
    const int hd = c.plane2d ? 1 : -1;
    Box dom0(IntVect(0), IntVect(c.n-1, ny-1, c.n-1));
    const Real dx = Real(1.0)/c.n;
    const Real ly = c.plane2d ? (c.ly > 0 ? Real(c.ly) : Real(0.7)*dx) : Real(1.0);
    RealBox rb({0.,0.,0.}, {1.,ly,1.});
    Array<int,3> isp{c.bcd[0] == pb::BC::Periodic, c.bcd[1] == pb::BC::Periodic, c.bcd[2] == pb::BC::Periodic};
    Box patch, patch2;                                   // in the index space of the coarser level (patch2: layout 2 only)
    Box prev_patch = dom0;                               // extent of the coarser level
    for (int l = 0; l < c.nlev; ++l) {
        IntVect r = (l == 0) ? IntVect(1) : IntVect(c.ratio);
        if (hd >= 0 && l > 0) { r[hd] = 1; }
        Box dl = (l == 0) ? dom0 : amrex::refine(h.geom[l-1].Domain(), r);
        Box lev;
        if (l == 0) { lev = dl; }
        else {
            if (l == 1 && c.full) { patch = prev_patch; }
            else if (l == 1 && c.layout == 1) {
                IntVect lo = prev_patch.smallEnd(), hi = prev_patch.bigEnd();
                for (int d = 0; d < 3; ++d) { if (d == hd) { continue; } hi[d] = lo[d] + (hi[d] - lo[d] + 1)/2 - 1; }
                patch = Box(lo, hi);
            }
            else if (l == 1 && c.layout == 2) {
                IntVect lo = prev_patch.smallEnd(), hi = prev_patch.bigEnd();
                const int w = hi[0] - lo[0] + 1;
                IntVect lo1 = lo, hi1 = hi, lo2 = lo, hi2 = hi;
                for (int d = 0; d < 3; ++d) {
                    if (d == hd) { continue; }
                    const int wd = hi[d] - lo[d] + 1;
                    lo1[d] = lo[d] + wd/8; hi1[d] = lo[d] + 3*wd/8 - 1;
                    lo2[d] = lo[d] + 5*wd/8; hi2[d] = lo[d] + 7*wd/8 - 1;
                }
                (void)w;
                patch = Box(lo1, hi1);
                patch2 = Box(lo2, hi2);
            }
            else {
                IntVect lo = prev_patch.smallEnd(), hi = prev_patch.bigEnd();
                for (int d = 0; d < 3; ++d) {
                    if (d == hd) { continue; }
                    const int w = hi[d] - lo[d] + 1;
                    lo[d] += w/4; hi[d] -= w/4;
                }
                patch = Box(lo, hi);
            }
            lev = amrex::refine(patch, r);
        }
        BoxList extra;
        if (l == 1 && c.layout == 2 && !c.full) { extra.push_back(amrex::refine(patch2, r)); }
        h.geom.emplace_back(dl, rb, CoordSys::cartesian, isp);
        BoxList bl(lev); for (auto const& b : extra) { bl.push_back(b); }
        BoxArray ba(bl);
        ba.maxSize(c.mgs);
        h.ba.push_back(ba);
        h.dm.push_back(make_dm(ba, dmkind));
        h.ratio.push_back(r);
        prev_patch = (l == 0) ? dom0 : lev;
        // next level's extent (in this level's index space) is `lev`
    }
    for (int l = 0; l < c.nlev; ++l) {
        h.rhs.push_back(std::make_unique<MultiFab>(h.ba[l], h.dm[l], 1, 0));
        h.phi.push_back(std::make_unique<MultiFab>(h.ba[l], h.dm[l], 1, 1));
        h.phi[l]->setVal(0.0);
        auto const dxa = h.geom[l].CellSizeArray();
        auto const lo = h.geom[l].ProbLoArray();
        for (MFIter mfi(*h.rhs[l]); mfi.isValid(); ++mfi) {
            auto const& a = h.rhs[l]->array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                a(i,j,k) = rhs_fn(c, lo[0] + (i+0.5)*dxa[0], c.plane2d ? Real(0.5) : lo[1] + (j+0.5)*dxa[1], lo[2] + (k+0.5)*dxa[2]);
            });
        }
        // independent covered mask
        auto m = std::make_unique<iMultiFab>(h.ba[l], h.dm[l], 1, 0);
        m->setVal(1);
        if (l + 1 < c.nlev) {
            BoxArray cba = h.ba[l+1]; cba.coarsen(h.ratio[l+1]);
            for (MFIter mfi(*m); mfi.isValid(); ++mfi) {
                auto const& a = m->array(mfi);
                for (auto const& is : cba.intersections(mfi.validbox())) {
                    amrex::LoopOnCpu(is.second, [&] (int i, int j, int k) { a(i,j,k) = 0; });
                }
            }
        }
        h.unc.push_back(std::move(m));
    }
    return h;
}

pb::PressureProblem make_problem (Cfg const& c, Hier& h)
{
    pb::PressureProblem p;
    for (int d = 0; d < 3; ++d) { p.bc[pb::face_index(d,0)] = c.bcd[d]; p.bc[pb::face_index(d,1)] = c.bcd[d]; }
    for (int l = 0; l < c.nlev; ++l) {
        pb::PressureLevel L;
        L.ba = h.ba[l]; L.dm = h.dm[l]; L.geom = h.geom[l]; L.ref_ratio = h.ratio[l];
        L.rhs = h.rhs[l].get(); L.phi = h.phi[l].get();
        p.levels.push_back(L);
    }
    return p;
}

// The same one-level layout through the single-level API (FFT/MLMG selector path).
pb::PressureProblem make_single_problem (Cfg const& c, Hier& h)
{
    pb::PressureProblem p;
    for (int d = 0; d < 3; ++d) { p.bc[pb::face_index(d,0)] = c.bcd[d]; p.bc[pb::face_index(d,1)] = c.bcd[d]; }
    p.ba = h.ba[0]; p.dm = h.dm[0]; p.geom = h.geom[0]; p.rhs = h.rhs[0].get(); p.phi = h.phi[0].get();
    return p;
}

std::vector<double> hier_gather (Hier const& h, int nlev)
{
    std::vector<double> all;
    for (int l = 0; l < nlev; ++l) {
        auto v = gather_box(*h.phi[l], h.geom[l].Domain());
        all.insert(all.end(), v.begin(), v.end());
    }
    return all;
}

struct Err { double l2 = 0.0, linf = 0.0, mean = 0.0; };

// Error of `phi` against `ref` over the uncovered cells (volume weighted L2, max). Removes the volume-weighted
// mean of the difference first when `remove_const`.
template <class F>
Err hier_error (Cfg const& c, Hier const& h, std::vector<MultiFab*> const& phi, F ref, bool remove_const, int only_level = -1)
{
    double sv = 0.0, se = 0.0;
    for (int l = 0; l < c.nlev; ++l) {
        if (only_level >= 0 && l != only_level) { continue; }
        auto const dxa = h.geom[l].CellSizeArray();
        const double v = dxa[0]*dxa[1]*dxa[2];
        auto const lo = h.geom[l].ProbLoArray();
        for (MFIter mfi(*phi[l]); mfi.isValid(); ++mfi) {
            auto const& a = phi[l]->const_array(mfi);
            auto const& u = h.unc[l]->const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                if (!u(i,j,k)) { return; }
                sv += v; se += v*(a(i,j,k) - ref(l, i, j, k, lo, dxa));
            });
        }
    }
    ParallelDescriptor::ReduceRealSum(sv); ParallelDescriptor::ReduceRealSum(se);
    Err e;
    e.mean = remove_const ? se/sv : 0.0;
    double s2 = 0.0, mx = 0.0;
    for (int l = 0; l < c.nlev; ++l) {
        if (only_level >= 0 && l != only_level) { continue; }
        auto const dxa = h.geom[l].CellSizeArray();
        const double v = dxa[0]*dxa[1]*dxa[2];
        auto const lo = h.geom[l].ProbLoArray();
        for (MFIter mfi(*phi[l]); mfi.isValid(); ++mfi) {
            auto const& a = phi[l]->const_array(mfi);
            auto const& u = h.unc[l]->const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                if (!u(i,j,k)) { return; }
                const double d = a(i,j,k) - ref(l, i, j, k, lo, dxa) - e.mean;
                s2 += v*d*d; mx = std::max(mx, std::abs(d));
            });
        }
    }
    ParallelDescriptor::ReduceRealSum(s2); ParallelDescriptor::ReduceRealMax(mx);
    e.l2 = std::sqrt(s2/sv); e.linf = mx;
    return e;
}

// --------------------------------------------------------------------------------------------------------------
// Face-gradient checks
// --------------------------------------------------------------------------------------------------------------
struct GradOut {
    std::vector<std::array<std::unique_ptr<MultiFab>,3>> store;
    std::vector<std::array<MultiFab*,3>> ptr;
};

GradOut make_grad (Cfg const& c, Hier const& h)
{
    GradOut g;
    g.store.resize(c.nlev); g.ptr.resize(c.nlev);
    for (int l = 0; l < c.nlev; ++l) {
        for (int d = 0; d < 3; ++d) {
            g.store[l][d] = std::make_unique<MultiFab>(amrex::convert(h.ba[l], IntVect::TheDimensionVector(d)), h.dm[l], 1, 0);
            g.store[l][d]->setVal(0.0);
            g.ptr[l][d] = g.store[l][d].get();
        }
    }
    return g;
}

// Run all face-gradient checks (e). Returns nothing; reports through ccheck / GRAD lines.
void gradient_checks (Cfg const& c, Hier& h, pb::PressureProblem const& p)
{
    GradOut g = make_grad(c, h);
    pb::PressureResult gr = pb::face_gradient_composite(p, g.ptr);
    ccheck(gr.status == pb::Status::Ok, "face_gradient_composite status Ok");
    if (gr.status != pb::Status::Ok) { return; }
    const int hd = c.plane2d ? 1 : -1;
    // (1) coarse faces under the fine level equal the plain average of the fine faces (independent loops on gathered arrays)
    double worst_avg = 0.0;
    for (int l = 1; l < c.nlev; ++l) {
        for (int d = 0; d < 3; ++d) {
            if (d == hd) { continue; }
            Box fdom = amrex::surroundingNodes(h.geom[l].Domain(), d);
            Box cdom = amrex::surroundingNodes(h.geom[l-1].Domain(), d);
            auto fv = gather_box(*g.ptr[l][d], fdom);
            auto cv = gather_box(*g.ptr[l-1][d], cdom);
            // fine faces present
            MultiFab ones(amrex::convert(h.ba[l], IntVect::TheDimensionVector(d)), h.dm[l], 1, 0); ones.setVal(1.0);
            auto fm = gather_box(ones, fdom);
            if (ParallelDescriptor::IOProcessor()) {
                const IntVect r = h.ratio[l];
                const Long fnx = fdom.length(0), fny = fdom.length(1);
                const Long cnx = cdom.length(0), cny = cdom.length(1);
                for (int k = cdom.smallEnd(2); k <= cdom.bigEnd(2); ++k)
                for (int j = cdom.smallEnd(1); j <= cdom.bigEnd(1); ++j)
                for (int i = cdom.smallEnd(0); i <= cdom.bigEnd(0); ++i) {
                    // fine faces covering this coarse face: along d the single fine face at r[d]*index, tangentially r cells
                    IntVect lo(i*(d==0 ? r[0] : r[0]), j*r[1], k*r[2]);
                    IntVect hi = lo;
                    for (int t = 0; t < 3; ++t) {
                        if (t == d) { hi[t] = lo[t]; }
                        else { hi[t] = lo[t] + r[t] - 1; }
                    }
                    double s = 0.0; int cnt = 0; bool all = true;
                    for (int kk = lo[2]; kk <= hi[2]; ++kk) for (int jj = lo[1]; jj <= hi[1]; ++jj) for (int ii = lo[0]; ii <= hi[0]; ++ii) {
                        if (!fdom.contains(ii,jj,kk)) { all = false; continue; }
                        const Long fi = (ii - fdom.smallEnd(0)) + fnx*((jj - fdom.smallEnd(1)) + fny*Long(kk - fdom.smallEnd(2)));
                        if (fm[fi] == 0.0) { all = false; continue; }
                        s += fv[fi]; ++cnt;
                    }
                    if (!all || cnt == 0) { continue; }
                    const Long ci = (i - cdom.smallEnd(0)) + cnx*((j - cdom.smallEnd(1)) + cny*Long(k - cdom.smallEnd(2)));
                    worst_avg = std::max(worst_avg, std::abs(cv[ci] - s/cnt));
                }
            }
        }
    }
    ParallelDescriptor::ReduceRealMax(worst_avg);
    ccheck(worst_avg <= 1e-13, "coarse faces under a finer level equal the average of the fine faces (max dev " + std::to_string(worst_avg) + ")");

    // (2) divergence consistency: sum_d (g_hi - g_lo)/dx = rhs (mean-removed) on uncovered cells.
    //     The mean removed in the solve is subtracted from rhs here (it is a constant shift on uncovered cells).
    double bsum = 0.0, vsum = 0.0;
    for (int l = 0; l < c.nlev; ++l) {
        auto const dxa = h.geom[l].CellSizeArray(); const double v = dxa[0]*dxa[1]*dxa[2];
        for (MFIter mfi(*h.rhs[l]); mfi.isValid(); ++mfi) {
            auto const& a = h.rhs[l]->const_array(mfi); auto const& u = h.unc[l]->const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { if (u(i,j,k)) { bsum += v*a(i,j,k); vsum += v; } });
        }
    }
    ParallelDescriptor::ReduceRealSum(bsum); ParallelDescriptor::ReduceRealSum(vsum);
    const bool singular = (c.bc != pb::BC::Dirichlet);
    const double shift = singular ? bsum/vsum : 0.0;
    auto divergence_error = [&] (bool raw_coarse, double& rel2, double& relmax) {
        double r2 = 0.0, b2 = 0.0, rm = 0.0, bm = 0.0;
        for (int l = 0; l < c.nlev; ++l) {
            auto const dxa = h.geom[l].CellSizeArray(); const double v = dxa[0]*dxa[1]*dxa[2];
            // optionally replace the C/F-face values of the coarse gradient by the plain coarse difference (negative control)
            std::array<std::unique_ptr<MultiFab>,3> alt;
            if (raw_coarse && l + 1 < c.nlev) {
                MultiFab ph(h.ba[l], h.dm[l], 1, 1); ph.setVal(0.0); MultiFab::Copy(ph, *h.phi[l], 0, 0, 1, 0);
                ph.FillBoundary(h.geom[l].periodicity());
                for (int d = 0; d < 3; ++d) {
                    alt[d] = std::make_unique<MultiFab>(amrex::convert(h.ba[l], IntVect::TheDimensionVector(d)), h.dm[l], 1, 0);
                    for (MFIter mfi(*alt[d]); mfi.isValid(); ++mfi) {
                        auto const& a = alt[d]->array(mfi); auto const& q = ph.const_array(mfi);
                        const Box dom = h.geom[l].Domain();
                        const IntVect e = IntVect::TheDimensionVector(d);
                        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                            const IntVect iv(i,j,k);
                            // face value outside the domain: Neumann 0; otherwise (phi(iv) - phi(iv-e))/dx (ghost cells filled)
                            if (!dom.contains(iv) && !h.geom[l].isPeriodic(d)) { a(i,j,k) = 0.0; return; }
                            a(i,j,k) = (q(i,j,k) - q(i-e[0],j-e[1],k-e[2]))/dxa[d];
                        });
                    }
                }
            }
            for (MFIter mfi(*h.rhs[l]); mfi.isValid(); ++mfi) {
                auto const& b = h.rhs[l]->const_array(mfi); auto const& u = h.unc[l]->const_array(mfi);
                Array4<Real const> gx = (raw_coarse && alt[0]) ? alt[0]->const_array(mfi) : g.ptr[l][0]->const_array(mfi);
                Array4<Real const> gy = (raw_coarse && alt[1]) ? alt[1]->const_array(mfi) : g.ptr[l][1]->const_array(mfi);
                Array4<Real const> gz = (raw_coarse && alt[2]) ? alt[2]->const_array(mfi) : g.ptr[l][2]->const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    if (!u(i,j,k)) { return; }
                    double div = (gx(i+1,j,k) - gx(i,j,k))/dxa[0] + (gz(i,j,k+1) - gz(i,j,k))/dxa[2];
                    if (hd != 1) { div += (gy(i,j+1,k) - gy(i,j,k))/dxa[1]; }
                    const double bb = b(i,j,k) - shift;
                    r2 += v*(bb-div)*(bb-div); b2 += v*bb*bb; rm = std::max(rm, std::abs(bb-div)); bm = std::max(bm, std::abs(bb));
                });
            }
        }
        ParallelDescriptor::ReduceRealSum(r2); ParallelDescriptor::ReduceRealSum(b2);
        ParallelDescriptor::ReduceRealMax(rm); ParallelDescriptor::ReduceRealMax(bm);
        rel2 = std::sqrt(r2/b2); relmax = rm/bm;
    };
    double r2, rm;
    divergence_error(false, r2, rm);
    Print() << std::setprecision(6) << "GRAD div_consistency rel2=" << r2 << " relmax=" << rm << "\n";
    ccheck(r2 <= 1e-8, "face-gradient divergence (with fine-flux averaged coarse faces) equals the RHS on uncovered cells, rel2 " + std::to_string(r2));
    if (c.nlev > 1 && !c.full) {
        double r2c, rmc;
        divergence_error(true, r2c, rmc);
        Print() << std::setprecision(6) << "GRAD div_consistency_negative_control rel2=" << r2c << " relmax=" << rmc << "\n";
        ccheck(r2c > 100*std::max(r2, 1e-12), "negative control: plain coarse differences at the C/F faces violate the divergence consistency (test is sensitive)");
    }
    // (3) accuracy of the gradients against the analytic ones (second order)
    double gerr = 0.0, gmax = 0.0; int gloc[5] = {-1,-1,0,0,0};
    for (int l = 0; l < c.nlev; ++l) {
        auto const dxa = h.geom[l].CellSizeArray(); auto const lo = h.geom[l].ProbLoArray();
        const int ks[3] = {c.kx, c.ky, c.kz};
        for (int d = 0; d < 3; ++d) {
            if (d == hd) { continue; }
            for (MFIter mfi(*g.ptr[l][d]); mfi.isValid(); ++mfi) {
                auto const& a = g.ptr[l][d]->const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    const Real x = lo[0] + (d==0 ? i : i+0.5)*dxa[0];
                    const Real y = c.plane2d ? Real(0.5) : lo[1] + (d==1 ? j : j+0.5)*dxa[1];
                    const Real z = lo[2] + (d==2 ? k : k+0.5)*dxa[2];
                    // analytic derivative by a centred difference of the exact function (spacing small)
                    const Real e = Real(1e-5);
                    const Real xp = x + (d==0 ? e : 0), xm = x - (d==0 ? e : 0);
                    const Real yp = y + (d==1 ? e : 0), ym = y - (d==1 ? e : 0);
                    const Real zp = z + (d==2 ? e : 0), zm = z - (d==2 ? e : 0);
                    const double ex = (exact_fn(c, xp, yp, zp) - exact_fn(c, xm, ym, zm))/(2*e);
                    // Dirichlet/Neumann domain faces: the stencil gradient at a domain face is one-sided; skip faces on the boundary
                    const Box dom = h.geom[l].Domain();
                    const int idx = (d==0 ? i : (d==1 ? j : k));
                    if (!h.geom[l].isPeriodic(d) && (idx == dom.smallEnd(d) || idx == dom.bigEnd(d)+1)) { return; }
                    if (std::abs(a(i,j,k) - ex) > gerr) { gerr = std::abs(a(i,j,k) - ex); gloc[0] = l; gloc[1] = d; gloc[2] = i; gloc[3] = j; gloc[4] = k; }
                    gmax = std::max(gmax, std::abs(ex));
                });
            }
        }
        (void)ks;
    }
    ParallelDescriptor::ReduceRealMax(gerr); ParallelDescriptor::ReduceRealMax(gmax);
    Print() << std::setprecision(6) << "GRAD accuracy(rank-local loc l,d,i,j,k=" << gloc[0] << "," << gloc[1] << "," << gloc[2] << "," << gloc[3] << "," << gloc[4] << ") max_abs_err=" << gerr << " max_grad=" << gmax << " rel=" << gerr/gmax << "\n";
    // Sanity bound, not an order claim: gradients on the fine side of a C/F face come from the interpolated fine ghost value
    // (error ~ h_c^2/h_f), so the bound grows with the ratio; everywhere else the gradient is second order.
    const double gtol = 0.06*std::pow(c.ratio/2.0, 2);
    ccheck(gerr/gmax <= gtol, "face gradient agrees with the analytic gradient (rel max error " + std::to_string(gerr/gmax) + " <= " + std::to_string(gtol) + ")");
}

// --------------------------------------------------------------------------------------------------------------
// mode=comp
// --------------------------------------------------------------------------------------------------------------
void run_comp (ParmParse& pp)
{
    Cfg c;
    pp.query("n", c.n); pp.query("nlev", c.nlev); pp.query("ratio", c.ratio);
    int full = 0, plane = 0; pp.query("full", full); pp.query("plane2d", plane);
    c.full = full != 0; c.plane2d = plane != 0;
    std::string bcs = "neumann"; pp.query("bc", bcs); c.set_bc(cparse_bc(bcs));
    {   // optional per-direction mix of neumann and periodic: bcs3 = "periodic neumann periodic"
        std::vector<std::string> b3;
        if (pp.queryarr("bcs3", b3) && b3.size() == 3) { for (int d = 0; d < 3; ++d) { c.bcd[d] = cparse_bc(b3[d]); } bcs = b3[0] + "-" + b3[1] + "-" + b3[2]; c.bc = pb::BC::Neumann; }
    }
    pp.query("mgs", c.mgs); pp.query("kx", c.kx); pp.query("ky", c.ky); pp.query("kz", c.kz);
    if (c.plane2d) { c.ky = 0; }
    pp.query("ly", c.ly); pp.query("layout", c.layout);
    int dmkind = 0; pp.query("dmkind", dmkind);
    int do_uniform = 1; pp.query("uniform", do_uniform);
    int do_grad = 0; pp.query("grad", do_grad);
    int do_repeat = 0; pp.query("repeat", do_repeat);
    int mgs2 = 0; pp.query("mgs2", mgs2);
    int dmkind2 = 1; pp.query("dmkind2", dmkind2);
    std::string out; pp.query("out", out);
    double tol_rel = 1e-12; pp.query("tol_rel", tol_rel);
    double eps_H = 1e-8; pp.query("eps_H", eps_H);

    Hier h = build_hier(c, dmkind);
    pb::PressureProblem p = make_problem(c, h);
    pb::PressureOptions o;
    o.tol_rel = tol_rel; o.removed_mean_warn = 1.0; o.residual_tol = eps_H; o.verbose = 1; pp.query("verbose", o.verbose);
    pb::PressureWorkspace ws;
    pb::PressureResult r = pb::solve_pressure(p, o, &ws);
    Long ncells = 0; for (int l = 0; l < c.nlev; ++l) { ncells += h.ba[l].numPts(); }
    std::vector<MultiFab*> phis; for (auto& m : h.phi) { phis.push_back(m.get()); }
    const bool remove_const = (c.bc != pb::BC::Dirichlet);
    auto exact_ref = [&] (int l, int i, int j, int k, GpuArray<Real,3> const& lo, GpuArray<Real,3> const& dxa) {
        (void)l;
        return double(exact_fn(c, lo[0] + (i+0.5)*dxa[0], c.plane2d ? Real(0.5) : lo[1] + (j+0.5)*dxa[1], lo[2] + (k+0.5)*dxa[2]));
    };
    Err ee = hier_error(c, h, phis, exact_ref, remove_const);
    auto hv = hier_gather(h, c.nlev);
    Print() << std::setprecision(6) << "COMP n=" << c.n << " nlev=" << c.nlev << " ratio=" << c.ratio << " full=" << full
            << " plane2d=" << plane << " bc=" << bcs << " nranks=" << ParallelDescriptor::NProcs() << " mgs=" << c.mgs
            << " ncells=" << ncells << " status=" << pb::to_string(r.status) << " iters=" << r.backend_status.iterations
            << " own_res=" << r.backend_status.own_residual << " true_rel2=" << r.residual_rel2
            << " true_relmax=" << r.residual_relmax << " residual_ok=" << (r.residual_ok ? 1 : 0)
            << " nwarn=" << r.warnings.size() << " err_l2=" << ee.l2 << " err_linf=" << ee.linf
            << " removed_rel=" << (r.components.empty() ? 0.0 : r.components[0].removed_rel)
            << " nunc=" << r.ncells_uncovered << " hash=" << hx(fnv(hv)) << "\n";
    ccheck(r.status == pb::Status::Ok, "composite solve status Ok");
    ccheck(r.residual_checked && r.residual_rel2 <= eps_H, "composite true residual " + std::to_string(r.residual_rel2) + " <= eps_H");

    // covered coarse cells hold the average-down of the fine solution
    {
        double dev = 0.0;
        for (int l = c.nlev - 1; l >= 1; --l) {
            MultiFab tmp(h.ba[l-1], h.dm[l-1], 1, 0);
            tmp.setVal(0.0);
            amrex::average_down(*h.phi[l], tmp, 0, 1, h.ratio[l]);
            for (MFIter mfi(tmp); mfi.isValid(); ++mfi) {
                auto const& t = tmp.const_array(mfi); auto const& a = h.phi[l-1]->const_array(mfi); auto const& u = h.unc[l-1]->const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    if (!u(i,j,k)) { dev = std::max(dev, std::abs(double(t(i,j,k) - a(i,j,k)))); }
                });
            }
        }
        ParallelDescriptor::ReduceRealMax(dev);
        ccheck(dev == 0.0, "covered coarse cells equal the average-down of the next finer level (max dev " + std::to_string(dev) + ")");
    }

    if (!out.empty() && ParallelDescriptor::IOProcessor()) {
        std::ofstream ofs(out + "_phi.bin", std::ios::binary);
        ofs.write(reinterpret_cast<const char*>(hv.data()), static_cast<std::streamsize>(hv.size()*sizeof(double)));
    }

    // Uniform single-level reference at the finest resolution (FFT).
    if (do_uniform) {
        int fr = 1; for (int l = 1; l < c.nlev; ++l) { fr *= c.ratio; }
        Cfg u = c; u.nlev = 1; u.n = c.n*fr; u.mgs = c.mgs; u.full = false;
        Hier hu = build_hier(u, 0);
        pb::PressureProblem pu = make_single_problem(u, hu);
        pb::PressureOptions ou = o; ou.tol_rel = 1e-12; ou.residual_tol = 1e-12;
        pb::PressureResult ru = pb::solve_pressure(pu, ou);
        ccheck(ru.status == pb::Status::Ok, "uniform fine reference solve status Ok (backend " + ru.backend + ")");
        std::vector<MultiFab*> pus{hu.phi[0].get()};
        Err eu = hier_error(u, hu, pus, exact_ref, remove_const);
        Err ef = hier_error(c, h, phis, exact_ref, remove_const, c.nlev - 1);
        // composite finest level against the uniform solution, on the finest level's cells
        const int lf = c.nlev - 1;
        MultiFab ref(h.ba[lf], h.dm[lf], 1, 0);
        ref.setVal(0.0);
        ref.ParallelCopy(*hu.phi[0], 0, 0, 1);
        MultiFab d(h.ba[lf], h.dm[lf], 1, 0);
        MultiFab::Copy(d, *h.phi[lf], 0, 0, 1, 0);
        MultiFab::Subtract(d, ref, 0, 0, 1, 0);
        double dm_ = 0.0;
        if (remove_const) { dm_ = d.sum(0)/double(h.ba[lf].numPts()); d.plus(-dm_, 0, 1, 0); }
        const double dl2 = d.norm2(0)/std::sqrt(double(h.ba[lf].numPts())), dinf = d.norm0(0);
        const double rl2 = ref.norm2(0)/std::sqrt(double(h.ba[lf].numPts()));
        // composite error and the uniform-coarse error for context
        Cfg cc = c; cc.nlev = 1; cc.full = false;
        Hier hc = build_hier(cc, 0);
        pb::PressureProblem pc = make_single_problem(cc, hc);
        pb::PressureResult rc = pb::solve_pressure(pc, ou);
        std::vector<MultiFab*> pcs{hc.phi[0].get()};
        Err ec = hier_error(cc, hc, pcs, exact_ref, remove_const);
        Print() << std::setprecision(6) << "UNIFORM n=" << c.n << " nlev=" << c.nlev << " ratio=" << c.ratio << " full=" << full << " bc=" << bcs
                << " uni_fine_err_l2=" << eu.l2 << " uni_coarse_err_l2=" << ec.l2 << " comp_err_l2=" << ee.l2 << " comp_fine_err_l2=" << ef.l2 << " comp_fine_err_linf=" << ef.linf << " uni_backend=" << ru.backend
                << " diff_fine_l2_abs=" << dl2 << " diff_fine_l2_rel=" << dl2/rl2 << " diff_fine_linf=" << dinf
                << " mean_const=" << dm_ << " coarse_status=" << pb::to_string(rc.status) << "\n";
        ccheck(ef.l2 <= ec.l2, "composite error on the finest level is below the uniform coarse error");
    }

    // Decomposition independence in process: other box size, other distribution mapping.
    if (mgs2 > 0) {
        Cfg c2 = c; c2.mgs = mgs2;
        Hier h2 = build_hier(c2, dmkind2);
        pb::PressureProblem p2 = make_problem(c2, h2);
        pb::PressureResult r2 = pb::solve_pressure(p2, o);
        ccheck(r2.status == pb::Status::Ok, "second decomposition solve status Ok");
        double num = 0.0, den = 0.0;
        auto hv2 = hier_gather(h2, c.nlev);
        if (ParallelDescriptor::IOProcessor()) {
            for (std::size_t i = 0; i < hv.size(); ++i) { num += (hv[i]-hv2[i])*(hv[i]-hv2[i]); den += hv[i]*hv[i]; }
        }
        ParallelDescriptor::ReduceRealSum(num); ParallelDescriptor::ReduceRealSum(den);
        const double rel = std::sqrt(num/den);
        Print() << std::setprecision(6) << "DECOMP mgs=" << c.mgs << " mgs2=" << mgs2 << " nboxes=" << h.ba[0].size() << "," << h2.ba[0].size()
                << " rel_l2=" << rel << " mean_removed_bitwise=" << ((r.components[0].removed_mean == r2.components[0].removed_mean) ? 1 : 0)
                << " gauge_bitwise=" << ((r.components[0].gauge_shift == r2.components[0].gauge_shift) ? 1 : 0) << " hash2=" << hx(fnv(hv2)) << "\n";
        ccheck(rel <= eps_H, "decomposition independence (box split and distribution mapping) rel L2 " + std::to_string(rel) + " <= eps_H");
        ccheck(r.components[0].removed_mean == r2.components[0].removed_mean, "removed mean bitwise identical across decompositions");
    }

    // Workspace reuse and run-to-run repeatability.
    if (do_repeat) {
        pb::PressureResult r3 = pb::solve_pressure(p, o);                  // fresh build, discard
        auto hv3 = hier_gather(h, c.nlev);
        pb::PressureResult r4 = pb::solve_pressure(p, o, &ws);             // reused workspace
        auto hv4 = hier_gather(h, c.nlev);
        Print() << "REPEAT hash_first=" << hx(fnv(hv)) << " hash_fresh=" << hx(fnv(hv3)) << " hash_reuse=" << hx(fnv(hv4))
                << " rebuilt_flags=" << (r3.workspace_rebuilt ? 1 : 0) << (r4.workspace_rebuilt ? 1 : 0) << "\n";
        ccheck(r3.status == pb::Status::Ok && r4.status == pb::Status::Ok, "repeat solves Ok");
        ccheck(hv == hv3, "run-to-run bitwise repeatability (fresh setup)");
        ccheck(hv == hv4, "solve with a reused workspace bitwise equals the first solve");
        ccheck(!r4.workspace_rebuilt, "a matching workspace is not rebuilt");
    }

    if (do_grad) { gradient_checks(c, h, p); }
}

// --------------------------------------------------------------------------------------------------------------
// mode=comp_ns2d
// --------------------------------------------------------------------------------------------------------------
void run_comp_ns2d (ParmParse& pp)
{
    // ns2d_16_int_1to2: level 0 = 16x1x16 periodic over [0,2pi]x[-0.05,0.05]x[0,2pi]; level 1 = the mesh
    // IJK=16,1,16, XB=pi/2..3pi/2 in x and z, i.e. the central 8x8 coarse cells, ratio (2,1,2).
    const Real twopi = Real(6.28318530718);
    RealBox rb({0.0, -0.05, 0.0}, {twopi, 0.05, twopi});
    Array<int,3> isp{1, 1, 1};
    Box d0(IntVect(0,0,0), IntVect(15,0,15));
    IntVect r(2,1,2);
    Box d1 = amrex::refine(d0, r);
    int mgs = 8; pp.query("mgs", mgs);
    int nranks_dm = 0; pp.query("dmkind", nranks_dm);
    Geometry g0(d0, rb, CoordSys::cartesian, isp), g1(d1, rb, CoordSys::cartesian, isp);
    BoxArray ba0(d0); ba0.maxSize(mgs);
    BoxArray ba1(amrex::refine(Box(IntVect(4,0,4), IntVect(11,0,11)), r)); ba1.maxSize(mgs);
    DistributionMapping dm0(ba0), dm1(ba1);
    MultiFab rhs0(ba0, dm0, 1, 0), rhs1(ba1, dm1, 1, 0), phi0(ba0, dm0, 1, 1), phi1(ba1, dm1, 1, 1);
    phi0.setVal(0.0); phi1.setVal(0.0);
    // Taylor-Green like velocity divergence source: smooth, mean-free, with structure in the refined centre
    auto src = [&] (Real x, Real z) { return std::cos(x)*std::sin(2*z) + Real(0.3)*std::sin(2*x)*std::cos(z)*std::exp(-((x-kPi)*(x-kPi) + (z-kPi)*(z-kPi))/Real(2.0)); };
    for (int l = 0; l < 2; ++l) {
        MultiFab& m = (l == 0) ? rhs0 : rhs1;
        auto const dxa = (l == 0 ? g0 : g1).CellSizeArray();
        for (MFIter mfi(m); mfi.isValid(); ++mfi) {
            auto const& a = m.array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { a(i,j,k) = src((i+0.5)*dxa[0], (k+0.5)*dxa[2]); });
        }
    }
    pb::PressureProblem p;
    p.bc.fill(pb::BC::Periodic);
    pb::PressureLevel L0; L0.ba = ba0; L0.dm = dm0; L0.geom = g0; L0.rhs = &rhs0; L0.phi = &phi0;
    pb::PressureLevel L1; L1.ba = ba1; L1.dm = dm1; L1.geom = g1; L1.ref_ratio = r; L1.rhs = &rhs1; L1.phi = &phi1;
    p.levels = {L0, L1};
    pb::PressureOptions o; o.residual_tol = 1e-8; o.removed_mean_warn = 1.0; pp.query("verbose", o.verbose);
    pb::PressureWorkspace ws;
    pb::PressureResult res = pb::solve_pressure(p, o, &ws);
    Print() << std::setprecision(6) << "NS2D status=" << pb::to_string(res.status) << " iters=" << res.backend_status.iterations
            << " true_rel2=" << res.residual_rel2 << " true_relmax=" << res.residual_relmax << " nunc=" << res.ncells_uncovered
            << " removed_rel=" << res.components[0].removed_rel << " nranks=" << ParallelDescriptor::NProcs() << "\n";
    ccheck(res.status == pb::Status::Ok, "ns2d_16 two-level composite solve Ok");
    ccheck(res.residual_rel2 <= 1e-8, "ns2d_16 two-level composite true residual <= 1e-8");
    ccheck(res.ncells_uncovered == 16*16 - 8*8 + 16*16, "uncovered cell count (192 coarse + 256 fine)");

    // y-direction: the hierarchy with a one-cell y must give the same answer as the same hierarchy with the
    // direction made trivial in the reference operator: check the divergence of the composite face gradient
    // (the y gradient is zero, the y term of the operator vanishes under periodic and Neumann).
    std::vector<std::array<MultiFab*,3>> gp(2);
    std::array<std::unique_ptr<MultiFab>,3> g0s, g1s;
    for (int d = 0; d < 3; ++d) {
        g0s[d] = std::make_unique<MultiFab>(amrex::convert(ba0, IntVect::TheDimensionVector(d)), dm0, 1, 0);
        g1s[d] = std::make_unique<MultiFab>(amrex::convert(ba1, IntVect::TheDimensionVector(d)), dm1, 1, 0);
        gp[0][d] = g0s[d].get(); gp[1][d] = g1s[d].get();
    }
    pb::PressureResult gr = pb::face_gradient_composite(p, gp);
    ccheck(gr.status == pb::Status::Ok, "ns2d_16 face gradient Ok");
    ccheck(gp[0][1]->norm0(0) == 0.0 && gp[1][1]->norm0(0) == 0.0, "ns2d_16 y gradient is exactly zero");
    // divergence of the (x,z) gradients vs the mean-removed RHS on uncovered cells
    double bsum = 0.0, vsum = 0.0;
    std::array<MultiFab*,2> rh{&rhs0, &rhs1};
    std::array<Geometry,2> gg{g0, g1};
    std::array<BoxArray,2> bb{ba0, ba1};
    iMultiFab m0(ba0, dm0, 1, 0); m0.setVal(1);
    {
        BoxArray cba = ba1; cba.coarsen(r);
        for (MFIter mfi(m0); mfi.isValid(); ++mfi) {
            auto const& a = m0.array(mfi);
            for (auto const& is : cba.intersections(mfi.validbox())) { amrex::LoopOnCpu(is.second, [&] (int i, int j, int k) { a(i,j,k) = 0; }); }
        }
    }
    for (int l = 0; l < 2; ++l) {
        auto const dxa = gg[l].CellSizeArray(); const double v = dxa[0]*dxa[1]*dxa[2];
        for (MFIter mfi(*rh[l]); mfi.isValid(); ++mfi) {
            auto const& a = rh[l]->const_array(mfi);
            Array4<int const> u; if (l == 0) { u = m0.const_array(mfi); }
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { if (l == 1 || u(i,j,k)) { bsum += v*a(i,j,k); vsum += v; } });
        }
    }
    ParallelDescriptor::ReduceRealSum(bsum); ParallelDescriptor::ReduceRealSum(vsum);
    const double shift = bsum/vsum;
    double r2 = 0.0, b2 = 0.0;
    for (int l = 0; l < 2; ++l) {
        auto const dxa = gg[l].CellSizeArray(); const double v = dxa[0]*dxa[1]*dxa[2];
        for (MFIter mfi(*rh[l]); mfi.isValid(); ++mfi) {
            auto const& a = rh[l]->const_array(mfi);
            auto const& gx = gp[l][0]->const_array(mfi); auto const& gz = gp[l][2]->const_array(mfi);
            Array4<int const> u; if (l == 0) { u = m0.const_array(mfi); }
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                if (l == 0 && !u(i,j,k)) { return; }
                const double div = (gx(i+1,j,k) - gx(i,j,k))/dxa[0] + (gz(i,j,k+1) - gz(i,j,k))/dxa[2];
                const double b_ = a(i,j,k) - shift;
                r2 += v*(b_-div)*(b_-div); b2 += v*b_*b_;
            });
        }
    }
    ParallelDescriptor::ReduceRealSum(r2); ParallelDescriptor::ReduceRealSum(b2);
    Print() << std::setprecision(6) << "NS2D div_gradient_vs_rhs rel2=" << std::sqrt(r2/b2) << "\n";
    ccheck(std::sqrt(r2/b2) <= 1e-8, "ns2d_16 divergence of the composite face gradient equals the mean-removed RHS (x,z only)");
    // Equivalence of the 2-D hierarchy to a 3-D one with several y cells and the same data in y (Neumann/periodic y):
    {
        const int ny = 4;
        RealBox rb3({0.0, 0.0, 0.0}, {twopi, 4*twopi/16, twopi});   // isotropic cells: dy = dx (no anisotropy for the smoother)
        Box e0(IntVect(0,0,0), IntVect(15,ny-1,15));
        Box e1 = amrex::refine(e0, IntVect(2,2,2));
        Geometry h0(e0, rb3, CoordSys::cartesian, isp), h1(e1, rb3, CoordSys::cartesian, isp);
        BoxArray a0(e0); a0.maxSize(mgs);
        BoxArray a1(amrex::refine(Box(IntVect(4,0,4), IntVect(11,ny-1,11)), IntVect(2,2,2))); a1.maxSize(mgs);
        DistributionMapping m0_(a0), m1_(a1);
        MultiFab s0(a0, m0_, 1, 0), s1(a1, m1_, 1, 0), q0(a0, m0_, 1, 1), q1(a1, m1_, 1, 1);
        q0.setVal(0.0); q1.setVal(0.0);
        for (int l = 0; l < 2; ++l) {
            MultiFab& m = (l == 0) ? s0 : s1;
            auto const dxa = (l == 0 ? h0 : h1).CellSizeArray();
            for (MFIter mfi(m); mfi.isValid(); ++mfi) {
                auto const& a = m.array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { a(i,j,k) = src((i+0.5)*dxa[0], (k+0.5)*dxa[2]); });
            }
        }
        pb::PressureProblem p3;
        p3.bc.fill(pb::BC::Periodic);
        pb::PressureLevel A0; A0.ba = a0; A0.dm = m0_; A0.geom = h0; A0.rhs = &s0; A0.phi = &q0;
        pb::PressureLevel A1; A1.ba = a1; A1.dm = m1_; A1.geom = h1; A1.ref_ratio = IntVect(2,2,2); A1.rhs = &s1; A1.phi = &q1;
        p3.levels = {A0, A1};
        pb::PressureResult r3 = pb::solve_pressure(p3, o);
        Print() << "NS2D variant3d status=" << pb::to_string(r3.status) << " iters=" << r3.backend_status.iterations << " true_rel2=" << r3.residual_rel2 << " msg=" << r3.message << "\n";
        ccheck(r3.status == pb::Status::Ok, "3-D four-cell-y variant with y-independent data solves Ok");
        // compare the j = 0 plane of the 3-D solution on level 0 uncovered cells and level 1 with the 2-D one; the grids differ in
        // y-resolution only, so x-z truncation errors are identical per level: compare coarse level only on uncovered cells
        double num = 0.0, den = 0.0;
        for (int l = 0; l < 2; ++l) {
            MultiFab& a2 = (l == 0) ? phi0 : phi1;
            MultiFab& a3 = (l == 0) ? q0 : q1;
            BoxArray b2d = (l == 0) ? ba0 : ba1;
            // project the 3-D solution (j = 0 plane) to the 2-D layout
            MultiFab t(b2d, (l == 0) ? dm0 : dm1, 1, 0);
            t.setVal(0.0);
            Box plane = (l == 0) ? Box(IntVect(0,0,0), IntVect(15,0,15)) : Box(IntVect(0,0,0), IntVect(31,0,31));
            auto v3 = gather_box(a3, (l == 0) ? e0 : e1);
            auto v2 = gather_box(a2, plane);
            if (ParallelDescriptor::IOProcessor()) {
                const Box big = (l == 0) ? e0 : e1;
                for (int k = plane.smallEnd(2); k <= plane.bigEnd(2); ++k) for (int i = plane.smallEnd(0); i <= plane.bigEnd(0); ++i) {
                    const Long i3 = i + Long(big.length(0))*(0 + Long(big.length(1))*k);
                    const Long i2 = i + Long(plane.length(0))*(0 + Long(plane.length(1))*k);
                    num += (v3[i3] - v2[i2])*(v3[i3] - v2[i2]); den += v2[i2]*v2[i2];
                }
            }
        }
        ParallelDescriptor::ReduceRealSum(num); ParallelDescriptor::ReduceRealSum(den);
        // The two hierarchies are the same discrete (x,z) problem when the y-cells carry identical data and y is periodic
        // (the y second difference of a y-independent field is 0). Level-1 y is refined to 8 cells in the 3-D variant.
        Print() << std::setprecision(6) << "NS2D y_equivalence rel_l2_vs_3d=" << std::sqrt(num/den) << "\n";
        ccheck(std::sqrt(num/den) <= 1e-8, "one-cell-y composite equals the same data on a 4-cell periodic y (rel L2 " + std::to_string(std::sqrt(num/den)) + ")");
    }
}

// --------------------------------------------------------------------------------------------------------------
// mode=comp_sel
// --------------------------------------------------------------------------------------------------------------
void run_comp_sel ()
{
    Cfg c; c.n = 16; c.nlev = 2; c.ratio = 2; c.mgs = 8; c.set_bc(pb::BC::Neumann);
    pb::PressureOptions o; o.verbose = 0;
    auto expect = [&] (std::string const& what, Cfg cfg, pb::Status st, std::function<void(pb::PressureProblem&)> mod, pb::BackendKind req, std::string const& frag) {
        Hier h = build_hier(cfg, 0);
        pb::PressureProblem p = make_problem(cfg, h);
        if (mod) { mod(p); }
        for (int l = 0; l < static_cast<int>(p.levels.size()); ++l) { p.levels[l].phi->setVal(7.0); }
        o.backend = req;
        pb::PressureResult r = pb::solve_pressure(p, o);
        bool ok = (r.status == st);
        if (st != pb::Status::Ok) {
            ok = ok && !r.message.empty() && (frag.empty() || r.message.find(frag) != std::string::npos);
            for (int l = 0; l < static_cast<int>(p.levels.size()); ++l) { ok = ok && p.levels[l].phi->min(0) == 7.0 && p.levels[l].phi->max(0) == 7.0; }
        }
        Print() << "  " << what << ": status=" << pb::to_string(r.status) << (r.message.empty() ? "" : " msg=\"" + r.message + "\"") << "\n";
        ccheck(ok, "composite selector: " + what);
    };
    using S = pb::Status; using K = pb::BackendKind;
    expect("two levels Ok (Auto -> MLMG)", c, S::Ok, nullptr, K::Auto, "");
    expect("two levels Ok (explicit MLMG)", c, S::Ok, nullptr, K::MLMG, "");
    expect("composite on the FFT backend not built", c, S::NotBuilt, nullptr, K::FFT, "FFT");
    expect("cylindrical not built", c, S::NotBuilt, [] (pb::PressureProblem& p) { p.cylindrical = true; }, K::Auto, "cylindrical");
    expect("mixed open/closed faces not built", c, S::NotBuilt, [] (pb::PressureProblem& p) { p.bc[pb::face_index(0,1)] = pb::BC::Dirichlet; }, K::Auto, "mixed");
    expect("non-uniform widths not built", c, S::NotBuilt, [] (pb::PressureProblem& p) { for (int d = 0; d < 3; ++d) { p.cell_width[d].assign(16, Real(1./16)); } p.cell_width[2][3] = Real(0.1); }, K::Auto, "non-uniform");
    {
        iMultiFab cls;
        expect("masked cells not built", c, S::NotBuilt, [&] (pb::PressureProblem& p) {
            cls.define(p.levels[0].ba, p.levels[0].dm, 1, 0); cls.setVal(1); p.cell_class = &cls; }, K::Auto, "masked");
    }
    {
        MultiFab a;
        expect("variable coefficient not built", c, S::NotBuilt, [&] (pb::PressureProblem& p) {
            a.define(p.levels[0].ba, p.levels[0].dm, 1, 0); a.setVal(1.0); p.cell_coef_a = &a; }, K::Auto, "coefficient");
    }
    {
        MultiFab a;
        expect("single-level gauge_weight with levels not built (use PressureLevel fields)", c, S::NotBuilt, [&] (pb::PressureProblem& p) {
            a.define(p.levels[0].ba, p.levels[0].dm, 1, 0); a.setVal(1.0); p.gauge_weight = &a; }, K::Auto, "PressureLevel");
    }
    {
        std::vector<std::unique_ptr<MultiFab>> g;
        auto field = [&g] (pb::PressureProblem& p, int l) { g.push_back(std::make_unique<MultiFab>(p.levels[l].ba, p.levels[l].dm, 1, 0)); g.back()->setVal(1.0); return g.back().get(); };
        expect("per-level gauge fields on every level Ok", c, S::Ok, [&] (pb::PressureProblem& p) {
            for (int l = 0; l < 2; ++l) { p.levels[l].gauge_weight = field(p, l); p.levels[l].gauge_offset = field(p, l); } }, K::Auto, "");
        expect("gauge_weight on one level only invalid", c, S::InvalidInput, [&] (pb::PressureProblem& p) {
            p.levels[0].gauge_weight = field(p, 0); }, K::Auto, "every level");
        expect("gauge_offset on the fine level only invalid", c, S::InvalidInput, [&] (pb::PressureProblem& p) {
            p.levels[1].gauge_offset = field(p, 1); }, K::Auto, "every level");
        expect("gauge field with another BoxArray invalid", c, S::InvalidInput, [&] (pb::PressureProblem& p) {
            for (int l = 0; l < 2; ++l) { p.levels[l].gauge_weight = field(p, l); }
            BoxArray b2(p.levels[0].ba.minimalBox()); b2.maxSize(4);
            g.push_back(std::make_unique<MultiFab>(b2, DistributionMapping(b2), 1, 0)); g.back()->setVal(1.0);
            p.levels[0].gauge_weight = g.back().get(); }, K::Auto, "differ from the level");
    }
    expect("ratio 3 not built", c, S::NotBuilt, [] (pb::PressureProblem& p) { p.levels[1].ref_ratio = IntVect(3); }, K::Auto, "ratio");
    expect("anisotropic ratio not built", c, S::NotBuilt, [] (pb::PressureProblem& p) { p.levels[1].ref_ratio = IntVect(2,2,4); }, K::Auto, "anisotropic");
    expect("dirichlet two levels Ok", [&] { Cfg d = c; d.set_bc(pb::BC::Dirichlet); return d; }(), S::Ok, nullptr, K::Auto, "");
    expect("periodic three levels Ok", [&] { Cfg d = c; d.set_bc(pb::BC::Periodic); d.nlev = 3; d.n = 32; return d; }(), S::Ok, nullptr, K::Auto, "");
    expect("one-level composite Ok", [&] { Cfg d = c; d.nlev = 1; return d; }(), S::Ok, nullptr, K::Auto, "");
    expect("plane2d Neumann Ok", [&] { Cfg d = c; d.plane2d = true; return d; }(), S::Ok, nullptr, K::Auto, "");
    // InvalidInput
    std::vector<std::unique_ptr<MultiFab>> keep;
    auto set_ba = [&keep] (pb::PressureProblem& p, int l, BoxArray const& b) {
        p.levels[l].ba = b; p.levels[l].dm = DistributionMapping(b);
        keep.push_back(std::make_unique<MultiFab>(b, p.levels[l].dm, 1, 0)); keep.back()->setVal(0.0); p.levels[l].rhs = keep.back().get();
        keep.push_back(std::make_unique<MultiFab>(b, p.levels[l].dm, 1, 1)); p.levels[l].phi = keep.back().get();
    };
    expect("level 0 not covering the domain invalid", c, S::InvalidInput, [&] (pb::PressureProblem& p) {
        set_ba(p, 0, BoxArray(Box(IntVect(0), IntVect(7)))); }, K::Auto, "cover the domain exactly");
    expect("fine boxes not aligned to the ratio invalid", c, S::InvalidInput, [&] (pb::PressureProblem& p) {
        set_ba(p, 1, BoxArray(Box(IntVect(5,5,5), IntVect(12,12,12)))); }, K::Auto, "aligned");
    expect("fine level outside the coarse BoxArray invalid", [&] { Cfg d = c; d.nlev = 3; d.n = 32; return d; }(), S::InvalidInput, [&] (pb::PressureProblem& p) {
        // level 2 box (level-2 indices) entirely outside the coarsened range of level 1
        set_ba(p, 2, BoxArray(Box(IntVect(0,0,0), IntVect(7,7,7)))); }, K::Auto, "nested");
    expect("no buffer cell around level 2 invalid (3 levels)", [&] { Cfg d = c; d.nlev = 3; d.n = 32; return d; }(), S::InvalidInput, [&] (pb::PressureProblem& p) {
        // level 1 covers [16,47]; a level 2 box flush with its lower edge has no coarse neighbour cell
        set_ba(p, 2, BoxArray(Box(IntVect(32,32,32), IntVect(63,63,63)))); }, K::Auto, "nested");
    expect("domain of level 1 not the refined domain invalid", c, S::InvalidInput, [] (pb::PressureProblem& p) {
        p.levels[1].ref_ratio = IntVect(4); }, K::Auto, "refined by ref_ratio");
    expect("periodicity not matching the BC invalid", c, S::InvalidInput, [] (pb::PressureProblem& p) {
        p.bc.fill(pb::BC::Periodic); }, K::Auto, "periodicity");
    {
        MultiFab other;
        expect("rhs layout mismatch invalid", c, S::InvalidInput, [&] (pb::PressureProblem& p) {
            BoxArray b(p.levels[1].ba.minimalBox()); b.maxSize(4); other.define(b, DistributionMapping(b), 1, 0); p.levels[1].rhs = &other; }, K::Auto, "rhs/phi BoxArray");
    }
    {   // a fine level flush with a non-periodic domain boundary is fine (no buffer needed there); periodic wrap is fine too
        Cfg d = c; d.full = true;
        expect("fine level over the whole domain Ok (non-periodic)", d, S::Ok, nullptr, K::Auto, "");
        d.set_bc(pb::BC::Periodic);
        expect("fine level over the whole domain Ok (periodic)", d, S::Ok, nullptr, K::Auto, "");
    }
    // legacy: nlevels > 1 on the single-level fields stays NotBuilt
    {
        Cfg d = c; d.nlev = 1;
        Hier h = build_hier(d, 0);
        pb::PressureProblem p = make_problem(d, h);
        pb::PressureProblem q;
        q.ba = p.levels[0].ba; q.dm = p.levels[0].dm; q.geom = p.levels[0].geom; q.rhs = p.levels[0].rhs; q.phi = p.levels[0].phi; q.bc = p.bc; q.nlevels = 2;
        pb::PressureResult r = pb::solve_pressure(q, o);
        ccheck(r.status == pb::Status::NotBuilt, "single-level fields with nlevels = 2 still NotBuilt");
    }
}

// --------------------------------------------------------------------------------------------------------------
// mode=comp_ws: D-058 workspace rebuild
// --------------------------------------------------------------------------------------------------------------
void run_comp_ws ()
{
    Cfg c; c.n = 16; c.nlev = 2; c.ratio = 2; c.mgs = 8; c.set_bc(pb::BC::Periodic);
    pb::PressureOptions o; o.verbose = 0; o.removed_mean_warn = 1.0; o.residual_tol = 1e-8;
    pb::PressureWorkspace ws;
    ccheck(!ws.built() && ws.num_levels() == 0, "empty workspace is not built");
    {
        Hier h = build_hier(c, 0);
        pb::PressureProblem p = make_problem(c, h);
        std::string msg;
        ccheck(ws.rebuild(p, &msg) == pb::Status::Ok && ws.built() && ws.num_levels() == 2 && ws.matches(p), "rebuild(hierarchy A)");
        pb::PressureResult r = pb::solve_pressure(p, o, &ws);
        ccheck(r.status == pb::Status::Ok && !r.workspace_rebuilt && r.residual_rel2 <= 1e-8, "solve with the workspace, hierarchy A");
    }   // hierarchy A (rhs, phi, BoxArrays) destroyed: the workspace must not touch them again
    {
        Cfg c2 = c; c2.mgs = 4;     // different decomposition: stale workspace
        Hier h2 = build_hier(c2, 0);
        pb::PressureProblem p2 = make_problem(c2, h2);
        ccheck(!ws.matches(p2), "workspace of A does not match the new decomposition");
        pb::PressureResult r = pb::solve_pressure(p2, o, &ws);
        ccheck(r.status == pb::Status::Ok && r.workspace_rebuilt && r.residual_rel2 <= 1e-8, "stale workspace is rebuilt automatically and flagged");
        ccheck(ws.matches(p2), "workspace matches after the automatic rebuild");
        // explicit regrid: a different patch
        Cfg c3 = c; c3.full = true;
        Hier h3 = build_hier(c3, 0);
        pb::PressureProblem p3 = make_problem(c3, h3);
        ccheck(ws.rebuild(p3) == pb::Status::Ok && ws.matches(p3), "rebuild(hierarchy C: fine level over the whole domain)");
        pb::PressureResult r3 = pb::solve_pressure(p3, o, &ws);
        ccheck(r3.status == pb::Status::Ok && !r3.workspace_rebuilt && r3.residual_rel2 <= 1e-8, "solve after explicit rebuild");
        // rebuild with an unsupported problem leaves the workspace empty
        pb::PressureProblem bad = p3; bad.cylindrical = true;
        std::string msg;
        ccheck(ws.rebuild(bad, &msg) == pb::Status::NotBuilt && !ws.built() && !msg.empty(), "failed rebuild empties the workspace and reports the reason");
    }
}


// --------------------------------------------------------------------------------------------------------------
// mode=comp_trigger: FR-039 trigger points on a hierarchy (full checks versus the cheap path)
// --------------------------------------------------------------------------------------------------------------
void run_comp_trigger ()
{
    for (pb::BC bc : {pb::BC::Neumann, pb::BC::Periodic}) {
        Cfg c; c.n = 16; c.nlev = 2; c.ratio = 2; c.mgs = 8; c.set_bc(bc);
        const std::string tag = (bc == pb::BC::Neumann ? "neumann: " : "periodic: ");
        Hier h = build_hier(c, 0);
        pb::PressureProblem p = make_problem(c, h);
        auto gather_all = [&] () {
            std::vector<double> all;
            for (int l = 0; l < c.nlev; ++l) {
                std::vector<double> v = gather_box(*h.phi[l], h.geom[l].Domain());
                all.insert(all.end(), v.begin(), v.end());
            }
            return all;
        };
        auto run = [&] (pb::PressureOptions const& oo, pb::PressureWorkspace* ws, std::vector<double>& out) {
            for (int l = 0; l < c.nlev; ++l) { h.phi[l]->setVal(0.0); }
            pb::PressureResult r = pb::solve_pressure(p, oo, ws);
            out = gather_all();
            return r;
        };
        pb::PressureOptions o; o.verbose = 0; o.removed_mean_warn = 1.0; o.residual_tol = 1e-8;
        std::vector<double> ref, v;
        pb::PressureResult r0 = run(o, nullptr, ref);
        ccheck(r0.status == pb::Status::Ok && r0.full_checks && r0.residual_checked && r0.triggers == pb::SolveDebug && r0.residual_rel2 <= 1e-8, tag + "default options: full checks on a hierarchy");
        pb::PressureOptions oc = o; oc.trigger = pb::SolveRoutine;
        pb::PressureResult r1 = run(oc, nullptr, v);
        ccheck(r1.status == pb::Status::Ok && !r1.full_checks && !r1.residual_checked && v == ref, tag + "Routine, no workspace: cheap path, solution bitwise equal");
        ccheck(r1.components.size() == r0.components.size() && r1.components[0].removed_mean == r0.components[0].removed_mean && r1.components[0].gauge_shift == r0.components[0].gauge_shift, tag + "cheap path: same removed mean and gauge constant (bitwise)");
        pb::PressureWorkspace ws;
        pb::PressureResult a = run(oc, &ws, v);
        ccheck(a.full_checks && (a.triggers & pb::SolveFirst) && a.residual_checked && v == ref, tag + "workspace, first solve: full checks");
        pb::PressureResult b = run(oc, &ws, v);
        ccheck(!b.full_checks && !b.residual_checked && v == ref, tag + "workspace, second solve: cheap path, bitwise equal");
        std::string msg;
        ccheck(ws.rebuild(p, &msg) == pb::Status::Ok && ws.regrid_pending(), tag + "rebuild() after a regrid marks the next solve (" + msg + ")");
        pb::PressureResult c1 = run(oc, &ws, v);
        ccheck(c1.full_checks && (c1.triggers & pb::SolveFirstAfterRegrid) && c1.residual_checked && v == ref, tag + "first solve after rebuild: full checks");
        ccheck(!run(oc, &ws, v).full_checks, tag + "following solve: cheap again");
        pb::PressureOptions od = oc; od.trigger = pb::SolveDebug;
        ccheck(run(od, &ws, v).full_checks, tag + "Debug flag: full checks");
        // a changed hierarchy seen by solve_pressure counts as a regrid
        {
            Cfg c2 = c; c2.mgs = 4;
            Hier h2 = build_hier(c2, 0);
            pb::PressureProblem p2 = make_problem(c2, h2);
            pb::PressureResult d = pb::solve_pressure(p2, oc, &ws);
            ccheck(d.status == pb::Status::Ok && d.workspace_rebuilt && d.full_checks && (d.triggers & pb::SolveFirstAfterRegrid), tag + "new decomposition in a solve: detected as regrid, full checks");
        }
        // the cheap path still reports non-convergence and unsupported input
        {
            pb::PressureProblem bad = p; bad.cylindrical = true;
            pb::PressureResult e = pb::solve_pressure(bad, oc, &ws);
            ccheck(e.status == pb::Status::NotBuilt && !e.message.empty(), tag + "cheap path: NotBuilt still reported");
            pb::PressureOptions ot = oc; ot.max_iter = 1; ot.tol_rel = 1e-14;
            pb::PressureResult nc = run(ot, nullptr, v);
            ccheck(nc.status == pb::Status::NotConverged, tag + "cheap path: NotConverged still reported");
        }
    }
}

// --------------------------------------------------------------------------------------------------------------
// mode=comp_gauge: D-067 mean removal (Volume / ScaledArithmetic) and gauge (rho, KRES) on a hierarchy.
//   Solves the same hierarchy with (A) Volume, no gauge fields, (B) Volume with per-level rho and KRES, (C) the FDS
//   parity switch ScaledArithmetic without gauge fields. RHS = manufactured RHS + rhs_offset (so that a mean is
//   removed). Writes level-concatenated full-domain dumps for the independent numpy check in tests/pb_test.py.
// --------------------------------------------------------------------------------------------------------------
Real gauge_rho_fn (Real x, Real y, Real z)
{
    return Real(1.2)*(Real(1.0) + Real(0.3)*std::sin(Real(2.)*kPi*x)*std::cos(kPi*y) + Real(0.2)*z);
}
Real gauge_kres_fn (Real x, Real y, Real z) { return Real(0.05)*(x*x + Real(2.)*y - z); }
Real probe_phi_fn (Real x, Real y, Real z) { return std::cos(Real(3.)*x + y) + z*z + Real(0.7) + Real(1.0e-3)*std::sin(Real(40.)*x*y); }

std::string hexd_c (double x) { std::uint64_t u; std::memcpy(&u, &x, 8); return hx(u); }

void run_comp_gauge (ParmParse& pp)
{
    Cfg c;
    pp.query("n", c.n); pp.query("nlev", c.nlev); pp.query("ratio", c.ratio);
    int plane = 0; pp.query("plane2d", plane); c.plane2d = plane != 0;
    std::string bcs = "neumann"; pp.query("bc", bcs); c.set_bc(cparse_bc(bcs));
    pp.query("mgs", c.mgs); pp.query("layout", c.layout);
    c.kx = 2; c.ky = c.plane2d ? 0 : 1; c.kz = 2;
    int dmkind = 0; pp.query("dmkind", dmkind);
    double offset = 0.37; pp.query("rhs_offset", offset);
    double tol_rel = 1e-13; pp.query("tol_rel", tol_rel);
    std::string out; pp.query("out", out);

    Hier h = build_hier(c, dmkind);
    for (auto& r : h.rhs) { r->plus(Real(offset), 0, 1, 0); }
    // per-level gauge fields and an analytic probe field
    std::vector<std::unique_ptr<MultiFab>> rho, kres, probe;
    for (int l = 0; l < c.nlev; ++l) {
        rho.push_back(std::make_unique<MultiFab>(h.ba[l], h.dm[l], 1, 0));
        kres.push_back(std::make_unique<MultiFab>(h.ba[l], h.dm[l], 1, 0));
        probe.push_back(std::make_unique<MultiFab>(h.ba[l], h.dm[l], 1, 0));
        auto const dxa = h.geom[l].CellSizeArray(); auto const lo = h.geom[l].ProbLoArray();
        for (MFIter mfi(*rho[l]); mfi.isValid(); ++mfi) {
            auto const& r = rho[l]->array(mfi); auto const& k_ = kres[l]->array(mfi); auto const& q = probe[l]->array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const Real x = lo[0] + (i+0.5)*dxa[0], y = lo[1] + (j+0.5)*dxa[1], z = lo[2] + (k+0.5)*dxa[2];
                r(i,j,k) = gauge_rho_fn(x, y, z); k_(i,j,k) = gauge_kres_fn(x, y, z); q(i,j,k) = probe_phi_fn(x, y, z);
            });
        }
    }
    pb::PressureOptions o; o.tol_rel = tol_rel; o.removed_mean_warn = 1.0e3; o.residual_tol = 1e-8; o.verbose = 1;
    auto solve = [&] (pb::MeanKind mk, bool gauge, std::vector<double>& dump) {
        pb::PressureProblem p = make_problem(c, h);
        p.mean_kind = mk;
        if (gauge) { for (int l = 0; l < c.nlev; ++l) { p.levels[l].gauge_weight = rho[l].get(); p.levels[l].gauge_offset = kres[l].get(); } }
        for (auto& f : h.phi) { f->setVal(0.0); }
        pb::PressureResult r = pb::solve_pressure(p, o);
        dump = hier_gather(h, c.nlev);
        return r;
    };
    std::vector<double> phiA, phiB, phiC;
    pb::PressureResult rA = solve(pb::MeanKind::Volume, false, phiA);
    pb::PressureResult rB = solve(pb::MeanKind::Volume, true, phiB);
    pb::PressureResult rC = solve(pb::MeanKind::ScaledArithmetic, false, phiC);
    ccheck(rA.status == pb::Status::Ok && rB.status == pb::Status::Ok && rC.status == pb::Status::Ok, "all three composite solves Ok");
    ccheck(rA.residual_ok && rB.residual_ok && rC.residual_ok, "true composite residual within eps_H in all three");
    auto const& cA = rA.components[0]; auto const& cB = rB.components[0]; auto const& cC = rC.components[0];
    ccheck(cA.removed_mean == cB.removed_mean, "gauge fields do not change the removed mean (bitwise)");
    if (cA.singular) {
        ccheck(std::abs(cA.gauge_shift) > 0.0 && cA.gauge_shift != cB.gauge_shift, "weighted gauge constant differs from the plain one (test is sensitive)");
    } else {
        // A component with an open face has a unique solution: no mean removal, no gauge, gauge fields are ignored.
        ccheck(cA.removed_mean == 0.0 && cC.removed_mean == 0.0, "non-singular component: nothing removed from the RHS");
        ccheck(cA.gauge_shift == 0.0 && cB.gauge_shift == 0.0 && cC.gauge_shift == 0.0, "non-singular component: no gauge shift");
        ccheck(phiA == phiB, "non-singular component: gauge fields leave the solution bitwise unchanged");
    }
    Print() << std::setprecision(17) << "GAUGE bc=" << bcs << " nranks=" << ParallelDescriptor::NProcs() << " mgs=" << c.mgs << " dmkind=" << dmkind
            << " nlev=" << c.nlev << " ratio=" << c.ratio << " plane2d=" << plane << " offset=" << offset
            << " removedV=" << hexd_c(cA.removed_mean) << " removedS=" << hexd_c(cC.removed_mean)
            << " removedV_val=" << cA.removed_mean << " removedS_val=" << cC.removed_mean
            << " removed_relV=" << cA.removed_rel << " removed_relS=" << cC.removed_rel
            << " shiftA=" << cA.gauge_shift << " shiftB=" << cB.gauge_shift << " shiftC=" << cC.gauge_shift
            << " nunc=" << rA.ncells_uncovered << " iters=" << rA.backend_status.iterations;
    for (int l = 0; l < c.nlev; ++l) {
        auto const dxa = h.geom[l].CellSizeArray();
        Print() << std::setprecision(17) << " vol" << l << "=" << double(dxa[0])*double(dxa[1])*double(dxa[2]);
    }
    Print() << "\n";

    // Building block, decomposition independent by construction: the gauge numerator and denominator of a fixed analytic
    // field with the same exact sums the solver uses. Printed as hex for bitwise comparison across rank counts and layouts.
    {
        std::vector<pb::SumTerm> num, den;
        std::vector<std::unique_ptr<MultiFab>> X;
        for (int l = 0; l < c.nlev; ++l) {
            X.push_back(std::make_unique<MultiFab>(h.ba[l], h.dm[l], 1, 0));
            MultiFab::Copy(*X[l], *probe[l], 0, 0, 1, 0);
            MultiFab::Subtract(*X[l], *kres[l], 0, 0, 1, 0);
            auto const dxa = h.geom[l].CellSizeArray();
            const double v = double(dxa[0])*double(dxa[1])*double(dxa[2]);
            pb::SumTerm tn; tn.mf = X[l].get(); tn.weight = v; tn.uncovered = h.unc[l].get(); tn.wfield = rho[l].get(); num.push_back(tn);
            pb::SumTerm td; td.mf = rho[l].get(); td.weight = v; td.uncovered = h.unc[l].get(); den.push_back(td);
        }
        const double sn = pb::exact_sum_multi(num, 1).sum[0], sd = pb::exact_sum_multi(den, 1).sum[0];
        Print() << std::setprecision(17) << "GAUGE_EXACT num=" << hexd_c(sn) << " den=" << hexd_c(sd) << " shift=" << hexd_c(sn/sd) << " shift_val=" << sn/sd << "\n";
    }

    // Independent divergence check of the ScaledArithmetic solution: div(grad phi_C) must equal b - mean(F)/v_l on the
    // uncovered cells, with mean(F) computed here from the dumped RHS (not from the library). The Volume zero mode is the
    // negative control: it must NOT match (volumes differ between levels).
    if (c.nlev > 1 && cA.singular) {
        for (int l = 0; l < c.nlev; ++l) { h.phi[l]->setVal(0.0); }
        pb::PressureProblem p = make_problem(c, h);
        p.mean_kind = pb::MeanKind::ScaledArithmetic;
        pb::PressureResult r = pb::solve_pressure(p, o);
        ccheck(r.status == pb::Status::Ok, "ScaledArithmetic solve for the divergence check Ok");
        GradOut g = make_grad(c, h);
        pb::PressureResult gr = pb::face_gradient_composite(p, g.ptr);
        ccheck(gr.status == pb::Status::Ok, "face_gradient_composite status Ok");
        double sF = 0.0, sV = 0.0, sVb = 0.0; double cnt = 0.0;
        for (int l = 0; l < c.nlev; ++l) {
            auto const dxa = h.geom[l].CellSizeArray(); const double v = double(dxa[0])*double(dxa[1])*double(dxa[2]);
            for (MFIter mfi(*h.rhs[l]); mfi.isValid(); ++mfi) {
                auto const& b = h.rhs[l]->const_array(mfi); auto const& u = h.unc[l]->const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { if (u(i,j,k)) { sF += v*b(i,j,k); sV += v; cnt += 1.0; (void)sVb; } });
            }
        }
        ParallelDescriptor::ReduceRealSum(sF); ParallelDescriptor::ReduceRealSum(sV); ParallelDescriptor::ReduceRealSum(cnt);
        const double meanF = sF/cnt, meanV = sF/sV;
        const int hd = c.plane2d ? 1 : -1;
        auto div_err = [&] (bool scaled) {
            double r2 = 0.0, b2 = 0.0;
            for (int l = 0; l < c.nlev; ++l) {
                auto const dxa = h.geom[l].CellSizeArray(); const double v = double(dxa[0])*double(dxa[1])*double(dxa[2]);
                const double shift = scaled ? meanF/v : meanV;
                for (MFIter mfi(*h.rhs[l]); mfi.isValid(); ++mfi) {
                    auto const& b = h.rhs[l]->const_array(mfi); auto const& u = h.unc[l]->const_array(mfi);
                    auto const& gx = g.ptr[l][0]->const_array(mfi); auto const& gy = g.ptr[l][1]->const_array(mfi); auto const& gz = g.ptr[l][2]->const_array(mfi);
                    amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                        if (!u(i,j,k)) { return; }
                        double div = (gx(i+1,j,k) - gx(i,j,k))/dxa[0] + (gz(i,j,k+1) - gz(i,j,k))/dxa[2];
                        if (hd != 1) { div += (gy(i,j+1,k) - gy(i,j,k))/dxa[1]; }
                        const double bb = b(i,j,k) - shift;
                        r2 += v*(bb-div)*(bb-div); b2 += v*bb*bb;
                    });
                }
            }
            ParallelDescriptor::ReduceRealSum(r2); ParallelDescriptor::ReduceRealSum(b2);
            return std::sqrt(r2/b2);
        };
        const double es = div_err(true), ev = div_err(false);
        Print() << std::setprecision(6) << "GAUGE_DIV scaled_zero_mode=" << es << " volume_zero_mode_control=" << ev << " meanF=" << meanF << " meanV=" << meanV << "\n";
        ccheck(es <= 1e-8, "ScaledArithmetic solution: div(grad phi) equals b - mean(v*b)/v on uncovered cells (independent formula), rel2 " + std::to_string(es));
        ccheck(ev > 100*std::max(es, 1e-12), "negative control: the Volume zero mode does not match the ScaledArithmetic solution (test is sensitive)");
    }

    if (!out.empty()) {
        auto w = [&] (std::string const& name, std::vector<double> const& v) {
            if (ParallelDescriptor::IOProcessor()) {
                std::ofstream ofs(out + name, std::ios::binary);
                ofs.write(reinterpret_cast<const char*>(v.data()), static_cast<std::streamsize>(v.size()*sizeof(double)));
            }
        };
        auto gatherl = [&] (std::vector<std::unique_ptr<MultiFab>> const& f) {
            std::vector<double> all;
            for (int l = 0; l < c.nlev; ++l) { auto v = gather_box(*f[l], h.geom[l].Domain()); all.insert(all.end(), v.begin(), v.end()); }
            return all;
        };
        std::vector<std::unique_ptr<MultiFab>> uncd;
        for (int l = 0; l < c.nlev; ++l) {
            uncd.push_back(std::make_unique<MultiFab>(h.ba[l], h.dm[l], 1, 0));
            for (MFIter mfi(*uncd[l]); mfi.isValid(); ++mfi) {
                auto const& a = uncd[l]->array(mfi); auto const& u = h.unc[l]->const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { a(i,j,k) = double(u(i,j,k)); });
            }
        }
        w("_phiA.bin", phiA); w("_phiB.bin", phiB); w("_phiC.bin", phiC);
        w("_rho.bin", gatherl(rho)); w("_kres.bin", gatherl(kres)); w("_unc.bin", gatherl(uncd));
        std::vector<std::unique_ptr<MultiFab>> rr;
        for (int l = 0; l < c.nlev; ++l) { rr.push_back(std::make_unique<MultiFab>(h.ba[l], h.dm[l], 1, 0)); MultiFab::Copy(*rr[l], *h.rhs[l], 0, 0, 1, 0); }
        w("_rhs.bin", gatherl(rr));
    }
}

} // anonymous namespace

int run_composite_mode (std::string const& mode, ParmParse& pp)
{
    if (mode == "comp") { run_comp(pp); }
    else if (mode == "comp_ns2d") { run_comp_ns2d(pp); }
    else if (mode == "comp_sel") { run_comp_sel(); }
    else if (mode == "comp_ws") { run_comp_ws(); }
    else if (mode == "comp_gauge") { run_comp_gauge(pp); }
    else if (mode == "comp_trigger") { run_comp_trigger(); }
    else { return -1; }
    return g_cfail;
}
