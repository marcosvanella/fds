// test_facetransfer.cpp (R4, D-060): the single face-velocity prolongation path (prolong_faces_level = prolong_faces_normal_linear on whole levels).
//  * every fine cell keeps the divergence of its coarse parent to round-off (random coarse field with divergence, and a discretely divergence-free one),
//    for ratio 2 in 3-D, ratio (4,2,4) mixed, and the 2-D case with a single cell in y (ratio 1 in y);
//  * fine faces that lie on a coarse face (the interface faces of a patch included) carry exactly the coarse value;
//  * max |div u - D| over the new cells is reported (D = 0 here; and D = the parent divergence, to show the reporting function);
//  * old fine data is copied bitwise over the overlap; the fine -> coarse face mean returns the coarse value;
//  * negative control: plain injection of the coarse face value breaks the parent divergence.
// Runs on 1 and 4 ranks (box-to-rank independent: values are a function of the face index only).
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>

#include "FaceTransfer.H"
#include "check.H"

using namespace fdsrt;

namespace {

double hrand(int i, int j, int k, int d, int salt)
{
    uint64_t h = 1469598103934665603ULL;
    for (long v : {long(i) + 1000, long(j) + 1000, long(k) + 1000, long(d), long(salt)}) { h ^= static_cast<uint64_t>(v); h *= 1099511628211ULL; h ^= h >> 29; h *= 0x9e3779b97f4a7c15ULL; }
    return static_cast<double>(h >> 11) / 9007199254740992.0 - 0.5;
}

struct Setup {
    amrex::Box cdom;
    amrex::IntVect ratio;
    amrex::Geometry cg, fg;
    amrex::BoxArray cba, fba;                 // cell-centred
    amrex::DistributionMapping cdm, fdm;
    std::array<amrex::MultiFab, 3> C, F;      // coarse / fine faces
};

void make(Setup& s, amrex::IntVect cn, amrex::IntVect ratio, amrex::Box fregion_c, int split_dir)
{
    s.cdom = amrex::Box(amrex::IntVect(0), cn - 1);
    s.ratio = ratio;
    amrex::RealBox rb({0, 0, 0}, {1, 1, 1});
    amrex::Array<int, 3> per{1, 1, 1};
    s.cg.define(s.cdom, rb, 0, per);
    s.fg = amrex::refine(s.cg, ratio);
    s.cba.define(s.cdom);
    s.cba.maxSize(8);
    s.cdm.define(s.cba);
    // fine patch = refined fregion_c, cut into two boxes along split_dir at a coarse-cell boundary
    amrex::Box fr = amrex::refine(fregion_c, ratio);
    amrex::Box a = fr, b = fr;
    const int cut = fregion_c.smallEnd(split_dir) + fregion_c.length(split_dir) / 2;   // first coarse cell of the second box
    a.setBig(split_dir, cut * ratio[split_dir] - 1);
    b.setSmall(split_dir, cut * ratio[split_dir]);
    amrex::BoxList bl; bl.push_back(a); bl.push_back(b);
    s.fba.define(bl);
    s.fdm.define(s.fba);
    for (int d = 0; d < 3; ++d) {
        s.C[d].define(amrex::convert(s.cba, amrex::IntVect::TheDimensionVector(d)), s.cdm, 1, 1);
        s.F[d].define(amrex::convert(s.fba, amrex::IntVect::TheDimensionVector(d)), s.fdm, 1, 1);
        s.F[d].setVal(-999.0);
    }
}

// coarse face values: random with divergence (kind 0), or the discrete curl of a random edge potential (kind 1; exactly div-free per cell, also in the 2-D case)
void fill_coarse(Setup& s, int kind)
{
    const amrex::IntVect n = s.cdom.length();
    auto wrap = [&](int v, int d) { return ((v % n[d]) + n[d]) % n[d]; };
    for (int d = 0; d < 3; ++d)
        for (amrex::MFIter mfi(s.C[d]); mfi.isValid(); ++mfi) {
            auto a = s.C[d].array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                if (kind == 0) { a(i, j, k) = hrand(wrap(i, 0), wrap(j, 1), wrap(k, 2), d, 1) * 2.0; return; }
                // potential A_e(i,j,k) on edges; face value u_d = d_e1 A_e2 - d_e2 A_e1 (centred on the face): use cell-edge potentials psi_e at integer indices
                auto A = [&](int e, int ii, int jj, int kk) { return hrand(wrap(ii, 0), wrap(jj, 1), wrap(kk, 2), e, 7); };
                // face of direction d at index (i,j,k): curl component d from the 4 edges around it: standard Yee/MAC discrete curl
                const int e1 = (d + 1) % 3, e2 = (d + 2) % 3;
                int p[3] = {i, j, k}, q1[3] = {i, j, k}, q2[3] = {i, j, k};
                // A_e2 at (p) and (p + e1), A_e1 at (p) and (p + e2): u = (A_e2(p+e1) - A_e2(p)) - (A_e1(p+e2) - A_e1(p))
                q1[e1] += 1; q2[e2] += 1;
                a(i, j, k) = (A(e2, q1[0], q1[1], q1[2]) - A(e2, p[0], p[1], p[2])) - (A(e1, q2[0], q2[1], q2[2]) - A(e1, p[0], p[1], p[2]));
            });
        }
    // the 2-D case (one cell in y): V must have the same value on both y faces for a zero divergence contribution; kind 1 only
    if (kind == 1 && n[1] == 1) {
        // with a single cell in y the wrap makes A(.., j+1) = A(.., j): the curl terms in y cancel by construction (periodic image), nothing to do
    }
}

struct Result { double max_div_diff = 0; double max_iface_diff = 0; double max_div_abs = 0; long nfaces_m0 = 0; long ncells = 0; };

Result evaluate(Setup& s)
{
    Result r;
    const amrex::IntVect R = s.ratio;
    // coarse faces under the fine boxes, gathered on the fine distribution (collective calls, so outside the box loop)
    amrex::BoxArray ccells = amrex::convert(s.fba, amrex::IntVect::TheCellVector());
    ccells.coarsen(R);
    std::array<amrex::MultiFab, 3> T;
    for (int d = 0; d < 3; ++d) {
        T[d].define(amrex::convert(ccells, amrex::IntVect::TheDimensionVector(d)), s.fdm, 1, 0);
        T[d].setVal(-1e30);
        T[d].ParallelCopy(s.C[d], 0, 0, 1, 0, 0, s.cg.periodicity());
    }
    for (amrex::MFIter mfi(s.F[0]); mfi.isValid(); ++mfi) {
        const amrex::Box cells = amrex::enclosedCells(mfi.validbox());
        auto u = s.F[0].const_array(mfi), v = s.F[1].const_array(mfi), w = s.F[2].const_array(mfi);
        auto cua = T[0].const_array(mfi), cva = T[1].const_array(mfi), cwa = T[2].const_array(mfi);
        const double rdxc[3] = {1.0 / s.cg.CellSize(0), 1.0 / s.cg.CellSize(1), 1.0 / s.cg.CellSize(2)};
        const double rdxf[3] = {1.0 / s.fg.CellSize(0), 1.0 / s.fg.CellSize(1), 1.0 / s.fg.CellSize(2)};
        amrex::LoopOnCpu(cells, [&](int i, int j, int k) {
            const amrex::IntVect c = amrex::coarsen(amrex::IntVect(i, j, k), R);
            const double dc = (cua(c[0] + 1, c[1], c[2]) - cua(c[0], c[1], c[2])) * rdxc[0] + (cva(c[0], c[1] + 1, c[2]) - cva(c[0], c[1], c[2])) * rdxc[1] +
                              (cwa(c[0], c[1], c[2] + 1) - cwa(c[0], c[1], c[2])) * rdxc[2];
            const double df = (u(i + 1, j, k) - u(i, j, k)) * rdxf[0] + (v(i, j + 1, k) - v(i, j, k)) * rdxf[1] + (w(i, j, k + 1) - w(i, j, k)) * rdxf[2];
            r.max_div_diff = std::max(r.max_div_diff, std::abs(df - dc));
            r.max_div_abs = std::max(r.max_div_abs, std::abs(df));
            ++r.ncells;
        });
        // faces on a coarse face: fine index multiple of the ratio in the normal direction; the value is the coarse face of the containing transverse cell
        for (int d = 0; d < 3; ++d) {
            const amrex::Box fb = amrex::convert(cells, amrex::IntVect::TheDimensionVector(d));
            auto f = (d == 0 ? u : d == 1 ? v : w);
            auto ca = (d == 0 ? cua : d == 1 ? cva : cwa);
            amrex::LoopOnCpu(fb, [&](int i, int j, int k) {
                const int fi[3] = {i, j, k};
                if (fi[d] % R[d] != 0) return;
                const amrex::IntVect c = amrex::coarsen(amrex::IntVect(i, j, k), R);
                int ci[3] = {c[0], c[1], c[2]};
                ci[d] = fi[d] / R[d];
                r.max_iface_diff = std::max(r.max_iface_diff, std::abs(f(i, j, k) - ca(ci[0], ci[1], ci[2])));
                ++r.nfaces_m0;
            });
        }
    }
    double x[4] = {r.max_div_diff, r.max_iface_diff, r.max_div_abs, 0.0};
    amrex::ParallelDescriptor::ReduceRealMax(x, 3);
    r.max_div_diff = x[0]; r.max_iface_diff = x[1]; r.max_div_abs = x[2];
    long nn[2] = {r.ncells, r.nfaces_m0};
    amrex::ParallelDescriptor::ReduceLongSum(nn, 2);
    r.ncells = nn[0]; r.nfaces_m0 = nn[1];
    return r;
}

void run_case(const std::string& name, amrex::IntVect cn, amrex::IntVect ratio, amrex::Box fregion_c, int split_dir, double scale_tol)
{
    for (int kind = 0; kind < 2; ++kind) {
        Setup s;
        make(s, cn, ratio, fregion_c, split_dir);
        fill_coarse(s, kind);
        for (int d = 0; d < 3; ++d) s.C[d].FillBoundary(s.cg.periodicity());
        const std::array<amrex::MultiFab*, 3> fine{&s.F[0], &s.F[1], &s.F[2]};
        const std::array<const amrex::MultiFab*, 3> crse{&s.C[0], &s.C[1], &s.C[2]};
        prolong_faces_level(fine, crse, s.cg, ratio);
        const Result r = evaluate(s);
        const std::string tag = name + (kind == 0 ? ", coarse field with divergence" : ", discretely divergence-free coarse field");
        // scale of the divergence: u / dx_coarse
        const double scale = 1.0 / s.cg.CellSize(0);
        if (amrex::ParallelDescriptor::IOProcessor())
            std::printf("  %-52s cells %ld  max|div_child - div_parent| %.2e  max|div_child| %.2e (scale %.1f)  interface/coarse-aligned faces %ld, max diff to coarse %.1e\n", tag.c_str(), r.ncells, r.max_div_diff, r.max_div_abs, scale,
                        r.nfaces_m0, r.max_iface_diff);
        CHECK_MSG(r.max_div_diff <= scale_tol * scale, tag + ": each fine cell has the divergence of its parent, got " + std::to_string(r.max_div_diff));
        CHECK_MSG(r.max_iface_diff == 0.0, tag + ": faces on a coarse face (interface faces included) carry exactly the coarse value");
        if (kind == 1) CHECK_MSG(r.max_div_abs <= scale_tol * scale, tag + ": divergence-free coarse gives divergence-free fine cells, got " + std::to_string(r.max_div_abs));
        // the reporting function: |div u - D| with D = 0 equals max|div_child| (kind 1: round-off)
        double rep = max_divergence_error_local({&s.F[0], &s.F[1], &s.F[2]}, s.fg, nullptr);
        amrex::ParallelDescriptor::ReduceRealMax(rep);
        CHECK_MSG(std::abs(rep - r.max_div_abs) <= 1e-12 * scale, tag + ": max_divergence_error_local agrees with the direct evaluation");

        // old fine data copied bitwise over the overlap: set a marker on the valid faces of the first fine box, prolong again with it as old data
        {
            std::array<amrex::MultiFab, 3> O;
            amrex::BoxList bl; bl.push_back(s.fba[0]);   // old level = first box only
            amrex::BoxArray oba(bl);
            amrex::DistributionMapping odm(oba);
            for (int d = 0; d < 3; ++d) {
                O[d].define(amrex::convert(oba, amrex::IntVect::TheDimensionVector(d)), odm, 1, 0);
                for (amrex::MFIter mfi(O[d]); mfi.isValid(); ++mfi) { auto a = O[d].array(mfi); amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { a(i, j, k) = 1000.0 + i + 31.0 * j + 977.0 * k + d; }); }
            }
            for (int d = 0; d < 3; ++d) s.F[d].setVal(-999.0);
            prolong_faces_level(fine, crse, s.cg, ratio, {&O[0], &O[1], &O[2]});
            long bad_old = 0, bad_new = 0;
            for (int d = 0; d < 3; ++d) {
                amrex::BoxArray ofaces = amrex::convert(oba, amrex::IntVect::TheDimensionVector(d));
                for (amrex::MFIter mfi(s.F[d]); mfi.isValid(); ++mfi) {
                    auto a = s.F[d].const_array(mfi);
                    amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                        const bool in_old = ofaces.contains(amrex::IntVect(i, j, k));
                        const double marker = 1000.0 + i + 31.0 * j + 977.0 * k + d;
                        if (in_old) { if (a(i, j, k) != marker) ++bad_old; }
                        else if (a(i, j, k) == -999.0 || (a(i, j, k) >= 900.0)) ++bad_new;
                    });
                }
            }
            amrex::ParallelDescriptor::ReduceLongSum(bad_old); amrex::ParallelDescriptor::ReduceLongSum(bad_new);
            CHECK_MSG(bad_old == 0 && bad_new == 0, tag + ": old fine faces copied bitwise over the overlap, new faces prolonged (bad " + std::to_string(bad_old) + "/" + std::to_string(bad_new) + ")");
        }

        // fine -> coarse mean: perturb the fine faces, the coarse face gets the mean of the fine faces lying on it
        {
            for (int d = 0; d < 3; ++d)
                for (amrex::MFIter mfi(s.F[d]); mfi.isValid(); ++mfi) { auto a = s.F[d].array(mfi); amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { a(i, j, k) = hrand(i, j, k, d, 3); }); }
            std::array<amrex::MultiFab, 3> C2;
            for (int d = 0; d < 3; ++d) { C2[d].define(s.C[d].boxArray(), s.cdm, 1, 0); C2[d].setVal(-5.0); }
            average_down_faces({&s.F[0], &s.F[1], &s.F[2]}, {&C2[0], &C2[1], &C2[2]}, ratio);
            long bad = 0, touched = 0;
            for (int d = 0; d < 3; ++d) {
                // reference: serial over a gathered copy of the fine faces on the owning rank of each coarse face is awkward; instead recompute the mean from F via ParallelCopy to a single-box array
                amrex::Box fb_all = amrex::convert(s.fg.Domain(), amrex::IntVect::TheDimensionVector(d));
                amrex::BoxArray one(fb_all);
                amrex::DistributionMapping odm(amrex::Vector<int>{amrex::ParallelDescriptor::IOProcessorNumber()});
                amrex::MultiFab G(one, odm, 1, 0);
                G.setVal(-77.0);
                G.ParallelCopy(s.F[d], 0, 0, 1, 0, 0);
                amrex::MultiFab Gc(amrex::convert(s.cba, amrex::IntVect::TheDimensionVector(d)), s.cdm, 1, 0);
                // compare on rank of IO: gather C2 as well
                amrex::BoxArray cone(amrex::convert(s.cdom, amrex::IntVect::TheDimensionVector(d)));
                amrex::MultiFab H(cone, odm, 1, 0);
                H.setVal(-88.0);
                H.ParallelCopy(C2[d], 0, 0, 1, 0, 0);
                if (amrex::ParallelDescriptor::IOProcessor()) {
                    auto g = G[0].const_array(); auto h = H[0].const_array();
                    amrex::LoopOnCpu(H[0].box(), [&](int i, int j, int k) {
                        if (h(i, j, k) == -5.0 || h(i, j, k) == -88.0) return;   // coarse face with no fine faces on it: untouched
                        ++touched;
                        const int ci[3] = {i, j, k};
                        int lo[3], hi[3];
                        for (int e = 0; e < 3; ++e) { lo[e] = ci[e] * ratio[e]; hi[e] = (e == d) ? lo[e] : lo[e] + ratio[e] - 1; }
                        double sum = 0; int n = 0;
                        for (int kk = lo[2]; kk <= hi[2]; ++kk) for (int jj = lo[1]; jj <= hi[1]; ++jj) for (int ii = lo[0]; ii <= hi[0]; ++ii) { sum += g(ii, jj, kk); ++n; }
                        if (std::abs(h(i, j, k) - sum / n) > 1e-15) ++bad;
                    });
                }
            }
            amrex::ParallelDescriptor::ReduceLongSum(bad); amrex::ParallelDescriptor::ReduceLongSum(touched);
            CHECK_MSG(touched > 0 && bad == 0, tag + ": average_down_faces gives the mean of the fine faces on each coarse face (faces touched " + std::to_string(touched) + ", bad " + std::to_string(bad) + ")");
        }

        // negative control: injection (every fine face = the coarse face value of the cell that contains it, no normal interpolation) breaks the parent divergence
        if (kind == 0) {
            for (int d = 0; d < 3; ++d) {
                std::array<amrex::MultiFab, 1> t;
                amrex::BoxArray cb = amrex::convert(amrex::BoxArray(s.fba).coarsen(ratio), amrex::IntVect::TheDimensionVector(d));
                (void)cb;
            }
            Setup bad;
            make(bad, cn, ratio, fregion_c, split_dir);
            fill_coarse(bad, 0);
            for (int d = 0; d < 3; ++d) bad.C[d].FillBoundary(bad.cg.periodicity());
            for (int d = 0; d < 3; ++d)
                for (amrex::MFIter mfi(bad.F[d]); mfi.isValid(); ++mfi) {
                    auto a = bad.F[d].array(mfi);
                    amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                        const amrex::IntVect c = amrex::coarsen(amrex::IntVect(i, j, k), ratio);
                        int ci[3] = {c[0], c[1], c[2]};
                        const int fi[3] = {i, j, k};
                        ci[d] = fi[d] / ratio[d];   // the coarse face at or below the fine face: injection
                        a(i, j, k) = hrand(((ci[0] % cn[0]) + cn[0]) % cn[0], ((ci[1] % cn[1]) + cn[1]) % cn[1], ((ci[2] % cn[2]) + cn[2]) % cn[2], d, 1) * 2.0;
                    });
                }
            const Result rb = evaluate(bad);
            CHECK_MSG(rb.max_div_diff > 1e-3 * scale, name + ": negative control, injection of the coarse face value breaks the parent divergence (" + std::to_string(rb.max_div_diff) + ")");
        }
    }
}

}  // namespace

int main(int, char** argv)
{
    int one = 1;
    amrex::Initialize(one, argv);
    {
        // 3-D, ratio 2: fine patch = coarse cells (2..9, 3..10, 4..11), two boxes split in x
        run_case("3-D ratio 2", amrex::IntVect(16), amrex::IntVect(2), amrex::Box(amrex::IntVect(2, 3, 4), amrex::IntVect(9, 10, 11)), 0, 1e-12);
        // mixed ratio (4,2,4), the patch touches the low x edge of the periodic domain
        run_case("3-D ratio (4,2,4)", amrex::IntVect(16), amrex::IntVect(4, 2, 4), amrex::Box(amrex::IntVect(0, 4, 4), amrex::IntVect(7, 11, 11)), 2, 1e-12);
        // 2-D: one cell in y (ratio 1 in y), ratio 2 in x and z
        run_case("2-D (1 cell in y) ratio (2,1,2)", amrex::IntVect(16, 1, 16), amrex::IntVect(2, 1, 2), amrex::Box(amrex::IntVect(3, 0, 3), amrex::IntVect(10, 0, 10)), 0, 1e-12);
        // 2-D ratio (4,1,4)
        run_case("2-D (1 cell in y) ratio (4,1,4)", amrex::IntVect(16, 1, 16), amrex::IntVect(4, 1, 4), amrex::Box(amrex::IntVect(4, 0, 2), amrex::IntVect(11, 0, 9)), 2, 1e-12);
    }
    const long nfail = fdstest::report("test_facetransfer");
    amrex::Finalize();
    return nfail == 0 ? 0 : 1;
}
