// probe_facedivfree.cpp (R4 step 1): what does the installed AMReX FaceDivFree interpolator do for the cases of the R4 design note?
// It calls face_divfree_interp.interp_arr directly on face-centred FArrayBoxes (one coarse region, one fine patch inside it), so no AmrCore,
// no registry and no driver code is involved. Numbers are printed; checks fail only if a documented expectation is violated.
//
//   probe_facedivfree                       all supported cases (3-D ratios 2 and 4, anisotropic 2x4x2, injection override, fallback)
//   probe_facedivfree --ratio-1-hidden      ratio (2,1,2) on a one-cell-thick domain: expected to ABORT inside AMReX (ctest WILL_FAIL)
//
// Index conventions: cell (i,j,k); the x-face (i,j,k) lies on the low side of cell (i,j,k); a coarse region of nc^3 cells, fine patch =
// refine(coarse cells [lo..hi]). The coarse fields come from a discrete curl of a hashed vector potential (exactly divergence-free) or from
// hashed face values (nonzero divergence).
#include <AMReX.H>
#include <AMReX_Interpolater.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>

#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>

#include "check.H"
#include "FaceTransfer.H"

using namespace amrex;

namespace {

double hash01(int i, int j, int k, int c)
{
    unsigned long long h = 1469598103934665603ULL;
    for (long long v : {static_cast<long long>(i), static_cast<long long>(j), static_cast<long long>(k), static_cast<long long>(c)}) {
        h ^= static_cast<unsigned long long>(v + 1000003); h *= 1099511628211ULL; h ^= h >> 29;
    }
    return static_cast<double>(h % 1000003ULL) / 1000003.0 - 0.5;
}

struct Result {
    double max_u = 0;            // max |face value| of the coarse field (scale)
    double div_fine_max = 0;     // max over fine cells of |div|
    double div_dev_max = 0;      // max over fine cells of |div(child) - div(parent)|
    double div_crse_max = 0;     // max over coarse cells of |div|
    double face_avg_dev = 0;     // max over coarse faces of |mean of fine faces on it - coarse value|
    double iface_dev = 0;        // max over the fine faces lying on the patch boundary of |fine - coarse value|
    int n_nan = 0;
};

enum class Kind { DivFree, Random };
enum class Method { AmrexAsIs, AmrexInjectInterface, NormalLinear };

// Face value of the coarse field. Divergence-free: discrete curl of the potential A (a hashed value at the edge of its component).
// hidden = one-cell direction y: potential Ay only, a function of (i,k), so v = 0 and nothing depends on j.
double crse_face(Kind kind, bool hidden, int d, int i, int j, int k, const Real* dx)
{
    if (kind == Kind::Random) return hash01(i, hidden ? 0 : j, k, 10 + d) + 0.3 * std::sin(0.4 * i + 0.2 * k);
    auto A = [&](int c, int a, int b, int e) { return hash01(a, hidden ? 0 : b, e, c) * (hidden ? (c == 1 ? 1.0 : 0.0) : 1.0); };
    if (d == 0) return (A(2, i, j + 1, k) - A(2, i, j, k)) / dx[1] - (A(1, i, j, k + 1) - A(1, i, j, k)) / dx[2];
    if (d == 1) return (A(0, i, j, k + 1) - A(0, i, j, k)) / dx[2] - (A(2, i + 1, j, k) - A(2, i, j, k)) / dx[0];
    return (A(1, i + 1, j, k) - A(1, i, j, k)) / dx[0] - (A(0, i, j + 1, k) - A(0, i, j, k)) / dx[1];
}

Result run(const IntVect& ratio, Kind kind, Method method, bool hidden)
{
    // coarse cells 0..n-1 in each direction (hidden: y has one cell); fine patch = coarse cells lo..hi, refined
    const IntVect nc = hidden ? IntVect(12, 1, 12) : IntVect(12, 12, 12);
    const Box cdom(IntVect(0), nc - 1);
    const Box cpatch(IntVect(3, 0, 3), hidden ? IntVect(8, 0, 8) : IntVect(8, 8, 8));
    const Box fpatch = refine(cpatch, ratio);
    const Real h = 1.0 / 12.0;
    const RealBox rb({0., 0., 0.}, {h * nc[0], h * nc[1], h * nc[2]});
    Geometry cgeom(cdom, rb, 0, {0, 0, 0});
    Geometry fgeom = amrex::refine(cgeom, ratio);
    const Real cdx[3] = {h, h, h};
    const Real fdx[3] = {h / ratio[0], h / ratio[1], h / ratio[2]};

    const Box cregion = grow(cpatch, 1);   // FaceDivFree::CoarseBox
    Array<FArrayBox, 3> cf, ff;
    Array<IArrayBox, 3> mk;
    Array<FArrayBox*, 3> cp, fp;
    Array<IArrayBox*, 3> mp{nullptr, nullptr, nullptr};
    Result r;
    for (int d = 0; d < 3; ++d) {
        const Box cb = convert(cregion, IntVect::TheDimensionVector(d));
        const Box fb = convert(fpatch, IntVect::TheDimensionVector(d));
        cf[d].resize(cb, 1); ff[d].resize(fb, 1); ff[d].setVal<RunOn::Host>(std::nan(""));
        auto a = cf[d].array();
        LoopOnCpu(cb, [&](int i, int j, int k) { a(i, j, k) = crse_face(kind, hidden, d, i, j, k, cdx); r.max_u = std::max(r.max_u, std::abs(a(i, j, k))); });
        cp[d] = &cf[d]; fp[d] = &ff[d];
    }
    if (method == Method::NormalLinear) {
        std::array<const FArrayBox*, 3> c{&cf[0], &cf[1], &cf[2]};
        std::array<FArrayBox*, 3> f{&ff[0], &ff[1], &ff[2]};
        fdsrt::prolong_faces_normal_linear(c, f, fpatch, ratio);
    } else {
        if (method == Method::AmrexInjectInterface) {
            // interface faces (on the boundary of the patch): fine value = coarse value, set before the call; mask = 0 there so AMReX skips them
            for (int d = 0; d < 3; ++d) {
                const Box cb = convert(cpatch, IntVect::TheDimensionVector(d));
                mk[d].resize(cb, 1); mk[d].setVal<RunOn::Host>(1);
                auto m = mk[d].array(); auto a = ff[d].array(); auto c = cf[d].const_array();
                LoopOnCpu(cb, [&](int i, int j, int k) {
                    const int idx[3] = {i, j, k};
                    const bool on_boundary = (idx[d] == cpatch.smallEnd(d)) || (idx[d] == cpatch.bigEnd(d) + 1);
                    if (!on_boundary) return;
                    m(i, j, k) = 0;
                    int lo[3] = {i * ratio[0], j * ratio[1], k * ratio[2]}; int hi[3] = {lo[0] + (d == 0 ? 0 : ratio[0] - 1), lo[1] + (d == 1 ? 0 : ratio[1] - 1), lo[2] + (d == 2 ? 0 : ratio[2] - 1)};
                    for (int kk = lo[2]; kk <= hi[2]; ++kk) for (int jj = lo[1]; jj <= hi[1]; ++jj) for (int ii = lo[0]; ii <= hi[0]; ++ii) a(ii, jj, kk) = c(i, j, k);
                });
                mp[d] = &mk[d];
            }
        }
        Vector<Array<BCRec, 3>> bcr;
        face_divfree_interp.interp_arr(cp, 0, fp, 0, 1, fpatch, ratio, mp, cgeom, fgeom, bcr, 0, 0, RunOn::Host);
    }
    for (int d = 0; d < 3; ++d) { const auto a = ff[d].const_array(); LoopOnCpu(ff[d].box(), [&](int i, int j, int k) { if (std::isnan(a(i, j, k))) ++r.n_nan; }); }

    auto fa = std::array<Array4<const Real>, 3>{ff[0].const_array(), ff[1].const_array(), ff[2].const_array()};
    auto ca = std::array<Array4<const Real>, 3>{cf[0].const_array(), cf[1].const_array(), cf[2].const_array()};
    auto divc = [&](int i, int j, int k) {
        return (ca[0](i + 1, j, k) - ca[0](i, j, k)) / cdx[0] + (ca[1](i, j + 1, k) - ca[1](i, j, k)) / cdx[1] + (ca[2](i, j, k + 1) - ca[2](i, j, k)) / cdx[2];
    };
    auto divf = [&](int i, int j, int k) {
        return (fa[0](i + 1, j, k) - fa[0](i, j, k)) / fdx[0] + (fa[1](i, j + 1, k) - fa[1](i, j, k)) / fdx[1] + (fa[2](i, j, k + 1) - fa[2](i, j, k)) / fdx[2];
    };
    LoopOnCpu(cpatch, [&](int i, int j, int k) { r.div_crse_max = std::max(r.div_crse_max, std::abs(divc(i, j, k))); });
    LoopOnCpu(fpatch, [&](int i, int j, int k) {
        const double df = divf(i, j, k);
        const double dc = divc(i / ratio[0], j / ratio[1], k / ratio[2]);   // fine indices are non-negative here
        r.div_fine_max = std::max(r.div_fine_max, std::abs(df));
        r.div_dev_max = std::max(r.div_dev_max, std::abs(df - dc));
    });
    // mean of the fine faces on each coarse face of the patch (including the boundary faces)
    for (int d = 0; d < 3; ++d) {
        LoopOnCpu(convert(cpatch, IntVect::TheDimensionVector(d)), [&](int i, int j, int k) {
            const int c[3] = {i, j, k};
            int lo[3], hi[3];
            for (int e = 0; e < 3; ++e) { lo[e] = c[e] * ratio[e]; hi[e] = (e == d) ? lo[e] : lo[e] + ratio[e] - 1; }
            double s = 0; int n = 0;
            for (int kk = lo[2]; kk <= hi[2]; ++kk) for (int jj = lo[1]; jj <= hi[1]; ++jj) for (int ii = lo[0]; ii <= hi[0]; ++ii) { s += fa[d](ii, jj, kk); ++n; }
            r.face_avg_dev = std::max(r.face_avg_dev, std::abs(s / n - ca[d](i, j, k)));
            if (c[d] == cpatch.smallEnd(d) || c[d] == cpatch.bigEnd(d) + 1)
                for (int kk = lo[2]; kk <= hi[2]; ++kk) for (int jj = lo[1]; jj <= hi[1]; ++jj) for (int ii = lo[0]; ii <= hi[0]; ++ii)
                    r.iface_dev = std::max(r.iface_dev, std::abs(fa[d](ii, jj, kk) - ca[d](i, j, k)));
        });
    }
    return r;
}

const char* mname(Method m) { return m == Method::AmrexAsIs ? "FaceDivFree as is" : m == Method::AmrexInjectInterface ? "FaceDivFree + interface injection" : "normal-linear (fallback)"; }

void report(const std::string& name, const IntVect& ratio, Kind kind, Method m, bool hidden, bool expect_fine_div_zero, bool expect_div_dev_zero, bool expect_iface_equal)
{
    const Result r = run(ratio, kind, m, hidden);
    const double s = r.max_u / (1.0 / 12.0);   // divergence scale |u|/dx
    if (amrex::ParallelDescriptor::IOProcessor())
        std::printf("  %-34s ratio (%d,%d,%d) %-9s %-7s: max|div c|/s %.2e  max|div f|/s %.2e  max|div f - div parent|/s %.2e  face-mean dev %.2e  interface dev %.2e  nan %d\n",
                    mname(m), ratio[0], ratio[1], ratio[2], kind == Kind::DivFree ? "div-free" : "random", hidden ? "hidden-y" : "3-D", r.div_crse_max / s, r.div_fine_max / s, r.div_dev_max / s, r.face_avg_dev, r.iface_dev, r.n_nan);
    CHECK_MSG(r.n_nan == 0, name + ": no NaN (all fine faces written)");
    CHECK_MSG(r.face_avg_dev < 1e-12 * r.max_u, name + ": mean of fine faces on a coarse face = coarse value");
    if (expect_fine_div_zero) CHECK_MSG(r.div_fine_max < 1e-10 * s, name + ": divergence-free coarse gives divergence-free fine");
    if (expect_div_dev_zero) CHECK_MSG(r.div_dev_max < 1e-10 * s, name + ": every child keeps the divergence of its parent");
    if (expect_iface_equal) CHECK_MSG(r.iface_dev < 1e-14 * r.max_u, name + ": interface faces equal the coarse value");
}

}  // namespace

int main(int argc, char** argv)
{
    int one = 1;
    amrex::Initialize(one, argv);
    const bool unsupported = argc > 1 && std::strcmp(argv[1], "--ratio-1-hidden") == 0;
    {
        if (unsupported) {
            // expected to abort in AMReX_Interpolater.cpp ("Only refinement ratio of 2 or 4 is supported")
            (void)run(IntVect(2, 1, 2), Kind::DivFree, Method::AmrexAsIs, true);
            if (amrex::ParallelDescriptor::IOProcessor()) std::printf("UNEXPECTED: ratio (2,1,2) did not abort\n");
            fdstest::counter().failures++;
        } else {
            for (const IntVect& r : {IntVect(2), IntVect(4), IntVect(2, 4, 2)}) {
                const std::string t = "ratio " + std::to_string(r[0]) + std::to_string(r[1]) + std::to_string(r[2]);
                report(t + " div-free", r, Kind::DivFree, Method::AmrexAsIs, false, true, true, false);
                report(t + " random", r, Kind::Random, Method::AmrexAsIs, false, false, true, false);
                report(t + " div-free + injection", r, Kind::DivFree, Method::AmrexInjectInterface, false, true, true, true);
                report(t + " random + injection", r, Kind::Random, Method::AmrexInjectInterface, false, false, true, true);
            }
            // fallback for the cases AMReX cannot do: coarse values on the faces, linear in the normal direction between them
            for (const IntVect& r : {IntVect(2), IntVect(4), IntVect(2, 4, 2)}) {
                const std::string t = "normal-linear ratio " + std::to_string(r[0]) + std::to_string(r[1]) + std::to_string(r[2]);
                report(t + " div-free", r, Kind::DivFree, Method::NormalLinear, false, true, true, true);
                report(t + " random", r, Kind::Random, Method::NormalLinear, false, false, true, true);
            }
            report("normal-linear hidden-y (2,1,2) div-free", IntVect(2, 1, 2), Kind::DivFree, Method::NormalLinear, true, true, true, true);
            report("normal-linear hidden-y (4,1,4) random", IntVect(4, 1, 4), Kind::Random, Method::NormalLinear, true, false, true, true);
        }
    }
    const long nfail = fdstest::report(unsupported ? "probe_facedivfree --ratio-1-hidden" : "probe_facedivfree");
    amrex::Finalize();
    return nfail == 0 ? 0 : 1;
}
