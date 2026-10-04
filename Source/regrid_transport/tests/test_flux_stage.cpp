// test_flux_stage.cpp (R2b): the interface flux overwrite end to end on the mock transport code (tests/mock_transport.H), run on 1 and 4 ranks.
//  A  static two-level 2-D case (ns2d_16_int_1to2 layout: 16x1x16 coarse cells, a refined square patch in the middle, ratio (2,1,2), the patch cut into 4 boxes so that box-box
//     faces on the fine level and coarse faces shared by coarse boxes are exercised). A blob that straddles the interface is advected (divergence-free rotation, face velocity from a
//     stream function so that the coarse face velocity equals the mean of the fine ones) and diffused. Composite mass of every species conserved to round-off with the overwrite ON;
//     NEGATIVE CONTROL: with it OFF the composite mass drifts by many orders of magnitude more.
//  B  uniform state stays uniform across the interface (round-off), overwrite ON.
//  C  three levels (ratio 2 and 2), finest first ordering: conservation to round-off, and a control with the overwrite OFF.
//  D  3-D case with a patch that has corners, constant translation: conservation, control.
//  E  full-domain level 1 = no interface faces: the run is bitwise the single-level run at the fine resolution (empty override set changes nothing); the number of entries is 0.
//  F  D-059 shared ghost cells: the count (printed) agrees with the conflicts counted by the real FDS-rule fill (LevelOps), for a square patch (4 corners), an L shape, a thin patch.
// The mock has no second ghost layer and no FDS ghost rules; those are checked by fr016_ghost_check and by the driver tests (README).
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>

#include "FluxOverrideOps.H"
#include "FluxStageRunner.H"
#include "GhostShare.H"
#include "LevelOps.H"
#include "check.H"
#include "mock_transport.H"

using namespace fdsrt;

namespace {

constexpr double PI = 3.14159265358979323846;

std::string sci(double v) { char b[32]; std::snprintf(b, sizeof b, "%.2e", v); return b; }

amrex::Geometry geom_of(const amrex::Box& dom, double L) { return amrex::Geometry(dom, amrex::RealBox({0., 0., 0.}, {L, L, L}), 0, {1, 1, 1}); }

// stream function of the 2-D rotation (x, z plane); periodic on [0,1]^2
double psi(double x, double z) { return std::sin(2 * PI * x) * std::sin(2 * PI * z) / (2 * PI); }

// face velocity from the stream function at the corners of the face (exactly divergence free on each level, coarse face = mean of the fine faces)
std::array<double, 3> psi_vel(int, int d, double x, double y, double z, const double* dx)
{
    // (x, y, z) is the CENTRE of the face; the face spans dz in z (for d = 0) or dx in x (for d = 2)
    (void)y;
    std::array<double, 3> v{0, 0, 0};
    if (d == 0) v[0] = (psi(x, z + 0.5 * dx[2]) - psi(x, z - 0.5 * dx[2])) / dx[2];
    if (d == 2) v[2] = -(psi(x + 0.5 * dx[0], z) - psi(x - 0.5 * dx[0], z)) / dx[0];
    return v;
}

double blob(double x, double y, double z, int n)
{
    (void)y;
    const double cx = 0.5 - 0.25 + 0.03 * n, cz = 0.5, s = 0.09;   // centred near the left patch edge (x = 0.25): straddles the interface
    return 1.0 + (1.0 + 0.5 * n) * std::exp(-((x - cx) * (x - cx) + (z - cz) * (z - cz)) / (s * s));
}

struct Hier2D {
    std::vector<mock::Level> lv;
};

// 2-D hierarchy: level 0 = 16x1x16 on [0,1]^3, maxSize 8 (4 boxes); level 1 = ratio (2,1,2) over coarse cells [4..11]x[4..11] cut into boxes of 8 (4 boxes); optional level 2 over
// level-1 cells [10..21]x[10..21] (inside the level-1 patch, which spans 8..23) with ratio (2,1,2), cut into boxes of 6.
void build_2d(mock::Transport& T, int nlev, bool full_domain)
{
    amrex::Box d0(amrex::IntVect(0, 0, 0), amrex::IntVect(15, 0, 15));
    T.lv.clear(); T.lv.resize(nlev);
    T.lv[0].geom = geom_of(d0, 1.0);
    T.lv[0].geom = amrex::Geometry(d0, amrex::RealBox({0., 0., 0.}, {1., 1. / 16, 1.}), 0, {1, 1, 1});
    T.lv[0].ba = amrex::BoxArray(d0);
    T.lv[0].ba.maxSize(amrex::IntVect(8, 1, 8));
    for (int l = 1; l < nlev; ++l) {
        const amrex::IntVect r(2, 1, 2);
        mock::Level& L = T.lv[l];
        L.ratio = r;
        L.geom = amrex::refine(T.lv[l - 1].geom, r);
        amrex::Box cb;
        if (l == 1) cb = full_domain ? d0 : amrex::Box(amrex::IntVect(4, 0, 4), amrex::IntVect(11, 0, 11));
        else cb = amrex::Box(amrex::IntVect(10, 0, 10), amrex::IntVect(21, 0, 21));
        L.ba = amrex::BoxArray(amrex::refine(cb, r));
        L.ba.maxSize(l == 1 ? amrex::IntVect(8, 1, 8) : amrex::IntVect(6, 1, 6));
    }
}

struct Result {
    double drift[4] = {0, 0, 0, 0};   // max over species of |M - M0| / M0 at the end, per species
    double maxdrift = 0;
    long entries = 0;
};

Result run(mock::Transport& T, int nc, bool overwrite, int nsteps, double D, const std::function<double(double, double, double, int)>& q0,
           const std::function<std::array<double, 3>(int, int, double, double, double, const double*)>& vel, double cfl_dt)
{
    std::vector<fdsrt::StageLevel> sl;
    for (auto& L : T.lv) { fdsrt::StageLevel s; s.geom = L.geom; s.ba = L.ba; s.dm = L.dm.empty() ? amrex::DistributionMapping(L.ba) : L.dm; s.ratio_from_parent = L.ratio; sl.push_back(s); }
    for (size_t l = 0; l < T.lv.size(); ++l) if (T.lv[l].dm.empty()) T.lv[l].dm = sl[l].dm;
    fdsrt::FluxStageRunner runner(sl, overwrite);
    T.runner = &runner;
    T.init(nc, D, q0, vel);
    long double m0[8];
    for (int n = 0; n < nc; ++n) m0[n] = T.composite(n);
    for (int s = 0; s < nsteps; ++s) T.step(cfl_dt);
    Result r;
    for (int n = 0; n < nc; ++n) {
        const double dm = static_cast<double>(std::abs(T.composite(n) - m0[n]) / std::abs(m0[n]));
        if (n < 4) r.drift[n] = dm;
        r.maxdrift = std::max(r.maxdrift, dm);
    }
    r.entries = runner.stats().adv_entries + runner.stats().dif_entries;
    amrex::ParallelAllReduce::Sum(r.entries, amrex::ParallelContext::CommunicatorSub());
    T.runner = nullptr;
    return r;
}

void test_A_B()
{
    const double dt = 0.2 * (1.0 / 32) / 1.0;   // fine dx = 1/32, |u| <= 1/(2pi) * 2pi = 1 ... stream function amplitude gives |u| <= 1
    {
        mock::Transport T; build_2d(T, 2, false);
        const Result on = run(T, 3, true, 25, 0.01, blob, psi_vel, dt);
        mock::Transport T2; build_2d(T2, 2, false);
        const Result off = run(T2, 3, false, 25, 0.01, blob, psi_vel, dt);
        if (amrex::ParallelDescriptor::IOProcessor())
            std::printf("  A two-level blob, 25 steps: composite mass drift overwrite ON %s (entries %ld), OFF %s\n", sci(on.maxdrift).c_str(), on.entries, sci(off.maxdrift).c_str());
        CHECK_MSG(on.maxdrift < 1e-13, sci(on.maxdrift));
        CHECK_MSG(on.entries > 0, "no override entries");
        CHECK_MSG(off.maxdrift > 1e3 * std::max(on.maxdrift, 1e-16) && off.maxdrift > 1e-7, sci(off.maxdrift));   // negative control: the test can see a failure
    }
    {   // B: uniform state
        mock::Transport T; build_2d(T, 2, false);
        auto uni = [](double, double, double, int n) { return 1.0 + 0.25 * n; };
        std::vector<fdsrt::StageLevel> sl;
        const Result r = run(T, 3, true, 25, 0.01, uni, psi_vel, dt);
        double mx = 0;
        // run() reset the runner; read the state
        for (int n = 0; n < 3; ++n) mx = std::max(mx, T.max_abs_diff_from(1.0 + 0.25 * n, n));
        if (amrex::ParallelDescriptor::IOProcessor()) std::printf("  B uniform state after 25 steps: max deviation %s (conservation drift %s)\n", sci(mx).c_str(), sci(r.maxdrift).c_str());
        CHECK_MSG(mx < 1e-13, sci(mx));
        CHECK_MSG(r.entries > 0, "no interface entries in the uniform case");
    }
}

void test_C()
{
    const double dt = 0.2 * (1.0 / 64);
    mock::Transport T; build_2d(T, 3, false);
    auto q0 = [](double x, double y, double z, int n) { return blob(x, y, z, n) + 0.5 * std::exp(-((x - 0.5) * (x - 0.5) + (z - 0.5) * (z - 0.5)) / 0.002); };
    const Result on = run(T, 2, true, 20, 0.01, q0, psi_vel, dt);
    mock::Transport T2; build_2d(T2, 3, false);
    const Result off = run(T2, 2, false, 20, 0.01, q0, psi_vel, dt);
    if (amrex::ParallelDescriptor::IOProcessor()) std::printf("  C three levels, 20 steps: drift ON %s (entries %ld), OFF %s\n", sci(on.maxdrift).c_str(), on.entries, sci(off.maxdrift).c_str());
    CHECK_MSG(on.maxdrift < 1e-13, sci(on.maxdrift));
    CHECK_MSG(off.maxdrift > 1e-7, sci(off.maxdrift));
}

void test_D()
{
    // 3-D: coarse 12^3 on [0,1]^3, patch over coarse cells [3..8]^3 with an extra tilted piece so that edges and corners exist; constant velocity
    mock::Transport T, T2;
    for (mock::Transport* P : {&T, &T2}) {
        amrex::Box d0(amrex::IntVect(0), amrex::IntVect(11));
        P->lv.clear(); P->lv.resize(2);
        P->lv[0].geom = geom_of(d0, 1.0);
        P->lv[0].ba = amrex::BoxArray(d0);
        P->lv[0].ba.maxSize(6);
        mock::Level& L = P->lv[1];
        L.ratio = amrex::IntVect(2);
        L.geom = amrex::refine(P->lv[0].geom, L.ratio);
        amrex::BoxList bl;
        bl.push_back(amrex::refine(amrex::Box(amrex::IntVect(3), amrex::IntVect(8)), L.ratio));
        bl.push_back(amrex::refine(amrex::Box(amrex::IntVect(9, 3, 3), amrex::IntVect(10, 5, 5)), L.ratio));   // a lug: more edges and corners
        L.ba = amrex::BoxArray(bl);
        L.ba.maxSize(6);
    }
    auto q0 = [](double x, double y, double z, int n) {
        return 1.0 + (1 + n) * std::exp(-((x - 0.55) * (x - 0.55) + (y - 0.5) * (y - 0.5) + (z - 0.5) * (z - 0.5)) / (0.1 * 0.1));
    };
    auto vel = [](int, int d, double, double, double, const double*) { std::array<double, 3> v{0, 0, 0}; v[d] = (d == 0 ? 0.7 : d == 1 ? -0.4 : 0.3); return v; };
    const double dt = 0.2 * (1.0 / 24) / 0.7;
    const Result on = run(T, 2, true, 12, 0.005, q0, vel, dt);
    const Result off = run(T2, 2, false, 12, 0.005, q0, vel, dt);
    if (amrex::ParallelDescriptor::IOProcessor()) std::printf("  D 3-D patch with lug, 12 steps: drift ON %s (entries %ld), OFF %s\n", sci(on.maxdrift).c_str(), on.entries, sci(off.maxdrift).c_str());
    CHECK_MSG(on.maxdrift < 1e-13, sci(on.maxdrift));
    CHECK_MSG(off.maxdrift > 1e-8, sci(off.maxdrift));
}

void test_E()
{
    // level 1 covers the whole domain: no interface; compare with a single-level run at the fine resolution, bitwise
    const double dt = 0.2 * (1.0 / 32);
    mock::Transport A; build_2d(A, 2, true);
    const Result ra = run(A, 3, true, 10, 0.01, blob, psi_vel, dt);
    mock::Transport B;
    {
        amrex::Box d1(amrex::IntVect(0, 0, 0), amrex::IntVect(31, 0, 31));
        B.lv.clear(); B.lv.resize(1);
        B.lv[0].geom = amrex::Geometry(d1, amrex::RealBox({0., 0., 0.}, {1., 1. / 16, 1.}), 0, {1, 1, 1});
        B.lv[0].ba = A.lv[1].ba;   // same boxes and mapping as the fine level of A
        B.lv[0].dm = A.lv[1].dm;
    }
    // psi_vel is called with the level number; it does not depend on it
    const Result rb = run(B, 3, true, 10, 0.01, blob, psi_vel, dt);
    bool same = true;
    for (amrex::MFIter mfi(A.lv[1].S[0]); mfi.isValid(); ++mfi) {
        auto a = A.lv[1].S[0].const_array(mfi);
        auto b = B.lv[0].S[0].const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { for (int n = 0; n < 3; ++n) if (a(i, j, k, n) != b(i, j, k, n)) same = false; });
    }
    int s = same ? 1 : 0;
    amrex::ParallelAllReduce::Min(s, amrex::ParallelContext::CommunicatorSub());
    if (amrex::ParallelDescriptor::IOProcessor()) std::printf("  E full-domain level 1: override entries %ld, fine-level state bitwise equal to the single-level fine run: %s (drift %s vs %s)\n", ra.entries, s ? "yes" : "NO", sci(ra.maxdrift).c_str(), sci(rb.maxdrift).c_str());
    CHECK_MSG(ra.entries == 0, "interface faces found with a full-domain level 1");
    CHECK_MSG(s == 1, "level 1 differs from the single-level run");
}

// G/H: patch touching periodic domain edges (interface faces at both the low and the high side of the domain), and the species-sum property of the override
void test_G_H()
{
    const double dt = 0.2 * (1.0 / 32);
    for (int use_overwrite = 1; use_overwrite >= 0; --use_overwrite) {
        mock::Transport T; build_2d(T, 2, false);
        // replace the patch by one that touches the low x edge and the high z edge of the periodic domain: coarse cells [0..3] x [12..15]
        T.lv[1].ba = amrex::BoxArray(amrex::refine(amrex::Box(amrex::IntVect(0, 0, 12), amrex::IntVect(3, 0, 15)), amrex::IntVect(2, 1, 2)));
        T.lv[1].ba.maxSize(amrex::IntVect(4, 1, 4));
        auto q0 = [](double x, double y, double z, int n) {   // blob across the periodic corner (x = 0, z = 1)
            (void)y;
            auto wrap = [](double d) { return d - std::round(d); };
            const double dx = wrap(x - 0.0), dz = wrap(z - 1.0);
            return 1.0 + (1.0 + n) * std::exp(-(dx * dx + dz * dz) / (0.1 * 0.1));
        };
        const Result r = run(T, 2, use_overwrite == 1, 20, 0.01, q0, psi_vel, dt);
        if (amrex::ParallelDescriptor::IOProcessor()) std::printf("  H patch on periodic edges, overwrite %s: drift %s (entries %ld)\n", use_overwrite ? "ON " : "OFF", sci(r.maxdrift).c_str(), r.entries);
        if (use_overwrite) { CHECK_MSG(r.maxdrift < 1e-13, sci(r.maxdrift)); CHECK(r.entries > 0); }
        else CHECK_MSG(r.maxdrift > 1e-7, sci(r.maxdrift));
    }
    {   // G: the override of a face whose fine fluxes sum to zero over the species sums to zero (the species-sum fix survives the area mean), and equals the plain mean
        amrex::Box d0(amrex::IntVect(0, 0, 0), amrex::IntVect(7, 0, 7));
        amrex::Geometry cg = amrex::Geometry(d0, amrex::RealBox({0., 0., 0.}, {1., 0.125, 1.}), 0, {1, 1, 1});
        amrex::BoxArray cba(d0);
        cba.maxSize(amrex::IntVect(4, 1, 4));
        amrex::DistributionMapping cdm(cba);
        const amrex::IntVect r(2, 1, 2);
        amrex::BoxArray fba(amrex::refine(amrex::Box(amrex::IntVect(2, 0, 2), amrex::IntVect(5, 0, 5)), r));
        fba.maxSize(amrex::IntVect(4, 1, 4));
        amrex::DistributionMapping fdm(fba);
        amrex::MultiFab ff[3];
        for (int d = 0; d < 3; ++d) {
            amrex::IntVect nodal(0); nodal[d] = 1;
            ff[d].define(amrex::convert(fba, nodal), fdm, 3, 0);
            for (amrex::MFIter mfi(ff[d]); mfi.isValid(); ++mfi) {
                auto a = ff[d].array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                    a(i, j, k, 0) = std::sin(1.3 * i + 0.7 * k + d);
                    a(i, j, k, 1) = std::cos(0.9 * i - 0.4 * k + 2 * d);
                    a(i, j, k, 2) = -(a(i, j, k, 0) + a(i, j, k, 1));
                });
            }
        }
        const amrex::MultiFab* fl[3] = {&ff[0], &ff[1], &ff[2]};
        OverrideStats st;
        auto lists = build_flux_overrides(cba, cdm, cg, fba, fdm, r, fl, &st);
        double worst_sum = 0, worst_mean = 0;
        long nface = 0;
        int ib = 0;
        for (amrex::MFIter mfi(cba, cdm); mfi.isValid(); ++mfi, ++ib)
            for (const FluxOverride& o : lists[ib])
                for (size_t f = 0; f < o.face.size(); ++f) {
                    ++nface;
                    worst_sum = std::max(worst_sum, std::abs(o.value[f * 3] + o.value[f * 3 + 1] + o.value[f * 3 + 2]));
                    // independent mean of the fine faces tiling this coarse face (face index a is the low face of cell a; fine face a*r_d in the normal direction)
                    const int d = o.dir;
                    int lo[3], hi[3];
                    const int cc[3] = {o.face[f][0], o.face[f][1], o.face[f][2]};
                    for (int e = 0; e < 3; ++e) { lo[e] = cc[e] * r[e]; hi[e] = lo[e] + r[e] - 1; }
                    lo[d] = hi[d] = cc[d] * r[d];
                    double m0 = 0; int cnt = 0;
                    for (int k = lo[2]; k <= hi[2]; ++k) for (int j = lo[1]; j <= hi[1]; ++j) for (int i = lo[0]; i <= hi[0]; ++i) {
                        // find the fine value: any fine box holding the face (global lookup through the analytic definition used above)
                        m0 += std::sin(1.3 * i + 0.7 * k + d); ++cnt;
                    }
                    worst_mean = std::max(worst_mean, std::abs(o.value[f * 3] - m0 / cnt));
                }
        long tot = nface;
        amrex::ParallelAllReduce::Sum(tot, amrex::ParallelContext::CommunicatorSub());
        amrex::ParallelAllReduce::Max(worst_sum, amrex::ParallelContext::CommunicatorSub());
        amrex::ParallelAllReduce::Max(worst_mean, amrex::ParallelContext::CommunicatorSub());
        if (amrex::ParallelDescriptor::IOProcessor()) std::printf("  G override faces %ld: max |species sum| %s, max deviation from the independent mean %s\n", tot, sci(worst_sum).c_str(), sci(worst_mean).c_str());
        CHECK(tot > 0);
        CHECK_MSG(worst_sum < 1e-15, sci(worst_sum));
        CHECK_MSG(worst_mean < 1e-15, sci(worst_mean));
    }
}

// D-059: shared ghost cells
long conflicts_from_levelops(const amrex::BoxArray& fine_ba, const amrex::IntVect& r, const amrex::Geometry& cg, const amrex::BoxArray& cba)
{
    amrex::DistributionMapping cdm(cba), fdm(fine_ba);
    amrex::Geometry fg = amrex::refine(cg, r);
    amrex::MultiFab crho(cba, cdm, 1, 2), frho(fine_ba, fdm, 1, 2);
    crho.setVal(1.0); frho.setVal(1.0);
    ScalarStage cs, fs;
    cs.rho = &crho; fs.rho = &frho;
    const CfStats st = fill_covered_ghosts_fds(cs, fs, fg, cg, r, nullptr);
    long c = st.conflicts;
    amrex::ParallelAllReduce::Sum(c, amrex::ParallelContext::CommunicatorSub());
    return c;
}

void test_F()
{
    amrex::Box d0(amrex::IntVect(0), amrex::IntVect(15));
    amrex::Geometry cg = geom_of(d0, 1.0);
    amrex::BoxArray cba(d0);
    cba.maxSize(8);
    const amrex::IntVect r(2);
    struct Case { const char* name; std::vector<amrex::Box> boxes; long expect_two_faces; };   // boxes in COARSE cells
    std::vector<Case> cases;
    cases.push_back({"square 8^3 patch", {amrex::Box(amrex::IntVect(4), amrex::IntVect(11))}, -1});
    cases.push_back({"L shape", {amrex::Box(amrex::IntVect(2, 2, 4), amrex::IntVect(9, 9, 11)), amrex::Box(amrex::IntVect(10, 2, 4), amrex::IntVect(13, 5, 11))}, -1});
    cases.push_back({"thin patch (3 coarse cells)", {amrex::Box(amrex::IntVect(5, 5, 5), amrex::IntVect(7, 10, 10))}, -1});
    for (const Case& c : cases) {
        amrex::BoxList bl;
        for (const amrex::Box& b : c.boxes) bl.push_back(amrex::refine(b, r));
        amrex::BoxArray fba(bl);
        fba.maxSize(8);
        const GhostShareCount n = count_shared_ghost_cells(fba, r, cg);
        const long lo = conflicts_from_levelops(fba, r, cg, cba);
        if (amrex::ParallelDescriptor::IOProcessor()) std::printf("  F D-059 %-28s %s | conflicts counted by the FDS-rule fill: %ld\n", c.name, n.to_string().c_str(), lo);
        CHECK_MSG(n.serving_several == lo, std::string(c.name) + ": " + std::to_string(n.serving_several) + " vs " + std::to_string(lo));
        CHECK(n.serving_several > 0);
    }
    // a patch of 8^3 coarse cells in 3-D: the covered cells serving more than one face are the edge and corner cells of the outer two layers; the count is deterministic
    amrex::BoxArray fba(amrex::refine(amrex::Box(amrex::IntVect(4), amrex::IntVect(11)), r));
    const GhostShareCount sq = count_shared_ghost_cells(fba, r, cg);
    CHECK(sq.serving_two_faces > 0);
}

}  // namespace

int main(int argc, char** argv)
{
    amrex::Initialize(argc, argv);
    {
        test_A_B();
        test_C();
        test_D();
        test_E();
        test_G_H();
        test_F();
    }
    const long nfail = fdstest::report("flux_stage");
    amrex::Finalize();
    return nfail == 0 ? 0 : 1;
}
