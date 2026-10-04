// loops_guards.cpp: pb_fds_loops_guards, guard and edge-index tests of the host loop translations (FdsPressureLoops.cpp).
//   C1  periodic wrap of the H boundary fill (L1220-L1222) needs the box to be the whole domain extent (throws otherwise);
//   C2  TUNNEL_PRECONDITIONER is refused (NotBuilt) by the H fill and by L1209;
//   C3  L1209 wind arrays U_WIND/V_WIND/W_WIND must cover 0..KBP1; z walls have KK=0 and KK=KBP1 and read the wind at those indices.
// Every check prints PASS/FAIL; exit 0 only if all pass. Mutants (PB_FDSLOOPS_MUTANT 19, 20, 21 in FdsPressureLoops.cpp: no whole-extent
// guard, the former wind guard 1..KBAR, no tunnel refusal) must make this program fail (tests/m5_fds_loops.cmake).
// Built also with -fsanitize=address (pb_fds_loops_guards_asan): the C3 cases then also prove that no read leaves the arrays.
// No AMReX, no MPI.
#include "FdsPressureLoops.H"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

using namespace pressure_backend::fdsloops;

namespace {
int g_fail = 0, g_n = 0;
void check (bool ok, const char* what) { ++g_n; if (!ok) ++g_fail; std::printf("  %s: %s\n", ok ? "PASS" : "FAIL", what); }

bool same_bits (double a, double b) { return std::memcmp(&a, &b, 8) == 0; }

template <class Ex, class Fn> bool throws_as (Fn fn) { try { fn(); } catch (Ex const&) { return true; } catch (...) { return false; } return false; }

struct Arr3 {   // owned Fortran-ordered array with explicit lower bounds
    std::vector<double> v; F3 f;
    Arr3 (int l0, int h0, int l1, int h1, int l2, int h2, double fill = 0.0) {
        f.lo[0] = l0; f.lo[1] = l1; f.lo[2] = l2; f.n[0] = h0 - l0 + 1; f.n[1] = h1 - l1 + 1; f.n[2] = h2 - l2 + 1;
        v.assign(static_cast<std::size_t>(f.n[0]) * f.n[1] * f.n[2], fill); f.p = v.data();
    }
    void fill_pattern (double s) { for (std::size_t q = 0; q < v.size(); ++q) v[q] = s * (1.0 + 0.37 * std::sin(0.9 * static_cast<double>(q) + 0.3)); }
};
struct Arr2 { std::vector<double> v; F2 f; Arr2 (int n0, int n1, double fill = 0.0) { v.assign(static_cast<std::size_t>(n0) * n1, fill); f.p = v.data(); f.lo[0] = 1; f.lo[1] = 1; f.n[0] = n0; f.n[1] = n1; } };
struct Arr1 { std::vector<double> v; F1 f; Arr1 (int lo, int hi, double fill = 0.0) { v.assign(hi - lo + 1, fill); f.p = v.data(); f.lo = lo; f.n = hi - lo + 1; } };

// ---------------------------------------------------------------------------------------------------------------- C1
void test_c1 ()
{
    std::printf("C1 periodic wrap needs the whole domain extent\n");
    const int N = 6;
    const Box3 dom = fds_interior(N, N, N);
    Arr2 bxs(N, N), bxf(N, N), bys(N, N), byf(N, N), bzs(N, N), bzf(N, N);
    for (int dir = 0; dir < 3; ++dir) {
        const char* nm = dir == 0 ? "x" : dir == 1 ? "y" : "z";
        HFillOptions opt; opt.domain = dom;
        auto fill = [&](Box3 const& b, int code, Arr3& h) {
            if (dir == 0) pres_h_bc_x(b, opt, code, 0.5, bxs.f, bxf.f, h.f);
            else if (dir == 1) pres_h_bc_y(b, opt, code, 0.5, bys.f, byf.f, h.f);
            else pres_h_bc_z(b, opt, code, 0.5, bzs.f, bzf.f, h.f);
        };
        Arr3 h(0, N + 1, 0, N + 1, 0, N + 1); h.fill_pattern(2.0);
        // whole extent: wraps; the ghost layers are the opposite interior layers
        fill(dom, 0, h);
        bool wrapped = true;
        for (int a = 1; a <= N; ++a) for (int c = 1; c <= N; ++c) {
            auto at = [&](int s, int p, int q) { return dir == 0 ? h.f(s, p, q) : dir == 1 ? h.f(p, s, q) : h.f(p, q, s); };
            wrapped = wrapped && same_bits(at(0, a, c), at(N, a, c)) && same_bits(at(N + 1, a, c), at(1, a, c));
        }
        check(wrapped, (std::string("whole extent wraps in ") + nm).c_str());
        // a box short of the domain at the high end, at the low end, and a box larger than the domain: all refused, nothing written
        Box3 hi_short = dom; hi_short.hi[dir] = N - 1;
        Box3 lo_short = dom; lo_short.lo[dir] = 2;
        Arr3 h2(0, N + 1, 0, N + 1, 0, N + 1); h2.fill_pattern(2.0);
        const std::vector<double> before = h2.v;
        check(throws_as<std::invalid_argument>([&] { fill(hi_short, 0, h2); (void)0; return 0; }) , (std::string("sub-box (high end cut) refused in ") + nm).c_str());
        check(throws_as<std::invalid_argument>([&] { fill(lo_short, 0, h2); return 0; }), (std::string("sub-box (low end cut) refused in ") + nm).c_str());
        check(h2.v == before, (std::string("refused call wrote nothing in ") + nm).c_str());
        // a non-periodic code on a sub-box is not a wrap and stays allowed (physical-face fill of a piece of a face)
        bool ok = true; try { fill(hi_short, 3, h2); } catch (...) { ok = false; }
        check(ok, (std::string("non-periodic code on a sub-box is allowed in ") + nm).c_str());
    }
}

// ---------------------------------------------------------------------------------------------------------------- C2
// A small valid L1209 context (closed Neumann walls only) for the tunnel check and the wind cases.
struct Ctx {
    int ib, jb, kb;
    Arr3 hp, kres, fvx, fvy, fvz, uu, vv, ww;
    Arr1 hx, hy, hz, dx, dy, dz, rdxn, rdyn, rdzn, uw, vw, ww_;
    Arr2 bxs, bxf, bys, byf, bzs, bzf;
    PoissonBcContext c;
    PoissonVent vent;
    Ctx (int i, int j, int k, int wind_lo = 0, int wind_hi = -1)
        : ib(i), jb(j), kb(k),
          hp(0, i + 1, 0, j + 1, 0, k + 1), kres(0, i + 1, 0, j + 1, 0, k + 1), fvx(0, i + 1, 0, j + 1, 0, k + 1), fvy(0, i + 1, 0, j + 1, 0, k + 1), fvz(0, i + 1, 0, j + 1, 0, k + 1),
          uu(-1, i + 1, 0, j + 1, 0, k + 1), vv(0, i + 1, -1, j + 1, 0, k + 1), ww(0, i + 1, 0, j + 1, -1, k + 1),
          hx(0, i + 1), hy(0, j + 1), hz(0, k + 1), dx(1, i), dy(1, j), dz(1, k), rdxn(0, i), rdyn(0, j), rdzn(0, k),
          uw(wind_lo, wind_hi < 0 ? k + 1 : wind_hi), vw(wind_lo, wind_hi < 0 ? k + 1 : wind_hi), ww_(wind_lo, wind_hi < 0 ? k + 1 : wind_hi),
          bxs(j, k), bxf(j, k), bys(i, k), byf(i, k), bzs(i, j), bzf(i, j)
    {
        hp.fill_pattern(3.0); kres.fill_pattern(0.5); fvx.fill_pattern(0.2); fvy.fill_pattern(0.3); fvz.fill_pattern(0.4);
        uu.fill_pattern(1.0); vv.fill_pattern(1.1); ww.fill_pattern(1.2);
        for (int q = 0; q <= i + 1; ++q) hx.v[q] = 0.9 + 0.01 * q;
        for (int q = 0; q <= j + 1; ++q) hy.v[q] = 0.8 + 0.01 * q;
        for (int q = 0; q <= k + 1; ++q) hz.v[q] = 0.7 + 0.01 * q;
        for (int q = 0; q < i; ++q) dx.v[q] = 0.1 + 0.001 * q;
        for (int q = 0; q < j; ++q) dy.v[q] = 0.12 + 0.001 * q;
        for (int q = 0; q < k; ++q) dz.v[q] = 0.14 + 0.001 * q;
        for (int q = 0; q <= i; ++q) rdxn.v[q] = 9.0 + 0.1 * q;
        for (int q = 0; q <= j; ++q) rdyn.v[q] = 8.0 + 0.1 * q;
        for (int q = 0; q <= k; ++q) rdzn.v[q] = 7.0 + 0.1 * q;
        for (std::size_t q = 0; q < uw.v.size(); ++q) { uw.v[q] = 0.11 * (q + 1); vw.v[q] = 0.23 * (q + 1); ww_.v[q] = 0.37 * (q + 1); }
        c.ibar = i; c.jbar = j; c.kbar = k;
        c.hp = hp.f; c.kres = kres.f; c.fvx = fvx.f; c.fvy = fvy.f; c.fvz = fvz.f; c.uu = uu.f; c.vv = vv.f; c.ww = ww.f;
        c.hx = hx.f; c.hy = hy.f; c.hz = hz.f; c.dx = dx.f; c.dy = dy.f; c.dz = dz.f; c.rdxn = rdxn.f; c.rdyn = rdyn.f; c.rdzn = rdzn.f;
        c.u_wind = uw.f; c.v_wind = vw.f; c.w_wind = ww_.f;
        c.u0 = 0.1; c.v0 = 0.2; c.w0 = 0.3; c.t = 2.0; c.dt = 0.01; c.t_begin = 0.0;
        c.bxs = bxs.f; c.bxf = bxf.f; c.bys = bys.f; c.byf = byf.f; c.bzs = bzs.f; c.bzf = bzf.f;
        c.evaluate_ramp = [](double x, int ri) { return ri < 1 ? 1.0 : (x * x + 0.75 * ri) / (1.0 + x * x); };
        vent.pressure_ramp_index = 1; vent.dynamic_pressure = 7.0; vent.ior = 3;
    }
};

void test_c2 ()
{
    std::printf("C2 TUNNEL_PRECONDITIONER is refused\n");
    const int N = 4;
    Arr2 b(N, N); Arr3 h(0, N + 1, 0, N + 1, 0, N + 1); h.fill_pattern(1.0);
    HFillOptions opt; opt.domain = fds_interior(N, N, N);
    const std::vector<double> before = h.v;
    bool ok = true; try { pres_h_bc_x(opt.domain, opt, 3, 0.5, b.f, b.f, h.f); pres_h_bc_y(opt.domain, opt, 3, 0.5, b.f, b.f, h.f); pres_h_bc_z(opt.domain, opt, 3, 0.5, b.f, b.f, h.f); } catch (...) { ok = false; }
    check(ok, "tunnel flag off: the H fill runs");
    h.v = before;
    opt.tunnel_preconditioner = true;
    check(throws_as<NotBuilt>([&] { pres_h_bc_x(opt.domain, opt, 3, 0.5, b.f, b.f, h.f); return 0; }), "x fill refuses (NotBuilt)");
    check(throws_as<NotBuilt>([&] { pres_h_bc_y(opt.domain, opt, 3, 0.5, b.f, b.f, h.f); return 0; }), "y fill refuses (NotBuilt)");
    check(throws_as<NotBuilt>([&] { pres_h_bc_z(opt.domain, opt, 0, 0.5, b.f, b.f, h.f); return 0; }), "z fill refuses (NotBuilt) also for a periodic code");
    check(h.v == before, "refused fills wrote nothing");
    try { pres_h_bc_x(opt.domain, opt, 3, 0.5, b.f, b.f, h.f); } catch (NotBuilt const& e) { check(std::strstr(e.what(), "TUNNEL_PRECONDITIONER") != nullptr, "message names TUNNEL_PRECONDITIONER"); }
    Ctx x(3, 3, 3);
    PoissonWall w; w.i = 1; w.j = 1; w.k = 1; w.ior = 1; w.pressure_bc_type = fdsconst::NEUMANN; w.boundary_type = fdsconst::SOLID_BOUNDARY;
    bool ok2 = true; try { pres_poisson_boundary_arrays(x.c, &w, 1); } catch (...) { ok2 = false; }
    check(ok2, "tunnel flag off: L1209 runs");
    x.c.tunnel_preconditioner = true;
    check(throws_as<NotBuilt>([&] { pres_poisson_boundary_arrays(x.c, &w, 1); return 0; }), "L1209 refuses (NotBuilt)");
}

// ---------------------------------------------------------------------------------------------------------------- C3
// One open z wall with the wind option on, at KK = 0 (IOR=3) or KK = KBP1 (IOR=-3): the value must be the one computed by hand from
// W_WIND(KK) (distinct values at every index), so a read one cell off is a different number.
void test_c3 ()
{
    std::printf("C3 L1209 wind guard and the ghost indices KK=0 and KK=KBP1\n");
    const int I = 3, J = 3, K = 4, KBP1 = K + 1;
    for (int zdir = 0; zdir < 2; ++zdir) {
        Ctx x(I, J, K);
        x.c.open_wind_boundary = true;
        PoissonWall w; w.i = 2; w.j = 3; w.ior = zdir == 0 ? 3 : -3; w.k = zdir == 0 ? 0 : KBP1;
        w.pressure_bc_type = fdsconst::DIRICHLET; w.boundary_type = fdsconst::OPEN_BOUNDARY; w.t_ign = 0.5; w.rho_f = 1.2; w.vent = &x.vent;
        // make the outflow test fail on purpose so that the wind H0 branch is the one used: WW(I,J,0) >= 0 for IOR 3, WW(I,J,KBAR) <= 0 for IOR -3
        x.ww.f(w.i, w.j, 0) = 0.25; x.ww.f(w.i, w.j, K) = -0.25;
        pres_poisson_boundary_arrays(x.c, &w, 1);
        const double ts = x.c.t - x.c.t_begin, rf = x.c.evaluate_ramp(ts, 1), pext = rf * 7.0;
        double want;
        const double eddy = 0.0;
        if (zdir == 0) {
            const double h0 = x.hp.f(w.i, w.j, 1) + 0.5 / (x.c.dt * x.rdzn.f(0)) * (x.ww_.f(w.k) + eddy - x.ww.f(w.i, w.j, 0));
            want = pext / w.rho_f + h0;   // WW(I,J,0) = 0.25 is not < 0
            check(same_bits(x.bzs.f(w.i, w.j), want), "IOR=3 wall at KK=0: BZS from W_WIND(0)");
            // a different W_WIND(1) must not matter, a different W_WIND(0) must
            const double keep = x.ww_.f(1); x.ww_.f(1) += 5.0; x.bzs.v.assign(x.bzs.v.size(), 0.0); pres_poisson_boundary_arrays(x.c, &w, 1);
            check(same_bits(x.bzs.f(w.i, w.j), want), "IOR=3: W_WIND(1) is not read");
            x.ww_.f(1) = keep; x.ww_.f(0) += 5.0; pres_poisson_boundary_arrays(x.c, &w, 1);
            check(!same_bits(x.bzs.f(w.i, w.j), want), "IOR=3: W_WIND(0) is read");
        } else {
            const double h0 = x.hp.f(w.i, w.j, K) - 0.5 / (x.c.dt * x.rdzn.f(K)) * (x.ww_.f(w.k) + eddy - x.ww.f(w.i, w.j, K));
            want = pext / w.rho_f + h0;   // WW(I,J,KBAR) = -0.25 is not > 0
            check(same_bits(x.bzf.f(w.i, w.j), want), "IOR=-3 wall at KK=KBP1: BZF from W_WIND(KBP1)");
            const double keep = x.ww_.f(K); x.ww_.f(K) += 5.0; x.bzf.v.assign(x.bzf.v.size(), 0.0); pres_poisson_boundary_arrays(x.c, &w, 1);
            check(same_bits(x.bzf.f(w.i, w.j), want), "IOR=-3: W_WIND(KBAR) is not read");
            x.ww_.f(K) = keep; x.ww_.f(KBP1) += 5.0; pres_poisson_boundary_arrays(x.c, &w, 1);
            check(!same_bits(x.bzf.f(w.i, w.j), want), "IOR=-3: W_WIND(KBP1) is read");
        }
    }
    // wind views that do not cover 0..KBP1 are refused, for each of the three arrays, whether or not the walls would read the missing index
    struct V { int lo, hi; const char* what; };
    const V views[] = {{1, K, "lower bound 1, extent KBAR (the former guard accepted it)"}, {1, KBP1, "lower bound 1, up to KBP1"}, {0, K, "0..KBAR (KBP1 missing)"}, {0, KBP1 - 2, "0..KBAR-1"}};
    for (V const& v : views) {
        for (int which = 0; which < 3; ++which) {
            // build with a full-size wind, then replace one view by the short one
            Ctx x(I, J, K);
            Arr1 shortw(v.lo, v.hi, 0.5);
            x.c.open_wind_boundary = true;
            (which == 0 ? x.c.u_wind : which == 1 ? x.c.v_wind : x.c.w_wind) = shortw.f;
            PoissonWall w; w.i = 2; w.j = 3; w.k = KBP1; w.ior = -3; w.pressure_bc_type = fdsconst::DIRICHLET; w.boundary_type = fdsconst::OPEN_BOUNDARY; w.vent = &x.vent;
            char buf[160]; std::snprintf(buf, sizeof buf, "%s_wind view %s throws", which == 0 ? "u" : which == 1 ? "v" : "w", v.what);
            check(throws_as<std::invalid_argument>([&] { pres_poisson_boundary_arrays(x.c, &w, 1); return 0; }), buf);
        }
    }
    // the same short wind is fine when the wind option is off (the arrays are not read)
    {
        Ctx x(I, J, K); Arr1 shortw(1, K, 0.5); x.c.u_wind = x.c.v_wind = x.c.w_wind = shortw.f;
        PoissonWall w; w.i = 2; w.j = 3; w.k = KBP1; w.ior = -3; w.pressure_bc_type = fdsconst::DIRICHLET; w.boundary_type = fdsconst::OPEN_BOUNDARY; w.vent = &x.vent;
        bool ok = true; try { pres_poisson_boundary_arrays(x.c, &w, 1); } catch (...) { ok = false; }
        check(ok, "wind option off: short wind views are not needed");
    }
}
}   // namespace

int main ()
{
    test_c1();
    test_c2();
    test_c3();
    std::printf("SUMMARY checks=%d failed=%d %s\n", g_n, g_fail, g_fail == 0 ? "PASS" : "FAIL");
    return g_fail == 0 ? 0 : 1;
}
