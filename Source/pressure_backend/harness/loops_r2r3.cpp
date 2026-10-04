// loops_r2r3.cpp: pb_fds_loops_r2r3, self-contained tests for the two call-site conditions of the sign-off that need no driver:
//   R2  PRHS from AMReX low-face-indexed flux arrays equals the FDS-index version, bitwise, on a box that does not start at 0
//       (pres_compute_rhs_div_faces / fds_face_view); negative controls: the nodal array wrapped "as is" differs, and mutant 22
//       (wrong shift direction) must fail this program.
//   R3  D-067 ordering on a singular all-Neumann component with non-uniform rho: gauge, then the H ghost fill (L1220-L1222), then
//       P = RHOP*(HP-KRES) (L1207) over the full box. Then sum_V P = 0 over the gas cells, and the ghost P equals RHOP_ghost*(H_ghost-KRES_ghost);
//       with the fill done BEFORE the gauge the ghost layer is off by exactly the gauge constant (P off by c*rho_ghost, summed c*sum rho).
// Acceptance of the AMR route against FDS stays on eps_H, not bitwise (the translations are bitwise only against FDS-compiled loops).
// No AMReX, no MPI.
#include "FdsPressureLoops.H"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

using namespace pressure_backend::fdsloops;

namespace {
int g_fail = 0, g_n = 0;
void check (bool ok, const char* what) { ++g_n; if (!ok) ++g_fail; std::printf("  %s: %s\n", ok ? "PASS" : "FAIL", what); }
bool same_bits (double a, double b) { return std::memcmp(&a, &b, 8) == 0; }

struct Arr3 {
    std::vector<double> v; F3 f;
    Arr3 (int l0, int h0, int l1, int h1, int l2, int h2) {
        f.lo[0] = l0; f.lo[1] = l1; f.lo[2] = l2; f.n[0] = h0 - l0 + 1; f.n[1] = h1 - l1 + 1; f.n[2] = h2 - l2 + 1;
        v.assign(static_cast<std::size_t>(f.n[0]) * f.n[1] * f.n[2], 0.0); f.p = v.data();
    }
    Arr3 (Arr3 const& o) : v(o.v), f(o.f) { f.p = v.data(); }   // a copy owns its data (the view must not alias the original)
    Arr3& operator= (Arr3 const&) = delete;
};
struct Arr1 { std::vector<double> v; F1 f; Arr1 (int lo, int hi) { v.assign(hi - lo + 1, 0.0); f.p = v.data(); f.lo = lo; f.n = hi - lo + 1; } };
struct Arr2 { std::vector<double> v; F2 f; Arr2 (int n0, int n1) { v.assign(static_cast<std::size_t>(n0) * n1, 0.0); f.p = v.data(); f.lo[0] = 1; f.lo[1] = 1; f.n[0] = n0; f.n[1] = n1; } };

struct Rng { unsigned long long s = 88172645463325252ULL; double next () { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return static_cast<double>(s >> 11) / 9007199254740992.0 - 0.5; } };

// ---------------------------------------------------------------------------------------------------------------- R2
void test_r2 ()
{
    std::printf("R2 PRHS from AMReX low-face-indexed arrays equals the FDS-index version (box not starting at 0)\n");
    struct B { Box3 b; const char* what; };
    const B boxes[] = {{{{1, 1, 1}, {5, 4, 6}}, "FDS box 1..N"}, {{{3, 5, -2}, {7, 9, 2}}, "box starting at (3,5,-2)"}, {{{-4, 0, 10}, {-1, 0, 13}}, "box with a one-cell y range at y=0 and negative x"}};
    Rng rng;
    for (B const& bx : boxes) {
        const Box3 b = bx.b;
        const int il = b.lo[0], ih = b.hi[0], jl = b.lo[1], jh = b.hi[1], kl = b.lo[2], kh = b.hi[2];
        // FDS-index data: FVX(il-1:ih,..) etc. over the box (+ the low neighbour), with a larger allocation so that bounds differ from the box
        Arr3 fvx(il - 2, ih + 1, jl - 1, jh + 1, kl - 1, kh + 1), fvy(il - 1, ih + 1, jl - 2, jh + 1, kl - 1, kh + 1), fvz(il - 1, ih + 1, jl - 1, jh + 1, kl - 2, kh + 1);
        Arr3 dddt(il - 1, ih + 1, jl - 1, jh + 1, kl - 1, kh + 1);
        for (auto* a : {&fvx, &fvy, &fvz, &dddt}) for (double& x : a->v) x = rng.next();
        Arr1 rdx(il - 1, ih + 1), rdy(jl - 1, jh + 1), rdz(kl - 1, kh + 1);
        for (auto* a : {&rdx, &rdy, &rdz}) for (double& x : a->v) x = 10.0 + 4.0 * rng.next();
        // reference: FDS-index call
        Arr3 prhs_fds(il, ih, jl, jh, kl, kh), prhs_amr(il, ih, jl, jh, kl, kh), prhs_nodal(il, ih, jl, jh, kl, kh);
        pres_compute_rhs_div(b, fvx.f, fvy.f, fvz.f, dddt.f, rdx.f, rdy.f, rdz.f, prhs_fds.f);
        // AMReX layout: nodal arrays over the faces; x-face i (low face of cell i) = FDS FVX(i-1). Own bounds, shifted by one against the FDS arrays.
        Arr3 ax(il - 1, ih + 2, jl - 1, jh + 1, kl - 1, kh + 1), ay(il - 1, ih + 1, jl - 1, jh + 2, kl - 1, kh + 1), az(il - 1, ih + 1, jl - 1, jh + 1, kl - 1, kh + 2);
        for (int k = kl - 1; k <= kh + 1; ++k) for (int j = jl - 1; j <= jh + 1; ++j) for (int i = il - 1; i <= ih + 2; ++i) ax.f(i, j, k) = (i - 1 >= il - 2 && i - 1 <= ih + 1) ? fvx.f(i - 1, j, k) : 1e300;
        for (int k = kl - 1; k <= kh + 1; ++k) for (int j = jl - 1; j <= jh + 2; ++j) for (int i = il - 1; i <= ih + 1; ++i) ay.f(i, j, k) = (j - 1 >= jl - 2 && j - 1 <= jh + 1) ? fvy.f(i, j - 1, k) : 1e300;
        for (int k = kl - 1; k <= kh + 2; ++k) for (int j = jl - 1; j <= jh + 1; ++j) for (int i = il - 1; i <= ih + 1; ++i) az.f(i, j, k) = (k - 1 >= kl - 2 && k - 1 <= kh + 1) ? fvz.f(i, j, k - 1) : 1e300;
        pres_compute_rhs_div_faces(b, ax.f, ay.f, az.f, dddt.f, rdx.f, rdy.f, rdz.f, prhs_amr.f);
        std::string w1 = std::string(bx.what) + ": face-indexed PRHS is bitwise the FDS-index PRHS";
        check(prhs_amr.v.size() == prhs_fds.v.size() && std::memcmp(prhs_amr.v.data(), prhs_fds.v.data(), 8 * prhs_fds.v.size()) == 0, w1.c_str());
        // the face arrays hold nothing but the right values: perturb the one face the box does not use, the result must not move
        ax.f(il - 1, jl - 1, kl - 1) += 1.0;
        Arr3 again(il, ih, jl, jh, kl, kh);
        pres_compute_rhs_div_faces(b, ax.f, ay.f, az.f, dddt.f, rdx.f, rdy.f, rdz.f, again.v.empty() ? prhs_amr.f : again.f);
        check(std::memcmp(again.v.data(), prhs_fds.v.data(), 8 * prhs_fds.v.size()) == 0, "a face outside the box stencil does not enter");
        // negative control: the nodal arrays wrapped as they are (no index shift) give a different PRHS; throws if the view does not cover
        bool differs = false, threw = false;
        try { pres_compute_rhs_div(b, ax.f, ay.f, az.f, dddt.f, rdx.f, rdy.f, rdz.f, prhs_nodal.f); differs = std::memcmp(prhs_nodal.v.data(), prhs_fds.v.data(), 8 * prhs_fds.v.size()) != 0; }
        catch (std::invalid_argument const&) { threw = true; }
        check(differs || threw, "negative control: the nodal arrays wrapped without the shift differ from the FDS-index result (or are refused)");
        // a face array that lacks the high face ih+1 is refused
        Arr3 short_x(il - 1, ih, jl - 1, jh + 1, kl - 1, kh + 1);
        bool refused = false;
        try { pres_compute_rhs_div_faces(b, short_x.f, ay.f, az.f, dddt.f, rdx.f, rdy.f, rdz.f, prhs_amr.f); } catch (std::invalid_argument const&) { refused = true; }
        check(refused, "x-face array without the face ih+1 is refused");
    }
}

// ---------------------------------------------------------------------------------------------------------------- R3
void test_r3 ()
{
    std::printf("R3 D-067 order: gauge, H ghost fill, then P (singular all-Neumann component, non-uniform rho)\n");
    const int nx = 5, ny = 4, nz = 6;
    Rng rng;
    Arr3 H(0, nx + 1, 0, ny + 1, 0, nz + 1), rho(-1, nx + 2, -1, ny + 2, -1, nz + 2), kres(0, nx + 1, 0, ny + 1, 0, nz + 1);
    Arr1 dx(1, nx), dy(1, ny), dz(1, nz);
    for (int q = 1; q <= nx; ++q) dx.f(q) = 0.1 + 0.03 * q; for (int q = 1; q <= ny; ++q) dy.f(q) = 0.2 - 0.02 * q; for (int q = 1; q <= nz; ++q) dz.f(q) = 0.15 + 0.01 * q * q;
    for (int k = -1; k <= nz + 2; ++k) for (int j = -1; j <= ny + 2; ++j) for (int i = -1; i <= nx + 2; ++i) rho.f(i, j, k) = 1.0 + 0.4 * (rng.next() + 0.5) + 0.3 * std::sin(0.7 * i + 0.4 * j - 0.2 * k);   // stratified-like, non-uniform
    for (int k = 0; k <= nz + 1; ++k) for (int j = 0; j <= ny + 1; ++j) for (int i = 0; i <= nx + 1; ++i) kres.f(i, j, k) = 0.5 * (1.0 + rng.next()) ;
    // a solution with an arbitrary additive constant (what a Poisson solver returns for a singular problem)
    for (int k = 1; k <= nz; ++k) for (int j = 1; j <= ny; ++j) for (int i = 1; i <= nx; ++i) H.f(i, j, k) = 37.5 + 0.8 * rng.next() + 0.1 * i;
    Arr2 bxs(ny, nz), bxf(ny, nz), bys(nx, nz), byf(nx, nz), bzs(nx, ny), bzf(nx, ny);
    for (auto* a : {&bxs, &bxf, &bys, &byf, &bzs, &bzf}) for (double& x : a->v) x = 0.3 * rng.next();   // Neumann data (H_n), all faces
    const Box3 box = fds_interior(nx, ny, nz);
    HFillOptions opt; opt.domain = box;
    const double dxi = 0.07, deta = 0.05, dzeta = 0.09;
    auto fill = [&](Arr3& h) {
        pres_h_bc_x(box, opt, 3, dxi, bxs.f, bxf.f, h.f); pres_h_bc_y(box, opt, 3, deta, bys.f, byf.f, h.f); pres_h_bc_z(box, opt, 3, dzeta, bzs.f, bzf.f, h.f);
    };
    // gauge of D-067: sum over gas cells of V*rho*(H-KRES) = 0 (V = dx dy dz, rho = RHOP, KRES = the FDS kinetic-energy field)
    auto gauge = [&](Arr3& h) {
        long double num = 0, den = 0;
        for (int k = 1; k <= nz; ++k) for (int j = 1; j <= ny; ++j) for (int i = 1; i <= nx; ++i) {
            const long double w = static_cast<long double>(dx.f(i)) * dy.f(j) * dz.f(k) * rho.f(i, j, k);
            num += w * (static_cast<long double>(h.f(i, j, k)) - kres.f(i, j, k)); den += w;
        }
        const double c = static_cast<double>(num / den);
        for (int k = 1; k <= nz; ++k) for (int j = 1; j <= ny; ++j) for (int i = 1; i <= nx; ++i) h.f(i, j, k) -= c;
        return c;
    };
    auto sumvp = [&](Arr3 const& p, long double& scale) {
        long double s = 0; scale = 0;
        for (int k = 1; k <= nz; ++k) for (int j = 1; j <= ny; ++j) for (int i = 1; i <= nx; ++i) {
            const long double w = static_cast<long double>(dx.f(i)) * dy.f(j) * dz.f(k) * rho.f(i, j, k); s += w * p.f(i, j, k) / rho.f(i, j, k); scale += std::fabs(w * p.f(i, j, k) / rho.f(i, j, k));
        }
        return s;
    };
    // sum_V P with P = rho*(H-KRES): weights V, not V*rho, because P already carries rho
    auto sumP = [&](Arr3 const& p) { long double s = 0, sc = 0; for (int k = 1; k <= nz; ++k) for (int j = 1; j <= ny; ++j) for (int i = 1; i <= nx; ++i) { const long double v = static_cast<long double>(dx.f(i)) * dy.f(j) * dz.f(k); s += v * p.f(i, j, k); sc += v * std::fabs(p.f(i, j, k)); } (void)sumvp; return std::make_pair(s, sc); };

    // correct order
    Arr3 Hc = H; const double c = gauge(Hc); fill(Hc);
    Arr3 Pc(0, nx + 1, 0, ny + 1, 0, nz + 1);
    pres_p_from_h(fds_with_ghosts(nx, ny, nz), rho.f, Hc.f, kres.f, Pc.f);
    // note: sum_V P here is sum V*rho*(H-KRES), which the gauge sets to zero
    auto sp = sumP(Pc);
    char buf[200];
    std::snprintf(buf, sizeof buf, "gauge constant c = %.6g was removed (large offset exercised)", c);
    check(std::fabs(c) > 30.0, buf);
    std::snprintf(buf, sizeof buf, "gauge, fill, P: sum_V P = %.3e, relative to sum_V|P| = %.2e", static_cast<double>(sp.first), static_cast<double>(sp.first / sp.second));
    check(std::fabs(static_cast<double>(sp.first / sp.second)) < 1e-14, buf);
    // ghost P = RHOP_ghost*(H_ghost-KRES_ghost), where H_ghost are the Neumann ghosts of the gauged interior
    bool ghost_ok = true; long ng = 0;
    for (int k = 0; k <= nz + 1; ++k) for (int j = 0; j <= ny + 1; ++j) for (int i = 0; i <= nx + 1; ++i) {
        const bool interior = i >= 1 && i <= nx && j >= 1 && j <= ny && k >= 1 && k <= nz;
        if (interior) continue;
        ++ng; ghost_ok = ghost_ok && same_bits(Pc.f(i, j, k), rho.f(i, j, k) * (Hc.f(i, j, k) - kres.f(i, j, k)));
    }
    std::snprintf(buf, sizeof buf, "ghost P (%ld ghost cells, edges and corners included) = rho_ghost*(H_ghost-KRES_ghost), bitwise", ng);
    check(ghost_ok, buf);
    // the face-normal ghost layers carry the gauged interior value: H(0,J,K) = H(1,J,K) - DXI*BXS(J,K) (all-Neumann), the same for the five other faces
    bool nm_ok = true;
    for (int k = 1; k <= nz; ++k) for (int j = 1; j <= ny; ++j) nm_ok = nm_ok && same_bits(Hc.f(0, j, k), Hc.f(1, j, k) - dxi * bxs.f(j, k)) && same_bits(Hc.f(nx + 1, j, k), Hc.f(nx, j, k) + dxi * bxf.f(j, k));
    for (int k = 1; k <= nz; ++k) for (int i = 1; i <= nx; ++i) nm_ok = nm_ok && same_bits(Hc.f(i, 0, k), Hc.f(i, 1, k) - deta * bys.f(i, k)) && same_bits(Hc.f(i, ny + 1, k), Hc.f(i, ny, k) + deta * byf.f(i, k));
    for (int j = 1; j <= ny; ++j) for (int i = 1; i <= nx; ++i) nm_ok = nm_ok && same_bits(Hc.f(i, j, 0), Hc.f(i, j, 1) - dzeta * bzs.f(i, j)) && same_bits(Hc.f(i, j, nz + 1), Hc.f(i, j, nz) + dzeta * bzf.f(i, j));
    check(nm_ok, "face ghost layers are the Neumann extrapolation of the gauged interior");

    // negative control: fill BEFORE the gauge. The interior is shifted by c afterwards, the ghost layer keeps the old offset.
    Arr3 Hw = H; fill(Hw); gauge(Hw);
    Arr3 Pw(0, nx + 1, 0, ny + 1, 0, nz + 1);
    pres_p_from_h(fds_with_ghosts(nx, ny, nz), rho.f, Hw.f, kres.f, Pw.f);
    long double dsum = 0, rsum = 0; double maxrel = 0; bool interior_same = true;
    for (int k = 0; k <= nz + 1; ++k) for (int j = 0; j <= ny + 1; ++j) for (int i = 0; i <= nx + 1; ++i) {
        const bool interior = i >= 1 && i <= nx && j >= 1 && j <= ny && k >= 1 && k <= nz;
        if (interior) { interior_same = interior_same && same_bits(Pw.f(i, j, k), Pc.f(i, j, k)); continue; }
        const bool face_ghost = (i >= 1 && i <= nx ? 0 : 1) + (j >= 1 && j <= ny ? 0 : 1) + (k >= 1 && k <= nz ? 0 : 1) == 1;
        if (!face_ghost) continue;   // edges and corners are filled by later statements from earlier ghosts; the face layers carry the control
        const double d = Pw.f(i, j, k) - Pc.f(i, j, k);
        dsum += d; rsum += rho.f(i, j, k);
        // the fill-before-gauge ghost of a Neumann face is H_old(1)+..., it keeps the old offset: H_ghost differs by exactly c
        maxrel = std::max(maxrel, std::fabs(d - c * rho.f(i, j, k)) / (std::fabs(c) * rho.f(i, j, k)));
    }
    check(interior_same, "fill-before-gauge: interior P unchanged (the gauge acts on the interior either way)");
    std::snprintf(buf, sizeof buf, "fill-before-gauge fails: face-ghost P off by c*rho_ghost (max rel deviation from that %.1e), summed %.6g = c*sum rho (%.6g)", maxrel, static_cast<double>(dsum), static_cast<double>(c * rsum));
    check(maxrel < 1e-12 && std::fabs(static_cast<double>(dsum - c * rsum)) < 1e-9 * std::fabs(static_cast<double>(c * rsum)) && std::fabs(static_cast<double>(dsum)) > 1.0, buf);
}
}   // namespace

int main ()
{
    try { test_r2(); test_r3(); }
    catch (std::exception const& e) { std::printf("  FAIL: exception %s\n", e.what()); ++g_fail; ++g_n; }
    std::printf("SUMMARY checks=%d failed=%d %s\n", g_n, g_fail, g_fail == 0 ? "PASS" : "FAIL");
    return g_fail == 0 ? 0 : 1;
}
