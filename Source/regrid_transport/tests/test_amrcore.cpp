// test_amrcore.cpp: R2a checks with AMReX (no FDS Fortran, no GPU): the AmrCore subclass builds the static two-level hierarchy of an
// ns2d_16_int_1to2_refinement-style input and notifies the driver's LevelRegistry; average-down conserves the integral to round-off; a uniform state
// stays uniform across the coarse-fine face after the ghost fill; ghost layers 1 and 2 hold the coarse value (also with ratio 4, periodic wrap and a
// neighbouring fine box). Usage: test_amrcore <cases dir>. Run on 1 rank (and 4 ranks by hand).
#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>

#include <cmath>
#include <array>
#include <cstdio>
#include <tuple>
#include <string>

#include "DriverAdapter.H"
#include "FluxOverrideOps.H"
#include "LevelOps.H"
#include "RegridAmrCore.H"
#include "check.H"
#include "mesh_text.H"

using namespace fdsrt;

// Test-only mirror of fdsamr::exchange_fields (Source/driver/GhostExchange.cpp, which also holds the FDS-linked BcStep and cannot be linked without the Fortran
// objects). If the driver list changes, this copy must follow; the sentinel checks below do not depend on the exact list for H/HS and RHO/ZZ.
namespace fdsamr {
std::vector<std::string> exchange_fields(int code, bool predictor)
{
    switch (code) {
    case 1: return {"RHOS", "ZZS", "MU", "KRES", "D", "TMP", "RSUM"};
    case 4: return {"RHO", "ZZ", "MU", "KRES", "DS", "TMP", "RSUM"};
    case 3: return {"US", "VS", "WS", "HS"};
    case 6: return {"U", "V", "W", "H"};
    case 5: return {"FVX", "FVY", "FVZ", predictor ? "H" : "HS"};
    default: return {};
    }
}
}  // namespace fdsamr

namespace {

double fval(int i, int j, int k, int n) { return 1.0 + 0.001 * i + 0.0173 * j + 0.31 * k + 7.0 * n + std::sin(0.37 * i + 0.11 * k); }
double gval(int i, int j, int k, int n) { return 100.0 + 0.5 * i - 0.25 * j + 0.125 * k + 3.0 * n; }   // coarse values: distinct per cell

void fill_all(amrex::MultiFab& mf, double (*f)(int, int, int, int), bool with_ghost)
{
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        auto a = mf.array(mfi);
        const amrex::Box b = with_ghost ? mfi.growntilebox() : mfi.validbox();
        for (int n = 0; n < mf.nComp(); ++n) amrex::LoopOnCpu(b, [&](int i, int j, int k) { a(i, j, k, n) = f(i, j, k, n); });
    }
}

// Direct test of the operators (no registry): coarse domain 16^3 periodic, fine boxes with given ratio.
void test_ops(int rr, bool two_boxes)
{
    const std::string tag = "ops ratio " + std::to_string(rr) + (two_boxes ? " two boxes" : " one box");
    amrex::Box cdom(amrex::IntVect(0), amrex::IntVect(15));
    amrex::RealBox rb({0, 0, 0}, {1, 1, 1});
    amrex::Array<int, 3> per{1, 1, 1};
    amrex::Geometry cg(cdom, rb, 0, per);
    amrex::IntVect ratio(rr);
    amrex::Geometry fg = amrex::refine(cg, ratio);
    amrex::BoxArray cba(cdom);
    cba.maxSize(8);
    amrex::DistributionMapping cdm(cba);
    // fine region = coarse cells (0..7, 4..11, 4..11): touches the low x edge, so the x ghost layers wrap through the periodic boundary
    amrex::BoxList bl;
    amrex::Box fregion = amrex::refine(amrex::Box(amrex::IntVect(0, 4, 4), amrex::IntVect(7, 11, 11)), ratio);
    if (two_boxes) {
        amrex::Box a = fregion, b = fregion;
        const int mid = fregion.smallEnd(0) + fregion.length(0) / 2 - 1;
        a.setBig(0, mid);
        b.setSmall(0, mid + 1);
        bl.push_back(a);
        bl.push_back(b);
    } else {
        bl.push_back(fregion);
    }
    amrex::BoxArray fba(bl);
    amrex::DistributionMapping fdm(fba);
    const int ncomp = 2, ng = 3;
    amrex::MultiFab C(cba, cdm, ncomp, 0), F(fba, fdm, ncomp, ng);
    fill_all(C, gval, false);
    F.setVal(-777.0);
    fill_all(F, fval, false);
    F.FillBoundary(fg.periodicity());   // same-level fill first (as the driver does)
    amrex::MultiFab F_before(fba, fdm, ncomp, ng);
    amrex::MultiFab::Copy(F_before, F, 0, 0, ncomp, ng);
    const long nf = fill_cf_ghosts_pc(F, C, fg, cg, ratio, 2, 0, ncomp);
    long bad_c = 0, bad_same = 0, bad_beyond = 0, n2 = 0;
    for (amrex::MFIter mfi(F); mfi.isValid(); ++mfi) {
        auto a = F.const_array(mfi);
        auto b0 = F_before.const_array(mfi);
        const amrex::Box vb = mfi.validbox();
        amrex::LoopOnCpu(mfi.fabbox(), [&](int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            if (vb.contains(iv)) return;
            const int layer = std::max({vb.smallEnd(0) - i, i - vb.bigEnd(0), vb.smallEnd(1) - j, j - vb.bigEnd(1), vb.smallEnd(2) - k, k - vb.bigEnd(2)});
            amrex::IntVect w = iv;   // wrapped position
            for (int d = 0; d < 3; ++d) { if (w[d] < 0) w[d] += fg.Domain().length(d); if (w[d] >= fg.Domain().length(d)) w[d] -= fg.Domain().length(d); }
            const bool in_fine = fba.contains(w);
            for (int n = 0; n < ncomp; ++n) {
                if (in_fine) { if (a(i, j, k, n) != b0(i, j, k, n)) ++bad_same; }
                else if (layer <= 2) {
                    const amrex::IntVect c = amrex::coarsen(iv, ratio);
                    amrex::IntVect cw = c;
                    for (int d = 0; d < 3; ++d) { if (cw[d] < 0) cw[d] += 16; if (cw[d] >= 16) cw[d] -= 16; }
                    if (a(i, j, k, n) != gval(cw[0], cw[1], cw[2], n)) ++bad_c;
                    if (n == 0) ++n2;
                } else if (a(i, j, k, n) != b0(i, j, k, n)) ++bad_beyond;
            }
        });
    }
    // layer 2 equals layer 1: the cell one step further out in a face-normal direction sits in the same coarse cell
    long bad_l12 = 0;
    for (amrex::MFIter mfi(F); mfi.isValid(); ++mfi) {
        auto a = F.const_array(mfi);
        const amrex::Box vb = mfi.validbox();
        const amrex::IntVect probes[3] = {{vb.smallEnd(0), vb.smallEnd(1) + 1, vb.smallEnd(2) + 1}, {vb.smallEnd(0) + 1, vb.smallEnd(1), vb.smallEnd(2) + 1},
                                          {vb.smallEnd(0) + 1, vb.smallEnd(1) + 1, vb.smallEnd(2)}};
        for (int d = 0; d < 3; ++d) {
            amrex::IntVect p1 = probes[d], p2 = probes[d];
            p1[d] = vb.smallEnd(d) - 1; p2[d] = vb.smallEnd(d) - 2;
            amrex::IntVect w = p1; for (int e = 0; e < 3; ++e) { if (w[e] < 0) w[e] += fg.Domain().length(e); }
            if (fba.contains(w)) continue;   // a neighbouring fine box holds it
            for (int n = 0; n < ncomp; ++n) if (a(p1[0], p1[1], p1[2], n) != a(p2[0], p2[1], p2[2], n)) ++bad_l12;
        }
    }
    amrex::ParallelDescriptor::ReduceLongSum(n2);   // ranks without a fine box have none
    CHECK_MSG(nf >= 0, tag);
    CHECK_MSG(n2 > 0, tag + ": some ghost cells were filled");
    CHECK_MSG(bad_c == 0, tag + ": ghost cells (layers 1,2) hold the value of the coarse cell that contains them, bad=" + std::to_string(bad_c));
    CHECK_MSG(bad_same == 0, tag + ": cells held by a neighbouring fine box are not overwritten");
    CHECK_MSG(bad_beyond == 0, tag + ": layer 3 is untouched");
    CHECK_MSG(bad_l12 == 0, tag + ": layer 2 equals layer 1");

    // average-down: covered coarse cells = mean of the fine cells, the rest untouched, integral conserved
    amrex::MultiFab C2(cba, cdm, ncomp, 0);
    fill_all(C2, gval, false);
    average_down_cells(F, C2, fg, cg, ratio, 0, ncomp);
    double sum_all = 0, sum_unc = 0, sum_fine = 0, maxdev = 0;
    long untouched_bad = 0, covered_n = 0;
    amrex::BoxArray covered = fba;
    covered.coarsen(ratio);
    for (amrex::MFIter mfi(C2); mfi.isValid(); ++mfi) {
        auto a = C2.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            const bool cov = covered.contains(amrex::IntVect(i, j, k));
            for (int n = 0; n < ncomp; ++n) {
                sum_all += a(i, j, k, n);
                if (!cov) { sum_unc += a(i, j, k, n); if (a(i, j, k, n) != gval(i, j, k, n)) ++untouched_bad; }
                else {
                    double m = 0;
                    for (int c = 0; c < rr * rr * rr; ++c) m += fval(i * rr + c % rr, j * rr + (c / rr) % rr, k * rr + c / (rr * rr), n);
                    m /= rr * rr * rr;
                    maxdev = std::max(maxdev, std::abs(a(i, j, k, n) - m) / std::abs(m));
                }
            }
            if (cov) ++covered_n;
        });
    }
    for (amrex::MFIter mfi(F); mfi.isValid(); ++mfi) {
        auto a = F.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { for (int n = 0; n < ncomp; ++n) sum_fine += a(i, j, k, n); });
    }
    amrex::ParallelDescriptor::ReduceRealSum(sum_all);
    amrex::ParallelDescriptor::ReduceRealSum(sum_unc);
    amrex::ParallelDescriptor::ReduceRealSum(sum_fine);
    amrex::ParallelDescriptor::ReduceRealMax(maxdev);
    amrex::ParallelDescriptor::ReduceLongSum(untouched_bad);
    amrex::ParallelDescriptor::ReduceLongSum(covered_n);
    const double vc = 1.0, vf = 1.0 / (rr * rr * rr);   // cell volumes in units of a coarse cell
    const double total_before = sum_unc * vc + sum_fine * vf, total_after = sum_all * vc;
    CHECK_MSG(covered_n == 8 * 8 * 8, tag + ": 512 covered coarse cells, got " + std::to_string(covered_n));
    CHECK_MSG(untouched_bad == 0, tag + ": uncovered coarse cells unchanged");
    CHECK_MSG(maxdev < 1e-14, tag + ": covered cell = mean of its fine cells, max rel dev " + std::to_string(maxdev));
    CHECK_MSG(std::abs(total_after - total_before) < 1e-13 * std::abs(total_before), tag + ": integral conserved, rel diff " + std::to_string(std::abs(total_after - total_before) / total_before));
}


// ---- FDS ghost rules (wall.f90 ASSIGN_GHOST_VALUE, velo.f90 VISCOSITY_BC) against an independent hand computation ----
// Coarse domain 16^3 (not periodic), fine patch = coarse cells 4..11 in each direction (so every face has 8 x 8 covered cells), cut into several fine boxes.
double rv(int salt, int i, int j, int k, int n)
{
    unsigned long long h = 1469598103934665603ULL ^ (unsigned long long)(salt * 7919 + n * 104729);
    for (int v : {i + 100, j + 100, k + 100}) { h ^= (unsigned long long)v; h *= 1099511628211ULL; }
    return static_cast<double>((h >> 11) & 0xFFFFFFFFFFFFFULL) / 4503599627370496.0;   // [0,1)
}
// initial values: salts 1-5 coarse (RHO, ZZ, TMP, RSUM, MU), 11-15 fine
double v_rho(int s, int i, int j, int k) { return 1.0 + rv(s, i, j, k, 0); }
double v_zz(int s, int i, int j, int k, int n) { return rv(s + 1, i, j, k, n); }
double v_tmp(int s, int i, int j, int k) { return 280.0 + 40.0 * rv(s + 2, i, j, k, 0); }
double v_rsum(int s, int i, int j, int k) { return 280.0 + 20.0 * rv(s + 3, i, j, k, 0); }
double v_mu(int s, int i, int j, int k) { return 1e-5 + 1e-5 * rv(s + 4, i, j, k, 0); }
double eos_rsum(const double* zz) { return 287.0 + 10.0 * zz[0] + 3.0 * zz[1]; }
double eos_pbar(const amrex::IntVect& c) { return 101325.0 + 11.0 * c[2] + 0.5 * c[0]; }
double clip01t(double x) { return std::min(1.0, std::max(0.0, x)); }
bool close_to(double a, double b) { return std::abs(a - b) <= 1e-13 * std::max(1e-5, std::abs(b)); }

void test_fds_rules(int rr, bool with_thermo)
{
    const std::string tag = "fds rules ratio " + std::to_string(rr) + (with_thermo ? " with EOS" : " without EOS");
    amrex::Box cdom(amrex::IntVect(0), amrex::IntVect(15));
    amrex::RealBox rb({0, 0, 0}, {1, 1, 1});
    amrex::Array<int, 3> per{0, 0, 0};
    amrex::Geometry cg(cdom, rb, 0, per);
    amrex::IntVect ratio(rr);
    amrex::Geometry fg = amrex::refine(cg, ratio);
    amrex::BoxArray cba(cdom);
    cba.maxSize(8);
    amrex::DistributionMapping cdm(cba);
    amrex::BoxArray fba(amrex::refine(amrex::Box(amrex::IntVect(4), amrex::IntVect(11)), ratio));
    fba.maxSize(4 * rr);   // 2 x 2 x 2 fine boxes: box-to-box faces inside the patch
    amrex::DistributionMapping fdm(fba);
    const int nzz = 2;
    amrex::MultiFab Crho(cba, cdm, 1, 3), Czz(cba, cdm, nzz, 2), Ctmp(cba, cdm, 1, 2), Crsum(cba, cdm, 1, 1), Cmu(cba, cdm, 1, 1);
    amrex::MultiFab Frho(fba, fdm, 1, 3), Fzz(fba, fdm, nzz, 2), Ftmp(fba, fdm, 1, 2), Frsum(fba, fdm, 1, 1), Fmu(fba, fdm, 1, 1);
    auto setf = [&](amrex::MultiFab& m, auto f) {
        m.setVal(-555.0);
        for (amrex::MFIter mfi(m); mfi.isValid(); ++mfi) {
            auto a = m.array(mfi);
            for (int n = 0; n < m.nComp(); ++n) amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { a(i, j, k, n) = f(i, j, k, n); });
        }
    };
    for (int lev = 0; lev < 2; ++lev) {
        const int s = lev == 0 ? 1 : 11;
        setf(lev == 0 ? Crho : Frho, [&](int i, int j, int k, int) { return v_rho(s, i, j, k); });
        setf(lev == 0 ? Czz : Fzz, [&](int i, int j, int k, int n) { return v_zz(s, i, j, k, n); });
        setf(lev == 0 ? Ctmp : Ftmp, [&](int i, int j, int k, int) { return v_tmp(s, i, j, k); });
        setf(lev == 0 ? Crsum : Frsum, [&](int i, int j, int k, int) { return v_rsum(s, i, j, k); });
        setf(lev == 0 ? Cmu : Fmu, [&](int i, int j, int k, int) { return v_mu(s, i, j, k); });
    }
    ScalarStage cs, fs;
    cs.rho = &Crho; cs.zz = &Czz; cs.tmp = &Ctmp; cs.rsum = &Crsum; cs.mean = {&Cmu};
    fs.rho = &Frho; fs.zz = &Fzz; fs.tmp = &Ftmp; fs.rsum = &Frsum; fs.mean = {&Fmu};
    Thermo th;
    th.n_tracked = 2; th.rsum = eos_rsum; th.pbar = eos_pbar;   // the same function is used on both levels, indexed by that level's cell
    const CfStats sc = fill_covered_ghosts_fds(cs, fs, fg, cg, ratio, with_thermo ? &th : nullptr);
    const CfStats sf = fill_fine_ghosts_fds(fs, cs, fg, cg, ratio, with_thermo ? &th : nullptr);

    // ---- coarse covered cells ----
    long n_checked = 0, bad = 0, bad_untouched = 0, conflicts_expected = 0;
    const double n_int = static_cast<double>(rr * rr);
    for (amrex::MFIter mfi(Crho); mfi.isValid(); ++mfi) {
        auto crho = Crho.const_array(mfi); auto czz = Czz.const_array(mfi); auto ctmp = Ctmp.const_array(mfi); auto crsum = Crsum.const_array(mfi); auto cmu = Cmu.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            const int p[3] = {i - 4, j - 4, k - 4};
            if (p[0] < 0 || p[0] > 7 || p[1] < 0 || p[1] > 7 || p[2] < 0 || p[2] > 7) {   // uncovered coarse cell: never written
                if (crho(i, j, k) != v_rho(1, i, j, k)) ++bad_untouched;
                return;
            }
            int ncand = 0, dsel = -1, ssel = -1, lay = 0;
            for (int layer = 1; layer <= 2; ++layer)
                for (int d = 0; d < 3; ++d)
                    for (int s = 0; s < 2; ++s)
                        if (p[d] == (s == 0 ? layer - 1 : 8 - layer)) { if (ncand++ == 0) { dsel = d; ssel = s; lay = layer; } }
            if (ncand > 1) ++conflicts_expected;
            if (ncand == 0) {   // deep covered cell: untouched by this function
                if (crho(i, j, k) != v_rho(1, i, j, k) || cmu(i, j, k) != v_mu(1, i, j, k)) ++bad_untouched;
                return;
            }
            if (ncand > 1 && lay == 1 && false) return;
            const int t1 = (dsel + 1) % 3, t2 = (dsel + 2) % 3;
            int l1[3] = {i, j, k};   // the layer-1 cell (the cell itself, or its neighbour towards the face for layer 2)
            if (lay == 2) l1[dsel] += (ssel == 0) ? -1 : 1;
            double s_rho = 0, s_rz[2] = {0, 0}, s_tmp = 0, s_rsum = 0, s_mu = 0;
            for (int a = 0; a < rr; ++a)
                for (int b = 0; b < rr; ++b) {
                    int f[3];
                    f[dsel] = (ssel == 0) ? l1[dsel] * rr : l1[dsel] * rr + rr - 1;   // the fine layer next to the face only, not all rr^3 children
                    f[t1] = l1[t1] * rr + a; f[t2] = l1[t2] * rr + b;
                    s_rho += v_rho(11, f[0], f[1], f[2]);
                    for (int n = 0; n < 2; ++n) s_rz[n] += v_rho(11, f[0], f[1], f[2]) * v_zz(11, f[0], f[1], f[2], n);
                    s_tmp += v_tmp(11, f[0], f[1], f[2]); s_rsum += v_rsum(11, f[0], f[1], f[2]); s_mu += v_mu(11, f[0], f[1], f[2]);
                }
            const double rho_e = s_rho / n_int;
            const double zz_e[2] = {clip01t(s_rz[0] / s_rho), clip01t(s_rz[1] / s_rho)};
            double tmp_e = s_tmp / n_int, rsum_e = s_rsum / n_int;
            if (with_thermo) { rsum_e = eos_rsum(zz_e); tmp_e = eos_pbar(amrex::IntVect(l1[0], l1[1], l1[2])) / (rsum_e * rho_e); }
            ++n_checked;
            bool ok = close_to(crho(i, j, k), rho_e) && close_to(czz(i, j, k, 0), zz_e[0]) && close_to(czz(i, j, k, 1), zz_e[1]) && close_to(ctmp(i, j, k), tmp_e);
            if (lay == 1) ok = ok && close_to(crsum(i, j, k), rsum_e) && close_to(cmu(i, j, k), s_mu / n_int);
            else ok = ok && crsum(i, j, k) == v_rsum(1, i, j, k) && cmu(i, j, k) == v_mu(1, i, j, k);   // layer 2: RSUM and MU have one ghost layer in FDS
            if (!ok && ncand == 1) ++bad;
        });
    }
    amrex::ParallelDescriptor::ReduceLongSum(n_checked);
    amrex::ParallelDescriptor::ReduceLongSum(conflicts_expected);
    CfStats scr = sc;
    amrex::ParallelDescriptor::ReduceLongSum(scr.conflicts);
    CHECK_MSG(n_checked > 0, tag + ": covered face cells were checked");
    CHECK_MSG(bad == 0, tag + ": covered ghost cells, layers 1 and 2, follow the FDS rule (face fine layer only, mass-weighted ZZ, layer 2 = layer 1), bad=" + std::to_string(bad));
    CHECK_MSG(bad_untouched == 0, tag + ": uncovered and deep covered cells are untouched");
    CHECK_MSG(scr.conflicts == conflicts_expected && conflicts_expected > 0, tag + ": edge/corner conflict cells counted, expected " + std::to_string(conflicts_expected) + " got " + std::to_string(scr.conflicts));

    // ---- fine ghost cells (face-adjacent ones) ----
    long nf = 0, bad_f = 0;
    const amrex::Box patch = amrex::refine(amrex::Box(amrex::IntVect(4), amrex::IntVect(11)), ratio);
    for (amrex::MFIter mfi(Frho); mfi.isValid(); ++mfi) {
        auto frho = Frho.const_array(mfi); auto fzz = Fzz.const_array(mfi); auto ftmp = Ftmp.const_array(mfi); auto frsum = Frsum.const_array(mfi); auto fmu = Fmu.const_array(mfi);
        const amrex::Box vb = mfi.validbox();
        amrex::LoopOnCpu(amrex::grow(vb, 3), [&](int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            if (patch.contains(iv)) return;   // inside the patch: another fine box
            int nout = 0, dn = 0, layer = 0;   // face-adjacent ghost cell of THIS box (exactly one direction outside the box), outside the whole patch
            for (int d = 0; d < 3; ++d) {
                const int dd = std::max({vb.smallEnd(d) - iv[d], iv[d] - vb.bigEnd(d), 0});
                if (dd > 0) { ++nout; dn = d; layer = dd; }
            }
            if (nout != 1) return;            // edge and corner ghost cells are not defined by the FDS rule
            const amrex::IntVect c = amrex::coarsen(iv, ratio);
            const double rho_c = v_rho(1, c[0], c[1], c[2]);
            const double zz_e[2] = {clip01t((rho_c * v_zz(1, c[0], c[1], c[2], 0)) / rho_c), clip01t((rho_c * v_zz(1, c[0], c[1], c[2], 1)) / rho_c)};
            amrex::IntVect l1 = iv;
            if (layer == 2) l1[dn] += (iv[dn] < patch.smallEnd(dn)) ? 1 : -1;
            if (layer > 2) { if (frho(i, j, k) != -555.0) ++bad_f; return; }   // layer 3 untouched
            ++nf;
            double tmp_e = v_tmp(1, c[0], c[1], c[2]), rsum_e = v_rsum(1, c[0], c[1], c[2]);
            if (with_thermo) { rsum_e = eos_rsum(zz_e); tmp_e = eos_pbar(l1) / (rsum_e * rho_c); }
            bool ok = close_to(frho(i, j, k), rho_c) && close_to(fzz(i, j, k, 0), zz_e[0]) && close_to(fzz(i, j, k, 1), zz_e[1]) && close_to(ftmp(i, j, k), tmp_e);
            if (layer == 1) ok = ok && close_to(frsum(i, j, k), rsum_e) && close_to(fmu(i, j, k), v_mu(1, c[0], c[1], c[2]));
            // layer 2: RSUM and MU have one ghost layer (the cell is outside their FABs): nothing to check
            if (!ok && ++bad_f < 4 && getenv("RT_DEBUG")) std::fprintf(stderr, "bad fine ghost %d %d %d layer %d: rho %g (exp %g) zz0 %g (%g) tmp %g (%g) rsum %g (%g) mu %g (%g)\n", i, j, k, layer, frho(i,j,k), rho_c, fzz(i,j,k,0), zz_e[0], ftmp(i,j,k), tmp_e, frsum(i,j,k), rsum_e, fmu(i,j,k), v_mu(1, c[0], c[1], c[2]));
        });
    }
    amrex::ParallelDescriptor::ReduceLongSum(nf);
    CHECK_MSG(nf > 0, tag + ": fine face ghost cells were checked");
    CHECK_MSG(bad_f == 0, tag + ": fine ghost cells follow the FDS rule (coarse value injected, ZZ clipped, TMP/RSUM rebuilt, layer 2 = layer 1), bad=" + std::to_string(bad_f));
    (void)sf;
}


// ---- interface flux overwrite: area-sum of fine face fluxes onto coarse faces (no driver hooks needed) ----
double ff(int d, int i, int j, int k, int n) { return 0.5 + rv(100 + d, i, j, k, n) - 0.25 * n; }   // fine face flux (any sign)

void test_flux_override(int rr, bool periodic_x)
{
    const std::string tag = "flux override ratio " + std::to_string(rr) + (periodic_x ? " periodic x, patch at the low edge" : " interior patch");
    amrex::Box cdom(amrex::IntVect(0), amrex::IntVect(15));
    amrex::RealBox rb({0, 0, 0}, {1, 1, 1});
    amrex::Array<int, 3> per{periodic_x ? 1 : 0, 0, 0};
    amrex::Geometry cg(cdom, rb, 0, per);
    amrex::IntVect ratio(rr);
    amrex::BoxArray cba(cdom);
    cba.maxSize(8);
    amrex::DistributionMapping cdm(cba);
    const amrex::Box pc = periodic_x ? amrex::Box(amrex::IntVect(0, 4, 4), amrex::IntVect(7, 11, 11)) : amrex::Box(amrex::IntVect(4), amrex::IntVect(11));   // patch in coarse cells
    amrex::BoxArray fba(amrex::refine(pc, ratio));
    fba.maxSize(4 * rr);
    amrex::DistributionMapping fdm(fba);
    const int nscal = 3;
    std::vector<amrex::MultiFab> F;
    F.reserve(3);
    for (int d = 0; d < 3; ++d) {
        amrex::IntVect nodal(0); nodal[d] = 1;
        F.emplace_back(amrex::convert(fba, nodal), fdm, nscal, 0);
        for (amrex::MFIter mfi(F[d]); mfi.isValid(); ++mfi) {
            auto a = F[d].array(mfi);
            for (int n = 0; n < nscal; ++n) amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { a(i, j, k, n) = ff(d, i, j, k, n); });
        }
    }
    const amrex::MultiFab* fp[3] = {&F[0], &F[1], &F[2]};
    OverrideStats st;
    auto out = build_flux_overrides(cba, cdm, cg, fba, fdm, ratio, fp, &st);

    auto covered = [&](amrex::IntVect p) {
        if (periodic_x) { if (p[0] < 0) p[0] += 16; if (p[0] > 15) p[0] -= 16; }
        for (int e = 0; e < 3; ++e) if (p[e] < 0 || p[e] > 15) return -1;
        return pc.contains(p) ? 1 : 0;
    };
    long bad = 0, unique_low = 0, n_entries = 0;
    int ib = 0;
    for (amrex::MFIter mfi(cba, cdm); mfi.isValid(); ++mfi, ++ib) {
        for (int d = 0; d < 3; ++d) {
            amrex::IntVect nodal(0); nodal[d] = 1;
            const amrex::Box nb = amrex::convert(mfi.validbox(), nodal);
            // expected list, sorted by (k,j,i)
            std::vector<std::array<int, 3>> exp;
            amrex::LoopOnCpu(nb, [&](int i, int j, int k) {
                amrex::IntVect lo(i, j, k), hi(i, j, k); lo[d] -= 1;
                const int a = covered(lo), b = covered(hi);
                if (a >= 0 && b >= 0 && a != b) exp.push_back({i, j, k});
            });
            const FluxOverride* ov = nullptr;
            for (const auto& o : out[ib]) if (o.dir == d) ov = &o;
            if (exp.empty()) { if (ov) ++bad; continue; }
            if (!ov || ov->face != exp || ov->nscal != nscal || ov->value.size() != exp.size() * nscal) { ++bad; continue; }
            const int t1 = (d + 1) % 3, t2 = (d + 2) % 3;
            for (std::size_t q = 0; q < exp.size(); ++q) {
                const auto& f = exp[q];
                if (q > 0) { const auto& g = exp[q - 1]; if (!(std::make_tuple(g[2], g[1], g[0]) < std::make_tuple(f[2], f[1], f[0]))) ++bad; }
                int a = f[d]; if (periodic_x && d == 0 && a == 16) a = 0;   // the periodic image of face 0 carries the same flux
                for (int n = 0; n < nscal; ++n) {
                    double sum = 0;
                    for (int u = 0; u < rr; ++u) for (int v = 0; v < rr; ++v) {
                        int g[3]; g[d] = a * rr; g[t1] = f[t1] * rr + u; g[t2] = f[t2] * rr + v;
                        sum += ff(d, g[0], g[1], g[2], n);
                    }
                    if (std::abs(ov->value[q * nscal + n] - sum / (rr * rr)) > 1e-14 * std::max(1.0, std::abs(sum))) ++bad;
                }
                ++n_entries;
                if (mfi.validbox().contains(amrex::IntVect(f[0], f[1], f[2]))) ++unique_low;   // the box that owns the cell above the face: each face counted once
            }
        }
    }
    amrex::ParallelDescriptor::ReduceLongSum(unique_low);
    amrex::ParallelDescriptor::ReduceLongSum(n_entries);
    long st_entries = st.entries;
    amrex::ParallelDescriptor::ReduceLongSum(st_entries);
    CHECK_MSG(bad == 0, tag + ": interface faces, sorted face lists, shared faces listed in every holding box, values = area-weighted sum of the fine faces, bad=" + std::to_string(bad));
    CHECK_MSG(st_entries == n_entries && n_entries > 0, tag + ": entry count");
    // faces on the outline of the patch, each counted once by the box that holds the cell above the face: 6 sides x 8 x 8. Periodic x with the patch on the low edge:
    // x faces 8 and 0 (the high-end image, face 16, is listed too but its cell above wraps) plus 4 y/z sides.
    CHECK_MSG(unique_low == 384, tag + ": 384 interface faces, got " + std::to_string(unique_low));

    // empty set: no fine level -> nothing listed, one empty entry per local box
    amrex::BoxArray none;
    auto empty = build_flux_overrides(cba, cdm, cg, none, fdm, ratio, fp, nullptr);
    long nonempty = 0;
    for (auto& v : empty) nonempty += static_cast<long>(v.size());
    CHECK_MSG(nonempty == 0 && static_cast<int>(empty.size()) == static_cast<int>(out.size()), tag + ": empty fine level gives empty lists (no-op)");
}

}  // namespace

int main(int argc, char** argv)
{
    int one = 1;   // our only argument is the cases directory: do not let AMReX read it as an inputs file
    amrex::Initialize(one, argv);
    const std::string dir = argc > 1 ? argv[1] : ".";
    {
        test_ops(2, false);
        test_ops(4, false);
        test_ops(2, true);
        test_ops(4, true);
        test_fds_rules(2, false);
        test_fds_rules(2, true);
        test_fds_rules(4, false);
        test_fds_rules(4, true);
        test_flux_override(2, false);
        test_flux_override(4, false);
        test_flux_override(2, true);

        // ---- the ns2d_16_int_1to2_refinement-style hierarchy through AmrCore and the driver's LevelRegistry ----
        const std::string text = fdsrt_test::read_file(dir + "/ns2d_16_int_1to2.fds");
        auto meshes = fdsrt_test::meshes_from_text(text);
        Report rep;
        AmrParams p = parse_amr_params(text, rep);
        Hierarchy h;
        const bool ok = build_hierarchy_from_meshes(meshes, p, {true, false, true}, h, rep);
        CHECK(ok && rep.ok());
        if (ok) {
            // level 0 as the driver would hold it: the level-0 boxes in order, FDS-style objects adopted by the registry
            amrex::Geometry g0 = make_level0_geometry(h);
            amrex::BoxList bl0;
            for (const GridBox& g : h.levels[0].grids) bl0.push_back(to_amrex(g.box));
            amrex::BoxArray ba0(bl0);
            amrex::DistributionMapping dm0(ba0);
            fdsamr::DomainInfo dom{};
            dom.periodic[0] = 1; dom.periodic[2] = 1;
            dom.n_tracked = 2; dom.n_total = 2; dom.nranks = amrex::ParallelDescriptor::NProcs();
            fdsamr::Level l0 = fdsamr::make_layout_level(0, g0, ba0, dm0, amrex::IntVect(1), dom);
            fdsamr::Fields F0(l0, 2);
            fdsamr::SideData sd0(l0, fdsamr::layout_cell_walls(l0));
            fdsamr::LevelRegistry reg(dom, 2);
            reg.adopt_level0(l0, F0, sd0);

            RegridAmrCore core(h, p);
            core.init_static(reg, false, RegridAmrCore::DmFn(), &dm0);
            CHECK(core.finestLevel() == 1);
            CHECK(reg.num_levels() == 2 && reg.has_level(1) && reg.n_make == 1);
            CHECK(core.boxArray(0) == ba0 && core.boxArray(1).size() == 1);
            CHECK(core.boxArray(1)[0] == amrex::Box(amrex::IntVect(8, 0, 8), amrex::IntVect(23, 0, 23)));
            CHECK(core.refRatio(0) == amrex::IntVect(2, 1, 2));
            CHECK(reg.level(1).ref_ratio_from_parent == amrex::IntVect(2, 1, 2));
            CHECK(std::abs(core.Geom(1).CellSize(0) - 0.5 * core.Geom(0).CellSize(0)) < 1e-15);
            CHECK(std::abs(core.Geom(1).CellSize(1) - core.Geom(0).CellSize(1)) < 1e-15);   // hidden direction is not refined
            CHECK(reg.fields(1)["RHO"].boxArray() == core.boxArray(1));
            const amrex::iMultiFab* cov = reg.covered_mask(0);
            long ncov = 0;
            if (cov) for (amrex::MFIter mfi(*cov); mfi.isValid(); ++mfi) { auto a = cov->const_array(mfi); amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { ncov += a(i, j, k); }); }
            amrex::ParallelDescriptor::ReduceLongSum(ncov);
            CHECK_MSG(ncov == 64, "covered coarse cells under the fine patch: 8 x 8, got " + std::to_string(ncov));

            const std::vector<std::string> names = {"RHO", "ZZ", "TMP", "MU"};
            // (1) uniform state stays uniform across the interface; fine level: valid cells 2.5, ghost cells start as the sentinel -777 (layer 3, H/HS must keep it)
            for (int l = 0; l < 2; ++l)
                for (const auto& n : reg.fields(l).names()) reg.fields(l)[n].setVal(-777.0);
            auto uval = [](const std::string& n) { return n == "ZZ" ? 0.25 : 2.5; };   // ZZ is a mass fraction: the FDS ghost rule clips it to [0,1]
            for (const auto& n : names) if (reg.fields(0).has(n)) reg.fields(0)[n].setVal(uval(n));   // coarse level: everywhere (valid and ghost cells)
            for (const auto& n : names) if (reg.fields(1).has(n)) { amrex::MultiFab& m = reg.fields(1)[n]; for (amrex::MFIter mfi(m); mfi.isValid(); ++mfi) m[mfi].setVal<amrex::RunOn::Host>(uval(n), mfi.validbox()); }
            CfHookStats hs;
            fdsamr::CfGhostHook hook = make_cf_ghost_hook(reg, ThermoProvider(), &hs);
            // code 4 (RHO, ZZ, ...) and code 3 (US.. and HS, which are not transferred)
            for (const auto& n : names) if (reg.fields(1).has(n)) reg.fields(1).fill_ghosts(n);
            fdsamr::CfGhostRequest rq{1, 4, true, fdsamr::exchange_fields(4, true)};
            hook(rq);
            fdsamr::CfGhostRequest rq3{1, 3, true, fdsamr::exchange_fields(3, true)};
            hook(rq3);
            long nonuni = 0, nlayer12 = 0, beyond_untouched = 1;
            for (const auto& n : {std::string("RHO"), std::string("ZZ")}) {
                const amrex::MultiFab& m = reg.fields(1)[n];
                const amrex::Box fdom = reg.level(1).geom.Domain();
                for (amrex::MFIter mfi(m); mfi.isValid(); ++mfi) {
                    auto a = m.const_array(mfi);
                    const amrex::Box vb = mfi.validbox();
                    amrex::LoopOnCpu(mfi.fabbox(), [&](int i, int j, int k) {
                        const amrex::IntVect iv(i, j, k);
                        if (vb.contains(iv)) return;
                        const int layer = std::max({vb.smallEnd(0) - i, i - vb.bigEnd(0), vb.smallEnd(1) - j, j - vb.bigEnd(1), vb.smallEnd(2) - k, k - vb.bigEnd(2)});
                        const bool inside = fdom.contains(iv);   // periodic wrap: the patch is away from the x/z domain edges, y is not periodic
                        for (int c = 0; c < m.nComp(); ++c) {
                            if (inside && layer <= 2) { ++nlayer12; if (a(i, j, k, c) != uval(n)) ++nonuni; }
                            else if (inside && layer == 3 && a(i, j, k, c) != -777.0) beyond_untouched = 0;
                            else if (!inside && a(i, j, k, c) != -777.0) beyond_untouched = 0;
                        }
                    });
                }
            }
            amrex::ParallelDescriptor::ReduceLongSum(nonuni);
            amrex::ParallelDescriptor::ReduceLongSum(nlayer12);
            CHECK_MSG(nlayer12 > 0 && nonuni == 0, "uniform state stays uniform in ghost layers 1 and 2 of RHO and ZZ across the coarse-fine faces, bad=" + std::to_string(nonuni));
            CHECK_MSG(beyond_untouched == 1, "layer 3 of RHO and cells outside the non-periodic y edge are not written");
            { long v[3] = {hs.calls, hs.fine_ghost_cells, hs.covered_cells}; amrex::ParallelDescriptor::ReduceLongSum(v, 3); hs.calls = v[0] / amrex::ParallelDescriptor::NProcs(); hs.fine_ghost_cells = v[1]; hs.covered_cells = v[2]; }
            CHECK_MSG(hs.calls == 1 && hs.fine_ghost_cells > 0 && hs.covered_cells > 0, "hook did its work for code 4; code 3 (US, VS, WS, HS) has no transferable field and does nothing");
            long hs_bad = 0;
            { const amrex::MultiFab& m = reg.fields(1)["HS"]; for (amrex::MFIter mfi(m); mfi.isValid(); ++mfi) { auto a = m.const_array(mfi); amrex::LoopOnCpu(mfi.fabbox(), [&](int i, int j, int k) { if (a(i, j, k) != -777.0) ++hs_bad; }); } }
            CHECK_MSG(hs_bad == 0, "H/HS ghost cells are not transferred by this hook");

            // (2) average-down through the registry: distinct fine values, coarse covered cells become the 2x1x2 mean
            for (int l = 0; l < 2; ++l) fill_all(reg.fields(l)["RHO"], l == 0 ? gval : fval, false);
            average_down_registry(reg, {"RHO"});
            double maxdev = 0;
            {
                const amrex::MultiFab& c = reg.fields(0)["RHO"];
                for (amrex::MFIter mfi(c); mfi.isValid(); ++mfi) {
                    auto a = c.const_array(mfi);
                    amrex::LoopOnCpu(mfi.validbox() & amrex::Box(amrex::IntVect(4, 0, 4), amrex::IntVect(11, 0, 11)), [&](int i, int j, int k) {
                        double m = 0;
                        for (int di = 0; di < 2; ++di) for (int dk = 0; dk < 2; ++dk) m += fval(2 * i + di, j, 2 * k + dk, 0);
                        m /= 4;
                        maxdev = std::max(maxdev, std::abs(a(i, j, k) - m) / std::abs(m));
                    });
                }
            }
            amrex::ParallelDescriptor::ReduceRealMax(maxdev);
            CHECK_MSG(maxdev < 1e-14, "registry average-down: covered cell = mean of its 4 fine cells, max rel dev " + std::to_string(maxdev));
        }
    }
    const long fails = fdstest::report("regrid_transport amrcore (R2a)");
    amrex::Finalize();
    return fails == 0 ? 0 : 1;
}
