// LevelOps.cpp: see LevelOps.H. Host loops; the grids of a level are small and this runs once per exchange, a device version comes with the GPU phase.
#include "LevelOps.H"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <array>

namespace fdsrt {

void average_down_cells(const amrex::MultiFab& fine, amrex::MultiFab& crse, const amrex::Geometry& fine_geom, const amrex::Geometry& crse_geom,
                        const amrex::IntVect& ratio, int scomp, int ncomp)
{
    amrex::average_down(fine, crse, fine_geom, crse_geom, scomp, ncomp, ratio);
}

void average_down_species(const amrex::MultiFab& rho_f, const amrex::MultiFab& zz_f, amrex::MultiFab& rho_c, amrex::MultiFab& zz_c, const amrex::iMultiFab& covered_c,
                          const amrex::Geometry& fine_geom, const amrex::Geometry& crse_geom, const amrex::IntVect& ratio)
{
    const int ns = zz_f.nComp();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(zz_c.nComp() == ns, "average_down_species: component counts differ");
    amrex::MultiFab rzf(zz_f.boxArray(), zz_f.DistributionMap(), ns, 0), rzc(zz_c.boxArray(), zz_c.DistributionMap(), ns, 0);
    for (amrex::MFIter mfi(rzf); mfi.isValid(); ++mfi) {
        auto r = rho_f.const_array(mfi); auto z = zz_f.const_array(mfi); auto o = rzf.array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { for (int n = 0; n < ns; ++n) o(i, j, k, n) = r(i, j, k) * z(i, j, k, n); });
    }
    for (amrex::MFIter mfi(rzc); mfi.isValid(); ++mfi) {
        auto r = rho_c.const_array(mfi); auto z = zz_c.const_array(mfi); auto o = rzc.array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { for (int n = 0; n < ns; ++n) o(i, j, k, n) = r(i, j, k) * z(i, j, k, n); });
    }
    average_down_cells(rzf, rzc, fine_geom, crse_geom, ratio, 0, ns);
    average_down_cells(rho_f, rho_c, fine_geom, crse_geom, ratio, 0, 1);
    for (amrex::MFIter mfi(rzc); mfi.isValid(); ++mfi) {
        auto q = rzc.const_array(mfi); auto m = covered_c.const_array(mfi);
        auto rho = rho_c.const_array(mfi); auto z = zz_c.array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            if (!m(i, j, k)) return;
            for (int n = 0; n < ns; ++n) z(i, j, k, n) = q(i, j, k, n) / rho(i, j, k);
        });
    }
}

long fill_cf_ghosts_pc(amrex::MultiFab& fine, const amrex::MultiFab& crse, const amrex::Geometry& fine_geom, const amrex::Geometry& crse_geom,
                       const amrex::IntVect& ratio, int nlayers, int scomp, int ncomp)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nlayers <= fine.nGrow(), "fill_cf_ghosts_pc: more layers than the MultiFab has");
    for (int d = 0; d < 3; ++d)
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ratio[d] == 1 || nlayers <= ratio[d], "fill_cf_ghosts_pc: layers beyond the ratio would reach a second coarse cell");
    const amrex::BoxArray& fba = fine.boxArray();
    const amrex::Box fdom = fine_geom.Domain();

    // coarse values under the grown fine boxes (including periodic images): a MultiFab on the coarsened grown fine boxes (they may overlap each other)
    amrex::BoxArray cba = fba;
    cba.grow(nlayers);
    cba.coarsen(ratio);
    amrex::MultiFab patch(cba, fine.DistributionMap(), ncomp, 0);
    patch.setVal(0.0);
    patch.ParallelCopy(crse, scomp, 0, ncomp, 0, 0, crse_geom.periodicity());

    long nfilled = 0;
    const amrex::IntVect dlen = fdom.size();
    for (amrex::MFIter mfi(fine); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        const amrex::Box gb = amrex::grow(vb, nlayers);
        auto fa = fine.array(mfi);
        auto pa = patch.const_array(mfi);
        amrex::LoopOnCpu(gb, [&](int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            if (vb.contains(iv)) return;
            // position of the cell inside the domain: wrap periodic directions, skip a cell outside a non-periodic edge
            amrex::IntVect w = iv;
            for (int d = 0; d < 3; ++d) {
                if (w[d] < fdom.smallEnd(d) || w[d] > fdom.bigEnd(d)) {
                    if (!fine_geom.isPeriodic(d)) return;
                    w[d] += (w[d] < fdom.smallEnd(d)) ? dlen[d] : -dlen[d];
                }
            }
            if (fba.contains(w)) return;   // a fine grid holds this cell: the same-level exchange filled it
            const amrex::IntVect c = amrex::coarsen(iv, ratio);
            for (int n = 0; n < ncomp; ++n) fa(i, j, k, scomp + n) = pa(c[0], c[1], c[2], n);
            ++nfilled;
        });
    }
    return nfilled;
}


// ---------------------------------------------------------------------------------------------------------------------------------------------------------
// FDS ghost rules (see LevelOps.H)
// ---------------------------------------------------------------------------------------------------------------------------------------------------------
namespace {

// Coarse values under the grown fine boxes (periodic images included), one MultiFab on the coarsened grown fine boxes.
amrex::MultiFab coarse_patch(const amrex::MultiFab& fine_layout, const amrex::MultiFab& crse, const amrex::Geometry& cg, const amrex::IntVect& ratio, int nlayers)
{
    amrex::BoxArray cba = fine_layout.boxArray();
    cba.grow(nlayers);
    cba.coarsen(ratio);
    amrex::MultiFab p(cba, fine_layout.DistributionMap(), crse.nComp(), 0);
    p.setVal(0.0);
    p.ParallelCopy(crse, 0, 0, crse.nComp(), 0, 0, cg.periodicity());
    return p;
}

// Wraps periodic directions; false when the cell lies outside a non-periodic domain edge.
bool wrap_cell(amrex::IntVect& w, const amrex::Geometry& g)
{
    const amrex::Box& dom = g.Domain();
    for (int d = 0; d < 3; ++d) {
        if (w[d] < dom.smallEnd(d)) { if (!g.isPeriodic(d)) return false; w[d] += dom.length(d); }
        else if (w[d] > dom.bigEnd(d)) { if (!g.isPeriodic(d)) return false; w[d] -= dom.length(d); }
    }
    return true;
}

double clip01(double x) { return std::max(0.0, std::min(1.0, x)); }

}  // namespace

CfStats fill_fine_ghosts_fds(const ScalarStage& fine, const ScalarStage& crse, const amrex::Geometry& fg, const amrex::Geometry& cg, const amrex::IntVect& ratio,
                             const Thermo* th)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fine.rho && crse.rho, "fill_fine_ghosts_fds: RHO/RHOS is required on both levels");
    for (int d = 0; d < 3; ++d) AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ratio[d] == 1 || ratio[d] >= 2, "fill_fine_ghosts_fds: ratio");
    const amrex::BoxArray& fba = fine.rho->boxArray();
    const int nzz = fine.zz ? fine.zz->nComp() : 0;
    amrex::MultiFab Prho = coarse_patch(*fine.rho, *crse.rho, cg, ratio, 2);
    amrex::MultiFab Pzz, Ptmp, Prsum;
    std::vector<amrex::MultiFab> Pmean;
    if (fine.zz) Pzz = coarse_patch(*fine.rho, *crse.zz, cg, ratio, 2);
    if (fine.tmp) Ptmp = coarse_patch(*fine.rho, *crse.tmp, cg, ratio, 2);
    if (fine.rsum) Prsum = coarse_patch(*fine.rho, *crse.rsum, cg, ratio, 1);
    for (std::size_t k = 0; k < fine.mean.size(); ++k) Pmean.push_back(coarse_patch(*fine.rho, *crse.mean[k], cg, ratio, 1));
    const bool use_th = th && th->rsum && th->pbar && fine.zz;

    CfStats st;
    std::vector<double> zg(std::max(nzz, 1));
    for (amrex::MFIter mfi(*fine.rho); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        const amrex::Box gb = amrex::grow(vb, 2);
        auto rho = fine.rho->array(mfi);
        auto prho = Prho.const_array(mfi);
        amrex::Array4<double> zz, tmp, rsum;
        amrex::Array4<const double> pzz, ptmp, prsum;
        if (fine.zz) { zz = fine.zz->array(mfi); pzz = Pzz.const_array(mfi); }
        if (fine.tmp) { tmp = fine.tmp->array(mfi); ptmp = Ptmp.const_array(mfi); }
        if (fine.rsum) { rsum = fine.rsum->array(mfi); prsum = Prsum.const_array(mfi); }
        const int ng_rho = fine.rho->nGrow(), ng_zz = fine.zz ? fine.zz->nGrow() : 0, ng_tmp = fine.tmp ? fine.tmp->nGrow() : 0;
        const int ng_rsum = fine.rsum ? fine.rsum->nGrow() : 0;
        amrex::LoopOnCpu(gb, [&](int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            if (vb.contains(iv)) return;
            amrex::IntVect w = iv;
            if (!wrap_cell(w, fg)) return;      // physical boundary: not ours
            if (fba.contains(w)) return;        // held by another fine box: same-level exchange
            int dist[3], layer = 0, dn = 0;
            for (int d = 0; d < 3; ++d) {
                dist[d] = std::max({vb.smallEnd(d) - iv[d], iv[d] - vb.bigEnd(d), 0});
                if (dist[d] > layer) { layer = dist[d]; dn = d; }
            }
            const amrex::IntVect c = amrex::coarsen(iv, ratio);
            // layer-1 cell of this ghost cell (for the background pressure of TMP: layer 2 copies layer 1)
            amrex::IntVect l1 = iv;
            if (layer == 2) l1[dn] += (iv[dn] < vb.smallEnd(dn)) ? 1 : -1;
            const double rho_c = prho(c[0], c[1], c[2]);
            const double rho_g = 0.0 + 1.0 * rho_c;      // RHO_OTHER = sum ARO*RHO, ARO = 1 (coarse cell is larger than the fine face)
            if (layer <= std::min(2, ng_rho)) { rho(i, j, k) = rho_g; ++st.ghost_cells; }
            for (int n = 0; n < nzz; ++n) zg[n] = clip01(((1.0 * rho_c) * pzz(c[0], c[1], c[2], n)) / rho_g);
            if (fine.zz && layer <= std::min(2, ng_zz)) for (int n = 0; n < nzz; ++n) zz(i, j, k, n) = zg[n];
            double rsum_g = 0, tmp_g = 0;
            if (use_th) {
                rsum_g = th->rsum(zg.data());
                tmp_g = th->pbar(l1) / (rsum_g * rho_g);
            }
            if (fine.tmp && layer <= std::min(2, ng_tmp)) tmp(i, j, k) = use_th ? tmp_g : ptmp(c[0], c[1], c[2]);
            if (fine.rsum && layer <= std::min(1, ng_rsum)) rsum(i, j, k) = use_th ? rsum_g : prsum(c[0], c[1], c[2]);
        });
        for (std::size_t m = 0; m < fine.mean.size(); ++m) {
            auto a = fine.mean[m]->array(mfi);
            auto pa = Pmean[m].const_array(mfi);
            const int ng = fine.mean[m]->nGrow();
            amrex::LoopOnCpu(amrex::grow(vb, std::min(1, ng)), [&](int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                if (vb.contains(iv)) return;
                amrex::IntVect w = iv;
                if (!wrap_cell(w, fg) || fba.contains(w)) return;
                const amrex::IntVect c = amrex::coarsen(iv, ratio);
                a(i, j, k) = pa(c[0], c[1], c[2]) / 1.0;   // N_INT_CELLS = 1: the single coarse cell next to the face
            });
        }
    }
    return st;
}

CfStats fill_covered_ghosts_fds(const ScalarStage& crse, const ScalarStage& fine, const amrex::Geometry& fg, const amrex::Geometry& cg, const amrex::IntVect& ratio,
                                const Thermo* th)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fine.rho && crse.rho, "fill_covered_ghosts_fds: RHO/RHOS is required on both levels");
    const amrex::BoxArray& fba = fine.rho->boxArray();
    const amrex::DistributionMapping& fdm = fine.rho->DistributionMap();
    amrex::BoxArray cfba = fba;
    cfba.coarsen(ratio);
    const int nzz = fine.zz ? fine.zz->nComp() : 0;
    // plain sums of the fine cells next to a face: [0] sum(ARO rho), [1..nzz] sum((ARO rho) ZZ), then the means (tmp, rsum, mean fields)
    std::vector<const amrex::MultiFab*> ml;
    if (fine.tmp) ml.push_back(fine.tmp);
    if (fine.rsum) ml.push_back(fine.rsum);
    for (auto* m : fine.mean) ml.push_back(m);
    const int nml = static_cast<int>(ml.size());
    const int nc = 1 + nzz + nml;
    const int i_tmp = fine.tmp ? 1 + nzz : -1;
    const int i_rsum = fine.rsum ? 1 + nzz + (fine.tmp ? 1 : 0) : -1;
    const int i_mean0 = 1 + nzz + (fine.tmp ? 1 : 0) + (fine.rsum ? 1 : 0);

    // aux[d][s]: for every covered coarse cell, the sums over the fine layer on its low (s = 0) or high (s = 1) face in direction d
    std::array<std::array<amrex::MultiFab, 2>, 3> A;
    const amrex::Real* fdx = fg.CellSize();
    const amrex::Real* cdx = cg.CellSize();
    for (int d = 0; d < 3; ++d)
        for (int s = 0; s < 2; ++s) {
            amrex::MultiFab aux(cfba, fdm, nc, 0);
            aux.setVal(0.0);
            const int t1 = (d + 1) % 3, t2 = (d + 2) % 3;
            const double aro = std::min(1.0, (fdx[t1] * fdx[t2]) / (cdx[t1] * cdx[t2]));
            for (amrex::MFIter mfi(aux); mfi.isValid(); ++mfi) {
                auto ar = aux.array(mfi);
                auto frho = fine.rho->const_array(mfi);
                amrex::Array4<const double> fzz;
                if (fine.zz) fzz = fine.zz->const_array(mfi);
                std::vector<amrex::Array4<const double>> fm;
                for (auto* m : ml) fm.push_back(m->const_array(mfi));
                amrex::LoopOnCpu(mfi.validbox(), [&](int ci, int cj, int ck) {
                    const int cc[3] = {ci, cj, ck};
                    int lo[3], hi[3];
                    for (int e = 0; e < 3; ++e) { lo[e] = cc[e] * ratio[e]; hi[e] = lo[e] + ratio[e] - 1; }
                    lo[d] = hi[d] = (s == 0) ? cc[d] * ratio[d] : cc[d] * ratio[d] + ratio[d] - 1;
                    double s_rho = 0.0;
                    std::vector<double> s_rz(std::max(nzz, 1), 0.0), s_m(std::max(nml, 1), 0.0);
                    for (int k = lo[2]; k <= hi[2]; ++k)
                        for (int j = lo[1]; j <= hi[1]; ++j)
                            for (int i = lo[0]; i <= hi[0]; ++i) {
                                s_rho = s_rho + aro * frho(i, j, k);
                                for (int n = 0; n < nzz; ++n) s_rz[n] = s_rz[n] + aro * frho(i, j, k) * fzz(i, j, k, n);
                                for (int m = 0; m < nml; ++m) s_m[m] = s_m[m] + fm[m](i, j, k);
                            }
                    ar(ci, cj, ck, 0) = s_rho;
                    for (int n = 0; n < nzz; ++n) ar(ci, cj, ck, 1 + n) = s_rz[n];
                    for (int m = 0; m < nml; ++m) ar(ci, cj, ck, 1 + nzz + m) = s_m[m];
                });
            }
            A[d][s].define(crse.rho->boxArray(), crse.rho->DistributionMap(), nc, 1);
            A[d][s].setVal(0.0);
            A[d][s].ParallelCopy(aux, 0, 0, nc, 0, 1, cg.periodicity());
        }

    CfStats st;
    const amrex::Box cdom = cg.Domain();
    auto covered = [&](const amrex::IntVect& p) {   // 1 covered, 0 not covered, -1 outside the domain
        amrex::IntVect w = p;
        if (!wrap_cell(w, cg)) return -1;
        return cfba.contains(w) ? 1 : 0;
    };
    struct Cand { int d, s, layer; amrex::IntVect src; };
    std::vector<double> zg(std::max(nzz, 1));
    for (amrex::MFIter mfi(*crse.rho); mfi.isValid(); ++mfi) {
        auto rho = crse.rho->array(mfi);
        amrex::Array4<double> zz, tmp, rsum;
        if (crse.zz) zz = crse.zz->array(mfi);
        if (crse.tmp) tmp = crse.tmp->array(mfi);
        if (crse.rsum) rsum = crse.rsum->array(mfi);
        std::vector<amrex::Array4<double>> cm;
        for (auto* m : crse.mean) cm.push_back(m->array(mfi));
        std::array<std::array<amrex::Array4<const double>, 2>, 3> a;
        for (int d = 0; d < 3; ++d) for (int s = 0; s < 2; ++s) a[d][s] = A[d][s].const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            const amrex::IntVect x(i, j, k);
            if (covered(x) != 1) return;
            std::vector<Cand> cand;
            for (int layer = 1; layer <= 2; ++layer)
                for (int d = 0; d < 3; ++d)
                    for (int s = 0; s < 2; ++s) {
                        const int nb = (s == 0) ? -1 : 1;   // direction towards the uncovered coarse cell
                        amrex::IntVect y1 = x, y2 = x;
                        y1[d] += nb; y2[d] += 2 * nb;
                        if (layer == 1) { if (covered(y1) == 0) cand.push_back({d, s, 1, x}); }
                        else if (covered(y1) == 1 && covered(y2) == 0) cand.push_back({d, s, 2, y1});
                    }
            if (cand.empty()) return;
            if (cand.size() > 1) ++st.conflicts;
            const Cand& c = cand[0];
            const auto& av = a[c.d][c.s];
            const int p = c.src[0], q = c.src[1], r = c.src[2];
            const double s_rho = av(p, q, r, 0);
            rho(i, j, k) = s_rho;
            for (int n = 0; n < nzz; ++n) zg[n] = clip01(av(p, q, r, 1 + n) / s_rho);
            if (crse.zz) for (int n = 0; n < nzz; ++n) zz(i, j, k, n) = zg[n];
            const bool use_th = th && th->rsum && th->pbar && crse.zz;
            const double n_int = static_cast<double>(ratio[(c.d + 1) % 3] * ratio[(c.d + 2) % 3]);
            double rsum_g = 0;
            if (use_th) rsum_g = th->rsum(zg.data());
            if (crse.tmp) tmp(i, j, k) = use_th ? th->pbar(c.src) / (rsum_g * s_rho) : av(p, q, r, i_tmp) / n_int;
            if (c.layer == 1) {
                if (crse.rsum) rsum(i, j, k) = use_th ? rsum_g : av(p, q, r, i_rsum) / n_int;
                for (std::size_t m = 0; m < cm.size(); ++m) cm[m](i, j, k) = av(p, q, r, i_mean0 + static_cast<int>(m)) / n_int;
            }
            ++st.ghost_cells;
        });
    }
    (void)cdom;
    return st;
}

}  // namespace fdsrt
