// FluxOverrideOps.cpp: see FluxOverrideOps.H.
#include "FluxOverrideOps.H"

#include <AMReX_MFIter.H>

#include <algorithm>

namespace fdsrt {

std::vector<std::vector<FluxOverride>> build_flux_overrides(const amrex::BoxArray& crse_ba, const amrex::DistributionMapping& crse_dm, const amrex::Geometry& cg,
                                                           const amrex::BoxArray& fine_ba, const amrex::DistributionMapping& fine_dm, const amrex::IntVect& ratio,
                                                           const amrex::MultiFab* const fine_flux[3], OverrideStats* stats)
{
    std::vector<std::vector<FluxOverride>> out;
    for (amrex::MFIter mfi(crse_ba, crse_dm); mfi.isValid(); ++mfi) out.emplace_back();
    if (fine_ba.empty()) return out;
    const int nscal = fine_flux[0]->nComp();
    amrex::BoxArray cfba = fine_ba;
    cfba.coarsen(ratio);
    const amrex::Box dom = cg.Domain();

    for (int d = 0; d < 3; ++d) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fine_flux[d] && fine_flux[d]->nComp() == nscal && fine_flux[d]->ixType().nodeCentered(d) && fine_flux[d]->nGrow() == 0,
                                         "build_flux_overrides: fine flux must be nodal in its direction, no ghost cells, same component count in all directions");
        const int t1 = (d + 1) % 3, t2 = (d + 2) % 3;
        const double aro = 1.0 / (static_cast<double>(ratio[t1]) * static_cast<double>(ratio[t2]));   // A_f / A_c, uniform cells
        amrex::IntVect nodal(0); nodal[d] = 1;
        // coarse-face values from the fine faces, on the coarsened fine boxes (nodal in d): the fine face a*r_d coincides with coarse face a
        amrex::BoxArray aux_ba = amrex::convert(cfba, nodal);
        amrex::MultiFab aux(aux_ba, fine_dm, nscal, 0);
        aux.setVal(0.0);
        for (amrex::MFIter mfi(aux); mfi.isValid(); ++mfi) {
            auto ar = aux.array(mfi);
            auto fl = fine_flux[d]->const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int ci, int cj, int ck) {
                const int cc[3] = {ci, cj, ck};
                int lo[3], hi[3];
                for (int e = 0; e < 3; ++e) { lo[e] = cc[e] * ratio[e]; hi[e] = lo[e] + ratio[e] - 1; }
                lo[d] = hi[d] = cc[d] * ratio[d];
                for (int n = 0; n < nscal; ++n) {
                    double sum = 0.0;
                    for (int k = lo[2]; k <= hi[2]; ++k)
                        for (int j = lo[1]; j <= hi[1]; ++j)
                            for (int i = lo[0]; i <= hi[0]; ++i) sum = sum + aro * fl(i, j, k, n);
                    ar(ci, cj, ck, n) = sum;
                }
            });
        }
        amrex::BoxArray cn_ba = amrex::convert(crse_ba, nodal);
        amrex::MultiFab cface(cn_ba, crse_dm, nscal, 0);
        cface.setVal(0.0);
        cface.ParallelCopy(aux, 0, 0, nscal, 0, 0, cg.periodicity());

        auto covered = [&](amrex::IntVect p) {   // 1 covered, 0 not, -1 outside a non-periodic domain edge
            for (int e = 0; e < 3; ++e) {
                if (p[e] < dom.smallEnd(e)) { if (!cg.isPeriodic(e)) return -1; p[e] += dom.length(e); }
                else if (p[e] > dom.bigEnd(e)) { if (!cg.isPeriodic(e)) return -1; p[e] -= dom.length(e); }
            }
            return cfba.contains(p) ? 1 : 0;
        };
        int ib = 0;
        for (amrex::MFIter mfi(crse_ba, crse_dm); mfi.isValid(); ++mfi, ++ib) {
            // the nodal box of this coarse box
            const amrex::Box nb = amrex::convert(mfi.validbox(), nodal);
            auto cv = cface.const_array(mfi);
            FluxOverride ov;
            ov.dir = d;
            ov.nscal = nscal;
            amrex::LoopOnCpu(nb, [&](int i, int j, int k) {   // k outer, j, i inner: sorted by (k,j,i)
                const amrex::IntVect hi_cell(i, j, k);
                amrex::IntVect lo_cell = hi_cell;
                lo_cell[d] -= 1;
                const int cl = covered(lo_cell), ch = covered(hi_cell);
                if (cl < 0 || ch < 0 || cl == ch) return;
                ov.face.push_back({i, j, k});
                for (int n = 0; n < nscal; ++n) ov.value.push_back(cv(i, j, k, n));
            });
            if (!ov.face.empty()) {
                if (stats) stats->entries += static_cast<long>(ov.face.size());
                out[ib].push_back(std::move(ov));
            }
        }
    }
    return out;
}

}  // namespace fdsrt
