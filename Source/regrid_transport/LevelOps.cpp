// LevelOps.cpp: see LevelOps.H. Host loops; the grids of a level are small and this runs once per exchange, a device version comes with the GPU phase.
#include "LevelOps.H"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

namespace fdsrt {

void average_down_cells(const amrex::MultiFab& fine, amrex::MultiFab& crse, const amrex::Geometry& fine_geom, const amrex::Geometry& crse_geom,
                        const amrex::IntVect& ratio, int scomp, int ncomp)
{
    amrex::average_down(fine, crse, fine_geom, crse_geom, scomp, ncomp, ratio);
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

}  // namespace fdsrt
