// CellTransfer.cpp: see CellTransfer.H.
#include "CellTransfer.H"

#include <AMReX_BoxArray.H>
#include <AMReX_ParallelDescriptor.H>

namespace fdsrt {

ProlongStats prolong_conserved(amrex::MultiFab& fine, const amrex::MultiFab& crse, const amrex::Geometry& crse_geom, const amrex::IntVect& ratio,
                               const ProlongOpts& /*opt*/, const amrex::MultiFab* old_fine)
{
    using namespace amrex;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(fine.nComp() == crse.nComp(), "prolong_conserved: component counts differ");
    const int nc = fine.nComp();
    BoxArray cba = fine.boxArray();
    for (int i = 0; i < static_cast<int>(cba.size()); ++i) {
        const Box b = cba[i];
        for (int d = 0; d < 3; ++d)
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(b.smallEnd(d) % ratio[d] == 0 && (b.bigEnd(d) + 1) % ratio[d] == 0, "prolong_conserved: fine box is not a union of whole coarse cells");
    }
    cba.coarsen(ratio);
    MultiFab tmp(cba, fine.DistributionMap(), nc, 0);
    tmp.setVal(0.0);
    tmp.ParallelCopy(crse, 0, 0, nc, 0, 0, crse_geom.periodicity());   // the parents (valid coarse cells, periodic images included)

    ProlongStats st;
    for (MFIter mfi(fine); mfi.isValid(); ++mfi) {
        auto q = tmp.const_array(mfi);
        auto f = fine.array(mfi);
        st.parents += tmp.box(mfi.index()).numPts();
        LoopOnCpu(mfi.validbox(), nc, [&](int i, int j, int k, int n) { f(i, j, k, n) = q(amrex::coarsen(i, ratio[0]), amrex::coarsen(j, ratio[1]), amrex::coarsen(k, ratio[2]), n); });
    }
    if (old_fine) fine.ParallelCopy(*old_fine, 0, 0, nc, 0, 0);   // bitwise copy of the cells that already had fine data
    long v[4] = {st.parents, st.limited, st.clips, st.unfixable};
    ParallelAllReduce::Sum(v, 4, ParallelContext::CommunicatorSub());
    st.parents = v[0]; st.limited = v[1]; st.clips = v[2]; st.unfixable = v[3];
    return st;
}

}  // namespace fdsrt
