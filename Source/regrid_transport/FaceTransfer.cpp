// FaceTransfer.cpp: see FaceTransfer.H.
#include "FaceTransfer.H"

#include <AMReX_Array4.H>
#include <AMReX_GpuLaunch.H>

namespace fdsrt {

void prolong_faces_normal_linear(const std::array<const amrex::FArrayBox*, 3>& crse, const std::array<amrex::FArrayBox*, 3>& fine,
                                 const amrex::Box& fine_cells, const amrex::IntVect& ratio)
{
    for (int d = 0; d < 3; ++d) {
        const amrex::Box fb = amrex::convert(fine_cells, amrex::IntVect::TheDimensionVector(d));
        auto c = crse[d]->const_array();
        auto f = fine[d]->array();
        const amrex::IntVect r = ratio;
        amrex::ParallelFor(fb, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const int fi[3] = {i, j, k};
            int ci[3];
            for (int e = 0; e < 3; ++e) ci[e] = amrex::coarsen(fi[e], r[e]);          // floor division, also for negative indices
            const int m = fi[d] - ci[d] * r[d];                                       // 0 .. r-1 within the coarse cell, 0 = on a coarse face
            const amrex::Real w = static_cast<amrex::Real>(m) / static_cast<amrex::Real>(r[d]);
            amrex::Real v = c(ci[0], ci[1], ci[2]);
            if (m != 0) {
                const int cj[3] = {ci[0] + (d == 0), ci[1] + (d == 1), ci[2] + (d == 2)};
                v = (amrex::Real(1) - w) * v + w * c(cj[0], cj[1], cj[2]);
            }
            f(i, j, k) = v;
        });
    }
}

}  // namespace fdsrt
