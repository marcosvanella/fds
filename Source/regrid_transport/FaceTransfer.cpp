// FaceTransfer.cpp: see FaceTransfer.H.
#include "FaceTransfer.H"

#include <AMReX_Array4.H>
#include <AMReX_BoxArray.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_MFIter.H>
#include <AMReX_Reduce.H>

#include <algorithm>
#include <cmath>

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

void prolong_faces_level(const std::array<amrex::MultiFab*, 3>& fine, const std::array<const amrex::MultiFab*, 3>& crse, const amrex::Geometry& crse_geom,
                         const amrex::IntVect& ratio, const std::array<const amrex::MultiFab*, 3>& old_fine)
{
    const amrex::BoxArray fcells = amrex::convert(fine[0]->boxArray(), amrex::IntVect::TheCellVector());
    for (int b = 0; b < static_cast<int>(fcells.size()); ++b) {
        const amrex::Box& bx = fcells[b];
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(amrex::refine(amrex::coarsen(bx, ratio), ratio) == bx, "prolong_faces_level: a fine box is not a union of whole coarse cells");
    }
    amrex::BoxArray ccells = fcells;
    ccells.coarsen(ratio);
    std::array<amrex::MultiFab, 3> tmp;
    for (int d = 0; d < 3; ++d) {
        amrex::BoxArray cb = amrex::convert(ccells, amrex::IntVect::TheDimensionVector(d));
        tmp[d].define(cb, fine[d]->DistributionMap(), 1, 0);
        tmp[d].setVal(0.0);
        tmp[d].ParallelCopy(*crse[d], 0, 0, 1, 0, 0, crse_geom.periodicity());
    }
    for (amrex::MFIter mfi(*fine[0]); mfi.isValid(); ++mfi) {
        const amrex::Box cells = amrex::enclosedCells(mfi.validbox());
        prolong_faces_normal_linear({&tmp[0][mfi], &tmp[1][mfi], &tmp[2][mfi]}, {&(*fine[0])[mfi], &(*fine[1])[mfi], &(*fine[2])[mfi]}, cells, ratio);
    }
    for (int d = 0; d < 3; ++d)
        if (old_fine[d]) fine[d]->ParallelCopy(*old_fine[d], 0, 0, 1, 0, 0);
}

void average_down_faces(const std::array<const amrex::MultiFab*, 3>& fine, const std::array<amrex::MultiFab*, 3>& crse, const amrex::IntVect& ratio)
{
    for (int d = 0; d < 3; ++d) {
        // coarse faces under the fine level: coarse-aligned temporary on the coarsened fine BoxArray, then copy into the real coarse array
        amrex::BoxArray fcells = amrex::convert(fine[d]->boxArray(), amrex::IntVect::TheCellVector());
        fcells.coarsen(ratio);
        amrex::BoxArray cb = amrex::convert(fcells, amrex::IntVect::TheDimensionVector(d));
        amrex::MultiFab tmp(cb, fine[d]->DistributionMap(), 1, 0);
        tmp.setVal(0.0);
        for (amrex::MFIter mfi(tmp); mfi.isValid(); ++mfi) {
            const amrex::Box cfaces = mfi.validbox();   // coarse faces nodal in d
            auto c = tmp.array(mfi);
            auto f = fine[d]->const_array(mfi);
            const amrex::IntVect r = ratio;
            amrex::LoopOnCpu(cfaces, [&](int i, int j, int k) {
                // fine faces on this coarse face: fine index r_d*c_d in the normal direction, r_e*c_e + 0..r_e-1 in the others
                const int ci[3] = {i, j, k};
                int lo[3], hi[3];
                for (int e = 0; e < 3; ++e) {
                    lo[e] = ci[e] * r[e];
                    hi[e] = (e == d) ? lo[e] : lo[e] + r[e] - 1;
                }
                double s = 0.0; int n = 0;
                for (int kk = lo[2]; kk <= hi[2]; ++kk) for (int jj = lo[1]; jj <= hi[1]; ++jj) for (int ii = lo[0]; ii <= hi[0]; ++ii) { s += f(ii, jj, kk); ++n; }
                c(i, j, k) = s / n;
            });
        }
        crse[d]->ParallelCopy(tmp, 0, 0, 1, 0, 0);
    }
}

double max_divergence_error_local(const std::array<const amrex::MultiFab*, 3>& vel, const amrex::Geometry& geom, const amrex::MultiFab* D)
{
    double worst = 0.0;
    const double rdx[3] = {1.0 / geom.CellSize(0), 1.0 / geom.CellSize(1), 1.0 / geom.CellSize(2)};
    for (amrex::MFIter mfi(*vel[0]); mfi.isValid(); ++mfi) {
        const amrex::Box cells = amrex::enclosedCells(mfi.validbox());
        auto u = vel[0]->const_array(mfi); auto v = vel[1]->const_array(mfi); auto w = vel[2]->const_array(mfi);
        amrex::Array4<const amrex::Real> dd;
        if (D) dd = D->const_array(mfi);
        amrex::LoopOnCpu(cells, [&](int i, int j, int k) {
            const double div = (u(i + 1, j, k) - u(i, j, k)) * rdx[0] + (v(i, j + 1, k) - v(i, j, k)) * rdx[1] + (w(i, j, k + 1) - w(i, j, k)) * rdx[2];
            worst = std::max(worst, std::abs(div - (D ? dd(i, j, k) : 0.0)));
        });
    }
    return worst;
}

}  // namespace fdsrt
