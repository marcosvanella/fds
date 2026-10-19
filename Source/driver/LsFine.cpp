#include "LsFine.H"

#include <AMReX_ParallelContext.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <vector>

extern "C" {
void fds_ls_prepare(int r, int* nx, int* ny, double* hmax);
void fds_ls_ref_heights(int nx, int ny, double* h);
void fds_ls_wind_from_pairs(const double* zlo, const double* ulo, const double* vlo, const double* zhi, const double* uhi, const double* vhi, int* nbad);
}

namespace fdsamr {

LsFineConfig& ls_fine_config()
{
    static LsFineConfig c;
    return c;
}

int ls_fine_gather_wind(const LevelRegistry& reg)
{
    const int R = ls_fine_config().ratio;
    int nx = 0, ny = 0;
    double hmax = 0.0;
    fds_ls_prepare(R, &nx, &ny, &hmax);
    if (nx != reg.level(0).geom.Domain().length(0) * R || ny != reg.level(0).geom.Domain().length(1) * R) amrex::Abort("ls_fine_gather_wind: the LS grid does not cover the level-0 domain footprint");
    // column pre-reduction: per LS column only the two leaf samples that bracket the reference height H of the column are kept (the highest below H, the lowest at or
    // above H): LS_REF_WIND uses no other sample. Ranks: MAX of the lower height, MIN of the upper height, then the u, v of the rank that holds each (SUM with zeros, exact).
    const std::size_t n = static_cast<std::size_t>(nx) * ny;
    std::vector<double> href(n, -1.0);
    fds_ls_ref_heights(nx, ny, href.data());
    const double BIG = 1e300;
    std::vector<double> zlo(n, -BIG), zhi(n, BIG), ulo(n, 0.0), vlo(n, 0.0), uhi(n, 0.0), vhi(n, 0.0);
    int cum = 1;
    for (int l = 0; l < reg.num_levels(); ++l) {
        if (l > 0) cum *= reg.level(l).ref_ratio_from_parent[0];
        if (!reg.has_level(l)) continue;
        if (R % cum != 0) amrex::Abort("ls_fine_gather_wind: the LS ratio is not a multiple of the refinement ratio of level " + std::to_string(l));
        const int q = R / cum;
        const amrex::Geometry& g = reg.level(l).geom;
        const double dzl = g.CellSize(2);
        const amrex::Box dom = g.Domain();
        const amrex::MultiFab& U = reg.fields(l)["U"];
        const amrex::MultiFab& V = reg.fields(l)["V"];
        const amrex::iMultiFab* cov = reg.covered_mask(l);
        for (amrex::MFIter mfi(U); mfi.isValid(); ++mfi) {
            const amrex::Box cb = amrex::enclosedCells(mfi.validbox());   // cells of the box (U is nodal in x)
            auto ua = U.const_array(mfi);
            auto va = V.const_array(mfi);
            amrex::Array4<const int> ca;
            if (cov) ca = cov->const_array(mfi);
            for (int k = cb.smallEnd(2); k <= cb.bigEnd(2); ++k) {
                const double zc = (k - dom.smallEnd(2) + 0.5) * dzl;
                for (int j = cb.smallEnd(1); j <= cb.bigEnd(1); ++j)
                    for (int i = cb.smallEnd(0); i <= cb.bigEnd(0); ++i) {
                        if (cov && ca(i, j, k) != 0) continue;
                        const double uc = 0.5 * (ua(i, j, k) + ua(i + 1, j, k)), vc = 0.5 * (va(i, j, k) + va(i, j + 1, k));
                        for (int jj = (j - dom.smallEnd(1)) * q; jj < (j - dom.smallEnd(1) + 1) * q; ++jj)
                            for (int ii = (i - dom.smallEnd(0)) * q; ii < (i - dom.smallEnd(0) + 1) * q; ++ii) {
                                const std::size_t idx = static_cast<std::size_t>(jj) * nx + ii;
                                const double h = href[idx];
                                if (h < 0.0) continue;
                                if (zc < h) { if (zc > zlo[idx]) { zlo[idx] = zc; ulo[idx] = uc; vlo[idx] = vc; } }
                                else if (zc < zhi[idx]) { zhi[idx] = zc; uhi[idx] = uc; vhi[idx] = vc; }
                            }
                    }
            }
        }
    }
    if (amrex::ParallelDescriptor::NProcs() > 1) {
        const MPI_Comm comm = amrex::ParallelContext::CommunicatorSub();
        std::vector<double> zl = zlo, zh = zhi;
        amrex::ParallelAllReduce::Max(zlo.data(), static_cast<int>(n), comm);
        amrex::ParallelAllReduce::Min(zhi.data(), static_cast<int>(n), comm);
        for (std::size_t m = 0; m < n; ++m) {
            if (zl[m] != zlo[m]) { ulo[m] = 0.0; vlo[m] = 0.0; }
            if (zh[m] != zhi[m]) { uhi[m] = 0.0; vhi[m] = 0.0; }
        }
        for (std::vector<double>* a : {&ulo, &vlo, &uhi, &vhi}) amrex::ParallelAllReduce::Sum(a->data(), static_cast<int>(n), comm);
    }
    int nbad = 0;
    fds_ls_wind_from_pairs(zlo.data(), ulo.data(), vlo.data(), zhi.data(), uhi.data(), vhi.data(), &nbad);
    return nbad;
}

}  // namespace fdsamr
