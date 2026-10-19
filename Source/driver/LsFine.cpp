#include "LsFine.H"

#include <AMReX_ParallelContext.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <vector>

extern "C" {
void fds_ls_prepare(int r, int* nx, int* ny, double* hmax);
void fds_ls_wind_from_slots(int nks, const double* zs, const double* us, const double* vs, const double* cnt, int* nbad);
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
    const double dz0 = reg.level(0).geom.CellSize(2);
    const double dzf = dz0 / R;
    const int nzf = reg.level(0).geom.Domain().length(2) * R;
    // slots: the leaf cells whose lower face is below hmax + dz0 can hold the two samples that bracket the reference height
    const int nks = std::min(nzf, static_cast<int>(std::ceil((hmax + dz0) / dzf - 1e-9)) + 1);
    const std::size_t n = static_cast<std::size_t>(nx) * ny * nks;
    std::vector<double> z(n, 0.0), u(n, 0.0), v(n, 0.0), cnt(n, 0.0);
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
            const amrex::Box vb = mfi.validbox();   // cells of the box (nodal MultiFab: validbox is the cell box grown to faces; convert)
            const amrex::Box cb = amrex::enclosedCells(vb);
            auto ua = U.const_array(mfi);
            auto va = V.const_array(mfi);
            amrex::Array4<const int> ca;
            if (cov) ca = cov->const_array(mfi);
            for (int k = cb.smallEnd(2); k <= cb.bigEnd(2); ++k) {
                const int ks = (k - dom.smallEnd(2)) * q;
                if (ks >= nks) continue;
                const double zc = (k - dom.smallEnd(2) + 0.5) * dzl;
                for (int j = cb.smallEnd(1); j <= cb.bigEnd(1); ++j)
                    for (int i = cb.smallEnd(0); i <= cb.bigEnd(0); ++i) {
                        if (cov && ca(i, j, k) != 0) continue;
                        const double uc = 0.5 * (ua(i, j, k) + ua(i + 1, j, k)), vc = 0.5 * (va(i, j, k) + va(i, j + 1, k));
                        for (int jj = (j - dom.smallEnd(1)) * q; jj < (j - dom.smallEnd(1) + 1) * q; ++jj)
                            for (int ii = (i - dom.smallEnd(0)) * q; ii < (i - dom.smallEnd(0) + 1) * q; ++ii) {
                                const std::size_t idx = (static_cast<std::size_t>(jj) * nx + ii) * nks + ks;
                                z[idx] = zc; u[idx] = uc; v[idx] = vc; cnt[idx] = 1.0;
                            }
                    }
            }
        }
    }
    if (amrex::ParallelDescriptor::NProcs() > 1) {
        const MPI_Comm comm = amrex::ParallelContext::CommunicatorSub();
        for (std::vector<double>* a : {&z, &u, &v, &cnt}) amrex::ParallelAllReduce::Sum(a->data(), static_cast<int>(n), comm);
    }
    int nbad = 0;
    fds_ls_wind_from_slots(nks, z.data(), u.data(), v.data(), cnt.data(), &nbad);
    return nbad;
}

}  // namespace fdsamr
