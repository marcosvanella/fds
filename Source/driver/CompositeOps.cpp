// CompositeOps.cpp: see CompositeOps.H.
#include "../pressure_backend/FaceSync.H"
#include "CompositeOps.H"

#include <AMReX_BoxArray.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <memory>

#include "FaceTransfer.H"
#include "../pressure_backend/ObstMask.H"
#include "../pressure_backend/ObstCoarseFaces.H"

namespace fdsamr {

void composite_zero_solid_faces(LevelRegistry& reg, int l0, int nl, const std::array<const char*, 3>& names, int sh)
{
    for (int l = l0; l < nl; ++l) {
        const amrex::iMultiFab& m = reg.side_data(l).mask();
        for (int d = 0; d < 3; ++d) {
            amrex::MultiFab& f = reg.fields(l)[names[d]];
            const amrex::IntVect e = amrex::IntVect::TheDimensionVector(d);
            for (amrex::MFIter mfi(f); mfi.isValid(); ++mfi) {
                const auto a = f.array(mfi); const auto s = m.const_array(mfi);
                amrex::LoopOnCpu(amrex::surroundingNodes(m.box(mfi.index()), d), [&](int i, int j, int k) {
                    if (s(i, j, k) | s(i - e[0], j - e[1], k - e[2])) a(i + sh * e[0], j + sh * e[1], k + sh * e[2]) = 0.0;
                });
            }
        }
        if (l >= 1) {   // fine box faces next to a SOLID cell of the covering coarser level: wall, no normal velocity / flux (the U/V/W MultiFabs are nodal in d: coarse_solid_faces needs the CELL box, else the high-side faces in d are taken one cell too far out and stay open)
            const amrex::IntVect rr = reg.level(l).ref_ratio_from_parent;
            for (int d = 0; d < 3; ++d) {
                amrex::MultiFab& f = reg.fields(l)[names[d]];
                const amrex::IntVect e = amrex::IntVect::TheDimensionVector(d);
                for (amrex::MFIter mfi(f); mfi.isValid(); ++mfi) {
                    auto a = f.array(mfi);
                    for (const pb::ThinFace& cf : pb::coarse_solid_faces(reg.level(l - 1).geom, l - 1, rr, reg.level(l).ba, amrex::enclosedCells(mfi.validbox()), reg.level(l).geom))
                        if (cf.dir == d) a(cf.i + sh * e[0], cf.j + sh * e[1], cf.k + sh * e[2]) = 0.0;
                }
            }
        }
        if (l >= 1) {   // zero-thickness OBST faces (pb::thin_faces): no normal velocity / flux through them
            for (const pb::ThinFace& tf : pb::thin_faces(l)) {
                const int d = tf.dir; const amrex::IntVect e = amrex::IntVect::TheDimensionVector(d), p(tf.i, tf.j, tf.k);
                amrex::MultiFab& f = reg.fields(l)[names[d]];
                for (amrex::MFIter mfi(f); mfi.isValid(); ++mfi)
                    if (amrex::surroundingNodes(mfi.validbox(), d).contains(p)) f.array(mfi)(p[0] + sh * e[0], p[1] + sh * e[1], p[2] + sh * e[2]) = 0.0;
            }
        }
    }
}

void composite_zero_solid_cells(LevelRegistry& reg, int l0, int nl, const std::vector<const char*>& names)
{
    for (int l = l0; l < nl; ++l) {
        const amrex::iMultiFab& m = reg.side_data(l).mask();
        for (const char* nm : names) {
            amrex::MultiFab& f = reg.fields(l)[nm];
            for (amrex::MFIter mfi(f); mfi.isValid(); ++mfi) {
                const auto a = f.array(mfi); const auto s = m.const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { if (s(i, j, k)) a(i, j, k) = 0.0; });
            }
        }
    }
}

void composite_average_down_faces(LevelRegistry& reg, int nl, const std::array<const char*, 3>& names)
{
    for (int l = nl - 1; l >= 1; --l) {
        Fields& Ff = reg.fields(l);
        Fields& Fc = reg.fields(l - 1);
        const amrex::IntVect ratio = reg.level(l).ref_ratio_from_parent;
        for (int d = 0; d < 3; ++d) pb::sync_faces_low_owner(Ff[names[d]], reg.level(l).geom);   // one copy per shared fine face before the average-down
        std::array<const amrex::MultiFab*, 3> fi{&Ff[names[0]], &Ff[names[1]], &Ff[names[2]]};
        std::array<amrex::MultiFab*, 3> co{&Fc[names[0]], &Fc[names[1]], &Fc[names[2]]};
        fdsrt::average_down_faces(fi, co, ratio, reg.level(l - 1).geom.periodicity());
    }
}

void composite_average_down_flux(LevelRegistry& reg, int nl)
{
    // FVX/FVY/FVZ are cell-shaped arrays holding the flux at the UPPER face of the cell (FDS index): the nodal face n of AMReX is FV(n-1). Shift into nodal temporaries, area-mean fine -> coarse, shift back.
    const char* nm[3] = {"FVX", "FVY", "FVZ"};
    for (int l = nl - 1; l >= 1; --l) {
        Fields& Ff = reg.fields(l);
        Fields& Fc = reg.fields(l - 1);
        const amrex::IntVect ratio = reg.level(l).ref_ratio_from_parent;
        std::array<std::unique_ptr<amrex::MultiFab>, 3> nf, nc;
        std::array<const amrex::MultiFab*, 3> fi{nullptr, nullptr, nullptr};
        std::array<amrex::MultiFab*, 3> co{nullptr, nullptr, nullptr};
        for (int d = 0; d < 3; ++d) {
            const amrex::MultiFab& sf = Ff[nm[d]]; const amrex::MultiFab& sc = Fc[nm[d]];
            amrex::BoxArray bf = sf.boxArray(); bf.surroundingNodes(d);
            amrex::BoxArray bc = sc.boxArray(); bc.surroundingNodes(d);
            nf[d].reset(new amrex::MultiFab(bf, sf.DistributionMap(), 1, 0));
            nc[d].reset(new amrex::MultiFab(bc, sc.DistributionMap(), 1, 0));
            const amrex::IntVect sh = amrex::IntVect::TheDimensionVector(d);
            for (amrex::MFIter mfi(*nf[d]); mfi.isValid(); ++mfi) {
                auto a = nf[d]->array(mfi); const auto s = sf.const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { a(i, j, k) = s(i - sh[0], j - sh[1], k - sh[2]); });
            }
            for (amrex::MFIter mfi(*nc[d]); mfi.isValid(); ++mfi) {
                auto a = nc[d]->array(mfi); const auto s = sc.const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { a(i, j, k) = s(i - sh[0], j - sh[1], k - sh[2]); });
            }
            fi[d] = nf[d].get(); co[d] = nc[d].get();
        }
        for (int d = 0; d < 3; ++d) pb::sync_faces_low_owner(*nf[d], reg.level(l).geom, true);   // the computed copy: FV(n) of the box below
        fdsrt::average_down_faces(fi, co, ratio, reg.level(l - 1).geom.periodicity());
        for (int d = 0; d < 3; ++d) {
            amrex::MultiFab& dc = Fc[nm[d]];
            const amrex::IntVect sh = amrex::IntVect::TheDimensionVector(d);
            for (amrex::MFIter mfi(dc); mfi.isValid(); ++mfi) {
                auto a = dc.array(mfi); const auto s = nc[d]->const_array(mfi); const amrex::Box vb = mfi.validbox();
                amrex::LoopOnCpu(vb, [&](int i, int j, int k) { a(i, j, k) = s(i + sh[0], j + sh[1], k + sh[2]); });
            }
            dc.FillBoundary(reg.level(l - 1).geom.periodicity());   // the low face of a box is its ghost FV(-1): refresh it from the averaged copy (periodic image included)
        }
    }
}

void composite_fill_fine_face_ghosts(LevelRegistry& reg, int l, const std::array<const char*, 3>& names)
{
    const Level& lf = reg.level(l);
    const Level& lc = reg.level(l - 1);
    Fields& Ff = reg.fields(l);
    Fields& Fc = reg.fields(l - 1);
    const amrex::IntVect ratio = lf.ref_ratio_from_parent;
    const amrex::Box dom = lf.geom.Domain();
    int ng = 1;
    for (int d = 0; d < 3; ++d) ng = std::max(ng, Ff[names[d]].nGrow());
    // grown cell boxes, a union of whole coarse cells (grow by a multiple of the ratio in the refined directions, clipped to the domain in a non-periodic direction)
    amrex::BoxList bl;
    for (int i = 0; i < static_cast<int>(lf.ba.size()); ++i) {
        amrex::Box b = lf.ba[i];
        for (int d = 0; d < 3; ++d) {
            if (ratio[d] <= 1) continue;
            const int g = ((ng + ratio[d] - 1) / ratio[d]) * ratio[d];
            int lo = b.smallEnd(d) - g, hi = b.bigEnd(d) + g;
            if (!lf.geom.isPeriodic(d)) { lo = std::max(lo, dom.smallEnd(d)); hi = std::min(hi, dom.bigEnd(d)); }
            b.setSmall(d, lo); b.setBig(d, hi);
        }
        bl.push_back(b);
    }
    amrex::BoxArray gba(bl);
    std::array<std::unique_ptr<amrex::MultiFab>, 3> tmp;
    std::array<amrex::MultiFab*, 3> tp{nullptr, nullptr, nullptr};
    std::array<const amrex::MultiFab*, 3> cp{nullptr, nullptr, nullptr};
    for (int d = 0; d < 3; ++d) {
        amrex::BoxArray fb = gba; fb.surroundingNodes(d);
        tmp[d].reset(new amrex::MultiFab(fb, lf.dm, 1, 0));
        tmp[d]->setVal(0.0);
        tp[d] = tmp[d].get();
        cp[d] = &Fc[names[d]];
    }
    fdsrt::prolong_faces_level(tp, cp, lc.geom, ratio);
    for (int d = 0; d < 3; ++d) {
        amrex::MultiFab& dst = Ff[names[d]];
        for (amrex::MFIter mfi(dst); mfi.isValid(); ++mfi) {
            const amrex::Box vb = amrex::surroundingNodes(mfi.validbox(), d);
            auto a = dst.array(mfi);
            const auto s = tmp[d]->const_array(mfi);
            const amrex::Box x = (*tmp[d])[mfi].box() & dst[mfi].box();
            amrex::LoopOnCpu(x, [&](int i, int j, int k) { if (!vb.contains(i, j, k)) a(i, j, k) = s(i, j, k); });
        }
        Ff.fill_ghosts(names[d]);   // one source per seam ghost
    }
}

void composite_fine_pressure_rhs(Fields& F, const Level& lev, amrex::MultiFab& rhs)
{
    const double rdx = 1.0 / lev.dx[0], rdy = 1.0 / lev.dx[1], rdz = 1.0 / lev.dx[2];
    const bool y_term = lev.geom.Domain().length(1) > 1;   // TWO_D: PRESSURE_SOLVER_COMPUTE_RHS has no y term
    for (amrex::MFIter mfi(rhs); mfi.isValid(); ++mfi) {
        auto r = rhs.array(mfi);
        const auto fx = F["FVX"].const_array(mfi);
        const auto fy = F["FVY"].const_array(mfi);
        const auto fz = F["FVZ"].const_array(mfi);
        const auto dd = F["DDDT"].const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            // FVX(I) is the flux at the upper face of cell I (FDS index, cell-shaped array): (FVX(I-1)-FVX(I))/dx
            const double t1 = (fx(i - 1, j, k) - fx(i, j, k)) * rdx;
            const double t2 = y_term ? (fy(i, j - 1, k) - fy(i, j, k)) * rdy : 0.0;
            const double t3 = (fz(i, j, k - 1) - fz(i, j, k)) * rdz;
            r(i, j, k) = t1 + t2 + t3 - dd(i, j, k);
        });
    }
}

long composite_set_fine_h_ghosts(const Level& lev, amrex::MultiFab& H, const std::array<amrex::MultiFab*, 3>& grad)
{
    long n = 0;
    const amrex::Box dom = lev.geom.Domain();
    auto owned = [&](amrex::IntVect iv) {
        for (int d = 0; d < 3; ++d) {
            if (lev.geom.isPeriodic(d)) {
                const int len = dom.length(d);
                iv[d] = dom.smallEnd(d) + ((iv[d] - dom.smallEnd(d)) % len + len) % len;
            } else if (iv[d] < dom.smallEnd(d) || iv[d] > dom.bigEnd(d)) {
                return true;   // a physical edge: not a coarse-fine ghost
            }
        }
        return lev.ba.contains(iv);
    };
    for (amrex::MFIter mfi(H); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        auto h = H.array(mfi);
        for (int d = 0; d < 3; ++d) {
            if (dom.length(d) == 1) continue;
            const auto g = grad[d]->const_array(mfi);
            for (int side = 0; side < 2; ++side) {
                amrex::Box face = vb;   // the layer of valid cells next to the ghost layer
                const int in = side == 0 ? vb.smallEnd(d) : vb.bigEnd(d);
                face.setSmall(d, in); face.setBig(d, in);
                amrex::LoopOnCpu(face, [&](int i, int j, int k) {
                    amrex::IntVect gc(i, j, k);
                    gc[d] += side == 0 ? -1 : 1;
                    if (owned(gc)) return;
                    amrex::IntVect fc(i, j, k);
                    fc[d] += side == 0 ? 0 : 1;   // face index: low face of the inside cell (low side), low face of the ghost cell (high side)
                    const double gr = g(fc[0], fc[1], fc[2]);
                    const double hin = h(i, j, k);
                    h(gc[0], gc[1], gc[2]) = side == 0 ? hin - gr * lev.dx[d] : hin + gr * lev.dx[d];
                    ++n;
                });
            }
        }
    }
    return n;
}

}  // namespace fdsamr
