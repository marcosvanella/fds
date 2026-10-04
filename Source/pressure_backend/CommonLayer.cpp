#include "CommonLayer.H"
#include "ExactSum.H"

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <cmath>
#include <limits>

namespace pb {

using namespace amrex;

const char* to_string (Status s)
{
    switch (s) {
    case Status::Ok: return "Ok";
    case Status::NotBuilt: return "NotBuilt";
    case Status::InvalidInput: return "InvalidInput";
    case Status::NotConverged: return "NotConverged";
    }
    return "?";
}

const char* to_string (BC b)
{
    switch (b) {
    case BC::Neumann: return "neumann";
    case BC::Periodic: return "periodic";
    case BC::Dirichlet: return "dirichlet";
    }
    return "?";
}

ComponentMap label_components (PressureProblem const& p)
{
    ComponentMap cm;
    cm.label.define(p.ba, p.dm, 1, 0);
    cm.label.setVal(0);
    ComponentInfo ci;
    ci.id = 0;
    bool open = false;
    for (int f = 0; f < 6; ++f) { open = open || (p.bc[f] == BC::Dirichlet); }
    ci.singular = !open;
    ci.ncells = p.ba.numPts();

    // Pin: lowest global index (x fastest) of the component, from the labelled cells.
    const Box dom = p.geom.Domain();
    const Long nx = dom.length(0), ny = dom.length(1);
    Long lowest = std::numeric_limits<Long>::max();
    for (MFIter mfi(cm.label); mfi.isValid(); ++mfi) {
        auto const& l = cm.label.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            if (l(i,j,k) == 0) {
                const Long idx = (i - dom.smallEnd(0)) + nx*((j - dom.smallEnd(1)) + ny*Long(k - dom.smallEnd(2)));
                lowest = std::min(lowest, idx);
            }
        });
    }
    ParallelDescriptor::ReduceLongMin(lowest);
    ci.pin = IntVect(AMREX_D_DECL(int(lowest % nx) + dom.smallEnd(0),
                                  int((lowest / nx) % ny) + dom.smallEnd(1),
                                  int(lowest / (nx*ny)) + dom.smallEnd(2)));
    cm.comps.push_back(ci);
    return cm;
}

namespace {
struct MeanInfo { std::vector<double> mean; std::vector<double> rms; std::vector<Long> n; };

// Exact volume-weighted mean per component (sum(v*x) / (v*count)) and rms of x.
MeanInfo exact_mean (MultiFab const& mf, ComponentMap const& cm, iMultiFab const* uncovered, Real vol)
{
    const int nc = static_cast<int>(cm.comps.size());
    ExactSumResult s = exact_sum(mf, 0, vol, nc, uncovered, &cm.label);
    MultiFab sq(mf.boxArray(), mf.DistributionMap(), 1, 0);
    for (MFIter mfi(sq); mfi.isValid(); ++mfi) {
        auto const& q = sq.array(mfi);
        auto const& a = mf.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { q(i,j,k) = a(i,j,k)*a(i,j,k); });
    }
    ExactSumResult s2 = exact_sum(sq, 0, 1.0, nc, uncovered, &cm.label);
    MeanInfo m;
    for (int c = 0; c < nc; ++c) {
        const double n = static_cast<double>(s.count[c]);
        m.n.push_back(s.count[c]);
        m.mean.push_back(n > 0 ? s.sum[c] / (vol*n) : 0.0);
        m.rms.push_back(n > 0 ? std::sqrt(s2.sum[c] / n) : 0.0);
    }
    return m;
}

void subtract_per_component (MultiFab& mf, ComponentMap const& cm, iMultiFab const* uncovered,
                             std::vector<double> const& shift)
{
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        auto const& a = mf.array(mfi);
        auto const& l = cm.label.const_array(mfi);
        Array4<int const> u;
        if (uncovered) { u = uncovered->const_array(mfi); }
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const int c = l(i,j,k);
            if (c < 0 || c >= static_cast<int>(shift.size())) { return; }
            if (u.contains(i,j,k) && u(i,j,k) == 0) { return; }
            a(i,j,k) -= shift[c];
        });
    }
}
}

void remove_mean (MultiFab& rhs, ComponentMap& cm, iMultiFab const* uncovered, Real vol,
                  MultiFab const* cell_volume, MeanKind kind)
{
    const int nc = static_cast<int>(cm.comps.size());
    if (!cell_volume) {
        // Uniform cells: volume-weighted mean of b and arithmetic mean of vol*b differ only by the constant vol.
        MeanInfo m = exact_mean(rhs, cm, uncovered, vol);
        const double bmax = rhs.norm0(0);
        const double floor_ = std::ldexp(bmax, -52);      // idempotence: below round-off of b itself, leave alone
        std::vector<double> shift(nc, 0.0);
        for (int c = 0; c < nc; ++c) {
            ComponentInfo& ci = cm.comps[c];
            if (!ci.singular) { continue; }
            if (std::abs(m.mean[c]) > floor_) { shift[c] = m.mean[c]; }
            ci.removed_mean = shift[c];
            ci.removed_rel = (m.rms[c] > 0.0) ? std::abs(m.mean[c]) / m.rms[c] : 0.0;   // measured, also when not removed
        }
        subtract_per_component(rhs, cm, uncovered, shift);
        return;
    }
    // Per-cell volumes v. F = v*b is the volume-scaled right-hand side of the finite-volume system.
    MultiFab F(rhs.boxArray(), rhs.DistributionMap(), 1, 0);
    MultiFab::Copy(F, rhs, 0, 0, 1, 0);
    MultiFab::Multiply(F, *cell_volume, 0, 0, 1, 0);
    ExactSumResult sF = exact_sum(F, 0, 1.0, nc, uncovered, &cm.label);
    std::vector<double> mean(nc, 0.0), rms(nc, 0.0);
    double scale = 0.0;                                   // largest |term| of the quantity whose mean is removed
    if (kind == MeanKind::Volume) {
        ExactSumResult sV = exact_sum(*cell_volume, 0, 1.0, nc, uncovered, &cm.label);
        MultiFab sq(rhs.boxArray(), rhs.DistributionMap(), 1, 0);
        MultiFab::Copy(sq, rhs, 0, 0, 1, 0); MultiFab::Multiply(sq, rhs, 0, 0, 1, 0);
        ExactSumResult s2 = exact_sum(sq, 0, 1.0, nc, uncovered, &cm.label);
        for (int c = 0; c < nc; ++c) {
            mean[c] = (sV.sum[c] > 0.0) ? sF.sum[c] / sV.sum[c] : 0.0;
            rms[c] = (sF.count[c] > 0) ? std::sqrt(s2.sum[c] / double(sF.count[c])) : 0.0;
        }
        scale = rhs.norm0(0);
    } else {
        MultiFab sq(rhs.boxArray(), rhs.DistributionMap(), 1, 0);
        MultiFab::Copy(sq, F, 0, 0, 1, 0); MultiFab::Multiply(sq, F, 0, 0, 1, 0);
        ExactSumResult s2 = exact_sum(sq, 0, 1.0, nc, uncovered, &cm.label);
        for (int c = 0; c < nc; ++c) {
            const double n = double(sF.count[c]);
            mean[c] = (n > 0) ? sF.sum[c] / n : 0.0;
            rms[c] = (n > 0) ? std::sqrt(s2.sum[c] / n) : 0.0;
        }
        scale = F.norm0(0);
    }
    const double floor_ = std::ldexp(scale, -52);
    std::vector<double> shift(nc, 0.0);
    for (int c = 0; c < nc; ++c) {
        ComponentInfo& ci = cm.comps[c];
        if (!ci.singular) { continue; }
        if (std::abs(mean[c]) > floor_) { shift[c] = mean[c]; }
        ci.removed_mean = shift[c];
        ci.removed_rel = (rms[c] > 0.0) ? std::abs(mean[c]) / rms[c] : 0.0;
    }
    if (kind == MeanKind::Volume) {
        subtract_per_component(rhs, cm, uncovered, shift);
    } else {
        for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {      // b_k -= mean(F)/v_k, so that sum(v*b) is zero
            auto const& a = rhs.array(mfi);
            auto const& l = cm.label.const_array(mfi);
            auto const& v = cell_volume->const_array(mfi);
            Array4<int const> u;
            if (uncovered) { u = uncovered->const_array(mfi); }
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const int c = l(i,j,k);
                if (c < 0 || c >= nc || shift[c] == 0.0) { return; }
                if (u.contains(i,j,k) && u(i,j,k) == 0) { return; }
                a(i,j,k) -= shift[c] / v(i,j,k);
            });
        }
    }
}

void apply_gauge (MultiFab& phi, ComponentMap& cm, iMultiFab const* uncovered, Real vol,
                  MultiFab const* cell_volume, MultiFab const* gauge_weight, MultiFab const* gauge_offset)
{
    const int nc = static_cast<int>(cm.comps.size());
    std::vector<double> shift(nc, 0.0);
    if (!cell_volume && !gauge_weight && !gauge_offset) {
        MeanInfo m = exact_mean(phi, cm, uncovered, vol);
        for (int c = 0; c < nc; ++c) {
            if (cm.comps[c].singular) { shift[c] = m.mean[c]; cm.comps[c].gauge_shift = shift[c]; }
        }
    } else {
        // shift = sum(W*(phi - g)) / sum(W), W = v*rho (v: cell volume, rho: gauge_weight), g: gauge_offset.
        MultiFab W(phi.boxArray(), phi.DistributionMap(), 1, 0);
        if (cell_volume) { MultiFab::Copy(W, *cell_volume, 0, 0, 1, 0); } else { W.setVal(vol); }
        if (gauge_weight) { MultiFab::Multiply(W, *gauge_weight, 0, 0, 1, 0); }
        MultiFab X(phi.boxArray(), phi.DistributionMap(), 1, 0);
        MultiFab::Copy(X, phi, 0, 0, 1, 0);
        if (gauge_offset) { MultiFab::Subtract(X, *gauge_offset, 0, 0, 1, 0); }
        ExactSumResult sx = exact_sum(X, 0, 1.0, nc, uncovered, &cm.label, &W);
        ExactSumResult sw = exact_sum(W, 0, 1.0, nc, uncovered, &cm.label);
        for (int c = 0; c < nc; ++c) {
            if (cm.comps[c].singular && sw.sum[c] > 0.0) { shift[c] = sx.sum[c] / sw.sum[c]; cm.comps[c].gauge_shift = shift[c]; }
        }
    }
    subtract_per_component(phi, cm, uncovered, shift);
}

void apply_operator (PressureProblem const& p, MultiFab& phi, MultiFab& out)
{
    AMREX_ALWAYS_ASSERT(phi.nGrow() >= 1);
    phi.FillBoundary(p.geom.periodicity());
    const Box dom = p.geom.Domain();
    const auto dx = p.geom.CellSizeArray();
    Real idx2[3] = {1.0/(dx[0]*dx[0]), 1.0/(dx[1]*dx[1]), 1.0/(dx[2]*dx[2])};
    std::array<BC,6> bc = p.bc;
    for (MFIter mfi(out); mfi.isValid(); ++mfi) {
        auto const& a = phi.const_array(mfi);
        auto const& o = out.array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const int iv[3] = {i, j, k};
            const Real c = a(i,j,k);
            Real r = 0.0;
            for (int d = 0; d < 3; ++d) {
                int lo[3] = {i, j, k}, hi[3] = {i, j, k};
                lo[d] -= 1; hi[d] += 1;
                Real vlo = a(lo[0],lo[1],lo[2]);   // ghost (periodic) or interior neighbour
                Real vhi = a(hi[0],hi[1],hi[2]);
                if (iv[d] == dom.smallEnd(d) && bc[face_index(d,0)] != BC::Periodic) {
                    vlo = (bc[face_index(d,0)] == BC::Neumann) ? c : -c;
                }
                if (iv[d] == dom.bigEnd(d) && bc[face_index(d,1)] != BC::Periodic) {
                    vhi = (bc[face_index(d,1)] == BC::Neumann) ? c : -c;
                }
                r += (vlo + vhi - Real(2.0)*c) * idx2[d];
            }
            o(i,j,k) = r;
        });
    }
}

ResidualNorms true_residual (PressureProblem const& p, MultiFab& phi, MultiFab const& rhs)
{
    MultiFab lphi(p.ba, p.dm, 1, 0);
    apply_operator(p, phi, lphi);
    MultiFab::Xpay(lphi, Real(-1.0), rhs, 0, 0, 1, 0);   // lphi = rhs - lphi
    ResidualNorms n;
    const double b2 = rhs.norm2(0), bm = rhs.norm0(0);
    n.rel2 = (b2 > 0.0) ? lphi.norm2(0) / b2 : lphi.norm2(0);
    n.relmax = (bm > 0.0) ? lphi.norm0(0) / bm : lphi.norm0(0);
    return n;
}

} // namespace pb
