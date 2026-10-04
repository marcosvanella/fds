#include "CommonLayer.H"
#include "ExactSum.H"

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelContext.H>
#include <AMReX_iMultiFab.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

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


bool is_masked (PressureProblem const& p)
{
    if (!p.cell_class) { return false; }
    return p.cell_class->max(0) != 0 || p.cell_class->min(0) != 0;
}

iMultiFab class_ghost (PressureProblem const& p)
{
    iMultiFab cls(p.ba, p.dm, 1, 1);
    cls.setVal(CellSolid);                                   // outside a non-periodic domain face: nothing to couple to
    if (p.cell_class) { iMultiFab::Copy(cls, *p.cell_class, 0, 0, 1, 0); } else { cls.setVal(CellGas, 0, 1, 0); }
    cls.FillBoundary(p.geom.periodicity());
    return cls;
}

MultiFab known_ghost (PressureProblem const& p)
{
    MultiFab kv(p.ba, p.dm, 1, 1);
    kv.setVal(Real(0));
    if (p.known_value) { MultiFab::Copy(kv, *p.known_value, 0, 0, 1, 0); }
    kv.FillBoundary(p.geom.periodicity());
    return kv;
}

namespace {
// Masked components: flood fill over gas-gas faces (frozen/masked-notes.md). Collective.
ComponentMap label_components_masked (PressureProblem const& p)
{
    const Box dom = p.geom.Domain();
    const Long nx = dom.length(0), ny = dom.length(1);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dom.numPts() < Long(std::numeric_limits<int>::max()), "masked labelling needs fewer than 2^31 cells");
    const int big = std::numeric_limits<int>::max();
    iMultiFab cls = class_ghost(p);
    iMultiFab lab(p.ba, p.dm, 1, 1);
    lab.setVal(big);
    auto gidx = [&] (int i, int j, int k) { return int((i - dom.smallEnd(0)) + nx*((j - dom.smallEnd(1)) + ny*Long(k - dom.smallEnd(2)))); };
    for (MFIter mfi(lab); mfi.isValid(); ++mfi) {
        auto const& l = lab.array(mfi);
        auto const& c = cls.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { if (c(i,j,k) == CellGas) { l(i,j,k) = gidx(i,j,k); } });
    }
    for (;;) {
        lab.FillBoundary(p.geom.periodicity());
        Long changed = 0;
        for (MFIter mfi(lab); mfi.isValid(); ++mfi) {
            auto const& l = lab.array(mfi);
            auto const& c = cls.const_array(mfi);
            const Box vb = mfi.validbox();
            const Dim3 lo = lbound(vb), hi = ubound(vb);
            auto relax = [&] (int i, int j, int k) {
                if (c(i,j,k) != CellGas) { return false; }
                int m = l(i,j,k);
                const int cand[6] = {c(i-1,j,k) == CellGas ? l(i-1,j,k) : big, c(i+1,j,k) == CellGas ? l(i+1,j,k) : big,
                                     c(i,j-1,k) == CellGas ? l(i,j-1,k) : big, c(i,j+1,k) == CellGas ? l(i,j+1,k) : big,
                                     c(i,j,k-1) == CellGas ? l(i,j,k-1) : big, c(i,j,k+1) == CellGas ? l(i,j,k+1) : big};
                for (int q = 0; q < 6; ++q) { m = std::min(m, cand[q]); }
                if (m < l(i,j,k)) { l(i,j,k) = m; return true; }
                return false;
            };
            bool again = true;
            while (again) {                                     // alternate sweeps inside the box until stable
                again = false;
                for (int k = lo.z; k <= hi.z; ++k) for (int j = lo.y; j <= hi.y; ++j) for (int i = lo.x; i <= hi.x; ++i) {
                    if (relax(i,j,k)) { again = true; ++changed; }
                }
                for (int k = hi.z; k >= lo.z; --k) for (int j = hi.y; j >= lo.y; --j) for (int i = hi.x; i >= lo.x; --i) {
                    if (relax(i,j,k)) { again = true; ++changed; }
                }
            }
        }
        ParallelDescriptor::ReduceLongSum(changed);
        if (changed == 0) { break; }
    }
    // Roots (cells whose label is their own index) -> sorted unique list on every rank.
    std::vector<int> roots;
    for (MFIter mfi(lab); mfi.isValid(); ++mfi) {
        auto const& l = lab.const_array(mfi);
        auto const& c = cls.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { if (c(i,j,k) == CellGas && l(i,j,k) == gidx(i,j,k)) { roots.push_back(l(i,j,k)); } });
    }
    std::vector<int> all_roots = roots;
#ifdef AMREX_USE_MPI
    {
        const int np = ParallelContext::NProcsSub();
        std::vector<int> cnt(np, 0), displ(np, 0);
        int mine = int(roots.size());
        MPI_Allgather(&mine, 1, MPI_INT, cnt.data(), 1, MPI_INT, ParallelContext::CommunicatorSub());
        int tot = 0;
        for (int r = 0; r < np; ++r) { displ[r] = tot; tot += cnt[r]; }
        all_roots.assign(tot, 0);
        MPI_Allgatherv(roots.data(), mine, MPI_INT, all_roots.data(), cnt.data(), displ.data(), MPI_INT, ParallelContext::CommunicatorSub());
    }
#endif
    std::sort(all_roots.begin(), all_roots.end());
    all_roots.erase(std::unique(all_roots.begin(), all_roots.end()), all_roots.end());
    const int nc = int(all_roots.size());

    ComponentMap cm;
    cm.label.define(p.ba, p.dm, 1, 0);
    cm.label.setVal(-1);
    std::vector<Long> count(nc, 0);
    std::vector<int> open(nc, 0);
    std::array<BC,6> const ebc = effective_bc(p.bc, dom);
    for (MFIter mfi(lab); mfi.isValid(); ++mfi) {
        auto const& l = lab.const_array(mfi);
        auto const& c = cls.const_array(mfi);
        auto const& o = cm.label.array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            if (c(i,j,k) != CellGas) { return; }
            const int id = int(std::lower_bound(all_roots.begin(), all_roots.end(), l(i,j,k)) - all_roots.begin());
            o(i,j,k) = id;
            count[id] += 1;
            const int iv[3] = {i, j, k};
            for (int d = 0; d < 3; ++d) {
                int a[3] = {i, j, k}, b[3] = {i, j, k};
                a[d] -= 1; b[d] += 1;
                if (iv[d] == dom.smallEnd(d) && ebc[face_index(d,0)] == BC::Dirichlet) { open[id] = 1; }
                if (iv[d] == dom.bigEnd(d) && ebc[face_index(d,1)] == BC::Dirichlet) { open[id] = 1; }
                if (c(a[0],a[1],a[2]) == CellKnown || c(b[0],b[1],b[2]) == CellKnown) { open[id] = 1; }
            }
        });
    }
    if (nc > 0) {
        ParallelDescriptor::ReduceLongSum(count.data(), nc);
        ParallelDescriptor::ReduceIntMax(open.data(), nc);
    }
    for (int q = 0; q < nc; ++q) {
        ComponentInfo ci;
        ci.id = q;
        ci.singular = (open[q] == 0);
        ci.ncells = count[q];
        const int r = all_roots[q];
        ci.pin = IntVect(AMREX_D_DECL(int(r % nx) + dom.smallEnd(0), int((r / nx) % ny) + dom.smallEnd(1), int(r / (nx*ny)) + dom.smallEnd(2)));
        cm.comps.push_back(ci);
    }
    return cm;
}
}

ComponentMap label_components (PressureProblem const& p, bool record_pin)
{
    if (is_masked(p)) { return label_components_masked(p); }
    ComponentMap cm;
    cm.label.define(p.ba, p.dm, 1, 0);
    cm.label.setVal(0);
    ComponentInfo ci;
    ci.id = 0;
    bool open = false;
    std::array<BC,6> const ebc = effective_bc(p.bc, p.geom.Domain());
    for (int f = 0; f < 6; ++f) { open = open || (ebc[f] == BC::Dirichlet); }
    ci.singular = !open;
    ci.ncells = p.ba.numPts();

    if (!record_pin) { cm.comps.push_back(ci); return cm; }
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
struct MeanInfo { std::vector<double> mean; std::vector<double> rms; std::vector<Long> n; std::vector<double> sum; };

// Exact volume-weighted mean per component (sum(v*x) / (v*count)) and rms of x.
MeanInfo exact_mean (MultiFab const& mf, ComponentMap const& cm, iMultiFab const* uncovered, Real vol, bool want_rms = true)
{
    const int nc = static_cast<int>(cm.comps.size());
    ExactSumResult s = exact_sum(mf, 0, vol, nc, uncovered, &cm.label);
    ExactSumResult s2;
    s2.sum.assign(nc, 0.0); s2.count.assign(nc, 0);
    if (want_rms) {
        MultiFab sq(mf.boxArray(), mf.DistributionMap(), 1, 0);
        for (MFIter mfi(sq); mfi.isValid(); ++mfi) {
            auto const& q = sq.array(mfi);
            auto const& a = mf.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { q(i,j,k) = a(i,j,k)*a(i,j,k); });
        }
        s2 = exact_sum(sq, 0, 1.0, nc, uncovered, &cm.label);
    }
    MeanInfo m;
    for (int c = 0; c < nc; ++c) {
        const double n = static_cast<double>(s.count[c]);
        m.n.push_back(s.count[c]);
        m.sum.push_back(s.sum[c]);
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
                  MultiFab const* cell_volume, MeanKind kind, bool measure)
{
    const int nc = static_cast<int>(cm.comps.size());
    if (!cell_volume) {
        // Uniform cells: volume-weighted mean of b and arithmetic mean of vol*b differ only by the constant vol.
        MeanInfo m = exact_mean(rhs, cm, uncovered, vol, measure);
        const double bmax = rhs.norm0(0);
        const double floor_ = std::ldexp(bmax, -52);      // idempotence: below round-off of b itself, leave alone
        std::vector<double> shift(nc, 0.0);
        for (int c = 0; c < nc; ++c) {
            ComponentInfo& ci = cm.comps[c];
            if (!ci.singular) { continue; }
            // ScaledArithmetic literally: arithmetic mean of F = vol*b, divided by vol to give the shift of b.
            double mean = m.mean[c];
            if (kind == MeanKind::ScaledArithmetic && m.n[c] > 0) { mean = (m.sum[c] / double(m.n[c])) / vol; }
            if (std::abs(mean) > floor_) { shift[c] = mean; }
            // removed_mean: the constant subtracted from b (Volume) or the mean of F = vol*b (ScaledArithmetic).
            ci.removed_mean = (kind == MeanKind::ScaledArithmetic) ? shift[c] * vol : shift[c];
            ci.removed_rel = (m.rms[c] > 0.0) ? std::abs(mean) / m.rms[c] : 0.0;   // measured, also when not removed
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
    std::array<BC,6> bc = effective_bc(p.bc, dom);
    if (is_masked(p)) {
        // Masked operator (frozen/masked-notes.md): gas rows only; a face to a Solid cell carries nothing, a face to a Known cell
        // carries 2 (g - phi)/dx^2, a domain face as in the unmasked operator. Zero on Solid and Known cells.
        iMultiFab cls = class_ghost(p);
        MultiFab kv = known_ghost(p);
        for (MFIter mfi(out); mfi.isValid(); ++mfi) {
            auto const& a = phi.const_array(mfi);
            auto const& o = out.array(mfi);
            auto const& cl = cls.const_array(mfi);
            auto const& g = kv.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                if (cl(i,j,k) != CellGas) { o(i,j,k) = Real(0); return; }
                const int iv[3] = {i, j, k};
                const Real c = a(i,j,k);
                Real r = 0.0;
                for (int d = 0; d < 3; ++d) {
                    for (int side = 0; side < 2; ++side) {
                        int nb[3] = {i, j, k};
                        nb[d] += (side == 0) ? -1 : 1;
                        const bool edge = (side == 0) ? (iv[d] == dom.smallEnd(d)) : (iv[d] == dom.bigEnd(d));
                        if (edge && bc[face_index(d,side)] != BC::Periodic) {
                            if (bc[face_index(d,side)] == BC::Dirichlet) { r += Real(-2.0)*c*idx2[d]; }
                            continue;
                        }
                        const int kc = cl(nb[0],nb[1],nb[2]);
                        if (kc == CellGas) { r += (a(nb[0],nb[1],nb[2]) - c) * idx2[d]; }
                        else if (kc == CellKnown) { r += Real(2.0) * (g(nb[0],nb[1],nb[2]) - c) * idx2[d]; }
                    }
                }
                o(i,j,k) = r;
            });
        }
        return;
    }
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

Status fold_boundary_data (PressureProblem const& p, BoundaryData const& bd, MultiFab& rhs, std::string* message)
{
    auto fail = [&] (Status st, std::string const& m) { if (message) { *message = m; } return st; };
    if (!p.levels.empty() || p.nlevels != 1) { return fail(Status::NotBuilt, "fold_boundary_data is built for single-level problems"); }
    if (rhs.boxArray() != p.ba || !(rhs.DistributionMap() == p.dm)) { return fail(Status::InvalidInput, "rhs BoxArray/DistributionMapping differ from the problem"); }
    const Box dom = p.geom.Domain();
    for (int d = 0; d < 3; ++d) {
        for (int side = 0; side < 2; ++side) {
            const int f = face_index(d, side);
            const bool has = (bd.value[f] != nullptr) || (bd.constant[f] != Real(0));
            if (p.bc[f] == BC::Periodic && has) { return fail(Status::InvalidInput, "boundary data given on a periodic face"); }
            if (bd.value[f] && bd.value[f]->nComp() < 1) { return fail(Status::InvalidInput, "boundary data MultiFab has no component"); }
            if (p.bc[f] == BC::Dirichlet && dom.length(d) == 1 && d != 1) { return fail(Status::NotBuilt, "a Dirichlet face in a one-cell x or z direction is not built"); }
        }
    }
    for (int d = 0; d < 3; ++d) {
        if (d == 1 && dom.length(1) == 1) { continue; }               // TWO_D: the y term does not exist
        const Real h = p.geom.CellSize(d);
        for (int side = 0; side < 2; ++side) {
            const int f = face_index(d, side);
            if (p.bc[f] == BC::Periodic) { continue; }
            if (!bd.value[f] && bd.constant[f] == Real(0)) { continue; }
            Real coef;                                                 // factor applied to the data value
            if (p.bc[f] == BC::Neumann) { coef = (side == 0 ? Real(1) : Real(-1)) / h; }
            else { coef = Real(-2) / (h*h); }
            Box layer = dom;
            if (side == 0) { layer.setBig(d, dom.smallEnd(d)); } else { layer.setSmall(d, dom.bigEnd(d)); }
            MultiFab tmp(p.ba, p.dm, 1, 0);
            tmp.setVal(Real(0));
            if (bd.value[f]) { tmp.ParallelCopy(*bd.value[f], 0, 0, 1); }
            for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
                const Box b = mfi.validbox() & layer;
                if (!b.ok()) { continue; }
                auto const& r = rhs.array(mfi);
                auto const& t = tmp.const_array(mfi);
                const Real cst = bd.value[f] ? Real(0) : bd.constant[f];
                const bool use_arr = (bd.value[f] != nullptr);
                amrex::LoopOnCpu(b, [&] (int i, int j, int k) { r(i,j,k) += coef * (use_arr ? t(i,j,k) : cst); });
            }
        }
    }
    if (message) { message->clear(); }
    return Status::Ok;
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

// One residual cell: sum over ranks of the value at iv (the cell belongs to exactly one box); optionally zeroes it.
static double take_cell (MultiFab& mf, IntVect const& iv, bool zero)
{
    double v = 0.0;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        if (mfi.validbox().contains(iv)) {
            auto const& a = mf.array(mfi);
            v += double(a(iv[0], iv[1], iv[2]));
            if (zero) { a(iv[0], iv[1], iv[2]) = Real(0); }
        }
    }
    ParallelDescriptor::ReduceRealSum(v);
    return v;
}

ResidualSums residual_sums (PressureProblem const& p, MultiFab& phi, MultiFab const& rhs, bool singular, IntVect const* pin)
{
    ResidualSums S;
    S.singular = singular;
    MultiFab r(p.ba, p.dm, 1, 0);
    apply_operator(p, phi, r);
    MultiFab::Xpay(r, Real(-1.0), rhs, 0, 0, 1, 0);          // r = rhs - L phi
    const double rn = r.norm2(0);
    S.r2 = rn*rn;
    S.rmax = r.norm0(0);
    S.sumr = exact_sum(r, 0, 1.0, 1, nullptr, nullptr).sum[0];
    S.sumw = double(p.ba.numPts());
    if (pin) {
        const double rp = take_cell(r, *pin, true);
        S.pin2 = rp*rp; S.pin_abs = std::abs(rp);
        const double rn2 = r.norm2(0);
        S.r2_nopin = rn2*rn2;
        S.rmax_nopin = r.norm0(0);
        S.pin_excluded = true;
    } else {
        S.r2_nopin = S.r2; S.rmax_nopin = S.rmax;
    }
    double pn = phi.norm2(0);
    if (is_masked(p)) {
        // Norms over gas cells only: r and rhs are zero on the other cells (the common layer zeroes b there); H is Known-valued/0 there.
        iMultiFab cls = class_ghost(p);
        MultiFab ph(p.ba, p.dm, 1, 0);
        MultiFab::Copy(ph, phi, 0, 0, 1, 0);
        Long ngas = 0;
        for (MFIter mfi(ph); mfi.isValid(); ++mfi) {
            auto const& a = ph.array(mfi);
            auto const& cl = cls.const_array(mfi);
            Long n = 0;
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { if (cl(i,j,k) != CellGas) { a(i,j,k) = Real(0); } else { ++n; } });
            ngas += n;
        }
        ParallelDescriptor::ReduceLongSum(ngas);
        pn = ph.norm2(0);
        S.sumw = double(ngas);
    }
    const double bn = rhs.norm2(0);
    S.b2 = bn*bn; S.phi2 = pn*pn; S.bmax = rhs.norm0(0);
    const auto dx = p.geom.CellSizeArray();
    const Box dom = p.geom.Domain();
    double a = 0.0;
    for (int d = 0; d < 3; ++d) { if (dom.length(d) > 1) { a += 4.0/(double(dx[d])*double(dx[d])); } }
    S.anorm = a;
    return S;
}

void evaluate_residual (PressureResult& R, PressureOptions const& o, ResidualSums const& S)
{
    auto rel = [&] (double sq) { return (S.b2 > 0.0) ? std::sqrt(sq / S.b2) : std::sqrt(sq); };
    R.residual_checked = true;
    R.residual_rel2 = rel(S.r2);
    R.residual_relmax = (S.bmax > 0.0) ? S.rmax / S.bmax : S.rmax;
    // Mean-removed: sum(w (r - m)^2) = sum(w r^2) - m^2 sum(w) with m the weighted mean (exact sum).
    const double m = (S.sumw > 0.0) ? S.sumr / S.sumw : 0.0;
    R.residual_rel2_mr = S.singular ? rel(std::max(0.0, S.r2 - m*m*S.sumw)) : R.residual_rel2;
    R.residual_rel2_nopin = rel(S.r2_nopin);
    R.residual_relmax_nopin = (S.bmax > 0.0) ? S.rmax_nopin / S.bmax : S.rmax_nopin;
    R.residual_pin_abs = S.pin_abs;
    R.residual_pin_rel = rel(S.pin2);
    const double bn = std::sqrt(S.b2), an_phi = S.anorm * std::sqrt(S.phi2);
    R.residual_backward = std::sqrt(S.r2_nopin) / (bn + an_phi > 0.0 ? bn + an_phi : 1.0);
    const double u = std::ldexp(1.0, -53);
    // Size-scaled limit (FR-031, docs/pressure/07 sections 14-16): the residual of a floating-point evaluation cannot be below
    // c u ||A|| ||H||_2 / ||b||_2 (c = kResidualRoundoff = 10), so the limit is max(residual_tol, that floor), for singular and
    // non-singular components alike. The compared value is the non-pin residual (equal to the full one when no pin was applied).
    R.residual_floor = kResidualRoundoff * u * an_phi / (bn > 0.0 ? bn : 1.0);
    R.residual_check = S.singular ? R.residual_rel2_nopin : R.residual_rel2;
    R.residual_limit = std::max(o.residual_tol, R.residual_floor);
    R.residual_ok = (R.residual_check <= R.residual_limit);
    if (!R.residual_ok) {
        std::ostringstream m2;
        if (S.singular) {
            m2 << "true " << S.what << "residual ||b-L*H||_2/||b||_2 = " << R.residual_check << " (pin cell excluded; raw "
               << R.residual_rel2 << ", pin row " << R.residual_pin_rel << ") exceeds " << R.residual_limit << " (residual_tol " << o.residual_tol
               << ", round-off floor " << R.residual_floor << ")";
        } else {
            m2 << "true " << S.what << "residual ||b-L*H||_2/||b||_2 = " << R.residual_rel2 << " exceeds " << R.residual_limit << " (residual_tol " << o.residual_tol
               << ", round-off floor " << R.residual_floor << ")";
        }
        R.warnings.push_back(m2.str());
    }
}

} // namespace pb
