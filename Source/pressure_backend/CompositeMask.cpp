// Composite masked labelling (D4) and Known-on-C/F refusal (Q2).
#include "CommonLayer.H"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <limits>
#include <set>

namespace pb {
using namespace amrex;

bool composite_is_masked (PressureProblem const& p)
{
    for (auto const& L : p.levels) {
        if (L.cell_class && (L.cell_class->max(0) != 0 || L.cell_class->min(0) != 0)) { return true; }
    }
    return false;
}

namespace {

Long gindex (Box const& dom, IntVect const& iv, Long off)
{
    const Long nx = dom.length(0), ny = dom.length(1);
    return off + Long(iv[0]-dom.smallEnd(0)) + nx*(Long(iv[1]-dom.smallEnd(1)) + ny*Long(iv[2]-dom.smallEnd(2)));
}

void copy_class (PressureLevel const& L, Geometry const& geom, iMultiFab& cls)
{
    cls.define(L.ba, L.dm, 1, 1);
    cls.setVal(CellGas);                  // C/F ghost cells are not filled below; never leave them uninitialised
    if (L.cell_class) {
        for (MFIter mfi(cls); mfi.isValid(); ++mfi) {
            auto const& o = cls.array(mfi);
            auto const& c = L.cell_class->const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { o(i,j,k) = c(i,j,k); });
        }
    } else {
        cls.setVal(CellGas);
    }
    cls.FillBoundary(geom.periodicity());
    Box const& dom = geom.Domain();
    for (MFIter mfi(cls); mfi.isValid(); ++mfi) {
        Box gbx = amrex::grow(mfi.validbox(), 1);
        auto const& o = cls.array(mfi);
        amrex::LoopOnCpu(gbx, [&] (int i, int j, int k) {
            IntVect iv(i,j,k);
            if (dom.contains(iv)) { return; }
            bool only_per = true;
            for (int d = 0; d < 3; ++d) {
                if (iv[d] < dom.smallEnd(d) || iv[d] > dom.bigEnd(d)) {
                    if (!geom.isPeriodic(d)) { only_per = false; }
                }
            }
            if (!only_per) { o(i,j,k) = CellSolid; }
        });
    }
}

bool wrap_iv (Geometry const& geom, IntVect& iv)
{
    Box const& dom = geom.Domain();
    for (int d = 0; d < 3; ++d) {
        if (iv[d] < dom.smallEnd(d) || iv[d] > dom.bigEnd(d)) {
            if (!geom.isPeriodic(d)) { return false; }
            const int len = dom.length(d);
            iv[d] = dom.smallEnd(d) + (((iv[d] - dom.smallEnd(d)) % len) + len) % len;
        }
    }
    return true;
}

} // namespace

std::string label_components_composite_masked (PressureProblem const& p,
    std::vector<std::unique_ptr<iMultiFab>> const& W_unc,
    std::vector<Geometry> const& geom,
    std::vector<IntVect> const& ratio,
    CompositeComponentMap& cm)
{
    const double t0 = amrex::second();
    MFIter::allowMultipleMFIters(true);
    const int nlev = int(p.levels.size());
    cm = CompositeComponentMap{};
    cm.label.resize(nlev);
    std::vector<iMultiFab> cls(nlev);
    std::array<BC,6> const ebc = effective_bc(p.bc, geom[0].Domain());

    Long nall = 0;
    for (int l = 0; l < nlev; ++l) { nall += geom[l].Domain().numPts(); }
    if (nall >= Long(std::numeric_limits<int>::max())) {
        return "masked composite labelling needs fewer than 2^31 cells across levels";
    }
    std::vector<Long> off(nlev, 0);
    for (int l = 1; l < nlev; ++l) { off[l] = off[l-1] + geom[l-1].Domain().numPts(); }

    for (int l = 0; l < nlev; ++l) {
        copy_class(p.levels[l], geom[l], cls[l]);
        cm.label[l] = std::make_unique<iMultiFab>(p.levels[l].ba, p.levels[l].dm, 1, 1);
        cm.label[l]->setVal(-1);
    }

    // Q2: Known on C/F.
    int known_cf = 0;
    for (int L = 1; L < nlev; ++L) {
        BoxArray cba = p.levels[L].ba; cba.coarsen(ratio[L]); cba.grow(1);
        iMultiFab ccls(cba, p.levels[L].dm, 1, 0);
        ccls.setVal(CellSolid);
        ccls.ParallelCopy(cls[L-1], 0, 0, 1);
        BoxArray const& fba = p.levels[L].ba;
        for (MFIter mfi(cls[L]); mfi.isValid(); ++mfi) {
            auto const& fl = cls[L].const_array(mfi);
            // ccls shares the fine DM: same MFIter index covers the coarsened+grown box.
            auto const& ca = ccls.const_array(mfi);
            Box const cbox = ccls[mfi].box();
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                for (int d = 0; d < 3; ++d) for (int s = 0; s < 2; ++s) {
                    IntVect n(i,j,k); n[d] += s ? 1 : -1;
                    if (fba.contains(n)) { continue; }
                    IntVect w = n;
                    if (!wrap_iv(geom[L], w)) { continue; }
                    if (fba.contains(w)) { continue; }
                    IntVect cc;
                    for (int tt = 0; tt < 3; ++tt) { cc[tt] = amrex::coarsen(w[tt], ratio[L][tt]); }
                    int kc = CellSolid;
                    if (cbox.contains(cc)) { kc = ca(cc[0],cc[1],cc[2]); }
                    if (fl(i,j,k) == CellKnown || kc == CellKnown) { known_cf = 1; }
                }
            });
        }
    }
    ParallelDescriptor::ReduceIntMax(known_cf);
    if (known_cf) {
        return "Known (OPEN) cell on a coarse/fine face is not built (Q2: vents stay inside one level; D3)";
    }

    // Seed.
    for (int l = 0; l < nlev; ++l) {
        Box const& dom = geom[l].Domain();
        for (MFIter mfi(*cm.label[l]); mfi.isValid(); ++mfi) {
            auto const& o = cm.label[l]->array(mfi);
            auto const& cl = cls[l].const_array(mfi);
            auto const& u = W_unc[l]->const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                if (u(i,j,k) != 0 && cl(i,j,k) == CellGas) {
                    o(i,j,k) = int(gindex(dom, IntVect(i,j,k), off[l]));
                }
            });
        }
        cm.label[l]->FillBoundary(geom[l].periodicity());
    }

    auto same_level_sweep = [&] (int l) {
        int ch = 0;
        Box const& dom = geom[l].Domain();
        for (MFIter mfi(*cm.label[l]); mfi.isValid(); ++mfi) {
            auto const& o = cm.label[l]->array(mfi);
            auto const& cl = cls[l].const_array(mfi);
            auto const& u = W_unc[l]->const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                if (u(i,j,k) == 0 || cl(i,j,k) != CellGas || o(i,j,k) < 0) { return; }
                int cur = o(i,j,k);
                const int iv[3] = {i,j,k};
                for (int d = 0; d < 3; ++d) for (int s = 0; s < 2; ++s) {
                    int nb[3] = {i,j,k}; nb[d] += s ? 1 : -1;
                    const bool edge = (s == 0) ? (iv[d] == dom.smallEnd(d)) : (iv[d] == dom.bigEnd(d));
                    if (edge && ebc[face_index(d,s)] != BC::Periodic) { continue; }
                    // periodic neighbours are read from the ghost cell (FillBoundary with periodicity): the wrapped index lies
                    // outside this box whenever the box does not span the domain
                    if (cl(nb[0],nb[1],nb[2]) == CellGas && o(nb[0],nb[1],nb[2]) >= 0) {
                        cur = std::min(cur, o(nb[0],nb[1],nb[2]));
                    }
                }
                if (cur < o(i,j,k)) { o(i,j,k) = cur; ch = 1; }
            });
        }
        cm.label[l]->FillBoundary(geom[l].periodicity());
        return ch;
    };

    int rounds = 0;
    for (;;) {
        int changed = 0;
        for (int l = 0; l < nlev; ++l) { changed |= same_level_sweep(l); }

        // C/F: couple uncovered coarse gas with uncovered fine gas across the interface.
        for (int L = 1; L < nlev; ++L) {
            BoxArray cba = p.levels[L].ba; cba.coarsen(ratio[L]); cba.grow(1);
            iMultiFab clab(cba, p.levels[L].dm, 1, 0); clab.setVal(-1);
            clab.ParallelCopy(*cm.label[L-1], 0, 0, 1);
            iMultiFab cunc(cba, p.levels[L].dm, 1, 0); cunc.setVal(0);
            cunc.ParallelCopy(*W_unc[L-1], 0, 0, 1);
            iMultiFab ccls(cba, p.levels[L].dm, 1, 0); ccls.setVal(CellSolid);
            ccls.ParallelCopy(cls[L-1], 0, 0, 1);

            // Fine takes min from coarse across C/F; collect fine->coarse mins in a Fab then ParallelCopy back.
            BoxArray cba0 = p.levels[L].ba; cba0.coarsen(ratio[L]);
            iMultiFab f2c(cba0, p.levels[L].dm, 1, 0);
            f2c.setVal(std::numeric_limits<int>::max());
            BoxArray const& fba = p.levels[L].ba;

            for (MFIter mfi(*cm.label[L]); mfi.isValid(); ++mfi) {
                auto const& o = cm.label[L]->array(mfi);
                auto const& cl = cls[L].const_array(mfi);
                auto const& u = W_unc[L]->const_array(mfi);
                auto const& la = clab.const_array(mfi);
                auto const& ua = cunc.const_array(mfi);
                auto const& ca = ccls.const_array(mfi);
                Box const cbox = clab[mfi].box();
                // f2c uses un-grown coarsened boxes; same index as fine.
                auto const& f2 = f2c.array(mfi);
                Box const f2box = f2c[mfi].box();
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    if (u(i,j,k) == 0 || cl(i,j,k) != CellGas || o(i,j,k) < 0) { return; }
                    for (int d = 0; d < 3; ++d) for (int s = 0; s < 2; ++s) {
                        IntVect n(i,j,k); n[d] += s ? 1 : -1;
                        if (fba.contains(n)) { continue; }
                        IntVect w = n;
                        if (!wrap_iv(geom[L], w) || fba.contains(w)) { continue; }
                        IntVect cc;
                        for (int tt = 0; tt < 3; ++tt) { cc[tt] = amrex::coarsen(w[tt], ratio[L][tt]); }
                        if (!cbox.contains(cc)) { continue; }
                        const int lc = la(cc[0],cc[1],cc[2]);
                        const int uc = ua(cc[0],cc[1],cc[2]);
                        const int kc = ca(cc[0],cc[1],cc[2]);
                        if (uc != 0 && kc == CellGas && lc >= 0) {
                            if (lc < o(i,j,k)) { o(i,j,k) = lc; changed = 1; }
                            if (f2box.contains(cc)) {
                                f2(cc[0],cc[1],cc[2]) = std::min(f2(cc[0],cc[1],cc[2]), o(i,j,k));
                            }
                        }
                    }
                });
            }
            cm.label[L]->FillBoundary(geom[L].periodicity());
            iMultiFab cmin(p.levels[L-1].ba, p.levels[L-1].dm, 1, 0);
            cmin.setVal(std::numeric_limits<int>::max());
            cmin.ParallelCopy(f2c, 0, 0, 1);
            for (MFIter mfi(*cm.label[L-1]); mfi.isValid(); ++mfi) {
                auto const& o = cm.label[L-1]->array(mfi);
                auto const& cl = cls[L-1].const_array(mfi);
                auto const& u = W_unc[L-1]->const_array(mfi);
                auto const& mn = cmin.const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    if (u(i,j,k) == 0 || cl(i,j,k) != CellGas || o(i,j,k) < 0) { return; }
                    if (mn(i,j,k) < o(i,j,k)) { o(i,j,k) = mn(i,j,k); changed = 1; }
                });
            }
            cm.label[L-1]->FillBoundary(geom[L-1].periodicity());
        }

        ParallelDescriptor::ReduceIntMax(changed);
        ++rounds;
        if (!changed) { break; }
        if (rounds > 100000) { return "masked composite labelling did not converge"; }
    }
    cm.rounds = rounds;

    // Unique roots: each rank contributes local roots; IO sorts unique list and Bcasts.
    std::vector<int> local_roots;
    for (int l = 0; l < nlev; ++l) {
        Box const& dom = geom[l].Domain();
        for (MFIter mfi(*cm.label[l]); mfi.isValid(); ++mfi) {
            auto const& o = cm.label[l]->const_array(mfi);
            auto const& cl = cls[l].const_array(mfi);
            auto const& u = W_unc[l]->const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                if (u(i,j,k) == 0 || cl(i,j,k) != CellGas || o(i,j,k) < 0) { return; }
                if (o(i,j,k) == int(gindex(dom, IntVect(i,j,k), off[l]))) { local_roots.push_back(o(i,j,k)); }
            });
        }
    }
    int nloc = int(local_roots.size());
    std::vector<int> rcounts(ParallelDescriptor::NProcs(), 0), displ(ParallelDescriptor::NProcs(), 0);
    rcounts[ParallelDescriptor::MyProc()] = nloc;
    ParallelDescriptor::ReduceIntSum(rcounts.data(), ParallelDescriptor::NProcs());
    int ntot = 0;
    for (int r = 0; r < ParallelDescriptor::NProcs(); ++r) { displ[r] = ntot; ntot += rcounts[r]; }
    std::vector<int> gathered(std::max(ntot, 1), 0);
    // MPI_Allgatherv via AMReX: use ParallelDescriptor::GatherLayoutDataToVector isn't right for ints.
    // Fallback: each root value is unique globally; OR all local roots into a bitset is impossible.
    // Use serialisation on rank 0 via Gatherv through Long vector + Bcast of unique sorted list.
    {
        std::vector<Long> send(nloc);
        for (int i = 0; i < nloc; ++i) { send[i] = local_roots[i]; }
        std::vector<Long> recv(std::max(ntot, 1));
        // Manual Gatherv using point-to-point is heavy; instead unique via set on each rank then iterative merge:
        std::set<int> uniq(local_roots.begin(), local_roots.end());
        // Exchange set sizes and contents with a ring reduce to rank 0
        int my_n = int(uniq.size());
        std::vector<int> my_v(uniq.begin(), uniq.end());
        if (ParallelDescriptor::IOProcessor()) {
            std::set<int> all(my_v.begin(), my_v.end());
            for (int src = 0; src < ParallelDescriptor::NProcs(); ++src) {
                if (src == ParallelDescriptor::MyProc()) { continue; }
                int sn = 0;
                ParallelDescriptor::Recv(&sn, size_t(1), src, 1701);
                std::vector<int> sv(sn);
                if (sn > 0) { ParallelDescriptor::Recv(sv.data(), size_t(sn), src, 1702); }
                all.insert(sv.begin(), sv.end());
            }
            gathered.assign(all.begin(), all.end());
            ntot = int(gathered.size());
        } else {
            ParallelDescriptor::Send(&my_n, size_t(1), ParallelDescriptor::IOProcessorNumber(), 1701);
            if (my_n > 0) { ParallelDescriptor::Send(my_v.data(), size_t(my_n), ParallelDescriptor::IOProcessorNumber(), 1702); }
        }
        ParallelDescriptor::Bcast(&ntot, 1, ParallelDescriptor::IOProcessorNumber());
        gathered.resize(ntot);
        if (ntot > 0) { ParallelDescriptor::Bcast(gathered.data(), ntot, ParallelDescriptor::IOProcessorNumber()); }
    }
    std::vector<int> const& all_roots = gathered;
    const int nc = int(all_roots.size());

    // Remap labels to component ids 0..nc-1; count cells; detect open (Dirichlet or Known nbr).
    std::vector<Long> count(nc, 0);
    std::vector<int> open(nc, 0);
    std::vector<int> pin_lev(nc, nlev), pin_lin(nc, std::numeric_limits<int>::max());
    // pin: lowest global index = the root value itself; decode level and cell from root.
    for (int l = 0; l < nlev; ++l) {
        Box const& dom = geom[l].Domain();
        for (MFIter mfi(*cm.label[l]); mfi.isValid(); ++mfi) {
            auto const& o = cm.label[l]->array(mfi);
            auto const& cl = cls[l].const_array(mfi);
            auto const& u = W_unc[l]->const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                if (u(i,j,k) == 0 || cl(i,j,k) != CellGas || o(i,j,k) < 0) { o(i,j,k) = -1; return; }
                auto it = std::lower_bound(all_roots.begin(), all_roots.end(), o(i,j,k));
                const int id = int(it - all_roots.begin());
                o(i,j,k) = id;
                count[id] += 1;
                const int iv[3] = {i,j,k};
                for (int d = 0; d < 3; ++d) {
                    if (iv[d] == dom.smallEnd(d) && ebc[face_index(d,0)] == BC::Dirichlet) { open[id] = 1; }
                    if (iv[d] == dom.bigEnd(d) && ebc[face_index(d,1)] == BC::Dirichlet) { open[id] = 1; }
                    for (int s = 0; s < 2; ++s) {
                        int nb[3] = {i,j,k}; nb[d] += s ? 1 : -1;
                        bool edge = (s == 0) ? (iv[d] == dom.smallEnd(d)) : (iv[d] == dom.bigEnd(d));
                        if (edge && ebc[face_index(d,s)] != BC::Periodic) { continue; }
                        if (cl(nb[0],nb[1],nb[2]) == CellKnown) { open[id] = 1; }     // periodic: ghost cell, as above
                    }
                }
            });
        }
        cm.label[l]->FillBoundary(geom[l].periodicity());
    }
    if (nc > 0) {
        ParallelDescriptor::ReduceLongSum(count.data(), nc);
        ParallelDescriptor::ReduceIntMax(open.data(), nc);
    }

    // Pin cell = decode root global index (lowest of the component).
    for (int q = 0; q < nc; ++q) {
        Long r = all_roots[q];
        int pl = 0;
        while (pl + 1 < nlev && r >= off[pl+1]) { ++pl; }
        Long loc = r - off[pl];
        Box const& dom = geom[pl].Domain();
        const Long nx = dom.length(0), ny = dom.length(1);
        IntVect iv(int(loc % nx) + dom.smallEnd(0),
                   int((loc / nx) % ny) + dom.smallEnd(1),
                   int(loc / (nx*ny)) + dom.smallEnd(2));
        ComponentInfo ci;
        ci.id = q;
        ci.singular = (open[q] == 0);
        ci.ncells = count[q];
        ci.pin = iv;
        cm.comps.push_back(ci);
        if (ci.singular) { cm.pins.push_back({pl, iv}); }
    }
    double secs = amrex::second() - t0;
    ParallelDescriptor::ReduceRealMax(secs);
    cm.seconds = secs;
    return {};
}

} // namespace pb
