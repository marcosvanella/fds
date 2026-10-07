// IlGhost.cpp: see IlGhost.H.
#include "IlGhost.H"

#include <AMReX_BoxArray.H>
#include <AMReX_Loop.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <cstring>
#include <numeric>
#include <vector>

namespace fdsrt {
namespace {

// Ghost-fastest layout: index = ll + nicr * (n + nra * ib). One pair slab = nplane contiguous segments of length (hi-lo).
void copy_slab(double* dst, const double* src, int lo, int hi, int nicr, int nplane)
{
    const int w = hi - lo;
    for (int s = 0; s < nplane; ++s) std::memcpy(dst + lo + s * nicr, src + lo + s * nicr, sizeof(double) * w);
}

void pack_slab(double* buf, const double* src, int lo, int hi, int nicr, int nplane)
{
    const int w = hi - lo;
    for (int s = 0; s < nplane; ++s) std::memcpy(buf + s * w, src + lo + s * nicr, sizeof(double) * w);
}

void unpack_slab(double* dst, const double* buf, int lo, int hi, int nicr, int nplane)
{
    const int w = hi - lo;
    for (int s = 0; s < nplane; ++s) std::memcpy(dst + lo + s * nicr, buf + s * w, sizeof(double) * w);
}

}  // namespace

IlGhostTable build_il_ghost_table(const std::vector<IlBox>& boxes, const std::vector<amrex::IntVect>& ratio, const std::vector<IlWall>& walls)
{
    int nlev = 0;
    for (const IlBox& b : boxes) nlev = std::max(nlev, b.level + 1);
    std::vector<amrex::IntVect> rr(nlev, amrex::IntVect(1));   // ratio from level 0
    std::vector<amrex::BoxList> bl(nlev);
    std::vector<std::vector<int>> id(nlev);                    // box index of each BoxArray entry
    for (int l = 1; l < nlev; ++l) rr[l] = rr[l - 1] * ratio[l - 1];
    for (int b = 0; b < static_cast<int>(boxes.size()); ++b) {
        bl[boxes[b].level].push_back(boxes[b].box);
        id[boxes[b].level].push_back(b);
    }
    std::vector<amrex::BoxArray> ba(bl.begin(), bl.end());

    // sender box and sender cells (region, sender level index space) of each wall; -1 when no box covers the ghost cell
    const int nw = static_cast<int>(walls.size());
    std::vector<int> src(nw, -1);
    std::vector<amrex::Box> reg(nw);
    for (int w = 0; w < nw; ++w) {
        const int lr = boxes[walls[w].box].level, d = std::abs(walls[w].ior) - 1;
        for (int k = 0; k < nlev && src[w] < 0; ++k) {
            const int ls = k <= lr ? lr - k : k;   // own level, coarser levels, finer levels
            amrex::Box r(walls[w].g, walls[w].g);
            if (ls < lr) r.coarsen(rr[lr] / rr[ls]);
            if (ls > lr) {
                r.refine(rr[ls] / rr[lr]);
                if (walls[w].ior > 0) r.setSmall(d, r.bigEnd(d)); else r.setBig(d, r.smallEnd(d));
            }
            const auto is = ba[ls].intersections(r);
            if (is.empty()) continue;
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(is.size() == 1 && is[0].second == r, "build_il_ghost_table: a wall face spans more than one box");
            src[w] = id[ls][is[0].first];
            reg[w] = r;
        }
    }

    std::vector<int> ord(nw);
    std::iota(ord.begin(), ord.end(), 0);
    std::stable_sort(ord.begin(), ord.end(), [&](int a, int b) { return std::make_pair(walls[a].box, src[a]) < std::make_pair(walls[b].box, src[b]); });
    IlGhostTable t;
    t.nic_min.assign(nw, 1);
    t.nic_max.assign(nw, 0);
    for (int i = 0; i < nw; ++i) {
        const int w = ord[i], s = src[w];
        if (s < 0) continue;
        if (t.pair_lo.empty() || walls[w].box != walls[ord[i - 1]].box || s != src[ord[i - 1]]) {
            t.pair_lo.push_back(t.nicr());
            t.pair_recv.push_back(walls[w].box);
        }
        const amrex::IntVect o = boxes[s].box.smallEnd() - 1;
        t.nic_min[w] = t.nicr() + 1;
        amrex::LoopOnCpu(reg[w], [&](int a, int b, int c) {
            t.ils_box.push_back(s + 1);
            t.iio.push_back(a - o[0]);
            t.jjo.push_back(b - o[1]);
            t.kko.push_back(c - o[2]);
        });
        t.nic_max[w] = t.nicr();
    }
    t.pair_lo.push_back(t.nicr());
    return t;
}

void exchange_il_ghost(const IlGhostTable& t, const int* box_rank, int nra, int nsb, double* ils, double* ilr)
{
    namespace PD = amrex::ParallelDescriptor;
    const int me = PD::MyProc(), nicr = t.nicr(), np = static_cast<int>(t.pair_lo.size()) - 1, nplane = nra * nsb;
    constexpr int tag0 = 1251;

    std::vector<std::vector<double>> sbuf, rbuf;
    std::vector<int> r_lo, r_hi;
    std::vector<PD::Message> msgs;

    // Same-rank slabs and post all recvs first.
    for (int p = 0; p < np; ++p) {
        const int lo = t.pair_lo[p], hi = t.pair_lo[p + 1];
        const int sr = box_rank[t.ils_box[lo] - 1], rr = box_rank[t.pair_recv[p]];
        if (sr == rr) {
            if (me == sr) copy_slab(ilr, ils, lo, hi, nicr, nplane);
            continue;
        }
        if (me != rr) continue;
        r_lo.push_back(lo);
        r_hi.push_back(hi);
        rbuf.emplace_back(static_cast<size_t>(hi - lo) * nplane);
        msgs.push_back(PD::Arecv(rbuf.back().data(), rbuf.back().size(), sr, tag0 + p));
    }
    // Pack and send.
    for (int p = 0; p < np; ++p) {
        const int lo = t.pair_lo[p], hi = t.pair_lo[p + 1];
        const int sr = box_rank[t.ils_box[lo] - 1], rr = box_rank[t.pair_recv[p]];
        if (sr == rr || me != sr) continue;
        sbuf.emplace_back(static_cast<size_t>(hi - lo) * nplane);
        pack_slab(sbuf.back().data(), ils, lo, hi, nicr, nplane);
        msgs.push_back(PD::Asend(sbuf.back().data(), sbuf.back().size(), rr, tag0 + p));
    }
    for (size_t i = 0; i < msgs.size(); ++i) msgs[i].wait();
    for (size_t i = 0; i < rbuf.size(); ++i) unpack_slab(ilr, rbuf[i].data(), r_lo[i], r_hi[i], nicr, nplane);
}

}  // namespace fdsrt
