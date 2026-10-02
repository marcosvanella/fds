// SideData.cpp: see SideData.H. Kernel-facing rules (M2a): (a) passive scalars are not involved (Fields.cpp); (b) only uniform Cartesian
// metrics are used.
#include "SideData.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>

namespace fdsamr {

namespace {
inline std::uint64_t mix(std::uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
}  // namespace

SideData::SideData(const Level0& l0, const CellWallProvider& provider)
    : m_mask(l0.ba, l0.dm, kMaskNComp, amrex::IntVect(kMaskNG))
{
    m_mask.setVal(0);
    const amrex::Box dom = l0.geom.Domain();
    for (amrex::MFIter mfi(m_mask); mfi.isValid(); ++mfi) {
        const int nm = mfi.index() + 1 + l0.fds_mesh_offset;   // FDS mesh number of the box (level 0: box index + 1)
        const amrex::Box vb = mfi.validbox();
        auto a = m_mask.array(mfi);
        provider(nm, vb, a);
        // Interface faces (mesh to mesh, same level) are open in the mask: FDS gives them a wall cell, single-mesh FDS does not.
        // Domain-boundary walls (including periodic faces) stay nonzero. A periodic domain edge carries the interface code when the neighbour is another mesh (or this mesh): it stays a wall.
        amrex::LoopOnCpu(vb, [&](int i, int j, int k) {
            for (int f = 0; f < 6; ++f) {
                int& c = a(i, j, k, 1 + f);
                if (c == kInterface) {
                    amrex::IntVect nb(i, j, k);
                    nb[f / 2] += (f % 2 == 0) ? -1 : 1;
                    c = dom.contains(nb) ? kNoWall : kWall;   // across a periodic domain edge FDS also codes the neighbour mesh as an interface: still a wall cell there
                } else if (c != kNoWall && c != kWall) {
                    amrex::Abort("SideData: invalid wall code from the provider");
                }
            }
        });
    }
    // valid+2 by exchange (full, edges and corners included: the clip reads them), periodic images included
    m_mask.FillBoundary(0, kMaskNComp - 1, m_mask.nGrowVect(), l0.geom.periodicity());
    // source flag from geometry: inside the domain only (a periodic image is not a source)
    for (amrex::MFIter mfi(m_mask); mfi.isValid(); ++mfi) {
        auto a = m_mask.array(mfi);
        amrex::LoopOnCpu(mfi.fabbox(), [&](int i, int j, int k) { a(i, j, k, 7) = dom.contains(amrex::IntVect(i, j, k)) ? 1 : 0; });
    }
}

std::uint64_t SideData::hash() const
{
    std::uint64_t h = 0;
    for (amrex::MFIter mfi(m_mask); mfi.isValid(); ++mfi) {
        auto a = m_mask.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            for (int c = 0; c < kMaskNComp; ++c) {
                const std::uint64_t key = mix((static_cast<std::uint64_t>(static_cast<std::uint32_t>(i)) << 40) ^
                                              (static_cast<std::uint64_t>(static_cast<std::uint32_t>(j)) << 20) ^
                                              static_cast<std::uint64_t>(static_cast<std::uint32_t>(k)) ^
                                              (static_cast<std::uint64_t>(c) << 60));
                h += key * static_cast<std::uint64_t>(a(i, j, k, c) + 1);   // wraps modulo 2^64 by definition
            }
        });
    }
    // exact integer sum over ranks (gather, then add in a fixed order: unsigned addition is associative anyway)
    const int np = amrex::ParallelDescriptor::NProcs();
    std::vector<std::uint64_t> all(np, 0);
    if (np > 1) {
        amrex::ParallelDescriptor::Gather(&h, 1, all.data(), 1, amrex::ParallelDescriptor::IOProcessorNumber());
        amrex::ParallelDescriptor::Bcast(all.data(), np, amrex::ParallelDescriptor::IOProcessorNumber());
        std::uint64_t s = 0;
        for (auto v : all) s += v;
        return s;
    }
    return h;
}

void SideData::counts(long out[3]) const
{
    long c[3] = {0, 0, 0};
    for (amrex::MFIter mfi(m_mask); mfi.isValid(); ++mfi) {
        auto a = m_mask.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
            c[0] += a(i, j, k, 0);
            for (int f = 1; f <= 6; ++f) c[1] += a(i, j, k, f);
            c[2] += a(i, j, k, 7);
        });
    }
    amrex::ParallelDescriptor::ReduceLongSum(c, 3);
    for (int i = 0; i < 3; ++i) out[i] = c[i];
}

}  // namespace fdsamr
