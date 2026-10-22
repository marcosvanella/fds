// FaceSync.H: one value per shared face of a face (nodal) MultiFab.
#pragma once
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <memory>

namespace pb {

// The owned faces of a face MultiFab: a copy on its BoxArray with each box trimmed to the faces it owns (see below), so every face is in exactly one box.
// An empty pointer for a cell-centred MultiFab. ParallelCopy from it is a single-source copy.
// Every high face below the domain top is trimmed whether or not a neighbour box has it: as a copy SOURCE this needs mf to cover the domain (proper nesting;
// level 0). sync_faces_low_owner below is safe on any level: a face no box owns keeps its value.
inline std::unique_ptr<amrex::MultiFab> owned_faces(const amrex::MultiFab& mf, const amrex::Geometry& geom, bool high_owner = false)
{
    const amrex::IntVect nod = mf.ixType().toIntVect();
    int d = -1; for (int e = 0; e < 3; ++e) if (nod[e]) d = e;
    if (d < 0 || nod.sum() != 1) return nullptr;
    const int top = geom.Domain().bigEnd(d) + 1, bot = geom.Domain().smallEnd(d);   // nodal indices of the domain's high and low faces
    amrex::BoxList bl;
    for (int b = 0; b < static_cast<int>(mf.size()); ++b) {
        amrex::Box x = mf.boxArray()[b];
        if (high_owner) { if (geom.isPeriodic(d) || x.smallEnd(d) > bot) x.growLo(d, -1); }
        else if (geom.isPeriodic(d) || x.bigEnd(d) < top) x.growHi(d, -1);
        bl.push_back(x);
    }
    const int nc = mf.nComp();
    auto own = std::make_unique<amrex::MultiFab>(amrex::BoxArray(std::move(bl)), mf.DistributionMap(), nc, 0);
    for (amrex::MFIter mfi(*own); mfi.isValid(); ++mfi) {
        auto o = own->array(mfi); const auto s = mf.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), nc, [&](int i, int j, int k, int n) { o(i, j, k, n) = s(i, j, k, n); });
    }
    return own;
}

// A face shared by two boxes (a seam, or the two copies of a periodic domain face) is valid in both, and each box may hold its own copy (equal only to round-off).
// Set every such face to the copy of the box that has it as its LOW face (the low-face owner); a non-periodic domain's high end keeps its only copy. Each face then
// has one source, so a later average-down or ParallelCopy does not depend on which copy is copied last (local before received: the rank layout).
// high_owner: the box that has the face as its HIGH face owns it instead (FVX/FVY/FVZ shifted to nodal faces: a box's low face there is its ghost FV(-1), the high face
// its computed FV(n)); a non-periodic domain's low end then keeps its only copy.
inline void sync_faces_low_owner(amrex::MultiFab& mf, const amrex::Geometry& geom, bool high_owner = false)
{
    const auto own = owned_faces(mf, geom, high_owner);
    if (own) mf.ParallelCopy(*own, 0, 0, mf.nComp(), 0, 0, geom.periodicity());
}

}  // namespace pb
