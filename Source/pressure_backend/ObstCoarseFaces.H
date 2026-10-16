// ObstCoarseFaces.H: fine-box faces whose neighbour cell is a SOLID cell of the covering coarser level (header-only, uses ObstMask.H).
// A fine box face with a coarse solid neighbour (the neighbour is not in the fine level's BoxArray, and the coarse cell holding it is an OBST cell of level lev-1) is a wall for the
// transport: its normal velocity is zero, ADV = DIF = 0, and it carries no coarse-fine override entry. Same rule for every level pair (l-1, l); level 0 is FDS's own OBST mask.
#pragma once
#include <AMReX_BaseFab.H>
#include <AMReX_BoxArray.H>
#include "ObstMask.H"

namespace pb {

// OBST cell (SOLID) of level `lev` at cell iv.
inline bool obst_solid_cell (amrex::Geometry const& g, int lev, amrex::IntVect const& iv)
{
    if (obst_list().empty()) return false;
    amrex::BaseFab<int> f(amrex::Box(iv, iv)); f.setVal<amrex::RunOn::Host>(0);
    fill_obst_mask(g, lev, f.box(), f.array());
    return f(iv) != 0;
}

// (also: a neighbour in another box of the same level that is an OBST cell of this level). Faces (dir, nodal fine index) of box vb of the fine level (ba = its BoxArray, ratio to the coarse level lev-1, geometry of lev-1) whose outside neighbour is a coarse SOLID cell.
inline std::vector<ThinFace> coarse_solid_faces (amrex::Geometry const& gcrs, int lcrs, amrex::IntVect const& ratio, amrex::BoxArray const& fba, amrex::Box const& vb, amrex::Geometry const& gfine)
{
    std::vector<ThinFace> out;
    if (obst_list().empty()) return out;
    const amrex::Box fdom = amrex::refine(gcrs.Domain(), ratio);
    const amrex::Box cb = amrex::coarsen(amrex::grow(vb, 1), ratio);
    amrex::BaseFab<int> m(cb); m.setVal<amrex::RunOn::Host>(0);
    fill_obst_mask(gcrs, lcrs, cb, m.array());
    amrex::BaseFab<int> mf(amrex::grow(vb, 1)); mf.setVal<amrex::RunOn::Host>(0);   // the level's own mask: a face between two boxes of this level with a solid cell across is a wall too
    fill_obst_mask(gfine, lcrs + 1, mf.box(), mf.array());
    for (int d = 0; d < 3; ++d)
        for (int side = 0; side < 2; ++side) {
            amrex::Box slab = vb; if (side == 0) slab.setBig(d, vb.smallEnd(d)); else slab.setSmall(d, vb.bigEnd(d));
            for (amrex::BoxIterator bi(slab); bi.ok(); ++bi) {
                amrex::IntVect n = bi(); n[d] += side == 0 ? -1 : 1;
                if (!fdom.contains(n)) continue;
                if (fba.contains(n) ? mf(n) == 0 : m(amrex::coarsen(n, ratio)) == 0) continue;
                amrex::IntVect f = bi(); if (side == 1) f[d] += 1;
                out.push_back(ThinFace{d, f[0], f[1], f[2]});
            }
        }
    return out;
}

}  // namespace pb
