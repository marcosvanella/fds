// IlGhost.H: L1251, the radiation intensity ghost table of AMR boxes (notes/radiation-il-ghost.md). Host only, built once per regrid.
//
// The AMR analogue of FDS's IL_S/IL_R (main.f90:2191-2254, radi.f90:4357-4361, 4831). The sweep (s5_rad_sweep, s5rad_sweep_k2.F90) uses it unchanged:
//   s5rad_ils   packs ILS(LL,N,IBND) = IL(IIO_S(LL),JJO_S(LL),KKO_S(LL),ILS_BOX(LL)) during the sweep of angle N, band IBND;
//   after each radiation update the caller copies ILS into ILR (whole array, the FDS lag: same-process IL_R = IL_S);
//   s5rad_wall1 reads at an INTERPOLATED_BOUNDARY wall IW the mean of ILR(W_NIC_MIN(IW):W_NIC_MAX(IW),N,IBND).
// So ILS and ILR share one layout ILR(NICR,NRA,NSB), position index fastest; this table fills no value and averages nothing.
//
// Positions: grouped by (receiver box, sender box) in ascending order, so each neighbour pair is one contiguous range (pair_lo); inside a pair the
// walls come in the caller's order and the sender cells of one wall are consecutive, K then J then I with I fastest (FDS IIO_R order).
// Sender cells of a wall with ghost cell g (receiver level lr, sender level ls, per-direction ratios):
//   ls == lr: g itself;   ls < lr: the coarse cell containing g (one entry);   ls > lr: the fine cells under the face, i.e. g refined and cut to the
//   layer next to the receiver (NIC = product of the tangential ratios; the sweep's mean is the area mean).
// Search order: the receiver's level, then coarser levels (nearest first), then finer levels. For a disjoint composite tiling (FDS-like meshes)
// exactly one box covers g; for full levels (coarse under fine) a box takes its ghosts from its own level first, then from the coarser one.
#ifndef FDSRT_ILGHOST_H_
#define FDSRT_ILGHOST_H_

#include <AMReX_Box.H>

#include <vector>

namespace fdsrt {

struct IlBox {
    int level;
    amrex::Box box;      // cells in the index space of `level`
};

struct IlWall {
    int box;             // receiver box (0-based index into the box list)
    int ior;             // FDS orientation: +1 is the low x face of the box (ghost i = lo-1), -1 the high x face, +-2 y, +-3 z
    amrex::IntVect g;    // ghost cell, index space of the receiver's level
};

struct IlGhostTable {
    std::vector<int> nic_min, nic_max;        // per wall, 1-based positions; nic_max = nic_min - 1 when no box covers the ghost cell (not an interpolated wall)
    std::vector<int> ils_box, iio, jjo, kko;  // per position: sender box (1-based) and sender cell in the sender's own 1-based index space
    std::vector<int> pair_lo;                 // pair p covers positions pair_lo[p] .. pair_lo[p+1]-1 (0-based); one more entry than pairs
    std::vector<int> pair_recv;               // 0-based receiver box of pair p (sender = ils_box[pair_lo[p]] - 1)
    int nicr() const { return static_cast<int>(ils_box.size()); }
};

// L1251 entry point. ratio[l] = refinement ratio from level l to level l+1 (per direction, 1 in a single-cell direction).
IlGhostTable build_il_ghost_table(const std::vector<IlBox>& boxes, const std::vector<amrex::IntVect>& ratio, const std::vector<IlWall>& walls);

// Multi-rank ILS→ILR exchange. Layout (ghost-fastest, matches test_il_ghost / FDS): index = ll + nicr * (n + nra * ibnd).
// Each pair_lo slab is one message (full NRA×NSB); direction sender-box rank → receiver-box rank. Same-rank pairs: memcpy the slabs.
// box_rank[b] = owning rank of box b (size = number of boxes). Pack of swept angles into ILS is the caller's job; this only copies ILS into ILR.
void exchange_il_ghost(const IlGhostTable& t, const int* box_rank, int nra, int nsb, double* ils, double* ilr);

}  // namespace fdsrt

#endif
