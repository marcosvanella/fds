// SideData.H: per-box side data rebuilt from the FDS set-up (M2a, S2): solid flag, face mask and source-cell flag on valid+2.
//
// Kernel-facing rules (M2a): (a) passive scalars are not involved: the mask is per cell, not per scalar (handled in Fields.cpp);
// (b) only uniform Cartesian metrics are used; this file has no metric data (cell metrics are constant per direction, Level0::dx).
//
// The mask follows the IR-007 face-mask rule (D-031): a face is nonzero wherever single-mesh FDS would give it a nonzero WALL_INDEX,
// i.e. it is zero only on faces between two boxes of the same level. It is NOT derived from WALL_INDEX of the box's own FDS mesh alone:
// FDS gives a mesh-to-mesh interface a wall cell (INTERPOLATED_BOUNDARY), which is zeroed here.
//
// iMultiFab with kNComp components and ng = 2 (valid+2), same BoxArray and DistributionMapping as level 0:
//   comp 0    SOLID cell flag (FDS CELL%SOLID)
//   comp 1..6 face mask of the cell's faces (-x, +x, -y, +y, -z, +z) = WALL_INDEX(-1), (1), (-2), (2), (-3), (3) as 0/1
//   comp 7    source flag: 1 if the cell lies inside the level-0 domain (valid cell or ghost of a neighbouring box). A periodic image
//             outside the domain and a physical-boundary ghost are 0: they never set clip flags.
#ifndef FDS_SIDEDATA_H_
#define FDS_SIDEDATA_H_

#include <AMReX_Array4.H>
#include <AMReX_iMultiFab.H>

#include <cstdint>
#include <functional>

#include "FdsAmr.H"

namespace fdsamr {

constexpr int kMaskNComp = 8;
constexpr int kMaskNG = 2;

// Wall code delivered by a provider for each face of a valid cell.
enum : int { kNoWall = 0, kWall = 1, kInterface = 2 };

// Fills comps 0..6 of the valid cells of the box `vbox` of FDS mesh `nm` (= box index + 1): comp 0 SOLID (0/1), comps 1..6 the wall code
// (kNoWall, kWall, kInterface) of the face. Called for local boxes only.
using CellWallProvider = std::function<void(int nm, const amrex::Box& vbox, amrex::Array4<int> const& a)>;

class SideData {
public:
    SideData(const Level0& l0, const CellWallProvider& provider);
    const amrex::iMultiFab& mask() const { return m_mask; }
    amrex::iMultiFab& mask() { return m_mask; }

    // Order-independent (exact, integer) hash over the valid cells of all components: equal for equal global fields, whatever the box
    // layout or rank count.
    std::uint64_t hash() const;
    // Counts of nonzero face flags over valid cells (sum over ranks): [0] solid cells, [1] wall faces, [2] cells flagged as source.
    void counts(long out[3]) const;

private:
    amrex::iMultiFab m_mask;
};

}  // namespace fdsamr

#endif
