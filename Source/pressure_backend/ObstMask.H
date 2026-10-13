// ObstMask.H: driver OBST -> pressure cell_class (Solid/Gas) per level, header-only (no Fortran dependency; the driver loads the list).
// List: the global OBST XB as read by FDS (after READ_OBST's clip to each mesh, input nudges kept), set once by the driver
// (driver/ObstList.H, from the bind(C) export fds_obst_count / fds_obst_xb). Snapping: FDS rounding at the level's dx
// (OB%I1 = NINT((XB1-XS)/DX), cells I1+1..I2 in FDS = I1..I2-1 here); an XB of zero thickness after snapping marks no cell.
#pragma once
#include <AMReX_Geometry.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Array4.H>
#include <array>
#include <cmath>
#include <vector>

namespace pb {

inline std::vector<std::array<double, 6>>& obst_list () { static std::vector<std::array<double, 6>> v; return v; }

// Cell-index box of one XB at the level of `geom` (FDS NINT rounding); empty when any direction snaps to zero thickness.
inline amrex::Box obst_cells (amrex::Geometry const& geom, std::array<double, 6> const& xb)
{
    amrex::IntVect lo, hi;
    for (int d = 0; d < 3; ++d) {
        const double plo = geom.ProbLo(d), dx = geom.CellSize(d);
        lo[d] = static_cast<int>(std::lround((xb[2*d] - plo) / dx));
        hi[d] = static_cast<int>(std::lround((xb[2*d+1] - plo) / dx)) - 1;
    }
    return amrex::Box(lo, hi);
}

// Cells of `bx` inside an OBST -> 1 (CellSolid), others -> 0 (CellGas). `lev` is informational (the snap uses geom's dx).
inline void fill_obst_mask (amrex::Geometry const& geom, int /*lev*/, amrex::Box const& bx, amrex::Array4<int> const& a)
{
    amrex::LoopOnCpu(bx, [&] (int i, int j, int k) { a(i, j, k) = 0; });
    for (auto const& xb : obst_list()) {
        const amrex::Box b = obst_cells(geom, xb) & bx;
        if (b.ok()) { amrex::LoopOnCpu(b, [&] (int i, int j, int k) { a(i, j, k) = 1; }); }
    }
}

// Thin per-level builder for PressureLevel::cell_class (valid cells, no ghost).
inline std::unique_ptr<amrex::iMultiFab> build_obst_class (amrex::Geometry const& geom, int lev, amrex::BoxArray const& ba, amrex::DistributionMapping const& dm)
{
    auto m = std::make_unique<amrex::iMultiFab>(ba, dm, 1, 0);
    for (amrex::MFIter mfi(*m); mfi.isValid(); ++mfi) { fill_obst_mask(geom, lev, mfi.validbox(), m->array(mfi)); }
    return m;
}

// Forced-tag box of one XB at the level of `geom`: rounded outward (floor/ceil; covers the NINT snap of this and every finer level), grown by `grow` cells.
inline amrex::Box obst_tag_box (amrex::Geometry const& geom, std::array<double, 6> const& xb, int grow)
{
    amrex::IntVect lo, hi;
    for (int d = 0; d < 3; ++d) {
        const double plo = geom.ProbLo(d), dx = geom.CellSize(d);
        lo[d] = static_cast<int>(std::floor((xb[2*d] - plo) / dx)) - grow;
        hi[d] = static_cast<int>(std::ceil((xb[2*d+1] - plo) / dx)) - 1 + grow;
    }
    return amrex::Box(lo, hi) & geom.Domain();
}

}  // namespace pb
