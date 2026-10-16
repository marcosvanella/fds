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
// D-010: XB of the stateful walls (burner VENTs etc.) that stay on level 0; RegridAmrCore clears the tags around them
inline std::vector<std::array<double, 6>>& nofine_list () { static std::vector<std::array<double, 6>> v; return v; }

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


// ---- thin OBST (owner rule): an OBST that vanishes at level l (zero thickness after the NINT snap in exactly one direction, thick in the other two) has no solid cells on
// level l or any finer level; there it is a zero-thickness wall FACE set. The coarsest vanishing level l* wins: the plane is the level-l* snap plane (x refined by the ratio) and the
// transverse extent is the level-l* snap refined (nested: every coarse face is r*r fine faces). Rods/points (two or more zero directions) keep the per-level snap (dropped when empty).
// The level geometries are set by the driver (level 0 first); the rule needs the cell size of every level up to the asked one.
inline std::vector<amrex::Geometry>& level_geoms () { static std::vector<amrex::Geometry> v; return v; }
inline void set_level_geoms (std::vector<amrex::Geometry> const& g) { level_geoms() = g; }
inline void note_level_geom (int lev, amrex::Geometry const& g) { auto& v = level_geoms(); if (static_cast<int>(v.size()) <= lev) v.resize(lev + 1, g); v[lev] = g; }   // a level's geometry, noted when its cell walls are built

struct ThinFace { int dir, i, j, k; };   // AMReX nodal index in the level's index space: the face between cell (i,j,k)-e_dir and cell (i,j,k) (low face of cell (i,j,k) in direction dir)

// Coarsest level l* <= lev at which `xb` is thin (-1: none). Fills the level-l* snap (lo/hi per direction) and the thin direction.
inline int thin_level (std::array<double, 6> const& xb, int lev, int& dir, int plo_[3], int phi_[3])
{
    auto const& G = level_geoms();
    if (G.empty()) return -1;
    for (int k = 0; k <= lev && k < static_cast<int>(G.size()); ++k) {
        int nz = 0, dz = 0; bool ok = true; int lo[3], hi[3];
        for (int d = 0; d < 3; ++d) {
            lo[d] = static_cast<int>(std::lround((xb[2*d] - G[0].ProbLo(d)) / G[k].CellSize(d)));
            hi[d] = static_cast<int>(std::lround((xb[2*d+1] - G[0].ProbLo(d)) / G[k].CellSize(d)));
            if (hi[d] == lo[d]) { ++nz; dz = d; } else if (hi[d] < lo[d]) ok = false;
        }
        if (ok && nz == 1) { dir = dz; for (int d = 0; d < 3; ++d) { plo_[d] = lo[d]; phi_[d] = hi[d]; } return k; }
    }
    return -1;
}

// Thin faces of level lev (>= 1) of every OBST, the same set on every box layout (geometry only); faces on a domain boundary plane are dropped (the domain wall exists anyway).
inline std::vector<ThinFace> thin_faces (int lev)
{
    std::vector<ThinFace> out;
    auto const& G = level_geoms();
    if (lev < 1 || lev >= static_cast<int>(G.size())) return out;
    for (auto const& xb : obst_list()) {
        int dir = 0, lo[3], hi[3];
        const int ls = thin_level(xb, lev, dir, lo, hi);
        if (ls < 0) continue;
        const int r = static_cast<int>(std::lround(G[ls].CellSize(dir) / G[lev].CellSize(dir)));   // refinement of the vanishing level (uniform ratio per direction)
        const int plane = lo[dir] * r;
        const amrex::Box dom = G[lev].Domain();
        if (plane <= dom.smallEnd(dir) || plane > dom.bigEnd(dir)) continue;
        const int t1 = (dir + 1) % 3, t2 = (dir + 2) % 3;
        const int rt1 = static_cast<int>(std::lround(G[ls].CellSize(t1) / G[lev].CellSize(t1))), rt2 = static_cast<int>(std::lround(G[ls].CellSize(t2) / G[lev].CellSize(t2)));
        for (int b = lo[t2] * rt2; b < hi[t2] * rt2; ++b)
            for (int a = lo[t1] * rt1; a < hi[t1] * rt1; ++a) {
                ThinFace f; f.dir = dir; int c[3]; c[dir] = plane; c[t1] = a; c[t2] = b; f.i = c[0]; f.j = c[1]; f.k = c[2]; out.push_back(f);
            }
    }
    return out;
}

// Cells of `bx` inside an OBST -> 1 (CellSolid), others -> 0 (CellGas). `lev` is informational (the snap uses geom's dx).
inline void fill_obst_mask (amrex::Geometry const& geom, int lev, amrex::Box const& bx, amrex::Array4<int> const& a)
{
    amrex::LoopOnCpu(bx, [&] (int i, int j, int k) { a(i, j, k) = 0; });
    for (auto const& xb : obst_list()) {
        if (lev >= 1) { int td = 0, tl[3], th[3]; if (thin_level(xb, lev, td, tl, th) >= 0) continue; }   // vanished at this or a coarser level: a thin face set (thin_faces), no solid cells
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
