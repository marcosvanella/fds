// OutputPlan.H: ADR-004 / FR-072 Step C skeleton (notes/output-plan.md). Static disjoint output meshes at L_out inside the
// refinable region and level-0 &MESH outside; writer ranks = contiguous blocks balanced by cell count (DistributionMapping(Vector<int>)).
// Fill: ParallelCopy COPY from level 0 + average_down from finer (host); OverrideSync on face data before a face fill.
// Fortran OUT_MESH binding (FDS_FINE_B_SET_VIEW / BUILD_FINE_BOX device) waits on N_OUT_MESHES > NMESHES validation (Mapper: tables-only until then).
#ifndef FDS_OUTPUTPLAN_H_
#define FDS_OUTPUTPLAN_H_

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

#include <memory>
#include <string>
#include <vector>

#include "ExactSum.H"
#include "LevelRegistry.H"

namespace fdsamr {

struct OutputPlanOptions {
    amrex::Box refinable;          // level-0 index space; empty = no refinable region (uniform: output = level-0 boxes)
    int L_out = 1;                 // output level inside refinable (ratio^L_out from level 0)
    amrex::IntVect ratio{2, 2, 2}; // refinement ratio per level (uniform)
    int blocking = 8;
    int maxsize = 32;
    std::vector<std::string> cell_names{"RHO", "TMP", "H", "D", "MU", "KRES", "RSUM"};
    std::vector<std::string> face_names{"U", "V", "W"};   // OverrideSync before fill
};

// One resolution group of the output mesh set (all boxes at the same cell size).
struct OutputGroup {
    int level = 0;                 // 0 = level-0 resolution; L_out = fine
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm; // writer layout
    std::unique_ptr<amrex::MultiFab> cell;   // ncomp = cell_names.size(), no ghost
    std::unique_ptr<amrex::MultiFab> face[3]; // one each, or null if face_names empty
    std::unique_ptr<amrex::MultiFab> amr_level; // D2 AMR_LEVEL quantity (which compute level filled the cell)
};

class OutputPlan {
public:
    // Build the static layout from the registry's level-0 domain and the options. Does not fill.
    void make(const LevelRegistry& reg, const OutputPlanOptions& o);
    // Fill cell (and face) MultiFabs from the composite hierarchy at this output time.
    void fill(const LevelRegistry& reg);
    // True when Fortran OUT_MESH objects can be bound (N_OUT_MESHES validated). Always false until that gate opens.
    bool bind_ready() const { return false; }
    const std::vector<OutputGroup>& groups() const { return m_g; }
    const OutputPlanOptions& options() const { return m_o; }
    // Writer rank of output-mesh index `m` in group `g` (contiguous blocks).
    int writer_of(int g, int m) const;

private:
    OutputPlanOptions m_o;
    std::vector<OutputGroup> m_g;
    static amrex::DistributionMapping writers_by_cells(const amrex::BoxArray& ba);
};

// D6: exact composite mass (RHO and RHO*ZZ_n) over uncovered cells of levels 0..nl-1. Uses ExactSum hierarchy / uncovered product.
std::vector<double> composite_mass_exact(const LevelRegistry& reg, int nl, int ns);

}  // namespace fdsamr

#endif
