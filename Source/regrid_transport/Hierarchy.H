// Hierarchy.H: IR-002 grouping of &MESH entries into AMR levels, input-time checks (FR-010, FR-013) and the static hierarchy (R0, R1).
// No AMReX dependency: everything here works on integer index boxes (IBox.H) and plain structs, so it is unit-tested without the
// AMReX library. RegridAmrCore.H turns a Hierarchy into an AMReX AmrCore.
//
// Terms: "level 0" is the coarsest grid; level l+1 is level l refined by REF_RATIO(l). "Hidden direction" is a direction in which
// every mesh has a single cell (a 2-D input in the 3-D code, FR-010); it has ratio 1 on every level and is left out of all checks.
// "Proper nesting" means that a fine grid, shrunk to coarse cells and grown by N_PROPER coarse cells, lies inside the next coarser level.
#ifndef FDSRT_HIERARCHY_H_
#define FDSRT_HIERARCHY_H_

#include <array>
#include <iosfwd>
#include <string>
#include <vector>

#include "AmrInput.H"
#include "IBox.H"

namespace fdsrt {

// One FDS &MESH (after MULT expansion). Index i of the vector is FDS mesh number i+1.
struct MeshInput {
    int ijk[3] = {1, 1, 1};
    double xb[6] = {0, 0, 0, 0, 0, 0};
    int rank = 0;
};

// Result of the IR-002 grouping: which level each mesh belongs to and its index box on that level.
struct Grouping {
    bool ok = false;
    std::array<bool, 3> hidden{false, false, false};
    std::array<bool, 3> periodic{false, false, false};
    double dom_lo[3] = {0, 0, 0};
    double dom_hi[3] = {0, 0, 0};
    double dx0[3] = {0, 0, 0};          // level-0 cell size
    IVec n0{1, 1, 1};                   // level-0 cell counts
    std::vector<int> level_of_mesh;     // -1 when the mesh was rejected
    std::vector<IBox> box_of_mesh;      // index box on the mesh's own level, origin at the level-0 domain corner
    std::vector<std::vector<int>> meshes_of_level;  // mesh indices (0-based) per level
    std::vector<IVec> cum_ratio;        // cell-size ratio level 0 / level l, per direction
    std::vector<IVec> ref_ratio;        // ratio from level l-1 to l (level 0: 1)
    std::vector<IBox> gap_boxes;        // level-0 cells covered by no mesh (non-box level 0, not supported before M2)
    int top_level() const { return static_cast<int>(meshes_of_level.size()) - 1; }
    // Level-0 meshes in input order (what assemble_level0 receives).
    std::vector<int> level0_meshes() const { return meshes_of_level.empty() ? std::vector<int>() : meshes_of_level[0]; }
};

// IR-002: groups the meshes by cell size, assigns levels using REF_RATIO, and runs the ratio checks (FR-010): a mesh pair that
// touches (also across a periodic boundary) with a ratio other than 1, 2 or 4, or a different ratio per direction, is an error that
// names both meshes (1-based) and the ratio. Hidden directions are ignored. Also rejects overlapping meshes, finer meshes whose
// edges are not on the parent level's cell faces, and finer meshes outside the level-0 domain.
Grouping group_meshes(const std::vector<MeshInput>& meshes, const AmrParams& p, const std::array<bool, 3>& periodic, Report& rep);

// A grid of the hierarchy. `mesh` is the 0-based input mesh it comes from (or part of), -1 for a grid added to fill the level-0
// domain under finer meshes or to satisfy proper nesting.
struct GridBox {
    IBox box;
    int mesh = -1;
};

// A box of the refinable region (IR-008) on one level, with its origin: a finer &MESH footprint (`mesh` >= 0, 0-based) or a declared
// &AMR_REGION (`region` >= 0, 0-based). Declared boxes are stored after snapping outward to the blocking factor (the effective region).
struct RegionBox {
    IBox box;
    int mesh = -1;
    int region = -1;
};

struct Level {
    int level = 0;
    IBox domain;
    IVec ref_from_parent{1, 1, 1};
    IVec cum_ratio{1, 1, 1};
    IVec blocking_factor{1, 1, 1};
    IVec max_grid_size{1, 1, 1};
    double dx[3] = {0, 0, 0};
    std::vector<GridBox> grids;
    std::vector<RegionBox> taggable;  // effective refinable region on this level (IR-008): cells that may be tagged for refinement to level+1
};

struct Hierarchy {
    int max_level = 0;                    // finest level that may exist (MAX_LEVEL)
    int top = 0;                          // finest level that has static grids (top <= max_level)
    std::vector<Level> levels;            // levels 0..max_level; levels above `top` have geometry but no grids yet
    double dom_lo[3] = {0, 0, 0};
    double dom_hi[3] = {0, 0, 0};
    std::array<bool, 3> periodic{false, false, false};
    std::array<bool, 3> hidden{false, false, false};
    int n_proper = 1;
    int top_level() const { return top; }
};

// Builds the static hierarchy (R1) from a successful grouping: level-0 boxes are the level-0 meshes (unchanged, in input order)
// plus the part of the domain that lies under finer meshes; level l >= 1 boxes are the meshes of that level plus the boxes needed
// for proper nesting of level l+1 (snapped outward to the blocking factor), cut to MAX_GRID_SIZE. Runs the blocking-factor checks
// (level-0 domain divisible by BLOCKING_FACTOR; level l >= 1 mesh edges aligned with it) and check_nesting(). Level-0 boxes are the
// meshes as given and need not be multiples of the blocking factor. A box domain is never padded (ruling N1a applies to non-box level 0
// only, which is rejected before M2). Declared &AMR_REGION boxes are snapped outward to the blocking factor of the level they refine
// into, with a warning; the snapped boxes are the effective refinable region and are what the dump reports.
bool build_static_hierarchy(const Grouping& g, const AmrParams& p, Hierarchy& out, Report& rep);

// FR-013 assertions on a hierarchy: grids inside the domain and disjoint, level 0 tiles the domain, level l >= 1 grids aligned with
// the blocking factor and no larger than MAX_GRID_SIZE, and proper nesting. Appends one error per violation (first few per level).
void check_nesting(const Hierarchy& h, Report& rep);

// Everything in one call: grouping, checks, hierarchy. Returns false (with errors in rep) if the input must be rejected.
bool build_hierarchy_from_meshes(const std::vector<MeshInput>& meshes, const AmrParams& p, const std::array<bool, 3>& periodic,
                                 Hierarchy& out, Report& rep);

// WARNING (not error) of FR-010 when the level-0 cell counts cannot be halved at least 3 times in every non-hidden direction;
// returns the number of halvings possible in all directions and sets coarsest[] to the resulting coarsest MLMG grid size.
int mlmg_coarsening(const Grouping& g, IVec& coarsest, Report& rep);

// Deterministic text dump of the hierarchy (FR-010 acceptance: the dump equals the input boxes after nesting/blocking).
void dump_hierarchy(std::ostream& os, const Hierarchy& h);

// IR-008: tags outside the refinable region are discarded and counted. clip_tags() removes tags (level-l cells) outside
// Level::taggable and returns how many were removed; discarded(l) is the running total per level; first_discard(l) is true once
// per level so that the log message is written once.
class TagClipper {
public:
    explicit TagClipper(const Hierarchy& h) : h_(h), count_(h.levels.size(), 0), reported_(h.levels.size(), false) {}
    long long clip_tags(int lev, std::vector<IVec>& tags);
    long long discarded(int lev) const { return count_[lev]; }
    bool first_discard(int lev);

private:
    const Hierarchy& h_;
    std::vector<long long> count_;
    std::vector<bool> reported_;
};

}  // namespace fdsrt

#endif
