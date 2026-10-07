// TagOps.H: C++ side of the tagging kernels (R4). The kernels are Fortran (rt_tag_kernels.F90, K2 default of ADR-001); this file only
// loops over the boxes of a level and passes flat pointers plus index bounds. No K1 kernel was needed (reasons in notes/r4-tagging-kernels.md).
//
// A TagBoxArray is AMReX's per-level array of one-byte tags (the cells flagged "refine here"). The functions below only SET tags.
//   tag_cells        : one criterion on one cell field (threshold, or undivided difference to the neighbours), with TAG_KEEP hysteresis
//   tag_boxes        : force-tag index boxes (user boxes, static finer &MESH footprints)
//   clip_tags_to_boxes: remove tags outside the allowed boxes (the refinable region of IR-008) and return how many were removed
//   count_tags       : number of tags set (valid cells, summed over ranks)
//   make_covered_mask: 1 on the cells of a level that lie under the next finer level (input of the hysteresis)
// The tag buffer (N_ERROR_BUF) is done here (buffer_tags: ghost exchange then Chebyshev dilation) so it is complete across box
// boundaries (D-080); AmrCore is given n_error_buf = 0 so its own buffer does not double-grow. Proper nesting stays with AmrMesh.
#ifndef FDSRT_TAGOPS_H_
#define FDSRT_TAGOPS_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_TagBox.H>

#include "AmrInput.H"

#include <array>
#include <vector>

namespace fdsrt {

enum class TagMode { Above = 0, Diff = 1 };

struct TagCriterion {
    TagMode mode = TagMode::Diff;
    double thr = 0.0;        // threshold (Above: on v - base; Diff: on the largest undivided difference, relative if `relative`)
    double base = 0.0;       // Above only: reference value (e.g. ambient temperature), the test is v - base > thr
    double keepfac = 1.0;    // TAG_KEEP: cells under the finer level use thr*keepfac (0 < keepfac <= 1); 1 = no hysteresis
    bool relative = false;   // Diff only: divide each difference by max(|v|, |v_nb|)
    int comp = 0;            // component of the field (the numerator when `den` is given)
    std::array<int, 3> dirs{1, 1, 1};   // 0 for a single-cell direction (not read)
};

// Sets the tag in every valid cell of `tags` where the criterion holds on field `q` (same BoxArray and DistributionMapping as `tags`).
// With `den` (same layout as q), the value is q(comp)/den(den_comp), e.g. a mass fraction from rho*Z and rho.
// `covered` (iMultiFab, one component) enables the hysteresis; nullptr = none. q (and den) need one valid ghost layer in the active directions.
void tag_cells(amrex::TagBoxArray& tags, const amrex::MultiFab& q, const amrex::MultiFab* den, int den_comp, const amrex::iMultiFab* covered,
               const TagCriterion& c);

// Sets all tags inside each box (index space of the level of `tags`).
void tag_boxes(amrex::TagBoxArray& tags, const std::vector<amrex::Box>& boxes);

// Keeps only the tags inside the union of `allowed`; returns the number of tags removed (sum over ranks).
long clip_tags_to_boxes(amrex::TagBoxArray& tags, const std::vector<amrex::Box>& allowed);

long count_tags(const amrex::TagBoxArray& tags);

// Grow every TagBox::SET cell by `nbuf` cells (Chebyshev, all directions including edges and corners) into TagBox::BUF.
// Uses a temporary TagBoxArray with `nbuf` ghost cells (FillBoundary with `geom` periodicity) and marks the CLEAR valid cells of
// `tags` that have a SET within `nbuf`; nbuf >= 0 (0: no-op). Matches the single-box behaviour of AMReX TagBoxArray::buffer and
// also covers cells across box and rank boundaries (D-080 / D-053).
void buffer_tags(amrex::TagBoxArray& tags, int nbuf, const amrex::Geometry& geom);

// Mask on (ba, dm): 1 where the cell is covered by `fine_ba` coarsened by `ratio`, else 0.
amrex::iMultiFab make_covered_mask(const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, const amrex::BoxArray& fine_ba, const amrex::IntVect& ratio);

// Fields for apply_tag_criteria (one level). Null pointers skip that criterion. zz is mass fraction when zz_den is null, else rho*Z with den = RHO.
struct TagFields {
    const amrex::MultiFab* tmp = nullptr;
    const amrex::MultiFab* rho = nullptr;
    const amrex::MultiFab* zz = nullptr;
    const amrex::MultiFab* zz_den = nullptr;
    int zz_comp = 0;
    const amrex::MultiFab* hrrpuv = nullptr;
    const amrex::MultiFab* vort = nullptr;   // cell-centred |ω| for TAG_VORT
    double tmp_ambient = 293.15;             // base for TAG_TMP_ABOVE (TMPA)
    std::array<int, 3> dirs{1, 1, 1};       // 0 in a single-cell direction
};

// OR of every set &AMR TAG_* criterion into `tags` (caller fills ghosts). keepfac from TAG_KEEP.
void apply_tag_criteria(amrex::TagBoxArray& tags, const AmrParams& p, const TagFields& f, const amrex::iMultiFab* covered);

}  // namespace fdsrt

#endif
