// Fields.H: MultiFab registry of the FDS mesh fields and the FDS <-> AMReX index maps (M2a, S2; IR-005).
//
// Kernel-facing rules (M2a):
//  (a) Passive scalars: N_TOTAL_SCALARS counts tracked species followed by passive scalars. They are handled HERE: the species arrays
//      ZZ and ZZS (rho*Z at the kernel stages, "RHO_ZZ") get ncomp = N_TOTAL_SCALARS components (component c <-> FDS index c+1), so a
//      passive scalar is one more component with the same ghost width and the same exchange. Nothing else in the driver special-cases it.
//      FX/FY/FZ and the other per-scalar scratch arrays are NOT registered (see is_per_box_scratch).
//  (b) Only uniform Cartesian metrics are used. R(I)/RRN(I) are dropped: the registry carries no metric arrays, and the index maps assume
//      a constant cell size per direction (Level0 rejects CYLINDRICAL and TRN* meshes, IR-002).
//
// Layout (D-031, inventory/mesh_fields.csv), one MultiFab per field on the level-0 BoxArray (box = FDS mesh):
//   cell-centred          RHO, RHOS ng=3 (D-031 native layout); TMP, ZZ, ZZS ng=2; H, HS, D, DS, ... ng=1 (as allocated by FDS)
//   face-centred (nodal)  U, V, W, US, VS, WS ng=1, nodal in their own direction
// FDS array bounds are reproduced exactly (ng == ng_fds) except RHO/RHOS, which carry one extra ghost layer (ng_fds = 2 is the FDS
// window inside the native FAB; the kernel binding of S3 exposes it by lower-bound remapping).
//
// Index map. Box valid cells: FDS I = 1..IBAR <-> AMReX lo..hi (I = a - lo + 1). Face-centred in direction d: FDS U(I) is the HIGH face
// of cell I, AMReX face a is the LOW face of cell a, so FDS U(I) <-> AMReX face lo+I (offset +1). FDS U(0) is the low face of the
// box (AMReX lo), U(IBAR) the high face (AMReX hi+1). FAB lower bound in FDS index: 1-ng (cell), -ng (nodal direction);
// extent: n + 2 ng (+1 in the nodal direction). With ng=1 these are exactly U(-1:IBP1,0:JBP1,0:KBP1) etc.
#ifndef FDS_FIELDS_H_
#define FDS_FIELDS_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "FdsAmr.H"

namespace fdsamr {

enum class Stag : int { Cell = 0, FaceX = 1, FaceY = 2, FaceZ = 3 };

struct FieldSpec {
    const char* name;
    Stag stag;
    int ng;            // ghost width of the native MultiFab
    int ng_fds;        // ghost width of the FDS allocation (init.f90); differs from ng only for RHO, RHOS
    bool per_scalar;   // ncomp = N_TOTAL_SCALARS (else 1)
    const char* fds_bounds;  // FDS declared bounds (documentation and inventory cross-check)
    bool nodal(int d) const { return static_cast<int>(stag) == d + 1; }
    int ncomp(int nscalars) const { return per_scalar ? nscalars : 1; }
};

const std::vector<FieldSpec>& field_table();
const FieldSpec* find_field(const std::string& name);       // nullptr if not registered
const std::vector<std::string>& default_field_names();      // fields allocated by default for the M2a step
bool is_per_box_scratch(const std::string& name);           // FX, FY, FZ and the other per-scalar scratch arrays

// ---- FDS <-> AMReX index maps (one constant offset per direction and box) ----
inline int to_amrex(const FieldSpec& s, int d, int lo, int I) { return lo + I - 1 + (s.nodal(d) ? 1 : 0); }
inline int to_fds(const FieldSpec& s, int d, int lo, int a) { return a - lo + 1 - (s.nodal(d) ? 1 : 0); }

struct FdsBounds {
    int lb[3];   // FDS lower bounds
    int ext[3];  // extents
};
// FDS bounds of the whole native FAB (ng ghost layers). cell_box = the box's valid CELL box (the FDS mesh)
inline FdsBounds fds_bounds(const FieldSpec& s, const amrex::Box& cell_box)
{
    FdsBounds b{};
    for (int d = 0; d < 3; ++d) {
        const int n = cell_box.length(d);
        b.lb[d] = s.nodal(d) ? -s.ng : 1 - s.ng;
        b.ext[d] = n + 2 * s.ng + (s.nodal(d) ? 1 : 0);
    }
    return b;
}
// FDS bounds of the array as FDS allocates it (ng_fds ghost layers): the window the unmodified kernels index. Equal to fds_bounds()
// except for RHO/RHOS, where the window is a strided sub-block of the native FAB.
inline FdsBounds fds_window(const FieldSpec& s, const amrex::Box& cell_box)
{
    FieldSpec t = s;
    t.ng = s.ng_fds;
    return fds_bounds(t, cell_box);
}
// FDS index range of the valid data (cells 1..n, faces 0..n)
inline void fds_valid_range(const FieldSpec& s, const amrex::Box& cell_box, int lo_i[3], int hi_i[3])
{
    for (int d = 0; d < 3; ++d) {
        lo_i[d] = s.nodal(d) ? 0 : 1;
        hi_i[d] = cell_box.length(d);
    }
}

// Column-major view of a FAB with FDS lower bounds: v(I,J,K,n), n = 1..ncomp (FDS scalar index). `base` points at element
// (lb, n=1); strides are in elements, so a window smaller than the FAB (RHO/RHOS) is a strided view.
struct FdsView {
    double* base = nullptr;
    int lb[3] = {0, 0, 0};
    int ext[3] = {0, 0, 0};
    long stride[3] = {1, 0, 0};
    long cstride = 0;
    int ncomp = 1;
    bool contiguous = true;
    double& operator()(int I, int J, int K, int n = 1) const
    {
        return base[(I - lb[0]) * stride[0] + (J - lb[1]) * stride[1] + (K - lb[2]) * stride[2] + (n - 1) * cstride];
    }
};
// window = true: bounds as FDS allocates the array (what the kernels see); false: the whole native FAB.
FdsView make_fds_view(const FieldSpec& s, amrex::FArrayBox& fab, const amrex::Box& cell_box, bool window = true);

// ---- the registry ----
class Fields {
public:
    // names: fields to allocate (default_field_names() if empty). nscalars = N_TOTAL_SCALARS.
    Fields(const Level0& l0, int nscalars, const std::vector<std::string>& names = {});

    bool has(const std::string& name) const { return m_mf.count(name) > 0; }
    amrex::MultiFab& operator[](const std::string& name);
    const amrex::MultiFab& operator[](const std::string& name) const;
    const FieldSpec& spec(const std::string& name) const;

    // Ghost fill between boxes of this level (neighbouring boxes and periodic images), faces, edges and corners: layout-independent.
    // Physical-boundary ghosts are not touched (wall BCs fill them). A face-only mode is deliberately not offered: with AMReX 26.09 the
    // FillBoundary `cross` flag is applied to remote (MPI) copies only, not to same-rank copies, so its result would depend on the
    // box-to-rank map (see README, "Ghost fill"). The D-031 redundant valid+1 clip may read edge and
    // corner ghosts of RHO/RHOS/ZZ, which the full fill provides.
    void fill_ghosts(const std::string& name) const;

    int nscalars() const { return m_ns; }
    long bytes() const;                       // allocated bytes summed over ranks' local FABs (this rank only)
    const Level0& level0() const { return m_l0; }
    std::vector<std::string> names() const;

private:
    const Level0& m_l0;
    int m_ns;
    std::map<std::string, std::unique_ptr<amrex::MultiFab>> m_mf;
};

}  // namespace fdsamr

#endif
