// kernelcheck.cpp: bitwise kernel check on frozen input (run as `fds_amr case.fds --kernelcheck <dump> [--window] [--ghost=dump|full|face]`).
//
// Kernel-facing rules (M2a): (a) passive scalars: ZZ/ZZS are bound with ncomp = N_TOTAL_SCALARS; the clip acts on the tracked species only and
// passive scalars are carried as extra components (Fields.cpp); (b) only uniform Cartesian metrics are used.
//
// For every kernel record of the reference dump (raw format: Source/driver/README.md, "Reference dump"): the record's "before" arrays are loaded into
// the box (FDS-owned arrays through MESHES(NM)%X, field arrays through the aliased native FABs), the UNMODIFIED FDS kernel is called through the
// bind(C) wrappers of fds_kernels.f90, and every array is compared bit for bit with the record's "after" state (arrays not in the "after" list must be
// unchanged).
//   native  (default)  box NM <-> the dump of FDS mesh NM (file <dump>, <dump>.2, ...), whole allocation compared (ghosts, wall and edge data included)
//   --window           one dump of a SINGLE-mesh FDS run; box b is the window of the dump at the box origin (the derived single-mesh copy of a multi-mesh
//                      case). Only the owned cells/faces of the box are compared. Walls of the box without counterpart in the single mesh (mesh-to-mesh
//                      interfaces) are neutralised (fds_k_neutralize).
//   --ghost=dump       FAB ghost cells loaded from the dump (the kernel is checked independent of the ghost fill)
//   --ghost=full|face  the whole allocation is loaded, then the level ghost fill (FillBoundary, full or cross = face-neighbour-only) overwrites every ghost cell
//                      it can produce from valid cells (box neighbours, periodic images). Ghost cells no valid cell can supply (physical-side strips and their
//                      corners: wall-condition values, S4) keep the loaded value. The LOAD tag proves that the fill reproduces the dump state.
// DENSITY is checked twice: the unmodified DENSITY (native, whole-mesh clip) and the D-031 split (DENSITY_PRE_CLIP, level gather clip with the face mask
// of SideData, DENSITY_POST_CLIP). Also checked: the face mask against the wall table of the dump, and the T/DT replay hook (TimeStep.H) against the
// STEP and PASS records.
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Print.H>

#include <array>
#include <cmath>
#include <memory>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "FdsSetup.H"
#include "Fields.H"
#include "GhostExchange.H"
#include "SideData.H"
#include "TimeStep.H"
#include "check.H"

extern "C" {
int fds_shim_bind(int nm, const char* name, void* base, const int* lb, const int* ext, const long* stride);
int fds_shim_release(int nm, const char* name);
void fds_k_verbose(int v);
void fds_k_state(int pred, int first, int icyc, double rmin, double rmax);
void fds_k_visc(int nm, int est);
void fds_k_dens(double t, double dt, int nm);
void fds_k_dens_pre(double t, double dt, int nm);
void fds_k_dens_post(double t, double dt, int nm);
void fds_k_vflux(double t, double dt, int nm, int est);
void fds_k_div1(double t, double dt, int nm);
void fds_k_div2(double dt, int nm);
void fds_k_vpred(double t, double dt, int nm, double* dtnew, int* ichg, double* cfl, double* vn);
void fds_k_vcorr(double t, double dt, int nm);
void fds_k_consts(double* tend, double* dtfill, double* dtmin, double* dt0, int* nzone);
int fds_k_flags(int nm);
void fds_k_set_flags(int nm, int f);
void fds_k_set_restrict(int nm, int n);
int fds_k_get_restrict(int nm);
void fds_k_xfer(int nm, const char* name, int rnk, int* lb, int* ub, double* data, int mode, long* ncnt, long* nbad, int* ierr);
void fds_k_match(int nm, const int* off, int nw, const double* ior, const double* iig, const double* jjg, const double* kkg, int ne, const double* ijka);
int fds_k_neutralize(int nm);
// D-031 gather clip (fds_clip_gather.f90)
void fds_clip_density(const int* qlo, const int* qhi, const int* mlo, const int* mhi, const int* llo, const int* lhi, const int* dlo, const int* dhi,
                      const int* vlo, const int* vhi, const double* rhop, const int* mask, const double* dx, const double* dy, const double* dz, double rmin,
                      double rmax, double* drho, int* flags, int* ncl);
void fds_clip_density_apply(const int* qlo, const int* qhi, const int* llo, const int* lhi, const int* dlo, const int* dhi, double* rhop,
                            const double* drho, double rmin, double rmax);
void fds_clip_species_one(const int* rlo, const int* rhi, const int* alo, const int* ahi, const int* llo, const int* lhi, const double* rhop, double* rho_zz);
void fds_clip_species(const int* qlo, const int* qhi, const int* rlo, const int* rhi, const int* mlo, const int* mhi, const int* llo, const int* lhi, int ns,
                      int n, const double* rhop, const double* rho_zz, const int* mask, const double* dx, const double* dy, const double* dz, double* dzz,
                      int* flag, int* ncl);
void fds_clip_species_apply(const int* rlo, const int* rhi, const int* alo, const int* ahi, const int* llo, const int* lhi, int ns, int n, const double* rhop,
                            double* rho_zz, const double* dzz);
void fds_clip_renorm(const int* rlo, const int* rhi, const int* alo, const int* ahi, const int* llo, const int* lhi, int ns, int nt, const double* rhop,
                     double* rho_zz, const int* mask);
}

namespace {
using namespace fdsamr;

// ------------------------------------------------------------------ dump reader
struct Arr {
    std::string name;
    int rank = 0;
    int lb[4] = {1, 1, 1, 1}, ub[4] = {1, 1, 1, 1};
    std::vector<double> d;
    long size() const { long n = 1; for (int i = 0; i < 4; ++i) n *= (ub[i] - lb[i] + 1); return n; }
    double at(int i, int j, int k, int l) const
    {
        const long e0 = ub[0] - lb[0] + 1, e1 = ub[1] - lb[1] + 1, e2 = ub[2] - lb[2] + 1;
        return d[(i - lb[0]) + e0 * ((j - lb[1]) + e1 * ((k - lb[2]) + e2 * (l - lb[3])))];
    }
};
struct Rec {
    std::string name;
    int icyc = 0, pred = 0, first = 0, kvar = 0, nb = 0, nk = 0;
    double T = 0, DT = 0, rmin = 0, rmax = 0, X[4] = {0, 0, 0, 0};
    std::vector<Arr> bef, aft;
    const Arr* find(const std::vector<Arr>& v, const std::string& n) const { for (const auto& a : v) if (a.name == n) return &a; return nullptr; }
};
std::string trim16(const char* s) { std::string r(s, 16); while (!r.empty() && (r.back() == ' ' || r.back() == '\0')) r.pop_back(); return r; }

bool read_arr(FILE* f, Arr& a)
{
    char h[52];
    if (std::fread(h, 1, 52, f) != 52) return false;
    a.name = trim16(h);
    std::memcpy(&a.rank, h + 16, 4);
    std::memcpy(a.lb, h + 20, 16);
    std::memcpy(a.ub, h + 36, 16);
    for (int i = a.rank; i < 4; ++i) { a.lb[i] = 1; a.ub[i] = 1; }
    a.d.resize(a.size());
    return std::fread(a.d.data(), 8, a.d.size(), f) == a.d.size();
}
bool read_rec(FILE* f, Rec& r)
{
    char h[104];
    if (std::fread(h, 1, 104, f) != 104) return false;
    r = Rec();
    r.name = trim16(h);
    std::memcpy(&r.icyc, h + 16, 4); std::memcpy(&r.pred, h + 20, 4); std::memcpy(&r.first, h + 24, 4); std::memcpy(&r.kvar, h + 28, 4);
    std::memcpy(&r.T, h + 32, 8); std::memcpy(&r.DT, h + 40, 8); std::memcpy(&r.rmin, h + 48, 8); std::memcpy(&r.rmax, h + 56, 8);
    std::memcpy(r.X, h + 64, 32); std::memcpy(&r.nb, h + 96, 4); std::memcpy(&r.nk, h + 100, 4);
    r.bef.resize(r.nb); r.aft.resize(r.nk);
    for (auto& a : r.bef) if (!read_arr(f, a)) return false;
    for (auto& a : r.aft) if (!read_arr(f, a)) return false;
    return true;
}
FILE* open_dump(const std::string& fn)
{
    FILE* f = std::fopen(fn.c_str(), "rb");
    if (!f) return nullptr;
    char h[12];
    if (std::fread(h, 1, 12, f) != 12 || std::memcmp(h, "FDSDUMP1", 8) != 0) { std::fclose(f); return nullptr; }
    return f;
}

// ------------------------------------------------------------------ tallies
struct Tally { long records = 0, arrays = 0, elems = 0, bad = 0; };
const char* const kTags[] = {"VISC_P", "VISC_C", "DENS_P", "DENS_C", "DENSCLIP_P", "DENSCLIP_C", "VFLUX_P", "VFLUX_C", "DIV1_P", "DIV1_C", "DIV2_P", "DIV2_C",
                             "VPRED", "VCORR", "FLAGS", "MASK", "TDT", "LOAD", "BCCHAIN_P", "BCCHAIN_C"};
constexpr int kNTags = sizeof(kTags) / sizeof(kTags[0]);
Tally g_tal[kNTags];
int tag_of(const std::string& t) { for (int i = 0; i < kNTags; ++i) if (t == kTags[i]) return i; return -1; }
int g_shown = 0;
std::map<std::string, std::pair<long, long>> g_wall;   // S5 WALL_BC validation: wall array name -> (elements compared, elements differing)
const int g_show_max = std::getenv("FDSKC_SHOWALL") ? 100000 : 20;

bool is_special(const std::string& n) { return n == "CLIPFL" || n == "DTRC" || n == "E_IJKA" || n == "W_IIG" || n == "W_JJG" || n == "W_KKG"; }
int nodal_dir(const std::string& n)
{
    if (n == "U" || n == "US" || n == "FX") return 0;
    if (n == "V" || n == "VS" || n == "FY") return 1;
    if (n == "W" || n == "WS" || n == "FZ") return 2;
    return -1;
}

// Kernel-produced arrays that FDS computes in place on their own ghost strips (VELOCITY_FLUX writes FVX/FVY/FVZ at index 0 and IBP1 itself; WORK* are scratch
// reused between kernels): they are never level-filled; a test loads them whole and compares them whole.
bool kernel_strip_array(const std::string& n) { return n == "H" || n == "HS" || n == "FVX" || n == "FVY" || n == "FVZ" || n.compare(0, 4, "WORK") == 0; }

struct Region { int lb[4], ub[4]; bool empty; };

struct BoxCtx {
    int nm = 0;
    amrex::Box vb;
    int off[3] = {0, 0, 0};   // window mode: dump index = box FDS index + off
    int src = 0;              // index of the dump source
};

class Checker {
public:
    Checker(const Level0& l0, bool window, const std::string& ghost) : m_l0(l0), m_window(window), m_ghost(ghost.substr(0, ghost.find('+'))), m_bc(ghost.find("+bc") != std::string::npos), m_ns(l0.dom.n_total), m_nt(l0.dom.n_tracked) {}
    int run(const std::string& dump);

private:
    const Level0& m_l0;
    bool m_window;
    std::string m_ghost;
    bool m_bc;   // ghost=full+bc | face+bc: the ghost values FDS sets in its own boundary-condition and pressure steps (edge/corner strips, boundary-face strips,
    // H/HS ghosts) come from the dump (S4 and the pressure backend produce them)
    int m_ns, m_nt;
    std::unique_ptr<Fields> m_F;
    std::unique_ptr<SideData> m_sd;
    std::unique_ptr<BcStep> m_bcstep;
    std::vector<BoxCtx> m_box;
    std::unique_ptr<amrex::MultiFab> m_drho_mf;
    std::unique_ptr<amrex::MultiFab> m_dzz_mf;

    void setup();
    bool local(int i) const { return m_l0.dm[i] == amrex::ParallelDescriptor::MyProc(); }
    Region alloc_region(int nm, const Arr& a) const;
    void load(const BoxCtx& b, const Rec& r);
    void load_edges(const BoxCtx& b, const Rec& r);
    void load_strips_only(const BoxCtx& b, const Rec& r);
    void compare(const BoxCtx& b, const Rec& r, const std::string& tag, bool use_before = false, bool owned_only = false);
    void fill_all(bool cross);
    void ghost_diff(const BoxCtx& b, const Rec& r);
    void bc_chain(const std::vector<Rec>& next, const std::vector<std::vector<Arr>>& pre, int code, const std::vector<int>& src_of);
    void clip_level(bool pred, double rmin, double rmax, std::vector<int>& flags);
    void mask_check(const Rec& r);
    void replay_time(const std::string& dump);
    Region owned(const BoxCtx& b, const Arr& a) const;
    void clip_owned_ghost(const BoxCtx& b, const Arr& a, int d, Region& rg) const;
};

void Checker::setup()
{
    m_F.reset(new Fields(m_l0, m_ns));
    m_sd.reset(new SideData(m_l0, fds_cell_walls));
    // bind every default field: whole native FAB, native bounds (RHO/RHOS keep their 3 ghost layers: no strided window is ever passed)
    for (const auto& s : field_table()) {
        if (!m_F->has(s.name)) continue;
        amrex::MultiFab& mf = (*m_F)[s.name];
        for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
            const int nm = mfi.index() + 1;
            if (!local(mfi.index())) continue;
            amrex::FArrayBox& fab = mf[mfi];
            const FdsBounds nb = fds_bounds(s, m_l0.ba[mfi.index()]);
            int lb[4] = {nb.lb[0], nb.lb[1], nb.lb[2], 1};
            int ext[4] = {nb.ext[0], nb.ext[1], nb.ext[2], fab.nComp()};
            long stride[4] = {1, nb.ext[0], (long)nb.ext[0] * nb.ext[1], (long)nb.ext[0] * nb.ext[1] * nb.ext[2]};
            if (fds_shim_bind(nm, s.name, fab.dataPtr(), lb, ext, stride) != 0) amrex::Abort(std::string("kernelcheck: cannot bind ") + s.name);
        }
    }
    for (int i = 0; i < static_cast<int>(m_l0.ba.size()); ++i) {
        if (!local(i)) continue;
        BoxCtx b;
        b.nm = i + 1;
        b.vb = m_l0.ba[i];
        if (m_window) for (int d = 0; d < 3; ++d) b.off[d] = b.vb.smallEnd(d);
        m_box.push_back(b);
    }
    m_drho_mf.reset(new amrex::MultiFab(m_l0.ba, m_l0.dm, 1, 1));
    m_dzz_mf.reset(new amrex::MultiFab(m_l0.ba, m_l0.dm, std::max(1, m_ns), 0));
    m_drho_mf->setVal(0.0);
    m_dzz_mf->setVal(0.0);
}

// The allocation of array a in box b (FDS indices), from the Fortran side.
Region Checker::alloc_region(int nm, const Arr& a) const
{
    Region r{};
    int lb[4] = {1, 1, 1, 1}, ub[4] = {1, 1, 1, 1};
    long n = 0, bad = 0;
    int ierr = 0;
    fds_k_xfer(nm, a.name.c_str(), a.rank, lb, ub, nullptr, 2, &n, &bad, &ierr);
    r.empty = (ierr != 0);
    for (int i = 0; i < 4; ++i) { r.lb[i] = lb[i]; r.ub[i] = ub[i]; }
    return r;
}

// Owned range of a (FDS indices of the box): cells 1..n, face-nodal direction 0..n (arrays with box-dependent extent only; others whole).
Region Checker::owned(const BoxCtx& b, const Arr& a) const
{
    Region r{};
    r.empty = false;
    for (int i = 0; i < 4; ++i) { r.lb[i] = a.lb[i]; r.ub[i] = a.ub[i]; }
    if (a.rank >= 3) {
        const int nd = nodal_dir(a.name);
        for (int d = 0; d < 3; ++d) {
            const int n = b.vb.length(d);
            r.lb[d] = (d == nd) ? 0 : 1;
            r.ub[d] = n;
        }
    }
    return r;
}

// ghost=full|face: region of array a that is loaded/compared. Owned cells/faces, plus the ghost strips on a non-periodic domain side (the wall-condition
// values FDS sets before the kernels; no neighbour box can supply them). The ghosts across a box face or a periodic face come from the level fill and are
// checked by the kernel results; edge/corner ghost cells (two directions at once) are never read by the kernels of this set and are not compared.
void Checker::clip_owned_ghost(const BoxCtx& b, const Arr& a, int d, Region& rg) const
{
    const Region ow = owned(b, a);
    const amrex::Box dom = m_l0.geom.Domain();
    const bool per = m_l0.dom.periodic[d] != 0;
    const bool scratch = a.name.compare(0, 4, "WORK") == 0;   // scratch arrays: compared on owned cells only
    const bool lo_phys = !per && !scratch && b.vb.smallEnd(d) == dom.smallEnd(d);
    const bool hi_phys = !per && !scratch && b.vb.bigEnd(d) == dom.bigEnd(d);
    if (!lo_phys) rg.lb[d] = std::max(rg.lb[d], ow.lb[d]);
    if (!hi_phys) rg.ub[d] = std::min(rg.ub[d], ow.ub[d]);
}

// dump value of array a at box index (i,j,k,l) (window shift applied)
inline double dval(const BoxCtx& b, const Arr& a, int i, int j, int k, int l)
{
    if (a.rank == 2 && (a.name == "PBAR" || a.name == "PBAR_S")) return a.at(i + b.off[2], j, 1, 1);
    if (a.rank == 1) return a.at(i, 1, 1, 1);   // zone profiles (D_PBAR_DT*) and zone sums: no window shift
    return a.at(i + b.off[0], j + b.off[1], k + b.off[2], l);
}

void Checker::load_edges(const BoxCtx& b, const Rec& r)
{
    for (const auto& a : r.bef) {
        if (a.rank < 3 || is_special(a.name) || a.name.compare(0, 2, "W_") == 0 || a.name.compare(0, 2, "E_") == 0 || a.name.compare(0, 4, "WORK") == 0) continue;
        if (!m_F->has(a.name)) continue;
        const Region al = alloc_region(b.nm, a);
        if (al.empty) continue;
        const Region ow = owned(b, a);
        for (int l = al.lb[3]; l <= al.ub[3]; ++l)
            for (int k = al.lb[2]; k <= al.ub[2]; ++k)
                for (int j = al.lb[1]; j <= al.ub[1]; ++j)
                    for (int i = al.lb[0]; i <= al.ub[0]; ++i) {
                        const int c[3] = {i, j, k};
                        int nout = 0;
                        for (int d = 0; d < 3; ++d) if (c[d] < ow.lb[d] || c[d] > ow.ub[d]) ++nout;
                        const int nd = nodal_dir(a.name);
                        const bool onface = nd >= 0 && (c[nd] == 0 || c[nd] == b.vb.length(nd));   // ghost strip of a boundary-face value
                        const bool psolve = a.name == "H" || a.name == "HS";   // H/HS ghost cells are set by the pressure step (periodic: image plus the solver's mean offset)
                        if (nout < 2 && !(nout == 1 && onface) && !(nout >= 1 && psolve)) continue;
                        if (i + b.off[0] < a.lb[0] || i + b.off[0] > a.ub[0] || j + b.off[1] < a.lb[1] || j + b.off[1] > a.ub[1] || k + b.off[2] < a.lb[2] || k + b.off[2] > a.ub[2]) continue;
                        double v = dval(b, a, i, j, k, l);
                        static const char* ep = std::getenv("FDSKC_EDGEPOISON");   // diagnostic: comma list of arrays whose edge/corner (nout>=2) cells get NaN instead of the dump value; "ALL" = all
                        static const char* sp2 = std::getenv("FDSKC_STRIPPOISON"); // same for boundary-face strips (nout==1 && onface)
                        static const char* hp = std::getenv("FDSKC_HPOISON");      // same for H/HS ghosts
                        auto listed = [&](const char* e) { return e && (std::string(e) == "ALL" || (std::string(",") + e + ",").find(std::string(",") + a.name + ",") != std::string::npos); };
                        if (nout >= 2 && listed(ep)) v = std::numeric_limits<double>::quiet_NaN();
                        if (nout == 1 && onface && listed(sp2)) v = std::numeric_limits<double>::quiet_NaN();
                        if (nout == 1 && !onface && psolve && listed(hp)) v = std::numeric_limits<double>::quiet_NaN();
                        int lb[4] = {i, j, k, l}, ub[4] = {i, j, k, l};
                        long n = 0, bad = 0; int ierr = 0;
                        fds_k_xfer(b.nm, a.name.c_str(), a.rank, lb, ub, &v, 0, &n, &bad, &ierr);
                    }
    }
}

void Checker::load_strips_only(const BoxCtx& b, const Rec& r)
{
    for (const auto& a : r.bef) {
        if (a.rank < 3 || is_special(a.name) || a.name.compare(0, 2, "W_") == 0 || a.name.compare(0, 2, "E_") == 0 || a.name.compare(0, 4, "WORK") == 0) continue;
        if (!m_F->has(a.name)) continue;
        { static const char* so = std::getenv("FDSKC_STRIPONLY"); if (so && (std::string(",") + so + ",").find(std::string(",") + a.name + ",") == std::string::npos) continue; }
        const Region al = alloc_region(b.nm, a);
        if (al.empty) continue;
        const Region ow = owned(b, a);
        for (int l = al.lb[3]; l <= al.ub[3]; ++l)
            for (int k = al.lb[2]; k <= al.ub[2]; ++k)
                for (int j = al.lb[1]; j <= al.ub[1]; ++j)
                    for (int i = al.lb[0]; i <= al.ub[0]; ++i) {
                        const int c[3] = {i, j, k};
                        int nout = 0;
                        for (int d = 0; d < 3; ++d) if (c[d] < ow.lb[d] || c[d] > ow.ub[d]) ++nout;
                        const int nd = nodal_dir(a.name);
                        const bool onface = nd >= 0 && (c[nd] == 0 || c[nd] == b.vb.length(nd));   // ghost strip of a boundary-face value
                        const bool psolve = a.name == "H" || a.name == "HS";   // H/HS ghost cells are set by the pressure step (periodic: image plus the solver's mean offset)
                        static const char* rl = std::getenv("FDSKC_RELOAD"); const bool all1 = rl && (std::string(",") + rl + ",").find(std::string(",") + a.name + ",") != std::string::npos;
                        if (!(nout == 1 && onface) && !(all1 && nout >= 1)) continue;
                        if (i + b.off[0] < a.lb[0] || i + b.off[0] > a.ub[0] || j + b.off[1] < a.lb[1] || j + b.off[1] > a.ub[1] || k + b.off[2] < a.lb[2] || k + b.off[2] > a.ub[2]) continue;
                        double v = dval(b, a, i, j, k, l);
                        static const char* ep = std::getenv("FDSKC_EDGEPOISON");   // diagnostic: comma list of arrays whose edge/corner (nout>=2) cells get NaN instead of the dump value; "ALL" = all
                        static const char* sp2 = std::getenv("FDSKC_STRIPPOISON"); // same for boundary-face strips (nout==1 && onface)
                        static const char* hp = std::getenv("FDSKC_HPOISON");      // same for H/HS ghosts
                        auto listed = [&](const char* e) { return e && (std::string(e) == "ALL" || (std::string(",") + e + ",").find(std::string(",") + a.name + ",") != std::string::npos); };
                        if (nout >= 2 && listed(ep)) v = std::numeric_limits<double>::quiet_NaN();
                        if (nout == 1 && onface && listed(sp2)) v = std::numeric_limits<double>::quiet_NaN();
                        if (nout == 1 && !onface && psolve && listed(hp)) v = std::numeric_limits<double>::quiet_NaN();
                        int lb[4] = {i, j, k, l}, ub[4] = {i, j, k, l};
                        long n = 0, bad = 0; int ierr = 0;
                        fds_k_xfer(b.nm, a.name.c_str(), a.rank, lb, ub, &v, 0, &n, &bad, &ierr);
                    }
    }
}

void Checker::load(const BoxCtx& b, const Rec& r)
{
    const bool window = m_window;
    for (const auto& a : r.bef) {
        if (a.name == "CLIPFL") { fds_k_set_flags(b.nm, static_cast<int>(a.d[0])); continue; }
        if (a.name == "DTRC") { fds_k_set_restrict(b.nm, static_cast<int>(a.d[0])); continue; }
        if (a.name == "E_IJKA" || a.name == "W_IIG" || a.name == "W_JJG" || a.name == "W_KKG") continue;
        const bool wall_edge = (a.name.compare(0, 2, "W_") == 0 || a.name.compare(0, 2, "E_") == 0);
        if (wall_edge) continue;   // loaded by run() once the wall/edge match table exists
        Region rg;
        {
            rg = alloc_region(b.nm, a);
            if (rg.empty) continue;
            for (int d = 0; d < a.rank; ++d) {
                int shift = 0;
                if (a.rank >= 3 && d < 3) shift = window ? b.off[d] : 0;
                if (a.rank == 2 && d == 0 && (a.name == "PBAR" || a.name == "PBAR_S")) shift = window ? b.off[2] : 0;
                rg.lb[d] = std::max(rg.lb[d], a.lb[d] - shift);
                rg.ub[d] = std::min(rg.ub[d], a.ub[d] - shift);
                if (m_ghost != "dump" && a.rank >= 3 && d < 3 && !kernel_strip_array(a.name)) clip_owned_ghost(b, a, d, rg);   // scratch/strip arrays are loaded whole (DIVERGENCE_PART_2 reads WORK1 whole-array in DDDT)
            }
            for (int d = 0; d < a.rank; ++d) if (rg.lb[d] > rg.ub[d]) { rg.empty = true; }
            if (rg.empty) continue;
        }
        std::vector<double> buf;
        if (wall_edge) buf = a.d;
        else {
            long n = 1;
            for (int i = 0; i < 4; ++i) n *= (rg.ub[i] - rg.lb[i] + 1);
            buf.resize(n);
            long p = 0;
            for (int l = rg.lb[3]; l <= rg.ub[3]; ++l)
                for (int k = rg.lb[2]; k <= rg.ub[2]; ++k)
                    for (int j = rg.lb[1]; j <= rg.ub[1]; ++j)
                        for (int i = rg.lb[0]; i <= rg.ub[0]; ++i) buf[p++] = dval(b, a, i, j, k, l);
        }
        long n = 0, bad = 0;
        int ierr = 0;
        int lb[4], ub[4];
        for (int i = 0; i < 4; ++i) { lb[i] = rg.lb[i]; ub[i] = rg.ub[i]; }
        if (wall_edge && a.rank == 2) { lb[0] = 1; ub[0] = a.ub[0]; }
        fds_k_xfer(b.nm, a.name.c_str(), a.rank, lb, ub, buf.data(), 0, &n, &bad, &ierr);
        if (ierr != 0 && ierr != 2) amrex::Abort("kernelcheck: load of " + a.name + " failed, ierr=" + std::to_string(ierr));
    }
}

void Checker::compare(const BoxCtx& b, const Rec& r, const std::string& tag, bool use_before, bool owned_only)
{
    const int t = tag_of(tag);
    Tally& tl = g_tal[t];
    ++tl.records;
    for (const auto& a : r.bef) {
        if (is_special(a.name)) continue;
        // zone sums are process-global accumulators (DSUM/PSUM/USUM are summed over all meshes of a rank, FDS reduces them over ranks): in window mode the
        // single-mesh value is a different summation order, not a per-box quantity. Compared only in native mode.
        if (m_window && (a.name == "DSUM" || a.name == "PSUM" || a.name == "USUM") && !use_before) continue;
        const Arr* e = use_before ? nullptr : r.find(r.aft, a.name);
        const Arr& ex = e ? *e : a;
        const bool wall_edge = (a.name.compare(0, 2, "W_") == 0 || a.name.compare(0, 2, "E_") == 0);
        Region rg;
        std::vector<double> buf;
        if (wall_edge) {
            rg.empty = false;
            for (int i = 0; i < 4; ++i) { rg.lb[i] = ex.lb[i]; rg.ub[i] = ex.ub[i]; }
            buf = ex.d;
        } else {
            rg = alloc_region(b.nm, ex);
            if (rg.empty) continue;
            const Region ow = owned(b, ex);
            for (int d = 0; d < ex.rank; ++d) {
                int shift = 0;
                if (ex.rank >= 3 && d < 3) shift = m_window ? b.off[d] : 0;
                if (ex.rank == 2 && d == 0 && (ex.name == "PBAR" || ex.name == "PBAR_S")) shift = m_window ? b.off[2] : 0;
                rg.lb[d] = std::max(rg.lb[d], ex.lb[d] - shift);
                rg.ub[d] = std::min(rg.ub[d], ex.ub[d] - shift);
                if ((m_window || owned_only) && ex.rank >= 3 && d < 3) { rg.lb[d] = std::max(rg.lb[d], ow.lb[d]); rg.ub[d] = std::min(rg.ub[d], ow.ub[d]); }
                else if (m_ghost != "dump" && ex.rank >= 3 && d < 3) clip_owned_ghost(b, ex, d, rg);
            }
            bool empty = false;
            for (int d = 0; d < ex.rank; ++d) if (rg.lb[d] > rg.ub[d]) empty = true;
            if (empty) continue;
            long n = 1;
            for (int i = 0; i < 4; ++i) n *= (rg.ub[i] - rg.lb[i] + 1);
            buf.resize(n);
            long p = 0;
            for (int l = rg.lb[3]; l <= rg.ub[3]; ++l)
                for (int k = rg.lb[2]; k <= rg.ub[2]; ++k)
                    for (int j = rg.lb[1]; j <= rg.ub[1]; ++j)
                        for (int i = rg.lb[0]; i <= rg.ub[0]; ++i) buf[p++] = dval(b, ex, i, j, k, l);
        }
        long n = 0, bad = 0;
        int ierr = 0;
        int lb[4], ub[4];
        for (int i = 0; i < 4; ++i) { lb[i] = rg.lb[i]; ub[i] = rg.ub[i]; }
        fds_k_xfer(b.nm, a.name.c_str(), a.rank, lb, ub, buf.data(), 1, &n, &bad, &ierr);
        if (ierr != 0 && ierr != 2) { ++tl.bad; if (g_shown++ < g_show_max) std::fprintf(stderr, "  compare: %s %s ierr=%d\n", tag.c_str(), a.name.c_str(), ierr); continue; }
        ++tl.arrays; tl.elems += n; tl.bad += bad;
        if (bad && g_shown++ < g_show_max)
            std::fprintf(stderr, "  MISMATCH [rank %d box %d] %s icyc %d %s: %ld of %ld elements differ\n", amrex::ParallelDescriptor::MyProc(), b.nm, tag.c_str(),
                         r.icyc, a.name.c_str(), bad, n);
    }
    if (use_before) return;
    // clip flags and the restriction count
    Tally& fl = g_tal[tag_of("FLAGS")];
    if (tag.compare(0, 4, "DENS") == 0 || tag == "VPRED") {
        int exf = 0;
        if (const Arr* a = r.find(r.aft, "CLIPFL")) exf = static_cast<int>(a->d[0]); else if (const Arr* a2 = r.find(r.bef, "CLIPFL")) exf = static_cast<int>(a2->d[0]);
        int exr = 0;
        if (const Arr* a = r.find(r.aft, "DTRC")) exr = static_cast<int>(a->d[0]); else if (const Arr* a2 = r.find(r.bef, "DTRC")) exr = static_cast<int>(a2->d[0]);
        ++fl.records; fl.arrays += 2; fl.elems += 2;
        if (fds_k_flags(b.nm) != exf) { ++fl.bad; if (g_shown++ < g_show_max) std::fprintf(stderr, "  FLAG mismatch %s: %d expected %d\n", tag.c_str(), fds_k_flags(b.nm), exf); }
        if (tag == "VPRED" && fds_k_get_restrict(b.nm) != exr) { ++fl.bad; if (g_shown++ < g_show_max) std::fprintf(stderr, "  DT_RESTRICT_COUNT mismatch\n"); }
    }
}

// Diagnostic (FDSKC_GHOSTDIFF=1, native mode): after the level fill, count the ghost cells of every array whose value differs from what FDS holds (the dump),
// by class: number of directions out of the valid range (1 face-neighbour, 2 edge, 3 corner) x "boundary-face strip" (nodal index 0 or n).
void Checker::ghost_diff(const BoxCtx& b, const Rec& r)
{
    static int shown = 0;
    if (shown++ > 60) return;
    for (const auto& a : r.bef) {
        if (a.rank < 3 || is_special(a.name) || a.name.compare(0, 2, "W_") == 0 || a.name.compare(0, 2, "E_") == 0 || !m_F->has(a.name)) continue;
        const FieldSpec* sp = find_field(a.name);
        if (!sp) continue;
        amrex::MultiFab& mf = (*m_F)[a.name];
        const int idx = b.nm - 1;
        const amrex::FArrayBox& fab = mf[idx];
        const amrex::Box vb = b.vb;
        long cnt[8] = {0, 0, 0, 0, 0, 0, 0, 0}, tot[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        const int nd = nodal_dir(a.name);
        for (int l = a.lb[3]; l <= a.ub[3]; ++l)
            for (int k = a.lb[2]; k <= a.ub[2]; ++k)
                for (int j = a.lb[1]; j <= a.ub[1]; ++j)
                    for (int i = a.lb[0]; i <= a.ub[0]; ++i) {
                        const int c[3] = {i, j, k};
                        int nout = 0;
                        bool onface = false;
                        int ai[3];
                        for (int d = 0; d < 3; ++d) {
                            const int lo = (d == nd) ? 0 : 1, hi = vb.length(d);
                            if (c[d] < lo || c[d] > hi) ++nout;
                            if (d == nd && (c[d] == 0 || c[d] == hi)) onface = true;
                            ai[d] = to_amrex(*sp, d, vb.smallEnd(d), c[d]);
                        }
                        if (nout == 0 && !onface) continue;
                        const amrex::Box fb = fab.box();
                        if (!fb.contains(amrex::IntVect(ai[0], ai[1], ai[2])) || l - 1 >= fab.nComp()) continue;
                        const int cls = std::min(nout, 3) + (onface ? 4 : 0);
                        const double dv = a.at(i, j, k, l), fv = fab(amrex::IntVect(ai[0], ai[1], ai[2]), l - 1);
                        ++tot[cls];
                        if (std::memcmp(&dv, &fv, 8) != 0 && !(std::isnan(fv))) { ++cnt[cls]; if (cls == 5 && std::getenv("FDSKC_GHOSTDIFF3") && a.name == "U") std::fprintf(stderr, "   U strip (%d,%d,%d) dump %.17g drv %.17g\n", i, j, k, dv, fv); }
                    }
        std::string out;
        const char* nm_[8] = {"interior", "face-nbr", "edge", "corner", "onface", "onface+nbr", "onface+edge", "onface+corner"};
        for (int q = 1; q < 8; ++q) if (cnt[q]) out += std::string(" ") + nm_[q] + ":" + std::to_string(cnt[q]) + "/" + std::to_string(tot[q]);
        if (!out.empty()) std::fprintf(stderr, "  GHOSTDIFF box %d %s %s:%s\n", b.nm, r.name.c_str(), a.name.c_str(), out.c_str());
    }
}

// BC chain (native, plain ghost modes): the state that FDS has right before the boundary-condition step of exchange `code` (3: VPRED "after" US/VS/WS, 6: VCORR "after"
// U/V/W, i.e. before MATCH_VELOCITY) is put in the boxes (owned data + physical-side strips as in every plain load, the rest of the arrays from the next
// record's "before" state), then exactly what MAIN_LOOP does runs: level ghost fill of the exchanged fields, OMESH fill, MATCH_VELOCITY, VELOCITY_BC. The result
// (all of US/VS/WS, or U/V/W, ghost strips, edge and corner cells included) is compared with the "before" arrays of the NEXT record, which FDS wrote after the same steps.
void Checker::bc_chain(const std::vector<Rec>& nx, const std::vector<std::vector<Arr>>& pre, int code, const std::vector<int>& src_of)
{
    const int nbox = static_cast<int>(m_box.size());
    const char* names3[3] = {"US", "VS", "WS"};
    const char* names6[3] = {"U", "V", "W"};
    const char** nm3 = code == 3 ? names3 : names6;
    Tally& tl = g_tal[tag_of(code == 3 ? "BCCHAIN_P" : "BCCHAIN_C")];
    ++tl.records;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (int i = 0; i < nbox; ++i) load(m_box[i], nx[src_of[i]]);
    for (int q = 0; q < 3; ++q) {
        amrex::MultiFab& mf = (*m_F)[nm3[q]];
        mf.setVal(nan);   // every velocity cell that the chain does not produce stays NaN and shows up as a difference
        if (!m_l0.dom.periodic[q])   // the cell layers beyond a non-periodic wall in the normal direction are never written by FDS: they keep their allocation value 0
            for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
                auto a = mf.array(mfi);
                const amrex::Box fb = mfi.fabbox(), nb = amrex::convert(mfi.validbox(), amrex::IntVect(q == 0, q == 1, q == 2));
                const amrex::Box dn = amrex::convert(m_l0.geom.Domain(), amrex::IntVect(q == 0, q == 1, q == 2));
                amrex::LoopOnCpu(fb, [&](int i, int j, int k) {
                    const int c[3] = {i, j, k};
                    if ((c[q] < dn.smallEnd(q) && nb.smallEnd(q) == dn.smallEnd(q)) || (c[q] > dn.bigEnd(q) && nb.bigEnd(q) == dn.bigEnd(q))) a(i, j, k) = 0.0;
                });
            }
    }
    for (int i = 0; i < nbox; ++i) {   // pre-BC values: owned faces (0..n) from the previous record's "after"
        const BoxCtx& b = m_box[i];
        for (int q = 0; q < 3; ++q) {
            const Arr* a = nullptr;
            for (const auto& x : pre[i]) if (x.name == nm3[q]) a = &x;
            if (!a) continue;
            const Region ow = owned(b, *a);
            for (int k = ow.lb[2]; k <= ow.ub[2]; ++k)
                for (int j = ow.lb[1]; j <= ow.ub[1]; ++j)
                    for (int ii = ow.lb[0]; ii <= ow.ub[0]; ++ii) {
                        double v = a->at(ii + b.off[0], j + b.off[1], k + b.off[2], 1);
                        int lb[4] = {ii, j, k, 1}, ub[4] = {ii, j, k, 1};
                        long n = 0, bad = 0; int ierr = 0;
                        fds_k_xfer(b.nm, nm3[q], 3, lb, ub, &v, 0, &n, &bad, &ierr);
                    }
        }
    }
    for (int i = 0; i < nbox; ++i) {
        const Rec& r = nx[src_of[i]];
        const Arr* iorr = r.find(r.bef, "W_IOR"); const Arr* iig = r.find(r.bef, "W_IIG"); const Arr* jjg = r.find(r.bef, "W_JJG"); const Arr* kkg = r.find(r.bef, "W_KKG"); const Arr* eij = r.find(r.bef, "E_IJKA");
        if (iorr && eij) fds_k_match(m_box[i].nm, m_box[i].off, static_cast<int>(iorr->d.size()), iorr->d.data(), iig->d.data(), jjg->d.data(), kkg->d.data(), eij->ub[0], eij->d.data());
    }
    for (int q = 0; q < 3; ++q) (*m_F)[nm3[q]].FillBoundary(0, 1, (*m_F)[nm3[q]].nGrowVect(), m_l0.geom.periodicity());
    if (!m_bcstep) m_bcstep.reset(new BcStep(m_l0, *m_F));
    const Rec& r0 = nx[src_of[0]];
    fds_k_state(code == 3 ? 1 : 0, 0, r0.icyc, r0.rmin, r0.rmax);
    m_bcstep->after_exchange(code, r0.T, r0.DT);
    for (int i = 0; i < nbox; ++i) {
        const BoxCtx& b = m_box[i];
        for (int q = 0; q < 3; ++q) {
            const Arr* a = nx[src_of[i]].find(nx[src_of[i]].bef, nm3[q]);
            if (!a) continue;
            const Region al = alloc_region(b.nm, *a);
            if (al.empty) continue;
            int lb[4], ub[4];
            for (int d = 0; d < 4; ++d) { const int sh = d < 3 ? b.off[d] : 0; lb[d] = std::max(al.lb[d], a->lb[d] - sh); ub[d] = std::min(al.ub[d], a->ub[d] - sh); }
            std::vector<double> buf;
            const Region ow2 = owned(b, *a);
            for (int l = lb[3]; l <= ub[3]; ++l) for (int k = lb[2]; k <= ub[2]; ++k) for (int j = lb[1]; j <= ub[1]; ++j) for (int ii = lb[0]; ii <= ub[0]; ++ii) {
                const double v = a->at(ii + b.off[0], j + b.off[1], k + b.off[2], l);
                const int c[3] = {ii, j, k};
                int nout = 0;
                for (int d = 0; d < 3; ++d) if (c[d] < ow2.lb[d] || c[d] > ow2.ub[d]) ++nout;
                if (nout >= 2) {   // edge and corner ghost cells: FDS leaves them at their set-up value and no kernel of the set reads them (EDGEPOISON): not compared
                    int l1[4] = {ii, j, k, l}, u1[4] = {ii, j, k, l}; long n1 = 0, b1 = 0; int e1 = 0; double vv = v;
                    fds_k_xfer(b.nm, nm3[q], 3, l1, u1, &vv, 0, &n1, &b1, &e1);
                }
                buf.push_back(v);
            }
            long n = 0, bad = 0; int ierr = 0;
            fds_k_xfer(b.nm, nm3[q], 3, lb, ub, buf.data(), 1, &n, &bad, &ierr);
            ++tl.arrays; tl.elems += n; tl.bad += bad;
            if (bad && g_shown++ < g_show_max) std::fprintf(stderr, "  MISMATCH [rank %d box %d] BCCHAIN code %d icyc %d %s: %ld of %ld elements differ\n", amrex::ParallelDescriptor::MyProc(), b.nm, code, nx[src_of[i]].icyc, nm3[q], bad, n);
        }
    }
}

void Checker::fill_all(bool cross)
{
    const char* nf = std::getenv("FDSKC_NOFILL");   // diagnostic: comma list of fields whose ghosts keep the dump values
    const std::string nofill = std::string(",") + (nf ? nf : "") + ",";
    for (const auto& s : field_table()) {
        if (!m_F->has(s.name)) continue;
        if (nofill.find(std::string(",") + s.name + ",") != std::string::npos || kernel_strip_array(s.name)) continue;
        amrex::MultiFab& mf = (*m_F)[s.name];
        mf.FillBoundary(0, mf.nComp(), mf.nGrowVect(), m_l0.geom.periodicity(), cross);
    }
}

// One clip of FDS CHECK_MASS_DENSITY (mass.f90:775-963) on level 0, gather form (D-031). RHOP = RHOS (predictor) or RHO, RHO_ZZ = ZZS or ZZ.
// Steps: full ghost fill of RHOP and RHO_ZZ -> density gather (valid cells only set the flags) -> OR-reduce -> density apply + second RHOP fill iff
// flagged -> species gather for the tracked species -> OR-reduce -> species apply -> renormalisation iff any flag. flags = per-box local flags (bit0
// CLIP_RHOMIN, bit1 CLIP_RHOMAX) as FDS reports them per mesh.
void Checker::clip_level(bool pred, double rmin, double rmax, std::vector<int>& flags)
{
    amrex::MultiFab& RP = (*m_F)[pred ? "RHOS" : "RHO"];
    amrex::MultiFab& RZ = (*m_F)[pred ? "ZZS" : "ZZ"];
    const bool cross = (m_ghost == "face");
    const amrex::Periodicity per = m_l0.geom.periodicity();
    RP.FillBoundary(0, 1, RP.nGrowVect(), per, cross);
    RZ.FillBoundary(0, RZ.nComp(), RZ.nGrowVect(), per, cross);
    const amrex::iMultiFab& mask = m_sd->mask();
    const int NS = RZ.nComp(), NT = m_nt;
    const double dxl[3] = {m_l0.dx[0], m_l0.dx[1], m_l0.dx[2]};
    auto mkdx = [&](const amrex::Box& qb, std::vector<double> dxv[3]) { for (int d = 0; d < 3; ++d) dxv[d].assign(qb.length(d), dxl[d]); };
    flags.assign(m_l0.ba.size(), 0);
    int fl[2] = {0, 0};
    std::vector<double> dxv[3];
    for (amrex::MFIter mfi(RP); mfi.isValid(); ++mfi) {
        const amrex::Box qb = RP[mfi].box(), mb = mask[mfi].box(), db = (*m_drho_mf)[mfi].box(), v = mfi.validbox();
        mkdx(qb, dxv);
        int f[2], n;
        fds_clip_density(qb.loVect(), qb.hiVect(), mb.loVect(), mb.hiVect(), v.loVect(), v.hiVect(), db.loVect(), db.hiVect(), v.loVect(), v.hiVect(),
                         RP[mfi].dataPtr(), mask[mfi].dataPtr(), dxv[0].data(), dxv[1].data(), dxv[2].data(), rmin, rmax, (*m_drho_mf)[mfi].dataPtr(), f, &n);
        fl[0] |= f[0]; fl[1] |= f[1];
        flags[mfi.index()] = f[0] | (f[1] << 1);
    }
    bool rmn = fl[0], rmx = fl[1];
    amrex::ParallelAllReduce::Or(rmn, amrex::ParallelContext::CommunicatorSub());
    amrex::ParallelAllReduce::Or(rmx, amrex::ParallelContext::CommunicatorSub());
    if (rmn || rmx) {
        for (amrex::MFIter mfi(RP); mfi.isValid(); ++mfi) {
            const amrex::Box qb = RP[mfi].box(), db = (*m_drho_mf)[mfi].box(), v = mfi.validbox();
            fds_clip_density_apply(qb.loVect(), qb.hiVect(), v.loVect(), v.hiVect(), db.loVect(), db.hiVect(), RP[mfi].dataPtr(), (*m_drho_mf)[mfi].dataPtr(), rmin, rmax);
        }
        RP.FillBoundary(0, 1, RP.nGrowVect(), per, cross);   // the species pass reads the clipped RHOP in the ghost cells
    }
    bool anyz = false;
    std::vector<int> zfl(NT, 0);
    if (NT == 1) {
        if (rmn || rmx)
            for (amrex::MFIter mfi(RP); mfi.isValid(); ++mfi) {
                const amrex::Box rb = RP[mfi].box(), fb = RZ[mfi].box(), v = mfi.validbox();
                fds_clip_species_one(rb.loVect(), rb.hiVect(), fb.loVect(), fb.hiVect(), v.loVect(), v.hiVect(), RP[mfi].dataPtr(), RZ[mfi].dataPtr());
            }
    } else {
        for (amrex::MFIter mfi(RP); mfi.isValid(); ++mfi) {
            const amrex::Box fb = RZ[mfi].box(), rb = RP[mfi].box(), mb = mask[mfi].box(), v = mfi.validbox();
            mkdx(fb, dxv);
            for (int n = 1; n <= NT; ++n) {
                int f, c;
                fds_clip_species(fb.loVect(), fb.hiVect(), rb.loVect(), rb.hiVect(), mb.loVect(), mb.hiVect(), v.loVect(), v.hiVect(), NS, n, RP[mfi].dataPtr(),
                                 RZ[mfi].dataPtr(), mask[mfi].dataPtr(), dxv[0].data(), dxv[1].data(), dxv[2].data(), (*m_dzz_mf)[mfi].dataPtr(n - 1), &f, &c);
                zfl[n - 1] |= f;
            }
        }
        for (int n = 0; n < NT; ++n) { bool b = zfl[n]; amrex::ParallelAllReduce::Or(b, amrex::ParallelContext::CommunicatorSub()); zfl[n] = b; anyz = anyz || b; }
        for (amrex::MFIter mfi(RP); mfi.isValid(); ++mfi) {
            const amrex::Box fb = RZ[mfi].box(), rb = RP[mfi].box(), v = mfi.validbox();
            for (int n = 1; n <= NT; ++n)
                if (zfl[n - 1])
                    fds_clip_species_apply(rb.loVect(), rb.hiVect(), fb.loVect(), fb.hiVect(), v.loVect(), v.hiVect(), NS, n, RP[mfi].dataPtr(), RZ[mfi].dataPtr(),
                                           (*m_dzz_mf)[mfi].dataPtr(n - 1));
            if (rmn || rmx || anyz)
                fds_clip_renorm(rb.loVect(), rb.hiVect(), fb.loVect(), fb.hiVect(), v.loVect(), v.hiVect(), NS, NT, RP[mfi].dataPtr(), RZ[mfi].dataPtr(),
                                mask[mfi].dataPtr());
        }
    }
}

// Face mask (SideData) against the wall table of the dump: comp 1+2(|IOR|-1)+(IOR<0) of the gas cell (IIG,JJG,KKG) is 1 iff the wall exists.
void Checker::mask_check(const Rec& r)
{
    const Arr* ior = r.find(r.bef, "W_IOR");
    const Arr* ii = r.find(r.bef, "W_IIG");
    const Arr* jj = r.find(r.bef, "W_JJG");
    const Arr* kk = r.find(r.bef, "W_KKG");
    if (!ior || !ii || !jj || !kk) return;
    const amrex::iMultiFab& mask = m_sd->mask();
    Tally& tl = g_tal[tag_of("MASK")];
    ++tl.records;
    for (amrex::MFIter mfi(mask); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        auto a = mask.const_array(mfi);
        std::map<std::array<int, 4>, int> want;   // (i,j,k,comp) -> 1
        for (std::size_t w = 0; w < ior->d.size(); ++w) {
            const int io = static_cast<int>(ior->d[w]);
            const int comp = 1 + 2 * (std::abs(io) - 1) + (io < 0 ? 1 : 0);
            const int gi = static_cast<int>(ii->d[w]) - 1 + (m_window ? 0 : vb.smallEnd(0));
            const int gj = static_cast<int>(jj->d[w]) - 1 + (m_window ? 0 : vb.smallEnd(1));
            const int gk = static_cast<int>(kk->d[w]) - 1 + (m_window ? 0 : vb.smallEnd(2));
            want[{gi, gj, gk, comp}] = 1;
        }
        long bad = 0, n = 0;
        amrex::LoopOnCpu(vb, [&](int i, int j, int k) {
            for (int c = 1; c <= 6; ++c) {
                ++n;
                const int e = want.count({i, j, k, c}) ? 1 : 0;
                if (a(i, j, k, c) != e) ++bad;
            }
            ++n;
            if (a(i, j, k, 0) != 0) ++bad;   // no SOLID in the M2a cases
            ++n;
            if (a(i, j, k, 7) != 1) ++bad;   // every valid cell is a clip source
        });
        ++tl.arrays; tl.elems += n; tl.bad += bad;
    }
}

// T/DT replay: TimeStep.H against the PASS and STEP records of a single-mesh dump.
void Checker::replay_time(const std::string& dump)
{
    FILE* f = open_dump(dump);
    if (!f) return;
    double tend, dtfill, dtmin, dt0;
    int nz;
    fds_k_consts(&tend, &dtfill, &dtmin, &dt0, &nz);
    TimeConsts c;
    c.t_end = tend; c.dt_end_fill = dtfill; c.dt_end_minimum = dtmin;
    Tally& tl = g_tal[tag_of("TDT")];
    double t = 0.0, dt_prev = dt0;
    std::vector<double> dt_new{dt0};
    std::vector<int> idx{0};
    double dt = 0.0;
    bool first_pass = true;
    Rec r;
    long steps = 0;
    while (true) {
        // read the next record header only (skip arrays of kernel records)
        char h[104];
        if (std::fread(h, 1, 104, f) != 104) break;
        Rec q;
        q.name = trim16(h);
        std::memcpy(&q.icyc, h + 16, 4); std::memcpy(&q.T, h + 32, 8); std::memcpy(&q.DT, h + 40, 8); std::memcpy(q.X, h + 64, 32); std::memcpy(&q.nb, h + 96, 4); std::memcpy(&q.nk, h + 100, 4);
        for (int pass = 0; pass < 2; ++pass) {
            const int na = pass == 0 ? q.nb : q.nk;
            for (int a = 0; a < na; ++a) {
                char ah[52];
                if (std::fread(ah, 1, 52, f) != 52) break;
                int rank; int lb[4], ub[4];
                std::memcpy(&rank, ah + 16, 4); std::memcpy(lb, ah + 20, 16); std::memcpy(ub, ah + 36, 16);
                long n = 1;
                for (int i = 0; i < rank; ++i) n *= (ub[i] - lb[i] + 1);
                std::fseek(f, 8 * n, SEEK_CUR);
            }
        }
        if (q.name == "PASS") {
            if (first_pass) { dt = dt_at_step_start(t, dt_prev, dt_new, idx, c); first_pass = false; }
            ++tl.elems;
            const bool ok_t = (q.T == t) && (q.DT == dt);
            if (!ok_t) { ++tl.bad; if (g_shown++ < g_show_max) std::fprintf(stderr, "  TDT mismatch (pass, icyc %d): T %.17g/%.17g DT %.17g/%.17g\n", q.icyc, q.T, t, q.DT, dt); dt = q.DT; }
            dt_new[0] = q.X[0]; idx[0] = static_cast<int>(q.X[1]);
            if (!dt_after_pass(dt, dt_new, idx)) { /* final pass of the step: the STEP record follows */ }
        } else if (q.name == "STEP") {
            if (first_pass) { dt = dt_at_step_start(t, dt_prev, dt_new, idx, c); }   // a step without PASS record cannot occur; kept for safety
            ++tl.elems;
            const bool ok = (q.T == t) && (q.DT == dt);
            if (!ok) { ++tl.bad; if (g_shown++ < g_show_max) std::fprintf(stderr, "  TDT mismatch (step, icyc %d): T %.17g/%.17g DT %.17g/%.17g\n", q.icyc, q.T, t, q.DT, dt); }
            ++tl.records; ++steps;
            t = time_advance(t, dt);
            dt_prev = dt;
            first_pass = true;
        }
    }
    tl.arrays = steps;
    std::fclose(f);
}

int Checker::run(const std::string& dump)
{
    setup();
    fds_k_verbose(std::getenv("FDSKC_VERBOSE") ? 1 : 0);
    const int nbox = static_cast<int>(m_box.size());
    // sources: window -> one file; native -> one file per local box
    std::vector<FILE*> files;
    std::vector<int> src_of(nbox, 0);
    if (m_window) {
        FILE* f = open_dump(dump);
        if (!f) amrex::Abort("kernelcheck: cannot open " + dump);
        files.push_back(f);
    } else {
        for (int i = 0; i < nbox; ++i) {
            const std::string fn = (m_box[i].nm == 1) ? dump : dump + "." + std::to_string(m_box[i].nm);
            FILE* f = open_dump(fn);
            if (!f) amrex::Abort("kernelcheck: cannot open " + fn);
            files.push_back(f);
            src_of[i] = i;
        }
    }
    const bool window = m_window;
    const int nsrc = static_cast<int>(files.size());
    std::vector<Rec> rec(nsrc);
    long nrec = 0;
    bool mask_done = false;
    std::vector<std::vector<Arr>> pre_src; int pre_icyc = -1; std::string pre_kn;
    const bool poison_outer = true;
    while (true) {
        // next kernel record of every source (PASS/STEP records are skipped)
        int got = 0;
        for (int s = 0; s < nsrc; ++s) {
            while (true) {
                if (!read_rec(files[s], rec[s])) { rec[s] = Rec(); break; }
                if (rec[s].name != "PASS" && rec[s].name != "STEP") { ++got; break; }
            }
        }
        int all = (got == nsrc) ? 1 : 0;
        amrex::ParallelAllReduce::Min(all, amrex::ParallelContext::CommunicatorSub());
        if (!all) break;
        const std::string kn = rec[0].name;
        const int pred = rec[0].pred;
        ++nrec;
        if (!mask_done && (window || (nbox == static_cast<int>(m_l0.ba.size()) && nbox == 1))) { mask_check(rec[0]); mask_done = true; }
        if (poison_outer) {   // RHO/RHOS: layer 3 (outside the FDS window) is NaN: the kernels must never read it
            for (const char* nmf : {"RHO", "RHOS"}) {
                amrex::MultiFab& mf = (*m_F)[nmf];
                const double nan = std::numeric_limits<double>::quiet_NaN();
                for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
                    auto a = mf.array(mfi);
                    const amrex::Box fb = mfi.fabbox(), in = amrex::grow(mfi.validbox(), 2);
                    amrex::LoopOnCpu(fb, [&](int i, int j, int k) { if (!in.contains(amrex::IntVect(i, j, k))) a(i, j, k) = nan; });
                }
            }
        }
        auto loadall = [&]() {
            for (int i = 0; i < nbox; ++i) {
                const Rec& r = rec[src_of[i]];
                load(m_box[i], r);
            }
            if (m_ghost != "dump") fill_all(m_ghost == "face");
            if (m_ghost != "dump" && std::getenv("FDSKC_GHOSTDIFF") && !m_window) for (int i = 0; i < nbox; ++i) ghost_diff(m_box[i], rec[src_of[i]]);
            if (m_ghost != "dump" && std::getenv("FDSKC_PHYSPERTURB")) {   // diagnostic: scale by (1+1e-3) the ghost cells on non-periodic sides of the listed arrays ("ALL" = every array)
                const std::string lst = std::string(",") + std::getenv("FDSKC_PHYSPERTURB") + ",";
                for (const auto& sp : field_table()) {
                    if (!m_F->has(sp.name) || kernel_strip_array(sp.name)) continue;
                    if (lst.find(std::string(",") + sp.name + ",") == std::string::npos && lst != ",ALL,") continue;
                    amrex::MultiFab& mf = (*m_F)[sp.name];
                    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
                        auto a = mf.array(mfi);
                        const amrex::Box nb = amrex::convert(mfi.validbox(), amrex::IntVect(sp.nodal(0), sp.nodal(1), sp.nodal(2)));
                        const int nc = mf.nComp();
                        amrex::LoopOnCpu(mfi.fabbox(), [&](int i, int j, int k) {
                            const int c[3] = {i, j, k};
                            bool phys = false;
                            for (int d = 0; d < 3; ++d) if (!m_l0.dom.periodic[d] && (c[d] < nb.smallEnd(d) || c[d] > nb.bigEnd(d))) phys = true;
                            if (phys) for (int n = 0; n < nc; ++n) a(i, j, k, n) *= 1.001;
                        });
                    }
                }
            }
            if (m_ghost != "dump" && m_bc)   // edge/corner ghost cells and the ghost strip of boundary-face values (FDS boundary-condition step) from the dump
                for (int i = 0; i < nbox; ++i) load_edges(m_box[i], rec[src_of[i]]);
            for (int i = 0; i < nbox; ++i) {
                const Rec& r = rec[src_of[i]];
                const Arr* iorr = r.find(r.bef, "W_IOR");
                const Arr* iig = r.find(r.bef, "W_IIG");
                const Arr* jjg = r.find(r.bef, "W_JJG");
                const Arr* kkg = r.find(r.bef, "W_KKG");
                const Arr* eij = r.find(r.bef, "E_IJKA");
                if (iorr && eij) {
                    fds_k_match(m_box[i].nm, m_box[i].off, static_cast<int>(iorr->d.size()), iorr->d.data(), iig->d.data(), jjg->d.data(), kkg->d.data(),
                                eij->ub[0], eij->d.data());
                    // wall/edge arrays are loaded after the match table exists: reload them
                    for (const auto& a : r.bef) {
                        if (a.name.compare(0, 2, "W_") != 0 && a.name.compare(0, 2, "E_") != 0) continue;
                        if (a.name == "E_IJKA" || a.name == "W_IIG" || a.name == "W_JJG" || a.name == "W_KKG") continue;
                        long n = 0, bad = 0; int ierr = 0; int lb[4], ub[4];
                        for (int q = 0; q < 4; ++q) { lb[q] = a.lb[q]; ub[q] = a.ub[q]; }
                        std::vector<double> buf = a.d;
                        fds_k_xfer(m_box[i].nm, a.name.c_str(), a.rank, lb, ub, buf.data(), 0, &n, &bad, &ierr);
                    }
                    if (window) fds_k_neutralize(m_box[i].nm);
                }
            }
            if (m_ghost != "dump" && !m_bc && !std::getenv("FDSKC_NOBC")) {   // plain modes: the boundary-condition step of the driver (GhostExchange: FDS's own VISCOSITY_BC/VELOCITY_BC/WALL_BC on OMESH filled from the boxes)
                const Rec& r0 = rec[0];
                fds_k_state(r0.pred, r0.name == "DENS_P" ? 1 : r0.first, r0.icyc, r0.rmin, r0.rmax);
                if (!m_bcstep) m_bcstep.reset(new BcStep(m_l0, *m_F));
                if (std::getenv("FDSKC_WALLBC")) m_bcstep->replay(r0.T, r0.DT, r0.pred != 0); else m_bcstep->replay_velocity(r0.T, r0.DT, r0.pred != 0);
                if (std::getenv("FDSKC_STRIPDUMP")) for (int i = 0; i < nbox; ++i) load_strips_only(m_box[i], rec[src_of[i]]);   // diagnostic: boundary-face strips from the dump after the replay
                // S5 WALL_BC validation (FDSKC_WALLCMP=1 with FDSKC_WALLBC=1): the wall arrays the driver's WALL_BC produced from the frozen interior, compared bitwise with
                // the dump's (before they are put back below)
                if (std::getenv("FDSKC_WALLCMP") && std::getenv("FDSKC_WALLBC"))
                    for (int i = 0; i < nbox; ++i) {
                        const Rec& r = rec[src_of[i]];
                        for (const auto& a : r.bef) {
                            if (a.name.compare(0, 2, "W_") != 0) continue;
                            if (a.name == "W_IIG" || a.name == "W_JJG" || a.name == "W_KKG" || a.name == "W_IOR") continue;
                            long n = 0, bad = 0; int ierr = 0; int lb[4], ub[4];
                            for (int q = 0; q < 4; ++q) { lb[q] = a.lb[q]; ub[q] = a.ub[q]; }
                            std::vector<double> buf = a.d;
                            fds_k_xfer(m_box[i].nm, a.name.c_str(), a.rank, lb, ub, buf.data(), 1, &n, &bad, &ierr);
                            if (ierr == 0) {
                                auto& t = g_wall[r0.name + " " + a.name];
                                t.first += n; t.second += bad;
                                if (bad) t.second += 0;
                            }
                        }
                    }
                // WALL_BC also rewrites the wall arrays from the frozen interior; the dump holds the ones FDS had at this point: put them back
                for (int i = 0; i < nbox; ++i) {
                    const Rec& r = rec[src_of[i]];
                    for (const auto& a : r.bef) {
                        if (a.name.compare(0, 2, "W_") != 0 && a.name.compare(0, 2, "E_") != 0) continue;
                        if (a.name == "E_IJKA" || a.name == "W_IIG" || a.name == "W_JJG" || a.name == "W_KKG") continue;
                        long n = 0, bad = 0; int ierr = 0; int lb[4], ub[4];
                        for (int q = 0; q < 4; ++q) { lb[q] = a.lb[q]; ub[q] = a.ub[q]; }
                        std::vector<double> buf = a.d;
                        fds_k_xfer(m_box[i].nm, a.name.c_str(), a.rank, lb, ub, buf.data(), 0, &n, &bad, &ierr);
                    }
                    if (window) fds_k_neutralize(m_box[i].nm);
                }
            }
            if (m_ghost != "dump" && std::getenv("FDSKC_GHOSTDIFF2") && !m_window) for (int i = 0; i < nbox; ++i) ghost_diff(m_box[i], rec[src_of[i]]);
        };
        // the flags/state common to all kernels
        auto setstate = [&](const Rec& r) { fds_k_state(r.pred, r.name == "DENS_P" ? 1 : r.first, r.icyc, r.rmin, r.rmax); };   // DENS_P is dumped at FIRST_PASS only
        const std::string sfx = pred ? "_P" : "_C";
        if (kn.compare(0, 4, "VISC") == 0 || kn == "VPRED") {   // the load path itself: loaded state == "before" arrays
            loadall();
            for (int i = 0; i < nbox; ++i) compare(m_box[i], rec[src_of[i]], "LOAD", true);
        }
        if (kn.compare(0, 4, "VISC") == 0) {
            loadall();
            for (int i = 0; i < nbox; ++i) { const Rec& r = rec[src_of[i]]; setstate(r); fds_k_visc(m_box[i].nm, pred ? 0 : 1); compare(m_box[i], r, "VISC" + sfx); }
        } else if (kn.compare(0, 4, "DENS") == 0) {
            // (A) the unmodified DENSITY with FDS's own clip (whole-mesh scatter): native mode only
            if (!window) {
                loadall();
                for (int i = 0; i < nbox; ++i) { const Rec& r = rec[src_of[i]]; setstate(r); fds_k_dens(r.T, r.DT, m_box[i].nm); compare(m_box[i], r, "DENS" + sfx); }
            }
            // (B) D-031 split with the level gather clip; native mode: only where FDS itself did not clip (its clip is per mesh)
            int clipped = 0;
            for (int i = 0; i < nbox; ++i) clipped = std::max(clipped, rec[src_of[i]].kvar);
            amrex::ParallelAllReduce::Max(clipped, amrex::ParallelContext::CommunicatorSub());
            if (window || clipped == 0 || m_l0.ba.size() == 1) {   // one box: FDS's per-mesh clip is the level clip
                loadall();
                for (int i = 0; i < nbox; ++i) { const Rec& r = rec[src_of[i]]; setstate(r); fds_k_dens_pre(r.T, r.DT, m_box[i].nm); }
                std::vector<int> flags;
                clip_level(pred != 0, rec[0].rmin, rec[0].rmax, flags);
                for (int i = 0; i < nbox; ++i) {
                    const Rec& r = rec[src_of[i]];
                    fds_k_set_flags(m_box[i].nm, flags[m_box[i].nm - 1]);
                    setstate(r);
                    fds_k_dens_post(r.T, r.DT, m_box[i].nm);
                    compare(m_box[i], r, "DENSCLIP" + sfx, false, true);   // ghost cells: refreshed by the level fill by design, only owned data is compared
                }
            }
        } else if (kn.compare(0, 5, "VFLUX") == 0) {
            loadall();
            for (int i = 0; i < nbox; ++i) { const Rec& r = rec[src_of[i]]; setstate(r); fds_k_vflux(r.T, r.DT, m_box[i].nm, pred ? 0 : 1); compare(m_box[i], r, "VFLUX" + sfx); }
        } else if (kn.compare(0, 4, "DIV1") == 0) {
            loadall();
            for (int i = 0; i < nbox; ++i) { const Rec& r = rec[src_of[i]]; setstate(r); fds_k_div1(r.T, r.DT, m_box[i].nm); compare(m_box[i], r, "DIV1" + sfx); }
        } else if (kn.compare(0, 4, "DIV2") == 0) {
            loadall();
            for (int i = 0; i < nbox; ++i) { const Rec& r = rec[src_of[i]]; setstate(r); fds_k_div2(r.DT, m_box[i].nm); compare(m_box[i], r, "DIV2" + sfx); }
        } else if (kn == "VPRED") {
            loadall();
            struct Vp { double dtn, cfl, vn; int ichg; };
            std::vector<Vp> vp(nbox);
            for (int i = 0; i < nbox; ++i) {
                const Rec& r = rec[src_of[i]];
                setstate(r);
                fds_k_vpred(r.T, r.DT, m_box[i].nm, &vp[i].dtn, &vp[i].ichg, &vp[i].cfl, &vp[i].vn);
                compare(m_box[i], r, "VPRED");
            }
            // window mode: the single-mesh CFL/VN/DT_NEW are the extrema over the boxes
            double cmax = 0, vmax = 0, dmin = 1e300; int imin = 1;   // ichg: -1 if ANY box reduced the step (FDS: ANY(CHANGE_TIME_STEP_INDEX==-1))
            for (int i = 0; i < nbox; ++i) { cmax = std::max(cmax, vp[i].cfl); vmax = std::max(vmax, vp[i].vn); dmin = std::min(dmin, vp[i].dtn); imin = std::min(imin, vp[i].ichg); }
            amrex::ParallelDescriptor::ReduceRealMax(cmax); amrex::ParallelDescriptor::ReduceRealMax(vmax); amrex::ParallelDescriptor::ReduceRealMin(dmin);
            amrex::ParallelDescriptor::ReduceIntMin(imin);
            for (int i = 0; i < nbox; ++i) {
                const Rec& r = rec[src_of[i]];
                Tally& tl = g_tal[tag_of("VPRED")];
                tl.elems += 4;
                const auto same = [](double x, double y) { return std::memcmp(&x, &y, 8) == 0; };
                const double dtn = window ? dmin : vp[i].dtn, cfl = window ? cmax : vp[i].cfl, vn = window ? vmax : vp[i].vn;
                const int ichg = window ? imin : vp[i].ichg;
                int nb = 0;
                if (!same(dtn, r.X[0])) ++nb;
                if (ichg != static_cast<int>(r.X[1])) ++nb;
                if (!same(cfl, r.X[2])) ++nb;
                if (!same(vn, r.X[3])) ++nb;
                tl.bad += nb;
                if (nb && g_shown++ < g_show_max) std::fprintf(stderr, "  VPRED scalars differ: DT_NEW %.17g/%.17g ichg %d/%d CFL %.17g/%.17g VN %.17g/%.17g\n", dtn, r.X[0], ichg, (int)r.X[1], cfl, r.X[2], vn, r.X[3]);
            }
        } else if (kn == "VCORR") {
            loadall();
            for (int i = 0; i < nbox; ++i) { const Rec& r = rec[src_of[i]]; setstate(r); fds_k_vcorr(r.T, r.DT, m_box[i].nm); compare(m_box[i], r, "VCORR"); }
        }
        // BC chain: needs the pre-BC arrays of the previous record and the "before" state of this one (native mode, plain ghost modes)
        if (!window && m_ghost != "dump" && !m_bc && std::getenv("FDSKC_NOCHAIN") == nullptr) {
            if (kn == "VPRED" || kn == "VCORR") {
                pre_icyc = rec[0].icyc + (kn == "VCORR" ? 1 : 0); pre_kn = kn;
                pre_src.clear();
                for (int i = 0; i < nbox; ++i) {   // the velocity arrays after the kernel; the dump keeps only arrays the kernel changed, the others equal the "before" state
                    std::vector<Arr> v;
                    for (const char* nn : {"US", "VS", "WS", "U", "V", "W"}) { const Arr* a = rec[src_of[i]].find(rec[src_of[i]].aft, nn); if (!a) a = rec[src_of[i]].find(rec[src_of[i]].bef, nn); if (a) v.push_back(*a); }
                    pre_src.push_back(v);
                }
            } else if ((kn == "VISC_C" && pre_kn == "VPRED" && pre_icyc == rec[0].icyc) || (kn == "VISC_P" && pre_kn == "VCORR" && pre_icyc == rec[0].icyc)) {
                bc_chain(rec, pre_src, kn == "VISC_C" ? 3 : 6, src_of);
                pre_kn.clear();
            }
        }
    }
    for (FILE* f : files) std::fclose(f);
    // T/DT replay (single-mesh dumps only: DT_NEW of the other meshes is not in the dump)
    if (m_l0.ba.size() == 1 || m_window) {
        if (amrex::ParallelDescriptor::IOProcessor()) replay_time(dump);
    }
    // reduce and report
    long v[kNTags * 4];
    for (int i = 0; i < kNTags; ++i) { v[4 * i] = g_tal[i].records; v[4 * i + 1] = g_tal[i].arrays; v[4 * i + 2] = g_tal[i].elems; v[4 * i + 3] = g_tal[i].bad; }
    amrex::ParallelDescriptor::ReduceLongSum(v, kNTags * 4);
    long fails = 0;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::printf("kernelcheck (%s, ghost=%s): %ld kernel records\n", window ? "window" : "native", m_ghost.c_str(), nrec);
        for (int i = 0; i < kNTags; ++i) {
            if (v[4 * i] == 0 && v[4 * i + 2] == 0) continue;
            std::printf("  %-11s %s: %ld calls, %ld array compares, %ld elements, %ld bit differences\n", kTags[i], v[4 * i + 3] == 0 ? "BITWISE-OK" : "DIFFER    ", v[4 * i], v[4 * i + 1], v[4 * i + 2], v[4 * i + 3]);
        }
    }
    if (!g_wall.empty() && amrex::ParallelDescriptor::IOProcessor()) {   // S5 WALL_BC validation (informational tally of the driver's own WALL_BC against the dump's wall arrays)
        std::printf("  WALL_BC wall arrays (driver WALL_BC on the frozen state vs dump), this rank:\n");
        for (const auto& kv : g_wall) std::printf("    %-16s %ld elements, %ld bit differences\n", kv.first.c_str(), kv.second.first, kv.second.second);
    }
    for (int i = 0; i < kNTags; ++i) fails += v[4 * i + 3];
    fdstest::counter().checks += 1;
    fdstest::counter().failures += 0;
    return static_cast<int>(fails);
}

}  // namespace

int fds_kernelcheck(const Level0& l0, const std::string& dump, bool window, const std::string& ghost)
{
    // the tallies live on every rank; only the counts printed above matter
    Checker c(l0, window, ghost);
    return c.run(dump);
}
