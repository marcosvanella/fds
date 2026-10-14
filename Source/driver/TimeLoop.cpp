// TimeLoop.cpp: see TimeLoop.H.
//
// Kernel-facing rules (M2a):
//  (a) Passive scalars: N_TOTAL_SCALARS beyond the tracked species are handled by Fields.cpp (ZZ/ZZS carry ncomp = N_TOTAL_SCALARS); every copy, exchange and clip
//      here loops over all components, the D-031 clip special-cases only the tracked species.
//  (b) Only uniform Cartesian metrics are used (R(I)/RRN(I) dropped); CYLINDRICAL and TRN* meshes are rejected (IR-002), the pressure call-through assumes IPS=0.
#include "ObstList.H"
#include "TimeLoop.H"
#include "LevelRegistry.H"
#include "PressureBcMap.H"
#include "CompositeOps.H"
#include "FluxStageRunner.H"
#include "DriverAdapter.H"
#include "PostRegridProjection.H"
#include "IlGhost.H"

#include <AMReX_BLassert.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelContext.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <algorithm>
#include <vector>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/resource.h>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <initializer_list>

#include "ExactSum.H"
#include "FdsSetup.H"
#include "PressureIface.H"

extern "C" {
int fds_shim_bind(int nm, const char* name, void* base, const int* lb, const int* ext, const long* stride);
void fds_k_state(int pred, int first, int icyc, double rmin, double rmax);
void fds_k_visc(int nm, int est);
void fds_k_xfer(int nm, const char* name, int rnk, int* lb, int* ub, double* data, int mode, long* ncnt, long* nbad, int* ierr);
void fds_k_dens_post(double t, double dt, int nm);
void fds_k_vflux(double t, double dt, int nm, int est);
void fds_k_div1(double t, double dt, int nm);
void fds_k_rad(double t, int nm);   // COMPUTE_RADIATION (+ COMBUSTION_BC), any level
void fds_k_rte_factor();
void fds_rad_ghost_set(int nicr, int nw, int nmmax, const int* w_nm, const int* w_ijk, const int* w_nic, const int* s_nm, const int* s_ijk, const int* s_loc, double** ils, double** ilr, int* nra, int* nsb);   // fds_kernels.f90
void fds_k_div2(double dt, int nm);
void fds_k_comb(double t, double dt, int nm);
void fds_k_vpred(double t, double dt, int nm, double* dtnew, int* ichg, double* cfl, double* vn);
void fds_k_vcorr(double t, double dt, int nm);
void fds_k_consts(double* tend, double* dtfill, double* dtmin, double* dt0, int* nzone);
void fds_k_set_flags(int nm, int f);
void fds_g_match_flux(int nm);
// end-of-step FDS outputs (patch 0006: fds_setup(mode=3) -> UPDATE_GLOBAL_OUTPUTS, EXCHANGE/DUMP_GLOBAL_OUTPUTS, WRITE_DIAGNOSTICS) and the hook module
void fds_setup(int mode, const char* fname, double* dt_out);
int fds_hook_step_outputs();
void fds_hook_set_step(double t, double dt, int icyc);
void fds_hook_set_mesh_dumps(int flag);  // patch 0011: per-mesh dump loop in fds_setup(mode=3), default off
// D-031 gather clip (fds_clip_gather.f90)
void fds_clip_density(const int* qlo, const int* qhi, const int* mlo, const int* mhi, const int* llo, const int* lhi, const int* dlo, const int* dhi,
                      const int* vlo, const int* vhi, const double* rhop, const int* mask, const double* dx, const double* dy, const double* dz, double rmin,
                      double rmax, double* drho, int* flags, int* ncl);
void fds_clip_density_apply(const int* qlo, const int* qhi, const int* llo, const int* lhi, const int* dlo, const int* dhi, double* rhop, const double* drho,
                            double rmin, double rmax);
void fds_clip_species_one(const int* rlo, const int* rhi, const int* alo, const int* ahi, const int* llo, const int* lhi, const double* rhop, double* rho_zz);
void fds_clip_species(const int* qlo, const int* qhi, const int* rlo, const int* rhi, const int* mlo, const int* mhi, const int* llo, const int* lhi, int ns,
                      int n, const double* rhop, const double* rho_zz, const int* mask, const double* dx, const double* dy, const double* dz, double* dzz,
                      int* flag, int* ncl);
void fds_clip_species_apply(const int* rlo, const int* rhi, const int* alo, const int* ahi, const int* llo, const int* lhi, int ns, int n, const double* rhop,
                            double* rho_zz, const double* dzz);
void fds_clip_renorm(const int* rlo, const int* rhi, const int* alo, const int* ahi, const int* llo, const int* lhi, int ns, int nt, const double* rhop,
                     double* rho_zz, const int* mask);
// fds_step.f90
void fds_p_params(int* ip, double* rp);
void fds_p_mesh_info(int nm, int* ip, double* rp);
void fds_p_get(int nm, const char* name, double* buf, long ntot, int* ierr);
void fds_p_mfd(int nm);
int fds_p_dens_skip();   // fds_step.f90: 1 when DENSITY returns at once (ICYC<=1 etc.): no level clip then
void fds_p_zero_dot();
void fds_p_esums(int nm, const int* mask, double* s);
void fds_p_add_qdot(int n, double v);
double fds_p_enthalpy_sum(int nm);
void fds_p_get_qr(int nm, double* buf, long ntot);
void fds_p_dens_pre(double t, double dt, int nm);
void fds_p_init_div();
void fds_p_zone_get(int n, double* d, double* p, double* u);
void fds_p_zone_set(int n, const double* d, const double* p, const double* u);
void fds_p_wall_dump(int nm, int ofx, int ofz, int kg0);
void fds_p_wall_table_dump(int nm, int ofx, int ofz);   // SCRATCH PATCH
void fds_hook_set_flag(int flag);
void fds_p_edge_dump(int nm, const char* fn, int nc);
void fds_p_iface_walls(int nm, int mode, const int* edge, int all_edges);
void fds_p_zone_terms(int nm, int* czone, double* dterm, double* pterm, int* wzone, double* wterm);
void fds_p_zone_volume_terms(int nm, int* czone, double* vterm);
double fds_p_zone_volume(int iz);
void fds_p_zone_volume_set(int iz, double v);
int fds_p_vent_count(int nm);
void fds_p_vent_get(int nm, int nv, int* ip, double* rp);
void fds_get_face_bc(int nm, int* code);   // fds_mesh_query.f90 (patch 0015)
void fds_get_face_wall_types(int nm, int f, int n1, int n2, int* t);   // fds_mesh_query.f90 (mixed face)
void fds_p_iter_init();
void fds_p_iter_inc();
int fds_p_iter_baro();
void fds_p_set_iter_baro(int f);
void fds_p_clear_attached(int nm);
void fds_p_baroclinic(double t, int nm);

void fds_p_noflux(double dt, int nm, int zero);
void fds_p_rhs(double t, double dt, int nm);
void fds_p_get_prhs(int nm, double* buf);
void fds_p_bmax(int nm, const int* flags, double* v);
void fds_p_get_bcface(int nm, int face, double* buf, long nexpect, int* ierr);
void fds_p_h_ghost(int nm, int pred, const int* flags);
void fds_p_resid(int nm);
void fds_p_velerr(double dt, int nm);
void fds_p_get_err(int nm, double* pe, double* ve, double* ptb);
void fds_p_set_wall_counter(int n);
// fds_flux_hooks.f90 (notes/flux-hooks-design.md)
void fds_flux_reserve(int ntot);
int fds_fine_b_set_view(int l, int ib, int which, const int* lb, const int* ext, void* p);   // fds_box_obj.f90: 0 ok, 3 no fine-level support in this tree
#ifdef FDS_HAVE_FINE_LEVEL
int fds_fine_level_create(int lev, int nbox, const int* ijk, const double* xb, const int* local, int* nm0);   // fds_fine_level.f90
void fds_fine_set_solid(int lev, int ib, const int* s);   // fds_fine_level.f90
int fds_fine_level_destroy(int lev);
#endif
int fds_wseam_refresh(int nm);
void fds_wseam_counts(int nm, int* next, int* nint);
int fds_wseam_upload(int nm);
int fds_wseam_check(int nm);
int fds_wseam_poke(int nm);
int fds_flux_register(int nm, int kind, int dir, const int* lb, const int* ext, double* p);
int fds_flux_set_mode(int nm, int kind, int mode);
int fds_flux_set_override(int nm, int kind, int n, const int* dir, const int* idx, int nsc, const double* val);
int fds_p_stop_status();
}

namespace fdsamr {

namespace {
inline bool is_local(const Level0& l0, int i) { return l0.dm[i] == amrex::ParallelDescriptor::MyProc(); }
template <class T> void all_max(T& v) { amrex::ParallelAllReduce::Max(v, amrex::ParallelContext::CommunicatorSub()); }
void die(const std::string& msg) { amrex::Abort("time loop: " + msg); }
}  // namespace

struct TimeLoop::Impl {
    TimeLoop& L;
    const Level0& l0;
    int IP[32];
    double RP[8];
    TimeConsts tc;
    int nbox = 0, ns = 0, nt = 0, nzone = 0;
    std::vector<int> mip;           // per box: 16 ints of fds_p_mesh_info (local boxes only)
    std::vector<double> dt_new;     // per box
    std::vector<int> chg;           // per box
    std::unique_ptr<amrex::MultiFab> drho, dzz, rhs, phi;
    int wall_counter = 0;
    bool first_pass = true;
    bool mms_done = false;
    int it_last_pred = 0, it_last_corr = 0;
    double perr = 0, verr = 0;
    double zone_rel_step = 0;
    // P3-C11/C12: per-step clip budget [L*nt + (N-1)]; pred+corr accumulate; reset at advance() start
    std::vector<double> clip_mass_acc;
    // Species budget of the D-031 clip (rho*Z of each tracked species moved by clip_level: density clip redistribution and the renormalisation into the largest
    // species are not species-conservative, as in FDS). Uncovered cells; the predictor's change enters the corrector with weight 1/2 (only its last pass counts).
    std::vector<double> clip_species_acc, clip_species_pred;
    std::vector<double> clip_species_sums(bool pred)
    {
        const amrex::MultiFab& RZ = (*c->F)[pred ? "ZZS" : "ZZ"];
        const int lev = static_cast<int>(c - lc.data());
        const amrex::iMultiFab* cov = L.m_reg ? L.m_reg->covered_mask(lev) : nullptr;
        const double vol = c->lev->dx[0] * c->lev->dx[1] * c->lev->dx[2];
        std::vector<double> m(nt, 0.0);
        for (amrex::MFIter mfi(RZ); mfi.isValid(); ++mfi) {
            const auto z = RZ.const_array(mfi);
            amrex::Array4<const int> cv;
            if (cov) cv = cov->const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { if (cov && cv(i, j, k) != 0) return; for (int n = 0; n < nt; ++n) m[n] += z(i, j, k, n) * vol; });
        }
        amrex::ParallelAllReduce::Sum(m.data(), nt, amrex::ParallelContext::CommunicatorSub());
        return m;
    }
    std::vector<long> clip_cells_acc;
    void ensure_clip_budget()
    {
        const int nl = std::max(1, static_cast<int>(lc.size()));
        const int nsp = std::max(1, nt);
        const std::size_t need = static_cast<std::size_t>(nl) * static_cast<std::size_t>(nsp);
        if (clip_mass_acc.size() < need) { clip_mass_acc.resize(need, 0.0); clip_cells_acc.resize(need, 0); }
    }
    void reset_clip_budget() { ensure_clip_budget(); std::fill(clip_mass_acc.begin(), clip_mass_acc.end(), 0.0); std::fill(clip_cells_acc.begin(), clip_cells_acc.end(), 0); }
    int clip_slot(int lev, int n0) const { return lev * nt + n0; }
    int passes = 0;
    int bc_type[6] = {0, 0, 0, 0, 0, 0};   // pb::BC per domain face (side*3+dir)
    FaceBc vel_bc{};                        // wall type per domain face (2*dir+side) for the fine-level velocity ghosts (BcStep::vel_bc, patch 0015)
    std::array<std::vector<int>, 6> face_wtype;   // mixed domain face: level-0 wall type per face cell (BcStep::face_wtype)
    bool solve_dir[3] = {true, true, true};
    std::string dir;
    std::vector<std::string> log_steps;
    std::vector<std::string> log_mass;
    // M2a measurements (NFR-030, A-22): wall seconds of this rank; the reported values are the maximum over ranks
    double t_loop = 0, t_pres = 0, t_solve = 0, t_out = 0;
    bool fds_outputs = false;       // FDS's own CHID_devc/_hrr/_mass/_steps/.out writers are driven (main.f90 patch 0006 present)

    explicit Impl(TimeLoop& l) : L(l), l0(l.m_l0) { std::memset(IP, 0, sizeof IP); std::memset(RP, 0, sizeof RP); }

    // ---- per-level context (S9). The stage functions below act on the CURRENT level c (select(lev)); level 0 is the default and the only bound level today. The code that is
    // level-0 only (set-up, pressure scheme, zone sums, outputs) keeps using l0 / L.m_F directly and is documented as such (notes/level-interface.md section 2).
    struct LevelCtx {
        const Level* lev = nullptr;
        Fields* F = nullptr;
        BcStep* bc = nullptr;
        SideData* sd = nullptr;
        amrex::MultiFab* drho = nullptr;   // scratch of the D-031 clip
        amrex::MultiFab* dzz = nullptr;
        int nbox = 0;
        int fds0 = 0;                      // FDS mesh number of box i is fds0 + i + 1
        int off = 0;                       // first slot of this level in the global per-box vectors (dt_new, chg)
    };
    std::vector<LevelCtx> lc;
    LevelCtx* c = nullptr;
    void select(int lev)
    {
        if (lev < 0 || lev >= static_cast<int>(lc.size())) {
            const bool known = L.m_reg && L.m_reg->has_level(lev);
            die("level " + std::to_string(lev) + (known ? " exists in the registry (layout, Fields, SideData) but has no FDS binding: " : " is not a bound level: ") +
                "the kernels take an FDS mesh number and a fine-level FDS MESH_TYPE object (metrics RDX..RDZN, wall tables, CELL_INDEX, POINT_TO_BOX route of patch 0005 DRAFT) does not exist yet (notes/level-interface.md section 3)");
        }
        c = &lc[lev];
    }
    int bi(int nm) const { return nm - 1 - c->fds0; }   // box index of FDS mesh number nm on the current level
    template <class Fn> void for_levels(Fn f)            // every bound level, coarse to fine; the current level is level 0 again afterwards
    {
        for (int l = 0; l < static_cast<int>(lc.size()); ++l) { select(l); f(l); }
        select(0);
    }
    int total_boxes() const { return lc.empty() ? 0 : lc.back().off + lc.back().nbox; }

    // ---- S12 level binding (TimeLoop::bind_level): a registry level > 0 gets FDS fine-box mesh objects (fds_fine_level.f90, patch 0007), views onto its Fields, a BcStep and the
    // scratch MultiFabs of the D-031 clip, and becomes level `lev` of lc. Only the top bound level can be unbound / rebound (mesh numbers of the levels above would move).
    struct OwnLevel { std::unique_ptr<amrex::MultiFab> drho, dzz; std::unique_ptr<BcStep> bc; };
    std::vector<OwnLevel> own;                                   // own[lev-1]
    static int view_code(const std::string& n)
    {
        static const char* names[23] = {"U", "V", "W", "US", "VS", "WS", "D", "DS", "H", "HS", "KRES", "FVX", "FVY", "FVZ", "RHO", "RHOS", "MU", "TMP", "Q", "RSUM", "ZZ", "ZZS", "DDDT"};
        for (int i = 0; i < 23; ++i) if (n == names[i]) return i + 1;
        return 0;
    }
    void bind_level(int lev)
    {
#ifndef FDS_HAVE_FINE_LEVEL
        die("bind_level(" + std::to_string(lev) + "): this build has no fine-level mesh objects (patch 0007 is not applied to Source/mesh.f90; see notes/level-interface.md section 3)");
#else
        if (lev < 1) die("bind_level: level 0 is bound at construction");
        if (!L.m_reg || !L.m_reg->has_level(lev)) die("bind_level(" + std::to_string(lev) + "): the level is not in the registry (make_level first)");
        if (lev > static_cast<int>(lc.size())) die("bind_level(" + std::to_string(lev) + "): levels are bound in order, the next one is " + std::to_string(lc.size()));
        if (lev < static_cast<int>(lc.size()) - 1) die("bind_level(" + std::to_string(lev) + "): only the top bound level can be rebound");
        if (lev == static_cast<int>(lc.size()) - 1) unbind_level(lev);
        const Level& Lv = L.m_reg->level(lev);
        Fields& F = L.m_reg->fields(lev);
        const int nb = static_cast<int>(Lv.ba.size());
        const int me = amrex::ParallelDescriptor::MyProc();
        std::vector<int> ijk(3 * nb), loc(nb);
        std::vector<double> xb(6 * nb);
        for (int i = 0; i < nb; ++i) {
            const amrex::Box b = Lv.ba[i];
            for (int d = 0; d < 3; ++d) {
                ijk[3 * i + d] = b.length(d);
                xb[6 * i + 2 * d] = Lv.geom.ProbLo(d) + b.smallEnd(d) * Lv.dx[d];
                xb[6 * i + 2 * d + 1] = Lv.geom.ProbLo(d) + (b.bigEnd(d) + 1) * Lv.dx[d];
            }
            loc[i] = (Lv.dm[i] == me) ? 1 : 0;
        }
        int nm0 = 0;
        {
            const int ierr = fds_fine_level_create(lev, nb, ijk.data(), xb.data(), loc.data(), &nm0);
            if (ierr != 0) die("bind_level(" + std::to_string(lev) + "): fds_fine_level_create failed (code " + std::to_string(ierr) + (ierr == 3 ? ": the level-0 meshes of this rank have more than one pressure zone" : "") + ")");
        }
        for (const auto& sp : field_table()) {
            const int which = view_code(sp.name);
            if (which == 0 || !F.has(sp.name)) continue;
            amrex::MultiFab& mf = F[sp.name];
            for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
                const amrex::Box cb = Lv.ba[mfi.index()];
                const FdsBounds nbd = fds_bounds(sp, cb);
                const int lb[4] = {nbd.lb[0], nbd.lb[1], nbd.lb[2], 1};
                const int ext[4] = {nbd.ext[0], nbd.ext[1], nbd.ext[2], mf[mfi].nComp()};
                if (fds_fine_b_set_view(lev, mfi.index() + 1, which, lb, ext, mf[mfi].dataPtr()) != 0)
                    die(std::string("bind_level: cannot bind the view of ") + sp.name + " on level " + std::to_string(lev));
            }
        }
        {   // CELL%SOLID of the fine cells from SideData comp 0
            const amrex::iMultiFab& sm = L.m_reg->side_data(lev).mask();
            for (amrex::MFIter mfi(sm); mfi.isValid(); ++mfi) fds_fine_set_solid(lev, mfi.index() + 1, sm[mfi].dataPtr());
        }
        own.resize(lev);
        OwnLevel& o = own[lev - 1];
        o.drho.reset(new amrex::MultiFab(Lv.ba, Lv.dm, 1, 1));
        o.dzz.reset(new amrex::MultiFab(Lv.ba, Lv.dm, std::max(1, ns), 0));
        o.drho->setVal(0.0); o.dzz->setVal(0.0);
        o.bc.reset(new BcStep(Lv, F));
        if (!std::getenv("FDSTL_LEGACY_BCSTEP_OFFSET")) o.bc->set_mesh_offset(nm0);   // (FDSTL_LEGACY_BCSTEP_OFFSET=1: offset 0 as before S14.4, reproducer of the level-0 ghost overwrite, notes/fine-velocity-bc.md)
        // BcStep calls the FDS routines with the FDS mesh number of the bound box (the registry Level carries offset 0)
        o.bc->ext_ghost = true;   // requirement of levels > 0 (BcStep::exchange): the driver saves the face velocities, MATCH_VELOCITY is not run
        o.bc->vel_bc = vel_bc;    // fine boxes have no wall cells: the domain-face wall type of level 0 (patch 0015)
        o.bc->face_wtype = face_wtype; o.bc->l0_len = l0.geom.Domain().length();   // mixed face: per level-0 face cell
        LevelCtx x;
        x.lev = &Lv; x.F = &F; x.bc = o.bc.get(); x.sd = &L.m_reg->side_data(lev); x.drho = o.drho.get(); x.dzz = o.dzz.get();
        x.nbox = nb; x.fds0 = nm0; x.off = total_boxes();
        lc.push_back(x);      // may move lc: c is reset below
        ++layout_gen;
        dt_new.resize(total_boxes(), L.m_dt);
        chg.resize(total_boxes(), 0);
        select(0);
#endif
    }
    void unbind_level(int lev)
    {
#ifdef FDS_HAVE_FINE_LEVEL
        if (lev < 1 || lev != static_cast<int>(lc.size()) - 1) die("unbind_level(" + std::to_string(lev) + "): only the top bound level (> 0) can be unbound");
        if (fds_fine_level_destroy(lev) != 0) die("unbind_level: fds_fine_level_destroy failed");
        lc.pop_back();
        ++layout_gen;
        own.resize(lev - 1);
        if (static_cast<int>(sfx.size()) > lev) sfx.resize(lev);   // the stage flux arrays of the level are re-registered with the next bind
        if (lc.size() == 1) { select(0); for (int k = 0; k < 2; ++k) { flux_set_override(k, {}); flux_set_mode(k, 0); } }   // no fine level left: level 0 must not keep the old interface override lists
        dt_new.resize(total_boxes()); chg.resize(total_boxes());
        select(0);
#else
        die("unbind_level: no fine-level support in this build");
#endif
    }
    bool local(int i) const { return c ? c->lev->dm[i] == amrex::ParallelDescriptor::MyProc() : is_local(l0, i); }

    // ------------------------------------------------------------ setup
    void check_scope()
    {
        std::vector<std::string> bad;
        auto chk = [&](bool cond, const char* what) { if (cond) bad.push_back(what); };
        chk(IP[7], "CC_IBM"); chk(IP[8], "HVAC"); chk(IP[11] != 0, "LEVEL_SET_MODE");
        chk(IP[12], "TUNNEL_PRECONDITIONER"); chk(IP[13] != 1, "pressure solver other than FFT"); chk(IP[14], "CYLINDRICAL"); chk(IP[15], "particle drag");
        // D-089: FREEZE_VELOCITY (IP[18]) allowed in AMR TimeLoop for cold-flow / passive-scalar adaptive demos; reacting-flow production must not depend on freeze.
        chk(IP[16], "SOLID_HEAT_TRANSFER_3D"); chk(IP[17] > 0, "Lagrangian particles"); chk(IP[19], "SOLID_PHASE_ONLY");
        chk(IP[22], "MPI particle exchange"); chk(IP[23], "obstruction mass exchange"); chk(IP[24] > 0, "agglomeration");
        // D-090: INIT_HRRPUV (IP[25]) allowed for manufactured-HRR / TAG_HRRPUV gates (R4-T2); reactions (IP[9]) allowed on level 0 (COMBUSTION in the corrector).
        chk(IP[26], "TERRAIN_CASE"); chk(IP[27], "READ_EXTERNAL"); chk(IP[28], "SYNTHETIC_EDDY_METHOD");
        if (!bad.empty()) {
            std::string s;
            for (auto& b : bad) s += " " + b;
            die("case outside the S5 scope (S6 and later):" + s);
        }
    }

    void setup()
    {
        fds_p_params(IP, RP);
        check_scope();
        nbox = static_cast<int>(l0.ba.size());
        ns = l0.dom.n_total;
        nt = l0.dom.n_tracked;
        double t_end, dtfill, dtmin, dt0; int nz;
        fds_k_consts(&t_end, &dtfill, &dtmin, &dt0, &nz);
        tc.t_end = t_end; tc.dt_end_fill = dtfill; tc.dt_end_minimum = dtmin;
        nzone = IP[5];
        L.m_rmin = RP[4]; L.m_rmax = RP[5];
        dt_new.assign(nbox, L.m_dt);
        chg.assign(nbox, 0);
        mip.assign(static_cast<std::size_t>(nbox) * 16, 0);
        for (int i = 0; i < nbox; ++i) {
            if (!local(i)) continue;
            double r[4];
            fds_p_mesh_info(i + 1, &mip[16 * static_cast<std::size_t>(i)], r);
            if (mip[16 * i + 6] != 0) die("IPS /= 0 (stretched or transposed Poisson direction) is outside the S5 scope");
        }
        L.m_F.reset(new Fields(l0, ns));
        L.m_sd.reset(new SideData(l0, fds_cell_walls));
        L.m_bc.reset(new BcStep(l0, *L.m_F));
        drho.reset(new amrex::MultiFab(l0.ba, l0.dm, 1, 1));
        dzz.reset(new amrex::MultiFab(l0.ba, l0.dm, std::max(1, ns), 0));
        drho->setVal(0.0); dzz->setVal(0.0);
        rhs.reset(new amrex::MultiFab(l0.ba, l0.dm, 1, 0));
        phi.reset(new amrex::MultiFab(l0.ba, l0.dm, 1, 1));
        rhs->setVal(0.0); phi->setVal(0.0);
        lc.assign(1, LevelCtx());
        lc[0].lev = &l0; lc[0].F = L.m_F.get(); lc[0].bc = L.m_bc.get(); lc[0].sd = L.m_sd.get(); lc[0].drho = drho.get(); lc[0].dzz = dzz.get();
        lc[0].nbox = nbox; lc[0].fds0 = l0.fds_mesh_offset; lc[0].off = 0;
        c = &lc[0];
        L.m_reg.reset(new LevelRegistry(l0.dom, ns));
        L.m_reg->adopt_level0(l0, *L.m_F, *L.m_sd);
        L.bind_fields(true);
        pressure_bc_map();
        zone_setup();
        fds_outputs = fds_hook_step_outputs() != 0 && std::getenv("FDSTL_NO_FDS_OUTPUTS") == nullptr;
        amrex::Print() << "FDS-AMReX: FDS output writers (devc/hrr/mass/steps/out) " << (fds_outputs ? "driven through fds_setup(mode=3)" : "not available (main.f90 patch 0006 absent or FDSTL_NO_FDS_OUTPUTS set)") << "\n";
        if (fds_outputs && std::getenv("FDSTL_MESH_DUMPS") != nullptr) {  // patch 0011 (per-mesh dumps: slices, boundary files, Plot3D, ...)
            fds_hook_set_mesh_dumps(1);
            amrex::Print() << "FDS-AMReX: per-mesh FDS dumps on (FDSTL_MESH_DUMPS, needs main.f90 patch 0011)\n";
        }
        if (const char* e = std::getenv("FDSTL_FLUXCHK")) fluxchk = std::atoi(e);
        if (const char* e = std::getenv("FDSTL_FLUXCHK_KINDS")) fluxchk_kinds = std::atoi(e);
        if (const char* e = std::getenv("FDSTL_WSEAM")) wseam = std::atoi(e);
        if (wseam) {
            long te = 0, ti = 0;
            each_local([&](int nm) { int a = 0, b = 0; fds_wseam_counts(nm, &a, &b); te += a; ti += b; });
            amrex::Print() << "FDS-AMReX: wall seam check FDSTL_WSEAM=" << wseam << " (W1 lists built, rank 0: " << te << " external and " << ti << " internal wall cells in its boxes; 2 = negative control)\n";
        }
        if (fluxchk) amrex::Print() << "FDS-AMReX: interface flux hook check FDSTL_FLUXCHK=" << fluxchk << " (1 empty override, 2 no-op override with the read-out values, 3 scaled values: negative control)\n";
        if (const char* e = std::getenv("FDSTL_IFACE")) iface_mode = std::atoi(e);
        L.m_bc->iface_hook = [this](bool on) { iface(on); };
        if (const char* e = std::getenv("FDSTL_EXTGHOST")) ext_ghost = std::atoi(e) != 0;
        fds_hook_set_flag(ext_ghost ? 1 : 0); L.m_bc->ext_ghost = ext_ghost;
        amrex::Print() << "FDS-AMReX: EXTERNAL_GHOSTS_FILLED=" << (ext_ghost ? 1 : 0) << ", interface walls as regular faces (FDSTL_IFACE)=" << iface_mode << "\n";
        if (iface_mode) {
            // FDS ends its own initialisation with MATCH_VELOCITY, VISCOSITY_BC, VELOCITY_BC (main.f90 ~445-455) on the multi-mesh state, which leaves the interface-edge
            // arrays (EDGE%OMEGA, EDGE%TAU) from the interpolation branch. Redo it here with the interfaces as regular faces, so that the first step starts from the
            // edge data that a single mesh has (no interface edge data).
            state(false, true);
            each_local([&](int nm) { fds_k_visc(nm, 0); });
            ghost_exchange(*L.m_F, 4);
            L.m_bc->after_exchange(4, L.m_t, L.m_dt);
            ghost_exchange(*L.m_F, 6);
            L.m_bc->after_exchange(6, L.m_t, L.m_dt);
            // and FDS's initial DIVERGENCE_PART_1 (main.f90 ~632, D of the first step) with the interfaces as regular faces
            fds_p_init_div();
            iface(true);
            each_local([&](int nm) { fds_k_div1(L.m_t, L.m_dt, nm); });
            iface(false);
        }
    }

    // FDS Poisson boundary codes of the domain faces -> pb::BC (S5: homogeneous data only; the data is checked at every solve)
    void pressure_bc_map()
    {
        int code[6] = {-1, -1, -1, -1, -1, -1}, nn[3] = {0, 0, 0};
        int whi[6] = {-1, -1, -1, -1, -1, -1}, wlo[6] = {-99, -99, -99, -99, -99, -99};   // max and -min of the wall types of the meshes on each domain face (patch 0015)
        const amrex::Box dom = l0.geom.Domain();
        for (int i = 0; i < nbox; ++i) {
            if (!local(i)) continue;
            const amrex::Box vb = l0.ba[i];
            int fc[6];
            fds_get_face_bc(i + 1, fc);
            for (int f = 0; f < 6; ++f)
                if (fc[f] >= 0 && (f % 2 ? vb.bigEnd(f / 2) == dom.bigEnd(f / 2) : vb.smallEnd(f / 2) == dom.smallEnd(f / 2))) { whi[f] = std::max(whi[f], fc[f]); wlo[f] = std::max(wlo[f], -fc[f]); }
            const int c3[3] = {mip[16 * i + 3], mip[16 * i + 4], mip[16 * i + 5]};
            for (int d = 0; d < 3; ++d) {
                nn[d] = dom.length(d);
                if (vb.smallEnd(d) == dom.smallEnd(d)) code[d] = c3[d];
                if (vb.bigEnd(d) == dom.bigEnd(d)) code[3 + d] = c3[d];
            }
        }
        if (std::getenv("FDSTL_PBV")) std::fprintf(stderr, "[rank %d] codes before reduce %d %d %d %d %d %d\n", amrex::ParallelDescriptor::MyProc(), code[0], code[1], code[2], code[3], code[4], code[5]);
        for (int f = 0; f < 6; ++f) all_max(code[f]);
        for (int f = 0; f < 6; ++f) { all_max(whi[f]); all_max(wlo[f]); vel_bc[f] = whi[f] == -wlo[f] ? whi[f] : 0; }   // one type over the whole face, else 0 (plain mirror)
        for (int f = 0; f < 6; ++f) {   // mixed face: the level-0 wall type per face cell (BcStep::face_wtype)
            face_wtype[f].clear();
            if (vel_bc[f] != 0 || whi[f] < 0) continue;
            const int d = f / 2, t1 = d == 0 ? 1 : 0, t2 = d == 2 ? 1 : 2, n1 = dom.length(t1);
            std::vector<int> m(static_cast<std::size_t>(n1) * dom.length(t2), -1);
            for (int i = 0; i < nbox; ++i) {
                const amrex::Box vb = l0.ba[i];
                if (!local(i) || (f % 2 ? vb.bigEnd(d) != dom.bigEnd(d) : vb.smallEnd(d) != dom.smallEnd(d))) continue;
                const int l1 = vb.length(t1), l2 = vb.length(t2);
                std::vector<int> t(static_cast<std::size_t>(l1) * l2);
                fds_get_face_wall_types(i + 1, f + 1, l1, l2, t.data());
                for (int b = 0; b < l2; ++b) for (int a = 0; a < l1; ++a) { int& x = m[(vb.smallEnd(t1) - dom.smallEnd(t1) + a) + n1 * (vb.smallEnd(t2) - dom.smallEnd(t2) + b)]; x = std::max(x, t[a + l1 * b]); }
            }
            amrex::ParallelAllReduce::Max(m.data(), static_cast<int>(m.size()), amrex::ParallelContext::CommunicatorSub());
            face_wtype[f] = m;
        }
        if (!L.m_o.quiet) amrex::Print() << "driver set-up: fine-level wall types -x..+z (1 INERT, 3 MIRROR, else plain mirror) " << vel_bc[0] << " " << vel_bc[1] << " " << vel_bc[2] << " " << vel_bc[3] << " " << vel_bc[4] << " " << vel_bc[5] << "\n";
        if (std::getenv("FDSTL_PBV")) std::fprintf(stderr, "[rank %d] codes %d %d %d %d %d %d\n", amrex::ParallelDescriptor::MyProc(), code[0], code[1], code[2], code[3], code[4], code[5]);
        DirBc dm[3];
        for (int d = 0; d < 3; ++d) {
            const int n = dom.length(d);
            solve_dir[d] = !(d == 1 && IP[20] && n == 1);   // TWO_D: no y solve
            dm[d] = map_pressure_bc_direction(d, code[d], code[3 + d], n, l0.dom.periodic[d] != 0, !solve_dir[d]);   // PressureBcMap.H (thin-direction rules, S9)
            if (!dm[d].error.empty()) die(dm[d].error);
        }
        if (!solve_dir[1]) ignored_direction_follows_open(dm[1], dm[0], dm[2]);
        for (int d = 0; d < 3; ++d) {
            if (!dm[d].note.empty() && !L.m_o.quiet) amrex::Print() << "driver set-up: " << dm[d].note << "\n";
            bc_type[d] = static_cast<int>(dm[d].lo);
            bc_type[3 + d] = static_cast<int>(dm[d].hi);
        }
    }

    // D-028 set-up sums: pressure-zone volumes (replace the FDS values) and vent areas (compared with FDS)
    void zone_setup()
    {
        if (nzone > 0) {
            std::vector<int> g; std::vector<double> v;
            for (int i = 0; i < nbox; ++i) {
                if (!local(i)) continue;
                const long n = static_cast<long>(mip[16 * i]) * mip[16 * i + 1] * mip[16 * i + 2];
                std::vector<int> cz(n); std::vector<double> vt(n);
                fds_p_zone_volume_terms(i + 1, cz.data(), vt.data());
                for (long c = 0; c < n; ++c) if (cz[c] > 0) { g.push_back(cz[c] - 1); v.push_back(vt[c]); }
            }
            const std::vector<double> vol = exact_group_sums(nzone, g, v);
            double rel = 0.0;
            for (int z = 0; z < nzone; ++z) {
                const double f = fds_p_zone_volume(z + 1);
                rel = std::max(rel, std::abs(vol[z] - f) / std::max(std::abs(vol[z]), 1e-300));
                fds_p_zone_volume_set(z + 1, vol[z]);
            }
            if (!L.m_o.quiet) amrex::Print() << "driver set-up: " << nzone << " pressure zone volume(s) by exact fixed-point sum; max relative difference to FDS accumulation " << rel << "\n";
        }
        const int nvt = IP[31];
        if (nvt > 0) {
            std::vector<int> g; std::vector<double> v;
            std::vector<double> fdsv(nvt, -1.0);
            for (int i = 0; i < nbox; ++i) {
                if (!local(i)) continue;
                const int nv = fds_p_vent_count(i + 1);
                for (int k = 1; k <= nv; ++k) {
                    int ip[9]; double rp[2];
                    fds_p_vent_get(i + 1, k, ip, rp);
                    if (ip[8]) continue;   // circular vent: FDS value kept
                    const int grp = ip[7] - 1, ior = std::abs(ip[6]);
                    if (grp < 0 || grp >= nvt) continue;
                    fdsv[grp] = std::max(fdsv[grp], rp[1]);
                    const double a = (ior == 1) ? l0.dx[1] * l0.dx[2] : (ior == 2) ? l0.dx[0] * l0.dx[2] : l0.dx[0] * l0.dx[1];
                    const long nc = (ior == 1) ? static_cast<long>(ip[3] - ip[2]) * (ip[5] - ip[4]) : (ior == 2) ? static_cast<long>(ip[1] - ip[0]) * (ip[5] - ip[4])
                                                                                                        : static_cast<long>(ip[1] - ip[0]) * (ip[3] - ip[2]);
                    for (long c = 0; c < nc; ++c) { g.push_back(grp); v.push_back(a); }
                }
            }
            const std::vector<double> ar = exact_group_sums(nvt, g, v);
            for (int q = 0; q < nvt; ++q) all_max(fdsv[q]);
            double rel = 0.0;
            for (int q = 0; q < nvt; ++q) if (fdsv[q] >= 0.0) rel = std::max(rel, std::abs(ar[q] - fdsv[q]) / std::max(std::abs(ar[q]), 1e-300));
            vent_rel = rel;
            if (!L.m_o.quiet) amrex::Print() << "driver set-up: " << nvt << " vent total area(s) by exact fixed-point sum; max relative difference to FDS TOTAL_FDS_AREA " << rel << "\n";
        }
    }
    double vent_rel = 0.0;

    // ------------------------------------------------------------ helpers

    // Interface walls as no-wall faces around the kernels that would otherwise add wall-face terms for them (env FDSTL_IFACE=1: interfaces inside the level, =2 also the
    // periodic images at the domain edge; default 1, FDSTL_IFACE=0 restores the S5 behaviour)
    int iface_mode = 1;
    bool ext_ghost = true;   // EXTERNAL_GHOSTS_FILLED (patches 0003/0004); FDSTL_EXTGHOST=0 returns to the OMESH-average route of FDS
    void iface(bool on, bool all = false)
    {
        if (!iface_mode || lev_of_c() > 0) return;   // a fine box has no FDS WALL cells and no box-interface walls (its faces are filled by the ghost hook)
        const amrex::Box dom = c->lev->geom.Domain();
        each_local([&](int nm) {
            const amrex::Box b = c->lev->ba[bi(nm)];
            int e[6] = {b.smallEnd(0) == dom.smallEnd(0), b.bigEnd(0) == dom.bigEnd(0), b.smallEnd(1) == dom.smallEnd(1), b.bigEnd(1) == dom.bigEnd(1), b.smallEnd(2) == dom.smallEnd(2), b.bigEnd(2) == dom.bigEnd(2)};
            fds_p_iface_walls(nm, on ? 1 : 0, e, (iface_mode == 2 || all) ? 1 : 0);
        });
    }

    void state(bool pred, bool first) { fds_k_state(pred ? 1 : 0, first ? 1 : 0, L.m_icyc, L.m_rmin, L.m_rmax); }
    template <class F> void each_local(F f) { for (int i = 0; i < c->nbox; ++i) if (local(i)) f(c->fds0 + i + 1); }

    // ------------------------------------------------------------ S10.2 interface flux hooks (notes/flux-hooks-design.md)
    // ADV (kind 0) = the stage product FX*UU read by the density update, DIF (kind 1) = RHO_D_DZDX/Y/Z as the divergence reads them; nodal MultiFabs per level, ncomp = N_TOTAL_SCALARS, no ghost layer.
    // FDS face I (high face of cell I) = AMReX face lo+I in the normal direction; tangentially FDS J = b-lo+1. Registered with the Fortran hook module per local box; the hooks do nothing
    // for a box whose mode is 0, so the unhooked step is the original code.
    struct StageFlux { std::unique_ptr<amrex::MultiFab> mf[2][3]; std::vector<int> mode = {0, 0}; };
    std::vector<StageFlux> sfx;        // by level
    // ---- S10.3 wall-state seam (ADR-001 W1/W2; notes/wall-seam-design.md). FDSTL_WSEAM=1: before the density update and before WALL_BC of every stage the four host-produced wall arrays
    // (UVW_SAVE, U_GHOST, V_GHOST, W_GHOST) are staged over the box's WLIST_EXT (the stand-in of the device upload) with a checksum assertion, and after the stage the host arrays
    // are asserted unchanged (the wall kernels only consume them). 2 = negative control: a host value is perturbed after the upload, the assertion must stop the run.
    int wseam = 0;
    void wseam_upload(const char* where)
    {
        if (!wseam) return;
        each_local([&](int nm) {
            if (fds_wseam_upload(nm) != 0) die(std::string("wall seam: the staged copy differs from the host arrays right after the upload (box ") + std::to_string(nm) + ", " + where + ")");
            if (wseam == 2) fds_wseam_poke(nm);
        });
    }
    void wseam_check(const char* where)
    {
        if (!wseam) return;
        each_local([&](int nm) {
            if (fds_wseam_check(nm) != 0)
                die(std::string("wall seam: UVW_SAVE/U_GHOST/V_GHOST/W_GHOST of box ") + std::to_string(nm) + " changed on the host after the upload (checksum assertion of ADR-001 W2, stage " + where + ")");
        });
    }
    int fluxchk = 0;                   // env FDSTL_FLUXCHK: 1 bracket every density / divergence-1 call with a read-out and an empty override pass; 2 override with the read-out values
                                       // of every face (a no-op: result must be bitwise the unhooked one); 3 as 2 with every value scaled by 1+1e-3 (negative control: result must differ)
    int fluxchk_kinds = 3;             // env FDSTL_FLUXCHK_KINDS: bit 0 ADV, bit 1 DIF (default both)
    long fluxchk_calls[2] = {0, 0}, fluxchk_faces[2] = {0, 0};
    double fluxchk_max[2] = {0, 0};
    int lev_of_c() const { return static_cast<int>(c - lc.data()); }
    StageFlux& stage_flux_arrays()
    {
        const int lev = lev_of_c();
        if (static_cast<int>(sfx.size()) <= lev) sfx.resize(lev + 1);
        StageFlux& S = sfx[lev];
        if (S.mf[0][0]) return S;
        fds_flux_reserve(c->fds0 + c->nbox);
        for (int k = 0; k < 2; ++k)
            for (int d = 0; d < 3; ++d) {
                const amrex::BoxArray nba = amrex::convert(c->lev->ba, amrex::IntVect::TheDimensionVector(d));
                S.mf[k][d].reset(new amrex::MultiFab(nba, c->lev->dm, std::max(1, ns), 0));
                S.mf[k][d]->setVal(0.0);
                for (amrex::MFIter mfi(*S.mf[k][d]); mfi.isValid(); ++mfi) {
                    const amrex::Box b = mfi.validbox();
                    const int lb[4] = {d == 0 ? 0 : 1, d == 1 ? 0 : 1, d == 2 ? 0 : 1, 1};
                    const int ext[4] = {b.length(0), b.length(1), b.length(2), std::max(1, ns)};
                    if (fds_flux_register(c->fds0 + mfi.index() + 1, k, d, lb, ext, (*S.mf[k][d])[mfi].dataPtr()) != 0) die("fds_flux_register failed");
                }
            }
        return S;
    }
    void flux_set_mode(int kind, int mode)
    {
        stage_flux_arrays().mode[kind] = mode;
        each_local([&](int nm) { if (fds_flux_set_mode(nm, kind, mode) != 0) die("fds_flux_set_mode failed"); });
    }
    // override lists (FluxOverride: AMReX face indices) -> the Fortran hook table; lists are in the order of the local boxes of the level (MFIter order). Mode becomes 2 (empty lists: nothing changes).
    void flux_set_override(int kind, const std::vector<std::vector<fdsrt::FluxOverride>>& per_box)
    {
        StageFlux& S = stage_flux_arrays();
        int k = 0;
        for (amrex::MFIter mfi(*S.mf[kind][0]); mfi.isValid(); ++mfi, ++k) {
            const int nm = c->fds0 + mfi.index() + 1;
            const amrex::Box vbx = c->lev->ba[mfi.index()]; const amrex::IntVect lo = vbx.smallEnd();
            std::vector<int> dir, idx; std::vector<double> val; int n = 0;
            if (k < static_cast<int>(per_box.size()))
                for (const fdsrt::FluxOverride& fo : per_box[k]) {
                    if (fo.nscal != ns) die("flux override: nscal differs from N_TOTAL_SCALARS");
                    if (fo.dir < 0 || fo.dir > 2) die("flux override: bad direction");
                    for (std::size_t f = 0; f < fo.face.size(); ++f) {
                        dir.push_back(fo.dir);
                        for (int e = 0; e < 3; ++e) idx.push_back(fo.face[f][e] - lo[e] + (e == fo.dir ? 0 : 1));
                        for (int q = 0; q < ns; ++q) val.push_back(fo.value[f * ns + q]);
                        ++n;
                    }
                }
            const int r = fds_flux_set_override(nm, kind, n, dir.data(), idx.data(), ns, val.data());
            if (r != 0) die(std::string("flux override rejected (") + (r == 2 ? "a face outside the valid range of the box" : "bad argument") + ") on box " + std::to_string(mfi.index()));
        }
        flux_set_mode(kind, 2);
    }
    // one list per local box and direction with every face of the read-out array (scaled): the whole-box no-op override
    std::vector<std::vector<fdsrt::FluxOverride>> flux_lists_from_readout(int kind, double scale)
    {
        StageFlux& S = stage_flux_arrays();
        std::vector<std::vector<fdsrt::FluxOverride>> per_box;
        for (amrex::MFIter mfi(*S.mf[kind][0]); mfi.isValid(); ++mfi) {
            per_box.emplace_back();
            for (int d = 0; d < 3; ++d) {
                fdsrt::FluxOverride fo; fo.dir = d; fo.nscal = ns;
                const amrex::Box b = (*S.mf[kind][d])[mfi].box();
                const amrex::Array4<const double> a = S.mf[kind][d]->const_array(mfi);
                for (int kk = b.smallEnd(2); kk <= b.bigEnd(2); ++kk)
                    for (int jj = b.smallEnd(1); jj <= b.bigEnd(1); ++jj)
                        for (int ii = b.smallEnd(0); ii <= b.bigEnd(0); ++ii) {
                            fo.face.push_back({ii, jj, kk});
                            for (int q = 0; q < ns; ++q) { fo.value.push_back(a(ii, jj, kk, q) * scale); fluxchk_max[kind] = std::max(fluxchk_max[kind], std::abs(a(ii, jj, kk, q))); }
                            ++fluxchk_faces[kind];
                        }
                per_box.back().push_back(std::move(fo));
            }
        }
        return per_box;
    }
    // ADV read-out of the stage (needs FX: call after the viscosity/MFD stage): the density update is skipped for the box by the hook, nothing else changes
    void flux_readout_adv(double t, double dt)
    {
        StageFlux& S = stage_flux_arrays(); const int old = S.mode[0];
        flux_set_mode(0, 1);
        each_local([&](int nm) { fds_p_dens_pre(t, dt, nm); });
        flux_set_mode(0, old);
    }
    std::vector<double> zsave;
    void zone_save() { if (nzone <= 0) return; zsave.assign(3 * nzone, 0.0); fds_p_zone_get(nzone, zsave.data(), zsave.data() + nzone, zsave.data() + 2 * nzone); }
    void zone_restore() { if (nzone <= 0) return; fds_p_zone_set(nzone, zsave.data(), zsave.data() + nzone, zsave.data() + 2 * nzone); }
    void div1_plain() { iface(true); each_local([&](int nm) { fds_k_div1(L.m_t, L.m_dt, nm); }); iface(false); }
    // DIF read-out: the stage divergence-1 itself with the hook copying the fluxes (the kernel result is the normal one). The zone sums are saved first so that the override pass can redo them.
    void flux_readout_dif()
    {
        StageFlux& S = stage_flux_arrays(); const int old = S.mode[1];
        flux_set_mode(1, 1); zone_save(); div1_plain(); flux_set_mode(1, old);
    }
    // DIF override pass: DIVERGENCE_PART_1 once more with the listed faces replaced (zone sums restored first; the kernel is re-runnable with identical inputs)
    void flux_apply_dif()
    {
        StageFlux& S = stage_flux_arrays();
        if (S.mode[1] != 2) return;
        zone_restore(); div1_plain();
    }
    // FDSTL_FLUXCHK: bracket the stage calls (see fluxchk above)
    void fluxchk_prepare(int kind)
    {
        ++fluxchk_calls[kind];
        if (fluxchk == 1) flux_set_override(kind, {});
        else flux_set_override(kind, flux_lists_from_readout(kind, fluxchk == 3 ? 1.001 : 1.0));
    }
    void fluxchk_density(bool pred, double t, double dt)
    {
        (void)pred;
        flux_readout_adv(t, dt);
        fluxchk_prepare(0);
    }
    void fluxchk_div1()
    {
        flux_readout_dif();
        fluxchk_prepare(1);
        flux_apply_dif();
        flux_set_override(1, {}); flux_set_mode(1, 0);
    }

    // ------------------------------------------------------------ D-031 level clip (port of the kernel check's clip_level; gather form)
    void clip_level(bool pred, std::vector<int>& flags)
    {
        Fields& F = *c->F;
        amrex::MultiFab& RP_ = F[pred ? "RHOS" : "RHO"];
        amrex::MultiFab& RZ = F[pred ? "ZZS" : "ZZ"];
        const amrex::Periodicity per = c->lev->geom.periodicity();
        RP_.FillBoundary(0, 1, RP_.nGrowVect(), per);
        RZ.FillBoundary(0, RZ.nComp(), RZ.nGrowVect(), per);
        const amrex::iMultiFab& mask = c->sd->mask();
        const int NS = RZ.nComp(), NT = nt;
        const double rmin = L.m_rmin, rmax = L.m_rmax;
        auto mkdx = [&](const amrex::Box& qb, std::vector<double> dxv[3]) { for (int d = 0; d < 3; ++d) dxv[d].assign(qb.length(d), c->lev->dx[d]); };
        flags.assign(c->nbox, 0);
        int fl[2] = {0, 0};
        std::vector<double> dxv[3];
        for (amrex::MFIter mfi(RP_); mfi.isValid(); ++mfi) {
            const amrex::Box qb = RP_[mfi].box(), mb = mask[mfi].box(), db = (*c->drho)[mfi].box(), v = mfi.validbox();
            mkdx(qb, dxv);
            int f[2], n;
            fds_clip_density(qb.loVect(), qb.hiVect(), mb.loVect(), mb.hiVect(), v.loVect(), v.hiVect(), db.loVect(), db.hiVect(), v.loVect(), v.hiVect(),
                             RP_[mfi].dataPtr(), mask[mfi].dataPtr(), dxv[0].data(), dxv[1].data(), dxv[2].data(), rmin, rmax, (*c->drho)[mfi].dataPtr(), f, &n);
            fl[0] |= f[0]; fl[1] |= f[1];
            flags[mfi.index()] = f[0] | (f[1] << 1);
        }
        bool rmn = fl[0], rmx = fl[1];
        amrex::ParallelAllReduce::Or(rmn, amrex::ParallelContext::CommunicatorSub());
        amrex::ParallelAllReduce::Or(rmx, amrex::ParallelContext::CommunicatorSub());
        if (rmn || rmx) {
            for (amrex::MFIter mfi(RP_); mfi.isValid(); ++mfi) {
                const amrex::Box qb = RP_[mfi].box(), db = (*c->drho)[mfi].box(), v = mfi.validbox();
                fds_clip_density_apply(qb.loVect(), qb.hiVect(), v.loVect(), v.hiVect(), db.loVect(), db.hiVect(), RP_[mfi].dataPtr(), (*c->drho)[mfi].dataPtr(), rmin, rmax);
            }
            RP_.FillBoundary(0, 1, RP_.nGrowVect(), per);
        }
        bool anyz = false;
        std::vector<int> zfl(NT, 0);
        if (NT == 1) {
            if (rmn || rmx)
                for (amrex::MFIter mfi(RP_); mfi.isValid(); ++mfi) {
                    const amrex::Box rb = RP_[mfi].box(), fb = RZ[mfi].box(), v = mfi.validbox();
                    fds_clip_species_one(rb.loVect(), rb.hiVect(), fb.loVect(), fb.hiVect(), v.loVect(), v.hiVect(), RP_[mfi].dataPtr(), RZ[mfi].dataPtr());
                }
        } else {
            std::vector<long> cclip(NT, 0);
            for (amrex::MFIter mfi(RP_); mfi.isValid(); ++mfi) {
                const amrex::Box fb = RZ[mfi].box(), rb = RP_[mfi].box(), mb = mask[mfi].box(), v = mfi.validbox();
                mkdx(fb, dxv);
                for (int n = 1; n <= NT; ++n) {
                    int f, ccnt;
                    fds_clip_species(fb.loVect(), fb.hiVect(), rb.loVect(), rb.hiVect(), mb.loVect(), mb.hiVect(), v.loVect(), v.hiVect(), NS, n, RP_[mfi].dataPtr(),
                                     RZ[mfi].dataPtr(), mask[mfi].dataPtr(), dxv[0].data(), dxv[1].data(), dxv[2].data(), (*c->dzz)[mfi].dataPtr(n - 1), &f, &ccnt);
                    zfl[n - 1] |= f;
                    cclip[n - 1] += ccnt;
                }
            }
            for (int n = 0; n < NT; ++n) { bool b = zfl[n]; amrex::ParallelAllReduce::Or(b, amrex::ParallelContext::CommunicatorSub()); zfl[n] = b; anyz = anyz || b; }
            // CLIP_MASS/CLIP_CELLS: D-028 sum of source-only Δ(ρZ)·V (not net DELTA); NCLIP cell counts. Density-only path (NT==1) leaves zeros.
            {
                const int lev = static_cast<int>(c - lc.data());
                if (lev >= 0 && nt > 0) {
                    ensure_clip_budget();
                    if (anyz) {
                        const double vol = c->lev->dx[0] * c->lev->dx[1] * c->lev->dx[2];
                        const amrex::iMultiFab* cov = L.m_reg ? L.m_reg->covered_mask(lev) : nullptr;
                        for (int n = 0; n < NT; ++n) {
                            long cc = cclip[n];
                            amrex::ParallelAllReduce::Sum(cc, amrex::ParallelContext::CommunicatorSub());
                            clip_cells_acc[clip_slot(lev, n)] += cc;
                        }
                        for (int n = 1; n <= NT; ++n) {
                            if (!zfl[n - 1]) continue;
                            std::vector<int> g;
                            std::vector<double> vsum;
                            for (amrex::MFIter mfi(RP_); mfi.isValid(); ++mfi) {
                                const auto rz = RZ.const_array(mfi);
                                const auto rp = RP_.const_array(mfi);
                                const auto dz = c->dzz->const_array(mfi);
                                const auto mk = mask.const_array(mfi);
                                amrex::Array4<const int> cv;
                                if (cov) cv = cov->const_array(mfi);
                                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                                    if (cov && cv(i, j, k) != 0) return;
                                    if (mk(i, j, k, 0) != 0) return;
                                    const double q = rz(i, j, k, n - 1), r = rp(i, j, k);
                                    if (q >= 0.0 && q <= r) return;
                                    g.push_back(0);
                                    vsum.push_back(dz(i, j, k, n - 1) * vol);
                                });
                            }
                            clip_mass_acc[clip_slot(lev, n - 1)] += exact_group_sums(1, g, vsum)[0];
                        }
                    }
                }
            }
            for (amrex::MFIter mfi(RP_); mfi.isValid(); ++mfi) {
                const amrex::Box fb = RZ[mfi].box(), rb = RP_[mfi].box(), v = mfi.validbox();
                for (int n = 1; n <= NT; ++n)
                    if (zfl[n - 1])
                        fds_clip_species_apply(rb.loVect(), rb.hiVect(), fb.loVect(), fb.hiVect(), v.loVect(), v.hiVect(), NS, n, RP_[mfi].dataPtr(), RZ[mfi].dataPtr(),
                                               (*c->dzz)[mfi].dataPtr(n - 1));
                if (rmn || rmx || anyz)
                    fds_clip_renorm(rb.loVect(), rb.hiVect(), fb.loVect(), fb.hiVect(), v.loVect(), v.hiVect(), NS, NT, RP_[mfi].dataPtr(), RZ[mfi].dataPtr(),
                                    mask[mfi].dataPtr());
            }
        }
    }

    // DENSITY of every box (the D-031 split form: pre-clip half, level gather clip, post-clip half)
    void density(bool pred, double t, double dt)
    {
        wseam_upload("DENSITY");
        const bool fx = fluxchk && (fluxchk_kinds & 1);
        if (fx) fluxchk_density(pred, t, dt);
        each_local([&](int nm) { fds_p_dens_pre(t, dt, nm); });
        const char* pcs = std::getenv("FDSTL_PRECLIP");     // SCRATCH PATCH: dump RHOS/ZZS before the level clip and DELTA_RHO after it
        const bool pcd = pcs && std::atoi(pcs) == L.m_icyc;
        auto dumpfab = [&](const char* tag, amrex::MultiFab& mf) {
            for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
                const amrex::FArrayBox& fab = mf[mfi]; const amrex::Box b = fab.box(); const amrex::Box v = mfi.validbox();
                std::FILE* fp = std::fopen((dir + "/raw_" + tag + "_b" + std::to_string(mfi.index()) + ".bin").c_str(), "wb"); if (!fp) continue;
                int hdr[16] = {b.smallEnd(0), b.smallEnd(1), b.smallEnd(2), b.bigEnd(0), b.bigEnd(1), b.bigEnd(2), fab.nComp(), 0, v.smallEnd(0), v.smallEnd(1), v.smallEnd(2), v.bigEnd(0), v.bigEnd(1), v.bigEnd(2), 0, 0};
                std::fwrite(hdr, sizeof(int), 16, fp); std::fwrite(fab.dataPtr(), sizeof(double), static_cast<std::size_t>(b.numPts()) * fab.nComp(), fp); std::fclose(fp);
            }
        };
        if (pcd) { dumpfab(pred ? "preclip_p_RHOS" : "preclip_c_RHO", (*c->F)[pred ? "RHOS" : "RHO"]); dumpfab(pred ? "preclip_p_ZZS" : "preclip_c_ZZ", (*c->F)[pred ? "ZZS" : "ZZ"]); }
        std::vector<int> flags;
        const std::vector<double> zb = clip_species_sums(pred);
        // DENSITY returns at once on the first step (ICYC<=1, mass.f90 guard): ZZ/ZZS still hold mass fractions, so the level clip (which treats them as RHO*ZZ and
        // renormalises them to RHO) must not run either; it turned Y_air=1 into RHO in an &INIT hot region (rho = rho^2 at step 2).
        if (fds_p_dens_skip() != 0) flags.assign(c->nbox, 0);
        else clip_level(pred, flags);
        {
            const std::vector<double> za = clip_species_sums(pred);
            const std::size_t need = static_cast<std::size_t>(nt) * lc.size(), off = static_cast<std::size_t>(nt) * static_cast<std::size_t>(c - lc.data());
            if (clip_species_acc.size() < static_cast<std::size_t>(nt)) clip_species_acc.resize(nt, 0.0);
            if (clip_species_pred.size() < need) clip_species_pred.resize(need, 0.0);
            for (int n = 0; n < nt; ++n) {
                if (pred) clip_species_pred[off + n] = 0.5 * (za[n] - zb[n]);
                else { clip_species_acc[n] += clip_species_pred[off + n] + (za[n] - zb[n]); clip_species_pred[off + n] = 0.0; }
            }
        }
        if (pcd) {
            double dmax = 0.0; long nz = 0;
            for (amrex::MFIter mfi(*c->drho); mfi.isValid(); ++mfi) { const amrex::FArrayBox& fab = (*c->drho)[mfi]; const double* p = fab.dataPtr(); const long np = fab.box().numPts(); for (long q = 0; q < np; ++q) { if (p[q] != 0.0) ++nz; dmax = std::max(dmax, std::abs(p[q])); } }
            std::fprintf(stderr, "PRECLIP icyc %d %s rmin %.17g rmax %.17g clip flags (rhomin|rhomax<<1) per box:", L.m_icyc, pred ? "pred" : "corr", L.m_rmin, L.m_rmax);
            for (int f : flags) std::fprintf(stderr, " %d", f);
            std::fprintf(stderr, "  DELTA_RHO nonzero %ld max|.| %.6g\n", nz, dmax);
            dumpfab(pred ? "preclip_p_DELTA_RHO" : "preclip_c_DELTA_RHO", *c->drho);
        }
        each_local([&](int nm) { fds_k_set_flags(nm, flags[bi(nm)]); fds_k_dens_post(t, dt, nm); });
        if (fx) { flux_set_override(0, {}); flux_set_mode(0, 0); }
        wseam_check("DENSITY");
    }

    // ------------------------------------------------------------ D-028: zone integrals (DSUM, PSUM, USUM) as exact fixed-point sums
    void zone_sums(bool pred)
    {
        if (nzone <= 0) return;
        const bool exact = L.m_o.exact_zone_sums || std::getenv("FDSTL_EXACT_ZONES") != nullptr;
        const bool diag = exact || std::getenv("FDSTL_ZONE_DIAG") != nullptr;
        std::vector<double> fds_loc(3 * nzone, 0.0);
        {
            std::vector<double> d(nzone), p(nzone), u(nzone);
            fds_p_zone_get(nzone, d.data(), p.data(), u.data());
            for (int z = 0; z < nzone; ++z) { fds_loc[z] = d[z]; fds_loc[nzone + z] = p[z]; fds_loc[2 * nzone + z] = u[z]; }
        }
        // FDS's own order: the kernels accumulated DSUM/PSUM/USUM per rank in box order; then the MPI_SUM reduction of INITIALIZE_DIVERGENCE_INTEGRALS' partner (main.f90 ~2111)
        amrex::ParallelAllReduce::Sum(fds_loc.data(), 3 * nzone, amrex::ParallelContext::CommunicatorSub());
        if (!diag) {
            if (std::getenv("FDSTL_ZONES")) for (int z = 0; z < nzone; ++z) amrex::Print() << "  [zone] icyc " << L.m_icyc << (pred ? " P " : " C ") << "FDS order DSUM " << std::setprecision(17) << fds_loc[z] << " PSUM " << fds_loc[nzone + z] << " USUM " << fds_loc[2 * nzone + z] << "\n";
            fds_p_zone_set(nzone, fds_loc.data(), fds_loc.data() + nzone, fds_loc.data() + 2 * nzone);
            return;
        }
        std::vector<int> g; std::vector<double> v;
        each_local([&](int nm) {
            const int i = nm - 1;
            const long nc = static_cast<long>(mip[16 * i]) * mip[16 * i + 1] * mip[16 * i + 2];
            const int nw = mip[16 * i + 7];
            std::vector<int> cz(nc), wz(std::max(nw, 1)); std::vector<double> dt_(nc), pt(nc), wt(std::max(nw, 1));
            fds_p_zone_terms(nm, cz.data(), dt_.data(), pt.data(), wz.data(), wt.data());
            for (long c = 0; c < nc; ++c) if (cz[c] > 0) {
                g.push_back(cz[c] - 1); v.push_back(dt_[c]);
                g.push_back(nzone + cz[c] - 1); v.push_back(pt[c]);
            }
            for (int w = 0; w < nw; ++w) if (wz[w] > 0) { g.push_back(2 * nzone + wz[w] - 1); v.push_back(wt[w]); }
        });
        const std::vector<double> s = exact_group_sums(3 * nzone, g, v);
        for (int q = 0; q < 3 * nzone; ++q) zone_rel_step = std::max(zone_rel_step, std::abs(s[q] - fds_loc[q]) / std::max(std::abs(s[q]), 1e-30));
        if (std::getenv("FDSTL_ZONES")) for (int z = 0; z < nzone; ++z) amrex::Print() << "  [zone] icyc " << L.m_icyc << (pred ? " P " : " C ") << "DSUM " << std::setprecision(17) << s[z] << " PSUM " << s[nzone + z] << " USUM " << s[2 * nzone + z] << " (FDS order: " << fds_loc[z] << " " << fds_loc[nzone + z] << " " << fds_loc[2 * nzone + z] << ")\n";
        if (exact) fds_p_zone_set(nzone, s.data(), s.data() + nzone, s.data() + 2 * nzone);
        else fds_p_zone_set(nzone, fds_loc.data(), fds_loc.data() + nzone, fds_loc.data() + 2 * nzone);
    }

    // ------------------------------------------------------------ pressure (PRESSURE_ITERATION_SCHEME, main.f90:1649)
    void physical_flags(int i, int fl[6]) const
    {
        const amrex::Box dom = l0.geom.Domain(), vb = l0.ba[i];
        for (int d = 0; d < 3; ++d) {
            fl[2 * d] = (!l0.dom.periodic[d] && vb.smallEnd(d) == dom.smallEnd(d)) ? 1 : 0;       // FDS order: xlo, xhi, ylo, yhi, zlo, zhi
            fl[2 * d + 1] = (!l0.dom.periodic[d] && vb.bigEnd(d) == dom.bigEnd(d)) ? 1 : 0;
        }
    }

    void solve_poisson(bool pred)
    {
        const double ts0 = amrex::second();
        solve_poisson_impl(pred);
        t_solve += amrex::second() - ts0;
    }
    void solve_poisson_impl(bool pred)
    {
        Fields& F = *L.m_F;
        level0_rhs_and_check();
        solve_level0(pred, F);
    }

    // FDS face order 0..5 (xlo,xhi,ylo,yhi,zlo,zhi) -> pb::face_index(dir,side) = side*3+dir
    static int fds_face_to_pb(int f) { return (f % 2) * 3 + (f / 2); }

    // Pack FDS BXS..BZF of level-0 meshes into 6 cell-centred MultiFabs on lv's BA (domain-face layers).
    // Periodic pb faces stay zero (fold rejects data on Periodic).
    void pack_l0_boundary_mfs(std::array<std::unique_ptr<amrex::MultiFab>, 6>& out)
    {
        for (int f = 0; f < 6; ++f) {
            out[f].reset(new amrex::MultiFab(l0.ba, l0.dm, 1, 0));
            out[f]->setVal(0.0);
        }
        for (int bi = 0; bi < nbox; ++bi) {
            if (!local(bi)) continue;
            const int nm = bi + 1;
            int fl[6];
            physical_flags(bi, fl);
            const int ibar = mip[16 * bi + 0], jbar = mip[16 * bi + 1], kbar = mip[16 * bi + 2];
            const amrex::Box vb = l0.ba[bi];
            const int i0 = vb.smallEnd(0), j0 = vb.smallEnd(1), k0 = vb.smallEnd(2);
            for (int ff = 0; ff < 6; ++ff) {
                if (!fl[ff]) continue;
                const int pb_f = fds_face_to_pb(ff);
                if (static_cast<pb::BC>(bc_type[pb_f]) == pb::BC::Periodic) continue;
                const int d = ff / 2;
                long nexpect = 0;
                if (d == 0) nexpect = (long)jbar * (long)kbar;
                else if (d == 1) nexpect = (long)ibar * (long)kbar;
                else nexpect = (long)ibar * (long)jbar;
                std::vector<double> buf((std::size_t)nexpect);
                int ierr = 0;
                fds_p_get_bcface(nm, ff, buf.data(), nexpect, &ierr);
                if (ierr != 0) die(std::string("fds_p_get_bcface face=") + std::to_string(ff) + " ierr=" + std::to_string(ierr));
                // write into the fab that owns this box
                for (amrex::MFIter mfi(*out[pb_f]); mfi.isValid(); ++mfi) {
                    if (mfi.index() != bi) continue;
                    auto a = out[pb_f]->array(mfi);
                    int n = 0;
                    if (d == 0) {
                        const int i = (ff % 2 == 0) ? vb.smallEnd(0) : vb.bigEnd(0);
                        for (int kk = 0; kk < kbar; ++kk)
                            for (int jj = 0; jj < jbar; ++jj) a(i, j0 + jj, k0 + kk) = buf[n++];
                    } else if (d == 1) {
                        const int j = (ff % 2 == 0) ? vb.smallEnd(1) : vb.bigEnd(1);
                        for (int kk = 0; kk < kbar; ++kk)
                            for (int ii = 0; ii < ibar; ++ii) a(i0 + ii, j, k0 + kk) = buf[n++];
                    } else {
                        const int k = (ff % 2 == 0) ? vb.smallEnd(2) : vb.bigEnd(2);
                        for (int jj = 0; jj < jbar; ++jj)
                            for (int ii = 0; ii < ibar; ++ii) a(i0 + ii, j0 + jj, k) = buf[n++];
                    }
                }
            }
        }
    }

    // Inject L0 face MultiFab onto a fine-level MultiFab (piecewise-constant via coarsened ParallelCopy).
    void inject_bd_to_level(const amrex::MultiFab& l0mf, amrex::MultiFab& fine, const amrex::IntVect& ratio)
    {
        fine.setVal(0.0);
        if (ratio == amrex::IntVect(1)) {
            fine.ParallelCopy(l0mf);
            return;
        }
        amrex::BoxArray cba = fine.boxArray();
        cba.coarsen(ratio);
        amrex::MultiFab ctmp(cba, fine.DistributionMap(), 1, 0);
        ctmp.setVal(0.0);
        ctmp.ParallelCopy(l0mf);
        for (amrex::MFIter mfi(fine); mfi.isValid(); ++mfi) {
            auto fa = fine.array(mfi);
            auto ca = ctmp.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                fa(i, j, k) = ca(amrex::coarsen(i, ratio[0]), amrex::coarsen(j, ratio[1]), amrex::coarsen(k, ratio[2]));
            });
        }
    }

    // Build BoundaryData (+ owned MultiFabs) for one level from L0-packed face MFs; skip Periodic faces.
    pb::BoundaryData make_bd_for_level(const Level& lv, const std::array<std::unique_ptr<amrex::MultiFab>, 6>& l0bd,
                                       std::array<std::unique_ptr<amrex::MultiFab>, 6>& owned)
    {
        pb::BoundaryData bd;
        // For levels > 1, ratio_from_parent is only relative to parent; inject stepwise from L0 index space:
        // use cumulative ref from geom Domain vs l0 Domain length.
        amrex::IntVect rL0(1);
        if (&lv != &l0) {
            const amrex::Box d0 = l0.geom.Domain(), df = lv.geom.Domain();
            for (int d = 0; d < 3; ++d) {
                const int n0 = d0.length(d), nf = df.length(d);
                rL0[d] = (n0 > 0) ? (nf / n0) : 1;
            }
        }
        for (int f = 0; f < 6; ++f) {
            if (static_cast<pb::BC>(bc_type[f]) == pb::BC::Periodic) {
                bd.value[f] = nullptr;
                bd.constant[f] = 0.0;
                continue;
            }
            owned[f].reset(new amrex::MultiFab(lv.ba, lv.dm, 1, 0));
            inject_bd_to_level(*l0bd[f], *owned[f], rL0);
            bd.value[f] = owned[f].get();
            bd.constant[f] = 0.0;
        }
        return bd;
    }

    // Mixed Dirichlet face: Neumann over the SOLID level-0 face cells in the composite solve (pb face_neumann; U and FV are 0 there, BcStep::zero_mixed_wall).
    void set_face_neumann(pb::PressureProblem& p) const
    {
        for (int d = 0; d < 3; ++d)
            for (int side = 0; side < 2; ++side) {
                const std::vector<int>& w = face_wtype[2 * d + side];
                const int pf = pb::face_index(d, side);
                if (w.empty() || p.bc[pf] != pb::BC::Dirichlet || std::find(w.begin(), w.end(), 1) == w.end()) continue;
                p.face_neumann[pf].resize(w.size());
                for (std::size_t c = 0; c < w.size(); ++c) p.face_neumann[pf][c] = w[c] == 1 ? 1 : 0;   // SOLID_BOUNDARY
            }
    }

    // Dirichlet (OPEN) domain faces of a fine level: H ghost = 2 g - H(inside), g = the level-0 BXS..BZF that fold_composite folded into the fine RHS.
    // The even mirror above is a Neumann ghost: the fine velocity correction then disagreed with the Dirichlet solve (div error along the face).
    void dirichlet_fine_h_ghosts(const Level& lv, amrex::MultiFab& H)
    {
        std::array<std::unique_ptr<amrex::MultiFab>, 6> l0bd, own;
        pack_l0_boundary_mfs(l0bd);
        const pb::BoundaryData bd = make_bd_for_level(lv, l0bd, own);
        const amrex::Box dom = lv.geom.Domain();
        for (int f = 0; f < 6; ++f) {
            if (static_cast<pb::BC>(bc_type[f]) != pb::BC::Dirichlet) continue;
            const int d = f % 3, side = f / 3, e = side ? dom.bigEnd(d) : dom.smallEnd(d);
            amrex::IntVect s(0); s[d] = side ? 1 : -1;
            for (amrex::MFIter mfi(H); mfi.isValid(); ++mfi) {
                amrex::Box b = mfi.validbox();
                if ((side ? b.bigEnd(d) : b.smallEnd(d)) != e) continue;
                b.setRange(d, e);
                auto h = H.array(mfi); const auto g = bd.value[f]->const_array(mfi);
                const std::vector<int>& w = face_wtype[2 * d + side];   // mixed face: even mirror (Neumann) over the SOLID level-0 cells
                const amrex::Box d0 = l0.geom.Domain();
                const int t1 = d == 0 ? 1 : 0, t2 = d == 2 ? 1 : 2, r1 = dom.length(t1) / d0.length(t1), r2 = dom.length(t2) / d0.length(t2);
                amrex::LoopOnCpu(b, [&](int i, int j, int k) {
                    const int c[3] = {i, j, k};
                    const bool wall = !w.empty() && w[amrex::coarsen(c[t1], r1) - d0.smallEnd(t1) + d0.length(t1) * (amrex::coarsen(c[t2], r2) - d0.smallEnd(t2))] == 1;
                    h(i + s[0], j + s[1], k + s[2]) = wall ? h(i, j, k) : 2.0 * g(i, j, k) - h(i, j, k);
                });
            }
        }
    }

    void fold_then_solve_level0(pb::PressureProblem& p, pb::PressureOptions& /*o*/)
    {
        std::array<std::unique_ptr<amrex::MultiFab>, 6> l0bd;
        pack_l0_boundary_mfs(l0bd);
        std::array<std::unique_ptr<amrex::MultiFab>, 6> owned;
        pb::BoundaryData bd = make_bd_for_level(l0, l0bd, owned);
        std::string msg;
        const pb::Status fs = pb::fold_boundary_data(p, bd, *rhs, &msg);
        if (fs != pb::Status::Ok) die(std::string("fold_boundary_data(single): ") + pb::to_string(fs) + ": " + msg);
    }

    void fold_composite(pb::PressureProblem& p, std::vector<amrex::MultiFab*>& rhs_levels)
    {
        std::array<std::unique_ptr<amrex::MultiFab>, 6> l0bd;
        pack_l0_boundary_mfs(l0bd);
        const int nl = static_cast<int>(p.levels.size());
        std::vector<pb::BoundaryData> bds(nl);
        // keep owned MFs alive for the fold call
        std::vector<std::array<std::unique_ptr<amrex::MultiFab>, 6>> owned(nl);
        for (int l = 0; l < nl; ++l) bds[l] = make_bd_for_level(*lc[l].lev, l0bd, owned[l]);
        std::string msg;
        const pb::Status fs = pb::fold_boundary_data(p, bds, rhs_levels, &msg);
        if (fs != pb::Status::Ok) die(std::string("fold_boundary_data(composite): ") + pb::to_string(fs) + ": " + msg);
        (void)owned;
    }

    // PRHS of the level-0 boxes (as FDS computed it, fds_p_rhs) into `rhs`, and the check that the Poisson boundary data are homogeneous (every solve)
    void level0_rhs_and_check()
    {
        for (amrex::MFIter mfi(*rhs); mfi.isValid(); ++mfi) {
            const int nm = mfi.index() + 1;
            const amrex::Box vb = mfi.validbox();
            std::vector<double> buf(static_cast<std::size_t>(vb.numPts()));
            fds_p_get_prhs(nm, buf.data());
            auto a = rhs->array(mfi);
            const int nx = vb.length(0), ny = vb.length(1);
            for (int k = 0; k < vb.length(2); ++k)
                for (int j = 0; j < ny; ++j)
                    for (int i = 0; i < nx; ++i) a(vb.smallEnd(0) + i, vb.smallEnd(1) + j, vb.smallEnd(2) + k) = buf[i + nx * (j + static_cast<std::size_t>(ny) * k)];
        }
        // S6: inhomogeneous BXS..BZF are folded into rhs before solve_pressure (fold_boundary_data). No die.
        if (std::getenv("FDSTL_PBV")) {
            double bm = 0.0;
            for (int i = 0; i < nbox; ++i) {
                if (!local(i)) continue;
                int fl[6], sf[6];
                physical_flags(i, fl);
                for (int f = 0; f < 6; ++f) sf[f] = (fl[f] && solve_dir[f / 2] && bc_type[(f % 2) * 3 + f / 2] != static_cast<int>(pb::BC::Periodic)) ? 1 : 0;
                double v[6];
                fds_p_bmax(i + 1, sf, v);
                for (double x : v) bm = std::max(bm, x);
            }
            all_max(bm);
            if (bm > 0.0) amrex::Print() << "  [pb] inhomogeneous BXS..BZF max|·|=" << bm << " (fold before solve)\n";
        }
    }
    // P2 solver study (RunOptions::pressure_*): backend, Krylov method and tolerance of every solve; one CSV row per solve (rank 0).
    // pressure/13 rules: shared true relative L2 residual check at the tolerance on every solve, setup frozen between regrids (workspaces), zero initial
    // guess unless pressure.warm_start=1 (then the previous solution: `warm` keeps phi between solves).
    std::FILE* pcsv = nullptr;
    bool study = false, warm = false;
    int nsetup = 0;               // operator set-ups seen so far: a later one (regrid, changed workspace) is a rebuild
    pb::PressureWorkspace pws0;   // level-0 solves of the study
    void study_options(pb::PressureOptions& o)
    {
        const RunOptions& ro = L.m_o;
        study = !(ro.pressure_backend.empty() && ro.krylov.empty() && ro.pressure_csv.empty());
        if (!study) return;
        const std::string& b = ro.pressure_backend;
        const std::string& k = ro.krylov;
        if (b == "mlmg") o.backend = pb::BackendKind::MLMG;
        else if (b == "hypre") o.backend = pb::BackendKind::HYPRE;
        else if (!b.empty() && b != "auto") die("--pressure-backend " + b + ": use auto, mlmg or hypre");
        using K = pb::HypreKrylov;
        o.mlmg_krylov = o.hypre.krylov = k == "none" ? K::None : k == "gmres" ? K::GMRES : k == "bicgstab" ? K::BiCGSTAB : k == "cg" ? K::PCG : K::Auto;
        if (o.mlmg_krylov == K::Auto && !k.empty()) die("--krylov " + k + ": use none, gmres, bicgstab or cg");
        o.tol_rel = ro.pressure_tol;
        pb::query_options(o);
        o.check_residual = true; o.trigger = pb::SolveDebug; o.residual_tol = o.tol_rel;   // full path: true residual on every solve
        warm = o.use_initial_guess;
        if (pcsv || !amrex::ParallelDescriptor::IOProcessor()) return;
        pcsv = std::fopen((ro.pressure_csv.empty() ? dir + "/" + ro.chid + "_pressure.csv" : ro.pressure_csv).c_str(), "w");
        if (!pcsv) die("cannot open the pressure CSV");
        std::setvbuf(pcsv, nullptr, _IOLBF, 0);
        std::fprintf(pcsv, "step,stage,level,backend,krylov_used,guess,first,rebuilt,iterations,own_residual,true_rel2,true_ok,setup_s,solve_s,rss_max_mb,rss_sum_mb,div_max\n");
    }
    void study_row(bool pred, const char* level, const pb::PressureResult& r)
    {
        if (!study) return;
        rusage ru; getrusage(RUSAGE_SELF, &ru);
        long rmax = ru.ru_maxrss, rsum = ru.ru_maxrss;   // KiB, peak resident set of this rank
        amrex::ParallelDescriptor::ReduceLongMax(rmax); amrex::ParallelDescriptor::ReduceLongSum(rsum);
        const pb::BackendStatus& s = r.backend_status;
        const bool first = !s.plan_reused, rebuilt = first && nsetup++ > 0;   // first: set up in this solve; rebuilt: a set-up after an earlier one
        if (pcsv) std::fprintf(pcsv, "%d,%s,%s,%s,%s,%s,%d,%d,%d,%.3e,%.3e,%d,%.3e,%.3e,%.1f,%.1f,\n", L.m_icyc, pred ? "pred" : "corr", level, r.backend.c_str(),
                               s.method.empty() ? "direct" : s.method.c_str(), warm ? "warm" : "zero", first ? 1 : 0, (rebuilt || r.workspace_rebuilt) ? 1 : 0, s.iterations,
                               s.own_residual, r.residual_rel2, r.residual_ok ? 1 : 0, s.setup_seconds, s.solve_seconds, rmax / 1024.0, rsum / 1024.0);
    }
    void solve_level0(bool pred, Fields& F)
    {
        if (std::getenv("FDSTL_PBV")) std::fprintf(stderr, "[rank %d] bc_type %d %d %d %d %d %d\n", amrex::ParallelDescriptor::MyProc(), bc_type[0], bc_type[1], bc_type[2], bc_type[3], bc_type[4], bc_type[5]);
        pb::PressureProblem p;
        p.ba = l0.ba; p.dm = l0.dm; p.geom = l0.geom;
        for (int f = 0; f < 6; ++f) p.bc[f] = static_cast<pb::BC>(bc_type[f]);
        p.rhs = rhs.get(); p.phi = phi.get();
        pb::PressureOptions o;
        o.backend = pb::BackendKind::FFT;
        o.verbose = std::getenv("FDSTL_PBV") ? 1 : 0;
        o.check_residual = std::getenv("FDSTL_PBV") != nullptr;
        if (!warm) phi->setVal(0.0);
        fold_then_solve_level0(p, o);
        study_options(o);
        const pb::PressureResult r = pb::solve_pressure(p, o, study ? &pws0 : nullptr);
        if (r.status != pb::Status::Ok && r.status != pb::Status::NotConverged) die(std::string("pressure interface: ") + pb::to_string(r.status) + ": " + r.message);
        study_row(pred, "0", r);
        if (std::getenv("FDSTL_PBV")) { amrex::Print() << "  [pb] backend=" << r.backend << " res_rel2=" << r.residual_rel2 << " relmax=" << r.residual_relmax << " residual_check=" << r.residual_check << " residual_limit=" << r.residual_limit << " residual_floor=" << r.residual_floor << " ok=" << r.residual_ok << " rhs max=" << rhs->norminf(0) << " phi max=" << phi->norminf(0); for (auto& c : r.components) amrex::Print() << " removed_mean=" << c.removed_mean << " rel=" << c.removed_rel; for (auto& w : r.warnings) amrex::Print() << " W:" << w; amrex::Print() << "\n"; }
        if (const char* e = std::getenv("FDSTL_PDUMP")) {
            // hand-over for the one-solve check against another Poisson solver (Role 2): the right-hand side PRHS (FDS IPS=0 layout, valid cells, I fastest) and the level solution
            // that becomes H/HS, both in the cell layout of the level (global, written by rank 0): <outdir>/pdump_<icyc>_<P|C>_{rhs,phi}.bin (doubles), dx dy dz in pdump_dx.txt
            if (std::atoi(e) == L.m_icyc) {
                const FieldSpec& cs = F.spec("D");
                const std::vector<double> a = gather_mf(*rhs, cs, 0), b = gather_mf(*phi, cs, 0);
                if (amrex::ParallelDescriptor::MyProc() == 0) {
                    const std::string base = dir + "/pdump_" + std::to_string(L.m_icyc) + (pred ? "_P_" : "_C_");
                    if (std::FILE* fp = std::fopen((base + "rhs.bin").c_str(), "wb")) { std::fwrite(a.data(), sizeof(double), a.size(), fp); std::fclose(fp); }
                    if (std::FILE* fp = std::fopen((base + "phi.bin").c_str(), "wb")) { std::fwrite(b.data(), sizeof(double), b.size(), fp); std::fclose(fp); }
                    if (std::FILE* fp = std::fopen((dir + "/pdump_dx.txt").c_str(), "w")) { std::fprintf(fp, "%.17g %.17g %.17g\n", l0.dx[0], l0.dx[1], l0.dx[2]); std::fclose(fp); }
                }
            }
        }
        const char* hn = pred ? "H" : "HS";
        amrex::MultiFab& H = F[hn];
        for (amrex::MFIter mfi(H); mfi.isValid(); ++mfi) {
            const amrex::Box vb = mfi.validbox();
            auto h = H.array(mfi);
            auto q = phi->const_array(mfi);
            amrex::LoopOnCpu(vb, [&](int i, int j, int k) { h(i, j, k) = q(i, j, k); });
        }
        H.FillBoundary(0, 1, H.nGrowVect(), l0.geom.periodicity());   // box-to-box and periodic ghosts
        each_local([&](int nm) { int fl[6]; physical_flags(nm - 1, fl); fds_p_h_ghost(nm, pred ? 1 : 0, fl); });   // physical faces: FDS's own H boundary values
        if (std::getenv("FDSTL_PBV")) {
            double mx = 0.0, mr = 0.0;
            for (amrex::MFIter mfi(H); mfi.isValid(); ++mfi) {
                auto h = H.const_array(mfi); auto b = rhs->const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) {
                    const double L = (h(i + 1, j, k) - 2 * h(i, j, k) + h(i - 1, j, k)) / (l0.dx[0] * l0.dx[0]) + (h(i, j + 1, k) - 2 * h(i, j, k) + h(i, j - 1, k)) / (l0.dx[1] * l0.dx[1]) +
                                     (h(i, j, k + 1) - 2 * h(i, j, k) + h(i, j, k - 1)) / (l0.dx[2] * l0.dx[2]);
                    mx = std::max(mx, std::abs(L - b(i, j, k))); mr = std::max(mr, std::abs(b(i, j, k)));
                });
            }
            double mw = 0.0, dphi = 0.0;
            for (amrex::MFIter mfi(H); mfi.isValid(); ++mfi) {
                auto h = H.const_array(mfi); auto b = rhs->const_array(mfi); auto q = phi->const_array(mfi);
                const amrex::Box vb = mfi.validbox();
                const int nx = vb.length(0), nz = vb.length(2);
                amrex::LoopOnCpu(vb, [&](int i, int j, int k) {
                    const int ip = (i + 1 - vb.smallEnd(0)) % nx + vb.smallEnd(0), im = (i - 1 - vb.smallEnd(0) + nx) % nx + vb.smallEnd(0);
                    const int kp = (k + 1 - vb.smallEnd(2)) % nz + vb.smallEnd(2), km = (k - 1 - vb.smallEnd(2) + nz) % nz + vb.smallEnd(2);
                    const double L = (h(ip, j, k) - 2 * h(i, j, k) + h(im, j, k)) / (l0.dx[0] * l0.dx[0]) + (h(i, j, kp) - 2 * h(i, j, k) + h(i, j, km)) / (l0.dx[2] * l0.dx[2]);
                    mw = std::max(mw, std::abs(L - b(i, j, k)));
                    dphi = std::max(dphi, std::abs(h(i, j, k) - q(i, j, k)));
                });
            }
            for (amrex::MFIter mfi(H); mfi.isValid(); ++mfi) {
                auto h = H.const_array(mfi); const amrex::Box vb = mfi.validbox();
                const int j = vb.smallEnd(1), k = 5;
                amrex::Print() << "  [pb] ghost probe: fab box " << H[mfi].box() << " h(lo-1)=" << h(vb.smallEnd(0) - 1, j, k) << " h(hi)=" << h(vb.bigEnd(0), j, k) << " h(hi+1)=" << h(vb.bigEnd(0) + 1, j, k)
                               << " h(lo)=" << h(vb.smallEnd(0), j, k) << " y: " << h(3, j - 1, k) << " " << h(3, j, k) << " " << h(3, j + 1, k) << " z: " << h(3, j, vb.smallEnd(2) - 1) << " " << h(3, j, vb.bigEnd(2)) << "\n";
            }
            amrex::Print() << "  [pb] residual from the H fab with ghosts: " << mx << " wrapped: " << mw << " |H-phi| " << dphi << " (rhs max " << mr << ")\n";
        }
    }

    // ------------------------------------------------------------ S14 composite (several bound levels)
    // One time step over a hierarchy: the stage pieces run per level (for_levels), the pressure is ONE composite Poisson solve over all levels (pressure_backend, Role 2), the
    // velocity of a fine face overwrites the coarse faces under it, the fine ghost faces come from the coarse level, and (comp_overwrite) the interface flux overwrite of the
    // scalar fluxes (D-050) runs between the read-out and the cell update. Same dt on every level, no subcycling (D-050). Without a second bound level none of this runs.
    bool multi() const { return lc.size() > 1; }
    bool comp_overwrite = true;                         // interface flux overwrite (FDSTL_OVERWRITE=0: negative control)
    int layout_gen = 0, comp_gen = -1;
    std::unique_ptr<pb::PressureWorkspace> pws;
    std::unique_ptr<fdsrt::FluxStageRunner> runner;
    std::vector<std::unique_ptr<amrex::MultiFab>> crhs, cphi;   // by level (index 0 unused: level 0 has rhs / phi)
    std::vector<std::unique_ptr<amrex::iMultiFab>> ccls;        // OBST cell_class by level (ObstMask.H); empty without OBST
    TimeLoop::CompositeReport crep;
    void ensure_comp()
    {
        if (comp_gen == layout_gen && !crhs.empty()) return;
        const int nl = static_cast<int>(lc.size());
        crhs.clear(); cphi.clear();
        crhs.resize(nl); cphi.resize(nl);
        std::vector<fdsrt::StageLevel> sl;
        for (int l = 0; l < nl; ++l) {
            const Level& lv = *lc[l].lev;
            if (l > 0) {
                crhs[l].reset(new amrex::MultiFab(lv.ba, lv.dm, 1, 0));
                cphi[l].reset(new amrex::MultiFab(lv.ba, lv.dm, 1, 0));
                crhs[l]->setVal(0.0); cphi[l]->setVal(0.0);
            }
            fdsrt::StageLevel x; x.geom = lv.geom; x.ba = lv.ba; x.dm = lv.dm; x.ratio_from_parent = l == 0 ? amrex::IntVect(1) : lv.ref_ratio_from_parent;
            sl.push_back(x);
        }
        if (const char* e = std::getenv("FDSTL_OVERWRITE")) comp_overwrite = std::atoi(e) != 0;
        runner.reset(new fdsrt::FluxStageRunner(sl, comp_overwrite));
        pws.reset(new pb::PressureWorkspace());
        load_obst_list();
        ccls.clear();
        // Level 0 is an FDS mesh: FDS treats its OBSTs itself (solve through solids, wall normal velocity, e.g. a MASS_FLUX burner, via NO_FLUX FV),
        // as with --levels 1. Only levels >= 1, where FDS has no OBST, are masked; a Solid L0 cell is under L1 (--tag-obst) or FDS's own.
        if (!pb::obst_list().empty()) for (int l = 0; l < nl; ++l) ccls.push_back(l == 0 ? nullptr : pb::build_obst_class(lc[l].lev->geom, l, lc[l].lev->ba, lc[l].lev->dm));
        comp_gen = layout_gen;
    }
    // Composite pressure solve of the stage: PRHS of every level, one solve, H (pred) / HS (corrector) of every level, the face gradient at the interface into the fine ghost cells.
    void solve_composite(bool pred)
    {
        ensure_comp();
        LevelRegistry& reg = *L.m_reg;
        const int nl = static_cast<int>(lc.size());
        composite_zero_solid_faces(reg, 1, nl, {"FVX", "FVY", "FVZ"}, -1);   // levels >= 1 only: level 0 keeps FDS's NO_FLUX FV (wall velocity, e.g. MASS_FLUX)   // no flux through a solid face: with U = 0 there, the face stays 0 (beta = 0 in the solve)
        level0_rhs_and_check();
        for (int l = 1; l < nl; ++l) composite_fine_pressure_rhs(*lc[l].F, *lc[l].lev, *crhs[l]);
        if (std::getenv("FDSTL_GHOSTDBG") && nl > 1) {   // composite right-hand side by level: its sum must vanish for a periodic domain
            for (int l = 0; l < nl; ++l) {
                const amrex::MultiFab& r = l == 0 ? *rhs : *crhs[l];
                const amrex::iMultiFab* cov = reg.covered_mask(l);
                const Level& lv = *lc[l].lev; const double vol = lv.dx[0] * lv.dx[1] * lv.dx[2];
                double s = 0.0, sa = 0.0, sc = 0.0;
                for (amrex::MFIter mfi(r); mfi.isValid(); ++mfi) {
                    auto a = r.const_array(mfi); amrex::Array4<const int> ca; if (cov) ca = cov->const_array(mfi);
                    amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { if (cov && ca(i, j, k) != 0) { sc += a(i, j, k) * vol; return; } s += a(i, j, k) * vol; sa += std::abs(a(i, j, k)) * vol; });
                }
                amrex::Print() << "RHSSUM " << (pred ? "pred" : "corr") << " L" << l << " sum(rhs*vol) uncovered " << s << " covered " << sc << " all " << s + sc << " sum|rhs|vol " << sa << "\n";
            }
        }
        if (std::getenv("FDSTL_GHOSTDBG") && nl > 1) {   // net outflow of FV through the boundary of the fine box, from the fine and from the coarse array (one fine box, x-z plane)
            const Level& lf = *lc[1].lev; const Level& l0v = *lc[0].lev;
            for (amrex::MFIter mfi((*lc[1].F)["FVX"]); mfi.isValid(); ++mfi) {
                const amrex::Box fb = mfi.validbox(); amrex::Box cb = fb; cb.coarsen(lf.ref_ratio_from_parent);
                auto fx = (*lc[1].F)["FVX"].const_array(mfi); auto fz = (*lc[1].F)["FVZ"].const_array(mfi);
                double bf = 0.0, nlo = 0.0, nhi = 0.0;
                for (int k = fb.smallEnd(2); k <= fb.bigEnd(2); ++k) { bf += (fx(fb.bigEnd(0), 0, k) - fx(fb.smallEnd(0) - 1, 0, k)) * lf.dx[2] * lf.dx[1]; nlo += fx(fb.smallEnd(0) - 1, 0, k); nhi += fx(fb.bigEnd(0), 0, k); }
                for (int i = fb.smallEnd(0); i <= fb.bigEnd(0); ++i) bf += (fz(i, 0, fb.bigEnd(2)) - fz(i, 0, fb.smallEnd(2) - 1)) * lf.dx[0] * lf.dx[1];
                amrex::Print() << "FLUXBND fine: net outflow " << bf << " (mean FVX at low face " << nlo / fb.length(2) << ", high face " << nhi / fb.length(2) << ")\n";
            }
            amrex::Box cb = lf.ba[0]; cb.coarsen(lf.ref_ratio_from_parent);
            for (amrex::MFIter m0((*lc[0].F)["FVX"]); m0.isValid(); ++m0) {
                auto cx = (*lc[0].F)["FVX"].const_array(m0); auto cz = (*lc[0].F)["FVZ"].const_array(m0);
                double bc = 0.0, clo = 0.0, chi = 0.0;
                for (int k = cb.smallEnd(2); k <= cb.bigEnd(2); ++k) { bc += (cx(cb.bigEnd(0), 0, k) - cx(cb.smallEnd(0) - 1, 0, k)) * l0v.dx[2] * l0v.dx[1]; clo += cx(cb.smallEnd(0) - 1, 0, k); chi += cx(cb.bigEnd(0), 0, k); }
                for (int i = cb.smallEnd(0); i <= cb.bigEnd(0); ++i) bc += (cz(i, 0, cb.bigEnd(2)) - cz(i, 0, cb.smallEnd(2) - 1)) * l0v.dx[0] * l0v.dx[1];
                {   // periodic wrap of level 0: FV at the low and high domain face
                    const amrex::Box dom = l0v.geom.Domain(); double a = 0, b = 0, c = 0, d = 0;
                    for (int k = dom.smallEnd(2); k <= dom.bigEnd(2); ++k) { a += cx(dom.smallEnd(0) - 1, 0, k); b += cx(dom.bigEnd(0), 0, k); }
                    for (int i = dom.smallEnd(0); i <= dom.bigEnd(0); ++i) { c += cz(i, 0, dom.smallEnd(2) - 1); d += cz(i, 0, dom.bigEnd(2)); }
                    amrex::Print() << "WRAPCHK level 0: sum FVX low " << a << " high " << b << "; sum FVZ low " << c << " high " << d << "\n";
                }
                amrex::Print() << "FLUXBND coarse: net outflow " << bc << " (mean FVX at low face " << clo / cb.length(2) << ", high face " << chi / cb.length(2) << ")\n";
            }
        }
        pb::PressureProblem p;
        for (int f = 0; f < 6; ++f) p.bc[f] = static_cast<pb::BC>(bc_type[f]);
        set_face_neumann(p);
        p.levels.resize(nl);
        for (int l = 0; l < nl; ++l) {
            pb::PressureLevel& pl = p.levels[l];
            const Level& lv = *lc[l].lev;
            pl.ba = lv.ba; pl.dm = lv.dm; pl.geom = lv.geom;
            pl.ref_ratio = l == 0 ? amrex::IntVect(1) : lv.ref_ratio_from_parent;
            pl.rhs = l == 0 ? rhs.get() : crhs[l].get();
            pl.phi = l == 0 ? phi.get() : cphi[l].get();
            if (!warm) pl.phi->setVal(0.0);
            if (l < static_cast<int>(ccls.size())) pl.cell_class = ccls[l].get();
        }
        pb::PressureOptions o;
        o.backend = pb::BackendKind::Auto;
        o.verbose = std::getenv("FDSTL_PBV") ? 1 : 0;
        o.check_residual = std::getenv("FDSTL_PBV") != nullptr || std::getenv("FDSTL_PBCHECK") != nullptr;   // the size-scaled residual check of Role 2 (residual_limit = max(residual_tol, residual_floor))
        o.trigger = pb::SolveRoutine;
        {
            std::vector<amrex::MultiFab*> rhs_levels(nl);
            for (int l = 0; l < nl; ++l) rhs_levels[l] = const_cast<amrex::MultiFab*>(p.levels[l].rhs);
            fold_composite(p, rhs_levels);
        }
        study_options(o);
        const pb::PressureResult r = pb::solve_pressure(p, o, pws.get());
        if (r.status != pb::Status::Ok && r.status != pb::Status::NotConverged) die(std::string("composite pressure solve: ") + pb::to_string(r.status) + ": " + r.message);
        study_row(pred, "all", r);
        ++crep.solves;
        for (const auto& ci : r.components) { crep.removed_mean = std::max(crep.removed_mean, std::abs(ci.removed_mean)); crep.removed_rel = std::max(crep.removed_rel, ci.removed_rel); }
        crep.iterations += r.backend_status.iterations;
        if (r.residual_checked) {
            ++crep.res_checked;
            crep.res_check_max = std::max(crep.res_check_max, r.residual_check);
            crep.res_limit_max = std::max(crep.res_limit_max, r.residual_limit);
            crep.res_floor_max = std::max(crep.res_floor_max, r.residual_floor);
            if (!r.residual_ok) ++crep.res_failed;
        }
        crep.backend = r.backend;
        std::vector<std::array<amrex::MultiFab*, 3>> grad(nl);
        std::vector<std::unique_ptr<amrex::MultiFab>> gown;
        for (int l = 0; l < nl; ++l)
            for (int d = 0; d < 3; ++d) {
                amrex::BoxArray fb = lc[l].lev->ba; fb.surroundingNodes(d);
                gown.emplace_back(new amrex::MultiFab(fb, lc[l].lev->dm, 1, 0));
                gown.back()->setVal(0.0);
                grad[l][d] = gown.back().get();
            }
        const pb::PressureResult rg = pb::face_gradient_composite(p, grad);
        if (rg.status != pb::Status::Ok) die(std::string("composite face gradient: ") + pb::to_string(rg.status) + ": " + rg.message);
        const char* hn = pred ? "H" : "HS";
        for (int l = 0; l < nl; ++l) {
            amrex::MultiFab& H = (*lc[l].F)[hn];
            const amrex::MultiFab& q = *p.levels[l].phi;
            for (amrex::MFIter mfi(H); mfi.isValid(); ++mfi) {
                auto h = H.array(mfi); auto a = q.const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { h(i, j, k) = a(i, j, k); });
            }
            H.FillBoundary(0, 1, H.nGrowVect(), lc[l].lev->geom.periodicity());
            if (l == 0) {
                select(0);
                each_local([&](int nm) { int fl[6]; physical_flags(nm - 1, fl); fds_p_h_ghost(nm, pred ? 1 : 0, fl); });
            } else {
                crep.h_ghost_cells += composite_set_fine_h_ghosts(*lc[l].lev, H, grad[l]);
                mirror_domain_edges(*lc[l].lev, H, -1);
                dirichlet_fine_h_ghosts(*lc[l].lev, H);
            }
        }
        (void)reg;
    }
    int pressure_scheme_composite(bool pred, double t, double dt, double& perr_out, double& verr_out)
    {
        // IP[0] (PRESSURE_ITERATIONS) is not iterated: the composite solve is one Poisson solve per stage (the FFT/MLMG solve is exact for the periodic and Neumann faces of Phase 3)
        LevelRegistry& reg = *L.m_reg;
        const int nl = static_cast<int>(lc.size());
        perr_out = 0.0; verr_out = 0.0;
        fds_p_iter_init();
        fds_p_iter_inc();
        select(0);
        if (IP[4]) for_levels([&](int) { each_local([&](int nm) { fds_p_baroclinic(t, nm); }); });   // BAROCLINIC_CORRECTION of every level (reads H/HS and RHO with the ghost layers of the level)
        L.m_bc->exchange_om();
        wrap_check("composite before match");
        each_local([&](int nm) { fds_g_match_flux(nm); });
        wrap_check("composite after match");
        for_levels([&](int) { each_local([&](int nm) { fds_p_noflux(dt, nm, 1); }); });
        for (int l = 1; l < nl; ++l) lc[l].bc->zero_wall_flux();   // NO_FLUX of a fine box at an INERT/MIRROR domain face (no wall cell there): FV = 0 on the wall face (patch 0017)
        wrap_check("composite after noflux");
        { ProfScope prof_ad(4); composite_average_down_flux(reg, nl); }
        wrap_check("composite after avg-down");
        select(0);
        each_local([&](int nm) { fds_p_rhs(t, dt, nm); });
        const double ts0 = amrex::second();
        solve_composite(pred);
        t_solve += amrex::second() - ts0;
        return 1;
    }
    // the density pieces of a stage over all levels: with the flux overwrite (several levels) the ADV read-out of every level, the override lists, then the cell update finest first
    void stage_density_all(bool pred, double t, double dt)
    {
        if (!multi() || !comp_overwrite) { for_levels([&](int) { density(pred, t, dt); }); return; }
        ensure_comp();
        FluxStages fs(L);
        const int nl = static_cast<int>(lc.size());
        {
            ProfScope prof_ow(6);
            for (int l = 0; l < nl; ++l) fs.compute_stage_fluxes(l, pred);
            runner->set_overrides(fs, fdsrt::FluxKind::Adv);
        }
        for (int l = nl - 1; l >= 0; --l) fs.apply_flux_divergence(l, pred);
    }
    // DIVERGENCE_PART_1 over all levels; with the flux overwrite: the DIF read-out of every level, the override lists, then the re-run where a list was set
    void divergence1_all()
    {
        if (!multi() || !comp_overwrite) { for_levels([&](int) { s_div1(); }); return; }
        ensure_comp();
        FluxStages fs(L);
        const int nl = static_cast<int>(lc.size());
        {
            ProfScope prof_ow(6);
            for (int l = 0; l < nl; ++l) fs.readout_dif(l);
            runner->set_overrides(fs, fdsrt::FluxKind::Dif);
        }
        for (int l = nl - 1; l >= 0; --l) fs.run_divergence_part1(l);
    }
    void avg_down_vel(bool pred)
    {
        if (!multi()) return;
        ProfScope prof_ad(4);
        composite_zero_solid_faces(*L.m_reg, 1, static_cast<int>(lc.size()), pred ? std::array<const char*, 3>{"US", "VS", "WS"} : std::array<const char*, 3>{"U", "V", "W"}, 0);
        composite_zero_solid_cells(*L.m_reg, 1, static_cast<int>(lc.size()), {"D", "DS", "DDDT"});   // solid fine cells: D/DS/DDDT stay 0 (div u = 0 there)
        if (pred) composite_average_down_faces(*L.m_reg, static_cast<int>(lc.size()), {"US", "VS", "WS"});
        else composite_average_down_faces(*L.m_reg, static_cast<int>(lc.size()), {"U", "V", "W"});
    }
    // Cause A (Role 2): an IC with nonzero wall-normal velocity on an INERT/MIRROR face has a net volume flux.
    // The Neumann projection removes the mean of that flux (removed_rel ~ 1) and the continuity update then grows mass at a
    // constant rate. Zero the normal velocity on Neumann faces (U/US, V/VS, W/WS) before the initial projection. Periodic and
    // Dirichlet (OPEN) faces are left alone. Level 0 VELOCITY_BC repeats this each step; fine boxes have no wall cells.
    void zero_closed_wall_normal()
    {
        const char* vel[3] = {"U", "V", "W"};
        const char* star[3] = {"US", "VS", "WS"};
        const char* divn[2] = {"D", "DS"};
        const int nl = static_cast<int>(lc.size());
        for (int l = 0; l < nl; ++l) {
            const amrex::Box dom = lc[l].lev->geom.Domain();
            const double rdx[3] = {1.0 / lc[l].lev->dx[0], 1.0 / lc[l].lev->dx[1], 1.0 / lc[l].lev->dx[2]};
            Fields& F = *lc[l].F;
            for (int d = 0; d < 3; ++d) {
                for (int side = 0; side < 2; ++side) {
                    if (bc_type[side * 3 + d] != static_cast<int>(pb::BC::Neumann)) continue;
                    const int face = side == 0 ? dom.smallEnd(d) : dom.bigEnd(d) + 1;
                    const char* names[2] = {vel[d], star[d]};
                    for (int which = 0; which < 2; ++which) {
                        if (!F.has(names[which]) || !F.has(divn[which])) continue;
                        amrex::MultiFab& mf = F[names[which]];
                        amrex::MultiFab& div = F[divn[which]];
                        for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
                            const amrex::Box vb = mfi.validbox();
                            if (side == 0 ? vb.smallEnd(d) != face : vb.bigEnd(d) != face) continue;
                            amrex::Box slab = vb;
                            slab.setSmall(d, face);
                            slab.setBig(d, face);
                            auto a = mf.array(mfi);
                            auto dd = div.array(mfi);
                            const double r = rdx[d];
                            amrex::LoopOnCpu(slab, [&](int i, int j, int k) {
                                const double wold = a(i, j, k);
                                if (wold == 0.0) return;
                                a(i, j, k) = 0.0;
                                // div uses (high - low)/dx. Dropping the low face raises div by w/dx; dropping the high face lowers it.
                                int c[3] = {i, j, k};
                                c[d] += side == 0 ? 0 : -1;
                                dd(c[0], c[1], c[2]) += (side == 0 ? wold : -wold) * r;
                            });
                        }
                    }
                }
            }
        }
    }
    // One corrector-form pass without the density update on all levels: the D, DS, KRES, MU, ... of a previous stage that a newly made level does not have (notes/level-binding.md)
    // FDSTL_GHOSTDBG: min/max of the scalars over the valid cells and over the whole FAB (ghost layers included), per level
    void ghost_report(const char* tag)
    {
        if (!std::getenv("FDSTL_GHOSTDBG")) return;
        nan_fields(tag);
        field_ranges(tag);
        for (int l = 0; l < static_cast<int>(lc.size()); ++l)
            for (const char* nm : {"RHO", "RHOS", "TMP", "RSUM", "ZZ", "ZZS"}) {
                const amrex::MultiFab& mf = (*lc[l].F)[nm];
                amrex::Print() << "GHOSTDBG " << tag << " L" << l << " " << nm << " valid [" << mf.min(0, 0) << ", " << mf.max(0, 0) << "] with ghosts [" << mf.min(0, mf.nGrow()) << ", " << mf.max(0, mf.nGrow()) << "]\n";
            }
    }
    // FDSTL_GHOSTDBG: max |div u - D*| over the uncovered cells by level, for the predicted (US.., DS) or the corrected (U.., D) velocity
    void div_report(const char* tag, bool pred)
    {
        if (!std::getenv("FDSTL_GHOSTDBG") || !multi()) return;
        (void)0;
        std::vector<fdsrt::ProjectionLevel> lv;
        for (int l = 0; l < static_cast<int>(lc.size()); ++l) {
            Fields& F = *lc[l].F;
            fdsrt::ProjectionLevel pl;
            pl.geom = L.m_reg->level(l).geom;
            pl.ref_ratio = l > 0 ? L.m_reg->level(l).ref_ratio_from_parent : amrex::IntVect(1);
            pl.vel = pred ? std::array<amrex::MultiFab*, 3>{&F["US"], &F["VS"], &F["WS"]} : std::array<amrex::MultiFab*, 3>{&F["U"], &F["V"], &F["W"]};
            pl.D = pred ? &F["DS"] : &F["D"];
            pl.covered = L.m_reg->covered_mask(l);
            lv.push_back(pl);
        }
        std::vector<double> bl;
        const double e = fdsrt::composite_divergence_error(lv, &bl);
        amrex::Print() << "DIVERR " << tag << " " << e << " (L0 " << bl[0] << ", L1 " << (bl.size() > 1 ? bl[1] : 0.0) << ")\n";
        if (std::atoi(std::getenv("FDSTL_GHOSTDBG")) >= 3)
            for (int l = 0; l < static_cast<int>(lc.size()); ++l) {   // where is the worst cell (x-z plane, one box per rank)
                Fields& F = *lc[l].F;
                const amrex::MultiFab& u = F[pred ? "US" : "U"]; const amrex::MultiFab& w = F[pred ? "WS" : "W"]; const amrex::MultiFab& d = F[pred ? "DS" : "D"];
                const amrex::Geometry& g = L.m_reg->level(l).geom;
                const amrex::iMultiFab* cov = L.m_reg->covered_mask(l);
                for (amrex::MFIter mfi(d); mfi.isValid(); ++mfi) {
                    auto ua = u.const_array(mfi); auto wa = w.const_array(mfi); auto da = d.const_array(mfi);
                    const amrex::Box b = mfi.validbox(); double worst = -1.0; int wi = 0, wk = 0;
                    amrex::Array4<const int> ca; if (cov) ca = cov->const_array(mfi);
                    amrex::LoopOnCpu(b, [&](int i, int j, int k) {
                        if (cov && ca(i, j, k) != 0) return;
                        const double dv = (ua(i + 1, j, k) - ua(i, j, k)) / g.CellSize(0) + (wa(i, j, k + 1) - wa(i, j, k)) / g.CellSize(2) - da(i, j, k);
                        if (std::abs(dv) > worst) { worst = std::abs(dv); wi = i; wk = k; }
                    });
                    if (pred) {   // does the stored pressure satisfy lap(H) = -div(FV) - DDDT on the cells away from the box edge?
                        const amrex::MultiFab& h = F["H"]; const amrex::MultiFab& fx = F["FVX"]; const amrex::MultiFab& fz = F["FVZ"]; const amrex::MultiFab& dd = F["DDDT"];
                        auto ha = h.const_array(mfi); auto fxa = fx.const_array(mfi); auto fza = fz.const_array(mfi); auto dda = dd.const_array(mfi);
                        double rw = -1.0, ri = -1.0; int ri_i = 0, ri_k = 0; const double dx = g.CellSize(0), dz = g.CellSize(2);
                        amrex::LoopOnCpu(b, [&](int i, int j, int k) {
                            if (cov && ca(i, j, k) != 0) return;
                            const double lap = (ha(i + 1, j, k) - 2 * ha(i, j, k) + ha(i - 1, j, k)) / (dx * dx) + (ha(i, j, k + 1) - 2 * ha(i, j, k) + ha(i, j, k - 1)) / (dz * dz);
                            const double res = lap + (fxa(i, j, k) - fxa(i - 1, j, k)) / dx + (fza(i, j, k) - fza(i, j, k - 1)) / dz + dda(i, j, k);
                            (void)res;
                            const double r2 = std::abs(lap + (fxa(i + 1, j, k) - fxa(i, j, k)) / dx + (fza(i, j, k + 1) - fza(i, j, k)) / dz + dda(i, j, k));
                            const double r1 = std::abs(res);
                            const bool inner = i > b.smallEnd(0) + 1 && i < b.bigEnd(0) - 1 && k > b.smallEnd(2) + 1 && k < b.bigEnd(2) - 1;
                            if (inner && std::min(r1, r2) > ri) { ri = std::min(r1, r2); ri_i = i; ri_k = k; }
                            if (!inner && std::min(r1, r2) > rw) rw = std::min(r1, r2);
                        });
                        amrex::Print() << "PRESRES " << tag << " L" << l << " interior worst " << ri << " at i " << ri_i << " k " << ri_k << ", box-edge cells worst " << rw << " (min over the two FV index conventions)\n";
                    }
                    if (pred) {
                        auto uo = F["U"].const_array(mfi); auto wo = F["W"].const_array(mfi); auto ha = F["H"].const_array(mfi); auto fxa = F["FVX"].const_array(mfi); auto fza = F["FVZ"].const_array(mfi);
                        const double dx = g.CellSize(0), dz = g.CellSize(2); const int i = wi, k = wk;
                        const double divU = (uo(i + 1, 0, k) - uo(i, 0, k)) / dx + (wo(i, 0, k + 1) - wo(i, 0, k)) / dz;
                        const double lap = (ha(i + 1, 0, k) - 2 * ha(i, 0, k) + ha(i - 1, 0, k)) / (dx * dx) + (ha(i, 0, k + 1) - 2 * ha(i, 0, k) + ha(i, 0, k - 1)) / (dz * dz);
                        const double dfA = (fxa(i + 1, 0, k) - fxa(i, 0, k)) / dx + (fza(i, 0, k + 1) - fza(i, 0, k)) / dz, dfB = (fxa(i, 0, k) - fxa(i - 1, 0, k)) / dx + (fza(i, 0, k) - fza(i, 0, k - 1)) / dz;
                        amrex::Print() << "CELLCHK " << tag << " L" << l << " (" << i << "," << k << ") divU " << divU << " lapH " << lap << " divFV(lo-face conv) " << dfA << " divFV(hi-face conv) " << dfB << " dt " << L.m_dt << " DS " << da(i, 0, k)
                                       << " => divUS expected " << divU - L.m_dt * (dfA + lap) << " / " << divU - L.m_dt * (dfB + lap) << "\n";
                    }
                    amrex::Print() << "DIVWORST " << tag << " L" << l << " box " << b << " worst " << worst << " at i " << wi << " k " << wk << "\n";
                }
            }
    }
    // FDSTL_GHOSTDBG=2: valid-cell min/max of the velocity, flux, pressure and divergence fields per level
    void field_ranges(const char* tag)
    {
        const char* e = std::getenv("FDSTL_GHOSTDBG");
        if (!e || std::atoi(e) < 2) return;
        for (int l = 0; l < static_cast<int>(lc.size()); ++l)
            for (const char* nm : {"U", "W", "US", "WS", "FVX", "FVZ", "H", "HS", "D", "DS", "DDDT", "MU", "KRES"}) {
                if (!lc[l].F->has(nm)) continue;
                const amrex::MultiFab& mf = (*lc[l].F)[nm];
                amrex::Print() << "RANGE " << tag << " L" << l << " " << nm << " [" << mf.min(0, 0) << ", " << mf.max(0, 0) << "]\n";
            }
    }
    void nan_fields(const char* tag)
    {
        for (int l = 0; l < static_cast<int>(lc.size()); ++l)
            for (const char* nm : {"RHO", "RHOS", "ZZ", "ZZS", "TMP", "D", "DS", "KRES", "MU", "RSUM", "Q", "U", "V", "W", "US", "VS", "WS", "H", "HS"}) {
                if (!lc[l].F->has(nm)) continue;
                const amrex::MultiFab& mf = (*lc[l].F)[nm];
                if (mf.contains_nan(0, mf.nComp(), 0)) amrex::Print() << "PRIMEDBG " << tag << ": level " << l << " " << nm << " has NaN in its valid cells\n";
            }
    }
    void nan_report(const char* tag)
    {
        if (!std::getenv("FDSTL_PRIMEDBG")) return;
        for (int l = 0; l < static_cast<int>(lc.size()); ++l) {   // the Fortran-side scratch arrays (fds_k_xfer mode 3)
            select(l);
            for (const auto& nr : std::initializer_list<std::pair<const char*, int>>{{"WORK1", 3}, {"WORK2", 3}, {"WORK3", 3}, {"WORK4", 3}, {"WORK5", 3}, {"WORK6", 3}, {"WORK7", 3}, {"WORK8", 3}, {"WORK9", 3},
                                                                                        {"SWORK1", 4}, {"SWORK2", 4}, {"SWORK3", 4}, {"SWORK4", 4}, {"DEL_RHO_D_DEL_Z", 4}, {"FX", 4}, {"FY", 4}, {"FZ", 4}})
                each_local([&](int nm) {
                    int lb[4] = {1, 1, 1, 1}, ub[4] = {1, 1, 1, 1}, ierr = 0; long n = 0, bad = 0;
                    fds_k_xfer(nm, nr.first, nr.second, lb, ub, nullptr, 2, &n, &bad, &ierr);
                    if (ierr != 0) return;
                    std::size_t cnt = 1;
                    bool sane = true;
                    for (int d = 0; d < nr.second; ++d) { const long ext = static_cast<long>(ub[d]) - lb[d] + 1; if (ext < 1 || ext > 100000) sane = false; else cnt *= static_cast<std::size_t>(ext); }
                    if (!sane || cnt > 100000000u) return;
                    std::vector<double> buf(cnt);
                    fds_k_xfer(nm, nr.first, nr.second, lb, ub, buf.data(), 3, &n, &bad, &ierr);
                    if (ierr != 0) return;
                    long nn = 0; for (double x : buf) if (std::isnan(x)) ++nn;
                    if (nn) amrex::Print() << "PRIMEDBG " << tag << ": level " << l << " scratch " << nr.first << " has " << nn << " NaN of " << cnt << "\n";
                });
        }
        select(0);
        nan_fields(tag);
    }
    // The state a newly bound level lacks before its first predictor: the ghost layers of the scalars (coarse-fine hook, domain-edge mirror) and of the velocities. The scalar and velocity
    // ghosts are those of the corrector end of a previous step (exchange 6 and 4). `full`: also the corrector-form stage pieces without the density update (viscosity, velocity flux, wall BC,
    // DIVERGENCE_PART_1/2) that fill D, DS, KRES, MU (FDSTL_PRIME_FULL=1; kept for the D_PBAR_DT binding investigation, notes/level-binding.md).
    void prime_levels()
    {
        const int nl = static_cast<int>(lc.size());
        const bool full = std::getenv("FDSTL_PRIME_FULL") != nullptr;
        auto probe = [&](const char* tag) {
            if (!std::getenv("FDSTL_TLDIAG") || nl < 2) return;
            const amrex::MultiFab& u = (*lc[0].F)["U"];
            for (amrex::MFIter mfi(u); mfi.isValid(); ++mfi) { auto a = u.const_array(mfi); amrex::Print() << "PROBE prime " << tag << ": U(5,0,-1) " << a(5, 0, -1) << " U(5,0,16) " << a(5, 0, 16) << "\n"; }
        };
        probe("entry");
        state(false, false);
        probe("after state");
        nan_report("start");
        for (int l = 1; l < nl; ++l) {
            composite_fill_fine_face_ghosts(*L.m_reg, l, {"U", "V", "W"});
            composite_fill_fine_face_ghosts(*L.m_reg, l, {"US", "VS", "WS"});
        }
        probe("after fine face ghosts");
        // level 0 keeps the ghost layers of its own initialisation: the corrector-form exchange/after_exchange 6 and 4 at t=0 (ICYC=0) overwrote the periodic ghost layer of U/W with wall-extrapolated values
        // (first predictor FVX/FVZ wrong at the periodic faces: the PRHS did not telescope, removed mean 1.6); only the new levels get them
        for (int l = 1; l < nl; ++l) { select(l); s_exchange(6, false); }
        probe("after exchange 6");
        for (int l = 1; l < nl; ++l) { select(l); s_after(6); }
        select(0);
        probe("after after-6");
        nan_report("exchange 6");
        if (full) for_levels([&](int) { s_visc_mfd(false); });
        nan_report("viscosity");
        for (int l = 1; l < nl; ++l) { select(l); s_exchange(4, false); }
        probe("after exchange 4");
        for (int l = 1; l < nl; ++l) { select(l); s_after(4); }
        select(0);
        probe("after after-4");
        nan_report("exchange 4");
        // The next predictor's DENSITY reads DEL_RHO_D_DEL_Z of the last corrector DIVERGENCE_PART_1: on a new fine box it does not exist and on level 0 it was made with the old
        // coarse-fine faces (DIF override of the old hierarchy), so the first step after a regrid moved rho*Z between species. Redo it on the new hierarchy (with the DIF override); keep D, DS.
        if (nl > 1 && !full && !std::getenv("FDSTL_NO_PRIME_DIV1")) {
            std::vector<std::unique_ptr<amrex::MultiFab>> keep;
            for (int l = 0; l < nl; ++l)
                for (const char* nm : {"D", "DS"}) { const amrex::MultiFab& f = (*lc[l].F)[nm]; keep.emplace_back(new amrex::MultiFab(f.boxArray(), f.DistributionMap(), f.nComp(), f.nGrow())); amrex::MultiFab::Copy(*keep.back(), f, 0, 0, f.nComp(), f.nGrow()); }
            zone_save();
            fds_p_init_div();   // as before the stage DIVERGENCE_PART_1 (the full-prime path does the same)
            // DIVERGENCE_PART_1 reads the wall-cell values of WALL_BC (after the velocity flux): stale ones from the old hierarchy moved rho*Z between species
            for_levels([&](int) { s_vflux(false); });
            for_levels([&](int) { s_wall_bc(false); });
            divergence1_all();
            zone_restore();
            int q = 0;
            for (int l = 0; l < nl; ++l)
                for (const char* nm : {"D", "DS"}) { amrex::MultiFab& f = (*lc[l].F)[nm]; amrex::MultiFab::Copy(f, *keep[q], 0, 0, f.nComp(), f.nGrow()); ++q; }
            nan_report("div1 (DEL_RHO_D_DEL_Z)");
        }
        if (!full) return;
        for_levels([&](int) { s_vflux(false); });
        nan_report("vflux");
        fds_p_init_div();
        for_levels([&](int) { s_wall_bc(false); });
        nan_report("wall_bc");
        for_levels([&](int) { s_div1(); });
        nan_report("div1");
        for_levels([&](int) { s_div2(); });
        nan_report("div2");
    }

    // FDSTL_GHOSTDBG: FVX/FVZ at the periodic wrap of level 0 (sums over the domain faces; a periodic level needs low == high for the PRHS to telescope)
    void wrap_check(const char* tag)
    {
        if (!std::getenv("FDSTL_GHOSTDBG")) return;
        const amrex::Box dom = lc[0].lev->geom.Domain();
        for (amrex::MFIter m0((*lc[0].F)["FVX"]); m0.isValid(); ++m0) {
            auto cx = (*lc[0].F)["FVX"].const_array(m0); auto cz = (*lc[0].F)["FVZ"].const_array(m0);
            double a = 0, b = 0, c = 0, d = 0;
            for (int k = dom.smallEnd(2); k <= dom.bigEnd(2); ++k) { a += cx(dom.smallEnd(0) - 1, 0, k); b += cx(dom.bigEnd(0), 0, k); }
            for (int i = dom.smallEnd(0); i <= dom.bigEnd(0); ++i) { c += cz(i, 0, dom.smallEnd(2) - 1); d += cz(i, 0, dom.bigEnd(2)); }
            amrex::Print() << "WRAP " << tag << " icyc " << L.m_icyc << ": sum FVX low " << a << " high " << b << "; sum FVZ low " << c << " high " << d << "\n";
        }
    }
    int pressure_scheme(bool pred, double t, double dt, double& perr_out, double& verr_out)
    {
        wrap_check(pred ? "pred scheme in" : "corr scheme in");
        if (multi()) return pressure_scheme_composite(pred, t, dt, perr_out, verr_out);
        Fields& F = *L.m_F;
        BcStep& bc = *L.m_bc;
        const char* hn = pred ? "H" : "HS";
        (void)hn;
        int iter = 0;
        double verr_old = 0.0, perr_old = 0.0;
        const int maxit = pred ? IP[2] : IP[1];
        fds_p_iter_init();
        for (;;) {
            fds_p_iter_inc();
            ++iter;
            const bool fvd = !pred && iter == 1;   // FVX/FVY/FVZ dumps of the first corrector iteration (stage_raw is silent unless FDSTL_STAGE selects the cycle)
            if (fvd) stage_raw("c_fv_pres_in", {"FVX", "FVY", "FVZ"});
            if (fds_p_iter_baro() || iter == 1) {
                if (IP[4]) each_local([&](int nm) { fds_p_baroclinic(t, nm); });
                bc.exchange_om();   // MESH_EXCHANGE(5): FVX, FVY, FVZ (and H/HS) to OMESH
                each_local([&](int nm) { fds_g_match_flux(nm); });
            }
            if (fvd) stage_raw("c_fv_matched", {"FVX", "FVY", "FVZ"});
            wrap_check("after match");
            each_local([&](int nm) { fds_p_noflux(dt, nm, iter == 1 ? 1 : 0); fds_p_rhs(t, dt, nm); });
            if (fvd) stage_raw("c_fv_noflux", {"FVX", "FVY", "FVZ"});
            solve_poisson(pred);
            each_local([&](int nm) { fds_p_resid(nm); });
            if (!IP[0]) break;
            F.fill_ghosts(pred ? "H" : "HS");
            bc.exchange_om();
            each_local([&](int nm) { fds_p_velerr(dt, nm); });
            double pe = 0.0, ve = 0.0;
            each_local([&](int nm) { double a, b, c; fds_p_get_err(nm, &a, &b, &c); pe = std::max(pe, a); ve = std::max(ve, b); });
            all_max(pe); all_max(ve);
            if (std::getenv("FDSTL_DEBUG")) { double ptb = 0; each_local([&](int nm) { double a, b, c; fds_p_get_err(nm, &a, &b, &c); ptb = std::max(ptb, c); }); amrex::Print() << "  [dbg] icyc " << L.m_icyc << (pred ? " P" : " C") << " pass " << passes << " iter " << iter << " dt=" << dt << " pe=" << pe << " ve=" << ve << " pois_err=" << ptb << "\n"; }
            perr_out = pe; verr_out = ve;
            if (pe < RP[0]) fds_p_set_iter_baro(0);
            if (pred && iter >= IP[2]) break;
            if (!pred && iter >= IP[1]) break;
            (void)maxit;
            if (pe < RP[0] && ve < RP[1]) break;
            if (IP[3] && L.m_icyc > 10) {
                if (iter > 3 && ve > RP[6] * verr_old && pe > RP[6] * perr_old) break;
                verr_old = ve; perr_old = pe;
            }
        }
        return iter;
    }

    // ------------------------------------------------------------ per-level stage bodies (S9): the exact calls advance() makes, on the current level c
    void s_visc_mfd(bool pred) { each_local([&](int nm) { fds_k_visc(nm, pred ? 0 : 1); fds_p_mfd(nm); }); }
    void s_exchange(int code, bool pred)
    {
        c->bc->exchange(code, pred);
        const int l = static_cast<int>(c - lc.data());
        if (l > 0 && (code == 3 || code == 6)) {
            ProfScope prof_cf(5);
            if (code == 3) composite_fill_fine_face_ghosts(*L.m_reg, l, {"US", "VS", "WS"});
            else composite_fill_fine_face_ghosts(*L.m_reg, l, {"U", "V", "W"});
        }
    }
    void s_after(int code) { c->bc->after_exchange(code, L.m_t, L.m_dt); }
    void s_vflux(bool pred)
    {
        iface(true);
        each_local([&](int nm) { fds_p_clear_attached(nm); fds_k_vflux(L.m_t, L.m_dt, nm, pred ? 0 : 1); });
        iface(false);
    }
    // Radiation ghost intensity of the fine boxes (radi.f90 RAD_GHOST, patch 0022): FR-062/L1251 table of every face ghost cell of every fine box (built once per layout;
    // the sender is a box of the receiver's level if one covers the ghost, else the coarser level), then the ILS -> ILR copy across ranks after each update.
    fdsrt::IlGhostTable rg_t; std::vector<int> rg_rank; double* rg_ils = nullptr; double* rg_ilr = nullptr; int rg_gen = -1, rg_nra = 0, rg_nsb = 0;
    void rad_ghost_table()   // before the radiation sweeps: a regrid must not leave the old table (old boxes' ghost cells) in use for the first update
    {
        if (lc.size() > 1 && rg_gen != layout_gen) {
            std::vector<fdsrt::IlBox> boxes; std::vector<int> nm; std::vector<amrex::IntVect> ratio; std::vector<fdsrt::IlWall> walls; std::vector<int> wijk;
            rg_rank.clear();
            for (std::size_t l = 0; l < lc.size(); ++l) {
                if (l > 0) ratio.push_back(lc[l].lev->ref_ratio_from_parent);
                for (int i = 0; i < lc[l].nbox; ++i) {
                    const amrex::Box b = lc[l].lev->ba[i];
                    boxes.push_back({static_cast<int>(l), b}); nm.push_back(lc[l].fds0 + i + 1); rg_rank.push_back(lc[l].lev->dm[i]);
                    if (l == 0) continue;
                    for (int d = 0; d < 3; ++d) for (int side = 0; side < 2; ++side) {   // face ghost cells, receiver order (walls of one box consecutive)
                        amrex::Box f = side ? amrex::adjCellHi(b, d) : amrex::adjCellLo(b, d);
                        for (amrex::BoxIterator it(f); it.ok(); ++it) {
                            const amrex::IntVect g = it();
                            walls.push_back({static_cast<int>(boxes.size()) - 1, side ? -(d + 1) : d + 1, g});
                            for (int q = 0; q < 3; ++q) wijk.push_back(g[q] - b.smallEnd(q) + 1);
                        }
                    }
                }
            }
            rg_t = fdsrt::build_il_ghost_table(boxes, ratio, walls);
            std::vector<int> wnm, wnic, snm, sijk, sloc;
            for (std::size_t w = 0; w < walls.size(); ++w) { wnm.push_back(nm[walls[w].box]); wnic.push_back(rg_t.nic_min[w]); wnic.push_back(rg_t.nic_max[w]); }
            for (int p = 0; p < rg_t.nicr(); ++p) { snm.push_back(nm[rg_t.ils_box[p] - 1]); sijk.push_back(rg_t.iio[p]); sijk.push_back(rg_t.jjo[p]); sijk.push_back(rg_t.kko[p]); sloc.push_back(rg_rank[rg_t.ils_box[p] - 1] == amrex::ParallelDescriptor::MyProc()); }
            fds_rad_ghost_set(rg_t.nicr(), static_cast<int>(walls.size()), *std::max_element(nm.begin(), nm.end()), wnm.data(), wijk.data(), wnic.data(), snm.data(), sijk.data(), sloc.data(),
                              &rg_ils, &rg_ilr, &rg_nra, &rg_nsb);
            rg_gen = layout_gen;   // new table: ILS holds the senders' isotropic intensity; exchanged at once (no ambient ghost on the first update)
            rad_ghost_exchange();
        }
    }
    void rad_ghost_exchange() { fdsrt::exchange_il_ghost(rg_t, rg_rank.data(), rg_nra, rg_nsb, rg_ils, rg_ilr); }
    void s_wall_bc(bool pred) { wseam_upload("WALL_BC"); c->bc->wall_bc(pred ? 1 : 0, L.m_t, L.m_dt); wseam_check("WALL_BC"); }
    void s_div1() { if (fluxchk && (fluxchk_kinds & 2)) { fluxchk_div1(); return; } div1_plain(); }
    void s_div2() { each_local([&](int nm) { fds_k_div2(L.m_dt, nm); }); }
    // velocity predictor: writes this rank's boxes of the level into the global per-box vectors at slot off + box (other slots untouched)
    void s_vpred(std::vector<double>& dn, std::vector<double>& ci)
    {
        each_local([&](int nm) {
            double dtn, cfl, vn; int ic;
            fds_k_vpred(L.m_t + L.m_dt, L.m_dt, nm, &dtn, &ic, &cfl, &vn);
            dn[c->off + bi(nm)] = dtn; ci[c->off + bi(nm)] = ic;
        });
    }
    void s_vcorr() { each_local([&](int nm) { fds_k_vcorr(L.m_t, L.m_dt, nm); }); }

    // ------------------------------------------------------------ one MAIN_LOOP iteration
    bool advance()
    {
        select(0);
        StepRecord rec;
        zone_rel_step = 0.0;
        reset_clip_budget();
        ++L.m_icyc;
        rec.step = L.m_icyc;
        L.m_dt = dt_at_step_start(L.m_t, L.m_dt, dt_new, chg, tc);
        rec.t0 = L.m_t;
        double& t = L.m_t;
        double& dt = L.m_dt;

        // ---------------- predictor
        state(true, true);
        for_levels([&](int) { s_visc_mfd(true); });
        div_report("start of step", false);
        first_pass = true;
        passes = 0;
        int stop = 0;
        for (;;) {
            ++passes;
            state(true, first_pass);
            if (passes == 1) {
                stage_scratch("p1_mfd", {{"FX", 4}, {"FY", 4}, {"FZ", 4}});
                stage_scratch("static", {{"RDX", 1}, {"RDY", 1}, {"RDZ", 1}, {"RDXN", 1}, {"RDYN", 1}, {"RDZN", 1}, {"R", 1}, {"RRN", 1}});
                stage_solid("static");
                stage_scratch("static", {{"MU_RSQMW_Z", 2}, {"K_RSQMW_Z", 2}, {"CP_Z", 2}, {"H_SENS_Z", 2}, {"RSQ_MW_Z", 1}, {"MWR_Z", 1}, {"MW", 1}}, true);
            }
            stage_raw(passes == 1 ? "p1_pre_dens" : "p2_pre_dens", {"RHO", "ZZ", "TMP", "U", "V", "W", "H", "HS", "MU", "KRES"});
            ghost_report("p_predens");
            stage_density_all(true, t, dt);
            ghost_report("p_dens");
            stage(passes == 1 ? "p1_dens" : "p2_dens", {"RHOS", "ZZS"});
            stage_raw(passes == 1 ? "p1_a_dens" : "p2_a_dens", {"RHOS", "TMP"});
            for_levels([&](int) { s_exchange(1, true); });
            ghost_report("p_exch1");
            stage_raw(passes == 1 ? "p1_b_fill" : "p2_b_fill", {"RHOS", "TMP"});
            for_levels([&](int) { s_after(1); });
            ghost_report(passes == 1 ? "p1_after1" : "p2_after1");
            stage_raw(passes == 1 ? "p1_c_visc" : "p2_c_visc", {"RHOS", "TMP"});
            stage_raw(passes == 1 ? "p1_prevflux" : "p2_prevflux", {"RHO", "RHOS", "U", "V", "W", "MU", "KRES", "H", "HS", "ZZ", "TMP"});
            if (const char* e = std::getenv("FDSTL_EDGES")) if (std::atoi(e) == L.m_icyc && passes == 1) each_local([&](int nm) { const std::string fn = dir + "/edges_p1_b" + std::to_string(nm) + ".bin"; fds_p_edge_dump(nm, fn.c_str(), static_cast<int>(fn.size())); });
            for_levels([&](int) { s_vflux(true); });
            stage(passes == 1 ? "p1_vflux" : "p2_vflux", {"FVX", "FVZ", "MU"});
            ghost_report("pd_vflux");
            fds_p_init_div();
            for_levels([&](int) { s_wall_bc(true); });
            ghost_report("pd_wallbc");
            if (std::getenv("FDSTL_WALLS") && L.m_icyc == std::atoi(std::getenv("FDSTL_WALLS")) && passes == 1)
                {   // SCRATCH PATCH: print plane FDSTL_WALLS_K (default 5, <0 = all planes); FDSTL_WALLTAB=1 also writes walltab_b<box>.bin (all external + internal walls)
                    const char* kk = std::getenv("FDSTL_WALLS_K"); const int k0 = kk ? std::atoi(kk) : 5;
                    each_local([&](int nm) { const amrex::Box b = l0.ba[nm - 1]; fds_p_wall_dump(nm, b.smallEnd(0), b.smallEnd(2), k0);
                        if (std::getenv("FDSTL_WALLTAB")) fds_p_wall_table_dump(nm, b.smallEnd(0), b.smallEnd(2)); });
                }
            stage_raw(passes == 1 ? "p1_prediv" : "p2_prediv", {"RHOS", "ZZS", "TMP", "RSUM", "U", "V", "W", "MU", "KRES", "D"});
            divergence1_all();
            ghost_report("pd_div1"); nan_report("pd_div1");
            stage(passes == 1 ? "p1_div1" : "p2_div1", {"DS", "MU", "KRES", "TMP", "RSUM"});
            if (passes == 1) stage_scratch("p1_div1", {{"WORK1", 3}, {"WORK2", 3}, {"WORK3", 3}, {"WORK4", 3}, {"WORK5", 3}, {"WORK6", 3}, {"WORK7", 3}, {"WORK8", 3}, {"WORK9", 3},
                                                       {"SWORK1", 4}, {"SWORK2", 4}, {"SWORK3", 4}, {"SWORK4", 4}, {"DEL_RHO_D_DEL_Z", 4}});
            zone_sums(true);
            ghost_report("pd_zone");
            for_levels([&](int) { s_div2(); });
            ghost_report("pd_div2");
            stage(passes == 1 ? "p1_div" : "p2_div", {"D", "DS", "DDDT"});
            state(true, first_pass);
            bool took = L.m_hook && L.m_hook(true, passes);
            if (!took) { const double tp0 = amrex::second(); rec.it_pred = pressure_scheme(true, t, dt, rec.perr_pred, rec.verr_pred); t_pres += amrex::second() - tp0; }
            // velocity predictor of every box, then the global DT decision (MIN over boxes and ranks)
            stage(passes == 1 ? "p1_press" : "p2_press", {"H", "FVX", "FVZ"});
            ghost_report("p_press");
            // (S9: the vectors hold every bound level, slot off + box; every box has one owner, so the Sum reduction is a gather; dt_at_step_start / dt_after_pass take the
            // MIN over all of them: one global dt over levels and ranks, the same in both stages, D-050)
            const int ntot = total_boxes();
            std::vector<double> dn(ntot, 0.0);
            std::vector<double> ci(ntot, 0.0);
            for_levels([&](int) { s_vpred(dn, ci); });
            amrex::ParallelAllReduce::Sum(dn.data(), ntot, amrex::ParallelContext::CommunicatorSub());
            amrex::ParallelAllReduce::Sum(ci.data(), ntot, amrex::ParallelContext::CommunicatorSub());
            for (int i = 0; i < ntot; ++i) { dt_new[i] = dn[i]; chg[i] = static_cast<int>(ci[i]); }
            stop = fds_p_stop_status();
            all_max(stop);
            bool nonfinite = false;
            for (double x : dt_new) if (!std::isfinite(x)) nonfinite = true;
            if (stop != 0 || nonfinite) { rec.dt = dt; L.m_steps.push_back(rec); return false; }
            if (!dt_after_pass(dt, dt_new, chg)) break;
            first_pass = false;
        }
        rec.passes = passes;
        stage("pred_end", {"US", "WS", "RHOS"});
        ghost_report("pred_end");
        avg_down_vel(true);
        div_report("after predictor", true);
        for_levels([&](int) { s_exchange(3, true); });
        for_levels([&](int) { s_after(3); });
        stage("pred_match", {"US", "WS"});

        // ---------------- corrector
        t = time_advance(t, dt);
        fds_p_zero_dot();   // MAIN_LOOP: Q_DOT = 0, M_DOT = 0 after T = T + DT
        state(false, first_pass);
        for_levels([&](int) { s_visc_mfd(false); });
        stage_raw("c_pre", {"RHOS", "TMP", "RHO", "US", "VS", "WS", "H", "HS"});
        stage_density_all(false, t, dt);
        stage("c_dens", {"RHO", "ZZ"});
        for_levels([&](int) { s_exchange(4, false); });
        for_levels([&](int) { s_after(4); });
        ghost_report("c_after4");
        stage_raw("c_prevflux", {"RHOS", "US", "VS", "WS", "MU", "KRES", "H", "HS", "ZZS", "TMP", "RHO"});
        for_levels([&](int) { s_vflux(false); });
        stage("c_vflux", {"FVX", "FVZ"});
        if (IP[9] > 0) for_levels([&](int) { each_local([&](int nm) { fds_k_comb(t, dt, nm); }); });   // main.f90 ~1060: COMBUSTION_LOAD_BALANCED after VELOCITY_FLUX, before WALL_BC; every mesh of every level
        // WP1b / 7.13b: edge tables + vorticity/stress work arrays right after VELOCITY_FLUX (corrector); FDSTL_STAGES=<icyc>
        stage_scratch("c_vfluxE", {{"WORK1", 3}, {"WORK2", 3}, {"WORK3", 3}, {"WORK4", 3}, {"WORK5", 3}, {"WORK6", 3}, {"DS", 3}, {"FVX", 3}, {"FVY", 3}, {"FVZ", 3}, {"ED_OMEGA", 2}, {"ED_TAU", 2}, {"CE_EDGE_INDEX", 4}});
        ++wall_counter;
        fds_p_set_wall_counter(wall_counter);
        for_levels([&](int) { s_wall_bc(false); });
        if (IP[10]) rad_ghost_table();   // fine-box ghost intensity table for the current layout
        if (IP[9] > 0 || IP[10] || IP[25]) for_levels([&](int) {
            each_local([&](int nm) { fds_k_rad(L.m_t, nm); });
            if (IP[10] && lc.size() > 1) rad_ghost_exchange();   // ILS -> ILR (all ranks) after each level: level L+1 sees this update's level-L intensity of the same angle subset; same-level faces lag one update, as MESH_EXCHANGE(2)
        });   // main.f90 ~1121-1142: the one COMPUTE_RADIATION (QR) and COMBUSTION_BC per box, after WALL_BC, before DIVERGENCE_PART_1
        if (IP[10]) fds_k_rte_factor();   // main.f90 CALCULATE_RTE_SOURCE_CORRECTION_FACTOR
        if (IP[10]) { double d; fds_setup(4, "", &d); }   // main.f90 MESH_EXCHANGE(2): lagged intensity exchange between level-0 meshes (IL_S -> IL_R)
        if (wall_counter == IP[21]) wall_counter = 0;
        fds_p_set_wall_counter(wall_counter);
        fds_p_init_div();
        stage_raw("c_fv_wallbc", {"FVX", "FVY", "FVZ"});   // FVX/FVY/FVZ bracket: where do they change after c_vflux (WP1b)
        divergence1_all();
        zone_sums(false);
        stage_raw("c_fv_div1", {"FVX", "FVY", "FVZ"});
        for_levels([&](int) { s_div2(); });
        state(false, first_pass);
        stage_raw("c_fv_div2", {"FVX", "FVY", "FVZ"});
        {
            bool took = L.m_hook && L.m_hook(false, 1);
            if (!took) { const double tp0 = amrex::second(); rec.it_corr = pressure_scheme(false, t, dt, rec.perr_corr, rec.verr_corr); t_pres += amrex::second() - tp0; }
        }
        stage("c_press", {"HS", "DS"});
        stage_raw("c_prevcorr", {"U", "V", "W", "US", "VS", "WS", "H", "HS", "FVX", "FVY", "FVZ"});   // velocity state seen by VELOCITY_CORRECTOR
        ghost_report("before_vcorr");
        for_levels([&](int) { s_vcorr(); });
        ghost_report("after_vcorr");
        avg_down_vel(false);
        div_report("after corrector", false);
        stage("c_end", {"U", "W"});
        for_levels([&](int) { s_exchange(6, false); });
        for_levels([&](int) { s_after(6); });
        rec.dt = dt;
        rec.zone_rel = zone_rel_step;
        L.m_steps.push_back(rec);
        outputs(rec);
        return true;
    }

    // ------------------------------------------------------------ outputs
    // Global array of one (component of a) registered field on rank 0: valid cells 1..n in FDS order (I fastest). For face fields FDS index I is the high face.
    std::vector<double> gather(const std::string& name, int comp)
    {
        Fields& F = *L.m_F;
        return gather_mf(F[name], F.spec(name), comp);
    }
    std::vector<double> gather_mf(const amrex::MultiFab& mf, const FieldSpec& s, int comp)
    {
        amrex::MultiFab tmp(l0.ba, l0.dm, 1, 0);
        for (amrex::MFIter mfi(tmp); mfi.isValid(); ++mfi) {
            auto t = tmp.array(mfi);
            auto a = mf.const_array(mfi);
            const int sh[3] = {s.nodal(0) ? 1 : 0, s.nodal(1) ? 1 : 0, s.nodal(2) ? 1 : 0};
            amrex::LoopOnCpu(mfi.validbox(), [&](int i, int j, int k) { t(i, j, k) = a(i + sh[0], j + sh[1], k + sh[2], comp); });
        }
        const amrex::Box dom = l0.geom.Domain();
        amrex::BoxArray ba1(dom);
        amrex::DistributionMapping dm1(ba1, 1);
        amrex::MultiFab one(ba1, dm1, 1, 0);
        one.setVal(0.0);
        one.ParallelCopy(tmp, 0, 0, 1);
        std::vector<double> out;
        if (amrex::ParallelDescriptor::MyProc() == 0) {
            out.resize(static_cast<std::size_t>(dom.numPts()));
            for (amrex::MFIter mfi(one); mfi.isValid(); ++mfi) {
                auto a = one.const_array(mfi);
                std::size_t p = 0;
                for (int k = dom.smallEnd(2); k <= dom.bigEnd(2); ++k)
                    for (int j = dom.smallEnd(1); j <= dom.bigEnd(1); ++j)
                        for (int i = dom.smallEnd(0); i <= dom.bigEnd(0); ++i) out[p++] = a(i, j, k);
            }
        }
        return out;
    }


    // Diagnostic (FDSTL_STAGE=<icyc>): dump the level fields after each stage of that cycle, to compare decompositions stage by stage
    void stage(const char* tag, std::initializer_list<const char*> fl)
    {
        const char* e = std::getenv("FDSTL_STAGE");
        if (!e || std::atoi(e) != L.m_icyc) return;
        for (const char* f : fl) {
            const std::vector<double> a = gather(f, 0);
            if (amrex::ParallelDescriptor::MyProc() != 0) continue;
            std::FILE* fp = std::fopen((dir + "/stage_" + tag + "_" + f + ".bin").c_str(), "wb");
            if (fp) { std::fwrite(a.data(), sizeof(double), a.size(), fp); std::fclose(fp); }
        }
    }


    // Diagnostic (FDSTL_STAGEG=<icyc>): raw FABs (valid + ghost cells, FDS/AMReX index of the box) of the named fields, one file per box, for ghost-layer comparisons across decompositions
    void stage_raw(const char* tag, std::initializer_list<const char*> fl)
    {
        const char* e = std::getenv("FDSTL_STAGEG");
        if (!e || std::atoi(e) != L.m_icyc) return;
        Fields& F = *L.m_F;
        for (const char* f : fl) {
            const amrex::MultiFab& mf = F[f];
            for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
                const amrex::FArrayBox& fab = mf[mfi];
                const amrex::Box b = fab.box();
                std::FILE* fp = std::fopen((dir + "/raw_" + tag + "_" + f + "_b" + std::to_string(mfi.index()) + ".bin").c_str(), "wb");
                if (!fp) continue;
                const amrex::Box v = mfi.validbox();
                int hdr[16] = {b.smallEnd(0), b.smallEnd(1), b.smallEnd(2), b.bigEnd(0), b.bigEnd(1), b.bigEnd(2), fab.nComp(), 0, v.smallEnd(0), v.smallEnd(1), v.smallEnd(2), v.bigEnd(0), v.bigEnd(1), v.bigEnd(2), 0, 0};
                std::fwrite(hdr, sizeof(int), 16, fp);
                std::fwrite(fab.dataPtr(), sizeof(double), static_cast<std::size_t>(b.numPts()) * fab.nComp(), fp);
                std::fclose(fp);
            }
        }
    }

    // Diagnostic (FDSTL_STAGES=<icyc>): Fortran-side scratch arrays that are not registered (FX FY FZ, WORK1..9, SWORK1..4, DEL_RHO_D_DEL_Z, CE_EDGE_INDEX, ED_OMEGA, ED_TAU), the mesh metrics and the species tables,
    // read through fds_k_xfer (mode 3) with the box's own allocation bounds. One file per box and array: scr_<tag>_<name>_b<box>.bin = int32 header[16] {rank, lb[4], ub[4], 0..}
    // followed by float64 data (first index fastest); the tables and metrics are written once per run under the tag "static" (box 0 only for the tables).
    void stage_scratch(const char* tag, std::initializer_list<std::pair<const char*, int>> names, bool tables = false)
    {
        const char* e = std::getenv("FDSTL_STAGES");
        if (!e || std::atoi(e) != L.m_icyc) return;
        for (const auto& nr : names)
            each_local([&](int nm) {
                if (tables && nm != 1) return;
                int lb[4] = {1, 1, 1, 1}, ub[4] = {1, 1, 1, 1}, ierr = 0; long n = 0, bad = 0;
                fds_k_xfer(nm, nr.first, nr.second, lb, ub, nullptr, 2, &n, &bad, &ierr);
                if (ierr != 0) return;
                std::size_t cnt = 1;
                for (int d = 0; d < nr.second; ++d) cnt *= static_cast<std::size_t>(ub[d] - lb[d] + 1);
                std::vector<double> buf(cnt);
                fds_k_xfer(nm, nr.first, nr.second, lb, ub, buf.data(), 3, &n, &bad, &ierr);
                if (ierr != 0) return;
                int hdr[16] = {nr.second, lb[0], lb[1], lb[2], lb[3], ub[0], ub[1], ub[2], ub[3], 0, 0, 0, 0, 0, 0, 0};
                std::FILE* fp = std::fopen((dir + "/scr_" + tag + "_" + nr.first + "_b" + std::to_string(nm - 1) + ".bin").c_str(), "wb");
                if (!fp) return;
                std::fwrite(hdr, sizeof(int), 16, fp);
                std::fwrite(buf.data(), sizeof(double), cnt, fp);
                std::fclose(fp);
            });
    }

    // FDSTL_STAGES: CELL%SOLID of every local box (SideData component 0, built from the FDS cells at set-up) as scr_<tag>_SOLID_b<box>.bin in the same layout as stage_scratch:
    // header int32[16] = {3, lb[4], ub[4], 0...} with the FDS index range 0..IBP1, 0..JBP1, 0..KBP1 (the valid cells and one ghost layer), then float64 (0 or 1), I fastest.
    // Ghost layer: SideData's fill (neighbour box value or periodic image; 0 at a closed domain edge), not FDS's own CELL(0,..)%SOLID.
    void stage_solid(const char* tag)
    {
        const char* e = std::getenv("FDSTL_STAGES");
        if (!e || std::atoi(e) != L.m_icyc) return;
        const amrex::iMultiFab& m = L.m_sd->mask();
        for (amrex::MFIter mfi(m); mfi.isValid(); ++mfi) {
            const amrex::Box vb = mfi.validbox();
            const auto a = m.const_array(mfi);
            const int n[3] = {vb.length(0) + 2, vb.length(1) + 2, vb.length(2) + 2};
            std::vector<double> buf(static_cast<std::size_t>(n[0]) * n[1] * n[2]);
            std::size_t q = 0;
            for (int k = 0; k < n[2]; ++k)
                for (int j = 0; j < n[1]; ++j)
                    for (int i = 0; i < n[0]; ++i) buf[q++] = a(vb.smallEnd(0) - 1 + i, vb.smallEnd(1) - 1 + j, vb.smallEnd(2) - 1 + k, 0);
            int hdr[16] = {3, 0, 0, 0, 1, n[0] - 1, n[1] - 1, n[2] - 1, 1, 0, 0, 0, 0, 0, 0, 0};
            std::FILE* fp = std::fopen((dir + "/scr_" + tag + "_SOLID_b" + std::to_string(mfi.index()) + ".bin").c_str(), "wb");
            if (!fp) continue;
            std::fwrite(hdr, sizeof(int), 16, fp);
            std::fwrite(buf.data(), sizeof(double), buf.size(), fp);
            std::fclose(fp);
        }
    }

    void write_mms(double t)
    {
        const std::vector<double> rho = gather("RHO", 0), z = gather("ZZ", 1), u = gather("U", 0), w = gather("W", 0), h = gather("H", 0);
        if (amrex::ParallelDescriptor::MyProc() != 0) return;
        const amrex::Box dom = l0.geom.Domain();
        std::FILE* f = std::fopen((dir + "/" + L.m_o.chid + "_mms.csv").c_str(), "w");
        if (!f) return;
        std::fprintf(f, "%d,%d,%d,%d,%d,%d\n", 1, dom.length(0), 1, dom.length(1), 1, dom.length(2));
        std::fprintf(f, "%22.15E\n", t);
        for (std::size_t c = 0; c < rho.size(); ++c) std::fprintf(f, "%22.15E,%22.15E,%22.15E,%22.15E,%22.15E\n", rho[c], z[c], u[c], w[c], h[c]);
        std::fclose(f);
    }

    std::string mass_row(double t)
    {
        Fields& F = *L.m_F;
        const double dv = l0.dx[0] * l0.dx[1] * l0.dx[2];
        std::vector<double> m(nt);
        double tot = 0.0;
        for (int n = 0; n < nt; ++n) { m[n] = exact_sum_product(F["RHO"], 0, F["ZZ"], n, dv); tot += m[n]; }
        char b[64];
        std::string s;
        std::snprintf(b, sizeof b, "%.7E,%.15E", t, tot); s += b;
        for (int n = 0; n < nt; ++n) { std::snprintf(b, sizeof b, ",%.15E", m[n]); s += b; }
        return s;
    }

    // Composite Q_ENTH, Q_CONV, Q_COND, Q_DIFF of the HRR file: UPDATE_HRR (inside fds_setup(3)) sums them over the level-0 meshes only. The walls: add to Q_DOT the
    // sums over the finer cells not covered by a still finer level, minus the covered coarser cells, so each region is counted once. Q_ENTH: replace UPDATE_HRR's
    // d/dt of the level-0 integral of rho*h_s by d/dt of the composite one, E (energy_sums). After a regrid (prime_levels) e_old is E on the new hierarchy, so the
    // regrid change of E (fine-from-coarse interpolation, coarse under removed fine cells) is left out of Q_ENTH.
    double e_old = 0.0;
    bool have_e_old = false;
    // E (composite) and E0 (level 0 alone, as UPDATE_HRR); w = composite wall terms minus the level-0 ones (Q_DOT(4:6)).
    void energy_sums(double& e, double& e0, double w[3])
    {
        const int nl = static_cast<int>(lc.size());
        e = e0 = 0.0; w[0] = w[1] = w[2] = 0.0;
        for (int l = 0; l < nl; ++l) {
            const amrex::MultiFab& q = (*lc[l].F)["Q"];
            amrex::iMultiFab cov(q.boxArray(), q.DistributionMap(), 1, 0);
            if (l + 1 < nl) cov = amrex::makeFineMask(q.boxArray(), q.DistributionMap(), (*lc[l + 1].F)["Q"].boxArray(), L.m_reg->level(l + 1).ref_ratio_from_parent, 0, 1);
            else cov.setVal(0);
            for (amrex::MFIter mfi(cov); mfi.isValid(); ++mfi) {
                const amrex::Box vb = mfi.validbox();
                auto c = cov.const_array(mfi);
                const int nm = lc[l].fds0 + mfi.index() + 1;
                std::vector<int> m(static_cast<std::size_t>(vb.numPts()));
                std::size_t n = 0;
                double t[4];
                if (l == 0) {   // level 0: all cells (E0), then remove the covered ones
                    std::fill(m.begin(), m.end(), 1);
                    fds_p_esums(nm, m.data(), t);
                    e0 += t[0]; e += t[0];
                    amrex::LoopOnCpu(vb, [&](int i, int j, int k) { m[n++] = c(i, j, k) != 0 ? 1 : 0; });
                    fds_p_esums(nm, m.data(), t);
                    e -= t[0]; for (int a = 0; a < 3; ++a) w[a] -= t[1 + a];
                } else {        // finer: add the cells not covered by the next level
                    amrex::LoopOnCpu(vb, [&](int i, int j, int k) { m[n++] = c(i, j, k) != 0 ? 0 : 1; });
                    fds_p_esums(nm, m.data(), t);
                    e += t[0]; for (int a = 0; a < 3; ++a) w[a] += t[1 + a];
                }
            }
        }
    }
    void energy_snapshot()
    {
        if (lc.size() < 2 || !L.m_reg || !fds_outputs) return;
        double e0, w[3];
        energy_sums(e_old, e0, w);
        have_e_old = true;
    }
    void composite_energy(double dt)
    {
        double e, e0, w[3];
        energy_sums(e, e0, w);
        double es_old = 0.0;   // UPDATE_HRR's ENTHALPY_SUM_OLD (local level-0 meshes)
        for (amrex::MFIter mfi((*lc[0].F)["Q"]); mfi.isValid(); ++mfi) es_old += fds_p_enthalpy_sum(lc[0].fds0 + mfi.index() + 1);
        for (int a = 0; a < 3; ++a) fds_p_add_qdot(4 + a, w[a]);
        if (L.m_icyc > 0 && have_e_old) fds_p_add_qdot(9, ((e - e_old) - (e0 - es_old)) / dt);
        e_old = e; have_e_old = true;
    }

    void outputs(const StepRecord& r)
    {
        log_steps.push_back([&] { char b[256]; std::snprintf(b, sizeof b, "%d,%.17E,%.17E,%d,%d,%d,%.6E,%.6E,%.6E,%.6E,%.3E", r.step, L.m_t, r.dt, r.passes, r.it_pred, r.it_corr, r.perr_pred, r.verr_pred, r.perr_corr, r.verr_corr, r.zone_rel); return std::string(b); }());
        log_mass.push_back(mass_row(L.m_t));
        if (fds_outputs) {
            // FDS's own writers (CHID_devc.csv, _hrr.csv, _mass.csv, _steps.csv, _cpu.csv and the .out step block) on the box data; main.f90 patch 0006
            const double to0 = amrex::second();
            fds_hook_set_step(L.m_t, r.dt, L.m_icyc);
            // Composite HRR: UPDATE_HRR sums Q*VC over the level-0 meshes, whose covered cells carry the coarse combustion of the averaged-down state. For the writers
            // only, the covered level-0 Q is replaced by the volume average of the finer Q (uncovered L0 + all finer cells, each counted once) and restored afterwards.
            // Same for QR (Q_RADI): the covered level-0 QR is the coarse radiation of the averaged-down state, while HRR is the composite one.
            std::vector<amrex::MultiFab> qsave;
            const int nlq = static_cast<int>(lc.size());
            std::vector<const char*> qnames;
            if (nlq > 1 && L.m_reg) for (const char* nm : {"Q", "QR"}) if (lc[0].F->has(nm)) qnames.push_back(nm);
            if (!qnames.empty() && lc[0].F->has("QR"))   // a fine box's radiation writes its own QR array, not the MultiFab: copy it in first
                for (int l = 1; l < nlq; ++l)
                    for (amrex::MFIter mfi((*lc[l].F)["QR"]); mfi.isValid(); ++mfi) {
                        amrex::FArrayBox& f = (*lc[l].F)["QR"][mfi];
                        fds_p_get_qr(lc[l].fds0 + mfi.index() + 1, f.dataPtr(), static_cast<long>(f.box().numPts()));
                    }
            for (const char* nm : qnames)
                for (int l = 0; l < nlq - 1; ++l) { const amrex::MultiFab& q = (*lc[l].F)[nm]; qsave.emplace_back(q.boxArray(), q.DistributionMap(), 1, 0); amrex::MultiFab::Copy(qsave.back(), q, 0, 0, 1, 0); }
            for (const char* nm : qnames)
                for (int l = nlq - 2; l >= 0; --l) amrex::average_down((*lc[l + 1].F)[nm], (*lc[l].F)[nm], 0, 1, L.m_reg->level(l + 1).ref_ratio_from_parent);
            if (nlq > 1 && L.m_reg) composite_energy(r.dt);
            double dtd = 0.0;
            fds_setup(3, "", &dtd);
            std::size_t qi = 0;
            for (const char* nm : qnames)
                for (int l = 0; l < nlq - 1; ++l) amrex::MultiFab::Copy((*lc[l].F)[nm], qsave[qi++], 0, 0, 1, 0);
            t_out += amrex::second() - to0;
        }
        if (!mms_done && (IP[6] == 7 || IP[6] == 11) && L.m_t >= RP[2]) { write_mms(L.m_t); mms_done = true; }
        if (!L.m_o.quiet)
            amrex::Print() << "step " << r.step << " T=" << L.m_t << " DT=" << r.dt << " passes=" << r.passes << " iters=" << r.it_pred << "/" << r.it_corr << " perr=" << r.perr_corr
                           << " verr=" << r.verr_corr << " zone_rel=" << r.zone_rel << "\n";
    }

    // NFR-030 / A-22 / NFR-031 measurements, one file per run (rank 0): wall seconds are the maximum over ranks, memory is summed over ranks.
    static double proc_status_kb(const char* key)
    {
        std::ifstream f("/proc/self/status");
        std::string line;
        const std::string k = std::string(key) + ":";
        while (std::getline(f, line))
            if (line.compare(0, k.size(), k) == 0) return std::atof(line.c_str() + k.size());
        return -1.0;
    }
    void write_perf()
    {
        double mx[4] = {t_loop, t_pres, t_solve, t_out};
        amrex::ParallelAllReduce::Max(mx, 4, amrex::ParallelContext::CommunicatorSub());
        double pr[7];
        std::copy(g_prof, g_prof + 7, pr);
        amrex::ParallelAllReduce::Max(pr, 7, amrex::ParallelContext::CommunicatorSub());
        double mem[3] = {proc_status_kb("VmHWM"), proc_status_kb("VmRSS"), static_cast<double>(L.m_F->bytes())};
        amrex::ParallelAllReduce::Sum(mem, 3, amrex::ParallelContext::CommunicatorSub());
        double hwm_max = proc_status_kb("VmHWM");
        amrex::ParallelAllReduce::Max(hwm_max, amrex::ParallelContext::CommunicatorSub());
        if (amrex::ParallelDescriptor::MyProc() != 0) return;
        const amrex::Box dom = l0.geom.Domain();
        std::FILE* f = std::fopen((dir + "/" + L.m_o.chid + "_driver_perf.csv").c_str(), "w");
        if (!f) return;
        std::fprintf(f, "ranks,boxes,nx,ny,nz,steps,t_loop_s,t_pressure_s,t_poisson_solve_s,t_fds_outputs_s,f_pres,f_solve,peak_rss_sum_kb,peak_rss_max_kb,rss_end_sum_kb,field_bytes_sum\n");
        std::fprintf(f, "%d,%d,%d,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%.4f,%.4f,%.0f,%.0f,%.0f,%.0f\n", amrex::ParallelDescriptor::NProcs(), nbox, dom.length(0), dom.length(1),
                     dom.length(2), L.m_icyc, mx[0], mx[1], mx[2], mx[3], mx[0] > 0 ? mx[1] / mx[0] : 0.0, mx[0] > 0 ? mx[2] / mx[0] : 0.0, mem[0], hwm_max, mem[1], mem[2]);
        std::fclose(f);
        if (std::getenv("FDSTL_PROFILE")) {   // diagnostic: where the ghost-layer wall time goes (max over ranks)
            std::FILE* pf = std::fopen((dir + "/" + L.m_o.chid + "_driver_profile.txt").c_str(), "w");
            if (pf) {
                std::fprintf(pf, "t_loop %.6f\nt_pressure %.6f\nt_poisson_solve %.6f\nt_fds_outputs %.6f\nghost_fill_boundary %.6f\nfill_omesh %.6f\nmatch_periodic_faces %.6f\nfds_bc_routines %.6f\n"
                                 "avg_down %.6f\ncf_ghost_fill %.6f\nflux_overwrite %.6f\n",
                             mx[0], mx[1], mx[2], mx[3], pr[0], pr[1], pr[2], pr[3], pr[4], pr[5], pr[6]);
                std::fclose(pf);
            }
        }
        // manifest of the field dump (the _final_<FIELD>.bin files below)
        std::ofstream m(dir + "/" + L.m_o.chid + "_final_manifest.txt");
        m << "# field dump for comparison: <chid>_final_<NAME>.bin, float64 little-endian, valid cells of the whole level, FDS order (I fastest, then J, then K)\n"
          << "# face fields (U,V,W): index I is the high face of cell I; ZZ<n>: scalar n of ZZ (ZZ1 .. ZZ" << ns << ")\n"
          << "nx=" << dom.length(0) << "\nny=" << dom.length(1) << "\nnz=" << dom.length(2) << "\nT=" << std::setprecision(17) << L.m_t << "\nstep=" << L.m_icyc << "\n"
          << "fields=U,V,W,H,HS,US,VS,WS,D,DS,RHO,TMP";
        for (int c = 0; c < ns; ++c) m << ",ZZ" << (c + 1);
        m << "\n";
    }

    void finish()
    {
        if (std::getenv("FDSTL_DIVERR")) {   // D-058 diagnostic (FluxStages::max_divergence_error), level 0, no effect on the results
            FluxStages fs(L);
            double e = fs.max_divergence_error(0);
            all_max(e);
            const double dmax = L.m_reg->fields(0)["D"].norminf(0, 0, false);   // scale of D (valid cells)
            amrex::Print() << "DIVERR level 0: max |div u - D| over the valid cells = " << std::setprecision(6) << e << " (1/s), max |D| = " << dmax << "\n";
        }
        if (fluxchk) {
            double mx0 = fluxchk_max[0], mx1 = fluxchk_max[1];
            all_max(mx0); all_max(mx1);
            amrex::Print() << "FLUXCHK mode " << fluxchk << ": ADV " << fluxchk_calls[0] << " calls (this rank " << fluxchk_faces[0] << " faces listed, max |flux| " << mx0 << "), DIF " << fluxchk_calls[1]
                           << " calls (" << fluxchk_faces[1] << " faces, max |flux| " << mx1 << ")\n";
        }
        if (amrex::ParallelDescriptor::MyProc() == 0) {
            std::FILE* f = std::fopen((dir + "/" + L.m_o.chid + "_driver_steps.csv").c_str(), "w");
            if (f) {
                std::fprintf(f, "step,T,DT,passes,it_pred,it_corr,perr_pred,verr_pred,perr_corr,verr_corr,zone_rel\n");
                for (auto& s : log_steps) std::fprintf(f, "%s\n", s.c_str());
                std::fclose(f);
            }
            f = std::fopen((dir + "/" + L.m_o.chid + "_driver_mass.csv").c_str(), "w");
            if (f) {
                std::fprintf(f, "Time,Total");
                for (int n = 0; n < nt; ++n) std::fprintf(f, ",S%d", n + 1);
                std::fprintf(f, "\n");
                for (auto& s : log_mass) std::fprintf(f, "%s\n", s.c_str());
                std::fclose(f);
            }
        }
        write_perf();
        // final fields (valid cells, FDS order, float64, written by rank 0) for the comparison with the baseline restart file
        static const char* const names[] = {"U", "V", "W", "H", "HS", "US", "VS", "WS", "D", "DS", "RHO", "TMP"};
        for (const char* nme : names) {
            const std::vector<double> a = gather(nme, 0);
            if (amrex::ParallelDescriptor::MyProc() == 0) {
                std::FILE* f = std::fopen((dir + "/" + L.m_o.chid + "_final_" + nme + ".bin").c_str(), "wb");
                if (f) { std::fwrite(a.data(), sizeof(double), a.size(), f); std::fclose(f); }
            }
        }
        for (int c = 0; c < ns; ++c) {
            const std::vector<double> a = gather("ZZ", c);
            if (amrex::ParallelDescriptor::MyProc() == 0) {
                std::FILE* f = std::fopen((dir + "/" + L.m_o.chid + "_final_ZZ" + std::to_string(c + 1) + ".bin").c_str(), "wb");
                if (f) { std::fwrite(a.data(), sizeof(double), a.size(), f); std::fclose(f); }
            }
        }
    }
};

// ====================================================================== public interface
TimeLoop::TimeLoop(const Level0& l0, double dt_setup, const RunOptions& o) : m_l0(l0), m_o(o)
{
    m.reset(new Impl(*this));
    m_dt = dt_setup;
    m->dir = o.outdir;
    m->setup();
}
TimeLoop::~TimeLoop() = default;

void TimeLoop::set_state(double t, double dt, int icyc) { m_t = t; m_dt = dt; m_icyc = icyc; }

void TimeLoop::bind_fields(bool copy_setup_state)
{
    Fields& F = *m_F;
    for (const auto& s : field_table()) {
        if (!F.has(s.name)) continue;
        amrex::MultiFab& mf = F[s.name];
        for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
            const int nm = mfi.index() + 1;
            amrex::FArrayBox& fab = mf[mfi];
            const amrex::Box cb = m_l0.ba[mfi.index()];
            if (copy_setup_state) {
                const FdsBounds w = fds_window(s, cb);
                const int ncomp = fab.nComp();
                const long ntot = static_cast<long>(w.ext[0]) * w.ext[1] * w.ext[2] * ncomp;
                std::vector<double> buf(ntot);
                int ierr = 0;
                fds_p_get(nm, s.name, buf.data(), ntot, &ierr);
                if (ierr != 0) die(std::string("cannot copy the FDS set-up state of ") + s.name + " (ierr " + std::to_string(ierr) + ")");
                const FdsView v = make_fds_view(s, fab, cb, true);
                long p = 0;
                for (int n = 1; n <= ncomp; ++n)
                    for (int K = v.lb[2]; K < v.lb[2] + v.ext[2]; ++K)
                        for (int J = v.lb[1]; J < v.lb[1] + v.ext[1]; ++J)
                            for (int I = v.lb[0]; I < v.lb[0] + v.ext[0]; ++I) v(I, J, K, n) = buf[p++];
            }
            const FdsBounds nb = fds_bounds(s, cb);
            int lb[4] = {nb.lb[0], nb.lb[1], nb.lb[2], 1};
            int ext[4] = {nb.ext[0], nb.ext[1], nb.ext[2], fab.nComp()};
            long stride[4] = {1, nb.ext[0], static_cast<long>(nb.ext[0]) * nb.ext[1], static_cast<long>(nb.ext[0]) * nb.ext[1] * nb.ext[2]};
            if (fds_shim_bind(nm, s.name, fab.dataPtr(), lb, ext, stride) != 0) die(std::string("cannot bind ") + s.name);
        }
    }
}

// ---- S9 per-level entry points (thin wrappers on the Impl stage bodies; level 0 is the current level again on return)
int TimeLoop::num_levels() const { return static_cast<int>(m->lc.size()); }
void TimeLoop::stage_state(bool predictor, bool first_pass) { m->state(predictor, first_pass); }
void TimeLoop::stage_init_divergence() { fds_p_init_div(); }
void TimeLoop::bind_level(int lev) { m->bind_level(lev); }
void TimeLoop::unbind_level(int lev) { m->unbind_level(lev); }
bool TimeLoop::level_bound(int lev) const { return lev >= 0 && lev < static_cast<int>(m->lc.size()); }
#define FDSTL_LEVEL_STAGE(lev, body) do { m->select(lev); body; m->select(0); } while (0)
void TimeLoop::stage_viscosity(int lev, bool predictor) { FDSTL_LEVEL_STAGE(lev, m->s_visc_mfd(predictor)); }
void TimeLoop::stage_density(int lev, bool predictor) { FDSTL_LEVEL_STAGE(lev, m->density(predictor, m_t, m_dt)); }
void TimeLoop::stage_exchange(int lev, int code, bool predictor) { FDSTL_LEVEL_STAGE(lev, m->s_exchange(code, predictor)); }
void TimeLoop::stage_boundary(int lev, int code) { FDSTL_LEVEL_STAGE(lev, m->s_after(code)); }
void TimeLoop::stage_velocity_flux(int lev, bool predictor) { FDSTL_LEVEL_STAGE(lev, m->s_vflux(predictor)); }
void TimeLoop::stage_wall_bc(int lev, bool predictor) { FDSTL_LEVEL_STAGE(lev, m->s_wall_bc(predictor)); }
void TimeLoop::stage_divergence1(int lev) { FDSTL_LEVEL_STAGE(lev, m->s_div1()); }
void TimeLoop::stage_divergence2(int lev) { FDSTL_LEVEL_STAGE(lev, m->s_div2()); }
void TimeLoop::stage_velocity_correct(int lev) { FDSTL_LEVEL_STAGE(lev, m->s_vcorr()); }
// ---- S10.2 interface flux hooks (see TimeLoop.H and notes/flux-hooks-design.md)
void TimeLoop::flux_readout_adv(int lev) { FDSTL_LEVEL_STAGE(lev, m->flux_readout_adv(m_t, m_dt)); }
void TimeLoop::flux_readout_dif(int lev) { FDSTL_LEVEL_STAGE(lev, m->flux_readout_dif()); }
void TimeLoop::flux_apply_dif(int lev) { FDSTL_LEVEL_STAGE(lev, m->flux_apply_dif()); }
void TimeLoop::flux_set_override(int lev, int kind, const std::vector<std::vector<fdsrt::FluxOverride>>& per_local_box)
{
    if (kind < 0 || kind > 1) die("flux_set_override: kind must be 0 (ADV) or 1 (DIF)");
    FDSTL_LEVEL_STAGE(lev, m->flux_set_override(kind, per_local_box));
}
const amrex::MultiFab& TimeLoop::flux_array(int lev, int kind, int dir)
{
    if (kind < 0 || kind > 1 || dir < 0 || dir > 2) die("flux_array: bad kind or direction");
    m->select(lev);
    const amrex::MultiFab& r = *m->stage_flux_arrays().mf[kind][dir];
    m->select(0);
    return r;
}
void TimeLoop::stage_velocity_update(int lev, bool predictor, double* local_dt, int* change_time_step)
{
    if (!predictor) die("stage_velocity_update: the corrector velocity update is stage_velocity_correct");
    m->select(lev);
    std::vector<double> dn(m->total_boxes(), 1.0e300), ci(m->total_boxes(), 0.0);
    m->s_vpred(dn, ci);
    double mn = 1.0e300; int ch = 0;
    for (int i = 0; i < m->c->nbox; ++i) if (m->local(i)) { mn = std::min(mn, dn[m->c->off + i]); ch |= static_cast<int>(ci[m->c->off + i]); }
    *local_dt = mn; *change_time_step = ch;
    if (static_cast<int>(m_last_local_dt.size()) <= lev) m_last_local_dt.resize(lev + 1, 1.0e300);
    m_last_local_dt[lev] = mn;
    m->select(0);
}
#undef FDSTL_LEVEL_STAGE
double TimeLoop::global_dt(const std::vector<double>& local_dt_per_level) const
{
    double mn = 1.0e300;
    for (double d : local_dt_per_level) mn = std::min(mn, d);
    amrex::ParallelAllReduce::Min(mn, amrex::ParallelContext::CommunicatorSub());
    return mn;
}
void TimeLoop::set_cf_ghost_hook(int lev, CfGhostHook h)
{
    m->select(lev);
    m->c->bc->cf_ghost_hook = std::move(h);
    m->select(0);
}

const TimeLoop::CompositeReport& TimeLoop::composite_report() const { return m->crep; }
void TimeLoop::pressure_csv_div(int step, const std::vector<double>& div_by_level)
{
    for (std::size_t l = 0; m->pcsv && l < div_by_level.size(); ++l) std::fprintf(m->pcsv, "%d,end,%zu,,,,,,,,,,,,,,%.3e\n", step, l, div_by_level[l]);
}
double TimeLoop::clip_species(int n) const { return n >= 0 && n < static_cast<int>(m->clip_species_acc.size()) ? m->clip_species_acc[n] : 0.0; }
double TimeLoop::clip_mass(int n, int l) const
{
    if (n < 1 || l < 0 || m->nt <= 0) return 0.0;
    const int idx = l * m->nt + (n - 1);
    if (idx >= static_cast<int>(m->clip_mass_acc.size())) return 0.0;
    return m->clip_mass_acc[idx];
}
long TimeLoop::clip_cells(int n, int l) const
{
    if (n < 1 || l < 0 || m->nt <= 0) return 0;
    const int idx = l * m->nt + (n - 1);
    if (idx >= static_cast<int>(m->clip_cells_acc.size())) return 0;
    return m->clip_cells_acc[idx];
}
void TimeLoop::prime_levels() { m->prime_levels(); m->energy_snapshot(); }
void TimeLoop::zero_closed_wall_normal() { m->zero_closed_wall_normal(); }
const int* TimeLoop::pressure_bc() const { return m->bc_type; }
std::array<std::vector<int>, 6> TimeLoop::pressure_face_neumann() const
{
    pb::PressureProblem p;
    for (int f = 0; f < 6; ++f) p.bc[f] = static_cast<pb::BC>(m->bc_type[f]);
    m->set_face_neumann(p);
    return p.face_neumann;
}
void TimeLoop::set_flux_overwrite(bool on) { m->comp_overwrite = on; m->comp_gen = -1; }
const fdsrt::StageStats* TimeLoop::flux_match() const { return m->runner ? &m->runner->stats() : nullptr; }
void TimeLoop::reset_flux_match() { if (m->runner) m->runner->reset_match(); }
void TimeLoop::write_perf() { m->write_perf(); }

bool TimeLoop::advance() { const double t0 = amrex::second(); const bool ok = m->advance(); m->t_loop += amrex::second() - t0; return ok; }

int TimeLoop::run()
{
    m->log_mass.push_back(m->mass_row(m_t));
    int fails = 0;
    for (;;) {
        if (m_o.max_steps >= 0 && m_icyc >= m_o.max_steps) break;
        const double tl0 = amrex::second();
        const bool adv_ok = m->advance();
        m->t_loop += amrex::second() - tl0;
        if (!adv_ok) { ++fails; amrex::Print() << "time loop: STOP (instability or non-finite DT) at step " << m_icyc << "\n"; break; }
        if (m_t >= m->tc.t_end && m_icyc > 0) break;
    }
    m->finish();
    return fails;
}

// ---- FluxStages: fdsrt::FluxAccess (D-058: max_divergence_error) on a TimeLoop (phase mapping documented in TimeLoop.H)
void FluxStages::compute_stage_fluxes(int level, bool) { m_tl.flux_readout_adv(level); }
void FluxStages::apply_flux_divergence(int level, bool predictor)
{
    // Density update only (consumes the ADV override list when one is set). DIVERGENCE_PART_1 is NOT run here: in FDS order it follows the exchange of RHOS/ZZS, the boundary steps, the
    // velocity flux and WALL_BC (and on several levels the coarse-fine ghost hook has to have filled the fine ghosts of the stage first), so the DIF read-out is a separate phase:
    // readout_dif(level) on every level, runner.set_overrides(Dif), run_divergence_part1(level) on every level (notes/flux-stage-wiring.md steps 3 to 6; D-061).
    m_tl.stage_density(level, predictor);
    if (level < static_cast<int>(m_div1_pending.size())) m_div1_pending[level] = 0;
}
void FluxStages::readout_dif(int level)
{
    m_tl.flux_readout_dif(level);           // DIVERGENCE_PART_1 plus the DIF read-out; stage_flux(level, Dif, d) is valid afterwards
    if (static_cast<int>(m_div1_pending.size()) <= level) m_div1_pending.resize(level + 1, 0);
    m_div1_pending[level] = 1;
}
const amrex::MultiFab& FluxStages::stage_flux(int level, fdsrt::FluxKind kind, int dir) const
{
    return m_tl.flux_array(level, kind == fdsrt::FluxKind::Adv ? 0 : 1, dir);
}
void FluxStages::set_flux_override(int level, fdsrt::FluxKind kind, const std::vector<std::vector<fdsrt::FluxOverride>>& per_local_box)
{
    m_tl.flux_set_override(level, kind == fdsrt::FluxKind::Adv ? 0 : 1, per_local_box);
}
double FluxStages::local_dt(int level) const
{
    if (level < 0 || level >= static_cast<int>(m_tl.m_last_local_dt.size())) return 1.0e300;
    return m_tl.m_last_local_dt[level];
}
void FluxStages::run_divergence_part1(int level)
{
    if (level < static_cast<int>(m_div1_pending.size()) && m_div1_pending[level]) { m_div1_pending[level] = 0; m_tl.flux_apply_dif(level); }
    else m_tl.stage_divergence1(level);
}
double FluxStages::max_divergence_error(int level) const
{
    // D-058 diagnostic: max over the valid cells of `level` on this rank of |div u - D|, with u the face velocities U, V, W of the level's Fields (AMReX face a = low face of cell a, the
    // divergence is the uniform-Cartesian difference with the level's dx) and D the divergence array. After a full step D equals div u up to the pressure tolerance.
    if (!m_tl.m_reg || !m_tl.m_reg->has_level(level)) amrex::Abort("FluxStages::max_divergence_error: no such level");
    Fields& F = m_tl.m_reg->fields(level);
    const Level& lev = m_tl.m_reg->level(level);
    double mx = 0.0;
    for (amrex::MFIter mfi(F["D"]); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        const auto d = F["D"].const_array(mfi);
        const auto u = F["U"].const_array(mfi); const auto v = F["V"].const_array(mfi); const auto w = F["W"].const_array(mfi);
        for (int k = vb.smallEnd(2); k <= vb.bigEnd(2); ++k)
            for (int j = vb.smallEnd(1); j <= vb.bigEnd(1); ++j)
                for (int i = vb.smallEnd(0); i <= vb.bigEnd(0); ++i) {
                    const double div = (u(i + 1, j, k) - u(i, j, k)) / lev.dx[0] + (v(i, j + 1, k) - v(i, j, k)) / lev.dx[1] + (w(i, j, k + 1) - w(i, j, k)) / lev.dx[2];
                    mx = std::max(mx, std::abs(div - d(i, j, k)));
                }
    }
    return mx;
}


}  // namespace fdsamr
