// TimeLoop.cpp: see TimeLoop.H.
//
// Kernel-facing rules (M2a):
//  (a) Passive scalars: N_TOTAL_SCALARS beyond the tracked species are handled by Fields.cpp (ZZ/ZZS carry ncomp = N_TOTAL_SCALARS); every copy, exchange and clip
//      here loops over all components, the D-031 clip special-cases only the tracked species.
//  (b) Only uniform Cartesian metrics are used (R(I)/RRN(I) dropped); CYLINDRICAL and TRN* meshes are rejected (IR-002), the pressure call-through assumes IPS=0.
#include "TimeLoop.H"
#include "LevelRegistry.H"
#include "PressureBcMap.H"

#include <AMReX_BLassert.H>
#include <AMReX_ParallelContext.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
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
void fds_k_div2(double dt, int nm);
void fds_k_vpred(double t, double dt, int nm, double* dtnew, int* ichg, double* cfl, double* vn);
void fds_k_vcorr(double t, double dt, int nm);
void fds_k_consts(double* tend, double* dtfill, double* dtmin, double* dt0, int* nzone);
void fds_k_set_flags(int nm, int f);
void fds_g_match_flux(int nm);
// end-of-step FDS outputs (patch 0006: fds_setup(mode=3) -> UPDATE_GLOBAL_OUTPUTS, EXCHANGE/DUMP_GLOBAL_OUTPUTS, WRITE_DIAGNOSTICS) and the hook module
void fds_setup(int mode, const char* fname, double* dt_out);
int fds_hook_step_outputs();
void fds_hook_set_step(double t, double dt, int icyc);
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
void fds_p_zero_dot();
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

extern double g_prof[4];   // wall seconds of the ghost-layer machinery (Fields.cpp), written by write_perf when FDSTL_PROFILE is set

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
    int passes = 0;
    int bc_type[6] = {0, 0, 0, 0, 0, 0};   // pb::BC per domain face (side*3+dir)
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
        static const char* names[22] = {"U", "V", "W", "US", "VS", "WS", "D", "DS", "H", "HS", "KRES", "FVX", "FVY", "FVZ", "RHO", "RHOS", "MU", "TMP", "Q", "RSUM", "ZZ", "ZZS"};
        for (int i = 0; i < 22; ++i) if (n == names[i]) return i + 1;
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
        own.resize(lev);
        OwnLevel& o = own[lev - 1];
        o.drho.reset(new amrex::MultiFab(Lv.ba, Lv.dm, 1, 1));
        o.dzz.reset(new amrex::MultiFab(Lv.ba, Lv.dm, std::max(1, ns), 0));
        o.drho->setVal(0.0); o.dzz->setVal(0.0);
        o.bc.reset(new BcStep(Lv, F));
        o.bc->ext_ghost = true;   // requirement of levels > 0 (BcStep::exchange): the driver saves the face velocities, MATCH_VELOCITY is not run
        LevelCtx x;
        x.lev = &Lv; x.F = &F; x.bc = o.bc.get(); x.sd = &L.m_reg->side_data(lev); x.drho = o.drho.get(); x.dzz = o.dzz.get();
        x.nbox = nb; x.fds0 = nm0; x.off = total_boxes();
        lc.push_back(x);      // may move lc: c is reset below
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
        own.resize(lev - 1);
        if (static_cast<int>(sfx.size()) > lev) sfx.resize(lev);   // the stage flux arrays of the level are re-registered with the next bind
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
        chk(IP[7], "CC_IBM"); chk(IP[8], "HVAC"); chk(IP[9] > 0, "reactions"); chk(IP[10], "radiation"); chk(IP[11] != 0, "LEVEL_SET_MODE");
        chk(IP[12], "TUNNEL_PRECONDITIONER"); chk(IP[13] != 1, "pressure solver other than FFT"); chk(IP[14], "CYLINDRICAL"); chk(IP[15], "particle drag");
        chk(IP[16], "SOLID_HEAT_TRANSFER_3D"); chk(IP[17] > 0, "Lagrangian particles"); chk(IP[18], "FREEZE_VELOCITY"); chk(IP[19], "SOLID_PHASE_ONLY");
        chk(IP[22], "MPI particle exchange"); chk(IP[23], "obstruction mass exchange"); chk(IP[24] > 0, "agglomeration"); chk(IP[25], "INIT_HRRPUV");
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
        const amrex::Box dom = l0.geom.Domain();
        for (int i = 0; i < nbox; ++i) {
            if (!local(i)) continue;
            const amrex::Box vb = l0.ba[i];
            const int c3[3] = {mip[16 * i + 3], mip[16 * i + 4], mip[16 * i + 5]};
            for (int d = 0; d < 3; ++d) {
                nn[d] = dom.length(d);
                if (vb.smallEnd(d) == dom.smallEnd(d)) code[d] = c3[d];
                if (vb.bigEnd(d) == dom.bigEnd(d)) code[3 + d] = c3[d];
            }
        }
        if (std::getenv("FDSTL_PBV")) std::fprintf(stderr, "[rank %d] codes before reduce %d %d %d %d %d %d\n", amrex::ParallelDescriptor::MyProc(), code[0], code[1], code[2], code[3], code[4], code[5]);
        for (int f = 0; f < 6; ++f) all_max(code[f]);
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
            for (amrex::MFIter mfi(RP_); mfi.isValid(); ++mfi) {
                const amrex::Box fb = RZ[mfi].box(), rb = RP_[mfi].box(), mb = mask[mfi].box(), v = mfi.validbox();
                mkdx(fb, dxv);
                for (int n = 1; n <= NT; ++n) {
                    int f, ccnt;
                    fds_clip_species(fb.loVect(), fb.hiVect(), rb.loVect(), rb.hiVect(), mb.loVect(), mb.hiVect(), v.loVect(), v.hiVect(), NS, n, RP_[mfi].dataPtr(),
                                     RZ[mfi].dataPtr(), mask[mfi].dataPtr(), dxv[0].data(), dxv[1].data(), dxv[2].data(), (*c->dzz)[mfi].dataPtr(n - 1), &f, &ccnt);
                    zfl[n - 1] |= f;
                }
            }
            for (int n = 0; n < NT; ++n) { bool b = zfl[n]; amrex::ParallelAllReduce::Or(b, amrex::ParallelContext::CommunicatorSub()); zfl[n] = b; anyz = anyz || b; }
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
        clip_level(pred, flags);
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
        // boundary data of the Poisson problem must be homogeneous on the faces that are solved (S5 scope); checked every solve
        {
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
            if (bm > 0.0) die("inhomogeneous Poisson boundary data (BXS..BZF /= 0): needs the S6 pressure interface extension (see notes/pressure-iface-review.md)");
        }
        if (std::getenv("FDSTL_PBV")) std::fprintf(stderr, "[rank %d] bc_type %d %d %d %d %d %d\n", amrex::ParallelDescriptor::MyProc(), bc_type[0], bc_type[1], bc_type[2], bc_type[3], bc_type[4], bc_type[5]);
        pb::PressureProblem p;
        p.ba = l0.ba; p.dm = l0.dm; p.geom = l0.geom;
        for (int f = 0; f < 6; ++f) p.bc[f] = static_cast<pb::BC>(bc_type[f]);
        p.rhs = rhs.get(); p.phi = phi.get();
        pb::PressureOptions o;
        o.backend = pb::BackendKind::FFT;
        o.verbose = std::getenv("FDSTL_PBV") ? 1 : 0;
        o.check_residual = std::getenv("FDSTL_PBV") != nullptr;
        phi->setVal(0.0);
        const pb::PressureResult r = pb::solve_pressure(p, o);
        if (r.status != pb::Status::Ok && r.status != pb::Status::NotConverged) die(std::string("pressure interface: ") + pb::to_string(r.status) + ": " + r.message);
        if (std::getenv("FDSTL_PBV")) { amrex::Print() << "  [pb] backend=" << r.backend << " res_rel2=" << r.residual_rel2 << " relmax=" << r.residual_relmax << " rhs max=" << rhs->norminf(0) << " phi max=" << phi->norminf(0); for (auto& c : r.components) amrex::Print() << " removed_mean=" << c.removed_mean << " rel=" << c.removed_rel; for (auto& w : r.warnings) amrex::Print() << " W:" << w; amrex::Print() << "\n"; }
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

    int pressure_scheme(bool pred, double t, double dt, double& perr_out, double& verr_out)
    {
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
    void s_exchange(int code, bool pred) { c->bc->exchange(code, pred); }
    void s_after(int code) { c->bc->after_exchange(code, L.m_t, L.m_dt); }
    void s_vflux(bool pred)
    {
        iface(true);
        each_local([&](int nm) { fds_p_clear_attached(nm); fds_k_vflux(L.m_t, L.m_dt, nm, pred ? 0 : 1); });
        iface(false);
    }
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
        if (lc.size() > 1) die("advance(): level " + std::to_string(lc.size() - 1) + " is bound but the composite pressure solve across levels is not built (Role 2); drive the per-level stage entry points instead");
        select(0);
        StepRecord rec;
        zone_rel_step = 0.0;
        ++L.m_icyc;
        rec.step = L.m_icyc;
        L.m_dt = dt_at_step_start(L.m_t, L.m_dt, dt_new, chg, tc);
        rec.t0 = L.m_t;
        double& t = L.m_t;
        double& dt = L.m_dt;

        // ---------------- predictor
        state(true, true);
        for_levels([&](int) { s_visc_mfd(true); });
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
            for_levels([&](int) { density(true, t, dt); });
            stage(passes == 1 ? "p1_dens" : "p2_dens", {"RHOS", "ZZS"});
            stage_raw(passes == 1 ? "p1_a_dens" : "p2_a_dens", {"RHOS", "TMP"});
            for_levels([&](int) { s_exchange(1, true); });
            stage_raw(passes == 1 ? "p1_b_fill" : "p2_b_fill", {"RHOS", "TMP"});
            for_levels([&](int) { s_after(1); });
            stage_raw(passes == 1 ? "p1_c_visc" : "p2_c_visc", {"RHOS", "TMP"});
            stage_raw(passes == 1 ? "p1_prevflux" : "p2_prevflux", {"RHO", "RHOS", "U", "V", "W", "MU", "KRES", "H", "HS", "ZZ", "TMP"});
            if (const char* e = std::getenv("FDSTL_EDGES")) if (std::atoi(e) == L.m_icyc && passes == 1) each_local([&](int nm) { const std::string fn = dir + "/edges_p1_b" + std::to_string(nm) + ".bin"; fds_p_edge_dump(nm, fn.c_str(), static_cast<int>(fn.size())); });
            for_levels([&](int) { s_vflux(true); });
            stage(passes == 1 ? "p1_vflux" : "p2_vflux", {"FVX", "FVZ", "MU"});
            fds_p_init_div();
            for_levels([&](int) { s_wall_bc(true); });
            if (std::getenv("FDSTL_WALLS") && L.m_icyc == std::atoi(std::getenv("FDSTL_WALLS")) && passes == 1)
                {   // SCRATCH PATCH: print plane FDSTL_WALLS_K (default 5, <0 = all planes); FDSTL_WALLTAB=1 also writes walltab_b<box>.bin (all external + internal walls)
                    const char* kk = std::getenv("FDSTL_WALLS_K"); const int k0 = kk ? std::atoi(kk) : 5;
                    each_local([&](int nm) { const amrex::Box b = l0.ba[nm - 1]; fds_p_wall_dump(nm, b.smallEnd(0), b.smallEnd(2), k0);
                        if (std::getenv("FDSTL_WALLTAB")) fds_p_wall_table_dump(nm, b.smallEnd(0), b.smallEnd(2)); });
                }
            stage_raw(passes == 1 ? "p1_prediv" : "p2_prediv", {"RHOS", "ZZS", "TMP", "RSUM", "U", "V", "W", "MU", "KRES", "D"});
            for_levels([&](int) { s_div1(); });
            stage(passes == 1 ? "p1_div1" : "p2_div1", {"DS", "MU", "KRES", "TMP", "RSUM"});
            if (passes == 1) stage_scratch("p1_div1", {{"WORK1", 3}, {"WORK2", 3}, {"WORK3", 3}, {"WORK4", 3}, {"WORK5", 3}, {"WORK6", 3}, {"WORK7", 3}, {"WORK8", 3}, {"WORK9", 3},
                                                       {"SWORK1", 4}, {"SWORK2", 4}, {"SWORK3", 4}, {"SWORK4", 4}, {"DEL_RHO_D_DEL_Z", 4}});
            zone_sums(true);
            for_levels([&](int) { s_div2(); });
            stage(passes == 1 ? "p1_div" : "p2_div", {"D", "DS", "DDDT"});
            state(true, first_pass);
            bool took = L.m_hook && L.m_hook(true, passes);
            if (!took) { const double tp0 = amrex::second(); rec.it_pred = pressure_scheme(true, t, dt, rec.perr_pred, rec.verr_pred); t_pres += amrex::second() - tp0; }
            // velocity predictor of every box, then the global DT decision (MIN over boxes and ranks)
            stage(passes == 1 ? "p1_press" : "p2_press", {"H", "FVX", "FVZ"});
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
        for_levels([&](int) { s_exchange(3, true); });
        for_levels([&](int) { s_after(3); });
        stage("pred_match", {"US", "WS"});

        // ---------------- corrector
        t = time_advance(t, dt);
        fds_p_zero_dot();   // MAIN_LOOP: Q_DOT = 0, M_DOT = 0 after T = T + DT
        state(false, first_pass);
        for_levels([&](int) { s_visc_mfd(false); });
        stage_raw("c_pre", {"RHOS", "TMP", "RHO", "US", "VS", "WS", "H", "HS"});
        for_levels([&](int) { density(false, t, dt); });
        stage("c_dens", {"RHO", "ZZ"});
        for_levels([&](int) { s_exchange(4, false); });
        for_levels([&](int) { s_after(4); });
        stage_raw("c_prevflux", {"RHOS", "US", "VS", "WS", "MU", "KRES", "H", "HS", "ZZS", "TMP", "RHO"});
        for_levels([&](int) { s_vflux(false); });
        stage("c_vflux", {"FVX", "FVZ"});
        ++wall_counter;
        fds_p_set_wall_counter(wall_counter);
        for_levels([&](int) { s_wall_bc(false); });
        if (wall_counter == IP[21]) wall_counter = 0;
        fds_p_set_wall_counter(wall_counter);
        fds_p_init_div();
        stage_raw("c_fv_wallbc", {"FVX", "FVY", "FVZ"});   // FVX/FVY/FVZ bracket: where do they change after c_vflux (WP1b)
        for_levels([&](int) { s_div1(); });
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
        for_levels([&](int) { s_vcorr(); });
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

    // Diagnostic (FDSTL_STAGES=<icyc>): Fortran-side scratch arrays that are not registered (FX FY FZ, WORK1..9, SWORK1..4, DEL_RHO_D_DEL_Z), the mesh metrics and the species tables,
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

    void outputs(const StepRecord& r)
    {
        log_steps.push_back([&] { char b[256]; std::snprintf(b, sizeof b, "%d,%.17E,%.17E,%d,%d,%d,%.6E,%.6E,%.6E,%.6E,%.3E", r.step, L.m_t, r.dt, r.passes, r.it_pred, r.it_corr, r.perr_pred, r.verr_pred, r.perr_corr, r.verr_corr, r.zone_rel); return std::string(b); }());
        log_mass.push_back(mass_row(L.m_t));
        if (fds_outputs) {
            // FDS's own writers (CHID_devc.csv, _hrr.csv, _mass.csv, _steps.csv, _cpu.csv and the .out step block) on the box data; main.f90 patch 0006
            const double to0 = amrex::second();
            fds_hook_set_step(L.m_t, r.dt, L.m_icyc);
            double dtd = 0.0;
            fds_setup(3, "", &dtd);
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
        double pr[4] = {g_prof[0], g_prof[1], g_prof[2], g_prof[3]};
        amrex::ParallelAllReduce::Max(pr, 4, amrex::ParallelContext::CommunicatorSub());
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
                std::fprintf(pf, "t_loop %.6f\nt_pressure %.6f\nt_poisson_solve %.6f\nt_fds_outputs %.6f\nghost_fill_boundary %.6f\nfill_omesh %.6f\nmatch_periodic_faces %.6f\nfds_bc_routines %.6f\n",
                             mx[0], mx[1], mx[2], mx[3], pr[0], pr[1], pr[2], pr[3]);
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
