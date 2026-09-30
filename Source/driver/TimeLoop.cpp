// TimeLoop.cpp: see TimeLoop.H.
//
// Kernel-facing rules (M2a):
//  (a) Passive scalars: N_TOTAL_SCALARS beyond the tracked species are handled by Fields.cpp (ZZ/ZZS carry ncomp = N_TOTAL_SCALARS); every copy, exchange and clip
//      here loops over all components, the D-031 clip special-cases only the tracked species.
//  (b) Only uniform Cartesian metrics are used (R(I)/RRN(I) dropped); CYLINDRICAL and TRN* meshes are rejected (IR-002), the pressure call-through assumes IPS=0.
#include "TimeLoop.H"

#include <AMReX_BLassert.H>
#include <AMReX_ParallelContext.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <initializer_list>

#include "ExactSum.H"
#include "FdsSetup.H"
#include "PressureIface.H"

extern "C" {
int fds_shim_bind(int nm, const char* name, void* base, const int* lb, const int* ext, const long* stride);
void fds_k_state(int pred, int first, int icyc, double rmin, double rmax);
void fds_k_visc(int nm, int est);
void fds_k_dens_post(double t, double dt, int nm);
void fds_k_vflux(double t, double dt, int nm, int est);
void fds_k_div1(double t, double dt, int nm);
void fds_k_div2(double dt, int nm);
void fds_k_vpred(double t, double dt, int nm, double* dtnew, int* ichg, double* cfl, double* vn);
void fds_k_vcorr(double t, double dt, int nm);
void fds_k_consts(double* tend, double* dtfill, double* dtmin, double* dt0, int* nzone);
void fds_k_set_flags(int nm, int f);
void fds_g_match_flux(int nm);
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
void fds_p_dens_pre(double t, double dt, int nm);
void fds_p_init_div();
void fds_p_zone_get(int n, double* d, double* p, double* u);
void fds_p_zone_set(int n, const double* d, const double* p, const double* u);
void fds_p_wall_dump(int nm, int ofx, int ofz, int kg0);
void fds_hook_set_flag(int flag);
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
    int passes = 0;
    int bc_type[6] = {0, 0, 0, 0, 0, 0};   // pb::BC per domain face (side*3+dir)
    bool solve_dir[3] = {true, true, true};
    std::string dir;
    std::vector<std::string> log_steps;
    std::vector<std::string> log_mass;

    explicit Impl(TimeLoop& l) : L(l), l0(l.m_l0) { std::memset(IP, 0, sizeof IP); std::memset(RP, 0, sizeof RP); }
    bool local(int i) const { return is_local(l0, i); }

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
        L.bind_fields(true);
        pressure_bc_map();
        zone_setup();
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
        for (int d = 0; d < 3; ++d) {
            const int n = dom.length(d);
            solve_dir[d] = !(d == 1 && IP[20] && n == 1);   // TWO_D: no y solve
            const bool per = l0.dom.periodic[d] != 0;
            const int lo = code[d], hi = code[3 + d];
            if (lo < 0 || hi < 0) die("no box found at a domain face (pressure boundary map)");
            if (per) {
                // a periodic level direction: the level-wide Poisson problem is periodic whatever the per-mesh FDS code is (single mesh: code 0; several meshes: the
                // mesh-to-mesh interpolated boundaries that the level solve replaces)
                bc_type[d] = bc_type[3 + d] = static_cast<int>(pb::BC::Periodic);
                continue;
            }
            if (lo == 0 || hi == 0) {
                if (n != 1) die("FDS pressure code 0 (periodic) on a non-periodic domain direction with more than one cell (inconsistent, S6)");
                bc_type[d] = bc_type[3 + d] = static_cast<int>(pb::BC::Neumann);   // one cell: the Neumann and periodic operators coincide
                continue;
            }
            if (lo > 4 || hi > 4) die("pressure boundary code 5/6 (cylindrical axis) is outside the S5 scope");
            bc_type[d] = static_cast<int>((lo == 1 || lo == 2) ? pb::BC::Dirichlet : pb::BC::Neumann);
            bc_type[3 + d] = static_cast<int>((hi == 1 || hi == 4) ? pb::BC::Dirichlet : pb::BC::Neumann);
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
        if (!iface_mode) return;
        const amrex::Box dom = l0.geom.Domain();
        each_local([&](int nm) {
            const amrex::Box b = l0.ba[nm - 1];
            int e[6] = {b.smallEnd(0) == dom.smallEnd(0), b.bigEnd(0) == dom.bigEnd(0), b.smallEnd(1) == dom.smallEnd(1), b.bigEnd(1) == dom.bigEnd(1), b.smallEnd(2) == dom.smallEnd(2), b.bigEnd(2) == dom.bigEnd(2)};
            fds_p_iface_walls(nm, on ? 1 : 0, e, (iface_mode == 2 || all) ? 1 : 0);
        });
    }

    void state(bool pred, bool first) { fds_k_state(pred ? 1 : 0, first ? 1 : 0, L.m_icyc, L.m_rmin, L.m_rmax); }
    template <class F> void each_local(F f) { for (int i = 0; i < nbox; ++i) if (local(i)) f(i + 1); }

    // ------------------------------------------------------------ D-031 level clip (port of the kernel check's clip_level; gather form)
    void clip_level(bool pred, std::vector<int>& flags)
    {
        Fields& F = *L.m_F;
        amrex::MultiFab& RP_ = F[pred ? "RHOS" : "RHO"];
        amrex::MultiFab& RZ = F[pred ? "ZZS" : "ZZ"];
        const amrex::Periodicity per = l0.geom.periodicity();
        RP_.FillBoundary(0, 1, RP_.nGrowVect(), per);
        RZ.FillBoundary(0, RZ.nComp(), RZ.nGrowVect(), per);
        const amrex::iMultiFab& mask = L.m_sd->mask();
        const int NS = RZ.nComp(), NT = nt;
        const double rmin = L.m_rmin, rmax = L.m_rmax;
        auto mkdx = [&](const amrex::Box& qb, std::vector<double> dxv[3]) { for (int d = 0; d < 3; ++d) dxv[d].assign(qb.length(d), l0.dx[d]); };
        flags.assign(nbox, 0);
        int fl[2] = {0, 0};
        std::vector<double> dxv[3];
        for (amrex::MFIter mfi(RP_); mfi.isValid(); ++mfi) {
            const amrex::Box qb = RP_[mfi].box(), mb = mask[mfi].box(), db = (*drho)[mfi].box(), v = mfi.validbox();
            mkdx(qb, dxv);
            int f[2], n;
            fds_clip_density(qb.loVect(), qb.hiVect(), mb.loVect(), mb.hiVect(), v.loVect(), v.hiVect(), db.loVect(), db.hiVect(), v.loVect(), v.hiVect(),
                             RP_[mfi].dataPtr(), mask[mfi].dataPtr(), dxv[0].data(), dxv[1].data(), dxv[2].data(), rmin, rmax, (*drho)[mfi].dataPtr(), f, &n);
            fl[0] |= f[0]; fl[1] |= f[1];
            flags[mfi.index()] = f[0] | (f[1] << 1);
        }
        bool rmn = fl[0], rmx = fl[1];
        amrex::ParallelAllReduce::Or(rmn, amrex::ParallelContext::CommunicatorSub());
        amrex::ParallelAllReduce::Or(rmx, amrex::ParallelContext::CommunicatorSub());
        if (rmn || rmx) {
            for (amrex::MFIter mfi(RP_); mfi.isValid(); ++mfi) {
                const amrex::Box qb = RP_[mfi].box(), db = (*drho)[mfi].box(), v = mfi.validbox();
                fds_clip_density_apply(qb.loVect(), qb.hiVect(), v.loVect(), v.hiVect(), db.loVect(), db.hiVect(), RP_[mfi].dataPtr(), (*drho)[mfi].dataPtr(), rmin, rmax);
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
                    int f, c;
                    fds_clip_species(fb.loVect(), fb.hiVect(), rb.loVect(), rb.hiVect(), mb.loVect(), mb.hiVect(), v.loVect(), v.hiVect(), NS, n, RP_[mfi].dataPtr(),
                                     RZ[mfi].dataPtr(), mask[mfi].dataPtr(), dxv[0].data(), dxv[1].data(), dxv[2].data(), (*dzz)[mfi].dataPtr(n - 1), &f, &c);
                    zfl[n - 1] |= f;
                }
            }
            for (int n = 0; n < NT; ++n) { bool b = zfl[n]; amrex::ParallelAllReduce::Or(b, amrex::ParallelContext::CommunicatorSub()); zfl[n] = b; anyz = anyz || b; }
            for (amrex::MFIter mfi(RP_); mfi.isValid(); ++mfi) {
                const amrex::Box fb = RZ[mfi].box(), rb = RP_[mfi].box(), v = mfi.validbox();
                for (int n = 1; n <= NT; ++n)
                    if (zfl[n - 1])
                        fds_clip_species_apply(rb.loVect(), rb.hiVect(), fb.loVect(), fb.hiVect(), v.loVect(), v.hiVect(), NS, n, RP_[mfi].dataPtr(), RZ[mfi].dataPtr(),
                                               (*dzz)[mfi].dataPtr(n - 1));
                if (rmn || rmx || anyz)
                    fds_clip_renorm(rb.loVect(), rb.hiVect(), fb.loVect(), fb.hiVect(), v.loVect(), v.hiVect(), NS, NT, RP_[mfi].dataPtr(), RZ[mfi].dataPtr(),
                                    mask[mfi].dataPtr());
            }
        }
    }

    // DENSITY of every box (the D-031 split form: pre-clip half, level gather clip, post-clip half)
    void density(bool pred, double t, double dt)
    {
        each_local([&](int nm) { fds_p_dens_pre(t, dt, nm); });
        std::vector<int> flags;
        clip_level(pred, flags);
        each_local([&](int nm) { fds_k_set_flags(nm, flags[nm - 1]); fds_k_dens_post(t, dt, nm); });
    }

    // ------------------------------------------------------------ D-028: zone integrals (DSUM, PSUM, USUM) as exact fixed-point sums
    void zone_sums(bool pred)
    {
        if (nzone <= 0) return;
        std::vector<int> g; std::vector<double> v;
        std::vector<double> fds_loc(3 * nzone, 0.0);
        {
            std::vector<double> d(nzone), p(nzone), u(nzone);
            fds_p_zone_get(nzone, d.data(), p.data(), u.data());
            for (int z = 0; z < nzone; ++z) { fds_loc[z] = d[z]; fds_loc[nzone + z] = p[z]; fds_loc[2 * nzone + z] = u[z]; }
        }
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
        amrex::ParallelAllReduce::Sum(fds_loc.data(), 3 * nzone, amrex::ParallelContext::CommunicatorSub());   // FDS's own order (per rank, then over ranks)
        for (int q = 0; q < 3 * nzone; ++q) zone_rel_step = std::max(zone_rel_step, std::abs(s[q] - fds_loc[q]) / std::max(std::abs(s[q]), 1e-30));
        if (std::getenv("FDSTL_ZONES")) for (int z = 0; z < nzone; ++z) amrex::Print() << "  [zone] icyc " << L.m_icyc << (pred ? " P " : " C ") << "DSUM " << std::setprecision(17) << s[z] << " PSUM " << s[nzone + z] << " USUM " << s[2 * nzone + z] << "\n";
        fds_p_zone_set(nzone, s.data(), s.data() + nzone, s.data() + 2 * nzone);
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
            if (fds_p_iter_baro() || iter == 1) {
                if (IP[4]) each_local([&](int nm) { fds_p_baroclinic(t, nm); });
                bc.exchange_om();   // MESH_EXCHANGE(5): FVX, FVY, FVZ (and H/HS) to OMESH
                each_local([&](int nm) { fds_g_match_flux(nm); });
            }
            each_local([&](int nm) { fds_p_noflux(dt, nm, iter == 1 ? 1 : 0); fds_p_rhs(t, dt, nm); });
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

    // ------------------------------------------------------------ one MAIN_LOOP iteration
    bool advance()
    {
        Fields& F = *L.m_F;
        BcStep& bc = *L.m_bc;
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
        each_local([&](int nm) { fds_k_visc(nm, 0); fds_p_mfd(nm); });
        first_pass = true;
        passes = 0;
        int stop = 0;
        for (;;) {
            ++passes;
            state(true, first_pass);
            density(true, t, dt);
            stage(passes == 1 ? "p1_dens" : "p2_dens", {"RHOS", "ZZS"});
            stage_raw(passes == 1 ? "p1_a_dens" : "p2_a_dens", {"RHOS", "TMP"});
            ghost_exchange(F, 1);
            stage_raw(passes == 1 ? "p1_b_fill" : "p2_b_fill", {"RHOS", "TMP"});
            bc.after_exchange(1, t, dt);
            stage_raw(passes == 1 ? "p1_c_visc" : "p2_c_visc", {"RHOS", "TMP"});
            stage_raw(passes == 1 ? "p1_prevflux" : "p2_prevflux", {"RHO", "RHOS", "U", "V", "W", "MU", "KRES", "H", "HS", "ZZ", "TMP"});
            iface(true);
            each_local([&](int nm) { fds_p_clear_attached(nm); fds_k_vflux(t, dt, nm, 0); });
            iface(false);
            stage(passes == 1 ? "p1_vflux" : "p2_vflux", {"FVX", "FVZ", "MU"});
            fds_p_init_div();
            bc.wall_bc(1, t, dt);
            if (std::getenv("FDSTL_WALLS") && L.m_icyc == std::atoi(std::getenv("FDSTL_WALLS")) && passes == 1)
                each_local([&](int nm) { const amrex::Box b = l0.ba[nm - 1]; fds_p_wall_dump(nm, b.smallEnd(0), b.smallEnd(2), 5); });
            stage_raw(passes == 1 ? "p1_prediv" : "p2_prediv", {"RHOS", "ZZS", "TMP", "RSUM", "U", "V", "W", "MU", "KRES", "D"});
            iface(true);
            each_local([&](int nm) { fds_k_div1(t, dt, nm); });
            iface(false);
            stage(passes == 1 ? "p1_div1" : "p2_div1", {"DS", "MU", "KRES", "TMP", "RSUM"});
            zone_sums(true);
            each_local([&](int nm) { fds_k_div2(dt, nm); });
            stage(passes == 1 ? "p1_div" : "p2_div", {"D", "DS", "DDDT"});
            state(true, first_pass);
            bool took = L.m_hook && L.m_hook(true, passes);
            if (!took) rec.it_pred = pressure_scheme(true, t, dt, rec.perr_pred, rec.verr_pred);
            // velocity predictor of every box, then the global DT decision (MIN over boxes and ranks)
            stage(passes == 1 ? "p1_press" : "p2_press", {"H", "FVX", "FVZ"});
            std::vector<double> dn(nbox, 0.0);
            std::vector<double> ci(nbox, 0.0);
            each_local([&](int nm) {
                double dtn, cfl, vn; int ic;
                fds_k_vpred(t + dt, dt, nm, &dtn, &ic, &cfl, &vn);
                dn[nm - 1] = dtn; ci[nm - 1] = ic;
            });
            amrex::ParallelAllReduce::Sum(dn.data(), nbox, amrex::ParallelContext::CommunicatorSub());
            amrex::ParallelAllReduce::Sum(ci.data(), nbox, amrex::ParallelContext::CommunicatorSub());
            for (int i = 0; i < nbox; ++i) { dt_new[i] = dn[i]; chg[i] = static_cast<int>(ci[i]); }
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
        ghost_exchange(F, 3);
        bc.after_exchange(3, t, dt);

        // ---------------- corrector
        t = time_advance(t, dt);
        state(false, first_pass);
        each_local([&](int nm) { fds_k_visc(nm, 1); fds_p_mfd(nm); });
        stage_raw("c_pre", {"RHOS", "TMP", "RHO"});
        density(false, t, dt);
        stage("c_dens", {"RHO", "ZZ"});
        ghost_exchange(F, 4);
        bc.after_exchange(4, t, dt);
        iface(true);
        each_local([&](int nm) { fds_p_clear_attached(nm); fds_k_vflux(t, dt, nm, 1); });
        iface(false);
        ++wall_counter;
        fds_p_set_wall_counter(wall_counter);
        bc.wall_bc(0, t, dt);
        if (wall_counter == IP[21]) wall_counter = 0;
        fds_p_set_wall_counter(wall_counter);
        fds_p_init_div();
        iface(true);
        each_local([&](int nm) { fds_k_div1(t, dt, nm); });
        iface(false);
        zone_sums(false);
        each_local([&](int nm) { fds_k_div2(dt, nm); });
        state(false, first_pass);
        {
            bool took = L.m_hook && L.m_hook(false, 1);
            if (!took) rec.it_corr = pressure_scheme(false, t, dt, rec.perr_corr, rec.verr_corr);
        }
        stage("c_press", {"HS", "DS"});
        each_local([&](int nm) { fds_k_vcorr(t, dt, nm); });
        stage("c_end", {"U", "W"});
        ghost_exchange(F, 6);
        bc.after_exchange(6, t, dt);
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
        if (!mms_done && (IP[6] == 7 || IP[6] == 11) && L.m_t >= RP[2]) { write_mms(L.m_t); mms_done = true; }
        if (!L.m_o.quiet)
            amrex::Print() << "step " << r.step << " T=" << L.m_t << " DT=" << r.dt << " passes=" << r.passes << " iters=" << r.it_pred << "/" << r.it_corr << " perr=" << r.perr_corr
                           << " verr=" << r.verr_corr << " zone_rel=" << r.zone_rel << "\n";
    }

    void finish()
    {
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

bool TimeLoop::advance() { return m->advance(); }

int TimeLoop::run()
{
    m->log_mass.push_back(m->mass_row(m_t));
    int fails = 0;
    for (;;) {
        if (m_o.max_steps >= 0 && m_icyc >= m_o.max_steps) break;
        if (!m->advance()) { ++fails; amrex::Print() << "time loop: STOP (instability or non-finite DT) at step " << m_icyc << "\n"; break; }
        if (m_t >= m->tc.t_end && m_icyc > 0) break;
    }
    m->finish();
    return fails;
}

}  // namespace fdsamr
