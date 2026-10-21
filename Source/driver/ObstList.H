// ObstList.H: loads pb::obst_list() (pressure_backend/ObstMask.H) once from the FDS OBST export (fds_mesh_query.f90: fds_obst_count, fds_obst_xb),
// gathered over all ranks so every rank holds the global list. Collective on first call.
#pragma once
#include "../pressure_backend/ObstMask.H"
#include <AMReX_ParallelDescriptor.H>

extern "C" {
int fds_obst_count (); void fds_obst_xb (int n, double* xb); int fds_stateful_vents (int n, double* xb, int ex); int fds_bed_xb (int n, double* xb);
void fds_obst_surf (int n, int* s); int fds_vent_list (int n, double* xb, int* v);                                          // fds_mesh_query.f90
void fds_fine_set_lists (int no, const double* oxb, const int* osf, int nv, const double* vxb, const int* vi);          // fds_fine_level.f90
}

namespace fdsamr {

// local flat list -> global flat list on every rank
inline std::vector<double> gather_flat (std::vector<double> const& loc)
{
    const int np = amrex::ParallelDescriptor::NProcs();
    std::vector<int> cnt(np, 0), off(np, 0);
    int m = static_cast<int>(loc.size());
#ifdef AMREX_USE_MPI
    MPI_Allgather(&m, 1, MPI_INT, cnt.data(), 1, MPI_INT, amrex::ParallelDescriptor::Communicator());
#else
    cnt[0] = m;
#endif
    int tot = 0;
    for (int r = 0; r < np; ++r) { off[r] = tot; tot += cnt[r]; }
    std::vector<double> all(tot);
#ifdef AMREX_USE_MPI
    MPI_Allgatherv(loc.data(), m, MPI_DOUBLE, all.data(), cnt.data(), off.data(), MPI_DOUBLE, amrex::ParallelDescriptor::Communicator());
#else
    all = loc;
#endif
    return all;
}

// local XB list (6 per entry) -> global list on every rank
inline void gather_xb (std::vector<double> const& loc, std::vector<std::array<double, 6>>& L)
{
    const std::vector<double> all = gather_flat(loc);
    const int tot = static_cast<int>(all.size());
    L.clear();
    for (int i = 0; i + 5 < tot; i += 6) { L.push_back({all[i], all[i+1], all[i+2], all[i+3], all[i+4], all[i+5]}); }
}

// pb::nofine_list(): stateful VENT XB (D-010); ex = 1 (bed runs) drops the analytic burners from it
inline void load_nofine (int ex)
{
    std::vector<double> loc(6 * static_cast<size_t>(fds_stateful_vents(0, nullptr, ex)), 0.0);
    if (!loc.empty()) { fds_stateful_vents(static_cast<int>(loc.size() / 6), loc.data(), ex); }
    gather_xb(loc, pb::nofine_list());
}

// pb::obst_list() (OBST XB) and pb::nofine_list() (stateful VENT XB, D-010), once
inline void load_obst_list ()
{
    static bool done = false;
    if (done) { return; }
    done = true;
    const int n = fds_obst_count();
    std::vector<double> loc(6 * static_cast<size_t>(n));
    if (n > 0) { fds_obst_xb(n, loc.data()); }
    gather_xb(loc, pb::obst_list());
    load_nofine(0);
    const int nb = fds_bed_xb(0, nullptr);   // global INIT list: no gather
    loc.assign(6 * static_cast<size_t>(nb), 0.0);
    if (nb > 0) { fds_bed_xb(nb, loc.data()); }
    for (int i = 0; i < nb; ++i) { pb::bed_list().push_back({loc[6*i], loc[6*i+1], loc[6*i+2], loc[6*i+3], loc[6*i+4], loc[6*i+5]}); }
    // per-level fine walls (fds_fine_level.f90): OBST face surfaces and all VENTs, global lists, handed to Fortran once
    std::vector<double> osf(7 * static_cast<size_t>(n));
    if (n > 0) { std::vector<int> si(7 * static_cast<size_t>(n)); fds_obst_surf(n, si.data()); osf.assign(si.begin(), si.end()); }
    const std::vector<double> gosf = gather_flat(osf);
    const int nv = fds_vent_list(0, nullptr, nullptr);
    std::vector<double> vxb(6 * static_cast<size_t>(nv)), vi_d(3 * static_cast<size_t>(nv));
    if (nv > 0) { std::vector<int> vi(3 * static_cast<size_t>(nv)); fds_vent_list(nv, vxb.data(), vi.data()); vi_d.assign(vi.begin(), vi.end()); }
    const std::vector<double> gvxb = gather_flat(vxb), gvi = gather_flat(vi_d);
    const std::vector<int> iosf(gosf.begin(), gosf.end()), ivi(gvi.begin(), gvi.end());
    std::vector<double> oxb;
    for (auto const& b : pb::obst_list()) oxb.insert(oxb.end(), b.begin(), b.end());
    fds_fine_set_lists(static_cast<int>(pb::obst_list().size()), oxb.data(), iosf.data(), static_cast<int>(gvi.size() / 3), gvxb.data(), ivi.data());
}

}  // namespace fdsamr
