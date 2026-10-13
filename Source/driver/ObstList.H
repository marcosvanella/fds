// ObstList.H: loads pb::obst_list() (pressure_backend/ObstMask.H) once from the FDS OBST export (fds_mesh_query.f90: fds_obst_count, fds_obst_xb),
// gathered over all ranks so every rank holds the global list. Collective on first call.
#pragma once
#include "../pressure_backend/ObstMask.H"
#include <AMReX_ParallelDescriptor.H>

extern "C" { int fds_obst_count (); void fds_obst_xb (int n, double* xb); }

namespace fdsamr {

inline void load_obst_list ()
{
    static bool done = false;
    if (done) { return; }
    done = true;
    const int n = fds_obst_count();
    std::vector<double> loc(6 * static_cast<size_t>(n));
    if (n > 0) { fds_obst_xb(n, loc.data()); }
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
    auto& L = pb::obst_list();
    L.clear();
    for (int i = 0; i + 5 < tot; i += 6) { L.push_back({all[i], all[i+1], all[i+2], all[i+3], all[i+4], all[i+5]}); }
}

}  // namespace fdsamr
