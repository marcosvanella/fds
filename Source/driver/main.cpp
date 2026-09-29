// main.cpp: C++ entry point of the FDS-AMReX driver (M2a, S1 skeleton).
//
// Order: MPI_Init_thread(FUNNELED) -> amrex::Initialize -> fds_setup(0) [unchanged FDS set-up, main.f90 patch 0001]
//        -> build the level-0 layout -> (S2 onward: fields, time loop) -> amrex::Finalize -> fds_setup(mode) [FDS end, calls
//        MPI_Finalize and STOP, so it must come last].
//
// Kernel-facing rules (M2a): (a) passive scalars are handled by Fields.cpp (S2); (b) only uniform Cartesian metrics are used.
#include <AMReX.H>
#include <AMReX_Print.H>

#include <mpi.h>

#include <cstdio>
#include <cstdlib>

#include "FdsAmr.H"

extern "C" void fds_setup(int mode, const char* fname, double* dt_out);

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::fprintf(stderr, "usage: fds_amr case.fds\n");
        return 1;
    }

    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

    // Exit code 0 only through the Fortran STOP in fds_setup(mode>0).
    {
        // build_parm_parse = false: argv[1] is the FDS input file, not an AMReX inputs file.
        amrex::Initialize(argc, argv, false, MPI_COMM_WORLD);
        int mode_end = 2;  // 2 = FDS setup-only stop message (S1 has no time loop yet)
        {
            double dt = 0.;
            fds_setup(0, argv[1], &dt);

            fdsamr::Level0 l0 = fdsamr::build_level0();
            fdsamr::print_level0(l0);
            amrex::Print() << "FDS-AMReX: initial dt from FDS set-up = " << dt << "\n";
        }
        amrex::Finalize();
        fds_setup(mode_end, argv[1], nullptr);  // does not return
    }
    return 0;
}
