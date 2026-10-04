// DriverModes.H: driver-level end-to-end modes of Role 3 (R2b), run INSIDE fds_amr on the real TimeLoop. Entered from the driver's main() as
//   fds_amr <case.fds> --rt-e2e <test> [options]        (hook: DriverHook.patch; sources: DriverSources.cmake)
// Nothing here is part of a normal run; a driver without the hook has no such mode and tests/run_e2e_driver.sh then skips (exit 77).
// The fine level runs the FDS kernels through TimeLoop::bind_level (patches 0007-0009, validated; run with gfortran and oneAPI).
#ifndef FDSRT_DRIVERMODES_H_
#define FDSRT_DRIVERMODES_H_

#include "FdsAmr.H"
#include "TimeLoop.H"

namespace fdsrt {

// Returns the number of failed checks (0 = pass). argv[2] is "--rt-e2e", argv[3] the test name ("transport", "ghost", "list"), the rest options (see DriverModes.cpp).
int driver_mode(int argc, char** argv, fdsamr::TimeLoop& loop, const fdsamr::Level0& l0, double dt_setup);

}  // namespace fdsrt

#endif
