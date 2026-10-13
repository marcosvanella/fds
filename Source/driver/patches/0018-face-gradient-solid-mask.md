# 0018 face_gradient_composite with a Solid-only mask (Source/pressure_backend/CompositeSolve.cpp)
Status: **apply** (not committed; Role 2 file).
Without it every composite step with OBST aborts ("face_gradient_composite with masked cells is not built"): TimeLoop::solve_composite and
PressureBackendSolver (post-regrid projection) both call it once cell_class is set. Fix: the gradient of the unmasked operator, then zero on
every face touching a Solid cell (beta = 0 there, the masked flux). C/F faces never touch a Solid cell (the solve refuses that), so they are
unchanged; zeroing happens before the D-032 average-down. Known cells stay NotBuilt.
