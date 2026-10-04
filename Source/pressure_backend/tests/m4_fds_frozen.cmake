# CTest registrations for the M4 FDS-derived frozen cases (included from harness/CMakeLists.txt; no other file needs to change).
# pb_fds_frozen (harness/m4_modes.cpp) rebuilds the problem from the FDS dump metadata, solves the dumped PRHS with FFT, MLMG and
# HYPRE and compares H with FDS's own H (volume-weighted mean removed on both sides) and between backends, eps_H = 1e-8.
# Data and provenance: frozen/fds_csmag32_periodic/README.md.
find_program(PB_MPIEXEC NAMES mpiexec mpirun REQUIRED)
set(PB_FDS_FROZEN_DIR ${CMAKE_CURRENT_LIST_DIR}/../frozen/fds_csmag32_periodic)
foreach(stage P C)
  foreach(np 1 2)
    add_test(NAME pb_fds_frozen_csmag32_periodic_${stage}_np${np}
      COMMAND ${PB_MPIEXEC} --oversubscribe --bind-to none -np ${np} $<TARGET_FILE:pb_fds_frozen>
              mode=fds_frozen prefix=${PB_FDS_FROZEN_DIR}/fds_csmag32_pdump_n000003_${stage}_m001_ max_grid_size=16 eps_H=1e-8)
    set_tests_properties(pb_fds_frozen_csmag32_periodic_${stage}_np${np} PROPERTIES
      ENVIRONMENT "OMP_NUM_THREADS=1;OMP_DYNAMIC=false" TIMEOUT 600
      FAIL_REGULAR_EXPRESSION "CMP [^\n]* FAIL;CHECK FAIL")
  endforeach()
endforeach()
