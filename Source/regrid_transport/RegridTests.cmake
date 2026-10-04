# Tests of the regrid data path (Role 3, R4 part 2). Included from CMakeLists.txt inside the AMReX test block; uses DRV, MPIEXEC_*.
# R4 D-060: the single face-velocity prolongation path (parent divergence kept, interface faces = coarse value, max|div u - D| report), 1 and 4 ranks
add_executable(test_facetransfer tests/test_facetransfer.cpp)
target_link_libraries(test_facetransfer PRIVATE fds_regrid_transport_amrex)
target_include_directories(test_facetransfer PRIVATE tests ${DRV}/tests)
add_test(NAME regrid_transport_facetransfer COMMAND test_facetransfer)
add_test(NAME regrid_transport_facetransfer_np4 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_facetransfer>)

set_tests_properties(regrid_transport_facetransfer regrid_transport_facetransfer_np4 PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
