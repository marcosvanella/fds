# Tests of the regrid data path (Role 3, R4 part 2). Included from CMakeLists.txt inside the AMReX test block; uses DRV, MPIEXEC_*.
# R4 D-060: the single face-velocity prolongation path (parent divergence kept, interface faces = coarse value, max|div u - D| report), 1 and 4 ranks
add_executable(test_facetransfer tests/test_facetransfer.cpp)
target_link_libraries(test_facetransfer PRIVATE fds_regrid_transport_amrex)
target_include_directories(test_facetransfer PRIVATE tests ${DRV}/tests)
add_test(NAME regrid_transport_facetransfer COMMAND test_facetransfer)
add_test(NAME regrid_transport_facetransfer_np4 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_facetransfer>)

set_tests_properties(regrid_transport_facetransfer regrid_transport_facetransfer_np4 PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")

# Restriction of species is mass weighted (rho and rho*Z), negative control with a linear Z average, agreement with the driver's RegistryTransfer (Role 1)
add_executable(test_species_avgdown tests/test_species_avgdown.cpp ${DRV}/RegistryTransfer.cpp)
target_link_libraries(test_species_avgdown PRIVATE fds_regrid_transport_amrex)
target_include_directories(test_species_avgdown PRIVATE tests ${DRV}/tests)
add_test(NAME regrid_transport_species_avgdown COMMAND test_species_avgdown)
add_test(NAME regrid_transport_species_avgdown_np4 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_species_avgdown>)
set_tests_properties(regrid_transport_species_avgdown regrid_transport_species_avgdown_np4 PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")

# Moving blob end to end: RegridAmrCore + LevelRegistry + the driver's RegistryTransfer, prescribed velocity, conservation across regrids, uniform-fine comparison, controls
add_executable(test_blob_registry tests/test_blob_registry.cpp ${DRV}/RegistryTransfer.cpp)
target_link_libraries(test_blob_registry PRIVATE fds_regrid_transport_amrex)
target_include_directories(test_blob_registry PRIVATE tests ${DRV}/tests)
if (RT_HAVE_PRESSURE_BACKEND)
    target_link_libraries(test_blob_registry PRIVATE fds_rt_pressure_solver)
    target_compile_definitions(test_blob_registry PRIVATE RT_HAVE_PRESSURE_BACKEND)
endif()
add_test(NAME regrid_transport_blob_registry COMMAND test_blob_registry)
add_test(NAME regrid_transport_blob_registry_np4 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_blob_registry>)
add_test(NAME regrid_transport_blob_registry_ranks COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_regrid_rank_check.sh $<TARGET_FILE:test_blob_registry> ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG})
# Phase 3 moving-blob rows (V&V plan 5.12: B01 to B05, C07, F03, F04) on the registry with the mock conservative transport (the real FDS stages need a driver build: tests/run_e2e_driver.sh)
add_test(NAME regrid_transport_moving_blob_p3 COMMAND test_blob_registry --gates)
add_test(NAME regrid_transport_moving_blob_p3_np4 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_blob_registry> --gates)
set_tests_properties(regrid_transport_moving_blob_p3 regrid_transport_moving_blob_p3_np4 PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1" LABELS "phase3_gate_mock_transport")
set_tests_properties(regrid_transport_blob_registry regrid_transport_blob_registry_np4 regrid_transport_blob_registry_ranks PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
# The post-regrid projection part of these tests runs on a MOCK composite solver (amrex::MLMG): not a Phase 3 gate until Role 2's composite solver replaces it (notes/r4-part2-results.md)
set_tests_properties(regrid_transport_blob_registry regrid_transport_blob_registry_np4 regrid_transport_blob_registry_ranks PROPERTIES LABELS "projection_real_solver")

# R4 GPU path of the tagging kernels (K2 Fortran OpenMP target): the same kernel source built twice, host default and offload source (RT_OFFLOAD),
# each against an independent reference, and their outputs compared bitwise. Without an accelerator the offload source runs its target regions on
# the host (this checks the directive path and the has_device_addr/is_device_ptr clauses, not device memory). On a GPU machine configure with
# -DRT_OFFLOAD_FLAGS="-mp=gpu;-gpu=mem:managed,nofma" (nvfortran; no fast-math flags, BF-02 and BF-03 of docs/tools) to run on the device.
if (CMAKE_Fortran_COMPILER AND OpenMP_Fortran_FOUND)
    set(RT_OFFLOAD_FLAGS "" CACHE STRING "extra Fortran compile+link flags of the offload build of the tagging kernels (default: OpenMP only, host fallback)")
    set(RT_TK_SRC ${CMAKE_CURRENT_SOURCE_DIR}/rt_tag_kernels.F90 ${CMAKE_CURRENT_SOURCE_DIR}/tests/tag_kernel_check.F90)
    set_source_files_properties(${RT_TK_SRC} PROPERTIES Fortran_PREPROCESS ON)
    add_executable(tag_kernel_check_host ${RT_TK_SRC})
    set_target_properties(tag_kernel_check_host PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_tk_host)
    target_link_libraries(tag_kernel_check_host PRIVATE OpenMP::OpenMP_Fortran)
    add_executable(tag_kernel_check_offload ${RT_TK_SRC})
    set_target_properties(tag_kernel_check_offload PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_tk_off)
    target_compile_definitions(tag_kernel_check_offload PRIVATE RT_OFFLOAD)
    target_compile_options(tag_kernel_check_offload PRIVATE ${RT_OFFLOAD_FLAGS})
    target_link_options(tag_kernel_check_offload PRIVATE ${RT_OFFLOAD_FLAGS})
    target_link_libraries(tag_kernel_check_offload PRIVATE OpenMP::OpenMP_Fortran)
    add_test(NAME regrid_transport_tagkernel_host COMMAND tag_kernel_check_host)
    add_test(NAME regrid_transport_tagkernel_offload COMMAND tag_kernel_check_offload)
    add_test(NAME regrid_transport_tagkernel_host_vs_offload
             COMMAND ${CMAKE_COMMAND} -DHOST=$<TARGET_FILE:tag_kernel_check_host> -DOFFLOAD=$<TARGET_FILE:tag_kernel_check_offload>
                     -DWORKDIR=${CMAKE_CURRENT_BINARY_DIR} -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/compare_tag_kernel_outputs.cmake)
    set_tests_properties(regrid_transport_tagkernel_host regrid_transport_tagkernel_offload regrid_transport_tagkernel_host_vs_offload
                         PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=2" PASS_REGULAR_EXPRESSION "TAGKERNEL PASS|bitwise identical")
else()
    message(STATUS "regrid_transport: no Fortran/OpenMP, tagging-kernel host/offload check not built")
    add_test(NAME regrid_transport_tagkernel_host_vs_offload COMMAND ${CMAKE_COMMAND} -E echo "SKIP: no Fortran OpenMP compiler")
    set_tests_properties(regrid_transport_tagkernel_host_vs_offload PROPERTIES SKIP_REGULAR_EXPRESSION "SKIP:")
endif()
