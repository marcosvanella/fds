# R2b tests and sources of the interface flux overwrite (Role 3). include()d from CMakeLists.txt inside the AMReX + BUILD_TESTING block.
# FluxStageRunner / GhostShare are added to the AMReX library here so that this file is the only place that names them.
target_sources(fds_regrid_transport_amrex PRIVATE ${CMAKE_CURRENT_LIST_DIR}/FluxStageRunner.cpp ${CMAKE_CURRENT_LIST_DIR}/GhostShare.cpp)

# Interface flux overwrite end to end on the mock transport (override lists, finest-first order, conservation at round-off, negative controls, D-059 counts), 1 and 4 ranks
add_executable(test_flux_stage tests/test_flux_stage.cpp)
target_link_libraries(test_flux_stage PRIVATE fds_regrid_transport_amrex)
target_include_directories(test_flux_stage PRIVATE tests ${DRV}/tests)
add_test(NAME regrid_transport_flux_stage COMMAND test_flux_stage)
add_test(NAME regrid_transport_flux_stage_np4 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_flux_stage>)
set_tests_properties(regrid_transport_flux_stage regrid_transport_flux_stage_np4 PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")

# Driver-level checks need a driver build (RT_DRIVER_BUILD = directory with fds_amr). Without it they SKIP (exit code 77) with a message.
set(RT_DRIVER_BUILD "" CACHE PATH "Build directory of the AMR driver (contains fds_amr); enables the driver-level regrid_transport tests")
add_test(NAME regrid_transport_single_level_equiv COMMAND bash ${CMAKE_CURRENT_LIST_DIR}/tests/run_single_level_equiv.sh "${RT_DRIVER_BUILD}")
add_test(NAME regrid_transport_e2e_driver COMMAND bash ${CMAKE_CURRENT_LIST_DIR}/tests/run_e2e_driver.sh "${RT_DRIVER_BUILD}")
set_tests_properties(regrid_transport_single_level_equiv regrid_transport_e2e_driver PROPERTIES SKIP_RETURN_CODE 77)
