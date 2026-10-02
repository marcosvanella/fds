# Role 3 library target and tests (regrid_transport).
#
# Stand-alone use (no FDS, no AMReX needed for the parser tests):
#   cmake -S Source/regrid_transport -B build-rt && cmake --build build-rt && ctest --test-dir build-rt
# Inside the FDS build the top-level file include()s this directory later (patch list in notes/); not yet included there.
if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    cmake_minimum_required(VERSION 3.16)
    project(fds_regrid_transport_standalone LANGUAGES CXX)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    include(CTest)
    # the installed AMReX was configured with Fortran interfaces: its package file needs the language enabled
    include(CheckLanguage)
    check_language(Fortran)
    if (CMAKE_Fortran_COMPILER)
        enable_language(Fortran)
    endif()
endif()

# AMReX-independent part: &AMR parser (R0), IR-002 grouping, input checks, static hierarchy (R1).
add_library(fds_regrid_transport_core STATIC AmrInput.cpp Hierarchy.cpp)
target_include_directories(fds_regrid_transport_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(fds_regrid_transport_core PRIVATE -Wall -Wextra)

# AMReX part (R2a): AmrCore subclass, level operators, glue to the driver's LevelRegistry. Built when AMReX is found (same install as Source/driver).
set(AMREX_INSTALL_DIR "/workspace/amrex-install" CACHE PATH "AMReX install prefix for the regrid_transport AMReX part")
if (NOT TARGET AMReX::amrex)
    find_package(AMReX 26.09 CONFIG QUIET PATHS ${AMREX_INSTALL_DIR} NO_DEFAULT_PATH)
    find_package(MPI QUIET COMPONENTS CXX)
endif()
if (TARGET AMReX::amrex AND TARGET MPI::MPI_CXX)
    set(DRV ${CMAKE_CURRENT_SOURCE_DIR}/../driver)   # read-only use of the driver's per-level objects
    # Driver objects the tests need (Fields, SideData, LevelRegistry, ...): AMReX only, no FDS Fortran.
    add_library(fds_rt_driver_objs STATIC ${DRV}/FdsAmr.cpp ${DRV}/Fields.cpp ${DRV}/SideData.cpp ${DRV}/ExactSum.cpp ${DRV}/LevelRegistry.cpp ${DRV}/GhostExchange.cpp)
    target_include_directories(fds_rt_driver_objs PUBLIC ${DRV})
    target_link_libraries(fds_rt_driver_objs PUBLIC AMReX::amrex MPI::MPI_CXX fds_regrid_transport_core)
    add_library(fds_regrid_transport_amrex STATIC RegridAmrCore.cpp LevelOps.cpp FluxOverrideOps.cpp DriverAdapter.cpp FaceTransfer.cpp TagOps.cpp CellTransfer.cpp rt_tag_kernels.F90)
    set_source_files_properties(rt_tag_kernels.F90 PROPERTIES Fortran_PREPROCESS ON)   # tagging kernels: Fortran, K2 style (ADR-001); the macros of the directives need the preprocessor
    target_include_directories(fds_regrid_transport_amrex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(fds_regrid_transport_amrex PUBLIC fds_regrid_transport_core fds_rt_driver_objs AMReX::amrex MPI::MPI_CXX)
    target_compile_options(fds_regrid_transport_amrex PRIVATE -Wall -Wextra)
    find_package(OpenMP QUIET COMPONENTS CXX)
    if (OpenMP_CXX_FOUND)
        target_link_libraries(fds_regrid_transport_amrex PUBLIC OpenMP::OpenMP_CXX)
    endif()
    find_package(OpenMP QUIET COMPONENTS Fortran)
    if (OpenMP_Fortran_FOUND)
        target_link_libraries(fds_regrid_transport_amrex PUBLIC OpenMP::OpenMP_Fortran)
    endif()
    # Compile check of the wiring into the driver's TimeLoop (not linked: TimeLoop needs the FDS Fortran objects)
    add_library(fds_rt_wiring_check OBJECT TimeLoopWiring.cpp)
    target_link_libraries(fds_rt_wiring_check PRIVATE fds_regrid_transport_amrex)
    set(RT_HAVE_AMREX ON)
endif()

# Interface target kept from the stub; consumers link this.
add_library(fds_regrid_transport INTERFACE)
target_include_directories(fds_regrid_transport INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(fds_regrid_transport INTERFACE fds_regrid_transport_core)

if (BUILD_TESTING)
    add_executable(test_amr_input tests/test_amr_input.cpp)
    target_link_libraries(test_amr_input PRIVATE fds_regrid_transport_core)
    target_compile_options(test_amr_input PRIVATE -Wall -Wextra)
    add_test(NAME regrid_transport_amr_input COMMAND test_amr_input)

    add_executable(test_hierarchy tests/test_hierarchy.cpp)
    target_link_libraries(test_hierarchy PRIVATE fds_regrid_transport_core)
    target_include_directories(test_hierarchy PRIVATE tests)
    target_compile_options(test_hierarchy PRIVATE -Wall -Wextra)
    add_test(NAME regrid_transport_hierarchy COMMAND test_hierarchy ${CMAKE_CURRENT_SOURCE_DIR}/tests/cases)

    if (RT_HAVE_AMREX)
        add_executable(test_amrcore tests/test_amrcore.cpp tests/fds_stubs.cpp)
        target_link_libraries(test_amrcore PRIVATE fds_regrid_transport_amrex)
        target_include_directories(test_amrcore PRIVATE tests ${DRV}/tests)
        target_compile_options(test_amrcore PRIVATE -Wall -Wextra)
        add_test(NAME regrid_transport_amrcore COMMAND test_amrcore ${CMAKE_CURRENT_SOURCE_DIR}/tests/cases)

        # R4: what AMReX FaceDivFree can and cannot do (notes/r4-facedivfree-findings.md); the ratio-1 case must abort inside AMReX (WILL_FAIL)
        add_executable(probe_facedivfree tests/probe_facedivfree.cpp)
        target_link_libraries(probe_facedivfree PRIVATE fds_regrid_transport_amrex)
        target_include_directories(probe_facedivfree PRIVATE tests ${DRV}/tests)
        add_test(NAME regrid_transport_facedivfree_probe COMMAND probe_facedivfree)
        add_test(NAME regrid_transport_facedivfree_ratio1_aborts COMMAND probe_facedivfree --ratio-1-hidden)
        set_tests_properties(regrid_transport_facedivfree_ratio1_aborts PROPERTIES WILL_FAIL TRUE)

        # R4: tagging kernels against an independent cell-by-cell evaluation (FR-011), run on 1 and 4 ranks
        add_executable(test_tagging tests/test_tagging.cpp)
        target_link_libraries(test_tagging PRIVATE fds_regrid_transport_amrex)
        target_include_directories(test_tagging PRIVATE tests ${DRV}/tests)
        add_test(NAME regrid_transport_tagging COMMAND test_tagging)
        add_test(NAME regrid_transport_tagging_np4 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_tagging>)

        # R4: conservative limited prolongation of rho*Z, parent-sum clip, regrid budget (FR-012), run on 1 and 4 ranks
        add_executable(test_celltransfer tests/test_celltransfer.cpp)
        target_link_libraries(test_celltransfer PRIVATE fds_regrid_transport_amrex)
        target_include_directories(test_celltransfer PRIVATE tests ${DRV}/tests)
        add_test(NAME regrid_transport_celltransfer COMMAND test_celltransfer)
        add_test(NAME regrid_transport_celltransfer_np4 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_celltransfer>)

        # FR-016 ghost-value comparison against the dumps of the instrumented FDS (A-09b): needs the dump files, so it is a test only when RT_FR016_DUMP is given
        add_executable(fr016_ghost_check tests/fr016_ghost_check.cpp tests/fds_stubs.cpp)
        target_link_libraries(fr016_ghost_check PRIVATE fds_regrid_transport_amrex)
        target_include_directories(fr016_ghost_check PRIVATE tests)
        if (RT_FR016_DUMP)
            add_test(NAME regrid_transport_fr016_ghosts COMMAND fr016_ghost_check ${CMAKE_CURRENT_SOURCE_DIR}/tests/cases ${RT_FR016_DUMP})
        endif()
    endif()
endif()
