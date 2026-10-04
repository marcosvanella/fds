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
add_library(fds_regrid_transport_core STATIC AmrInput.cpp Hierarchy.cpp InputConverter.cpp)
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
    add_library(fds_regrid_transport_amrex STATIC RegridAmrCore.cpp LevelOps.cpp FluxOverrideOps.cpp DriverAdapter.cpp FaceTransfer.cpp TagOps.cpp CellTransfer.cpp PostRegridProjection.cpp rt_tag_kernels.F90)
    set_source_files_properties(rt_tag_kernels.F90 PROPERTIES Fortran_PREPROCESS ON)   # tagging kernels: Fortran, K2 style (ADR-001); the macros of the directives need the preprocessor
    target_include_directories(fds_regrid_transport_amrex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
    # Build-time choice of the tagging kernels (same entry points either way): host loops (default) or the OpenMP target offload variant of rt_tag_kernels.F90.
    # Device use needs AMReX built with a GPU backend (TagBox / FArrayBox data are then device pointers); RT_TAG_OFFLOAD_FLAGS carries the compiler flags (nvfortran: -mp=gpu;-gpu=mem:managed,nofma).
    option(RT_TAG_OFFLOAD "build the tagging kernels with the OpenMP target offload variant (RT_OFFLOAD)" OFF)
    set(RT_TAG_OFFLOAD_FLAGS "" CACHE STRING "extra Fortran compile and link flags of the offload variant of the tagging kernels")
    if (RT_TAG_OFFLOAD)
        set_source_files_properties(rt_tag_kernels.F90 PROPERTIES COMPILE_DEFINITIONS RT_OFFLOAD COMPILE_OPTIONS "${RT_TAG_OFFLOAD_FLAGS}")
        if (RT_TAG_OFFLOAD_FLAGS)
            target_link_options(fds_regrid_transport_amrex PUBLIC ${RT_TAG_OFFLOAD_FLAGS})
        endif()
        if (NOT AMReX_GPU_BACKEND STREQUAL "CUDA" AND NOT AMReX_GPU_BACKEND STREQUAL "HIP" AND NOT AMReX_GPU_BACKEND STREQUAL "SYCL")
            message(WARNING "regrid_transport: RT_TAG_OFFLOAD with an AMReX that has no GPU backend: the target regions run on the host or on managed memory (bring-up only); the tag arrays are not device pointers")
        endif()
        message(STATUS "regrid_transport: tagging kernels built as the OpenMP target offload variant")
    endif()
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
    # Binding of the post-regrid projection to Role 2's composite MLMG (Source/pressure_backend; needs AMReX with FFT and linear solvers, as that directory does)
    set(PB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../pressure_backend)
    if (EXISTS ${PB_DIR}/CompositeSolve.cpp AND AMReX_FFT AND AMReX_LSOLVERS)
        if (NOT TARGET pressure_backend)
            add_library(pressure_backend STATIC ${PB_DIR}/CommonLayer.cpp ${PB_DIR}/PressureIface.cpp ${PB_DIR}/FFTBackend.cpp ${PB_DIR}/MLMGBackend.cpp ${PB_DIR}/CompositeSolve.cpp)
            # Role 2's HYPRE backend: needs an AMReX built with HYPRE; otherwise the driver's link stub (Source/driver/HypreStub.cpp) provides the symbols when it exists
            if (EXISTS ${PB_DIR}/HypreBackend.cpp)
                if (AMReX_HYPRE_FOUND)
                    target_sources(pressure_backend PRIVATE ${PB_DIR}/HypreBackend.cpp)
                elseif (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../driver/HypreStub.cpp)
                    target_sources(pressure_backend PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../driver/HypreStub.cpp)
                else()
                    message(WARNING "regrid_transport: the pressure backend needs HYPRE symbols and neither AMReX_HYPRE_FOUND nor a stub exists; link of the real-solver tests will fail")
                endif()
            endif()
            target_include_directories(pressure_backend PRIVATE ${PB_DIR})   # private: its ExactSum.H must not shadow the driver's header of the same name
            target_link_libraries(pressure_backend PUBLIC AMReX::amrex)
        endif()
        add_library(fds_rt_pressure_solver STATIC PressureBackendSolver.cpp)
        target_include_directories(fds_rt_pressure_solver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} PRIVATE ${PB_DIR})
        target_link_libraries(fds_rt_pressure_solver PUBLIC fds_regrid_transport_amrex pressure_backend)
        target_compile_options(fds_rt_pressure_solver PRIVATE -Wall -Wextra)
        set(RT_HAVE_PRESSURE_BACKEND ON)
    else()
        message(STATUS "regrid_transport: pressure_backend or AMReX FFT/linear solvers not available; the projection tests use the AMReX-MLMG mock only")
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

    # D-076 input converter: classification, level-0-only input text, hierarchy, round trip, error cases (pure C++), and its command line tool
    add_executable(fds_amr_convert_input ConvertInputTool.cpp)
    target_link_libraries(fds_amr_convert_input PRIVATE fds_regrid_transport_core)
    target_compile_options(fds_amr_convert_input PRIVATE -Wall -Wextra)
    add_executable(test_input_converter tests/test_input_converter.cpp)
    target_link_libraries(test_input_converter PRIVATE fds_regrid_transport_core)
    target_include_directories(test_input_converter PRIVATE tests)
    target_compile_options(test_input_converter PRIVATE -Wall -Wextra)
    add_test(NAME regrid_transport_input_converter COMMAND test_input_converter ${CMAKE_CURRENT_SOURCE_DIR}/tests/cases)
    # the tool on the FR-016 gate mesh lines: convert, then convert its own output (a level-0-only input converts to itself)
    add_test(NAME regrid_transport_convert_input_tool
             COMMAND sh -c "$<TARGET_FILE:fds_amr_convert_input> ${CMAKE_CURRENT_SOURCE_DIR}/tests/cases/ns2d_16_int_1to2.fds -o ${CMAKE_CURRENT_BINARY_DIR}/ns2d_level0.fds && $<TARGET_FILE:fds_amr_convert_input> ${CMAKE_CURRENT_BINARY_DIR}/ns2d_level0.fds -o ${CMAKE_CURRENT_BINARY_DIR}/ns2d_level0_again.fds && cmp ${CMAKE_CURRENT_BINARY_DIR}/ns2d_level0.fds ${CMAKE_CURRENT_BINARY_DIR}/ns2d_level0_again.fds")
    # negative control: without the &AMR line the same tool run must fail
    add_test(NAME regrid_transport_convert_input_tool_no_amr
             COMMAND sh -c "grep -v '^&AMR' ${CMAKE_CURRENT_SOURCE_DIR}/tests/cases/ns2d_16_int_1to2.fds > ${CMAKE_CURRENT_BINARY_DIR}/ns2d_noamr.fds; ! $<TARGET_FILE:fds_amr_convert_input> ${CMAKE_CURRENT_BINARY_DIR}/ns2d_noamr.fds")

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

        # R4: dynamic hierarchy through RegridAmrCore (tags, region clip, make/remake/clear, transfer, conservation across regrids), 1 and 4 ranks
        add_executable(test_regrid_core tests/test_regrid_core.cpp tests/fds_stubs.cpp)
        target_link_libraries(test_regrid_core PRIVATE fds_regrid_transport_amrex)
        target_include_directories(test_regrid_core PRIVATE tests ${DRV}/tests)
        add_test(NAME regrid_transport_regrid_core COMMAND test_regrid_core)
        add_test(NAME regrid_transport_regrid_core_np4 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_regrid_core>)
        add_executable(test_tag_buffer_large tests/test_tag_buffer_large.cpp tests/fds_stubs.cpp)
        target_link_libraries(test_tag_buffer_large PRIVATE fds_regrid_transport_amrex)
        target_include_directories(test_tag_buffer_large PRIVATE tests ${DRV}/tests)
        add_test(NAME regrid_transport_tag_buffer_large COMMAND test_tag_buffer_large)
        add_test(NAME regrid_transport_tag_buffer_large_np4 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_tag_buffer_large>)
        add_test(NAME regrid_transport_regrid_core_ranks COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_regrid_rank_check.sh $<TARGET_FILE:test_regrid_core> ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG})
        set_tests_properties(regrid_transport_regrid_core regrid_transport_regrid_core_np4 regrid_transport_regrid_core_ranks PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")   # one thread per rank: the test is small and oversubscribed threads only slow it

        # R4 part 2 tests (face prolongation, registry transfer, dynamic blob end to end): RegridTests.cmake
        include(${CMAKE_CURRENT_SOURCE_DIR}/RegridTests.cmake)

        # R2b: interface flux overwrite (runner, shared-ghost count, mock transport tests): FluxTests.cmake
        include(${CMAKE_CURRENT_SOURCE_DIR}/FluxTests.cmake)

        # FR-016 ghost-value comparison against the dumps of the instrumented FDS (A-09b): needs the dump files, so it is a test only when RT_FR016_DUMP is given
        add_executable(fr016_ghost_check tests/fr016_ghost_check.cpp tests/fds_stubs.cpp)
        target_link_libraries(fr016_ghost_check PRIVATE fds_regrid_transport_amrex)
        target_include_directories(fr016_ghost_check PRIVATE tests)
        if (RT_FR016_DUMP)
            add_test(NAME regrid_transport_fr016_ghosts COMMAND fr016_ghost_check ${CMAKE_CURRENT_SOURCE_DIR}/tests/cases ${RT_FR016_DUMP})
        endif()
    endif()
endif()
