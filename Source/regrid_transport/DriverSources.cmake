# DriverSources.cmake: adds Role 3's driver-level end-to-end modes (DriverModes.cpp, `fds_amr <case> --rt-e2e <test>`) and the sources they need to the driver target `fds`.
# Included from Source/driver/CMakeLists.txt after the target's own sources (one line, see DriverHook.patch). Defines FDSRT_DRIVER_MODES for the driver's main.cpp.
set(_rt ${CMAKE_CURRENT_LIST_DIR})
set(_drv ${_rt}/../driver)
target_sources(fds PRIVATE
    ${_rt}/DriverModes.cpp
    ${_rt}/TimeLoopWiring.cpp
    ${_rt}/DriverAdapter.cpp
    ${_rt}/LevelOps.cpp
    ${_rt}/CellTransfer.cpp
    ${_rt}/FaceTransfer.cpp
    ${_rt}/FluxOverrideOps.cpp
    ${_rt}/FluxStageRunner.cpp
    ${_rt}/GhostShare.cpp
    # DriverAdapter.cpp also holds install_post_regrid_projection, which names the regrid core and the projection: they have to be linked even though these modes do not call them
    ${_rt}/PostRegridProjection.cpp
    ${_rt}/RegridAmrCore.cpp
    ${_rt}/AmrInput.cpp
    ${_rt}/Hierarchy.cpp
    ${_rt}/TagOps.cpp
    ${_rt}/rt_tag_kernels.F90
    ${_drv}/RegistryTransfer.cpp)
target_compile_definitions(fds PRIVATE FDSRT_DRIVER_MODES)
set_source_files_properties(${_rt}/rt_tag_kernels.F90 PROPERTIES Fortran_PREPROCESS ON)
