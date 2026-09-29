# CTest registrations (included from harness/CMakeLists.txt). All tests run the harness through pb_test.py.
find_package(Python3 REQUIRED COMPONENTS Interpreter)
find_program(PB_MPIEXEC NAMES mpiexec mpirun REQUIRED)
set(PB_ENV "OMP_NUM_THREADS=1;OMP_DYNAMIC=false")

function(pb_add_test name)
  add_test(NAME ${name}
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_LIST_DIR}/pb_test.py
            --harness $<TARGET_FILE:pb_harness> --mpiexec ${PB_MPIEXEC}
            --work ${CMAKE_CURRENT_BINARY_DIR}/work_${name} ${ARGN})
  set_tests_properties(${name} PROPERTIES ENVIRONMENT "${PB_ENV}" TIMEOUT 600)
endfunction()

pb_add_test(pb_selector      selector)
pb_add_test(pb_exactsum      exactsum)
foreach(bc neumann periodic)
  pb_add_test(pb_fft_vs_mlmg_${bc}_64     fftmlmg --bc ${bc} --n "64 64 64")
  pb_add_test(pb_fft_vs_mlmg_${bc}_34x18x32 fftmlmg --bc ${bc} --n "34 18 32")
  pb_add_test(pb_frozen_${bc}_34x18x32    frozen  --bc ${bc} --n "34 18 32")
  pb_add_test(pb_decomp_${bc}_64          decomp  --bc ${bc} --n "64 64 64")
  pb_add_test(pb_repeat_${bc}_64          repeat  --bc ${bc} --n "64 64 64")
endforeach()
pb_add_test(pb_singular_hook singular)
