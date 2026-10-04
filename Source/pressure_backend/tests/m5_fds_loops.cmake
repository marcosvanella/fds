# CTest registrations for the host translations of FDS pressure-area loops (FdsPressureLoops.cpp; included from harness/CMakeLists.txt).
# pb_fds_loops (harness/loops_modes.cpp) compares the C++ functions bitwise with the verbatim upstream loops (pres.f90 of FireX 36975d7)
# compiled by gfortran (tests/fds_loops/ref_loops.f90.in, tests/fds_loops_test.py). Needs only a Fortran compiler; no AMReX, no MPI.
# Mutation check: PB_FDSLOOPS_MUTANT builds change one operand, sign or index and the test must fail on each.
find_package(Python3 REQUIRED COMPONENTS Interpreter)
set(PB_LOOPS_PRES ${CMAKE_CURRENT_LIST_DIR}/../../pres.f90 CACHE FILEPATH "pres.f90 (FireX 36975d7 text) the verbatim loops are copied from")
set(PB_LOOPS_GITREPO "" CACHE PATH "optional git repository to try 'git show 36975d7:Source/pres.f90' in when PB_LOOPS_PRES drifted")
set(PB_LOOPS_WORK ${CMAKE_CURRENT_BINARY_DIR}/fds_loops_work)
set(PB_LOOPS_TEST ${CMAKE_CURRENT_LIST_DIR}/fds_loops_test.py)
set(PB_LOOPS_NMUT 18)
set(PB_LOOPS_GIT_ARGS "")
if(PB_LOOPS_GITREPO)
  set(PB_LOOPS_GIT_ARGS --git-repo ${PB_LOOPS_GITREPO})
endif()

function(pb_loops_exe name)
  add_executable(${name} ${CMAKE_CURRENT_LIST_DIR}/../harness/loops_modes.cpp ${CMAKE_CURRENT_LIST_DIR}/../FdsPressureLoops.cpp)
  target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/..)
  target_compile_options(${name} PRIVATE -ffp-contract=off)
  target_compile_features(${name} PRIVATE cxx_std_17)
  foreach(d ${ARGN})
    target_compile_definitions(${name} PRIVATE ${d})
  endforeach()
endfunction()

pb_loops_exe(pb_fds_loops)
foreach(m RANGE 1 ${PB_LOOPS_NMUT})
  pb_loops_exe(pb_fds_loops_mut${m} PB_FDSLOOPS_MUTANT=${m})
endforeach()

add_test(NAME pb_fdsloops_ref
  COMMAND ${Python3_EXECUTABLE} ${PB_LOOPS_TEST} make-ref --pres ${PB_LOOPS_PRES} --fortran ${CMAKE_Fortran_COMPILER} --work ${PB_LOOPS_WORK} ${PB_LOOPS_GIT_ARGS})
set_tests_properties(pb_fdsloops_ref PROPERTIES FIXTURES_SETUP fdsloops_ref TIMEOUT 300)

add_test(NAME pb_fdsloops_bitwise
  COMMAND ${Python3_EXECUTABLE} ${PB_LOOPS_TEST} check --work ${PB_LOOPS_WORK} --exe $<TARGET_FILE:pb_fds_loops>)
set_tests_properties(pb_fdsloops_bitwise PROPERTIES FIXTURES_REQUIRED fdsloops_ref TIMEOUT 300)

foreach(m RANGE 1 ${PB_LOOPS_NMUT})
  add_test(NAME pb_fdsloops_mutant_${m}
    COMMAND ${Python3_EXECUTABLE} ${PB_LOOPS_TEST} check --expect-fail --work ${PB_LOOPS_WORK} --exe $<TARGET_FILE:pb_fds_loops_mut${m}>)
  set_tests_properties(pb_fdsloops_mutant_${m} PROPERTIES FIXTURES_REQUIRED fdsloops_ref TIMEOUT 300)
endforeach()

add_test(NAME pb_fdsloops_drift
  COMMAND ${Python3_EXECUTABLE} ${PB_LOOPS_TEST} drift --pres ${PB_LOOPS_PRES} --work ${CMAKE_CURRENT_BINARY_DIR}/fds_loops_drift)

# Against arrays written by a real FDS run (frozen/fds_loops_cases, made with frozen/fds_loops_dump_hook.py on a scratch copy of pres.f90): L1211,
# L1207, L1209 (including the OPEN-boundary walls) and L1220-L1222 are bitwise equal to what FDS computed. The velocity arrays of FDS start at
# index -1 in their own direction (init.f90); pb_fdsloops_fds_uu_lb0 reads them with lower bound 0 and must fail (negative control for the
# former "open walls differ" finding, frozen/fds-loops-notes.md).
add_test(NAME pb_fdsloops_fds
  COMMAND ${Python3_EXECUTABLE} ${PB_LOOPS_TEST} check-fds --work ${PB_LOOPS_WORK} --exe $<TARGET_FILE:pb_fds_loops>
          --archive ${CMAKE_CURRENT_LIST_DIR}/../frozen/fds_loops_cases/fds_loops_cases.tar)
set_tests_properties(pb_fdsloops_fds PROPERTIES TIMEOUT 300)
add_test(NAME pb_fdsloops_fds_uu_lb0
  COMMAND ${Python3_EXECUTABLE} ${PB_LOOPS_TEST} check-fds --expect-fail --exe-arg uu_lb0=1 --work ${PB_LOOPS_WORK}_lb0 --exe $<TARGET_FILE:pb_fds_loops>
          --archive ${CMAKE_CURRENT_LIST_DIR}/../frozen/fds_loops_cases/fds_loops_cases.tar)
set_tests_properties(pb_fdsloops_fds_uu_lb0 PROPERTIES TIMEOUT 300)

# Guards and edge indices (harness/loops_guards.cpp): C1 whole-extent guard of the periodic H fill, C2 TUNNEL_PRECONDITIONER refusal,
# C3 L1209 wind arrays must cover 0..KBP1 and the z-wall indices KK=0 / KK=KBP1. Mutants 19 (no whole-extent guard), 20 (the former
# wind guard 1..KBAR) and 21 (no tunnel refusal) must fail it. The AddressSanitizer build proves that no read leaves an array; mutant 20
# under the sanitizer reproduces the former heap-buffer-overflow at the W_WIND(K) read (tests expects it to fail).
function(pb_loops_guards_exe name)
  add_executable(${name} ${CMAKE_CURRENT_LIST_DIR}/../harness/loops_guards.cpp ${CMAKE_CURRENT_LIST_DIR}/../FdsPressureLoops.cpp)
  target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/..)
  target_compile_options(${name} PRIVATE -ffp-contract=off)
  target_compile_features(${name} PRIVATE cxx_std_17)
  foreach(d ${ARGN})
    target_compile_definitions(${name} PRIVATE ${d})
  endforeach()
endfunction()
pb_loops_guards_exe(pb_fds_loops_guards)
add_test(NAME pb_fdsloops_guards COMMAND pb_fds_loops_guards)
foreach(m 19 20 21)
  pb_loops_guards_exe(pb_fds_loops_guards_mut${m} PB_FDSLOOPS_MUTANT=${m})
  add_test(NAME pb_fdsloops_guards_mutant_${m} COMMAND pb_fds_loops_guards_mut${m})
  set_tests_properties(pb_fdsloops_guards_mutant_${m} PROPERTIES WILL_FAIL TRUE)
endforeach()
include(CheckCXXSourceCompiles)
set(CMAKE_REQUIRED_FLAGS "-fsanitize=address")
set(CMAKE_REQUIRED_LINK_OPTIONS "-fsanitize=address")
check_cxx_source_compiles("int main() { return 0; }" PB_LOOPS_HAVE_ASAN)
unset(CMAKE_REQUIRED_FLAGS)
unset(CMAKE_REQUIRED_LINK_OPTIONS)
if(PB_LOOPS_HAVE_ASAN)
  pb_loops_guards_exe(pb_fds_loops_guards_asan)
  target_compile_options(pb_fds_loops_guards_asan PRIVATE -fsanitize=address -fno-omit-frame-pointer -g)
  target_link_options(pb_fds_loops_guards_asan PRIVATE -fsanitize=address)
  add_test(NAME pb_fdsloops_guards_asan COMMAND pb_fds_loops_guards_asan)
  pb_loops_guards_exe(pb_fds_loops_guards_asan_mut20 PB_FDSLOOPS_MUTANT=20)
  target_compile_options(pb_fds_loops_guards_asan_mut20 PRIVATE -fsanitize=address -fno-omit-frame-pointer -g)
  target_link_options(pb_fds_loops_guards_asan_mut20 PRIVATE -fsanitize=address)
  add_test(NAME pb_fdsloops_guards_asan_mutant_20 COMMAND pb_fds_loops_guards_asan_mut20)
  set_tests_properties(pb_fdsloops_guards_asan_mutant_20 PROPERTIES WILL_FAIL TRUE)
endif()

# R2 (AMReX low-face index mapping of PRHS) and R3 (D-067 order gauge -> H ghost fill -> P on a singular Neumann component), harness/loops_r2r3.cpp.
# Mutant 22 (wrong direction of the face-index shift) must fail it. Acceptance against FDS for the AMR route stays on eps_H, not bitwise.
function(pb_loops_r2r3_exe name)
  add_executable(${name} ${CMAKE_CURRENT_LIST_DIR}/../harness/loops_r2r3.cpp ${CMAKE_CURRENT_LIST_DIR}/../FdsPressureLoops.cpp)
  target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/..)
  target_compile_options(${name} PRIVATE -ffp-contract=off)
  target_compile_features(${name} PRIVATE cxx_std_17)
  foreach(d ${ARGN})
    target_compile_definitions(${name} PRIVATE ${d})
  endforeach()
endfunction()
pb_loops_r2r3_exe(pb_fds_loops_r2r3)
add_test(NAME pb_fdsloops_r2r3 COMMAND pb_fds_loops_r2r3)
pb_loops_r2r3_exe(pb_fds_loops_r2r3_mut22 PB_FDSLOOPS_MUTANT=22)
add_test(NAME pb_fdsloops_r2r3_mutant_22 COMMAND pb_fds_loops_r2r3_mut22)
set_tests_properties(pb_fdsloops_r2r3_mutant_22 PROPERTIES WILL_FAIL TRUE)
