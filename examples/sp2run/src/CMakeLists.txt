cmake_minimum_required(VERSION 2.8)
project(main Fortran)

enable_language(C)
enable_language(CXX)

include(FindOpenMP)

include(FindPkgConfig)
pkg_check_modules(BML bml REQUIRED)

include_directories(${BML_INCLUDEDIR})

add_executable(main
  main.F90)
set_target_properties(main
  PROPERTIES
  COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
  LINK_FLAGS ${OpenMP_Fortran_FLAGS})

target_link_libraries(main ${MY_FLAGS})
