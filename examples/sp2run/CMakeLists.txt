cmake_minimum_required(VERSION 2.8)
project(main Fortran)

include(FindPkgConfig)
pkg_check_modules(BML bml REQUIRED)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}) 

set(MY_FLAGS " " CACHE STRING "Linking with this")

message("linking with this ${MY_FLAGS}")
add_subdirectory(src)
