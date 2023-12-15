cmake_minimum_required(VERSION 3.16)
project(appendix C CXX Fortran)

set(dir ${CMAKE_CURRENT_SOURCE_DIR}/build/)
set(CMAKE_BUILD_DIRECTORY ${dir})
#set(CMAKE_INSTALL_BINDIR ${dir})

include(FindPkgConfig)

find_package(PROGRESS CONFIG QUIET)
pkg_check_modules(PROGRESS REQUIRED progress)
message(STATUS "Found progress: ${PROGRESS_LDFLAGS}")
list(APPEND LINK_LIBRARIES ${PROGRESS_LDFLAGS})

find_package(BML CONFIG QUIET)
pkg_check_modules(BML REQUIRED bml)
list(APPEND LINK_LIBRARIES BML::bml)
list(APPEND LINK_LIBRARIES ${BML_LDFLAGS})
message(STATUS "Found bml: ${BML_LDFLAGS}")

if(PROGRESS_MPI)
  message(STATUS "Will build with MPI")
  add_definitions(-DDO_MPI)
endif()

message(STATUS "Extra FC Flags ${EXTRA_FCFLAGS}")

if(DEFINED EXTRA_FCFLAGS)
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTRA_FCFLAGS}")
endif()

list(APPEND LINK_LIBRARIES ${EXTRA_LINK_FLAGS})
find_library(FOUND_METIS metis)
list(APPEND LINK_LIBRARIES ${FOUND_METIS})

message(STATUS "Project sources = " ${PROJECT_SOURCE_DIR} )
include_directories(${PROJECT_SOURCE_DIR}/)
include_directories(${CMAKE_BINARY_DIR}/)
include_directories(${BML_INCLUDEDIR})
include_directories(${PROGRESS_INCLUDEDIR})

function(progress_appendix myappendix main_and_srcs)
list(GET main_and_srcs 0 main)
include_directories(${PROGRESS_INCLUDEDIR})
add_executable(${myappendix} ${main})
target_sources(${myappendix} PRIVATE ${ARGN})
target_link_libraries(${myappendix} PUBLIC
    ${LINK_LIBRARIES})
  set_target_properties(${myappendix}
    PROPERTIES
    LINK_FLAGS "--coverage")
endfunction(progress_appendix)


progress_appendix(part1 part1.F90
                            aux_mod.F90
                            )

progress_appendix(part2 part2.F90
                            aux_mod.F90
                            )

progress_appendix(dmconstruction_graphBased dmconstruction_graphBased.F90
                            aux_mod.F90
                            )
install(TARGETS part1 part2 DESTINATION ${CMAKE_INSTALL_BINDIR})
