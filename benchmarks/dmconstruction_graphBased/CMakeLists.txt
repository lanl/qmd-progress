cmake_minimum_required(VERSION 3.23)
project(appendix C CXX Fortran)

set(dir ${CMAKE_CURRENT_SOURCE_DIR}/build/)
set(CMAKE_BUILD_DIRECTORY ${dir})
set(CMAKE_CURRENT_BINARY_DIR ${dir})

include(FindPkgConfig)
if(NOT TARGET PROGRESS::progress)
  find_package(PROGRESS CONFIG QUIET)
endif()
if(TARGET PROGRESS::progress)
  message(STATUS "Found progress: PROGRESS::progress")
  list(APPEND LINK_LIBRARIES PROGRESS::progress)
else()
  pkg_check_modules(PROGRESS REQUIRED progress)
  message(STATUS "Found progress: ${PROGRESS_LDFLAGS}")
  list(APPEND LINK_LIBRARIES ${PROGRESS_LDFLAGS})
endif()

if(NOT TARGET BML::bml)
  find_package(BML CONFIG QUIET)
endif()
if(TARGET BML::bml)
  message(STATUS "Found bml: BML::bml")
  list(APPEND LINK_LIBRARIES BML::bml)
else()
  pkg_check_modules(BML REQUIRED bml)
  message(STATUS "Found bml: ${BML_LDFLAGS}")
  list(APPEND LINK_LIBRARIES ${BML_LDFLAGS})
endif()

find_library(FOUND_METIS metis)
list(APPEND LINK_LIBRARIES ${FOUND_METIS})

message(STATUS "Project sources = " ${PROJECT_SOURCE_DIR} )
include_directories(${PROJECT_SOURCE_DIR}/)
include_directories(${CMAKE_BINARY_DIR}/)
include_directories(${BML_INCLUDEDIR})
include_directories(${PROGRESS_INCLUDEDIR})

function(progress_appendix myappendix main_and_srcs)
list(GET main_and_srcs 0 main)
add_executable(${myappendix} ${main})
target_sources(${myappendix} PRIVATE ${ARGN})
target_link_libraries(${myappendix} PUBLIC
    ${LINK_LIBRARIES})
  set_target_properties(${myappendix}
    PROPERTIES
    LINK_FLAGS "--coverage")
endfunction(progress_appendix)

progress_appendix(main dmconstruction_graphBased_bio.F90
                            aux_mod.F90
                            )

SET(DESTINATION ${PROJECT_BINARY_DIR})
install(FILES ${CMAKE_BINARY_DIR}/progress.pc
  DESTINATION ${PROJECT_BINARY_DIR})

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
