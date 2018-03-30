cmake_minimum_required(VERSION 2.8.7)
project(progress C CXX Fortran)

# The library version is versioned off the major version. If the API
# changes, the library version should be bumped.
set(PROJECT_VERSION_MAJOR "1")
set(PROJECT_VERSION_MINOR "1")
set(PROJECT_VERSION_PATCH "0")

set(PROJECT_VERSION
  "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")
set(PROJECT_DESCRIPTION
  "Progress Library "
  "A library for solvers used in quantum chemistry packages")

include(GNUInstallDirs)

set(CMAKE_INSTALL_PKG_CONFIG_DIR ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH TRUE)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
  message(STATUS "No user specified build type, using default...")
endif()
message(STATUS "Build type is ${CMAKE_BUILD_TYPE}")

string(TOUPPER ${CMAKE_BUILD_TYPE} TEMP)
if(TEMP STREQUAL DEBUG)
  add_definitions(-DGLOBAL_DEBUG=PROGRESS_LOG_DEBUG)
else()
  add_definitions(-DGLOBAL_DEBUG=PROGRESS_LOG_INFO)
endif()

set(GNU_C_FLAGS_DEBUG -O0 -g -Wall -Wimplicit -Wno-unknown-pragmas --coverage
  -save-temps -std=c99)
set(GNU_C_FLAGS_RELEASE -O3 -g -std=c99)
set(GNU_Fortran_FLAGS_DEBUG -O0 -g --coverage
  -ffree-line-length-none)
set(GNU_Fortran_FLAGS_RELEASE -O3 -g -ffree-line-length-none)

set(Intel_C_FLAGS_DEBUG -O0 -g -std=c99 -check=conversions,stack,uninit
  -traceback -qopt-report=5)
set(Intel_C_FLAGS_RELEASE -O2 -g -std=c99 -qopt-report=5)
set(Intel_Fortran_FLAGS_DEBUG -O0 -g -check all -assume realloc_lhs
  -traceback -qopt-report=5)
set(Intel_Fortran_FLAGS_RELEASE -O2 -g -assume realloc_lhs -qopt-report=5)

if(CMAKE_BUILD_TYPE)
  string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
  if(BUILD_TYPE STREQUAL "DEBUG" OR BUILD_TYPE STREQUAL "RELEASE")
    if(CMAKE_Fortran_FLAGS STREQUAL "")
      if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU"
          OR CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
        set(CMAKE_Fortran_FLAGS_${BUILD_TYPE}
          ${${CMAKE_Fortran_COMPILER_ID}_Fortran_FLAGS_${BUILD_TYPE}})
        string(REPLACE ";" " "
          CMAKE_Fortran_FLAGS_${BUILD_TYPE}
          "${CMAKE_Fortran_FLAGS_${BUILD_TYPE}}")
	if(DEFINED EXTRA_FCFLAGS)
          set(CMAKE_Fortran_FLAGS_${BUILD_TYPE} "${CMAKE_Fortran_FLAGS_${BUILD_TYPE}} ${EXTRA_FCFLAGS}")
	endif()
        message(STATUS "Setting Fortran compiler flags to "
          "${CMAKE_Fortran_FLAGS_${BUILD_TYPE}}")
      else()
        message(STATUS "Unknown Fortran compiler ${CMAKE_Fortran_COMPILER_ID}")
      endif()
    endif()
  endif()
endif()

if(CMAKE_Fortran_FLAGS STRGREATER "")
  message(STATUS "Using user supplied Fortran compiler flags ${CMAKE_Fortran_FLAGS}")
endif()

include(FindPkgConfig)

pkg_check_modules(BML bml REQUIRED)
message(STATUS "Found bml: ${BML_LDFLAGS}")
list(APPEND LINK_LIBRARIES ${BML_LDFLAGS})

set(BLAS_VENDOR ""
  CACHE STRING "If set, the preferred BLAS/LAPACK vendor. Possible choices: {Intel,MKL,ACML}")

set(PROGRESS_GRAPHLIB FALSE CACHE BOOL "Compile with externel graph libraries")
if(PROGRESS_GRAPHLIB)
  find_library(FOUND_METIS metis)
  if(NOT FOUND_METIS)
    message(FATAL_ERROR "Could not find metis library")
  endif()
  message(STATUS "Found metis: ${FOUND_METIS}")
  add_definitions(-DDO_GRAPHLIB)
  list(APPEND LINK_LIBRARIES ${FOUND_METIS})

  # Get some more information on how metis was built.
  find_file(FOUND_METIS_H "metis.h")

  file(STRINGS ${FOUND_METIS_H} IDXTYPEWIDTH REGEX "define IDXTYPEWIDTH")
  string(REPLACE " " ";" TEMP ${IDXTYPEWIDTH})
  list(GET TEMP 2 IDXTYPEWIDTH)
  message(STATUS "metis IDXTYPEWIDTH: ${IDXTYPEWIDTH}")

  if(IDXTYPEWIDTH EQUAL 32)
    add_definitions(-DMETIS_INDEX_KIND=4)
  else()
    add_definitions(-DMETIS_INDEX_KIND=8)
  endif()

  file(STRINGS ${FOUND_METIS_H} REALTYPEWIDTH REGEX "define REALTYPEWIDTH")
  string(REPLACE " " ";" TEMP ${REALTYPEWIDTH})
  list(GET TEMP 2 REALTYPEWIDTH)
  message(STATUS "metis REALTYPEWIDTH: ${REALTYPEWIDTH}")

  if(REALTYPEWIDTH EQUAL 32)
    add_definitions(-DMETIS_REAL_KIND=0e0)
  else()
    add_definitions(-DMETIS_REAL_KIND=0d0)
  endif()
else()
    add_definitions(-DMETIS_INDEX_KIND=4)
    add_definitions(-DMETIS_REAL_KIND=0d0)
endif()

set(BML_OPENMP TRUE CACHE BOOL "Whether the bml was built with OpenMP")
set(PROGRESS_OPENMP TRUE CACHE BOOL "Whether to compile with OpenMP support")
if(PROGRESS_OPENMP OR BML_OPENMP)
  include(FindOpenMP)
  if(OPENMP_FOUND)
    if(CMAKE_VERSION VERSION_LESS 3.1)
      # There is not Fortran support in FindOpenMP before version 3.1. We will
      # have to hack this.
      message(WARNING "Versions of <cmake-3.1 can not detect Fortran "
        "compiler OpenMP flags. We will assume that the flags found for "
        "the C compiler also work for the Fortran compiler. If this is "
        "not the case, please send email to <nbock@lanl.gov>.")
      set(OpenMP_Fortran_FLAGS ${OpenMP_C_FLAGS})
    endif()
  else()
    message(WARNING "Could not get the compilers to use OpenMP. "
      "Will pretend that this never happened and compile the library "
      "without OpenMP.")
    unset(OpenMP_Fortran_FLAGS)
  endif()
else()
  message(STATUS "Will not build with OpenMP")
endif()

set(PROGRESS_MPI FALSE CACHE BOOL "Compile with MPI support")
if(PROGRESS_MPI)
  message(STATUS "Will build with MPI")
  add_definitions(-DDO_MPI)
#  add_definitions(-DDO_MPI_BLOCK)
endif()

set(SANITY_CHECK FALSE CACHE BOOL "Additional sanity checks")
if(SANITY_CHECK)
  message(STATUS "Will add sanity checks")
  add_definitions(-DSANITY_CHECK)
endif()

add_subdirectory(src)

find_package(Doxygen)

if(DOXYGEN_FOUND)
  set(DOXYGEN_INPUT
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/tests/src)
  string(REPLACE ";" " " DOXYGEN_INPUT "${DOXYGEN_INPUT}")
  configure_file(docs/Doxyfile.in Doxyfile)
  add_custom_target(docs COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile)
endif(DOXYGEN_FOUND)

set(PROGRESS_TESTING FALSE
  CACHE BOOL "Whether to build the test suite.")

if(PROGRESS_TESTING)
  message(STATUS "Setting up test suite")
  find_program(VALGRIND valgrind
    DOC "valgrind - a suite of tools for debugging and profiling programs")
  if(VALGRIND)
    message(STATUS "Found valgrind, will test memory")
  endif()
  enable_testing()
  add_subdirectory(tests)
endif()

set(PROGRESS_EXAMPLES FALSE
  CACHE BOOL "Whether to build example programs.")

if(PROGRESS_EXAMPLES)
  message(STATUS "Setting up to build example programs")
  add_subdirectory(examples)
endif()

if(OPENMP_FOUND)
  message(STATUS "Additional Fortran compiler link flags: '${OpenMP_Fortran_FLAGS}'")
endif()

configure_file(src/progress.pc.in progress.pc)
install(FILES ${CMAKE_BINARY_DIR}/progress.pc
  DESTINATION ${CMAKE_INSTALL_PKG_CONFIG_DIR})
