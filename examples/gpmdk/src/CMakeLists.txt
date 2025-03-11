cmake_minimum_required(VERSION 3.10.0)
project(appendix C CXX Fortran)

set(dir ${CMAKE_CURRENT_SOURCE_DIR}/build/)
set(CMAKE_BUILD_DIRECTORY ${dir})
set(CMAKE_CURRENT_BINARY_DIR ${dir})

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

find_library(FOUND_METIS metis)
  if(NOT FOUND_METIS)
    message(FATAL_ERROR "Could not find metis library")
  endif()
  message(STATUS "Found metis: ${FOUND_METIS}")
  add_definitions(-DDO_GRAPHLIB)
  list(APPEND LINK_LIBRARIES ${FOUND_METIS})

  get_filename_component(METIS_LIB ${FOUND_METIS} DIRECTORY)
  list(APPEND LINK_LIBRARIES "-L/${METIS_LIB} -lmetis")

if(PROGRESS_MPI)
  message(STATUS "Will build with MPI")
  add_definitions(-DDO_MPI)
endif()

if(SANITY_CHECK)
  message(STATUS "Will build with SANITY CHECK")
  add_definitions(-DSANITY_CHECK)
endif()
	
message(STATUS "Extra FC Flags ${EXTRA_FCFLAGS}")

if(DEFINED EXTRA_FCFLAGS)
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTRA_FCFLAGS}")
endif()


if(LIB)
  list(APPEND SHAREDLIB "-fPIC -shared")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTRA_FCFLAGS} ${SHAREDLIB}")
  add_definitions(-DLIBON)
endif()

message(STATUS "Project sources = " ${PROJECT_SOURCE_DIR} )
include_directories(${CMAKE_BINARY_DIR}/)
include_directories(${BML_INCLUDEDIR})
include_directories(${PROGRESS_INCLUDEDIR})

#add_subdirectory(hamiltonian)
#set(dir "${CMAKE_CURRENT_SOURCE_DIR}/hamiltonian")
if(LATTE)
  set(extras)
  add_definitions(-DUSE_LATTE)
  message(STATUS "Will build with LATTE modules in" "${LATTE_DIR}")
  if(${LATTE_PREC} STREQUAL "SINGLE")
    add_definitions(-DSINGLEPREC)
  else()
    add_definitions(-DDOUBLEPREC)
  endif()
  set(dir "${LATTE_DIR}/src/")

list(
  APPEND extras
  "${dir}/setuparray.F90"
  "${dir}/myprecision.F90"
  "${dir}/constants_mod.F90"     
  "${dir}/readtb.F90"
 "${dir}/spinarray.F90"
 "${dir}/restartarray.F90"
 "${dir}/kspacearray.F90"
 "${dir}/neblistarray.F90"
 "${dir}/readcr.F90"
 "${dir}/xboarray.F90"
 "${dir}/nonoarray.F90"
 "${dir}/univarray.F90"
 "${dir}/ppotarray.F90"
 "${dir}/virialarray.F90"
 "${dir}/mdarray.F90"
 "${dir}/univscaling.F90"
 "${dir}/univscaling_function.F90"
 "${dir}/wignerd.F90"
 "${dir}/wignerarray.F90"  
 "${dir}/univtailcoef.F90"
 "${dir}/latteparser_latte_mod.F90"
 "${dir}/wignerd.F90"
 "${dir}/wignerarray.F90"  
 "${dir}/univtailcoef.F90"
 "${dir}/latteparser_latte_mod.F90"
 "${dir}/coulombarray.F90"
 "${dir}/sparsearray.F90"
 "${dir}/relaxcommon.F90"
 "${dir}/openfiles_mod.F90"
 "${dir}/kernelparser_mod.F90"
 "${dir}/fermicommon.F90"
 "${dir}/errors.F90"
 "${dir}/timer_mod.F90"
 "${dir}/gradH.F90"
 "${dir}/dfdx.F90"
 "${dir}/am.F90"
 "${dir}/bm.F90"
 "${dir}/dunivscaling.F90"
 "${dir}/gethdim.F90"
 "${dir}/gencutofflist.F90"
 "${dir}/build_integral_map.F90"
 "${dir}/bldnewH.F90"
 "${dir}/genX.F90"
 "${dir}/angfactor.F90"
 "${dir}/genXprogress.F90"
 "${dir}/allocatenono.F90"
 "${dir}/genorbitallist.F90"
 "${dir}/getmatindlist.F90"
 "${dir}/neblists.F90"
 "${dir}/pbc.F90"
 "${dir}/genHonsite.F90"
 "${dir}/allocatenebarrays.F90"
 "${dir}/deallocatenono.F90"
 "${dir}/deallocatenebarrays.F90"
)
endif()


list(
  APPEND common_sources 
"gpmdcov_allocation.F90"
"gpmdcov_buildz.F90"
"gpmdcov_constraints.F90"
"gpmdcov_diagonalize.F90"
"gpmdcov_dm_min.F90"
"gpmdcov_dump.F90"
"gpmdcov_energandforces.F90"
"gpmdcov_finalize.F90"
"gpmdcov_firstcharges.F90"
"gpmdcov_graphpart.F90"
"gpmdcov_init.F90"
"gpmdcov_initparts.F90"
"gpmdcov_kernel.F90"
"gpmdcov_langevin.F90"
"gpmdcov_mdloop.F90"
"gpmdcov_assert.F90"
"gpmdcov_mod.F90"
"gpmdcov_part.F90"
"gpmdcov_preparemd.F90"
"gpmdcov_reshuffle.F90"
"gpmdcov_restart.F90"
"gpmdcov_rhosolver.F90"
"gpmdcov_test.F90"
"gpmdcov_vars.F90"
"gpmdcov_writeout.F90"
"gpmdcov_neighbor.F90"
"gpmdcov_parser.F90"
"gpmdcov_latte.F90"
"gpmdcov_rhosolver.F90"
"gpmdcov_nonequilibrium.F90"
"gpmdcov_constraints.F90"
"gpmdcov_highorder.F90"
"gpmdcov_field.F90"
"gpmdcov_nvtx.F90"
"gpmdcov_response.F90"
"gpmdcov_dos.F90"
"gpmdcov_dispersion.F90"
                            )


function(progress_appendix myappendix main_and_srcs)
list(GET main_and_srcs 0 main)
include_directories(${PROGRESS_INCLUDEDIR})
add_executable(${myappendix} ${main} ${common_sources} ${extras})
target_sources(${myappendix} PRIVATE ${ARGN})
target_link_libraries(${myappendix} PUBLIC
    ${LINK_LIBRARIES})
  set_target_properties(${myappendix}
    PROPERTIES
    LINK_FLAGS "")
  #add_subdirectory(hamiltonian)
endfunction(progress_appendix)

function(progress_appendix_library myappendix main_and_srcs)
list(GET main_and_srcs 0 main)
include_directories(${PROGRESS_INCLUDEDIR})
add_library(${myappendix} SHARED ${main} ${common_sources} ${extras})
target_sources(${myappendix} PRIVATE ${ARGN})
target_link_libraries(${myappendix} PUBLIC
    ${LINK_LIBRARIES})
  set_target_properties(${myappendix}
    PROPERTIES
    LINK_FLAGS "")
endfunction(progress_appendix_library)

if(GPMDK_NVTX)
  add_definitions(-DUSE_NVTX)
endif()

if(LIB)
  progress_appendix_library( gpmdk gpmdcov_lib.F90 gpmdcov_lib_mod.F90 )
  #progress_appendix( gpmd_a gpmd_secuential.F90 gpmdcov_lib_mod.F90)
  file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/gpmd.py
    DESTINATION ${PROJECT_BINARY_DIR})
else()
  progress_appendix( gpmdk gpmdk.F90 )
  #progress_appendix( gpmd_a gpmd_secuential.F90 gpmdcov_lib_mod.F90)
endif()

SET(DESTINATION ${PROJECT_BINARY_DIR})
install(FILES ${CMAKE_BINARY_DIR}/progress.pc
  DESTINATION ${PROJECT_BINARY_DIR})

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
