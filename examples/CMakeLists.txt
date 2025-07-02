include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${CMAKE_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/tests/src)
include_directories(${BML_INCLUDEDIR})

function(progress_example myexample path_to_example)
  add_executable(${myexample} ${path_to_example})
  target_link_libraries(${myexample} PUBLIC
    progress
    ${LINK_LIBRARIES})
  set_target_properties(${myexample}
    PROPERTIES
    LINK_FLAGS "--coverage")
  if(OPENMP_FOUND)
    set_target_properties(${myexample}
      PROPERTIES
      COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
      LINK_FLAGS ${OpenMP_Fortran_FLAGS})
  endif()
  if(MPI_Fortran_FOUND AND MPI_Fortran_COMPILE_FLAGS)
    set_target_properties(${myexample}
      PROPERTIES
      COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS}
      LINK_FLAGS ${MPI_Fortran_LINK_FLAGS})
  endif()
endfunction(progress_example)

function(progress_gpmdk_example myexample main_and_srcs)
list(GET main_and_srcs 0 main)
list(
  APPEND common_sources 
"gpmdk/src/gpmdcov_allocation.F90"
"gpmdk/src/gpmdcov_buildz.F90"
"gpmdk/src/gpmdcov_constraints.F90"
"gpmdk/src/gpmdcov_diagonalize.F90"
"gpmdk/src/gpmdcov_dm_min.F90"
"gpmdk/src/gpmdcov_dump.F90"
"gpmdk/src/gpmdcov_energandforces.F90"
"gpmdk/src/gpmdcov_finalize.F90"
"gpmdk/src/gpmdcov_firstcharges.F90"
"gpmdk/src/gpmdcov_graphpart.F90"
"gpmdk/src/gpmdcov_init.F90"
"gpmdk/src/gpmdcov_initparts.F90"
"gpmdk/src/gpmdcov_kernel.F90"
"gpmdk/src/gpmdcov_langevin.F90"
"gpmdk/src/gpmdcov_mdloop.F90"
"gpmdk/src/gpmdcov_assert.F90"
"gpmdk/src/gpmdcov_mod.F90"
"gpmdk/src/gpmdcov_part.F90"
"gpmdk/src/gpmdcov_preparemd.F90"
"gpmdk/src/gpmdcov_reshuffle.F90"
"gpmdk/src/gpmdcov_restart.F90"
"gpmdk/src/gpmdcov_rhosolver.F90"
"gpmdk/src/gpmdcov_test.F90"
"gpmdk/src/gpmdcov_vars.F90"
"gpmdk/src/gpmdcov_writeout.F90"
"gpmdk/src/gpmdcov_neighbor.F90"
"gpmdk/src/gpmdcov_parser.F90"
"gpmdk/src/gpmdcov_latte.F90"
"gpmdk/src/gpmdcov_rhosolver.F90"
"gpmdk/src/gpmdcov_nonequilibrium.F90"
"gpmdk/src/gpmdcov_constraints.F90"
"gpmdk/src/gpmdcov_highorder.F90"
"gpmdk/src/gpmdcov_field.F90"
"gpmdk/src/gpmdcov_response.F90"
"gpmdk/src/gpmdcov_dos.F90"
"gpmdk/src/gpmdcov_dispersion.F90"

)

  if(PROGRESS_GPMDK_LIB)
    add_library(${myexample} SHARED ${main} ${common_sources})
    set_target_properties(${myexample}
      PROPERTIES
      COMPILE_FLAGS "-fPIC -shared"
      LINK_FLAGS "-fPIC -shared")
    add_definitions(-DLIBON)
    file(COPY gpmdk/src/gpmd.py
         DESTINATION ${PROJECT_BINARY_DIR})
  else()
    add_executable(${myexample} ${main} ${common_sources})
  endif()
  target_link_libraries(${myexample} PUBLIC
    progress
    ${LINK_LIBRARIES})
  set_target_properties(${myexample}
    PROPERTIES
    LINK_FLAGS "--coverage")
  if(OPENMP_FOUND)
    set_target_properties(${myexample}
      PROPERTIES
      COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
      LINK_FLAGS ${OpenMP_Fortran_FLAGS})
  endif()
  if(MPI_Fortran_FOUND AND MPI_Fortran_COMPILE_FLAGS)
    set_target_properties(${myexample}
      PROPERTIES
      COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS}
      LINK_FLAGS ${MPI_Fortran_LINK_FLAGS})
  endif()
endfunction(progress_gpmdk_example)


if(PROGRESS_CINTERFACE)
  function(progress_c_example myexample path_to_example)
    add_executable(${myexample} ${path_to_example})
    target_link_libraries(${myexample} PUBLIC
      progress
      ${LINK_LIBRARIES})
    set_target_properties(${myexample}
      PROPERTIES
      LINK_FLAGS "--coverage")
    if(OPENMP_FOUND)
      set_target_properties(${myexample}
        PROPERTIES
        COMPILE_FLAGS ${OpenMP_C_FLAGS}
        LINK_FLAGS ${OpenMP_C_FLAGS})
    endif()
    if(MPI_C_FOUND AND MPI_C_COMPILE_FLAGS)
      set_target_properties(${myexample}
        PROPERTIES
        COMPILE_FLAGS ${MPI_C_COMPILE_FLAGS}
        LINK_FLAGS ${MPI_C_LINK_FLAGS})
    endif()
  endfunction(progress_c_example)
endif()


progress_example(changecoords changecoords/changecoords.F90)
progress_example(getdihedral getdihedral/getdihedral.F90)
progress_example(gpscf_dist gpdist/gpscf_dist.F90)
#MEW Comment out as doesn't build on crusher
#progress_example(gpmd_dist gpdist/gpmd_dist.F90)
progress_example(gploop gploop/gploop.F90)
progress_example(gpmd gpmd/gpmd.F90)
progress_example(gpscf gpscf/gpscf.F90)
progress_example(gptest gptest/gptest.F90)
progress_example(gpmdcov gpmdcov/gpmdcov.F90)
progress_example(gpscfcov gpscfcov/gpscfcov.F90)
progress_example(metis_partition Mtest/metis_partition.F90)
progress_example(moveatoms moveatoms/moveatoms.F90)
progress_example(replicate replicate/replicate.F90)
progress_example(sp2run sp2run/src/main.F90)
progress_example(scf scf/src/main.F90)
progress_example(mdresponse mdresponse/mdresponse.F90)
if(PROGRESS_GPMDK)
  if(PROGRESS_GPMDK_LIB)
    progress_gpmdk_example(gpmdk gpmdk/src/gpmdcov_lib_mod.F90 gpmdk/src/gpmdcov_lib.F90)
  else()
    progress_gpmdk_example(gpmdk gpmdk/src/gpmdk.F90)
  endif()
endif()

if(PROGRESS_CINTERFACE)
  progress_c_example(cinterface c_interface/progress.c)
endif()

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
