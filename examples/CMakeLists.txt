include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${CMAKE_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/tests/src)
include_directories(${BML_INCLUDEDIR})

function(progress_example myexample main_and_srcs)
  list(GET main_and_srcs 0 main)
  add_executable(${myexample} ${main} )
  target_sources(${myexample} PRIVATE ${ARGN})
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


progress_example(changecoords changecoords/changecoords.F90)
progress_example(getdihedral getdihedral/getdihedral.F90)
progress_example(gpscf_dist gpdist/gpscf_dist.F90)
progress_example(gpmd_dist gpdist/gpmd_dist.F90)
progress_example(gploop gploop/gploop.F90)
progress_example(gpmd gpmd/gpmd.F90)
progress_example(gpscf gpscf/gpscf.F90)
progress_example(gptest gptest/gptest.F90)
progress_example(gpmdcov gpmdcov/gpmdcov.F90 
		gpmdcov/gpmdcov_mod.F90
		gpmdcov/gpmdcov_vars.F90
		gpmdcov/gpmdcov_init.F90
		gpmdcov/gpmdcov_restart.F90
		gpmdcov/gpmdcov_part.F90
		gpmdcov/gpmdcov_reshuffle.F90
		gpmdcov/gpmdcov_graphpart.F90
		gpmdcov/gpmdcov_writeout.F90
		gpmdcov/gpmdcov_initparts.F90
		gpmdcov/gpmdcov_buildz.F90
		gpmdcov/gpmdcov_rhosolver.F90
		gpmdcov/gpmdcov_firstcharges.F90
		gpmdcov/gpmdcov_dm_min.F90
		gpmdcov/gpmdcov_energandforces.F90
                gpmdcov/gpmdcov_preparemd.F90
                gpmdcov/gpmdcov_mdloop.F90
                gpmdcov/gpmdcov_dump.F90
                gpmdcov/gpmdcov_finalize.F90)
progress_example(gpscfcov gpscfcov/gpscfcov.F90)
progress_example(gpmdcov_original gpmdcov/gpmdcov_original.F90
		gpmdcov/gpmdcov_vars.F90
		gpmdcov/gpmdcov_rhosolver.F90)
progress_example(metis_partition Mtest/metis_partition.F90)
progress_example(moveatoms moveatoms/moveatoms.F90)
progress_example(replicate replicate/replicate.F90)
progress_example(sp2run sp2run/src/main.F90)
progress_example(scf scf/src/main.F90)
progress_example(mdresponse mdresponse/mdresponse.F90)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
