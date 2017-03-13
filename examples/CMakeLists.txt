include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/tests/src)
include_directories(${CMAKE_SOURCE_DIR}/build/src)
include_directories(${BML_INCLUDEDIR})

function(progress_example myexample path_to_example)
  add_executable(${myexample} ${path_to_example})
  target_link_libraries(${myexample} progress ${LINK_LIBRARIES})
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
progress_example(gpscf_dist gpdist/gpscf_dist.F90)
progress_example(gpmd_dist gpdist/gpmd_dist.F90)
progress_example(gploop gploop/gploop.F90)
progress_example(gpmd gpmd/gpmd.F90)
progress_example(gpscf gpscf/gpscf.F90)
progress_example(gptest gptest/gptest.F90)
progress_example(gpmdcov gpmdcov/gpmdcov.F90)
progress_example(gpscfcov gpscfcov/gpscfcov.F90)
progress_example(metis_partition Mtest/metis_partition.F90)
progress_example(sp2run sp2run/src/main.F90)
progress_example(scf scf/src/main.F90)
progress_example(mdresponse mdresponse/mdresponse.F90)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
