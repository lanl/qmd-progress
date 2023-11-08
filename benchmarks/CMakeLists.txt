include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${CMAKE_BINARY_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/tests/src)
include_directories(${BML_INCLUDEDIR})

function(progress_benchmark mybenchmark)
  add_executable(${mybenchmark} ${ARGN})
  target_link_libraries(${mybenchmark} PUBLIC
    progress
    ${LINK_LIBRARIES})
  set_target_properties(${mybenchmark}
    PROPERTIES
    LINK_FLAGS "--coverage")
  if(OPENMP_FOUND)
    set_target_properties(${mybenchmark}
      PROPERTIES
      COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
      LINK_FLAGS ${OpenMP_Fortran_FLAGS})
  endif()
  if(MPI_Fortran_FOUND AND MPI_Fortran_COMPILE_FLAGS)
    set_target_properties(${mybenchmark}
      PROPERTIES
      COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS}
      LINK_FLAGS ${MPI_Fortran_LINK_FLAGS})
  endif()
endfunction(progress_benchmark)

progress_benchmark(dmconstruction dmconstruction/dmconstruction.F90)

progress_benchmark(dmconstruction3d dmconstruction3d/dmconstruction3d.F90)

progress_benchmark(dmconstruction_gp dmconstruction_graphBased/dmconstruction_graphBased.F90)

progress_benchmark(dmconstruction_gpbio dmconstruction_graphBased/dmconstruction_graphBased_bio.F90 dmconstruction_graphBased/aux_mod.F90)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

install(TARGETS dmconstruction
        DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
