include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/tests/src)
include_directories(${CMAKE_SOURCE_DIR}/build/src)
include_directories(${BML_INCLUDEDIR})

list(APPEND LINK_LIBRARIES "-L${BML_LIBDIR} -l${BML_LIBRARIES}")

add_executable(changecoords changecoords/changecoords.F90)
add_executable(gpscf_dist gpdist/gpscf_dist.F90)
add_executable(gpmd_dist gpdist/gpmd_dist.F90)
add_executable(gploop gploop/gploop.F90)
add_executable(gpmd gpmd/gpmd.F90)
add_executable(gpscf gpscf/gpscf.F90)
add_executable(gptest gptest/gptest.F90)
add_executable(gpmdcov gpmdcov/gpmdcov.F90)
add_executable(gpscfcov gpscfcov/gpscfcov.F90)
#add_executable(Mtest Mtest/Mtest.F90)
add_executable(metis_partition Mtest/metis_partition.F90)
add_executable(sp2run sp2run/src/main.F90)
add_executable(scf scf/src/main.F90)

target_link_libraries(changecoords progress ${LINK_LIBRARIES})
target_link_libraries(gpscf_dist progress ${LINK_LIBRARIES})
target_link_libraries(gpmd_dist progress ${LINK_LIBRARIES})
target_link_libraries(gploop progress ${LINK_LIBRARIES})
target_link_libraries(gpmd progress ${LINK_LIBRARIES})
target_link_libraries(gpscf progress ${LINK_LIBRARIES})
target_link_libraries(gptest progress ${LINK_LIBRARIES})
target_link_libraries(gpmdcov progress ${LINK_LIBRARIES})
target_link_libraries(gpscfcov progress ${LINK_LIBRARIES})
#target_link_libraries(Mtest progress ${LINK_LIBRARIES})
target_link_libraries(metis_partition progress ${LINK_LIBRARIES})
target_link_libraries(sp2run progress ${LINK_LIBRARIES})
target_link_libraries(scf progress ${LINK_LIBRARIES})

set_target_properties(changecoords
  PROPERTIES
  LINK_FLAGS "--coverage")
set_target_properties(gpscf_dist
  PROPERTIES
  LINK_FLAGS "--coverage")
set_target_properties(gpmd_dist
  PROPERTIES
  LINK_FLAGS "--coverage")
set_target_properties(gploop
  PROPERTIES
  LINK_FLAGS "--coverage")
set_target_properties(gpmd
  PROPERTIES
  LINK_FLAGS "--coverage")
set_target_properties(gpscf
  PROPERTIES
  LINK_FLAGS "--coverage")
set_target_properties(gptest
  PROPERTIES
  LINK_FLAGS "--coverage")
#  set_target_properties(Mtest
set_target_properties(gpmdcov
  PROPERTIES
  LINK_FLAGS "--coverage")
set_target_properties(gpscfcov
  PROPERTIES
  LINK_FLAGS "--coverage")
#  PROPERTIES
#  LINK_FLAGS "--coverage")
  set_target_properties(metis_partition
  PROPERTIES
  LINK_FLAGS "--coverage")
set_target_properties(sp2run
  PROPERTIES
  LINK_FLAGS "--coverage")
set_target_properties(scf
  PROPERTIES
  LINK_FLAGS "--coverage")

if(OPENMP_FOUND)
  set_target_properties(changecoords
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
    LINK_FLAGS ${OpenMP_Fortran_FLAGS})
  set_target_properties(gpscf_dist
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
    LINK_FLAGS ${OpenMP_Fortran_FLAGS})
  set_target_properties(gpmd_dist
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
    LINK_FLAGS ${OpenMP_Fortran_FLAGS})
  set_target_properties(gploop
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
    LINK_FLAGS ${OpenMP_Fortran_FLAGS})
  set_target_properties(gpmd
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
    LINK_FLAGS ${OpenMP_Fortran_FLAGS})
  set_target_properties(gpscf
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
    LINK_FLAGS ${OpenMP_Fortran_FLAGS})
  set_target_properties(gptest
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
    LINK_FLAGS ${OpenMP_Fortran_FLAGS})
  set_target_properties(gpmdcov
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
    LINK_FLAGS ${OpenMP_Fortran_FLAGS})
  set_target_properties(gpscfcov
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
    LINK_FLAGS ${OpenMP_Fortran_FLAGS})
  #  set_target_properties(Mtest
  #  PROPERTIES
  #  COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
  #  LINK_FLAGS ${OpenMP_Fortran_FLAGS})
  set_target_properties(metis_partition
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
    LINK_FLAGS ${OpenMP_Fortran_FLAGS})
  set_target_properties(sp2run
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
    LINK_FLAGS ${OpenMP_Fortran_FLAGS})
  set_target_properties(scf
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
    LINK_FLAGS ${OpenMP_Fortran_FLAGS})
endif()

if(MPI_Fortran_FOUND AND MPI_Fortran_COMPILE_FLAGS)
  set_target_properties(changecoords
    PROPERTIES
    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS}
    LINK_FLAGS ${MPI_Fortran_LINK_FLAGS})
  set_target_properties(gpscf_dist
    PROPERTIES
    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS}
    LINK_FLAGS ${MPI_Fortran_LINK_FLAGS})
  set_target_properties(gpmd_dist
    PROPERTIES
    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS}
    LINK_FLAGS ${MPI_Fortran_LINK_FLAGS})
  set_target_properties(gploop
    PROPERTIES
    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS}
    LINK_FLAGS ${MPI_Fortran_LINK_FLAGS})
  set_target_properties(gpmd
    PROPERTIES
    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS}
    LINK_FLAGS ${MPI_Fortran_LINK_FLAGS})
  set_target_properties(gpscf
    PROPERTIES
    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS}
    LINK_FLAGS ${MPI_Fortran_LINK_FLAGS})
  set_target_properties(gptest
    PROPERTIES
    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS}
    LINK_FLAGS ${MPI_Fortran_LINK_FLAGS})
  set_target_properties(gpmdcov
    PROPERTIES
    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS}
    LINK_FLAGS ${MPI_Fortran_LINK_FLAGS})
  set_target_properties(gpscfcov
    PROPERTIES
    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS}
    LINK_FLAGS ${MPI_Fortran_LINK_FLAGS})
  # set_target_properties(Mtest
  #  PROPERTIES
  #  COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS}
  #  LINK_FLAGS ${MPI_Fortran_LINK_FLAGS})
  set_target_properties(metis_partition
    PROPERTIES
    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS}
    LINK_FLAGS ${MPI_Fortran_LINK_FLAGS})
  set_target_properties(sp2run
    PROPERTIES
    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS}
    LINK_FLAGS ${MPI_Fortran_LINK_FLAGS})
  set_target_properties(scf
    PROPERTIES
    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS}
    LINK_FLAGS ${MPI_Fortran_LINK_FLAGS})
endif()

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
