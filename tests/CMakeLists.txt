include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/tests/src)
include_directories(${CMAKE_SOURCE_DIR}/build/src)
include_directories(${BML_INCLUDEDIR})

add_executable(main
  src/prec.F90
  src/main.F90
  src/test_subgraphloop.F90)

target_link_libraries(main progress ${LINK_LIBRARIES})

set_target_properties(main
  PROPERTIES
  LINK_FLAGS "--coverage")
if(OPENMP_FOUND)
  set_target_properties(main
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
    LINK_FLAGS ${OpenMP_Fortran_FLAGS})
endif()
if(MPI_Fortran_FOUND AND MPI_Fortran_COMPILE_FLAGS)
  set_target_properties(main
    PROPERTIES
    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS}
    LINK_FLAGS ${MPI_Fortran_LINK_FLAGS})
endif()

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

enable_testing()

function(progress_test mytest)
  #if (PROGRESS_MPI)
  #  add_test(NAME ${mytest}
  #   COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/main  ${mytest}
  #    message(${msg}))
  #else()
    add_test(NAME ${mytest}
      COMMAND ${CMAKE_BINARY_DIR}/main  ${mytest}
      message(${msg}))
  #endif()
endfunction(progress_test)

function(copy_data data_folder)
  set(S ${CMAKE_SOURCE_DIR}/tests/${data_folder})
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${S}  ${CMAKE_BINARY_DIR}/tests/)
endfunction(copy_data)

copy_data(tests_data)

# PROGRESS modules
progress_test(prg_density_cheb_fermi)
progress_test(prg_density)
progress_test(prg_density_T)
progress_test(prg_density_T_fermi)
progress_test(prg_sp2_basic)
progress_test(prg_implicit_fermi)
progress_test(prg_sp2_alg1_dense)
progress_test(prg_sp2_alg2_dense)
progress_test(prg_sp2_alg1_ellpack)
progress_test(prg_sp2_alg2_ellpack)
progress_test(prg_sp2_alg2_ellpack_poly)
progress_test(prg_sp2_alg1_seq_dense)
progress_test(prg_sp2_alg2_seq_dense)
progress_test(prg_sp2_alg1_seq_ellpack)
progress_test(prg_sp2_alg2_seq_ellpack)
progress_test(prg_sp2_alg1_seq_inplace_dense)
progress_test(prg_sp2_alg2_seq_inplace_dense)
progress_test(prg_sp2_alg1_seq_inplace_ellpack)
progress_test(prg_sp2_alg2_seq_inplace_ellpack)
progress_test(prg_sp2_fermi_dense)
progress_test(prg_equal_partition)
progress_test(prg_file_partition)
progress_test(prg_subgraphsp2_equal)
progress_test(prg_deorthogonalize_dense)
progress_test(prg_orthogonalize_dense)
progress_test(prg_buildzdiag)
progress_test(prg_buildzsparse)

progress_test(prg_system_parse_write_xyz)
progress_test(prg_system_parse_write_pdb)
progress_test(prg_system_parse_write_dat)

# LATTE modules
progress_test(load_tbparms_latte)
progress_test(load_bintTBparamsH)

# EXTERNAL modules
#progress_test(get_hshuckel)
