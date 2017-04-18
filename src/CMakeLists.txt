include_directories(${BML_INCLUDEDIR})

add_library(progress
  prg_charges_mod.F90
  prg_densitymatrix_mod.F90
  prg_dos_mod.F90
  prg_extras_mod.F90
  prg_genz_mod.F90
  prg_graph_mod.F90
  prg_graphsp2parser_mod.F90
  prg_homolumo_mod.F90
  prg_initmatrices_mod.F90
  prg_kernelparser_mod.F90
  prg_memory_consumption.c
  prg_nonortho_mod.F90
  prg_normalize_mod.F90
  prg_openfiles_mod.F90
  prg_parallel_mod.F90
  prg_partition_mod.F90
  prg_progress_mod.F90
  prg_ptable_mod.F90
  prg_pulaycomponent_mod.F90
  prg_pulaymixer_mod.F90
  prg_response_mod.F90
  prg_sp2_fermi_mod.F90
  prg_sp2_mod.F90
  prg_sp2parser_mod.F90
  prg_subgraphloop_mod.F90
  prg_system_mod.F90
  prg_timer_mod.F90
  prg_xlbo_mod.F90
  # LATTE modules
  ./latte_mods/coulomb_latte_mod.F90
  ./latte_mods/ham_latte_mod.F90
  ./latte_mods/hsderivative_latte_mod.F90
  ./latte_mods/huckel_latte_mod.F90
  ./latte_mods/latteparser_latte_mod.F90
  ./latte_mods/md_latte_mod.F90
  ./latte_mods/neighborlist_latte_mod.F90
  ./latte_mods/nonorthocoulombforces_latte_mod.F90
  ./latte_mods/physmathconst_latte_mod.F90
  ./latte_mods/ppot_latte_mod.F90
  ./latte_mods/slaterkosterforce_latte_mod.F90
  ./latte_mods/tbparams_latte_mod.F90)
  # Wrappers 

if(OPENMP_FOUND)
  set_target_properties(progress
    PROPERTIES COMPILE_FLAGS ${OpenMP_Fortran_FLAGS})
endif()

if(EXTRA_LINK_FLAGS)
  message(STATUS "Linking with these extra flags: ${EXTRA_LINK_FLAGS}")
  target_link_libraries(progress ${EXTRA_LINK_FLAGS})
endif()

install(TARGETS progress
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/prg_charges_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_densitymatrix_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_dos_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_extras_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_genz_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_graph_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_graphsp2parser_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_homolumo_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_initmatrices_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_kernelparser_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_nonortho_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_normalize_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_openfiles_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_parallel_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_partition_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_progress_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_ptable_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_pulaycomponent_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_pulaymixer_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_response_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_sp2_fermi_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_sp2_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_sp2parser_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_subgraphloop_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_system_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_timer_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/prg_xlbo_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/coulomb_latte_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/latteparser_latte_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/ham_latte_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/hsderivative_latte_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/huckel_latte_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/md_latte_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/neighborlist_latte_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/nonorthocoulombforces_latte_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/slaterkosterforce_latte_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/physmathconst_latte_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/ppot_latte_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/tbparams_latte_mod.mod
  ./wrappers/prg_rhosolver_inc.F90

  DESTINATION include)
