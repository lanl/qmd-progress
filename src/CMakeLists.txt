include_directories(${BML_INCLUDEDIR})
include_directories(${BML_INCLUDEDIR}/latte_mods)
add_library(progress  
  charges_mod.F90 
  densitymatrix_mod.F90
  dos_mod.F90 
  extras_mod.F90 
  genz_mod.F90
  graph_mod.F90
  graph_sp2parser_mod.F90
  homolumo_mod.F90
  initmatrices_mod.F90
  kernelparser_mod.F90
  nonortho_mod.F90
  normalize_mod.F90
  openfiles_mod.F90
  parallel_mod.F90
  partition_mod.F90
  progress_mod.F90
  pulaycomponent_mod.F90
  pulaymixer_mod.F90  
  ptable_mod.F90
  response_mod.F90
  sp2_mod.F90
  sp2parser_mod.F90
  subgraphloop_mod.F90
  system_mod.F90
  timer_mod.F90
  xlbo_mod.F90
  # LATTE modes
  ./latte_mods/coulomb_latte_mod.F90
  ./latte_mods/latteparser_latte_mod.F90
  ./latte_mods/nonorthocoulombforces_latte_mod.F90
  ./latte_mods/ham_latte_mod.F90
  ./latte_mods/hsderivative_latte_mod.F90
  ./latte_mods/huckel_latte_mod.F90
  ./latte_mods/md_latte_mod.F90
  ./latte_mods/neighborlist_latte_mod.F90
  ./latte_mods/slaterkosterforce_latte_mod.F90
  ./latte_mods/physmathconst_latte_mod.F90
  ./latte_mods/ppot_latte_mod.F90
  ./latte_mods/tbparams_latte_mod.F90
  )

set_target_properties(progress
  PROPERTIES COMPILE_FLAGS ${OpenMP_Fortran_FLAGS})

if(EXTRA_LINK_FLAGS)
  message(STATUS "Linking with these extra flags: ${EXTRA_LINK_FLAGS}")
  target_link_libraries(progress ${EXTRA_LINK_FLAGS})
endif()

install(TARGETS progress
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} 
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/charges_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/densitymatrix_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/dos_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/extras_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/genz_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/graph_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/graph_sp2parser_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/homolumo_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/initmatrices_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/nonortho_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/normalize_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/kernelparser_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/openfiles_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/parallel_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/partition_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/pulaycomponent_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/pulaymixer_mod.mod  
  ${CMAKE_CURRENT_BINARY_DIR}/progress_mod.mod  
  ${CMAKE_CURRENT_BINARY_DIR}/ptable_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/response_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/sp2_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/sp2parser_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/subgraphloop_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/system_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/timer_mod.mod
  ${CMAKE_CURRENT_BINARY_DIR}/xlbo_mod.mod
  # LATTE modes
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
  DESTINATION include)


