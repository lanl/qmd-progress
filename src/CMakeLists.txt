include_directories(${BML_INCLUDEDIR})

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
  ptable_mod.F90
  pulaycomponent_mod.F90
  pulaymixer_mod.F90
  response_mod.F90
  sp2_fermi_mod.F90
  sp2_mod.F90
  sp2parser_mod.F90
  subgraphloop_mod.F90
  system_mod.F90
  timer_mod.F90
  xlbo_mod.F90
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
  charges_mod.mod
  densitymatrix_mod.mod
  dos_mod.mod
  extras_mod.mod
  genz_mod.mod
  graph_mod.mod
  graph_sp2parser_mod.mod
  homolumo_mod.mod
  initmatrices_mod.mod
  kernelparser_mod.mod
  nonortho_mod.mod
  normalize_mod.mod
  openfiles_mod.mod
  parallel_mod.mod
  partition_mod.mod
  progress_mod.mod
  ptable_mod.mod
  pulaycomponent_mod.mod
  pulaymixer_mod.mod
  response_mod.mod
  sp2_fermi_mod.mod
  sp2_mod.mod
  sp2parser_mod.mod
  subgraphloop_mod.mod
  system_mod.mod
  timer_mod.mod
  xlbo_mod.mod
  # LATTE modes
  coulomb_latte_mod.mod
  latteparser_latte_mod.mod
  ham_latte_mod.mod
  hsderivative_latte_mod.mod
  huckel_latte_mod.mod
  md_latte_mod.mod
  neighborlist_latte_mod.mod
  nonorthocoulombforces_latte_mod.mod
  slaterkosterforce_latte_mod.mod
  physmathconst_latte_mod.mod
  ppot_latte_mod.mod
  tbparams_latte_mod.mod
  DESTINATION include)
