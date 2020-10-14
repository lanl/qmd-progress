include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/sp2_tensorcore)
include_directories(${PROJECT_SOURCE_DIR}/src/sp2_tensorcore/src)
include_directories(${PROJECT_SOURCE_DIR}/build/src/sp2_tensorcore/src)
include_directories(${BML_INCLUDEDIR})

add_executable(DNN_run
  ./src/sp2_tensorcore/src/sp2_dnn_tcore.cu  
  ./src/sp2_tensorcore/src/tcore_hp_emulator.cu)
  #sp2_dnn_tcore.cu  
  #tcore_hp_emulator.cu)

#target_link_libraries(DNN_run progress ${LINK_LIBRARIES})

#set_target_properties(DNN_run
#  PROPERTIES
#  LINK_FLAGS "--coverage")
#if(OPENMP_FOUND)
#  set_target_properties(main
#    PROPERTIES
#    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
#    LINK_FLAGS ${OpenMP_Fortran_FLAGS})
#endif()
#if(MPI_Fortran_FOUND AND MPI_Fortran_COMPILE_FLAGS)
#  set_target_properties(main
#    PROPERTIES
#    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS}
#    LINK_FLAGS ${MPI_Fortran_LINK_FLAGS})
#endif()

#SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

