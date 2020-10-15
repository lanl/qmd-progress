include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/sp2_tensorcore)
include_directories(${PROJECT_SOURCE_DIR}/src/sp2_tensorcore/src)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

message(STATUS "Compiling SP2 tensor core")

find_package(CUDA)

cuda_add_executable(DNN_run
  ${PROJECT_SOURCE_DIR}/src/sp2_tensorcore/src/sp2_dnn_tcore.cu  
  ${PROJECT_SOURCE_DIR}/src/sp2_tensorcore/src/tcore_hp_emulator.cu
  ${PROJECT_SOURCE_DIR}/src/sp2_tensorcore/src/linalg_tools.cu)

set(TCSP2_INCLUDE "${PROJECT_SOURCE_DIR}/src/sp2_tensorcore/include") 

message(STATUS "Including the following directories" ${TCSP2_INCLUDE})

target_include_directories(DNN_run PUBLIC ${TCSP2_INCLUDE} ${CUDA_INCLUDE_DIRS}  ${CUDA_CUBLAS_DIRS})

set(CUDA_TOOLKIT_ROOT_DIR $ENV{CUDA_LIB})

target_link_libraries(DNN_run -L${CUDA_TOOLKIT_ROOT_DIR} ${BLAS_LIBRARIES} ${CUDA_LIBRARIES} -lcublas)


