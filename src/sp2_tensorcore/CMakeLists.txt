include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/sp2_tensorcore)
include_directories(${PROJECT_SOURCE_DIR}/src/sp2_tensorcore/src)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})
#include_directories(${PROJECT_SOURCE_DIR}/src/sp2_tensorcore/lib)
#include_directories(${PROJECT_SOURCE_DIR}/build/src/sp2_tensorcore/src)
#include_directories(${BML_INCLUDEDIR})

message(STATUS "Compiling SP2 tensor core")
find_package(CUDA)
#findCUDAToolkit()
cuda_add_executable(DNN_run
  ${PROJECT_SOURCE_DIR}/src/sp2_tensorcore/src/sp2_dnn_tcore.cu  
  ${PROJECT_SOURCE_DIR}/src/sp2_tensorcore/src/tcore_hp_emulator.cu
  ${PROJECT_SOURCE_DIR}/src/sp2_tensorcore/src/linalg_tools.cu)

set(TCSP2_LIBRARIES "-L/usr/local/cuda-11.0/lib64/") 
set(TCSP2_LIBRARIES1 "${PROJECT_SOURCE_DIR}/src/sp2_tensorcore/src/linalg_tools.cu  ${PROJECT_SOURCE_DIR}/src/sp2_tensorcore/src/tcore_hp_emulator.cu")


set(TCSP2_INCLUDE "${PROJECT_SOURCE_DIR}/src/sp2_tensorcore/include") 

message(STATUS "Including the following directories")
message(STATUS ${TCSP2_INCLUDE})
message(STATUS "Linking the following directories" ${TCSP2_LIBRARIES})

target_include_directories(DNN_run PUBLIC ${TCSP2_INCLUDE} ${CUDA_INCLUDE_DIRS}  ${CUDA_CUBLAS_DIRS})

set(CUDA_TOOLKIT_ROOT_DIR $ENV{CUDA_LIB})
#target_link_libraries(DNN_run ${TCSP2_LIBRARIES} ${BLAS_LIBRARIES} ${CUDA_LIBRARIES} ${TCSP2_LIBRARIES1} -lcublas)
target_link_libraries(DNN_run -L${CUDA_TOOLKIT_ROOT_DIR} ${BLAS_LIBRARIES} ${CUDA_LIBRARIES} -lcublas)


