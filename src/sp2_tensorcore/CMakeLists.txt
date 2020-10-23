include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/sp2_tensorcore)
include_directories(${PROJECT_SOURCE_DIR}/src/sp2_tensorcore/src)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

message(STATUS "Compiling SP2 tensor core")

find_package(CUDA)

cuda_add_library(prg_sp2_tc 
  ${PROJECT_SOURCE_DIR}/src/sp2_tensorcore/src/prg_sp2_tensorcore.cu  
  ${PROJECT_SOURCE_DIR}/src/sp2_tensorcore/src/tcore_hp_emulator.cu
  ${PROJECT_SOURCE_DIR}/src/sp2_tensorcore/src/linalg_tools.cu)

set(TCSP2_INCLUDE "${PROJECT_SOURCE_DIR}/src/sp2_tensorcore/include") 

message(STATUS "Including the following directories" ${TCSP2_INCLUDE})

target_include_directories(prg_sp2_tc PUBLIC ${TCSP2_INCLUDE} ${CUDA_INCLUDE_DIRS}  ${CUDA_CUBLAS_DIRS})
message(STATUS "Cuda include directories" ${CUDA_INCLUDE_DIRS})

set(CUDA_DIR $ENV{CUDADIR})
set(CUBLAS_LIBRARIES "-L${CUDA_DIR}/lib64/ -L${CUDA_DIR}/lib64/")

target_link_libraries(prg_sp2_tc ${CUBLAS_LIBRARIES} ${BLAS_LIBRARIES} ${CUDA_LIBRARIES} -lcublas -lcudadevrt -lcudart)

set(TC_INCLUDEDIR ${PROJECT_SOURCE_DIR}/install/include)
set(TC_INSTALLDIR ${PROJECT_SOURCE_DIR}/install/lib64)

install(TARGETS prg_sp2_tc
	OBJECTS ARCHIVE DESTINATION ${TC_INSTALLDIR})

install(FILES
	${TCSP2_INCLUDE}/tcore_hp_emulator.cuh
	${TCSP2_INCLUDE}/linalg_tools.cuh
	${TCSP2_INCLUDE}/prg_sp2_tensorcore.cuh
	DESTINATION ${TC_INCLUDEDIR})
