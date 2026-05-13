include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/mlsp2_tensorcore)
include_directories(${PROJECT_SOURCE_DIR}/src/mlsp2_tensorcore/src)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

message(STATUS "Compiling ML-SP2 on Nvidia Tensor Cores!!")

find_package(CUDA REQUIRED)

# Compiling cuda files to make a C++ library
cuda_add_library(prg_mlsp2_tc
  ${PROJECT_SOURCE_DIR}/src/mlsp2_tensorcore/src/prg_mlsp2_tensorcore.cu  
  ${PROJECT_SOURCE_DIR}/src/mlsp2_tensorcore/src/tcore_hp_emulator.cu
  ${PROJECT_SOURCE_DIR}/src/mlsp2_tensorcore/src/linalg_tools.cu)

# Gathering cuda objects. Instead we are manually uning the objects
# that were generated when prg_mlsp2_tc library was added.
set(MLSP2_TC_O 
"${PROJECT_BINARY_DIR}/src/mlsp2_tensorcore/CMakeFiles/prg_mlsp2_tc.dir/src/prg_mlsp2_tc_generated_prg_mlsp2_tensorcore.cu.o"
"${PROJECT_BINARY_DIR}/src/mlsp2_tensorcore/CMakeFiles/prg_mlsp2_tc.dir/src/prg_mlsp2_tc_generated_linalg_tools.cu.o" 
"${PROJECT_BINARY_DIR}/src/mlsp2_tensorcore/CMakeFiles/prg_mlsp2_tc.dir/src/prg_mlsp2_tc_generated_tcore_hp_emulator.cu.o" 
)

# Making the Fortran lib out of the precompiled cuda objects.
add_library(prg_mlsp2_tc_fortran
  ${PROJECT_SOURCE_DIR}/src/mlsp2_tensorcore/src/prg_mlsp2_tensorcore_mod.F90
  ${MLSP2_TC_O})

set_target_properties(prg_mlsp2_tc_fortran PROPERTIES
  Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/src
)

set(MLSP2TC_INCLUDE "${PROJECT_SOURCE_DIR}/src/mlsp2_tensorcore/src") 

message(STATUS "Including the following directories" ${MLSP2TC_INCLUDE})

target_include_directories(prg_mlsp2_tc PUBLIC ${MLSP2TC_INCLUDE} ${CUDA_INCLUDE_DIRS} ${CUDA_CUBLAS_DIRS})
message(STATUS "Cuda include directories" ${CUDA_INCLUDE_DIRS})

set(CUDA_DIR $ENV{CUDADIR})
#set(CUBLAS_LIBRARIES "-L${CUDA_DIR}/lib64/ -L${CUDA_DIR}/lib64/")

target_link_libraries(prg_mlsp2_tc ${BLAS_LIBRARIES} ${CUDA_LIBRARIES} -lcublas -lcudadevrt -lcudart -lnvToolsExt)
target_link_libraries(prg_mlsp2_tc_fortran ${BLAS_LIBRARIES} ${CUDA_LIBRARIES} -lcublas -lcudadevrt -lcudart -lnvToolsExt)

set(TC_INCLUDEDIR ${PROJECT_SOURCE_DIR}/install/include)
set(TC_INSTALLDIR ${PROJCT_SOURCE_DIR}/install/lib64)

install(TARGETS prg_mlsp2_tc
	OBJECTS ARCHIVE DESTINATION ${TC_INSTALLDIR})
install(TARGETS prg_mlsp2_tc_fortran
	OBJECTS ARCHIVE DESTINATION ${TC_INSTALLDIR})

install(FILES
        ${PROJECT_BINARY_DIR}/src/prg_mlsp2_tensorcore_mod.mod
        ${MLSP2_TC_O}
	${MLSP2TC_INCLUDE}/tcore_hp_emulator.cuh
	${MLSP2TC_INCLUDE}/linalg_tools.cuh
	${MLSP2TC_INCLUDE}/prg_mlsp2_tensorcore.cuh
	DESTINATION ${TC_INCLUDEDIR})


