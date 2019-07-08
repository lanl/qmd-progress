#!/bin/bash

# Make sure all the paths are correct

rm -r build
rm -r install

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MEMBERWORK/csc304/magma-2.5.0/lib

MY_PATH=$(pwd)
BML_LIB="${MY_PATH}/../bml/install"
MAGMA_PATH="${MY_PATH}/../magma-2.5.0/"
export ESSL_DIR=${ESSL_DIR:="${OLCF_ESSL_ROOT}"}
export CC=${CC:=xlC_r}
export FC=${FC:=xlf90_r}
#export FC=${FC:=mpif90}
export CXX=${CXX:=xlC_r}
export BLAS_VENDOR=${BLAS_VENDOR:=None}
export PROGRESS_OPENMP=${PROGRESS_OPENMP:=yes}
#export PROGRESS_MPI=${PROGRESS_MPI:=yes}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/install"}
export PROGRESS_TESTING=${PROGRESS_TESTING:=yes}
export PROGRESS_EXAMPLES=${PROGRESS_EXAMPLES:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Release}
export EXTRA_FCFLAGS=${EXTRA_FCFLAGS:="-I${MAGMA_PATH}/include -I${OLCF_CUDA_ROOT}/include -I${BML_LIB}/include/ -qessl -qstrict=all -qsmp=omp -O2 -qextname -qxlf2003=polymorphic -qthreaded"}
export EXTRA_LINK_FLAGS=${EXTRA_LINK_FLAGS:="-L${BML_LIB}/lib64/ -lbml_fortran -lbml -L${OLCF_ESSL_ROOT}/lib64/ -qsmp=omp -lessl -lesslsmp  -qextname -lxlopt -lxlf90_r -lxlfmath -lxl -lxlsmp -L${MAGMA_PATH}/lib/ -lmagma -L${OLCF_CUDA_ROOT}/lib64/ -lcublas -lcudart"}
export CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH:=$BML_LIB}



./build.sh configure

                                                                                                                                                                                              
                                    
