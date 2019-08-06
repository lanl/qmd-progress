#!/bin/bash

# Make sure all the paths are correct

rm -r build
rm -r install

MY_PATH=$(pwd)
BML_LIB="${MY_PATH}/../bml/install"
export ESSL_DIR=${ESSL_DIR:="${OLCF_ESSL_ROOT}"}
export MAGMA_PATH=${MAGMA_PATH:="$MEMBERWORK/csc304/magma-2.5.0"}
export CC=${CC:=gcc}
export FC=${FC:=gfortran}
#export FC=${FC:=mpif90}
export CXX=${CXX:=g++}
export BLAS_VENDOR=${BLAS_VENDOR:=Auto}
export PROGRESS_OPENMP=${PROGRESS_OPENMP:=yes}
#export PROGRESS_MPI=${PROGRESS_MPI:=yes}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/install"}
export PROGRESS_TESTING=${PROGRESS_TESTING:=yes}
export PROGRESS_EXAMPLES=${PROGRESS_EXAMPLES:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Release}


export EXTRA_FCFLAGS=${EXTRA_FCFLAGS:=" -I${MAGMA_PATH}/include -I${BML_LIB}/include/ -ffree-line-length-none -fopenmp -lpthread -I${OLCF_OPENBLAS_ROOT}/include  "}
export EXTRA_LINK_FLAGS=${EXTRA_LINK_FLAGS:="-fopenmp -L${BML_LIB}/lib64/ -lbml_fortran -lbml -L${MAGMA_PATH}/lib/ -lmagma  -fopenmp -fopenmp -lpthread -L${OLCF_OPENBLAS_ROOT}/lib64 -lopenblas "}
export CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH:=$BML_LIB}




./build.sh configure

                                                                                                                                                                                              
                                    
