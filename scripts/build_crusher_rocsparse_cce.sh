#!/bin/bash

rm -r build
rm -r install

# Set BML Library location
MY_PATH=`pwd`
export BML_DIR=${MY_PATH}/../bml/install
# Configuring PROGRESS with OpenMP
export CC=${CC:=cc}
export FC=${FC:=ftn}
export CXX=${CXX:=CC}
export BLAS_VENDOR=${BLAS_VENDOR:=OpenBLAS}
export BML_OPENMP=yes
export PROGRESS_OPENMP=${PROGRESS_OPENMP:=yes}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/install"}
export PROGRESS_GRAPHLIB=${PROGRESS_GRAPHLIB:=no}
export PROGRESS_TESTING=${PROGRESS_TESTING:=no}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=RelWithDebInfo}
export PROGRESS_EXAMPLES=${PROGRESS_EXAMPLES:=no}
export PROGRESS_BENCHMARKS=${PROGRESS_BENCHMARKS:=yes}
export EXTRA_FCFLAGS="-hsystem_alloc"
export EXTRA_LINK_FLAGS=""
./build.sh configure

