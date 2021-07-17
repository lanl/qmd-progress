#!/bin/bash
module load cmake
module load cuda
module load gcc/8.1.1
module load netlib-lapack
module load openblas
module load magma

rm -r build
rm -r install

# Set BML Library location
BML_LIB=/gpfs/alpine/csc304/proj-shared/bml

MY_PATH=`pwd`

# Configuring PROGRESS with OpenMP
export MAGMA_PATH=${MAGMA_PATH:=${OLCF_MAGMA_ROOT}}
export CC=${CC:=gcc}
export FC=${FC:=gfortran}
export CXX=${CXX:=g++}
export PROGRESS_OPENMP=${PROGRESS_OPENMP:=yes}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/install"}
export PROGRESS_GRAPHLIB=${PROGRESS_GRAPHLIB:=no}
export PROGRESS_TESTING=${PROGRESS_TESTING:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=RelWithDebInfo}
export PROGRESS_EXAMPLES=${PROGRESS_EXAMPLES:=yes}
export PKG_CONFIG_PATH=${BML_LIB}/lib64/pkgconfig

./build.sh install
