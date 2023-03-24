#!/bin/bash
module load cmake
module load gcc/9.1.0

rm -r build
rm -r install

# Set BML Library location
BML_LIB=/autofs/nccs-svm1_proj/csc304/bml/summit_gcc_magma

MY_PATH=`pwd`

#get jsrun with full path
JSRUN=$(which jsrun)
echo ${JSRUN}

# Configuring PROGRESS with OpenMP
export CC=${CC:=gcc}
export FC=${FC:=gfortran}
export CXX=${CXX:=g++}
export PROGRESS_OPENMP=${PROGRESS_OPENMP:=yes}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/install"}
export PROGRESS_GRAPHLIB=${PROGRESS_GRAPHLIB:=no}
export PROGRESS_TESTING=${PROGRESS_TESTING:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=RelWithDebInfo}
export PROGRESS_EXAMPLES=${PROGRESS_EXAMPLES:=yes}
export PROGRESS_BENCHMARKS=${PROGRESS_BENCHMARKS:=yes}
export BML_PREFIX_PATH=${BML_PREFIX_PATH:=$BML_LIB}

export PROGRESS_NONMPI_PRECOMMAND=${PROGRESS_NONMPI_PRECOMMAND:=${JSRUN}}
export PROGRESS_NONMPI_PRECOMMAND_ARGS=${PROGRESS_NONMPI_PRECOMMAND_ARGS:="-n1;-a1;-g1;-c7;--smpiargs=off"}

./build.sh install
