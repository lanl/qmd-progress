#!/bin/bash
rm -r build
rm -r install

# Set BML Library location
MY_PATH=`pwd`
export BML_DIR=/autofs/nccs-svm1_proj/csc304/bml/crusher/install

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
export BML_PREFIX_PATH=${BML_PREFIX_PATH:=$BML_DIR}

./build.sh configure

