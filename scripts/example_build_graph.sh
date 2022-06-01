#!/bin/bash

# Make a copy of this script and modify the METIS and BML library locations

rm -r build
rm -r install

# Set METIS and BML Library locations
METIS_LIB="/tmp/metis-5.1.0/build/Linux-x86_64/libmetis"
METIS_INCLUDE="/tmp/metis-5.1.0/include"
BML_LIB="/tmp/bml/install"

MY_PATH=`pwd`

# Configuring PROGRESS with OpenMP
FC=mpif90
export CC=${CC:=gcc}
export FC=${FC:=mpif90}
export CXX=${CXX:=g++}
export PROGRESS_OPENMP=${PROGRESS_OPENMP:=yes}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/install"}
export PROGRESS_GRAPHLIB=${PROGRESS_GRAPHLIB:=yes}
export PROGRESS_TESTING=${PROGRESS_TESTING:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Release}
export PROGRESS_EXAMPLES=${PROGRESS_EXAMPLES:=yes}
export PKG_CONFIG_PATH="$BML_LIB/lib/pkgconfig:$BML_LIB/lib64/pkgconfig"
export CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH:="$BML_LIB:$METIS_LIB:$METIS_INCLUDE"}
export EXTRA_FCFLAGS=${EXTRA_FCFLAGS:=""}
export EXTRA_LINK_FLAGS=${EXTRA_LINK_FLAGS:=""}
./build.sh configure

# Make PROGRESS library and examples after running this script:
#   cd build
#   make
