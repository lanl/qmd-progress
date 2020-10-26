#!/bin/bash

# Make a copy of this script and modify the METIS and BML library locations

rm -r build
rm -r install

# Set METIS and BML Library locations
METIS_LIB="$HOME/metis-5.1.0/build/Linux-x86_64/libmetis"
BML_LIB="$HOME/bml/install"

MY_PATH=`pwd`

# Configuring PROGRESS with OpenMP
export CC=${CC:=gcc}
export FC=${FC:=gfortran}
export CXX=${CXX:=g++}
export PROGRESS_OPENMP=${PROGRESS_OPENMP:=yes}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/install"}
export PROGRESS_GRAPHLIB=${PROGRESS_GRAPHLIB:=no}
export PROGRESS_SP2TC=${PROGRESS_SP2TC:=Fortran}
#export PROGRESS_SP2TC=${PROGRESS_SP2TC:=C++}
export PROGRESS_TESTING=${PROGRESS_TESTING:=no}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Debug}
export PROGRESS_EXAMPLES=${PROGRESS_EXAMPLES:=no}
export CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH:=$BML_LIB}
export EXTRA_FCFLAGS=${EXTRA_FCFLAGS:=""}
export EXTRA_LINK_FLAGS=${EXTRA_LINK_FLAGS:=""}
./build.sh configure 

# Make PROGRESS library and examples after running this script:
#   cd build
#   make
