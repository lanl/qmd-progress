#!/bin/bash

# Make a copy of this script and modify the METIS and BML library locations

rm -r build
rm -r install

# Set METIS and BML Library locations
METIS_LIB="$HOME/metis-5.1.0/build/Linux-x86_64/libmetis"
BML_LIB="$HOME/bml/install"

MY_PATH=`pwd`

# Configuring PROGRESS with OpenMP
export CC=${CC:=xlC_r}
export FC=${FC:=xlf90_r}
export CXX=${CXX:=xlC_r}
export PROGRESS_OPENMP=${PROGRESS_OPENMP:=yes}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/install"}
export PROGRESS_GRAPHLIB=${PROGRESS_GRAPHLIB:=no}
export PROGRESS_TESTING=${PROGRESS_TESTING:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Release}
export PROGRESS_EXAMPLES=${PROGRESS_EXAMPLES:=yes}
export PROGRESS_EXAMPLES=${PROGRESS_EXAMPLES:=yes}
export CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH:=$BML_LIB}
export EXTRA_FCFLAGS=${EXTRA_FCFLAGS:="-I$BML_LIB/include"}
export EXTRA_LINK_FLAGS=${EXTRA_LINK_FLAGS:="-L$BML_LIB/lib64/ -lbml_fortran -lbml"}
./build.sh configure

# Make PROGRESS library and examples after running this script:
#   cd build
#   make
