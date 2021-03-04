#!/bin/bash

# Make a copy of this script and modify the METIS and BML library locations

rm -r build
rm -r install

# Set METIS and BML Library locations

METIS_LIB="$HOME/BESGraph/metis-5.1.0/"
BML_LIB="$HOME/BESGraph/bml/install"

MY_PATH=`pwd`

# Configuring PROGRESS with OpenMP
export CC=${CC:=gcc}
export FC=${FC:=mpif90}
export CXX=${CXX:=g++}
export PROGRESS_OPENMP=${PROGRESS_OPENMP:=yes}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/install"}
export PROGRESS_GRAPHLIB=${PROGRESS_GRAPHLIB:=yes}
export PROGRESS_MPI=${PROGRESS_MPI:=yes}
export PROGRESS_TESTING=${PROGRESS_TESTING:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Release}
export PROGRESS_EXAMPLES=${PROGRESS_EXAMPLES:=yes}
export CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH:=$BML_LIB}
#export EXTRA_FCFLAGS=${EXTRA_FCFLAGS:="-I$MPI_INCLUDE"}
#export EXTRA_LINK_FLAGS=${EXTRA_LINK_FLAGS:="-L$MPI_LIB"}
export CMAKE_INCLUDE_PATH=${CMAKE_INCLUDE_PATH:="$METIS_LIB/include"}
export CMAKE_LIBRARY_PATH=${CMAKE_LIBRARY_PATH:="$METIS_LIB/lib"}
./build.sh configure

# Make PROGRESS library and examples after running this script:
#   cd build
#   make
