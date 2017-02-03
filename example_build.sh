#!/bin/bash

# Make a copy of this script and modify the METIS and BML library locations

rm -r build
rm -r install 

# Set METIS and BML Library locations
METIS_LIB="/usr/local/lib"
BML_LIB="$HOME/bml/lib"

MY_PATH=`pwd`

# Configuring PROGRESS with OpenMP and METIS Graph Library
CC=mpicc FC=mpif90 BLAS_VENDOR=MKL PKG_CONFIG_PATH=$BML_LIB/pkgconfig PROGRESS_OPENMP=yes INSTALL_DIR="$MY_PATH/install" PROGRESS_GRAPHLIB=yes EXTRA_LINK_FLAGS="-L$METIS_LIB -lmetis" PROGRESS_TESTING=yes CMAKE_BUILD_TYPE=Release PROGRESS_EXAMPLES=yes ./build.sh configure

# Make PROGRESS library and examples
# cd build
# make
# make test
