#!/bin/bash

# Make a copy of this script and modify the METIS and BML library locations

rm -r build
rm -r install 

# Set METIS and BML Library locations
METIS_LIB="$HOME/metis-5.1.0/build/Linux-x86_64/libmetis"
BML_LIB="$HOME/bml/install/lib"

MY_PATH=`pwd`

# Configuring PROGRESS with OpenMP
CC=gcc FC=gfortran BLAS_VENDOR=GNU PKG_CONFIG_PATH=$BML_LIB/pkgconfig PROGRESS_OPENMP=yes INSTALL_DIR="$MY_PATH/install" PROGRESS_GRAPHLIB=no PROGRESS_TESTING=yes CMAKE_BUILD_TYPE=Release PROGRESS_EXAMPLES=yes ./build.sh configure


# Configuring PROGRESS with OpenMP, MPI and METIS Graph Library
#CC=mpicc FC=mpifort BLAS_VENDOR=GNU PKG_CONFIG_PATH=$BML_LIB/pkgconfig PROGRESS_OPENMP=yes PROGRESS_MPI=yes INSTALL_DIR="$MY_PATH/install" PROGRESS_GRAPHLIB=yes EXTRA_LINK_FLAGS="-L$METIS_LIB -lmetis" PROGRESS_TESTING=yes CMAKE_BUILD_TYPE=Release PROGRESS_EXAMPLES=yes ./build.sh configure

# Make PROGRESS library and examples
# cd build
# make
# make test
