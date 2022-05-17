#!/bin/bash

# Make a copy of this script and modify the PROGRESS_SP2TC flag as needed.
# Tested on V100 and A100 GPUs. See https://doi.org/10.1021/acs.jctc.1c00057
# for details.
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
export PROGRESS_SP2TC=${PROGRESS_SP2TC:=Fortran} # Choose Fortran or C++
export PROGRESS_TESTING=${PROGRESS_TESTING:=no}
export PROGRESS_EXAMPLES=${PROGRESS_EXAMPLES:=yes}
export CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH:=$BML_LIB}
export EXTRA_FCFLAGS=${EXTRA_FCFLAGS:="-llapack -lblas -fopenmp"}
export EXTRA_LINK_FLAGS=${EXTRA_LINK_FLAGS:="-llapack -lblas -fopenmp"}
./build.sh configure 

