#!/bin/bash

# Replace USER by your user name. 

rm -r build

module load gcc/5.4.0
module load openmpi/1.10.3-gcc_5.4.0
module load cmake
module load mkl
source /projects/opt/intel/parallel_studio_xe_2016/mkl/bin/mklvars.sh intel64

CC=gcc \
FC=gfortran \
BLAS_VENDOR=MKL \
PKG_CONFIG_PATH=/home/USER/bml/install/lib64/pkgconfig PROGRESS_OPENMP=yes \
PROGRESS_GRAPHLIB=yes \
EXTRA_LINK_FLAGS="-L/home/USER/metis-4.1.0/build/Linux-x86_64/libmetis
-lmetis" \
PROGRESS_TESTING=yes CMAKE_BUILD_TYPE=Release FC=gfortran \
PROGRESS_EXAMPLES=yes \
./build.sh configure


