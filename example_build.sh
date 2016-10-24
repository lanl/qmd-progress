#!/bin/bash

# Check metis and bml path. 

rm -r build
rm -r install 

MY_PATH=`pwd`
cd ../
HOME_PATH=`pwd`
cd $MY_PATH

CC=gcc FC=gfortran BLAS_VENDOR=MKL PKG_CONFIG_PATH=$HOME_PATH/bml/install/lib/pkgconfig PROGRESS_OPENMP=yes INSTALL_DIR="$MY_PATH/install" PROGRESS_GRAPHLIB=yes EXTRA_LINK_FLAGS="-L/usr/local/lib -lmetis" PROGRESS_TESTING=yes CMAKE_BUILD_TYPE=Release FC=gfortran PROGRESS_EXAMPLES=yes ./build.sh configure


