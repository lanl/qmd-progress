#!/bin/bash

set -ev

cd ~
git clone --depth=1 https://github.com/lanl/bml.git;
cd bml
CC=gcc-6 FC=gfortran-6 CXX=g++-6 \
  CMAKE_BUILD_TYPE=Release BLAS_VENDOR=GNU \
  BML_OPENMP=yes BML_TESTING=no \
  ./build.sh install
