#!/bin/bash

set -e -x

cd ~
git clone https://github.com/lanl/bml.git

cd bml
echo "Installing bml version $(git describe)"
CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:-Release} \
  BLAS_VENDOR=GNU \
  BML_OPENMP=yes \
  BML_TESTING=no \
  ./build.sh install
