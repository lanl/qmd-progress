#!/bin/bash

set -e -x

cd ~
git clone --depth=1 https://github.com/lanl/bml.git

cd bml
CMAKE_BUILD_TYPE=Release \
  BLAS_VENDOR=GNU \
  BML_OPENMP=yes \
  BML_TESTING=no \
  ./build.sh install
