#!/bin/bash

set -e -u -x

if [[ -v TEST_SETTINGS ]]; then
  [[ -f ${TEST_SETTINGS} ]] && source ${TEST_SETTINGS}
fi

cd ~
git clone https://github.com/lanl/bml.git

cd bml
echo "Installing bml version $(git describe)"
CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:-Release} \
  BLAS_VENDOR=Generic \
  BML_OPENMP=yes \
  BML_TESTING=no \
  ./build.sh install
