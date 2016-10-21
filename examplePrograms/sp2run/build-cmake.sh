#!/bin/bash

# To build a testing program for the progress library.
# Example using this script: 
#         ./buid-cmake.sh ifort compile.
#                        or
#         ./buid-cmake.sh gfortran compile.


TOP_DIR="${PWD}"
BUILD_DIR="${BUILD_DIR:=${PWD}/build}"
COMP=$1
TOR=$2
BML_DIR=/home/christian/bml/install
PROGRESS_DIR=/home/christian/progress/build/src

prepare() {
  read -p $BUILD_DIR" is going to be removed. Are you sure? (press ctrl+c to abort)"
  rm -rf $BUILD_DIR
  mkdir -v -p "${BUILD_DIR}" || exit
}

configure() {
  pushd "${BUILD_DIR}"
  set -v
  if [[ $COMP == "ifort" ]]
  then

    cmake ../ \
      -DCMAKE_Fortran_COMPILER="${FC:=ifort}" \
      -DCMAKE_Fortran_FLAGS="${CMAKE_Fortran_FLAGS:=-openmp \
      -I/$BML_DIR/include -I/$PROGRESS_DIR/}" \
      -DMY_FLAGS="${MY_FLAGS:=-lmkl_lapack95_lp64 -lmkl_intel_lp64 \
      -lmkl_intel_thread -lmkl_core -lmkl_intel_thread -lpthread -openmp -O0 -g\
      -L/$PROGRESS_DIR/ -lprogress \
      -L/$BML_DIR/lib64  -lbml}"

  elif [[ $COMP == "gfortran" ]]
  then

    cmake ../ \
      -DCMAKE_Fortran_COMPILER="gfortran" \
      -DCMAKE_Fortran_FLAGS="${CMAKE_Fortran_FLAGS:=-fopenmp -lgcov --coverage -I/$BML_DIR/include -I/$PROGRESS_DIR/}" \
      -DMY_FLAGS="${MY_FLAGS:=-Wl,--no-as-needed  -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core \
      -lmkl_gnu_thread -ldl -lpthread -lm -fopenmp \
      -L/$PROGRESS_DIR/ -lprogress \
      -L/$BML_DIR/lib/ -lbml}"

  else
    echo "ERROR: no valid compiler selected (ifort or gfortran)"
    exit
  fi

  $* || exit
  set +v
  popd
}

compile() {
  make -C "${BUILD_DIR}" VERBOSE=1
}

prepare
configure

if [[ $TOR == "compile" ]]
then
  compile 
fi  


