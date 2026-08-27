#!/bin/bash

[[ -e build ]] && rm -rfv build
[[ -e install ]] && rm -rfv install

export PROGRESS_MPI=${PROGRESS_MPI:=yes}

if [[ ${PROGRESS_MPI} = "yes" ]]; then
    export CC=${CC:=mpicc}
    export FC=${FC:=mpifort}
    export CXX=${CXX:=mpic++}
else
    export CC=${CC:=gcc}
    export FC=${FC:=gfortran}
    export CXX=${CXX:=g++}
fi

#export FFLAGS="-I$CONDA_PREFIX/lib"
#export FCFLAGS="-I$CONDA_PREFIX/lib"
#export LD_LIBRARY_FLAGS=$LD_LIBRARY_FLAGS:"$CONDA_PREFIX/lib"

export PROGRESS_MPI=no
export BML_DIR=${BML_DIR:=/Users/mewall/packages/gpmd/bml/install}
export PROGRESS_OPENMP=${PROGRESS_OPENMP:=yes}
export PROGRESS_GRAPHLIB=${PROGRESS_GRAPHLIB:=no}
export PROGRESS_TESTING=${PROGRESS_TESTING:=yes}
export PROGRESS_EXAMPLES=${PROGRESS_EXAMPLES:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Release}
export CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}:${BML_DIR}
export VERBOSE_MAKEFILE=${VERBOSE_MAKEFILE:=yes}
export COMMAND=${1:-compile}

./build.sh ${COMMAND}
