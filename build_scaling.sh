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

export PROGRESS_OPENMP=${PROGRESS_OPENMP:=yes}
export PROGRESS_GRAPHLIB=${PROGRESS_GRAPHLIB:=yes}
export PROGRESS_TESTING=${PROGRESS_TESTING:=yes}
export PROGRESS_EXAMPLES=${PROGRESS_EXAMPLES:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Release}
export VERBOSE_MAKEFILE=${VERBOSE_MAKEFILE:=yes}
export COMMAND=${1:-compile}

./build.sh ${COMMAND}
