#!/bin/bash

rm -r build
rm -r install

module load cmake
module load gcc/9.3.0
module load spectrum-mpi

# Set BML Library location
BML_LIB="/gpfs/alpine/mat190/scratch/jeanluc/GIT/bml/install/lib64"

MY_PATH=`pwd`

#get jsrun with full path
JSRUN=$(which jsrun)
echo ${JSRUN}

# Configuring PROGRESS with OpenMP
export CC=${CC:=mpicc}
export FC=${FC:=mpifort}

export PROGRESS_OPENMP=${PROGRESS_OPENMP:=yes}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/install"}
export PROGRESS_GRAPHLIB=${PROGRESS_GRAPHLIB:=no}
export PROGRESS_TESTING=${PROGRESS_TESTING:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Release}
export PROGRESS_EXAMPLES=${PROGRESS_EXAMPLES:=yes}
export PROGRESS_MPI=${PROGRESS_MPI:=yes}
export PROGRESS_BENCHMARKS=${PROGRESS_BENCHMARKS:=yes}

export PROGRESS_NONMPI_PRECOMMAND=${PROGRESS_NONMPI_PRECOMMAND:=${JSRUN}}
export PROGRESS_NONMPI_PRECOMMAND_ARGS=${PROGRESS_NONMPI_PRECOMMAND_ARGS:="-n1;-a1;-g1;-c7;--smpiargs=off"}

PKG_CONFIG_PATH=$BML_LIB/pkgconfig

./build.sh configure

# Make PROGRESS library and examples after running this script:
#   cd build
#   make
