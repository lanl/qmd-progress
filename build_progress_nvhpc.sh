#!/bin/bash

# Make a copy of this script and modify the METIS and BML library locations

rm -r build
rm -r install_nvhpc

# Set METIS and BML Library locations
METIS_LIB="/usr/projects/icapt/mewall/venado/packages/metis-5.1.0/install"
BML_LIB="/usr/projects/icapt/mewall/venado/packages/qmd-progress/bml/install_nvhpc"
#MPI_LIB="/projects/darwin-nv/rhel9/aarch64/packages/nvhpc/Linux_aarch64/24.3/comm_libs/openmpi/openmpi-3.1.5/lib"

MY_PATH=`pwd`

# Configuring PROGRESS with OpenMP
export CC=$(which mpicc)
export FC=$(which mpifort)
export CXX=$(which mpic++)
#export CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS:="-ffixed-line-length-none"}
#export CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS:="-ef -DCRAY_SDK"}
export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:$BML_LIB/lib/pkgconfig:$BML_LIB/lib64/pkgconfig"
export PROGRESS_OPENMP=${PROGRESS_OPENMP:=yes}
export PROGRESS_MPI=${PROGRESS_MPI:=yes}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/install_nvhpc"}
export PROGRESS_GRAPHLIB=${PROGRESS_GRAPHLIB:=yes}
export PROGRESS_TESTING=${PROGRESS_TESTING:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Release}
export PROGRESS_EXAMPLES=${PROGRESS_EXAMPLES:=yes}
#export PROGRESS_BENCHMARKS=yes
export CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH:=$BML_LIB/:$METIS_LIB/}
#export EXTRA_FCFLAGS=${EXTRA_FCFLAGS:="-ffree-line-length-none  -Wno-pedantic -fallow-argument-mismatch -fdefault-integer-8"}
#export EXTRA_FCFLAGS=${EXTRA_FCFLAGS:="-I${MPI_LIB} -g -ffree-line-length-none  -Wno-pedantic -fallow-argument-mismatch"}
#export EXTRA_LINK_FLAGS=${EXTRA_LINK_FLAGS:="-L${MPI_LIB} -lmpi -lmpi_mpifh -ffree-line-length-none  -Wno-pedantic -fallow-argument-mismatch"}
#export EXTRA_FCFLAGS=${EXTRA_FCFLAGS:="-I${MPI_LIB} -g"}
#export EXTRA_LINK_FLAGS=${EXTRA_LINK_FLAGS:="-L${MPI_LIB} -lmpi -lmpi_mpifh"}
export EXTRA_FCFLAGS=${EXTRA_FCFLAGS:=""}
export EXTRA_LINK_FLAGS=${EXTRA_LINK_FLAGS:=""}
export CMAKE_INCLUDE_PATH=${CMAKE_INCLUDE_PATH:="$METIS_LIB/include"}
export CMAKE_LIBRARY_PATH=${CMAKE_LIBRARY_PATH:="$METIS_LIB/lib"}
#export BLAS_LIBRARIES="-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl"
./build.sh configure

# Make PROGRESS library and examples after running this script:
#   cd build
#   make
