#!/bin/bash

# Make a copy of this script and modify the METIS and BML library locations

rm -r build
#rm -r install_magma_2.7.2
rm -r install_hackathon

# Set METIS and BML Library locations
METIS_LIB="/usr/projects/icapt/mewall/packages/gpmd/metis-5.1.0/"
export BML_LIB=${BML_LIB:="/usr/projects/icapt/mewall/packages/gpmd/qmd-progress/bml/install_hackathon"}
#BML_LIB="/usr/projects/icapt/mewall/packages/gpmd/bml/install_magma_2.7.2"

MY_PATH=`pwd`

# Configuring PROGRESS with OpenMP
export CC=${CC:=cc}
export FC=${FC:=ftn}
export CXX=${CXX:=CC}
export CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS:="-ffixed-line-length-none"}
#export CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS:="-ef -DCRAY_SDK"}
export PKG_CONFIG_PATH="$BML_LIB/lib/pkgconfig:$BML_LIB/lib64/pkgconfig:$PKG_CONFIG_PATH"
export PROGRESS_OPENMP=${PROGRESS_OPENMP:=yes}
export PROGRESS_MPI=${PROGRESS_MPI:=yes}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/install_hackathon"}
export PROGRESS_GRAPHLIB=${PROGRESS_GRAPHLIB:=yes}
export PROGRESS_TESTING=${PROGRESS_TESTING:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Release}
export PROGRESS_EXAMPLES=${PROGRESS_EXAMPLES:=yes}
#export PROGRESS_BENCHMARKS=yes
export CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH:=$METIS_LIB/}
#export EXTRA_FCFLAGS=${EXTRA_FCFLAGS:="-ffree-line-length-none  -Wno-pedantic -fallow-argument-mismatch -fdefault-integer-8"}
export EXTRA_FCFLAGS=${EXTRA_FCFLAGS:="-g -O2 -ffree-line-length-none  -Wno-pedantic -fallow-argument-mismatch"}
export EXTRA_LINK_FLAGS=${EXTRA_LINK_FLAGS:="-g -O2 -ffree-line-length-none  -Wno-pedantic -fallow-argument-mismatch"}
export CMAKE_INCLUDE_PATH=${CMAKE_INCLUDE_PATH:="$METIS_LIB/include"}
export CMAKE_LIBRARY_PATH=${CMAKE_LIBRARY_PATH:="$METIS_LIB/lib"}
#export BLAS_LIBRARIES="-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl"
./build.sh configure

# Make PROGRESS library and examples after running this script:
#   cd build
#   make
