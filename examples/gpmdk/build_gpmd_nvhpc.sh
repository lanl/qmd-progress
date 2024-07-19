WORKDIR=/usr/projects/icapt/mewall/venado/packages
rm -rf build_nvhpc
mkdir build_nvhpc
cd build_nvhpc
cmake  -DCMAKE_Fortran_COMPILER="mpifort" -DPROGRESS_MPI="yes" -DGPMDK_NVTX="yes" \
-DEXTRA_FCFLAGS="-g -mp -L${NVHPC_ROOT}/cuda/lib64 -lnvToolsExt" \
-DCMAKE_PREFIX_PATH="$WORKDIR/qmd-progress/install_nvhpc/;$WORKDIR/qmd-progress/bml/install_nvhpc;$WORKDIR/metis-5.1.0/install"  ../src/
make

