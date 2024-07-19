WORKDIR=/usr/projects/icapt/mewall/venado/packages
rm -rf build_nvhpc
mkdir build_nvhpc
cd build_nvhpc
#cmake .. -DCMAKE_INSTALL_PREFIX=$WORKDIR/magma-2.7.2/install_nvhpc -DGPU_TARGET=sm_90 -DCMAKE_CUDA_ARCHITECTURES=90 -DLAPACK_LIBRARIES="-L/projects/darwin-nv/rhel9/aarch64/packages/nvpl/23.11/lib -lnvpl_lapack_lp64_gomp -lnvpl_blas_lp64_gomp"
CMAKE_PREFIX_PATH=$NVHPC_ROOT/math_libs/lib64 cmake .. -DCMAKE_INSTALL_PREFIX=$WORKDIR/magma-2.7.2/install_nvhpc -DGPU_TARGET=sm_90 -DCMAKE_CUDA_ARCHITECTURES=90 
make -j16 install
cd -
