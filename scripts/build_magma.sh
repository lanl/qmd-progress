rm -rf build
mkdir build
cd build
CMAKE_PREFIX_PATH=$CUDA_HOME/../../math_libs/$CRAY_CUDATOOLKIT_VERSION cmake .. -DCMAKE_INSTALL_PREFIX=/usr/projects/icapt/mewall/venado/packages/magma-2.7.2/install -DGPU_TARGET=sm_90 
