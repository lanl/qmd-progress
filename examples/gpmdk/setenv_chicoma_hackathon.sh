module switch PrgEnv-cray PrgEnv-gnu
module swap gcc gcc/11.2.0
export PKG_CONFIG_PATH=/opt/cray/xpmem/2.5.2-2.4_3.45__gd0f7936.shasta/lib64/pkgconfig:$PKG_CONFIG_PATH
module load cmake
module load cudatoolkit
module load intel-mkl
