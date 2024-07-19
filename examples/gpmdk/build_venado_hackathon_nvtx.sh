WORKDIR=/usr/projects/icapt/mewall/venado/packages
rm -rf build_hackathon
mkdir build_hackathon
cd build_hackathon
cmake  -DCMAKE_Fortran_COMPILER="ftn" -DPROGRESS_MPI="yes" -DLIB="no" -DGPMDK_NVTX="yes" \
-DEXTRA_FCFLAGS="-g -O2 -Wall -Wunused -fopenmp -ffree-line-length-none -ffpe-trap=invalid,overflow,zero -lnvToolsExt" \
-DCMAKE_PREFIX_PATH="$WORKDIR/qmd-progress/install_hackathon/;$WORKDIR/qmd-progress/bml/install_hackathon;$WORKDIR/metis-5.1.0/install"  ../src/
#-DEXTRA_FCFLAGS="-g -O2 -DHPCTOOLKIT_PROFILE -Wall -Wunused -fopenmp -ffree-line-length-none -ffpe-trap=invalid,overflow,zero -L$HPCTOOLKIT/lib/hpctoolkit -lhpctoolkit" \
make -j1

