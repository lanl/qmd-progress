WORKDIR=/usr/projects/icapt/mewall/packages/gpmd
rm -rf build_hackathon
mkdir build_hackathon
cd build_hackathon
cmake  -DCMAKE_Fortran_COMPILER="ftn" -DPROGRESS_MPI="yes" \
-DEXTRA_FCFLAGS="-g -O2 -Wall -Wunused -fopenmp -ffree-line-length-none -ffpe-trap=invalid,overflow,zero" \
-DCMAKE_PREFIX_PATH="$WORKDIR/qmd-progress/install_hackathon/;$WORKDIR/qmd-progress/bml/install_hackathon;$WORKDIR/metis-5.1.0/"  ../src/
#-DEXTRA_FCFLAGS="-g -O2 -DHPCTOOLKIT_PROFILE -Wall -Wunused -fopenmp -ffree-line-length-none -ffpe-trap=invalid,overflow,zero -L$HPCTOOLKIT/lib/hpctoolkit -lhpctoolkit" \
make
