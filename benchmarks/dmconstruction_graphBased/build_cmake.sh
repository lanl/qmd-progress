mkdir build
cd build

METIS_LIB="$HOME/metis-5.1.0/lib"
METIS_INCLUDE="$HOME/metis-5.1.0/include"
PROGRESS="$HOME/qmd-progress/install/"
BML="$HOME/bml/install"
BML_INCLUDE="$HOME/bml/install/include"
FC="/usr/bin/mpif90"
cmake  -DCMAKE_PREFIX_PATH="$BML;$BML_INCLUDE;$PROGRESS;$METIS_LIB;$METIS_INCLUDE" ../
make  
