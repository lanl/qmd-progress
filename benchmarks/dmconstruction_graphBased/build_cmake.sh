mkdir build
cd build

METIS_LIB="$HOME/metis-5.1.0/lib"
METIS_INCLUDE="$HOME/metis-5.1.0/include"
PROGRESS="$HOME/qmd-progress/install/"
BML="$HOME/bml/install"

cmake  -DCMAKE_PREFIX_PATH="$BML;$PROGRESS;$METIS_LIB;$METIS_INCLUDE" ../
make  
