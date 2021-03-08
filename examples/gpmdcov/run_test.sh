#!/bin/bash

export OMP_NUM_THREADS=20

RUN="$HOME/BESGraph/qmd-progress/build/gpmdcov"
#rm *.tmp
for name in DiagEf BlockPart; do
  cp ./tests/input_$name'.in' input.tmp
  mpirun -np 2 $RUN input.tmp &> out.tmp 
  ./get_energy.py out.tmp &> energy.tmp
  result=`./test-energy.py --reference ./tests/ref.$name'.out' --current energy.tmp `
  if [ "$result" != "Error" ]
  then
  echo "ERROR at "$name
  exit -1
  else
  echo "Test for "$name" ... OK"
  fi
done

