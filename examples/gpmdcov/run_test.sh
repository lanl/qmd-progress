#!/bin/bash

export OMP_NUM_THREADS=20

RUN="$HOME/BESGraph/qmd-progress/build/gpmdcov"
rm *.tmp
for name in gpmdcov_Init gpmdcov_Part gpmdcov_InitParts gpmdcov_FirstCharges \
gpmdcov_DM_Min gpmdcov_PrepareMD gpmdcov_MDloop; do
  sed -e s/"gpmdcov_"/$name/g ./tests/input_test.in > input.tmp
  mpirun -np 2 $RUN input_tmp.in &> /dev/null
  error=`diff log.gpmdcov ./tests/ref.$name`
  if [ "$error" != "" ]
  then
  echo "ERROR at "$name
  exit -1
  else
  echo "Test for "$name" ... OK"
  fi
done

