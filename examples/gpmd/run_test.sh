#!/bin/bash

# Script to run this test this program (GPMD).
# Written by C. F. A. Negre Jul. 2016 Los Alamos Nat. Lab.

set -e                                          # This will exit the script if there is any error
set -u
set -x
MY_PATH=`pwd`                                   # Capturing the local path of the folder where we are running.

#RUN="mpirun -np 1 ../../build/gpmd"  #GPMD program
RUN="../../build/gpmd"  #GPMD program
cd ../latteTBparams/
PARAMS="$( pwd | sed 's_/_\\/_g' )"

cd $MY_PATH

for name in ch4 sucrose ; do

  INFILE="input_"$name".in"
  REF="energy_"$name".out"
  COORDS="coords_"$name".dat"
  STR="PAR"
  STRR='"'$PARAMS'"'

  cp  ./tests/$INFILE .
  sed -e s/"PAR"/$STRR/g $INFILE  >  input_tmp.in

  cp  ./tests/$REF .
  cp  ./tests/$COORDS coords_tmp.dat

  echo -e "\nTesting for "$name" \n"
  grep -A 20 DESCRIPTION $INFILE


  time $RUN input_tmp.in > out
  echo ""
  grep -e "Energy Total \[eV\] =" out | sed -e 's/Energy Total \[eV\]/ /g' | awk 'NF>1{print $2}' > energy.out
#  python get_energy.py out > energy.out
  python test-energy.py --reference $REF --current energy.out --reltol 0.000001
  rm $INFILE
  rm $REF
  rm coords_tmp.dat
  rm input_tmp.in
  rm out

done

rm Test_*

echo -e "\nEnd of run and test"
