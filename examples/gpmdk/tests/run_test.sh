#!/bin/bash

export OMP_NUM_THREADS=20
source ../vars
RUN="$HOME/gpmd/build/gpmdcov"
#RUN="/home/cnegre/qmd-progress/build/gpmdcov"
#rm *.tmp

declare -A Ref 
Ref["KerFullBuildAlways"]=16
Ref["KerFullBuildOnce"]=17
Ref["KerFullBuildOnceBlock"]=20
Ref["KerFullBuildAlwaysBlock"]=20
Ref["KerFullBuildOnceMetis"]=23
Ref["KerFullBuildAlwaysMetis"]=23
Ref["KerByPartsBuildOnceBlockUpdate"]=21
Ref["KerByPartsBuildOnceMetisUpdate"]=20
Ref["KerByPartsBuildOnceBlock"]=18
Ref["KerByPartsBuildAlwaysMetis"]=18

if [ 1 == 1 ]
then 
for name in KerFullBuildAlways KerFullBuildOnce KerFullBuildOnceBlock \
	    KerFullBuildAlwaysBlock KerFullBuildOnceMetis KerFullBuildAlwaysMetis \
	    KerByPartsBuildOnceBlockUpdate KerByPartsBuildOnceMetisUpdate ; do

  cp ./input_$name'.in' input.tmp
  if [ "$1" = "mpi" ]
  then
   echo "Testing with MPI ..."
   #mpirun -np 2 $RUN input.tmp &> out.tmp
   mpirun -np 2 ../build/gpmdcov input.tmp &> out.tmp
  else
   $RUN input.tmp &> out.tmp
  fi
  nscf=`grep "SCF converged" out.tmp | awk 'NF>1{print $(NF-2)}'`
  echo $nscf
  if [ "$nscf" -le "${Ref[$name]}" ]  
  then
  echo "Test for "$name" ... OK"
  else
  echo "ERROR at "$name 
  exit -1
  fi

done
fi


for name in DiagEf BlockPart KerUpdate; do
  cp ./input_$name'.in' input.tmp
  if [ "$1" = "mpi" ]
  then
   mpirun -np 2 $RUN input.tmp &> out.tmp 
  else
   $RUN input.tmp &> out.tmp 
  fi
  ./get_energy.py out.tmp &> energy.tmp
  result=`./test-energy.py --reference ./ref.$name --current energy.tmp `
  if [ "$result" != "Error" ]
  then
  echo "ERROR at "$name
  exit -1
  else
  echo "Test for "$name" ... OK"
  fi
done

