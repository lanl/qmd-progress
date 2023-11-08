#!/bin/bash

#source vars !Sourcing environmental variables

export OMP_NUM_THREADS=8
run=../../build/dmconstruction3d
nreps=2

device="myCPU" #Architecture name
alg="sp2_alg1" #Name of the algorithm 
tag="prg_sp2_alg1" #Tag for naming output files

for format in Ellpack Dense
do
  for system in semicond
  do
    fileout="times_${system}_${alg}_${device}_${format}.dat"
    rm $fileout
    for i in 1024 2000 3456 8192 11664 16000
    do
      echo "Format, System, Size:" $format"," $system"," $i
      sed 's/NOrbs=.*/NOrbs= '$i'/g' input.in.$system > tmp
      sed 's/BMLType=.*/BMLType= '$format'/g' tmp > input.in
      #jsrun -n1 -a1 -g1 -c21 -bpacked:21 ./main input.in $nreps > out$i$device$alg$format$system
      $run  input.in $nreps > out$i$device$alg$format$system
      time=`grep $tag out$i$device$alg$format$system | awk 'NF>1{print $NF}'`
      echo $i $time  >> $fileout
    done
  done
done

