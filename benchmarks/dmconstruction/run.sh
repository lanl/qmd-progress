#!/bin/bash

#source vars !Sourcing environmental variables

export OMP_NUM_THREADS=12
run=../../build/dmconstruction
nreps=2

device="myCPU" #Architecture name
alg="sp2_alg1" #Name of the algorithm 
tag="prg_sp2_alg1" #Tag for naming output files

for format in Ellpack CSR Ellblock Dense
do
  for system in softmatt semicond metal
  do
    fileout="times_${system}_${alg}_${device}_${format}.dat"
    rm $fileout
    #for i in 1000 2000 4000 6000 8000 10000 12000 14000 16000 18000 20000
    for i in 100 200 400
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

