#!/bin/bash

export OMP_NUM_THREADS=4
run=../../build/dmconstruction_bio
nreps=2

device="myCPU" #Architecture name

for format in Ellpack Dense
do
  for system in bio
  do
    fileout="times_${device}_${format}.dat"
    rm $fileout
    for ((i=1; i<=3; i++))
    do
      for ((j=1; j<=$i; j++))
      do
        for ((k=1; k<=$j; k++))
        do
          echo "Format, Replicas:" $format"," $i, $j, $k
          sed 's/ReplicateX=.*/ReplicateX= '$i'/g' input.in.$system > tmp
          sed 's/ReplicateY=.*/ReplicateY= '$j'/g' tmp > input.in
          sed 's/ReplicateZ=.*/ReplicateZ= '$k'/g' input.in > tmp
          sed 's/BMLType=.*/BMLType= '$format'/g' tmp > input.in
          #jsrun -n1 -a1 -g1 -c21 -bpacked:21 $run input.in $nreps > R$i$j$k$device$format$system.out
          $run  input.in $nreps > R$i$j$k$device$format$system.out
          norbs=`grep orbitals R$i$j$k$device$format$system.out | awk 'NF>1{print $NF}'`
          time1=`grep SP2 R$i$j$k$device$format$system.out | awk 'NF>1{print $NF}'`
          time2=`grep Diagonalization R$i$j$k$device$format$system.out | awk 'NF>1{print $NF}'`
          echo $norbs $time1 $time2
          echo $norbs $time1 $time2  >> $fileout
        done
      done
    done
  done
done

