#!/bin/bash

export OMP_NUM_THREADS=18

RUN="$HOME/BESGraph/qmd-progress/build/gpmdcov"
#OMP_NUM_THREADS=32 mpirun -np 2 $RUN input.in | tee out  
#OMP_NUM_THREADS=40 mpirun -np 1 $RUN input.in | tee out  
OMP_NUM_THREADS=64 mpirun -np 1 $RUN input.in | tee out  
#OMP_NUM_THREADS=20 $RUN input.in | tee out  
#$RUN input.in | tee out 

