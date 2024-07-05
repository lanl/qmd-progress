#!/bin/bash


#RUN="$PWD/../build_hackathon/gpmdk"
#OMP_NUM_THREADS=32 mpiexec -np 3 $RUN input.in | tee out  
#OMP_NUM_THREADS=36 mpirun -np 1 $RUN input.in | tee out
#OMP_NUM_THREADS=36 mpiexec -np 1 $RUN input.in | tee out

#OMP_NUM_THREADS=36 mpirun -np 3 --map-by socket:pe=18 --bind-to core:overload-allowed  $RUN input.in | tee out
OMP_NUM_THREADS=32 srun -n 4 -c 32 bash wrapper.sh   
#OMP_NUM_THREADS=64 $RUN input.in | tee out  
#OMP_NUM_THREADS=36 $RUN input.in | tee out 
#$RUN input.in | tee out 

