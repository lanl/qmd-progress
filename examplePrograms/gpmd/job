#!/bin/bash
#SBATCH -n 40 #Number of cores 

#SBATCH --ntasks-per-node=40
#
#SBATCH -p ccs6 --time=168:00:00
## OMP_NUM_THREADS controls the number of threads your application use
## This variable cat be set by the following command :

ulimit -s unlimited
module load gcc/5.4.0
module load openmpi/1.10.3-gcc_5.4.0
module load cmake
module load mkl

source /projects/opt/intel/parallel_studio_xe_2016/mkl/bin/mklvars.sh intel64
export OMP_NUM_THREADS=40
export KMP_STACKSIZE=3200M


## Copy files to work directory:
#cp $SUBMITDIR/YourDatafile $SCRATCH

## Mark outfiles for automatic copying to $SUBMITDIR:
#chkfile YourOutputfile

## Run command
./run.sh 
