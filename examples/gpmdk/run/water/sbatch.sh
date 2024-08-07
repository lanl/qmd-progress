#!/bin/bash
##SBATCH -p gpu
##SBATCH -p gpu_debug
##SBATCH --reservation=gpu_debug
##SBATCH --qos=debug 
#SBATCH --time 2:00:00
#SBATCH -N 4
##SBATCH -A ichelp_g
##SBATCH -A w23_macroqmd_g

cd ${PWD}
source /usr/projects/icapt/mewall/packages/gpmd/gpmd/setenv_gpu.sh

#export MPICH_ALLREDUCE_NO_SMP=1
export MPICH_SMP_SINGLE_COPY_MODE=NONE
#export MPICH_OPT_THREAD_SYNC=0

OMP_NUM_THREADS=32 srun -n 8 --ntasks-per-node=4 --cpus-per-task=32 bash wrapper.sh   

