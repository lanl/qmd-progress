#!/bin/bash
#SBATCH -p gpu
##SBATCH -p gpu_debug
##SBATCH --reservation=gpu_debug
##SBATCH --qos=debug 
#SBATCH --time 0:15:00
#SBATCH -N 16
#SBATCH -A w25_proteinqmd_g
##SBATCH -A w23_macroqmd_g

cd ${PWD}
source /usr/projects/icapt/libraries/qmd-progress/scripts/setenv_chicoma_nvhpc.sh

#export MPICH_ALLREDUCE_NO_SMP=1
#export MPICH_SMP_SINGLE_COPY_MODE=NONE
#export MPICH_OPT_THREAD_SYNC=0

OMP_NUM_THREADS=32 srun -n 64 --ntasks-per-node=4 --cpus-per-task=32 bash wrapper.sh   

