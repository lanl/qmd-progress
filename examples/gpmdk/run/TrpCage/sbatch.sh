#!/bin/bash
##SBATCH -p gpu
##SBATCH -p gpu_debug
##SBATCH --reservation=gpu_debug
##SBATCH --qos=debug 
#SBATCH --time 2:00:00
#SBATCH -N 64
##SBATCH -A ichelp_g
##SBATCH -A w23_macroqmd_g

cd ${PWD}
source /usr/projects/icapt/mewall/venado/packages/qmd-progress/scripts/setenv_venado.sh

#export MPICH_ALLREDUCE_NO_SMP=1
export MPICH_SMP_SINGLE_COPY_MODE=NONE
#export MPICH_OPT_THREAD_SYNC=0

OMP_NUM_THREADS=72 srun -n 64 --ntasks-per-node=4 --cpus-per-task=72 bash wrapper.sh   

