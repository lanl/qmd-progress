#!/bin/bash
export CUDA_VISIBLE_DEVICES=$(echo "${SLURM_LOCALID}%4" |bc)

#RUN="nsys profile /usr/projects/icapt/mewall/packages/gpmd/gpmd/build_gpu/gpmdcov"
#RUN="map --profile /usr/projects/icapt/mewall/packages/gpmd/gpmd/build_forge/gpmdcov"
RUN="/usr/projects/icapt/mewall/venado/packages/qmd-progress/examples/gpmdk/build_hackathon/gpmdk"
#stdbuf -o0 $RUN input_ICH.in |& tee out_ICH
if [ $SLURM_PROCID -eq 0 ]; then
#	nsys profile -e NSYS_MPI_STORE_TEAMS_PER_RANK=1 -t mpi,openmp,cuda,nvtx --mpi-impl=mpich --delay=90 $RUN input.in >& out_0
        nsys profile -e NSYS_MPI_STORE_TEAMS_PER_RANK=1 -c cudaProfilerApi  -t mpi,openmp,cuda,nvtx --mpi-impl=mpich $RUN input.in >& out_$SLURM_PROCID
	#nsys profile -t mpi,openmp,cuda,nvtx --mpi-impl=mpich --delay=90 -o report_%q{SLURM_PROCID} $RUN input.in >& out_$SLURM_PROCID
else
	$RUN input.in >& /dev/null
	#nsys profile --start-later -t mpi,openmp,cuda,nvtx --mpi-impl=mpich -o report_%q{SLURM_PROCID} $RUN input.in >& /dev/null
	#nsys profile -t mpi,openmp,cuda,nvtx --mpi-impl=mpich --delay=90 -o report_%q{SLURM_PROCID} $RUN input.in >& /dev/null
fi
