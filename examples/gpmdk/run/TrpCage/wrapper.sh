#!/bin/bash
export CUDA_VISIBLE_DEVICES=$(echo "${SLURM_LOCALID}%4" |bc)

#RUN="nsys profile /usr/projects/icapt/mewall/packages/gpmd/gpmd/build_gpu/gpmdcov"
#RUN="map --profile /usr/projects/icapt/mewall/packages/gpmd/gpmd/build_forge/gpmdcov"
RUN="/usr/projects/icapt/mewall/packages/gpmd/qmd-progress/examples/gpmdk/build_hackathon/gpmdk"
#stdbuf -o0 $RUN input_ICH.in |& tee out_ICH
if [ $SLURM_PROCID -eq 0 ]; then
	stdbuf -o0 $RUN input.in >& out_0
else
	$RUN input.in >& /dev/null
fi
