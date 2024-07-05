#!/bin/bash
export CUDA_VISIBLE_DEVICES=$(echo "${SLURM_LOCALID}%4" |bc)

RUN="$PWD/../../build_hackathon/gpmdk"
#stdbuf -o0 $RUN input_ICH.in |& tee out_ICH
if [ $SLURM_PROCID -eq 0 ]; then
       	stdbuf -o0 $RUN input.in >& out
else
	$RUN input.in >& /dev/null
fi
