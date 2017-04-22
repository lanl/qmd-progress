set -u
set -x

#export OMPI_MCA_opal_paffinity_alone=0
#export OMPI_MCA_plm_rsh_no_tree_spawn=1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
: ${MPIRUN:=mpirun}
NODES=${1:-1}

if ${MPIRUN} --version > /dev/null 2>&1; then
    ${MPIRUN} \
        -np ${NODES} \
        --map-by node \
        --hostfile hostfile \
        --mca plm_rsh_no_tree_spawn 1 \
        -x OMP_NUM_THREADS \
        -x LD_LIBRARY_PATH \
        ../../build/gpmdcov input.in \
        | tee out
else
    ../../build/gpmdcov input.in | tee  out
fi
