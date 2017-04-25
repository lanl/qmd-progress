set -u
set -x

#export OMPI_MCA_opal_paffinity_alone=0
#export OMPI_MCA_plm_rsh_no_tree_spawn=1
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
: ${MPIRUN:=mpirun}
NODES=${1:-1}
: ${MAP:=core} # Map (node, core)

if ${MPIRUN} --version > /dev/null 2>&1; then
    ${MPIRUN} \
        -np ${NODES} \
        --map-by ${MAP} \
        --hostfile ~/hostfile \
        --mca plm_rsh_no_tree_spawn 1 \
        --mca orte_base_help_aggregate 0 \
        -x OMP_NUM_THREADS \
        -x LD_LIBRARY_PATH \
        ../../build/gpmdcov input.in \
        2>&1 | tee out
else
    ../../build/gpmdcov input.in | tee  out
fi
