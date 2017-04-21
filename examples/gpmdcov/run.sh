set -u

export OMPI_MCA_opal_paffinity_alone=0
export OMP_NUM_THREADS=${OMP_NUM_THREADS:=8}
: ${MPIRUN:=mpirun}
NODES=${1:-1}

if ${MPIRUN} --version > /dev/null 2>&1; then
    ${MPIRUN} \
        -np ${NODES} \
        --map-by node \
        --host 10.209.225.48 \
        -x OMP_NUM_THREADS \
        ../../build/gpmdcov ./input.in \
        | tee out
else
    ../../build/gpmdcov ./input.in | tee  out
fi
