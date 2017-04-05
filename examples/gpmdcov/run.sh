set -u

export OMPI_MCA_opal_paffinity_alone=0
export OMP_NUM_THREADS=${OMP_NUM_THREADS:=8}
: ${MPIRUN:=mpirun}
NODES=${1:-1}

if [[ -x ${MPIRUN} ]]; then
    ${MPIRUN} \
        -np ${NODES} \
        --map-by node \
        -x OMP_NUM_THREADS \
        ../../build/gpmdcov ./input.in \
        | tee out
else
    ../../build/gpmdcov ./input.in | tee  out
fi
#/usr/local/bin/mpiexec -n 4 -bycore -x 1  ../../build/gpmdcov ./input.in | tee out
