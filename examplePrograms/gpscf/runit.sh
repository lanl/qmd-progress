export OMPI_MCA_opal_paffinity_alone=0
export OMP_NUM_THREADS=$2
mpirun -np $1 -bynode -x OMP_NUM_THREADS ../../build/gpscf
