export OMPI_MCA_opal_paffinity_alone=0
export OMP_NUM_THREADS=8
#mpirun -np $1 -bynode -x OMP_NUM_THREADS ../../build/gpscf
#../../build/gpmdcov ./input.in | tee  out
/usr/local/bin/mpiexec -n 4 -bycore -x 1  ../../build/gpmdcov ./input.in | tee out
