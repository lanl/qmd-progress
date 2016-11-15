export OMPI_MCA_opal_paffinity_alone=0
export OMP_NUM_THREADS=$2
#mpirun -np $1 --map-by node -x OMP_NUM_THREADS ../../build/gpmd_dist ./input.in
mpirun -np $1 --map-by node -x OMP_NUM_THREADS ../../build/gpmd_dist ./input_300.in
#mpirun -np $1 --map-by node -x OMP_NUM_THREADS ../../build/gpmd_dist ./input_3000.in
#mpirun -np $1 --map-by node -x OMP_NUM_THREADS ../../build/gpmd_dist ./input_6000.in
#mpirun -np $1 --map-by node -x OMP_NUM_THREADS ../../build/gpmd_dist ./input_20000.in
