export OMPI_MCA_opal_paffinity_alone=0
export OMP_NUM_THREADS=$1
../../build/gpmd_dist ./input.in
