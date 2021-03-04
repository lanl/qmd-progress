export OMP_NUM_THREADS=20
mpirun -np 2 $HOME/BESGraph/qmd-progress/build/gpmdcov input.in
#mpiexec -np 1  $HOME/BESGraph/qmd-progress/build/gpmdcov input.in
#$HOME/BESGraph/qmd-progress/build/gpmdcov input.in
