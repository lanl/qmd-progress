#gfortran -c -g -I$CUDA_INCLUDE -I$HOME/qmd-progress/install/include prg_sp2_tensorcore_mod.F90 -L$HOME/qmd-progress/install/lib64 -lprogress -lprg_sp2_tc -L$HOME/bml/install/lib64 -lbml_fortran -lbml -L/usr/local/cuda-10.2/lib64 -lcublas -lcudadevrt -lcudart -lstdc++ -lgfortran -lcuda 

gfortran   -I$HOME/bml/install/include  -I$HOME/qmd-progress/install/include -g -o main main.F90 -L$HOME/qmd-progress/install/lib64/ -lprogress -lprg_sp2_tc_fortran -L/usr/local/cuda-10.2/lib64 -lcudart -lcuda -lcublas -lstdc++

