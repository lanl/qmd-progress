
gfortran   -I$HOME/bml/install/include  -I$HOME/qmd-progress/install/include -g -o main main.F90 -L$HOME/qmd-progress/install/lib64/ -lprogress -lprg_sp2_tc_fortran -L/usr/local/cuda-11.2/lib64 -lcudart -lcuda -lcublas -lstdc++

