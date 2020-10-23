gfortran -c -I$CUDA_INCLUDE -I$HOME/qmd-progress/install/include -I$HOME/qmd-progress/install/include prg_sp2_tensorcore_mod.F90 -L/$HOME/qmd-progress/install/lib64 -lprg_sp2_tc -L/usr/local/cuda-10.2/lib64 -lcublas -L$HOME/qmd-progress/install/lib64 -lprogress -L$HOME/bml/install/lib64 -lbml_fortran -lbml -lcudart -lstdc++

nvcc   -I$HOME/bml/install/include  -I$HOME/qmd-progress/install/include -o main main.cu  -L/$HOME/qmd-progress/install/lib64/ -lprg_sp2_tc -L/usr/local/cuda-10.2/lib64 -lcublas   -L$HOME/qmd-progress/install/lib64 -lprogress -L$HOME/bml/install/lib64 -lbml_fortran -lbml -lcudart -lstdc++

