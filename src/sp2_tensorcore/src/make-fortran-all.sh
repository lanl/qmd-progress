nvcc -c -g -I$CUDA_INCLUDE linalg_tools.cu -o linalg_tools.o -lcublas -lcudadevrt -lcudart -lstdc++ -lgfortran -lcuda 
nvcc -c -g -I$CUDA_INCLUDE tcore_hp_emulator.cu -o tcore_hp_emulator.o -lcublas -lcudadevrt -lcudart -lstdc++ -lgfortran -lcuda
nvcc -c -g -I$CUDA_INCLUDE prg_sp2_tensorcore.cu -o prg_sp2_tensorcore.o -lcublas -lcudadevrt -lcudart -lstdc++ -lgfortran -lcuda 

gfortran -c -I$CUDA_INCLUDE prg_sp2_tensorcore_mod.F90 -o prg_sp2_tensorcore_mod.o -L/usr/local/cuda-11.0/lib64 -lcublas -lcudadevrt -lcudart -lstdc++ -lgfortran -lcuda
echo "interface compilation complete"

ar -cvq prg_sp2_tensorcore_mod.a prg_sp2_tensorcore_mod.o prg_sp2_tensorcore.o tcore_hp_emulator.o linalg_tools.o 

#gfortran -I$CUDA_INCLUDE -o run test_main.f90 test_mod.F90 test.o -L/usr/local/cuda-11.0/lib64 -lcudart -lcuda -lcublas -lstdc++
gfortran -I$CUDA_INCLUDE -o run main.F90 prg_sp2_tensorcore_mod.a -L/usr/local/cuda-11.0/lib64 -lcudart -lcuda -lcublas -lstdc++
echo "done"

