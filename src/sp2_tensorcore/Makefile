## USER SPECIFIC DIRECTORIES ##
#
# ##CUDA directory:
#CUDA_ROOT_DIR=. #/usr/local/cuda
#
# ##########################################################
#
# ## CC COMPILER OPTIONS ##
#
# ## CC compiler options:
# CC=gcc
# CC_FLAGS=
# CC_LIBS=
#
# ##########################################################
#
# ## NVCC COMPILER OPTIONS ##
#
# # NVCC compiler options:
NVCC=nvcc
#NVCC_FLAGS= 
#NVCC_LIBS=

# # CUDA library directory:
CUDA_LIB_DIR= -Llib
# # CUDA include directory:
CUDA_INC_DIR= -Iinclude
# # CUDA linking libraries:
CUDA_LINK_LIBS= -lcublas src/prg_sp2_tensorcore.cu src/linalg_tools.cu src/tcore_hp_emulator.cu
#
# ##########################################################
#
# ## Project file structure ##
#
# # Source file directory:
SRC_DIR = src
#
# # Object file directory:
OBJ_DIR = ./
#
# # Include header file diretory:
INC_DIR = include
#
# ##########################################################
#
# ## Make variables ##
#
# # Target executable name:
EXE = main
#
# # Object files:
OBJS = main.o 
OBJS =+ prg_sp2_tensorcore.o
# ##########################################################
#
# ## Compile ##
#
# # Link CUDA compiled object files to target executable:
#prg_sp2_tensorcore.o : src/prg_sp2_tensorcore.cu
#	nvcc -c -std=c++11 src/prg_sp2_tensorcore.cu -o prg_sp2_tensorcore.o -Iinclude -lcublas
main : main.o
	nvcc main.o -o main src/prg_sp2_tensorcore.cu src/linalg_tools.cu src/tcore_hp_emulator.cu -Iinclude -lcublas
main.o : src/main.cu
	nvcc -c -std=c++11 src/main.cu -o main.o -Iinclude


#Clean objects in object directory.
clean:
	$(RM) bin/* *.o
