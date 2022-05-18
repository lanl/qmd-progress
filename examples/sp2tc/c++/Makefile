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
CUDA_LIB_DIR= -L$(CUDA_ROOT)/lib64
# # CUDA include directory:
CUDA_INC_DIR= -I$(HOME)/qmd-progress-new/install/include
#CUDA_INC_DIR+= -I$(HOME)/bml/install/include

# # CUDA linking libraries:
CUDA_LINK_LIBS= -lcublas -lcusolver -L$(HOME)/qmd-progress-new/install/lib64 -lprogress -lprg_sp2_tc 
#CUDA_LINK_LIBS+= -L$(HOME)/bml/install/lib64 -lbml

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

main : main.o
	nvcc main.o -o main $(CUDA_LINK_LIBS)
main.o : main.cu
	nvcc -c -std=c++11 main.cu $(CUDA_INC_DIR) -o main.o 


#Clean objects in object directory.
clean:
	$(RM) bin/* *.o
