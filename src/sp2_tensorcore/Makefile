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
CUDA_LINK_LIBS= -lcublas lib/linalg_tools.cu lib/tcore_hp_emulator.cu
#
# ##########################################################
#
# ## Project file structure ##
#
# # Source file directory:
SRC_DIR = src
#
# # Object file directory:
OBJ_DIR = bin
#
# # Include header file diretory:
INC_DIR = include
#
# ##########################################################
#
# ## Make variables ##
#
# # Target executable name:
EXE = DNN_run
#
# # Object files:
OBJS = $(OBJ_DIR)/sp2_dnn_tcore.o
# ##########################################################
#
# ## Compile ##
#
# # Link CUDA compiled object files to target executable:
DNN_run : bin/sp2_dnn_tcore.o
	nvcc bin/sp2_dnn_tcore.o -o DNN_run -lcublas lib/linalg_tools.cu lib/tcore_hp_emulator.cu -Iinclude
bin/sp2_dnn_tcore.o : src/sp2_dnn_tcore.cu
	nvcc -c -std=c++11 src/sp2_dnn_tcore.cu -o bin/sp2_dnn_tcore.o


#Clean objects in object directory.
clean:
	$(RM) bin/* *.o
