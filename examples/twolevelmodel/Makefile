# MAKEFILE FOR twolevelmodel

MAGMA_PATH = ${HOME}/magma-2.5.1

MAGMA_INCLUDE =  -I${MAGMA_PATH}/include
MAGMA_FLAGS= -L${MAGMA_PATH}/lib/ -lmagma  
ESSL_FLAGS= -fopenmp -lpthread -L${OLCF_ESSL_ROOT}/lib64 -lesslsmp


progress_PATH = ${HOME}/qmd-progress
bml_PATH = ${HOME}/bml

F90      = gfortran
COMPFLAGS += -fopenmp -I/$(bml_PATH)/install/include -I/$(progress_PATH)/install/include $(MAGMA_INCLUDE)

LINKFLAGS+= -L$(progress_PATH)/install/lib64 -lprogress 
LINKFLAGS+= -L$(bml_PATH)/install/lib64  -lbml_fortran -lbml
LINKFLAGS+= $(MAGMA_FLAGS) 
LINKFLAGS+= $(ESSL_FLAGS) 

all:

	$(F90) $(COMPFLAGS) -o  main main.F90 $(LINKFLAGS)

clean:
	rm -f *.o *.mod 
	rm -f main
	
