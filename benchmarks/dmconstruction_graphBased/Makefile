# MAKEFILE FOR twolevelmodel

progress_PATH = ${HOME}/Copa/qmd-progress
bml_PATH = ${HOME}/bml
metis=${HOME}/metis-5.1.0/
MAGMA_PATH=${HOME}/magma

F90      = mpif90 
COMPFLAGS += -fopenmp -I/$(bml_PATH)/install/include -I/$(progress_PATH)/install/include -I/$(metis)/include
COMPFLAGS += -I${MAGMA_PATH}/include -ffree-line-length-none -fopenmp -lpthread -I${CUDA_INCLUDE_DIRS}

LINKFLAGS+= -L$(progress_PATH)/install/lib -lprogress 
LINKFLAGS+= -L$(bml_PATH)/install/lib  -lbml_fortran -lbml
LINKFLAGS+= -L$(metis)/lib -lmetis
LINKFLAGS+= -llapack -lblas
LINKFLAGS+= -fopenmp  -L${MAGMA_PATH}/lib/ -lmagma  -L${CUDA_CUDART_LIBRARY} -lcublas -lcudart -lcusolver



all:
	(make bml)
	(make progress)
	(make gp)

clean:
	rm -f *.o *.mod 
	rm -f main

bml: 
	(cd $(bml_PATH)/build ; make -j; make install)

progress:
	(cd $(progress_PATH)/build ; make -j; make install)

gp:
	$(F90) $(COMPFLAGS) -o  dmconstruction_graphBased dmconstruction_graphBased.F90 $(LINKFLAGS)


