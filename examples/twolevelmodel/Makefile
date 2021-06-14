# MAKEFILE FOR twolevelmodel

progress_PATH = ${HOME}/qmd-progress
bml_PATH = ${HOME}/bml

F90      = gfortran
COMPFLAGS += -fopenmp -I/$(bml_PATH)/install/include -I/$(progress_PATH)/install/include

LINKFLAGS+= -L$(progress_PATH)/install/lib64 -lprogress
LINKFLAGS+= -L$(bml_PATH)/install/lib64  -lbml_fortran -lbml
LINKFLAGS+= -llapack -lblas  -fopenmp 

all:

	$(F90) $(COMPFLAGS) -o  main main.F90 $(LINKFLAGS)

clean:
	rm -f *.o *.mod 
	rm -f main
	
