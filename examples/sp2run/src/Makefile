# MAKEFILE FOR SP2run

progress_PATH = ${HOME}/qmd-progress
bml_PATH = ${HOME}/bml

F90      = gfortran
COMPFLAGS += -fopenmp -I/$(bml_PATH)/install/include -I/$(progress_PATH)/install/include

LINKFLAGS+= -L$(progress_PATH)/install/lib -lprogress
LINKFLAGS+= -L$(bml_PATH)/install/lib  -lbml_fortran -lbml
LINKFLAGS+= -llapack -lblas

all:

	$(F90) $(COMPFLAGS) -o  main main.F90 $(LINKFLAGS)

clean:
	rm -f *.o *.mod 
	rm -f main
	
