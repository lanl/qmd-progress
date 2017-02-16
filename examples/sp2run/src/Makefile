# MAKEFILE FOR SP2run

progress_PATH = ${HOME}/qmd-progress
bml_PATH = ${HOME}/bml

F90      = gfortran
F90FLAGS += -fopenmp -I/$(bml_PATH)/install/include -I/$(progress_PATH)/install/include

LINKFLAGS+= -L$(progress_PATH)/install/lib -lprogress
LINKFLAGS+= -L$(bml_PATH)/install/lib -lbml
LINKFLAGS+= -llapack -lblas

all:

	$(F90) $(F90FLAGS) -o  main main.F90 $(LINKFLAGS)

clean:
	rm -f *.o *.mod 
	rm -f main
	
