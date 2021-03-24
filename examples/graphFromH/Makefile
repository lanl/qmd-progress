# MAKEFILE FOR twolevelmodel

progress_PATH = ${HOME}/Copa/qmd-progress
bml_PATH = ${HOME}/bml
metis=${HOME}/metis-5.1.0/

F90      = mpif90 
COMPFLAGS += -fopenmp -I/$(bml_PATH)/install/include -I/$(progress_PATH)/install/include -I/$(metis)/include

LINKFLAGS+= -L$(progress_PATH)/install/lib -lprogress 
LINKFLAGS+= -L$(bml_PATH)/install/lib  -lbml_fortran -lbml
LINKFLAGS+= -L$(metis)/lib -lmetis
LINKFLAGS+= -llapack -lblas

all:
	(make bml)
	(make progress)
	(make gp)
	(make sp2)

clean:
	rm -f *.o *.mod 
	rm -f main

bml: 
	(cd $(bml_PATH)/build ; make install)

progress:
	(cd $(progress_PATH)/build ; make install)

gp:
	$(F90) $(COMPFLAGS) -o  main main.F90 $(LINKFLAGS)

sp2:
	$(F90) $(COMPFLAGS) -o  main_sp2 main_sp2graph.F90 $(LINKFLAGS)
	$(F90) $(COMPFLAGS) -o  main_gather main_gather.F90 $(LINKFLAGS)

