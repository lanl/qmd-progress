# MAKEFILE FOR THIS PROGRAM

progress_PATH = ${HOME}/qmd-progress
bml_PATH = ${HOME}/bml

F90      = gfortran
COMPFLAGS += -fopenmp -I/$(bml_PATH)/install/include -I/$(progress_PATH)/install/include

LINKFLAGS+= -L$(progress_PATH)/install/lib64 -lprogress
LINKFLAGS+= -L$(bml_PATH)/install/lib64 -lbml_fortran -lbml -Wall
#LINKFLAGS+= -llapack -lblas

LINKFLAGS+= -lm -Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_lapack95_lp64 -lmkl_gf_lp64 \
            -lmkl_gnu_thread -lmkl_core -lmkl_gnu_thread -lmkl_core -ldl -lpthread -lm

all:

	$(F90) $(COMPFLAGS) -o  main main.F90 $(LINKFLAGS)

clean:
	rm -f *.o *.mod 
	rm -f main
	
