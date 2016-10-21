# MAKEFILE FOR SP2run

BML_DIR=/home/christian/bml/install
PROGRESS_DIR=/home/christian/progress/build/src

#F90      = ifort  
F90      = gfortran
#F90FLAGS += -fast -openmp 
F90FLAGS += -fopenmp -I/$(BML_DIR)/include -I/$(PROGRESS_DIR)/

LIBS+= -Wl,--no-as-needed  -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm -fopenmp 
LIBS+= -L/$(PROGRESS_DIR)/ -lprogress
LIBS+= -L/$(BML_DIR)/lib/ -lbml


# how to produce objects from source files

%.o:%.F90
	$(F90) $(F90FLAGS) -c $<

# General vars

TARGET = ../build/main

all: $(TARGET)

PROGRAM = main.o 

#MODULES = add_all_module_here.o

OBJS = $(PROGRAM) $(MODULES)

# dependencies

main.o : 
#module1.o : module2.o	

$(MODULES) :
$(PROGRAM) : $(MODULES)

# rules

$(TARGET):	$(MODULES) $(PROGRAM) 
	$(F90) -o  $(TARGET) $(LINKFLAGS) $(OBJS) $(LIBS) $(LIBSP)  

clean:
	rm -f *.o *.mod 
	rm -f ../build/main
	
