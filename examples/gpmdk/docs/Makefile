# Makefile for the GPMD documentation

DIR=$(shell pwd)

all:
	(make doxy)	

clean:
	rm -f *.backup

doxy:
	(doxygen ./source-doxygen/gpmd_doxyfile)	
	(cd ./build-doxygen/latex ; make ; cp refman.pdf ../../gpmd-doxygen.pdf)	

sphinx:
	(cp ../README.rst ./source-sphinx/README.rst)
	(cd ./source-sphinx/; make latex)
	(cd ./build-sphinx/latex/; make; cp gpmd.pdf ../../gpmd-sphinx.pdf)

