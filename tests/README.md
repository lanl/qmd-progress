# Testing the Progress library {#test}

Testing program for the progress library 
========================================

## To run the tests: ##

Go into the build folder and type:

      make test

To run the tests in verbose mode: 

      make test ARGS="-V"

## To run a single test: ##

To run a test on its own (in build) we just need to type: 

    ctest -R <test_name> --verbose

, where "test_name" is the name of the test we want to run. 
Right now the keywords (test_name) we can pass are the following:

- density : Tests the diagonalization routine to build the density.
- sp2_short : Tests the first version of sp2
- sp2_alg1 :  Algorithm 1 for sp2
- sp2_alg2 :  Algorithm 2 for sp2
- sp2_alg2_ellpack : Algorithm 2 for sp2 with ellpack
- sp2_alg1_seq : See sp2_mod.F90 source file 
- sp2_alg2_seq : See sp2_mod.F90 source file
- deorthogonalize_dense: See nonortho.F90 source file
- orthogonalize_dense: See nonortho.F90 source file
- buildzdiag: See genz_mod.F90 source file

## To add a test: ##

- add the corresponding name of the test in /progress/tests/CMakeLists.txt

- add the corresponding keyword and test in /progress/tests/src/main.F90

- Copy any file that is necessary to run (data) in /progress/tests/tests_data/

- reconfigure and recompile

