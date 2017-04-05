---
title: progress
---

# A library for quantum chemistry solvers.

PROGRESS: Parallel, Rapid _O(N)_ and Graph-based Recursive Electronic
Structure Solver. **LA-CC-16-068**

- This library is focused on the development of general solvers that are
  commonly used in _quantum chemistry packages_.

- This library has to be compiled with the [_Basic Matrix Library_
  (BML)](https://qmmd.github.io/bml/).

## Authors

(in alphabetical order)

- Anders M. N. Niklasson <<amn@lanl.gov>>
- Christian F. A. Negre <<cnegre@lanl.gov>>
- Marc J. Cawkwell <<cawkwell@lanl.gov>>
- Nicolas Bock <<nicolasbock@gmail.com>>
- Susan M. Mniszewski <<smm@lanl.gov>>
- Michael E. Wall <<mewall@lanl.gov>>

## Build Dependencies

- `>=OpenMP-3.1`
- `>=metis-5.0` if building with `PROGRESS_GRAPHLIB`

(On some distributions, metis is available as a package. Make sure you install
the `-dev` package. For example, Ubuntu requires `libmetis-dev`.)

## How to build

    $ CMAKE_PREFIX_PATH=<BML install path> ./build.sh

## How to install

    $ cd build
    $ sudo make install

To specify the Intel Fortran compiler:

    $ FC=ifort PKG_CONFIG_PATH=<BML install path>/lib/pkgconfig ./build.sh

To build with the gfortran compiler and OpenMP:

    $ CC=gcc FC=gfortran \
        CMAKE_BUILD_TYPE=Release \
        PROGRESS_OPENMP=yes \
        CMAKE_PREFIX_PATH=<BML install path> \
        CMAKE_INSTALL_PREFIX=<PROGRESS install path> \
        ./build.sh configure

To build with OpenMP, MPI and testing enabled:

    $ CC=mpicc FC=mpif90 \
        CMAKE_BUILD_TYPE=Release \
        PROGRESS_OPENMP=yes \
        PROGRESS_MPI=yes \
        PROGRESS_TESTING=yes \
        CMAKE_PREFIX_PATH=<BML install path> \
        CMAKE_INSTALL_PREFIX=<PROGRESS install path> \
        ./build.sh configure

To build with OpenMP, MPI, testing enabled and example programs built:

	$ CC=mpicc FC=mpif90 \
	    CMAKE_BUILD_TYPE=Release \
	    PROGRESS_OPENMP=yes \
	    PROGRESS_MPI=yes \
	    PROGRESS_TESTING=yes \
	    PROGRESS_EXAMPLES=yes \
	    CMAKE_PREFIX_PATH=<BML install path> \
	    CMAKE_INSTALL_PREFIX=<PROGRESS install path> \
	    ./build.sh configure

To build with OpenMP and MPI and testing enabled and example programs built
and the METIS graph partitioning library:

	$ CC=mpicc FC=mpif90 \
	    CMAKE_BUILD_TYPE=Release \
	    PROGRESS_OPENMP=yes \
	    PROGRESS_MPI=yes \
	    PROGRESS_GRAPHLIB=yes \
	    PROGRESS_TESTING=yes \
	    PROGRESS_EXAMPLES=yes \
	    CMAKE_PREFIX_PATH=<BML install path> \
	    CMAKE_INSTALL_PREFIX=<PROGRESS install path> \
	    ./build.sh configure
