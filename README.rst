This website is intended to provide some guidance on how to get and install the
**PROGRESS** library. LA-UR number **LA-UR-17-27372**

.. list-table::
  :header-rows: 1

  * - Issues
    - Pull Requests
    - CI
  * - .. image:: https://img.shields.io/github/issues/lanl/qmd-progress.svg
        :alt: GitHub issues
        :target: https://github.com/lanl/qmd-progress/issues
    - .. image:: https://img.shields.io/github/issues-pr/lanl/qmd-progress.svg
        :alt: GitHub pull requests
        :target: https://github.com/lanl/qmd-progress/pulls
    - .. image:: https://github.com/lanl/qmd-progress/workflows/CI/badge.svg
        :alt: GitHub Actions
        :target: https://github.com/lanl/qmd-progress/actions

A library for quantum chemistry solvers
=======================================

**PROGRESS**: Parallel, Rapid **O(N)** and Graph-based Recursive Electronic
Structure Solver. **LA-CC-16-068**

- This library is focused on the development of general solvers that are
  commonly used in **quantum chemistry packages**.

- This library has to be compiled with the `Basic Matrix Library (BML)
  <https://basic-matrix-library.readthedocs.io/en/latest/>`_.

- Our webpage can be found at https://qmd-progress.readthedocs.io/

Authors
=======

(in alphabetical order)

- Anders M. N. Niklasson <amn@lanl.gov>
- Christian F. A. Negre <cnegre@lanl.gov>
- Marc J. Cawkwell <cawkwell@lanl.gov>
- Nicolas Bock <nicolasbock@gmail.com>
- Susan M. Mniszewski <smm@lanl.gov>
- Michael E. Wall <mewall@lanl.gov>

Contributors
============

- Jesse Grindstaff <grindstaff@lanl.gov>
- Alicia Welden <welden@umich.edu>
- Nestor Aguirre <nfaguirrec@lanl.gov>
- Jean-Luc Fattebert <fattebertj@ornl.gov>

Build Dependencies
==================

- **>=OpenMP-3.1**
- **>=metis-5.0** if building with **PROGRESS_GRAPHLIB**

Note that on some distributions, metis is available as a package. Make sure you
install the **-dev** package. For example, Ubuntu requires **libmetis-dev**.

Testing in our CI container
===========================

We are switching our CI tests from Travis-CI to GitHub Actions because Travis-CI
is `limiting the number of builds for open source projects
<https://blog.travis-ci.com/2020-11-02-travis-ci-new-billing>`_. Our workflow
uses a
`custom Docker image <https://hub.docker.com/r/nicolasbock/qmd-progress>`_ which
comes with the necessary compiler tool chain and a pre-installed `bml` library
to build and test the **qmd-progress** library. Using **docker** is a convenient
and quick way to develop, build, and test the **qmd-progress** library.

.. code-block:: console

    ./scripts/run-local-docker-container.sh

Inside the container:

.. code-block:: console

    ./build.sh compile

Alternatively, you can run one of the CI tests by executing e.g.

.. code-block:: console

    ./scripts/ci-with-graphlib-debug.sh

Build and Install Instructions
==============================

How to build
------------

.. code-block:: console

    CMAKE_PREFIX_PATH=<BML install path> ./build.sh

How to install
--------------

.. code-block:: console

    cd build
    sudo make install

To specify the Intel Fortran compiler:

.. code-block:: console

    FC=ifort PKG_CONFIG_PATH=<BML install path>/lib/pkgconfig ./build.sh

To build with the gfortran compiler and OpenMP:

.. code-block:: console

    CC=gcc FC=gfortran \
        CMAKE_BUILD_TYPE=Release \
        PROGRESS_OPENMP=yes \
        CMAKE_PREFIX_PATH=<BML install path> \
        CMAKE_INSTALL_PREFIX=<PROGRESS install path> \
        ./build.sh configure

To build with OpenMP, MPI and testing enabled:

.. code-block:: console

    CC=mpicc FC=mpif90 \
        CMAKE_BUILD_TYPE=Release \
        PROGRESS_OPENMP=yes \
        PROGRESS_MPI=yes \
        PROGRESS_TESTING=yes \
        CMAKE_PREFIX_PATH=<BML install path> \
        CMAKE_INSTALL_PREFIX=<PROGRESS install path> \
        ./build.sh configure

To build with OpenMP, MPI, testing enabled and example programs built:

.. code-block:: console

    CC=mpicc FC=mpif90 \
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

.. code-block:: console

    CC=mpicc FC=mpif90 \
	    CMAKE_BUILD_TYPE=Release \
	    PROGRESS_OPENMP=yes \
	    PROGRESS_MPI=yes \
	    PROGRESS_GRAPHLIB=yes \
	    PROGRESS_TESTING=yes \
	    PROGRESS_EXAMPLES=yes \
	    CMAKE_PREFIX_PATH=<BML install path> \
	    CMAKE_INSTALL_PREFIX=<PROGRESS install path> \
	    ./build.sh configure

Citing
======

.. code-block:: bibtex

    @misc{2016progress,
        title={\textrm{PROGRESS} Version 1.0},
        author={Niklasson, Anders M. and
                Mniszewski, Susan M and
                Negre, Christian F. A. and
                Wall, Michael E. and
                Cawkwell, Marc J., and
                Nicolas Bock},
        year={2016},
        url = {https://github.com/lanl/qmd-progress},
        institution={Los Alamos National Laboratory (LANL), Los Alamos, NM (United States)}
    }

Support acknowledges
====================

This development is currently supported by the Exascale Computing Project
(17-SC-20-SC), a collaborative effort of two U.S. Department of Energy
organizations (Office of Science and the National Nuclear Security
Administration) responsible for the planning and preparation of a capable
exascale ecosystem, including software, applications, hardware, advanced system
engineering, and early testbed platforms, in support of the nation’s exascale
computing imperative.

Basic Energy Sciences (LANL2014E8AN) and the Laboratory Directed Research and
Development Program of Los Alamos National Laboratory. To tests these
developments we used resources provided by the Los Alamos National Laboratory
Institutional Computing Program, which is supported by the U.S. Department of
Energy National Nuclear Security Administration
