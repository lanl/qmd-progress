Graph-partition Quantum based Molecular Dynamics with Kernel (GPMDK)
=======================================================

About
********
The following research code is based on techniques that were recently published in  `J. Chem. Phys. 158, 074108 (2023) <https://pubs.aip.org/aip/jcp/article/158/7/074108/2877017/Graph-based-quantum-response-theory-and-shadow>`_. 
Briefly, the graph-based electronic structure theory is applied by using the atomic graph and Quantum Molecular Dynamics is performed on a completely memory distributed way. Moreover, this code make use of two LANL developed libraries, namely `PROGRESS <https://qmd-progress.readthedocs.io/en/latest/>`_ and `BML <https://basic-matrix-library.readthedocs.io/en/stable/>`_  which need to be properly installed in order to build the GPMD code. More about this libraries can be found in the following references  `arXiv:2401.13772 <https://arxiv.org/abs/2401.13772>`_ and `J.Supercomput. 74 (11): 6201–19. <https://link.springer.com/article/10.1007/s11227-018-2533-0>`_. 


Getting started
******************
We assume that you have access to the LANL Institutional Computing (IC) as well as access to Darwin. In order to get a Darwin account, please visit the following 
`lik <https://int.lanl.gov/org/ddste/aldsc/ccs/applied-computer-science-ccs-7/access-request-process.shtml>`_.
For Institutional computing and getting access to Chicoma, please, ask the project manager (whoever got 
time to run on IC machines by submitting a proposal) to include you in the specific 
project. Then verify that you are properly included in `here <https://hpcaccounts.lanl.gov/projects/request/xd_g_turq>`_. 

This material can be cloned or downloaded from the PROGRESS github repository (`github link <https://github.com/lanl/qmd-progress.git>`_).
The repository can be cloned as follows:

 .. code-block:: bash
 
    git clone https://github.com/lanl/qmd-progress.git

Look in the examples/gpmdk subdirectory.
    
Further build instructions for gpmdk to be written

**Building required packages METIS and Magma**

Download METIS (latest tested version: 5.1.0):
http://glaros.dtc.umn.edu/gkhome/metis/metis/download
SCP and unpack using: 
    .. code-block:: bash
        scp metis-5.1.0.tar.gz username@ch-fe1:.
        tar -xvf metis-5.1.0.tar.gz

Build  METIS:
// Move into METIS directory
    .. code-block:: bash
        cd metis-5.1.0/
// -cc is the C compiler to use
//-prefix is where to build the metis libraries and can be changed
    .. code-block:: bash
        make config -shared=1 -cc=cc -prefix=~/local
        make install

Download Magma (currently using version 2.5.4) from MAGMA website:
https://icl.utk.edu/magma/downloads/

Secure copy magma tar ball to chicoma:
// Replace “username” with your user name
    .. code-block:: bash
        scp magma-2.5.4.tar.gz username@ch-fe1:.

Untar and unzip the tarball and rename the directory:
    .. code-block:: bash
        tar -xvf magma-2.5.4.tar.gz
        mv magma-2.5.4/ magma

Build magma using the “build_magma.sh” script provided in gpmd/chicomaGPU:
// Start outside of the magma directory then type the following command:
    .. code-block:: bash
        bash build_magma.sh

**Building BML and QMD-PROGRESS**

Clone into the BML repository:
    .. code-block:: bash
        git clone https://github.com/lanl/bml.git
Build BML by running the build_bml.sh script provided in gpmd/chicomaGPU in the directory immediately above bml/ :
    .. code-block:: bash
	bash build_bml.sh
BML requires a BLAS library (basic linear algebra subprograms) which should be available via the intel-mkl module. If you run into issues with BLAS not being found, check to make sure the intel-mkl module is loaded correctly (this should have been done by the setup-envs.sh script). You can check which modules are loaded using the following command:
    .. code-block:: bash
        module list
Clone into the QMD-PROGRESS repository:
    .. code-block:: bash
        git clone https://github.com/lanl/qmd-progress.git
Build QMD-PROGRESS by running the build_progress.sh script provided in gpmd/chicomaGPU in the directory immediately above qmd-progress:
// ** Note: QMD-PROGRESS requires BML, so be sure to build BML first
    .. code-block:: bash
        bash build_progress.sh

**Building GPMDK**

Build GPMDK after building all other dependencies (METIS, Magma, BML, and QMD-PROGRESS) by running the build_cmake.sh script provided with GPMDK :
    .. code-block:: bash
        bash build_cmake.sh

Running GPMDK
*******************

Go into the ``run`` folder. There is an input file ``input.in`` providing most of the variables that are needed by the code. There 
are also several ``pdb`` files with chemical systems that one can use as examples.


