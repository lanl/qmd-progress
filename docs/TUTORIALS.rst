
Tutorials
===========

In this section, we will provide some examples of how to effectively use the PROGRESS library. To follow along, it is important to be already familiar with the `BML <https://basic-matrix-library.readthedocs.io/en/stable/>`_ library and have it successfully compiled and installed. All the PROGRESS libraries API calls details can be found at `PROGRESS <_static/doxy/group__PROGRESS.html>`_. Here we will not provide details on the API calls. This documentation has been approved for unlimited release with **LA-UR-23-28151**.

Sample Code
-----------

To begin, let's create a sample code that can be easily compiled by utilizing the necessary "flags" and environment variables from both PROGRESS and BML libraries. In this example, we assume that these libraries are located in the ``$HOME`` directory.

1. Open your terminal and navigate to your home directory, and create a new directory called  "exampleCode":: 
   
    cd ~
    mkdir exampleCode

2. Change into the "exampleCode" directory, and create a folder called "src" to store our source code files::
   
    cd exampleCode
    mkdir src

3. Now you are ready to start writing a sample code within the "src" directory. This structure will help keep the code organized and maintainable throughout the development process.

Here is an example of a sample code file (e.g., ``main.F90``) that you can create inside the "src" directory:

.. code-block:: fortran
   
   program main

    use bml
    use prg_sp2_mod
    use myMod

    call mySub()
  
  end program main

Because most certainly the code will have modules. Let's also create a module (e.g., ``myMod.F90``) for our code with an example 
subroutine:

.. code-block:: fortran
   
   module myMod
   contains

   subroutine mySub 

    write(*,*)"Hello from Sub!"

   end subroutine mySub   

  end module myMod
  
4. Next, we will write a simple ``CMakeLists.txt`` inside ``/src`` which will allow us to automatically 
compile our code using all the FORTRAN flags that were used for BML and PROGRESS. We will also add directives to 
link to both BML and PROGRESS libraries. The ``CMakeLists.txt`` will look like: 

.. code-block:: cmake

  cmake_minimum_required(VERSION 3.10.0)
  project(main VERSION 1.0.0 LANGUAGES Fortran)

  set(dir ${CMAKE_CURRENT_SOURCE_DIR}/build/)
  set(CMAKE_BUILD_DIRECTORY ${dir})

  include(FindPkgConfig)

  find_package(PROGRESS CONFIG QUIET)
  pkg_check_modules(PROGRESS REQUIRED progress)
  message(STATUS "Found progress: ${PROGRESS_LDFLAGS}")
  list(APPEND LINK_LIBRARIES ${PROGRESS_LDFLAGS})

  find_package(BML CONFIG QUIET)
  pkg_check_modules(BML REQUIRED bml)
  list(APPEND LINK_LIBRARIES BML::bml)
  list(APPEND LINK_LIBRARIES ${BML_LDFLAGS})
  message(STATUS "Found bml: ${BML_LDFLAGS}")

  if(PROGRESS_MPI)
    message(STATUS "Will build with MPI")
    add_definitions(-DDO_MPI)
  endif()

  message(STATUS "Project sources = " ${PROJECT_SOURCE_DIR} )
  include_directories(${PROJECT_SOURCE_DIR}/)
  include_directories(${CMAKE_BINARY_DIR}/)
  include_directories(${BML_INCLUDEDIR})
  include_directories(${PROGRESS_INCLUDEDIR})

  function(progress_appendix myappendix main_and_srcs)
  list(GET main_and_srcs 0 main)
  include_directories(${PROGRESS_INCLUDEDIR})
  add_executable(${myappendix} ${main})
  target_sources(${myappendix} PRIVATE ${ARGN})
  target_link_libraries(${myappendix} PUBLIC ${LINK_LIBRARIES})
  set_target_properties(${myappendix} PROPERTIES LINK_FLAGS "")
  endfunction(progress_appendix)

  progress_appendix(main main.F90
                            myMod.F90
                            )

  install(TARGETS main DESTINATION ${CMAKE_INSTALL_BINDIR})
                                                        

Feel free to modify the code according to your requirements and desired functionality. More modules can be easily added 
in the ``CMakeLists.txt`` file. 
Once you have completed your sample code, you can proceed with compiling it as follows::
  
    mkdir build ; cd build 
    cmake -DCMAKE_PREFIX_PATH="$HOME/qmd-progress/install/;$HOME/bml/install" ../src/
    make 

Remember to refer to the documentation of the PROGRESS and BML libraries for further details on how to utilize their features effectively. In order to run the code we just need to type::

    ./main

Let's now build a sample Hamiltonian matrix according to reference [Finkelstein]_. Details on the parameters 
and how to use this API call ca be found at: `Model Hamiltonian <_static/doxy/namespaceprg__modelham__mod.html#ae10c14620b7d6a3b001a3ca0eb785fff>`_. The code will need to be changed as follows:

.. code-block:: fortran
  
   module myMod
    use bml
    use prg_modelham_mod
    contains

    subroutine mySub
     implicit none
     real(8) :: ea, eb, dab, daiaj, dbibj, dec, rcoeff
     integer :: norbs, prec, seed, verbose
     logical :: reshuffle
     type(bml_matrix_t) ::  ham_bml

     norbs=100
     prec = kind(1.0d0)
     call bml_zero_matrix("dense",bml_element_real,prec,norbs,norbs,ham_bml)

     ea = 0.0d0; eb = 0.0d0; dab = -2.0d0; daiaj = 0.0d0 ; dbibj = -1.0d0
     dec = -1000.0d0; rcoeff = 0.0d0; reshuffle = .false. ; seed = 123; verbose = 1
     call prg_twolevel_model(ea, eb, dab, daiaj, dbibj, &
       dec, rcoeff, reshuffle, seed, ham_bml, verbose)

   end subroutine mySub

  end module myMod

Running this code will produce a 100x100 Model Hamiltonian Matrix that one can use to test any PROGRESS algorithm.  The output will only show part of the matrix: 

.. code-block:: bash

    h_bml
        0.000 -2.000  0.000 -2.020  
        -2.000  0.000 -2.000 -1.000 
        0.000 -2.000  0.000 -2.000  

Building a Density Matrix
-------------------------
One of the most important bottlenecks in computational chemistry is the calculation of the density matrix (DM). Usually this is calculated by direct application of the Fermi function. The method involves performing a matrix diagonalization in which all the computational effort is concentrated. Here we will use a PROGRESS library call to build the density matrix from the Hamiltonian using different methods.


Direct Fermi function application
#################################

Follow the steps provided on the section before to obtain a Hamiltonian matrix to work with. Then, add the density matrix
module in the scope section on the ``myMod.F90`` module file as follows:

.. code-block:: fortran

   module myMod
    use bml
    use prg_modelham_mod
    use prg_densitymatrix_mod !Density matrix module


Add the following lines to the scope of the subroutine 
in order to define the necessary variables: 

.. code-block:: fortran

     real(8), allocatable :: eigenvalues(:)
     real(8) :: bndfil
     integer, parameter :: dp = 8
     real(8) :: threshold
     type(bml_matrix_t) ::  rho_bml


Then, add the following lines 
after the Hamiltonian is constructed: 

.. code-block:: fortran
   
   allocate(eigenvalues(norbs))
   call bml_zero_matrix(bml_type,bml_element_real,dp,norbs,norbs,rho_bml)
   threshold = 1.0D-5 !Threshold value to eliminate small elements
   bndfil = 0.5 !Electronic filling factor (half of the states will be filled)
   !Computing the density matrix with diagonalization
   call prg_build_density_T0(ham_bml, rho_bml, threshold, bndfil, eigenvalues)
   call bml_print_matrix("rho_bml",rho_bml,0,10,0,10)

This will construct the DM with a direct application of the Fermi function. For a theoretical explanation on this
see [Koskinen]_ and [Niklasson]_ . One can use the output eigenvalues to plot the Density Of States (DOS) by Adding the following line in the scope of the subroutine:

.. code-block:: fortran

    use prg_dos_mod 

and the following code block after the DM is constructed:

.. code-block:: fortran

    !Computing the Fermi Level/Chemical potential
    ef = (eigenvalues(int(norbs/2)+1) + eigenvalues(int(norbs/2)))/2
    eigenvalues = eigenvalues - ef

    !Writting the total DOS
    call prg_write_tdos(eigenvalues, 0.05d0, 10000, -20.0d0, 20.0d0, "tdos.dat")

One can then plot the data from ``tdos.dat`` using `xmgrace <https://plasma-gate.weizmann.ac.il/Grace/>`_ or any other plotting tool. To know more about the parametes used in the ``prg_write_tdos`` subroutine, reffer to `prg_dos_mod <_static/doxy/namespaceprg__dos__mod.html>`_.

SP2 Algorithm
###############

In this section we will apply the "Second order spectral purification method," or SP2 algorithm. This algorithm consists of a series of matrix multiplications that attempt to "purify" the spectrum of the Hamiltonian matrix, resulting in a matrix with eigenvalues 0 or 1 depending on whether the initial eigenvalue of the Hamiltonian was above or below the Fermi level. We will hence replace the codeblock above by the following one:

.. code-block:: fortran

   call bml_zero_matrix(bml_type,bml_element_real,dp,norbs,norbs,rho_bml)
   threshold = 1.0D-5 !Threshold value to eliminate small elements
   bndfil = 0.5 !Electronic filling factor (half of the states will be filled)
   
   call prg_sp2_alg1(ham_bml,rho_bml,threshold,bndfil,15,100 &
         ,"Rel",1.0D-10,20)


This will solve for DM using the SP2 method.

Congruence transfomation 
--------------------------------------

We will construct the congruence transformation from the overlap matrix. For this, we will use a proxy overlap where orbitals i amd j are overlapping with 
a function :math:`S_{ij} = \exp(-|j - i|)`. Note that typically the overlap matrix is computed from the chemical system and further details about this could be found in [Negre2016]_. 
We will start adding the following lines to the module scope: ``use prg_genz_mod; use prg_nonortho_mod``. 
The following are the heading lines to ba added to the scope of the routine: 

.. code-block:: fortran

     type(bml_matrix_t) ::  smat_bml
     real(8), allocatable :: smat(:,:)
     integer :: i,j


The condeblock to be added to generate the overlap matrix ``smat_bml`` is the following: 

.. code-block:: fortran

     allocate(smat(norbs,norbs))
     do i = 1,norbs
      do j = 1,norbs
        smat(i,j) = exp(-1.0*real(abs(j-i),8))
      enddo
     enddo

     call bml_import_from_dense(bml_type,smat,smat_bml,threshold,norbs)

To obtain a congruence transformation matrix ``zmat_bml`` we will add the following lines:

.. code-block:: fortran

     call bml_zero_matrix(bml_type,bml_element_real,dp,norbs,norbs,zmat_bml)
     call prg_buildzdiag(smat_bml,zmat_bml,threshold,norbs,bml_type)
     call bml_print_matrix("zmat_bml",zmat_bml,0,10,0,10)

Other linear scaling algorithms can be also used in combination with sparse bml matrix types. 
This can be seen in: `Congruence transformation <_static/doxy/namespaceprg__genz__mod.html>`_.
Once the matrix zmat_bml is obtained one can "orthogonalize" the Hamiltonian matrix using routines 
in `Orthogonalization/deorthogonalization  <_static/doxy/namespaceprg__nonortho__mod.html>`_.

Handling chemical system
------------------------

Although this is not the main purpose of the progress library, several tools are in place to handle chemical systems. For instance, one can read and write a ``pdb``, ``xyz``, ``dat`` (LATTE input), and ``lmp`` (lammps input) file by calling a routine. The module to be used is the ``prg_system_mod``.
The system derived type is then used to access all the systems information, including coordinates and atomic types. An example follows. Lets create a coordinate ``coords.xyz`` file as follows::

    3 
    h2o initial system
    O 0.0 0.0 0.0 
    H 0.0 0.0 1.0
    H 0.0 1.0 0.0 

This system can be read/parsed as follows: 

.. code-block:: fortran

    call prg_parse_system(sy,"coords.xyz")

Details about the system type can be found at: `System type <_static/doxy/structprg__system__mod_1_1system__type.html>`_.


Referenece
----------

.. [Koskinen] Koskinen, Pekka, and Ville Mäkinen. 2009. “Density-Functional Tight-Binding for Beginners.” Computational Materials Science 47 (1): 237–53.
.. [Niklasson] Niklasson, Anders M. N., and Matt Challacombe. 2004. “Density Matrix Perturbation Theory.” Physical Review Letters 92 (19): 193001.
.. [Finkelstein] J. Finkelstein, C. Negre, J-L. Fattebert. 2023. `"A fast, dense Chebyshev solver for electronic structure on GPUs" <https://arxiv.org/abs/2306.12616>`_.
.. [Negre2016] Negre, Christian F. A., Susan M. Mniszewski, Marc J. Cawkwell, Nicolas Bock, Michael E. Wall, and Anders M. N. Niklasson. 2016. “Recursive Factorization of the Inverse Overlap Matrix in Linear-Scaling Quantum Molecular Dynamics Simulations.” Journal of Chemical Theory and Computation 12 (7): 3063–73.
