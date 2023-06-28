
Tutorials
=========

In this tutorial, we will provide some examples of how to effectively use the library. To follow along, it is important that you are already familiar with the xx library and have successfully compiled it.

Sample Code
-----------

To begin, let's create a sample code that can be easily compiled by utilizing the necessary "flags" and environment variables from both xx and x libraries. In this example, we assume that the x and xx libraries are located in the `$HOME` directory.

1. Open your terminal and navigate to your home directory::

   cd ~

2. Create a new directory called "exampleCode" to store our code::

   mkdir exampleCode

3. Change into the "exampleCode" directory::

   cd exampleCode

4. Create a folder called "src" to hold our source code files::

   mkdir src

Now you are ready to start writing your sample code within the "src" directory. This structure will help keep your code organized and maintainable throughout the development process.

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
  
Next, we will write a symple ``CMakeList.txt`` inside ``/scr`` which will allow us to automatically 
compile our code using all the fortran flags as for BML and PROGRESS and will also lint to these 
libraries. The ``CMakeList.txt`` will look like: 

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
in the ``CMakeList.txt`` file. 
Once you have completed writing your sample code, you can proceed with compiling it as follows::
  
    mkdir build ; cd build 
    cmake -DCMAKE_PREFIX_PATH="$HOME/qmd-progress/install/;$HOME/bml/install" ../src/

Remember to refer to the documentation of the xx and x libraries for further details on how to utilize their features effectively. Let's now build a sample Hamiltonian according to reference J. Finkelstein, C. Negre, J-L. Fattebert;  https://arxiv.org/abs/2306.12616; 2023. For this we will change the module as follows:

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
     dec = 0.01d0; rcoeff = 0.0d0; reshuffle = .false. ; seed = 123; verbose = 1
     call prg_twolevel_model(ea, eb, dab, daiaj, dbibj, &
       dec, rcoeff, reshuffle, seed, ham_bml, verbose)

   end subroutine mySub

  end module myMod


Running this code will produce a 100x100 Model Hamiltonian Matrix that one can use to test any PROGRESS algorithm. 

Building a Density Matrix
-------------------------

Direct Fermi function application
#################################

Follow the steps provided on the section before to obtain a Hamiltonian matrix to work with. Add the following lines 
after the Hamiltonian is constructed. 

.. code-block:: fortran
   
   allocate(eigenvalues(norbs))
   call bml_zero_matrix(mham%bml_type,bml_element_real,dp,norbs,norbs,rho_bml)
   threshold = 1.0D-5 !Threshold value to eliminate small elements
   bndfil = 0.5 !Electronic filling factor (half of the states will be filled)
   !Computing the density matrix with diagonalization
   call prg_build_density_T0(ham_bml, rho_bml, threshold, bndfil, eigenvalues)

This will construct the DM with a direct application of the Fermi function. For a theoretical explanation on this
see [Koskinen]_ and [Niklasson]_ .

# SP2 Algorithm #

We will apply the SP2 algorithm ...

Congruence transfomation 
########################

We will construct the congruence transform 


Handling chemical system
------------------------


Running on HPC machines 
-----------------------

Adding a new PROGRESS fun


Referenece
----------

.. [Koskinen] Koskinen, Pekka, and Ville Mäkinen. 2009. “Density-Functional Tight-Binding for Beginners.” Computational Materials Science 47 (1): 237–53.
.. [Niklasson] Niklasson, Anders M. N., and Matt Challacombe. 2004. “Density Matrix Perturbation Theory.” Physical Review Letters 92 (19): 193001.


