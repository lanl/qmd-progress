!> High-level program to perform GRAPH-BASED QMD.
!!
!! \ingroup PROGRAMS
!! \brief  High-level program to perform GRAPH-BASED QMD with a DFTB Hamiltonian
!!  using a full parallel Graph-based approach.
!! \note To test this program with the \verbatim run_test \endverbatim script.
!!
!! \author C. F. A. Negre
!! (cnegre@lanl.gov)
!!
!! \author S. Mniszewski
!! (smn@lanl.gov)
!!
!! \author A. M. N. Niklasson
!! (amn@lanl.gov)
!!
!! \author M. E. Wall
!! (mewall@lanl.gov)
!!
!! Verbose levels:
!!   >= 0- Only print tags.
!!   >= 1- Print useful scalar physical data (e.g., Total Energy)
!!   >= 2- Print single file output data (e.g., 0-SCF Charges)
!!   >= 3- Write out trajectory data.
!!   >= 4- Write out physically meaningful 1D arrays (e.g., Charges for the species)
!!   >= 5- Write out physically meaningful 2D arrays (e.g., H matrix)
!!
program gpmd
 
  use gpmdcov_vars 

!!!!!!!!!!!!!!!!!!!!!!!!
  !> Main program driver
!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize the program variables and parse input files.
  call gpmdcov_Init()

  !We give a first guess of the Fermi level.
  Ef =  lt%efermi

  !The inverse of the electronic temperature.
  beta = 1.0_dp/lt%kbt 

  !> Initial partition of the system based on the covalency graph.
  !! This will need to be replace by a first SP2 algorithm to compute a
  !! first density matrix.
  call gpmdcov_Part()

  !> Initialize partitions.
  call gpmdcov_InitParts()

  !> Comput first charges.
  call gpmdcov_FirstCharges()

  !> First SCF loop up to maxscf.
  call gpmdcov_DM_Min(lt%maxscf,sy%net_charge,.true.)

  !> First calculation of energies and forces.
  call gpmdcov_EnergAndForces(sy%net_charge)

  !> Setup the Molecular Dynamics (MD) calculation.
  call gpmdcov_PrepareMD()

  !> Perform the MD simulation.
  call gpmdcov_MDloop()

  !> Finalize the program.
  call gpmdcov_Finalize()

end program gpmd
