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

  ! Local modules
  use gpmdcov_vars
  use gpmdcov_dm_min_mod
  use gpmdcov_energandforces_mod
  use gpmdcov_prepareMD_mod
  use gpmdcov_mdloop_mod
  use gpmdcov_init_mod
  use gpmdcov_part_mod
  use gpmdcov_assert_mod
  use gpmdcov_mod
  use gpmdcov_diagonalize_mod
  use gpmdcov_kernel_mod

  ! Latte modules 
#ifdef USE_LATTE 
  use gpmdcov_latte_mod
  use setuparray
  use prg_system_mod
  use latteparser
  use constants_mod
  use bml
  use nonoarray
  use NEBLISTARRAY 
#endif

  !integer :: ipreMD
!!!!!!!!!!!!!!!!!!!!!!!!
  !> Main program driver
!!!!!!!!!!!!!!!!!!!!!!!!
  
  !> Initialize the program variables and parse input files.
  call gpmdcov_Init(.false.)
  
#ifdef USE_LATTE
  call init_latte()
#endif

  call gpmdcov_assert_input(myRank)
 
  !We give a first guess of the Fermi level.
  Ef =  lt%efermi

  !The inverse of the electronic temperature.
  beta = 1.0_dp/lt%kbt

  !> Initial partition of the system based on the covalency graph.
  call gpmdcov_Part(1)

  if(lt%stopAt == "gpmdcov_Part") stop
  
  !> Initialize partitions.
  call gpmdcov_InitParts()
  if(lt%stopAt == "gpmdcov_InitParts") stop
  
  !> Comput first charges.
  if(lt%method == "DiagEfFull") eig = .false.
  if(.not. eig) then
        call gpmdcov_Diagonalize_H0()
        if(lt%MuCalcType == "FromParts" .or. lt%MuCalcType == "Combined")then
                call gpmdcov_muFromParts()
                if(err_status)return
        endif
  endif
  call gpmdcov_FirstCharges(eig)
  if(lt%stopAt == "gpmdcov_FirstCharges") stop


  !> First SCF loop up to maxscf.
  if(eig)then 
        call gpmdcov_DM_Min(lt%maxscf,sy%net_charge,.true.)
  else 
        call gpmdcov_DM_Min_Eig(lt%maxscf,sy%net_charge,.true.,.false.)
  endif

  if(lt%stopAt == "gpmdcov_DM_Min") stop
  
  !> First calculation of energies and forces.
  call gpmdcov_EnergAndForces(sy%net_charge)
  if(lt%stopAt == "gpmdcov_Energ") stop

  !> Setup the Molecular Dynamics (MD) calculation.
  call gpmdcov_PrepareMD(gpmdt%temp0)
  if(lt%stopAt == "gpmdcov_PrepareMD") stop
  
  !> Perform the MD simulation.
  call gpmdcov_MDloop()
  if(lt%stopAt == "gpmdcov_MDloop") stop

  !> Finalize the program.
  call gpmdcov_Finalize()

end program gpmd
