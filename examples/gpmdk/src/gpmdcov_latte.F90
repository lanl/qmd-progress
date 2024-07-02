!> Module for calling original latte routines
!! \brief This module will be used to call latte routines 
!! using the origial code
!!
!! @ingroup GPMDCOV
!!
!!
module gpmdcov_latte_mod

  use gpmdcov_vars
  use prg_system_mod
  ! Latte modules 
#ifdef USE_LATTE 
  use setuparray
  use prg_system_mod
  use latteparser
  use constants_mod
  use bml
  use nonoarray
  use NEBLISTARRAY
  use univarray
#endif

public :: init_latte

contains

  !> Initialize LATTE variables by passing
  !! the total system.
  !! \brief This will read the latte input file and the TB 
  !! parametere. !!!WARNING!!! The parameters (cutoff, threshold, etc.)
  !! need to be consistent with the gpmd input files.
  !!
  subroutine init_latte()
    type(system_type) :: mysy
#ifdef USE_LATTE
    !> Read the LATTE control file
    call parse_control("latte.in")
    call readtb()
#endif
  end subroutine init_latte
  

  !> Get the H and S based on the coordinates and symbols
  !! \brief Constructing the latte Hamiltonian and Overlap matrices. 
  !! This is done with the information contained in the latte input 
  !! file and the TB parameters files.
  !! \param myham0_bml Non-scf Hamiltonian bml format.
  !! \param myover_bml Overlap matrix in bml format.
  !! \param mycoords Coordinates of the system/sub-system.
  !! coords(1,10): x-coordinate of atom 10.
  !! \param mysymbols Symbol for every atom.
  !! \param mylattice_vectors Lattice vectors. mylattice_vectors(1,:): First 
  !! lattice vector.
  !! \param thresh Threshold value used to compute the matrices. 
  !!
  subroutine get_hsmat_latte(myham0_bml,myover_bml,mycoords,mysymbols,&
        mylattice_vectors,thresh)
    implicit none
    type(bml_matrix_t) :: myham0_bml,myover_bml
    real(dp), allocatable :: mycoords(:,:)
    character(2), allocatable :: mysymbols(:)
    real(dp), allocatable :: mylattice_vectors(:,:)
    real(dp) :: thresh
 
#ifdef USE_LATTE
  
    nats = size(mycoords,dim=2)
    if(allocated(cr)) deallocate(cr)
    allocate(cr(3,nats))  
    cr = mycoords  
    box = mylattice_vectors
    if(allocated(atele)) deallocate(atele) 
    allocate(atele(nats))
    do i=1,nats
      atele(i) = trim(adjustl(mysymbols(i)))
    enddo
    call readcr() 

    !Get the dimension of the Hamiltonian
    call gethdim()
    IF (basistype .EQ. "NONORTHO")then 
      call allocatenono
    endif

    !Get cutoff list. This might also be precomputed and passed 
    !for only the subsystem.
    call gencutofflist()

    !!!WARNING
    !Build integral mapp. This needs to be computed only once
    !for the sull system and should be passed for the partial 
    !system through this routine. The integral map could become 
    !a memory problem if computed for the full sytem!  
    call build_integral_map()

    !Construct the Hamiltonian and Overlap
    call genhonsite()
    call genorbitallist() 
    call getmatindlist() 
    call allocatenebarrays()

    !!!WARNING
    !This could be avoided since H and S for the part only need 
    !the tranlation vectors of the full system but not necessary the 
    !NL. At this point the NL was already applied.
    call neblists(0)
    
    call bldnewhs()

    call bml_zero_matrix("dense",bml_element_real,latteprec,hdim,hdim,myham0_bml)
    call bml_import_from_dense("dense", H, myham0_bml, thresh, hdim)
    call bml_zero_matrix("dense",bml_element_real,latteprec,hdim,hdim,myover_bml)
    call bml_import_from_dense("dense",smat , myover_bml, thresh, hdim)
    
    if(lt%verbose >= 1 .and. myRank == 1)then 
      call bml_print_matrix("ham_bml",myham0_bml,0,10,0,10)
      call bml_print_matrix("over_bml",myover_bml,0,10,0,10)
    endif 

    call deallocatenono()
    call deallocatenebarrays()

    !The following are variables which are deallocated 
    !in the deallocateall latte module which implies 
    !the use of many other modules not necessary for gpmd.
    deallocate(f, fpp, ftot)
    deallocate(deltaq, mycharge)
    deallocate(elempointer)
    deallocate(lcnshift)
    deallocate(fpul, fscoul) 
    if(allocated(fsspin)) deallocate(fsspin)
    if(allocated(respchi))deallocate(respchi)
    deallocate(qlist)
    deallocate(h,hdiag,h0,bo,h_onsite)
    deallocate(cutoff_list)
    deallocate(igl_map,orbital_list)


    !stop

#endif
  end subroutine get_hsmat_latte



end module gpmdcov_latte_mod


