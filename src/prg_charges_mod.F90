!> A module to compute the Mulliken charges of a chemical system.
!! \brief This module contains routines that compute properties related to charges.
!! @ingroup PROGRESS
!!
module prg_charges_mod

  use prg_openfiles_mod
  use bml
  use prg_parallel_mod
  use prg_graph_mod
  use prg_system_mod

  implicit none     

  private 

  integer, parameter :: dp = kind(1.0d0)

  public :: prg_get_charges, prg_get_hscf

contains

  !> Constructs the charges from the density matrix.
  !! \param rho_bml Density matrix in bml format.
  !! \param over_bml Overlap matrix in bml format.
  !! \param hindex Start and end index for every atom in the system. 
  !! \param charges Output parameter that gives the vectorized charges.
  !! \param threshold Threshold value for matrix elements.
  !!
  subroutine prg_get_charges(rho_bml,over_bml,hindex,charges,numel,spindex,mdimin,threshold)
    implicit none
    character(20)                        ::  bml_type, bml_mode
    integer                              ::  i, j, nats, norb
    integer                              ::  mdim
    integer, intent(in)                  ::  hindex(:,:), mdimin, spindex(:)
    real(dp)                             ::  znuc
    real(dp), allocatable                ::  rho_diag(:)
    real(dp), allocatable, intent(inout) ::  charges(:)
    real(dp), intent(in)                 ::  numel(:), threshold
    type(bml_matrix_t)                   ::  aux_bml
    type(bml_matrix_t), intent(inout)    ::  over_bml, rho_bml

    norb = bml_get_N(rho_bml)
    bml_type = bml_get_type(rho_bml)
    bml_mode = bml_get_distribution_mode(rho_bml)
    nats = size(hindex,dim=2)

    if(mdimin.lt.0)then 
      mdim = norb
    else
      mdim = mdimin
    endif  

    if(.not.allocated(charges)) allocate(charges(nats)) 
    allocate(rho_diag(norb))

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,mdim,aux_bml, &
      bml_mode) 

    call bml_multiply(rho_bml,over_bml,aux_bml,1.0_dp,0.0_dp,threshold)

#ifdef DO_MPI
    if (getNRanks() .gt. 1 .and. &
      bml_get_distribution_mode(aux_bml) == BML_DMODE_DISTRIBUTED) then
      call prg_allGatherParallel(aux_bml)
    endif
#endif
   
    call bml_get_diagonal(aux_bml,rho_diag)
    
    do i = 1,nats
      znuc = numel(spindex(i))
      charges(i)=0.0_dp
      do j = hindex(1,i),hindex(2,i)
        charges(i) = charges(i) + rho_diag(j)
      enddo
      charges(i) = charges(i) - znuc    
    enddo    
    
    deallocate(rho_diag)
    call bml_deallocate(aux_bml)

  end subroutine prg_get_charges 

  !> Constructs the SCF hamiltonian given H0, HubbardU and charges.
  !! This routine does: 
  !! \f$ H = \sum_i U_i q_i + V_i; \f$, where \f$ U \f$ is the Hubbard parameter for every atom i.
  !! \f$ V \f$ is the coulombic potential for every atom i. 
  !! \param ham_bml Hamiltonian in bml format.
  !! \param over_bml Overlap in bml format.
  !! \param hindex Start and end index for every atom in the system. 
  !! \param hubbardu Hubbard parameter for every atom. 
  !! \param charges Charges for every atom.
  !! \param coulomb_pot_r Coulombic potential (r contribution)
  !! \param coulomb_pot_k Coulombic potential (k contribution)
  !! \param mdim Maximum nonzeroes elements per row for every row.
  !! \param threshold Threshold value for matrix elements.
  !!
  subroutine prg_get_hscf(ham0_bml,over_bml,ham_bml,spindex,hindex,hubbardu,charges,&
      coulomb_pot_r,coulomb_pot_k,mdimin,threshold)
    implicit none
    character(20)                      ::  bml_type
    integer                            ::  i, j, nats, norb, mdim
    integer, intent(in)                ::  hindex(:,:), mdimin, spindex(:)
    real(dp), allocatable              ::  coulomb_pot(:), diagonal(:)
    real(dp), intent(in)               ::  charges(:), coulomb_pot_k(:), coulomb_pot_r(:), hubbardu(:)
    real(dp), intent(in)               ::  threshold
    type(bml_matrix_t)                 ::  aux_bml
    type(bml_matrix_t), intent(in)     ::  ham0_bml, over_bml
    type(bml_matrix_t), intent(inout)  ::  ham_bml

    nats = size(coulomb_pot_r)
    allocate(coulomb_pot(nats))

    norb = bml_get_N(ham0_bml) 

    if(mdimin.lt.0)then 
      mdim = norb
    else
      mdim = mdimin
    endif

    allocate(diagonal(norb))

    call bml_copy_new(ham0_bml,ham_bml)
        
    bml_type = bml_get_type(ham_bml)

    coulomb_pot = coulomb_pot_k + coulomb_pot_r

    call bml_zero_matrix(bml_type,bml_element_real,dp,mdim,norb,aux_bml) 

    do i = 1,nats
      do j = hindex(1,i),hindex(2,i) 
        diagonal(j) = hubbardu(spindex(i))*charges(i) + coulomb_pot(i)
      enddo
    enddo
    
    call bml_set_diagonal(aux_bml,diagonal)    
                   
    call bml_multiply(over_bml,aux_bml,ham_bml,0.5_dp,1.0_dp,threshold) !  h = h + 0.5*s*h1
      
    call bml_multiply(aux_bml,over_bml,ham_bml,0.5_dp,1.0_dp,threshold) !  h = h + 0.5*h1*s     

    deallocate(diagonal)
    deallocate(coulomb_pot)
    call bml_deallocate(aux_bml)

  end subroutine prg_get_hscf

end module prg_charges_mod    
