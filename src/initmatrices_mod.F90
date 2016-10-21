!> Initialization module. 
!! \ingroup PROGRESS
!! \brief Routines in this module are used to initialize several matrices that will be used in the 
!! code. 
!! 
module initmatrices_mod

  use bml 

  implicit none

  private 

  integer, parameter :: dp = kind(1.0d0)

  public :: init_hsmat,init_pzmat,init_ortho

contains 

  !> Initialize Hamiltonian and Overlap Matrix.
  !! \brief Allocation of the Hamiltonian and Overlap matrix into bml formats. 
  !! \param ham_bml Hamiltonian in bml format. 
  !! \param over_bml Overlap in bml format. 
  !! \param threshold Threshold value for matrix elements. 
  !! \param mdim Max nonzero elements per row for every row see \cite Mniszewski2015 .
  !! \param norb Total number of orbitals.
  subroutine init_hsmat(ham_bml,over_bml,bml_type,mdim,norb)
    implicit none
    type(bml_matrix_t), intent (inout) :: ham_bml,over_bml 
    integer, intent(in) :: norb
    integer, intent(inout) :: mdim
    character(20) :: bml_type

    if(mdim < 0)mdim = norb 
    !Allocate bml's  
    call bml_zero_matrix(bml_type,bml_element_real,dp,mdim,norb,ham_bml)                 
    call bml_zero_matrix(bml_type,bml_element_real,dp,mdim,norb,over_bml)      

  end subroutine init_hsmat  

  !> Initialize Density matrix and Inverse square root Overlap.
  !! \brief Allocation of the Density matrix and Inverse square root Overlap matrix into bml formats. 
  !! \param rho_bml Density matrix in bml format. 
  !! \param zmat_bml Inverse square root Overlap in bml format. 
  !! \param threshold Threshold value for matrix elements. 
  !! \param mdim Max nonzero elements per row for every row see \cite Mniszewski2015 .
  !! \param norb Total number of orbitals.
  subroutine init_pzmat(rho_bml,zmat_bml,bml_type,mdim,norb)
    implicit none
    type(bml_matrix_t), intent (inout) :: rho_bml,zmat_bml 
    integer, intent(in) :: norb
    integer, intent(inout) :: mdim
    character(20) :: bml_type

    if(mdim < 0)mdim = norb 
    !Allocate bml's  
    call bml_zero_matrix(bml_type,bml_element_real,dp,mdim,norb,rho_bml)                 
    call bml_zero_matrix(bml_type,bml_element_real,dp,mdim,norb,zmat_bml)      

  end subroutine init_pzmat

  !> Initialize The orthogonal versions of Hamiltonian and Density Matrix.
  !! \brief Allocation of the orthogonal Hamiltonian and Density matrix into bml formats. 
  !! \param orthoh_bml Orthogonal Hamiltonian in bml format. 
  !! \param orthop_bml Orthogonal Density Matrix in bml format. 
  !! \param threshold Threshold value for matrix elements. 
  !! \param mdim Max nonzero elements per row for every row see \cite Mniszewski2015 .
  !! \param norb Total number of orbitals.
  subroutine init_ortho(orthoh_bml,orthop_bml,bml_type,mdim,norb)
    implicit none
    type(bml_matrix_t), intent (inout) :: orthoh_bml,orthop_bml 
    integer, intent(in) :: norb
    integer, intent(inout) :: mdim
    character(20) :: bml_type

    if(mdim < 0)mdim = norb 
    !Allocate bml's  
    call bml_noinit_matrix(bml_type,bml_element_real,dp,mdim,norb,orthoh_bml)                 
    call bml_noinit_matrix(bml_type,bml_element_real,dp,mdim,norb,orthop_bml)      

  end subroutine init_ortho

end module initmatrices_mod
