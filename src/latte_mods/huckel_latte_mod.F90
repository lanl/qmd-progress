!> A module to compute Extended huckel Hamiltonian and Overlap matrices.
!! \brief This module will compute H ans S from the EH parameter.
!! @ingroup EXTERNAL
!!
module huckel_latte_mod

  use bml

  private
  public :: get_hshuckel

  integer, parameter :: dp = kind(1.0d0)

  type :: single_orb  
    integer :: orb_of_e
    real(dp) :: VSIP
    real(dp) :: expo
    real(dp) :: exp2(2)
    real(dp) :: coef2(2)
  end type 

  type :: atom_parameter 
    character(3) :: symbol
    integer :: valence_electron
    type(single_orb) :: orb(4)
  end type 

  type :: atom
    integer atomtype;
    real(dp) x;
    real(dp) y;
    real(dp) z;
  end type 

contains

  subroutine get_hshuckel(ham_bml,over_bml,coordinates,spindex,spatnum&
      ,parampath,bml_type,mdim,threshold& 
      ,nsp,splist,basis,numel,onsite_energ,&
      norbi,hubbardu)

    character(len=*), intent(in) :: parampath
    type(atom_parameter), allocatable :: period(:)
    type(bml_matrix_t), intent(inout) :: ham_bml
    type(bml_matrix_t), intent(inout) :: over_bml
    integer, intent(inout) :: nsp
    integer, intent(in) :: spindex(:)
    integer, intent(in) :: spatnum(:)
    character(2), allocatable, intent(inout) :: splist(:)
    character(4),allocatable, intent(inout) :: basis(:)
    real(dp),allocatable, intent(inout) :: numel(:)
    real(dp),allocatable, intent(inout) :: onsite_energ(:,:)
    integer,allocatable, intent(inout) :: norbi(:)
    real(dp),allocatable, intent(inout) :: hubbardu(:)
    real(dp), intent(in) :: coordinates(:,:)
    type(atom), allocatable :: molecule(:)
    integer :: no_atoms, i,j, indexi, indexj
    integer :: noorb_i, atom_row, no_orbitals
    integer :: atom_col, noorb_j, ii, jj
    real(dp) :: delx, dely, delz, S(16,16),H(16,16)
    character(len=*), intent(in) :: bml_type
    integer, intent(inout) :: mdim
    real(dp), intent(in) :: threshold
    real(dp), allocatable :: row(:)

  end subroutine get_hshuckel

end module huckel_latte_mod
