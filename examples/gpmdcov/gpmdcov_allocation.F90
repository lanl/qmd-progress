module gpmdcov_allocation_mod

  use bml
 
  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  public :: gpmdcov_reallocate_intVect
  public :: gpmdcov_reallocate_realVect
  public :: gpmdcov_reallocate_realMat
  public :: gpmdcov_reallocate_denseBmlRealMat

contains

  !> To reallocate a vector
  !!
  subroutine gpmdcov_reallocate_realVect(myVect,mySize)
    implicit none 
    real(dp),allocatable :: myVect(:)
    integer :: mySize

    if(allocated(myVect))then
      deallocate(myVect)
    endif
    allocate(myVect(mySize))

  end subroutine gpmdcov_reallocate_realVect

  !> To reallocate a matrix
  !!
  subroutine gpmdcov_reallocate_realMat(myMat,myRows,myCols)
    implicit none
    real(dp),allocatable :: myMat(:,:)
    integer :: myRows,myCols

    if(allocated(myMat))then
      deallocate(myMat)
    endif
    allocate(myMat(myRows,myCols))

  end subroutine gpmdcov_reallocate_realMat


  !> To reallocate a vector of integers
  !!
  subroutine gpmdcov_reallocate_intVect(myVect,mySize)
    implicit none
    integer, allocatable :: myVect(:)
    integer :: mySize

    if(allocated(myVect))then
      deallocate(myVect)
    endif
    allocate(myVect(mySize))

  end subroutine gpmdcov_reallocate_intVect


 !> To reallocate a bml matrix
 !!
  subroutine gpmdcov_reallocate_denseBmlRealMat(myMat,mySize)

    implicit none
    integer :: mySize
    type(bml_matrix_t) :: myMat

    if(bml_allocated(myMat))then
      call bml_deallocate(myMat)
    endif
    call bml_zero_matrix("dense",bml_element_real,dp,mySize,mySize,myMat)

  end subroutine gpmdcov_reallocate_denseBmlRealMat



end module gpmdcov_allocation_mod
