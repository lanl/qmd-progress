module gpmdcov_allocation_mod

  use bml
 
  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  public :: gpmdcov_reallocate_intVect
  public :: gpmdcov_reallocate_realVect
  public :: gpmdcov_reallocate_charVect
  public :: gpmdcov_reallocate_realMat
  public :: gpmdcov_reallocate_denseBmlRealMat
  public :: gpmdcov_vect2MatInt
  public :: gpmdcov_mat2VectInt
  public :: gpmdcov_bml_set_N
  public :: gpmdcov_bml_allocated
  
contains

  function gpmdcov_bml_allocated(myMat)
    implicit none
    type(bml_matrix_t) :: myMat
    logical :: gpmdcov_bml_allocated

    gpmdcov_bml_allocated = bml_allocated(myMat)
  end function gpmdcov_bml_allocated
  
  subroutine gpmdcov_bml_set_N(myMat,mySize)
    implicit none
    integer :: mySize
    type(bml_matrix_t) :: myMat

    call bml_set_N_dense(myMat,mySize)

  end subroutine gpmdcov_bml_set_N
    
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


  !> To reallocate a vector of integers
  !!
  subroutine gpmdcov_reallocate_charVect(myVect,mySize)
    implicit none
    character(*), allocatable :: myVect(:)
    integer :: mySize

    if(allocated(myVect))then
      deallocate(myVect)
    endif
    allocate(myVect(mySize))

  end subroutine gpmdcov_reallocate_charVect


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

  !> To reformat an integer matrix into a vector
  !!
  subroutine gpmdcov_mat2VectInt(myMat,myVect,rows,cols)
    implicit none
    integer, allocatable :: myVect(:)
    integer, intent(in) :: myMat(:,:)
    integer, intent(in) :: cols,rows
    integer :: i,j

    if(allocated(myVect))then
      deallocate(myVect)
    endif
    allocate(myVect(cols*rows))

    do i=1,rows
        do j=1,cols
                myVect(j + (i-1)*cols) = myMat(j,i)
         enddo
    enddo

  end subroutine gpmdcov_mat2VectInt

  !> To reformat an integer vector into a matrix 
  !!
  subroutine gpmdcov_vect2MatInt(myVect,myMat,rows,cols)
    implicit none
    integer, allocatable :: myMat(:,:)
    integer, intent(in) :: myVect(:)
    integer, intent(in) :: cols,rows
    integer :: i,j

    if(allocated(myMat))then
      deallocate(myMat)
    endif
    allocate(myMat(cols,rows))
    
    do i=1,rows
        do j=1,cols
                myMat(j,i) = myVect(j + (i-1)*cols)
         enddo
    enddo

  end subroutine gpmdcov_vect2MatInt





end module gpmdcov_allocation_mod
