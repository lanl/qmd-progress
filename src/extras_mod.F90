!> Extra routines: 
!! @ingroup PROGRESS 
!! \brief A module to add any extra routine considered necessary but which is NOT 
!! essential for any other PROGRESS routines.
!!      
module extras_mod

  use openfiles_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  public :: print_matrix, delta, mls, get_mem

contains

  !> To write a dense matrix to screen.
  !! \param matname Matrix name.
  !! \param amat Matrix to be printed. 
  !! \param i1 Print from row i1.
  !! \param i2 Print up to from row i2.
  !! \param j1 Print from column j1.
  !! \param j2 Print up to column j2.  
  !! 
  subroutine print_matrix(matname,amat,i1,i2,j1,j2)
    implicit none
    integer :: ndim, i, j
    integer, intent (in) :: i1,i2,j1,j2
    integer :: ii2,jj2
    real(dp), intent (in) :: amat(:,:)
    character(len=*) :: matname

    if(i1 > i2)stop "Error at print_matrix, i1 > i2"
    if(j1 > j2)stop "Error at print_matrix, j1 > j2"      

    ndim = size(amat,dim=1)
    if(i2 > ndim)then 
      ii2=ndim
    else
      ii2=i2
    endif   

    if(j2 > ndim)then 
      jj2=ndim
    else
      jj2=j2
    endif   

    write(*,*)""
    write(*,*)" ============================================== "
    write(*,*)matname
    do i = i1, ii2
      WRITE(*,'(10F15.10)') (AMAT(i,j), j = j1,jj2)
    end do
    write(*,*)" ============================================== "
    write(*,*)""

  end subroutine print_matrix


  !> To get the actual time in milliseconds. 
  !! \param mls Output value with the machine time in milliseconds. 
  !!
  function mls()
    real(dp) :: mls
    integer :: timevector(8)

    mls = 0.0_dp
    call date_and_time(values=timevector)
    mls=timevector(5)*60.0_dp*60.0_dp*1000.0_dp + timevector(6)*60.0_dp*1000.0_dp &
      + timevector(7)*1000.0_dp + timevector(8)

  end function mls

  !> Delta function ||X^tSX - I||. CFAN, March 2015.
  !! 
  subroutine delta(x,s,nn,dta)
    implicit none
    integer :: i, j, nn
    real(dp) :: x(nn,nn),s(nn,nn),temp1(nn,nn),temp2(nn,nn), dta
    real(dp) :: identity(nn,nn)

    identity=0.0

    do j = 1, nn
      identity(j,j)=1.0
    enddo

    temp1=matmul(transpose(x),s)
    temp2=matmul(temp1,x)

    do j = 1, nn
      identity(j,j)=1.0
    enddo

    temp1=0.0
    do i = 1, nn
      do j = 1, nn
        temp1(i,j) = identity(i,j)-temp2(i,j)
      enddo
    enddo

    !Take the max absolute value of the leading eigenvectors.
    call twonorm(temp1,nn,dta)

  end subroutine delta

  
  !> Get proc memory
  !! \param procname Process name to get the mem usage.
  !! \param tag Tag to pprint the processor mem usage. 
  !! 
  subroutine get_mem(procname,tag)
    implicit none
    character(*), intent(in) :: procname
    character(*), intent(in) :: tag
    character(200) :: command    

    command = "echo 'Used mem "//tag//"=' $(top -n 1 -b | grep  "//procname//" | head -n 1 | awk '{print $10;}')"
    
    call system(command)    

  end subroutine get_mem
    
  ! norm2. CFAN, March 2015.
  subroutine twonorm(a,nn,norm2)
    implicit none
    integer :: info, nn
    real(dp) :: a(nn,nn), norm2
    integer :: tmp_lwork
    real(dp) :: utmp(nn,nn), tmp_evals(nn)
    real(dp), allocatable :: tmp_work(:)

    tmp_lwork=3*nn -1
    allocate(tmp_work(tmp_lwork))

    utmp=a

    call dsyev("v", "u", nn, utmp, nn, tmp_evals, tmp_work, &
    tmp_lwork,  info) 

    norm2=max(abs(tmp_evals(1)),abs(tmp_evals(nn)))

    deallocate(tmp_work)

  end subroutine twonorm

end module extras_mod   
