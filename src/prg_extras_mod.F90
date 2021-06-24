!> Extra routines.
!! @ingroup PROGRESS
!! \brief A module to add any extra routine considered necessary but which is NOT
!! essential for any other PROGRESS routine.
!!
module prg_extras_mod

  use prg_openfiles_mod

  implicit none

  private

  interface
    subroutine prg_memory_consumption(vm_peak, vm_size, pid, ppid) &
         bind(C, name = "prg_memory_consumption")
      use, intrinsic :: iso_C_binding
      integer(C_LONG_LONG), intent(inout) :: vm_peak
      integer(C_LONG_LONG), intent(inout) :: vm_size
      integer(C_LONG_LONG), intent(inout) :: pid
      integer(C_LONG_LONG), intent(inout) :: ppid
    end subroutine prg_memory_consumption
  end interface

  interface to_string
    module procedure to_string_integer
    module procedure to_string_long_long
    module procedure to_string_double
  end interface to_string

  integer, parameter :: dp = kind(1.0d0)

  public :: mls
  public :: prg_delta
  public :: prg_get_mem
  public :: prg_print_matrix
  public :: to_string
  public :: prg_norm2

contains

  !> Convert integer to string.
  !! \param i The integer
  !! \return The string
  !!
  function to_string_integer(i)

    character(len=:), allocatable :: to_string_integer
    integer, intent(in) :: i
    character(len=20) :: buffer

    write(buffer, "(I20)") i
    allocate(character(len_trim(adjustl(buffer))) :: to_string_integer)
    to_string_integer = trim(adjustl(buffer))

  end function to_string_integer

  !> Convert integer to string.
  !! \param i The integer
  !! \return The string
  !!
  function to_string_long_long(i)

    use, intrinsic :: iso_C_binding

    character(len=:), allocatable :: to_string_long_long
    integer(kind=C_LONG_LONG), intent(in) :: i
    character(len=30) :: buffer

    write(buffer, "(I30)") i
    allocate(character(len_trim(adjustl(buffer))) :: to_string_long_long)
    to_string_long_long = trim(adjustl(buffer))

  end function to_string_long_long

  !> Convert double to string.
  !! \param x The double
  !! \return The string
  !!
  function to_string_double(x)

    character(len=:), allocatable :: to_string_double
    double precision, intent(in) :: x
    character(len=20) :: buffer

    write(buffer, "(ES20.8)") x
    allocate(character(len_trim(adjustl(buffer))) :: to_string_double)
    to_string_double = trim(adjustl(buffer))

  end function to_string_double

  !> To write a dense matrix to screen.
  !! \param matname Matrix name.
  !! \param amat Matrix to be printed.
  !! \param i1 Print from row i1.
  !! \param i2 Print up to from row i2.
  !! \param j1 Print from column j1.
  !! \param j2 Print up to column j2.
  !!
  subroutine prg_print_matrix(matname,amat,i1,i2,j1,j2)

    integer :: ndim, i, j
    integer, intent (in) :: i1,i2,j1,j2
    integer :: ii2,jj2
    real(dp), intent (in) :: amat(:,:)
    character(len=*) :: matname

    if(i1 > i2)stop "Error at prg_print_matrix, i1 > i2"
    if(j1 > j2)stop "Error at prg_print_matrix, j1 > j2"

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
      write(*,'(10F15.10)') (AMAT(i,j), j = j1,jj2)
    end do
    write(*,*)" ============================================== "
    write(*,*)""

  end subroutine prg_print_matrix


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

  !> Delta function ||X^tSX - I||.
  !! \param x input matrix.
  !! \param s overlap matrix.
  !! \param dta Delta output value.
  !!
  subroutine prg_delta(x,s,nn,dta)

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
    call prg_twonorm(temp1,nn,dta)

  end subroutine prg_delta


  !> Get proc memory
  !! \param procname Process name to get the mem usage.
  !! \param tag Tag to pprint the processor mem usage.
  !!
  subroutine prg_get_mem(procname, tag)

    use, intrinsic :: iso_C_binding

    character(*), intent(in) :: procname
    character(*), intent(in) :: tag
    character(200) :: command
    integer(kind=C_LONG_LONG) :: vm_peak, vm_size, pid, ppid

    call prg_memory_consumption(vm_peak, vm_size, pid, ppid)

    write(*, *) "Used mem "//trim(tag) &
         //" (pid "//to_string(pid)//", " &
         //" ppid "//to_string(ppid)//") = " &
         //trim(to_string(vm_size))//" MiB (" &
         //trim(to_string(vm_peak))//" MiB)"

  end subroutine prg_get_mem

  !> Gets the norm2 of a square matrix.
  !! \param a Square matrix.
  !! \param nn Matrix size.
  !! \param norm2 Two-norm of matrix a.
  !!
  subroutine prg_twonorm(a,nn,norm2)

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

  end subroutine prg_twonorm

  !> Gets the norm2 of a vector.
  !! \param a Vector.
  !!
  real(dp) function prg_norm2(a)

    integer :: nn, i
    real(dp), intent(in) :: a(:)

    nn = size(a,dim=1)

#ifdef NORM2
    prg_norm2 = norm2(a)
#else
    do i = 1, nn
      prg_norm2 = prg_norm2 + a(i)*a(i)
    enddo
    prg_norm2 = sqrt(prg_norm2)
#endif

  end function prg_norm2

end module prg_extras_mod
