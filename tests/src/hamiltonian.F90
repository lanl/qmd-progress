
module hamiltonian_mod

  use prg_openfiles_mod
  use bml

  implicit none

  private

  integer, parameter, public :: dp = kind(1.0d0)
  public :: h_read, read_matrix

contains

  ! This reads the Hamiltonian from a file
  !
  subroutine h_read(ham,hdim)

    implicit none
    character(20)          ::  dummy1, dummy2
    integer                ::  Naux, cont, i, ii
    integer                ::  io, j, jj, nnz
    integer, intent(out)   ::  hdim
    real(dp)               ::  val
    real(dp), allocatable  ::  ham(:,:)

    call prg_open_file_to_read(io,'hamiltonian_ortho.mtx')

    read(io,*)dummy1
    read(io,*)hdim, hdim, nnz
    read(io,*)dummy2

    allocate(ham(hdim,hdim))

    ham=0.0_dp

    cont = 0

    do i=1,nnz
       read(io,*)ii,jj,val
       ham(ii,jj) = val
    enddo

    do i=1,hdim
       do j=1,hdim
          if(j.gt.i)ham(i,j)=ham(j,i) !Enforced symmetrization
       enddo
    enddo

    close(io)

  end subroutine h_read


  ! This reads a square matrix from a file
  !
  subroutine read_matrix(mat,hdim,filename)

    implicit none
    character(20)          ::  dummy1, dummy2
    character(len=*)       ::  filename
    integer                ::  Naux, cont, i, ii
    integer                ::  io, j, jj, nnz
    integer, intent(out)   ::  hdim
    real(dp)               ::  val
    real(dp), allocatable  ::  mat(:,:)

    call prg_open_file_to_read(io,trim(filename))

    read(io,*)dummy1
    read(io,*)hdim, hdim, nnz
    read(io,*)dummy2

    allocate(mat(hdim,hdim))

    mat=0.0_dp

    cont = 0

    do i=1,nnz
       read(io,*)ii,jj,val
       mat(ii,jj) = val
    enddo

    do i=1,hdim
       do j=1,hdim
          if(j.gt.i)mat(i,j)=mat(j,i) !Enforced symmetrization
       enddo
    enddo

    close(io)

  end subroutine read_matrix

end module hamiltonian_mod
