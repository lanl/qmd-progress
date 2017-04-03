
module hamiltonian_mod

  use accuracy_mod 
  use prg_openfiles_mod
  use bml

  implicit none

  private 

  public :: h_read, read_matrix

contains

  ! This reads the hamiltonian from a file          
  ! Hamiltonian. Units must be in eV 	          
  subroutine h_read(ham,hdim)

    implicit none
    integer :: i,j,cont,Naux,io
    integer :: ii, jj, nnz
    integer, intent(out) :: hdim
    real(dp) :: prom, conv_units, val
    real(dp), allocatable :: ham(:,:)
    character(20) :: dummy1, dummy2    

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
        if(j.GT.i)ham(i,j)=ham(j,i) !Enforced symmetrization
      enddo
    enddo

    close(io)

  end subroutine h_read

  
  !! This reads a quare matrix from a file          
  !!
  subroutine read_matrix(mat,hdim,filename)

    implicit none
    integer :: i,j,cont,Naux,io
    integer :: ii, jj, nnz
    integer, intent(out) :: hdim
    real(dp) :: prom, conv_units, val
    real(dp), allocatable :: mat(:,:)
    character(20) :: dummy1, dummy2   
    character(len=*) :: filename  

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
        if(j.GT.i)mat(i,j)=mat(j,i) !Enforced symmetrization
      enddo
    enddo

    close(io)

  end subroutine read_matrix
    
end module hamiltonian_mod
