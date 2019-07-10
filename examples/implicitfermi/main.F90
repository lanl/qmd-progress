
!!
program main

  use bml

  !progress lib modes
  use prg_implicit_fermi_mod
  !use hamiltonian_mod

  implicit none

  integer, parameter :: dp = 8
  integer :: norb, mdim, verbose, nsteps,osteps 
  type(bml_matrix_t) :: ham_bml, p_bml, xi0_bml
  character(20) :: bml_type
  real(dp) :: threshold, mu, beta, occErrLimit, nocc  
  integer :: i,j
  real,allocatable :: vector(:)
  !real :: bml_get 

  

  !Some parameters that can be changed depending on the test.
  bml_type = "dense"
  threshold = 1.0d-9
  mdim = -1
  verbose = 1
  norb = 384 
  nsteps = 8
  osteps = 1
  occErrLimit = 1.0d-6
  mu = 0.2
  beta = 0.25 
  nocc = 300

  allocate(vector(384)) 

  call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,ham_bml)
  
  !The following Hamiltonian belongs to a water box structure
  !which was precalculated with dftb+
  call bml_read_matrix(ham_bml,'hamiltonian_ortho.mtx')
  !call bml_print_matrix("ham",ham_bml,1,384,1,384) 
  
  open(1, file = 'ham384.txt', status = 'new')
  do j=1,384
        call bml_get_row(ham_bml,j,vector)
        do i=1,384
                write(1,*) vector(i)  
        end do
  end do 

  !Allocate the density matrix
  call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,p_bml)

  call prg_implicit_fermi(ham_bml, xi0_bml, p_bml, nsteps, nocc,  mu, beta, osteps, occErrLimit, threshold)

  !call bml_print_matrix("density matrix",p_bml,0,10,0,10) 

  call bml_deallocate(ham_bml)
  call bml_deallocate(p_bml)
end program main
