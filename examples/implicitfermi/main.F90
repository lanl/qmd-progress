
!!
program main

  use bml

  !progress lib modes
  use prg_implicit_fermi_mod
  !use hamiltonian_mod

  implicit none

  integer, parameter :: dp = 8
  integer :: norb, mdim, verbose, nsteps,osteps 
  type(bml_matrix_t) :: ham_bml, p_bml, p1_bml
  character(20) :: bml_type
  real(dp) :: threshold, mu, beta, occErrLimit, nocc, error  
  !integer :: i,j
  !real,allocatable :: vector(:)
  !real :: bml_get 

  

  !Some parameters that can be changed depending on the test.
  bml_type = "ellpack"
  threshold = 0.000000001
  mdim = -1
  verbose = 1
  norb = 384
  nsteps = 10
  osteps = 1
  occErrLimit = 1.0
  mu = 0.0
  beta = 100 
  nocc = 200

  call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,ham_bml)
  
  !The following Hamiltonian belongs to a water box structure
  !which was precalculated with dftb+
  call bml_read_matrix(ham_bml,'hamiltonian_ortho.mtx')
  call bml_print_matrix("hamiltonian",ham_bml,0,10,0,10)   

  !Allocate the density matrix
  call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,p_bml)
  call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,p1_bml)

  call prg_implicit_fermi(ham_bml, p_bml, nsteps, nocc,  mu, beta, osteps, occErrLimit, threshold)

  call bml_print_matrix("density matrix",p_bml,0,10,0,10) 
  write(*,*) "mu = ", mu 
  mu = 0.0
  call prg_test_density_matrix(ham_bml, p1_bml, beta, mu, nocc, osteps, occErrLimit, threshold)

  call bml_print_matrix("density matrix",p1_bml,0,10,0,10)
  write(*,*) "mu = ", mu
  call bml_add(p_bml, p1_bml, 1.0_dp, -1.0_dp, threshold)

  error = bml_fnorm(p_bml)

  write(*,*) "error =", error

  call bml_deallocate(ham_bml)
  call bml_deallocate(p1_bml)
  call bml_deallocate(p_bml)
end program main
