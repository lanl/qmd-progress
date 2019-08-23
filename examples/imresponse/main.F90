
!!
program main

  ! Program to test prg_implicit_fermi_response for calculating perturbations in
  ! the density matrix from perturbations in the hamiltonian with implicit fermi 

  use bml

  !progress lib modes
  use prg_implicit_fermi_mod
  !use hamiltonian_mod

  implicit none

  integer, parameter :: dp = 8
  integer :: N, M, order, verbose, nsteps,osteps 
  type(bml_matrix_t) :: H0_bml, H1_bml, H2_bml, H3_bml, P0_bml, P1_bml, P2_bml, P3_bml
  type(bml_matrix_t), allocatable :: H_bml(:) 
  character(20) :: bml_type
  real(dp) :: threshold, mu, beta, occErrLimit, nocc, error, lin_tol, h, lambda
  real(dp), allocatable :: mu_list(:)  

  

  !Some parameters that can be changed depending on the test.
  bml_type = "ellpack"
  threshold = 0.000000000
  verbose = 1
  N = 450
  M = 450
  nsteps = 18
  osteps = 1
  occErrLimit = 0.001
  mu = -0.18 
  beta = 50
  nocc = 100.75
  order = 1
  lin_tol = 0.00000000001
  h = 0.0001
  lambda = 0.001
  allocate(mu_list(order))
  allocate(H_bml(3))
 
  call bml_zero_matrix(bml_type,bml_element_real,dp,N,M,H0_bml)
 
        call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, H1_bml) 
        call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, H2_bml) 
        call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, H3_bml) 
        call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, P0_bml)
        call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, P1_bml)
        call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, P2_bml)
        call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, P3_bml)
 
  call bml_read_matrix(H0_bml,'F_ort_0064.mtx')
  call bml_print_matrix("hamiltonian",H0_bml,0,10,0,10)   

        call bml_copy(H0_bml, H1_bml)
        call bml_scale(0.1_dp, H1_bml)
        call bml_copy(H0_bml, H2_bml)
        call bml_scale(0.1_dp, H2_bml)
        call bml_copy(H0_bml, H3_bml)
        call bml_scale(0.1_dp, H3_bml)

  call prg_implicit_fermi_response(H0_bml, H1_bml, H2_bml, H3_bml, P0_bml, P1_bml, P2_bml, P3_bml, & 
                                   nsteps, mu, mu_list, beta, nocc, occErrLimit, lin_tol, order, threshold)

       call bml_scale(1000.0_dp, P1_bml)
       call bml_scale(1000.0_dp, P2_bml)
       call bml_scale(1000.0_dp, P3_bml)
       call bml_print_matrix("P1", P1_bml, 0,10,0,10)
       call bml_print_matrix("P2", P2_bml, 0,10,0,10)
       call bml_print_matrix("P3", P3_bml, 0,10,0,10)

  H_bml(1) = H1_bml
  H_bml(2) = H2_bml
  H_bml(3) = H3_bml
  call prg_finite_diff(H0_bml, H_bml, mu, mu_list, beta, order, lambda, h, threshold)

  call bml_deallocate(H0_bml)
  call bml_deallocate(H1_bml)
  call bml_deallocate(H2_bml)
  call bml_deallocate(H3_bml)
  call bml_deallocate(P0_bml)
  call bml_deallocate(P1_bml)
  call bml_deallocate(P2_bml)
  call bml_deallocate(P3_bml)
  deallocate(mu_list) 
end program main
