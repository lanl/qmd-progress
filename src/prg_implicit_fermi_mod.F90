! The Implicit Recursive Fermi O(N) module.
!! \ingroup PROGRESS
!! \brief Here are subroutines implementing Niklasson's implicit recursive fermi dirac exact 
!! density matrix purification algorithm.
!!
module prg_implicit_fermi_mod

  use bml   
  use prg_normalize_mod
  use prg_densitymatrix_mod 
  use prg_timer_mod
  use prg_parallel_mod

  implicit none

  private  !Everything is private by default

  integer, parameter :: dp = kind(1.0d0)

  public :: prg_implicit_fermi
  public :: prg_implicit_fermi_zero
  public :: prg_test_density_matrix
  public :: prg_implicit_fermi_response 
  public :: prg_finite_diff  

contains

  !> Recursive Implicit Fermi Dirac for finite temperature.
  !! \param h_bml Input Hamiltonian matrix.
  !! \param p_bml Output density matrix.
  !! \param nsteps Number of recursion steps.
  !! \param k Expansion order 
  !! \param nocc Number of occupied states.
  !! \param mu Shifted chemical potential
  !! \param beta Input inverse temperature.
  !! \param method 0 - conjugate gradient, 1 - newton-schultz
  !! \param osteps Outer loop steps to converge chemical potential
  !! \param occErrLimit Occupation error limit.
  !! \param threshold Threshold for multiplication.
  !! \param tol Tolerance for linear system solver
  subroutine prg_implicit_fermi(h_bml, p_bml, nsteps, k, nocc, &
       mu, beta, method, osteps, occErrLimit, threshold, tol)

    implicit none

    type(bml_matrix_t), intent(in) :: h_bml
    type(bml_matrix_t), intent(inout) :: p_bml
    integer, intent(in) :: osteps, nsteps, method, k
    real(dp), intent(in) :: nocc, threshold
    real(dp), intent(in) :: tol
    real(dp), intent(in) :: occErrLimit, beta
    real(dp), intent(inout) :: mu

    type(bml_matrix_t) :: w_bml, y_bml, d_bml, p2_bml, aux1_bml, aux2_bml, I_bml, ai_bml
    real(dp) :: trdPdmu, trP0, occErr
    real(dp) :: cnst, ofactor
    real(dp), allocatable :: trace(:), gbnd(:)
    character(20) :: bml_type
    integer :: N, M, i, iter, exp_order

    bml_type = bml_get_type(h_bml)
    N = bml_get_N(h_bml)
    M = bml_get_M(h_bml)

    allocate(trace(2))   
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, p2_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, d_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, w_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, y_bml)
    if (k .gt. 2) then 
       call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, aux1_bml)
       call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, aux2_bml)
    endif
    if (method .eq. 1) then
       call bml_identity_matrix(bml_type, bml_element_real, dp, N, M, I_bml)
       call bml_identity_matrix(bml_type, bml_element_real, dp, N, M, ai_bml)
    endif

    occErr = 10000.0_dp
    iter = 0
    exp_order = k**nsteps
    cnst = beta/(4*exp_order)


    do while ((osteps .eq. 0 .and. occErr .gt. occErrLimit) .or. &
         (osteps .gt. 0 .and. iter .lt. osteps))
       iter = iter + 1

       ! Normalization 
       ! P0 = 0.5*II - cnst*(H0-mu0*II)
       call bml_copy(h_bml, p_bml)
       call prg_normalize_implicit_fermi(p_bml, cnst, mu)

       if (method .eq. 0) then   
          write(*,*) "Doing CG"
          do i = 1, nsteps

             if (k .eq. 2) then 
                call bml_multiply_x2(p_bml, p2_bml, threshold, trace)

                ! Y = 2*(P2-P) + II
                call bml_copy(p2_bml, y_bml)
                call bml_add(y_bml, p_bml, 1.0_dp, -1.0_dp, threshold)
                call bml_scale_add_identity(y_bml, 2.0_dp, 1.0_dp, threshold)
             else 
                call prg_setup_linsys(p_bml, y_bml, p2_bml, d_bml, w_bml, aux1_bml, aux2_bml, k, threshold) 
             end if
             call prg_conjgrad(y_bml, p_bml, p2_bml, d_bml, w_bml, tol, threshold)
          enddo
       else 
          write(*,*) "Doing NS"
          do i = 1, nsteps

             if (k .eq. 2) then 
                call bml_multiply_x2(p_bml, p2_bml, threshold, trace)

                ! Y = 2*(P2-P) + II
                call bml_copy(p2_bml, y_bml)
                call bml_add(y_bml, p_bml, 1.0_dp, -1.0_dp, threshold)
                call bml_scale_add_identity(y_bml, 2.0_dp, 1.0_dp, threshold)
             else
                call prg_setup_linsys(p_bml, y_bml, p2_bml, d_bml, w_bml, aux1_bml, aux2_bml, k, threshold)
             end if
             if (i .eq. 1) then 
                call prg_conjgrad(y_bml, ai_bml, I_bml, d_bml, w_bml, 0.9_dp, threshold)
             end if
             call prg_newtonschulz(y_bml, ai_bml, d_bml, w_bml, tol, threshold)
             call bml_multiply(ai_bml, p2_bml, p_bml, 1.0_dp, 0.0_dp, threshold)
          enddo

       end if
       call bml_print_matrix("implicit",p_bml,0,10,0,10)
       trdPdmu = bml_trace(p_bml)
       trP0 = trdPdmu
       trdPdmu = trdPdmu - bml_sum_squares(p_bml) ! sum p(i,j)**2
       trdPdmu = beta * trdPdmu
       mu = mu + (nocc - trP0)/trdPdmu
       occErr = abs(trP0 - nocc)
    enddo

    ! Adjust occupation
    ! X = II-P0
    call bml_copy(p_bml, d_bml)
    call bml_scale_add_identity(d_bml, -1.0_dp, 1.0_dp, threshold)

    call bml_multiply(p_bml, d_bml, w_bml, 1.0_dp, 0.0_dp, threshold)
    ofactor = ((nocc - trP0)/trdPdmu) * beta
    call bml_add(p_bml, w_bml, 1.0_dp, ofactor, threshold)
    !call bml_print_matrix("density matrix",p_bml,0,10,0,10)

    deallocate(trace)

    call bml_deallocate(p2_bml)
    call bml_deallocate(w_bml)
    call bml_deallocate(d_bml)
    call bml_deallocate(y_bml)
    if (k .gt. 2) then 
       call bml_deallocate(aux1_bml)
       call bml_deallocate(aux2_bml) 
    endif
    if (method .eq. 1) then
       call bml_deallocate(ai_bml)
       call bml_deallocate(I_bml)
    endif

  end subroutine prg_implicit_fermi

  !> Recursive Implicit Fermi Dirac for zero temperature.
  !! \param h_bml Input Hamiltonian matrix.
  !! \param p_bml Output density matrix.
  !! \param nsteps Number of recursion steps.
  !! \param mu Shifted chemical potential
  !! \param beta Input inverse temperature.
  !! \param method 0 - conjugate gradient, 1 - newton-schultz
  !! \param threshold Threshold for multiplication.
  !! \param tol Tolerance for linear system solver
  subroutine prg_implicit_fermi_zero(h_bml, p_bml, nsteps, mu, method, threshold, tol)

    implicit none

    type(bml_matrix_t), intent(in) :: h_bml
    type(bml_matrix_t), intent(inout) :: p_bml
    integer, intent(in) :: nsteps, method
    real(dp), intent(in) :: mu, threshold
    real(dp), intent(inout), optional :: tol

    type(bml_matrix_t) :: w_bml, y_bml, d_bml, p2_bml, aux1_bml, aux2_bml, I_bml, ai_bml
    real(dp) :: cnst 
    real(dp), allocatable :: trace(:), gbnd(:)
    character(20) :: bml_type
    integer :: N, M, i 

    bml_type = bml_get_type(h_bml)
    N = bml_get_N(h_bml)
    M = bml_get_M(h_bml)

    allocate(trace(2))
    allocate(gbnd(2))
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, p2_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, d_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, w_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, y_bml)
    if (method .eq. 1) then
       call bml_identity_matrix(bml_type, bml_element_real, dp, N, M, I_bml)
       call bml_identity_matrix(bml_type, bml_element_real, dp, N, M, ai_bml)
    endif

    call bml_copy(h_bml, p_bml)
    call bml_gershgorin(p_bml, gbnd) 
    cnst = 0.5*min(1/(mu-gbnd(1)),1/(gbnd(2)-mu))
    call prg_normalize_implicit_fermi(p_bml, cnst, mu)

    if (method .eq. 0) then
       write(*,*) "Doing CG"
       do i = 1, nsteps

          call bml_multiply_x2(p_bml, p2_bml, threshold, trace)

          ! Y = 2*(P2-P) + II
          call bml_copy(p2_bml, y_bml)
          call bml_add(y_bml, p_bml, 1.0_dp, -1.0_dp, threshold)
          call bml_scale_add_identity(y_bml, 2.0_dp, 1.0_dp, threshold)
          call prg_conjgrad(y_bml, p_bml, p2_bml, d_bml, w_bml, tol, threshold)
       enddo
    else
       write(*,*) "Doing NS"
       do i = 1, nsteps

          ! Y = 2*(P2-P) + II
          call bml_copy(p2_bml, y_bml)
          call bml_add(y_bml, p_bml, 1.0_dp, -1.0_dp, threshold)
          call bml_scale_add_identity(y_bml, 2.0_dp, 1.0_dp, threshold)
          if (i .eq. 1) then
             call prg_conjgrad(y_bml, ai_bml, I_bml, d_bml, w_bml, 0.9_dp, threshold)
          end if
          call prg_newtonschulz(y_bml, ai_bml, d_bml, w_bml, tol, threshold)
          call bml_multiply(ai_bml, p2_bml, p_bml, 1.0_dp, 0.0_dp, threshold)
       enddo
    endif

    deallocate(gbnd)
    deallocate(trace)

    call bml_deallocate(p2_bml)
    call bml_deallocate(w_bml)
    call bml_deallocate(d_bml)
    call bml_deallocate(y_bml)
    if (method .eq. 1) then
       call bml_deallocate(ai_bml)
       call bml_deallocate(I_bml)
    endif

  end subroutine prg_implicit_fermi_zero

  !> Recursive Implicit Fermi Dirac for finite temperature.
  !! \param H0_bml Input Hamiltonian matrix.
  !! \param H1_bml, H2_bml, H3_bml Input First to third order perturbations of H0.
  !! \param P0_bml Output density matrix.
  !! \param P1_bml, P2_bml, P3_bml Output First to third order density matrix response.   
  !! \param nsteps Number of recursion steps.
  !! \param mu0 Shifted chemical potential.
  !! \param mu Pre-allocated array of length order.
  !! \param beta Input inverse temperature.  
  !! \param nocc Number of occupied states.
  !! \param occ_tol Occupation error tolerance.
  !! \param lin_tol Linear solver tolerance.
  !! \param order Calculate response up to this order.
  !! \param threshold Threshold for matrix algebra. 
  subroutine prg_implicit_fermi_response(H0_bml, H1_bml, H2_bml, H3_bml, P0_bml, P1_bml, P2_bml, P3_bml, & 
       nsteps, mu0, mu, beta, nocc, occ_tol, lin_tol, order, threshold)

    implicit none 

    type(bml_matrix_t), intent(in) :: H0_bml, H1_bml, H2_bml, H3_bml
    type(bml_matrix_t), intent(inout) :: P0_bml, P1_bml, P2_bml, P3_bml
    real(dp), intent(inout) :: mu0 
    real(dp), allocatable, intent(inout) :: mu(:)
    real(dp), intent(in) :: beta, occ_tol, lin_tol, nocc
    integer, intent(in) :: nsteps
    type(bml_matrix_t) :: I_bml, tmp1_bml, tmp2_bml, C0_bml, T_bml, Ti_bml
    type(bml_matrix_t), allocatable :: B_bml(:), P_bml(:), C_bml(:), H_bml(:) 
    real(dp), allocatable :: p_trace(:), trace(:)
    character(20) :: bml_type
    real(dp) :: occ_err, p0_trace, pmu_trace, cnst, threshold, tol, lambda, h  
    integer :: N, M, order, i, j, k

    k = 0
    occ_err = 10000.0
    allocate(p_trace(order))
    allocate(B_bml(order))
    allocate(C_bml(order))
    allocate(P_bml(order))
    allocate(H_bml(order))

    bml_type = bml_get_type(H0_bml)
    N = bml_get_N(H0_bml)
    M = bml_get_M(H0_bml) 

    do i = 1, order 
       mu(i) = 0.0_dp
       call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, B_bml(i))
       call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, C_bml(i))
    end do

    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, tmp1_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, tmp2_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, C0_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, T_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, Ti_bml)
    call bml_identity_matrix(bml_type, bml_element_real, dp, N, M, I_bml)

    H_bml(1) = H1_bml
    P_bml(1) = P1_bml
    if (order .gt. 1) then
       H_bml(2) = H2_bml
       P_bml(2) = P2_bml
    end if
    if (order .gt. 2) then
       H_bml(3) = H3_bml
       P_bml(3) = P3_bml
    end if

    cnst = beta/(2**(2+nsteps))

    do while (occ_err .gt. occ_tol)
       k = k + 1
       call bml_copy(H0_bml, P0_bml)
       call prg_normalize_implicit_fermi(P0_bml, cnst, mu0)

       do j = 1, order
          call bml_copy(H_bml(j), P_bml(j))
          call prg_normalize_implicit_fermi(P_bml(j), cnst, mu(j))
          call bml_scale_add_identity(P_bml(j), 1.0_dp, -0.5_dp, threshold) 
       end do

       do i = 1, nsteps

          ! Calculate coefficient matrices 

          ! C0 = P0^2
          call bml_multiply(P0_bml, P0_bml, C0_bml, 1.0_dp, 0.0_dp, threshold)
          ! C1 = P0*P1+P1*P0, B1 = 2(P1 - C1)
          call bml_multiply(P0_bml, P_bml(1), C_bml(1), 1.0_dp, 0.0_dp, threshold)
          call bml_multiply(P_bml(1), P0_bml, C_bml(1), 1.0_dp, 1.0_dp, threshold)
          call bml_copy(P_bml(1), B_bml(1)) 
          call bml_add(B_bml(1), C_bml(1), 2.0_dp, -2.0_dp, threshold)                
          if (order > 1) then
             ! C2 = P1^2 + P0*P2 + P2*P0, B2 = 2(P2 - C2)
             call bml_multiply(P_bml(1), P_bml(1), C_bml(2), 1.0_dp, 0.0_dp, threshold)
             call bml_multiply(P0_bml, P_bml(2), C_bml(2), 1.0_dp, 1.0_dp, threshold)
             call bml_multiply(P_bml(2), P0_bml, C_bml(2), 1.0_dp, 1.0_dp, threshold)
             call bml_copy(P_bml(2), B_bml(2))
             call bml_add(B_bml(2), C_bml(2), 2.0_dp, -2.0_dp, threshold)
          end if
          if (order > 2) then 
             ! C3 = P1*P2 + P2+P1 + P0*P3 + P3*P0, B3 = 2(P3 - C3) 
             call bml_multiply(P_bml(1), P_bml(2), C_bml(3), 1.0_dp, 0.0_dp, threshold)
             call bml_multiply(P_bml(2), P_bml(1), C_bml(3), 1.0_dp, 1.0_dp, threshold)
             call bml_multiply(P0_bml, P_bml(3), C_bml(3), 1.0_dp, 1.0_dp, threshold)
             call bml_multiply(P_bml(3), P0_bml, C_bml(3), 1.0_dp, 1.0_dp, threshold)
             call bml_copy(P_bml(3), B_bml(3))
             call bml_add(B_bml(3), C_bml(3), 2.0_dp, -2.0_dp, threshold)
          endif
          ! T = 2P0^2 - 2P0 + I
          call bml_copy(C0_bml, T_bml)
          call bml_add(T_bml, P0_bml, 1.0_dp, -1.0_dp, threshold)
          call bml_scale_add_identity(T_bml, 2.0_dp, 1.0_dp, threshold)
          ! Find T-inverse 
          if (i .eq. 1) then
             call prg_conjgrad(T_bml, Ti_bml, I_bml, tmp1_bml, tmp2_bml, 0.01_dp, threshold)
             call bml_identity_matrix(bml_type, bml_element_real, dp, N, M, I_bml)
          end if
          call prg_newtonschulz(T_bml, Ti_bml, tmp1_bml, tmp2_bml, lin_tol, threshold)
          ! Get next P0
          call bml_multiply(Ti_bml, C0_bml, P0_bml, 1.0_dp, 0.0_dp, threshold)
          ! Get next P1
          call bml_multiply(B_bml(1), P0_bml, C_bml(1), 1.0_dp, 1.0_dp, threshold)
          call bml_multiply(Ti_bml, C_bml(1), P_bml(1), 1.0_dp, 0.0_dp, threshold)
          if (order > 1) then 
             ! Get next P2 
             call bml_multiply(B_bml(2), P0_bml, C_bml(2), 1.0_dp, 1.0_dp, threshold)
             call bml_multiply(B_bml(1), P_bml(1), C_bml(2), 1.0_dp, 1.0_dp, threshold)
             call bml_multiply(Ti_bml, C_bml(2), P_bml(2), 1.0_dp, 0.0_dp, threshold)
          end if
          if (order > 2) then 
             ! Get next P3 
             call bml_multiply(B_bml(3), P0_bml, C_bml(3), 1.0_dp, 1.0_dp, threshold)
             call bml_multiply(B_bml(2), P_bml(1), C_bml(3), 1.0_dp, 1.0_dp, threshold)
             call bml_multiply(B_bml(1), P_bml(2), C_bml(3), 1.0_dp, 1.0_dp, threshold)
             call bml_multiply(Ti_bml, C_bml(3), P_bml(3), 1.0_dp, 0.0_dp, threshold)
          endif
       enddo

       ! Pmu = beta*P0(I-P0)
       call bml_copy(P0_bml, tmp1_bml)
       call bml_scale_add_identity(tmp1_bml, -1.0_dp, 1.0_dp, threshold)
       call bml_multiply(P0_bml, tmp1_bml, tmp2_bml, beta, 0.0_dp, threshold)

       pmu_trace = bml_trace(tmp2_bml)
       p0_trace = bml_trace(P0_bml)
       occ_err = abs(p0_trace-nocc)
       mu0 = mu0 + (nocc - p0_trace)/pmu_trace
       do i = 1, order
          p_trace(i) = bml_trace(P_bml(i))
          mu(i) = mu(i) - p_trace(i)/pmu_trace
          occ_err = occ_err + abs(p_trace(i))
       enddo

       write(*,*) "occ_err =", occ_err
       if (k .gt. 50) then 
          write(*,*) "Chemical potential is not converging"
          exit 
       endif

    enddo

    do i = 1, order
       call bml_deallocate(B_bml(i))
       call bml_deallocate(C_bml(i))
    end do

    call bml_deallocate(I_bml)
    call bml_deallocate(tmp1_bml)
    call bml_deallocate(tmp2_bml)
    call bml_deallocate(T_bml)
    call bml_deallocate(Ti_bml)
    deallocate(p_trace)
    deallocate(B_bml)
    deallocate(C_bml)

  end subroutine prg_implicit_fermi_response

  !> Calculate density matrix response from perturbations in the Hamiltonian
  !using finite differences.
  !! \param H0_bml Input Hamiltonian matrix.
  !! \param H_list Input List of one to third order Hamiltonian perturbations 
  !! \param mu0 Shifted chemical potential.
  !! \param mu List of first to third order perturbations in the chemical
  !!  potential.
  !! \param beta Input inverse temperature.  
  !! \param order Calculate response up to this order.
  !! \param lambda Perturbation parameter
  !! \param h Finite difference step size 
  !! \param threshold Threshold for matrix algebra. 
  subroutine prg_finite_diff(H0_bml, H_list, mu0, mu_list, beta, order, lambda, h, threshold)

    implicit none

    type(bml_matrix_t), intent(in) :: H0_bml
    real(dp), intent(in) :: mu0
    type(bml_matrix_t), allocatable, intent(in) :: H_list(:)
    real(dp), allocatable, intent(in) :: mu_list(:)
    real(dp), intent(in) :: lambda, beta, threshold, h 
    integer, intent(in) :: order 
    character(20) :: bml_type
    real(dp) :: mu_1minus, mu_1plus, mu_2minus, mu_2plus, mu_central
    real(dp) :: lambda_f, lambda_b, lambda_2f, lambda_2b
    integer :: N, M, i 
    type(bml_matrix_t) :: D0_bml, D1minus_bml, D1plus_bml, D2plus_bml, D2minus_bml, tmp1_bml 

    bml_type = bml_get_type(H0_bml)
    N = bml_get_N(H0_bml)
    M = bml_get_M(H0_bml)

    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, tmp1_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, D0_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, D1minus_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, D1plus_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, D2plus_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, D2minus_bml)

    lambda_f = lambda+h
    lambda_b = lambda-h
    lambda_2f = lambda+2*h
    lambda_2b = lambda-2*h

    call bml_copy(H0_bml, tmp1_bml) 
    mu_central = mu0
    do i = 1, order  
       call bml_add(tmp1_bml, H_list(i), 1.0_dp, lambda**i, threshold)
       mu_central = mu_central + lambda**i*mu_list(i) 
    end do
    call prg_get_density_matrix(tmp1_bml, D0_bml, beta, mu_central, threshold)

    call bml_copy(H0_bml, tmp1_bml)
    mu_1plus = mu0
    do i = 1, order 
       call bml_add(tmp1_bml, H_list(i), 1.0_dp, lambda_f**i, threshold)
       mu_1plus = mu_1plus + lambda_f**i*mu_list(i) 
    end do
    call prg_get_density_matrix(tmp1_bml, D1plus_bml, beta, mu_1plus, threshold)
    call bml_copy(D1plus_bml, tmp1_bml)
    call bml_add(tmp1_bml, D0_bml, 1.0_dp/h, -1.0_dp/h, threshold)

    call bml_scale(1000.0_dp, tmp1_bml)
    call bml_print_matrix("Finite diff - Order 1 * 1000", tmp1_bml, 0,10,0,10)

    if (order .gt. 1) then 
       call bml_copy(H0_bml, tmp1_bml)
       mu_1minus = mu0
       do i = 1, order
          call bml_add(tmp1_bml, H_list(i), 1.0_dp, lambda_b**i, threshold)
          mu_1minus = mu_1minus + lambda_b**i*mu_list(i)
       end do
       call prg_get_density_matrix(tmp1_bml, D1minus_bml, beta, mu_1minus, threshold)
       call bml_copy(D0_bml, tmp1_bml)
       call bml_add(tmp1_bml, D1minus_bml, -1.0_dp/(h*h), 0.5_dp/(h*h), threshold)
       call bml_add(tmp1_bml, D1plus_bml, 1.0_dp, 0.5_dp/(h*h))

       call bml_scale(1000.0_dp, tmp1_bml)
       call bml_print_matrix("Finite diff - Order 2 * 1000", tmp1_bml, 0,10,0,10)          
    end if

    if (order .gt. 2) then 
       call bml_copy(H0_bml, tmp1_bml)
       mu_2minus = mu0
       do i = 1, order
          call bml_add(tmp1_bml, H_list(i), 1.0_dp, lambda_2b**i, threshold)
          mu_2minus = mu_2minus + lambda_2b**i*mu_list(i)
       end do
       call prg_get_density_matrix(tmp1_bml, D2minus_bml, beta, mu_2minus, threshold)
       call bml_copy(H0_bml, tmp1_bml)
       mu_2plus = mu0
       do i = 1, order
          call bml_add(tmp1_bml, H_list(i), 1.0_dp, lambda_2f**i, threshold)
          mu_2plus = mu_2plus + lambda_2f**i*mu_list(i)
       end do
       call prg_get_density_matrix(tmp1_bml, D2plus_bml, beta, mu_2plus, threshold)
       call bml_copy(D2plus_bml, tmp1_bml)
       call bml_add(tmp1_bml, D2minus_bml, 1.0_dp/(12.0_dp*h**3), -1.0_dp/(12.0_dp*h**3), threshold)
       call bml_add(tmp1_bml, D1plus_bml, 1.0_dp, -1.0/(6.0_dp*h**3), threshold)
       call bml_add(tmp1_bml, D1minus_bml, 1.0_dp, 1.0/(6.0_dp*h**3), threshold)

       call bml_scale(1000.0_dp, tmp1_bml)
       call bml_print_matrix("Finite diff - Order 3 * 1000", tmp1_bml, 0,10,0,10)
    end if

    call bml_deallocate(tmp1_bml)
    call bml_deallocate(D0_bml)
    call bml_deallocate(D1minus_bml)
    call bml_deallocate(D1plus_bml)
    call bml_deallocate(D2plus_bml)
    call bml_deallocate(D2minus_bml)    

  end subroutine prg_finite_diff

  subroutine prg_setup_linsys(p_bml, A_bml, b_bml, p2_bml, y_bml, aux_bml, &
       aux1_bml, k, threshold) 

    implicit none 

    type(bml_matrix_t), intent(inout) :: A_bml, b_bml, p2_bml, y_bml, aux_bml, aux1_bml
    type(bml_matrix_t), intent(in) :: p_bml
    real(dp), intent(in) :: threshold
    integer, intent(in) :: k 
    character(20) :: bml_type
    integer :: M, N, i

    if (k .eq. 2) then 
       call bml_multiply(p_bml, p_bml, b_bml, 1.0_dp, 0.0_dp, threshold)
       call bml_copy(b_bml, A_bml)
       call bml_add(A_bml, p_bml, 2.0_dp, -2.0_dp, threshold)
       call bml_scale_add_identity(A_bml, 1.0_dp, 1.0_dp, threshold) 

    else 
       call bml_multiply(p_bml, p_bml, p2_bml, 1.0_dp, 0.0_dp, threshold)
       call bml_copy(p2_bml, y_bml)
       call bml_add(y_bml, p_bml, 1.0_dp, -2.0_dp, threshold)
       call bml_scale_add_identity(y_bml, 1.0_dp, 1.0_dp, threshold)

       call bml_copy(p2_bml, b_bml)
       call bml_copy(y_bml, A_bml) 
       do i = 1,(k-2)/2 
          call bml_multiply(b_bml, p2_bml, aux_bml, 1.0_dp, 0.0_dp, threshold)
          call bml_multiply(A_bml, y_bml, aux1_bml, 1.0_dp, 0.0_dp, threshold) 
          call bml_copy(aux_bml, b_bml)
          call bml_copy(aux1_bml, A_bml)
       enddo
       call bml_add(A_bml, b_bml, 1.0_dp, 1.0_dp, threshold)
    end if

  end subroutine prg_setup_linsys

  !> Find the inverse of the matrix A with Newton-Schulz iteration
  !! \param a_bml Input matrix A
  !! \param ai_bml Input starting guess and output inverse 
  !! \param r_bml Auxillary matrix 
  !! \param tmp_bml Auxillary matrix 
  !! \param tol Convergence criterion (Frobenius norm of residual matrix) 
  !! \param threshold Threshold for matrix algebra 

  subroutine prg_newtonschulz(a_bml, ai_bml, r_bml, tmp_bml, tol, threshold)

    implicit none

    type(bml_matrix_t), intent(inout) :: ai_bml, r_bml, tmp_bml
    type(bml_matrix_t), intent(in) :: a_bml
    real(dp), intent(in) :: threshold, tol
    real(dp) :: norm
    integer :: i     

    norm = 1.0
    i = 0
    do while(norm > tol)
       call bml_copy(ai_bml, tmp_bml)
       call bml_multiply(a_bml, ai_bml, r_bml, 1.0_dp, 0.0_dp, threshold)
       call bml_scale_add_identity(r_bml, -1.0_dp, 1.0_dp, threshold)
       norm = bml_fnorm(r_bml)
       !   write(*,*) "norm = ", norm
       if (norm < tol) then
          exit 
       end if
       call bml_multiply(tmp_bml, r_bml, ai_bml, 1.0_dp, 1.0_dp, threshold)
       i = i + 1
    enddo
    ! write(*,*) "Number of NS iterations:", i
  end subroutine prg_newtonschulz

  ! Preconditioned CG, preconditioner inverse diagonal of A
  !> Solve the system AX = B with conjugate gradient 
  !! \param A_bml Coefficient matrix A 
  !! \param p_bml Output solution X
  !! \param p2_bml Right side matrix B
  !! \param d_bml Auxillary matrix 
  !! \param w_bml Auxillary matrix 
  !! \param cg_tol Convergence condition (squared Frobenius norm of residual
  !! matrix)
  !! \param threshold Threshold for matrix algebra 
  subroutine prg_pcg(A_bml, p_bml, p2_bml, d_bml, wtmp_bml, cg_tol, threshold)

    implicit none

    type(bml_matrix_t), intent(in) :: A_bml
    type(bml_matrix_t), intent(inout) :: p_bml, p2_bml, d_bml, wtmp_bml
    real(dp), intent(in) :: cg_tol, threshold

    type(bml_matrix_t) :: M_bml, z_bml
    real(dp), allocatable :: diagonal(:)
    real(dp) :: alpha, beta
    character(20) :: bml_type
    integer :: k,N,M
    real(dp) :: r_norm_old, r_norm_new

    bml_type = bml_get_type(p_bml)
    N = bml_get_N(p_bml)
    M = bml_get_M(p_bml)

    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, z_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, M_bml)

    allocate(diagonal(N))
    call bml_get_diagonal(A_bml, diagonal)
    do k = 1,N
       diagonal(k) = 1.0_dp/diagonal(k)
    enddo
    call bml_set_diagonal(M_bml, diagonal)

    call bml_multiply(A_bml, p_bml, p2_bml, -1.0_dp, 1.0_dp, threshold)
    call bml_multiply(M_bml, p2_bml, z_bml, 1.0_dp, 0.0_dp, threshold) 
    r_norm_new = bml_trace_mult(z_bml, p2_bml)
    call bml_copy(z_bml, d_bml)
    k = 0

    do while (bml_sum_squares(p2_bml) .gt. cg_tol)

       write(*,*) "r_norm", bml_sum_squares(p2_bml)
       k = k + 1
       if (k .ne. 1) then 
          beta = r_norm_new/r_norm_old
          call bml_add(d_bml, z_bml, beta, 1.0_dp, threshold)
       endif

       call bml_multiply(A_bml, d_bml, wtmp_bml, 1.0_dp, 0.0_dp, threshold)
       alpha = bml_trace_mult(p2_bml,z_bml)/bml_trace_mult(d_bml, wtmp_bml)
       call bml_add(p_bml, d_bml, 1.0_dp, alpha, threshold)
       call bml_add(p2_bml, wtmp_bml, 1.0_dp, -alpha, threshold)
       call bml_multiply(M_bml, p2_bml, z_bml, 1.0_dp, 0.0_dp, threshold)
       r_norm_old = r_norm_new
       r_norm_new = bml_trace_mult(p2_bml,z_bml)
       if (k .gt. 100) then
          write(*,*) "PCG is not converging"
          stop
       endif
    enddo
    write(*,*) "Number of iterations:", k

    call bml_deallocate(z_bml)
    call bml_deallocate(M_bml)
    deallocate(diagonal)

  end subroutine prg_pcg

  !> Solve the system AX = B with conjugate gradient 
  !! \param A_bml Coefficient matrix A 
  !! \param p_bml Output solution X
  !! \param p2_bml Right side matrix B
  !! \param d_bml Auxillary matrix 
  !! \param w_bml Auxillary matrix 
  !! \param cg_tol Convergence condition (squared Frobenius norm of residual matrix)
  !! \param threshold Threshold for matrix algebra 
  subroutine prg_conjgrad(A_bml, p_bml, p2_bml, d_bml, w_bml, cg_tol, threshold)

    implicit none

    type(bml_matrix_t), intent(in) :: A_bml
    type(bml_matrix_t), intent(inout) :: p_bml, p2_bml, d_bml, w_bml
    real(dp), intent(in) :: cg_tol, threshold

    real(dp) :: alpha, beta
    integer :: k
    real(dp) :: r_norm_old, r_norm_new

    call bml_multiply(A_bml, p_bml, p2_bml, -1.0_dp, 1.0_dp, threshold)
    r_norm_new = bml_sum_squares(p2_bml)
    k = 0

    do while (r_norm_new .gt. cg_tol) 

       !   write(*,*) r_norm_new
       k = k + 1
       if (k .eq. 1) then 
          write(*,*) r_norm_new 
          call bml_copy(p2_bml, d_bml)
       else 
          beta = r_norm_new/r_norm_old
          call bml_add(d_bml, p2_bml, beta, 1.0_dp, threshold)
       endif

       call bml_multiply(A_bml, d_bml, w_bml, 1.0_dp, 0.0_dp, threshold)
       alpha = r_norm_new/bml_trace_mult(d_bml, w_bml)

       call bml_add(p_bml, d_bml, 1.0_dp, alpha, threshold)
       call bml_add(p2_bml, w_bml, 1.0_dp, -alpha, threshold)
       r_norm_old = r_norm_new
       r_norm_new = bml_sum_squares(p2_bml)
       if (k .gt. 50) then
          write(*,*) "Conjugate gradient is not converging"
          stop
       endif
    enddo
    write(*,*) "Number of CG-iterations:", k

  end subroutine prg_conjgrad

  subroutine prg_get_density_matrix(ham_bml, p_bml, beta, mu, threshold)

    implicit none

    type(bml_matrix_t), intent(in) :: ham_bml
    type(bml_matrix_t), intent(inout) :: p_bml
    real(dp), intent(in) ::  beta, threshold
    real(dp), intent(in) :: mu
    character(20) :: bml_type
    integer :: N, M, i
    real(dp), allocatable ::  eigenvalues(:)
    type(bml_matrix_t) :: eigenvectors_bml,occupation_bml,aux_bml,aux1_bml,i_bml

    bml_type = bml_get_type(p_bml)
    N = bml_get_N(p_bml)
    M = bml_get_M(p_bml)

    allocate(eigenvalues(N))

    call bml_zero_matrix(bml_type,bml_element_real,dp,N,M,eigenvectors_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,N,M,occupation_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,N,M,aux_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,N,M,aux1_bml)
    call bml_identity_matrix(bml_type,bml_element_real,dp,N,M,i_bml)

    call bml_diagonalize(ham_bml,eigenvalues,eigenvectors_bml)

    do i=1,N
       eigenvalues(i) = fermi(eigenvalues(i),mu,beta)
    enddo

    call bml_set_diagonal(occupation_bml, eigenvalues)
    call bml_multiply(eigenvectors_bml, occupation_bml, aux_bml, 1.0_dp, 0.0_dp,threshold)
    call bml_transpose(eigenvectors_bml, aux1_bml)
    call bml_multiply(aux_bml, aux1_bml, p_bml, 1.0_dp, 0.0_dp, threshold)

    call bml_deallocate(eigenvectors_bml)
    call bml_deallocate(occupation_bml)
    call bml_deallocate(aux_bml)
    call bml_deallocate(aux1_bml)
    call bml_deallocate(i_bml)

    deallocate(eigenvalues)

  end subroutine prg_get_density_matrix

  !> Calculate the density matrix with diagonalization
  !! \param ham_bml Input hamiltonian
  !! \param p_bml Output density matrix
  !! \param beta Inverse temperature 
  !! \param mu Chemical potential 
  !! \param nocc Number of occupied states
  !! \param osteps Outer loop steps to converge chemical potential
  !! \param occErrLimit Occupation error limit.
  !! \param threshold Threshold for matrix algebra 
  subroutine prg_test_density_matrix(ham_bml, p_bml, beta, mu, nocc, osteps, occErrLimit, threshold)

    implicit none 

    type(bml_matrix_t), intent(in) :: ham_bml
    type(bml_matrix_t), intent(inout) :: p_bml
    real(dp), intent(in) ::  beta, nocc, occErrLimit, threshold
    real(dp), intent(inout) :: mu
    integer, intent(in) :: osteps
    character(20) :: bml_type
    integer :: N, M, i, iter
    real(dp) :: trdPdmu, trP0, ofactor, occErr
    real(dp), allocatable ::  eigenvalues(:)
    type(bml_matrix_t) :: eigenvectors_bml,occupation_bml,aux_bml,aux1_bml,i_bml

    bml_type = bml_get_type(p_bml)
    N = bml_get_N(p_bml)
    M = bml_get_M(p_bml)

    allocate(eigenvalues(N))

    call bml_zero_matrix(bml_type,bml_element_real,dp,N,M,eigenvectors_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,N,M,occupation_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,N,M,aux_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,N,M,aux1_bml)
    call bml_identity_matrix(bml_type,bml_element_real,dp,N,M,i_bml)

    occErr = 1000.0_dp
    iter = 0

    do while ((osteps .eq. 0 .and. occErr .gt. occErrLimit) .or. &
         (osteps .gt. 0 .and. iter .lt. osteps))
       iter = iter + 1

       call bml_diagonalize(ham_bml,eigenvalues,eigenvectors_bml)

       do i=1,N   
          eigenvalues(i) = fermi(eigenvalues(i),mu,beta)
       enddo

       call bml_set_diagonal(occupation_bml, eigenvalues) 
       call bml_multiply(eigenvectors_bml, occupation_bml, aux_bml, 1.0_dp, 0.0_dp,threshold)
       call bml_transpose(eigenvectors_bml, aux1_bml)
       call bml_multiply(aux_bml, aux1_bml, p_bml, 1.0_dp, 0.0_dp, threshold)
       trdPdmu = bml_trace(p_bml)
       trP0 = trdPdmu
       trdPdmu = trdPdmu - bml_sum_squares(p_bml) ! sum p(i,j)**2
       trdPdmu = beta * trdPdmu
       mu = mu + (nocc - trP0)/trdPdmu
       occErr = abs(trP0 - nocc)
    enddo

    ! Adjust occupation
    ! X = II-P0
    call bml_copy(p_bml, aux_bml)
    call bml_scale_add_identity(aux_bml, -1.0_dp, 1.0_dp, threshold)

    call bml_multiply(p_bml, aux_bml, aux1_bml, 1.0_dp, 0.0_dp, threshold)
    ofactor = ((nocc - trP0)/trdPdmu) * beta
    call bml_add(p_bml, aux1_bml, 1.0_dp, ofactor, threshold)
    !call bml_print_matrix("density matrix",p_bml,0,10,0,10)

    call bml_deallocate(eigenvectors_bml)
    call bml_deallocate(occupation_bml)
    call bml_deallocate(aux_bml)
    call bml_deallocate(aux1_bml)
    call bml_deallocate(i_bml)

    deallocate(eigenvalues)

  end subroutine prg_test_density_matrix

  !> Gives the Fermi distribution value for energy e.
  !! \param e Energy.
  !! \param mu Fermi energy.
  !!
  real(dp) function fermi(e,mu,beta)

    real(dp), intent(in) :: e, mu, beta

    fermi = 1.0_dp/(1.0_dp+exp(beta*(e-mu)))

  end function fermi

end module prg_implicit_fermi_mod
