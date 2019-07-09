! The Implicit Recursive Fermi O(N) module.
!! \ingroup PROGRESS
!! \brief This subroutine implements Niklasson's implicit recursive fermi dirac exact 
!! density matrix purification algorithm.
!!
module prg_implicit_fermi_mod

  use bml   
  use prg_normalize_mod
  use prg_timer_mod
  use prg_parallel_mod

  implicit none

  private  !Everything is private by default

  integer, parameter :: dp = kind(1.0d0)

  public :: prg_implicit_fermi

contains

  !> Recursive Implicit Fermi Dirac.
  !! \param h_bml Input Hamiltonian matrix.
  !! \param xi0_bml Initial guess of first inverse.
  !! \param p_bml Output density matrix.
  !! \param nsteps Number of sp2 iterations.
  !! \param nocc Number of occupied states.
  !! \param mu Shifted chemical potential
  !! \param beta Input inverse temperature.
  !! \param osteps Outer loop steps to converge chemical potential
  !! \param occErrLimit Occupation error limit.
  !! \param threshold Threshold for multiplication.
  !!
  subroutine prg_implicit_fermi(h_bml, xi0_bml, p_bml, nsteps, nocc, &
       mu, beta, osteps, occErrLimit, threshold)

    implicit none

    type(bml_matrix_t), intent(in) :: h_bml
    type(bml_matrix_t), intent(inout) :: xi0_bml, p_bml
    integer, intent(in) :: osteps, nsteps
    real(dp), intent(in) :: nocc, threshold
    real(dp), intent(in) :: occErrLimit, beta
    real(dp), intent(inout) :: mu

    type(bml_matrix_t) :: p2_bml, i_bml, xi_bml
    type(bml_matrix_t) :: xtmp_bml, x_bml, y_bml
    real(dp) :: trdPdmu, trP0, occErr
    real(dp) :: cnst, ofactor
    real(dp), allocatable :: trace(:), gbnd(:)
    character(20) :: bml_type
    integer :: N, M, i, j, iter, NrSchultz
    logical :: firstTime

    bml_type = bml_get_type(h_bml)
    N = bml_get_N(h_bml)
    M = bml_get_M(h_bml)

    allocate(trace(2))

    call bml_identity_matrix(bml_type, bml_element_real, dp, N, M, i_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, xtmp_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, p2_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, xi_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, x_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, y_bml)

    occErr = 1.0_dp
    firstTime = .false.
    iter = 0
    cnst = beta/(2**(2+nsteps))

    do while ((osteps .eq. 0 .and. occErr .gt. occErrLimit) .or. &
         (osteps .gt. 0 .and. iter .lt. osteps))
       iter = iter + 1

       ! Normalization 
       ! P0 = 0.5*II - cnst*(H0-mu0*II)
       call bml_copy(h_bml, p_bml)
       call prg_normalize_implicit_fermi(p_bml, cnst, mu) 

       if (.not. bml_allocated(xi0_bml)) then
          call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, xi0_bml)
          firstTime = .true.      
       endif

       call bml_copy(xi0_bml, xi_bml) !XI = XI0

       do i = 1, nsteps
          call bml_multiply_x2(p_bml, p2_bml, threshold, trace)

          ! Y = 2*(P02-P0) + II
          call bml_copy(p2_bml, y_bml)
          call bml_add(y_bml, p_bml, 1.0_dp, -1.0_dp, threshold)
          call bml_add(y_bml, i_bml, 2.0_dp, 1.0_dp, threshold)

          !! First time do a full inverse
          if (firstTime .eqv. .true. ) then
             call bml_inverse(y_bml, xi_bml)
             firstTime = .false.
          else
             NrSchultz = 4
             do j = 1, NrSchultz
                call bml_copy(xi_bml, xtmp_bml)
                call bml_multiply(xtmp_bml, y_bml, x_bml, -1.0_dp, 0.0_dp)
                call bml_multiply(x_bml, xtmp_bml, xi_bml, 1.0_dp, 2.0_dp) 
             enddo
          endif

          if (i == 1) then
             call bml_copy(xi_bml, xi0_bml)
          endif

          call bml_multiply(xi_bml, p2_bml, p_bml, 1.0_dp, 0.0_dp, threshold)
       enddo

       trdPdmu = bml_trace(p_bml)
       trP0 = trdPdmu
       trdPdmu = trdPdmu - bml_sum_squares(p_bml) ! sum p(i,j)**2

       trdPdmu = beta * trdPdmu
       mu = mu + (nocc - trP0)/trdPdmu
       occErr = abs(trP0 - nocc)
    enddo

    ! Adjust occupation
    !! X = II-P0
    call bml_copy(i_bml, x_bml)
    call bml_add(x_bml, p_bml, 1.0_dp, -1.0_dp, threshold)

    call bml_multiply(p_bml, x_bml, xi_bml, 1.0_dp, 0.0_dp, threshold)
    ofactor = ((nocc - trP0)/trdPdmu) * beta
    call bml_add(p_bml, xi_bml, 1.0_dp, ofactor, threshold)

    deallocate(trace)

    call bml_deallocate(xi_bml)
    call bml_deallocate(i_bml)
    call bml_deallocate(p2_bml)
    call bml_deallocate(xtmp_bml)
    call bml_deallocate(x_bml)
    call bml_deallocate(y_bml)

  end subroutine prg_implicit_fermi




end module prg_implicit_fermi_mod
