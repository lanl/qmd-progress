!> The SP2 Recursive Fermi O(N) module.
!! \ingroup PROGRESS
    !
    ! This subroutine implements Niklasson's SP2 recursive fermi dirac exact 
    ! density matrix purification algorithm.
    !

module prg_sp2_fermirec_mod

  use bml   
  use prg_normalize_mod
  use prg_timer_mod
  use prg_parallel_mod
  
  implicit none
 
  private  !Everything is private by default
  
  integer, parameter :: dp = kind(1.0d0)

  public :: prg_sp2_fermirec

contains
 
  !! What is MaxIt for?

  !> Recursive SP2 Fermi.
  !! \param h_bml Input Hamiltonian matrix.
  !! \param xio_bml Input ? matrix.
  !! \param nsteps Number of sp2 iterations.
  !! \param nocc Number of occupied states.
  !! \param threshold Threshold for multiplication.
  !! \param occErrLimit Occupation error limit.
  !! \param traceLimit Trace limit.
  !! \param x_bml Output initial matrix.
  !! \param mu Shifted chemical potential
  !! \param beta Output inverse temperature.
  subroutine prg_sp2_fermirec(h_bml, xio_bml, nsteps, nocc, threshold, occErrLimit, &
      x_bml, mu, beta)

    implicit none

    type(bml_matrix_t), intent(in) :: h_bml
    type(bml_matrix_t), intent(inout) :: x_bml
    integer, intent(in) :: nsteps
    real(dp), intent(in) :: nocc, threshold
    real(dp), intent(in) :: occErrLimit
    real(dp), intent(inout) :: mu, beta

    type(bml_matrix_t) :: x1_bml, x2_bml, i_bml, xi_bml
    type(bml_matrix_t) :: y_bml
    real(dp) :: trdPdmu, trP0, occErr
    real(dp) :: cnst, ofactor
    real(dp), allocatable :: trace(:)
    character(20) :: bml_type
    integer :: N, M, i, cnt, NrSchultz
    logical :: firstTime

    bml_type = bml_get_type(h_bml)
    N = bml_get_N(h_bml)
    M = bml_get_M(h_bml)

    allocate(trace(2))

    call bml_identity_matrix(bml_type, bml_element_real, dp, N, M, i_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, x1_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, x2_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, xi_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, y_bml)

    occErr = 1.0_dp
    cnt = 0
    firstTime = .true.

    do while (occErr .gt. occErrLimit)

      ! calculate beta before entering this subroutine?
      !beta = 1.0_dp/(kB*T) ! Temp in K

      ! Is this a normalization? 
      ! yes, make as separate routine in normalize_mod
      ! P0 = 0.5*II - cnst*(H0-mu0*II)
      cnst = beta/(2**(2+nsteps))
      call bml_copy(h_bml, x_bml)
      call bml_add_deprecated(1.0_dp, x_bml, -mu, i_bml, threshold)
      call bml_add_deprecated(-cnst, x_bml, 0.5_dp, i_bml, threshold)
      call bml_copy(xio_bml, xi_bml) !XI = XIO

      do i = 1, nsteps
        call bml_multiply_x2(x_bml, x2_bml, threshold, trace)
        ! Y = 2*(P02-P0) + II
        call bml_copy(x2_bml, y_bml)
        call bml_add_deprecated(1.0_dp, y_bml, -1.0_dp, x_bml, threshold)
        call bml_add_deprecated(2.0_dp, y_bml, 1.0_dp, i_bml, threshold)

        if (firstTime .eqv. .true.) then
          call bml_inverse(y_bml, xi_bml)
        else
          NrSchulz = 4
          do i = 1, NrSchulz
            call bml_copy(xi_bml, xtmp_bml)
            call bml_multiply(xtmp_bml, y_bml, x_bml, &
              -1.0_dp, 0.0_dp)
            call bml_multiply(x_bml, xtmp_bml, xi_bml, &
               1.0_dp, 2.0_dp) 
          enddo
        endif
        firsttime = .false.


        if (j == 1) then
          call bml_copy(xi_bml, xio_bml)
        endif

        call bml_multiply(xi_bml, x2_bml, x_bml, 1.0_dp, 0.0_dp, threshold)
      enddo

      ! Need to write sumDiagonal routine?
      trdPdmu = bml_sumDiagonal(x_bml)
      trP0 = trdPdmu
      ! Need to write sumX2 routine?
      trdPdmu = trdPdmu - bml_sumX2(x_bml) ! sum x(i,j)**2

      trdPdmu = beta * trdPdmu
      mu = mu + (nocc - trP0)/trdPdmu
      occErr = abs(trP0 - nocc)
      cnt = cnt + 1
      if (cnt .ge. MaxIt) then
        occErr = 0.0_dp
      endif
    enddo

    ! Adjust occupation
    call bml_copy(i_bml, x1_bml)
    call bml_add_deprecated(1.0_dp, x1_bml, -1.0_dp, x_bml, threshold)
    call bml_multiply(x_bml, x1_bml, xi_bml, 1.0_dp, 0.0_dp, threshold)
    ofactor = ((nocc - trP0)/trdPmu) * beta
    call bml_add_deprecated(1.0_dp, x_bml, ofactor, xi_bml, threshold)

    deallocate(trace)

    call bml_deallocate(xi_bml)
    call bml_deallocate(i_bml)
    call bml_deallocate(x1_bml)
    call bml_deallocate(x2_bml)
    call bml_deallocate(y_bml)

  end subroutine prg_sp2_fermirec

  end module prg_sp2_fermirec_mod
