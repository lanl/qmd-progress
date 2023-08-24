!> The SP2 Fermi module.
!! \ingroup PROGRESS
!!
!! \brief This subroutine implements Niklasson's truncated SP2 density matrix
!!  purification algorithm.
!!
module prg_sp2_fermi_mod

  use bml
  use prg_normalize_mod
  use prg_timer_mod
  use prg_parallel_mod

  implicit none

  private  !Everything is private by default

  integer, parameter :: dp = kind(1.0d0)

  public :: prg_sp2_fermi_init
  public :: prg_sp2_fermi_init_norecs
  public :: prg_sp2_fermi
  public :: prg_sp2_entropy_function
  public :: sp2_entropy_ts
  public :: sp2_inverse

contains

  !> Truncated SP2 prg_initialization.
  !! \param h_bml Input Hamiltonian matrix.
  !! \param nsteps Number of sp2 iterations.
  !! \param nocc Number of occupied states.
  !! \param tscale Temperature rescaling factor.
  !! \param threshold Threshold for multiplication.
  !! \param occErrLimit Occupation error limit.
  !! \param traceLimit Trace limit.
  !! \param x_bml Output prg_initial matrix.
  !! \param mu Shifted chemical potential
  !! \param beta Output inverse temperature.
  !! \param h1 Output temperature-scaled minimum gershgorin bound.
  !! \param hN Output temperature-scaled maximum gershgorin bound.
  !! \param sgnlist SP2 sequence
  subroutine prg_sp2_fermi_init(h_bml, nsteps, nocc, tscale, threshold, &
       occErrLimit, traceLimit, x_bml, mu, beta, h1, hN, sgnlist)

    implicit none

    type(bml_matrix_t), intent(in) :: h_bml
    type(bml_matrix_t), intent(inout) :: x_bml
    integer, intent(in) :: nsteps
    integer, intent(inout) :: sgnlist(:)
    real(dp), intent(in) :: nocc, tscale, threshold
    real(dp), intent(in) :: occErrLimit, traceLimit
    real(dp), intent(inout) :: mu, beta, h1, hN

    type(bml_matrix_t) :: x1_bml, x2_bml, tmp_bml, i_bml
    real(dp) :: lambda, occErr, sfactor
    real(dp) :: traceX0, traceX1, traceX2, traceX
    real(dp), allocatable :: gbnd(:), trace(:)
    logical :: firstTime
    character(20) :: bml_type
    integer :: N, M, i

    bml_type = bml_get_type(h_bml)
    N = bml_get_N(h_bml)
    M = bml_get_M(h_bml)

    allocate(gbnd(2))

    !> Calculate Gershgorin bounds and rescale
    call bml_gershgorin(h_bml, gbnd)
    mu = 0.5 * (gbnd(2) + gbnd(1))
    h1 = tscale * gbnd(1)
    hN = tscale * gbnd(2)

    deallocate(gbnd)

    lambda = 0.0_dp
    occErr = 1.0_dp
    firstTime = .true.

    allocate(trace(2))

    call bml_identity_matrix(bml_type, bml_element_real, dp, N, M, i_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, x1_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, x2_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, tmp_bml)

    do while (occErr .gt. occErrLimit)

      call bml_copy(h_bml, x_bml)
      call prg_normalize_fermi(x_bml, h1, hN, mu)

      ! X1 = -I/(hN-h1)
      call bml_copy(i_bml, x1_bml)
      sfactor = -1.0_dp / (hN - h1)
      call bml_scale(sfactor, x1_bml)

      do i = 1, nsteps
        call bml_multiply_x2(x_bml, x2_bml, threshold, trace)
        traceX0 = trace(1)
        traceX2 = trace(2)

        !> Determine sequence branching first time through
        if (firstTime .eqv. .true.) then
          if (abs(traceX2-nocc) .lt. abs(2.0_dp*traceX0-traceX2-nocc)) then
            sgnlist(i) = -1
          else
            sgnlist(i) = 1
          end if
        end if

        ! X1 = X1 + sgnlist(i)*(X1 - X0*X1 - X1*X0)
        ! if sgnlist == 1, X1 = 2 * X1 - (X0*X1 + X1*X0)
        ! if sgnlist == -1, X1 = X0*X1 + X1*X0
        !
        ! tmp = X0*X1 + X1*X0
        call bml_multiply(x_bml, x1_bml, tmp_bml, 1.0_dp, 0.0_dp, threshold)
        call bml_multiply(x1_bml, x_bml, tmp_bml, 1.0_dp, 1.0_dp, threshold)

        if (sgnlist(i) .eq. 1) then
          ! X1 = 2 * X1 - tmp
          call bml_add_deprecated(2.0_dp, x1_bml, -1.0_dp, tmp_bml, threshold)
        else
          call bml_copy(tmp_bml, x1_bml)
        endif

        ! X0 = X0 + sgnlist(i)*(X0 - X0_2)
        ! if sgnlist == 1, X0 = 2.0*X0 - X0_2
        ! if sgnlist == -1, X0 = X0_2
        !
        if (sgnlist(i) .eq. 1) then
          ! X0 = 2 * X0 - X2
          call bml_add_deprecated(2.0_dp, x_bml, -1.0_dp, x2_bml, threshold)
        else
          call bml_copy(x2_bml, x_bml)
        endif

      end do

      firstTime = .false.
      traceX0 = bml_trace(x_bml)
      traceX1 = bml_trace(x1_bml)
      occErr = abs(nocc - traceX0)
      ! Newton=Rhapson step to correct for occupation
      if (abs(traceX1) .gt. traceLimit) then
        lambda = (nocc - traceX0) / traceX1
      else
        lambda = 0.0_dp
      end if
      mu = mu + lambda
    end do

    deallocate(trace)
    call bml_deallocate(x2_bml)

    ! X0*(I-X0)
    ! I = I - X0
    call bml_add_deprecated(1.0_dp, i_bml, -1.0_dp, x_bml, threshold)
    ! tmp = X0*I
    call bml_multiply(x_bml, i_bml, tmp_bml, 1.0_dp, 0.0_dp, threshold)
    traceX = bml_trace(tmp_bml)
    traceX1 = bml_trace(x1_bml)
    if (abs(traceX) .gt. traceLimit) then
      beta = -traceX1 / traceX
    else
      beta = -1000.0_dp
    end if

    ! X = 2 * X
    !call bml_scale(2.0_dp, x_bml)

    call bml_deallocate(tmp_bml)
    call bml_deallocate(i_bml)
    call bml_deallocate(x1_bml)

  end subroutine prg_sp2_fermi_init


  !> Truncated SP2 prg_initialization.
  !! This routine also gives back the Number of SP2 recursive steps
  !! that gets a Pseudo-Fermi distribution with a temperature close to
  !! the target temperature which is entered using parameter beta  = (1/KbT).
  !! \param h_bml Input Hamiltonian matrix.
  !! \param nsteps Output number of sp2 iterations.
  !! \param nocc Number of occupied states.
  !! \param tscale Temperature rescaling factor.
  !! \param threshold Threshold for multiplication.
  !! \param occErrLimit Occupation error limit.
  !! \param traceLimit Trace limit.
  !! \param x_bml Output prg_initial matrix.
  !! \param mu Shifted chemical potential
  !! \param beta Input guess and output inverse temperature.
  !! \param h1 Output temperature-scaled minimum gershgorin bound.
  !! \param hN Output temperature-scaled maximum gershgorin bound.
  !! \param sgnlist SP2 sequence
  !! \param verbose Optional parameter for verbosity.
  subroutine prg_sp2_fermi_init_norecs(h_bml, nsteps, nocc, tscale, threshold, &
       occErrLimit, traceLimit, x_bml, mu, beta, h1, hN, sgnlist, verbose)

    implicit none

    type(bml_matrix_t), intent(in) :: h_bml
    type(bml_matrix_t), intent(inout) :: x_bml
    integer, intent(inout) :: nsteps
    integer, intent(inout) :: sgnlist(:)
    real(dp), intent(in) :: nocc, tscale, threshold
    real(dp), intent(in) :: occErrLimit, traceLimit
    real(dp), intent(inout) :: mu, beta, h1, hN

    type(bml_matrix_t) :: x1_bml, x2_bml, tmp_bml, i_bml
    real(dp) :: lambda, occErr, sfactor, maxder
    real(dp) :: traceX0, traceX1, traceX2, traceX, beta0
    real(dp), allocatable :: gbnd(:), trace(:), probe(:), probe_2(:)
    logical :: firstTime
    character(20) :: bml_type
    integer :: N, M, i
    integer, optional :: verbose

    if(present(verbose))then
      if(verbose >= 1) write(*,*)"Inside prg_sp2_fermi_init_norecs ..."
    endif

    bml_type = bml_get_type(h_bml)
    N = bml_get_N(h_bml)
    M = bml_get_M(h_bml)
    beta0 = beta

    allocate(gbnd(2))

    !> Calculate Gershgorin bounds and rescale
    call bml_gershgorin(h_bml, gbnd)
    mu = 0.5 * (gbnd(2) + gbnd(1))
    h1 = tscale * gbnd(1)
    hN = tscale * gbnd(2)

    lambda = 0.0_dp
    occErr = 1.0_dp
    firstTime = .true.

    allocate(trace(2))

    call bml_identity_matrix(bml_type, bml_element_real, dp, N, M, i_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, x1_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, x2_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, tmp_bml)

    allocate(probe(1000),probe_2(1000))


    do while (occErr .gt. occErrLimit)

      call bml_copy(h_bml, x_bml)
      call prg_normalize_fermi(x_bml, h1, hN, mu)

      ! Probe function to compute the derivative at mu.
      do i = 1,1000
        probe(i) = 1.0_dp - i*0.001_dp
      enddo

      ! X1 = -I/(hN-h1)
      call bml_copy(i_bml, x1_bml)
      sfactor = -1.0_dp / (hN - h1)
      call bml_scale(sfactor, x1_bml)

      do i = 1, nsteps
        call bml_multiply_x2(x_bml, x2_bml, threshold, trace)

        probe_2(:) = probe(:)*probe(:)

        traceX0 = trace(1)
        traceX2 = trace(2)

        !> Determine sequence branching first time through
        if (firstTime .eqv. .true.) then
          if (abs(traceX2-nocc) .lt. abs(2.0_dp*traceX0-traceX2-nocc)) then
            sgnlist(i) = -1
          else
            sgnlist(i) = 1
          end if
        end if

        ! X1 = X1 + sgnlist(i)*(X1 - X0*X1 - X1*X0)
        ! if sgnlist == 1, X1 = 2 * X1 - (X0*X1 + X1*X0)
        ! if sgnlist == -1, X1 = X0*X1 + X1*X0
        !
        ! tmp = X0*X1 + X1*X0
        call bml_multiply(x_bml, x1_bml, tmp_bml, 1.0_dp, 0.0_dp, threshold)
        call bml_multiply(x1_bml, x_bml, tmp_bml, 1.0_dp, 1.0_dp, threshold)

        if (sgnlist(i) .eq. 1) then
          ! X1 = 2 * X1 - tmp
          call bml_add_deprecated(2.0_dp, x1_bml, -1.0_dp, tmp_bml, threshold)
        else
          call bml_copy(tmp_bml, x1_bml)
        endif

        ! X0 = X0 + sgnlist(i)*(X0 - X0_2)
        ! if sgnlist == 1, X0 = 2.0*X0 - X0_2
        ! if sgnlist == -1, X0 = X0_2
        !
        if (sgnlist(i) .eq. 1) then
          ! X0 = 2 * X0 - X2
          call bml_add_deprecated(2.0_dp, x_bml, -1.0_dp, x2_bml, threshold)
          probe = 2.0_dp*probe - probe_2
        else
          call bml_copy(x2_bml, x_bml)
          probe = probe_2
        endif

        maxder = absmaxderivative(probe,0.001_dp)
        beta = (4.0_dp*maxder)/(gbnd(2)-gbnd(1))
        if(beta > beta0) then
          nsteps = i
          exit
        endif

      end do

      ! Write probe function into a file (only for debugging purposes)
      ! do i = 1,1000
      !   write(1000,*)gbnd(1) + ((gbnd(2)-gbnd(1))/1000.0_dp)*i,probe(i)
      ! enddo

      if(present(verbose))then
        if(verbose >= 1) then
          write(*,*)"beta = ",beta
          write(*,*)"kbT = ",1.0_dp/beta
        endif
      endif

      firstTime = .false.
      traceX0 = bml_trace(x_bml)
      traceX1 = bml_trace(x1_bml)
      occErr = abs(nocc - traceX0)
      ! Newton=Rhapson step to correct for occupation
      if (abs(traceX1) .gt. traceLimit) then
        lambda = (nocc - traceX0) / traceX1
      else
        lambda = 0.0_dp
      end if
      mu = mu + lambda
    end do

    deallocate(trace)
    call bml_deallocate(x2_bml)

    ! X0*(I-X0)
    ! I = I - X0
    call bml_add_deprecated(1.0_dp, i_bml, -1.0_dp, x_bml, threshold)
    ! tmp = X0*I
    call bml_print_matrix("x_bml",x_bml,1,4,1,4)
    call bml_print_matrix("i_bml",i_bml,1,4,1,4)
    call bml_multiply(x_bml, i_bml, tmp_bml, 1.0_dp, 0.0_dp, threshold)

    call bml_print_matrix("tmp_bml",tmp_bml,1,4,1,4)

    traceX = bml_trace(tmp_bml)
    traceX1 = bml_trace(x1_bml)
    if (abs(traceX) .gt. traceLimit) then
      beta = -traceX1 / traceX
    else
      beta = -1000.0_dp
    end if

    ! X = 2 * X
    !call bml_scale(2.0_dp, x_bml)
    deallocate(gbnd)
    call bml_deallocate(tmp_bml)
    call bml_deallocate(i_bml)
    call bml_deallocate(x1_bml)

  end subroutine prg_sp2_fermi_init_norecs

  !> Calculate Truncated SP2.
  !! \param h_bml Hamiltonian matrix
  !! \param osteps Outer loop steps
  !! \param nsteps Number of sequence branches
  !! \param nocc Number of occupation states
  !! \param mu Shifted chemical potential
  !! \param beta Inverse temperature
  !! \param h1 Minimum scaled Gershgorin bound.
  !! \param hN Maximum scaled Gershgorin bound.
  !! \param sgnlist SP2 sequence
  !! \param threshold Threshold for multiplies
  !! \param eps Occupation error limit
  !! \param traceLimit Trace limit
  !! \param x_bml Output density matrix
  subroutine prg_sp2_fermi(h_bml, osteps, nsteps, nocc, mu, beta, h1, hN, sgnlist, &
       threshold, eps, traceLimit, x_bml)

    implicit none

    type(bml_matrix_t), intent(in) :: h_bml
    type(bml_matrix_t), intent(inout) :: x_bml
    integer, intent(in) :: osteps, nsteps
    integer, intent(in) :: sgnlist(:)
    real(dp), intent(in) :: nocc, threshold, eps, traceLimit
    real(dp), intent(inout) :: beta, h1, hN
    real(dp), intent(inout) :: mu

    type(bml_matrix_t) :: x2_bml, dx_bml, i_bml
    real(dp), allocatable :: trace(:), gbnd(:)
    real(dp) :: sfactor, occErr, traceX0, traceX2, traceDX, lambda
    real(dp) :: a, b
    integer :: iter, i, N, M
    character(20) :: bml_type

    bml_type = bml_get_type(h_bml)
    N = bml_get_N(h_bml)
    M = bml_get_M(h_bml)

    allocate(trace(2))

    call bml_identity_matrix(bml_type, bml_element_real, dp, N, M, i_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, dx_bml)
    call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, x2_bml)

    occErr = 1.0_dp + eps
    iter = 0
    do while ((osteps .eq. 0 .and. occErr .gt. eps) .or. &
         (osteps .gt. 0 .and. iter .lt. osteps))
      iter = iter + 1
      call bml_copy(h_bml, x_bml)
      call prg_normalize_fermi(x_bml, h1, hN, mu)

      do i = 1, nsteps
        call bml_multiply_x2(x_bml, x2_bml, threshold, trace)
        traceX0 = trace(1)
        traceX2 = trace(2)

        ! X0 = X0 + sgnlist(i)*(X0 - X0_2)
        if (sgnlist(i) .eq. 1) then
          call bml_add_deprecated(2.0_dp, x_bml, -1.0_dp, x2_bml, threshold)
        else
          call bml_copy(x2_bml, x_bml)
        endif

      end do

      traceX0 = bml_trace(x_bml)
      occErr = abs(nocc - traceX0)

      ! DX = -beta*X0*(I-X0)
      call bml_copy(i_bml, x2_bml)
      call bml_add_deprecated(1.0_dp, x2_bml, -1.0_dp, x_bml, threshold)
      call bml_multiply(x_bml, x2_bml, dx_bml, -beta, 0.0_dp, threshold)
      traceDX = bml_trace(dx_bml)

      ! Newton-Rhapson step to correct for occupation
      if (abs(traceDX) .gt. traceLimit) then
        lambda = (nocc - traceX0) / traceDX
      else
        lambda = 0.0_dp
      end if
      mu = mu + lambda

    end do

    ! Correction for occupation
    call bml_add_deprecated(1.0_dp, x_bml, lambda, dx_bml, threshold)

    ! X = 2*X
    !call bml_scale(2.0_dp, x_bml)

    deallocate(trace)
    call bml_deallocate(i_bml)
    call bml_deallocate(x2_bml)
    call bml_deallocate(dx_bml)

  end subroutine prg_sp2_fermi

  !> Calculate SP2 entropy function using gaussian quadrature.
  !! Note that GG and ee are allocated and returned
  !! from this routine.
  !! \param mu Shifted chemical potential
  !! \param h1 Minimum scaled Gershgorin bound
  !! \param hN Maximum scaled Gershgorin bound
  !! \param nsteps Number of SP2 steps
  !! \param sgnlist SP2 sequence
  !! \param GG Entropy function
  !! \param ee 1D mesh
  subroutine prg_sp2_entropy_function(mu, h1, hN, nsteps, sgnlist, GG, ee)

    implicit none

    real(dp), intent(in) :: mu, h1, hN
    integer, intent(in) :: nsteps, sgnlist(:)
    real(dp), intent(inout), allocatable :: GG(:), ee(:)

    real(dp) :: dh, c1, c2, x1, x2, iint
    real(dp) :: a, b, c_1, c_2, x_1, x_2
    integer :: i, k, N

    dh = 0.0001_dp
    c1 = 1.0_dp
    c2 = 1.0_dp
    x1 = 0.5773502691896257_dp
    x2 = -x1
    N = 10001 ! number of elements in ee

    if (allocated(GG)) then
      deallocate(GG)
    end if
    if (allocated(ee)) then
      deallocate(ee)
    end if

    allocate(GG(N))
    allocate(ee(N))

    ! Set ee = 0.0 to 1.0 by 0.0001
    do i = 0, N-1
      ee(i+1) = i * dh
    end do

    GG(1) = 0.0_dp
    iint = 0.0_dp
    do i = 1, N-1
      a = ee(i)
      b = ee(i+1)
      c_1 = c1*(b-a)/2.0_dp
      c_2 = c2*(b-a)/2.0_dp;
      x_1 = ((b-a)*x1+(b+a))/2.0_dp
      x_2 = ((b-a)*x2+(b+a))/2.0_dp
      iint = iint + c_1*sp2_inverse(x_1,mu,h1,hN,nsteps,sgnlist)
      iint = iint + c_2*sp2_inverse(x_2,mu,h1,hN,nsteps,sgnlist)
      GG(i+1) = iint
    end do

    GG = GG-GG(N)*ee

  end subroutine prg_sp2_entropy_function

  !> Test SP2 entropy.
  !! Get the entropy contribution TS to the total free energy.
  !! \param D0_bml BML matrix
  !! \param GG Entropy function
  !! \param ee 1D mesh
  !! \param TS Energy contribution
  function sp2_entropy_ts(D0_bml, GG, ee) result(TS)

    implicit none

    type(bml_matrix_t), intent(in) :: D0_bml
    real(dp), intent(in) :: GG(*), ee(*)
    real(dp) :: TS

    type(bml_matrix_t) :: aux_bml
    real(dp), allocatable :: hh(:)
    real(dp) :: hs, threshold
    integer :: s, j, N
    character(20) :: bml_type, bml_dmode

    N = bml_get_N(D0_bml)
    bml_type = bml_get_type(D0_bml)
    bml_dmode = bml_get_distribution_mode(D0_bml)

    allocate(hh(N))

    ! Diagonalize D0
    call bml_zero_matrix(bml_type, bml_element_real, &
         dp, N, N, aux_bml)

    call bml_diagonalize(D0_bml, hh, aux_bml)

    call bml_deallocate(aux_bml)

    TS = 0.0_dp

    do s = 1,N
      hs = abs(hh(s))
      j = floor(hs/0.0001_dp + 0.00000_dp) + 1

      if (j .gt. 0 .and. j .le. 10001) then
        TS = TS + ((hs-ee(j))*GG(j+1) + &
             (ee(j+1)-hs)*GG(j))/(ee(j+1)-ee(j))
      end if
    end do

    deallocate(hh)

  end function sp2_entropy_ts

  !> Calculate the SP2 inverse.
  !! \param f Occupation factor
  !! \param mu Shifted chemical potential
  !! \param h1 Minimum scaled Gershgorin bound
  !! \param hN Maximum scaled Gershgorin bound
  !! \param nsteps Numbers of SP2 iterations
  !! \param sgnlist SP2 sequence
  !! \param ee Energy value
  function sp2_inverse(f, mu, h1, hN, nsteps, sgnlist) result(ee)

    implicit none

    real(dp), intent(in) :: f, mu, h1, hN
    integer, intent(in) :: nsteps, sgnlist(:)

    real(dp) :: sgn, c, ee
    integer :: i

    ee = f
    do i = nsteps, 1, -1
      sgn = float(sgnlist(i))
      c = (1.0_dp+sgn)/2.0_dp
      ee = c - sgn*sqrt(abs(c - sgn*ee))
    end do

    ee = hN-mu-ee*(hN-h1)

  end function sp2_inverse

  !> Gets the absolute maximum of the derivative of a function.
  !! \param func.
  !! \param de Energy step.
  !!
  real(dp) function absmaxderivative(func,de)

    implicit none
    real(dp), intent(in) :: func(:), de
    integer :: j

    absmaxderivative = -10000.0d0

    do j=1,size(func, dim=1)-1
      if(abs(func(j+1) - func(j))/de > absmaxderivative) &
           absmaxderivative = abs(func(j+1) - func(j))/de
    enddo

  end function absmaxderivative


end module prg_sp2_fermi_mod
