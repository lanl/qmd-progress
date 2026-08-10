!> A module to perform XLBO integration.
!! @ingroup PROGRESS
!!
!! \brief This module will be used to compute integrate the dynamical variable "n" in xlbo.
!!
module prg_xlbo_mod

  use prg_openfiles_mod
  use bml
  use prg_kernelparser_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  !> Coefficients for K=5 (6-point history) XLBO dissipation - uniform timestep
  real(dp), parameter :: C0 = -6.0_dp
  real(dp), parameter :: C1 = 14.0_dp
  real(dp), parameter :: C2 = -8.0_dp
  real(dp), parameter :: C3 = -3.0_dp
  real(dp), parameter :: C4 = 4.0_dp
  real(dp), parameter :: C5 = -1.0_dp
  real(dp), parameter :: kappa = 1.82_dp
  real(dp), parameter :: alpha = 0.018_dp

  !> K=5 Variable Timestep Coefficient Lookup Tables
  !! Fixed normalization: c_4 = 4.0, c_5 = -1.0 for all patterns
  !! Indexed by 5-bit pattern: bit 0 = most recent dt, bit 4 = oldest dt
  !! Bit value: 0 = half step (dt/2), 1 = full step (dt)
  real(dp), parameter :: XLBO_K5_C0(0:31) = [ &
    -6.00000000_dp, -0.85714286_dp, -1.92857143_dp, -0.32539683_dp, &
    -1.33333333_dp,  0.08571429_dp,  0.87500000_dp,  0.40625000_dp, &
    22.00000000_dp,  4.92857143_dp, 16.00000000_dp,  4.23809524_dp, &
    44.33333333_dp,  9.57142857_dp, 28.37500000_dp,  7.47656250_dp, &
   -62.00000000_dp,-10.50000000_dp,-29.42857143_dp, -6.58730159_dp, &
   -63.00000000_dp,-11.34285714_dp,-30.75000000_dp, -7.11718750_dp, &
   -98.00000000_dp,-14.71428571_dp,-39.00000000_dp, -7.83333333_dp, &
   -75.66666667_dp,-11.80000000_dp,-30.00000000_dp, -6.00000000_dp]

  real(dp), parameter :: XLBO_K5_C1(0:31) = [ &
    14.00000000_dp,  3.00000000_dp,  3.00000000_dp,  0.52380952_dp, &
     1.33333333_dp, -1.94285714_dp, -2.12500000_dp, -1.64062500_dp, &
   -64.00000000_dp,-31.00000000_dp,-31.00000000_dp,-14.28571429_dp, &
  -110.33333333_dp,-46.28571429_dp,-50.62500000_dp,-21.72265625_dp, &
   160.00000000_dp, 53.00000000_dp, 53.00000000_dp, 19.23809524_dp, &
   147.00000000_dp, 47.77142857_dp, 52.25000000_dp, 18.63671875_dp, &
   245.00000000_dp, 68.00000000_dp, 68.00000000_dp, 21.00000000_dp, &
   170.66666667_dp, 44.80000000_dp, 49.00000000_dp, 14.00000000_dp]

  real(dp), parameter :: XLBO_K5_C2(0:31) = [ &
    -8.00000000_dp, -0.57142857_dp,  0.71428571_dp,  2.05555556_dp, &
     1.66666667_dp,  3.70000000_dp,  3.06250000_dp,  3.06250000_dp, &
    67.00000000_dp, 47.28571429_dp, 34.00000000_dp, 26.66666667_dp, &
    79.33333333_dp, 48.00000000_dp, 32.81250000_dp, 23.51562500_dp, &
  -133.00000000_dp,-63.00000000_dp,-40.28571429_dp,-23.77777778_dp, &
   -94.00000000_dp,-42.80000000_dp,-27.12500000_dp,-15.42187500_dp, &
  -192.00000000_dp,-73.14285714_dp,-44.00000000_dp,-19.83333333_dp, &
  -102.66666667_dp,-35.70000000_dp,-21.00000000_dp, -8.00000000_dp]

  real(dp), parameter :: XLBO_K5_C3(0:31) = [ &
    -3.00000000_dp, -4.57142857_dp, -4.78571429_dp, -5.25396825_dp, &
    -4.66666667_dp, -4.84285714_dp, -4.81250000_dp, -4.82812500_dp, &
   -28.00000000_dp,-24.21428571_dp,-22.00000000_dp,-19.61904762_dp, &
   -16.33333333_dp,-14.28571429_dp,-13.56250000_dp,-12.26953125_dp, &
    32.00000000_dp, 17.50000000_dp, 13.71428571_dp,  8.12698413_dp, &
     7.00000000_dp,  3.37142857_dp,  2.62500000_dp,  0.90234375_dp, &
    42.00000000_dp, 16.85714286_dp, 12.00000000_dp,  3.66666667_dp, &
     4.66666667_dp, -0.30000000_dp, -1.00000000_dp, -3.00000000_dp]

  real(dp), parameter :: XLBO_K5_C4(0:31) = [ &
     4.0_dp,  4.0_dp,  4.0_dp,  4.0_dp,  4.0_dp,  4.0_dp,  4.0_dp,  4.0_dp, &
     4.0_dp,  4.0_dp,  4.0_dp,  4.0_dp,  4.0_dp,  4.0_dp,  4.0_dp,  4.0_dp, &
     4.0_dp,  4.0_dp,  4.0_dp,  4.0_dp,  4.0_dp,  4.0_dp,  4.0_dp,  4.0_dp, &
     4.0_dp,  4.0_dp,  4.0_dp,  4.0_dp,  4.0_dp,  4.0_dp,  4.0_dp,  4.0_dp]

  real(dp), parameter :: XLBO_K5_C5(0:31) = [ &
    -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, &
    -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, &
    -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, &
    -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp]

  real(dp), parameter :: XLBO_K5_dK(0:31) = [ &
     0.75000000_dp,  0.28571429_dp,  0.39285714_dp,  0.17063492_dp, &
     0.33333333_dp,  0.06785714_dp,  0.01562500_dp,  0.07812500_dp, &
     2.00000000_dp,  1.14285714_dp,  2.25000000_dp,  1.38095238_dp, &
     5.08333333_dp,  2.71428571_dp,  4.70312500_dp,  2.83203125_dp, &
     7.00000000_dp,  3.00000000_dp,  4.89285714_dp,  2.53968254_dp, &
     8.25000000_dp,  3.72857143_dp,  5.78125000_dp,  3.08984375_dp, &
    11.75000000_dp,  4.57142857_dp,  7.00000000_dp,  3.33333333_dp, &
    10.66666667_dp,  4.32500000_dp,  6.25000000_dp,  3.00000000_dp]

  !> Pattern-specific alpha values for K=5 variable timesteps.
  !! Computed as alpha = 0.018 * 3.0 / d_K, capped at alpha_max/2, so every
  !! pattern stays within the stability bound while maximizing dissipation.
  !! Indexed by 5-bit pattern: bit 0 = most recent dt, bit 4 = oldest dt.
  real(dp), parameter :: XLBO_K5_alpha(0:31) = [ &
        0.072000_dp,     0.186207_dp,     0.138461_dp,     0.129496_dp, &
        0.163636_dp,     0.126162_dp,     0.152756_dp,     0.121484_dp, &
        0.008830_dp,     0.012298_dp,     0.024000_dp,     0.015600_dp, &
        0.004534_dp,     0.005701_dp,     0.011489_dp,     0.019081_dp, &
        0.006900_dp,     0.018000_dp,     0.011040_dp,     0.021261_dp, &
        0.006546_dp,     0.014477_dp,     0.009343_dp,     0.017476_dp, &
        0.004596_dp,     0.011817_dp,     0.007714_dp,     0.016216_dp, &
        0.005061_dp,     0.012420_dp,     0.008640_dp,     0.018000_dp]

  !> Coefficients for K=10 (11-point history) XLBO dissipation
  !> From Niklasson et al. JCP 2009 Table I extended
  real(dp), parameter :: C0_K10 = -858.0_dp
  real(dp), parameter :: C1_K10 = 2652.0_dp
  real(dp), parameter :: C2_K10 = -3094.0_dp
  real(dp), parameter :: C3_K10 = 1496.0_dp
  real(dp), parameter :: C4_K10 = 272.0_dp
  real(dp), parameter :: C5_K10 = -952.0_dp
  real(dp), parameter :: C6_K10 = 731.0_dp
  real(dp), parameter :: C7_K10 = -322.0_dp
  real(dp), parameter :: C8_K10 = 88.0_dp
  real(dp), parameter :: C9_K10 = -14.0_dp
  real(dp), parameter :: C10_K10 = 1.0_dp
  real(dp), parameter :: kappa_K10 = 1.88_dp
  real(dp), parameter :: alpha_K10 = 0.036e-3_dp

  real(dp), parameter :: cc = 0.9_dp;        ! Scaled prg_delta kernel

  !> General xlbo solver type
  !!
  type, public :: xlbo_type

    character(20) :: jobname

    integer :: verbose

    !> Max SCF iterations at every XLBO MD step.
    integer :: maxscfiter

    !> Max SCF iterations for the first minit steps.
    integer :: maxscfinititer

    real(dp) :: threshold

    !> Use SCF the first M_prg_init MD steps
    integer :: minit

    !> Scaled prg_delta Kernel
    real(dp) :: cc

    !> Timestep history for adaptive time step
    !> Size 10 to support K=10 (11-point history); only first 5 used for K=5
    real(dp) :: dt_history(10)
    integer :: nsteps_taken

  end type xlbo_type

  public :: prg_parse_xlbo, prg_xlbo_nint, prg_xlbo_nint_kernel, prg_xlbo_fcoulupdate
  public :: prg_xlbo_nint_kernelTimesRes

contains

  !> The parser for XLBO parser.
  !!
  subroutine prg_parse_xlbo(xlbo,filename)

    implicit none
    type(xlbo_type), intent(inout) :: xlbo
    integer, parameter :: nkey_char = 1, nkey_int = 4, nkey_re = 2, nkey_log = 1
    character(len=*) :: filename

    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=100) :: &
         'JobName=']
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
         'MyJob']

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
         'Verbose=','Mprg_init=','MaxSCFIter=', 'MaxSCFInitIter=']
    integer :: valvector_int(nkey_int) = (/ &
         0, 6, 0, 4 /)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
         'NumThresh=', 'ScaledDeltaKernel=' ]
    real(dp) :: valvector_re(nkey_re) = (/&
         0.0, 0.99 /)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=100) :: &
         'Log1=']
    logical :: valvector_log(nkey_log) = (/&
         .false. /)

    !Start and stop characters
    character(len=50), parameter :: startstop(2) = [character(len=50) :: &
         'XLBO{', '}']

    call prg_parsing_kernel(keyvector_char,valvector_char&
         ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
         keyvector_log,valvector_log,trim(filename),startstop)

    !Characters
    xlbo%JobName = valvector_char(1)

    !Reals
    xlbo%threshold = valvector_re(1)
    xlbo%cc = valvector_re(2)

    !Integers
    xlbo%verbose = valvector_int(1)
    xlbo%minit = valvector_int(2)
    xlbo%maxscfiter = valvector_int(3)
    xlbo%maxscfinititer = valvector_int(4)

    !Initialize timestep history
    xlbo%dt_history = 0.0_dp
    xlbo%nsteps_taken = 0

  end subroutine prg_parse_xlbo



  !> Compute K=5 variable timestep lookup index from dt_history
  !! \brief Converts dt_history into a 5-bit integer for coefficient lookup
  !! \param dt_history Timestep history (most recent first, 5 elements)
  !! \return Bit pattern: 0 = half step (dt/2), 1 = full step (dt)
  function get_K5_history_index(dt_history) result(index)
    implicit none
    real(dp), intent(in) :: dt_history(5)
    integer :: index
    integer :: k

    index = 0
    do k = 1, 5
      ! If timestep is close to 1.0 (full step), set bit k-1
      if (abs(dt_history(k) - 1.0_dp) < 0.1_dp) then
        index = ibset(index, k-1)
      endif
    end do
  end function get_K5_history_index

  !> Compute the per-step XLBO coefficients shared by all integration routines.
  !! \brief Reproduces the uniform, warmup, and adaptive-timestep paths in one place.
  !!  - Without dt (uniform timestep): standard Verlet coefficients and the fixed
  !!    K=5 constants C0..C5 with the constant alpha.
  !!  - With dt but fewer than 6 steps of history (warmup): variable-timestep Verlet
  !!    coefficients, fixed K=5 constants, and the pattern-specific alpha.
  !!  - With dt and full history (adaptive): variable-timestep Verlet coefficients,
  !!    the pattern-specific K=5 coefficient table, and the pattern-specific alpha.
  !! The returned Cn_use / alpha_use / P_* / ka_scale let each caller integrate with a
  !! single expression regardless of path.
  subroutine xlbo_compute_coeffs(xl,P_n,P_n1,ka_scale, &
       C0u,C1u,C2u,C3u,C4u,C5u,alpha_use,dt)
    implicit none
    type(xlbo_type), intent(in) :: xl
    real(dp), intent(out) :: P_n, P_n1, ka_scale
    real(dp), intent(out) :: C0u, C1u, C2u, C3u, C4u, C5u, alpha_use
    real(dp), intent(in), optional :: dt
    logical :: allow_adaptive_timestep
    integer :: hist_idx
    real(dp) :: dt_n, dt_prev, r

    allow_adaptive_timestep = present(dt) .and. xl%nsteps_taken >= 6

    ! Dissipation coefficients default to the fixed K=5 constants.
    C0u = C0; C1u = C1; C2u = C2; C3u = C3; C4u = C4; C5u = C5
    alpha_use = alpha

    ! Variable/uniform Verlet coefficients and the kappa/alpha scaling factor.
    if (present(dt)) then
      dt_n = dt                    ! Current timestep (ratio: 1.0 full, 0.5 half)
      dt_prev = xl%dt_history(1)   ! Previous timestep
      if (dt_prev > 1.0e-12_dp) then
        r = dt_n / dt_prev
        P_n = 1.0_dp + r
        P_n1 = r
        ka_scale = 0.5_dp * dt_n * (dt_n + dt_prev)
      else
        ! First step or no history yet - uniform coefficients, dt^2 scaling.
        P_n = 2.0_dp
        P_n1 = 1.0_dp
        ka_scale = dt_n * dt_n
      endif
      ! Pattern-specific alpha applies during warmup and adaptive phases.
      alpha_use = XLBO_K5_alpha(get_K5_history_index(xl%dt_history(1:5)))
    else
      ! Uniform timestep: standard Verlet coefficients, no extra scaling.
      P_n = 2.0_dp
      P_n1 = 1.0_dp
      ka_scale = 1.0_dp
    endif

    ! Full history available: use the pattern-specific coefficient table.
    if (allow_adaptive_timestep) then
      hist_idx = get_K5_history_index(xl%dt_history(1:5))
      C0u = XLBO_K5_C0(hist_idx)
      C1u = XLBO_K5_C1(hist_idx)
      C2u = XLBO_K5_C2(hist_idx)
      C3u = XLBO_K5_C3(hist_idx)
      C4u = XLBO_K5_C4(hist_idx)
      C5u = XLBO_K5_C5(hist_idx)
      alpha_use = XLBO_K5_alpha(hist_idx)
    endif

  end subroutine xlbo_compute_coeffs

  !> Push the current timestep ratio onto the history stack.
  !! \brief dt is a ratio (1.0 = full user timestep, 0.5 = half). Only called when
  !! adaptive timestepping is active, so the caller guards on present(dt).
  subroutine xlbo_push_dt_history(xl,dt)
    implicit none
    type(xlbo_type), intent(inout) :: xl
    real(dp), intent(in) :: dt

    xl%dt_history(5) = xl%dt_history(4)
    xl%dt_history(4) = xl%dt_history(3)
    xl%dt_history(3) = xl%dt_history(2)
    xl%dt_history(2) = xl%dt_history(1)
    xl%dt_history(1) = dt
    xl%nsteps_taken = xl%nsteps_taken + 1

  end subroutine xlbo_push_dt_history

  !> This routine integrates the dynamical variable "n"
  !! \param charges
  subroutine prg_xlbo_nint(charges,n,n_0,n_1,n_2,n_3,n_4,n_5,mdstep,xl,dt)
    implicit none
    real(dp), allocatable, intent(inout) :: n(:), n_0(:), n_1(:), n_2(:), n_3(:), n_4(:), n_5(:)
    real(dp), allocatable, intent(in) :: charges(:)
    type(xlbo_type), intent(inout) :: xl
    integer, intent(in) :: mdstep
    real(dp), intent(in), optional :: dt
    integer :: nats
    real(dp) :: alpha_use
    real(dp) :: C0u, C1u, C2u, C3u, C4u, C5u
    real(dp) :: P_n, P_n1, ka_scale

    nats = size(charges,dim=1)

    if(.not.allocated(n))then
      allocate(n(nats))
      allocate(n_0(nats))
      allocate(n_1(nats))
      allocate(n_2(nats))
      allocate(n_3(nats))
      allocate(n_4(nats))
      allocate(n_5(nats))
    endif

    if(mdstep.le.1)then
      n = charges;
      n_0 = charges;
      n_1 = charges;
      n_2 = charges;
      n_3 = charges;
      n_4 = charges;
      n_5 = charges;
      xl%dt_history = 0.0_dp
      xl%nsteps_taken = 0
    endif

    call xlbo_compute_coeffs(xl,P_n,P_n1,ka_scale,C0u,C1u,C2u,C3u,C4u,C5u,alpha_use,dt)

    n = P_n*n_0 - P_n1*n_1 + xl%cc*ka_scale*kappa*(charges-n) &
         + alpha_use*(C0u*n_0+C1u*n_1+C2u*n_2+C3u*n_3+C4u*n_4+C5u*n_5)

    ! Shift history arrays
    n_5 = n_4; n_4 = n_3; n_3 = n_2; n_2 = n_1; n_1 = n_0; n_0 = n

    if (present(dt)) call xlbo_push_dt_history(xl,dt)

  end subroutine prg_xlbo_nint


  !> This routine integrates the dynamical variable "n"
  !! \param charges
  subroutine prg_xlbo_nint_kernel(charges,n,n_0,n_1,n_2,n_3,n_4,n_5,mdstep,kernel,xl,dt)
    implicit none
    real(dp), allocatable, intent(inout) :: n(:), n_0(:), n_1(:), n_2(:), n_3(:), n_4(:), n_5(:)
    real(dp), allocatable, intent(in) :: charges(:)
    real(dp), allocatable, intent(in) :: kernel(:,:)
    type(xlbo_type), intent(inout) :: xl
    integer, intent(in) :: mdstep
    real(dp), intent(in), optional :: dt
    integer :: nats
    real(dp) :: alpha_use
    real(dp) :: C0u, C1u, C2u, C3u, C4u, C5u
    real(dp) :: P_n, P_n1, ka_scale

    nats = size(charges,dim=1)

    if(.not.allocated(n))then
      allocate(n(nats))
      allocate(n_0(nats))
      allocate(n_1(nats))
      allocate(n_2(nats))
      allocate(n_3(nats))
      allocate(n_4(nats))
      allocate(n_5(nats))
    endif

    if(mdstep.le.1)then
      n = charges;
      n_0 = charges;
      n_1 = charges;
      n_2 = charges;
      n_3 = charges;
      n_4 = charges;
      n_5 = charges;
      xl%dt_history = 0.0_dp
      xl%nsteps_taken = 0
    endif

    ! From developper's code
    !   dn2dt2 = -MATMUL(KK0,(q-n))
    !   n = 2*n_0 - n_1 + kappa*dn2dt2 +
    !   alpha*(C0*n_0+C1*n_1+C2*n_2+C3*n_3+C4*n_4+C5*n_5+C6*n_6)
    !   n_6 = n_5; n_5 = n_4; n_4 = n_3; n_3 = n_2; n_2 = n_1; n_1 = n_0; n_0 = n

    call xlbo_compute_coeffs(xl,P_n,P_n1,ka_scale,C0u,C1u,C2u,C3u,C4u,C5u,alpha_use,dt)

    n = P_n*n_0 - P_n1*n_1 - ka_scale*kappa*matmul(kernel,(charges-n)) &
         + alpha_use*(C0u*n_0+C1u*n_1+C2u*n_2+C3u*n_3+C4u*n_4+C5u*n_5)

    ! Shift history arrays
    n_5 = n_4; n_4 = n_3; n_3 = n_2; n_2 = n_1; n_1 = n_0; n_0 = n

    if (present(dt)) call xlbo_push_dt_history(xl,dt)

  end subroutine prg_xlbo_nint_kernel


  !> This routine integrates the dynamical variable "n"
  !! \brief In this case we are passing a premultiplied ressidue x kernel
  !! tis is done to avoid rank-specific multiplication within this routine.
  !! \param charges
  subroutine prg_xlbo_nint_kernelTimesRes(charges,n,n_0,n_1,n_2,n_3,n_4,n_5,mdstep,kernelTimesRes,xl,dt)
    implicit none
    real(dp), allocatable, intent(inout) :: n(:), n_0(:), n_1(:), n_2(:), n_3(:), n_4(:), n_5(:)
    real(dp), allocatable, intent(in) :: charges(:)
    real(dp), allocatable, intent(in) :: kernelTimesRes(:)
    type(xlbo_type), intent(inout) :: xl
    integer, intent(in) :: mdstep
    real(dp), intent(in), optional :: dt
    integer :: nats
    real(dp) :: alpha_use
    real(dp) :: C0u, C1u, C2u, C3u, C4u, C5u
    real(dp) :: P_n, P_n1, ka_scale

    nats = size(charges,dim=1)

    if(.not.allocated(n))then
      allocate(n(nats))
      allocate(n_0(nats))
      allocate(n_1(nats))
      allocate(n_2(nats))
      allocate(n_3(nats))
      allocate(n_4(nats))
      allocate(n_5(nats))
    endif

    if(mdstep.le.1)then
      n = charges;
      n_0 = charges;
      n_1 = charges;
      n_2 = charges;
      n_3 = charges;
      n_4 = charges;
      n_5 = charges;
      xl%dt_history = 0.0_dp
      xl%nsteps_taken = 0
    endif

    call xlbo_compute_coeffs(xl,P_n,P_n1,ka_scale,C0u,C1u,C2u,C3u,C4u,C5u,alpha_use,dt)

    n = P_n*n_0 - P_n1*n_1 - ka_scale*kappa*kernelTimesRes &
         & + alpha_use*(C0u*n_0+C1u*n_1+C2u*n_2+C3u*n_3+C4u*n_4+C5u*n_5)

    ! Shift history arrays
    n_5 = n_4; n_4 = n_3; n_3 = n_2; n_2 = n_1; n_1 = n_0; n_0 = n

    if (present(dt)) call xlbo_push_dt_history(xl,dt)

  end subroutine prg_xlbo_nint_kernelTimesRes



  !> Adjust forces for the linearized XLBOMD functional
  !! \param charges
  subroutine prg_xlbo_fcoulupdate(fcoul,charges,n)
    implicit none
    real(dp), intent(inout) :: fcoul(:,:),charges(:)
    real(dp), intent(inout) :: n(:)

    fcoul(1,:) = (2.0_dp*charges(:)-n(:))*fcoul(1,:)/n(:);
    fcoul(2,:) = (2.0_dp*charges(:)-n(:))*fcoul(2,:)/n(:);
    fcoul(3,:) = (2.0_dp*charges(:)-n(:))*fcoul(3,:)/n(:);

  end subroutine prg_xlbo_fcoulupdate

end module prg_xlbo_mod
