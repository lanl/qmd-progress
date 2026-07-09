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

  !> Pattern-specific alpha values for K=5 XLBO dissipation
  !! Conservative scaling: alpha = 0.000796 × 3.0 / d_K
  !! All values ≤ α_max/2, ensuring stability for all 32 patterns
  !! Limiting pattern: 6 (d_K=0.015625, exactly at α_max/2)
  !! Average safety margin: 20× above usage (suitable for large systems)
  !! Indexed by 5-bit pattern: bit 0 = most recent dt, bit 4 = oldest dt
  real(dp), parameter :: XLBO_K5_alpha(0:31) = [ &
        0.003184_dp,     0.008358_dp,     0.006079_dp,     0.013995_dp, &
        0.007164_dp,     0.035192_dp,     0.152832_dp,     0.030566_dp, &
        0.001194_dp,     0.002090_dp,     0.001061_dp,     0.001729_dp, &
        0.000470_dp,     0.000880_dp,     0.000508_dp,     0.000843_dp, &
        0.000341_dp,     0.000796_dp,     0.000488_dp,     0.000940_dp, &
        0.000289_dp,     0.000640_dp,     0.000413_dp,     0.000773_dp, &
        0.000203_dp,     0.000522_dp,     0.000341_dp,     0.000716_dp, &
        0.000224_dp,     0.000552_dp,     0.000382_dp,     0.000796_dp]

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

    !> Timestep history for interpolation-based integration
    !> Size 10 to support K=10 (11-point history); only first 5 used for K=5
    real(dp) :: dt_history(10)
    integer :: nsteps_taken

    !> Use extended history (K=10, 11-point) instead of default (K=5, 6-point)
    logical :: extended_history

    !> Use variable timestep coefficients (alternative to interpolation for K=5)
    logical :: use_variable_coeffs

  end type xlbo_type

  public :: prg_parse_xlbo, prg_xlbo_nint, prg_xlbo_nint_kernel, prg_xlbo_fcoulupdate
  public :: prg_xlbo_nint_kernelTimesRes

contains

  !> The parser for XLBO parser.
  !!
  subroutine prg_parse_xlbo(xlbo,filename)

    implicit none
    type(xlbo_type), intent(inout) :: xlbo
    integer, parameter :: nkey_char = 1, nkey_int = 4, nkey_re = 2, nkey_log = 2
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
         'Log1=', 'ExtendedHistory=']
    logical :: valvector_log(nkey_log) = (/&
         .false., .false. /)

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

    !Logicals
    xlbo%extended_history = valvector_log(2)

    !Integers
    xlbo%verbose = valvector_int(1)
    xlbo%minit = valvector_int(2)
    xlbo%maxscfiter = valvector_int(3)
    xlbo%maxscfinititer = valvector_int(4)

    !Initialize timestep history
    xlbo%dt_history = 0.0_dp
    xlbo%nsteps_taken = 0
    xlbo%use_variable_coeffs = .false.

  end subroutine prg_parse_xlbo


  !> Compute cubic spline second derivatives using tridiagonal solver
  !! \param x Array of x values (must be monotonic)
  !! \param y Array of y values
  !! \param n Number of points
  !! \param y2 Output: second derivatives at each point (natural boundary conditions)
  subroutine cubic_spline_coeffs(x, y, n, y2)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: x(0:n-1), y(0:n-1)
    real(dp), intent(out) :: y2(0:n-1)
    real(dp) :: u(0:n-1), sig, p
    integer :: i

    ! Natural spline: second derivative = 0 at endpoints
    y2(0) = 0.0_dp
    u(0) = 0.0_dp

    ! Forward sweep of tridiagonal solver
    do i = 1, n-2
      sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
      p = sig * y2(i-1) + 2.0_dp
      y2(i) = (sig - 1.0_dp) / p
      u(i) = (y(i+1) - y(i)) / (x(i+1) - x(i)) - (y(i) - y(i-1)) / (x(i) - x(i-1))
      u(i) = (6.0_dp * u(i) / (x(i+1) - x(i-1)) - sig * u(i-1)) / p
    enddo

    ! Natural spline: second derivative = 0 at right endpoint
    y2(n-1) = 0.0_dp

    ! Back substitution
    do i = n-2, 0, -1
      y2(i) = y2(i) * y2(i+1) + u(i)
    enddo

  end subroutine cubic_spline_coeffs


  !> Evaluate cubic spline at a given point
  !! \param xa Array of x values
  !! \param ya Array of y values
  !! \param y2a Array of second derivatives (from cubic_spline_coeffs)
  !! \param n Number of points
  !! \param x Point at which to evaluate
  !! \param y Output: interpolated value
  subroutine cubic_spline_eval(xa, ya, y2a, n, x, y)
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: xa(0:n-1), ya(0:n-1), y2a(0:n-1), x
    real(dp), intent(out) :: y
    integer :: klo, khi, k
    real(dp) :: h, a, b

    ! Binary search for bracketing interval
    ! Handle descending arrays (time values are negative and decreasing)
    klo = 0
    khi = n - 1
    if (xa(0) > xa(n-1)) then
      ! Descending array
      do while (khi - klo > 1)
        k = (khi + klo) / 2
        if (xa(k) < x) then
          khi = k
        else
          klo = k
        endif
      enddo
    else
      ! Ascending array
      do while (khi - klo > 1)
        k = (khi + klo) / 2
        if (xa(k) > x) then
          khi = k
        else
          klo = k
        endif
      enddo
    endif

    ! Evaluate cubic polynomial in this interval
    h = xa(khi) - xa(klo)
    if (abs(h) < 1.0e-12_dp) then
      ! Degenerate case: coincident points
      y = ya(klo)
      return
    endif

    a = (xa(khi) - x) / h
    b = (x - xa(klo)) / h
    y = a * ya(klo) + b * ya(khi) + &
        ((a**3 - a) * y2a(klo) + (b**3 - b) * y2a(khi)) * (h**2) / 6.0_dp

  end subroutine cubic_spline_eval


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

  !> Interpolate charges from non-uniform to uniform time grid using cubic spline interpolation
  !! \brief Given historical charges at non-uniform time points, interpolate to uniform grid
  !! \param dt_history Timestep history (most recent first)
  !! \param n_0, n_1, n_2, n_3, n_4, n_5 Charge arrays at non-uniform times
  !! \param ni_0, ni_1, ni_2, ni_3, ni_4, ni_5 Output: interpolated charges at uniform times
  !! \param nats Number of atoms
  subroutine prg_xlbo_interpolate_charges(dt_history, n_0, n_1, n_2, n_3, n_4, n_5, &
                                           ni_0, ni_1, ni_2, ni_3, ni_4, ni_5, nats)
    implicit none
    real(dp), intent(in) :: dt_history(5)
    real(dp), intent(in) :: n_0(:), n_1(:), n_2(:), n_3(:), n_4(:), n_5(:)
    real(dp), intent(out) :: ni_0(:), ni_1(:), ni_2(:), ni_3(:), ni_4(:), ni_5(:)
    integer, intent(in) :: nats

    real(dp) :: t(0:5), t_uniform(0:5)
    real(dp) :: dt_uniform
    integer :: i, j, iat
    real(dp) :: y(0:5), y2(0:5)  ! Function values and second derivatives for one atom

    ! Build non-uniform source grid (where charges are stored)
    t(0) = 0.0_dp
    t(1) = -dt_history(1)
    t(2) = t(1) - dt_history(2)
    t(3) = t(2) - dt_history(3)
    t(4) = t(3) - dt_history(4)
    t(5) = t(4) - dt_history(5)

    ! Build uniform target grid using fixed half-step spacing (dt/2)
    ! For adaptive timestepping, this provides a consistent interpolation target
    dt_uniform = 0.5_dp
    do i = 0, 5
      t_uniform(i) = -i * dt_uniform
    enddo

    ! Debug output for first atom to verify interpolation accuracy
    if (nats > 0) then
      write(*,*) "XLBO K=5 Interpolation Debug:"
      write(*,*) "  dt_history:", dt_history
      write(*,*) "  Source grid t:", t
      write(*,*) "  Target grid t_uniform:", t_uniform
      write(*,*) "  dt_uniform:", dt_uniform
    endif

    ! Interpolate each atom independently using cubic splines
    do iat = 1, nats
      ! Gather charges for this atom
      y(0) = n_0(iat)
      y(1) = n_1(iat)
      y(2) = n_2(iat)
      y(3) = n_3(iat)
      y(4) = n_4(iat)
      y(5) = n_5(iat)

      ! Compute spline second derivatives (natural boundary conditions)
      call cubic_spline_coeffs(t, y, 6, y2)

      ! Evaluate spline at uniform target points
      call cubic_spline_eval(t, y, y2, 6, t_uniform(0), ni_0(iat))
      call cubic_spline_eval(t, y, y2, 6, t_uniform(1), ni_1(iat))
      call cubic_spline_eval(t, y, y2, 6, t_uniform(2), ni_2(iat))
      call cubic_spline_eval(t, y, y2, 6, t_uniform(3), ni_3(iat))
      call cubic_spline_eval(t, y, y2, 6, t_uniform(4), ni_4(iat))
      call cubic_spline_eval(t, y, y2, 6, t_uniform(5), ni_5(iat))

      ! Debug output for first atom: check if coincident points match
      if (iat == 1) then
        write(*,*) "  First atom source charges:", y
        write(*,*) "  First atom interpolated charges:", ni_0(iat), ni_1(iat), ni_2(iat), &
                   ni_3(iat), ni_4(iat), ni_5(iat)
        do i = 0, 5
          do j = 0, 5
            if (abs(t_uniform(i) - t(j)) < 1.0e-10_dp) then
              write(*,*) "  Coincident point: t_uniform(",i,")=", t_uniform(i), &
                         " matches t(",j,")=", t(j)
              write(*,*) "    Expected value:", y(j), " Interpolated:", &
                         merge(ni_0(iat), merge(ni_1(iat), merge(ni_2(iat), merge(ni_3(iat), &
                         merge(ni_4(iat), ni_5(iat), i==5), i==4), i==3), i==2), i==1)
            endif
          enddo
        enddo
      endif
    enddo

  end subroutine prg_xlbo_interpolate_charges


  !> Interpolate charges from non-uniform to uniform time grid using cubic spline interpolation (11-point version for K=10)
  !! \brief Given historical charges at non-uniform time points, interpolate to uniform grid
  !! \param dt_history Timestep history (most recent first) - 10 elements
  !! \param n_0..n_10 Charge arrays at non-uniform times (11 points)
  !! \param ni_0..ni_10 Output: interpolated charges at uniform times
  !! \param nats Number of atoms
  subroutine prg_xlbo_interpolate_charges_K10(dt_history, n_0, n_1, n_2, n_3, n_4, n_5, &
                                               n_6, n_7, n_8, n_9, n_10, &
                                               ni_0, ni_1, ni_2, ni_3, ni_4, ni_5, &
                                               ni_6, ni_7, ni_8, ni_9, ni_10, nats)
    implicit none
    real(dp), intent(in) :: dt_history(10)
    real(dp), intent(in) :: n_0(:), n_1(:), n_2(:), n_3(:), n_4(:), n_5(:)
    real(dp), intent(in) :: n_6(:), n_7(:), n_8(:), n_9(:), n_10(:)
    real(dp), intent(out) :: ni_0(:), ni_1(:), ni_2(:), ni_3(:), ni_4(:), ni_5(:)
    real(dp), intent(out) :: ni_6(:), ni_7(:), ni_8(:), ni_9(:), ni_10(:)
    integer, intent(in) :: nats

    real(dp) :: t(0:10), t_uniform(0:10)
    real(dp) :: dt_uniform
    integer :: i, j, iat
    real(dp) :: y(0:10), y2(0:10)  ! Function values and second derivatives for one atom

    ! Build non-uniform source grid (where charges are stored)
    t(0) = 0.0_dp
    t(1) = -dt_history(1)
    do i = 2, 10
      t(i) = t(i-1) - dt_history(i)
    enddo

    ! Build uniform target grid using fixed half-step spacing (dt/2)
    ! For adaptive timestepping, this provides a consistent interpolation target
    dt_uniform = 0.5_dp
    do i = 0, 10
      t_uniform(i) = -i * dt_uniform
    enddo

    ! Debug output to verify interpolation accuracy
    if (nats > 0) then
      write(*,*) "XLBO K=10 Interpolation Debug:"
      write(*,*) "  dt_history:", dt_history
      write(*,*) "  Source grid t:", t
      write(*,*) "  Target grid t_uniform:", t_uniform
      write(*,*) "  dt_uniform:", dt_uniform
    endif

    ! Interpolate each atom independently using cubic splines
    do iat = 1, nats
      ! Gather charges for this atom
      y(0) = n_0(iat)
      y(1) = n_1(iat)
      y(2) = n_2(iat)
      y(3) = n_3(iat)
      y(4) = n_4(iat)
      y(5) = n_5(iat)
      y(6) = n_6(iat)
      y(7) = n_7(iat)
      y(8) = n_8(iat)
      y(9) = n_9(iat)
      y(10) = n_10(iat)

      ! Compute spline second derivatives (natural boundary conditions)
      call cubic_spline_coeffs(t, y, 11, y2)

      ! Evaluate spline at uniform target points
      call cubic_spline_eval(t, y, y2, 11, t_uniform(0), ni_0(iat))
      call cubic_spline_eval(t, y, y2, 11, t_uniform(1), ni_1(iat))
      call cubic_spline_eval(t, y, y2, 11, t_uniform(2), ni_2(iat))
      call cubic_spline_eval(t, y, y2, 11, t_uniform(3), ni_3(iat))
      call cubic_spline_eval(t, y, y2, 11, t_uniform(4), ni_4(iat))
      call cubic_spline_eval(t, y, y2, 11, t_uniform(5), ni_5(iat))
      call cubic_spline_eval(t, y, y2, 11, t_uniform(6), ni_6(iat))
      call cubic_spline_eval(t, y, y2, 11, t_uniform(7), ni_7(iat))
      call cubic_spline_eval(t, y, y2, 11, t_uniform(8), ni_8(iat))
      call cubic_spline_eval(t, y, y2, 11, t_uniform(9), ni_9(iat))
      call cubic_spline_eval(t, y, y2, 11, t_uniform(10), ni_10(iat))

      ! Debug output for first atom: check if coincident points match
      if (iat == 1) then
        write(*,*) "  First atom source charges:", y
        write(*,*) "  First atom interpolated charges:"
        write(*,*) "    ", ni_0(iat), ni_1(iat), ni_2(iat), ni_3(iat), ni_4(iat), ni_5(iat)
        write(*,*) "    ", ni_6(iat), ni_7(iat), ni_8(iat), ni_9(iat), ni_10(iat)
        do i = 0, 10
          do j = 0, 10
            if (abs(t_uniform(i) - t(j)) < 1.0e-10_dp) then
              write(*,'(A,I2,A,F8.4,A,I2,A,F8.4)') "  Coincident: t_uniform(", i, ")=", &
                      t_uniform(i), " matches t(", j, ")=", t(j)
              if (i == 0) write(*,'(A,F12.8,A,F12.8)') "    Expected:", y(j), " Got:", ni_0(iat)
              if (i == 1) write(*,'(A,F12.8,A,F12.8)') "    Expected:", y(j), " Got:", ni_1(iat)
              if (i == 2) write(*,'(A,F12.8,A,F12.8)') "    Expected:", y(j), " Got:", ni_2(iat)
              if (i == 3) write(*,'(A,F12.8,A,F12.8)') "    Expected:", y(j), " Got:", ni_3(iat)
              if (i == 4) write(*,'(A,F12.8,A,F12.8)') "    Expected:", y(j), " Got:", ni_4(iat)
              if (i == 5) write(*,'(A,F12.8,A,F12.8)') "    Expected:", y(j), " Got:", ni_5(iat)
              if (i == 6) write(*,'(A,F12.8,A,F12.8)') "    Expected:", y(j), " Got:", ni_6(iat)
              if (i == 7) write(*,'(A,F12.8,A,F12.8)') "    Expected:", y(j), " Got:", ni_7(iat)
              if (i == 8) write(*,'(A,F12.8,A,F12.8)') "    Expected:", y(j), " Got:", ni_8(iat)
              if (i == 9) write(*,'(A,F12.8,A,F12.8)') "    Expected:", y(j), " Got:", ni_9(iat)
              if (i == 10) write(*,'(A,F12.8,A,F12.8)') "    Expected:", y(j), " Got:", ni_10(iat)
            endif
          enddo
        enddo
      endif
    enddo

  end subroutine prg_xlbo_interpolate_charges_K10


  !> This routine integrates the dynamical variable "n"
  !! \param charges
  subroutine prg_xlbo_nint(charges,n,n_0,n_1,n_2,n_3,n_4,n_5,mdstep,xl,dt, &
                           n_6,n_7,n_8,n_9,n_10)
    implicit none
    real(dp), allocatable, intent(inout) :: n(:), n_0(:), n_1(:), n_2(:), n_3(:), n_4(:), n_5(:)
    real(dp), allocatable, intent(inout), optional :: n_6(:), n_7(:), n_8(:), n_9(:), n_10(:)
    real(dp), allocatable, intent(in) :: charges(:)
    type(xlbo_type), intent(inout) :: xl
    integer, intent(in) :: mdstep
    real(dp), intent(in), optional :: dt
    integer :: nats
    real(dp) :: kappa_use, alpha_use
    real(dp), allocatable :: ni_0(:), ni_1(:), ni_2(:), ni_3(:), ni_4(:), ni_5(:)
    real(dp), allocatable :: ni_6(:), ni_7(:), ni_8(:), ni_9(:), ni_10(:)
    logical :: use_interpolation, use_K10
    integer :: hist_idx
    real(dp) :: C0_use, C1_use, C2_use, C3_use, C4_use, C5_use, d_K_use
    real(dp) :: dt_n, dt_prev, r, P_n_coeff, P_n1_coeff, kappa_alpha_scale

    nats = size(charges,dim=1)

    ! Determine if we should use K=10
    use_K10 = xl%extended_history .and. present(n_6) .and. present(n_7) .and. &
              present(n_8) .and. present(n_9) .and. present(n_10)

    if(.not.allocated(n))then
      allocate(n(nats))
      allocate(n_0(nats))
      allocate(n_1(nats))
      allocate(n_2(nats))
      allocate(n_3(nats))
      allocate(n_4(nats))
      allocate(n_5(nats))
      if (use_K10) then
        allocate(n_6(nats))
        allocate(n_7(nats))
        allocate(n_8(nats))
        allocate(n_9(nats))
        allocate(n_10(nats))
      endif
    endif

    if(mdstep.le.1)then
      n = charges;
      n_0 = charges;
      n_1 = charges;
      n_2 = charges;
      n_3 = charges;
      n_4 = charges;
      n_5 = charges;
      if (use_K10) then
        n_6 = charges;
        n_7 = charges;
        n_8 = charges;
        n_9 = charges;
        n_10 = charges;
      endif
      xl%dt_history = 0.0_dp
      xl%nsteps_taken = 0
    endif

    ! Determine if we should use interpolation
    if (use_K10) then
      use_interpolation = present(dt) .and. xl%nsteps_taken >= 11
    else
      use_interpolation = present(dt) .and. xl%nsteps_taken >= 6
    endif

    ! Select parameters based on K value
    ! Use base kappa and alpha without dt scaling (scaling is in kappa_alpha_scale)
    if (use_K10) then
      kappa_use = kappa_K10
      alpha_use = alpha_K10
    else
      kappa_use = kappa
      alpha_use = alpha
    endif

    ! Compute variable timestep Verlet coefficients
    if (present(dt)) then
      ! Get timestep history for last two steps
      dt_n = xl%dt_history(1)      ! Most recent timestep
      dt_prev = xl%dt_history(2)   ! Previous timestep

      ! Check if we have valid history (dt_prev > 0)
      if (dt_prev > 1.0e-12_dp) then
        ! Compute Verlet coefficients for variable timesteps
        r = dt_n / dt_prev  ! Timestep ratio
        P_n_coeff = 1.0_dp + r
        P_n1_coeff = r

        ! Compute kappa and alpha scaling factor
        ! This comes from: 0.5 * dt_n * (dt_n + dt_{n-1})
        kappa_alpha_scale = 0.5_dp * dt_n * (dt_n + dt_prev)
      else
        ! First step or dt_prev not set yet - use uniform coefficients
        P_n_coeff = 2.0_dp
        P_n1_coeff = 1.0_dp
        kappa_alpha_scale = dt_n * dt_n  ! dt^2 for first step
      endif
    else
      ! Uniform timestep: standard Verlet coefficients
      P_n_coeff = 2.0_dp
      P_n1_coeff = 1.0_dp
      kappa_alpha_scale = 1.0_dp
    endif

    if (use_interpolation) then
      if (use_K10) then
        ! Allocate interpolated charge arrays for K=10
        allocate(ni_0(nats), ni_1(nats), ni_2(nats), ni_3(nats), ni_4(nats), ni_5(nats))
        allocate(ni_6(nats), ni_7(nats), ni_8(nats), ni_9(nats), ni_10(nats))

        ! Interpolate historical charges to uniform grid (11 points)
        call prg_xlbo_interpolate_charges_K10(xl%dt_history, n_0, n_1, n_2, n_3, n_4, n_5, &
                                               n_6, n_7, n_8, n_9, n_10, &
                                               ni_0, ni_1, ni_2, ni_3, ni_4, ni_5, &
                                               ni_6, ni_7, ni_8, ni_9, ni_10, nats)

        ! Integration using interpolated charges (K=10) with variable timestep Verlet
        n = P_n_coeff*ni_0 - P_n1_coeff*ni_1 + xl%cc*kappa_alpha_scale*kappa_use*(charges-n) &
             + kappa_alpha_scale*alpha_use*(C0_K10*ni_0+C1_K10*ni_1+C2_K10*ni_2+C3_K10*ni_3+C4_K10*ni_4+C5_K10*ni_5 &
                        +C6_K10*ni_6+C7_K10*ni_7+C8_K10*ni_8+C9_K10*ni_9+C10_K10*ni_10)

        deallocate(ni_0, ni_1, ni_2, ni_3, ni_4, ni_5)
        deallocate(ni_6, ni_7, ni_8, ni_9, ni_10)
      else
        ! K=5 with variable timesteps: choose method
        if (xl%use_variable_coeffs) then
          ! New method: Use variable timestep coefficients directly (no interpolation)

          ! Get coefficient index from timestep history
          hist_idx = get_K5_history_index(xl%dt_history(1:5))

          ! Lookup coefficients for this specific history pattern
          C0_use = XLBO_K5_C0(hist_idx)
          C1_use = XLBO_K5_C1(hist_idx)
          C2_use = XLBO_K5_C2(hist_idx)
          C3_use = XLBO_K5_C3(hist_idx)
          C4_use = C4  ! Fixed by normalization
          C5_use = C5  ! Fixed by normalization
          d_K_use = XLBO_K5_dK(hist_idx)

          ! Use pattern-specific alpha value (hybrid strategy: all 32 patterns stable)
          alpha_use = XLBO_K5_alpha(hist_idx)

          ! Integration using raw charges with variable coefficients and variable timestep Verlet
          n = P_n_coeff*n_0 - P_n1_coeff*n_1 + xl%cc*kappa_alpha_scale*kappa_use*(charges-n) &
               + kappa_alpha_scale*alpha_use*(C0_use*n_0+C1_use*n_1+C2_use*n_2+C3_use*n_3+C4_use*n_4+C5_use*n_5)

        else
          ! Old method: Interpolate to uniform grid, use standard coefficients
          allocate(ni_0(nats), ni_1(nats), ni_2(nats), ni_3(nats), ni_4(nats), ni_5(nats))

          ! Interpolate historical charges to uniform grid (6 points)
          call prg_xlbo_interpolate_charges(xl%dt_history, n_0, n_1, n_2, n_3, n_4, n_5, &
                                             ni_0, ni_1, ni_2, ni_3, ni_4, ni_5, nats)

          ! Integration using interpolated charges (K=5) with variable timestep Verlet
          n = P_n_coeff*ni_0 - P_n1_coeff*ni_1 + xl%cc*kappa_alpha_scale*kappa_use*(charges-n) &
               + kappa_alpha_scale*alpha_use*(C0*ni_0+C1*ni_1+C2*ni_2+C3*ni_3+C4*ni_4+C5*ni_5)

          deallocate(ni_0, ni_1, ni_2, ni_3, ni_4, ni_5)
        endif
      endif
    else
      ! Integration using raw charges with variable timestep Verlet
      if (use_K10) then
        n = P_n_coeff*n_0 - P_n1_coeff*n_1 + xl%cc*kappa_alpha_scale*kappa_use*(charges-n) &
             + kappa_alpha_scale*alpha_use*(C0_K10*n_0+C1_K10*n_1+C2_K10*n_2+C3_K10*n_3+C4_K10*n_4+C5_K10*n_5 &
                        +C6_K10*n_6+C7_K10*n_7+C8_K10*n_8+C9_K10*n_9+C10_K10*n_10)
      else
        n = P_n_coeff*n_0 - P_n1_coeff*n_1 + xl%cc*kappa_alpha_scale*kappa_use*(charges-n) &
             + kappa_alpha_scale*alpha_use*(C0*n_0+C1*n_1+C2*n_2+C3*n_3+C4*n_4+C5*n_5)
      endif
    endif

    ! Shift history arrays
    if (use_K10) then
      n_10 = n_9; n_9 = n_8; n_8 = n_7; n_7 = n_6; n_6 = n_5
    endif
    n_5 = n_4; n_4 = n_3; n_3 = n_2; n_2 = n_1; n_1 = n_0; n_0 = n

    ! Update timestep history if dt provided (dt is ratio: 1.0 for full step, 0.5 for half step)
    ! History stores ratios that indicate spacing relative to user timestep
    ! Interpolation always maps to fixed dt/2 uniform grid for stability
    if (present(dt)) then
      if (use_K10) then
        xl%dt_history(10) = xl%dt_history(9)
        xl%dt_history(9) = xl%dt_history(8)
        xl%dt_history(8) = xl%dt_history(7)
        xl%dt_history(7) = xl%dt_history(6)
        xl%dt_history(6) = xl%dt_history(5)
      endif
      xl%dt_history(5) = xl%dt_history(4)
      xl%dt_history(4) = xl%dt_history(3)
      xl%dt_history(3) = xl%dt_history(2)
      xl%dt_history(2) = xl%dt_history(1)
      xl%dt_history(1) = dt  ! Store ratio (1.0 = full user timestep, 0.5 = half user timestep)
      xl%nsteps_taken = xl%nsteps_taken + 1
    endif

  end subroutine prg_xlbo_nint


  !> This routine integrates the dynamical variable "n"
  !! \param charges
  subroutine prg_xlbo_nint_kernel(charges,n,n_0,n_1,n_2,n_3,n_4,n_5,mdstep,kernel,xl,dt, &
                                   n_6,n_7,n_8,n_9,n_10)
    implicit none
    real(dp), allocatable, intent(inout) :: n(:), n_0(:), n_1(:), n_2(:), n_3(:), n_4(:), n_5(:)
    real(dp), allocatable, intent(inout), optional :: n_6(:), n_7(:), n_8(:), n_9(:), n_10(:)
    real(dp), allocatable, intent(in) :: charges(:)
    real(dp), allocatable, intent(in) :: kernel(:,:)
    type(xlbo_type), intent(inout) :: xl
    integer, intent(in) :: mdstep
    real(dp), intent(in), optional :: dt
    integer :: nats
    real(dp) :: kappa_use, alpha_use
    real(dp), allocatable :: ni_0(:), ni_1(:), ni_2(:), ni_3(:), ni_4(:), ni_5(:)
    real(dp), allocatable :: ni_6(:), ni_7(:), ni_8(:), ni_9(:), ni_10(:)
    real(dp), allocatable :: KK0n(:)
    logical :: use_interpolation, use_K10
    integer :: hist_idx
    real(dp) :: C0_use, C1_use, C2_use, C3_use, C4_use, C5_use, d_K_use
    real(dp) :: dt_n, dt_prev, r, P_n_coeff, P_n1_coeff, kappa_alpha_scale

    nats = size(charges,dim=1)

    ! Determine if we should use K=10
    use_K10 = xl%extended_history .and. present(n_6) .and. present(n_7) .and. &
              present(n_8) .and. present(n_9) .and. present(n_10)

    if(.not.allocated(n))then
      allocate(n(nats))
      allocate(n_0(nats))
      allocate(n_1(nats))
      allocate(n_2(nats))
      allocate(n_3(nats))
      allocate(n_4(nats))
      allocate(n_5(nats))
      if (use_K10) then
        allocate(n_6(nats))
        allocate(n_7(nats))
        allocate(n_8(nats))
        allocate(n_9(nats))
        allocate(n_10(nats))
      endif
    endif

    if(mdstep.le.1)then
      n = charges;
      n_0 = charges;
      n_1 = charges;
      n_2 = charges;
      n_3 = charges;
      n_4 = charges;
      n_5 = charges;
      if (use_K10) then
        n_6 = charges;
        n_7 = charges;
        n_8 = charges;
        n_9 = charges;
        n_10 = charges;
      endif
      xl%dt_history = 0.0_dp
      xl%nsteps_taken = 0
    endif

    ! Determine if we should use interpolation
    if (use_K10) then
      use_interpolation = present(dt) .and. xl%nsteps_taken >= 11
    else
      use_interpolation = present(dt) .and. xl%nsteps_taken >= 6
    endif

    ! Select parameters based on K value
    ! Use base kappa and alpha without dt scaling (scaling is in kappa_alpha_scale)
    if (use_K10) then
      kappa_use = kappa_K10
      alpha_use = alpha_K10
    else
      kappa_use = kappa
      alpha_use = alpha
    endif

    ! Compute variable timestep Verlet coefficients
    if (present(dt)) then
      ! Get timestep history for last two steps
      dt_n = xl%dt_history(1)      ! Most recent timestep
      dt_prev = xl%dt_history(2)   ! Previous timestep

      ! Check if we have valid history (dt_prev > 0)
      if (dt_prev > 1.0e-12_dp) then
        ! Compute Verlet coefficients for variable timesteps
        r = dt_n / dt_prev  ! Timestep ratio
        P_n_coeff = 1.0_dp + r
        P_n1_coeff = r

        ! Compute kappa and alpha scaling factor
        ! This comes from: 0.5 * dt_n * (dt_n + dt_{n-1})
        kappa_alpha_scale = 0.5_dp * dt_n * (dt_n + dt_prev)
      else
        ! First step or dt_prev not set yet - use uniform coefficients
        P_n_coeff = 2.0_dp
        P_n1_coeff = 1.0_dp
        kappa_alpha_scale = dt_n * dt_n  ! dt^2 for first step
      endif
    else
      ! Uniform timestep: standard Verlet coefficients
      P_n_coeff = 2.0_dp
      P_n1_coeff = 1.0_dp
      kappa_alpha_scale = 1.0_dp
    endif

    ! From developper's code
    !   dn2dt2 = -MATMUL(KK0,(q-n))
    !   n = 2*n_0 - n_1 + kappa*dn2dt2 +
    !   alpha*(C0*n_0+C1*n_1+C2*n_2+C3*n_3+C4*n_4+C5*n_5+C6*n_6)
    !   n_6 = n_5; n_5 = n_4; n_4 = n_3; n_3 = n_2; n_2 = n_1; n_1 = n_0; n_0 = n

    if (use_interpolation) then
      if (use_K10) then
        ! Allocate interpolated charge arrays for K=10
        allocate(ni_0(nats), ni_1(nats), ni_2(nats), ni_3(nats), ni_4(nats), ni_5(nats))
        allocate(ni_6(nats), ni_7(nats), ni_8(nats), ni_9(nats), ni_10(nats))

        ! Interpolate historical charges to uniform grid (11 points)
        call prg_xlbo_interpolate_charges_K10(xl%dt_history, n_0, n_1, n_2, n_3, n_4, n_5, &
                                               n_6, n_7, n_8, n_9, n_10, &
                                               ni_0, ni_1, ni_2, ni_3, ni_4, ni_5, &
                                               ni_6, ni_7, ni_8, ni_9, ni_10, nats)

        ! Integration using interpolated charges (K=10) with variable timestep Verlet
        n = P_n_coeff*ni_0 - P_n1_coeff*ni_1 - kappa_alpha_scale*kappa_use*matmul(kernel,(charges-n)) &
             + kappa_alpha_scale*alpha_use*(C0_K10*ni_0+C1_K10*ni_1+C2_K10*ni_2+C3_K10*ni_3+C4_K10*ni_4+C5_K10*ni_5 &
                        +C6_K10*ni_6+C7_K10*ni_7+C8_K10*ni_8+C9_K10*ni_9+C10_K10*ni_10)

        deallocate(ni_0, ni_1, ni_2, ni_3, ni_4, ni_5)
        deallocate(ni_6, ni_7, ni_8, ni_9, ni_10)
      else
        ! K=5 with variable timesteps: choose method
        if (xl%use_variable_coeffs) then
          ! New method: Use variable timestep coefficients directly (no interpolation)

          ! Get coefficient index from timestep history
          hist_idx = get_K5_history_index(xl%dt_history(1:5))

          ! Lookup coefficients for this specific history pattern
          C0_use = XLBO_K5_C0(hist_idx)
          C1_use = XLBO_K5_C1(hist_idx)
          C2_use = XLBO_K5_C2(hist_idx)
          C3_use = XLBO_K5_C3(hist_idx)
          C4_use = XLBO_K5_C4(hist_idx)
          C5_use = XLBO_K5_C5(hist_idx)
          d_K_use = XLBO_K5_dK(hist_idx)

          ! Use pattern-specific alpha value (hybrid strategy: all 32 patterns stable)
          alpha_use = XLBO_K5_alpha(hist_idx)

          ! Integration using raw charges with variable coefficients and variable timestep Verlet
          n = P_n_coeff*n_0 - P_n1_coeff*n_1 - kappa_alpha_scale*kappa_use*matmul(kernel,(charges-n)) &
               + kappa_alpha_scale*alpha_use*(C0_use*n_0+C1_use*n_1+C2_use*n_2+C3_use*n_3+C4_use*n_4+C5_use*n_5)

        else
          ! Old method: Interpolate to uniform grid, use standard coefficients
          allocate(ni_0(nats), ni_1(nats), ni_2(nats), ni_3(nats), ni_4(nats), ni_5(nats))

          ! Interpolate historical charges to uniform grid (6 points)
          call prg_xlbo_interpolate_charges(xl%dt_history, n_0, n_1, n_2, n_3, n_4, n_5, &
                                             ni_0, ni_1, ni_2, ni_3, ni_4, ni_5, nats)

          ! Integration using interpolated charges (K=5) with variable timestep Verlet
          n = P_n_coeff*ni_0 - P_n1_coeff*ni_1 - kappa_alpha_scale*kappa_use*matmul(kernel,(charges-n)) &
               + kappa_alpha_scale*alpha_use*(C0*ni_0+C1*ni_1+C2*ni_2+C3*ni_3+C4*ni_4+C5*ni_5)

          deallocate(ni_0, ni_1, ni_2, ni_3, ni_4, ni_5)
        endif
      endif
    else
      ! Integration using raw charges with variable timestep Verlet
      if (use_K10) then
        n = P_n_coeff*n_0 - P_n1_coeff*n_1 - kappa_alpha_scale*kappa_use*matmul(kernel,(charges-n)) &
             + kappa_alpha_scale*alpha_use*(C0_K10*n_0+C1_K10*n_1+C2_K10*n_2+C3_K10*n_3+C4_K10*n_4+C5_K10*n_5 &
                        +C6_K10*n_6+C7_K10*n_7+C8_K10*n_8+C9_K10*n_9+C10_K10*n_10)
      else
        n = P_n_coeff*n_0 - P_n1_coeff*n_1 - kappa_alpha_scale*kappa_use*matmul(kernel,(charges-n)) &
             + kappa_alpha_scale*alpha_use*(C0*n_0+C1*n_1+C2*n_2+C3*n_3+C4*n_4+C5*n_5)
      endif
    endif

    ! Shift history arrays
    if (use_K10) then
      n_10 = n_9; n_9 = n_8; n_8 = n_7; n_7 = n_6; n_6 = n_5
    endif
    n_5 = n_4; n_4 = n_3; n_3 = n_2; n_2 = n_1; n_1 = n_0; n_0 = n

    ! Update timestep history if dt provided (dt is ratio: 1.0 for full step, 0.5 for half step)
    ! History stores ratios that indicate spacing relative to user timestep
    ! Interpolation always maps to fixed dt/2 uniform grid for stability
    if (present(dt)) then
      if (use_K10) then
        xl%dt_history(10) = xl%dt_history(9)
        xl%dt_history(9) = xl%dt_history(8)
        xl%dt_history(8) = xl%dt_history(7)
        xl%dt_history(7) = xl%dt_history(6)
        xl%dt_history(6) = xl%dt_history(5)
      endif
      xl%dt_history(5) = xl%dt_history(4)
      xl%dt_history(4) = xl%dt_history(3)
      xl%dt_history(3) = xl%dt_history(2)
      xl%dt_history(2) = xl%dt_history(1)
      xl%dt_history(1) = dt  ! Store ratio (1.0 = full user timestep, 0.5 = half user timestep)
      xl%nsteps_taken = xl%nsteps_taken + 1
    endif

  end subroutine prg_xlbo_nint_kernel


  !> This routine integrates the dynamical variable "n"
  !! \brief In this case we are passing a premultiplied ressidue x kernel
  !! tis is done to avoid rank-specific multiplication within this routine.
  !! \param charges
  subroutine prg_xlbo_nint_kernelTimesRes(charges,n,n_0,n_1,n_2,n_3,n_4,n_5,mdstep,kernelTimesRes,xl,dt, &
                                          n_6,n_7,n_8,n_9,n_10)
    implicit none
    real(dp), allocatable, intent(inout) :: n(:), n_0(:), n_1(:), n_2(:), n_3(:), n_4(:), n_5(:)
    real(dp), allocatable, intent(inout), optional :: n_6(:), n_7(:), n_8(:), n_9(:), n_10(:)
    real(dp), allocatable, intent(in) :: charges(:)
    real(dp), allocatable, intent(in) :: kernelTimesRes(:)
    type(xlbo_type), intent(inout) :: xl
    integer, intent(in) :: mdstep
    real(dp), intent(in), optional :: dt
    integer :: nats
    real(dp) :: kappa_use, alpha_use
    real(dp), allocatable :: ni_0(:), ni_1(:), ni_2(:), ni_3(:), ni_4(:), ni_5(:)
    real(dp), allocatable :: ni_6(:), ni_7(:), ni_8(:), ni_9(:), ni_10(:)
    logical :: use_interpolation, use_K10
    integer :: hist_idx
    real(dp) :: C0_use, C1_use, C2_use, C3_use, C4_use, C5_use, d_K_use
    real(dp) :: dt_n, dt_prev, r, P_n_coeff, P_n1_coeff, kappa_alpha_scale

    nats = size(charges,dim=1)

    ! Determine if we should use K=10
    use_K10 = xl%extended_history .and. present(n_6) .and. present(n_7) .and. &
              present(n_8) .and. present(n_9) .and. present(n_10)

    if(.not.allocated(n))then
      allocate(n(nats))
      allocate(n_0(nats))
      allocate(n_1(nats))
      allocate(n_2(nats))
      allocate(n_3(nats))
      allocate(n_4(nats))
      allocate(n_5(nats))
      if (use_K10) then
        allocate(n_6(nats))
        allocate(n_7(nats))
        allocate(n_8(nats))
        allocate(n_9(nats))
        allocate(n_10(nats))
      endif
    endif

    if(mdstep.le.1)then
      n = charges;
      n_0 = charges;
      n_1 = charges;
      n_2 = charges;
      n_3 = charges;
      n_4 = charges;
      n_5 = charges;
      if (use_K10) then
        n_6 = charges;
        n_7 = charges;
        n_8 = charges;
        n_9 = charges;
        n_10 = charges;
      endif
      xl%dt_history = 0.0_dp
      xl%nsteps_taken = 0
    endif

    ! Determine if we should use interpolation
    if (use_K10) then
      use_interpolation = present(dt) .and. xl%nsteps_taken >= 11
    else
      use_interpolation = present(dt) .and. xl%nsteps_taken >= 6
    endif

    ! Select parameters based on K value
    ! Use base kappa and alpha without dt scaling (scaling is in kappa_alpha_scale)
    if (use_K10) then
      kappa_use = kappa_K10
      alpha_use = alpha_K10
    else
      kappa_use = kappa
      alpha_use = alpha
    endif

    ! Compute variable timestep Verlet coefficients
    if (present(dt)) then
      ! Get timestep history for last two steps
      dt_n = xl%dt_history(1)      ! Most recent timestep
      dt_prev = xl%dt_history(2)   ! Previous timestep

      ! Check if we have valid history (dt_prev > 0)
      if (dt_prev > 1.0e-12_dp) then
        ! Compute Verlet coefficients for variable timesteps
        r = dt_n / dt_prev  ! Timestep ratio
        P_n_coeff = 1.0_dp + r
        P_n1_coeff = r

        ! Compute kappa and alpha scaling factor
        ! This comes from: 0.5 * dt_n * (dt_n + dt_{n-1})
        kappa_alpha_scale = 0.5_dp * dt_n * (dt_n + dt_prev)
      else
        ! First step or dt_prev not set yet - use uniform coefficients
        P_n_coeff = 2.0_dp
        P_n1_coeff = 1.0_dp
        kappa_alpha_scale = dt_n * dt_n  ! dt^2 for first step
      endif
    else
      ! Uniform timestep: standard Verlet coefficients
      P_n_coeff = 2.0_dp
      P_n1_coeff = 1.0_dp
      kappa_alpha_scale = 1.0_dp
    endif

    if (use_interpolation) then
      if (use_K10) then
        ! Allocate interpolated charge arrays for K=10
        allocate(ni_0(nats), ni_1(nats), ni_2(nats), ni_3(nats), ni_4(nats), ni_5(nats))
        allocate(ni_6(nats), ni_7(nats), ni_8(nats), ni_9(nats), ni_10(nats))

        ! Interpolate historical charges to uniform grid (11 points)
        call prg_xlbo_interpolate_charges_K10(xl%dt_history, n_0, n_1, n_2, n_3, n_4, n_5, &
                                               n_6, n_7, n_8, n_9, n_10, &
                                               ni_0, ni_1, ni_2, ni_3, ni_4, ni_5, &
                                               ni_6, ni_7, ni_8, ni_9, ni_10, nats)

        ! Integration using interpolated charges (K=10) with variable timestep Verlet
        n = P_n_coeff*ni_0 - P_n1_coeff*ni_1 - kappa_alpha_scale*kappa_use*kernelTimesRes &
             & + kappa_alpha_scale*alpha_use*(C0_K10*ni_0+C1_K10*ni_1+C2_K10*ni_2+C3_K10*ni_3+C4_K10*ni_4+C5_K10*ni_5 &
                          +C6_K10*ni_6+C7_K10*ni_7+C8_K10*ni_8+C9_K10*ni_9+C10_K10*ni_10)

        deallocate(ni_0, ni_1, ni_2, ni_3, ni_4, ni_5)
        deallocate(ni_6, ni_7, ni_8, ni_9, ni_10)
      else
        ! K=5 with variable timesteps: choose method
        if (xl%use_variable_coeffs) then
          ! New method: Use variable timestep coefficients directly (no interpolation)

          ! Get coefficient index from timestep history
          hist_idx = get_K5_history_index(xl%dt_history(1:5))

          ! Lookup coefficients for this specific history pattern
          C0_use = XLBO_K5_C0(hist_idx)
          C1_use = XLBO_K5_C1(hist_idx)
          C2_use = XLBO_K5_C2(hist_idx)
          C3_use = XLBO_K5_C3(hist_idx)
          d_K_use = XLBO_K5_dK(hist_idx)

          ! Use pattern-specific alpha value (hybrid strategy: all 32 patterns stable)
          alpha_use = XLBO_K5_alpha(hist_idx)

          ! Integration using raw charges with variable coefficients and variable timestep Verlet
          n = P_n_coeff*n_0 - P_n1_coeff*n_1 - kappa_alpha_scale*kappa_use*kernelTimesRes &
               & + kappa_alpha_scale*alpha_use*(C0_use*n_0+C1_use*n_1+C2_use*n_2+C3_use*n_3+C4*n_4+C5*n_5)

        else
          ! Old method: Interpolate to uniform grid, use standard coefficients
          allocate(ni_0(nats), ni_1(nats), ni_2(nats), ni_3(nats), ni_4(nats), ni_5(nats))

          ! Interpolate historical charges to uniform grid (6 points)
          call prg_xlbo_interpolate_charges(xl%dt_history, n_0, n_1, n_2, n_3, n_4, n_5, &
                                             ni_0, ni_1, ni_2, ni_3, ni_4, ni_5, nats)

          ! Integration using interpolated charges (K=5) with variable timestep Verlet
          n = P_n_coeff*ni_0 - P_n1_coeff*ni_1 - kappa_alpha_scale*kappa_use*kernelTimesRes &
               & + kappa_alpha_scale*alpha_use*(C0*ni_0+C1*ni_1+C2*ni_2+C3*ni_3+C4*ni_4+C5*ni_5)

          deallocate(ni_0, ni_1, ni_2, ni_3, ni_4, ni_5)
        endif
      endif
    else
      ! Integration using raw charges with variable timestep Verlet
      if (use_K10) then
        n = P_n_coeff*n_0 - P_n1_coeff*n_1 - kappa_alpha_scale*kappa_use*kernelTimesRes &
             & + kappa_alpha_scale*alpha_use*(C0_K10*n_0+C1_K10*n_1+C2_K10*n_2+C3_K10*n_3+C4_K10*n_4+C5_K10*n_5 &
                          +C6_K10*n_6+C7_K10*n_7+C8_K10*n_8+C9_K10*n_9+C10_K10*n_10)
      else
        n = P_n_coeff*n_0 - P_n1_coeff*n_1 - kappa_alpha_scale*kappa_use*kernelTimesRes &
             & + kappa_alpha_scale*alpha_use*(C0*n_0+C1*n_1+C2*n_2+C3*n_3+C4*n_4+C5*n_5)
      endif
    endif

    ! Shift history arrays
    if (use_K10) then
      n_10 = n_9; n_9 = n_8; n_8 = n_7; n_7 = n_6; n_6 = n_5
    endif
    n_5 = n_4; n_4 = n_3; n_3 = n_2; n_2 = n_1; n_1 = n_0; n_0 = n

    ! Update timestep history if dt provided (dt is ratio: 1.0 for full step, 0.5 for half step)
    ! History stores ratios that indicate spacing relative to user timestep
    ! Interpolation always maps to fixed dt/2 uniform grid for stability
    if (present(dt)) then
      if (use_K10) then
        xl%dt_history(10) = xl%dt_history(9)
        xl%dt_history(9) = xl%dt_history(8)
        xl%dt_history(8) = xl%dt_history(7)
        xl%dt_history(7) = xl%dt_history(6)
        xl%dt_history(6) = xl%dt_history(5)
      endif
      xl%dt_history(5) = xl%dt_history(4)
      xl%dt_history(4) = xl%dt_history(3)
      xl%dt_history(3) = xl%dt_history(2)
      xl%dt_history(2) = xl%dt_history(1)
      xl%dt_history(1) = dt  ! Store ratio (1.0 = full user timestep, 0.5 = half user timestep)
      xl%nsteps_taken = xl%nsteps_taken + 1
    endif

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
