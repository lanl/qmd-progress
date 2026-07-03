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

  !> Coefficients for K=5 dissipation (6 history points)
  real(dp), parameter :: C0_K5 = -6.0_dp
  real(dp), parameter :: C1_K5 = 14.0_dp
  real(dp), parameter :: C2_K5 = -8.0_dp
  real(dp), parameter :: C3_K5 = -3.0_dp
  real(dp), parameter :: C4_K5 = 4.0_dp
  real(dp), parameter :: C5_K5 = -1.0_dp
  real(dp), parameter :: kappa_K5 = 1.82_dp
  real(dp), parameter :: alpha_K5 = 0.018_dp

  !> Coefficients for K=9 dissipation (10 history points)
  real(dp), parameter :: C0_K9 = -286.0_dp
  real(dp), parameter :: C1_K9 = 858.0_dp
  real(dp), parameter :: C2_K9 = -936.0_dp
  real(dp), parameter :: C3_K9 = 364.0_dp
  real(dp), parameter :: C4_K9 = 168.0_dp
  real(dp), parameter :: C5_K9 = -300.0_dp
  real(dp), parameter :: C6_K9 = 184.0_dp
  real(dp), parameter :: C7_K9 = -63.0_dp
  real(dp), parameter :: C8_K9 = 12.0_dp
  real(dp), parameter :: C9_K9 = -1.0_dp
  real(dp), parameter :: kappa_K9 = 1.89_dp
  real(dp), parameter :: alpha_K9 = 0.00012_dp

  !> Legacy parameters (for backward compatibility)
  real(dp), parameter :: kappa = 1.82_dp
  real(dp), parameter :: alpha = 0.018_dp
  real(dp), parameter :: cc = 0.9_dp;        ! Scaled prg_delta kernel

  !> Adaptive coefficient type for variable timesteps
  type :: xlbo_coeff_set
    real(dp) :: c0, c1, c2, c3, c4, c5
    real(dp) :: d5_tilde
  end type xlbo_coeff_set

  !> Lookup table for all 32 timestep patterns (K=5)
  !! Pattern encoding: 5-bit binary where bit i = 0 for δt/2, 1 for δt
  !! All patterns rigorously derived from Niklasson et al., JCP 130, 214109 (2009)
  type(xlbo_coeff_set), parameter :: coeff_table(0:31) = [ &
    xlbo_coeff_set(-6.0_dp, 14.0_dp, -8.0_dp, -3.0_dp, 4.0_dp, -1.0_dp, -0.750_dp), & ! 0: [½,½,½,½,½]
    xlbo_coeff_set(-0.9_dp, 3.0_dp, -0.6_dp, -4.6_dp, 4.0_dp, -1.0_dp, -0.286_dp), & ! 1: [1,½,½,½,½]
    xlbo_coeff_set(-1.9_dp, 3.0_dp, 0.7_dp, -4.8_dp, 4.0_dp, -1.0_dp, -0.393_dp), & ! 2: [½,1,½,½,½]
    xlbo_coeff_set(-0.3_dp, 0.5_dp, 2.1_dp, -5.3_dp, 4.0_dp, -1.0_dp, -0.171_dp), & ! 3: [1,1,½,½,½]
    xlbo_coeff_set(-1.3_dp, 1.3_dp, 1.7_dp, -4.7_dp, 4.0_dp, -1.0_dp, -0.333_dp), & ! 4: [½,½,1,½,½]
    xlbo_coeff_set(0.1_dp, -1.9_dp, 3.7_dp, -4.8_dp, 4.0_dp, -1.0_dp, -0.068_dp), & ! 5: [1,½,1,½,½]
    xlbo_coeff_set(0.9_dp, -2.1_dp, 3.1_dp, -4.8_dp, 4.0_dp, -1.0_dp, 0.016_dp), & ! 6: [½,1,1,½,½]
    xlbo_coeff_set(0.4_dp, -1.6_dp, 3.1_dp, -4.8_dp, 4.0_dp, -1.0_dp, 0.078_dp), & ! 7: [1,1,1,½,½]
    xlbo_coeff_set(22.0_dp, -64.0_dp, 67.0_dp, -28.0_dp, 4.0_dp, -1.0_dp, 2.000_dp), & ! 8: [½,½,½,1,½]
    xlbo_coeff_set(4.9_dp, -31.0_dp, 47.3_dp, -24.2_dp, 4.0_dp, -1.0_dp, 1.143_dp), & ! 9: [1,½,½,1,½]
    xlbo_coeff_set(16.0_dp, -31.0_dp, 34.0_dp, -22.0_dp, 4.0_dp, -1.0_dp, 2.250_dp), & ! 10: [½,1,½,1,½]
    xlbo_coeff_set(4.2_dp, -14.3_dp, 26.7_dp, -19.6_dp, 4.0_dp, -1.0_dp, 1.381_dp), & ! 11: [1,1,½,1,½]
    xlbo_coeff_set(44.3_dp, -110.3_dp, 79.3_dp, -16.3_dp, 4.0_dp, -1.0_dp, 5.083_dp), & ! 12: [½,½,1,1,½]
    xlbo_coeff_set(9.6_dp, -46.3_dp, 48.0_dp, -14.3_dp, 4.0_dp, -1.0_dp, 2.714_dp), & ! 13: [1,½,1,1,½]
    xlbo_coeff_set(28.4_dp, -50.6_dp, 32.8_dp, -13.6_dp, 4.0_dp, -1.0_dp, 4.703_dp), & ! 14: [½,1,1,1,½]
    xlbo_coeff_set(7.5_dp, -21.7_dp, 23.5_dp, -12.3_dp, 4.0_dp, -1.0_dp, 2.832_dp), & ! 15: [1,1,1,1,½]
    xlbo_coeff_set(-62.0_dp, 160.0_dp, -133.0_dp, 32.0_dp, 4.0_dp, -1.0_dp, -7.000_dp), & ! 16: [½,½,½,½,1]
    xlbo_coeff_set(-10.5_dp, 53.0_dp, -63.0_dp, 17.5_dp, 4.0_dp, -1.0_dp, -3.000_dp), & ! 17: [1,½,½,½,1]
    xlbo_coeff_set(-29.4_dp, 53.0_dp, -40.3_dp, 13.7_dp, 4.0_dp, -1.0_dp, -4.893_dp), & ! 18: [½,1,½,½,1]
    xlbo_coeff_set(-6.6_dp, 19.2_dp, -23.8_dp, 8.1_dp, 4.0_dp, -1.0_dp, -2.540_dp), & ! 19: [1,1,½,½,1]
    xlbo_coeff_set(-63.0_dp, 147.0_dp, -94.0_dp, 7.0_dp, 4.0_dp, -1.0_dp, -8.250_dp), & ! 20: [½,½,1,½,1]
    xlbo_coeff_set(-11.3_dp, 47.8_dp, -42.8_dp, 3.4_dp, 4.0_dp, -1.0_dp, -3.729_dp), & ! 21: [1,½,1,½,1]
    xlbo_coeff_set(-30.7_dp, 52.2_dp, -27.1_dp, 2.6_dp, 4.0_dp, -1.0_dp, -5.781_dp), & ! 22: [½,1,1,½,1]
    xlbo_coeff_set(-7.1_dp, 18.6_dp, -15.4_dp, 0.9_dp, 4.0_dp, -1.0_dp, -3.090_dp), & ! 23: [1,1,1,½,1]
    xlbo_coeff_set(-98.0_dp, 245.0_dp, -192.0_dp, 42.0_dp, 4.0_dp, -1.0_dp, -11.750_dp), & ! 24: [½,½,½,1,1]
    xlbo_coeff_set(-14.7_dp, 68.0_dp, -73.1_dp, 16.9_dp, 4.0_dp, -1.0_dp, -4.571_dp), & ! 25: [1,½,½,1,1]
    xlbo_coeff_set(-39.0_dp, 68.0_dp, -44.0_dp, 12.0_dp, 4.0_dp, -1.0_dp, -7.000_dp), & ! 26: [½,1,½,1,1]
    xlbo_coeff_set(-7.8_dp, 21.0_dp, -19.8_dp, 3.7_dp, 4.0_dp, -1.0_dp, -3.333_dp), & ! 27: [1,1,½,1,1]
    xlbo_coeff_set(-75.7_dp, 170.7_dp, -102.7_dp, 4.7_dp, 4.0_dp, -1.0_dp, -10.667_dp), & ! 28: [½,½,1,1,1]
    xlbo_coeff_set(-11.8_dp, 44.8_dp, -35.7_dp, -0.3_dp, 4.0_dp, -1.0_dp, -4.325_dp), & ! 29: [1,½,1,1,1]
    xlbo_coeff_set(-30.0_dp, 49.0_dp, -21.0_dp, -1.0_dp, 4.0_dp, -1.0_dp, -6.250_dp), & ! 30: [½,1,1,1,1]
    xlbo_coeff_set(-6.0_dp, 14.0_dp, -8.0_dp, -3.0_dp, 4.0_dp, -1.0_dp, -3.000_dp) ]    ! 31: [1,1,1,1,1]

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
    real(dp) :: dt_history(5)
    integer :: nsteps_taken

  end type xlbo_type

  public :: prg_parse_xlbo, prg_xlbo_nint, prg_xlbo_nint_kernel, prg_xlbo_fcoulupdate
  public :: prg_xlbo_nint_kernelTimesRes

contains

  !> Encode timestep history into pattern index (0-31)
  subroutine encode_timestep_pattern(dt_history, dt_base, pattern_idx, threshold)
    implicit none
    real(dp), intent(in) :: dt_history(5)
    real(dp), intent(in) :: dt_base
    integer, intent(out) :: pattern_idx
    real(dp), intent(in), optional :: threshold
    real(dp) :: tol
    integer :: i, bit_value
    real(dp) :: dist_full, dist_half

    if (present(threshold)) then
      tol = threshold
    else
      tol = 0.25_dp * dt_base
    endif

    pattern_idx = 0

    do i = 0, 4
      dist_full = abs(dt_history(i+1) - dt_base)
      dist_half = abs(dt_history(i+1) - 0.5_dp*dt_base)

      if (dist_full < tol) then
        bit_value = 1
      else if (dist_half < tol) then
        bit_value = 0
      else
        write(*,'(A,I1,A,ES14.6,A,ES14.6,A,ES14.6)') &
          'ERROR: dt_history(', i+1, ') = ', dt_history(i+1), &
          ' is neither dt_base (', dt_base, ') nor dt_base/2 (', 0.5_dp*dt_base, ')'
        write(*,'(A,ES14.6)') '  Distance to dt_base:   ', dist_full
        write(*,'(A,ES14.6)') '  Distance to dt_base/2: ', dist_half
        write(*,'(A,ES14.6)') '  Tolerance:             ', tol
        stop 'XLBO: invalid timestep pattern'
      endif

      pattern_idx = pattern_idx + bit_value * 2**i
    enddo

    if (pattern_idx < 0 .or. pattern_idx > 31) then
      write(*,'(A,I0)') 'ERROR: pattern_idx = ', pattern_idx
      stop 'XLBO: pattern index out of range'
    endif

  end subroutine encode_timestep_pattern


  !> Get adaptive coefficients for timestep pattern
  subroutine get_xlbo_coefficients_adaptive(dt_history, dt_base, &
                                             c0, c1, c2, c3, c4, c5, &
                                             alpha_out, alpha_base)
    implicit none
    real(dp), intent(in) :: dt_history(5)
    real(dp), intent(in) :: dt_base
    real(dp), intent(out) :: c0, c1, c2, c3, c4, c5
    real(dp), intent(out) :: alpha_out
    real(dp), intent(in) :: alpha_base

    integer :: pattern_idx
    type(xlbo_coeff_set) :: coeffs

    call encode_timestep_pattern(dt_history, dt_base, pattern_idx)

    coeffs = coeff_table(pattern_idx)

    c0 = coeffs%c0
    c1 = coeffs%c1
    c2 = coeffs%c2
    c3 = coeffs%c3
    c4 = coeffs%c4
    c5 = coeffs%c5

    alpha_out = alpha_base

  end subroutine get_xlbo_coefficients_adaptive

  !> Select dissipation coefficients based on available history
  !!
  subroutine prg_xlbo_select_coeffs(nsteps, use_K9, kappa_out, alpha_out, &
       C0, C1, C2, C3, C4, C5, C6, C7, C8, C9)
    implicit none
    integer, intent(in) :: nsteps
    logical, intent(out) :: use_K9
    real(dp), intent(out) :: kappa_out, alpha_out
    real(dp), intent(out) :: C0, C1, C2, C3, C4, C5
    real(dp), intent(out), optional :: C6, C7, C8, C9

    ! Use K=9 coefficients if we have 10+ history points available
    if (nsteps >= 10) then
      use_K9 = .true.
      kappa_out = kappa_K9
      alpha_out = alpha_K9
      C0 = C0_K9
      C1 = C1_K9
      C2 = C2_K9
      C3 = C3_K9
      C4 = C4_K9
      C5 = C5_K9
      if (present(C6)) C6 = C6_K9
      if (present(C7)) C7 = C7_K9
      if (present(C8)) C8 = C8_K9
      if (present(C9)) C9 = C9_K9
    else
      ! Use K=5 coefficients for early steps
      use_K9 = .false.
      kappa_out = kappa_K5
      alpha_out = alpha_K5
      C0 = C0_K5
      C1 = C1_K5
      C2 = C2_K5
      C3 = C3_K5
      C4 = C4_K5
      C5 = C5_K5
      if (present(C6)) C6 = 0.0_dp
      if (present(C7)) C7 = 0.0_dp
      if (present(C8)) C8 = 0.0_dp
      if (present(C9)) C9 = 0.0_dp
    endif
  end subroutine prg_xlbo_select_coeffs

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

    !Logicals

    !Integers
    xlbo%verbose = valvector_int(1)
    xlbo%minit = valvector_int(2)
    xlbo%maxscfiter = valvector_int(3)
    xlbo%maxscfinititer = valvector_int(4)

    !Initialize timestep history
    xlbo%dt_history = 0.0_dp
    xlbo%nsteps_taken = 0

  end subroutine prg_parse_xlbo


  !> Interpolate charges from non-uniform to uniform time grid using Lagrange interpolation
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
    real(dp) :: L(0:5,0:5)  ! Lagrange basis functions: L(k,i) = L_i(t_uniform(k))
    real(dp) :: dt_uniform
    integer :: i, j, k
    real(dp) :: numer, denom

    ! Build non-uniform source grid (where charges are stored)
    t(0) = 0.0_dp
    t(1) = -dt_history(1)
    t(2) = t(1) - dt_history(2)
    t(3) = t(2) - dt_history(3)
    t(4) = t(3) - dt_history(4)
    t(5) = t(4) - dt_history(5)

    ! Build uniform target grid using most recent timestep
    dt_uniform = dt_history(1)
    do i = 0, 5
      t_uniform(i) = -i * dt_uniform
    enddo

    ! Compute Lagrange basis functions L_i(t_uniform(k)) for all i,k
    ! L_i(t) = product over j≠i of [(t - t(j)) / (t(i) - t(j))]
    do k = 0, 5  ! For each target time
      do i = 0, 5  ! For each basis function
        L(k,i) = 1.0_dp
        do j = 0, 5
          if (j /= i) then
            numer = t_uniform(k) - t(j)
            denom = t(i) - t(j)
            L(k,i) = L(k,i) * (numer / denom)
          endif
        enddo
      enddo
    enddo

    ! Interpolate charges for each atom using the precomputed basis
    ! ni(k) = sum over i of [n_i * L(k,i)]
    ni_0 = L(0,0)*n_0 + L(0,1)*n_1 + L(0,2)*n_2 + L(0,3)*n_3 + L(0,4)*n_4 + L(0,5)*n_5
    ni_1 = L(1,0)*n_0 + L(1,1)*n_1 + L(1,2)*n_2 + L(1,3)*n_3 + L(1,4)*n_4 + L(1,5)*n_5
    ni_2 = L(2,0)*n_0 + L(2,1)*n_1 + L(2,2)*n_2 + L(2,3)*n_3 + L(2,4)*n_4 + L(2,5)*n_5
    ni_3 = L(3,0)*n_0 + L(3,1)*n_1 + L(3,2)*n_2 + L(3,3)*n_3 + L(3,4)*n_4 + L(3,5)*n_5
    ni_4 = L(4,0)*n_0 + L(4,1)*n_1 + L(4,2)*n_2 + L(4,3)*n_3 + L(4,4)*n_4 + L(4,5)*n_5
    ni_5 = L(5,0)*n_0 + L(5,1)*n_1 + L(5,2)*n_2 + L(5,3)*n_3 + L(5,4)*n_4 + L(5,5)*n_5

  end subroutine prg_xlbo_interpolate_charges


  !> This routine integrates the dynamical variable "n"
  !! \param charges
  subroutine prg_xlbo_nint(charges,n,n_0,n_1,n_2,n_3,n_4,n_5,mdstep,xl,dt,dt_base_in)
    implicit none
    real(dp), allocatable, intent(inout) :: n(:), n_0(:), n_1(:), n_2(:), n_3(:), n_4(:), n_5(:)
    real(dp), allocatable, intent(in) :: charges(:)
    type(xlbo_type), intent(inout) :: xl
    integer, intent(in) :: mdstep
    real(dp), intent(in), optional :: dt
    real(dp), intent(in), optional :: dt_base_in
    integer :: nats
    real(dp) :: kappa_use, alpha_use
    real(dp) :: c0_use, c1_use, c2_use, c3_use, c4_use, c5_use
    real(dp), save :: dt_base = -1.0_dp
    logical :: use_adaptive

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
      dt_base = -1.0_dp  ! Reset base timestep for new simulation
    endif

    ! Set base timestep from optional argument or first dt call
    if (present(dt_base_in)) then
      dt_base = dt_base_in
    else if (present(dt) .and. dt_base < 0.0_dp) then
      dt_base = dt
    endif

    ! Scale kappa with (dt/dt_base)^2 if dt provided, otherwise use fixed kappa
    if (present(dt) .and. dt_base > 0.0_dp) then
      kappa_use = kappa * (dt / dt_base)**2
    else
      kappa_use = kappa
    endif

    ! Determine if we should use adaptive coefficients
    use_adaptive = present(dt) .and. xl%nsteps_taken >= 5

    if (use_adaptive) then
      ! Get adaptive coefficients for the timestep pattern
      call get_xlbo_coefficients_adaptive(xl%dt_history, dt_base, &
                                           c0_use, c1_use, c2_use, c3_use, c4_use, c5_use, &
                                           alpha_use, alpha)

      ! Integration using adaptive coefficients (applied directly to history)
      n = 2.0_dp*n_0 - n_1 + xl%cc*kappa_use*(charges-n) &
           + alpha_use*(c0_use*n_0 + c1_use*n_1 + c2_use*n_2 + c3_use*n_3 + c4_use*n_4 + c5_use*n_5)
    else
      ! Integration using uniform coefficients (standard behavior)
      n = 2.0_dp*n_0 - n_1 + xl%cc*kappa_use*(charges-n) &
           + alpha*(C0_K5*n_0 + C1_K5*n_1 + C2_K5*n_2 + C3_K5*n_3 + C4_K5*n_4 + C5_K5*n_5)
    endif

    n_5 = n_4; n_4 = n_3; n_3 = n_2; n_2 = n_1; n_1 = n_0; n_0 = n

    ! Update timestep history if dt provided
    if (present(dt)) then
      xl%dt_history(5) = xl%dt_history(4)
      xl%dt_history(4) = xl%dt_history(3)
      xl%dt_history(3) = xl%dt_history(2)
      xl%dt_history(2) = xl%dt_history(1)
      xl%dt_history(1) = dt
      xl%nsteps_taken = xl%nsteps_taken + 1
    endif

  end subroutine prg_xlbo_nint


  !> This routine integrates the dynamical variable "n"
  !! \param charges
  subroutine prg_xlbo_nint_kernel(charges,n,n_0,n_1,n_2,n_3,n_4,n_5,mdstep,kernel,xl,dt,dt_base_in)
    implicit none
    real(dp), allocatable, intent(inout) :: n(:), n_0(:), n_1(:), n_2(:), n_3(:), n_4(:), n_5(:)
    real(dp), allocatable, intent(in) :: charges(:)
    real(dp), allocatable, intent(in) :: kernel(:,:)
    type(xlbo_type), intent(inout) :: xl
    integer, intent(in) :: mdstep
    real(dp), intent(in), optional :: dt
    real(dp), intent(in), optional :: dt_base_in
    integer :: nats
    real(dp) :: kappa_use, alpha_use
    real(dp) :: c0_use, c1_use, c2_use, c3_use, c4_use, c5_use
    real(dp), save :: dt_base = -1.0_dp
    logical :: use_adaptive

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
      dt_base = -1.0_dp  ! Reset base timestep for new simulation
    endif

    ! Set base timestep from optional argument or first dt call
    if (present(dt_base_in)) then
      dt_base = dt_base_in
    else if (present(dt) .and. dt_base < 0.0_dp) then
      dt_base = dt
    endif

    ! Scale kappa with (dt/dt_base)^2 if dt provided, otherwise use fixed kappa
    if (present(dt) .and. dt_base > 0.0_dp) then
      kappa_use = kappa * (dt / dt_base)**2
    else
      kappa_use = kappa
    endif

    ! Determine if we should use adaptive coefficients
    use_adaptive = present(dt) .and. xl%nsteps_taken >= 5

    ! From developper's code
    !   dn2dt2 = -MATMUL(KK0,(q-n))
    !   n = 2*n_0 - n_1 + kappa*dn2dt2 +
    !   alpha*(C0*n_0+C1*n_1+C2*n_2+C3*n_3+C4*n_4+C5*n_5+C6*n_6)
    !   n_6 = n_5; n_5 = n_4; n_4 = n_3; n_3 = n_2; n_2 = n_1; n_1 = n_0; n_0 = n

    if (use_adaptive) then
      ! Get adaptive coefficients for the timestep pattern
      call get_xlbo_coefficients_adaptive(xl%dt_history, dt_base, &
                                           c0_use, c1_use, c2_use, c3_use, c4_use, c5_use, &
                                           alpha_use, alpha)

      ! Integration using adaptive coefficients (applied directly to history)
      n = 2.0_dp*n_0 - n_1 - 1.0_dp*kappa_use*matmul(kernel,(charges-n)) &
           + alpha_use*(c0_use*n_0 + c1_use*n_1 + c2_use*n_2 + c3_use*n_3 + c4_use*n_4 + c5_use*n_5)
    else
      ! Integration using uniform coefficients (standard behavior)
      !call bml_print_matrix("ker",kernel,1,10,1,10)
      !write(*,*)matmul(kernel,(charges-n))
      !n = 2.0_dp*n_0 - n_1 + xl%cc*kappa*(charges-n) &
      n = 2.0_dp*n_0 - n_1 - 1.0_dp*kappa_use*matmul(kernel,(charges-n)) &
           + alpha*(C0_K5*n_0 + C1_K5*n_1 + C2_K5*n_2 + C3_K5*n_3 + C4_K5*n_4 + C5_K5*n_5)
    endif

    n_5 = n_4; n_4 = n_3; n_3 = n_2; n_2 = n_1; n_1 = n_0; n_0 = n

    ! Update timestep history if dt provided
    if (present(dt)) then
      xl%dt_history(5) = xl%dt_history(4)
      xl%dt_history(4) = xl%dt_history(3)
      xl%dt_history(3) = xl%dt_history(2)
      xl%dt_history(2) = xl%dt_history(1)
      xl%dt_history(1) = dt
      xl%nsteps_taken = xl%nsteps_taken + 1
    endif

  end subroutine prg_xlbo_nint_kernel


  !> This routine integrates the dynamical variable "n"
  !! \brief In this case we are passing a premultiplied ressidue x kernel
  !! tis is done to avoid rank-specific multiplication within this routine.
  !! \param charges
  subroutine prg_xlbo_nint_kernelTimesRes(charges,n,n_0,n_1,n_2,n_3,n_4,n_5,mdstep,kernelTimesRes,xl,dt,dt_base_in)
    implicit none
    real(dp), allocatable, intent(inout) :: n(:), n_0(:), n_1(:), n_2(:), n_3(:), n_4(:), n_5(:)
    real(dp), allocatable, intent(in) :: charges(:)
    real(dp), allocatable, intent(in) :: kernelTimesRes(:)
    type(xlbo_type), intent(inout) :: xl
    integer, intent(in) :: mdstep
    real(dp), intent(in), optional :: dt
    real(dp), intent(in), optional :: dt_base_in
    integer :: nats
    real(dp) :: kappa_use, alpha_use
    real(dp) :: c0_use, c1_use, c2_use, c3_use, c4_use, c5_use
    real(dp), save :: dt_base = -1.0_dp
    logical :: use_adaptive

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
      dt_base = -1.0_dp  ! Reset base timestep for new simulation
    endif

    ! Set base timestep from optional argument or first dt call
    if (present(dt_base_in)) then
      dt_base = dt_base_in
    else if (present(dt) .and. dt_base < 0.0_dp) then
      dt_base = dt
    endif

    ! Scale kappa with (dt/dt_base)^2 if dt provided, otherwise use fixed kappa
    if (present(dt) .and. dt_base > 0.0_dp) then
      kappa_use = kappa * (dt / dt_base)**2
    else
      kappa_use = kappa
    endif

    ! Determine if we should use adaptive coefficients
    use_adaptive = present(dt) .and. xl%nsteps_taken >= 5

    if (use_adaptive) then
      ! Get adaptive coefficients for the timestep pattern
      call get_xlbo_coefficients_adaptive(xl%dt_history, dt_base, &
                                           c0_use, c1_use, c2_use, c3_use, c4_use, c5_use, &
                                           alpha_use, alpha)

      ! Integration using adaptive coefficients (applied directly to history)
      n = 2.0_dp*n_0 - n_1 - 1.0_dp*kappa_use*kernelTimesRes &
           & + alpha_use*(c0_use*n_0 + c1_use*n_1 + c2_use*n_2 + c3_use*n_3 + c4_use*n_4 + c5_use*n_5)
    else
      ! Integration using uniform coefficients (standard behavior)
      n = 2.0_dp*n_0 - n_1 - 1.0_dp*kappa_use*kernelTimesRes &
           & + alpha*(C0_K5*n_0 + C1_K5*n_1 + C2_K5*n_2 + C3_K5*n_3 + C4_K5*n_4 + C5_K5*n_5)
    endif

    n_5 = n_4; n_4 = n_3; n_3 = n_2; n_2 = n_1; n_1 = n_0; n_0 = n

    ! Update timestep history if dt provided
    if (present(dt)) then
      xl%dt_history(5) = xl%dt_history(4)
      xl%dt_history(4) = xl%dt_history(3)
      xl%dt_history(3) = xl%dt_history(2)
      xl%dt_history(2) = xl%dt_history(1)
      xl%dt_history(1) = dt
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
