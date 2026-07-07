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

  !> Coefficients for K=5 (6-point history) XLBO dissipation
  real(dp), parameter :: C0 = -6.0_dp
  real(dp), parameter :: C1 = 14.0_dp
  real(dp), parameter :: C2 = -8.0_dp
  real(dp), parameter :: C3 = -3.0_dp
  real(dp), parameter :: C4 = 4.0_dp
  real(dp), parameter :: C5 = -1.0_dp
  real(dp), parameter :: kappa = 1.82_dp
  real(dp), parameter :: alpha = 0.018_dp

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


  !> Interpolate charges from non-uniform to uniform time grid using Lagrange interpolation (11-point version for K=10)
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
    real(dp) :: L(0:10,0:10)  ! Lagrange basis functions: L(k,i) = L_i(t_uniform(k))
    real(dp) :: dt_uniform
    integer :: i, j, k
    real(dp) :: numer, denom

    ! Build non-uniform source grid (where charges are stored)
    t(0) = 0.0_dp
    t(1) = -dt_history(1)
    do i = 2, 10
      t(i) = t(i-1) - dt_history(i)
    enddo

    ! Build uniform target grid using most recent timestep
    dt_uniform = dt_history(1)
    do i = 0, 10
      t_uniform(i) = -i * dt_uniform
    enddo

    ! Compute Lagrange basis functions L_i(t_uniform(k)) for all i,k
    ! L_i(t) = product over j≠i of [(t - t(j)) / (t(i) - t(j))]
    do k = 0, 10  ! For each target time
      do i = 0, 10  ! For each basis function
        L(k,i) = 1.0_dp
        do j = 0, 10
          if (j /= i) then
            numer = t_uniform(k) - t(j)
            denom = t(i) - t(j)
            L(k,i) = L(k,i) * (numer / denom)
          endif
        enddo
      enddo
    enddo

    ! Interpolate charges for each atom using the precomputed basis
    ni_0 = L(0,0)*n_0 + L(0,1)*n_1 + L(0,2)*n_2 + L(0,3)*n_3 + L(0,4)*n_4 + L(0,5)*n_5 + &
           L(0,6)*n_6 + L(0,7)*n_7 + L(0,8)*n_8 + L(0,9)*n_9 + L(0,10)*n_10
    ni_1 = L(1,0)*n_0 + L(1,1)*n_1 + L(1,2)*n_2 + L(1,3)*n_3 + L(1,4)*n_4 + L(1,5)*n_5 + &
           L(1,6)*n_6 + L(1,7)*n_7 + L(1,8)*n_8 + L(1,9)*n_9 + L(1,10)*n_10
    ni_2 = L(2,0)*n_0 + L(2,1)*n_1 + L(2,2)*n_2 + L(2,3)*n_3 + L(2,4)*n_4 + L(2,5)*n_5 + &
           L(2,6)*n_6 + L(2,7)*n_7 + L(2,8)*n_8 + L(2,9)*n_9 + L(2,10)*n_10
    ni_3 = L(3,0)*n_0 + L(3,1)*n_1 + L(3,2)*n_2 + L(3,3)*n_3 + L(3,4)*n_4 + L(3,5)*n_5 + &
           L(3,6)*n_6 + L(3,7)*n_7 + L(3,8)*n_8 + L(3,9)*n_9 + L(3,10)*n_10
    ni_4 = L(4,0)*n_0 + L(4,1)*n_1 + L(4,2)*n_2 + L(4,3)*n_3 + L(4,4)*n_4 + L(4,5)*n_5 + &
           L(4,6)*n_6 + L(4,7)*n_7 + L(4,8)*n_8 + L(4,9)*n_9 + L(4,10)*n_10
    ni_5 = L(5,0)*n_0 + L(5,1)*n_1 + L(5,2)*n_2 + L(5,3)*n_3 + L(5,4)*n_4 + L(5,5)*n_5 + &
           L(5,6)*n_6 + L(5,7)*n_7 + L(5,8)*n_8 + L(5,9)*n_9 + L(5,10)*n_10
    ni_6 = L(6,0)*n_0 + L(6,1)*n_1 + L(6,2)*n_2 + L(6,3)*n_3 + L(6,4)*n_4 + L(6,5)*n_5 + &
           L(6,6)*n_6 + L(6,7)*n_7 + L(6,8)*n_8 + L(6,9)*n_9 + L(6,10)*n_10
    ni_7 = L(7,0)*n_0 + L(7,1)*n_1 + L(7,2)*n_2 + L(7,3)*n_3 + L(7,4)*n_4 + L(7,5)*n_5 + &
           L(7,6)*n_6 + L(7,7)*n_7 + L(7,8)*n_8 + L(7,9)*n_9 + L(7,10)*n_10
    ni_8 = L(8,0)*n_0 + L(8,1)*n_1 + L(8,2)*n_2 + L(8,3)*n_3 + L(8,4)*n_4 + L(8,5)*n_5 + &
           L(8,6)*n_6 + L(8,7)*n_7 + L(8,8)*n_8 + L(8,9)*n_9 + L(8,10)*n_10
    ni_9 = L(9,0)*n_0 + L(9,1)*n_1 + L(9,2)*n_2 + L(9,3)*n_3 + L(9,4)*n_4 + L(9,5)*n_5 + &
           L(9,6)*n_6 + L(9,7)*n_7 + L(9,8)*n_8 + L(9,9)*n_9 + L(9,10)*n_10
    ni_10 = L(10,0)*n_0 + L(10,1)*n_1 + L(10,2)*n_2 + L(10,3)*n_3 + L(10,4)*n_4 + L(10,5)*n_5 + &
            L(10,6)*n_6 + L(10,7)*n_7 + L(10,8)*n_8 + L(10,9)*n_9 + L(10,10)*n_10

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
    ! Scale kappa with dt^2 when dt is present, regardless of interpolation
    ! The kappa term couples to the current SCF charges and must scale with actual timestep
    if (use_K10) then
      if (present(dt)) then
        kappa_use = kappa_K10 * dt**2
      else
        kappa_use = kappa_K10
      endif
      alpha_use = alpha_K10
    else
      if (present(dt)) then
        kappa_use = kappa * dt**2
      else
        kappa_use = kappa
      endif
      alpha_use = alpha
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

        ! Integration using interpolated charges (K=10)
        n = 2.0_dp*ni_0 - ni_1 + xl%cc*kappa_use*(charges-n) &
             + alpha_use*(C0_K10*ni_0+C1_K10*ni_1+C2_K10*ni_2+C3_K10*ni_3+C4_K10*ni_4+C5_K10*ni_5 &
                        +C6_K10*ni_6+C7_K10*ni_7+C8_K10*ni_8+C9_K10*ni_9+C10_K10*ni_10)

        deallocate(ni_0, ni_1, ni_2, ni_3, ni_4, ni_5)
        deallocate(ni_6, ni_7, ni_8, ni_9, ni_10)
      else
        ! Allocate interpolated charge arrays for K=5
        allocate(ni_0(nats), ni_1(nats), ni_2(nats), ni_3(nats), ni_4(nats), ni_5(nats))

        ! Interpolate historical charges to uniform grid (6 points)
        call prg_xlbo_interpolate_charges(xl%dt_history, n_0, n_1, n_2, n_3, n_4, n_5, &
                                           ni_0, ni_1, ni_2, ni_3, ni_4, ni_5, nats)

        ! Integration using interpolated charges (K=5)
        n = 2.0_dp*ni_0 - ni_1 + xl%cc*kappa_use*(charges-n) &
             + alpha_use*(C0*ni_0+C1*ni_1+C2*ni_2+C3*ni_3+C4*ni_4+C5*ni_5)

        deallocate(ni_0, ni_1, ni_2, ni_3, ni_4, ni_5)
      endif
    else
      ! Integration using raw charges (standard behavior)
      if (use_K10) then
        n = 2.0_dp*n_0 - n_1 + xl%cc*kappa_use*(charges-n) &
             + alpha_use*(C0_K10*n_0+C1_K10*n_1+C2_K10*n_2+C3_K10*n_3+C4_K10*n_4+C5_K10*n_5 &
                        +C6_K10*n_6+C7_K10*n_7+C8_K10*n_8+C9_K10*n_9+C10_K10*n_10)
      else
        n = 2.0_dp*n_0 - n_1 + xl%cc*kappa_use*(charges-n) &
             + alpha_use*(C0*n_0+C1*n_1+C2*n_2+C3*n_3+C4*n_4+C5*n_5)
      endif
    endif

    ! Shift history arrays
    if (use_K10) then
      n_10 = n_9; n_9 = n_8; n_8 = n_7; n_7 = n_6; n_6 = n_5
    endif
    n_5 = n_4; n_4 = n_3; n_3 = n_2; n_2 = n_1; n_1 = n_0; n_0 = n

    ! Update timestep history if dt provided (dt is ratio: 1.0 or 0.5)
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
      xl%dt_history(1) = dt  ! Store ratio, not absolute timestep
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
    ! Scale kappa with dt^2 when dt is present, regardless of interpolation
    ! The kappa term couples to the current SCF charges and must scale with actual timestep
    if (use_K10) then
      if (present(dt)) then
        kappa_use = kappa_K10 * dt**2
      else
        kappa_use = kappa_K10
      endif
      alpha_use = alpha_K10
    else
      if (present(dt)) then
        kappa_use = kappa * dt**2
      else
        kappa_use = kappa
      endif
      alpha_use = alpha
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

        ! Integration using interpolated charges (K=10)
        n = 2.0_dp*ni_0 - ni_1 - 1.0_dp*kappa_use*matmul(kernel,(charges-n)) &
             + alpha_use*(C0_K10*ni_0+C1_K10*ni_1+C2_K10*ni_2+C3_K10*ni_3+C4_K10*ni_4+C5_K10*ni_5 &
                        +C6_K10*ni_6+C7_K10*ni_7+C8_K10*ni_8+C9_K10*ni_9+C10_K10*ni_10)

        deallocate(ni_0, ni_1, ni_2, ni_3, ni_4, ni_5)
        deallocate(ni_6, ni_7, ni_8, ni_9, ni_10)
      else
        ! Allocate interpolated charge arrays for K=5
        allocate(ni_0(nats), ni_1(nats), ni_2(nats), ni_3(nats), ni_4(nats), ni_5(nats))

        ! Interpolate historical charges to uniform grid (6 points)
        call prg_xlbo_interpolate_charges(xl%dt_history, n_0, n_1, n_2, n_3, n_4, n_5, &
                                           ni_0, ni_1, ni_2, ni_3, ni_4, ni_5, nats)

        ! Integration using interpolated charges (K=5)
        n = 2.0_dp*ni_0 - ni_1 - 1.0_dp*kappa_use*matmul(kernel,(charges-n)) &
             + alpha_use*(C0*ni_0+C1*ni_1+C2*ni_2+C3*ni_3+C4*ni_4+C5*ni_5)

        deallocate(ni_0, ni_1, ni_2, ni_3, ni_4, ni_5)
      endif
    else
      ! Integration using raw charges (standard behavior)
      if (use_K10) then
        n = 2.0_dp*n_0 - n_1 - 1.0_dp*kappa_use*matmul(kernel,(charges-n)) &
             + alpha_use*(C0_K10*n_0+C1_K10*n_1+C2_K10*n_2+C3_K10*n_3+C4_K10*n_4+C5_K10*n_5 &
                        +C6_K10*n_6+C7_K10*n_7+C8_K10*n_8+C9_K10*n_9+C10_K10*n_10)
      else
        n = 2.0_dp*n_0 - n_1 - 1.0_dp*kappa_use*matmul(kernel,(charges-n)) &
             + alpha_use*(C0*n_0+C1*n_1+C2*n_2+C3*n_3+C4*n_4+C5*n_5)
      endif
    endif

    ! Shift history arrays
    if (use_K10) then
      n_10 = n_9; n_9 = n_8; n_8 = n_7; n_7 = n_6; n_6 = n_5
    endif
    n_5 = n_4; n_4 = n_3; n_3 = n_2; n_2 = n_1; n_1 = n_0; n_0 = n

    ! Update timestep history if dt provided (dt is ratio: 1.0 or 0.5)
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
      xl%dt_history(1) = dt  ! Store ratio, not absolute timestep
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
    ! Scale kappa with dt^2 when dt is present, regardless of interpolation
    ! The kappa term couples to the current SCF charges and must scale with actual timestep
    if (use_K10) then
      if (present(dt)) then
        kappa_use = kappa_K10 * dt**2
      else
        kappa_use = kappa_K10
      endif
      alpha_use = alpha_K10
    else
      if (present(dt)) then
        kappa_use = kappa * dt**2
      else
        kappa_use = kappa
      endif
      alpha_use = alpha
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

        ! Integration using interpolated charges (K=10)
        n = 2.0_dp*ni_0 - ni_1 - 1.0_dp*kappa_use*kernelTimesRes &
             & + alpha_use*(C0_K10*ni_0+C1_K10*ni_1+C2_K10*ni_2+C3_K10*ni_3+C4_K10*ni_4+C5_K10*ni_5 &
                          +C6_K10*ni_6+C7_K10*ni_7+C8_K10*ni_8+C9_K10*ni_9+C10_K10*ni_10)

        deallocate(ni_0, ni_1, ni_2, ni_3, ni_4, ni_5)
        deallocate(ni_6, ni_7, ni_8, ni_9, ni_10)
      else
        ! Allocate interpolated charge arrays for K=5
        allocate(ni_0(nats), ni_1(nats), ni_2(nats), ni_3(nats), ni_4(nats), ni_5(nats))

        ! Interpolate historical charges to uniform grid (6 points)
        call prg_xlbo_interpolate_charges(xl%dt_history, n_0, n_1, n_2, n_3, n_4, n_5, &
                                           ni_0, ni_1, ni_2, ni_3, ni_4, ni_5, nats)

        ! Integration using interpolated charges (K=5)
        n = 2.0_dp*ni_0 - ni_1 - 1.0_dp*kappa_use*kernelTimesRes &
             & + alpha_use*(C0*ni_0+C1*ni_1+C2*ni_2+C3*ni_3+C4*ni_4+C5*ni_5)

        deallocate(ni_0, ni_1, ni_2, ni_3, ni_4, ni_5)
      endif
    else
      ! Integration using raw charges (standard behavior)
      if (use_K10) then
        n = 2.0_dp*n_0 - n_1 - 1.0_dp*kappa_use*kernelTimesRes &
             & + alpha_use*(C0_K10*n_0+C1_K10*n_1+C2_K10*n_2+C3_K10*n_3+C4_K10*n_4+C5_K10*n_5 &
                          +C6_K10*n_6+C7_K10*n_7+C8_K10*n_8+C9_K10*n_9+C10_K10*n_10)
      else
        n = 2.0_dp*n_0 - n_1 - 1.0_dp*kappa_use*kernelTimesRes &
             & + alpha_use*(C0*n_0+C1*n_1+C2*n_2+C3*n_3+C4*n_4+C5*n_5)
      endif
    endif

    ! Shift history arrays
    if (use_K10) then
      n_10 = n_9; n_9 = n_8; n_8 = n_7; n_7 = n_6; n_6 = n_5
    endif
    n_5 = n_4; n_4 = n_3; n_3 = n_2; n_2 = n_1; n_1 = n_0; n_0 = n

    ! Update timestep history if dt provided (dt is ratio: 1.0 or 0.5)
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
      xl%dt_history(1) = dt  ! Store ratio, not absolute timestep
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
