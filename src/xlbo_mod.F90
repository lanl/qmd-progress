!> A module to perform XLBO integration.
!! @ingroup PROGRESS 
!!
!! \brief This module will be used to compute integrate the dynamical variable "n" in xlbo. 
!!      
module xlbo_mod

  use openfiles_mod
  use bml   
  use kernelparser_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)
  
  !> Coefficients for modified Verlet integration
  real(dp), parameter :: C0 = -6.0_dp
  real(dp), parameter :: C1 = 14.0_dp
  real(dp), parameter :: C2 = -8.0_dp
  real(dp), parameter :: C3 = -3.0_dp
  real(dp), parameter :: C4 = 4.0_dp
  real(dp), parameter :: C5 = -1.0_dp; 

  !> Coefficients for modified Verlet integration
  real(dp), parameter :: kappa = 1.82_dp; 
  real(dp), parameter :: alpha = 0.018_dp;                         
  real(dp), parameter :: cc = 0.9_dp;        ! Scaled delta kernel

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

     !> Use SCF the first M_init MD steps
     integer :: minit

     !> Scaled delta Kernel
     real(dp) :: cc 

  end type xlbo_type 

  public :: parse_xlbo, xlbo_nint, xlbo_fcoulupdate

contains 

  !> The parser for XLBO parser.
  !!  
  subroutine parse_xlbo(xlbo,filename)

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
    'Verbose=','Minit=','MaxSCFIter=', 'MaxSCFInitIter=']                                   
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
     
    call parsing_kernel(keyvector_char,valvector_char&
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

  end subroutine parse_xlbo


  !> This routine integrates the dynamical variable "n"
  !! \param charges 
  subroutine xlbo_nint(charges,n,n_0,n_1,n_2,n_3,n_4,n_5,mdstep,xl)
  	implicit none 
  	real(dp), allocatable, intent(inout) :: n(:), n_0(:), n_1(:), n_2(:), n_3(:), n_4(:), n_5(:)
  	real(dp), allocatable, intent(in) :: charges(:)
    type(xlbo_type), intent(in) :: xl

    integer, intent(in) :: mdstep
    integer :: nats

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
    endif

    n = 2.0_dp*n_0 - n_1 + xl%cc*kappa*(charges-n) &
      + alpha*(C0*n_0+C1*n_1+C2*n_2+C3*n_3+C4*n_4+C5*n_5);
    n_5 = n_4; n_4 = n_3; n_3 = n_2; n_2 = n_1; n_1 = n_0; n_0 = n;

  end subroutine xlbo_nint


  !> Adjust forces for the linearized XLBOMD functional
  !! \param charges 
  subroutine xlbo_fcoulupdate(fcoul,charges,n)
    implicit none
    real(dp), intent(inout) :: fcoul(:,:),charges(:)
    real(dp), intent(inout) :: n(:)

        fcoul(1,:) = (2.0_dp*charges(:)-n(:))*fcoul(1,:)/n(:); 
        fcoul(2,:) = (2.0_dp*charges(:)-n(:))*fcoul(2,:)/n(:); 
        fcoul(3,:) = (2.0_dp*charges(:)-n(:))*fcoul(3,:)/n(:); 
  
  end subroutine xlbo_fcoulupdate

end module xlbo_mod

