!> LATTE parser.
!! \ingroup LATTE
!! This module is used to parse all the necessary input variables for a LATTE TB run (SCF/OPT/MD)
!! Adding a new input keyword to the parser:
!! - If the variable is real, we have to increase nkey_re.
!! - Add the keyword (character type) in the keyvector_re vector.
!! - Add a default value (real type) in the valvector_re.
!! - Define a new variable int he latte type and pass the value through valvector_re(num)
!! where num is the position of the new keyword in the vector.
!!
module latteparser_latte_mod

  use prg_openfiles_mod
  use prg_kernelparser_mod
  use bml

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  !> General latte input variables type.
  !!
  type, public :: latte_type

    !> Name of the current job.
    character(20) :: jobname

    !> Verbosity level.
    integer :: verbose

    !> Threshold values for matrix elements.
    real(dp) :: threshold

    !> Max nonzero elements per row for every row see \cite Mniszewski2015 .
    integer :: mdim

    !> Matrix format (Dense or Ellpack).
    character(20) :: bml_type

    !> Distribution mode (sequential, distributed, or graph_distributed).
    character(20) :: bml_dmode

    !> Coulomb Accuracy.
    real(dp) :: coul_acc

    !> Pulay mixing coefficient.
    real(dp) :: pulaycoeff

    !> Linear mixing coefficient.
    real(dp) :: mixcoeff

    !> Coulomb Accuracy.
    integer :: mpulay

    !> Maximum SCF iterations.
    integer :: maxscf

    !> SCF tolerance.
    real(dp) :: scftol

    !> Z Matrix calculation type.
    character(20) :: ZMat

    !> Solver method
    character(20) :: method

    !> Estimated ration between real & k space time efficiency.
    real(dp) :: timeratio

    !> Total number of steps for MD simulation.
    integer :: mdsteps

    !> Total number of steps for MD simulation.
    real(dp) :: timestep

    !> Total number of steps for MD simulation.
    character(100) :: parampath

    !> File containing coordinates.
    character(100) :: coordsfile

    !> File containing coordinates.
    integer :: nlisteach

    !> Restart calculation.
    logical :: restart

    !> Chemical potential initial guess of value.
    real(dp) :: efermi
    
    !> Electronic temperature kbT (in eV)
    real(dp) :: kbt

    !> Logical variable. If set to T efermi will be adjusted dynamically
    logical :: mumd

  end type latte_type

  public :: parse_latte

contains

  !> The parser for Latte General input variables.
  !!
  subroutine parse_latte(latte,filename)

    implicit none
    type(latte_type) :: latte
    integer, parameter :: nkey_char = 7, nkey_int = 6, nkey_re = 9, nkey_log = 2
    character(len=*) :: filename

    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=100) :: &
         'JobName=', 'BMLType=','ZMat=','Method=','ParamPath=','CoordsFile=', &
         'BMLDistributionType=']
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
         'MyJob'   , 'Dense'   ,'Diag','Diag','/home/name/','coords.dat', &
         'Sequential']

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
         'MDim=', 'Verbose=', 'MPulay=', 'MaxSCFIter=', 'MDSteps=', 'NlistEach=']
    integer :: valvector_int(nkey_int) = (/ &
         -1   ,     0    ,      5       ,  100, 100, 1 /)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
         'Threshold=','CoulAcc=','PulayCoeff=','SCFTol=','TimeRatio=','MixCoeff=','TimeStep=', &
          'EFermi=','kbt=']
    real(dp) :: valvector_re(nkey_re) = (/&
         0.00001    ,   0.00001    ,0.01    ,   0.001 ,  10.0, 0.5, 0.5, -1.0,  
         0.0 /)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=100) :: &
         'Restart=', 'MuMD=']
    logical :: valvector_log(nkey_log) = (/&
         .false., .false./)

    !Start and stop characters
    character(len=50), parameter :: startstop(2) = [character(len=50) :: &
         'Latte{', '}']

    call prg_parsing_kernel(keyvector_char,valvector_char&
         ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
         keyvector_log,valvector_log,trim(filename),startstop)

    !Characters
    latte%JobName = valvector_char(1)
    latte%ZMat = valvector_char(3)
    latte%method = valvector_char(4)
    latte%parampath = valvector_char(5)
    latte%coordsfile = valvector_char(6)

    if(valvector_char(2) == "Dense")then
      latte%bml_type = BML_MATRIX_DENSE
    elseif(valvector_char(2) == "Ellpack")then
      latte%bml_type = BML_MATRIX_ELLPACK
    endif

    if(valvector_char(7) == "Distributed")then
      latte%bml_dmode = BML_DMODE_DISTRIBUTED
    elseif(valvector_char(7) == "Sequential")then
      latte%bml_dmode = BML_DMODE_SEQUENTIAL
    endif

    !Reals
    latte%threshold = valvector_re(1)
    latte%coul_acc = valvector_re(2)
    latte%pulaycoeff = valvector_re(3)
    latte%scftol = valvector_re(4)
    latte%timeratio = valvector_re(5)
    latte%mixcoeff = valvector_re(6)
    latte%timestep = valvector_re(7)
    latte%efermi = valvector_re(8)
    latte%kbt = valvector_re(9)

    !Logicals
    latte%restart = valvector_log(1)
    latte%mumd = valvector_log(2)

    !Integers
    latte%mdim = valvector_int(1)
    latte%verbose = valvector_int(2)
    latte%mpulay = valvector_int(3)
    latte%maxscf = valvector_int(4)
    latte%mdsteps = valvector_int(5)
    latte%nlisteach = valvector_int(6)

  end subroutine parse_latte

end module latteparser_latte_mod
