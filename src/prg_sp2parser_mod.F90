!> SP2 parser.
!! \ingroup PROGRESS
!! This module is used to parse all the input variables for the SP2 method
!! electronic structure solver.
!! Adding a new input keyword to the parser:
!! - If the variable is real, we have to increase nkey_re.
!! - Add the keyword (character type) in the keyvector_re vector.
!! - Add a default value (real type) in the valvector_re.
!! - Define a new variable and pass the value through valvector_re(num)
!! where num is the position of the new keyword in the vector.
!!
module prg_sp2parser_mod

  use prg_openfiles_mod
  use prg_kernelparser_mod
  use bml

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  !> General SP2 solver type
  !!
  type, public :: sp2data_type
     character(20) :: jobname
     integer :: verbose
     integer :: minsp2iter
     integer :: maxsp2iter
     real(dp) :: sp2tol
     real(dp) :: threshold
     real(dp) :: bndfil
     integer :: mdim
     integer :: ndim
     character :: sdim(3)
     real(dp) :: pdim(3)
     character(20) :: bml_type
     character(10) :: sp2conv
     character(10) :: flavor
  end type sp2data_type

  public :: prg_parse_sp2

contains

  !> The parser for SP2 solver.
  !!
  subroutine prg_parse_sp2(sp2data,filename)

    implicit none
    type(sp2data_type), intent(inout) :: sp2data
    integer, parameter :: nkey_char = 4, nkey_int = 6, nkey_re = 3, nkey_log = 1
    character(len=*) :: filename

    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=100) :: &
      'JobName=', 'BMLType=','SP2Conv=','Flavor=' ]
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
      'MyJob'   , 'Dense'   ,'REL', 'Alg2' ]

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
    'MDim=', 'VarInt=', 'MinSP2Iter=', 'MaxSP2Iter=','Ndim=','Verbose=']
    integer :: valvector_int(nkey_int) = (/ &
       -1   ,     0    ,      10       ,      100 , 1, 0 /)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
      'NumThresh=','SP2Tol=','BndFil=' ]
    real(dp) :: valvector_re(nkey_re) = (/&
         0.0      ,   0.00000001    ,0.0 /)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=100) :: &
      'DUMMY=']
    logical :: valvector_log(nkey_log) = (/&
     .false./)

    !Start and stop characters
    character(len=50), parameter :: startstop(2) = [character(len=50) :: &
      'SP2{', '}']

    call prg_parsing_kernel(keyvector_char,valvector_char&
    ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
    keyvector_log,valvector_log,trim(filename),startstop)

    !Characters
    sp2data%JobName = valvector_char(1)

    if(valvector_char(2) == "Dense")then
       sp2data%bml_type = BML_MATRIX_DENSE
    elseif(valvector_char(2) == "Ellpack")then
       sp2data%bml_type = BML_MATRIX_ELLPACK
    endif
    sp2data%sp2conv = valvector_char(3)
    sp2data%flavor = valvector_char(4)

    !Reals
    sp2data%threshold = valvector_re(1)
    sp2data%sp2tol = valvector_re(2)
    sp2data%BndFil = valvector_re(3)

    !Logicals

    !Integers
    sp2data%mdim = valvector_int(1)
    sp2data%minsp2iter = valvector_int(3)
    sp2data%maxsp2iter = valvector_int(4)
    sp2data%ndim = valvector_int(5)
    sp2data%verbose = valvector_int(6)

  end subroutine prg_parse_sp2

end module prg_sp2parser_mod
