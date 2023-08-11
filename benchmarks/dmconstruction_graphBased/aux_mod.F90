!> Auxiliary module module.
!! \brief This module contains auxiliary routines.
!!
module aux_mod

  use bml
  use prg_kernelparser_mod

  implicit none

  private  !Everything is private by default

  integer, parameter :: dp = kind(1.0d0)

  !> General hamiltonian/system parameters type
  type, public :: bioham_type
    integer :: replicatex
    integer :: replicatey
    integer :: replicatez
    character(100) :: jobname
    character(100) :: bml_type
    character(50) :: systemfilename
    character(100) :: parampath
    real(dp) :: threshold
    logical :: reshuffle
    integer :: mdim
    integer :: verbose
  end type bioham_type

  !> Graph partitioning input parameters
  type, public :: gppar_type
    character(100) :: jobname
    integer :: numparts
    real(dp) :: threshold
  end type gppar_type

  public :: prg_parse_bioham, prg_parse_gppar

contains

  !> Model Ham parse.
  subroutine prg_parse_bioham(bioham,filename)

    implicit none
    type(bioham_type), intent(inout) :: bioham
    integer, parameter :: nkey_char = 4, nkey_int = 6, nkey_re = 2, nkey_log = 1
    character(len=*) :: filename

    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=50) :: &
         'JobName=', 'BMLType=', 'SystemFileName=', 'ParamPath=' ]
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
         'GetModelBioHam', 'Dense', 'prot.pdb', './latteTBparams' ]

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
         'ReplicateX=','ReplicateY=','ReplicateZ=', 'Mdim', 'Verbosity=', 'Seed=']
    integer :: valvector_int(nkey_int) = (/ &
         1, 1, 1, 0, 1, 100   /)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
         'Threshold=', 'Aux=']
    real(dp) :: valvector_re(nkey_re) = (/&
         0.0, 0.0 /)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=50) :: &
         'Reshuffle=']
    logical :: valvector_log(nkey_log) = (/&
         .false./)

    !Start and stop characters
    character(len=50), parameter :: startstop(2) = [character(len=50) :: &
         'BIOHAM{', '}']

    call prg_parsing_kernel(keyvector_char,valvector_char&
         ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
         keyvector_log,valvector_log,trim(filename),startstop)

    !Characters
    bioham%JobName = valvector_char(1)

    if(valvector_char(2) == "Dense")then
      bioham%bml_type = BML_MATRIX_DENSE
    elseif(valvector_char(2) == "Ellpack")then
      bioham%bml_type = BML_MATRIX_ELLPACK
    elseif(valvector_char(2) == "CSR")then
      bioham%bml_type = BML_MATRIX_CSR
    elseif(valvector_char(2) == "Ellblock")then
      bioham%bml_type = BML_MATRIX_ELLBLOCK
    endif

    bioham%systemfilename = valvector_char(3)
    bioham%parampath = valvector_char(4)

    !Integers
    bioham%replicatex = valvector_int(1)
    bioham%replicatey = valvector_int(2)
    bioham%replicatez = valvector_int(3)
    bioham%mdim = valvector_int(4)
    bioham%verbose = valvector_int(5)

    !Reals
    bioham%threshold = valvector_re(1)

    !Logicals
    bioham%reshuffle = valvector_log(1)

  end subroutine prg_parse_bioham

  !> Graph partitioning parameters input parser
  subroutine prg_parse_gppar(gppar,filename)

    implicit none
    type(gppar_type), intent(inout) :: gppar
    integer, parameter :: nkey_char = 1, nkey_int = 1, nkey_re = 2, nkey_log = 1
    character(len=*) :: filename

    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=50) :: &
         'JobName=' ]
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
         'DoGraphPart' ]

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
         'NumberOfParts=']
    integer :: valvector_int(nkey_int) = (/ &
         2  /)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
         'Threshold=', 'Aux=']
    real(dp) :: valvector_re(nkey_re) = (/&
         0.0, 0.0 /)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=50) :: &
         'algo=']
    logical :: valvector_log(nkey_log) = (/&
         .false./)

    !Start and stop characters
    character(len=50), parameter :: startstop(2) = [character(len=50) :: &
         'GPART{', '}']

    call prg_parsing_kernel(keyvector_char,valvector_char&
         ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
         keyvector_log,valvector_log,trim(filename),startstop)

    !Characters
    gppar%JobName = valvector_char(1)

    !Reals
    gppar%threshold = valvector_re(1)

    !Integers
    gppar%numparts = valvector_int(1)

  end subroutine prg_parse_gppar

end module aux_mod
