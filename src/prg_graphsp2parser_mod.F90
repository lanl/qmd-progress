!> Graph partitioning SP2 parser.
!! \ingroup PROGRESS
!! \brief This module is used to parse all the neccesary input variables for
!! graph-based SP2 electronic structure solver.
!! Adding a new input keyword to the parser:
!! - If the variable is real, we have to increase nkey_re.
!! - Add the keyword (character type) in the keyvector_re vector.
!! - Add a default value (real type) in the valvector_re.
!! - Define a new variable and pass the value through valvector_re(num)
!! where num is the position of the new keyword in the vector.
!!
module prg_graphsp2parser_mod

  use prg_openfiles_mod
  use prg_kernelparser_mod
  use bml

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  !> General SP2 solver type
  !!
  type, public :: gsp2data_type
    character(20) :: jobname
    character(50) :: hamfile
    integer :: verbose
    integer :: minsp2iter
    integer :: maxsp2iter
    integer :: nodesPerPart
    integer :: natoms
    integer :: partition_count
    
    !> For partitioning in x,y,z direction
    !! This is only used in case of "domain/spatial" type
    !! of partitioning
    integer :: nx, ny, nz

    real(dp) :: sp2tol
    real(dp) :: threshold
    real(dp) :: bndfil
    real(dp) :: gthreshold
    real(dp) :: errlimit
    integer :: mdim
    integer :: ndim
    character :: sdim(3)
    real(dp) :: pdim(3)
    character(20) :: bml_type
    character(10) :: sp2conv
    character(10) :: graph_element
    character(10) :: partition_type
    character(10) :: partition_refinement
    logical :: double_jump
    real(dp) :: covgfact !Factor for tuning the extension of the covalency
    real(dp) :: nlgcut   !Radius cutoff for the hmiltonian (distance) based graph
    integer :: parteach !Do the partition each PartEach mdsteps
  end type gsp2data_type

  public :: prg_parse_gsp2

contains

  !> The parser for SP2 solver.
  !!
  subroutine prg_parse_gsp2(gsp2data,filename)

    implicit none
    type(gsp2data_type), intent(inout) :: gsp2data
    integer, parameter :: nkey_char = 7, nkey_int = 11, nkey_re = 7, nkey_log = 2
    character(len=*) :: filename

    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=100) :: &
         'JobName=', 'BMLType=','SP2Conv=', 'HamFile=', 'GraphElement=', &
         'PartitionType=', 'PartitionRefinement=']
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
         'MyJob', 'Dense','REL', 'text.mtx', 'Atom', 'Block', 'None']

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
         'Mdim=', 'MinSP2Iter=', 'MaxSP2Iter=','Ndim=', 'NodesPerPart=', 'NAtoms=', &
         'PartitionCount=', 'PartEach=', 'PartitionCountX=','PartitionCountY=',&
         &'PartitionCountZ=']
    integer :: valvector_int(nkey_int) = (/ &
         -1, 10, 100, 1, 16, 1, 1,1,0,0,0 /)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
         'MatrixThreshold=','SP2Tol=','BndFil=', 'GraphThreshold=', 'ErrLimit=', 'CovGraphFact=', 'NLGraphCut=' ]
    real(dp) :: valvector_re(nkey_re) = (/&
         0.00001,     0.00000001, 0.0,  0.00000000001,    0.0, 2.5, 2.5 /)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=100) :: &
         'DoubleJump=', 'Log2=']
    logical :: valvector_log(nkey_log) = (/&
         .true., .false./)

    !Start and stop characters
    character(len=50), parameter :: startstop(2) = [character(len=50) :: &
         'GSP2{', '}']

    call prg_parsing_kernel(keyvector_char,valvector_char&
         ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
         keyvector_log,valvector_log,trim(filename),startstop)

    !Characters
    gsp2data%JobName = valvector_char(1)

    if(valvector_char(2) == "Dense")then
      gsp2data%bml_type = BML_MATRIX_DENSE
    elseif(valvector_char(2) == "Ellpack")then
      gsp2data%bml_type = BML_MATRIX_ELLPACK
    elseif(valvector_char(2) == "Ellblock")then
      gsp2data%bml_type = BML_MATRIX_ELLBLOCK
    endif
    gsp2data%sp2conv = valvector_char(3)
    gsp2data%hamfile = valvector_char(4)
    gsp2data%graph_element = valvector_char(5)

    gsp2data%partition_type = valvector_char(6)
    gsp2data%partition_refinement = valvector_char(7)

    !Reals
    gsp2data%threshold = valvector_re(1)
    gsp2data%sp2tol = valvector_re(2)
    gsp2data%bndFil = valvector_re(3)
    gsp2data%gthreshold = valvector_re(4)
    gsp2data%errlimit = valvector_re(5)
    gsp2data%covgfact = valvector_re(6)
    gsp2data%nlgcut = valvector_re(7)

    !Logicals
    gsp2data%double_jump = valvector_log(1)

    !Integers
    gsp2data%mdim = valvector_int(1)
    gsp2data%minsp2iter = valvector_int(2)
    gsp2data%maxsp2iter = valvector_int(3)
    gsp2data%ndim = valvector_int(4)
    gsp2data%nodesPerPart = valvector_int(5)
    gsp2data%natoms = valvector_int(6)
    gsp2data%partition_count= valvector_int(7)
    gsp2data%parteach= valvector_int(8)
    gsp2data%nx = valvector_int(9)
    gsp2data%ny = valvector_int(10)
    gsp2data%nz = valvector_int(11)

  end subroutine prg_parse_gsp2

end module prg_graphsp2parser_mod
