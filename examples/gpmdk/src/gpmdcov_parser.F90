!> Module for parsing related functions.
!! \brief This module will be used to have parsing routines.
!!
!! @ingroup GPMDCOV
!!
!!
module gpmdcov_parser_mod

  use prg_kernelparser_mod

  integer, parameter :: dp = kind(1.0d0)

  !> gpmd type
  type, public :: gpmd_type

    !> Job name.
    character(len=100) :: job_name

    !> Output file name
    character(len=20) :: output_file_name

    !> Optimize with dumped MD
    logical :: dovelresc

    !> Initial temperature
    real(dp) :: temp0

    !> Velocity rescaling factor
    real(dp) :: velresc_fact

    !> SMD Force Constant for Start of Switch
    real(dp) :: smdforceconststart

    !> SMD Force Constant for End of Switch
    real(dp) :: smdforceconstend

    !> SMD optimal separation distance
    real(dp) :: smdr0

    !> SMD minimum separation distance
    real(dp) :: smdminsep

    !> SMD maximum separation distance
    real(dp) :: smdmaxsep

    !Output control 
    logical :: writetraj

    !Restart file control 
    logical :: writerestart

    !Write residue in trajectory
    logical :: writeResTraj

    !Restart from dump file
    logical :: restartfromdump

    !> Write trajectory coordinates each 
    integer :: writetreach

    !> Write trajectory coordinates each 
    integer :: writerestarteach

    !> Replicae system in x direction 
    integer :: replicatex
    
    !> Replicae system in y direction 
    integer :: replicatey

    !> Replicae system in z direction 
    integer :: replicatez

    !> Track the rectivity  
    logical :: trackreactivity

    !> Track parts from partition. If 0, no part is tracked
    integer :: tracknparts

    !> Track parts from partition
    integer, allocatable :: trackparts(:)

    !> Dump into a file each dumpeach steps
    integer :: dumpeach 

    !> Number of steered pairs for SMD
    integer :: smdnumpairs

    !> SMD atom 1 indicies
    integer, allocatable :: smdatomind1(:)

    !> SMD atom 2 indicies
    integer, allocatable :: smdatomind2(:) 

    !> Use LATTE code to compute Hamiltonians
    logical :: useLatte

    !> Use deuterium mass for H atoms
    logical :: htod

    !> Use Steered MD (SMD)
    logical :: usesmd

    !> Output trajectory format
    character(len=100) :: traj_format

    !> Use Langevin integration
    logical :: langevin

    !> Langevin method
    character(len=100) :: langevin_method

    !> Langevin decay constant
    real(dp) :: langevin_gamma

    !> Number of initial minimization steps before dynamics
    integer :: minimization_steps

    !> Compute currents option 
    logical :: compcurr

    !> Currents Threshold 
    real(dp) :: currthr

    !> Translate and fold to box
    logical :: trfl

    !> Use Vectorized SKBlock method in PROGRESS
    logical :: usevectsk
    
    !> Apply a voltage on selceted atomic sites
    logical :: applyv
    
    !> Voltage filename
    character(100) :: voltagef
    
  end type gpmd_type

  private

  public :: gpmdcov_parse

contains

  !> Gpmdcov parser.
  !! \brief This module is used to parse all the input variables for this program.
  !! Adding a new input keyword to the parser:
  !! - If the variable is real, we have to increase nkey_re.
  !! - Add the keyword (character type) in the keyvector_re vector.
  !! - Add a default value (real type) in the valvector_re.
  !! - Define a new variable and pass the value through valvector_re(num)
  !! where num is the position of the new keyword in the vector.
  !! \param filename File name for the input.
  !! \param gpmd type. 
  !!
  subroutine gpmdcov_parse(filename,gpmdt)
    implicit none 
    character(len=*), intent(in) :: filename
    type(gpmd_type), intent(inout) :: gpmdt
    integer, parameter :: nkey_char = 5, nkey_int = 9, nkey_re = 7, nkey_log = 13
    integer :: i
    real(dp) :: realtmp
    character(20) :: dummyc
    logical :: inGPMD
     
    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=50) :: &
         & 'JobName=', 'Var2C=', 'TrajectoryFormat=', 'LangevinMethod=', 'VoltageFile=']
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
         &'MyMol', 'Def2', 'PDB', 'Goga', 'None']

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
         & 'WriteCoordsEach=',"Var2I=","ReplicateX=","ReplicateY=","ReplicateZ=","PartsToTrack=",&
         & "DumpEach=","MinimizationSteps=","SMDNumPairs="]
    integer :: valvector_int(nkey_int) = (/ &
         & 1, 1, 0, 0, 0, 0, 0, 0, 0/)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
         & 'VRFactor=','InitialTemperature=','LangevinGamma=','SMDForceConstantStart=',&
         & 'SMDForceConstantEnd=','SMDR0=','CurrentThreshold=']
    real(dp) :: valvector_re(nkey_re) = (/&
         & 0.0_dp, 0.0_dp, 0.01_dp, 0.0_dp,0.2_dp,2.0_dp,0.1_dp/)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=50) :: &
         &'DoVelocityRescale=','WriteResidueInTrajectory=','WriteTrajectory=','TrackReactivity=',&
         &'RestartFromDump=','UseLATTE=','HtoD=','LangevinDynamics=','UseSMD=', &
         &'ComputeCurrents=', 'TranslateAndFoldToBox=', 'UseVectSKBlock=', 'ApplyVoltage=']
    logical :: valvector_log(nkey_log) = (/&
         &.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false., &
         &.false.,.True.,.false.,.false./)

    !Start and stop characters
    character(len=50), parameter :: startstop(2) = [character(len=50) :: &
         'GPMD{', '}']

    call prg_parsing_kernel(keyvector_char,valvector_char&
         ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
         keyvector_log,valvector_log,trim(filename),startstop)

    !Characters
    gpmdt%job_name = valvector_char(1)
    gpmdt%traj_format = valvector_char(3)
    gpmdt%langevin_method = valvector_char(4)
    
    !Integer
    gpmdt%writetreach = valvector_int(1)
    gpmdt%replicatex = valvector_int(3)
    gpmdt%replicatey = valvector_int(4)
    gpmdt%replicatez = valvector_int(5)
    gpmdt%tracknparts = valvector_int(6)
    inGPMD = .false.
    if(gpmdt%tracknparts > 0)then
      write(*,*)"Reading parts to track ..."
      open(1, file=trim(filename))
      do i = 1,10000
        read(1,*) dummyc
        write(*,*) trim(adjustl(dummyc))
        if(trim(adjustl(dummyc)) == "Parts[")then
          exit
        end if
        if(trim(adjustl(dummyc)) == "GPMD{") inGPMD = .true.
                
        if(trim(dummyc) == "}" .and. inGPMD)then
          write(*,*)'ERROR: No part numbers defined'
          write(*,*)"Here is an example block you should add to the"
          write(*,*)"input file"
          write(*,*)""
          write(*,*)"Parts["
          write(*,*)"   1"
          write(*,*)"   100"
          write(*,*)"]"
          stop
        end if
      end do
      allocate(gpmdt%trackparts(gpmdt%tracknparts))
      do i = 1,gpmdt%tracknparts
        read(1,*)gpmdt%trackparts(i)
      end do
      write(*,*)""
      close(1)
    endif

    gpmdt%dumpeach = valvector_int(7)
    gpmdt%minimization_steps = valvector_int(8)
    gpmdt%smdnumpairs = valvector_int(9)
    if(gpmdt%smdnumpairs > 0) then
      write(*,*) "Reading SMD pair atom indicies"
      open(1, file=trim(filename))
      do i = 1,10000
        read(1,*) dummyc
        !write(*,*) trim(adjustl(dummyc))
        if(trim(adjustl(dummyc)) == "SMDPairs[")then
          exit
        end if
        if(trim(adjustl(dummyc)) == "GPMD{") inGPMD = .true.
                
        if(trim(dummyc) == "}" .and. inGPMD)then
          write(*,*)'ERROR: No SMD Pairs defined'
          write(*,*)"Here is an example block you should add to the"
          write(*,*)"input file"
          write(*,*)""
          write(*,*)"SMDPairs["
          write(*,*)"   1 3"
          write(*,*)"   100 102"
          write(*,*)"]"
          stop
        end if
      end do
      allocate(gpmdt%smdatomind1(gpmdt%smdnumpairs))
      allocate(gpmdt%smdatomind2(gpmdt%smdnumpairs))
      do i = 1,gpmdt%smdnumpairs
        read(1,*) gpmdt%smdatomind1(i), gpmdt%smdatomind2(i)
      end do
      write(*,*)""
      close(1)
    endif
    do i = 1,gpmdt%smdnumpairs
      write(*,*) "SMD Pairs Atom 1 ",gpmdt%smdatomind1(i)," SMD Pairs Atom 2 ",gpmdt%smdatomind2(i)
    enddo

    !Reals
    gpmdt%velresc_fact = valvector_re(1)
    gpmdt%temp0 = valvector_re(2)
    gpmdt%langevin_gamma = valvector_re(3)
    gpmdt%smdforceconststart = valvector_re(4)
    gpmdt%smdforceconstend = valvector_re(5)
    gpmdt%smdr0 = valvector_re(6)
    gpmdt%currthr = valvector_re(7)

    !Logs
    gpmdt%dovelresc = valvector_log(1)
    gpmdt%writeResTraj = valvector_log(2)
    gpmdt%writeTraj = valvector_log(3)
    gpmdt%trackreactivity = valvector_log(4)
    gpmdt%restartfromdump = valvector_log(5)
    gpmdt%useLatte = valvector_log(6)
    gpmdt%htod = valvector_log(7)
    gpmdt%langevin = valvector_log(8)
    gpmdt%usesmd = valvector_log(9)
    gpmdt%compcurr = valvector_log(10)
    gpmdt%trfl = valvector_log(11)
    gpmdt%usevectsk = valvector_log(12)
    gpmdt%applyv = valvector_log(13)

    if(gpmdt%applyv)then 
        gpmdt%voltagef = valvector_char(5)
    endif

  end subroutine gpmdcov_parse

end module gpmdcov_parser_mod
