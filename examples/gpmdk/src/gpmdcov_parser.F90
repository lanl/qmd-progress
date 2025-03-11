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

    !> XLBO ON/OFF option
    logical :: xlboON
   
    !> Coarse QMD 
    logical :: coarseqmd

    !> Do fine tol scf each x number of steps
    integer :: finetoleach

    !> Fine tol for scf 
    real(dp) :: finetol

    !> Coarse tol for scf 
    real(dp) :: coarsetol

    !> MD step to start profiling
    integer :: profile_start_step

    !> MD step to stop profiling
    integer :: profile_stop_step

    !> Use London dispersion forces
    logical :: disp

    !> Freeze atoms 
    logical :: freeze

    !> File containing the atoms to freeze
    character(100) :: freezef 
    
  end type gpmd_type

  !> electrontic structure output type
  !! Controls output of electronic structure information
  !! as collected across QMD simulation
  type, public :: estructout_type
          !> TDOS Output file name
          character (len=100) :: tdos_output_filename

          !> PDOS Output file name
          character (len=100) :: pdos_output_filename

          !> Write out Density of State
          logical :: write_tdos

          !> Compute Projected Density of State
          logical :: compute_pdos

          !> Write out Projected Density of State
          logical :: write_pdos

          !> Left boundary of DOS energy window
          real(dp) :: tdos_emin

          !> Right boundary of DOS energy window
          real(dp) :: tdos_emax

          !> TDOS sigma value
          real(dp) :: tdos_sigma

          !> TDOS Number of Points
          integer :: tdos_num_points

          !> PDOS Number of Atoms in Subsystem
          integer :: pdos_num_atoms

          !> PDOS Atom Numbers
          integer, allocatable :: pdos_atoms(:)

  end type estructout_type

  private

  public :: gpmdcov_parse, gpmdcov_estructout_parse

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
    integer, parameter :: nkey_char = 6, nkey_int = 12, nkey_re = 9, nkey_log = 17
    integer :: i
    real(dp) :: realtmp
    character(20) :: dummyc
    logical :: inGPMD
     
    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=50) :: &
         & 'JobName=', 'Var2C=', 'TrajectoryFormat=', 'LangevinMethod=', 'VoltageFile='&
         & ,'FreezeFile=']
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
         &'MyMol', 'Def2', 'PDB', 'Goga', 'None','None']

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
         & 'WriteCoordsEach=',"Var2I=","ReplicateX=","ReplicateY=","ReplicateZ=","PartsToTrack=",&
         & "DumpEach=","MinimizationSteps=","SMDNumPairs=","FineTolEach=","ProfileStartStep=","ProfileStopStep="]
    integer :: valvector_int(nkey_int) = (/ &
         & 1, 1, 0, 0, 0, 0, 0, 0, 0, 5, -1, -1/)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
         & 'VRFactor=','InitialTemperature=','LangevinGamma=','SMDForceConstantStart=',&
         & 'SMDForceConstantEnd=','SMDR0=','CurrentThreshold=',"FineTol=","CoarseTol="]
    real(dp) :: valvector_re(nkey_re) = (/&
         & 0.0_dp, 0.0_dp, 0.01_dp, 0.0_dp,0.2_dp,2.0_dp,0.1_dp,1.0d-5,0.01_dp/)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=50) :: &
         &'DoVelocityRescale=','WriteResidueInTrajectory=','WriteTrajectory=','TrackReactivity=',&
         &'RestartFromDump=','UseLATTE=','HtoD=','LangevinDynamics=','UseSMD=', &
         &'ComputeCurrents=', 'TranslateAndFoldToBox=', 'UseVectSKBlock=', 'ApplyVoltage=','XLBO=','CoarseQMD=',&
         &'UseDispersion=','UseFreeze=']
    logical :: valvector_log(nkey_log) = (/&
         &.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false., &
         &.false.,.True.,.false.,.false.,.true.,.false.,.false.,.false./)

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
    gpmdt%finetoleach = valvector_int(10)
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
    gpmdt%profile_start_step = valvector_int(10)
    gpmdt%profile_stop_step = valvector_int(11)
    
    !Reals
    gpmdt%velresc_fact = valvector_re(1)
    gpmdt%temp0 = valvector_re(2)
    gpmdt%langevin_gamma = valvector_re(3)
    gpmdt%smdforceconststart = valvector_re(4)
    gpmdt%smdforceconstend = valvector_re(5)
    gpmdt%smdr0 = valvector_re(6)
    gpmdt%currthr = valvector_re(7)
    gpmdt%finetol = valvector_re(8)
    gpmdt%coarsetol = valvector_re(9)

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
    gpmdt%xlboON = valvector_log(14)
    gpmdt%coarseqmd = valvector_log(15)
    gpmdt%disp = valvector_log(16)
    gpmdt%freeze = valvector_log(17)

    if(gpmdt%applyv)then 
        gpmdt%voltagef = valvector_char(5)
    endif

    if(gpmdt%freeze)then
        gpmdt%freezef = valvector_char(6)
    endif


  end subroutine gpmdcov_parse


  !> Electronic sturcture output control parser.
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
  subroutine gpmdcov_estructout_parse(filename,estrout)
    implicit none 
    character(len=*), intent(in) :: filename
    type(estructout_type), intent(inout) :: estrout
    integer, parameter :: nkey_char = 2, nkey_int = 2, nkey_re = 3, nkey_log = 3
    integer :: i, counter, ind, j, startatom, endatom, single
    character(20) :: dummyc
    character(20) :: oneline
    logical :: inESTR
     
    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=50) :: &
         & 'TDOSOutputFileName=', 'PDOSOutputFileName=']
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
         &'tdos_output', 'pdos_output']

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
         & 'TDOSNumPoints=', 'PDOSNumAtoms=']
    integer :: valvector_int(nkey_int) = (/ &
         & 1000, 0/)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
         & 'TDOSEmin=', 'TDOSEmax=', 'TDOSSigma=']
    real(dp) :: valvector_re(nkey_re) = (/&
         & -20.0_dp, 20.0_dp, 0.5_dp/)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=50) :: &
         &'WriteTDOS=', 'ComputePDOS=', 'WritePDOS=']
    logical :: valvector_log(nkey_log) = (/&
         &.false., .false., .false./)

    !Start and stop characters
    character(len=50), parameter :: startstop(2) = [character(len=50) :: &
         'ESTRUCTOUT{', '}']

    call prg_parsing_kernel(keyvector_char,valvector_char&
         ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
         keyvector_log,valvector_log,trim(filename),startstop)

    !Characters
    estrout%tdos_output_filename = valvector_char(1)
    estrout%pdos_output_filename = valvector_char(2)
    
    !Integer
    estrout%tdos_num_points = valvector_int(1)
    estrout%pdos_num_atoms = valvector_int(2)
    inESTR = .false.
    if(estrout%pdos_num_atoms > 0) then
      write(*,*) "Reading PDOS atoms"
      open(1, file=trim(filename))
      do i = 1,10000
        read(1,*) dummyc
        !write(*,*) trim(adjustl(dummyc))
        if(trim(adjustl(dummyc)) == "PDOSAtoms[")then
          exit
        end if
        if(trim(adjustl(dummyc)) == "ESTRUCTOUT{") inESTR = .true.
                
        if(trim(dummyc) == "}" .and. inESTR)then
          write(*,*)'ERROR: No PDOS Atoms defined'
          write(*,*)"Here is an example block you should add to the"
          write(*,*)"input file"
          write(*,*)""
          write(*,*)"PDOSAtoms["
          write(*,*)"   1:10"
          write(*,*)"   12"
          write(*,*)"   15:30"
          write(*,*)"]"
          write(*,*)""
          write(*,*)"This will assign atoms 1 through 10, 12, and 15" 
          write(*,*)"through 30 to the PDOS subsystem"
          stop
        end if
      end do
      allocate(estrout%pdos_atoms(estrout%pdos_num_atoms))
      counter=1
      do i = 1,estrout%pdos_num_atoms
        !! Read through PDOS atom lines: maximum number of lines is number of atoms
        !! if each atom number is on its own line
        read(1,*) oneline
        !! Exit if the end of the PDOSAtoms bracket is reached
        !! This will happen if one or more atom groups is given
        if(trim(adjustl(oneline)) == "]") then
                exit
        endif
        !! If an atom group is given, there will be a colon separating the 
        !! first and last atom number in the group
        ind = index(oneline,":")
        if(ind .gt. 0) then
            read(oneline(1:ind-1), '(i5)') startatom
            read(oneline(ind+1:len(oneline)), '(i5)') endatom
            do j = startatom, endatom
              !write(*,*) "Adding atom to array ",j
              estrout%pdos_atoms(counter) = j
              counter = counter + 1
            enddo
        else
                !write(*,*) "Adding single atom to array ",oneline
                read(oneline,'(i5)') single
                estrout%pdos_atoms(counter) = single
                counter = counter + 1
        endif
      end do
      write(*,*)""
      close(1)
    endif
    write(*,*) "PDOS Atoms ",estrout%pdos_atoms
    !call gpmdcov_msI("gpmdcov_parser","PDOS Atom Array " &
    !     & // to_string(estrout%pdos_atoms),lt%verbose,myRank)


    !Reals
    estrout%tdos_emin = valvector_re(1)
    estrout%tdos_emax = valvector_re(2)
    estrout%tdos_sigma = valvector_re(3)

    !Logs
    estrout%write_tdos = valvector_log(1)
    estrout%compute_pdos = valvector_log(2)
    estrout%write_pdos = valvector_log(3)

  end subroutine gpmdcov_estructout_parse

end module gpmdcov_parser_mod
