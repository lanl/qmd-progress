!> High-level program to move two set of atoms away from each other.
!! \brief This program can be used to move a set of atoms from the rest.
!! Takes .pdb, .xyz, .dat and .gen as inputs
!!
!! \ingroup PROGRAMS
!!
!! Example using this code:
!!
!!     \verbatim  moveatoms coords.xyz coords_out.pdb input.in \endverbatim
!!
!!
program moveatoms

  !PROGRESS lib modes.
  use prg_system_mod
  use prg_kernelparser_mod
  ! use misc
  integer, parameter                ::  dp = kind(1.0d0)
  character(30)                     ::  jobname, typeofmove
  real(dp)                          ::  dir(3), dist, point1(3), point2(3)
  integer                           ::  pos1, pos2, move_atoms(2)
  type(system_type)                 ::  sy
  integer                           ::  i
  character(30)                     ::  filein, fileout, inputfile
  character(30)                     ::  namein, nameout
  character(3)                      ::  extin, extout
  character(2)                      ::  flag
  character(1), allocatable         ::  tempc(:)
  character(len=30)                 ::  tempcflex
  integer                           ::  l, lenc
  real(dp)                          ::  gc(3)
  real(dp), allocatable             ::  origin(:)

  call getarg(1, filein)
  call getarg(2, fileout)
  call getarg(3, inputfile)

  if(filein == "")then
     write(*,*)""
     write(*,*)"Usage:"
     write(*,*)""
     write(*,*)"  $ moveatoms <filein> <fileout> <inputfile>"
     write(*,*)""
     write(*,*)"<filein>:  Input coordinates file "
     write(*,*)"<fileout>: Output coordinates file "
     write(*,*)"<inputfile>: File containing the input values  "
     write(*,*)""
     stop
  endif

  lenc=len(adjustl(trim(filein)))
  if(.not.allocated(tempc))allocate(tempc(lenc))
  tempcflex = adjustl(trim(filein))
  namein = adjustl(trim(tempcflex(1:lenc-4)))
  extin = adjustl(trim(tempcflex(lenc-2:lenc+1)))

  lenc=len(adjustl(trim(fileout)))
  if(.not.allocated(tempc))allocate(tempc(lenc))
  tempcflex = adjustl(trim(fileout))
  nameout = adjustl(trim(tempcflex(1:lenc-4)))
  extout = adjustl(trim(tempcflex(lenc-2:lenc+1)))

  write(*,*)extin,extout

  call parse_move(jobname, typeofmove, dir, dist, pos1, pos2, point1, &
       point2, move_atoms,adjustl(trim(inputfile))) !Reads the system coordinate.
  call prg_parse_system(sy,adjustl(trim(namein)),extin) !Reads the system coordinate.


  select case (typeofmove)
  case ("DirectionAndDistance")
     d=sqrt( dir(1)**2 +  dir(2)**2 +  dir(3)**2)

     dir(1)= dir(1)/d
     dir(2)= dir(2)/d
     dir(3)= dir(3)/d

     dir(1)= dir(1)* dist
     dir(2)= dir(2)* dist
     dir(3)= dir(3)* dist

  case ("AtomToAtom")
     dir(1)=sy%coordinate(1, pos2)-sy%coordinate(1, pos1)
     dir(2)=sy%coordinate(2, pos2)-sy%coordinate(2, pos1)
     dir(3)=sy%coordinate(3, pos2)-sy%coordinate(3, pos1)

   case ("PointToPoint")
      dir(1)=point2(1)-point1(1)
      dir(2)=point2(2)-point1(2)
      dir(3)=point2(3)-point1(3)

  case ("PointToPointAndDistance")
     dir(1)=point2(1)-point1(1)
     dir(2)=point2(2)-point1(2)
     dir(3)=point2(3)-point1(3)

     d=sqrt( dir(1)**2 +  dir(2)**2 +  dir(3)**2)

     dir(1)= dir(1)/d
     dir(2)= dir(2)/d
     dir(3)= dir(3)/d

     dir(1)= dir(1)* dist
     dir(2)= dir(2)* dist
     dir(3)= dir(3)* dist

  case default
     stop "No TypeOfMove defined"
  end select

  !move atoms
  if( move_atoms(2).gt.sy%nats)  move_atoms(2)=sy%nats

  do i= move_atoms(1), move_atoms(2)
     sy%coordinate(1,i) = sy%coordinate(1,i) +  dir(1)
     sy%coordinate(2,i) = sy%coordinate(2,i) +  dir(2)
     sy%coordinate(3,i) = sy%coordinate(3,i) +  dir(3)
  enddo

  call prg_write_system(sy,adjustl(trim(nameout)),extout)

end program moveatoms


!> The parser for this program
!!
subroutine parse_move(jobname, typeofmove, dir, dist, pos1, pos2, point1, &
     point2, move_atoms, filename)

  use prg_kernelparser_mod
  implicit none
  integer, parameter                ::  dp = kind(1.0d0)
  integer, parameter                ::  nkey_char = 2, nkey_int = 4, nkey_re = 10, nkey_log = 1
  character(len=*)                  ::  filename, typeofmove
  character(20)                     ::  jobname
  real(dp)                          ::  dir(3), dist, point1(3), point2(3)
  integer                           ::  pos1, pos2, move_atoms(2)

  !Library of keywords with the respective defaults.
  character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=100) :: &
       'Jobname=','TypeOfMove=']
  character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
       'MyJob','DirectionAndDistance']

  character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
       'Atom1=','Atom2=','ListToMove1=','ListToMove2=']
  integer :: valvector_int(nkey_int) = (/ &
       0,0,0,0/)

  character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
       'X=','Y=','Z=','Distance=','X1=','Y1=','Z1=','X2=','Y2=','Z2=']
  real(dp) :: valvector_re(nkey_re) = (/&
       0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)

  character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=100) :: &
       'DUMMY=']
  logical :: valvector_log(nkey_log) = (/&
       .false./)

  !Start and stop characters
  character(len=50), parameter :: startstop(2) = [character(len=50) :: &
       'MOVE{', '}']

  call prg_parsing_kernel(keyvector_char,valvector_char&
       ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
       keyvector_log,valvector_log,trim(filename),startstop)

  dir(1)= valvector_re(1);  dir(2)= valvector_re(2);  dir(3)= valvector_re(3)
  dist= valvector_re(4)
  point1(1)= valvector_re(5);  point1(2)= valvector_re(6);  point1(3)= valvector_re(7)
  point2(1)= valvector_re(8);  point2(2)= valvector_re(9);  point2(3)= valvector_re(10)
  jobname = valvector_char(1); typeofmove = valvector_char(2)
  pos1= valvector_int(1);  pos2= valvector_int(2)
  move_atoms(1)= valvector_int(3);  move_atoms(2)= valvector_int(4)

end subroutine parse_move
