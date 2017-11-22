!> Program to compute the dihedral angle given four atom
!! indices and an input file.
!!
!! \ingroup PROGRAMS
!!
!! Example using this program:
!!
!!     \verbatim  getdihedral coords.xyz 1 2 3 4 \endverbatim
!!
!!
program getdihedral

  !PROGRESS lib modes.
  use prg_system_mod
  integer, parameter                ::  dp = kind(1.0d0)
  real(dp)                          ::  dihedral, mv1, mv2, v1(3), v2(3)
  real(dp)                          ::  dotprod, cosdir, v2xv20(3), v1xv10(3)
  real(dp)                          ::  v10(3),v20(3), cprod(3), normcprod, sindir
  type(system_type)                 ::  sy
  integer                           ::  i, id1,id2,id3,id4
  character(30)                     ::  filein
  character(30)                     ::  namein
  character(3)                      ::  extin
  character(2)                      ::  index1, index2, index3, index4
  character(1), allocatable         ::  tempc(:)
  character(len=30)                 ::  tempcflex
  integer                           ::  l, lenc

  call getarg(1, filein)
  call getarg(2, index1)
  call getarg(3, index2)
  call getarg(4, index3)
  call getarg(5, index4)

  if(filein == "")then
     write(*,*)""
     write(*,*)"Usage:"
     write(*,*)""
     write(*,*)"  $ getdihedreal <filein> <index1> <index2> <index3> <index4>"
     write(*,*)""
     write(*,*)"<filein>:  Input coordinates file "
     write(*,*)"<index 1-4>: Indexes to determine the dihedral angle"
     write(*,*)""
     stop
  endif

  read(index1,*) id1
  read(index2,*) id2
  read(index3,*) id3
  read(index4,*) id4

  lenc=len(adjustl(trim(filein)))
  if(.not.allocated(tempc))allocate(tempc(lenc))
  tempcflex = adjustl(trim(filein))
  namein = adjustl(trim(tempcflex(1:lenc-4)))
  extin = adjustl(trim(tempcflex(lenc-2:lenc+1)))

  call prg_parse_system(sy,adjustl(trim(namein)),extin) !Reads the system coordinate.

  call prg_get_dihedral(sy%coordinate,id1,id2,id3,id4,dihedral)

  write(*,*)"Dihedral =",dihedral

end program getdihedral
