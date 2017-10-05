!> Program to compute the dihedral angle given four atom
!! indexes and an input file.
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
  real(dp)                          ::  v10(3),v20(3)
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

  open(1,file="tmp")
  write(1,*)index1," ",index2," ",index3," ",index4
  close(1)
  open(1,file="tmp")
  read(1,*)id1,id2,id3,id4
  close(1)
  call system("rm tmp")

  lenc=len(adjustl(trim(filein)))
  if(.not.allocated(tempc))allocate(tempc(lenc))
  tempcflex = adjustl(trim(filein))
  namein = adjustl(trim(tempcflex(1:lenc-4)))
  extin = adjustl(trim(tempcflex(lenc-2:lenc+1)))

  call prg_parse_system(sy,adjustl(trim(namein)),extin) !Reads the system coordinate.

  v1=sy%coordinate(:,id4) - sy%coordinate(:,id3)
  v10=sy%coordinate(:,id2) - sy%coordinate(:,id3)
  v2=sy%coordinate(:,id1) - sy%coordinate(:,id2)
  v20=sy%coordinate(:,id3) - sy%coordinate(:,id2)

  v1xv10(1)=v1(2)*v10(3)-v1(3)*v10(2)
  v1xv10(2)=-(v1(1)*v10(3)-v1(3)*v10(1))
  v1xv10(3)=v1(1)*v10(2)-v1(2)*v10(1)

  v2xv20(1)=v2(2)*v20(3)-v2(3)*v20(2)
  v2xv20(2)=-(v2(1)*v20(3)-v2(3)*v20(1))
  v2xv20(3)=v2(1)*v20(2)-v2(2)*v20(1)

  dotprod = v1xv10(1)*v2xv20(1) + v1xv10(2)*v2xv20(2) + v1xv10(3)*v2xv20(3)
  mv1= sqrt(v1xv10(1)*v1xv10(1) + v1xv10(2)*v1xv10(2) + v1xv10(3)*v1xv10(3))
  mv2= sqrt(v2xv20(1)*v2xv20(1) + v2xv20(2)*v2xv20(2) + v2xv20(3)*v2xv20(3))

  cosdir = dotprod/(mv1*mv2)

  dihedral=acos(-cosdir)

  write(*,*)"Dihedral =",360*dihedral/(2.0*3.14159265359)

end program getdihedral
