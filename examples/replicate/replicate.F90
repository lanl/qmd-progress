!> Program to replicate the system coordinates along the lattice vectors.
!!
!! \ingroup PROGRAMS
!!
!! Example using this program:
!!
!!     \verbatim  replicate input.xyz output.xyz 2 2 2 \endverbatim
!!
program replicate

  !PROGRESS lib modes.
  use prg_system_mod
  integer, parameter                ::  dp = kind(1.0d0)
  type(system_type)                 ::  sy
  integer                           ::  nx, ny, nz
  character(30)                     ::  filein, fileout, namein, nameout
  character(3)                      ::  extin, extout
  character(2)                      ::  nxc, nyc, nzc
  character(1), allocatable         ::  tempc(:)
  character(len=30)                 ::  tempcflex
  integer                           ::  lenc

  call getarg(1, filein)
  call getarg(2, fileout)
  call getarg(3, nxc)
  call getarg(4, nyc)
  call getarg(5, nzc)

  if(filein == "")then
     write(*,*)""
     write(*,*)"Usage:"
     write(*,*)""
     write(*,*)"  $ replicate <filein> <fileout> <nx> <ny> <nz>"
     write(*,*)""
     write(*,*)"<filein>:  Input coordinates file "
     write(*,*)"<filein>:  Output coordinates file "
     write(*,*)"<n* >: Integer lattice translation"
     write(*,*)""
     stop
  endif

  read(nxc,*) nx
  read(nyc,*) ny
  read(nzc,*) nz

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

  call prg_parse_system(sy,adjustl(trim(namein)),extin)
  call prg_replicate(sy%coordinate,sy%symbol,sy%lattice_vector,nx,ny,nz)
  call prg_centeratbox(sy%coordinate,sy%lattice_vector,0)

  sy%nats = size(sy%coordinate, dim=2)

  call prg_write_system(sy,adjustl(trim(nameout)),extout)

end program replicate
