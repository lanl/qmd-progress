!> High-level program to change coordinates formats
!! \brief This program can be used to change coordinate formats
!! such as .pdb, .xyz, .dat and .gen
!!
!! \ingroup PROGRAMS
!!
!! Example using this code:
!!
!!     \verbatim  chancoords coords.xyz coords.pdb  \endverbatim
!!
!!
program changecoords

  !PROGRESS lib modes.
  use system_mod

  implicit none
  integer, parameter                ::  dp = kind(1.0d0)
  type(system_type)                 ::  system
  integer                           ::  i
  character(30)                     ::  filein, fileout
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
  call getarg(3, flag)

  if(filein == "")then
    write(*,*)""
    write(*,*)"Usage:"
    write(*,*)""
    write(*,*)"  $ changecoords <filein> <fileout> <flag>"
    write(*,*)""
    write(*,*)"<filein>:  Input coordinates file "
    write(*,*)"<fileout>: Output coordinates file "
    write(*,*)"<flag>:    -c: (center), -f (fold to box)"
    write(*,*)""
    stop
  endif

  write(*,*)"Changing from  ",trim(filein)," to ", trim(fileout)

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

  call parse_system(system,adjustl(trim(namein)),extin) !Reads the system coordinate.

  !Displace the geometric center to the center of the box
  if(flag.EQ."-c")then

    gc= 0.0d0

    do i=1,system%nats
      gc=gc + system%coordinate(:,i)
    enddo
    gc=gc/real(system%nats,dp)

    do i=1,system%nats
      system%coordinate(:,i) = system%coordinate(:,i) - gc
    enddo

    do i=1,system%nats
      system%coordinate(1,i) = system%coordinate(1,i) + system%lattice_vector(1,1)/2.0d0
      system%coordinate(2,i) = system%coordinate(2,i) + system%lattice_vector(2,2)/2.0d0
      system%coordinate(3,i) = system%coordinate(3,i) + system%lattice_vector(3,3)/2.0d0
    enddo

  endif


  if(flag.EQ."-f")then

    gc= 0.0d0

    call translateandfoldtobox(system%coordinate,system%lattice_vector,origin)

  endif

  call write_system(system,adjustl(trim(nameout)),extout) !Reads the system coordinate.

end program changecoords
