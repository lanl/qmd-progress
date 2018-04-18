!> High-level program to change coordinates formats.
!! \brief This program can be used to change coordinate formats.
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
  use prg_system_mod

  implicit none
  integer, parameter                ::  dp = kind(1.0d0)
  type(system_type)                 ::  sy
  integer                           ::  lenc, indexi
  character(30)                     ::  filein, fileout
  character(30)                     ::  namein, nameout
  character(3)                      ::  extin, extout
  character(20)                     ::  indexc
  character(2)                      ::  flag
  character(1), allocatable         ::  tempc(:)
  character(len=30)                 ::  tempcflex
  real(dp), allocatable             ::  origin(:)

  call getarg(1, filein)
  call getarg(2, fileout)
  call getarg(3, flag)
  call getarg(4, indexc)

  if(filein == "")then
    write(*,*)""
    write(*,*)"Usage:"
    write(*,*)""
    write(*,*)"  $ changecoords <filein> <fileout> <flag>"
    write(*,*)""
    write(*,*)"<filein>:  Input coordinates file "
    write(*,*)"<fileout>: Output coordinates file "
    write(*,*)"<flag>:    -c: (center), -f (fold to box), -w(wrap around) <index>"
    write(*,*)"            <index>: Atom index to be wrapped with the system using PBC"
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

  call prg_parse_system(sy,adjustl(trim(namein)),extin)
  
  select case(flag)
    case("-c")
      call prg_centeratbox(sy%coordinate,sy%lattice_vector,1)
    case("-f")
      call prg_translateandfoldtobox(sy%coordinate,sy%lattice_vector,origin,1)
    case("-w")
      if(indexc == "")then
        write(*,*)"ERROR: Please provide the atom index for flag -w"
        stop
      endif
      read(indexc,*)indexi
      call prg_wraparound(sy%coordinate,sy%lattice_vector,indexi,1)
    case("")
        write(*,*)"Proceeding without transformation ..."
    case default
      write(*,*) "Invalid flag ",flag
      stop
  end select

  call prg_write_system(sy,adjustl(trim(nameout)),extout)

end program changecoords
