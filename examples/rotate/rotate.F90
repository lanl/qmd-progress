!> High-level program to rotate a set of atoms within a molecular structure.
!! \brief This program can be used to rotate a set of atoms respect to the rest.
!! Takes .pdb, .xyz, .dat and .gen as inputs
!!
!! \ingroup PROGRAMS
!!
!! Example using this code:
!!
!!     \verbatim  rotate coords.xyz coords_out.pdb input.in \endverbatim
!!
!!
program rotate

  !PROGRESS lib modes.
  use prg_system_mod
  use prg_syrotation_mod
  use prg_kernelparser_mod
  
  integer, parameter                ::  dp = kind(1.0d0)
  type(system_type)                 ::  sy
  character(40)                     ::  filein, fileout, inputfile
  character(30)                     ::  namein, nameout
  character(3)                      ::  extin, extout
  character(1), allocatable         ::  tempc(:)
  character(len=30)                 ::  tempcflex
  integer                           ::  lenc
  type(rotation_type)               ::  rot

  call getarg(1, filein)
  call getarg(2, fileout)
  call getarg(3, inputfile)

  if(filein == "")then
     write(*,*)""
     write(*,*)"Usage:"
     write(*,*)""
     write(*,*)"  $ rotate <filein> <fileout> <inputfile>"
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

  call prg_parse_rotation(rot,adjustl(trim(inputfile))) 
  call prg_parse_system(sy,adjustl(trim(namein)),extin) 
  call prg_rotate(rot,sy%coordinate,1)
  call prg_write_system(sy,"out",'pdb')

end program rotate
