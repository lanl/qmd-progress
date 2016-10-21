!> Molecular dynamics module.
!! @ingroup LATTE 
!!
!! \brief This module will be used to perform operations related with MD.
!!      
module md_latte_mod

  use openfiles_mod
  use bml   
  use kernelparser_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)
  
  real(dp), parameter :: F2V = 0.01602176487_dp/1.660548782_dp

  !> General md type
  !!
  type, public :: md_type
     
     character(20) :: integrator
     
     integer :: verbose
     
     !> MD time step in fs
     real(dp) :: timestep

     !> Max SCF iterations for the first minit steps.
     integer :: mdSteps    

     !> Max SCF iterations for the first minit steps.
     integer :: writeeach        

  end type md_type 

  public :: parse_md, halfVerlet, updatecoords

contains 

  !> The parser for md solver.
  !!  
  subroutine parse_md(md,filename)

    implicit none
    type(md_type), intent(inout) :: md
    integer, parameter :: nkey_char = 1, nkey_int = 3, nkey_re = 1, nkey_log = 1
    character(len=*) :: filename    
    
    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=100) :: &
      'JobName=']
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
      'MyJob']

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
    'Verbose=','MDSteps=','WriteEach=']                                   
    integer :: valvector_int(nkey_int) = (/ &
       0, 100, 10 /)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
      'TimeStep=' ]
    real(dp) :: valvector_re(nkey_re) = (/&
         0.25  /)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=100) :: &
      'Log1=']
    logical :: valvector_log(nkey_log) = (/&
     .false. /)

    !Start and stop characters
    character(len=50), parameter :: startstop(2) = [character(len=50) :: &
      'MD{', '}']
     
    call parsing_kernel(keyvector_char,valvector_char&
    ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
    keyvector_log,valvector_log,trim(filename),startstop)

    !Characters 
    md%integrator = valvector_char(1)

    !Reals         
    md%timestep = valvector_re(1)

    !Logicals    

    !Integers
    md%verbose = valvector_int(1)
    md%mdsteps = valvector_int(2)
    md%writeeach = valvector_int(3)

  end subroutine parse_md

  !> Half Verlet integration.
  !! \param mass Mass of every atom in the system.
  !! \param FTOT Total force on every atom.
  !! \param timestep Time step for Verlet integration. 
  !! \param VX Velocities for the x direction. 
  !! \param VY Velocities for the y direction. 
  !! \param VZ Velocities for the z direction. 
  subroutine halfVerlet(mass,FTOT,timestep,VX,VY,VZ)
    implicit none
    real(dp), intent(in) :: mass(:)
    real(dp), intent(in) :: FTOT(:,:)
	  real(dp), intent(in) :: timestep
	  real(dp), intent(inout) :: VX(:), VY(:), VZ(:)

      VX(:) = VX(:) + 0.5_dp*timestep*(F2V*FTOT(1,:)/mass(:))
      VY(:) = VY(:) + 0.5_dp*timestep*(F2V*FTOT(2,:)/mass(:))      
      VZ(:) = VZ(:) + 0.5_dp*timestep*(F2V*FTOT(3,:)/mass(:))
    
  end subroutine halfVerlet

  !> Update atomic positions.
  !! \param origin Coordinate origin.
  !! \param lattice_vector Lattice vectors of the system.
  !! \param timestep Time step for Verlet integration.   
  !! \param VX Velocities for the x direction. 
  !! \param VY Velocities for the y direction. 
  !! \param VZ Velocities for the z direction. 
  !! \param coordinate Coordinates of every atom in the system.
  subroutine updatecoords(origin,lattice_vector,timestep,VX,VY,VZ,coordinate)
    implicit none
    real(dp), intent(in) :: origin(:),lattice_vector(:,:)
    real(dp), intent(in) :: VX(:), VY(:), VZ(:), timestep
    real(dp), intent(inout) ::  coordinate(:,:)
    integer :: i

    do i=1,size(coordinate,dim=2)
      coordinate(1,i) = coordinate(1,i) + timestep*VX(i)
      if(coordinate(1,i)-origin(1).gt.lattice_vector(1,1))then
        coordinate(1,i) = coordinate(1,i) - lattice_vector(1,1)
      endif
      if(coordinate(1,i)-origin(1).lt.0.0_dp)then
        coordinate(1,i) = coordinate(1,i) + lattice_vector(1,1)  
      endif
      coordinate(2,i) = coordinate(2,i) + timestep*VY(i)
      if(coordinate(2,i)-origin(2).gt.lattice_vector(2,2))then
        coordinate(2,i) = coordinate(2,i) - lattice_vector(2,2)   
      endif 
      if(coordinate(2,i)-origin(2).lt.0.0_dp)then
        coordinate(2,i) = coordinate(2,i) + lattice_vector(2,2)   
      endif
      coordinate(3,i) = coordinate(3,i) + timestep*VZ(i)
      if(coordinate(3,i)-origin(3).gt.lattice_vector(3,3))then
        coordinate(3,i) = coordinate(3,i) - lattice_vector(3,3)   
      endif
      if(coordinate(3,i)-origin(3).lt.0.0_dp)then
        coordinate(3,i) = coordinate(3,i) + lattice_vector(3,3)
      endif
    enddo          

  end subroutine updatecoords      

end module md_latte_mod

