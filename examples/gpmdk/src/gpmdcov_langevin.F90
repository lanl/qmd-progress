module gpmdcov_langevin_mod
  private
  
  integer, parameter :: dp = kind(1.0d0)

  real(dp), parameter :: F2V = 0.01602176487_dp/1.660548782_dp

  real(dp), parameter :: KE2T = 1.0_dp/0.000086173435_dp

  real(dp), parameter :: MVV2KE = 166.0538782_dp/1.602176487_dp

  real(dp), parameter :: twopi = 8.0_dp*atan(1.0_dp)
  
  public :: gpmdcov_uniform_to_normal, LangevinDVGoga, LangevinDVSivaOne, LangevinDVSivaTwo, LangevinDXGoga
  
contains

  !> Convert uniform randomly distributed rands to standard-normally distributed rands
  !> Requires an input 3-d array of rands with length 2 in last dim
  !> Returns the normally distributed rands in the first position of the last dim
  !! \param rands 3-d array of rands with length 2 in last dim
  
  subroutine gpmdcov_uniform_to_normal(rands)

    real, intent(inout) :: rands(:,:,:)
    real, allocatable :: amp(:,:), tinies(:)
    logical, allocatable :: tinymask(:,:)
    real :: smallest

    if(.not.allocated(amp))then
       allocate(amp(size(rands,dim=1),size(rands,dim=2)))
       allocate(tinymask(size(amp,dim=1),size(amp,dim=2)))
    endif

    !***Important***
    !Replace rands that are too small for the log()
    !If this is left out then eventually the program will
    !   crash
    smallest = tiny(smallest)
    tinymask = rands(:,:,1).lt.smallest
    if(any(tinymask))then
       tinies = pack(rands(:,:,1),tinymask)
       do while(any(tinies.lt.smallest))
          call random_number(tinies)
       end do
       rands(:,:,1) = unpack(tinies,tinymask,rands(:,:,1))
    endif

    !Generate normally distributed rands using Box-Muller
    amp = log(rands(:,:,1))
    amp = sqrt(-2.0 * amp)

    rands(:,:,1) = amp(:,:)*cos(real(twopi)*rands(:,:,2))
    rands(:,:,2) = amp(:,:)*sin(real(twopi)*rands(:,:,2))

  end subroutine gpmdcov_uniform_to_normal
    
  !> Half Verlet Langevin integration. Step one, before coordinate update.
  !! \param mass Mass of every atom in the system.
  !! \param FTOT Total force on every atom.
  !! \param timestep Time step for Verlet integration (fs).
  !! \param VX Velocities for the x direction.
  !! \param VY Velocities for the y direction.
  !! \param VZ Velocities for the z direction.
  !! \param f Friction constant
  !! \param T Temperature (K)
  subroutine LangevinDVGoga(mass,FTOT,timestep,V,DV,gamma,T,rands)

    real(dp), intent(in) :: mass(:)
    real(dp), intent(in) :: FTOT(:,:)
    real(dp), intent(in) :: timestep
    real(dp), intent(inout) :: V(:,:)
    real(dp), intent(out) :: DV(:,:)
    real(dp), intent(in) :: gamma
    real(dp), intent(in) :: T
    real, intent(in) :: rands(:,:,:)
    integer :: i
    real(dp) :: f

    f = 1._dp - exp(-gamma*timestep)

    
    do i = 1,3
       V(i,:) = V(i,:) + timestep*(F2V*FTOT(i,:)/mass(:))
       DV(i,:) = -f*V(i,:) + sqrt(T/KE2T/MVV2KE*f*(2.0_dp-f)/mass(:))*rands(:,i,1)
    end do
    !DVX(:) = -f*VX(:) + sqrt(T/KE2T*f*(2.0_dp-f)/mass(:))*rands(:,1,1)
    !DVY(:) = -f*VY(:) + sqrt(T/KE2T*f*(2.0_dp-f)/mass(:))*rands(:,2,1)
    !DVZ(:) = -f*VZ(:) + sqrt(T/KE2T*f*(2.0_dp-f)/mass(:))*rands(:,3,1)
    !VX(:) = VX(:) + 0.5_dp*timestep*(F2V*FTOT(1,:)/mass(:))
    !VY(:) = VY(:) + 0.5_dp*timestep*(F2V*FTOT(2,:)/mass(:))
    !VZ(:) = VZ(:) + 0.5_dp*timestep*(F2V*FTOT(3,:)/mass(:))

  end subroutine LangevinDVGoga

  subroutine LangevinDVSivaOne(mass,FTOT,timestep,V,gamma,T,rands)

    real(dp), intent(in) :: mass(:)
    real(dp), intent(in) :: FTOT(:,:)
    real(dp), intent(in) :: timestep
    real(dp), intent(inout) :: V(:,:)
    real(dp), intent(in) :: gamma
    real(dp), intent(in) :: T
    real, intent(in) :: rands(:,:,:)
    real(dp) :: a, b, c
    integer :: i

    ! real(dp) :: alpha, beta, norm_alpha, norm_beta, cdf_alpha, cdf_beta, sigfac
    ! real, allocatable :: thresholded_rands(:,:)
    ! logical, allocatable :: mask1(:,:), mask2(:,:)

    ! if(.not.allocated(thresholded_rands))then
    !    allocate(thresholded_rands(size(rands,dim=1),size(rands,dim=2)))
    !    allocate(mask1(size(rands,dim=1),size(rands,dim=2)))
    !    allocate(mask2(size(rands,dim=1),size(rands,dim=2)))
    ! endif
    
    ! alpha = -6.0_dp
    ! beta = 6.0_dp

    ! mask1 = rands(:,:,1).lt.alpha.or.rands(:,:,1).gt.beta
    ! mask2 = rands(:,:,2).lt.alpha.or.rands(:,:,2).gt.beta

    ! if(ANY(mask1.or.mask2))then
    !    write(*,*)"Langevin random velocity and correction are both out of range"
    !    stop
    ! endif

    ! thresholded_rands = merge(rands(:,:,2),rands(:,:,1),mask1)
    
    ! norm_alpha = exp(-alpha*alphs/2.0_dp)/sqrt(twopi)
    ! norm_beta = exp(-beta*beta/2.0_dp)/sqrt(twopi)
    ! cdf_alpha = (1. + erf(alpha/sqrt(2.0_dp)))/2.0_dp
    ! cdf_beta = (1. + erf(beta/sqrt(2.0_dp)))/2.0_dp
    ! sigfac = sqrt(1.0_dp - (beta*norm_beta-alpha*norm_alpha)/(cdf_beta - cdf_alpha) - ((norm_beta-norm_alpha)/(cdf_beta-cdf_alpha)))
    
    a = exp(-gamma*timestep)
    c = 2.0_dp/gamma/timestep
    b = sqrt(c*tanh(1.0_dp/c))
    
    do i = 1,3
       !V(i,:) = sqrt(a)*V(i,:) + sqrt((1.0_dp-a)*T/KE2T/MVV2KE/mass(:))*thresholded_rands(:,i)/sigfac
!       do j = 1,size(mass)
!          write(*,*)mass(j),rands(j,i,1)
!          V(i,j) = sqrt(a)*V(i,j) + sqrt((1.0_dp-a)*T/KE2T/MVV2KE/mass(j))*rands(j,i,1)
!       enddo
       V(i,:) = sqrt(a)*V(i,:) + sqrt((1.0_dp-a)*T/KE2T/MVV2KE/mass(:))*rands(:,i,1)
       V(i,:) = V(i,:) + b*timestep/2.0_dp*(F2V*FTOT(i,:)/mass(:))
    end do
    
    !DVX(:) = -f*VX(:) + sqrt(T/KE2T*f*(2.0_dp-f)/mass(:))*rands(:,1,1)
    !DVY(:) = -f*VY(:) + sqrt(T/KE2T*f*(2.0_dp-f)/mass(:))*rands(:,2,1)
    !DVZ(:) = -f*VZ(:) + sqrt(T/KE2T*f*(2.0_dp-f)/mass(:))*rands(:,3,1)
    !VX(:) = VX(:) + 0.5_dp*timestep*(F2V*FTOT(1,:)/mass(:))
    !VY(:) = VY(:) + 0.5_dp*timestep*(F2V*FTOT(2,:)/mass(:))
    !VZ(:) = VZ(:) + 0.5_dp*timestep*(F2V*FTOT(3,:)/mass(:))

  end subroutine LangevinDVSivaOne

  subroutine LangevinDVSivaTwo(mass,FTOT,timestep,V,gamma,T,rands)

    real(dp), intent(in) :: mass(:)
    real(dp), intent(in) :: FTOT(:,:)
    real(dp), intent(in) :: timestep
    real(dp), intent(inout) :: V(:,:)
    real(dp), intent(in) :: gamma
    real(dp), intent(in) :: T
    real, intent(in) :: rands(:,:,:)
    real(dp) :: a, b, c
    integer :: i
    ! real(dp) :: alpha, beta, norm_alpha, norm_beta, cdf_alpha, cdf_beta, sigfac
    ! real, allocatable :: thresholded_rands(:,:)
    ! logical, allocatable :: mask1(:,:), mask2(:,:)

    ! if(.not.allocated(thresholded_rands))then
    !    allocate(thresholded_rands(size(rands,dim=1),size(rands,dim=2)))
    !    allocate(mask1(size(rands,dim=1),size(rands,dim=2)))
    !    allocate(mask2(size(rands,dim=1),size(rands,dim=2)))
    ! endif
    
    ! alpha = -6.0_dp
    ! beta = 6.0_dp

    ! mask1 = rands(:,:,1).lt.alpha.or.rands(:,:,1).gt.beta
    ! mask2 = rands(:,:,2).lt.alpha.or.rands(:,:,2).gt.beta

    ! if(ANY(mask1.or.mask2))then
    !    write(*,*)"Langevin random velocity and correction are both out of range"
    !    stop
    ! endif

    ! thresholded_rands = merge(rands(:,:,2),rands(:,:,1),mask1)
        
    ! norm_alpha = exp(-alpha*alphs/2.0_dp)/sqrt(twopi)
    ! norm_beta = exp(-beta*beta/2.0_dp)/sqrt(twopi)
    ! cdf_alpha = (1. + erf(alpha/sqrt(2.0_dp)))/2.0_dp
    ! cdf_beta = (1. + erf(beta/sqrt(2.0_dp)))/2.0_dp
    ! sigfac = sqrt(1.0_dp - (beta*norm_beta-alpha*norm_alpha)/(cdf_beta - cdf_alpha) - ((norm_beta-norm_alpha)/(cdf_beta-cdf_alpha)))

    a = exp(-gamma*timestep)
    c = 2.0_dp/gamma/timestep
    b = sqrt(c*tanh(1.0_dp/c))
    
    do i = 1,3
       V(i,:) = V(i,:) + b*timestep/2.0_dp*(F2V*FTOT(i,:)/mass(:))
       !V(i,:) = sqrt(a)*V(i,:) + sqrt((1.0_dp-a)*T/KE2T/MVV2KE/mass(:))*thresholded_rands(:,i)/sigfac
       V(i,:) = sqrt(a)*V(i,:) + sqrt((1.0_dp-a)*T/KE2T/MVV2KE/mass(:))*rands(:,i,1)
    end do
    
    !DVX(:) = -f*VX(:) + sqrt(T/KE2T*f*(2.0_dp-f)/mass(:))*rands(:,1,1)
    !DVY(:) = -f*VY(:) + sqrt(T/KE2T*f*(2.0_dp-f)/mass(:))*rands(:,2,1)
    !DVZ(:) = -f*VZ(:) + sqrt(T/KE2T*f*(2.0_dp-f)/mass(:))*rands(:,3,1)
    !VX(:) = VX(:) + 0.5_dp*timestep*(F2V*FTOT(1,:)/mass(:))
    !VY(:) = VY(:) + 0.5_dp*timestep*(F2V*FTOT(2,:)/mass(:))
    !VZ(:) = VZ(:) + 0.5_dp*timestep*(F2V*FTOT(3,:)/mass(:))

  end subroutine LangevinDVSivaTwo

  !> Update atomic positions.
  !! \param origin Coordinate origin.
  !! \param lattice_vector Lattice vectors of the system.
  !! \param timestep Time step for Verlet integration.
  !! \param VX Velocities for the x direction.
  !! \param VY Velocities for the y direction.
  !! \param VZ Velocities for the z direction.
  !! \param coordinate Coordinates of every atom in the system.
  subroutine LangevinDXGoga(origin,lattice_vector,timestep,VX,VY,VZ,DVX,DVY,DVZ,coordinate)
    implicit none
    real(dp), intent(in) :: origin(:),lattice_vector(:,:)
    real(dp), intent(in) :: VX(:), VY(:), VZ(:), timestep
    real(dp), intent(in) :: DVX(:), DVY(:), DVZ(:)
    real(dp), intent(inout) ::  coordinate(:,:)
    integer :: i

    do i=1,size(coordinate,dim=2)
      coordinate(1,i) = coordinate(1,i) + timestep*(VX(i)+0.5_dp*DVX(i))
      if(coordinate(1,i)-origin(1).gt.lattice_vector(1,1))then
        coordinate(1,i) = coordinate(1,i) - lattice_vector(1,1)
      endif
      if(coordinate(1,i)-origin(1).lt.0.0_dp)then
        coordinate(1,i) = coordinate(1,i) + lattice_vector(1,1)
      endif
      coordinate(2,i) = coordinate(2,i) + timestep*(VY(i)+0.5_dp*DVY(i))
      if(coordinate(2,i)-origin(2).gt.lattice_vector(2,2))then
        coordinate(2,i) = coordinate(2,i) - lattice_vector(2,2)
      endif
      if(coordinate(2,i)-origin(2).lt.0.0_dp)then
        coordinate(2,i) = coordinate(2,i) + lattice_vector(2,2)
      endif
      coordinate(3,i) = coordinate(3,i) + timestep*(VZ(i)+0.5_dp*DVZ(i))
      if(coordinate(3,i)-origin(3).gt.lattice_vector(3,3))then
        coordinate(3,i) = coordinate(3,i) - lattice_vector(3,3)
      endif
      if(coordinate(3,i)-origin(3).lt.0.0_dp)then
        coordinate(3,i) = coordinate(3,i) + lattice_vector(3,3)
      endif
    enddo

  end subroutine LangevinDXGoga

end module gpmdcov_langevin_mod
