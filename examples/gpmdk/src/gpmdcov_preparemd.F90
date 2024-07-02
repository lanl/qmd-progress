module gpmdcov_prepareMD_mod

contains 

  !>  Preparing for MD
  !!
  subroutine gpmdcov_prepareMD(temp0)
    use gpmdcov_vars
    use gpmdcov_writeout_mod
    Implicit none
    real(dp), intent(in) :: temp0

    call gpmdcov_msI("gpmdcov_prepareMD","In gpmdcov_prepareMD...",lt%verbose,myRank)

    !> Initialize velocities
    if(.not.allocated(sy%velocity))then
      allocate(sy%velocity(3,sy%nats))
      sy%velocity(1,:) = 0.0_dp; sy%velocity(2,:) = 0.0_dp; sy%velocity(3,:) = 0.0_dp;    ! Initialize velocities
    endif

    !> Kinetic energy in eV (MVV2KE: unit conversion)
    MVV2KE = 166.0538782_dp/1.602176487_dp

    KE2T = 1.0_dp/0.000086173435_dp

    !> Pressure in bar from eV / Angstrom^3

    EVOVERV2P = 1.602176487_dp*1000.0_dp
    
    if(temp0 > 1.0E-10)then
       if(.not.gpmdt%restartfromdump) call gpmdcov_addVelocity(temp0,sy%velocity,sy%mass)
    endif

  end subroutine gpmdcov_prepareMD

  !>  Adding random velocity 
  !!
  subroutine gpmdcov_addVelocity(temp0,velocity,mass)
     use gpmdcov_vars
     implicit none 
     real(dp), intent(inout) :: velocity(:,:)
     real(dp), intent(in) :: mass(:)
     real(dp), intent(in) :: temp0
     real(dp) :: kinE,iVel,ran
     integer :: nats,myi,myseed,ssize
     integer, allocatable :: seedin(:) 

     nats = size(velocity,dim=2)

     kinE = temp0 * (1.5_dp)*real(nats,dp)/KE2T
     kinE = kinE/(0.5_dp*MVV2KE)
     kinE = kinE/real(nats,dp)

     myseed = 12345
     call random_seed()
     call random_seed(size=ssize)
     allocate(seedin(ssize))
     seedin = myseed
     call random_seed(PUT=seedin)

     !call random_seed(myseed)

     !> Distribute kin E 
     do myi = 1,nats
      iVel = sqrt(kinE/mass(myi))
      call random_number(ran)
      write(*,*)
      velocity(1,myi) = (2.0_dp*ran - 1.0_dp)
      call random_number(ran)
      velocity(2,myi) = (2.0_dp*ran - 1.0_dp)
      call random_number(ran)
      velocity(3,myi) = (2.0_dp*ran - 1.0_dp)

      velocity(:,myi) = velocity(:,myi)/(norm2(velocity(:,myi)))
      velocity(:,myi) = iVel*sy%velocity(:,myi)
     enddo

#ifdef DO_MPI
      if (numRanks .gt. 1) then  !THIS IS VERY IMPORTANT

        call prg_sumRealReduceN(velocity(1,:), nats)
        call prg_sumRealReduceN(velocity(2,:), nats)
        call prg_sumRealReduceN(velocity(3,:), nats)

        velocity = velocity/real(numRanks,dp)
      endif
#endif



  end subroutine gpmdcov_addVelocity



end module gpmdcov_prepareMD_mod
