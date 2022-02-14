module gpmdcov_MDloop_mod

contains

  !>  Main MD loop
  !!  This routine performs the MD loops up to "ls%mdsteps"
  !!
  subroutine gpmdcov_MDloop()

    use gpmdcov_vars
    use gpmdcov_dm_min_mod
    use gpmdcov_energandforces_mod
    use gpmdcov_part_mod
    use gpmdcov_writeout_mod
    use gpmdcov_kernel_mod

    real(dp) :: mls_ii, resnorm
    real(dp), allocatable :: kernelTimesRes(:)
    real(dp), allocatable :: KK0Res(:)


    call gpmdcov_msI("gpmdcov_MDloop","In gpmdcov_MDloop ...",lt%verbose,myRank)

    do mdstep = 1,lt%mdsteps

      mls_ii = mls()

      if(myRank == 1)then
        write(*,*)""
        write(*,*)"         #######################"
        write(*,*)"           MDStep =",mdstep
        write(*,*)"         #######################"
        write(*,*)""
      endif

      !> Get Kinetic energy
      EKIN = 0.0_dp
      do i=1,sy%nats
        EKIN = EKIN + &
       & sy%mass(i)*(sy%velocity(1,i)**2+sy%velocity(2,i)**2+sy%velocity(3,i)**2)
      enddo
      EKIN = 0.5_dp*MVV2KE*EKIN

      !! Statistical temperature in Kelvin
      Temp = (2.0_dp/3.0_dp)*KE2T*EKIN/real(sy%nats,dp);
      !! Total Energy in eV
      Energy = EKIN + EPOT;
      !! Time in fs
      Time = mdstep*lt%timestep;

      if(myRank == 1)then
        write(*,*)"Time [fs] = ",Time
        write(*,*)"Energy Kinetic [eV] = ",EKIN
        write(*,*)"Energy Potential [eV] = ",EPOT
        write(*,*)"Energy Total [eV] = ",Energy
        write(*,*)"Temperature [K] = ",Temp
      endif

      call gpmdcov_msI("gpmdcov_MDloop","Time for Preliminars "//to_string(mls() - mls_ii)//" ms",lt%verbose,myRank)

      !> First 1/2 of Leapfrog step
      if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_start(dyn_timer,"Half Verlet")
      call halfVerlet(sy%mass,sy%force,lt%timestep,sy%velocity(1,:),sy%velocity(2,:),sy%velocity(3,:))
      if(lt%verbose >= 1) call prg_timer_stop(dyn_timer,1)

      if(myRank == 1 .and. lt%verbose.GE.5)then
        write(*,*)"Velocities"
        do i = 1,sy%nats
          write(*,*)i,sy%velocity(1,i),sy%velocity(2,i),sy%velocity(3,i)
        enddo
      endif

      !> Update positions
      if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_start(dyn_timer,"Update positions")
      call updatecoords(origin,sy%lattice_vector,lt%timestep,sy%velocity(1,:),sy%velocity(2,:),sy%velocity(3,:),sy%coordinate)
      if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_stop(dyn_timer,1)

#ifdef DO_MPI
      if (getNRanks() .gt. 1) then
        call prg_sumRealReduceN(sy%coordinate(1,:), sy%nats)
        call prg_sumRealReduceN(sy%coordinate(2,:), sy%nats)
        call prg_sumRealReduceN(sy%coordinate(3,:), sy%nats)
        sy%coordinate = sy%coordinate/real(getNRanks(),dp)
      endif
#endif

      mls_i = mls()
      if(lt%doKernel)then
        if(kernel%kernelType == "ByParts")then
          allocate(kernelTimesRes(sy%nats))
          if(mdstep.le.1)then 
                n = sy%net_charge
          endif
          if(mdstep > 1 .and. kernel%rankNUpdate > 0 .and. &
                & mod(mdstep,kernel%updateEach) == 0)then
           mls_ii = mls()
                call gpmdcov_rankN_update_byParts(sy%net_charge,n,syprt,syprtk,kernel%rankNUpdate,KK0Res)
                call gpmdcov_msI("gpmdcov_MDloop","Time for gpmdcov_rankN_update_byParts"//to_string(mls() - mls_ii)//" ms",lt%verbose,myRank)
          else
                mls_ii = mls()
                call gpmdcov_applyKernel(sy%net_charge,n,syprtk,KK0Res)
                call gpmdcov_msI("gpmdcov_MDloop","Time for gpmdcov_applyKernel"//to_string(mls() - mls_ii)//" ms",lt%verbose,myRank)
          endif
          write(*,*)"After Kernel Apply"
          write(*,*)"Before ninit"
          call prg_xlbo_nint_kernelTimesRes(sy%net_charge,n,n_0,&
          &n_1,n_2,n_3,n_4,n_5,mdstep,KK0Res,xl)
          write(*,*)"After ninit"
          deallocate(kernelTimesRes)
        else
          call prg_xlbo_nint_kernel(sy%net_charge,n,n_0,n_1,n_2,n_3,n_4,n_5,mdstep,Ker,xl)
        endif
      else
        call prg_xlbo_nint(sy%net_charge,n,n_0,n_1,n_2,n_3,n_4,n_5,mdstep,xl)
      endif
      call gpmdcov_msI("gpmdcov_MDloop","Time for prg_xlbo_nint"//to_string(mls() - mls_i)//" ms",lt%verbose,myRank)

      !> Update neighbor list (Actialized every nlisteach times steps)
      if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_start(dyn_timer,"Build Nlist")
      if(mod(mdstep,lt%nlisteach) == 0 .or. mdstep == 0 .or. mdstep == 1)then
        call destroy_nlist(nl)
        call build_nlist_int(sy%coordinate,sy%lattice_vector,coulcut,nl,lt%verbose)
      endif
      if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_stop(dyn_timer,1)

      !> Repartition.
      ! This builds the new graph.
      mls_i = mls()
      call gpmdcov_Part()
      call gpmdcov_msI("gpmdcov_MDloop","Time for gpmdcov_Part&
      &"//to_string(mls() - mls_i)//" ms",lt%verbose,myRank)
      !> Reprg_initialize parts.
      mls_i = mls()
      call gpmdcov_InitParts()
      call gpmdcov_msI("gpmdcov_MDloop","Time for gpmdcov_InitParts &
      &"//to_string(mls() - mls_i)//" ms",lt%verbose,myRank)


      if(mdstep >= 2) resnorm =  norm2(sy%net_charge - n)/sqrt(dble(sy%nats))
      call gpmdcov_msI("gpmdcov_MDloop","ResNorm = "//to_string(resnorm),lt%verbose,myRank)
      
      mls_i = mls()
      Nr_SCF_It = xl%maxscfiter;
      !> Use SCF the first M_prg_init MD steps
      if(mdstep < xl%minit)then
        Nr_SCF_It = xl%maxscfInitIter
      else
        Nr_SCF_It = xl%maxscfiter
      endif

      !> SCF loop

      if(newPart) Nr_SCF_It = xl%maxscfInitIter

      if(Nr_SCF_It.ne.0)then 
        if(eig)then 
                call gpmdcov_dm_min(Nr_SCF_It,n,.true.)
        else
                call gpmdcov_dm_min_Eig(Nr_SCF_It,n,.true.)
        endif
      endif

      sy%net_charge = n

      call gpmdcov_msI("gpmdcov_MDloop","Time for gpmdcov_DM_Min_1 &
      &"//to_string(mls() - mls_i)//" ms",lt%verbose,myRank)

      mls_i = mls()
      
      if(eig)then 
        call gpmdcov_DM_Min(1,sy%net_charge,.false.)
      else
        call gpmdcov_DM_Min_Eig(1,sy%net_charge,.false.)
      endif
      
      call gpmdcov_msI("gpmdcov_MDloop","Time for gpmdcov_DM_Min_2 &
      &"//to_string(mls() - mls_i)//" ms",lt%verbose,myRank)

      mls_i = mls()
      call gpmdcov_EnergAndForces(n)
      call gpmdcov_msI("gpmdcov_MDloop","Time for gpmd_EnergAndForces &
      &"//to_string(mls() - mls_i)//" ms",lt%verbose,myRank)

      mls_i = mls()
      !> Adjust forces for the linearized XLBOMD functional
      call prg_xlbo_fcoulupdate(Coul_Forces,sy%net_charge,n)

      !> Total XLBOMD force
      !       sy%force = SKForce + PairForces + FPUL + Coul_Forces +
      !       FSCOUL;
      sy%force = collectedforce + PairForces + Coul_Forces

      !> Integrate second 1/2 of leapfrog step
      call halfVerlet(sy%mass,sy%force,lt%timestep,sy%velocity(1,:),sy%velocity(2,:),sy%velocity(3,:))

      if(lt%verbose >= 3 .and. myRank == 1)then
        call prg_write_trajectory(sy,mdstep,5,lt%timestep,"trajectory","pdb")
      endif

      call gpmdcov_msI("gpmdcov_MDloop","Time for MD iter &
      &"//to_string(mls() - mls_ii)//" ms",lt%verbose,myRank)

      ! Save MD state each 120 steps
      if(mod(mdstep,120) == 0)call gpmdcov_dump()

    enddo
    ! End of MD loop.

  end subroutine gpmdcov_MDloop

end module gpmdcov_MDloop_mod

