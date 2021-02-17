
  !>  Main MD loop
  !!  This routine performs the MD loops up to "ls%mdsteps"
  !!
  subroutine gpmdcov_MDloop()

    use gpmdcov_vars

    real(dp) :: mls_ii

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
        EKIN = EKIN + sy%mass(i)*(sy%velocity(1,i)**2+sy%velocity(2,i)**2+sy%velocity(3,i)**2)
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
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Preliminars "//to_string(mls() - mls_ii)//" ms"
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul1 "//to_string(mls() - mls_ii)//" ms"

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
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul2 "//to_string(mls() - mls_ii)//" ms"

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
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul3 "//to_string(mls() - mls_ii)//" ms"
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
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for gpmdcov_Part "//to_string(mls() - mls_i)//" ms"
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul4 "//to_string(mls() - mls_ii)//" ms"
      !> Reprg_initialize parts.
      mls_i = mls()
      call gpmdcov_InitParts()
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for gpmdcov_InitParts "//to_string(mls() - mls_i)//" ms"
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul5 "//to_string(mls() - mls_ii)//" ms"

      mls_i = mls()
      call prg_xlbo_nint(sy%net_charge,n,n_0,n_1,n_2,n_3,n_4,n_5,mdstep,xl)
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for prg_xlbo_nint "//to_string(mls() - mls_i)//" ms"
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul5 "//to_string(mls() - mls_ii)//" ms"

      mls_i = mls()
      Nr_SCF_It = xl%maxscfiter;

      !> Use SCF the first M_prg_init MD steps
      if(mdstep < xl%minit)then
        Nr_SCF_It = xl%maxscfInitIter
      else
        Nr_SCF_It = xl%maxscfiter
      endif

      !> SCF loop
      if(Nr_SCF_It.ne.0)call gpmdcov_dm_min(Nr_SCF_It,n,.true.)

      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul6 "//to_string(mls() - mls_ii)//" ms"

      sy%net_charge = n

      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for gpmdcov_DM_Min_1 "//to_string(mls() - mls_i)//" ms"

      mls_i = mls()
      write(*,*)"Aditional DM construction ..."
      call gpmdcov_DM_Min(1,sy%net_charge,.false.)
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for gpmdcov_DM_Min_2 "//to_string(mls() - mls_i)//" ms"
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul7 "//to_string(mls() - mls_ii)//" ms"

      mls_i = mls()
      call gpmdcov_EnergAndForces(n)
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for gpmd_EnergAndForces "//to_string(mls() - mls_i)//" ms"
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul8 "//to_string(mls() - mls_ii)//" ms"

      mls_i = mls()
      !> Adjust forces for the linearized XLBOMD functional
      call prg_xlbo_fcoulupdate(Coul_Forces,sy%net_charge,n)
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul9 "//to_string(mls() - mls_ii)//" ms"

      !> Total XLBOMD force
      !       sy%force = SKForce + PairForces + FPUL + Coul_Forces + FSCOUL;
      sy%force = collectedforce + PairForces + Coul_Forces

      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul10 "//to_string(mls() - mls_ii)//" ms"
      !> Integrate second 1/2 of leapfrog step
      call halfVerlet(sy%mass,sy%force,lt%timestep,sy%velocity(1,:),sy%velocity(2,:),sy%velocity(3,:))

      if(lt%verbose >= 3 .and. myRank == 1)then
        call prg_write_trajectory(sy,mdstep,5,lt%timestep,"trajectory","pdb")
      endif
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Rest "//to_string(mls() - mls_i)//" ms"
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for Cumul11 "//to_string(mls() - mls_ii)//" ms"

      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for MD iter "//to_string(mls() - mls_ii)//" ms"

      ! Save MD state each 120 steps
      if(mod(mdstep,150) == 0)call gpmdcov_dump()

    enddo
    ! End of MD loop.

  end subroutine gpmdcov_MDloop


