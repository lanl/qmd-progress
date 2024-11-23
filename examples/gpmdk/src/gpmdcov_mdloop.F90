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
    use gpmdcov_neighbor_mod
    use gpmdcov_langevin_mod
    use gpmdcov_preparemd_mod

#ifdef USE_NVTX
    use gpmdcov_nvtx_mod
#endif

    real(dp) :: mls_md, mls_md1, mls_md2, resnorm
    real(dp), allocatable :: kernelTimesRes(:), n1(:)
    real(dp), allocatable :: KK0Res(:)
    type(system_type) :: syaux
    integer, allocatable :: seedin(:)
    integer :: myseed,ssize
    real, allocatable :: langevin_rands(:,:,:)
    real(dp), allocatable :: DV(:,:)
    integer :: maxv_atom_axis(2)
    integer :: nanv_atom_axis(2)
    real(dp) :: virial(3,3)
    real(dp) :: ke_tensor(3,3)
    real(dp) :: pressure_tensor(3,3)
    integer :: total_steps
    integer :: cuda_error
    
    interface
       integer(c_int) function cudaProfilerStart() bind(c,name="cudaProfilerStart")
         use iso_c_binding
       end function cudaProfilerStart
    end interface

    interface
       integer(c_int) function cudaProfilerStop() bind(c,name="cudaProfilerStop")
         use iso_c_binding
       end function cudaProfilerStop
    end interface
                                                     
    ! Prepare for Langevin dynamics if needed
    if(gpmdt%langevin)then
       
    !Initialize random number generator
       myseed = 12345
       call random_seed()
       call random_seed(size=ssize)
       allocate(seedin(ssize))
       seedin = myseed
       call random_seed(PUT=seedin)
    !Allocate array of random numbers for Langevin integration
       allocate(langevin_rands(sy%nats,3,2))
       !Allocate delta_V arrays
       allocate(DV(3,sy%nats))
    endif
    
    call gpmdcov_msI("gpmdcov_MDloop","In gpmdcov_MDloop ...",lt%verbose,myRank)
    savets = lt%timestep
    !do mdstep = -1,lt%mdsteps
    if(.not.gpmdt%restartfromdump.and.gpmdt%minimization_steps.ne.0)then
       sy%velocity = 0.0_dp
    endif
    ! Compute box volume
    call gpmdcov_get_vol(sy%lattice_vector,sy%volr)

    total_steps = lt%mdsteps + gpmdt%minimization_steps
    
    do mdstep = 1,total_steps
      !    if(mdstep < 0)then
      !            savets = lt%timestep
      !            lt%timestep = 0
      !    else
      !            lt%timestep = savets
      !    endif


      mls_md = mls()
#ifdef USE_NVTX
      if (mdstep == 10) then
              cuda_error = cudaProfilerStart()
      endif
      if (mdstep == 15) then
              cuda_error = cudaProfilerStop()
      endif     
      call nvtxStartRange("MD_iter",1)
#endif
      mls_md1 = mls()

      if(myRank == 1)then
        write(*,*)""
        write(*,*)"         #######################"
        if(mdstep.le.gpmdt%minimization_steps)then
           write(*,*)"           Min Step =",mdstep
        else
           write(*,*)"           MDStep =",mdstep-gpmdt%minimization_steps
        endif
        write(*,*)"         #######################"
        write(*,*)""
      endif

      maxv_atom_axis = MAXLOC(ABS(sy%velocity))
      call gpmdcov_msI("gpmdcov_MDloop","Maximum Velocity "//to_string(MAXVAL(ABS(sy%velocity)))//" &
        &for (atom,axis) = ("//to_string(maxv_atom_axis(2))//","//to_string(maxv_atom_axis(1))//")",lt%verbose,myRank)

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

      !! Statistical pressure
      do i = 1,3
         do j = 1,3
            ke_tensor(i,j) = MVV2KE*sum(sy%mass(:)*sy%velocity(i,:)*sy%velocity(j,:))
            virial(i,j) = MVV2KE*F2V*sum(sy%coordinate(i,:)*sy%force(j,:))
         enddo
      enddo

      pressure_tensor = EVOVERV2P*(ke_tensor + virial)/sy%volr
      
      if(myRank == 1)then
        write(*,*)"Time [fs] = ",Time
        write(*,*)"Energy Kinetic [eV] = ",EKIN
        write(*,*)"Energy Potential [eV] = ",EPOT
        write(*,*)"Energy Total [eV] = ",Energy
        write(*,*)"Temperature [K] = ",Temp
        write(*,*)"Pressure [bar] = ",pressure_tensor(1,1)+pressure_tensor(2,2)+pressure_tensor(3,3)
      endif

      call gpmdcov_msI("gpmdcov_MDloop","Time for Preliminaries "//to_string(mls() - mls_md1)//" ms",lt%verbose,myRank)
      mls_md1 = mls()

      if(.not.gpmdt%langevin)then

         !> First 1/2 of Leapfrog step
         call gpmdcov_msMem("gpmdcov_mdloop", "Before halfVerlet",lt%verbose,myRank)
         if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_start(dyn_timer,"Half Verlet")
         call halfVerlet(sy%mass,sy%force,lt%timestep,sy%velocity(1,:),sy%velocity(2,:),sy%velocity(3,:))
         if(lt%verbose >= 1) call prg_timer_stop(dyn_timer,1)
         call gpmdcov_msMem("gpmdcov_mdloop", "After halfVerlet",lt%verbose,myRank)

         if(myRank == 1 .and. lt%verbose.GE.5)then
            do i = 1,sy%nats
               write(*,*)i,sy%velocity(1,i),sy%velocity(2,i),sy%velocity(3,i)
            enddo
         endif
         !> Update positions
         call gpmdcov_msMem("gpmdcov_mdloop", "Before updatecoords",lt%verbose,myRank)
         if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_start(dyn_timer,"Update positions")
         call updatecoords(origin,sy%lattice_vector,lt%timestep,sy%velocity(1,:),sy%velocity(2,:),sy%velocity(3,:),sy%coordinate)
         if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_stop(dyn_timer,1)
         call gpmdcov_msMem("gpmdcov_mdloop", "After updatecoords",lt%verbose,myRank)
#ifdef DO_MPI
         if (numRanks .gt. 1) then  !THIS IS VERY IMPORTANT
            call prg_sumRealReduceN(sy%coordinate(1,:), sy%nats)
            call prg_sumRealReduceN(sy%coordinate(2,:), sy%nats)
            call prg_sumRealReduceN(sy%coordinate(3,:), sy%nats)
            
            call prg_sumRealReduceN(sy%velocity(1,:), sy%nats)
            call prg_sumRealReduceN(sy%velocity(2,:), sy%nats)
            call prg_sumRealReduceN(sy%velocity(3,:), sy%nats)
            
            sy%coordinate = sy%coordinate/real(numRanks,dp)
            sy%velocity = sy%velocity/real(numRanks,dp)
         endif
#endif

         call gpmdcov_msI("gpmdcov_MDloop","Time for halfVerlet, Update coords, and sumReduce "&
              &//to_string(mls() - mls_md1)//" ms",lt%verbose,myRank)

         mls_md1 = mls()
      endif

      !call prg_translateandfoldtobox(sy%coordinate,sy%lattice_vector,origin)
      !origin = 0.0_dp

      if(gpmdt%langevin.and.gpmdt%langevin_method.eq."Siva")then
         call random_number(langevin_rands)
         call gpmdcov_uniform_to_normal(langevin_rands)
         !if (myRank.eq.0)then
         !   call random_number(langevin_rands)
         !   call uniform_to_normal(langevin_rands)
         !endif
         !if (numranks.gt.1)then
         !   call prg_bcastRealParallel(langevin_rands,size(langevin_rands),0)
         !endif
         call gpmdcov_msI("gpmdcov_MDloop","Langevin Siva integration with (gamma,temp) = ("&
              &//to_string(gpmdt%langevin_gamma)//","//to_string(gpmdt%temp0)//")",lt%verbose,myRank)
         call LangevinDVSivaOne(sy%mass,sy%force,lt%timestep,sy%velocity,gpmdt%langevin_gamma,gpmdt%temp0,langevin_rands)
         call updatecoords(origin,sy%lattice_vector,lt%timestep,sy%velocity(1,:),sy%velocity(2,:),sy%velocity(3,:),sy%coordinate)
#ifdef DO_MPI
         if (numRanks .gt. 1) then  !THIS IS VERY IMPORTANT
            call prg_sumRealReduceN(sy%coordinate(1,:), sy%nats)
            call prg_sumRealReduceN(sy%coordinate(2,:), sy%nats)
            call prg_sumRealReduceN(sy%coordinate(3,:), sy%nats)
            
            call prg_sumRealReduceN(sy%velocity(1,:), sy%nats)
            call prg_sumRealReduceN(sy%velocity(2,:), sy%nats)
            call prg_sumRealReduceN(sy%velocity(3,:), sy%nats)
            
            sy%coordinate = sy%coordinate/real(numRanks,dp)
            sy%velocity = sy%velocity/real(numRanks,dp)
         endif
#endif
      endif

      if(mdstep >= 1)then
        if(lt%doKernel)then

          if(kernel%xlbolevel1)then
            call gpmdcov_msI("gpmdcov_MDloop","Doing XLBO level 1",lt%verbose,myRank)
            if(kernel%kernelType == "ByParts")then
              !propagate n
              allocate(kernelTimesRes(sy%nats))
              if(mdstep.le.1)then
                n = sy%net_charge
                call gpmdcov_applyKernel(sy%net_charge,n,syprtk,KK0Res)
                call prg_xlbo_nint_kernelTimesRes(sy%net_charge,n,n_0,&
                     &n_1,n_2,n_3,n_4,n_5,mdstep,KK0Res,xl)
              endif
              if(mdstep > 1 .and. kernel%rankNUpdate > 0 .and. &
                   & mod(mdstep,kernel%updateEach) == 0)then
                call gpmdcov_msI("gpmdcov_MDloop","Integrating n ...",lt%verbose,myRank)

                !call gpmdcov_applyKernel(sy%net_charge,n,syprtk,KK0Res)
                call prg_xlbo_nint_kernelTimesRes(sy%net_charge,n,n_0,&
                     &n_1,n_2,n_3,n_4,n_5,mdstep,KK0Res,xl)
                !Use n > H >  to get q_min
                ! call gpmdcov_DM_Min_Eig(1,sy%net_charge,.false.)
                !Compute KK0Res
                ! call gpmdcov_rankN_update_byParts(sy%net_charge,n,syprt,syprtk,kernel%rankNUpdate,KK0Res)
                !Use KK0Res to update n to and n1
                ! sy%net_charge = n - KK0Res
                !Compute H again and get q_min from n1 later in the code
              endif
              deallocate(kernelTimesRes)
            else
              STOP "XLBOLevel1 not implemented for other than kernelType= ByParts"
            endif

          else
            if(kernel%kernelType == "ByParts")then
              allocate(kernelTimesRes(sy%nats))
              if(mdstep.le.1)then
                n = sy%net_charge
              endif
              if(mdstep > 1 .and. kernel%rankNUpdate > 0 .and. &
                   & mod(mdstep,kernel%updateEach) == 0)then
                mls_md2 = mls()
                call gpmdcov_msMem("gpmdcov_mdloop", "Before gpmdcov_rankN_update_byParts",lt%verbose,myRank)
#ifdef USE_NVTX
                 call nvtxStartRange("RankN_update",2)
#endif
                call gpmdcov_rankN_update_byParts(sy%net_charge,n,syprt,syprtk,kernel%rankNUpdate,KK0Res)
#ifdef USE_NVTX
                call nvtxEndRange
#endif
                call gpmdcov_msMem("gpmdcov_mdloop", "After gpmdcov_rankN_update_byParts",lt%verbose,myRank)
                call gpmdcov_msI("gpmdcov_MDloop","Time for gpmdcov_rankN_update_byParts "&
                     &//to_string(mls() - mls_md2)//" ms",lt%verbose,myRank)
                if(gpmdt%trackreactivity)then
                  syaux%nats = sy%nats
                  allocate(syaux%symbol(sy%nats)); syaux%symbol = sy%symbol
                  allocate(syaux%coordinate(3,sy%nats)); syaux%coordinate = sy%coordinate
                  allocate(syaux%net_charge(sy%nats)); syaux%net_charge = KK0Res*100.0_dp
                  if(myRank == 1)then 
                   call prg_write_trajectory(syaux,mdstep,gpmdt%writetreach,lt%timestep,&
                           &adjustl(trim(lt%jobname))//"_reactivity","xyz")
                  endif
                  call prg_destroy_system(syaux)
                endif

              else
                mls_md2 = mls()
                call gpmdcov_msMem("gpmdcov_mdloop", "Before gpmdcov_applyKernel",lt%verbose,myRank)
                call gpmdcov_applyKernel(sy%net_charge,n,syprtk,KK0Res)
                call gpmdcov_msMem("gpmdcov_mdloop", "After gpmdcov_applyKernel",lt%verbose,myRank)
                call gpmdcov_msI("gpmdcov_MDloop","Time for gpmdcov_applyKernel "&
                     &//to_string(mls() - mls_md2)//" ms",lt%verbose,myRank)
              endif
              call gpmdcov_msMem("gpmdcov_mdloop", "Before prg_xlbo_nint_kernelTimesRes",lt%verbose,myRank)
              call prg_xlbo_nint_kernelTimesRes(sy%net_charge,n,n_0,&
                   &n_1,n_2,n_3,n_4,n_5,mdstep,KK0Res,xl)
              call gpmdcov_msMem("gpmdcov_mdloop", "After prg_xlbo_nint_kernelTimesRes",lt%verbose,myRank)
              deallocate(kernelTimesRes)
            else
              call gpmdcov_msMem("gpmdcov_mdloop", "Before prg_xlbo_nint_kernel",lt%verbose,myRank)
              call prg_xlbo_nint_kernel(sy%net_charge,n,n_0,n_1,n_2,n_3,n_4,n_5,mdstep,Ker,xl)
              call gpmdcov_msMem("gpmdcov_mdloop", "After prg_xlbo_nint_kernel",lt%verbose,myRank)
            endif
          endif
        else
          call gpmdcov_msMem("gpmdcov_mdloop", "Before prg_xlbo_nint",lt%verbose,myRank)
          
          if(gpmdt%xlboon)then
                call prg_xlbo_nint(sy%net_charge,n,n_0,n_1,n_2,n_3,n_4,n_5,mdstep,xl)
          else
                n = sy%net_charge
          endif

          if(gpmdt%coarseqmd)then 
                if(mod(mdstep,gpmdt%finetoleach) == 0)then 
                        lt%scftol = gpmdt%finetol !1.0d-5
                else
                        lt%scftol = gpmdt%coarsetol !0.01_dp
                endif
          endif

          call gpmdcov_msMem("gpmdcov_mdloop", "After prg_xlbo_nint",lt%verbose,myRank)
        endif
      else
        n = sy%net_charge
      endif
      call gpmdcov_msI("gpmdcov_MDloop","Time for prg_xlbo_nint "//to_string(mls() - mls_md1)//" ms",lt%verbose,myRank)

      !> Update neighbor list (Actialized every nlisteach times steps)
      mls_md1 = mls()
      if(mod(mdstep,lt%nlisteach) == 0 .or. mdstep == 0 .or. mdstep == 1)then
        call gpmdcov_msMem("gpmdcov_mdloop", "Before build_nlist_int",lt%verbose,myRank)
        call gpmdcov_destroy_nlist(nl,lt%verbose)
        !call destroy_nlist(nl)
        if(nlistSparse)then
#ifdef USE_NVTX
           call nvtxStartRange("build_nlist_sparse_v2",3)
#endif
           call gpmdcov_build_nlist_sparse_v2(sy%coordinate,sy%lattice_vector,coulcut,nl,lt%verbose,myRank,numRanks)
#ifdef USE_NVTX
           call nvtxEndRange
#endif
        else 
#ifdef USE_NVTX
           call nvtxStartRange("build_nlist_full",3)
#endif
           call gpmdcov_build_nlist_full(sy%coordinate,sy%lattice_vector,coulcut,nl,lt%verbose,myRank,numRanks)
#ifdef USE_NVTX
           call nvtxEndRange
#endif
        endif
        !LBox(1) =  sy%lattice_vector(1,1)
        !LBox(2) =  sy%lattice_vector(2,2)
        !LBox(3) =  sy%lattice_vector(3,3)

        !call gpmdcov_nearestneighborlist(nl%nrnnlist,nl%nndist,nl%nnRx,nl%nnRy,nl%nnRz,nl%nnType, &
        !          &sy%coordinate(1,:),sy%coordinate(2,:),sy%coordinate(2,:),LBox,coulcut,sy%nats,200)

        !call build_nlist_int(sy%coordinate,sy%lattice_vector,coulcut,nl,lt%verbose)
        call gpmdcov_msMem("gpmdcov_mdloop", "After build_nlist_int",lt%verbose,myRank)
      endif
      call gpmdcov_msI("gpmdcov_MDloop","Time for build_nlist_int "&
           &//to_string(mls() - mls_md1)//" ms",lt%verbose,myRank)
      !stop
      !> Repartition.
      ! This builds the new graph.
      mls_md1 = mls()
      call gpmdcov_msMem("gpmdcov_mdloop", "Before gpmdcov_Part",lt%verbose,myRank)
#ifdef USE_NVTX
           call nvtxStartRange("Part",4)
#endif

           call gpmdcov_Part(2)
#ifdef USE_NVTX
           call nvtxEndRange
#endif
      call gpmdcov_msMem("gpmdcov_mdloop", "After gpmdcov_Part",lt%verbose,myRank)
      call gpmdcov_msI("gpmdcov_MDloop","Time for gpmdcov_Part &
           &"//to_string(mls() - mls_i)//" ms",lt%verbose,myRank)
      !> Reprg_initialize parts.
      mls_i = mls()
      call gpmdcov_msMem("gpmdcov_mdloop", "Before gpmdcov_InitParts",lt%verbose,myRank)
#ifdef USE_NVTX
      call nvtxStartRange("InitParts",5)
#endif
       !if((mod(mdstep,lt%nlisteach) == 0 ) .or. (mod(mdstep,gsp2%parteach) == 0) &
       !        &.or. mdstep == 0 .or. mdstep == 1) call gpmdcov_InitParts()
       call gpmdcov_InitParts()
#ifdef USE_NVTX
      call nvtxEndRange
#endif
      call gpmdcov_msMem("gpmdcov_mdloop", "After gpmdcov_InitParts",lt%verbose,myRank)
      call gpmdcov_msI("gpmdcov_MDloop","Time for gpmdcov_InitParts &
           &"//to_string(mls() - mls_md1)//" ms",lt%verbose,myRank)

      mls_md1 = mls()
      resnorm = 0.0_dp

      if((mdstep >= 2) .and. (.not. kernel%xlbolevel1)) resnorm =  norm2(sy%net_charge - n)/sqrt(dble(sy%nats))

      Nr_SCF_It = xl%maxscfiter;
      !> Use SCF the first MD steps
      if(mdstep < xl%minit)then
        Nr_SCF_It = xl%maxscfInitIter
      else
        Nr_SCF_It = xl%maxscfiter
      endif

      !> SCF loop
      !if(newPart .and. (.not.lt%dokernel)) then
      !        Nr_SCF_It = xl%maxscfInitIter
      !        newPart = .false.
      !endif


      if(Nr_SCF_It .gt. 0)then
        if(eig)then
          call gpmdcov_msMem("gpmdcov_mdloop", "Before gpmdcov_dm_min",lt%verbose,myRank)
          call gpmdcov_dm_min(Nr_SCF_It,n,.true.)
          call gpmdcov_msMem("gpmdcov_mdloop", "After gpmdcov_dm_min",lt%verbose,myRank)
        else
          call gpmdcov_msMem("gpmdcov_mdloop", "Before gpmdcov_dm_min_Eig",lt%verbose,myRank)
          call gpmdcov_dm_min_Eig(Nr_SCF_It,n,.true.,.false.)
          call gpmdcov_msMem("gpmdcov_mdloop", "After gpmdcov_dm_min_Eig",lt%verbose,myRank)
        endif
      endif

      sy%net_charge = n

#ifdef USE_NVTX
           call nvtxStartRange("DM_min",6)
#endif

      if(gpmdt%xlboON)then
      if(eig)then
        call gpmdcov_msMem("gpmdcov_mdloop", "Before gpmdcov_dm_min",lt%verbose,myRank)
           call gpmdcov_DM_Min(1,sy%net_charge,.false.)
        call gpmdcov_msMem("gpmdcov_mdloop", "After gpmdcov_dm_min",lt%verbose,myRank)
      else
        call gpmdcov_msMem("gpmdcov_mdloop", "Before gpmdcov_dm_min_Eig",lt%verbose,myRank)
        call gpmdcov_DM_Min_Eig(1,sy%net_charge,.false.,.false.)
        call gpmdcov_msMem("gpmdcov_mdloop", "After gpmdcov_dm_min_Eig",lt%verbose,myRank)
      endif
      endif

      call gpmdcov_msI("gpmdcov_MDloop","Time for gpmdcov_DM_Min_1 &
           &"//to_string(mls() - mls_md1)//" ms",lt%verbose,myRank)

#ifdef USE_NVTX
           call nvtxEndRange
#endif
      if(kernel%xlbolevel1)then
        allocate(n1(sy%nats))
        if(mdstep > 1)then
          !sy%net_charge = n
           !Compute KK0Res
#ifdef USE_NVTX
           call nvtxStartRange("RankN_update",2)
#endif
           call gpmdcov_rankN_update_byParts(sy%net_charge,n,syprt,syprtk,kernel%rankNUpdate,KK0Res)
#ifdef USE_NVTX
           call nvtxEndRange
#endif
          !Use KK0Res to update n to and n1
          n1 = n - KK0Res
          sy%net_charge = n1
#ifdef USE_NVTX
          call nvtxStartRange("DM_min_eig",6)
#endif
          call gpmdcov_DM_Min_Eig(1,sy%net_charge,.false.,.false.)
#ifdef USE_NVTX
          call nvtxEndRange
#endif
          !Compute H again and get q_min from n1 later in the code
          resnorm =  norm2(sy%net_charge - n1)/sqrt(dble(sy%nats))
        endif
      endif


      mls_md1 = mls()
      call gpmdcov_msI("gpmdcov_MDloop","ResNorm = "//to_string(resnorm),lt%verbose,myRank)
      if(myRank == 1)then
         if(mdstep.le.gpmdt%minimization_steps)then
            write(*,'(A35,I15,A1,F18.5,A1,ES12.5,A1,ES12.5,A1,ES12.5)')"Minstep, Energy, Egap, Resnorm, Temp", &
                 &mdstep," ", Energy," ", egap_glob," ", resnorm," ", Temp
         else
            write(*,'(A35,I15,A1,F18.5,A1,ES12.5,A1,ES12.5,A1,ES12.5)')"Mdstep, Energy, Egap, Resnorm, Temp", &
                 &mdstep-gpmdt%minimization_steps," ", Energy," ", egap_glob," ", resnorm," ", Temp
         endif
        !write(*,*)"Step, Energy, EGap, Resnorm", mdstep, Energy, egap_glob, resnorm
      endif
#ifdef USE_NVTX
      call nvtxStartRange("EnergAndForces",7)
#endif
      call gpmdcov_msMem("gpmdcov_mdloop", "Before gpmdcov_EnergAndForces",lt%verbose,myRank)
      if(kernel%xlbolevel1)then
        if(mdstep <= 1) n1 = n
        call gpmdcov_EnergAndForces(n1)
        deallocate(n1)
      else
        call gpmdcov_EnergAndForces(n)
      endif
      call gpmdcov_msMem("gpmdcov_mdloop", "After gpmdcov_EnergAndForces",lt%verbose,myRank)
      call gpmdcov_msI("gpmdcov_MDloop","Time for gpmdcov_EnergAndForces &
           &"//to_string(mls() - mls_md1)//" ms",lt%verbose,myRank)
#ifdef USE_NVTX
      call nvtxEndRange
#endif
      mls_md1 = mls()
      !> Adjust forces for the linearized XLBOMD functional
      call gpmdcov_msMem("gpmdcov_mdloop", "Before prg_xlbo_fcoulupdate",lt%verbose,myRank)
      call prg_xlbo_fcoulupdate(Coul_Forces,sy%net_charge,n)
      call gpmdcov_msMem("gpmdcov_mdloop", "After prg_xlbo_fcoulupdate",lt%verbose,myRank)
      call gpmdcov_msI("gpmdcov_MDloop","Time for prg_xlbo_fcoulupdate &
           &"//to_string(mls() - mls_md1)//" ms",lt%verbose,myRank)

      mls_md1 = mls()

      !> Total XLBOMD force
      !       sy%force = SKForce + PairForces + FPUL + Coul_Forces +
      !       FSCOUL;
      sy%force = collectedforce + PairForces + Coul_Forces

      !> Integrate second 1/2 of leapfrog step
      if(gpmdt%dovelresc .eqv. .true.)then
         call gpmdcov_msI("gpmdcov_MDloop","Doing Velocity Rescale",lt%verbose,myRank)
         if(gpmdt%temp0 .ne. 0.0_dp)then
            sy%velocity = sqrt(gpmdt%temp0/Temp)*sy%velocity
         else
            sy%velocity = gpmdt%velresc_fact*sy%velocity
         endif
      endif

      
      if(.not.gpmdt%restartfromdump)then
         if(mdstep.lt.gpmdt%minimization_steps)then
            call gpmdcov_msI("gpmdcov_MDloop","Zeroing velocities during minimization",lt%verbose,myRank)
            sy%velocity = 0.0_dp
         elseif(mdstep.eq.gpmdt%minimization_steps)then
            if(.not.gpmdt%langevin.and.gpmdt%temp0.gt.1.0E-10)then
               call gpmdcov_addVelocity(gpmdt%temp0,sy%velocity,sy%mass)
            else
               sy%velocity = 0.0_dp
            endif
         endif
      endif

      call gpmdcov_msMem("gpmdcov_mdloop", "Before halfVerlet",lt%verbose,myRank)
      if(gpmdt%langevin.and.gpmdt%langevin_method.eq."Siva")then
         call random_number(langevin_rands)
         call gpmdcov_uniform_to_normal(langevin_rands)
!         if (myRank.eq.0)then
!            call random_number(langevin_rands)
!            call uniform_to_normal(langevin_rands)
!         endif
!         if (numranks.gt.1)then
!            call prg_bcastRealParallel(langevin_rands,size(langevin_rands),0)
!         endif
         call LangevinDVSivaTwo(sy%mass,sy%force,lt%timestep,sy%velocity,gpmdt%langevin_gamma,gpmdt%temp0,langevin_rands)
      endif

      if(gpmdt%langevin.and.gpmdt%langevin_method.eq."Goga")then
         call random_number(langevin_rands)
         call gpmdcov_uniform_to_normal(langevin_rands)
         call gpmdcov_msI("gpmdcov_MDloop","Langevin integration with (gamma,temp) = ("&
              &//to_string(gpmdt%langevin_gamma)//","//to_string(gpmdt%temp0)//")",lt%verbose,myRank)
         !First update the velocities using the force field and compute DV
         call LangevinDVGoga(sy%mass,sy%force,lt%timestep,sy%velocity,DV,gpmdt%langevin_gamma,gpmdt%temp0,langevin_rands)
         !Next update the coordinates using the updated velocities and DV
         call LangevinDXGoga(origin,sy%lattice_vector,lt%timestep,sy%velocity(1,:),sy%velocity(2,:),sy%velocity(3,:),&
           &DV(1,:),DV(2,:),DV(3,:),sy%coordinate)
         !Finally update the velocities using DV
         sy%velocity(:,:) = sy%velocity(:,:) + DV(:,:)
         call gpmdcov_msMem("gpmdcov_mdloop", "After Langevin",lt%verbose,myRank)
#ifdef DO_MPI
         if (numRanks .gt. 1) then  !THIS IS VERY IMPORTANT
            call prg_sumRealReduceN(sy%coordinate(1,:), sy%nats)
            call prg_sumRealReduceN(sy%coordinate(2,:), sy%nats)
            call prg_sumRealReduceN(sy%coordinate(3,:), sy%nats)
            
            call prg_sumRealReduceN(sy%velocity(1,:), sy%nats)
            call prg_sumRealReduceN(sy%velocity(2,:), sy%nats)
            call prg_sumRealReduceN(sy%velocity(3,:), sy%nats)
            
            sy%coordinate = sy%coordinate/real(numRanks,dp)
            sy%velocity = sy%velocity/real(numRanks,dp)
         endif
#endif

      else
         call halfVerlet(sy%mass,sy%force,lt%timestep,sy%velocity(1,:),sy%velocity(2,:),sy%velocity(3,:))
         call gpmdcov_msMem("gpmdcov_mdloop", "After halfVerlet",lt%verbose,myRank)
      endif

      if(gpmdt%writetraj .and. myRank == 1 .and. mdstep.ge.gpmdt%minimization_steps)then
        if(gpmdt%traj_format .eq. "XYZ")then
           call prg_write_trajectory(sy,mdstep-gpmdt%minimization_steps,gpmdt%writetreach,&
             &lt%timestep,adjustl(trim(lt%jobname))//"_trajectory","xyz")
        else
           call prg_write_trajectory(sy,mdstep-gpmdt%minimization_steps,gpmdt%writetreach,&
             &lt%timestep,adjustl(trim(lt%jobname))//"_trajectory","pdb")
        endif
      endif
      if(gpmdt%dumpeach .gt. 0)then
         if(mod(mdstep,gpmdt%dumpeach) .eq. 0)then
            call gpmdcov_dump()
         endif
      endif
      call gpmdcov_msI("gpmdcov_MDloop","Time for rest, inc halfVerlet &
           &"//to_string(mls() - mls_md1)//" ms",lt%verbose,myRank)

      call gpmdcov_msI("gpmdcov_MDloop","Time for MD iter &
           &"//to_string(mls() - mls_md)//" ms",lt%verbose,myRank)
#ifdef USE_NVTX
      call nvtxEndRange
#endif
      
      ! Save MD state each 120 steps
      if(gpmdt%dumpeach .gt. 0)then
         if(mod(mdstep,gpmdt%dumpeach) == 0)call gpmdcov_dump()
      endif
    enddo
    ! End of MD loop.

  end subroutine gpmdcov_MDloop

end module gpmdcov_MDloop_mod
