module gpmdcov_DM_Min_mod

contains

  !>  SCF loop
  !!
  subroutine gpmdcov_DM_Min(Nr_SCF,nguess,mix)

    use gpmdcov_vars
    use gpmdcov_mod
    use gpmdcov_rhosolver_mod
    use gpmdcov_writeout_mod
    use gpmdcov_kernel_mod
    use gpmdcov_diagonalize_mod

    integer, intent(in) :: Nr_SCF
    real(dp), allocatable :: nguess(:),residues(:)
    logical, intent(in) :: mix
    logical :: err
    real(dp) :: tch1

    converged = .false.
    charges_old = nguess
    
    call gpmdcov_msMem("gpmdcov_DM_Min","Before gpmd_DM_Min",lt%verbose,myRank)

    ! Beginning of the SCF loop.
    if(.not.allocated(auxcharge))allocate(auxcharge(sy%nats))

    !call gpmdcov_getKernel(sy%nats)
    !stop

    do iscf=1,Nr_SCF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! This is done for the whole system
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (myRank == 1) write(*,*)"SCF iter", iscf

      call gpmdcov_msI("gpmdcov_DM_Min","rank "//to_string(myRank)//" SCF iter "//to_string(iscf),lt%verbose,myRank)

      !> Real contribution to the Coul energy. The outputs are coul_forces_r,coul_pot_r.
      call gpmdcov_msI("gpmdcov_DM_Min","In real Coul ...",lt%verbose,myRank)

      if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_start(dyn_timer,"Real coul")
      !       call get_ewald_list_real(sy%spindex,sy%splist,sy%coordinate&
      !         ,nguess,tb%hubbardu,sy%lattice_vector,&
      !         sy%volr,lt%coul_acc,lt%timeratio,nl%nnRx,nl%nnRy,&
      !         nl%nnRz,nl%nrnnlist,nl%nnType,coul_forces_r,coul_pot_r);
      call get_ewald_list_real_dcalc(sy%spindex,sy%splist,sy%coordinate&
           ,nguess,tb%hubbardu,sy%lattice_vector,&
           sy%volr,lt%coul_acc,lt%timeratio,nl%nnIx,nl%nnIy,&
           nl%nnIz,nl%nrnnlist,nl%nnType,coul_forces_r,coul_pot_r);
      if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_stop(dyn_timer,1)


      !> Reciprocal contribution to the Coul energy. The outputs are coul_forces_k,coul_pot_k.
      call gpmdcov_msI("gpmdcov_DM_Min","In recip Coul ...",lt%verbose,myRank)

      if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_start(dyn_timer,"Recip coul")
      call get_ewald_recip(sy%spindex,sy%splist,sy%coordinate&
           ,nguess,tb%hubbardu,sy%lattice_vector,&
           sy%recip_vector,sy%volr,lt%coul_acc,coul_forces_k,coul_pot_k);
      if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_stop(dyn_timer,1)

      if(iscf == Nr_SCF) converged = .true.

!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Loop over parts
!!!!!!!!!!!!!!!!!!!!!!!!!!!

      auxcharge = 0.0_dp
      mls_i = mls()
#ifdef DO_MPI
      !    !do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)
      do iptt=1,partsInEachRank(myRank)
        ipt= reshuffle(iptt,myRank)
#else
        !     do ipt = 1,gpat%TotalParts
#endif
        norb = syprt(ipt)%estr%norbs

        if(.not.allocated(syprt(ipt)%estr%coul_pot_k))then
          allocate(syprt(ipt)%estr%coul_pot_k(syprt(ipt)%nats))
          allocate(syprt(ipt)%estr%coul_pot_r(syprt(ipt)%nats))
        endif

        syprt(ipt)%estr%coul_pot_k = 0.0_dp
        syprt(ipt)%estr%coul_pot_r = 0.0_dp
        syprt(ipt)%net_charge = 0.0_dp

        !> Get Coulombic potential and charges for the part.
        do j=1,gpat%sgraph(ipt)%lsize
          jj = gpat%sgraph(ipt)%core_halo_index(j)+1
          syprt(ipt)%estr%coul_pot_k(j) = coul_pot_k(jj)
          syprt(ipt)%estr%coul_pot_r(j) = coul_pot_r(jj)
          syprt(ipt)%net_charge(j) = nguess(jj)
        enddo

        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%ham)
        call bml_copy_new(syprt(ipt)%estr%ham0,syprt(ipt)%estr%ham)

        !> Get the scf hamiltonian. The output is ham_bml.
        call gpmdcov_msI("gpmdcov_DM_Min","In prg_get_hscf ...",lt%verbose,myRank)
        call prg_get_hscf(syprt(ipt)%estr%ham0,syprt(ipt)%estr%over,syprt(ipt)%estr%ham,syprt(ipt)%spindex,&
             syprt(ipt)%estr%hindex,tb%hubbardu,syprt(ipt)%net_charge,&
             syprt(ipt)%estr%coul_pot_r,syprt(ipt)%estr%coul_pot_k,lt%mdim,lt%threshold)

        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%oham)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%orho)

        !> Orthogonalize ham.
        if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_start(ortho_timer)
        call prg_orthogonalize(syprt(ipt)%estr%ham,syprt(ipt)%estr%zmat,syprt(ipt)%estr%oham,&
             lt%threshold,lt%bml_type,0)
        if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_stop(ortho_timer)

        !> Now solve for the desity matrix.
        call gpmdcov_RhoSolver(syprt(ipt)%estr%oham,syprt(ipt)%estr%orho,syprt(ipt)%estr%evects)
        norbsInEachCHAtRank(iptt) = size(syprt(ipt)%estr%evals,dim=1)

        !call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,rho_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%rho)

        !> Deprg_orthogonalize orthop_bml to get the density matrix rho_bml.
        if(printRank() == 1 .and. lt%verbose >= 1) call prg_timer_start(deortho_timer)
        call gpmdcov_msI("gpmdcov_DM_Min","Entering prg_deorthogonalize...",lt%verbose,myRank)
        call prg_deorthogonalize(syprt(ipt)%estr%orho,syprt(ipt)%estr%zmat,syprt(ipt)%estr%rho,&
             lt%threshold,lt%bml_type,0)
        call gpmdcov_msI("gpmdcov_DM_Min","Leaving prg_deorthogonalize...",lt%verbose,myRank)
        if(printRank() == 1 .and. lt%verbose >= 1) call prg_timer_stop(deortho_timer)

        !> Get the system charges from rho
        call gpmdcov_msI("gpmdcov_DM_Min","Getting system charges...",lt%verbose,myRank)
        call prg_get_charges(syprt(ipt)%estr%rho,syprt(ipt)%estr%over,syprt(ipt)%estr%hindex,syprt(ipt)%net_charge,tb%numel,&
             syprt(ipt)%spindex,lt%mdim,lt%threshold)
        call gpmdcov_msI("gpmdcov_DM_Min","Leaving Get charges...",lt%verbose,myRank)

        call gpmdcov_msIII("gpmdcov_DM_Min","Total charge of the part ="//to_string(printRank())// &
             & to_string(sum(syprt(ipt)%net_charge(:)))//to_string(size(syprt(ipt)%net_charge,dim=1)),lt%verbose,myRank)

        call gpmdcov_msIII("gpmdcov_DM_Min","Total charge of the core part = "//to_string(printRank())//&
             & to_string(sum(syprt(ipt)%net_charge(1:gpat%sgraph(ipt)%llsize))),lt%verbose,myRank)

        do ii=1,gpat%sgraph(ipt)%llsize
          j = gpat%sgraph(ipt)%core_halo_index(ii)+1
          auxcharge(j) = syprt(ipt)%net_charge(ii)
        enddo

      enddo

      call gpmdcov_msIII("gpmdcov_DM_Min","Time for get qs of all parts "//to_string(mls() - mls_i)//" ms",lt%verbose,myRank)

      mls_i = mls()


#ifdef DO_MPI
      if (getNRanks() > 1) then
        call prg_sumRealReduceN(auxcharge(:), sy%nats)
      endif
#endif
      call gpmdcov_msIII("gpmdcov_DM_Min","MPI rank finished prg_sumRealReduceN for qs "//to_string(mls() - mls_i)//" ms",lt%verbose,myRank)

      nguess = auxcharge

      scferror = norm2(nguess(:)-charges_old(:))/sqrt(real(sy%nats,dp))
      write(*,*)"Res Vect", nguess(1:5)-charges_old(1:5)
      if(mix)then
        if(kernel%kernelMixing .and. lt%dokernel)then
          if( scferror < 0.1_dp)then
            if(firstKernel)then
              call gpmdcov_getKernel(sy%nats)
              KSum = 0.0_dp
              !firstKernel = .false.
            endif
            ! scferror = norm2(nguess(:)-charges_old(:))/sqrt(real(sy%nats,dp))
            nguess = charges_old - MATMUL(Ker,(nguess-charges_old))

            call gpmdcov_msI("gpmdcov_DM_Min","SCF: Doing Kernel assisted mixing",lt%verbose,myRank)
          else
            ! scferror = norm2(nguess(:)-charges_old(:))/sqrt(real(sy%nats,dp))
            nguess = charges_old + 0.1_dp*(nguess-charges_old)
            call gpmdcov_msI("gpmdcov_DM_Min","SCF: Doing Kernel assisted mixing (linear part)",lt%verbose,myRank)
          endif
        else
          !LOOK at the res. If you have a real kernel it should behave much lower.
          call prg_qmixer(nguess,charges_old,dqin,&
               dqout,scferror,iscf,lt%pulaycoeff,lt%mpulay,0)
          call gpmdcov_msI("gpmdcov_DM_Min","SCF: Doing Pulay/DIIS mixing ",lt%verbose,myRank)
        endif
      else
        call prg_linearmixer(nguess,charges_old,scferror,lt%mixcoeff,lt%verbose)
      endif

      call gpmdcov_msI("gpmdcov_DM_Min","SCF Error "//to_string(scferror),lt%verbose,myRank)

      tch1 = sum(charges_old)
      charges_old = nguess
      tch = sum(nguess)

      call gpmdcov_msII("gpmdcov_dm_min", "Total charge ="//to_string(tch),lt%verbose,myRank)

      !Get the chemical potential. Feb 2021 implemetation
    !       if(lt%MuCalcType == "Dyn")then
    !         call gpmdcov_muDyn(nguess,Nr_SCF)
    !       elseif(lt%MuCalcType == "FromParts")then
    !         call gpmdcov_muFromParts()
    !       elseif(lt%MuCalcType == "Combined")then
    !         call gpmdcov_muDyn(nguess,Nr_SCF)
    !         call gpmdcov_muFromParts()
    !       else
    !         call gpmdcov_msI("gpmdcov_getmu","No Mu Calculation method. I will use &
    !              & a fixed mu instead ...",lt%verbose,myRank)
    !       endif

      if(converged)then ! To do a last extra step.
        exit
      else
        if(scferror < lt%scftol .and. iscf > 1) then
          if (myRank  ==  1) then
            write(*,*)""; write(*,*)"SCF converged within",iscf,"steps ..."
            write(*,*)"SCF error =",scferror
          endif
          converged = .true.
        endif
      endif
    enddo

    !KERNEL
   ! if(lt%doKernel)then
   !   if(mod(mdstep,kernel%updateEach) == 0)call gpmdcov_getKernel(sy%nats) 
   ! endif
    !END KERNEL


    call gpmdcov_msI("gpmdcov_dm_min", "Total charge ="//to_string(tch),lt%verbose,myRank)

    !> End of SCF loop.
    if (lt%verbose >= 6 .and. converged)then
      ipt = 1
      bndfil = real(0.5*nel/norb)
      allocate(eigenvalues(norb))
      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,aux_bml)
      call bml_write_matrix(syprt(ipt)%estr%ham,"H.mtx")

      call prg_build_density_T0(syprt(ipt)%estr%oham, aux_bml, lt%threshold, bndfil, eigenvalues)

      sparsity = bml_get_sparsity(aux_bml, 1.0D-5)
      write(*,*)"Sparsity Rho=",sparsity

      !Getting the fermi level
      ef = (eigenvalues(int(bndfil*norb)+1) + eigenvalues(int(bndfil*norb)))/2
      !eigenvalues = eigenvalues - ef
      write(*,*)"HOMO=",eigenvalues(int(norb/2))
      write(*,*)"LUMO=",eigenvalues(int(norb/2)+1)
      write(*,*)"Ef=", ef, int(bndfil*norb), norb
      !Writting the total DOS
      !call prg_write_tdos(eigenvalues, 0.05d0, 10000, -20.0d0, 20.0d0, "tdos.dat")

      deallocate(eigenvalues)
      call bml_deallocate(aux_bml)
      stop

    endif

    deallocate(auxcharge)

    if(myRank == 1 .and. lt%verbose >= 3)then
      write(*,*)""
      write(*,*)"Total charge of the system (After SCF) =",sum(nguess(:),size(nguess,dim=1))
      write(*,*)""
      call prg_write_system(sy,"charged_system","pdb")
      write(*,*)""; write(*,*)"System charges (After SCF):"
      do j=1,sy%nats
        write(*,*)sy%symbol(j),nguess(j)
      enddo
    endif


    call gpmdcov_msMem("gpmdcov_dm_min", "After gpmd_DM_Min",lt%verbose,myRank)

  end subroutine gpmdcov_DM_Min



  subroutine gpmdcov_DM_Min_Eig(Nr_SCF,nguess,mix)

    use gpmdcov_vars
    use gpmdcov_mod
    use gpmdcov_rhosolver_mod
    use gpmdcov_writeout_mod
    use gpmdcov_kernel_mod
    use gpmdcov_Diagonalize_mod

    integer, intent(in) :: Nr_SCF
    real(dp), allocatable :: nguess(:),residues(:),kernelTimesRes(:)
    logical, intent(in) :: mix
    logical :: err
    real(dp) :: tch1, KSum

    converged = .false.
    charges_old = nguess

    call gpmdcov_msMem("gpmdcov_dm_min_eig","Before gpmd_DM_Min_Eig",lt%verbose,myRank)

    ! Beginning of the SCF loop.
    if(.not.allocated(auxcharge))allocate(auxcharge(sy%nats))

    do iscf=1,Nr_SCF

      if (myRank == 1) write(*,*)"SCF iter", iscf

      call gpmdcov_msI("gpmdcov_DM_Min","rank "//to_string(myRank)//" SCF iter"//to_string(iscf),lt%verbose,myRank)

      !> Real contribution to the Coul energy. The outputs are
      !coul_forces_r,coul_pot_r.
      call gpmdcov_msI("gpmdcov_DM_Min","In real Coul ...",lt%verbose,myRank)

      if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_start(dyn_timer,"Real coul")
      call get_ewald_list_real_dcalc(sy%spindex,sy%splist,sy%coordinate&
           ,nguess,tb%hubbardu,sy%lattice_vector,&
           sy%volr,lt%coul_acc,lt%timeratio,nl%nnIx,nl%nnIy,&
           nl%nnIz,nl%nrnnlist,nl%nnType,coul_forces_r,coul_pot_r);
      if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_stop(dyn_timer,1)

      !> Reciprocal contribution to the Coul energy. The outputs are
      !coul_forces_k,coul_pot_k.
      call gpmdcov_msI("gpmdcov_DM_Min","In recip Coul ...",lt%verbose,myRank)

      if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_start(dyn_timer,"Recip coul")
      call get_ewald_recip(sy%spindex,sy%splist,sy%coordinate&
           &,nguess,tb%hubbardu,sy%lattice_vector,&
           &sy%recip_vector,sy%volr,lt%coul_acc,coul_forces_k,coul_pot_k);
      if(myRank == 1 .and. lt%verbose >= 1) call prg_timer_stop(dyn_timer,1)

      if(iscf == Nr_SCF) converged = .true.

      call gpmdcov_msMem("gpmdcov_dm_min_eig", "Before gpmd_diagonalize_H1",lt%verbose,myRank)
      call gpmdcov_diagonalize_H1(nguess)

      if(lt%MuCalcType == "FromParts" .or. lt%MuCalcType == "Combined")then
        call gpmdcov_muFromParts()
      endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Loop over parts
!!!!!!!!!!!!!!!!!!!!!!!!!!!

      auxcharge = 0.0_dp
      mls_i = mls()
#ifdef DO_MPI
      !    !do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)
      do iptt=1,partsInEachRank(myRank)
        ipt= reshuffle(iptt,myRank)
#else
        !     do ipt = 1,gpat%TotalParts
#endif
        norb = syprt(ipt)%estr%norbs

        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%orho)

        call prg_build_density_fromEvalsAndEvects(syprt(ipt)%estr%evects, syprt(ipt)%estr%evals,&
             & syprt(ipt)%estr%orho, lt%threshold, bndfil, lt%kbt, Ef, lt%verbose)

        norbsInEachCHAtRank(iptt) = size(syprt(ipt)%estr%evals,dim=1)

        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%rho)

        !> Deprg_orthogonalize orthop_bml to get the density matrix rho_bml.
        if(printRank() == 1 .and. lt%verbose >= 1) call prg_timer_start(deortho_timer)
        call gpmdcov_msI("gpmdcov_DM_Min","Entering prg_deorthogonalize...",lt%verbose,myRank)
        call prg_deorthogonalize(syprt(ipt)%estr%orho,syprt(ipt)%estr%zmat,syprt(ipt)%estr%rho,&
             lt%threshold,lt%bml_type,0)
        call gpmdcov_msI("gpmdcov_DM_Min","Leaving prg_deorthogonalize...",lt%verbose,myRank)
        if(printRank() == 1 .and. lt%verbose >= 1) call prg_timer_stop(deortho_timer)

        !> Get the system charges from rho
        call gpmdcov_msI("gpmdcov_DM_Min","Getting system charges...",lt%verbose,myRank)
        call prg_get_charges(syprt(ipt)%estr%rho,syprt(ipt)%estr%over,syprt(ipt)%estr%hindex,syprt(ipt)%net_charge,tb%numel,&
             &syprt(ipt)%spindex,lt%mdim,lt%threshold)
        call gpmdcov_msI("gpmdcov_DM_Min","Leaving Get charges...",lt%verbose,myRank)

        call gpmdcov_msIII("gpmdcov_DM_Min","Total charge of the part ="//to_string(printRank())// &
             & to_string(sum(syprt(ipt)%net_charge(:)))//to_string(size(syprt(ipt)%net_charge,dim=1)),lt%verbose,myRank)

        call gpmdcov_msIII("gpmdcov_DM_Min","Total charge of the core part = "//to_string(printRank())//&
             & to_string(sum(syprt(ipt)%net_charge(1:gpat%sgraph(ipt)%llsize))),lt%verbose,myRank)

        do ii=1,gpat%sgraph(ipt)%llsize
          j = gpat%sgraph(ipt)%core_halo_index(ii)+1
          auxcharge(j) = syprt(ipt)%net_charge(ii)
        enddo

      enddo

      mls_i = mls()
#ifdef DO_MPI
      if (getNRanks() > 1) then
        call prg_sumRealReduceN(auxcharge(:), sy%nats)
      endif
#endif
      call gpmdcov_msIII("gpmdcov_DM_Min","MPI rank finished prg_sumRealReduceN&
        &for qs "//to_string(mls() - mls_i)//" ms",lt%verbose,myRank)
      nguess = auxcharge

      scferror = norm2(nguess(:)-charges_old(:))/sqrt(real(sy%nats,dp))
      write(*,*)"Res Vect", nguess(1:5)-charges_old(1:5)
      if(mix)then
        if( kernel%kernelMixing .and. lt%dokernel)then
          if( scferror < 0.1_dp)then
            if(firstKernel)then
              if(kernel%kernelType == "Full")then 
                call gpmdcov_getKernel(sy%nats)
              elseif(kernel%kernelType == "ByBlocks")then
                call gpmdcov_getKernel_byBlocks(sy%nats)
              elseif(kernel%kernelType == "ByParts")then
                call gpmdcov_getKernel_byParts(sy%nats,syprt)
              else
                write(*,*)"The TypeOfKernel is not implemented"
                stop
              endif
              !firstKernel = .false.
            endif
            ! scferror = norm2(nguess(:)-charges_old(:))/sqrt(real(sy%nats,dp))
            if(kernel%kernelType == "ByParts")then  
               
               allocate(kernelTimesRes(sy%nats)) 
               call gpmdcov_applyKernel(nguess,charges_old,kernelTimesRes)
               nguess = charges_old - kernelTimesRes 
               deallocate(kernelTimesRes) 
                
            else
                nguess = charges_old - MATMUL(Ker,(nguess-charges_old))
            endif

            !nguess = charges_old - MATMUL(Ker,(nguess-charges_old))
            call gpmdcov_msI("gpmdcov_DM_Min","SCF: Doing Kernel assisted mixing",lt%verbose,myRank)
          else
            ! scferror = norm2(nguess(:)-charges_old(:))/sqrt(real(sy%nats,dp))
            nguess = charges_old + lt%pulaycoeff*(nguess-charges_old)
            call gpmdcov_msI("gpmdcov_DM_Min","SCF: Doing Kernel assisted mixing (linear part)",lt%verbose,myRank)
          endif
        else
          !LOOK at the res. If you have a real kernel it should behave much
          !lower.
          call prg_qmixer(nguess,charges_old,dqin,&
               dqout,scferror,iscf,lt%pulaycoeff,lt%mpulay,0)
          call gpmdcov_msI("gpmdcov_DM_Min","SCF: Doing Pulay/DIIS mixing",lt%verbose,myRank)
        endif
      else
        call prg_linearmixer(nguess,charges_old,scferror,lt%mixcoeff,lt%verbose)
      endif

      call gpmdcov_msI("gpmdcov_DM_Min","SCF Error = "//to_string(scferror),lt%verbose,myRank)

      tch1 = sum(charges_old)
      charges_old = nguess
      tch = sum(nguess)

      call gpmdcov_msII("gpmdcov_dm_min", "Total charge="//to_string(tch),lt%verbose,myRank)

      if(converged)then ! To do a last extra step.
        if(lt%dokernel)then
              if(kernel%kernelType == "Full")then
                call gpmdcov_getKernel(sy%nats)
              elseif(kernel%kernelType == "ByBlocks")then
                call gpmdcov_getKernel_byBlocks(sy%nats)
              elseif(kernel%kernelType == "ByParts")then
                call gpmdcov_getKernel_byParts(sy%nats,syprt)
              endif
        endif

        exit
      else
        if(scferror < lt%scftol .and. iscf > 1) then
          if (myRank  ==  1) then
            write(*,*)""; write(*,*)"SCF converged within",iscf,"steps ..."
            write(*,*)"SCF error =",scferror
          endif
          converged = .true.
        endif
      endif
    enddo
    deallocate(auxcharge)

    call gpmdcov_msMem("gpmdcov_dm_min", "After gpmd_DM_Min",lt%verbose,myRank)

  end subroutine gpmdcov_DM_Min_Eig


end module gpmdcov_DM_Min_mod
