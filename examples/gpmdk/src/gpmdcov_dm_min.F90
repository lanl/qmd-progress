module gpmdcov_DM_Min_mod
    use gpmdcov_nonequilibrium_mod

#ifdef USE_NVTX
    use gpmdcov_nvtx_mod
#endif

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
    use gpmdcov_field

    integer, intent(in) :: Nr_SCF
    real(dp), allocatable :: nguess(:)
    logical, intent(in) :: mix
    real(dp) :: tch1, mls_coul

    converged = .false.
    charges_old = nguess

    call gpmdcov_msMem("gpmdcov_DM_Min","Berofe gpmd_DM_Min",lt%verbose,myRank)

    if(.not.allocated(auxcharge))allocate(auxcharge(sy%nats))

    ! Beginning of the SCF loop.
    do iscf=1,Nr_SCF

      if (myRank == 1) write(*,*)"SCF iter", iscf

      call gpmdcov_msI("gpmdcov_DM_Min","rank "//to_string(myRank)//" SCF iter "//to_string(iscf),lt%verbose,myRank)

      !> Real contribution to the Coulomb energy. The outputs are coul_forces_r,coul_pot_r.
      call gpmdcov_msII("gpmdcov_DM_Min","In real Coul ...",lt%verbose,myRank)
      if(myRank == 1 .and. lt%verbose >= 1) mls_coul = mls()
      !       call get_ewald_list_real(sy%spindex,sy%splist,sy%coordinate&
      !         ,nguess,tb%hubbardu,sy%lattice_vector,&
      !         sy%volr,lt%coul_acc,lt%timeratio,nl%nnRx,nl%nnRy,&
      !         nl%nnRz,nl%nrnnlist,nl%nnType,coul_forces_r,coul_pot_r);

      call get_ewald_list_real_dcalc_vect(sy%spindex,sy%splist,sy%coordinate&
           ,nguess,tb%hubbardu,sy%lattice_vector,&
           sy%volr,lt%coul_acc,lt%timeratio,nl%nnIx,nl%nnIy,&
           nl%nnIz,nl%nrnnlist,nl%nnType,coul_forces_r,coul_pot_r);
      call gpmdcov_msII("gpmdcov_DM_Min","Time real coul "//to_string(mls() - mls_coul)//" ms",lt%verbose,myRank)
      if(myRank == 1 .and. lt%verbose >= 1) write(*,*)"Time for real coul",mls() - mls_coul


      !> Reciprocal contribution to the Coul energy. The outputs are coul_forces_k,coul_pot_k.
      call gpmdcov_msII("gpmdcov_DM_Min","In recip Coul ...",lt%verbose,myRank)
      if(myRank == 1 .and. lt%verbose >= 1) mls_coul = mls()
      call get_ewald_recip(sy%spindex,sy%splist,sy%coordinate&
           ,nguess,tb%hubbardu,sy%lattice_vector,&
           sy%recip_vector,sy%volr,lt%coul_acc,coul_forces_k,coul_pot_k);
      call gpmdcov_msII("gpmdcov_DM_Min","Time recip coul "//to_string(mls() - mls_coul)//" ms",lt%verbose,myRank)

      if(iscf == Nr_SCF) converged = .true.

!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Loop over parts
!!!!!!!!!!!!!!!!!!!!!!!!!!!

      auxcharge = 0.0_dp
      mls_i = mls()
#ifdef DO_MPI
     !do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)
      do iptt=1,partsInEachRank(myRank)
        ipt= reshuffle(iptt,myRank)
#else
      do ipt = 1,gpat%TotalParts
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
        call gpmdcov_msIII("gpmdcov_DM_Min","In prg_get_hscf ...",lt%verbose,myRank)
        call prg_get_hscf(syprt(ipt)%estr%ham0,syprt(ipt)%estr%over,syprt(ipt)%estr%ham,syprt(ipt)%spindex,&
             syprt(ipt)%estr%hindex,tb%hubbardu,syprt(ipt)%net_charge,&
             syprt(ipt)%estr%coul_pot_r,syprt(ipt)%estr%coul_pot_k,lt%mdim,lt%threshold)

        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%oham)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%orho)

        !> Orthogonalize ham.
        call prg_orthogonalize(syprt(ipt)%estr%ham,syprt(ipt)%estr%zmat,syprt(ipt)%estr%oham,&
             lt%threshold,lt%bml_type,0)

        !> Now solve for the desity matrix.
        call gpmdcov_RhoSolver(syprt(ipt)%estr%oham,syprt(ipt)%estr%orho,syprt(ipt)%estr%evects)
        norbsInEachCHAtRank(iptt) = size(syprt(ipt)%estr%evals,dim=1)

        !call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,rho_bml)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%rho)

        !> Deprg_orthogonalize orthop_bml to get the density matrix rho_bml.
        call gpmdcov_msIII("gpmdcov_DM_Min","Entering prg_deorthogonalize...",lt%verbose,myRank)
        call prg_deorthogonalize(syprt(ipt)%estr%orho,syprt(ipt)%estr%zmat,syprt(ipt)%estr%rho,&
             lt%threshold,lt%bml_type,0)
        call gpmdcov_msIII("gpmdcov_DM_Min","Leaving prg_deorthogonalize...",lt%verbose,myRank)

        !> Get the system charges from rho
        call gpmdcov_msIII("gpmdcov_DM_Min","Getting system charges...",lt%verbose,myRank)
        call prg_get_charges(syprt(ipt)%estr%rho,syprt(ipt)%estr%over,syprt(ipt)%estr%hindex,syprt(ipt)%net_charge,tb%numel,&
             syprt(ipt)%spindex,lt%mdim,lt%threshold)
        call gpmdcov_msIII("gpmdcov_DM_Min","Leaving Get charges...",lt%verbose,myRank)

        call gpmdcov_msIII("gpmdcov_DM_Min","Total charge of the part for Rank "//to_string(printRank())// &
             & "="//to_string(sum(syprt(ipt)%net_charge(:)))// &
             & "("//to_string(size(syprt(ipt)%net_charge,dim=1))//")",lt%verbose,myRank)

        call gpmdcov_msIII("gpmdcov_DM_Min","Total charge of the core part for Rank "//to_string(printRank())//&
             & "="//to_string(sum(syprt(ipt)%net_charge(1:gpat%sgraph(ipt)%llsize))),lt%verbose,myRank)

        do ii=1,gpat%sgraph(ipt)%llsize
          j = gpat%sgraph(ipt)%core_halo_index(ii)+1
          auxcharge(j) = syprt(ipt)%net_charge(ii)
        enddo

      enddo

      call gpmdcov_msII("gpmdcov_DM_Min","Time for get qs of all parts "//to_string(mls() - mls_i)//" ms",lt%verbose,myRank)

      call prg_wait()

      mls_i = mls()


#ifdef DO_MPI
      if (getNRanks() > 1) then
        call prg_sumRealReduceN(auxcharge(:), sy%nats)
      endif
#endif
      call gpmdcov_msII("gpmdcov_DM_Min","MPI rank finished prg_sumRealReduceN for qs "&
              &//to_string(mls() - mls_i)//" ms",lt%verbose,myRank)

      nguess = auxcharge

      scferror = norm2(nguess(:)-charges_old(:))/sqrt(real(sy%nats,dp))
      if(mix)then
        if(kernel%kernelMixing .and. lt%dokernel)then
          if( scferror < 0.1_dp)then
            if(firstKernel)then
              call gpmdcov_getKernel(sy%nats)
              KSum = 0.0_dp
            endif
            nguess = charges_old - MATMUL(Ker,(nguess-charges_old))
            call gpmdcov_msI("gpmdcov_DM_Min","SCF: Doing Kernel assisted mixing",lt%verbose,myRank)
          else
            nguess = charges_old + 0.1_dp*(nguess-charges_old)
            call gpmdcov_msI("gpmdcov_DM_Min","SCF: Doing Kernel assisted mixing (linear part)",lt%verbose,myRank)
          endif
        else
          call prg_qmixer(nguess,charges_old,dqin,&
               dqout,scferror,iscf,lt%pulaycoeff,lt%mpulay,0)
          call gpmdcov_msI("gpmdcov_DM_Min","SCF: Doing Pulay/DIIS mixing ",lt%verbose,myRank)
        endif
      else
        call prg_linearmixer(nguess,charges_old,scferror,lt%mixcoeff,lt%verbose)
      endif

      if(myRank == 1) write(*,*)"gpmdcov_DM_Min","SCF Error "//to_string(scferror)

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
      !         call gpmdcov_muFromParts()
      !         call gpmdcov_muDyn(nguess,Nr_SCF)
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

    call gpmdcov_msI("gpmdcov_dm_min", "Total charge ="//to_string(tch),lt%verbose,myRank)

    !> End of SCF loop.
    if (lt%verbose == 6 .and. converged)then
      ipt = 1
      bndfil = real(0.5*nel/norb)
      allocate(eigenvalues(norb))
      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,aux_bml)
      call bml_write_matrix(syprt(ipt)%estr%ham,"H.mtx")

      call prg_build_density_T0(syprt(ipt)%estr%oham, aux_bml, lt%threshold, bndfil, eigenvalues)
      call bml_write_matrix(aux_bml,"D.mtx")

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
      write(*,*)"Total charge of the system (After SCF) =",sum(nguess(:),dim=1),size(nguess(:),dim=1)
      write(*,*)""
      call prg_write_system(sy,"charged_system","pdb")
      write(*,*)""; write(*,*)"System charges (After SCF):"
      do j=1,sy%nats
        write(*,*)sy%symbol(j),nguess(j)
      enddo
    endif

    call gpmdcov_msMem("gpmdcov_dm_min", "After gpmd_DM_Min",lt%verbose,myRank)

  end subroutine gpmdcov_DM_Min


  subroutine gpmdcov_DM_Min_Eig(Nr_SCF,nguess,mix,applyField)

    use gpmdcov_vars
    use gpmdcov_mod
    use gpmdcov_rhosolver_mod
    use gpmdcov_writeout_mod
    use gpmdcov_kernel_mod
    use gpmdcov_Diagonalize_mod

    integer, intent(in) :: Nr_SCF
    real(dp), allocatable :: nguess(:),kernelTimesRes(:)
    logical, intent(in) :: mix,applyField
    real(dp) :: tch1
    real(8) :: mls_v, mls_coul, mls_mu, mls_red, mls_mix, mls_scfIter
    real(dp), allocatable :: KK0Res(:)

    converged = .false.
    if(.not.allocated(charges_old))allocate(charges_old(sy%nats))
    charges_old = nguess

    call gpmdcov_msMem("gpmdcov_dm_min_eig","Before gpmd_DM_Min_Eig",lt%verbose,myRank)

    ! Beginning of the SCF loop.
    if(.not.allocated(auxcharge))allocate(auxcharge(sy%nats))

    do iscf=1,Nr_SCF
      
      mls_scfIter = mls()
      
      call gpmdcov_msI("gpmdcov_DM_Min","rank "//to_string(myRank)//" SCF iter"//to_string(iscf),lt%verbose,myRank)

      !> Real contribution to the Coul energy. The outputs are
      !coul_forces_r,coul_pot_r.
      call gpmdcov_msI("gpmdcov_DM_Min","In real Coul ...",lt%verbose,myRank)
#ifdef DO_MPI
#ifdef USE_NVTX
        call nvtxStartRange("BarrierBeforeEwald",2)
        call prg_barrierParallel
        call nvtxEndRange
#endif
#endif

    call gpmdcov_msMem("gpmdcov_dm_min_eig","Before get_ewald_list_real_dcalc_vect",lt%verbose,myRank)
    if(myRank == 1 .and. lt%verbose >= 1) mls_coul = mls()
      call get_ewald_list_real_dcalc_vect(sy%spindex,sy%splist,sy%coordinate&
           ,nguess,tb%hubbardu,sy%lattice_vector,&
           sy%volr,lt%coul_acc,lt%timeratio,nl%nnIx,nl%nnIy,&
           nl%nnIz,nl%nrnnlist,nl%nnType,coul_forces_r,coul_pot_r);
    call gpmdcov_msII("gpmdcov_DM_Min","Time real coul "//to_string(mls() - mls_coul)//" ms",lt%verbose,myRank)
    call gpmdcov_msMem("gpmdcov_dm_min_eig","After get_ewald_list_real_dcalc_vect",lt%verbose,myRank)

      !> Reciprocal contribution to the Coul energy. The outputs are
      !coul_forces_k,coul_pot_k.
      call gpmdcov_msI("gpmdcov_DM_Min","In recip Coul ...",lt%verbose,myRank)

    call gpmdcov_msMem("gpmdcov_dm_min_eig","Before get_ewald_recip",lt%verbose,myRank)
    if(myRank == 1 .and. lt%verbose >= 1) mls_coul = mls()
      call get_ewald_recip(sy%spindex,sy%splist,sy%coordinate&
           &,nguess,tb%hubbardu,sy%lattice_vector,&
           &sy%recip_vector,sy%volr,lt%coul_acc,coul_forces_k,coul_pot_k);
    call gpmdcov_msII("gpmdcov_DM_Min","Time recip coul "//to_string(mls() - mls_coul)//" ms",lt%verbose,myRank)
    call gpmdcov_msMem("gpmdcov_dm_min_eig","After get_ewald_recip",lt%verbose,myRank)

#ifdef DO_MPI
#ifdef USE_NVTX
        call nvtxStartRange("BarrierAfterEwald",3)
        call prg_barrierParallel
        call nvtxEndRange
#endif
#endif
      if(iscf == Nr_SCF) converged = .true.

      call gpmdcov_msMem("gpmdcov_dm_min_eig", "Before gpmd_diagonalize_H1",lt%verbose,myRank)
      call gpmdcov_diagonalize_H1(nguess)
      call gpmdcov_msMem("gpmdcov_dm_min_eig", "After gpmd_diagonalize_H1",lt%verbose,myRank)

      if(lt%MuCalcType == "FromParts" .or. lt%MuCalcType == "Combined")then 
      call gpmdcov_msMem("gpmdcov_dm_min_eig", "Before gpmdcov_muFromParts",lt%verbose,myRank)

#ifdef DO_MPI
#ifdef USE_NVTX
        call nvtxStartRange("BarrierAfterH1",4)
        call prg_barrierParallel
        call nvtxEndRange
#endif
#endif
      if(myRank == 1 .and. lt%verbose >= 1) mls_mu = mls()
        call gpmdcov_muFromParts()
      call gpmdcov_msII("gpmdcov_DM_Min","Time for get Mu "//to_string(mls() - mls_mu)//" ms",lt%verbose,myRank)
      call gpmdcov_msMem("gpmdcov_dm_min_eig", "After gpmdcov_muFromParts",lt%verbose,myRank)
      else
        call gpmdcov_msI("gpmdcov_getmu","No Mu Calculation method. I will use &
             & a fixed mu instead ...",lt%verbose,myRank)
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Loop over parts
!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call gpmdcov_msMem("gpmdcov_dm_min_eig", "Before Loop over parts",lt%verbose,myRank)
      auxcharge = 0.0_dp
      mls_i = mls()
#ifdef DO_MPI
      !    !do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)
      do iptt=1,partsInEachRank(myRank)
        ipt= reshuffle(iptt,myRank)
#else
      do iptt = 1,gpat%TotalParts
        ipt = iptt
#endif
        norb = syprt(ipt)%estr%norbs

        if(bml_allocated(syprt(ipt)%estr%orho))call bml_deallocate(syprt(ipt)%estr%orho)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%orho)

        mls_v = mls()
        call prg_build_density_fromEvalsAndEvects(syprt(ipt)%estr%evects, syprt(ipt)%estr%evals,&
             & syprt(ipt)%estr%orho, lt%threshold, bndfil, lt%kbt, Ef, lt%verbose)
        call gpmdcov_msI("prg_build_density_fromEvalsAndEvects","Time build density "&
                &//to_string(mls() - mls_v)//" ms",lt%verbose,myRank)

        norbsInEachCHAtRank(iptt) = size(syprt(ipt)%estr%evals,dim=1)

        if(bml_allocated(syprt(ipt)%estr%rho))call bml_deallocate(syprt(ipt)%estr%rho)
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%rho)

        !> Deprg_orthogonalize orthop_bml to get the density matrix rho_bml.
        mls_v = mls()
        call gpmdcov_msIII("gpmdcov_DM_Min","Entering prg_deorthogonalize...",lt%verbose,myRank)
        call prg_deorthogonalize(syprt(ipt)%estr%orho,syprt(ipt)%estr%zmat,syprt(ipt)%estr%rho,&
             lt%threshold,lt%bml_type,0)
        call gpmdcov_msIII("gpmdcov_DM_Min","Leaving prg_deorthogonalize...",lt%verbose,myRank)
        call gpmdcov_msI("prg_build_density_fromEvalsAndEvects","Time for deorthogonalize "&
                &//to_string(mls() - mls_v)//" ms",lt%verbose,myRank)

        !> Get the system charges from rho
        call gpmdcov_msIII("gpmdcov_DM_Min","Getting system charges...",lt%verbose,myRank)
        mls_v = mls()
        call prg_get_charges(syprt(ipt)%estr%rho,syprt(ipt)%estr%over,syprt(ipt)%estr%hindex,syprt(ipt)%net_charge,tb%numel,&
             &syprt(ipt)%spindex,lt%mdim,lt%threshold)
        call gpmdcov_msIII("gpmdcov_DM_Min","Leaving Get charges...",lt%verbose,myRank)
        call gpmdcov_msI("prg_build_density_fromEvalsAndEvects","Time for get charges "&
                &//to_string(mls() - mls_v)//" ms",lt%verbose,myRank)

        call gpmdcov_msIII("gpmdcov_DM_Min","Total charge of the part for Rank "//to_string(printRank())// &
             & "="//to_string(sum(syprt(ipt)%net_charge(:)))// &
             & "("//to_string(size(syprt(ipt)%net_charge,dim=1))//")",lt%verbose,myRank)

        call gpmdcov_msIII("gpmdcov_DM_Min","Total charge of the core part for Rank "//to_string(printRank())//&
             & "="//to_string(sum(syprt(ipt)%net_charge(1:gpat%sgraph(ipt)%llsize))),lt%verbose,myRank)

        do ii=1,gpat%sgraph(ipt)%llsize
          j = gpat%sgraph(ipt)%core_halo_index(ii)+1
          auxcharge(j) = syprt(ipt)%net_charge(ii)
        enddo

        !> Get the system currents from rho and ham 
        if(gpmdt%compcurr)then
          if(MDStep >= 2 .and. (mod(MDstep,gsp2%parteach) .ne. 0))then
            norb_core = syprt(ipt)%estr%hindex(2,gpat%sgraph(ipt)%llsize)
            call gpmdcov_get_currents(syprt(ipt)%estr%ham,oldsyprt(ipt)%estr%rho,syprt(ipt)%estr%zmat,&
            &syprt(ipt)%estr%hindex,gpat%sgraph(ipt)%llsize,norb_core,gpat%sgraph(ipt)%core_halo_index,sy%symbol,&
            &gpmdt%currthr)
          endif
            
          !Storing the old subsystem density matrix
          if(bml_allocated(oldsyprt(ipt)%estr%rho))then
            call bml_deallocate(oldsyprt(ipt)%estr%rho)
          endif
          call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,oldsyprt(ipt)%estr%rho)
          call bml_copy(syprt(ipt)%estr%rho, oldsyprt(ipt)%estr%rho)
        endif

      enddo

      mls_red = mls()
#ifdef DO_MPI
      if (getNRanks() > 1) then
        call prg_sumRealReduceN(auxcharge(:), sy%nats)
      endif
#endif
      call gpmdcov_msMem("gpmdcov_dm_min_eig", "After Loop over parts",lt%verbose,myRank)
      call gpmdcov_msIII("gpmdcov_DM_Min","Time for prg_sumRealReduceN&
           &for qs "//to_string(mls() - mls_red)//" ms",lt%verbose,myRank)
      call gpmdcov_msMem("gpmdcov_dm_min_eig", "Before Kernel logic",lt%verbose,myRank)
      nguess = auxcharge
      scferror = 0.0d0
      scferror = norm2(nguess(:)-charges_old(:))/sqrt(real(sy%nats,dp))
!      if(myRank == 1) write(*,*)"ResVect", nguess(1:5)-charges_old(1:5)
!      if(myRank == 1) write(*,*)"nguess", nguess(1:5)
!      if(myRank == 1) write(*,*)"charges_old", charges_old(1:5)
!      if(myRank == 1) write(*,*)"SCF Error",iscf,scferror
!      if(myRank == 1) write(*,*)"mdstep,FirstKernel,NewPart,scferror,mix",mdstep,firstKernel,newPart,scferror,mix

      mls_mix = mls()
      if(mix)then
!        if(myRank == 1) write(*,*)mdstep,"InIf (mix)"
        if( kernel%kernelMixing .and. lt%dokernel)then
!          if(myRank == 1) write(*,*)mdstep,"InIf ( kernel%kernelMixing .and. lt%dokernel)"
          if( scferror < 0.01_dp .or. newPart)then
!            if(myRank == 1) write(*,*)mdstep,"InIf ( scferror < 0.01_dp .or. newPart)"
            if(firstKernel .or. kernel%buildAlways)then
!              if(myRank == 1) write(*,*)mdstep,"InIf (firstKernel .or. kernel%buildAlways)"
              if(kernel%kernelType == "Full")then
                call gpmdcov_getKernel(sy%nats)
              elseif(kernel%kernelType == "ByBlocks")then
                call gpmdcov_getKernel_byBlocks(sy%nats)
              elseif(kernel%kernelType == "ByParts")then
                call gpmdcov_getKernel_byParts(syprt,syprtk)
              else
                write(*,*)"The TypeOfKernel is not implemented"
                stop
              endif
              
              firstKernel = .false.
              newPart = .false.
              if((kernel%kernelType == "ByParts" ) .and. (.not. kernel%updateAfterBuild))then
                allocate(kernelTimesRes(sy%nats))
!                if(myRank == 1) write(*,*)mdstep,"InIf ( scferror < 0.01_dp .or. newPart), Applying K"
                !CHANGE Check with Anders
                !if(mdstep <= xl%minit)then 
                        call gpmdcov_applyKernel(nguess,charges_old,syprtk,kernelTimesRes)
                        nguess = charges_old - kernelTimesRes
                !endif
                !        call gpmdcov_rankN_update_byParts(nguess,charges_old,syprt,syprtk,kernel%rankNUpdate,KK0Res)
                !       nguess = charges_old - KK0Res
                       !nguess = nguess
                !endif

                deallocate(kernelTimesRes)
              elseif(kernel%kernelType == "ByParts" .and. kernel%updateAfterBuild)then 
                call gpmdcov_rankN_update_byParts(nguess,charges_old,syprt,syprtk,kernel%rankNUpdate,KK0Res)
                nguess = charges_old - KK0Res
              else
                nguess = charges_old - MATMUL(Ker,(nguess-charges_old))
              endif
            else
              !Update and apply
              if(myRank == 1) write(*,*)mdstep,"InIf elseOf( scferror < 0.01_dp .or. newPart)"
              if(kernel%rankNUpdate > 0)then 
                if(kernel%kernelType == "Full")then
                        write(*,*)"The Kernel rank N updates for the KernelType=Full is not implemented"
                        stop
                elseif(kernel%kernelType == "ByParts")then
                        call gpmdcov_rankN_update_byParts(nguess,charges_old,syprt,syprtk,kernel%rankNUpdate,KK0Res)
                        nguess = charges_old - KK0Res
                endif
              else !Just apply the old one
                if(kernel%kernelType == "Full")then
                        nguess = charges_old - MATMUL(Ker,(nguess-charges_old))
                elseif(kernel%kernelType == "ByParts")then
                        allocate(kernelTimesRes(sy%nats))
                        call gpmdcov_applyKernel(nguess,charges_old,syprtk,kernelTimesRes)
                        nguess = charges_old - kernelTimesRes
                        deallocate(kernelTimesRes)
                endif
              endif
              
            endif
            call gpmdcov_msI("gpmdcov_DM_Min","SCF: Doing Kernel assisted mixing",lt%verbose,myRank)
          else !if( scferror < 0.01_dp .or. newPart)
            if(kernel%initialMixingWith == "Lin")then
                nguess = charges_old + lt%pulaycoeff*(nguess-charges_old)
                call gpmdcov_msI("gpmdcov_DM_Min","SCF: Doing Initial Kernel mixing using linear mixing",lt%verbose,myRank)
            else
                call prg_qmixer(nguess,charges_old,dqin,&
                dqout,scferror,iscf,lt%pulaycoeff,lt%mpulay,0)
                call gpmdcov_msI("gpmdcov_DM_Min","SCF: Doing Initial Kernel mixing using Pulay/DIIS mixing",lt%verbose,myRank)
            endif
          endif
        else !if( kernel%kernelMixing .and. lt%dokernel)
          call prg_qmixer(nguess,charges_old,dqin,&
               dqout,scferror,iscf,lt%pulaycoeff,lt%mpulay,0)
          call gpmdcov_msI("gpmdcov_DM_Min","SCF: Doing Pulay/DIIS mixing",lt%verbose,myRank)
        endif
      else !if(mix)
        if(lt%dokernel .and. firstKernel)then
              if(kernel%kernelType == "Full")then
                call gpmdcov_getKernel(sy%nats)
              elseif(kernel%kernelType == "ByBlocks")then
                call gpmdcov_getKernel_byBlocks(sy%nats)
              elseif(kernel%kernelType == "ByParts")then
                call gpmdcov_getKernel_byParts(syprt,syprtk)
              else
                write(*,*)"The TypeOfKernel is not implemented"
                stop
              endif

              firstKernel = .false.
              !CHANGE
             ! if(kernel%kernelType == "ByParts")then
             !   allocate(kernel imesRes(sy%nats))
             !   call gpmdcov_applyKernel(nguess,charges_old,syprtk,kernelTimesRes)
             !   nguess = charges_old - kernelTimesRes
             !   deallocate(kernelTimesRes)
             ! else
             !   nguess = charges_old - MATMUL(Ker,(nguess-charges_old))
             ! endif
              
        else
              !CHANGE
              !call prg_linearmixer(nguess,charges_old,scferror,lt%mixcoeff,lt%verbose)
              nguess = nguess
        endif      
  
      endif
      call gpmdcov_msI("gpmdcov_DM_Min","Time for Mix &
           & "//to_string(mls() - mls_mix)//" ms",lt%verbose,myRank)

      call gpmdcov_msMem("gpmdcov_dm_min_eig", "After Kernel logic",lt%verbose,myRank)

      if(myRank == 1) write(*,*)"gpmdcov_DM_Min","SCF Error "//to_string(scferror)

      tch1 = sum(charges_old)
      charges_old = nguess
      tch = sum(nguess)
      call prg_wait
      if(lt%MuCalcType == "Combined" .or. lt%MuCalcType == "Dyn")then
        call gpmdcov_muDyn(nguess,Nr_SCF)
      endif

      call gpmdcov_msII("gpmdcov_dm_min", "Total charge="//to_string(tch),lt%verbose,myRank)
      call gpmdcov_msI("gpmdcov_DM_Min","Time for one SCF iter "&
           &//to_string(mls() - mls_scfIter)//" ms",lt%verbose,myRank)

      if(converged)then ! To do a last extra step.
      !  if(lt%dokernel)then
      !    if(kernel%kernelType == "Full")then
      !      call gpmdcov_getKernel(sy%nats)
      !    elseif(kernel%kernelType == "ByBlocks")then
      !      call gpmdcov_getKernel_byBlocks(sy%nats)
      !    elseif(kernel%kernelType == "ByParts")then
      !      call gpmdcov_getKernel_byParts(sy%nats,syprt,syprtk)
      !    endif
      !  endif

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
    newPart = .false.
    deallocate(auxcharge)

    call gpmdcov_msMem("gpmdcov_dm_min", "After gpmd_DM_Min",lt%verbose,myRank)

  end subroutine gpmdcov_DM_Min_Eig


end module gpmdcov_DM_Min_mod
