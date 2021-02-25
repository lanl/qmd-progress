module gpmdcov_EnergAndForces_mod

  contains

  subroutine gpmdcov_EnergAndForces(charges)
    
    use gpmdcov_vars
    use gpmdcov_writeout_mod

    Implicit none
    real(dp), intent(in) :: charges(:)
    real(dp), allocatable :: ebandvector(:)

    call gpmdcov_msMem("gpmdcov","Before gpmd_EnergAndForces",lt%verbose,myRank)

    if(.not.allocated(coul_forces)) allocate(coul_forces(3,sy%nats))
    if(.not.allocated(FPUL))allocate(FPUL(3,sy%nats))
    if(.not.allocated(FSCOUL))allocate(FSCOUL(3,sy%nats))
    if(.not.allocated(SKForce))allocate(SKForce(3,sy%nats))
    if(.not.allocated(collectedforce))allocate(collectedforce(3,sy%nats))
    if(.not.allocated(ebandvector))allocate(ebandvector(gpat%TotalParts))

    FPUL = 0.0_dp
    FSCOUL = 0.0_dp
    SKForce = 0.0_dp
    collectedforce = 0.0_dp
    ebandvector = 0.0_dp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Loop over all the parts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef DO_MPI
    !   !do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)
    do iptt=1,partsInEachRank(myRank)
      ipt= reshuffle(iptt,myRank)
#else
      !      do ipt = 1,gpat%TotalParts
#endif
      !Distribute the charges back to the parts.
      do j=1,gpat%sgraph(ipt)%lsize
        jj = gpat%sgraph(ipt)%core_halo_index(j)+1
        syprt(ipt)%net_charge(j) = charges(jj)
      enddo

      norb = syprt(ipt)%estr%norbs

      norb_core = syprt(ipt)%estr%hindex(2,gpat%sgraph(ipt)%llsize)

      if(bml_get_N(aux_bml).gt.0)then
        call bml_deallocate(aux_bml)
        call bml_deallocate(aux1_bml)
        deallocate(row)
      endif

      allocate(row(norb))

      !> Get Electronic energy
      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,nOrb,nOrb,aux1_bml)
      call bml_copy_new(syprt(ipt)%estr%rho,aux_bml)


      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,nOrb,nOrb,rhoat_bml)
      call prg_build_atomic_density(rhoat_bml,tb%numel,syprt(ipt)%estr%hindex,syprt(ipt)%spindex,norb,&
           lt%bml_type)


      call bml_add_deprecated(1.0_dp,aux_bml,-1.0_dp,rhoat_bml,lt%threshold)
      call bml_multiply(aux_bml,syprt(ipt)%estr%ham,aux1_bml,1.0d0, 0.0d0,lt%threshold)
      row=0.0_dp
      call bml_deallocate(rhoat_bml)
      call bml_get_diagonal(aux1_bml,row)

      TRRHOH = 0.0_dp
      do i=1,norb_core
        TRRHOH= TRRHOH+ row(i)
      enddo

      call gpmdcov_message("gpmdcov_EnergAndForces","Energy Band for part =&
      & "//to_string(ipt)//"= "//to_string(TRRHOH),lt%verbose,myRank)

      call bml_deallocate(aux_bml)
      call bml_deallocate(aux1_bml)
      call bml_deallocate(syprt(ipt)%estr%oham)
      deallocate(row)

      syprt(ipt)%estr%eband = TRRHOH
      ebandvector(ipt) = TRRHOH

      dx = 0.0001_dp;

      call get_dH(dx,syprt(ipt)%coordinate,syprt(ipt)%estr%hindex,syprt(ipt)%spindex,intPairsH,onsitesH,syprt(ipt)%symbol,&
           syprt(ipt)%lattice_vector, norb, tb%norbi, lt%bml_type, &
           lt%threshold, dH0x_bml,dH0y_bml,dH0z_bml)

      call get_dS(dx,syprt(ipt)%coordinate,syprt(ipt)%estr%hindex,syprt(ipt)%spindex,intPairsS,onsitesS,syprt(ipt)%symbol,&
           syprt(ipt)%lattice_vector, norb, tb%norbi, lt%bml_type, &
           lt%threshold, dSx_bml,dSy_bml,dSz_bml)

      if(printRank() == 1 .and. lt%verbose >= 10)then
        call bml_print_matrix("dH0x_bml",dH0x_bml,0,10,0,10)
        call bml_print_matrix("dH0y_bml",dH0y_bml,0,10,0,10)
        call bml_print_matrix("dH0z_bml",dH0z_bml,0,10,0,10)

        call bml_print_matrix("dSx_bml",dSx_bml,0,10,0,10)
        call bml_print_matrix("dSy_bml",dSy_bml,0,10,0,10)
        call bml_print_matrix("dSz_bml",dSz_bml,0,10,0,10)
      endif

      call get_skforce(syprt(ipt)%nats,syprt(ipt)%estr%rho,dH0x_bml,dH0y_bml,&
           dH0z_bml,syprt(ipt)%estr%hindex,syprt(ipt)%estr%SKForce,lt%threshold)

      call prg_get_pulayforce(syprt(ipt)%nats,syprt(ipt)%estr%zmat,syprt(ipt)%estr%ham,syprt(ipt)%estr%rho,&
           dSx_bml,dSy_bml,dSz_bml,syprt(ipt)%estr%hindex,syprt(ipt)%estr%FPUL,lt%threshold)

      call get_nonortho_coul_forces(syprt(ipt)%nats, norb, dSx_bml,dSy_bml,dSz_bml,&
           syprt(ipt)%estr%hindex,syprt(ipt)%spindex,syprt(ipt)%estr%rho,syprt(ipt)%net_charge,syprt(ipt)%estr%coul_pot_r,&
           syprt(ipt)%estr%coul_pot_k,tb%hubbardu,syprt(ipt)%estr%FSCOUL,lt%threshold)

      do i=1,gpat%sgraph(ipt)%llsize
        FPUL(:,gpat%sgraph(ipt)%core_halo_index(i)+1) = syprt(ipt)%estr%FPUL(:,i)
        FSCOUL(:,gpat%sgraph(ipt)%core_halo_index(i)+1) = syprt(ipt)%estr%FSCOUL(:,i)
        SKForce(:,gpat%sgraph(ipt)%core_halo_index(i)+1) = syprt(ipt)%estr%SKForce(:,i)
      enddo

      call bml_deallocate(dSx_bml)
      call bml_deallocate(dSy_bml)
      call bml_deallocate(dSz_bml)

      call bml_deallocate(dH0x_bml)
      call bml_deallocate(dH0y_bml)
      call bml_deallocate(dH0z_bml)

      call bml_deallocate(syprt(ipt)%estr%rho)
      call bml_deallocate(syprt(ipt)%estr%ham)
      call bml_deallocate(syprt(ipt)%estr%ham0)
      call bml_deallocate(syprt(ipt)%estr%over)
      call bml_deallocate(syprt(ipt)%estr%zmat)

    enddo

    collectedforce = FPUL + FSCOUL + SKForce

    mls_i = mls()

#ifdef DO_MPI
    if (getNRanks() .gt. 1) then
      !       call prg_sumRealReduceN(collectedforce(1,:), sy%nats)
      !       call prg_sumRealReduceN(collectedforce(2,:), sy%nats)
      !       call prg_sumRealReduceN(collectedforce(3,:), sy%nats)
      call prg_sumRealReduceN(collectedforce, sy%nats*3)
      call prg_sumRealReduceN(ebandvector, gpat%TotalParts)
    endif
#endif

    call gpmdcov_msI("gpmdcov_EnergAndForces","MPI rank finished prg_sumRealReduceN &
        &for Forces"//to_string(mls() - mls_i),lt%verbose,myRank)

    coul_forces =  coul_forces_r + coul_forces_k

    !> Get Repulsive energy and forces
    if(lt%verbose >= 1) call prg_timer_start(dyn_timer,"Get pair pot")
    !     call get_PairPot_contrib(sy%coordinate,sy%lattice_vector,sy%spindex,ppot,PairForces,ERep)
    call get_PairPot_contrib_int(sy%coordinate,sy%lattice_vector,nl%nnIx,nl%nnIy,&
         nl%nnIz,nl%nrnnlist,nl%nnType,sy%spindex,ppot,PairForces,ERep)
    if(lt%verbose >= 1) call prg_timer_stop(dyn_timer,1)

    !> Get Coulombic energy
    ECoul = 0.0;
    do i = 1,sy%nats
      ECoul = ECoul + charges(i)*(tb%hubbardu(sy%spindex(i))*charges(i) + coul_pot_r(i) + coul_pot_k(i) );
    enddo

    Etot = sum(ebandvector(:)) - 0.5_dp*ECoul  + ERep

    if(myRank == 1 .and. lt%verbose >= 1)then
      write(*,*)"Energy Coulomb = ", ECoul
      write(*,*)"Energy Band =", sum(ebandvector(:))
      write(*,*)"Energy Electronic =", Etot
      write(*,*)"Energy Repulsive = ", ERep
    endif

    EPOT = Etot;

    if(.not.allocated(sy%force))allocate(sy%force(3,sy%nats))

    sy%force =  collectedforce + coul_forces + PairForces

    if(myRank == 1 .and. lt%verbose >= 2)then
      write(*,*)""; write(*,*)"FPUL + FSCOUL + SKForce"
      do i = 1,sy%nats
        write(*,*)"Coul Force",i,collectedforce(1,i),collectedforce(2,i),collectedforce(3,i)
      enddo

      write(*,*)""; write(*,*)"Coulomb forces:"
      do i = 1,sy%nats
        write(*,*)"Coul Force",i,coul_forces(1,i),coul_forces(2,i),coul_forces(3,i)
      enddo

      write(*,*)""; write(*,*)"Repulsive forces:"
      do i = 1,sy%nats
        write(*,*)i,PairForces(1,i),PairForces(2,i),PairForces(3,i)
      enddo

      write(*,*)""; write(*,*) "Total forces"
      do i = 1,sy%nats
        write(*,*)i,sy%force(1,i),sy%force(2,i),sy%force(3,i)
      enddo
    endif

    deallocate(ebandvector)

    call gpmdcov_msMem("gpmdcov","After gpmd_EnergAndForces",lt%verbose,myRank)

  end subroutine gpmdcov_EnergAndForces

end module gpmdcov_EnergAndForces_mod
