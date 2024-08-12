module gpmdcov_EnergAndForces_mod

  use gpmdcov_mod

  contains

  subroutine gpmdcov_EnergAndForces(charges)
    
    use gpmdcov_vars
    use gpmdcov_writeout_mod
    use gpmdcov_kernel_mod
    use gpmdcov_constraints_mod

    Implicit none
    real(dp), intent(in) :: charges(:)
    real(dp), allocatable :: ebandvector(:)
    real(dp), allocatable ::  R1(:), R2(:)
    real(dp) :: dcoords(3), dist
    real(dp) :: smd_total_force(3), smd_total_energy
    real(dp) :: smd_total_energy_allpairs
    real(dp) :: smd_test_force(3), smd_test_energy
    real(dp) :: deltas(3)
    real(dp) :: delta_h, energy_plus, energy_minus, denergy, differ
    integer :: k
    logical :: testsmd

    call gpmdcov_msMem("gpmdcov","Before gpmd_EnergAndForces",lt%verbose,myRank)

    if(.not.allocated(coul_forces)) allocate(coul_forces(3,sy%nats))
    if(.not.allocated(GFPUL))allocate(GFPUL(3,sy%nats))
    if(.not.allocated(GFSCOUL))allocate(GFSCOUL(3,sy%nats))
    if(.not.allocated(SKForce))allocate(SKForce(3,sy%nats))
    if(.not.allocated(collectedforce))allocate(collectedforce(3,sy%nats))
    if(.not.allocated(ebandvector))allocate(ebandvector(gpat%TotalParts))

    GFPUL = 0.0_dp
    GFSCOUL = 0.0_dp
    SKForce = 0.0_dp
    collectedforce = 0.0_dp
    ebandvector = 0.0_dp
    smd_total_force(:) = 0.0_dp
    smd_total_energy = 0.0_dp
    smd_total_energy_allpairs = 0.0_dp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Loop over all the parts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef DO_MPI
    !do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)
    do iptt=1,partsInEachRank(myRank)
      ipt= reshuffle(iptt,myRank)
#else
    do ipt = 1,gpat%TotalParts
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

      if(gpmdt%usevectsk)then

         call get_dH_or_dS_vect(dx,syprt(ipt)%coordinate,syprt(ipt)%estr%hindex,&
              &syprt(ipt)%spindex,intPairsH,onsitesH,syprt(ipt)%symbol,&
              &syprt(ipt)%lattice_vector, norb, tb%norbi, lt%bml_type, &
              &lt%threshold, dH0x_bml,dH0y_bml,dH0z_bml)
         
         call get_dH_or_dS_vect(dx,syprt(ipt)%coordinate,syprt(ipt)%estr%hindex,&
              &syprt(ipt)%spindex,intPairsS,onsitesS,syprt(ipt)%symbol,&
              &syprt(ipt)%lattice_vector, norb, tb%norbi, lt%bml_type, &
              &lt%threshold, dSx_bml,dSy_bml,dSz_bml)
      else
         call get_dH(dx,syprt(ipt)%coordinate,syprt(ipt)%estr%hindex,&
              &syprt(ipt)%spindex,intPairsH,onsitesH,syprt(ipt)%symbol,&
              &syprt(ipt)%lattice_vector, norb, tb%norbi, lt%bml_type, &
              &lt%threshold, dH0x_bml,dH0y_bml,dH0z_bml)

         call get_dS(dx,syprt(ipt)%coordinate,syprt(ipt)%estr%hindex,&
              &syprt(ipt)%spindex,intPairsS,onsitesS,syprt(ipt)%symbol,&
              &syprt(ipt)%lattice_vector, norb, tb%norbi, lt%bml_type, &
              &lt%threshold, dSx_bml,dSy_bml,dSz_bml)
      endif

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

      !call prg_PulayComponentT(syprt(ipt)%estr%rho,syprt(ipt)%estr%ham,syprt(ipt)%estr%zmat,syprt(ipt)%estr%FPUL,lt%threshold &
      ! &,lt%mdim,lt%bml_type,lt%verbose)


      call get_nonortho_coul_forces(syprt(ipt)%nats, norb, dSx_bml,dSy_bml,dSz_bml,&
           syprt(ipt)%estr%hindex,syprt(ipt)%spindex,syprt(ipt)%estr%rho,syprt(ipt)%net_charge,syprt(ipt)%estr%coul_pot_r,&
           syprt(ipt)%estr%coul_pot_k,tb%hubbardu,syprt(ipt)%estr%FSCOUL,lt%threshold)

      do i=1,gpat%sgraph(ipt)%llsize
        GFPUL(:,gpat%sgraph(ipt)%core_halo_index(i)+1) = syprt(ipt)%estr%FPUL(:,i)
        GFSCOUL(:,gpat%sgraph(ipt)%core_halo_index(i)+1) = syprt(ipt)%estr%FSCOUL(:,i)
        SKForce(:,gpat%sgraph(ipt)%core_halo_index(i)+1) = syprt(ipt)%estr%SKForce(:,i)
      enddo

      call bml_deallocate(dSx_bml)
      call bml_deallocate(dSy_bml)
      call bml_deallocate(dSz_bml)

      call bml_deallocate(dH0x_bml)
      call bml_deallocate(dH0y_bml)
      call bml_deallocate(dH0z_bml)

      !if(.not.kernel%xlbolevel1)then
      call bml_deallocate(syprt(ipt)%estr%rho)
      call bml_deallocate(syprt(ipt)%estr%ham)
      call bml_deallocate(syprt(ipt)%estr%ham0)
      !endif
      !call bml_deallocate(syprt(ipt)%estr%over)
      !call bml_deallocate(syprt(ipt)%estr%zmat)

    enddo

    collectedforce = GFPUL + GFSCOUL + SKForce
    !collectedforce =  GFSCOUL + SKForce

    call gpmdcov_msMem("gpmdcov","Before steered MD (SMD) check",lt%verbose,myRank)
    
    !> Steered MD Force (is using SMD)
    if(gpmdt%usesmd) then

      !> Allocate SMD Arrays
      !! R1 and R2 for xyz coordinates for steered atoms in each pair
      !! directional_smd_force for xyz force components
      if(.not.allocated(R1))allocate(R1(3))
      if(.not.allocated(R2))allocate(R2(3))

      do i = 1,gpmdt%smdnumpairs
        !> Determine coordinates for steered atoms in each pair 
        R1 = sy%coordinate(:,gpmdt%smdatomind1(i))
        R2 = sy%coordinate(:,gpmdt%smdatomind2(i))

        write(*,*) "gpmdcov_EnergAndForces   SMD Pair Number ",i," Atoms ",gpmdt%smdatomind1(i),&
                &" and ",gpmdt%smdatomind2(i)
        !> Call constraints subroutine, harmonic to linear
        !! collectedforce will be updated for steered atoms
        call gpmdcov_constraint_harmonicToLinear(R1, R2, smd_total_force, smd_total_energy, lt%verbose)
     
        !> Update collectedforce to include SMD force for steered atoms
        collectedforce(:,gpmdt%smdatomind1(i)) = collectedforce(:,gpmdt%smdatomind1(i)) + smd_total_force(:)
        collectedforce(:,gpmdt%smdatomind2(i)) = collectedforce(:,gpmdt%smdatomind2(i)) - smd_total_force(:)

        !> Print velocity logging for steered atoms and total SMD energy/force terms
        if(MDstep .gt. 2) then
              call gpmdcov_msI("gpmdcov_EnergAndForces","Velocity of first steered atom " &
                      &//to_string(norm2(sy%velocity(:,gpmdt%smdatomind1(i)))), lt%verbose, myRank)
              call gpmdcov_msI("gpmdcov_EnergAndForces","Velocity of second steered atom " &
                      &//to_string(norm2(sy%velocity(:,gpmdt%smdatomind2(i)))), lt%verbose, myRank)
        endif
        call gpmdcov_msI("gpmdcov_EnergAndForces","SMD Force Magnitude " &
                & // to_string(norm2(smd_total_force)),lt%verbose,myRank)
        call gpmdcov_msI("gpmdcov_EnergAndForces","SMD Total Energy " &
             & // to_string(smd_total_energy),lt%verbose,myRank)
        do k = 1,3
           dcoords(k) = modulo(((R1(k) - R2(k)) + &
                   &0.5_dp*sy%lattice_vector(k,k)),sy%lattice_vector(k,k)) - &
                   &0.5_dp * sy%lattice_vector(k,k)
        enddo
        dist = norm2(dcoords)
        call gpmdcov_msI("gpmdcov_EnergAndForces","SMD distance " &
             & // to_string(dist),lt%verbose,myRank)

        !> Add energy for current steered atom pair to total steering energy
        !! This will then be added to the total energy of the system
        smd_total_energy_allpairs = smd_total_energy_allpairs + smd_total_energy
      enddo

      !> SMD Finite Difference Test
      !! Only test for last SMD atom pair
      !! Set testsmd to true to run finite difference derivative tests
      testsmd = .false.
      if(testsmd) then
              write(*,*) "gpmdcov_EnergAndForces testing SMD"
              call gpmdcov_msI("gpmdcov_EnergAndForces","Testing SMD Derivatives: SMD force (x-direction) " &
                      &//to_string(smd_total_force(1)), lt%verbose, myRank)
              !! Test with three different delta h values
              !! Current testing is only done in the x-direction by updating the x coordinate of R1
              deltas = (/0.001, 0.01, 0.1 /)
              do k=1,3
                  delta_h = deltas(k)
                  R1(1) = R1(1) + delta_h
                  call gpmdcov_constraint_harmonicToLinear(R1, R2, smd_test_force, smd_test_energy, &
                                                           lt%verbose)
                  energy_plus = smd_test_energy
                  R1(1) = R1(1) - 2.0_dp * delta_h
                  call gpmdcov_constraint_harmonicToLinear(R1, R2, smd_test_force, smd_test_energy, &
                                                           lt%verbose)
                  energy_minus = smd_test_energy

                  !! Compute finite difference derivative of energy to compare to output force
                  denergy = (energy_minus - energy_plus)/(2.0_dp*delta_h)
                  call gpmdcov_msI("gpmdcov_EnergAndForces","Testing SMD Derivatives: d energy (x-direction) " &
                          &//to_string(denergy), lt%verbose, myRank)
                  differ = (smd_total_force(1) - denergy)
                  call gpmdcov_msI("gpmdcov_EnergAndForces","Testing SMD Derivatives: difference (force - denergy)  " &
                          &//to_string(differ), lt%verbose, myRank)
                  call gpmdcov_msI("gpmdcov_EnergAndForces","Testing SMD Derivatives: delta h used " &
                          &//to_string(delta_h), lt%verbose, myRank)

                  !> Reset R1(1) for next test
                  R1(1) = R1(1) + delta_h
              enddo
      endif

      !> Deallocate SMD Arrays
      deallocate(R1)
      deallocate(R2)
    endif


    mls_i = mls()

#ifdef DO_MPI
    if (getNRanks() .gt. 1) then
      !       call prg_sumRealReduceN(collectedforce(1,:), sy%nats)
      !       call prg_sumRealReduceN(collectedforce(2,:), sy%nats)
       !       call prg_sumRealReduceN(collectedforce(3,:), sy%nats)
      call prg_barrierParallel
      call prg_sumRealReduceN(collectedforce, sy%nats*3)
      call prg_sumRealReduceN(ebandvector, gpat%TotalParts)
    endif
#endif

    call gpmdcov_msI("gpmdcov_EnergAndForces","MPI rank finished prg_sumRealReduceN &
        &for Forces"//to_string(mls() - mls_i),lt%verbose,myRank)

    coul_forces =  coul_forces_r + coul_forces_k

    !> Get Repulsive energy and forces
    !     call get_PairPot_contrib(sy%coordinate,sy%lattice_vector,sy%spindex,ppot,PairForces,ERep)
    call get_PairPot_contrib_int(sy%coordinate,sy%lattice_vector,nl%nnIx,nl%nnIy,&
         nl%nnIz,nl%nrnnlist,nl%nnType,sy%spindex,ppot,PairForces,ERep)

    !> Get Coulombic energy
    ECoul = 0.0;
    do i = 1,sy%nats
      ECoul = ECoul + charges(i)*(tb%hubbardu(sy%spindex(i))*charges(i) + coul_pot_r(i) + coul_pot_k(i) );
    enddo

    Etot = sum(ebandvector(:)) - 0.5_dp*ECoul  + ERep + smd_total_energy_allpairs

    entropy = 0.0_dp
    if(lt%Entropy)then
        if(lt%MuCalcType == "FromParts" .or. lt%MuCalcType == "Combined")then
          call gpmdcov_getentropy(evalsAll, dvalsAll, beta, Ef, myRank,entropy, lt%verbose)
        else
          write(*,*)"ERROR: Entropy calculation is only valid for MuCalcType= FromParts, or MuCalcType= Combined."
          stop
        endif
    endif
     
    EPOT = Etot + entropy

    if((myRank == 1) .and. (lt%verbose >= 2))then
      write(*,*)"Energy Coulomb = ", ECoul
      write(*,*)"Energy Band =", sum(ebandvector(:))
      write(*,*)"Energy Repulsive = ", ERep
      write(*,*)"Energy Entropic = ", entropy
      write(*,*)"Energy Electronic (Total) =", EPot
    endif

    if(.not.allocated(sy%force))allocate(sy%force(3,sy%nats))

    !TOTAL FORCES
    sy%force =  collectedforce +  PairForces + coul_forces
    !sy%force =  SKForce + GFSCOUL + GFPUL +  PairForces + coul_forces
    !sy%force =  SKForce + GFSCOUL + GFPUL   PairForces + coul_forces
    !sy%force =  coul_forces
    !write(*,*)"FORCESSS",sy%force
    !sy%force =  SKForce + GFSCOUL +  PairForces + coul_forces
    !sy%force =  GFSCOUL 
    !sy%force = SKForce + GFSCOUL + coul_forces + PairForces
    !write(*,*)"FORCES!!!",SKForce,GFSCOUL,coul_forces
    !sy%force =   collectedforce
    !sy%force =   coul_forces

    if(myRank == 1 .and. lt%verbose >= 3)then
      write(*,*)""; write(*,*)"FPUL + FSCOUL + SKForce"
      do i = 1,sy%nats
        write(*,*)"Collected Force",i,collectedforce(1,i),collectedforce(2,i),collectedforce(3,i)
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

  subroutine gpmdcov_getentropy(evals, dvals, beta, mu, rank,entropy, verbose)
    use gpmdcov_writeout_mod
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer                ::  i, norbs
    integer, intent(in)    ::  rank
    real(dp)  ::  fermi
    real(dp), intent(in)   ::  evals(:), dvals(:)
    real(dp), intent(in)   ::  beta
    real(dp), intent(inout)   ::  mu, entropy
    integer, optional, intent(in) :: verbose
    real(dp), allocatable :: fvals(:)

    if (.not.allocated(fvals))then
       allocate(fvals(size(evals)))
    endif
    
    call gpmdcov_msMem("gpmdcov_getentropy","Getting entropic energy contribution ...",verbose,rank)
    
    norbs = size(evals, dim = 1)
    
    entropy = 0.0_dp
    fermi = 0.0_dp
    call gpmdcov_fermifunction(beta,evals,mu,fvals)
    do i = 1,norbs
      fermi = fvals(i)
      if(abs(fermi) > 10d-9 .and. abs(fermi-1.0_dp) > 10d-9)then
        entropy = entropy + (2.0_dp/beta)*dvals(i)*(fermi*log(fermi) + (1.0_dp-fermi)*log(1.0_dp-fermi))
      endif
    enddo

  end subroutine gpmdcov_getentropy

end module gpmdcov_EnergAndForces_mod

