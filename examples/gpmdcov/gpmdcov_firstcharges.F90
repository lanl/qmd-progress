!> First computation of charges.
!! \brief Here we compute the first "non-scf charges" based on H0.
!!
subroutine gpmdcov_FirstCharges()

  use gpmdcov_vars

  if(lt%verbose >= 1 .and. myRank == 1)call prg_get_mem("gpmdcov", "Before gpmd_FirstCharges")

  if(.not.allocated(sy%net_charge))allocate(sy%net_charge(sy%nats))

!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Loop over parts
!!!!!!!!!!!!!!!!!!!!!!!!!!!

  sy%net_charge = 0.0_dp

#ifdef DO_MPI
  !!do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)

  allocate(norbsInEachCHAtRank(partsInEachRank(myRank)))

  do iptt=1,partsInEachRank(myRank)
    ipt= reshuffle(iptt,myRank)
#else
  do ipt = 1,gpat%TotalParts
#endif

    norb = syprt(ipt)%estr%norbs

    call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%rho)

    !> Initialize the orthogonal versions of ham and rho.
    call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%oham)
    call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%orho)

    !> Orthogonalize ham.
    if(lt%verbose >= 1 .and. myRank == 1) call prg_timer_start(ortho_timer)
    call prg_orthogonalize(syprt(ipt)%estr%ham0,syprt(ipt)%estr%zmat,syprt(ipt)%estr%oham,&
         lt%threshold,lt%bml_type,lt%verbose)
    if(lt%verbose >= 1 .and. myRank == 1) call prg_timer_stop(ortho_timer)

    call gpmdcov_RhoSolver(syprt(ipt)%estr%oham,syprt(ipt)%estr%orho)
    write(*,*)syprt(ipt)%estr%aux(3,:)
    syprt(ipt)%estr%norbsInCore = size(syprt(ipt)%estr%aux(1,:))
    norbsInEachCHAtRank(iptt) = syprt(ipt)%estr%norbsInCore

    !> Deorthogonalize rho.
    if(lt%verbose >= 1 .and. myRank == 1) call prg_timer_start(deortho_timer)
    call prg_deorthogonalize(syprt(ipt)%estr%orho,syprt(ipt)%estr%zmat,syprt(ipt)%estr%rho,&
         lt%threshold,lt%bml_type,lt%verbose)
    if(lt%verbose >= 1 .and. myRank == 1) call prg_timer_stop(deortho_timer)

    !> Get charges based on rho. rho_bml is the input and sy%net_charge is the outputs vector containing
    !! the charges.
    call prg_get_charges(syprt(ipt)%estr%rho, syprt(ipt)%estr%over, syprt(ipt)%estr%hindex, syprt(ipt)%net_charge, tb%numel,&
         syprt(ipt)%spindex, lt%mdim, lt%threshold)

    if(lt%verbose.GE.4 .and. myRank  ==  1)then
      write(*,*)""
      write(*,*)"Total charge of part "//to_string(ipt) &
           //" = "//to_string(sum(syprt(ipt)%net_charge(:),size(syprt(ipt)%net_charge,dim=1)))
      write(*,*)""
      write(*,*)""; write(*,*)"Part charges:"
      do j=1,syprt(ipt)%nats
        write(*,*)syprt(ipt)%symbol(j),syprt(ipt)%net_charge(j)
      enddo
    endif

    do j=1,gpat%sgraph(ipt)%llsize
      jj = gpat%sgraph(ipt)%core_halo_index(j)+1
      sy%net_charge(jj) = syprt(ipt)%net_charge(j)
    enddo

    call bml_deallocate(syprt(ipt)%estr%oham)
    call bml_deallocate(syprt(ipt)%estr%orho)
    call bml_deallocate(syprt(ipt)%estr%rho)

  enddo

  write(*,*)"norbsInEachCHAtRank",myRank,norbsInEachCHAtRank
  norbsInRank = sum(norbsInEachCHAtRank)
  write(*,*)"Number of orbitals in rank",myRank, "=", norbsInRank
  write(*,*)"Total number of parts", gpat%TotalParts
  allocate(evalsInRank(norbsInRank))
  allocate(fvalsInRank(norbsInRank))
  allocate(dvalsInRank(norbsInRank))

  shift = 0
  do iptt=1,partsInEachRank(myRank)
    ipt= reshuffle(iptt,myRank)
    do i = 1, norbsInEachCHAtRank(iptt)
      evalsInRank(i+shift) = syprt(ipt)%estr%aux(1,i)
      fvalsInRank(i+shift) = syprt(ipt)%estr%aux(2,i)
      dvalsInRank(i+shift) = syprt(ipt)%estr%aux(3,i)
    enddo
    shift = shift + norbsInEachCHAtRank(iptt)
  enddo

  call prg_wait()
  write(*,*)"evalsInRank",evalsInRank
  allocate(norbsInEachCH(gpat%TotalParts))
  norbsInEachCH = 0
  nRanks = getNRanks()
  write(*,*)"Total Number of Raks =",nRanks
  allocate(npartsVect(nRanks))
  allocate(displ(nRanks))
  allocate(norbsInEachRank(nRanks))
  do i = 1,nRanks
    npartsVect(i) = partsInEachRank(i)

    if(i == 1)then
      shift = 0
    else
      shift = shift + partsInEachRank(i-1)
    endif
    displ(i) = shift
  enddo

  write(*,*)"norbsInEachCH",norbsInEachCH
  write(*,*)"partsInEachRank",partsInEachRank
  write(*,*)"displ",displ

  call allGatherVIntParallel(norbsInEachCHAtRank, partsInEachRank(myRank),norbsInEachCH ,npartsVect, displ)

  write(*,*)"norbsInEachCH",norbsInEachCH
  kk = 0
  deallocate(displ); allocate(displ(nRanks))
  norbsInEachRank = 0.0_dp
  do i = 1,nRanks
    do j = 1,partsInEachRank(i)
      kk = kk + 1
      norbsInEachRank(i) = norbsInEachRank(i) + norbsInEachCH(kk)
    enddo
    if(i == 1)then
      shift = 0
    else
      shift = shift + norbsInEachRank(i-1)
    endif
    displ(i) = shift
  enddo
  totalNorbs = sum(norbsInEachRank)
  write(*,*)"norbsInRank",norbsInEachRank
  write(*,*)"norbsInEachCH",norbsInEachCH

  allocate(evalsAll(totalNorbs))
  allocate(fvalsAll(totalNorbs))
  allocate(dvalsAll(totalNorbs))
  num = getNRanks() + 220
  if(getNRanks() .EQ. 1)then
    do i = 1,norbsInRank
      write(220,*)evalsInRank(i)
    enddo
  else
    do i = 1,norbsInRank
      write(221,*)evalsInRank(i)
    enddo
  endif

  call prg_wait()

  write(*,*)"sizes",size(evalsInRank,dim=1),size(evalsAll,dim=1),size(norbsInEachRank,dim=1),size(displ,dim=1)
  write(*,*)"displ",displ
  call allGatherVRealParallel(evalsInRank, norbsInRank, evalsAll ,norbsInEachRank, displ)
  call allGatherVRealParallel(fvalsInRank, norbsInRank, fvalsAll ,norbsInEachRank, displ)
  call allGatherVRealParallel(dvalsInRank, norbsInRank, dvalsAll ,norbsInEachRank, displ)

  deallocate(displ)

  write(*,*)"norbsInEachCH",norbsInEachCH
  write(*,*)"dvalsAll",dvalsAll


  if (myRank .eq. 1) then
    do i = 1,size(evalsAll,dim=1)
      write(121,*)evalsAll(i)
      write(122,*)dvalsAll(i)
    enddo
  endif
  Ef = 0.5_dp*(maxval(evalsAll)-minval(evalsAll))
  nocc = bndfilTotal*real(sy%estr%norbs,dp)
  write(*,*)"nel, nocc", sy%estr%nel, nocc
  call gpmdcov_musearch(evalsAll, fvalsAll, dvalsAll, beta, nocc, 100, 10d-10, Ef, 1)

  write(*,*)"Mu",Ef

  stop

  ! End of loop over parts.

  mls_i = mls()

#ifdef DO_MPI
  if (getNRanks() .gt. 1) then
    call prg_sumRealReduceN(sy%net_charge, sy%nats)
  endif
#endif

  if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"MPI rank finished prg_sumIntReduceN for qs", mls() - mls_i

  !> Gather charges from all the parts.
  if(.not.allocated(charges_old))allocate(charges_old(sy%nats))

  if (printRank()  ==  1) then
    if(lt%verbose >= 2)then
      write(*,*)""
      write(*,*)"Total charge of the system = "//to_string(sum(sy%net_charge(:),size(sy%net_charge,dim=1)))
      write(*,*)""
      if(lt%verbose >= 5) call prg_write_system(sy,"charged_system","pdb")
      write(*,*)""; write(*,*)"Full System charges:"
      do j=1,sy%nats
        write(*,*)j,sy%symbol(j),sy%net_charge(j)
      enddo
    endif
  endif

  sy%net_charge(:) = sy%net_charge(:)- sum(sy%net_charge(:))/real(sy%nats)

  charges_old = sy%net_charge
  if(lt%verbose >= 1 .and. myRank == 1)call prg_get_mem("gpmdcov", "After gpmd_FirstCharges")
end subroutine gpmdcov_FirstCharges


