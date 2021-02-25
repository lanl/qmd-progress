!> First computation of charges.
!! \brief Here we compute the first "non-scf charges" based on H0.
!!
subroutine gpmdcov_FirstCharges()

  use gpmdcov_vars
  use gpmdcov_mod
  use gpmdcov_rhosolver_mod
  use gpmdcov_writeout_mod

  call gpmdcov_msMem("gpmdcov_FirstCharges","Before gpmd_FirstCharges",lt%verbose,myRank)

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
    !  do ipt = 1,gpat%TotalParts
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

    norbsInEachCHAtRank(iptt) = size(syprt(ipt)%estr%aux(1,:),dim=1)

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

  ! End of loop over parts.
  mls_i = mls()

  ! Get the chemical potential. Feb 2021 implemetation
  call gpmdcov_getmu()

#ifdef DO_MPI
  if (getNRanks() .gt. 1) then
    call prg_sumRealReduceN(sy%net_charge, sy%nats)
  endif
#endif

  call gpmdcov_msI("gpmdcov_FirstCharges","MPI rank finished prg_sumIntReduceN for qs ="//to_string(mls() - mls_i),lt%verbose,myRank)

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
  call gpmdcov_msMem("gpmdcov", "After gpmd_FirstCharges",lt%verbose,myRank)

end subroutine gpmdcov_FirstCharges


