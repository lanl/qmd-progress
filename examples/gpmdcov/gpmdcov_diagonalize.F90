!> Diagonalization in parallel with subgraphs.
!!
module gpmdcov_diagonalize_mod
  use gpmdcov_vars
  use gpmdcov_mod
  use gpmdcov_writeout_mod

contains

  subroutine gpmdcov_Diagonalize_H0

    call gpmdcov_msMem("gpmdcov_Diagonalize_H0","Before gpmd_Diagonalize",lt%verbose,myRank)

#ifdef DO_MPI
    allocate(norbsInEachCHAtRank(partsInEachRank(myRank)))
    do iptt=1,partsInEachRank(myRank)
      ipt= reshuffle(iptt,myRank)
#else
    allocate(norbsInEachCHAtRank(gpat%TotalParts))
    do iptt = 1,gpat%TotalParts
       ipt = iptt
#endif

      norb = syprt(ipt)%estr%norbs

      !> Initialize the orthogonal versions of ham and rho.
      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%oham)

      !> Orthogonalize ham.
      if(lt%verbose >= 1 .and. myRank == 1) call prg_timer_start(ortho_timer)
      call prg_orthogonalize(syprt(ipt)%estr%ham0,syprt(ipt)%estr%zmat,syprt(ipt)%estr%oham,&
           lt%threshold,lt%bml_type,0)
      if(lt%verbose >= 1 .and. myRank == 1) call prg_timer_stop(ortho_timer)

      !> Allocation
      if(allocated(syprt(ipt)%estr%evals))then 
        deallocate(syprt(ipt)%estr%evals)
        deallocate(syprt(ipt)%estr%dvals)
        call bml_deallocate(syprt(ipt)%estr%evects)
      endif
      allocate(syprt(ipt)%estr%evals(norb))
      allocate(syprt(ipt)%estr%dvals(norb))
      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%evects)

      !> Getting eivenvectors and eigenvalues for Core+Halos
      call prg_get_evalsDvalsEvects(syprt(ipt)%estr%oham, lt%threshold, syprt(ipt)%estr%hindex, &
           &  gpat%sgraph(ipt)%llsize, syprt(ipt)%estr%evals, syprt(ipt)%estr%dvals, &
           & syprt(ipt)%estr%evects, lt%verbose)

      call bml_print_matrix("ham",syprt(ipt)%estr%oham,0,10,0,10)

      norbsInEachCHAtRank(iptt) = size(syprt(ipt)%estr%evals,dim=1)

      call bml_deallocate(syprt(ipt)%estr%oham)
      write(*,*)"IPT",ipt,"syprt(ipt)%estr%evals(i)",syprt(ipt)%estr%evals

    enddo

    write(*,*)"norbsInEachCHAtRank",norbsInEachCHAtRank
    mls_i = mls()

    call gpmdcov_msMem("gpmdcov_Diagonalize", "After gpmd_Diagonalize",lt%verbose,myRank)
  end subroutine gpmdcov_Diagonalize_H0


  subroutine gpmdcov_Diagonalize_H1(nguess)
   real(dp), allocatable :: nguess(:)

    call gpmdcov_msMem("gpmdcov_Diagonalize_H1","Before gpmd_Diagonalize",lt%verbose,myRank)

    mls_i = mls()
#ifdef DO_MPI
    if(allocated(norbsInEachCHAtRank))then 
        deallocate(norbsInEachCHAtRank)
        allocate(norbsInEachCHAtRank(partsInEachRank(myRank)))
    endif
    do iptt=1,partsInEachRank(myRank)
      ipt= reshuffle(iptt,myRank)
#else
    if(allocated(norbsInEachCHAtRank))then 
        deallocate(norbsInEachCHAtRank) 
        allocate(norbsInEachCHAtRank(gpat%TotalParts))
    endif 
    do iptt = 1,gpat%TotalParts
       ipt = iptt
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
      call gpmdcov_msI("gpmdcov_DM_Min","In prg_get_hscf...",lt%verbose,myRank)
      call prg_get_hscf(syprt(ipt)%estr%ham0,syprt(ipt)%estr%over,syprt(ipt)%estr%ham,syprt(ipt)%spindex,&
           &syprt(ipt)%estr%hindex,tb%hubbardu,syprt(ipt)%net_charge,&
           &syprt(ipt)%estr%coul_pot_r,syprt(ipt)%estr%coul_pot_k,lt%mdim,lt%threshold)

      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%oham)

      !> Orthogonalize ham.
      call prg_orthogonalize(syprt(ipt)%estr%ham,syprt(ipt)%estr%zmat,syprt(ipt)%estr%oham,&
           &lt%threshold,lt%bml_type,0)

      if(allocated(syprt(ipt)%estr%evals))then 
        deallocate(syprt(ipt)%estr%evals)
        deallocate(syprt(ipt)%estr%dvals)
        call bml_deallocate(syprt(ipt)%estr%evects)
      endif
        allocate(syprt(ipt)%estr%evals(norb))
        allocate(syprt(ipt)%estr%dvals(norb))
        call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%evects)

      !> Diagonalize ham
      call prg_get_evalsDvalsEvects(syprt(ipt)%estr%oham, lt%threshold, syprt(ipt)%estr%hindex, &
           &  gpat%sgraph(ipt)%llsize, syprt(ipt)%estr%evals, syprt(ipt)%estr%dvals, &
           & syprt(ipt)%estr%evects, lt%verbose)

      norbsInEachCHAtRank(iptt) = size(syprt(ipt)%estr%evals,dim=1)
      
      call bml_deallocate(syprt(ipt)%estr%oham)

    enddo


  end subroutine gpmdcov_Diagonalize_H1


end module gpmdcov_Diagonalize_mod
