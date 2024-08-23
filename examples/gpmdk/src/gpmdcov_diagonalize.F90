!> Diagonalization in parallel with subgraphs.
!!
module gpmdcov_diagonalize_mod
  use gpmdcov_vars
  use gpmdcov_mod
  use gpmdcov_writeout_mod
  use prg_extras_mod
  
#ifdef USE_NVTX
    use gpmdcov_nvtx_mod
#endif

contains

  subroutine gpmdcov_Diagonalize_H0
    implicit none
    real(dp) :: mlsi

    if(myRank == 1) mlsi = mls()
    call gpmdcov_msMem("gpmdcov_Diagonalize_H0","Before gpmd_Diagonalize_H0 ",lt%verbose,myRank)

#ifdef DO_MPI
    if(allocated(norbsInEachCHAtRank))deallocate(norbsInEachCHAtRank)
    allocate(norbsInEachCHAtRank(partsInEachRank(myRank)))
    do iptt=1,partsInEachRank(myRank)
      ipt= reshuffle(iptt,myRank)
#else
    if(allocated(norbsInEachCHAtRank))deallocate(norbsInEachCHAtRank)
    allocate(norbsInEachCHAtRank(gpat%TotalParts))
    do iptt = 1,gpat%TotalParts
      ipt = iptt
#endif

      norb = syprt(ipt)%estr%norbs

      !> Initialize the orthogonal versions of ham and rho.
      if(bml_allocated(syprt(ipt)%estr%oham))call bml_deallocate(syprt(ipt)%estr%oham)
      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%oham)

      !> Orthogonalize ham.
      call prg_orthogonalize(syprt(ipt)%estr%ham0,syprt(ipt)%estr%zmat,syprt(ipt)%estr%oham,&
           lt%threshold,lt%bml_type,0)

      !> Allocation
      if(allocated(syprt(ipt)%estr%evals))deallocate(syprt(ipt)%estr%evals)
      if(allocated(syprt(ipt)%estr%dvals))deallocate(syprt(ipt)%estr%dvals)
      if(bml_allocated(syprt(ipt)%estr%evects)) call bml_deallocate(syprt(ipt)%estr%evects)

      allocate(syprt(ipt)%estr%evals(norb))
      allocate(syprt(ipt)%estr%dvals(norb))
      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%evects)

      !> Getting eivenvectors and eigenvalues for Core+Halos
      call prg_get_evalsDvalsEvects(syprt(ipt)%estr%oham, lt%threshold, syprt(ipt)%estr%hindex, &
           &  gpat%sgraph(ipt)%llsize, syprt(ipt)%estr%evals, syprt(ipt)%estr%dvals, &
           & syprt(ipt)%estr%evects, lt%verbose)

      if(lt%verbose >= 4 .and. myRank == 1) call bml_print_matrix("ham",syprt(ipt)%estr%oham,0,10,0,10)

      norbsInEachCHAtRank(iptt) = size(syprt(ipt)%estr%evals,dim=1)

      call bml_deallocate(syprt(ipt)%estr%oham)

      call gpmdcov_msVectRel("Evals for part "//to_string(ipt),syprt(ipt)%estr%evals,lt%verbose,4,myRank)

    enddo

    call gpmdcov_msVectInt("norbsInEachCHAtRank ",norbsInEachCHAtRank,lt%verbose,2,myRank)

    call gpmdcov_msMem("gpmdcov_Diagonalize_H0", "After gpmd_Diagonalize_H0",lt%verbose,myRank)
    call gpmdcov_msII("gpmdcov_Diagonalize_H0", "Time for gpmd_Diagonalize_H0 "//to_string(mls() - mlsi),lt%verbose,myRank)

  end subroutine gpmdcov_Diagonalize_H0


  subroutine gpmdcov_Diagonalize_H1(nguess)
    real(dp), allocatable :: nguess(:)
    real(dp) :: mlsi

    mlsi = mls()
#ifdef DO_MPI
    if(allocated(norbsInEachCHAtRank)) deallocate(norbsInEachCHAtRank)
    allocate(norbsInEachCHAtRank(partsInEachRank(myRank)))
    do iptt=1,partsInEachRank(myRank)
      ipt= reshuffle(iptt,myRank)
#else
    if(allocated(norbsInEachCHAtRank)) deallocate(norbsInEachCHAtRank)
    allocate(norbsInEachCHAtRank(gpat%TotalParts))
    do iptt = 1,gpat%TotalParts
      ipt = iptt
#endif
      norb = syprt(ipt)%estr%norbs

      call gpmdcov_msMem("gpmdcov_Diagonalize_H1", "Before pot alloc ",lt%verbose,myRank)
      if(allocated(syprt(ipt)%estr%coul_pot_k))deallocate(syprt(ipt)%estr%coul_pot_k)
      if(allocated(syprt(ipt)%estr%coul_pot_r))deallocate(syprt(ipt)%estr%coul_pot_r)
      allocate(syprt(ipt)%estr%coul_pot_k(syprt(ipt)%nats))
      allocate(syprt(ipt)%estr%coul_pot_r(syprt(ipt)%nats))

      syprt(ipt)%estr%coul_pot_k = 0.0_dp
      syprt(ipt)%estr%coul_pot_r = 0.0_dp
      syprt(ipt)%net_charge = 0.0_dp
      call gpmdcov_msMem("gpmdcov_Diagonalize_H1", "After pot alloc ",lt%verbose,myRank)

      !> Get Coulombic potential and charges for the part.
      call gpmdcov_msMem("gpmdcov_Diagonalize_H1", "Before pot asign ",lt%verbose,myRank)
      do j=1,gpat%sgraph(ipt)%lsize
        jj = gpat%sgraph(ipt)%core_halo_index(j)+1
        syprt(ipt)%estr%coul_pot_k(j) = coul_pot_k(jj)
        syprt(ipt)%estr%coul_pot_r(j) = coul_pot_r(jj)
        syprt(ipt)%net_charge(j) = nguess(jj)
      enddo
      call gpmdcov_msMem("gpmdcov_Diagonalize_H1", "After pot asign ",lt%verbose,myRank)

      call gpmdcov_msMem("gpmdcov_Diagonalize_H1", "Before hscf ",lt%verbose,myRank)
      if(bml_allocated(syprt(ipt)%estr%ham)) call bml_deallocate(syprt(ipt)%estr%ham)
      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%ham)
      !call bml_copy_new(syprt(ipt)%estr%ham0,syprt(ipt)%estr%ham)
      call bml_copy(syprt(ipt)%estr%ham0,syprt(ipt)%estr%ham)
      !> Get the scf hamiltonian. The output is ham_bml.
      call gpmdcov_msIII("gpmdcov_DM_Min","In prg_get_hscf...",lt%verbose,myRank)
#ifdef DO_MPI
#ifdef USE_NVTX
        call nvtxStartRange("BarrierBeforePrgGetHscf",2)
        call prg_barrierParallel
        call nvtxEndRange
#endif
#endif
      call prg_get_hscf(syprt(ipt)%estr%ham0,syprt(ipt)%estr%over,syprt(ipt)%estr%ham,syprt(ipt)%spindex,&
           &syprt(ipt)%estr%hindex,tb%hubbardu,syprt(ipt)%net_charge,&
           &syprt(ipt)%estr%coul_pot_r,syprt(ipt)%estr%coul_pot_k,lt%mdim,lt%threshold)
      call gpmdcov_msMem("gpmdcov_Diagonalize_H1", "After hscf ",lt%verbose,myRank)

      call gpmdcov_msMem("gpmdcov_Diagonalize_H1", "Before orth ",lt%verbose,myRank)
      if(bml_allocated(syprt(ipt)%estr%oham)) call bml_deallocate(syprt(ipt)%estr%oham)
      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%oham)
#ifdef DO_MPI
#ifdef USE_NVTX
        call nvtxStartRange("BarrierBeforeOrthogonalizeH1",3)
        call prg_barrierParallel
        call nvtxEndRange
#endif
#endif

      !> Orthogonalize ham.
      call prg_orthogonalize(syprt(ipt)%estr%ham,syprt(ipt)%estr%zmat,syprt(ipt)%estr%oham,&
           &lt%threshold,lt%bml_type,0)
      call gpmdcov_msMem("gpmdcov_Diagonalize_H1", "After orth",lt%verbose,myRank)

      call gpmdcov_msMem("gpmdcov_Diagonalize_H1", "Before evals realloc",lt%verbose,myRank)
      if(allocated(syprt(ipt)%estr%evals)) deallocate(syprt(ipt)%estr%evals)
      if(allocated(syprt(ipt)%estr%dvals)) deallocate(syprt(ipt)%estr%dvals)
      if(bml_allocated(syprt(ipt)%estr%evects)) call bml_deallocate(syprt(ipt)%estr%evects)

      allocate(syprt(ipt)%estr%evals(norb))
      allocate(syprt(ipt)%estr%dvals(norb))
      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%evects)
      call gpmdcov_msMem("gpmdcov_Diagonalize_H1", "After evals realloc",lt%verbose,myRank)

      !> Diagonalize ham
      call gpmdcov_msMem("gpmdcov_Diagonalize_H1", "Before diag",lt%verbose,myRank)
      call prg_get_evalsDvalsEvects(syprt(ipt)%estr%oham, lt%threshold, syprt(ipt)%estr%hindex, &
           &  gpat%sgraph(ipt)%llsize, syprt(ipt)%estr%evals, syprt(ipt)%estr%dvals, &
           & syprt(ipt)%estr%evects, lt%verbose)

#ifdef DO_MPI
#ifdef USE_NVTX
        call nvtxStartRange("BarrierAfterDiagonalizeH1",4)
        call prg_barrierParallel
        call nvtxEndRange
#endif
#endif
      norbsInEachCHAtRank(iptt) = size(syprt(ipt)%estr%evals,dim=1)

      call bml_deallocate(syprt(ipt)%estr%oham)
      call gpmdcov_msMem("gpmdcov_Diagonalize_H1", "After diag",lt%verbose,myRank)

    enddo
    call gpmdcov_msII("gpmdcov_Diagonalize_H1", "Time for gpmd_Diagonalize_H1 "//to_string(mls() - mlsi),lt%verbose,myRank)

  end subroutine gpmdcov_Diagonalize_H1

end module gpmdcov_Diagonalize_mod
