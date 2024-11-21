!>  Initialize the partition.
!!
subroutine gpmdcov_InitParts
  use gpmdcov_vars
  use gpmdcov_writeout_mod
  use gpmdcov_nonequilibrium_mod, only : gpmdcov_apply_voltage 
#ifdef USE_LATTE
  use gpmdcov_latte_mod
#endif
#ifdef USE_NVTX
    use gpmdcov_nvtx_mod
#endif

  

  implicit none 
  integer :: norbsCore

  call gpmdcov_msMem("gpmdcov","Before gpmd_InitParts",lt%verbose,myRank)

  !Vectorizing part sizes
  if(allocated(partsSizes)) deallocate(partsSizes)
  allocate(partsSizes(gpat%TotalParts))

  partsSizes=0.0_dp
  do i = 1,gpat%TotalParts
    partsSizes(i) = gpat%sgraph(i)%llsize
  enddo

#ifdef DO_MPI
  !do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)
  do iptt=1,partsInEachRank(myRank)
    ipt= reshuffle(iptt,myRank)
#else
  do ipt = 1,gpat%TotalParts
#endif

    if(lt%verbose >= 3)then
      write(*,*)""
      write(*,*)"         #####################################"
      write(*,*)"         (rank "//to_string(myRank)//") Initializing partition "//to_string(ipt)
      write(*,*)"         #####################################"
      write(*,*)""
    end if

    if(lt%verbose >= 3)then
      write(*,*)"Number of atoms in the core =      "//to_string(gpat%sgraph(ipt)%llsize)
      write(*,*)"Number of atoms in the core+halo = "//to_string(gpat%sgraph(ipt)%lsize)
      write(*,*)""
    end if

    !> Get the mapping of the Hamiltonian index with the atom index
    if(allocated(syprt(ipt)%estr%hindex))deallocate(syprt(ipt)%estr%hindex)
    allocate(syprt(ipt)%estr%hindex(2,syprt(ipt)%nats))

    !call get_hindex(syprt(ipt)%spindex,tb%norbi,syprt(ipt)%estr%hindex,norb,norbCores)

    !We use a new routine to get the number of orbitals in the core
    call get_hindex_coreHalo(syprt(ipt)%spindex,gpat%sgraph(ipt)%llsize,tb%norbi,syprt(ipt)%estr%hindex,norb,norbsCore,lt%verbose)
    syprt(ipt)%estr%norbs = norb
    syprt(ipt)%estr%norbsCore = norbsCore

    if(bml_allocated(syprt(ipt)%estr%ham0))then 
            call bml_deallocate(syprt(ipt)%estr%ham0)
            call bml_deallocate(syprt(ipt)%estr%over)
    endif



    call bml_noinit_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%ham0)
    call bml_noinit_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%over)
#ifdef DO_MPI
#ifdef USE_NVTX
        call nvtxStartRange("BarrierBeforeGetHS",2)
        call prg_barrierParallel
        call nvtxEndRange
#endif
#endif

    !Construction of the hamiltonian for every part
#ifdef USE_LATTE
    call get_hsmat_latte(syprt(ipt)%estr%ham0,syprt(ipt)%estr%over,syprt(ipt)%coordinate,syprt(ipt)%symbol,&
        syprt(ipt)%lattice_vector,lt%threshold)
#else
    if(gpmdt%usevectsk)then
       call get_hsmat_vect(syprt(ipt)%estr%ham0,syprt(ipt)%estr%over,syprt(ipt)%coordinate,&
            syprt(ipt)%lattice_vector,syprt(ipt)%spindex,&
            tb%norbi,syprt(ipt)%estr%hindex,onsitesH,onsitesS,intPairsH,intPairsS,lt%threshold)
    else
       call get_hsmat(syprt(ipt)%estr%ham0,syprt(ipt)%estr%over,syprt(ipt)%coordinate,&
            syprt(ipt)%lattice_vector,syprt(ipt)%spindex,&
            tb%norbi,syprt(ipt)%estr%hindex,onsitesH,onsitesS,intPairsH,intPairsS,lt%threshold)
    endif
#endif

    if(gpmdt%applyv)then
        call gpmdcov_apply_voltage(sy%nats,syprt(ipt)%nats,syprt(ipt)%estr%hindex,gpat%sgraph(ipt)%core_halo_index,&
            &syprt(ipt)%estr%ham0,syprt(ipt)%estr%over)
    endif 

    if (myRank  ==  1 .and. lt%verbose >= 5)then
      write(*,*)"H0 and S for part:"
      call bml_print_matrix("H0",syprt(ipt)%estr%ham0,0,6,0,6)
      call bml_print_matrix("S",syprt(ipt)%estr%over,0,6,0,6)
      if (lt%verbose >= 6)then
        call bml_write_matrix(syprt(ipt)%estr%ham0,"H0.mtx")
        call bml_write_matrix(syprt(ipt)%estr%over,"S.mtx")
    !    stop
      endif
    endif
#ifdef DO_MPI
#ifdef USE_NVTX
        call nvtxStartRange("BarrierAfterGetHS",3)
        call prg_barrierParallel
        call nvtxEndRange
#endif
#endif

    !> Get occupation based on last shell population.
    !  WARNING: This could change depending on the TB method being used.
    !nel = sum(element_numel(syprt(ipt)%atomic_number(:)),&
    !     size(syprt(ipt)%atomic_number,dim=1))
    nel = sum(element_numel(syprt(ipt)%atomic_number(:)))
    bndfil = nel/(2.0_dp*norb)

    !> Initialize the density matrix (rho_bml) and inverse overlap factor (zmat_bml).
    if(bml_allocated(syprt(ipt)%estr%zmat)) call bml_deallocate(syprt(ipt)%estr%zmat)
    call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,syprt(ipt)%estr%zmat)

    !> Get the Inverse square root overlap matrix.
    if(lt%verbose >= 3) call prg_timer_start(dyn_timer,"Build Z for part")
    call gpmdcov_buildz(syprt(ipt)%estr%over,syprt(ipt)%estr%zmat)
    if(lt%verbose >= 3) call prg_timer_stop(dyn_timer,1)

#ifdef DO_MPI
#ifdef USE_NVTX
        call nvtxStartRange("BarrierAfterGenZ",4)
        call prg_barrierParallel
        call nvtxEndRange
#endif
#endif
    if(myRank == 1 .and. lt%verbose >= 5)then
      write(*,*)"Z matrix for part:"
      call bml_print_matrix("Z",syprt(ipt)%estr%zmat,0,6,0,6)
    endif

  enddo
  call gpmdcov_msMem("gpmdcov","After gpmd_InitParts",lt%verbose,myRank)

end subroutine gpmdcov_InitParts

