module gpmdcov_Part_mod

contains

  !> Partition by systems
  !!
  subroutine gpmdcov_Part

    use gpmdcov_vars
    use gpmdcov_reshuffle_mod
    use gpmdcov_writeout_mod

    integer, allocatable :: graph_h(:,:)
    integer, allocatable :: graph_p(:,:)
    real(dp)             :: mls_ii
    integer              :: iipt

    if(mdstep < 1)then
      call gpmdcov_msMem("gpmdcov_Part", "Before prg_get_covgraph",lt%verbose,myRank)
      call prg_get_covgraph(sy,nl%nnStructMindist,nl%nnStruct,nl%nrnnstruct&
           ,gsp2%bml_type,gsp2%covgfact,g_bml,gsp2%mdim,lt%verbose)
      call gpmdcov_msMem("gpmdcov_Part", "After prg_get_covgraph",lt%verbose,myRank)
    else

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Anders' way of graph construction.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call gpmdcov_msI("gpmdcov_Part","In prg_get_covgraph_h ...",lt%verbose,myRank)
      mls_ii = mls()
      call prg_get_covgraph_h(sy,nl%nnStructMindist,nl%nnStruct,nl%nrnnstruct,gsp2%nlgcut,graph_h,gsp2%mdim,lt%verbose)
      call gpmdcov_msI("gpmdcov_Part","In prg_get_covgraph_h ..."//to_string(mls()-mls_ii)//" ms",lt%verbose,myRank)

#ifdef DO_MPI
      !do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)
      do iipt=1,partsInEachRank(myRank)
        ipt= reshuffle(iipt,myRank)
#else
      do ipt = 1,gpat%TotalParts
#endif

        call prg_collect_graph_p(syprt(ipt)%estr%orho,gpat%sgraph(ipt)%llsize,sy%nats,syprt(ipt)%estr%hindex,&
             gpat%sgraph(ipt)%core_halo_index,graph_p,gsp2%gthreshold,gsp2%mdim,lt%verbose)

        call bml_deallocate(syprt(ipt)%estr%orho)

      enddo

      mls_i = mls()

#ifdef DO_MPI
      if (getNRanks() > 1) then
        call prg_sumIntReduceN(graph_p, mdim*sy%nats)
      endif
#endif

      call gpmdcov_msI("gpmdcov_Part","Time for prg_sumIntReduceN for graph "//to_string(mls() - mls_i)//" ms",lt%verbose,myRank)

      call gpmdcov_msI("gpmdcov_Part","In prg_merge_graph ...",lt%verbose,myRank)
      mls_ii = mls()
      call prg_merge_graph(graph_p,graph_h)
      call gpmdcov_msI("gpmdcov_Part","Time for prg_merge_graph "//to_string(mls()-mls_ii)//" ms",lt%verbose,myRank)

      deallocate(graph_h)

      !Transform graph into bml format.
      if(mod(mdstep,gsp2%parteach)==0.or.mdstep == 0 .or.mdstep == 1)then
        if(bml_get_N(gcov_bml).GT.0) call bml_deallocate(g_bml)
        ! call bml_zero_matrix(gsp2%bml_type,bml_element_real,kind(1.0),sy%nats,mdim,g_bml,lt%bml_dmode)
        call bml_zero_matrix(gsp2%bml_type,bml_element_real,kind(1.0),sy%nats,mdim,g_bml)

        call gpmdcov_msMem("gpmdcov_Part","Before prg_graph2bml",lt%verbose,myRank)
        mls_ii = mls()
        call prg_graph2bml(graph_p,gsp2%bml_type,g_bml)
        call gpmdcov_msMem("gpmdcov_Part","After prg_graph2bml",lt%verbose,myRank)
        call gpmdcov_msI("gpmdcov_Part","Time for prg_graph2bml "//to_string(mls()-mls_ii)//" ms",lt%verbose,myRank)
      endif


    endif

    if(allocated(syprt))then
      do ipt=1,gpat%TotalParts
        call prg_destroy_subsystems(syprt(ipt),lt%verbose)
      enddo
      deallocate(syprt)
    endif

    if (myRank  ==  1 .and. lt%verbose >= 5) then
      call bml_print_matrix("gcov",g_bml,0,15,0,15)
    endif

    if(mod(mdstep,gsp2%parteach)==0.or.mdstep == 0 .or.mdstep == 1)then
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"In graph_part .."
      mls_ii = mls()
      call gpmd_graphpart()
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"Time for gpmd_graphpart "//to_string(mls()-mls_ii)//" ms"
      write(*,*)"MPI rank",myRank, "done with graph_part .."
    endif

    !To partition by molecule.
    ! write(*,*)"part by mol"
    ! call prg_molpartition(sy,nparts_cov,nl%nnStructMindist,nl%nnStruct,nl%nrnnstruct,"O ",gpat)

#ifdef SANITY_CHECK
    write(*, *) "sanity check before bml_matrix2submatrix_index"
    do ipt = 1,gpat%TotalParts
      do iipt = ipt+1,gpat%TotalParts
        do i = 1, gpat%sgraph(ipt)%llsize
          do j = 1, gpat%sgraph(iipt)%llsize
            if(gpat%sgraph(ipt)%core_halo_index(i) == gpat%sgraph(iipt)%core_halo_index(j))then
              write(*,*)"cores are repeated in partitions",mdstep
              write(*,*)ipt,gpat%sgraph(ipt)%core_halo_index(i),iipt,gpat%sgraph(ipt)%core_halo_index(j)
              write(*,*)i,j
              stop
            endif
          enddo
        enddo
      enddo
    enddo
#endif

    mls_ii = mls()
    do i=1,gpat%TotalParts
      call bml_matrix2submatrix_index(g_bml,&
           gpat%sgraph(i)%nodeInPart,gpat%nnodesInPart(i),&
           gpat%sgraph(i)%core_halo_index, &
           vsize,.true.)
      gpat%sgraph(i)%lsize = vsize(1)
      gpat%sgraph(i)%llsize = vsize(2)
    enddo

    call gpmdcov_msI("gpmdcov_Part","Time for bml_matrix2submatrix_index "//to_string(mls()-mls_ii)//" ms",lt%verbose,myRank)

    call gpmdcov_reshuffle()

#ifdef SANITY_CHECK
    write(*, *) "sanity check after bml_matrix2submatrix_index"
    do ipt = 1,gpat%TotalParts
      do iipt = ipt+1,gpat%TotalParts
        do i = 1, gpat%sgraph(ipt)%llsize
          do j = 1, gpat%sgraph(iipt)%llsize
            if(gpat%sgraph(ipt)%core_halo_index(i) == gpat%sgraph(iipt)%core_halo_index(j))then
              write(*,*)"cores are repeated in partitions",mdstep
              write(*,*)ipt,gpat%sgraph(ipt)%core_halo_index(i),iipt,gpat%sgraph(ipt)%core_halo_index(j)
              write(*,*)i,j
              stop
            endif
          enddo
        enddo
      enddo
    enddo
#endif

    if(allocated(syprt))deallocate(syprt)
    allocate(syprt(gpat%TotalParts))

    !> For every partition get the partial CH systems.
    call gpmdcov_msI("gpmdcov_Part","Getting CH subsystems ...",lt%verbose,myRank)
    call gpmdcov_msMem("gpmdcov_Part","Before prg_get_subsystem",lt%verbose,myRank)
    mls_ii = mls()

#ifdef DO_MPI
    !do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)
    do iipt=1,partsInEachRank(myRank)
      ipt= reshuffle(iipt,myRank)
#else
    do ipt = 1,gpat%TotalParts
#endif

      call prg_get_subsystem(sy,gpat%sgraph(ipt)%lsize,gpat%sgraph(ipt)%core_halo_index,syprt(ipt))
    enddo

    call gpmdcov_msI("gpmdcov_Part","Time for prg_get_subsystem "//to_string(mls()-mls_ii)//" ms",lt%verbose,myRank)
    call gpmdcov_msMem("gpmdcov","After prg_get_subsystem",lt%verbose,myRank)

    !To analyze partitions with VMD.
    if(myRank == 1)then
      if(mod(mdstep,20) == 0.or.mdstep == 2)then
        call gpmdcov_writeout()
      endif
    endif

  end subroutine gpmdcov_Part

end module gpmdcov_Part_mod
