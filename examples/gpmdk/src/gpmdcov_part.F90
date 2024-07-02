module gpmdcov_Part_mod

contains

  !> Partition by systems
  !!
  subroutine gpmdcov_Part(ipreMD)

    use gpmdcov_vars
    use gpmdcov_reshuffle_mod
    use gpmdcov_writeout_mod
    use gpmdcov_allocation_mod

    integer, allocatable :: graph_h(:,:)
    integer, allocatable :: graph_p(:,:)
    real(dp)             :: mls_ii
    real, allocatable :: onesMat(:,:)
    integer              :: iipt,ipreMD
    integer :: maxCoreHalo, minCoreHalo, averageCoreHalo
    integer :: maxCoreHaloLoc, maxCoreHaloRank
    integer :: coreHaloP1, coreP1
    integer :: myMdim

    if(gsp2%mdim < 0)then 
            myMdim = sy%nats
    elseif(gsp2%mdim > sy%nats)then 
            myMdim = sy%nats
    else 
            myMdim = gsp2%mdim
    endif

    if(ipreMD == 1)then
      call gpmdcov_msMem("gpmdcov_Part", "Before prg_get_covgraph",lt%verbose,myRank)
      call prg_get_covgraph(sy,nl%nnStruct,nl%nrnnstruct&
           ,gsp2%bml_type,gsp2%covgfact,g_bml,myMdim,lt%verbose)
      call gpmdcov_msMem("gpmdcov_Part", "After prg_get_covgraph",lt%verbose,myRank)
    else

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Anders' way of graph construction.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(gpat%TotalParts > 1)then
      call gpmdcov_msI("gpmdcov_Part","In prg_get_covgraph_h ...",lt%verbose,myRank)
      mls_ii = mls()
      call prg_get_covgraph_h(sy,nl%nnStruct,nl%nrnnstruct,gsp2%nlgcut,graph_h,myMdim,lt%verbose)
      call gpmdcov_msIII("gpmdcov_Part","In prg_get_covgraph_h ..."//to_string(mls()-mls_ii)//" ms",lt%verbose,myRank)

#ifdef DO_MPI
      !do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)
      do iipt=1,partsInEachRank(myRank)
        ipt= reshuffle(iipt,myRank)
#else
      do ipt = 1,gpat%TotalParts
#endif

        call prg_collect_graph_p(syprt(ipt)%estr%orho,gpat%sgraph(ipt)%llsize,sy%nats,syprt(ipt)%estr%hindex,&
             gpat%sgraph(ipt)%core_halo_index,graph_p,gsp2%gthreshold,myMdim,lt%verbose)

        call bml_deallocate(syprt(ipt)%estr%orho)

      enddo

      mls_i = mls()

!       call gpmdcov_mat2VectInt(graph_p,auxVectInt,sy%nats,myMdim)

#ifdef DO_MPI
      if (getNRanks() > 1) then
       call prg_sumIntReduceN(graph_p, myMdim*sy%nats)
 !      call prg_sumIntReduceN(auxVectInt, myMdim*sy%nats)
      endif
#endif

 !     call gpmdcov_vect2MatInt(auxVectInt,graph_p,sy%nats,myMdim)
 !     deallocate(auxVectInt)
!      write(*,*)graph_p
      call gpmdcov_msIII("gpmdcov_Part","Time for prg_sumIntReduceN for graph "//to_string(mls() - mls_i)//" ms",lt%verbose,myRank)

      call gpmdcov_msI("gpmdcov_Part","In prg_merge_graph ...",lt%verbose,myRank)
      mls_ii = mls()
      call prg_wait()
      call prg_merge_graph(graph_p,graph_h)
      call prg_wait()
      call gpmdcov_msIII("gpmdcov_Part","Time for prg_merge_graph "//to_string(mls()-mls_ii)//" ms",lt%verbose,myRank)

      !call prg_wait()
      !call prg_wait()

      !Transform graph into bml format.
  !    if(mod(mdstep,gsp2%parteach)==0 .or. mdstep <= 1)then !Only for debug purposes (Halos are updated every step)
        if(bml_allocated(g_bml)) call bml_deallocate(g_bml)
        !call bml_zero_matrix(gsp2%bml_type,bml_element_real,kind(1.0),sy%nats,myMdim,g_bml,lt%bml_dmode)
        call bml_zero_matrix(gsp2%bml_type,bml_element_real,kind(1.0),sy%nats,myMdim,g_bml)
        !call bml_zero_matrix(gsp2%bml_type,bml_element_real,dp,sy%nats,myMdim,g_bml)

        call gpmdcov_msMem("gpmdcov_Part","Before prg_graph2bml",lt%verbose,myRank)
        mls_ii = mls()
        call prg_graph2bml(graph_p,gsp2%bml_type,g_bml)
!        if(lt%verbose == 7 .and. myRank == 1)then 
!                call bml_write_matrix(g_bml,"g_bml.mtx")
!                stop
!        endif


        if(allocated(graph_p)) deallocate(graph_p)
        if(allocated(graph_h)) deallocate(graph_h)


        call gpmdcov_msMem("gpmdcov_Part","After prg_graph2bml",lt%verbose,myRank)
        call gpmdcov_msIII("gpmdcov_Part","Time for prg_graph2bml "//to_string(mls()-mls_ii)//" ms",lt%verbose,myRank)
     !endif
      else
        allocate(onesMat(sy%nats,sy%nats))
        onesMat = 1.0
        call bml_zero_matrix(gsp2%bml_type,bml_element_real,kind(1.0),sy%nats,myMdim,g_bml)
        call bml_import_from_dense(gsp2%bml_type, onesMat, g_bml, 0.0_dp, sy%nats)
        deallocate(onesMat)
      endif

    endif

    if(allocated(syprt))then
      do ipt=1,gpat%TotalParts
        call prg_destroy_estr(syprt(ipt)%estr)
      enddo
      
      do ipt=1,gpat%TotalParts
        call prg_destroy_subsystems(syprt(ipt),lt%verbose)
      enddo
      deallocate(syprt)
    endif

    if (myRank  ==  1 .and. lt%verbose >= 5) then
      call bml_print_matrix("gcov",g_bml,0,15,0,15)
    endif

    if(mod(mdstep,gsp2%parteach)==0 .or. mdstep <= 1)then
      if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"In graph_part .."
      mls_ii = mls()
      call gpmd_graphpart()
      firstKernel = .true.
      if(mdstep >= 1) newPart = .true.
      if(lt%verbose >= 3 .and. myRank == 1)write(*,*)"Time for gpmd_graphpart "//to_string(mls()-mls_ii)//" ms"
    endif

   !> \todo have the preconditioner construction independent of the partitioning
   !if(mod(mdstep,10)==0 .or. mdstep <= 1)then
   !   if(lt%verbose >= 1 .and. myRank == 1)write(*,*)"In graph_part .."
   !   firstKernel = .true.
   !endif



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
      if(myRank == 1 .and. lt%verbose == 3) write(*,*)"part",i,"cores, cores+halo",vsize(2),vsize(1)
    enddo

    if(myRank == 1)then
        maxCoreHalo = maxval(gpat%sgraph(:)%lsize)
        maxCoreHaloLoc = maxloc(gpat%sgraph(:)%lsize,dim=1)
        maxCoreHaloRank = 0
        do i = 1,gpat%totalProcs
           if (maxCoreHaloLoc.ge.gpat%localPartMin(i).and.maxCoreHaloLoc.le.gpat%localPartMax(i)) then
              maxCoreHaloRank = i
           endif
        enddo
        minCoreHalo = minval(gpat%sgraph(:)%lsize)
        if(gpat%TotalParts.eq.0)then
           averageCoreHalo = sum(gpat%sgraph(:)%lsize)
        else
           averageCoreHalo = sum(gpat%sgraph(:)%lsize)/gpat%TotalParts
        endif
        call gpmdcov_msI("gpmdcov_Part","Max and min core+halo "//to_string(maxCoreHalo)//" "//&
                &to_string(minCoreHalo),lt%verbose,myRank)
        call gpmdcov_msI("gpmdcov_Part","Max core+halo in (part rank) "//to_string(maxCoreHaloLoc)//" "//&
                &to_string(maxCoreHaloRank),lt%verbose,myRank)
        call gpmdcov_msI("gpmdcov_Part","Average core+halo "//to_string(averageCoreHalo)//" ",1,myRank)
        if(gpmdt%tracknparts > 0)then 
           do i = 1,gpmdt%tracknparts    
             coreHaloP1 = gpat%sgraph(gpmdt%trackparts(i))%lsize
             coreP1 = gpat%sgraph(gpmdt%trackparts(i))%llsize
             call gpmdcov_msI("gpmdcov_Part","Part "//to_string(gpmdt%trackparts(i))//&
                     &" core+halo "//to_string(coreHaloP1)//" ",1,myRank)
             call gpmdcov_msI("gpmdcov_Part","Part "//to_string(gpmdt%trackparts(i))//&
                     &" core "//to_string(coreP1)//" ",1,myRank)
           enddo
        endif
     endif

    call gpmdcov_msIII("gpmdcov_Part","Time for bml_matrix2submatrix_index "//to_string(mls()-mls_ii)//" ms",lt%verbose,myRank)

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
    !For currents
    if(gpmdt%compcurr)then 
      if(MDStep == 0 .or. mod(mdstep,gsp2%parteach)==0)then 
        if(allocated(oldsyprt))deallocate(oldsyprt)
        allocate(oldsyprt(gpat%TotalParts))
      endif
    endif

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

    call gpmdcov_msIII("gpmdcov_Part","Time for prg_get_subsystem "//to_string(mls()-mls_ii)//" ms",lt%verbose,myRank)
    call gpmdcov_msMem("gpmdcov","After prg_get_subsystem",lt%verbose,myRank)

    !To analyze partitions with VMD.
    if(lt%verbose >= 4)then
      if(mod(mdstep,gsp2%parteach) == 0 .or. mdstep == 1)then
        call gpmdcov_writepartitionout(sy,syprt,gpat,reshuffle,partsInEachRank,myRank)
      endif
    endif    


  end subroutine gpmdcov_Part

end module gpmdcov_Part_mod
