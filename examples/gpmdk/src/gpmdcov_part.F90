module gpmdcov_Part_mod

contains

  !> Partition by systems
  !!
  subroutine gpmdcov_Part(ipreMD)

    use gpmdcov_vars
    use gpmdcov_reshuffle_mod
    use gpmdcov_writeout_mod
    use gpmdcov_allocation_mod
    implicit none
    integer, allocatable :: graph_h(:,:)
    integer, allocatable, save :: graph_p(:,:), graph_p_flat(:)
    integer, allocatable, save :: graph_p_old(:,:)
    integer, allocatable, save :: G_added(:,:), G_removed(:,:)
    integer, allocatable, save :: N_added(:), N_removed(:), NNZ1(:), NNZ2(:), NNZ_updated(:)
    logical, allocatable, save :: v(:), v_check(:)
    integer :: n_atoms, max_updates, k, ktot_a, ktot_r
    real(dp)             :: mls_ii
    real, allocatable :: onesMat(:,:)
    integer                    ::iipt
    integer, intent(in)        :: ipreMD
    integer :: maxCoreHalo, minCoreHalo, averageCoreHalo
    integer :: maxCoreHaloLoc, maxCoreHaloRank
    integer :: coreHaloP1, coreP1
    integer :: myMdim
    logical :: check_chi,check_graph,graphs_end

    if(gsp2%mdim < 0)then 
       myMdim = sy%nats
    elseif(gsp2%mdim > sy%nats)then 
       myMdim = sy%nats
    else 
       myMdim = gsp2%mdim
    endif

    if(.not.allocated(graph_p))write(*,*)"GPMDCOV_PART: graph_p not allocated upon entry"
    if(ipreMD == 1)then
       write(*,*)"DEBUG: Doing ipreMD graph partitioning at mdstep ",mdstep

       call gpmdcov_msMem("gpmdcov_Part", "Before prg_get_covgraph",lt%verbose,myRank)
       call prg_get_covgraph(sy,nl%nnStruct,nl%nrnnstruct&
            ,gsp2%bml_type,gsp2%covgfact,g_bml,myMdim,lt%verbose)


       call gpmdcov_msMem("gpmdcov_Part", "After prg_get_covgraph",lt%verbose,myRank)
    else !ipreMD == 1
#ifdef DO_MPI
       n_atoms = sy%nats
       max_updates = 100
       if(.not.allocated(graph_p_old))then
          allocate(graph_p(myMdim,n_atoms))
          graph_p = 0
          allocate(graph_p_old(myMdim,n_atoms))
          graph_p_old = 0
          allocate(G_added(max_updates,n_atoms))
          allocate(G_removed(max_updates,n_atoms))
          allocate(N_added(n_atoms))
          allocate(N_removed(n_atoms))
          allocate(NNZ1(n_atoms))
          allocate(NNZ2(n_atoms))
          allocate(NNZ_updated(n_atoms))
          allocate(v(n_atoms))
          allocate(v_check(n_atoms))
       endif
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Anders' way of graph construction.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       if(gpat%TotalParts > 1)then
          if(.not.gsp2%small_subgraphs)then
             call gpmdcov_msI("gpmdcov_Part","In prg_get_covgraph_h ...",lt%verbose,myRank)
             mls_ii = mls()
             call prg_get_covgraph_h(sy,nl%nnStruct,nl%nrnnstruct,gsp2%nlgcut,graph_h,myMdim,lt%verbose)
             call gpmdcov_msII("gpmdcov_Part","In prg_get_covgraph_h ..."//to_string(mls()-mls_ii)//" ms",lt%verbose,myRank)
             graph_p = 0
          endif
#ifdef DO_MPI
      !do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)
          do iipt=1,partsInEachRank(myRank)
             ipt= reshuffle(iipt,myRank)
#else
          do ipt = 1,gpat%TotalParts
#endif
             if(gsp2%small_subgraphs)then
                call prg_collect_extended_graph_p(syprt(ipt)%estr%orho,gpat%sgraph(ipt)%llsize,sy%nats,syprt(ipt)%estr%hindex,&
                     gpat%sgraph(ipt)%core_halo_index,graph_p,gsp2%gthreshold,myMdim,gsp2%alpha,syprt(ipt)%coordinate,sy%coordinate,sy%lattice_vector,lt%verbose)
             else
                call prg_collect_graph_p(syprt(ipt)%estr%orho,gpat%sgraph(ipt)%llsize,sy%nats,syprt(ipt)%estr%hindex,&
                     gpat%sgraph(ipt)%core_halo_index,graph_p,gsp2%gthreshold,myMdim,lt%verbose)
             endif
        
             call bml_deallocate(syprt(ipt)%estr%orho)

          enddo
          
          mls_i = mls()

!       call gpmdcov_mat2VectInt(graph_p,auxVectInt,sy%nats,myMdim)
          if(.not.allocated(graph_p_flat))allocate(graph_p_flat(myMdim*sy%nats))
#ifdef DO_MPI
          if (getNRanks() > 1) then
             call prg_barrierParallel
             if((gsp2%parteach == 1) .or. (mod(mdstep,gsp2%parteach)==0) .or. (mdstep <= 1))then
                !graph_p_flat = RESHAPE(graph_p,shape(graph_p_flat))
                !call prg_sumIntReduceN(graph_p_flat, size(graph_p_flat))
                !graph_p = RESHAPE(graph_p_flat,shape(graph_p))
                call prg_sumIntReduceN(graph_p, myMdim*sy%nats)
                 write(*,*)"DEBUG: Doing full graph reduction at mdstep ",mdstep
                graph_p_old = graph_p
             else
                write(*,*)"DEBUG: Doing graph update reduction at mdstep ",mdstep

                ktot_a = 0
                NNZ1 = count(graph_p_old.ne.0,DIM=1)
                NNZ2 = count(graph_p.ne.0,DIM=1)
                
                G_added = 0
                N_added = 0
                v = .false.
                do iipt=1,partsInEachRank(myRank)
                   ipt= reshuffle(iipt,myRank)
                   do ii = 1,gpat%sgraph(ipt)%llsize
                      i = gpat%sgraph(ipt)%core_halo_index(ii) + 1
                      do j = 1,NNZ1(i)
                         v(graph_p_old(j,i)) = .true.
                      end do
                      k = 0;
                      do j = 1,NNZ2(i)
                         if (v(graph_p(j,i)) .eqv. .false.)then
                            k = k + 1
                            G_added(k,i) = graph_p(j,i)
                         endif
                      end do
                      N_added(i) = k  ! Number of added edges for each vertex i
                      if(ktot_a.lt.k)then
                         ktot_a = k
                      endif
                      v(graph_p_old(1:NNZ1(i),i)) = .false.
                      v(graph_p(1:NNZ2(i),i)) = .false.
                   end do
                enddo
                ! Removed edges
                ktot_r = 0
                G_removed = 0
                N_removed = 0
                v = .false.
                do iipt=1,partsInEachRank(myRank)
                   ipt= reshuffle(iipt,myRank)
                   do ii = 1,gpat%sgraph(ipt)%llsize
                      i = gpat%sgraph(ipt)%core_halo_index(ii) + 1
                      do j = 1,NNZ2(i)
                         v(graph_p(j,i)) = .true.
                      end do
                      k = 0;
                      do j = 1,NNZ1(i)
                         if (v(graph_p_old(j,i)) .eqv. .false.)then
                            k = k + 1
                            G_removed(k,i) = graph_p_old(j,i)
                         endif
                      end do
                      N_removed(i) = k  ! Number of added edges for each vertex i
                      if(ktot_r.lt.k)then
                         ktot_r = k
                      endif
                      v(graph_p_old(1:NNZ1(i),i)) = .false.
                      v(graph_p(1:NNZ2(i),i)) = .false.
                   end do
                enddo
                ! % Check NNZ_Updated: NNZ_Updated = NNZ1 + N_Added - N_Removed
                ! [NNZ_Updated, NNZ1 + N_added - N_removed]
                
                ! % Check NNZ_Updated: G_Updated = G2 But edges are not in the same order
                call prg_maxIntReduce2(ktot_a,ktot_r)
                write(*,*)"DEBUG: max ktot_a = ",ktot_a
                write(*,*)"DEBUG: max ktot_r = ",ktot_r
                if((ktot_a<=max_updates).and.(ktot_r<=max_updates))then
                   ! %% Use G_removed and G_added to update from G1 to G2
                   call prg_sumIntReduceN(G_added(1:ktot_a,:),n_atoms*ktot_a)
                   call prg_sumIntReduceN(G_removed(1:ktot_r,:),n_atoms*ktot_r)
                   NNZ_updated = 0
                   v = .false.
                   v_check = .false.
                   !$omp parallel do &
                   !$omp default(none) &
                   !$omp private(i,j,k) &
                   !$omp firstprivate(v,v_check) &
                   !$omp shared(G_added,G_removed,N_added,N_removed) &
                   !$omp shared(graph_p_old,graph_p,NNZ1,NNZ2,NNZ_updated,n_atoms,mymdim)
                   do i = 1,n_atoms
                      v = .false.
                      !$omp loop
                      do j = 1,NNZ1(i)
                         v(graph_p_old(j,i)) = .true.   
                      end do
                      !$omp loop
                      do j = 1,N_removed(i)
                         v(G_removed(j,i)) = .false.  ! % Remove edges
                      end do
                      k = 0
!!$omp loop reduction(+:k)
                      !$omp loop
                      do j = 1,mymdim
                         graph_p(j,i) = 0
                      enddo
                      do j = 1,NNZ1(i)
                         if (v(graph_p_old(j,i)) .eqv. .true.)then ! % Account only for the remaining edges  
                            k = k + 1;
                            graph_p(k,i) = graph_p_old(j,i);
                         end if
                      end do
                      NNZ_updated(i) = k + N_added(i)
                      !$omp loop
                      do j = k+1,NNZ_updated(i)
                         graph_p(j,i) = G_added(j-k,i) ! Add new edges at the end
                      end do
                      k = max(NNZ1(i),NNZ2(i))
                      graph_p_old(1:k,i) = graph_p(1:k,i)
                   end do
                   !$omp end parallel do
                else
                   write(*,*)"GPMDCOV_PART: WARNING: Number of changes exceeds max_updates. System might be unstable. Doing full reduction."
                   call prg_sumIntReduceN(graph_p, myMdim*sy%nats)
                   graph_p_old = graph_p
                endif
             endif
             !      call prg_sumIntReduceN(auxVectInt, myMdim*sy%nats)
          endif
#endif
          !     call gpmdcov_vect2MatInt(auxVectInt,graph_p,sy%nats,myMdim)
          !     deallocate(auxVectInt)
          !      write(*,*)graph_p
          call gpmdcov_msII("gpmdcov_Part","Time for prg_sumIntReduceN for graph "//to_string(mls() - mls_i)//" ms",lt%verbose,myRank)
          if(.not.gsp2%small_subgraphs)then
             call gpmdcov_msI("gpmdcov_Part","In prg_merge_graph ...",lt%verbose,myRank)
             mls_ii = mls()
             call prg_wait()
             call prg_merge_graph(graph_p,graph_h)
             call prg_wait()
             call gpmdcov_msII("gpmdcov_Part","Time for prg_merge_graph "//to_string(mls()-mls_ii)//" ms",lt%verbose,myRank)
          endif

          !call prg_wait()
          !call prg_wait()

      !Transform graph into bml format.
  !    if(mod(mdstep,gsp2%parteach)==0 .or. mdstep <= 1)then !Only for debug purposes (Halos are updated every step)
          !if(bml_allocated(g_bml)) call bml_deallocate(g_bml)
          !call bml_zero_matrix(gsp2%bml_type,bml_element_real,kind(1.0),sy%nats,myMdim,g_bml,lt%bml_dmode)
          if(.not.bml_allocated(g_bml))then
             call bml_zero_matrix(gsp2%bml_type,bml_element_real,kind(1.0),sy%nats,myMdim,g_bml)
          endif
          !call bml_zero_matrix(gsp2%bml_type,bml_element_real,dp,sy%nats,myMdim,g_bml)

          call gpmdcov_msMem("gpmdcov_Part","Before prg_graph2bml",lt%verbose,myRank)
          mls_ii = mls()
          call prg_graph2bml(graph_p,gsp2%bml_type,g_bml)
          !        if(lt%verbose == 7 .and. myRank == 1)then 
          !                call bml_write_matrix(g_bml,"g_bml.mtx")
          !                stop
          !        endif
          
          
          !if(allocated(graph_p)) deallocate(graph_p)
          if(allocated(graph_h)) deallocate(graph_h)
          
          
          call gpmdcov_msMem("gpmdcov_Part","After prg_graph2bml",lt%verbose,myRank)
          call gpmdcov_msII("gpmdcov_Part","Time for prg_graph2bml "//to_string(mls()-mls_ii)//" ms",lt%verbose,myRank)
          !endif
       else
          allocate(onesMat(sy%nats,sy%nats))
          onesMat = 1.0
          if(.not.bml_allocated(g_bml))then
             call bml_zero_matrix(gsp2%bml_type,bml_element_real,kind(1.0),sy%nats,myMdim,g_bml)
          endif
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
       if(lt%verbose >= 2 .and. myRank == 1)write(*,*)"Time for gpmd_graphpart "//to_string(mls()-mls_ii)//" ms"
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
            vsize,(ipremd==1).or.(.not.gsp2%small_subgraphs))
       gpat%sgraph(i)%lsize = vsize(1)
       gpat%sgraph(i)%llsize = vsize(2)
       if(myRank == 1 .and. lt%verbose == 3) write(*,*)"part",i,"cores, cores+halo",vsize(2),vsize(1)
       check_chi=.true.
       if(check_chi)then
          do j=1,vsize(1)
             if ((gpat%sgraph(i)%core_halo_index(j)+1).gt.sy%nats)then
                write(*,*)"GPMDCOV_PART: ERROR core_halo_index(",j,") + 1 exceeds nats = ",sy%nats
             endif
          enddo
       endif
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

    call gpmdcov_msII("gpmdcov_Part","Time for bml_matrix2submatrix_index "//to_string(mls()-mls_ii)//" ms",lt%verbose,myRank)
    
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

    if(.not.allocated(graph_p))write(*,*)"GPMDCOV_PART: graph_p not allocated on exit"
    
  end subroutine gpmdcov_Part
  
end module gpmdcov_Part_mod
