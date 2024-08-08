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
    integer, allocatable :: graph_p(:,:)
    integer, allocatable, save :: graph_p_old(:,:)
    integer, allocatable, save :: G_added(:,:), G_removed(:,:), G_updated(:,:)
    integer, allocatable, save :: N_added(:), N_removed(:), NNZ1(:), NNZ2(:), NNZ_updated(:)
    logical, allocatable, save :: v(:)
    integer :: na, usz, k, ktot
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
#ifdef DO_MPI
       na = sy%nats
       usz = 100
    if(.not.allocated(graph_p_old))then
       allocate(graph_p_old(myMdim,na))
       graph_p_old = 0
       allocate(v(na))
       allocate(G_added(usz,na))
       allocate(G_removed(usz,na))
       allocate(G_updated(usz,na))
       allocate(N_added(na))
       allocate(N_removed(na))
       allocate(NNZ1(na))
       allocate(NNZ2(na))
       allocate(NNZ_updated(na))
    endif
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Anders' way of graph construction.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(gpat%TotalParts > 1)then
      call gpmdcov_msI("gpmdcov_Part","In prg_get_covgraph_h ...",lt%verbose,myRank)
      mls_ii = mls()
      call prg_get_covgraph_h(sy,nl%nnStruct,nl%nrnnstruct,gsp2%nlgcut,graph_h,myMdim,lt%verbose)
      call gpmdcov_msII("gpmdcov_Part","In prg_get_covgraph_h ..."//to_string(mls()-mls_ii)//" ms",lt%verbose,myRank)

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
         call prg_barrierParallel
         if((gsp2%parteach == 1) .or. (mod(mdstep,gsp2%parteach)==1) .or. (mdstep <= 1))then 
            call prg_sumIntReduceN(graph_p, myMdim*sy%nats)
            write(*,*)"DEBUG: Doing full graph reduction at mdstep ",mdstep
            graph_p_old = graph_p
         else
            !MATLAB code
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! %% Detects new and removed edges between two graphs G1 and G2 %%
            ! %%         The routine does not require ordering              %%
            ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! clear;
            ! N = 10; M = 5;  % N nodes and max M edges

            ! % Create first graph randomly
            ! G1 = zeros(N,M); 
            ! NNZ1 = floor(1 + M*rand(10,1));  % [1,M] for each node
            ! TIO = [1:N]';
            ! for i = 1:N
            !   for k = 1:7  % reshuffle ordering
            !     a = floor(1 + N*rand(1));
            !     b = floor(1 + N*rand(1));
            !     tmp = TIO(a); TIO(a) = TIO(b); TIO(b) = tmp;
            !   end
            !   for j = 1:NNZ1(i)
            !     G1(i,j) = TIO(j);
            !   end
            ! end

            ! % Create second graph randomly
            ! G2 = zeros(N,M); 
            ! NNZ2 = floor(1 + M*rand(10,1));  % [1,M] for each node
            ! TIO = [1:N]';
            ! for i = 1:N
            !   for k = 1:7  % reshuffle ordering
            !     a = floor(1 + N*rand(1));
            !     b = floor(1 + N*rand(1));
            !     tmp = TIO(a); TIO(a) = TIO(b); TIO(b) = tmp;
            !   end
            !   for j = 1:NNZ2(i)
            !     G2(i,j) = TIO(j);
            !   end
            ! end

            ! % Analyze the difference from G1 to G2

            ! %  Added edges
            ! G_added = zeros(N,M); N_added = zeros(N,1);
            ! v = zeros(1,N);
            ! for i = 1:N
            !   for j = 1:NNZ1(i)
            !     v(G1(i,j)) = 1;
            !   end
            !   k = 0;
            !   for j = 1:NNZ2(i)
            !     if (v(G2(i,j)) == 0)
            !       k = k + 1;
            !       G_added(i,k) = G2(i,j);
            !     end
            !   end
            !   N_added(i) = k;  % Number of added edges for each vertex i
            !   v(G1(i,1:NNZ1(i))) = 0;
            !   v(G2(i,1:NNZ2(i))) = 0;
            ! end

            ! % Removed edges
            ! G_removed = zeros(N,M); N_removed = zeros(N,1);
            ! v = zeros(1,N);
            ! for i = 1:N
            !   for j = 1:NNZ2(i)
            !     v(G2(i,j)) = 1;
            !   end
            !   k = 0;
            !   for j = 1:NNZ1(i)
            !     if (v(G1(i,j)) == 0)
            !       k = k + 1;
            !       G_removed(i,k) = G1(i,j);
            !     end
            !   end
            !   N_removed(i) = k; % Number of removed edges for each vertex i
            !   v(G1(i,1:NNZ1(i))) = 0;
            !   v(G2(i,1:NNZ2(i))) = 0;
            ! end

            ! %% Use G_removed and G_added to update from G1 to G2
            ! G_Updated = zeros(N,M);
            ! NNZ_Updated = zeros(N,1);
            ! v = zeros(N,1); % Temporary vector that keeps track of elements that are there and then removed.
            ! for i = 1:N
            !   for j = 1:NNZ1(i)
            !     v(G1(i,j)) = 1;   
            !   end
            !   for j = 1:N_removed(i)
            !     v(G_removed(i,j)) = 0;  % Remove edges
            !   end
            !   cnt = 0;
            !   for j = 1:NNZ1(i)
            !     if v(G1(i,j)) > 0; % Account only for the remaining edges  
            !       cnt = cnt + 1;
            !       G_Updated(i,cnt) = G1(i,j);
            !     end
            !     NNZ_Updated(i) = cnt + N_added(i);
            !   end
            !   for j = cnt+1:NNZ_Updated(i)
            !     G_Updated(i,j) = G_added(i,j-cnt); % Add new edges at the end
            !   end
            ! end

            ! % Check NNZ_Updated: NNZ_Updated = NNZ1 + N_Added - N_Removed
            ! [NNZ_Updated, NNZ1 + N_added - N_removed]

            ! % Check NNZ_Updated: G_Updated = G2 But edges are not in the same order

            ! %  Added edges
            ktot = 0
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
                  ktot = ktot + k
                  v(graph_p_old(1:NNZ1(i),i)) = .false.
                  v(graph_p(1:NNZ2(i),i)) = .false.
               end do
            enddo
            ! Removed edges
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
                  ktot = ktot + k
                  v(graph_p_old(1:NNZ1(i),i)) = .false.
                  v(graph_p(1:NNZ2(i),i)) = .false.
               end do
            enddo
            write(*,*)"DEBUG: ktot = ",ktot
            call prg_sumIntReduceN(graph_p, myMdim*sy%nats)
            write(*,*)"DEBUG: Doing graph update reduction at mdstep ",mdstep
                        ! %% Use G_removed and G_added to update from G1 to G2
            G_updated = 0
            NNZ_updated = 0
            v = .false. ! % Temporary vector that keeps track of elements that are there and then removed.
            do iipt=1,partsInEachRank(myRank)
               ipt= reshuffle(iipt,myRank)
               do ii = 1,gpat%sgraph(ipt)%llsize
                  i = gpat%sgraph(ipt)%core_halo_index(ii) + 1
                  do j = 1,NNZ1(i)
                     v(graph_p_old(j,i)) = .true.   
                  end do
                  do j = 1,N_removed(i)
                     v(G_removed(j,i)) = .false.  ! % Remove edges
                  end do
                  k = 0
                  do j = 1,NNZ1(i)
                     if (v(graph_p_old(j,i)) .eqv. .true.)then ! % Account only for the remaining edges  
                        k = k + 1;
                        G_updated(k,i) = graph_p_old(j,i);
                     end if
                     NNZ_updated(i) = k + N_added(i);
                  end do
                  do j = k+1,NNZ_updated(i)
                     G_updated(j,i) = G_added(j-k,i) ! Add new edges at the end
                  end do
               end do
            enddo
            ! % Check NNZ_Updated: NNZ_Updated = NNZ1 + N_Added - N_Removed
            ! [NNZ_Updated, NNZ1 + N_added - N_removed]

            ! % Check NNZ_Updated: G_Updated = G2 But edges are not in the same order

            graph_p_old = graph_p
         endif
         !      call prg_sumIntReduceN(auxVectInt, myMdim*sy%nats)
      endif
#endif

 !     call gpmdcov_vect2MatInt(auxVectInt,graph_p,sy%nats,myMdim)
 !     deallocate(auxVectInt)
!      write(*,*)graph_p
      call gpmdcov_msII("gpmdcov_Part","Time for prg_sumIntReduceN for graph "//to_string(mls() - mls_i)//" ms",lt%verbose,myRank)

      call gpmdcov_msI("gpmdcov_Part","In prg_merge_graph ...",lt%verbose,myRank)
      mls_ii = mls()
      call prg_wait()
      call prg_merge_graph(graph_p,graph_h)
      call prg_wait()
      call gpmdcov_msII("gpmdcov_Part","Time for prg_merge_graph "//to_string(mls()-mls_ii)//" ms",lt%verbose,myRank)

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
        call gpmdcov_msII("gpmdcov_Part","Time for prg_graph2bml "//to_string(mls()-mls_ii)//" ms",lt%verbose,myRank)
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


  end subroutine gpmdcov_Part

end module gpmdcov_Part_mod
