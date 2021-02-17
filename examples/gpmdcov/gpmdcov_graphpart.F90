
  !> Graph partitioning subroutine
  subroutine gpmd_graphpart
    use gpmdcov_vars

    integer :: tnnz

    !     if(.not.allocated(xadj))then
    !     !> Symmetrize and Threshold the Matrix
    !     if(gsp2%mdim > 0)then
    !       mdim = gsp2%mdim
    !     else
    !       mdim = sy%nats
    !     endif

    !     call bml_write_matrix(g_bml,"g_bml_bef")
    !     call bml_zero_matrix(gsp2%bml_type,bml_element_real,kind(1.0),sy%nats,mdim,copy_g_bml)
    !     call bml_threshold(g_bml, gsp2%gthreshold)
    !     call bml_transpose(g_bml, copy_g_bml)
    !     call bml_add_deprecated(0.5_dp,g_bml,0.5_dp,copy_g_bml,0.0_dp)
    !     call bml_threshold(g_bml, gsp2%gthreshold)
    !     call bml_deallocate(copy_g_bml)
    !     call bml_write_matrix(g_bml,"g_bml_aft")
    ! #ifdef DO_MPI_BLOCK
    !     call prg_allGatherParallel(g_bml)
    ! #endif

    !     endif

    !> Create graph partitioning - Use Block or METIS or METIS+SA or METIS+KL
    if(lt%verbose >= 1) call prg_timer_start(part_timer)

    !> Block partitioning
    if (gsp2%partition_type == "Block") then
      !> Partition by orbital or atom
      if (gsp2%graph_element == "Orbital") then
        call prg_equalPartition(gpat, gsp2%nodesPerPart, sy%nats)
      else
        call prg_equalGroupPartition(gpat, hindex, nnodes, gsp2%nodesPerPart, sy%nats)
      endif

      !> METIS, METIS+SA, or METIS+KL partitioning
    else
      if(.not.allocated(xadj))then

        allocate(xadj(sy%nats+1))

        tnnz = 0  !This must be done at the bml level
        do i=1,sy%nats
          tnnz = tnnz + bml_get_row_bandwidth(g_bml,i)
        enddo

        !       allocate(adjncy(sy%nats*sy%nats)) !Old way
        allocate(adjncy(tnnz+1))

        call bml_adjacency(g_bml, xadj, adjncy, 1)

        !         call bml_deallocate(g_bml)

        !          call prg_sortadj(xadj, adjncy)  !Not needed anymore since the g_bml is sorted.

        nnodes = sy%nats
      else
        nnodes = sy%nats
      endif

#ifdef DO_GRAPHLIB
      nparts = gsp2%partition_count

      if (first_part) then
        allocate(part(nnodes))
        allocate(core_count(nparts))
      endif

      ! if(allocated(CH_count)) deallocate(CH_count)
      ! if(allocated(Halo_count)) deallocate(Halo_count)
      allocate(CH_count(nparts))
      allocate(Halo_count(nparts, nnodes))

      !> For METIS, if no refinement of first time, do full partitioning
      if (gsp2%partition_refinement == 'None' .or. first_part) then

        if(allocated(gpat%nnodesInPart))then
          call prg_destroyGraphPartitioning(gpat)
        endif
        !> Which METIS partitioning
        select case(gsp2%partition_type)
        case("METIS")
          call prg_metisPartition(gpat, sy%nats, sy%nats, xadj, adjncy, nparts, part, core_count,&
               CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)
        case("METIS+SA")
          call prg_metisPartition(gpat, sy%nats, sy%nats, xadj, adjncy, nparts, part, core_count,&
               CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)
          call prg_simAnnealing(gpat, xadj, adjncy, part, core_count, CH_count, &
               Halo_count, sumCubes, maxCH,smooth_maxCH,pnorm, niter, seed)
        case("METIS+KL")
          call prg_metisPartition(gpat, sy%nats, sy%nats, xadj, adjncy, nparts, part, core_count,&
               CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)
          call prg_KernLin2(gpat, xadj, adjncy, part, core_count, CH_count, &
               Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)
        case default
          write(*,*)"No METIS partitioning specified"
          stop ""
        end select

        first_part = .false.

        !> Not first time, do refinement
      else
        if (.not. first_part) then

          !> Which refinement
          select case(gsp2%partition_refinement)
          case("SA")
            call prg_simAnnealing(gpat, xadj, adjncy, part, core_count, CH_count, &
                 Halo_count, sumCubes, maxCH,smooth_maxCH,pnorm, niter, seed)
          case("KL")
            call prg_KernLin2(gpat, xadj, adjncy, part, core_count, CH_count, &
                 Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)
          case default
            write(*,*)"No refinement specified"
            stop ""
          end select

        endif
      endif

      deallocate(xadj)
      deallocate(adjncy)
      deallocate(CH_count)
      deallocate(Halo_count)

#endif
    endif

    if(lt%verbose >= 1) call prg_timer_stop(part_timer)

  end subroutine gpmd_graphpart

