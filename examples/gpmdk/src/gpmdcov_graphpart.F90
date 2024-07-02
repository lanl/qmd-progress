
  !> Graph partitioning subroutine
  subroutine gpmd_graphpart
    use gpmdcov_vars
    use gpmdcov_writeout_mod

    integer :: tnnz

    !> Create graph partitioning - Use Block or METIS or METIS+SA or METIS+KL
    if(lt%verbose >= 1) call prg_timer_start(part_timer)

    !> Block partitioning
    if (gsp2%partition_type == "Block") then
      !> Partition by orbital or atom
      if (gsp2%graph_element == "Atom") then
        if(gsp2%partition_count == 1)then 
          gpat%TotalParts = 1
          gsp2%nodesPerPart = sy%nats
        endif
        call prg_equalPartition(gpat, gsp2%nodesPerPart, sy%nats)
      else
        call prg_equalGroupPartition(gpat, hindex, nnodes, gsp2%nodesPerPart, sy%nats)
      endif

      !> METIS, METIS+SA, or METIS+KL partitioning
    else

#ifdef DO_GRAPHLIB

      if(.not.allocated(xadj))then

        allocate(xadj(sy%nats+1))

        tnnz = 0  !This must be done at the bml level
        do i=1,sy%nats
          tnnz = tnnz + bml_get_row_bandwidth(g_bml,i)
        enddo

        allocate(adjncy(tnnz+1))

        call bml_adjacency(g_bml, xadj, adjncy, 1)

        nnodes = sy%nats
      else
        nnodes = sy%nats
      endif

      nparts = gsp2%partition_count

      if (first_part) then
        allocate(part(nnodes))
        allocate(core_count(nparts))
      endif

      if(allocated(CH_count)) deallocate(CH_count)
      if(allocated(Halo_count)) deallocate(Halo_count)
      allocate(CH_count(nparts))
      allocate(Halo_count(nparts, nnodes))

      !> For METIS, if no refinement of first time, do full partitioning
      if (gsp2%partition_refinement == 'None' .or. first_part .or. mdstep <= 1) then
        call gpmdcov_message("graph_part","Will do METIS partition",lt%verbose,2,myrank)

        call prg_destroyGraphPartitioning(gpat)

        !> Which METIS partitioning
        select case(gsp2%partition_type)
       
        case("Block")
           if (gsp2%graph_element == "Atom") then
              call prg_equalPartition(gpat, gsp2%nodesPerPart, sy%nats)
           else
              call prg_equalGroupPartition(gpat, hindex, nnodes, gsp2%nodesPerPart, sy%nats)
           endif

        case("METIS")
          call prg_metisPartition(gpat, sy%nats, sy%nats, xadj, adjncy, nparts, part, core_count,&
               CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)

        case("METIS+SA")

          call prg_metisPartition(gpat, sy%nats, sy%nats, xadj, adjncy, nparts, part, core_count,&
               CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)

          seed = mdstep 
          call prg_simAnnealing(gpat, xadj, adjncy, part, core_count, CH_count, &
               Halo_count, sumCubes, maxCH,smooth_maxCH,pnorm,1,seed)

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
        call gpmdcov_message("graph_part","In refinment",lt%verbose,1,myrank)
        if (.not. first_part) then

          !> Which refinement
          select case(gsp2%partition_refinement)
          case("SA")
            call prg_simAnnealing(gpat, xadj, adjncy, part, core_count, CH_count, &
                 Halo_count, sumCubes, maxCH,smooth_maxCH,pnorm, 10, seed)
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
#else
      write(*,*)"ERROR: METIS option will only work if the code is compiled with GRAPHLIB"
#endif 

    endif

    !write(*,*)"PARTITION"
    !write(*,*)gpat%sgraph(1)%nodeInPart


    if(lt%verbose >= 1) call prg_timer_stop(part_timer)

  end subroutine gpmd_graphpart

