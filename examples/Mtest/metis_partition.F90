!> Standalone graph partitioning test program.
!!
program metisPartition_test

  use bml 
  use prg_progress_mod
  use prg_parallel_mod
  use prg_timer_mod
  use prg_graphsp2parser_mod
  use prg_sp2_mod
  use prg_graph_mod
  use prg_subgraphLoop_mod
  use prg_homolumo_mod
  use prg_graph_mod
  use prg_partition_mod

  implicit none     

  integer, parameter :: dp = kind(1.0d0)

  type(gsp2data_type)           :: gsp2
  type(bml_matrix_t)            :: h_bml
  type(graph_partitioning_t)    :: gp
  real(dp)                      :: sumCubes, maxCH, check_cost, smooth_maxCH, pnorm = 10
  integer                       :: N, M, i,j, neighbor, it, nparts=32
  integer, allocatable          :: xadj(:), adjncy(:), part(:)
  integer, allocatable          ::  CH_count(:), core_count(:)
  integer, allocatable          :: Halo_count(:,:), copy_Halo_count(:,:)
  integer                       :: niter=500, seed=1, node, rand_part, backup
  real                          :: u
  character(len=100) :: pname


  ! Start progress
  call progress_init()

  if (printRank() .eq. 1) then
      write(*,*) "Metis + SA Partition_test start ... "
  endif

  !> parsing input file.
   call parse_gsp2(gsp2,"input.in")
   
  !allocate
  call bml_zero_matrix(gsp2%bml_type, BML_ELEMENT_REAL, dp, gsp2%ndim, gsp2%mdim, h_bml)
  
  !read
  call bml_read_matrix(h_bml, gsp2%hamfile)
  
  N = bml_get_N(h_bml)
  M = bml_get_M(h_bml)
  
  !> allocate arrays 
  allocate(part(N))
  allocate(xadj(N+1))
  allocate(adjncy(N*M))
  allocate(CH_count(nparts))
  allocate(core_count(nparts))
  allocate(Halo_count(nparts,N))



  !! Threshold the graph
  call bml_threshold(h_bml, gsp2%gthreshold)

  !> Get graph data structures
  call bml_adjacency(h_bml, xadj, adjncy, 1);

  !partition with METIS
#ifdef DO_GRAPHLIB
  !call timer_start(graphsp2_timer, "TIme for METIS")
  call metisPartition(gp, N, N, xadj, adjncy, nparts, part, core_count, CH_count, Halo_count, sumCubes, maxCH,smooth_maxCH, pnorm)
  !call timer_stop(graphsp2_timer,1) 
#endif    
    
  !> compute cost of METIS partition
  
  call costPartition(gp, xadj, adjncy, part, core_count, CH_count, Halo_count, sumCubes, maxCH,smooth_maxCH, pnorm)
  write(*,*) "METIS_obj_value", maxCH
  
  !Improve partition with SA
  !call timer_start(dyn_timer, "TIme for Refinement")
  !call Kernlin_queue(gp, xadj, adjncy, part, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)
  !call KernLin(gp, xadj, adjncy, part, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm,10, seed)
  !call KernLin2(gp, xadj, adjncy, part, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)
  call simAnnealing_old(gp, xadj, adjncy, part, core_count, CH_count, Halo_count, sumCubes, maxCH,smooth_maxCH, pnorm, niter, seed)
  !call Kernlin_queue(gp, xadj, adjncy, part, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)
  !call simAnnealing(gp, xadj, adjncy, part, core_count, CH_count, Halo_count, sumCubes, maxCH,smooth_maxCH, pnorm, niter, seed)
  !call Kernlin_queue(gp, xadj, adjncy, part, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)
  !call costPartition(gp, xadj, adjncy, part, core_count, CH_count, Halo_count, sumCubes, maxCH,smooth_maxCH, pnorm)
  !call timer_stop(dyn_timer,1) 
  write(*,*) "METIS+Refine_obj_value", maxCH


  
  call bml_deallocate(h_bml)
  call destroyGraphPartitioning(gp)
  deallocate(xadj)
  deallocate(adjncy)
  deallocate(part)
  deallocate(CH_count)
  deallocate(Halo_count)
  
  call progress_shutdown()

  call exit(0)


end program metisPartition_test
