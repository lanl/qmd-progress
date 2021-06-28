!> The subgraphloop module.
!
!

module prg_subgraphloop_mod

  use prg_graph_mod
  use prg_sp2_mod
  use prg_timer_mod
  use prg_parallel_mod
  use bml
  use omp_lib

  implicit none

  private  !Everything is private by default

  integer, parameter :: dp = kind(1.0d0)

  public :: prg_subgraphSP2Loop
  public :: prg_balanceParts
  public :: prg_partOrdering
  public :: prg_getPartitionHalosFromGraph
  public :: prg_getGroupPartitionHalosFromGraph
  public :: prg_collectMatrixFromParts

contains

  !! Subgraph SP2 loop.
  !!
  !! \param h_bml Input Hamiltonian matrix
  !! \param g_bml Input Graph as matrix
  !! \param rho_bml Output density matrix
  !! \param gp Input/Output graph partitioning
  !! \param threshold threshold for matrix operations
  subroutine prg_subgraphSP2Loop(h_bml, g_bml, rho_bml, gp, threshold)

    type (bml_matrix_t), intent(in) :: h_bml, g_bml
    type (bml_matrix_t), intent(inout) :: rho_bml
    type (graph_partitioning_t), intent(inout) :: gp
    real(dp), intent(in) :: threshold

    integer :: i, j, k
    integer, allocatable :: vsize(:), vector(:)
    type (bml_matrix_t) :: x_bml
    real(dp) :: thresh0

    allocate(vsize(2))

    thresh0 = 0.0_dp

    call prg_timer_start(subind_timer)
    ! Determine elements for each subgraph
    !$omp parallel do default(none) &
    !$omp private(i, vsize, vector) &
    !$omp shared(gp, h_bml, g_bml) 
    do i = 1, gp%totalParts
  
      if(allocated(vector))deallocate(vector)
   allocate(vector(gp%sgraph(i)%lsize))
   vector(:) = gp%sgraph(i)%core_halo_index(1:gp%sgraph(i)%lsize)
       call bml_matrix2submatrix_index(g_bml, &
            gp%sgraph(i)%nodeInPart, gp%nnodesInPart(i), &
            gp%sgraph(i)%core_halo_index, &
            vsize, .true., h_bml)

       gp%sgraph(i)%lsize = vsize(1)
       gp%sgraph(i)%llsize = vsize(2)
       write(*,*)"CH",vector
    enddo
    !$omp end parallel do
    write(*,*)"AAA2"

    !    do i = 1, gp%totalParts
    !        write(*,*)"i = ", i, " core size = ", gp%sgraph(i)%llsize, &
    !                  " full size = ", gp%sgraph(i)%lsize
    !        write(*,*)"nnodes = ", gp%nnodesInPart(i), &
    !                  (gp%sgraph(i)%nodeInPart(k),k=1,gp%nnodesInPart(i))
    !    enddo

    call prg_timer_stop(subind_timer)

    deallocate(vsize)

    write(*,*)"AAA3"
    ! Balance parts by size of subgraph
    if (getNRanks() > 1) then
      call prg_balanceParts(gp)
      call prg_partOrdering(gp)
    endif
<<<<<<< HEAD
=======
    write(*,*)"AAA4"
>>>>>>> 62443b10f4fdae4b4415886b5128c81e376021a8
    ! Process each part one at a time
    !do i = 1, gp%nparts

    write(*,*)"Local parts", gp%localPartMin(gp%myRank+1), gp%localPartMax(gp%myRank+1)
    do i = gp%localPartMin(gp%myRank+1), gp%localPartMax(gp%myRank+1)

      call bml_zero_matrix(BML_MATRIX_DENSE, BML_ELEMENT_REAL, dp, &
           gp%sgraph(i)%lsize, gp%sgraph(i)%lsize, x_bml);

<<<<<<< HEAD
=======
    write(*,*)"AAA4",i

>>>>>>> 62443b10f4fdae4b4415886b5128c81e376021a8
         if(allocated(vector))deallocate(vector)
   allocate(vector(gp%sgraph(i)%lsize))
   vector(:) = gp%sgraph(i)%core_halo_index(1:gp%sgraph(i)%lsize)
   write(*,*)"vector",vector
       ! Extract subgraph and prg_normalize
       call prg_timer_start(subext_timer)
       call bml_matrix2submatrix(h_bml, x_bml, &
            vector, gp%sgraph(i)%lsize)
       call prg_timer_stop(subext_timer)

    write(*,*)"AAA5",thresh0,gp%maxIter,gp%mineval,gp%maxeval,gp%sgraph(i)%llsize
    write(*,*)"PP",gp%pp
    write(*,*)"VV",gp%vv
    call bml_print_matrix("rho_sp2",h_bml,0,10,0,10)
       ! Run SP2 on subgraph/submatrix
       call prg_timer_start(subsp2_timer)
       !call prg_sp2_submatrix_inplace(x_bml, threshold, gp%pp, &
       call prg_sp2_submatrix_inplace(x_bml, thresh0, gp%pp, &
            gp%maxIter, gp%vv, gp%mineval, gp%maxeval, &
            gp%sgraph(i)%llsize)
       call prg_timer_stop(subsp2_timer) 

<<<<<<< HEAD
=======
    write(*,*)"AAA6"
>>>>>>> 62443b10f4fdae4b4415886b5128c81e376021a8
    call bml_print_matrix("rho_sp2",x_bml,0,10,0,10)
       ! Reassemble into density matrix
       call prg_timer_start(suball_timer)
       call bml_submatrix2matrix(x_bml, rho_bml, &
            gp%sgraph(i)%core_halo_index, gp%sgraph(i)%lsize, &
            gp%sgraph(i)%llsize, threshold)
       call prg_timer_stop(suball_timer)

      call bml_deallocate(x_bml)

    enddo
<<<<<<< HEAD
=======
    write(*,*)"AAA5"
>>>>>>> 62443b10f4fdae4b4415886b5128c81e376021a8
    ! Fnorm
    call prg_fnormGraph(gp)

    ! Collect density matrix over local and distributed parts
    call prg_collectMatrixFromParts(gp, rho_bml)

  end subroutine prg_subgraphSP2Loop

  !> Collect distributed parts into same matrix.
  !!
  !! \param gp Graph partitioning
  !! \param rho_bml Matrix to be collected into
  subroutine prg_collectMatrixFromParts(gp, rho_bml)

    implicit none

    type (graph_partitioning_t), intent(inout) :: gp
    type (bml_matrix_t), intent(inout) :: rho_bml

    if (getNRanks() > 1) then
       write(*,*)getNRanks(),gp%nnodesInPart

       ! Save original domain
       call bml_save_domain(rho_bml)
       ! Update density matrix domain based on partitions of orbitals
       ! for gather
       call bml_update_domain(rho_bml, gp%localPartMin, gp%localPartMax, &
            gp%nnodesInPart) 
      !write(*,*)"gp%reorder",gp%reorder
      !stop
       call bml_reorder(rho_bml, gp%reorder);
       ! Exchange/gather density matrix across ranks
       call bml_allGatherVParallel(rho_bml)

      ! Reorder density matrix to match original
      call bml_reorder(rho_bml, gp%order)

      ! Restore to original domain
      call bml_restore_domain(rho_bml)

    endif

  end subroutine prg_collectMatrixFromParts

  !! Balance subgraphs by size.
  subroutine prg_balanceParts(gp)

    type (graph_partitioning_t), intent(inout) :: gp

    type (subgraph_t), allocatable :: temp_sgraph(:)
    type (subgraph_t) :: temp
    integer :: nranks, myRank, nParts, avgparts
    integer :: i, j, rid, sid, onum, pnode, ip

    integer :: tempp
    integer, allocatable :: tsize(:), torder(:)

    nranks = getNRanks()
    myRank = getMyRank()

    nParts = gp%totalParts
    avgparts = nParts / nranks
    ! Sort parts by core+halo size
    allocate(temp_sgraph(nParts))
    allocate(tsize(nParts))
    allocate(torder(nParts))

    do i = 1, nParts
      temp_sgraph(i) = gp%sgraph(i)
      torder(i) = i
      tsize(i) = temp_sgraph(i)%llsize
    enddo

    !! Sort parts by core+halo size
    do i = 1, nParts-1
      do j = i, nParts
        if (tsize(i) .lt. tsize(j)) then
          tempp = tsize(i)
          tsize(i) = tsize(j)
          tsize(j) = tempp

          tempp = torder(i)
          torder(i) = torder(j)
          torder(j) = tempp
        endif
      enddo
    enddo

    !! Print ordered subgraph sizes after sorting
    !   if (printRank() .eq. 1) then
    !      write(*,*) "after sort:"
    !      do i = 1, nParts
    !        write(*,*)"order ", i, " part = ", temp_sgraph(i)%part, " size = ", &
    !                  temp_sgraph(i)%lsize
    !      enddo
    !    endif

    !> Renumber parts
    !! Handle unbalanced numbers of parts.
    rid = 1
    do i = 1, nranks
      sid = sid + 1
      do j = i, nParts, nranks
        gp%sgraph(rid) = temp_sgraph(torder(j))
        !write(*,*) "rank = ", i, " part = "), rid, " from = ", j
        rid = rid + 1
        sid = sid + 1
      enddo
      if (sid .le. gp%localPartExtent(i)) then
        gp%sgraph(rid) = temp_sgraph(torder(nranks*avgparts+i))
        !write(*,*) "rank = ", i, " part = ", rid, " from = ", nranks*avgparts+i
        rid = rid + 1
      endif
    enddo
    !! Print ordered subgraph sizes after renumbering
    !    if (printRank() .eq. 1) then
    !      write(*,*) "after renumber:"
    !      do i = 1, nParts
    !        write(*,*)"order ", i, " part = ", gp%sgraph(i)%part, &
    !                  " core = ", gp%sgraph(i)%llsize, &
    !                  " size = ", gp%sgraph(i)%lsize, &
    !                  " min node = ", minval(gp%sgraph(i)%nodeInPart) + 1, &
    !                  " max node = ", maxval(gp%sgraph(i)%nodeInPart) + 1
    !      enddo
    !    endif

    !! Reset number of nodes per part
    if(allocated(gp%nnodesInPart))deallocate(gp%nnodesInPart)
    allocate(gp%nnodesInPart(nParts))
    if(allocated(gp%nnodesInPartAll))deallocate(gp%nnodesInPartAll)
    allocate(gp%nnodesInPartAll(nParts))
    do i = 1, nParts
      gp%nnodesInPart(i) = gp%sgraph(i)%llsize
      gp%nnodesInPartAll(i) = gp%sgraph(i)%llsize
    enddo

    deallocate(temp_sgraph)
    deallocate(tsize)
    deallocate(torder)

  end subroutine prg_balanceParts

  !> Set row ordering bases on parts
  !!
  !! \param gp Graph partitioning
  subroutine prg_partOrdering(gp)

    type (graph_partitioning_t), intent(inout) :: gp

    integer :: i, j, onum, pnode, nParts

    nParts = gp%totalParts

    !! Set ordering and reordering based on partitions
    onum = 0
    do i = 1, nParts
      !!write(*,*) i, " Part ", ip, " nnodes = ", gp%sgraph(i)%llsize
      do j = 1, gp%sgraph(i)%llsize
        pnode = gp%sgraph(i)%nodeInPart(j)
        !!write(*,*) "  ", j, "  Node ", pnode, " to ", onum
        gp%reorder(pnode+1) = onum
        gp%order(onum+1) = pnode
        onum = onum + 1
      enddo
    enddo
  end subroutine prg_partOrdering

  !> Get core+halo indeces for all partitions only using the graph.
  !!
  !! \param gp Graph partitioning
  !! \param g_bml Graph
  !! \param hnode Group start indeces
  !! \param djflg Double jump flag (true/false)
  subroutine prg_getGroupPartitionHalosFromGraph(gp, g_bml, hnode, djflag)

    implicit none

    type (graph_partitioning_t), intent(inout) :: gp
    type (bml_matrix_t), intent(in) :: g_bml
    integer, intent(in) :: hnode(*)
    logical, intent(in) :: djflag

    integer :: i, j, inode
    integer, allocatable :: vsize(:)

    allocate(vsize(2))

    !> Determine halo elements for each subgraph
    !$omp parallel do default(none) &
    !$omp private(i, j, inode, vsize) &
    !$omp shared(gp, g_bml, hnode, djflag)
    do i = 1, gp%totalParts

      do j = 1, gp%nnodesInPart(i)
        inode = gp%sgraph(i)%nodeInPart(j)
        gp%sgraph(i)%nodeInPart(j) = hnode(inode+1) - 1
      enddo

      call bml_matrix2submatrix_index(g_bml, &
           gp%sgraph(i)%nodeInPart, gp%nnodesInPart(i), &
           gp%sgraph(i)%core_halo_index, &
           vsize, djflag)

      gp%sgraph(i)%lsize = vsize(1)
      gp%sgraph(i)%llsize = vsize(2)

    enddo
    !$omp end parallel do

    deallocate(vsize)

  end subroutine prg_getGroupPartitionHalosFromGraph

  !> Get core+halo indeces for all partitions only using the graph.
  !!
  !! \param gp Graph partitioning
  !! \param g_bml Graph
  !! \param djflg Double jump flag (true/false)
  subroutine prg_getPartitionHalosFromGraph(gp, g_bml, djflag)

    implicit none

    type (graph_partitioning_t), intent(inout) :: gp
    type (bml_matrix_t), intent(in) :: g_bml
    logical, intent(in) :: djflag

    integer :: i
    integer, allocatable :: vsize(:)

    allocate(vsize(2))

    !> Determine halo elements for each subgraph
    !$omp parallel do default(none) &
    !$omp private(i, vsize) &
    !$omp shared(gp, g_bml, djflag)
    do i = 1, gp%totalParts

      call bml_matrix2submatrix_index(g_bml, &
           gp%sgraph(i)%nodeInPart, gp%nnodesInPart(i), &
           gp%sgraph(i)%core_halo_index, &
           vsize, djflag)

      gp%sgraph(i)%lsize = vsize(1)
      gp%sgraph(i)%llsize = vsize(2)

    enddo
    !$omp end parallel do

    deallocate(vsize)

  end subroutine prg_getPartitionHalosFromGraph

end module prg_subgraphloop_mod
