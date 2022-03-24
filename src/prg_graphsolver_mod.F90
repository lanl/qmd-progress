!> Module for graph-based solvers.
!!
!! \ingroup PROGRESS
!!
module prg_graphsolver_mod

  use bml
  use prg_parallel_mod
  use prg_progress_mod
  use prg_partition_mod
  use prg_graph_mod
  use prg_system_mod
  use prg_timer_mod
  use prg_extras_mod
  use prg_densitymatrix_mod
  use prg_subgraphloop_mod

  implicit none

  private  !Everything is private by default

  integer, parameter :: dp = kind(1.0d0)

  public :: prg_build_densityGP_T0

contains

  !> Builds the density matrix from \f$ H_0 \f$ using a graph-based approach.
  !! \param ham_bml Input Orthogonalized Hamiltonian matrix.
  !! \param g_bml Matrix to extract the graph from.
  !! \param rho_bml Density matrix.
  !! \param threshold Threshold for sparse matrix algebra.
  !! \param bndfil Filing factor.
  !!
  subroutine prg_build_densityGP_T0(ham_bml, g_bml, rho_bml, threshold, bndfil, &
       & Ef, nparts, verbose)

    integer, parameter ::  dp = kind(1.0d0)
    type(bml_matrix_t), intent(in) :: ham_bml, g_bml
    type(bml_matrix_t) :: aux_bml
    type(bml_matrix_t), intent(out) :: rho_bml
    real(dp) :: bndfil, threshold
    integer ::  myRank, norbs, nnodes, tnnz, i, vsize(2), inorbs
    integer, intent(inout) :: nparts
    real(dp), intent(inout) :: Ef
    character(20) :: bml_type
    integer, allocatable :: xadj(:), adjncy(:), vector(:)
    integer, allocatable  ::  part(:), core_count(:), Halo_count(:,:), CH_count(:)
    type(graph_partitioning_t) ::  gpat
    real(dp) :: maxCH, pnorm=6, smooth_maxCH, sumCubes, mlsi, mlsii
    type(system_type), allocatable    ::  syprt(:)
    integer, optional, intent(in) :: verbose

    ! Initialize progress MPI and get ranks
    if(getNRanks() == 0) call prg_progress_init()
    !call prg_progress_init()
    myRank = getMyRank() + 1
    write(6,*) 'myRank=', myRank
    bml_type = bml_get_type(ham_bml)

    ! Allocate bml matrices
    norbs = bml_get_N(ham_bml)
    if(.not. bml_allocated(rho_bml))then
      call bml_zero_matrix(bml_type,bml_element_real,dp,norbs,norbs,rho_bml)
    endif
    call bml_zero_matrix(bml_type,bml_element_real,dp,norbs,norbs,aux_bml)

    call bml_print_matrix("ham",ham_bml,0,10,0,10)

    ! Doing graph partitioning using Metis
    allocate(xadj(norbs+1)) ! Adjacency in csr format

    tnnz = 0
    do i=1,norbs
      tnnz = tnnz + bml_get_row_bandwidth(g_bml,i)
    enddo
    write(6,*) 'tnnz=', tnnz

    allocate(adjncy(tnnz+1))

    call bml_adjacency(g_bml, xadj, adjncy, 1)

    nnodes = norbs

    allocate(part(nnodes))
    allocate(core_count(nparts))
    allocate(CH_count(nparts))
    allocate(Halo_count(nparts, nnodes))

    call prg_metisPartition(gpat, nnodes, nnodes, xadj, adjncy, nparts, part, core_count,&
         CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)

    if(present(verbose) .and. (verbose > 1 .and. myRank == 1)) then
      write(*,*)"gpat = ",gpat%pname
      write(*,*)"totalParts = ",gpat%totalParts
      write(*,*)"totalNodes = ",gpat%totalNodes
      write(*,*)"partMin = ",gpat%localPartMin(1)
      write(*,*)"partMax = ",gpat%localPartMax(1)
    endif

    allocate(syprt(gpat%TotalParts))

    call prg_wait

    if(present(verbose) .and. (verbose > 1 .and. myRank == 1)) mlsi = mls()
    mlsi = mls()

    ! Extract core halo indices from partitions
    do i=1,gpat%TotalParts
      call bml_matrix2submatrix_index(g_bml,&
           gpat%sgraph(i)%nodeInPart,gpat%nnodesInPart(i),&
           gpat%sgraph(i)%core_halo_index, &
           vsize,.true.)
      gpat%sgraph(i)%lsize = vsize(1)
      gpat%sgraph(i)%llsize = vsize(2)
      write(6,*) 'nodeInpart size', size(gpat%sgraph(i)%nodeInPart)
      write(6,*) 'nodeInpart(i)', gpat%nnodesInPart(i)

      if(present(verbose) .and. (verbose > 1 .and. myRank == 1))&
           write(*,*)"Core and CH sizes:",vsize(2),vsize(1)
    enddo

    ! Construct DM for every part
    do i = gpat%localPartMin(myRank), gpat%localPartMax(myRank)
      !  do i=1,gpat%TotalParts ! If no MPI
      call bml_matrix2submatrix_index(g_bml,&
           gpat%sgraph(i)%nodeInPart,gpat%nnodesInPart(i),&
           gpat%sgraph(i)%core_halo_index, &
           vsize,.true.)
      gpat%sgraph(i)%lsize = vsize(1)
      gpat%sgraph(i)%llsize = vsize(2)
      if(myRank == 1)write(*,*)"Core and CH sizes:",gpat%sgraph(i)%llsize,gpat%sgraph(i)%lsize
      inorbs = vsize(1)
      call bml_zero_matrix("dense",bml_element_real,dp,inorbs,inorbs,syprt(i)%estr%ham)
      if(allocated(vector))deallocate(vector)
      allocate(vector(gpat%sgraph(i)%lsize))
      vector(:) = gpat%sgraph(i)%core_halo_index(1:gpat%sgraph(i)%lsize)
      call bml_matrix2submatrix(ham_bml, syprt(i)%estr%ham, vector, &
           & inorbs)

      call bml_print_matrix("ham_part_bml",syprt(i)%estr%ham,0,10,0,10)

      !Computing the density matrix with diagonalization
      mlsii = mls()
      call bml_zero_matrix("dense",bml_element_real,dp,inorbs,inorbs,syprt(i)%estr%rho)
      call prg_build_density_T_Fermi(syprt(i)%estr%ham,syprt(i)%estr%rho,bndfil, 0.1_dp, 0.0_dp)

      !call bml_print_matrix("rho_part_bml",syprt(i)%estr%rho,0,10,0,10)
      if(present(verbose) .and. verbose > 1 .and. myRank == 1)&
           write(*,*)"Time for prg_build_density_T0",mls()-mlsii
      call bml_submatrix2matrix(syprt(i)%estr%rho,rho_bml,vector,gpat%sgraph(i)%lsize,gpat%sgraph(i)%llsize,threshold)
    enddo

    ! Collect the small DM matrices into its full system corresponding
    call prg_partOrdering(gpat)
    call prg_collectMatrixFromParts(gpat, rho_bml)
    call prg_wait

    call bml_print_matrix("rho_bml",rho_bml,0,10,0,10)
    if(present(verbose) .and. (verbose > 1 .and. myRank == 1)) &
         write(*,*)"Total time for graph solver =",mls()-mlsi

    CALL BML_DEALLOCATE(AUX_BML)
    do i = gpat%localPartMin(myRank), gpat%localPartMax(myRank)
      call bml_deallocate(syprt(i)%estr%ham)
      call bml_deallocate(syprt(i)%estr%rho)
    enddo
    deallocate(syprt)

    !call prg_progress_shutdown
    if(getNRanks() == 0) call prg_progress_shutdown()
  end subroutine prg_build_densityGP_T0

end module prg_graphsolver_mod
