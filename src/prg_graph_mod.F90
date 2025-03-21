!> The graph module.
!! \ingroup PROGRESS
!
!

module prg_graph_mod

  use bml
  use prg_parallel_mod
  use omp_lib

  implicit none

  private  !Everything is private by default

  integer, parameter :: dp = kind(1.0d0)

  public :: subgraph_t
  public :: graph_partitioning_t
  public :: prg_initSubgraph
  public :: prg_destroySubgraph
  public :: prg_initGraphPartitioning
  public :: prg_destroyGraphPartitioning
  public :: prg_printGraphPartitioning
  public :: prg_equalPartition
  public :: prg_equalGroupPartition
  public :: prg_filePartition
  public :: prg_fnormGraph
  public :: prg_sedacsPartition

  !> Subgraph type
  type subgraph_t

    !> Partition number
    integer :: part

    !> Size of original matrix (h x h)
    integer :: hsize

    !> Size of full subgraph (l x l)
    integer :: lsize

    !> Size of core subgraph
    integer :: llsize

    !> Indeces from original matrix for subgraph core+halo extraction
    integer, allocatable :: core_halo_index(:)

    !> Nodes in this partition
    integer, allocatable :: nodeInPart(:)

    !> Trace per iteration
    !    real(dp) :: vvx(100)

  end type subgraph_t

  !> Graph partitioning type
  type graph_partitioning_t

    !> Partition name
    character(len=100) :: pname

    !> Local processor
    integer :: myRank

    !> Number of processors
    integer :: totalProcs

    !> Total number of global partitions
    integer :: totalParts

    !> Total number of global groups, nodes (or matrix rows)
    integer :: totalNodes

    !> Total number of global nodes (or matrix rows)
    integer :: totalNodes2

    !> Minimum global part number
    integer :: globalPartMin

    !> Maximum global part number
    integer :: globalPartMax

    !> Total global parts
    integer :: globalPartExtent

    !> Minimum part per processor
    integer, allocatable :: localPartMin(:)

    !> Maximum part per processor
    integer, allocatable :: localPartMax(:)

    !> Number of parts per processor
    integer, allocatable :: localPartExtent(:)

    !> Original ordering if required
    integer, allocatable :: order(:)

    !> Reordering if required
    integer, allocatable :: reorder(:)

    !> Total number of local partitions
    integer :: nparts

    !> For partitioning in x,y,z direction
    !! This is only used in case of "domain" type
    !! of partitioning
    integer :: nx, ny, nz

    !> Number of nodes in each local partition
    integer, allocatable :: nnodesInPart(:)

    !> Number of nodes in each partition
    integer, allocatable :: nnodesInPartAll(:)

    !> Sequence for SP2
    integer :: pp(100)

    !> Number of SP2 iterations
    integer :: maxIter

    !> Homo value
    real(dp) :: ehomo

    !> Lumo value
    real(dp) :: elumo

    !> Min eval for prg_normalize
    real(dp) :: mineval

    !> Max eval for prg_normalize
    real(dp) :: maxeval

    !> Trace per iteration
    real(dp) :: vv(100)

    !> Subgraph details
    type (subgraph_t), allocatable :: sgraph(:)

  end type graph_partitioning_t

contains

  !> Initialize subgraph.
  !! \param sg Subgraph
  !! \param pnum Part number
  !! \param hsize Size of full matrix
  subroutine prg_initSubgraph(sg, pnum, hsize)

    type (subgraph_t), intent(inout) :: sg
    integer, intent(in) :: pnum, hsize

    sg%part = pnum
    sg%hsize = hsize
    sg%lsize = 0
    sg%llsize = 0

    allocate(sg%core_halo_index(hsize))

  end subroutine prg_initSubgraph

  !> Destroy subgraph.
  !! \param sg Subgraph
  subroutine prg_destroySubgraph(sg)

    type (subgraph_t), intent(inout) :: sg

    if (allocated(sg%core_halo_index) .eqv. .true.) &
         deallocate(sg%core_halo_index)
    if (allocated(sg%nodeInPart) .eqv. .true.) deallocate(sg%nodeInPart)

  end subroutine prg_destroySubgraph

  !> Initialize graph partitioning.
  !! \param gp Graph partitioning
  !! \param pname Partitioning name
  !! \param np Number of partitions
  !! \param nnodes Number of groups/nodes
  !! \param nnodes2 Number of nodes
  subroutine prg_initGraphPartitioning(gp, pname, np, nnodes, nnodes2)

    type (graph_partitioning_t), intent(inout) :: gp
    character(len=*), intent(in) :: pname
    integer, intent(in) :: np, nnodes, nnodes2

    integer :: nprocs, avgparts
    integer :: i, it, last, nleft

    nprocs = getNRanks()
    gp%myRank = getMyRank()

    !! Global
    gp%pname = pname
    gp%totalProcs = nprocs
    gp%totalParts = np
    gp%totalNodes = nnodes
    gp%totalNodes2 = nnodes2

    !! Global bounds
    gp%globalPartMin = 1
    gp%globalPartMax = np
    gp%globalPartExtent = gp%globalPartMax - gp%globalPartMin + 1

    allocate(gp%localPartMin(nprocs))
    allocate(gp%localPartMax(nprocs))
    allocate(gp%localPartExtent(nprocs))

    ! Distribute parts evenly among ranks
    avgparts = gp%totalParts / nprocs
    do i = 1, nprocs
      gp%localPartExtent(i) = avgparts
    enddo

    nleft = gp%totalParts - nprocs * avgparts
    if (nleft .gt. 0) then
      do i = 1, nleft
        gp%localPartExtent(i) = gp%localPartExtent(i) + 1
      enddo
    endif

    gp%localPartMin(1) = 1
    gp%localPartMax(1) = gp%localPartExtent(1)
    do i = 2, nprocs
      gp%localPartMin(i) = gp%localPartMax(i-1) + 1
      gp%localPartMax(i) = gp%localPartMin(i) + gp%localPartExtent(i) - 1
    enddo

    gp%nparts = gp%localPartExtent(gp%myRank+1)

    allocate(gp%nnodesInPart(np))
    allocate(gp%nnodesInPartAll(np))
    allocate(gp%sgraph(np))

    !! For reordering
    allocate(gp%order(nnodes))
    allocate(gp%reorder(nnodes))

    gp%maxIter = 0
    gp%mineval = 0.0_dp
    gp%maxeval = 0.0_dp

    if (printRank() .eq. 1) then
      write(*,*)
      write(*,*) "total procs = ", gp%totalProcs
      write(*,*) "total nodes = ", gp%totalNodes
      write(*,*) "total nodes2 = ", gp%totalNodes2
      write(*,*) "total parts = ", gp%totalParts
      write(*,*) "local parts = ", gp%nparts
      write(*,*)
      write(*,*) "globalPartMin = ", gp%globalPartMin, &
           " globalPartMax = ", gp%globalPartMax, &
           " globalPartExtent = ", gp%globalPartExtent

      do i = 1, nprocs
        write(*,*) "rank = ", i-1, &
             " localPartMin = ", gp%localPartMin(i), &
             " localPartMax = ", gp%localPartMax(i), &
             " localPartExtent = ", gp%localPartExtent(i)
      enddo

      write(*,*)
    endif

  end subroutine prg_initGraphPartitioning

  !> Destroy graph partitioning
  !! \param sg Subgraph
  subroutine prg_destroyGraphPartitioning(gp)

    type (graph_partitioning_t), intent(inout) :: gp

    integer :: i

    if (allocated(gp%localPartMin))deallocate(gp%localPartMin)
    if (allocated(gp%localPartMax))deallocate(gp%localPartMax)
    if (allocated(gp%localPartExtent))deallocate(gp%localPartExtent)

    if (allocated(gp%order))deallocate(gp%order)
    if (allocated(gp%reorder))deallocate(gp%reorder)

    if(allocated(gp%nnodesInPart)) deallocate(gp%nnodesInPart)
    if(allocated(gp%nnodesInPartAll))deallocate(gp%nnodesInPartAll)

    if (allocated(gp%sgraph)) then
      do i = 1, gp%totalParts
        call prg_destroySubgraph(gp%sgraph(i))
      enddo
      deallocate(gp%sgraph)
    endif

  end subroutine prg_destroyGraphPartitioning

  !> Print graph partitioning structure data
  !! \param gp Graph partitioning
  subroutine prg_printGraphPartitioning(gp)

    type (graph_partitioning_t), intent(in) :: gp

    integer :: i, j

    if (gp%myRank .ne. 0) return

    ! Global data
    write(*,*) ""
    write(*,*) ""
    write(*,*) "Graph partitioning:"
    write(*,*) ""

    write(*,*) "name = ", gp%pname
    write(*,*) "totalProcs = ", gp%totalProcs
    write(*,*) "totalParts = ", gp%totalParts
    write(*,*) "totalNodes = ", gp%totalNodes
    write(*,*) "totalNodes2 = ", gp%totalNodes2
    write(*,*) ""
    write(*,*) "globalPartMin = ", gp%globalPartMin, &
         " globalPartMax = ", gp%globalPartMax, &
         " globalPartExtent = ", gp%globalPartExtent
    write(*,*) ""

    ! Local data
    write(*,*) "local parts = ", gp%nparts
    do i = 1, gp%totalProcs
      write(*,*) "rank = ", i-1, " localPartMin = ", gp%localPartMin(i), &
           " localPartMax = ", gp%localPartMax(i), &
           " localPartExtent = ", gp%localPartExtent(i)
    enddo
    write(*,*) ""

    ! SP2 data
    write(*,*) "Number of iterations = ", gp%maxIter
    write(*,*) "SP2 sequence = ", (gp%pp(i),i=1,gp%maxIter)
    write(*,*) "mineval = ", gp%mineval , "maxeval = ", gp%maxeval
    write(*,*) ""

    ! For each subgraph
    do i = 1, gp%nparts
      write(*,*) "part = ", i, " hsize = ", gp%sgraph(i)%hsize, &
           " lsize = ", gp%sgraph(i)%lsize, &
           " llsize = ", gp%sgraph(i)%llsize
      write(*,*) ""
      write(*,*) "Number of core nodes in part = ", gp%nnodesInPart(i)
      write(*,*) "nodeInPart = ", &
           (gp%sgraph(i)%nodeInPart(j), j=1,gp%nnodesInPart(i))
      write(*,*) ""

      write(*,*) "core_halo_index = ", &
           (gp%sgraph(i)%core_halo_index(j), j=1,gp%sgraph(i)%lsize)
      write(*,*) ""

    enddo

  end subroutine prg_printGraphPartitioning


  !> Create equal graph partitions, based on number of rows/orbitals
  !! \param gp Graph partitioning`
  !! \param nodesPerPart Number of core nodes per partition
  !! \param nnodes Total nodes in Hamiltonian matrix
  subroutine prg_equalPartition(gp, nodesPerPart, nnodes)

    type (graph_partitioning_t), intent(inout) :: gp
    integer, intent(in) :: nodesPerPart, nnodes

    integer :: i, j, it, np, psize
    character(len=100) :: pname

    !! Init graph partitioning
    if (nnodes.le.nodesPerPart)then
      np = 1
    else
      np = ceiling(real(nnodes) / real(nodesPerPart))
    endif
    pname = '("equalParts")'
    call prg_destroyGraphPartitioning(gp)
    call prg_initGraphPartitioning(gp, pname, np, nnodes, nnodes)

    ! Assign node ids (mapped to orbitals as rows) to each node in each
    ! partition
    !$omp parallel do default(none) private(i) &
    !$omp private(it,j,psize) &
    !$omp shared(gp,nnodes,nodesPerPart)
    do i = 1, gp%totalParts
      call prg_initSubgraph(gp%sgraph(i), i, nnodes)
      if ((i * nodesPerPart) .le. nnodes) then
        psize = nodesPerPart
      else
        psize = nnodes - (nodesPerPart * (i-1))
      endif
      allocate(gp%sgraph(i)%nodeInPart(psize))
      do j = 1, psize
        it = (i-1) * nodesPerPart + j-1;
        gp%sgraph(i)%nodeInPart(j) = it
      enddo
      gp%nnodesInPart(i) = psize
      gp%nnodesInPartAll(i) = psize
    enddo
    !$omp end parallel do

  end subroutine prg_equalPartition


  !> Create a partitioning based on node flips (as implemented in SEDACS - with several changes)
  !! \param gp Graph partitioning type
  !! \param whichParts_guess_saved Array indicating the "color" a node belongs to
  !! \param g_bml Adjacency matrix (BML type)
  !! \param nparts Number of total parts
  !! \param nnodes Total nodes in the graph (typically consistent with number of atoms)
  !! \param verb Verbosity level
  !!
  subroutine prg_sedacsPartition(gp,coords,whichParts_guess_saved,g_bml,nparts,nnodes,verb)
    implicit none
    type (graph_partitioning_t), intent(inout) :: gp
    integer, intent(in) :: nnodes, nparts
    integer, allocatable, intent(inout) :: whichParts_guess_saved(:)
    integer, optional, intent(in) :: verb
    integer :: i, cnt, part, j, it, psize, cut, cutOld, ac, nac
    integer ::  mdim, nodesPerPart,toBeDistributed,sumSizes,rem
    integer, allocatable :: whichParts(:)
    integer, allocatable :: graph(:,:),degs(:)
    real(kind(1.0)), allocatable :: row(:)
    real(dp) :: bal,sumDegs
    real(dp), allocatable, intent (in) :: coords(:,:)
    character(len=100) :: pname
    type (bml_matrix_t), intent(in) :: g_bml
    logical :: writeOut

    if(present(verb))then
      if(verb >= 1)then
        write(*,*)"Doing SEDACS type of partitioning ..."
        if(verb >= 2) writeOut = .true.
      endif
    endif

    !If there is no guess we will do a "chunk/block" partitioning
    !This initial pratitioning has all the nodes evenly distributed
    if(.not. allocated(whichParts_guess_saved))then
      nodesPerPart = int(nnodes/nparts)
      rem = nnodes - nodesPerPart*nparts
      toBeDistributed = rem

      pname = '("equalParts")'
      call prg_destroyGraphPartitioning(gp)
      call prg_initGraphPartitioning(gp, pname, nparts, nnodes, nnodes)

      sumSizes = 0
      do i = 1, gp%totalParts
        call prg_initSubgraph(gp%sgraph(i), i, nnodes)
        if (sumSizes .le. nnodes - nodesPerPart) then
          psize = nodesPerPart + min(toBeDistributed,1)
          toBeDistributed = max(toBeDistributed - 1,0)
        else
          psize = nnodes - sumSizes
        endif
        write(*,*)i,psize,sumSizes
        allocate(gp%sgraph(i)%nodeInPart(psize))
        do j = 1, psize
          it = sumSizes + j-1;
          gp%sgraph(i)%nodeInPart(j) = it
        enddo
        sumSizes = sumSizes + psize
        gp%nnodesInPart(i) = psize
        gp%nnodesInPartAll(i) = psize
      enddo

      call prg_get_parts_indices(gp,nnodes,whichParts)
      allocate(whichParts_guess_saved(nnodes))
    else
      allocate(whichParts(nnodes))
      whichParts = whichParts_guess_saved
    endif

    !Get degrees and build 2D-sedacs-like graph
    !In this case the degreess are saved in a separate array
    mdim = bml_get_m(g_bml)
    allocate(graph(nnodes,mdim))
    allocate(degs(nnodes))
    allocate(row(nnodes))
    row = 0.0
    degs = 0
    do i = 1,nnodes
      call bml_get_row(g_bml,i,row)
      do j = 1,nnodes
        if(row(j) .ge. 0.5)then
          degs(i) = degs(i) + 1
          graph(i,degs(i)) = j
        endif
      enddo
    enddo
    deallocate(row)

    sumDegs = sum(degs)
    call prg_get_graphCut(whichParts,degs,graph,cut) !Get the initial cut
    if(writeOut) write(*,*)"Iter, Cut, RelCut, and Bal",0,cut, cut/sumDegs,bal
    cutOld = cut
    do i = 1,200
      call do_flips_precomp(whichParts,coords,graph,g_bml,degs,nnodes,nparts,bal)
      call prg_get_graphCut(whichParts,degs,graph,cut)
      call prg_get_balancing(whichParts,nparts,bal)
      if(writeOut) write(*,*)"Iter, Cut, RelCut, and Bal",i,cut, cut/sumDegs,bal
      if(cut == cutOld)then
        exit
      else
        cutOld = cut
      endif
    enddo

    deallocate(graph)
    deallocate(degs)

    pname = '("sedacsPartition")'
    call prg_destroyGraphPartitioning(gp)
    call prg_initGraphPartitioning(gp, pname, nparts, nnodes, nnodes)

    !Assign node ids to each node in each part
    !$omp parallel do default(none) private(i) &
    !$omp private(cnt,it,j,psize,part) &
    !$omp shared(gp,nnodes,whichParts)
    do i = 1, gp%totalParts
      call prg_initSubgraph(gp%sgraph(i), i, nnodes)
      psize = int(sum(whichParts, mask=(whichParts==i))/i)
      allocate(gp%sgraph(i)%nodeInPart(psize))
      cnt = 0
      do j = 1,nnodes
        part = whichParts(j)
        if(part == i)then
          cnt = cnt + 1
          gp%sgraph(i)%nodeInPart(cnt) = j - 1
        endif
      enddo
      gp%nnodesInPart(i) = psize
      gp%nnodesInPartAll(i) = psize
    enddo
    !$omp end parallel do

    whichParts_guess_saved = whichParts
    deallocate(whichParts)

  end subroutine prg_sedacsPartition


  !! Do node partition flips with precomputed cuts.
  !! \brief This function is a special case of do_flips where
  !! the cuts around a node are precomputed for all possible
  !! part index that same node could have. This will differ from the do_flips
  !! since everytime there is a flip, there is no actualization of the cuts. The
  !! price to pay is the need of more iterations until convergence.
  !! \param whichPart partition indexing vector.
  !! \param graph Graph to be partition. graph[i,0] = degree of node i.
  !! graph[i,j>0] = the node conected to node i.
  !! \param nnodes Number of nodes.
  !! \param nparts Number of parts.
  !! \return whichPartNew New partition indexing verctor.
  !
  subroutine do_flips_precomp(whichParts,coords,graph,adj,degs,nnodes,nparts,bal)
    implicit none
    integer, allocatable :: cutsI(:,:)
    integer, intent(in) :: nnodes,nparts
    integer, allocatable, intent(inout) :: graph(:,:),degs(:)
    integer, allocatable, intent(inout) :: whichParts(:)
    real(dp) :: bal
    real(kind(1.0)), allocatable :: row(:)
    type(bml_matrix_t) :: adj
    integer :: i,j,deg,ind,ii,origCut,newCut, sumDeltaCut,ac,nac,newCut1
    integer :: origCutI,newCutI,origCutJ,newCutJ
    integer :: partIndexI, partIndexII,partIndexJ,k,cut,cutOld,newPartIndex,inDefect,discon
    integer, allocatable :: momentI(:),momentsI(:,:),partSizes(:)
    real(dp), allocatable, intent (in) :: coords(:,:)
    logical :: disconnected, conditionToFlip


    allocate(cutsI(nnodes,nparts))
    cutsI = 0
    sumDeltaCut = 0
    ac = 0
    nac = 0
    inDefect = 0

    !$omp parallel do default(none) private(i) &
    !$omp private(ii,deg,partIndexII,ind) &
    !$omp shared(degs,graph,nnodes,whichParts,cutsI)
    do i = 1,nnodes
      deg = degs(i)
      !Get the max cut a node could have
      cutsI(i,:) = deg
      !Lets look at every neighbor
      do ii = 1,deg
        ind = graph(i,ii)
        partIndexII = whichParts(ind)
        !Everytime there is a neighbor in a certain part
        !it will decrese the cut of I if I would be on that
        !same part.
        cutsI(i,partIndexII) = cutsI(i,partIndexII) - 1
      enddo
    enddo
    !$omp end parallel do

    !Get sizes
    allocate(row(nnodes))
    allocate(partSizes(nparts))
    partSizes = 0

    do i = 1,nnodes
      partIndexI = whichParts(i)
      partSizes(partIndexI) = partSizes(partIndexI) + 1
    enddo

    do i = 1,nnodes
      call bml_get_row(adj,i,row)
      if(cutsI(i,whichParts(i)) > 0)then
        do j = i+1, nnodes
          if((whichParts(i) .ne. whichParts(j)))then
            !conditionToFlip = .false.
            !disconnected = .false.

            !Look at their neighbors and count the cuts
            partIndexI = whichParts(i)
            partIndexJ = whichParts(j)
            origCutI = cutsI(i,partIndexI)
            origCutJ = cutsI(j,partIndexJ)

            newCutI = cutsI(i,partIndexJ)

            if(newCutI < degs(i))then
              newCutJ = cutsI(j,partIndexI)
              if(newCutJ < degs(j))then

                !Now we know the cut when I is in partIndexI and J
                origCut = origCutI + origCutJ
                newCut = newCutI + newCutJ

                if((newCut <= origCut))then

                  ac = ac + 1
                  whichParts(i) = partIndexJ
                  whichParts(j) = partIndexI

                  !Actualizing possible cuts of neighbors of i and j
                  do ii = 1,degs(i)
                    ind = graph(i,ii)
                    !The neighs of i "if now they are in the "new color
                    !of i", their cut "at that color" will be decreased
                    !by one.
                    cutsI(ind,partIndexJ) = cutsI(ind,partIndexJ) - 1
                    cutsI(ind,partIndexI) = cutsI(ind,partIndexI) + 1
                  enddo
                  do ii = 1,degs(j)
                    ind = graph(j,ii)
                    !Same for the neighs of j
                    cutsI(ind,partIndexI) = cutsI(ind,partIndexI) - 1
                    cutsI(ind,partIndexJ) = cutsI(ind,partIndexJ) + 1
                  enddo
                endif
              endif
            else
              nac = nac + 1
            endif

          elseif(whichParts(i) .eq. whichParts(j))then
            if(row(j) < 0.5)then
              origCut = 0
              newCut = 0
              partIndexI = whichParts(i)
              partIndexJ = whichParts(j)

              origCut = cutsI(i,partIndexI)
              origCut = origCut + cutsI(j,partIndexJ)

              ! Change color of I
              inDefect = minloc(partSizes,dim=1)

              newPartIndex = inDefect

              newCut = cutsI(i,newPartIndex)
              newCut = newCut + cutsI(j,partIndexJ)

              newCut1 = cutsI(i,partIndexI)
              newCut1 = newCut1 + cutsI(j,newPartIndex)

              if(newCut < origCut)then
                whichParts(i) = newPartIndex
                whichParts(j) = partIndexJ
                if(newPartIndex .ne.partIndexJ)then
                  partSizes(newPartIndex) = partSizes(newPartIndex) + 1
                  partSizes(partIndexJ) = partSizes(partIndexJ) - 1
                endif
                !sumDeltaCut = newCut - origCut + sumDeltaCut
              elseif(newCut1 < origCut)then
                whichParts(i) = partIndexI
                whichParts(j) = newPartIndex
                if(newPartIndex .ne. partIndexI)then
                  partSizes(newPartIndex) = partSizes(newPartIndex) + 1
                  partSizes(partIndexJ) = partSizes(partIndexJ) - 1
                endif
                !sumDeltaCut = newCut - origCut + sumDeltaCut
              else
                whichParts(i) = partIndexI
                whichParts(j) = partIndexJ
              endif
            endif
          else
            cycle
          endif
        enddo
      endif
    enddo

    deallocate(partSizes)
    deallocate(cutsI)
    deallocate(row)

  end subroutine do_flips_precomp



  !> Get balancing
  !! \param whichParts Vector containing the "color" (part)
  !! of each node.
  !! \param np Total number of parts
  !! \param bal Computed balancing
  !!
  subroutine prg_get_balancing(whichParts,np,bal)
    implicit none
    integer, allocatable, intent(in) :: whichParts(:)
    integer, allocatable :: partsSizes(:)
    integer :: np, i, nnodes
    real(dp), intent(out) :: bal

    allocate(partsSizes(np))
    nnodes = size(whichParts)
    partsSizes = 0
    do i = 1,nnodes
      !write(*,*)"whichParts",i,nnodes,whichParts(i),np
      if(whichParts(i) > np)then
        write(*,*)"!!!ERROR: Part index",whichParts(i),&
             "is larger than the total number of parts,",np
        stop
      endif
      partsSizes(whichParts(i)) = partsSizes(whichParts(i)) + 1
    enddo
    bal = real(maxval(partsSizes))/real(minval(partsSizes))
    deallocate(partsSizes)
  end subroutine prg_get_balancing

  !> Get part indices
  !! \param gp Graph partitioning type
  !! \param nnods Number of nodes
  !! \param whichParts Vector containing the "color" (part)
  !! of each node.
  !!
  subroutine prg_get_parts_indices(gp,nnodes,whichParts)
    type (graph_partitioning_t), intent(in) :: gp
    integer, allocatable,intent(out) :: whichParts(:)
    integer, intent(in) :: nnodes
    integer :: i,j,k,kk,partIndex

    if(.not. allocated(whichParts)) allocate(whichParts(nnodes))

    partIndex = 0
    whichParts = 0

    !$omp parallel do default(none) private(i) &
    !$omp private(k,kk) &
    !$omp shared(gp,nnodes,whichParts)
    do i = 1,gp%totalParts !Loop over parts
      do j = 1,gp%nnodesInPart(i)  !Loop over nodes in part
        k = gp%sgraph(i)%nodeInPart(j)
        kk = k + 1 !Indexig starts from 0 for all sgraphs (beacuse of historical reasons)
        whichParts(kk) = i
      enddo
    enddo
    !$omp end parallel do

  end subroutine prg_get_parts_indices


  !> Get graph cut
  !! \brief Get the total cut (edges crossing parts)
  !! \param whichParts Vector containing the "color" (part)
  !! of each node.
  !! \param degs Degrees of every node
  !! \param graph 2D array containing the information of every
  !! connection. graph(i,j): index of the jth neigbor of i.
  !! \param cut Computed cut
  !!
  subroutine prg_get_graphCut(whichParts,degs,graph,cut)
    integer, allocatable, intent(in) :: whichParts(:)
    integer, intent(out) :: cut
    integer :: nnodes,i,j,partIndexI,partIndexJ,k,kk
    integer :: ind
    integer, allocatable, intent(in) :: graph(:,:),degs(:)

    cut = 0
    nnodes = size(whichParts)

    do i = 1,nnodes
      !write(*,*)"graph",graph(i,1:degs(i))
      partIndexI = whichParts(i)
      !Look at the neighbors and check if they are in different part
      do j = 1,degs(i)
        ind = graph(i,j)
        partIndexJ = whichParts(ind)
        if(int(partIndexI - partIndexJ) .ne. 0)then
          cut = cut + 1
        endif
      enddo
    enddo
    cut = 0.5_dp*cut !Because of the double counting

  end subroutine prg_get_graphCut


  !> Create equal group graph partitions, based on number of atoms/groups
  !! \param gp Graph partitioning
  !! \param hindex Node indeces that represent ranges of atoms/groups
  !! \param ngroup Number of group nodes
  !! \param nodesPerPart Number of core nodes per partition
  !! \param nnodes Total nodes in Hamiltonian matrix
  subroutine prg_equalGroupPartition(gp, hindex, ngroup, nodesPerPart, nnodes)

    type (graph_partitioning_t), intent(inout) :: gp
    integer, intent(in) :: nodesPerPart, nnodes, ngroup
    integer, intent(in) :: hindex(2,ngroup)

    integer :: i, j, k, ll, np, psize, ind, ptotal
    character(len=100) :: pname

    !! Init graph partitioning
    np = ceiling(real(ngroup) / real(nodesPerPart))
    write(pname, '("equalGroupParts")')
    call prg_destroyGraphPartitioning(gp)
    call prg_initGraphPartitioning(gp, pname, np, ngroup, nnodes)

    !! Assign node ids (mapped to orbitals as rows) to each node in each
    !! partition
    !$omp parallel do default(none) &
    !$omp private(i, j, k, ll, ind, psize, ptotal) &
    !$omp shared(gp, hindex, nnodes, ngroup, nodesPerPart)
    do i = 1, gp%totalParts
      call prg_initSubgraph(gp%sgraph(i), i, nnodes)

      !! Figure out number of groups in part
      if ((i * nodesPerPart) .le. ngroup) then
        psize = nodesPerPart
      else
        psize = ngroup - (nodesPerPart * (i-1))
      endif

      !! Figure out total nodes/rows in part
      ind = (i-1)*nodesPerPart
      ptotal = 0
      do j = 1, psize
        ptotal = ptotal + hindex(2, ind+j) - hindex(1, ind+j) + 1
      enddo
      gp%nnodesInPart(i) = ptotal
      gp%nnodesInPartAll(i) = ptotal

      !! Enumerate all nodes in part
      allocate(gp%sgraph(i)%nodeInPart(ptotal))
      ll = 1
      do j = 1, psize
        do k = hindex(1, ind+j), hindex(2, ind+j)
          gp%sgraph(i)%nodeInPart(ll) = k-1
          ll = ll + 1
        enddo
      enddo
    enddo
    !$omp end parallel do

    !    do i = 1, gp%totalParts
    !      write(*,*) "part ", i, ": ", gp%nnodesInPart(i), " nodes"
    !      write(*,*) "     ", &
    !        (gp%sgraph(i)%nodeInPart(ll),ll = 1, gp%nnodesInPart(i))
    !    enddo

  end subroutine prg_equalGroupPartition

  !> Read graph partitions from a file, based on number of rows/orbitals
  !! \param partFile File containing core nodes for each partition
  !! \param gp Graph partitioning
  subroutine prg_filePartition(gp, partFile)

    type (graph_partitioning_t), intent(inout) :: gp
    character(len=*), intent(in) :: partFile

    call prg_readPart(gp, partFile)

  end subroutine prg_filePartition

  !> Read parts (core) from part file.
  !! \param gp Graph partitioning
  !! \param partFile Partition file
  subroutine prg_readPart(gp, partFile)

    character(len=*), intent(in) :: partFile
    type (graph_partitioning_t), intent(inout) :: gp

    integer :: pfile
    integer :: totalNodes, totalParts
    integer :: i, j, ip, pnode
    character(len=100) :: pname

    pfile = 10
    open(unit=pfile, status='old', FILE=partFile)

    read(pfile, *) pname
    read(pfile, *) totalNodes, totalParts

    call prg_initGraphPartitioning(gp, pname, totalParts, totalNodes, totalNodes)

    !! Read in part sizes
    read(pfile, *) (gp%nnodesInPartAll(i), i=1,gp%totalParts)
    do i = 1, gp%totalParts
      gp%nnodesInPart(i) = gp%nnodesInPartAll(i)
    enddo

    !! Read in nodes for each part
    do i = 1, gp%totalParts
      read(pfile, *) ip
      call prg_initSubgraph(gp%sgraph(i), i, totalNodes)
      allocate(gp%sgraph(i)%nodeInPart(gp%nnodesInPart(i)))
      read(pfile, *) (gp%sgraph(i)%nodeInPart(j),j=1,gp%nnodesInPart(i))
      do j = 1, gp%nnodesInPart(i)
        gp%sgraph(i)%nodeInPart(j) = gp%sgraph(i)%nodeInPart(j) + 1
      enddo
    enddo

    close(pfile)

  end subroutine prg_readPart

  !> Accumulate trace norm across all subgraphs
  !! \param gp Graph partitioning
  subroutine prg_fnormGraph(gp)

    type(graph_partitioning_t), intent(inout) :: gp

    integer :: i, j
#ifdef DO_MPI
    real(dp), allocatable :: sLocal(:), sGlobal(:)
#endif

#ifdef DO_MPI
    ! Sum traces from all parts on all ranks
    if (getNRanks() .gt. 1) then
      allocate(sLocal(gp%maxIter))
      allocate(sGlobal(gp%maxIter))
      do i = 1, gp%maxIter
        sLocal(i) = gp%vv(i)
      enddo
      call sumRealParallel(sLocal, sGlobal, gp%maxIter);
      do i = 1, gp%maxIter
        gp%vv(i) = sGlobal(i)
      enddo
      deallocate(sLocal, sGlobal)
    endif
#endif

    !! Take sqrt for fnorm per iter
    do i = 1, gp%maxIter
      gp%vv(i) = sqrt(gp%vv(i))
    enddo

    if (printRank() .eq. 1) then
      write(*,*)
      write(*,*) "prg_fnormGraph:"
      do i = 1, gp%maxIter
        write(*,*) "iter = ", i, " fnorm = ", gp%vv(i)
      enddo
    endif

  end subroutine prg_fnormGraph

end module prg_graph_mod
