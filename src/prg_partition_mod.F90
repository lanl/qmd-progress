!> The partition module.
!! \brief Contains different partitioning algorihms such as Metis, Simulated Annealing etc.
!!  Also contains optimization routines to improve upon existing partitioning,
!! such as simulated annealing, etc.
!!
module prg_partition_mod

  use bml
  use prg_graph_mod
  use prg_parallel_mod
  use prg_extras_mod
  use, intrinsic :: iso_c_binding

  implicit none

  private  !Everything is private by default

  integer, parameter :: dp = kind(1.0d0)

  !> From /usr/include/metis.h
  !!
  !! IDXTYPEWIDTH = 32 --> metis_index_kind = 4
  !! IDXTYPEWIDTH = 64 --> metis_index_kind = 8
  integer, parameter :: metis_index_kind = METIS_INDEX_KIND

  !> From /usr/include/metis.h
  !!
  !! REALTYPEWIDTH = 32 --> metis_real_kind = kind(0e0)
  !! REALTYPEWIDTH = 64 --> metis_real_kind = kind(0d0)
  integer, parameter :: metis_real_kind = kind(METIS_REAL_KIND)

  public :: prg_metisPartition
  public :: prg_costPartition
  public :: update_prg_costPartition
  public :: prg_simAnnealing
  public :: prg_check_arrays
  public :: prg_KernLin
  public :: prg_KernLin2
  public :: prg_Kernlin_queue
  public :: prg_update_gp
  public :: prg_simAnnealing_old

#ifdef DO_GRAPHLIB
  interface

    integer function METIS_SetDefaultOptions(options) &
        bind(C, name="METIS_SetDefaultOptions")

      import metis_index_kind

      integer(kind=metis_index_kind), intent(in) :: options(*)

    end function METIS_SetDefaultOptions

    integer function METIS_PartGraphKway(nvtxs, ncon, xadj, adjncy, vwgt, &
        vsize, adjwgt, nparts, tpwgts, ubvec, options, objval, part) &
        bind(C, name="METIS_PartGraphKway")

      import metis_index_kind
      import metis_real_kind

      integer(kind=metis_index_kind), intent(in) :: nvtxs(*)
      integer(kind=metis_index_kind), intent(in) :: ncon(*)
      integer(kind=metis_index_kind), intent(in) :: xadj(*)
      integer(kind=metis_index_kind), intent(in) :: adjncy(*)
      integer(kind=metis_index_kind), intent(in) :: vwgt(*)
      integer(kind=metis_index_kind), intent(in) :: vsize(*)
      integer(kind=metis_index_kind), intent(in) :: adjwgt(*)
      integer(kind=metis_index_kind), intent(in) :: nparts(*)
      real(kind=metis_real_kind), intent(in) :: tpwgts(*)
      real(kind=metis_real_kind), intent(in) :: ubvec(*)
      integer(kind=metis_index_kind), intent(in) :: options(*)
      integer(kind=metis_index_kind), intent(inout) :: objval(*)
      integer(kind=metis_index_kind), intent(inout) :: part(*)

    end function METIS_PartGraphKway

  end interface
#endif

contains

#ifdef DO_GRAPHLIB
  subroutine METIS_SetDefaultOptions_wrapper(options)

    integer(kind=metis_index_kind), intent(in) :: options(:)
    integer :: result

    result = METIS_SetDefaultOptions(options)
    if (result /= 1) then
      write(*, *) "error calling METIS_SetDefaultOptions"
      stop
    end if

  end subroutine METIS_SetDefaultOptions_wrapper

  subroutine METIS_PartGraphKway_wrapper(nvtxs, ncon, xadj, adjncy, vwgt, &
      vsize, adjwgt, nparts, tpwgts, ubvec, options, objval, part)

    integer, intent(in) :: nvtxs
    integer, intent(in) :: ncon
    integer, intent(in) :: xadj(:)
    integer, intent(in) :: adjncy(:)
    integer, pointer, intent(in) :: vwgt(:)
    integer, pointer, intent(in) :: vsize(:)
    integer, pointer, intent(in) :: adjwgt(:)
    integer, intent(in) :: nparts
    double precision, pointer, intent(in) :: tpwgts(:)
    double precision, pointer, intent(in) :: ubvec(:)
    integer(kind=metis_index_kind), intent(in) :: options(:)
    integer, intent(inout) :: objval
    integer, intent(inout) :: part(:)

    integer(kind=metis_index_kind) :: nvtxs_metis(1)
    integer(kind=metis_index_kind) :: ncon_metis(1)
    integer(kind=metis_index_kind), allocatable :: xadj_metis(:)
    integer(kind=metis_index_kind), allocatable :: adjncy_metis(:)
    integer(kind=metis_index_kind), pointer :: vwgt_metis(:) => null()
    integer(kind=metis_index_kind), pointer :: vsize_metis(:) => null()
    integer(kind=metis_index_kind), pointer :: adjwgt_metis(:) => null()
    integer(kind=metis_index_kind) :: nparts_metis(1)
    real(kind=metis_real_kind), pointer :: tpwgts_metis(:) => null()
    real(kind=metis_real_kind), pointer :: ubvec_metis(:) => null()
    integer(kind=metis_index_kind) :: objval_metis(1)
    integer(kind=metis_index_kind), allocatable :: part_metis(:)

    integer :: result

    nvtxs_metis(1) = nvtxs
    ncon_metis(1) = ncon

    allocate(xadj_metis(size(xadj)))
    xadj_metis = xadj

    allocate(adjncy_metis(size(adjncy)))
    adjncy_metis = adjncy

    if (associated(vwgt)) then
      allocate(vwgt_metis(size(vwgt)))
      vwgt_metis = vwgt
    end if

    if (associated(vsize)) then
      allocate(vsize_metis(size(vsize)))
      vsize_metis = vsize
    end if

    if (associated(adjwgt)) then
      allocate(adjwgt_metis(size(adjwgt)))
      adjwgt_metis = adjwgt
    end if

    nparts_metis(1) = nparts

    if (associated(tpwgts)) then
      allocate(tpwgts_metis(size(tpwgts)))
      tpwgts_metis = tpwgts
    end if

    if (associated(ubvec)) then
      allocate(ubvec_metis(size(ubvec)))
      ubvec_metis = ubvec
    end if

    objval_metis(1) = objval
    part_metis = part

    result = METIS_PartGraphKway(nvtxs_metis, ncon_metis, xadj_metis, adjncy_metis, vwgt_metis, vsize_metis, adjwgt_metis, &
      nparts_metis, tpwgts_metis, ubvec_metis, options, objval_metis, part_metis)
    if (result /= 1) then
      write(*, *) "error calling METIS_PartGraphKway"
      stop
    end if

    if (associated(vwgt_metis)) then
      deallocate(vwgt_metis)
    end if

    if (associated(vsize_metis)) then
      deallocate(vsize_metis)
    end if

    if (associated(adjwgt_metis)) then
      deallocate(adjwgt_metis)
    end if

    if (associated(tpwgts_metis)) then
      deallocate(tpwgts_metis)
    end if

    if (associated(ubvec_metis)) then
      deallocate(ubvec_metis)
    end if

    objval = objval_metis(1)
    part = part_metis

  end subroutine METIS_PartGraphKway_wrapper
#endif

  !> Create graph partitions minizing number of cut edges
  !! \param gp Graph partitioning`
  !! \param ngroups Number of groups/nodes
  !! \param nnodes Number of nodes
  !! \param xadj CSR array of graph nodes
  !! \param adjncy CSR array of graph neighbors
  !! \param nparts Number of Parts
  !! \param part Partition vector
  !! \param core_count Array: number of core vertices in each part
  !! \param CH_count Array: number of core+halo vertices in each part
  !! \param Halo_count 2D Array of size nparts by totalNodes: Halo_count(i,j) = k, node j is a halo of part i\
  !! with k connections
  !! \param sumCubes Sum of cubes objective value
  !! \param maxCh maximum core-halo part size obective value
  subroutine prg_metisPartition(gp, ngroups, nnodes, xadj, adjncy, nparts, part, core_count, CH_count, Halo_count, sumCubes, &
      maxCH, smooth_maxCH, pnorm)

    implicit none

    type (graph_partitioning_t), intent(inout) :: gp

    integer(kind=metis_index_kind), allocatable :: options(:)  ! options for metis (64-bit wide)
    integer, allocatable, intent(inout)         :: xadj(:), adjncy(:), part(:)
    integer, intent(inout)                      :: nparts
    integer                                     :: ncon, objval
    integer                                     :: i, j
    integer, target                             :: dummy_vwgt, dummy_vsize, dummy_adjwgt
    real(8), target                             :: dummy_tpwgts, dummy_ubvec
    real(dp), intent (inout)                    :: sumCubes, maxCH, smooth_maxCH, pnorm
    integer, intent (in)                        :: ngroups, nnodes
    integer, allocatable, intent(inout)         :: CH_count(:), core_count(:)
    integer, allocatable, intent(inout)         :: Halo_count(:,:)
    integer, allocatable                        :: copy_core_count(:)
    integer, pointer                            :: vwgt(:) => null(), vsize(:) => null(), adjwgt(:) => null()
    ! type(c_ptr)                                 :: vwgt, vsize, adjwgt
    ! type(c_ptr)                                 :: tpwgts, ubvec
    real(8), pointer                            :: tpwgts(:) => null(), ubvec(:) => null()
    character(len=100)                          :: pname

    allocate(options(40))
    allocate(copy_core_count(nparts))

    write(pname, '("metisParts")')

    options=0
#ifdef DO_GRAPHLIB
    call METIS_SetDefaultOptions_wrapper(options)
#endif

    ncon        = 1
    objval      = 1

    options(1)  = 1 !METIS_PTYPE_KWAY
    options(2)  = 0 !METIS_OBJTYPE_CUT
    options(9)  = 1 !METIS_OPTION_SEED
    options(18) = 1 !Fortran-style numbering is assumed that starts from 1

    !> prg_initialize
    Halo_count  = 0
    CH_count    = 0
    core_count  = 0

    if (printRank() .eq. 1) then
      write(*,*) "prg_metisPartition_test start ..."
    endif
    ! prg_initialize gp
    call prg_initGraphPartitioning(gp, pname, nparts, ngroups, nnodes)

    !> Partition graph into nparts'
    write(*,*) "The number of nodes in the graph is:", gp%totalNodes, &
      gp%totalNodes2, ncon, nparts, objval

#ifdef DO_GRAPHLIB
    call METIS_PartGraphKway_wrapper(gp%totalNodes, ncon, xadj, adjncy, vwgt, &
      vsize, adjwgt, nparts, tpwgts, ubvec, options, objval, part)
#endif

    !> Compute cost of partition
    call prg_costPartition(gp, xadj, adjncy, part, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)
    write(*,*) "Cost of METIS", sumCubes, maxCH, pnorm

    !prg_initialize and fill up subgraph structure
    !! Assign node ids (mapped to orbitals as rows) to each node in each
    do i = 1, nparts
      gp%nnodesInPartAll(i) = core_count(i)
      copy_core_count(i) = core_count(i)
      call prg_initSubgraph(gp%sgraph(i), i, nnodes)
      allocate(gp%sgraph(i)%nodeInPart(core_count(i)))
      gp%nnodesInPart(i) = core_count(i)
    enddo

    !Assign node ids to sgraph
    do i = 1, gp%totalNodes
      copy_core_count(part(i)) = copy_core_count(part(i)) - 1
      ! --> core_count((part(i))) - copy_core_count(part(i)) <--- = 1,2,3, ..., core_count( partnumber ) with respect to i

      ! NOTE: nodes in gp%sgraph()%nodeInPart() are currently 0 based!
      gp%sgraph(part(i))%nodeInPart(core_count((part(i))) - copy_core_count(part(i)) ) = i -1
    end do
    do i = 1, nparts
      do j = 1, core_count(i)
        if( part( gp%sgraph(i)%nodeInPart(j)+1 ) /= i) then
          write(*,*) "ERROR: subgraph struc incorrect!!", "node=",gp%sgraph(i)%nodeInPart(j)+1 , &
            "part=",i, "actual_part=", part(gp%sgraph(i)%nodeInPart(j)+1 )
          stop
        end if

      end do
    end do

  end subroutine prg_metisPartition

  !> Compute cost of a partition
  !!!
  !! \param gp Graph partitioning
  !! \param xadj CSR array of graph nodes
  !! \param adjncy CSR array of graph neighbors
  !! \param nparts Number of Parts
  !! \param partNumber Partition vector
  !! \param core_count Array: number of core vertices in each part
  !! \param CH_count Array: number of core+halo vertices in each part
  !! \param Halo_count 2D Array of size nparts by totalNodes: Halo_count(i,j) = k, node j is a halo of part i\
  !! with k connections
  !! \param sumCubes Sum of cubes objective value
  !! \param maxCh maximum core-halo part size obective value
  subroutine prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)

    type (graph_partitioning_t), intent(inout)  :: gp
    integer, allocatable, intent(inout)         :: xadj(:), adjncy(:)
    integer, allocatable, intent(in)            :: partNumber(:)
    integer, allocatable, intent(inout)         :: core_count(:)
    integer                                     :: totalParts, totalNodes,  i, j, neighbor
    real(dp), intent (inout)                    :: sumCubes, maxCH, smooth_maxCH, pnorm
    integer, allocatable, intent(inout)         ::  CH_count(:)
    integer, allocatable, intent(inout)         :: Halo_count(:,:)
    real(dp)                                    :: temp

    maxCH        = 0
    smooth_maxCH = 0
    sumCubes     = 0
    totalParts   = gp%totalParts
    totalNodes   = gp%totalNodes

    !> prg_initialize
    Halo_count = 0
    CH_count   = 0
    core_count = 0

    do i = 1, totalNodes
      CH_count(partNumber(i)) = CH_count(partNumber(i)) + 1 !core count
      core_count(partNumber(i)) = core_count(partNumber(i)) + 1 !core count
      do j = xadj(i), xadj(i + 1) - 1
        neighbor = adjncy(j)
        if (partNumber(i) /= partNumber(neighbor)) then
          if (Halo_count(partNumber(i) ,neighbor) == 0) then
            CH_count(partNumber(i)) = CH_count(partNumber(i)) + 1 !halo count
            Halo_count(partNumber(i), neighbor) = 1
          else
            Halo_count(partNumber(i), neighbor) = Halo_count(partNumber(i), neighbor) + 1
          end if
        end if
      end do
    end do

    do i = 1, totalParts
      if (core_count(i) <= 1) then
        print *, "core count <= 1 for partition "//to_string(i)//"!"
        stop
      end if
      temp = real(CH_count(i), dp)
      sumCubes = sumCubes+  temp*temp*temp
      smooth_maxCH = smooth_maxCH + temp**int(pnorm)
      if (CH_count(i) > maxCH) then
        maxCH = CH_count(i)
      end if
    end do
    smooth_maxCH = smooth_maxCH**(1/pnorm)

  end subroutine prg_costPartition

  !> Update cost of partition and the different parameters
  !> node is moves into new_part
  !> For each neighbor of node, the following cases hold:
  !> Case 1: neighbor is in old_part
  !> Case 2: neighbor is in new_part
  !> Case 3: neighbor is neither in old_ or new_part
  !! \param gp Graph partitioning
  !! \param xadj CSR array of graph nodes
  !! \param adjncy CSR array of 1043365660.0000000graph neighbors
  !! \param nparts Number of Parts
  !! \param partNumber Partition vector
  !! \param core_count Array: number of core vertices in each part
  !! \param CH_count Array: number of core+halo vertices in each part
  !! \param Halo_count 2D Array of size nparts by totalNodes: Halo_count(i,j) = k, node j is a halo of part i\
  !! with k connections
  !! \param sumCubes Sum of cubes objective value
  !! \param maxCh maximum core-halo part size obective value
  !! \param node Vertex that has moved to new_part
  !! \param new_part new part that node has moved to
  subroutine update_prg_costPartition(gp, xadj, adjncy, partNumber,core_count, CH_count, Halo_count, sumCubes, maxCH,smooth_maxCH, pnorm, node, new_part)

    type (graph_partitioning_t), intent(inout)  :: gp
    integer, allocatable, intent(inout)         :: xadj(:), adjncy(:)
    integer, allocatable, intent(inout)         :: partNumber(:), core_count(:)
    integer                                     :: totalParts, totalNodes,  i, j, neighbor
    real(dp), intent (inout)                    :: sumCubes, maxCH, smooth_maxCH, pnorm
    integer,  intent (in)                       :: node, new_part
    integer, allocatable, intent(inout)         :: CH_count(:)
    integer, allocatable, intent(inout)         :: Halo_count(:,:)
    integer                                     :: old_part
    real(dp)                                    :: temp

    totalParts = gp%totalParts
    totalNodes = gp%totalNodes
    old_part   = partNumber(node)

    if (old_part /= new_part) then
      core_count(new_part) = core_count(new_part)  + 1
      core_count(old_part) = core_count(old_part)  - 1
      CH_count(old_part) = CH_count(old_part) - 1 !core--
      CH_count(new_part) = CH_count(new_part) + 1 !core--
      do i=xadj(node), xadj(node+1) -1
        neighbor = adjncy(i)
        if (node /= neighbor) then
          if(partNumber(neighbor) == old_part) then !case 1
            Halo_count(old_part, node) = Halo_count(old_part, node) + 1
            if (Halo_count(old_part, node) == 1) then
              CH_count(old_part) = CH_count(old_part) + 1 !halo++
            end if
            Halo_count(new_part, neighbor) = Halo_count(new_part, neighbor) + 1
            if (Halo_count(new_part, neighbor) == 1) then
              CH_count(new_part) = CH_count(new_part) + 1 !halo++
            end if
          else if (partNumber(neighbor) == new_part) then !case 2
            Halo_count(old_part, neighbor) = Halo_count(old_part, neighbor) - 1
            if (Halo_count(old_part, neighbor) == 0) then
              CH_count(old_part) = CH_count(old_part) - 1 !halo--
            else if (Halo_count(old_part, neighbor) < 0) then
              write(*,*) "ERROR: Halo_count value cannot be negative, case 2i"
              write(*,*) "input matrix should be perfectly symmetric"
              stop
            end if
            Halo_count(new_part, node) = Halo_count(new_part, node) - 1
            if (Halo_count(new_part, node) == 0) then
              CH_count(new_part) = CH_count(new_part) - 1 !halo--, halo has become a core
            else if (Halo_count(new_part, node) < 0) then
              write(*,*) "ERROR: Halo_count value cannot be negative, case 2ii"
              write(*,*) "input matrix should be perfectly symmetric"
              stop
            end if
          else !case 3
            Halo_count(old_part, neighbor) = Halo_count(old_part, neighbor) - 1
            if (Halo_count(old_part, neighbor) == 0) then
              CH_count(old_part) = CH_count(old_part) - 1 !halo--
            else if (Halo_count(old_part, neighbor) < 0) then
              write(*,*) "ERROR: Halo_count value cannot be negative, case 3"
              write(*,*) "input matrix should be perfectly symmetric"
              stop
            end if
            Halo_count(new_part, neighbor) = Halo_count(new_part, neighbor) + 1
            if (Halo_count(new_part, neighbor) == 1) then
              CH_count(new_part) = CH_count(new_part) + 1 !halo++
            end if
          end if
        end if
      end do
      partNumber(node) = new_part
      sumCubes     = 0
      maxCH        = 0
      smooth_maxCH = 0
      do i=1, totalParts
        temp = real(CH_count(i), dp)
        sumCubes = sumCubes+  temp*temp*temp
        smooth_maxCH = smooth_maxCH + temp**int(pnorm)
        if (CH_count(i) > maxCH) then
          maxCH = real(CH_count(i),dp)
        end if
      end do
      smooth_maxCH = smooth_maxCH**(1/pnorm)
    end if

  end subroutine update_prg_costPartition

  !> Compute acceptance probability for simulated annealing
  !! \param it iteration
  !! \param prg_delta (new_obj_value - old_obj_value)
  !! \param r acceptance probability
  subroutine prg_accept_prob(it, prg_delta, r)
    integer, intent(in)     :: it
    real(dp), intent(in)    :: prg_delta
    real, intent(inout)     :: r
    real                    :: temp
    temp = 1.0/it
    r = 1
    if (prg_delta > 0) then
      r = exp(-(prg_delta/10.0)/temp)
    end if

  end subroutine prg_accept_prob

  !>Choose objective function to work with
  !! \param cost output according to chosen obj_fun
  !! \param sumCubes Sum of cubes obj value
  !! \param maxCH maximum core-halo part size obective value
  !! \param obj_fun 0=sumcubes, 1=maxCH
  subroutine prg_costIndex(cost, sumCubes, maxCH, smooth_maxCH, obj_fun)
    real(dp), intent (inout) :: cost, maxCH,smooth_maxCH, sumCubes
    integer, intent (inout) ::  obj_fun

    cost = -1

    if (obj_fun .eq. 0) then
      cost = sumCubes
    else if (obj_fun .eq. 1) then
      cost = maxCH
    else if (obj_fun .eq. 2) then
      cost = smooth_maxCH
    end if

  end subroutine prg_costIndex

  !> Pick a random node
  !! \param gp graph partitioning structure
  !! \param node output node
  !! \param seed random seed
  subroutine prg_rand_node(gp,node, seed)
    type (graph_partitioning_t), intent(inout) :: gp
    integer, intent(inout) :: node, seed
    integer ::  totalNodes,  i
    real :: u
    !call srand(seed)
    call random_number(u)
    node = FLOOR(gp%totalNodes*u) + 1

  end subroutine prg_rand_node

  !> Graph partitioning based on Simulated Annealing
  !! \param gp Graph partitioning
  !! \param xadj CSR array of graph nodes
  !! \param adjncy CSR array of graph neighbors
  !! \param nparts Number of Parts
  !! \param partNumber Partition vector
  !! \param core_count Array: number of core vertices in each part
  !! \param CH_count Array: number of core+halo vertices in each part
  !! \param Halo_count 2D Array of size nparts by totalNodes: Halo_count(i,j) = k, node j is a halo of part i\
  !! with k connections
  !! \param sumCubes Sum of cubes objective value
  !! \param maxCh maximum core-halo part size obective value
  !! \param niter Number of iterations
  !! \param seed Random seed
  subroutine prg_simAnnealing(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH,smooth_maxCH,pnorm, niter, seed)

    type (graph_partitioning_t), intent(inout)  :: gp
    integer, allocatable, intent(inout)         :: xadj(:), adjncy(:)
    integer, allocatable, intent(inout)         :: partNumber(:), core_count(:)
    integer                                     :: totalParts, totalNodes,  it, i,j,k, neighbor,  node, part_backup
    integer                                     :: totalNodes2
    real(dp)                                    :: cost, prev_cost, prg_delta, prev_maxCH
    real(dp), intent (inout)                    :: sumCubes, maxCH, smooth_maxCH, pnorm
    integer,  intent (in)                       :: niter
    integer, intent(inout)                      :: seed
    integer, allocatable, intent(inout)         :: CH_count(:)
    integer, allocatable, intent(inout)         :: Halo_count(:,:)
    integer, allocatable                        :: copy_core_count(:), empty_parts(:)
    integer                                     :: obj_fun = 2, min_CH_part, no_empty_parts,  new_part
    real                                        :: r, u
    character(len=100)                          :: pname

    totalParts = gp%totalParts
    totalNodes = gp%totalNodes
    totalNodes2 = gp%totalNodes2
    allocate(copy_core_count(totalParts))
    allocate(empty_parts(totalParts))
    !call srand(seed)
    write(*,*) "SA called..."
    if (totalNodes .lt. totalParts) then
      write(*,*) "ERROR: Number of parts cannot be greater than number of nodes."
      stop
    end if

    !> Compute current cost of partition
    call prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)

    !> Choose objective function to minimize
    call prg_costIndex(cost, sumCubes, maxCH,smooth_maxCH, obj_fun)
    prev_cost = cost

    !> Perform SA
    do it=1, niter
      call prg_rand_node(gp, node, seed)
      !> Find part with smalles size (should be included in update_prg_costPartition
      min_CH_part = 1
      do k =1, totalParts
        if (CH_count(min_CH_part) .gt. CH_count(k)) then
          min_CH_part = k
        end if
      end do

      !> if part(node) == max_ch_part, try to move node and it's neighbors to min_ch_part
      !> else move neighbors to part(node)
      if (CH_count( partNumber(node) ) .eq. maxCH) then
        do j = xadj(node), xadj(node+1)-1
          neighbor = adjncy(j)
          part_backup = partNumber(neighbor)
          call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, &
            smooth_maxCH, pnorm, neighbor, min_CH_part)
          call prg_costIndex(cost, sumCubes, maxCH, smooth_maxCH, obj_fun)
          prg_delta = cost - prev_cost
          call prg_accept_prob(it, prg_delta, r)
          call random_number(u)
          if (u > r) then ! reject
            call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, neighbor, part_backup)
          else
            prev_cost = cost
          end if
        end do
      else
        if (CH_count( min_CH_part) .eq. 0) then
          do j = xadj(node), xadj(node+1)-1
            neighbor = adjncy(j)
            part_backup = partNumber(neighbor)
            call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, neighbor, min_CH_part)
            call prg_costIndex(cost, sumCubes, maxCH, smooth_maxCH, obj_fun)
            prg_delta = cost - prev_cost
            call prg_accept_prob(it, prg_delta, r)
            call random_number(u)
            if (u > r) then ! reject
              call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm,  neighbor, part_backup)
            else
              prev_cost = cost
            end if
          end do
        else
          do j = xadj(node), xadj(node+1)-1
            neighbor = adjncy(j)
            part_backup = partNumber(neighbor)
            call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, neighbor, partNumber(node))
            call prg_costIndex(cost, sumCubes, maxCH, smooth_maxCH, obj_fun)
            prg_delta = cost - prev_cost
            call prg_accept_prob(it, prg_delta, r)
            call random_number(u)
            if (u > r) then ! reject
              call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, neighbor, part_backup)
            else
              prev_cost = cost
            end if
          end do
        end if
      end if
    end do

    !> Check empty part exist
    !> move nodes from maxpart to empty part
    no_empty_parts = 0
    do i=1,totalParts
      if (CH_count(i) .eq. 0) then
        no_empty_parts = no_empty_parts + 1
        empty_parts(no_empty_parts) = i
      end if
    end do
    prev_maxCH = maxCH
    do node=1,totalNodes
      if (no_empty_parts .le. 0) then
        exit
      end if

      if (CH_count(partNumber(node)) .eq. maxCH) then
        new_part = empty_parts(no_empty_parts)
        call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, node,new_part )
        prev_maxCH = maxCH
        no_empty_parts = no_empty_parts - 1

        !> move it neighbor in the max parts to the newpart
        do j = xadj(node), xadj(node+1)-1
          neighbor = adjncy(j)
          if (CH_count(partNumber(neighbor)) .eq. maxCH) then
            part_backup = partNumber(neighbor)
            call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH,smooth_maxCH, pnorm, neighbor, new_part)
            call prg_costIndex(cost, sumCubes, maxCH,smooth_maxCH, obj_fun)
            if (maxCH .ge. prev_maxCH) then ! reject
              call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, neighbor, part_backup )
              call prg_costIndex(cost, sumCubes, maxCH, smooth_maxCH, obj_fun)
            else
              prev_cost = cost
              prev_maxCH = maxCH
            end if
          end if
        end do
      end if
    end do

    write(*,*) "Cost of meTIS+SA", sumCubes, maxCH, smooth_maxCH



    !> Update graph structure
    !prg_initialize and fill up subgraph structure
    !! Assign node ids (mapped to orbitals as rows) to each node in each
    call prg_destroyGraphPartitioning(gp)
    write(pname, '("SAParts")')
    call prg_initGraphPartitioning(gp, pname, totalParts, totalNodes, totalNodes2)
    do i = 1, totalParts
      gp%nnodesInPartAll(i) = core_count(i)
      copy_core_count(i) = core_count(i)
      call prg_initSubgraph(gp%sgraph(i), i, gp%totalNodes2)
      allocate(gp%sgraph(i)%nodeInPart(core_count(i)))
      gp%nnodesInPart(i) = core_count(i)
    enddo

    !Assign node ids to sgraph
    do i=1, gp%totalNodes
      copy_core_count(partNumber(i)) =copy_core_count(partNumber(i)) - 1
      ! core_count((part(i))) - copy_core_count(part(i)) = node postion in array
      gp%sgraph(partNumber(i))%nodeInPart(core_count((partNumber(i))) - copy_core_count(partNumber(i)) ) = i -1 !NOTE: nodes in gp%sgraph()%nodeInPart() are currently 0 based!
    end do


    !> For debuging
    do i = 1, totalParts
      do j = 1, core_count(i)
        if( partNumber( gp%sgraph(i)%nodeInPart(j)+1 ) /= i) then
          write(*,*) "ERROR: subgraph struc incorrect!!", "node=",gp%sgraph(i)%nodeInPart(j)+1 , "part=",i, "actual_part=", partNumber(gp%sgraph(i)%nodeInPart(j)+1 )
          stop
        end if
      end do
    end do
    do i=1, totalParts
      write(*,*) "part=",i, "C=", core_count(i), "CH=", CH_count(i)
      if (CH_count(i) .eq. 0) then
        write(*,*) "ERROR: SA produced an empty part"
        stop
      end if
    end do
  end subroutine prg_simAnnealing


  !> Graph partitioning based on inspired by Kernighan-Lin
  !> Review METiS manual for description of k-way implementation of KL
  !> Pick a core together with its halos
  !> Place free vertices on a priority queue with (key, value) =(prg_delta, best_part),with prg_delta = change in obj_value
  !> Dequeue and allow hill climbing
  !! \param gp Graph partitioning
  !! \param xadj CSR array of graph nodes
  !! \param adjncy CSR array of graph neighbors
  !! \param nparts Number of Parts
  !! \param partNumber Partition vector
  !! \param core_count Array: number of core vertices in each part
  !! \param CH_count Array: number of core+halo vertices in each part
  !! \param Halo_count 2D Array of size nparts by totalNodes: Halo_count(i,j) = k, node j is a halo of part i\
  !! with k connections
  !! \param sumCubes Sum of cubes objective value
  !! \param maxCh maximum core-halo part size obective value
  !! \param nconverg number of before convergence
  !! \param seed random number generator seed


  subroutine prg_KernLin(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, nconverg, seed)

    type (graph_partitioning_t), intent(inout)    :: gp
    integer, allocatable, intent(inout)           :: xadj(:), adjncy(:)
    integer, allocatable, intent(inout)           :: partNumber(:), core_count(:)
    integer                                       :: totalParts, totalNodes, i, iit, j,k, neighbor,  node, part_backup, h, node2
    real(dp), intent(inout)                       :: sumCubes, maxCH, smooth_maxCH, pnorm
    real(dp)                                      :: cost, prev_cost, prev_iteration_cost, prev_maxCH, minCH
    integer, intent(inout)                        :: seed
    integer, intent(in)                           :: nconverg
    integer, allocatable, intent(inout)           :: CH_count(:)
    integer, allocatable, intent(inout)           :: Halo_count(:,:)
    integer, allocatable                          :: copy_core_count(:)
    integer                                       :: obj_fun = 2, counter, min_part, max_Climb = 1, climb_counter, temp, convg_counter, converge, no_locked_nodes, empty_counter, backup, best_node, best_part,no_empty_parts, new_part
    integer                                       :: totalNodes2
    real                                          :: r, u
    integer, allocatable                          :: vertex_locked(:), hedge_span(:), node_backup(:), node_part_backup(:), nodes(:), empty_parts(:)
    character(len=100)                            :: pname

    totalNodes = gp%totalNodes
    totalNodes2 = gp%totalNodes2
    totalParts = gp%totalParts

    !> Allocate arrays
    allocate(vertex_locked(totalNodes))
    allocate(hedge_span(totalParts))
    allocate(node_backup(totalNodes))
    allocate(node_part_backup(totalNodes))
    allocate(nodes(totalNodes))
    allocate(copy_core_count(totalParts))
    allocate(empty_parts(totalParts))

    !> Initialize variables
    vertex_locked   = 0
    hedge_span      = 0
    counter         = 0
    climb_counter   = 1
    converge        = 0
    convg_counter   = 0
    iit             = 0
    core_count      = 0
    CH_count        = 0
    Halo_count      = 0
    no_locked_nodes = 0

    !> Initialize array of nodes
    do  i = 1, totalNodes
      nodes(i) = i
    end do

    !> Randomize nodes
    call prg_rand_shuffle(nodes, seed)


    !>Compute current cost of partition
    call prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)

    !>Choose objective function to minimize
    call prg_costIndex(cost, sumCubes, maxCH, smooth_maxCH, obj_fun)
    prev_cost           = cost
    prev_iteration_cost = cost
    prev_maxCH          = maxCH

    !> iterate over the columns of the matrix, ie the hyperedges


    do while ( converge .eq. 0 .and. iit .lt. nconverg)
      iit = iit + 1
      vertex_locked=0 ! Free vertices
      no_locked_nodes = 0
      call prg_rand_shuffle(nodes, seed)

      !> KL iteration
      do i=1, gp%totalNodes
        h = nodes(i) !h represents a hyperedge

        !> let min_part be the smallest CH_part
        minCH = totalNodes + 1
        do j = 1, totalParts
          if (CH_count(j) .lt. minCH) then
            min_part = j
            minCH = CH_count(j)
          end if
        end do
        if (min_part .eq. -1) then
          min_part = partNumber(h) !hyperedge h contains node h
          do j = xadj(h), xadj(h+1)-1
            node = adjncy(j)
            if (hedge_span( partNumber(node))==0) then
              counter = counter + 1
              hedge_span( partNumber(node)) = CH_count(partNumber(node))
              if (CH_count(partNumber(node)) .le. CH_count(min_part)) then
                min_part = partNumber(node)
              end if
            end if
            if (counter == totalParts) then
              counter = 0
              exit
            end if
          end do
        end if
        hedge_span    = 0

        if(h .eq. 0) then
          write(*,*) "error h =0"
          stop
        end if
        !> Try and move free nodes to min_part
        climb_counter = 1
        do j = xadj(h), xadj(h+1)-1
          node = adjncy(j)

          if (vertex_locked(node) .eq. 0  ) then
            part_backup = partNumber(node)
            node_part_backup(climb_counter) = part_backup
            node_backup(climb_counter) = node
            call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH,smooth_maxCH, pnorm, node, min_part)
            call prg_costIndex(cost, sumCubes, maxCH,smooth_maxCH, obj_fun)
            if (cost .le. prev_cost) then !accept
              prev_cost = cost
              no_locked_nodes = no_locked_nodes + climb_counter
              !write(*,*) maxCH

              !> lock vertices
              !> (climb_counter) vertices have been accepted
              !> need to lock (climb_counter) vertices
              !! Last vertex to be moved is node_backup(climb_counter)

              temp = climb_counter
              do k=1, temp
                vertex_locked(node_backup(climb_counter)) = 1
                climb_counter = climb_counter - 1
              end do

              !> reset
              climb_counter = 1
              node_backup = -1 !for debug purposes
              node_part_backup = -1 !for debug purposes

            else
              if (climb_counter .lt. max_climb) then !climb
                climb_counter = climb_counter + 1
              else !undo climb
                do k=1, max_climb
                  call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, node_backup(climb_counter), node_part_backup(climb_counter))
                  call prg_costIndex(cost, sumCubes, maxCH, smooth_maxCH, obj_fun)

                  climb_counter = climb_counter - 1
                end do
                if (prev_cost .ne. cost) then
                  write(*,*) "ERROR: There was an error in undo step 2", node, cost, prev_cost, j
                  stop
                end if
                climb_counter = 1 !reset
              end if
            end if
          end if

          !            end if
        end do
        !end of hyperedge, undo any climbs
        if(prev_cost .ne. cost) then
          temp = climb_counter

          !if climb_counter = 1, we cannot have prev_cost /= cost
          do k=1, temp-1
            climb_counter = climb_counter - 1
            call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH,smooth_maxCH, pnorm, node_backup(climb_counter), node_part_backup(climb_counter))
            call prg_costIndex(cost, sumCubes, maxCH, smooth_maxCH, obj_fun)

          end do
          if (prev_cost .ne. cost) then
            write(*,*) "ERROR: Undo after hyperedge"
            stop
          end if
        end if

        !> If all vertices locked, go to next iteration
        if (no_locked_nodes .eq. gp%totalNodes) then
          exit
        end if


        !> If empty parts exit, place a vertex in max_part there
        empty_counter = 0
        do j=1,totalParts
          if (core_count(j) .eq. 0) then
            empty_counter = empty_counter +1
            empty_parts(empty_counter) = j
          end if
        end do
        if (empty_counter .gt. 0) then
          do j= 1,totalNodes
            if (CH_count(partNumber(j) ) .eq. maxCH) then
              do k = xadj(j), xadj(j+1)-1
                !> Place j and it's neighbors that are in the max part into the empty part
                node2 = adjncy(k)
                backup = partNumber(j)
                if (partNumber(node2) .eq. partNumber(j) ) then
                  call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH,smooth_maxCH, pnorm, node2,empty_parts(empty_counter) )
                  call prg_costIndex(cost, sumCubes, maxCH,smooth_maxCH,  obj_fun)
                  prev_cost = cost
                  !prev_maxCH = maxCH
                  if (prev_maxCH .lt. maxCH) then !undo
                    call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH,smooth_maxCH, pnorm, node2,backup )
                    call prg_costIndex(cost, sumCubes, maxCH, smooth_maxCH, obj_fun)
                    prev_cost = cost
                    prev_maxCH = maxCH
                  end if
                end if
              end do
              empty_counter = empty_counter - 1
              if (empty_counter .eq. 0) then
                exit
              end if
            end if
          end do
        end if
      end do !End KL iterations

      !> Check Convergence
      if (prev_iteration_cost .eq. cost) then
        convg_counter = convg_counter + 1
        if (convg_counter .eq. nconverg) then
          converge =1
        end if
      else
        prev_iteration_cost = cost
        convg_counter = 0
      end if

    end do !while loop

    !> Check empty part exist
    !> move nodes from maxpart to empty part
    no_empty_parts = 0

    do i=1,totalParts
      if (CH_count(i) .eq. 0) then
        no_empty_parts = no_empty_parts + 1
        empty_parts(no_empty_parts) = i
      end if
    end do
    prev_maxCH = maxCH
    do node=1,totalNodes
      if (no_empty_parts .le. 0) then
        exit
      end if

      if (CH_count(partNumber(node)) .eq. maxCH .and. core_count(partNumber(node)) .ne. 1 ) then
        new_part = empty_parts(no_empty_parts)
        call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, node,new_part )
        prev_maxCH = maxCH
        no_empty_parts = no_empty_parts - 1

        !> move it neighbor in the max parts to the newpart
        do j = xadj(node), xadj(node+1)-1
          neighbor = adjncy(j)
          if (CH_count(partNumber(neighbor)) .eq. maxCH) then
            part_backup = partNumber(neighbor)
            call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH,smooth_maxCH, pnorm, neighbor, new_part)
            call prg_costIndex(cost, sumCubes, maxCH,smooth_maxCH, obj_fun)
            if (maxCH .ge. prev_maxCH) then ! reject
              call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, neighbor, part_backup )
              call prg_costIndex(cost, sumCubes, maxCH, smooth_maxCH, obj_fun)

            else
              prev_cost = cost
              prev_maxCH = maxCH
            end if
          end if
        end do
      end if
    end do






    write(*,*) "Cost of KL:", cost, maxCH, "No iterations:", iit

    !> Update graph structure
    !prg_initialize and fill up subgraph structure
    !! Assign node ids (mapped to orbitals as rows) to each node in each
    call prg_destroyGraphPartitioning(gp)
    write(pname, '("KLParts")')
    call prg_initGraphPartitioning(gp, pname, totalParts, totalNodes, totalNodes2)

    !> Allocate subgraph structure
    do i = 1, totalParts
      gp%nnodesInPartAll(i) = core_count(i)
      copy_core_count(i) = core_count(i)
      call prg_initSubgraph(gp%sgraph(i), i, gp%totalNodes2)
      allocate(gp%sgraph(i)%nodeInPart(core_count(i)))
      gp%nnodesInPart(i) = core_count(i)
    enddo

    !> Assign node ids to sgraph
    do i=1, gp%totalNodes
      copy_core_count(partNumber(i)) =copy_core_count(partNumber(i)) - 1
      !!>core_count((part(i))) - copy_core_count(part(i)) = node postion in array
      gp%sgraph(partNumber(i))%nodeInPart(core_count((partNumber(i))) - copy_core_count(partNumber(i)) ) = i -1 !NOTE: nodes in gp%sgraph()%nodeInPart() are currently 0 based!
    end do


    do i=1, totalParts
      write(*,*) "part=",i, "C=", core_count(i), "CH=", CH_count(i)
      if (CH_count(i) .eq. 0) then
        write(*,*) "ERROR: KL produced an empty part"
        stop
      end if
    end do
    deallocate(vertex_locked)
    deallocate(hedge_span)
    deallocate(node_backup)
    deallocate(node_part_backup)
    deallocate(nodes)
    deallocate(copy_core_count)

  end subroutine prg_KernLin



  subroutine prg_update_gp(gp, partNumber, core_count)
    type (graph_partitioning_t), intent(inout)    :: gp
    integer, allocatable, intent(inout)           :: partNumber(:), core_count(:)
    integer                                       :: totalParts, totalNodes, i
    integer                                       :: totalNodes2
    integer, allocatable                          :: copy_core_count(:)
    character(len=100)                            :: pname

    totalNodes = gp%totalNodes
    totalNodes2 = gp%totalNodes2
    totalParts = gp%totalParts
    allocate(copy_core_count(totalParts))
    !> Update graph structure
    !prg_initialize and fill up subgraph structure
    !! Assign node ids (mapped to orbitals as rows) to each node in each
    call prg_destroyGraphPartitioning(gp)
    write(pname, '("Parts")')
    call prg_initGraphPartitioning(gp, pname, totalParts, totalNodes, totalNodes2)

    !> Allocate subgraph structure
    do i = 1, totalParts
      gp%nnodesInPartAll(i) = core_count(i)
      copy_core_count(i) = core_count(i)
      call prg_initSubgraph(gp%sgraph(i), i, gp%totalNodes2)
      allocate(gp%sgraph(i)%nodeInPart(core_count(i)))
      gp%nnodesInPart(i) = core_count(i)
    enddo

    !> Assign node ids to sgraph
    do i=1, gp%totalNodes
      copy_core_count(partNumber(i)) =copy_core_count(partNumber(i)) - 1
      !!>core_count((part(i))) - copy_core_count(part(i)) = node postion in array
      gp%sgraph(partNumber(i))%nodeInPart(core_count((partNumber(i))) - copy_core_count(partNumber(i)) ) = i -1 !NOTE: nodes in gp%sgraph()%nodeInPart() are currently 0 based!
    end do


    deallocate(copy_core_count)
  end subroutine prg_update_gp


  !> Randomly shuffle array
  subroutine prg_rand_shuffle(array, seed)

    integer, intent(inout)  :: array(:), seed
    integer                 :: i, randpos, temp
    real                    :: r

    !> Random seed
    call srand(seed)

    !> Shuffle array
    do i = size(array), 2, -1
      call random_number(r)
      randpos = int(r * i) + 1
      temp = array(randpos)
      array(randpos) = array(i)
      array(i) = temp
    end do

  end subroutine prg_rand_shuffle


  !> Error checking
  !> Checking that core_count, CH_count, Halo_count match
  subroutine prg_check_arrays(gp, core_count, CH_count, Halo_count)
    type (graph_partitioning_t), intent(inout) :: gp
    integer, allocatable, intent(inout)  :: core_count(:)
    integer, allocatable, intent(inout) ::  CH_count(:)
    integer, allocatable, intent(inout)  :: Halo_count(:,:)
    integer                              :: i,j, check

    do i=1,gp%totalParts
      check = core_count(i)
      do j=1, gp%totalNodes
        if (Halo_count(i,j) >0) then
          check = check + 1
        end if
      end do
      if (check /= CH_count(i)) then

        write(*,*) "ERROR: Halo_count is incorrect!"
        write(*,*) "check=", check, "CH_count(i)", CH_count(i)
        stop
      end if
    end do

    write(*,*) "prg_check_arrays PASSED!"
  end subroutine prg_check_arrays

  !> Greedy algorithm. At each step it chooses the (vertex, new_part) pair with highest gain
  !> Currently implementation is very slow
  subroutine prg_Kernlin_queue(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)

    type (graph_partitioning_t), intent(inout)    :: gp
    integer, allocatable, intent(inout)           :: xadj(:), adjncy(:)
    integer, allocatable, intent(inout)           :: partNumber(:), core_count(:)
    integer                                       :: totalParts, totalNodes, i, iit, j,k, neighbor,  node, part_backup, h, node2
    real(dp), intent(inout)                       :: sumCubes, maxCH, smooth_maxCH, pnorm
    real(dp)                                      :: cost, prev_cost, prev_iteration_cost, prev_maxCH, minCH, best_obj_val, current_cost
    integer, allocatable, intent(inout)           :: CH_count(:)
    integer, allocatable, intent(inout)           :: Halo_count(:,:)
    integer, allocatable                          :: copy_core_count(:)
    integer                                       :: backup, best_node, best_part
    real                                          :: r, u
    integer, allocatable                          :: vertex_locked(:), hedge_span(:), node_backup(:), node_part_backup(:), nodes(:), emptyParts(:)
    character(len=100)                            :: pname



    do i =1, 5
      call prg_find_best_move(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, best_node, best_part )
      if(best_node .eq. 0) then
        write(*,*) "error: node is 0"
        stop
      end if
      call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, best_node, best_part)


    end do
    write(*,*) "kl finished", gp%totalParts
    do i=1,gp%totalParts
      write(*,*) "part=", i, core_count(i), "ch=", CH_count(i)
    end do

  end subroutine prg_Kernlin_queue

  !> For kerlin_queue to find (vertex, new_part) pair with highest gain
  subroutine prg_find_best_move(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, best_node, best_part )

    type (graph_partitioning_t), intent(inout)    :: gp
    integer, allocatable, intent(inout)           :: xadj(:), adjncy(:)
    integer, allocatable, intent(inout)           :: partNumber(:), core_count(:)
    integer                                       :: totalParts, totalNodes, i, iit, j,k, neighbor,  node, part_backup, h, node2
    integer                                       :: totalNodes2
    real(dp), intent(inout)                       :: sumCubes, maxCH, smooth_maxCH, pnorm
    real(dp)                                      :: cost, prev_cost, prev_iteration_cost, prev_maxCH, minCH, best_obj_val, current_cost
    integer, intent(inout)                        :: best_node, best_part
    integer, allocatable, intent(inout)           :: CH_count(:)
    integer, allocatable, intent(inout)           :: Halo_count(:,:)
    integer, allocatable                          :: copy_core_count(:)
    integer                                       :: obj_fun = 2,  backup
    real                                          :: r, u
    integer, allocatable                          :: vertex_locked(:), hedge_span(:), node_backup(:), node_part_backup(:), nodes(:), emptyParts(:)
    character(len=100)                            :: pname
    !!  iterate through hyperedges
    !!  find min part and empty part
    !!  move vertices into part incident to h with smalles CH_count or empty part if it exists

    totalNodes = gp%totalNodes
    totalNodes2 = gp%totalNodes2
    totalParts = gp%totalParts

    best_node = 0
    best_part = 0
    call prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)
    call prg_costIndex(cost, sumCubes, maxCH, smooth_maxCH, obj_fun)
    best_obj_val = cost

    do i=1, totalNodes
      do j=1,totalParts
        backup = partNumber(i)
        call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, i, j)
        call prg_costIndex(cost, sumCubes, maxCH, smooth_maxCH, obj_fun)

        if (cost .le. best_obj_val) then
          best_node = i
          best_part = j
          best_obj_val = cost
        end if
        call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, i, backup)
      end do
    end do

  end subroutine prg_find_best_move

  subroutine prg_KernLin2(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)

    type (graph_partitioning_t), intent(inout)    :: gp
    integer, allocatable, intent(inout)           :: xadj(:), adjncy(:)
    integer, allocatable, intent(inout)           :: partNumber(:), core_count(:)
    integer                                       :: totalParts, totalNodes, i, iit, j,k, neighbor,  node, part_backup
    integer                                       :: totalNodes2
    real(dp), intent(inout)                       :: sumCubes, maxCH, smooth_maxCH, pnorm
    real(dp)                                      :: cost, prev_cost, prev_maxCH, minCH
    integer, allocatable, intent(inout)           :: CH_count(:)
    real                                          :: r
    integer, allocatable, intent(inout)           :: Halo_count(:,:)
    integer, allocatable                          :: copy_core_count(:),  empty_parts(:)
    integer                                       :: obj_fun = 2, largest_Hedge = -1, search_part, min_CH_part, new_part, no_empty_parts, seed =1
    character(len=100)                            :: pname
    !!  iterate through hyperedges
    !!  find min part and empty part
    !!  move vertices into part incident to h with smalles CH_count or empty part if it exists

    totalNodes = gp%totalNodes
    totalNodes2 = gp%totalNodes2
    totalParts = gp%totalParts

    !> Allocate arrays
    allocate(copy_core_count(totalParts))
    allocate(empty_parts(totalParts))

    call prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)
    call prg_costIndex(cost, sumCubes, maxCH, smooth_maxCH, obj_fun)
    prev_cost = cost
    do iit = 1, 40
      do i = 1, totalParts
        !> Pick hyperedge with largest size or random hyperedge with probability 0.5
        !> We wiil change it to pick hyperedge with highest priority, where priority will be defined according
        !> to different metrics
        call random_number(r)
        if (r .ge. 0.7) then
          call prg_get_largest_Hedge_in_part(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, i, largest_Hedge)
        else
          call prg_rand_node(gp, largest_Hedge, seed)

        end if
        !> Find part with smalles size (should be included in update_prg_costPartition
        min_CH_part = 1
        do k =1, totalParts
          if (CH_count(min_CH_part) .gt. CH_count(k)) then
            min_CH_part = k
          end if
        end do

        !> if current part is max, move to min_part
        !> then move subsets (neighbors)
        if (CH_count(partNumber(largest_Hedge)) .eq. maxCH) then
          !>Move hyperedges to minCH part
          new_part = min_CH_part
        else
          new_part = partNumber(largest_Hedge)
        end if

        !> Try and move intersecting hyperedges
        do j = xadj(largest_Hedge), xadj(largest_Hedge + 1)-1
          neighbor = adjncy(j)
          part_backup = partNumber(neighbor)
          call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, neighbor, new_part)
          call prg_costIndex(cost, sumCubes, maxCH, smooth_maxCH, obj_fun)
          if (cost .gt. prev_cost) then ! reject
            call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, neighbor, part_backup)
            call prg_costIndex(cost, sumCubes, maxCH, smooth_maxCH, obj_fun)
          else
            prev_cost = cost
          end if
        end do
      end do
    end do

    !> Move k number of vertices. k should be small i.e k <=20, k set in prg_Kernlin_queue
    !> Only use this for small systems
    !call prg_Kernlin_queue(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)

    !> Check empty part exist
    !> move nodes from maxpart to empty part
    no_empty_parts = 0

    do i=1,totalParts
      if (CH_count(i) .eq. 0) then
        no_empty_parts = no_empty_parts + 1
        empty_parts(no_empty_parts) = i
      end if
    end do

    prev_maxCH = maxCH
    do node=1,totalNodes
      if (no_empty_parts .le. 0) then
        exit
      end if

      if (CH_count(partNumber(node)) .eq. maxCH .and. core_count(partNumber(node)) .ne. 1 ) then
        new_part = empty_parts(no_empty_parts)
        call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, node,new_part )
        prev_maxCH = maxCH
        no_empty_parts = no_empty_parts - 1

        !> move it neighbor in the max parts to the newpart
        do j = xadj(node), xadj(node+1)-1
          neighbor = adjncy(j)
          if (CH_count(partNumber(neighbor)) .eq. maxCH) then
            part_backup = partNumber(neighbor)
            call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH,smooth_maxCH, pnorm, neighbor, new_part)
            call prg_costIndex(cost, sumCubes, maxCH,smooth_maxCH, obj_fun)
            if (maxCH .ge. prev_maxCH) then ! reject
              call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, neighbor, part_backup )
              call prg_costIndex(cost, sumCubes, maxCH, smooth_maxCH, obj_fun)

            else
              prev_cost = cost
              prev_maxCH = maxCH
            end if
          end if
        end do
      else if (core_count(partNumber(node)) .ne. 1 ) then
        new_part = empty_parts(no_empty_parts)
        call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, node,new_part )
        prev_maxCH = maxCH
        no_empty_parts = no_empty_parts - 1
      end if
    end do


    !> Update graph structure
    !prg_initialize and fill up subgraph structure
    !! Assign node ids (mapped to orbitals as rows) to each node in each
    call prg_destroyGraphPartitioning(gp)
    write(pname, '("KLParts")')
    call prg_initGraphPartitioning(gp, pname, totalParts, totalNodes, totalNodes2)

    !> Allocate subgraph structure
    do i = 1, totalParts
      gp%nnodesInPartAll(i) = core_count(i)
      copy_core_count(i) = core_count(i)
      call prg_initSubgraph(gp%sgraph(i), i, gp%totalNodes2)
      allocate(gp%sgraph(i)%nodeInPart(core_count(i)))
      gp%nnodesInPart(i) = core_count(i)
    enddo

    !> Assign node ids to sgraph
    do i=1, gp%totalNodes
      copy_core_count(partNumber(i)) =copy_core_count(partNumber(i)) - 1
      !!>core_count((part(i))) - copy_core_count(part(i)) = node postion in array
      gp%sgraph(partNumber(i))%nodeInPart(core_count((partNumber(i))) - copy_core_count(partNumber(i)) ) = i -1 !NOTE: nodes in gp%sgraph()%nodeInPart() are currently 0 based!
    end do


    do i=1, totalParts
      write(*,*) "part=",i, "C=", core_count(i), "CH=", CH_count(i)
      if (CH_count(i) .eq. 0) then
        write(*,*) "ERROR: KL produced an empty part"
        stop
      end if
    end do

  end subroutine prg_KernLin2


  subroutine prg_get_largest_Hedge_in_part(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, search_part, largest_Hedge)
    type (graph_partitioning_t), intent(inout)    :: gp
    integer, allocatable, intent(inout)           :: xadj(:), adjncy(:)
    integer, allocatable, intent(inout)           :: partNumber(:), core_count(:)
    integer                                       :: totalParts, totalNodes, i, iit, j,k, neighbor,  node
    integer                                       :: totalNodes2
    real(dp), intent(inout)                       :: sumCubes, maxCH, smooth_maxCH, pnorm
    integer, allocatable, intent(inout)           :: CH_count(:)
    integer, allocatable, intent(inout)           :: Halo_count(:,:)
    integer, allocatable                          :: copy_core_count(:)
    integer                                       :: obj_fun = 2, largest_Hsize
    integer, intent(inout)                        :: search_part, largest_Hedge
    real                                          :: r, u


    totalNodes = gp%totalNodes
    totalNodes2 = gp%totalNodes2
    totalParts = gp%totalParts
    !> i can be viewed as a hyperedge
    !> for all hyperedges in search_part, pick the one with largest size
    largest_Hsize = 0
    do i=1, totalNodes
      if (partNumber(i) .eq. search_part) then
        if ( xadj(i + 1) - xadj(i) .gt. largest_Hsize) then
          largest_Hsize = xadj(i + 1) - xadj(i)
          largest_Hedge = i
        end if
      end if
    end do

  end subroutine


  subroutine prg_simAnnealing_old(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH,smooth_maxCH,pnorm, niter, seed)

    type (graph_partitioning_t), intent(inout)  :: gp
    integer, allocatable, intent(inout)         :: xadj(:), adjncy(:)
    integer, allocatable, intent(inout)         :: partNumber(:), core_count(:)
    integer                                     :: totalParts, totalNodes,  it, i,j,k, neighbor,  node, part_backup, obj_fun=0
    integer                                     :: totalNodes2
    real(dp)                                    :: cost, prev_cost, prg_delta, prev_maxCH
    real(dp), intent (inout)                    :: sumCubes, maxCH, smooth_maxCH, pnorm
    integer,  intent (in)                       :: niter
    integer, intent(inout)                      :: seed
    integer, allocatable, intent(inout)         :: CH_count(:)
    integer, allocatable, intent(inout)         :: Halo_count(:,:)
    integer, allocatable                        :: copy_core_count(:), empty_parts(:)
    real                                        :: r, u
    character(len=100)                          :: pname

    totalParts = gp%totalParts
    totalNodes = gp%totalNodes
    totalNodes2 = gp%totalNodes2

    allocate(copy_core_count(totalParts))
    allocate(empty_parts(totalParts))
    !call srand(seed)
    write(*,*) "SA called..."
    if (totalNodes .lt. totalParts) then
      write(*,*) "ERROR: Number of parts cannot be greater than number of nodes."
      stop
    end if

    !> Compute current cost of partition
    call prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)

    !> Choose objective function to minimize
    call prg_costIndex(cost, sumCubes, maxCH,smooth_maxCH, obj_fun)
    prev_cost = cost

    !> Perform SA
    do it=1, niter
      call prg_rand_node(gp, node, seed)
      do j = xadj(node), xadj(node+1)-1
        neighbor = adjncy(j)
        part_backup = partNumber(neighbor)
        call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, neighbor, partNumber(node))
        call prg_costIndex(cost, sumCubes, maxCH, smooth_maxCH, obj_fun)
        prg_delta = cost - prev_cost
        call prg_accept_prob(it, prg_delta, r)
        call random_number(u)
        if (u > r) then ! reject
          call update_prg_costPartition(gp, xadj, adjncy, partNumber, core_count, CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm, neighbor, part_backup)
        else
          prev_cost = cost
          !write(*,*) maxCH
        end if
      end do

    end do



    write(*,*) "Cost of meTIS+SA", sumCubes, maxCH, smooth_maxCH



    !> Update graph structure
    !prg_initialize and fill up subgraph structure
    !! Assign node ids (mapped to orbitals as rows) to each node in each
    call prg_destroyGraphPartitioning(gp)
    write(pname, '("SAParts")')
    call prg_initGraphPartitioning(gp, pname, totalParts, totalNodes, totalNodes2)
    do i = 1, totalParts
      gp%nnodesInPartAll(i) = core_count(i)
      copy_core_count(i) = core_count(i)
      call prg_initSubgraph(gp%sgraph(i), i, gp%totalNodes2)
      allocate(gp%sgraph(i)%nodeInPart(core_count(i)))
      gp%nnodesInPart(i) = core_count(i)
    enddo

    !Assign node ids to sgraph
    do i=1, gp%totalNodes
      copy_core_count(partNumber(i)) =copy_core_count(partNumber(i)) - 1
      ! core_count((part(i))) - copy_core_count(part(i)) = node postion in array
      gp%sgraph(partNumber(i))%nodeInPart(core_count((partNumber(i))) - copy_core_count(partNumber(i)) ) = i -1 !NOTE: nodes in gp%sgraph()%nodeInPart() are currently 0 based!
    end do


    !> For debuging
    do i = 1, totalParts
      do j = 1, core_count(i)
        if( partNumber( gp%sgraph(i)%nodeInPart(j)+1 ) /= i) then
          write(*,*) "ERROR: subgraph struc incorrect!!", "node=",gp%sgraph(i)%nodeInPart(j)+1 , "part=",i, "actual_part=", partNumber(gp%sgraph(i)%nodeInPart(j)+1 )
          stop
        end if
      end do
    end do
    do i=1, totalParts
      write(*,*) "part=",i, "C=", core_count(i), "CH=", CH_count(i)
      if (CH_count(i) .eq. 0) then
        write(*,*) "ERROR: SA produced an empty part"

      end if
    end do
  end subroutine prg_simAnnealing_old


end module prg_partition_mod
