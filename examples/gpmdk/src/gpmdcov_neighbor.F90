module gpmdcov_neighbor_mod

  use bml
  use prg_progress_mod
  use prg_parallel_mod
  use gpmdcov_writeout_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: low = kind(100)

  !> NL type
  type, public :: neighlist_type  !< The neighbor list type

    !> Number of atoms of the system.
    integer :: nats

    !> Number of atoms within distance of Rcut from atom I including atoms in the skin
    integer, allocatable :: nrnnlist(:)

    !> Distance between atom I(in box) and J (including atoms in the skin)
    real(dp), allocatable :: nndist(:,:)

    !> x-coordinte of neighbor J to I within RCut (including atoms in the skin)
    real(dp), allocatable :: nnRx(:,:)
    !> y-coordinte of neighbor J to I within RCut (including atoms in the skin)
    real(dp), allocatable :: nnRy(:,:)
    !> z-coordinte of neighbor J to I within RCut (including atoms in the skin)
    real(dp), allocatable :: nnRz(:,:)

    !> x-integer translation of neighbor J to I within RCut (including atoms in the skin)
    integer, allocatable :: nnIx(:,:)
    !> y-integer translation of neighbor J to I within RCut (including atoms in the skin)
    integer, allocatable :: nnIy(:,:)
    !> z-integer translation of neighbor J to I within RCut (including atoms in the skin)
    integer, allocatable :: nnIz(:,:)

    !> The neighbor J of I corresponds to some translated atom number in the box that we need to keep track of.
    integer, allocatable :: nnType(:,:)

    !> The neigbors J to I within Rcut that are all within the box (not in the skin).
    integer, allocatable :: nnStruct(:,:)

    !> Number of neigbors to I within Rcut that are all within the box (not in the skin).
    integer, allocatable :: nrnnStruct(:)

    !> Minimum distance between neighbors. nnStructMindist(i,j) means the Minimum distance
    !! between nnStruct(i,j) and i considering translations.
    real(dp), allocatable :: nnStructMindist(:,:)

  end type neighlist_type

  public :: gpmdcov_build_nlist_full, gpmdcov_destroy_nlist, gpmdcov_build_nlist_sparse
  public :: gpmdcov_build_nlist_sparse_v2, gpmdcov_build_nlist_sparse_v3, gpmd_nearestneighborlist,gpmdcov_get_vol
  public ::  gpmdcov_get_nlist_box_indices

contains


  !> Destroy the neigbor list to recover memory.
  !! \param nl Neigbor list structure.
  !!
  subroutine gpmdcov_destroy_nlist(nl,verbose)
    implicit none
    type(neighlist_type), intent(inout) ::  nl
    integer, optional, intent(in) :: verbose

    if(present(verbose))then
      if(verbose >= 1) write(*,*)"At gpmdcov_destroy_nlist. Destroying neighbor list ..."
    endif

    if(allocated(nl%nrnnlist))deallocate(nl%nrnnlist)
    if(allocated(nl%nndist))deallocate(nl%nndist)
    if(allocated(nl%nnRx))deallocate(nl%nnRx)
    if(allocated(nl%nnRy))deallocate(nl%nnRy)
    if(allocated(nl%nnRz))deallocate(nl%nnRz)
    if(allocated(nl%nnIx))deallocate(nl%nnIx)
    if(allocated(nl%nnIy))deallocate(nl%nnIy)
    if(allocated(nl%nnIz))deallocate(nl%nnIz)
    if(allocated(nl%nnType))deallocate(nl%nnType)
    if(allocated(nl%nnStruct))deallocate(nl%nnStruct)
    if(allocated(nl%nrnnStruct))deallocate(nl%nrnnStruct)
    if(allocated(nl%nnStructMindist))deallocate(nl%nnStructMindist)

  end subroutine gpmdcov_destroy_nlist

  !>  Build the neighbor list with a "brute force method"
  !! \brief It will bild a neighbor list using an "all to all" approach
  !! \param coords System coordinates. coords(1,7): x-coordinate of atom 7.
  !! \param lattice_vectors. Lattice vectors of the system box. lattice_vectors(1,3): z-coordinate of vector 1.
  !! \param nl Neighbor list type.
  !! \param verbose Verbosity level.
  !! \param rank MPI rank
  subroutine gpmdcov_build_nlist_full(coords,lattice_vectors,rcut,nl,verbose,rank,numranks)
    implicit none
    integer                              ::  cnt, i, j, k, maxneigh
    integer                              ::  myNumranks, myrank, nats, natsPerRank
    integer, allocatable                 ::  vectNnIx(:), vectNnIy(:), vectNnIz(:), vectNnStruct(:)
    integer, allocatable                 ::  vectNnType(:), vectNrnnStruct(:), vectNrnnlist(:)
    integer, intent(in)                  ::  verbose
    integer, optional, intent(in)        ::  numranks, rank
    real(dp)                             ::  coordsNeigh(3), density, distance, translation(3)
    real(dp)                             ::  volBox, lvm(3), lvd(3,3)
    real(dp), allocatable                ::  fcoords(:,:), fdvarray(:,:), dvarray(:,:), darray(:)
    real(dp), allocatable, intent(in)    ::  coords(:,:), lattice_vectors(:,:)
    real(dp), intent(in)                 ::  rcut
    type(neighlist_type), intent(inout)  ::  nl
#ifdef DO_MPI
    integer, allocatable :: rankRange(:,:)
#endif

    if(present(rank).and.present(numranks))then
      myrank = rank
      myNumranks = numranks
    else
      myrank = 1
      myNumranks = 1
    endif

    call gpmdcov_msI("gpmdcov_build_neigborlist","Building neighbor list ...",verbose,myrank)

    nats = size(coords,dim=2) !Get the number of atoms

    natsPerRank = int(nats/myNumranks)

#ifdef DO_MPI
    allocate(rankRange(2,myNumranks))
    do i = 1,myNumranks
      rankRange(1,i) = (i-1)*natsPerRank + 1
      rankRange(2,i) = i*natsPerRank
    enddo
    rankRange(2,myNumranks) = rankRange(2,myNumranks) + (nats - myNumranks*natsPerRank)
#endif

    !We will have approximatly [pi * rcut^3 * atomic density] number of neighbors.
    !This will always be bounded by [(2*rcut)^3 * atomic density] which we will use as
    !our max number of neighbors.
    call gpmdcov_get_vol(lattice_vectors,volBox)
    density = nats/volBox
    maxneigh = int(floor(8.0_dp * density * rcut**3))
    allocate(vectNnType(maxneigh*nats))
    vectNnType = 0
    allocate(vectNnStruct(maxneigh*nats))
    vectNnStruct = 0
    allocate(vectNrnnStruct(nats))
    vectNrnnStruct = 0
    allocate(vectNrnnlist(nats))
    vectNrnnlist = 0
    
    allocate(fcoords(nats,3))
    allocate(fdvarray(nats,3))
    allocate(dvarray(nats,3))
    allocate(darray(nats))
    
    do i = 1,3
       lvm(i) = sqrt(lattice_vectors(i,1)*lattice_vectors(i,1) &
            + lattice_vectors(i,2)*lattice_vectors(i,2) &
            + lattice_vectors(i,3)*lattice_vectors(i,3))
       lvd(i,:) = lattice_vectors(i,:)/lvm(i)
       fcoords(:,i) = (lvd(i,1)*coords(1,:)+lvd(i,2)*coords(2,:)+lvd(i,3)*coords(3,:))/lvm(i)
    end do

    !$omp parallel do default(none) private(i) &
    !$omp private(cnt,j) &
    !$omp private(distance,translation,coordsNeigh) &
    !$omp private(dvarray, darray, fdvarray) &
    !$omp shared(fcoords,coords,rcut,vectNnType) &
    !$omp shared(verbose,vectNnIx,vectNnIy,vectNnIz) &
    !$omp shared(vectNnStruct,vectNrnnStruct,vectNrnnlist) &
    !$omp shared(lattice_vectors)&
    !$omp shared(maxneigh) &
#ifdef DO_MPI
    !$omp shared(nats,rankRange,myRank)
    do i = rankRange(1,myRank),rankRange(2,myRank) !For every atom in the rank range
#else
    !$omp shared(nats)
    do i = 1,nats !For every atom
#endif
      do k = 1,3
         fdvarray(:,k) = modulo((fcoords(:,k) - fcoords(i,k)) + 0.5,1.) - 0.5
      enddo

      dvarray = matmul(fdvarray,lattice_vectors)
      darray = norm2(dvarray,DIM=2)

      cnt = 0
      do j = 1,nats !For every of its possible neighbors
         if (darray(j) .lt. rcut .and. darray(j) .gt. 1d-12) then
            cnt = cnt + 1
            vectNntype((i-1)*maxneigh + cnt) = j ! j is a neighbor of i by some translation
            vectNnstruct((i-1)*maxneigh + cnt) = j ! j is a neighbor of i by some translation
         endif
      enddo
      vectNrnnStruct(i) = cnt
      vectNrnnlist(i) = cnt
    enddo
    !$omp end parallel do

    !We do a sum reduction on all the vectors
#ifdef DO_MPI
    call prg_sumIntReduceN(vectNntype,maxneigh*nats)
    call prg_sumIntReduceN(vectNnstruct,maxneigh*nats)
    call prg_sumIntReduceN(vectNrnnStruct,nats)
    call prg_sumIntReduceN(vectNrnnlist,nats)
#endif

    !We deallocate the auxiliary vectors and transform them to matrices.
    if(.not.allocated(nl%nnType))allocate(nl%nnType(maxneigh,nats))
    call vectorToMatrixInt(maxneigh,nats,vectNnType,nl%nnType)
    deallocate(vectNnType)

    if(.not.allocated(nl%nnStruct))allocate(nl%nnStruct(maxneigh,nats))
    call vectorToMatrixInt(maxneigh,nats,vectNnStruct,nl%nnStruct)
    deallocate(vectNnStruct)

    if(.not.allocated(nl%nrnnStruct))allocate(nl%nrnnStruct(nats))
    nl%nrnnStruct = vectNrnnStruct
    deallocate(vectNrnnStruct)

    if(.not.allocated(nl%nrnnlist))allocate(nl%nrnnlist(nats))
    nl%nrnnlist = vectNrnnlist
    deallocate(vectNrnnlist)

#ifdef DO_MPI
    deallocate(rankRange)
#endif

    deallocate(fcoords)
    deallocate(fdvarray)
    deallocate(darray)
    deallocate(dvarray)
    
  end subroutine gpmdcov_build_nlist_full


  !>  Build the neighbor list with a "brute force method"
  !! \brief It will bild a neighbor list using an "all to all" approach
  !! \param coords System coordinates. coords(1,7): x-coordinate of atom 7.
  !! \param lattice_vectors. Lattice vectors of the system box. lattice_vectors(1,3): z-coordinate of vector 1.
  !! \param nl Neighbor list type.
  !! \param verbose Verbosity level.
  !! \param rank MPI rank
  subroutine gpmdcov_build_nlist_sparse(coords,lattice_vectors,rcut,nl,verbose,rank,numranks)
    implicit none
    integer                              ::  NBox, cnt, i, ibox
    integer                              ::  ith, ix, iy, iz
    integer                              ::  j, jbox, jj, jxBox
    integer                              ::  jyBox, jzBox, maxInBox, maxNeigh
    integer                              ::  myNumranks, myrank, nats, natsPerRank
    integer                              ::  nx, ny, nz, tx
    integer                              ::  ty, tz
    integer, allocatable                 ::  boxOfI(:), inbox(:,:), ithFromXYZ(:,:,:)
    integer, allocatable                 ::  totPerBox(:), xBox(:), yBox(:), zBox(:)
    integer, allocatable                 ::  vectNnIx(:), vectNnIy(:), vectNnIz(:), vectNnStruct(:)
    integer, allocatable                 ::  vectNnType(:), vectNrnnStruct(:), vectNrnnlist(:)
    integer, intent(in)                  ::  verbose
    integer, optional, intent(in)        ::  numranks, rank
    real(dp)                             ::  coordsNeigh(3), density, distance, translation(3)
    real(dp)                             ::  volBox, minx, miny, minz, smallReal
    real(dp), allocatable, intent(in)    ::  coords(:,:), lattice_vectors(:,:)
    real(dp), intent(in)                 ::  rcut
    type(neighlist_type), intent(inout)  ::  nl
#ifdef DO_MPI
    integer, allocatable :: rankRange(:,:)
#endif

    if(present(rank).and.present(numranks))then
      myrank = rank
      myNumranks = numranks
    else
      myrank = 1
      myNumranks = 1
    endif

    call gpmdcov_msI("gpmdcov_build_neigborlist","Building neighbor list ...",verbose,myrank)

    nats = size(coords,dim=2) !Get the number of atoms
    natsPerRank = int(nats/myNumranks)

#ifdef DO_MPI
    allocate(rankRange(2,myNumranks))
    do i = 1,myNumranks
      rankRange(1,i) = (i-1)*natsPerRank + 1
      rankRange(2,i) = i*natsPerRank
    enddo
    rankRange(2,myNumranks) = rankRange(2,myNumranks) + (nats - myNumranks*natsPerRank)
#endif

    !We will have approximatly [(4/3)*pi * rcut^3 * atomic density] number of neighbors.
    !A very large atomic density could be 1 atom per (1.0 Ang)^3 = 1 atoms per Ang^3  
    call gpmdcov_get_vol(lattice_vectors,volBox)
    density = 1.0_dp 
    maxneigh = int(floor(3.14592_dp * 4.0_dp/3.0_dp * density * rcut**3))

    !We assume the box is orthogonal
    nx = 1 + floor(lattice_vectors(1,1)/rcut)
    ny = 1 + floor(lattice_vectors(2,2)/rcut)
    nz = 1 + floor(lattice_vectors(3,3)/rcut)

    NBox = nx*ny*nz
    maxInBox = int(density*rcut**3) !Upper boud for the max number of atoms per box

    allocate(inbox(NBox,maxInBox))
    inbox = 0
    allocate(totPerBox(Nbox))
    totPerBox = 0
    allocate(boxOfI(nats))
    boxOfI = 0
    allocate(xBox(Nbox))
    xBox = 0
    allocate(yBox(Nbox))
    yBox = 0
    allocate(zBox(Nbox))
    zBox = 0
    allocate(ithFromXYZ(nx,ny,nz))
    ithFromXYZ = 0

    allocate(vectNnType(maxneigh*nats))
    vectNnType = 0
    allocate(vectNnIx(maxneigh*nats))
    vectNnIx = 0
    allocate(vectNnIy(maxneigh*nats))
    vectNnIy = 0
    allocate(vectNnIz(maxneigh*nats))
    vectNnIz = 0
    allocate(vectNnStruct(maxneigh*nats))
    vectNnStruct = 0
    allocate(vectNrnnStruct(nats))
    vectNrnnStruct = 0
    allocate(vectNrnnlist(nats))
    vectNrnnlist = 0

    minx = 1.0d10
    miny = 1.0d10
    minz = 1.0d10
    do i = 1,nats 
       minx = min(minx,coords(1,i))
       miny = min(miny,coords(2,i))
       minz = min(minz,coords(3,i))
    enddo

    smallReal = 0.0_dp
    !Search for the box coordinate and index of every atom
    do i = 1,nats
      !Index every atom respect to the discretized position on the simulation box.
      !tranlation = coords(:,i) - origin !For the general case we need to make sure coords ar > 0 
      ix = 1+ int(floor((coords(1,i) - minx + smallReal)/rcut)) !small box x-index of atom i
      iy = 1+ int(floor((coords(2,i) - miny + smallReal)/rcut)) !small box y-index //
      iz = 1+ int(floor((coords(3,i) - minz + smallReal)/rcut)) !small box z-index //
      
      if(ix > nx .or. ix < 0)Stop "Error in box index"
      if(iy > ny .or. iy < 0)Stop "Error in box index"
      if(iz > nz .or. iz < 0)Stop "Error in box index"

      ith =  ix + (iy-1)*(nx) + (iz-1)*(nx)*(ny)  !Get small box index
      boxOfI(i) = ith

      !From index to box coordinates
      xBox(ith) = ix
      yBox(ith) = iy
      zBox(ith) = iz

      !From box coordinates to index
      ithFromXYZ(ix,iy,iz) = ith

      totPerBox(ith) = totPerBox(ith) + 1 !How many per box
      if(totPerBox(ith) > maxInBox) Stop "Exceeding the max in box allowed"
      inbox(ith,totPerBox(ith)) = i !Who is in both ith

    enddo

     !For each atom we will look around to see who are its neighbors
    !$omp parallel do default(none) private(i) &
    !$omp private(ibox,ix,iy,iz) &
    !$omp private(jxbox,jybox,jzbox,jbox) &
    !$omp private(distance,translation,coordsNeigh) &
    !$omp private(cnt,j,jj,tx,ty,tz) &
    !$omp shared(xbox,ybox,zbox,boxOfI) &
    !$omp shared(nx,ny,nz,inbox,totperbox) &
    !$omp shared(coords,rcut,vectNnType,ithfromxyz) &
    !$omp shared(vectNnIx,vectNnIy,vectNnIz) &
    !$omp shared(vectNnStruct,vectNrnnStruct,vectNrnnlist) &
    !$omp shared(lattice_vectors)&
    !$omp shared(maxneigh) &
#ifdef DO_MPI
    !$omp shared(rankRange,myRank,nats)
    do i = rankRange(1,myRank),rankRange(2,myRank) !For every atom in the rank range
#else
    !$omp shared(nats)
    do i = 1,nats !For every atom
#endif
      cnt = 0
      !Which box it beongs to
      ibox = boxOfI(i)
      !Look inside the box and the neighboring boxes
      do ix = -1,1
        do iy = -1,1
          do iz = -1,1
            !Get neigh box coordinate
            jxBox = xBox(ibox) + ix
            jyBox = yBox(ibox) + iy
            jzBox = zBox(ibox) + iz
            tx = 0.0_dp ; ty = 0.0_dp ; tz = 0.0_dp
            if(jxBox <= 0)then
              jxBox = nx
              tx = -1
            elseif(jxBox > nx)then
              jxBox = 1
              tx = 1
            endif
            if(jyBox <= 0)then
              jyBox = ny
              ty = -1
            elseif(jyBox > ny)then
              jyBox = 1
              ty = 1
            endif
            if(jzBox <= 0)then
              jzBox = nz
              tz = -1
            elseif(jzBox > nz)then
              jzBox = 1
              tz = 1
            endif

            !Get the neigh box index
            jbox = ithFromXYZ(jxBox,jyBox,jzBox)

            !Now loop over the atoms in the jbox
            do j = 1,totPerBox(jbox)
              jj = inbox(jbox,j) !Get atoms in box j
              translation = tx*lattice_vectors(1,:) + ty*lattice_vectors(2,:) + tz*lattice_vectors(3,:)
              coordsNeigh = coords(:,jj) + translation
              distance = norm2(coords(:,i) - coordsNeigh)
              if (distance .lt. rcut .and. distance .gt. 1d-12) then
                cnt = cnt + 1
                vectNntype((i-1)*maxneigh + cnt) = jj ! jj is a neighbor of i by some translation
                vectNnstruct((i-1)*maxneigh + cnt) = jj ! jj is a neighbor of i by some translation
                vectNnIx((i-1)*maxneigh + cnt) = tx
                vectNnIy((i-1)*maxneigh + cnt) = ty
                vectNnIz((i-1)*maxneigh + cnt) = tz

              endif
            enddo
          enddo
        enddo
      enddo
      vectNrnnStruct(i) = cnt
      vectNrnnlist(i) = cnt
    enddo
!$omp end parallel do

    deallocate(inbox)
    deallocate(totPerBox)
    deallocate(boxOfI)
    deallocate(xBox)
    deallocate(yBox)
    deallocate(zBox)
    deallocate(ithFromXYZ)


    !We do a sum reduction on all the vectors
#ifdef DO_MPI
    call prg_sumIntReduceN(vectNntype,maxneigh*nats)
    call prg_sumIntReduceN(vectNnstruct,maxneigh*nats)
    call prg_sumIntReduceN(vectNnIx,maxneigh*nats)
    call prg_sumIntReduceN(vectNnIy,maxneigh*nats)
    call prg_sumIntReduceN(vectNnIz,maxneigh*nats)
    call prg_sumIntReduceN(vectNrnnStruct,nats)
    call prg_sumIntReduceN(vectNrnnlist,nats)
#endif

    !We deallocate the auxiliary vectors and transform them to matrices.
    if(.not.allocated(nl%nnType))allocate(nl%nnType(maxneigh,nats))
    call vectorToMatrixInt(maxneigh,nats,vectNnType,nl%nnType)
    deallocate(vectNnType)

    if(.not.allocated(nl%nnIx))allocate(nl%nnIx(maxneigh,nats))
    call vectorToMatrixIntLow(maxneigh,nats,vectNnIx,nl%nnIx)
    deallocate(vectNnIx)

    if(.not.allocated(nl%nnIy))allocate(nl%nnIy(maxneigh,nats))
    call vectorToMatrixIntLow(maxneigh,nats,vectNnIy,nl%nnIy)
    deallocate(vectNnIy)

    if(.not.allocated(nl%nnIz))allocate(nl%nnIz(maxneigh,nats))
    call vectorToMatrixIntLow(maxneigh,nats,vectNnIz,nl%nnIz)
    deallocate(vectNnIz)

    if(.not.allocated(nl%nnStruct))allocate(nl%nnStruct(maxneigh,nats))
    call vectorToMatrixInt(maxneigh,nats,vectNnStruct,nl%nnStruct)
    deallocate(vectNnStruct)

    if(.not.allocated(nl%nrnnStruct))allocate(nl%nrnnStruct(nats))
    nl%nrnnStruct = vectNrnnStruct
    deallocate(vectNrnnStruct)

    if(.not.allocated(nl%nrnnlist))allocate(nl%nrnnlist(nats))
    nl%nrnnlist = vectNrnnlist
    deallocate(vectNrnnlist)

#ifdef DO_MPI
    deallocate(rankRange)
#endif

  end subroutine gpmdcov_build_nlist_sparse

  !> Convert a vector to a matrix.
  !! \param rows Number of rows.
  !! \param cols Number of columns.
  !! \param vect Vector input.
  !! \param mat Matix output.
  subroutine vectorToMatrixIntLow(rows,cols,vect,mat)
    implicit none
    integer, allocatable, intent(in) :: vect(:)
    integer(low), allocatable, intent(inout) :: mat(:,:)
    integer, intent(in) :: rows, cols
    integer :: i,j
    do j = 1,cols
      do i = 1,rows
        mat(i,j) = vect((j-1)*rows + i)
      enddo
    enddo
  end subroutine vectorToMatrixIntLow

  !> Convert a vector to a matrix.
  !! \param rows Number of rows.
  !! \param cols Number of columns.
  !! \param vect Vector input.
  !! \param mat Matix output.
  subroutine vectorToMatrixInt(rows,cols,vect,mat)
    implicit none
    integer, allocatable, intent(in) :: vect(:)
    integer, allocatable, intent(inout) :: mat(:,:)
    integer, intent(in) :: rows, cols
    integer :: i,j
    do j = 1,cols
      do i = 1,rows
        mat(i,j) = vect((j-1)*rows + i)
      enddo
    enddo
  end subroutine vectorToMatrixInt

  !> Gets the volume of the simulation box
  !! \brief Given an array of lattice vectors, it return the box volume
  !! \param lattice_vector Lattice vectors in an array. latice_vectors(1,3) means the z-coordinate
  !! of the first lattice vector.
  !! \param volBox Volume of the cell.
  subroutine gpmdcov_get_vol(lattice_vectors,volBox)
    implicit none
    real(dp)                             ::  a1xa2(3), a2xa3(3), a3xa1(3)
    real(dp)                             ::  pi
    real(dp), intent(in)                 ::  lattice_vectors(:,:)
    real(dp), intent(inout)              ::  volBox

    volBox=0.0_dp

    pi = 3.14159265358979323846264338327950_dp

    a1xa2(1) = lattice_vectors(1,2)*lattice_vectors(2,3) - lattice_vectors(1,3)*lattice_vectors(2,2)
    a1xa2(2) = -lattice_vectors(1,1)*lattice_vectors(2,3) + lattice_vectors(1,3)*lattice_vectors(2,1)
    a1xa2(3) =  lattice_vectors(1,1)*lattice_vectors(2,2) - lattice_vectors(1,2)*lattice_vectors(2,1)

    a2xa3(1) = lattice_vectors(2,2)*lattice_vectors(3,3) - lattice_vectors(2,3)*lattice_vectors(3,2)
    a2xa3(2) = -lattice_vectors(2,1)*lattice_vectors(3,3) + lattice_vectors(2,3)*lattice_vectors(3,1)
    a2xa3(3) =  lattice_vectors(2,1)*lattice_vectors(3,2) - lattice_vectors(2,2)*lattice_vectors(3,1)

    a3xa1(1) = lattice_vectors(3,2)*lattice_vectors(1,3) - lattice_vectors(3,3)*lattice_vectors(1,2)
    a3xa1(2) = -lattice_vectors(3,1)*lattice_vectors(1,3) + lattice_vectors(3,3)*lattice_vectors(1,1)
    a3xa1(3) =  lattice_vectors(3,1)*lattice_vectors(1,2) - lattice_vectors(3,2)*lattice_vectors(1,1)

    !Get the volume of the cell
    volBox = lattice_vectors(1,1)*a2xa3(1)+ lattice_vectors(1,2)*a2xa3(2)+lattice_vectors(1,3)*a2xa3(3)

  end subroutine gpmdcov_get_vol



  !! \brief It will bild a neighbor list using an "all to all" approach
  !! \param coords System coordinates. coords(1,7): x-coordinate of atom 7.
  !! \param lattice_vectors. Lattice vectors of the system box. lattice_vectors(1,3): z-coordinate of vector 1.
  !! \param nl Neighbor list type.
  !! \param verbose Verbosity level.
  !! \param rank MPI rank
  subroutine gpmdcov_get_nlist_box_indices(coords,boxOfI,lattice_vectors,nx,ny,nz,verbose,rank,numranks)
    implicit none
    integer                              ::  NBox, cnt, i, ibox
    integer                              ::  ith, ix, iy, iz
    integer                              ::  j, jbox, jj, jxBox
    integer                              ::  jyBox, jzBox, maxInBox, maxNeigh
    integer                              ::  myNumranks, myrank, nats, natsPerRank
    integer, intent(in)                  ::  nx, ny, nz
    integer :: tx
    integer                              ::  ty, tz,nx1,ny1,nz1
    integer, allocatable                 ::  inbox(:,:), ithFromXYZ(:,:,:)
    integer, allocatable, intent(out)     ::  boxOfI(:)
    integer, allocatable                 ::  totPerBox(:), xBox(:), yBox(:), zBox(:)
    integer, intent(in)                  ::  verbose
    integer, optional, intent(in)        ::  numranks, rank
    real(dp)                             ::  coordsNeigh(3), density, distance, translation(3)
    real(dp)                             ::  volBox, minx, miny, minz, smallReal, mlsnl
    real(dp)                             ::  maxx, maxy, maxz
    real(dp), allocatable, intent(in)    ::  coords(:,:), lattice_vectors(:,:)
    real(dp)                 ::  rcutx,rcuty,rcutz
#ifdef DO_MPI
    integer, allocatable :: rankRange(:,:)
#endif

    if(present(rank).and.present(numranks))then
      myrank = rank
      myNumranks = numranks
    else
      myrank = 1
      myNumranks = 1
    endif

    call gpmdcov_msI("gpmdcov_build_neigborlist","Building neighbor list ...",verbose,myrank)

    nats = size(coords,dim=2) !Get the number of atoms
    natsPerRank = int(nats/myNumranks)

    !We will have approximatly [(4/3)*pi * rcut^3 * atomic density] number of neighbors.
    !A very large atomic density could be 1 atom per (1.0 Ang)^3 = 1 atoms per Ang^3  
    call gpmdcov_get_vol(lattice_vectors,volBox)
    density = 1.0_dp
    maxneigh = int(floor(3.14592_dp * (4.0_dp/3.0_dp) * density * (rcutx*rcuty*rcutz)))

    minx = 1.0d10
    miny = 1.0d10
    minz = 1.0d10
    maxx = -1.0d10
    maxy = -1.0d10
    maxz = -1.0d10
    do i = 1,nats
       minx = min(minx,coords(1,i))
       miny = min(miny,coords(2,i))
       minz = min(minz,coords(3,i))
       maxx= max(maxx,coords(1,i))
       maxy = max(maxy,coords(2,i))
       maxz = max(maxz,coords(3,i))
    enddo

    !We assume the box is orthogona
    rcutx = (maxx - minx)/(real(nx)) 
    rcuty = (maxy - miny)/(real(ny)) 
    rcutz = (maxz - minz)/(real(nz)) 

    NBox = nx*ny*nz
    maxInBox = int(density*(rcutx*rcuty*rcutz)) !Upper boud for the max number of atoms per box
    mlsnl = mls()
    allocate(inbox(NBox,maxInBox))
    inbox = 0
    allocate(totPerBox(Nbox))
    totPerBox = 0
    allocate(boxOfI(nats))
    boxOfI = 0
    allocate(xBox(Nbox))
    xBox = 0
    allocate(yBox(Nbox))
    yBox = 0
    allocate(zBox(Nbox))
    zBox = 0
    allocate(ithFromXYZ(nx,ny,nz))
    ithFromXYZ = 0

    smallReal = -0.0000001_dp
    !Search for the box coordinate and index of every atom
    do i = 1,nats
      !Index every atom respect to the discretized position on the simulation box.
      !tranlation = coords(:,i) - origin !For the general case we need to make sure coords ar > 0 
      ix = 1 + int(floor(abs((coords(1,i) - minx + smallReal)/(rcutx)))) !small box x-index of atom i
      iy = 1 + int(floor(abs((coords(2,i) - miny + smallReal)/(rcuty)))) !small box y-index //
      iz = 1 + int(floor(abs((coords(3,i) - minz + smallReal)/(rcutz)))) !small box z-index //

      write(*,*)coords(:,i),ix,iy,iz,(coords(3,i) - minz + smallReal)/rcutz
      if(ix > nx .or. ix < 0)then 
              write(*,*)"ix",ix
              Stop "Error in box index"
      endif
      if(iy > ny .or. iy < 0)then 
              write(*,*)"iy",iy
              Stop "Error in box index"
      endif
      if(iz > nz .or. iz < 0)then 
              write(*,*)"iz",iz
              Stop "Error in box index"
      endif

      ith =  ix + (iy-1)*(nx) + (iz-1)*(nx)*(ny)  !Get small box index
      boxOfI(i) = ith

      !From index to box coordinates
      xBox(ith) = ix
      yBox(ith) = iy
      zBox(ith) = iz

      !From box coordinates to index
      ithFromXYZ(ix,iy,iz) = ith

      totPerBox(ith) = totPerBox(ith) + 1 !How many per box
      if(totPerBox(ith) > maxInBox) Stop "Exceeding the max in box allowed"
      inbox(ith,totPerBox(ith)) = i !Who is in box ith
      write(*,*)ix,iy,iz,boxOfI(i),coords(:,i)

    enddo


    end subroutine gpmdcov_get_nlist_box_indices
 

  !! \brief It will bild a neighbor list using an "all to all" approach
  !! \param coords System coordinates. coords(1,7): x-coordinate of atom 7.
  !! \param lattice_vectors. Lattice vectors of the system box. lattice_vectors(1,3): z-coordinate of vector 1.
  !! \param nl Neighbor list type.
  !! \param verbose Verbosity level.
  !! \param rank MPI rank
  subroutine gpmdcov_build_nlist_sparse_v2(coords,lattice_vectors,rcut,nl,verbose,rank,numranks)
    implicit none
    integer                              ::  NBox, cnt, i, ibox
    integer                              ::  ith, ix, iy, iz
    integer                              ::  j, jbox, jj, jxBox
    integer                              ::  jyBox, jzBox, maxInBox, maxNeigh
    integer                              ::  myNumranks, myrank, nats, natsPerRank
    integer                              ::  nx, ny, nz, tx
    integer                              ::  ty, tz
    integer, allocatable                 ::  boxOfI(:), inbox(:,:), ithFromXYZ(:,:,:)
    integer, allocatable                 ::  totPerBox(:), xBox(:), yBox(:), zBox(:)
    integer, intent(in)                  ::  verbose
    integer, optional, intent(in)        ::  numranks, rank
    real(dp)                             ::  coordsNeigh(3), density, distance, translation(3)
    real(dp)                             ::  volBox, minx, miny, minz, smallReal, mlsnl
    real(dp), allocatable, intent(in)    ::  coords(:,:), lattice_vectors(:,:)
    real(dp), intent(in)                 ::  rcut
    type(neighlist_type), intent(inout)  ::  nl
#ifdef DO_MPI
    integer, allocatable :: rankRange(:,:)
#endif

    if(present(rank).and.present(numranks))then
      myrank = rank
      myNumranks = numranks
    else
      myrank = 1
      myNumranks = 1
    endif

    call gpmdcov_msI("gpmdcov_build_neigborlist","Building neighbor list ...",verbose,myrank)

    nats = size(coords,dim=2) !Get the number of atoms
    natsPerRank = int(nats/myNumranks)

    !We will have approximatly [(4/3)*pi * rcut^3 * atomic density] number of neighbors.
    !A very large atomic density could be 1 atom per (1.0 Ang)^3 = 1 atoms per Ang^3  
    call gpmdcov_get_vol(lattice_vectors,volBox)
    density = 1.0_dp 
    maxneigh = int(floor(3.14592_dp * (4.0_dp/3.0_dp) * density * rcut**3))

    !We assume the box is orthogonal
    nx = 1 + floor(lattice_vectors(1,1)/(2.0_dp*rcut))
    ny = 1 + floor(lattice_vectors(2,2)/(2.0_dp*rcut))
    nz = 1 + floor(lattice_vectors(3,3)/(2.0_dp*rcut))
    !write(*,*)"nx,ny,nz",nx,ny,nz
    !stop
    NBox = nx*ny*nz
    maxInBox = int(density*rcut**3) !Upper boud for the max number of atoms per box
    mlsnl = mls()
    allocate(inbox(NBox,maxInBox))
    inbox = 0
    allocate(totPerBox(Nbox))
    totPerBox = 0
    allocate(boxOfI(nats))
    boxOfI = 0
    allocate(xBox(Nbox))
    xBox = 0
    allocate(yBox(Nbox))
    yBox = 0
    allocate(zBox(Nbox))
    zBox = 0
    allocate(ithFromXYZ(nx,ny,nz))
    ithFromXYZ = 0

    minx = 1.0d10
    miny = 1.0d10
    minz = 1.0d10
    do i = 1,nats 
       minx = min(minx,coords(1,i))
       miny = min(miny,coords(2,i))
       minz = min(minz,coords(3,i))
    enddo

    smallReal = 0.0_dp
    !Search for the box coordinate and index of every atom
    do i = 1,nats
      !Index every atom respect to the discretized position on the simulation box.
      !tranlation = coords(:,i) - origin !For the general case we need to make sure coords ar > 0 
      ix = 1+ int(floor((coords(1,i) - minx + smallReal)/(2.0_dp*rcut))) !small box x-index of atom i
      iy = 1+ int(floor((coords(2,i) - miny + smallReal)/(2.0_dp*rcut))) !small box y-index //
      iz = 1+ int(floor((coords(3,i) - minz + smallReal)/(2.0_dp*rcut))) !small box z-index //
      
      if(ix > nx .or. ix < 0)Stop "Error in box index"
      if(iy > ny .or. iy < 0)Stop "Error in box index"
      if(iz > nz .or. iz < 0)Stop "Error in box index"

      ith =  ix + (iy-1)*(nx) + (iz-1)*(nx)*(ny)  !Get small box index
      boxOfI(i) = ith

      !From index to box coordinates
      xBox(ith) = ix
      yBox(ith) = iy
      zBox(ith) = iz

      !From box coordinates to index
      ithFromXYZ(ix,iy,iz) = ith

      totPerBox(ith) = totPerBox(ith) + 1 !How many per box
      if(totPerBox(ith) > maxInBox) Stop "Exceeding the max in box allowed"
      inbox(ith,totPerBox(ith)) = i !Who is in both ith

    enddo

    if(.not.allocated(nl%nnType))allocate(nl%nnType(maxneigh,nats))
    if(.not.allocated(nl%nnIx))allocate(nl%nnIx(maxneigh,nats))
    if(.not.allocated(nl%nnIy))allocate(nl%nnIy(maxneigh,nats))
    if(.not.allocated(nl%nnIz))allocate(nl%nnIz(maxneigh,nats))
    if(.not.allocated(nl%nnStruct))allocate(nl%nnStruct(maxneigh,nats))
    if(.not.allocated(nl%nrnnStruct))allocate(nl%nrnnStruct(nats))
    if(.not.allocated(nl%nrnnlist))allocate(nl%nrnnlist(nats))


     !For each atom we will look around to see who are its neighbors
    !$omp parallel do default(none) private(i) &
    !$omp private(ibox,ix,iy,iz) &
    !$omp private(jxbox,jybox,jzbox,jbox) &
    !$omp private(distance,translation,coordsNeigh) &
    !$omp private(cnt,j,jj,tx,ty,tz) &
    !$omp shared(nx,ny,nz,boxOfI) &
    !$omp shared(xBox,yBox,zBox) &
    !$omp shared(coords,rcut,totPerBox) &
    !$omp shared(nl,inbox,ithFromXYZ) &
    !$omp shared(lattice_vectors)&
    !$omp shared(maxneigh) &
    !$omp shared(nats)
    do i = 1,nats !For every atom
!#endif
      cnt = 0
      !Which box it beongs to
      ibox = boxOfI(i)
      !Look inside the box and the neighboring boxes
      do ix = -1,1
        do iy = -1,1
          do iz = -1,1
            !Get neigh box coordinate
            jxBox = xBox(ibox) + ix
            jyBox = yBox(ibox) + iy
            jzBox = zBox(ibox) + iz
            tx = 0.0_dp ; ty = 0.0_dp ; tz = 0.0_dp
            if(jxBox <= 0)then
              jxBox = nx
              tx = -1
            elseif(jxBox > nx)then
              jxBox = 1
              tx = 1
            endif
            if(jyBox <= 0)then
              jyBox = ny
              ty = -1
            elseif(jyBox > ny)then
              jyBox = 1
              ty = 1
            endif
            if(jzBox <= 0)then
              jzBox = nz
              tz = -1
            elseif(jzBox > nz)then
              jzBox = 1
              tz = 1
            endif

            !Get the neigh box index
            jbox = ithFromXYZ(jxBox,jyBox,jzBox)

            !Now loop over the atoms in the jbox
            do j = 1,totPerBox(jbox)
              jj = inbox(jbox,j) !Get atoms in box j
              translation = tx*lattice_vectors(1,:) + ty*lattice_vectors(2,:) + tz*lattice_vectors(3,:)
              coordsNeigh = coords(:,jj) + translation
              distance = norm2(coords(:,i) - coordsNeigh)
              if (distance .lt. rcut .and. distance .gt. 1d-12) then
                cnt = cnt + 1
                nl%Nntype(cnt,i) = jj ! jj is a neighbor of i by some translation
                nl%Nnstruct(cnt,i) = jj ! jj is a neighbor of i by some translation
                nl%NnIx(cnt,i) = tx
                nl%NnIy(cnt,i) = ty
                nl%NnIz(cnt,i) = tz
              endif
            enddo
          enddo
        enddo
      enddo
      nl%NrnnStruct(i) = cnt
      nl%Nrnnlist(i) = cnt

    enddo
    !$omp end parallel do

    deallocate(inbox)
    deallocate(totPerBox)
    deallocate(boxOfI)
    deallocate(xBox)
    deallocate(yBox)
    deallocate(zBox)
    deallocate(ithFromXYZ)

  end subroutine gpmdcov_build_nlist_sparse_v2


  !>  Build the neighbor list with a "brute force method"
  !! \brief It will bild a neighbor list using an "all to all" approach
  !! \param coords System coordinates. coords(1,7): x-coordinate of atom 7.
  !! \param lattice_vectors. Lattice vectors of the system box. lattice_vectors(1,3): z-coordinate of vector 1.
  !! \param nl Neighbor list type.
  !! \param verbose Verbosity level.
  !! \param rank MPI rank
  subroutine gpmdcov_build_nlist_sparse_v3(coords,lattice_vectors,rcut,nl,verbose,rank,numranks)
    implicit none
    integer                              ::  NBox, cnt, i, ibox
    integer                              ::  ith, ix, iy, iz
    integer                              ::  j, jbox, jj, jxBox
    integer                              ::  jyBox, jzBox, maxInBox, maxNeigh
    integer                              ::  myNumranks, myrank, nats, natsPerRank
    integer                              ::  nx, ny, nz, tx
    integer                              ::  ty, tz
    integer, allocatable                 ::  boxOfI(:), inbox(:,:), ithFromXYZ(:,:,:)
    integer, allocatable                 ::  totPerBox(:), xBox(:), yBox(:), zBox(:)
    integer, intent(in)                  ::  verbose
    integer, optional, intent(in)        ::  numranks, rank
    real(dp)                             ::  coordsNeigh(3), density, distance, translation(3)
    real(dp)                             ::  volBox, minx, miny, minz, smallReal, mlsnl
    real(dp), allocatable, intent(in)    ::  coords(:,:), lattice_vectors(:,:)
    real(dp), intent(in)                 ::  rcut
    type(neighlist_type), intent(inout)  ::  nl
#ifdef DO_MPI
    integer, allocatable :: rankRange(:,:)
#endif

    if(present(rank).and.present(numranks))then
      myrank = rank
      myNumranks = numranks
    else
      myrank = 1
      myNumranks = 1
    endif

    call gpmdcov_msI("gpmdcov_build_neigborlist","Building neighbor list ...",verbose,myrank)

    nats = size(coords,dim=2) !Get the number of atoms
    natsPerRank = int(nats/myNumranks)

#ifdef DO_MPI
    allocate(rankRange(2,myNumranks))
    do i = 1,myNumranks
      rankRange(1,i) = (i-1)*natsPerRank + 1
      rankRange(2,i) = i*natsPerRank
    enddo
    rankRange(2,myNumranks) = rankRange(2,myNumranks) + (nats - myNumranks*natsPerRank)
#endif


    !We will have approximatly [(4/3)*pi * rcut^3 * atomic density] number of neighbors.
    !A very large atomic density could be 1 atom per (1.0 Ang)^3 = 1 atoms per Ang^3  
    call gpmdcov_get_vol(lattice_vectors,volBox)
    density = 1.0_dp 
    maxneigh = int(floor(3.14592_dp * 4.0_dp/3.0_dp * density * rcut**3))

    !We assume the box is orthogonal
    nx = 1 + floor(lattice_vectors(1,1)/rcut)
    ny = 1 + floor(lattice_vectors(2,2)/rcut)
    nz = 1 + floor(lattice_vectors(3,3)/rcut)

    NBox = nx*ny*nz
    maxInBox = int(density*rcut**3) !Upper boud for the max number of atoms per box
    mlsnl = mls()
    allocate(inbox(NBox,maxInBox))
    inbox = 0
    allocate(totPerBox(Nbox))
    totPerBox = 0
    allocate(boxOfI(nats))
    boxOfI = 0
    allocate(xBox(Nbox))
    xBox = 0
    allocate(yBox(Nbox))
    yBox = 0
    allocate(zBox(Nbox))
    zBox = 0
    allocate(ithFromXYZ(nx,ny,nz))
    ithFromXYZ = 0

    minx = 1.0d10
    miny = 1.0d10
    minz = 1.0d10
    do i = 1,nats 
       minx = min(minx,coords(1,i))
       miny = min(miny,coords(2,i))
       minz = min(minz,coords(3,i))
    enddo

    smallReal = 0.0_dp
    write(*,*)"nlist Time allocs",mls() - mlsnl
    mlsnl = mls()
    !Search for the box coordinate and index of every atom
    do i = 1,nats
      !Index every atom respect to the discretized position on the simulation box.
      !tranlation = coords(:,i) - origin !For the general case we need to make sure coords ar > 0 
      ix = 1+ int(floor((coords(1,i) - minx + smallReal)/rcut)) !small box x-index of atom i
      iy = 1+ int(floor((coords(2,i) - miny + smallReal)/rcut)) !small box y-index //
      iz = 1+ int(floor((coords(3,i) - minz + smallReal)/rcut)) !small box z-index //
      
      if(ix > nx .or. ix < 0)Stop "Error in box index"
      if(iy > ny .or. iy < 0)Stop "Error in box index"
      if(iz > nz .or. iz < 0)Stop "Error in box index"

      ith =  ix + (iy-1)*(nx) + (iz-1)*(nx)*(ny)  !Get small box index
      boxOfI(i) = ith

      !From index to box coordinates
      xBox(ith) = ix
      yBox(ith) = iy
      zBox(ith) = iz

      !From box coordinates to index
      ithFromXYZ(ix,iy,iz) = ith

      totPerBox(ith) = totPerBox(ith) + 1 !How many per box
      if(totPerBox(ith) > maxInBox) Stop "Exceeding the max in box allowed"
      inbox(ith,totPerBox(ith)) = i !Who is in both ith

    enddo
    write(*,*)"nlist Time seting up boxes",mls() - mlsnl
    mlsnl = mls()

    if(.not.allocated(nl%nnType))allocate(nl%nnType(maxneigh,nats))
    if(.not.allocated(nl%nnIx))allocate(nl%nnIx(maxneigh,nats))
    if(.not.allocated(nl%nnIy))allocate(nl%nnIy(maxneigh,nats))
    if(.not.allocated(nl%nnIz))allocate(nl%nnIz(maxneigh,nats))
    if(.not.allocated(nl%nnStruct))allocate(nl%nnStruct(maxneigh,nats))
    if(.not.allocated(nl%nrnnStruct))allocate(nl%nrnnStruct(nats))
    if(.not.allocated(nl%nrnnlist))allocate(nl%nrnnlist(nats))


     !For each atom we will look around to see who are its neighbors
    !$omp parallel do default(none) private(i) &
    !$omp private(ibox,ix,iy,iz) &
    !$omp private(jxbox,jybox,jzbox,jbox) &
    !$omp private(distance,translation,coordsNeigh) &
    !$omp private(cnt,j,jj,tx,ty,tz) &
    !$omp shared(nx,ny,nz,inbox) &
    !$omp shared(coords,rcut,ithFromXYZ,boxOfI,totPerBox) &
    !$omp shared(nl,xBox,yBox,zBox) &
    !$omp shared(lattice_vectors)&
    !$omp shared(maxneigh) &
#ifdef DO_MPI
    !$omp shared(rankRange,myRank,nats)
    do i = rankRange(1,myRank),rankRange(2,myRank) !For every atom in the rank range
#else
    !$omp shared(nats)
    do i = 1,nats !For every atom
#endif
      cnt = 0
      !Which box it beongs to
      ibox = boxOfI(i)
      !Look inside the box and the neighboring boxes
      do ix = -1,1
        do iy = -1,1
          do iz = -1,1
            !Get neigh box coordinate
            jxBox = xBox(ibox) + ix
            jyBox = yBox(ibox) + iy
            jzBox = zBox(ibox) + iz
            tx = 0.0_dp ; ty = 0.0_dp ; tz = 0.0_dp
            if(jxBox <= 0)then
              jxBox = nx
              tx = -1
            elseif(jxBox > nx)then
              jxBox = 1
              tx = 1
            endif
            if(jyBox <= 0)then
              jyBox = ny
              ty = -1
            elseif(jyBox > ny)then
              jyBox = 1
              ty = 1
            endif
            if(jzBox <= 0)then
              jzBox = nz
              tz = -1
            elseif(jzBox > nz)then
              jzBox = 1
              tz = 1
            endif

            !Get the neigh box index
            jbox = ithFromXYZ(jxBox,jyBox,jzBox)

            !Now loop over the atoms in the jbox
            do j = 1,totPerBox(jbox)
              jj = inbox(jbox,j) !Get atoms in box j
              translation = tx*lattice_vectors(1,:) + ty*lattice_vectors(2,:) + tz*lattice_vectors(3,:)
              coordsNeigh = coords(:,jj) + translation
              distance = norm2(coords(:,i) - coordsNeigh)
              if (distance .lt. rcut .and. distance .gt. 1d-12) then
                cnt = cnt + 1
                nl%Nntype(cnt,i) = jj ! jj is a neighbor of i by some translation
                nl%Nnstruct(cnt,i) = jj ! jj is a neighbor of i by some translation
                nl%NnIx(cnt,i) = tx
                nl%NnIy(cnt,i) = ty
                nl%NnIz(cnt,i) = tz
              endif
            enddo
          enddo
        enddo
      enddo
      nl%NrnnStruct(i) = cnt
      nl%Nrnnlist(i) = cnt

    enddo
    !$omp end parallel do
    write(*,*)"nlist Time nlist rsearch",mls() - mlsnl

    deallocate(inbox)
    deallocate(totPerBox)
    deallocate(boxOfI)
    deallocate(xBox)
    deallocate(yBox)
    deallocate(zBox)
    deallocate(ithFromXYZ)

        !We do a sum reduction on all the vectors
#ifdef DO_MPI
    call prg_sumIntReduceN(nl%Nntype,maxneigh*nats)
    call prg_sumIntReduceN(nl%Nnstruct,maxneigh*nats)
    call prg_sumIntReduceN(nl%NnIx,maxneigh*nats)
    call prg_sumIntReduceN(nl%NnIy,maxneigh*nats)
    call prg_sumIntReduceN(nl%NnIz,maxneigh*nats)
    call prg_sumIntReduceN(nl%NrnnStruct,nats)
    call prg_sumIntReduceN(nl%Nrnnlist,nats)
#endif

    write(*,*)"nlist Time for transfer",mls() - mlsnl

  end subroutine gpmdcov_build_nlist_sparse_v3


  subroutine gpmd_nearestneighborlist(nrnnlist,nndist,nnRx,nnRy,nnRz,nnType, &
                            &R_X,R_Y,R_Z,LBox,Rcut,Nr_atoms,Max_Nr_Neigh)
implicit none

integer,    parameter       ::  MSkin = 2
real(dp), parameter       :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0
integer,    intent(in)      :: Nr_atoms, Max_Nr_Neigh
real(dp), intent(in)      :: Rcut
real(dp), intent(in)      :: R_X(Nr_atoms), R_Y(Nr_atoms), R_Z(Nr_atoms), LBox(3)
real(dp)                  :: Rx(MSkin*Nr_atoms), Ry(MSkin*Nr_atoms), Rz(MSkin*Nr_atoms)
integer,    intent(out)     :: nrnnlist(Nr_atoms), nnType(Max_Nr_Neigh,Nr_atoms)
real(dp), intent(out)     :: nndist(Max_Nr_Neigh,Nr_atoms), nnRx(Max_Nr_Neigh,Nr_atoms)
real(dp), intent(out)     :: nnRy(Max_Nr_Neigh,Nr_atoms), nnRz(Max_Nr_Neigh,Nr_atoms)
integer       :: i,j,k,l,m,t,nx,ny,nz,cell,type(10*Nr_atoms),Nskin,N
integer       :: head(Nr_atoms*MSkin), list(MSkin*Nr_atoms), tmp(Nr_atoms), cnt, cct
real(dp)    :: Lx,Ly,Lz,Tx,Ty,Tz,RR(3),TT(3),dist,dLx,dLy,dLz

real(dp) :: start, finish

start = mls()

N = Nr_atoms
Lx = LBox(1)
Ly = LBox(2)
Lz = LBox(3) ! Dimensions of periodic BC
nx = floor(Lx/Rcut)
ny = floor(Ly/Rcut)
nz = floor(Lz/Rcut) ! Division into # cell boxes: nx, ny, nz
Rx(1:Nr_atoms) = R_X(1:Nr_atoms) !+Lx
Ry(1:Nr_atoms) = R_Y(1:Nr_atoms) !+Ly
Rz(1:Nr_atoms) = R_Z(1:Nr_atoms) !+Lz

if ((min(nx,ny,nz).lt.3).or.(Nr_atoms.lt.80)) then

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,l,m,RR,TT,Tx,Ty,Tz,dist,tmp,cnt)
 do i = 1,N
   cnt = 0
   tmp = ZERO
   RR(1) = Rx(i)
   RR(2) = Ry(i)
   RR(3) = Rz(i)
   do m = 1,N
     do j = -1,1
     do k = -1,1
     do l = -1,1
       Tx = Rx(m)+j*Lx  ! Search all neigbors within a single translation (multiple translations could be necessary for small systems!
       Ty = Ry(m)+k*Ly
       Tz = Rz(m)+l*Lz
       TT(1) = Tx
       TT(2) = Ty
       TT(3) = Tz
       dist = norm2(RR-TT)
       !if ((dist.lt.Rcut).and.(dist.gt.1e-12)) then ! Neighbors within Rcut inlcuidng translated atoms in the "skin"
       if ((dist < Rcut)) then
         cnt = cnt + 1
         nndist(cnt,i) = dist
         nnRx(cnt,i) = Tx
         nnRy(cnt,i) = Ty
         nnRz(cnt,i) = Tz
         nnType(cnt,i) = m  ! Neigbor is number of original ordering number m in the box that might have been stranslated to the skin
         tmp(m) = m
       endif
     enddo
     enddo
     enddo
   enddo
   nrnnlist(i) = cnt
 enddo
!$OMP END PARALLEL DO

else

 head = ZERO ! Linked list that keeps track of all atoms in the nx*ny*nz small boxes
 list = ZERO !
 do i = 1,N
  cell = 1 + floor(nx*Rx(i)/Lx) + floor(ny*Ry(i)/Ly)*nx + floor(nz*Rz(i)/Lz)*nx*ny
  list(i) = head(cell)
  head(cell) = i
  type(i) = i
 enddo

 !%%% And now add a skin or surface buffer to account for periodic BC, all 26 of them!
 cnt = 0
 do i = 1,nx*ny  ! All boxes in the first (z=0) layer
   t = head(i)
   do while (t.gt.0)    ! and all atoms of this first layer
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)
      Ry(N+cnt) = Ry(t)
      Rz(N+cnt) = Rz(t)+Lz ! Add skin atoms with coordinates translated by Lz
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = 1,nx*ny*nz,nx  ! All boxes in another (x=0) layer
   t = head(i)
   do while (t.gt.0)   ! and all atoms in that layer
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)+Lx ! Add skin atoms with coordinates translated by Lx
      Ry(N+cnt) = Ry(t)
      Rz(N+cnt) = Rz(t)
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = 1,nx*ny*nz,nx*ny  ! Continue ...
   do k = i,i+nx-1
     t = head(k)
     do while (t.gt.0)
        cnt = cnt + 1
        Rx(N+cnt) = Rx(t)
        Ry(N+cnt) = Ry(t)+Ly
        Rz(N+cnt) = Rz(t)
        type(N+cnt) = t
        t = list(t)
     enddo
   enddo
 enddo
 cct = 0
 do i = nx*ny*(nz-1)+1,nx*ny*nz
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      cct = cct + 1
      Rx(N+cnt) = Rx(t)
      Ry(N+cnt) = Ry(t)
      Rz(N+cnt) = Rz(t)-Lz
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = nx,nx*ny*nz,nx
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)-Lx
      Ry(N+cnt) = Ry(t)
      Rz(N+cnt) = Rz(t)
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = nx*(ny-1)+1,nx*ny*nz,nx*ny
   do k = i,i+nx-1
     t = head(k)
     do while (t.gt.0)
        cnt = cnt + 1
        Rx(N+cnt) = Rx(t)
        Ry(N+cnt) = Ry(t)-Ly
        Rz(N+cnt) = Rz(t)
        type(N+cnt) = t
        t = list(t)
     enddo
   enddo
 enddo
 do i = 1,nx
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)
      Ry(N+cnt) = Ry(t)+Ly
      Rz(N+cnt) = Rz(t)+Lz
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = nx*(ny-1)+1,nx*(ny-1)+nx
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)
      Ry(N+cnt) = Ry(t)-Ly
      Rz(N+cnt) = Rz(t)+Lz
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = nx*ny*(nz-1)+1,nx*ny*(nz-1)+nx
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)
      Ry(N+cnt) = Ry(t)+Ly
      Rz(N+cnt) = Rz(t)-Lz
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = nx*ny*(nz-1)+nx*(ny-1)+1,nx*ny*(nz-1)+nx*(ny-1)+nx
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)
      Ry(N+cnt) = Ry(t)-Ly
      Rz(N+cnt) = Rz(t)-Lz
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = 1,nx*ny,nx
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)+Lx
      Ry(N+cnt) = Ry(t)
      Rz(N+cnt) = Rz(t)+Lz
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = nx*ny*(nz-1)+1,nx*ny*nz,nx
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)+Lx
      Ry(N+cnt) = Ry(t)
      Rz(N+cnt) = Rz(t)-Lz
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = nx*ny*(nz-1)+nx,nx*ny*nz,nx
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)-Lx
      Ry(N+cnt) = Ry(t)
      Rz(N+cnt) = Rz(t)-Lz
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = nx,nx*ny,nx
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)-Lx
      Ry(N+cnt) = Ry(t)
      Rz(N+cnt) = Rz(t)+Lz
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = 1,nx*ny*nz,nx*ny
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)+Lx
      Ry(N+cnt) = Ry(t)+Ly
      Rz(N+cnt) = Rz(t)
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = nx*(ny-1)+1,nx*ny*nz,nx*ny
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)+Lx
      Ry(N+cnt) = Ry(t)-Ly
      Rz(N+cnt) = Rz(t)
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = nx*(ny-1)+nx,nx*ny*nz,nx*ny
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)-Lx
      Ry(N+cnt) = Ry(t)-Ly
      Rz(N+cnt) = Rz(t)
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 do i = nx,nx*ny*nz,nx*ny
   t = head(i)
   do while (t.gt.0)
      cnt = cnt + 1
      Rx(N+cnt) = Rx(t)-Lx
      Ry(N+cnt) = Ry(t)+Ly
      Rz(N+cnt) = Rz(t)
      type(N+cnt) = t
      t = list(t)
   enddo
 enddo
 t = head(1)
 do while (t.gt.0)
    cnt = cnt + 1
    Rx(N+cnt) = Rx(t)+Lx
    Ry(N+cnt) = Ry(t)+Ly
    Rz(N+cnt) = Rz(t)+Lz
    type(N+cnt) = t
    t = list(t)
 enddo
 t = head(nx)
 do while (t.gt.0)
    cnt = cnt + 1
    Rx(N+cnt) = Rx(t)-Lx
    Ry(N+cnt) = Ry(t)+Ly
    Rz(N+cnt) = Rz(t)+Lz
    type(N+cnt) = t
    t = list(t)
 enddo
 t = head(nx*(ny-1)+1)
 do while (t.gt.0)
    cnt = cnt + 1
    Rx(N+cnt) = Rx(t)+Lx
    Ry(N+cnt) = Ry(t)-Ly
    Rz(N+cnt) = Rz(t)+Lz
    type(N+cnt) = t
    t = list(t)
 enddo
 t = head(nx*(ny-1)+nx)
 do while (t.gt.0)
    cnt = cnt + 1
    Rx(N+cnt) = Rx(t)-Lx
    Ry(N+cnt) = Ry(t)-Ly
    Rz(N+cnt) = Rz(t)+Lz
    type(N+cnt) = t
    t = list(t)
 enddo
 t = head(nx*ny*(nz-1)+1)
 do while (t.gt.0)
    cnt = cnt + 1
    Rx(N+cnt) = Rx(t)+Lx
    Ry(N+cnt) = Ry(t)+Ly
    Rz(N+cnt) = Rz(t)-Lz
    type(N+cnt) = t
    t = list(t)
 enddo
 t = head(nx*ny*(nz-1)+nx)
 do while (t.gt.0)
    cnt = cnt + 1
    Rx(N+cnt) = Rx(t)-Lx
    Ry(N+cnt) = Ry(t)+Ly
    Rz(N+cnt) = Rz(t)-Lz
    type(N+cnt) = t
    t = list(t)
 enddo
 t = head(nx*ny*(nz-1)+nx*(ny-1)+1)
 do while (t.gt.0)
    cnt = cnt + 1
    Rx(N+cnt) = Rx(t)+Lx
    Ry(N+cnt) = Ry(t)-Ly
    Rz(N+cnt) = Rz(t)-Lz
    type(N+cnt) = t
    t = list(t)
 enddo
 t = head(nx*ny*nz)
 do while (t.gt.0)
    cnt = cnt + 1
    Rx(N+cnt) = Rx(t)-Lx
    Ry(N+cnt) = Ry(t)-Ly
    Rz(N+cnt) = Rz(t)-Lz
    type(N+cnt) = t
    t = list(t)
 enddo
 Nskin = cnt

 !! And now create a list for everything including the Skin/buffer layer
! [max(Rz),min(Rz)]
 dLx = Lx/nx
 dLy = Ly/ny
 dLz = Lz/nz
 Rx = Rx+dLx
 Ry = Ry+dLy
 Rz = Rz+dLz ! Shift to avoid negative coordinates
 Lx = Lx + 2*dLx
 Ly = Ly + 2*dLy
 Lz = Lz + 2*dLz
 nx = nx+2
 ny = ny+2
 nz = nz+2
 head = ZERO
 list = ZERO
 do i = 1,N+Nskin
  cell = 1 + floor(nx*Rx(i)/Lx) + floor(ny*Ry(i)/Ly)*nx + floor(nz*Rz(i)/Lz)*nx*ny
  list(i) = head(cell)
  head(cell) = i
 enddo

!! OPENMP LOOP, JUST LIKE IN CHRISTIAN'S CODE, ESSENTIALLY ALL TIME AND WORK IS SPENT HERE
 do i = 1,N ! Go through all atoms
   cnt = 0
   do j = -1,1  ! Translate position to neighboring small boxes
   do k = -1,1
   do l = -1,1
     Tx = Rx(i)+j*dLx
     Ty = Ry(i)+k*dLy
     Tz = Rz(i)+l*dLz
     cell = 1 + floor(nx*Tx/Lx) + floor(ny*Ty/Ly)*nx + floor(nz*Tz/Lz)*nx*ny ! and extract all atoms in those small
     t = head(cell)                                                          ! neighbor boxes ...
     do while (t.gt.0)
       RR(1) = Rx(i) - Rx(t)
       RR(2) = Ry(i) - Ry(t)
       RR(3) = Rz(i) - Rz(t)
       dist = norm2(RR)
       !if (dist.lt.Rcut) then                       ! All atoms including the skin within Rcut WITH ITSELF!
       if (dist.lt.Rcut .and. dist .gt. 1d-12) then  ! All atoms including the skin within Rcut WITHOUT ITSELF!
         cnt = cnt + 1
         nndist(cnt,i) = dist
         nnRx(cnt,i) = Rx(t)-dLx  ! Coordinates without the shift
         nnRy(cnt,i) = Ry(t)-dLy
         nnRz(cnt,i) = Rz(t)-dLz
         nnType(cnt,i) = type(t)
       endif
       t = list(t)
     enddo
   enddo
   enddo
   enddo
   nrnnlist(i) = cnt
 enddo
 finish = mls()
print '("Time for last nn_list = ",f6.3," seconds.")',finish-start
endif

end subroutine gpmd_nearestneighborlist


end module gpmdcov_neighbor_mod
