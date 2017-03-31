!> A module to read and handle a the nearest neighbor list.
!! \brief This module will construct a neigbor list for every atom.
!! @ingroup LATTE 
!!
!! \todo THIS ROUTINES NEEDS TO RUN PARALLEL
!! 
module neighborlist_latte_mod

  use prg_openfiles_mod
  use prg_ptable_mod

  implicit none     

  private 

  integer, parameter :: dp = kind(1.0d0)

  !> System type
  type, public :: neighlist_type  !< The molecular system type.

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
    integer(1), allocatable :: nnIx(:,:)
    !> y-integer translation of neighbor J to I within RCut (including atoms in the skin)
    integer(1), allocatable :: nnIy(:,:) 
    !> z-integer translation of neighbor J to I within RCut (including atoms in the skin)
    integer(1), allocatable :: nnIz(:,:)          

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

  public :: build_nlist, build_nlist_int, destroy_nlist

contains


  !> Destroy the neigbor list to recover memory. 
  !! \param nl Neigbor list structure.
  !!
  subroutine destroy_nlist(nl)
    implicit none
    type(neighlist_type), intent(inout)  ::  nl
    
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
     
  end subroutine destroy_nlist
  
    
  !> The neigbor list construction. 
  !! \param coords system coordinates for which neighbor list should be constructed.
  !! \param lattice_vectors lattice vectors of the system.
  !! \param rcut coulomb cut off radius.
  !! \param nl neighbor list type. 
  !! \param verbose verbosity level.
  !!
  !! WARNING: This list only works for orthogonal lattice_vectors
  !! \todo Generalize neighbor list construction to nonorthogonal lattice vectors.
  !! 
  !!
  subroutine build_nlist_int(coords,lattice_vectors,rcut,nl,verbose)
    implicit none
    integer                              ::  cell, cnt, cnt2, i
    integer                              ::  j, k, l, m
    integer                              ::  nats, natspblock, nx, ny
    integer                              ::  nz, Nskin, ss, t
    integer, intent(in)                  ::  verbose 
    integer, allocatable                 ::  ntype(:), tmp(:)
    real(dp)                             ::  Lx, Ly, Lz, Tx
    real(dp)                             ::  dLx, dLy, dLz
    real(dp)                             ::  Ty, Tz, calpha, coulcut
    real(dp)                             ::  coulvol, dist, pi, sqrtx
    real(dp), allocatable                ::  head(:), list(:), buffer(:,:), distvec(:), trtmp(:,:)
    real(dp), intent(in)                 ::  coords(:,:), lattice_vectors(:,:), rcut
    type(neighlist_type), intent(inout)  ::  nl
    logical(1)                           ::  found

    nats = size(coords,dim=2)
    pi = 3.14159265358979323846264338327950_dp

    Lx = lattice_vectors(1,1)
    Ly = lattice_vectors(2,2)
    Lz = lattice_vectors(3,3)

    ! Division into # cell boxes: nx, ny, nz
    nx = floor(Lx/Rcut)
    ny = floor(Ly/Rcut)
    nz = floor(Lz/Rcut)

    ! Asuming an upper bound average density of 1 atom/ang^3 
    natspblock = floor(0.5d0*Rcut**3)

    if(verbose >= 0)  write(*,*) "In build_nlist ..."
    if(verbose >= 1)  write(*,*) "Number of atoms per block = ",natspblock
    if(verbose >= 1)  write(*,*) "min(nx,ny,nz) =", min(nx,ny,nz)

    if(.not.(allocated(nl%nrnnlist)))then
      if(min(Lx,Ly,Lz)/2.0_dp < rcut)then 
        allocate(nl%nnIx(natspblock,nats));
        allocate(nl%nnIy(natspblock,nats));
        allocate(nl%nnIz(natspblock,nats));
      endif
      allocate(nl%nnType(natspblock,nats))
      allocate(nl%nnStruct(natspblock,nats))
      allocate(nl%nrnnStruct(nats))
      allocate(nl%nrnnlist(nats))      
    endif  

    if(min(nx,ny,nz) < 1000)then   ! Brute force for small systems!

      if(verbose >= 1)  write(*,*) "Performing brute force for small system ..."
      
      allocate(tmp(nats));          

      !$omp parallel do default(none) private(i) &
      !$omp private(cnt,cnt2,found,tmp,m,j,k,l,ss,Tx,Ty,Tz,dist) &
      !$omp shared(nats,coords,Lx,Ly,Lz,Rcut,nl)
      do i = 1,nats
        cnt = 0
        tmp = 0
        do m = 1,nats
          found = .false.
          do j = -1,1
            do k = -1,1
              do l = -1,1
                Tx = coords(1,m)+j*Lx;
                Ty = coords(2,m)+k*Ly;
                Tz = coords(3,m)+l*Lz;
                dist = (coords(1,i)-Tx)**2 + (coords(2,i)-Ty)**2 + (coords(3,i)-Tz)**2
                dist = sqrt(dist)
                if (dist .lt. Rcut .and. dist .gt. 1d-12) then 
                  cnt = cnt + 1
                  nl%nnType(cnt,i) = m
                  tmp(m) = m          
                  if(allocated(nl%nnIx))then
                    nl%nnIx(cnt,i) = j
                    nl%nnIy(cnt,i) = k
                    nl%nnIz(cnt,i) = l
                  else
                    found = .true.                    
                  endif                                              
                endif
                if(found .eqv. .true.)exit
              enddo
              if(found .eqv. .true.)exit
            enddo
            if(found .eqv. .true.)exit
          enddo

        enddo

        nl%nrnnlist(i) = cnt;
        cnt2 = 0;

        do ss = 1,nats
          if (tmp(ss) .gt. 0)then 
            cnt2 = cnt2 + 1
            nl%nnStruct(cnt2,i) = ss
          endif
        enddo
        nl%nrnnStruct(i) = cnt2
      enddo
      !$omp end parallel do

      deallocate(tmp)

    else ! Do the same but now with linked lists in O(N)

      allocate(head(nx*ny*nz));    
      allocate(list(nats));        
      allocate(ntype(10*nats));
      allocate(buffer(3,nats*26)) !Allocate max buffer atoms
      allocate(trtmp(3,nats*26)) !Allocate max buffer atoms

      buffer(:,1:nats) = coords

      head = 0 
      list = 0 
      do i = 1,nats  
        cell = 1 + floor(nx*coords(1,i)/Lx) + floor(ny*coords(2,i)/Ly)*nx &
          + floor(nz*coords(3,i)/Lz)*nx*ny;
        list(i) = head(cell);
        head(cell) = i;
        ntype(i) = i;        
      enddo

      !And now add a skin or surface buffer to account for periodic BC, all 26 of them!
      cnt = 0;
      do i = 1,nx*ny  !All boxes in the first (z=0) layer
        t = head(i);
        do while (t > 0) ! and all atoms of this first layer
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t);
          buffer(2,nats+cnt) = coords(2,t);
          buffer(3,nats+cnt) = coords(3,t)+Lz;
          ntype(nats+cnt) = t;
          trtmp(1,nats+cnt) = 0 ; trtmp(2,nats+cnt) = 0; trtmp(3,nats+cnt) = 1
          t = list(t);          
        enddo
      enddo

      do i = 1,nx*ny*nz,nx
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t)+Lx;
          buffer(2,nats+cnt) = coords(2,t);
          buffer(3,nats+cnt) = coords(3,t);
          ntype(nats+cnt) = t;
          trtmp(1,nats+cnt) = 1 ; trtmp(2,nats+cnt) = 0; trtmp(3,nats+cnt) = 0          
          t = list(t);
        enddo
      enddo

      do i = 1,nx*ny*nz,nx*ny
        do k = i,i+nx-1
          t = head(k);
          do while(t > 0)
            cnt = cnt + 1;
            buffer(1,nats+cnt) = coords(1,t);
            buffer(2,nats+cnt) = coords(2,t)+Ly;
            buffer(3,nats+cnt) = coords(3,t);
            ntype(nats+cnt) = t;
            trtmp(1,nats+cnt) = 0 ; trtmp(2,nats+cnt) = 1; trtmp(3,nats+cnt) = 0                      
            t = list(t);
          enddo
        enddo
      enddo

      do i = nx*ny*(nz-1)+1,nx*ny*nz
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t);
          buffer(2,nats+cnt) = coords(2,t);
          buffer(3,nats+cnt) = coords(3,t)-Lz;
          ntype(nats+cnt) = t;
          trtmp(1,nats+cnt) = 0 ; trtmp(2,nats+cnt) = 0; trtmp(3,nats+cnt) = -1                    
          t = list(t);
        enddo
      enddo

      do i = nx,nx*ny*nz,nx
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t)-Lx;
          buffer(2,nats+cnt) = coords(2,t);
          buffer(3,nats+cnt) = coords(3,t);
          ntype(nats+cnt) = t;
          trtmp(1,nats+cnt) = -1 ; trtmp(2,nats+cnt) = 0; trtmp(3,nats+cnt) = 0                    
          t = list(t);
        enddo
      enddo

      do i = nx*(ny-1)+1,nx*ny*nz,nx*ny
        do k = i,i+nx-1
          t = head(k);
          do while(t > 0)
            cnt = cnt + 1;
            buffer(1,nats+cnt) = coords(1,t);
            buffer(2,nats+cnt) = coords(2,t)-Ly;
            buffer(3,nats+cnt) = coords(3,t);
            ntype(nats+cnt) = t;
            trtmp(1,nats+cnt) = 0 ; trtmp(2,nats+cnt) = -1; trtmp(3,nats+cnt) = 0                      
            t = list(t);
          enddo
        enddo
      enddo

      do i = 1,nx
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t);
          buffer(2,nats+cnt) = coords(2,t)+Ly;
          buffer(3,nats+cnt) = coords(3,t)+Lz;
          ntype(nats+cnt) = t;
          trtmp(1,nats+cnt) = 0 ; trtmp(2,nats+cnt) = 1; trtmp(3,nats+cnt) = 1                    
          t = list(t);
        enddo
      enddo

      do i = nx*(ny-1)+1,nx*(ny-1)+nx
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t);
          buffer(2,nats+cnt) = coords(2,t)-Ly;
          buffer(3,nats+cnt) = coords(3,t)+Lz;
          ntype(nats+cnt) = t;
          trtmp(1,nats+cnt) = 0 ; trtmp(2,nats+cnt) = -1; trtmp(3,nats+cnt) = 1          
          t = list(t);
        enddo
      enddo

      do i = nx*ny*(nz-1)+1,nx*ny*(nz-1)+nx
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t);
          buffer(2,nats+cnt) = coords(2,t)+Ly;
          buffer(3,nats+cnt) = coords(3,t)-Lz;
          ntype(nats+cnt) = t;
          trtmp(1,nats+cnt) = 0 ; trtmp(2,nats+cnt) = 1; trtmp(3,nats+cnt) = -1                    
          t = list(t);
        enddo
      enddo

      do i = nx*ny*(nz-1)+nx*(ny-1)+1,nx*ny*(nz-1)+nx*(ny-1)+nx
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t);
          buffer(2,nats+cnt) = coords(2,t)-Ly;
          buffer(3,nats+cnt) = coords(3,t)-Lz;
          trtmp(1,nats+cnt) = 0 ; trtmp(2,nats+cnt) = -1; trtmp(3,nats+cnt) = -1                    
          ntype(nats+cnt) = t;
          t = list(t);
        enddo
      enddo

      do i = 1,nx*ny,nx
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t)+Lx;
          buffer(2,nats+cnt) = coords(2,t);
          buffer(3,nats+cnt) = coords(3,t)+Lz;
          ntype(nats+cnt) = t;
          trtmp(1,nats+cnt) = 1 ; trtmp(2,nats+cnt) = 0; trtmp(3,nats+cnt) = 1          
          t = list(t);
        enddo
      enddo

      do i = nx*ny*(nz-1)+1,nx*ny*nz,nx
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t)+Lx;
          buffer(2,nats+cnt) = coords(2,t);
          buffer(3,nats+cnt) = coords(3,t)-Lz;
          ntype(nats+cnt) = t;
          trtmp(1,nats+cnt) = 1 ; trtmp(2,nats+cnt) = 0; trtmp(3,nats+cnt) = -1                    
          t = list(t);
        enddo
      enddo

      do i = nx*ny*(nz-1)+nx,nx*ny*nz,nx
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t)-Lx;
          buffer(2,nats+cnt) = coords(2,t);
          buffer(3,nats+cnt) = coords(3,t)-Lz;
          trtmp(1,nats+cnt) = -1 ; trtmp(2,nats+cnt) = 0; trtmp(3,nats+cnt) = -1                    
          ntype(nats+cnt) = t;
          t = list(t);
        enddo
      enddo

      do i = nx,nx*ny,nx
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t)-Lx;
          buffer(2,nats+cnt) = coords(2,t);
          buffer(3,nats+cnt) = coords(3,t)+Lz;
          ntype(nats+cnt) = t;
          trtmp(1,nats+cnt) = -1 ; trtmp(2,nats+cnt) = 0; trtmp(3,nats+cnt) = 1                    
          t = list(t);
        enddo
      enddo

      do i = 1,nx*ny*nz,nx*ny
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t)+Lx;
          buffer(2,nats+cnt) = coords(2,t)+Ly;
          buffer(3,nats+cnt) = coords(3,t);
          ntype(nats+cnt) = t;
          trtmp(1,nats+cnt) = 1 ; trtmp(2,nats+cnt) = 1; trtmp(3,nats+cnt) = 0                    
          t = list(t);
        enddo
      enddo

      do i = nx*(ny-1)+1,nx*ny*nz,nx*ny
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t)+Lx;
          buffer(2,nats+cnt) = coords(2,t)-Ly;
          buffer(3,nats+cnt) = coords(3,t);
          trtmp(1,nats+cnt) = 1 ; trtmp(2,nats+cnt) = -1; trtmp(3,nats+cnt) = 0                    
          ntype(nats+cnt) = t;
          t = list(t);
        enddo
      enddo

      do i = nx*(ny-1)+nx,nx*ny*nz,nx*ny
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t)-Lx;
          buffer(2,nats+cnt) = coords(2,t)-Ly;
          buffer(3,nats+cnt) = coords(3,t);
          trtmp(1,nats+cnt) = -1 ; trtmp(2,nats+cnt) = -1; trtmp(3,nats+cnt) = 0                    
          ntype(nats+cnt) = t;
          t = list(t);
        enddo
      enddo

      do i = nx,nx*ny*nz,nx*ny
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t)-Lx;
          buffer(2,nats+cnt) = coords(2,t)+Ly;          
          buffer(3,nats+cnt) = coords(3,t);
          trtmp(1,nats+cnt) = -1 ; trtmp(2,nats+cnt) = 1; trtmp(3,nats+cnt) = 0                    
          ntype(nats+cnt) = t;
          t = list(t);
        enddo
      enddo

      t = head(1);
      do while(t > 0)
        cnt = cnt + 1;
        buffer(1,nats+cnt) = coords(1,t)+Lx;
        buffer(2,nats+cnt) = coords(2,t)+Ly;
        buffer(3,nats+cnt) = coords(3,t)+Lz;
        trtmp(1,nats+cnt) = 1 ; trtmp(2,nats+cnt) = 1; trtmp(3,nats+cnt) = 1                  
        ntype(nats+cnt) = t;
        t = list(t);
      enddo

      t = head(nx);
      do while(t > 0)
        cnt = cnt + 1;
        buffer(1,nats+cnt) = coords(1,t)-Lx;
        buffer(2,nats+cnt) = coords(2,t)+Ly;
        buffer(3,nats+cnt) = coords(3,t)+Lz;
        ntype(nats+cnt) = t;
        trtmp(1,nats+cnt) = -1 ; trtmp(2,nats+cnt) = 1; trtmp(3,nats+cnt) = 1                  
        t = list(t);
      enddo

      t = head(nx*(ny-1)+1);
      do while(t > 0)
        cnt = cnt + 1;
        buffer(1,nats+cnt) = coords(1,t)+Lx;
        buffer(2,nats+cnt) = coords(2,t)-Ly;
        buffer(3,nats+cnt) = coords(3,t)+Lz;
        ntype(nats+cnt) = t;
        trtmp(1,nats+cnt) = 1 ; trtmp(2,nats+cnt) = -1; trtmp(3,nats+cnt) = 1                  
        t = list(t);
      enddo

      t = head(nx*(ny-1)+nx);
      do while(t > 0)
        cnt = cnt + 1;
        buffer(1,nats+cnt) = coords(1,t)-Lx;
        buffer(2,nats+cnt) = coords(2,t)-Ly;
        buffer(3,nats+cnt) = coords(3,t)+Lz;
        ntype(nats+cnt) = t;
        trtmp(1,nats+cnt) = -1 ; trtmp(2,nats+cnt) = -1; trtmp(3,nats+cnt) = 1                  
        t = list(t);
      enddo

      t = head(nx*ny*(nz-1)+1);
      do while(t > 0)
        cnt = cnt + 1;
        buffer(1,nats+cnt) = coords(1,t)+Lx;
        buffer(2,nats+cnt) = coords(2,t)+Ly;
        buffer(3,nats+cnt) = coords(3,t)-Lz;
        ntype(nats+cnt) = t;
        trtmp(1,nats+cnt) = 1 ; trtmp(2,nats+cnt) = 1; trtmp(3,nats+cnt) = -1                  
        t = list(t);
      enddo

      t = head(nx*ny*(nz-1)+nx);
      do while(t > 0)
        cnt = cnt + 1;
        buffer(1,nats+cnt) = coords(1,t)-Lx;
        buffer(2,nats+cnt) = coords(2,t)+Ly;
        buffer(3,nats+cnt) = coords(3,t)-Lz;
        ntype(nats+cnt) = t;
        trtmp(1,nats+cnt) = -1 ; trtmp(2,nats+cnt) = 1; trtmp(3,nats+cnt) = -1                  
        t = list(t);
      enddo

      t = head(nx*ny*(nz-1)+nx*(ny-1)+1);
      do while(t > 0)
        cnt = cnt + 1;
        buffer(1,nats+cnt) = coords(1,t)+Lx;
        buffer(2,nats+cnt) = coords(2,t)-Ly;
        buffer(3,nats+cnt) = coords(3,t)-Lz;
        ntype(nats+cnt) = t;
        trtmp(1,nats+cnt) = 1 ; trtmp(2,nats+cnt) = -1; trtmp(3,nats+cnt) = -1                  
        t = list(t);
      enddo

      t = head(nx*ny*nz);
      do while(t > 0)
        cnt = cnt + 1;
        buffer(1,nats+cnt) = coords(1,t)-Lx;
        buffer(2,nats+cnt) = coords(2,t)-Ly;
        buffer(3,nats+cnt) = coords(3,t)-Lz;
        trtmp(1,nats+cnt) = -1 ; trtmp(2,nats+cnt) = -1; trtmp(3,nats+cnt) = -1                  
        ntype(nats+cnt) = t;
        t = list(t);
      enddo

      Nskin = cnt;

      ! And now create a list do everything including the Skin/buffer layer
      dLx = Lx/real(nx,dp); dLy = Ly/real(ny,dp); dLz = Lz/real(nz,dp);
      buffer(1,:) = buffer(1,:)+dLx; 
      buffer(2,:) = buffer(2,:)+dLy; 
      buffer(3,:) = buffer(3,:)+dLz; ! Shift to avoid negative coordinates
      Lx = Lx + 2*dLx; 
      Ly = Ly + 2*dLy; 
      Lz = Lz + 2*dLz;
      nx = nx+2; 
      ny = ny+2; 
      nz = nz+2;

      deallocate(head) 
      deallocate(list)
      allocate(head(nx*ny*nz))
      allocate(list(nats+Nskin))

      do i = 1,nats+Nskin
        cell = 1 + floor(nx*buffer(1,i)/Lx) + floor(ny*buffer(2,i)/Ly)*nx + floor(nz*buffer(3,i)/Lz)*nx*ny;
        list(i) = head(cell);
        head(cell) = i;
      enddo

      !$omp parallel do default(none) private(i) &
      !$omp private(cnt,j,k,l,Tx,Ty,Tz) &
      !$omp private(cell,dist,t) &      
      !$omp shared(nats,dLx,dLy,dLz,buffer,nx,ny,nz,Lx,Ly,Lz)&
      !$omp shared(head,Rcut,trtmp,ntype,list,nl)
      do i = 1,nats
        cnt = 0;
        do j = -1,1
          do k = -1,1
            do l = -1,1
              Tx = buffer(1,i)+j*dLx;
              Ty = buffer(2,i)+k*dLy;
              Tz = buffer(3,i)+l*dLz;
              cell = 1 + floor(nx*Tx/Lx) + floor(ny*Ty/Ly)*nx + floor(nz*Tz/Lz)*nx*ny;
              t = head(cell);
              do while(t > 0)

                dist = (buffer(1,i)-buffer(1,t))**2 + (buffer(2,i)-buffer(2,t))**2 &
                  + (buffer(3,i)-buffer(3,t))**2

                dist = sqrt(dist)                                 

                if (dist < Rcut)then
                  cnt = cnt + 1;
                  nl%nnIx(cnt,i) = trtmp(1,t)
                  nl%nnIy(cnt,i) = trtmp(2,t)
                  nl%nnIz(cnt,i) = trtmp(3,t)
                  nl%nnType(cnt,i) = ntype(t);
                  nl%nnStruct(cnt,i) = ntype(t);
                  nl%nnStructMindist(cnt,i) = dist

                  !                   distvec(ntype(t)) = min(distvec(ntype(t)),dist)

                endif

                t = list(t);
              enddo
            enddo
          enddo          
        enddo

        nl%nrnnlist(i) = cnt;
        nl%nrnnStruct(i) = cnt;

        !         do ss = 1,nats
        !           if (tmp(ss) .gt. 0)then 
        !             cnt2 = cnt2 + 1
        !             nl%nnStruct(cnt2,i) = ss
        !             nl%nnStructMindist(cnt2,i) = distvec(ss)
        !           endif
        !         enddo
        !         nl%nrnnStruct(i) = cnt2

      enddo
      !$omp end parallel do

      deallocate(ntype)
      deallocate(head)
      deallocate(list)    
      deallocate(buffer) !Allocate max buffer atoms
      deallocate(trtmp)

    endif



  end subroutine build_nlist_int


  !> The neigbor list construction.
  !! \param coords system coordinates for which neighbor list should be constructed.
  !! \param lattice_vectors lattice vectors of the system.
  !! \param rcut coulomb cut off radius.
  !! \param nl neighbor list type. 
  !! \param verbose verbosity level.
  !!
  !! WARNING: This list only works for orthogonal lattice_vectors
  !! \todo Generalize neighbor list construction to nonorthogonal lattice vectors.
  !! 
  !!
  subroutine build_nlist(coords,lattice_vectors,rcut,nl,verbose)
    implicit none

    integer                              ::  cell, cnt, cnt2, i
    integer                              ::  j, k, l, m
    integer                              ::  nats, natspblock, nx, ny
    integer                              ::  nz, Nskin, ss, t
    integer, intent(in)                  ::  verbose 
    integer, allocatable                 ::  ntype(:), tmp(:)
    real(dp)                             ::  Lx, Ly, Lz, Tx
    real(dp)                             ::  dLx, dLy, dLz
    real(dp)                             ::  Ty, Tz, calpha, coulcut
    real(dp)                             ::  coulvol, dist, pi, sqrtx
    real(dp), allocatable                ::  head(:), list(:), buffer(:,:), distvec(:)
    real(dp), intent(in)                 ::  coords(:,:), lattice_vectors(:,:), rcut
    type(neighlist_type), intent(inout)  ::  nl

    nats = size(coords,dim=2)
    pi = 3.14159265358979323846264338327950_dp

    Lx = lattice_vectors(1,1)
    Ly = lattice_vectors(2,2)
    Lz = lattice_vectors(3,3)

    ! Division into # cell boxes: nx, ny, nz
    nx = floor(Lx/Rcut)
    ny = floor(Ly/Rcut)
    nz = floor(Lz/Rcut)

    ! Asuming an upper bound average density of 1 atom/ang^3 
    natspblock = floor(1.0d0*Rcut**3)

    if(verbose >= 0)  write(*,*) "In build_nlist ..."
    if(verbose >= 1)  write(*,*) "Number of atoms per block = ",natspblock
    if(verbose >= 1)  write(*,*)  "min(nx,ny,nz) =", min(nx,ny,nz)

    if(.not.(allocated(nl%nrnnlist)))then
      allocate(nl%nndist(natspblock,nats))
      allocate(nl%nnRx(natspblock,nats));
      allocate(nl%nnRy(natspblock,nats));
      allocate(nl%nnRz(natspblock,nats));
      allocate(nl%nnType(natspblock,nats))
      allocate(nl%nnStruct(natspblock,nats))
      allocate(nl%nrnnStruct(nats))
      allocate(nl%nrnnlist(nats))
      allocate(nl%nnStructMindist(natspblock,nats))
    endif

    allocate(tmp(nats));    

    if(min(nx,ny,nz) < 3)then   ! Brute force for small systems!

      if(verbose.GT.1)  write(*,*) "Performing brute force for small system ..."

      allocate(distvec(nats))
      distvec=1.0d5

      !$omp parallel do default(none) private(i) &
      !$omp private(cnt,cnt2,tmp,distvec,m,j,k,l,ss,Tx,Ty,Tz,dist) &
      !$omp shared(nats,coords,Lx,Ly,Lz,Rcut,nl)
      do i = 1,nats
        cnt = 0
        tmp = 0
        distvec = 1.0d5
        do m = 1,nats
          do j = -1,1
            do k = -1,1
              do l = -1,1
                Tx = coords(1,m)+j*Lx;
                Ty = coords(2,m)+k*Ly;
                Tz = coords(3,m)+l*Lz;
                dist = (coords(1,i)-Tx)**2 + (coords(2,i)-Ty)**2 + (coords(3,i)-Tz)**2
                dist = sqrt(dist)
                if (dist .lt. Rcut .and. dist .gt. 1d-12) then 
                  cnt = cnt + 1
                  nl%nndist(cnt,i) = dist                  
                  nl%nnRx(cnt,i) = Tx
                  nl%nnRy(cnt,i) = Ty
                  nl%nnRz(cnt,i) = Tz
                  nl%nnType(cnt,i) = m
                  tmp(m) = m
                  distvec(m) = min(distvec(m),dist)
                  ! if(i.eq.5)then 
                  !   write(*,*)i,m,dist,cnt
                  ! endif
                endif
              enddo
            enddo
          enddo

        enddo

        nl%nrnnlist(i) = cnt;
        cnt2 = 0;

        do ss = 1,nats
          if (tmp(ss) .gt. 0)then 
            cnt2 = cnt2 + 1
            nl%nnStruct(cnt2,i) = ss
            nl%nnStructMindist(cnt2,i) = distvec(ss)
          endif
        enddo
        nl%nrnnStruct(i) = cnt2
      enddo
      !$omp end parallel do

      deallocate(distvec)      

    else ! Do the same but now with linked lists in O(N)

      allocate(head(nx*ny*nz));    
      allocate(list(nats));        
      allocate(ntype(10*nats));
      allocate(buffer(3,nats*26)) !Allocate max buffer atoms

      buffer(:,1:nats) = coords

      head = 0 
      list = 0 
      do i = 1,nats  
        cell = 1 + floor(nx*coords(1,i)/Lx) + floor(ny*coords(2,i)/Ly)*nx &
          + floor(nz*coords(3,i)/Lz)*nx*ny;
        list(i) = head(cell);
        head(cell) = i;
        ntype(i) = i;        
      enddo

      !And now add a skin or surface buffer to account for periodic BC, all 26 of them!
      cnt = 0;
      do i = 1,nx*ny  !All boxes in the first (z=0) layer
        t = head(i);
        do while (t > 0) ! and all atoms of this first layer
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t);
          buffer(2,nats+cnt) = coords(2,t);
          buffer(3,nats+cnt) = coords(3,t)+Lz;
          ntype(nats+cnt) = t;
          t = list(t);          
        enddo
      enddo

      do i = 1,nx*ny*nz,nx
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t)+Lx;
          buffer(2,nats+cnt) = coords(2,t);
          buffer(3,nats+cnt) = coords(3,t);
          ntype(nats+cnt) = t;
          t = list(t);
        enddo
      enddo

      do i = 1,nx*ny*nz,nx*ny
        do k = i,i+nx-1
          t = head(k);
          do while(t > 0)
            cnt = cnt + 1;
            buffer(1,nats+cnt) = coords(1,t);
            buffer(2,nats+cnt) = coords(2,t)+Ly;
            buffer(3,nats+cnt) = coords(3,t);
            ntype(nats+cnt) = t;
            t = list(t);
          enddo
        enddo
      enddo

      do i = nx*ny*(nz-1)+1,nx*ny*nz
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t);
          buffer(2,nats+cnt) = coords(2,t);
          buffer(3,nats+cnt) = coords(3,t)-Lz;
          ntype(nats+cnt) = t;
          t = list(t);
        enddo
      enddo

      do i = nx,nx*ny*nz,nx
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t)-Lx;
          buffer(2,nats+cnt) = coords(2,t);
          buffer(3,nats+cnt) = coords(3,t);
          ntype(nats+cnt) = t;
          t = list(t);
        enddo
      enddo

      do i = nx*(ny-1)+1,nx*ny*nz,nx*ny
        do k = i,i+nx-1
          t = head(k);
          do while(t > 0)
            cnt = cnt + 1;
            buffer(1,nats+cnt) = coords(1,t);
            buffer(2,nats+cnt) = coords(2,t)-Ly;
            buffer(3,nats+cnt) = coords(3,t);
            ntype(nats+cnt) = t;
            t = list(t);
          enddo
        enddo
      enddo

      do i = 1,nx
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t);
          buffer(2,nats+cnt) = coords(2,t)+Ly;
          buffer(3,nats+cnt) = coords(3,t)+Lz;
          ntype(nats+cnt) = t;
          t = list(t);
        enddo
      enddo

      do i = nx*(ny-1)+1,nx*(ny-1)+nx
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t);
          buffer(2,nats+cnt) = coords(2,t)-Ly;
          buffer(3,nats+cnt) = coords(3,t)+Lz;
          ntype(nats+cnt) = t;
          t = list(t);
        enddo
      enddo

      do i = nx*ny*(nz-1)+1,nx*ny*(nz-1)+nx
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t);
          buffer(2,nats+cnt) = coords(2,t)+Ly;
          buffer(3,nats+cnt) = coords(3,t)-Lz;
          ntype(nats+cnt) = t;
          t = list(t);
        enddo
      enddo

      do i = nx*ny*(nz-1)+nx*(ny-1)+1,nx*ny*(nz-1)+nx*(ny-1)+nx
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t);
          buffer(2,nats+cnt) = coords(2,t)-Ly;
          buffer(3,nats+cnt) = coords(3,t)-Lz;
          ntype(nats+cnt) = t;
          t = list(t);
        enddo
      enddo

      do i = 1,nx*ny,nx
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t)+Lx;
          buffer(2,nats+cnt) = coords(2,t);
          buffer(3,nats+cnt) = coords(3,t)+Lz;
          ntype(nats+cnt) = t;
          t = list(t);
        enddo
      enddo

      do i = nx*ny*(nz-1)+1,nx*ny*nz,nx
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t)+Lx;
          buffer(2,nats+cnt) = coords(2,t);
          buffer(3,nats+cnt) = coords(3,t)-Lz;
          ntype(nats+cnt) = t;
          t = list(t);
        enddo
      enddo

      do i = nx*ny*(nz-1)+nx,nx*ny*nz,nx
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t)-Lx;
          buffer(2,nats+cnt) = coords(2,t);
          buffer(3,nats+cnt) = coords(3,t)-Lz;
          ntype(nats+cnt) = t;
          t = list(t);
        enddo
      enddo

      do i = nx,nx*ny,nx
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t)-Lx;
          buffer(2,nats+cnt) = coords(2,t);
          buffer(3,nats+cnt) = coords(3,t)+Lz;
          ntype(nats+cnt) = t;
          t = list(t);
        enddo
      enddo

      do i = 1,nx*ny*nz,nx*ny
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t)+Lx;
          buffer(2,nats+cnt) = coords(2,t)+Ly;
          buffer(3,nats+cnt) = coords(3,t);
          ntype(nats+cnt) = t;
          t = list(t);
        enddo
      enddo

      do i = nx*(ny-1)+1,nx*ny*nz,nx*ny
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t)+Lx;
          buffer(2,nats+cnt) = coords(2,t)-Ly;
          buffer(3,nats+cnt) = coords(3,t);
          ntype(nats+cnt) = t;
          t = list(t);
        enddo
      enddo

      do i = nx*(ny-1)+nx,nx*ny*nz,nx*ny
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t)-Lx;
          buffer(2,nats+cnt) = coords(2,t)-Ly;
          buffer(3,nats+cnt) = coords(3,t);
          ntype(nats+cnt) = t;
          t = list(t);
        enddo
      enddo

      do i = nx,nx*ny*nz,nx*ny
        t = head(i);
        do while(t > 0)
          cnt = cnt + 1;
          buffer(1,nats+cnt) = coords(1,t)-Lx;
          buffer(2,nats+cnt) = coords(2,t)+Ly;
          buffer(3,nats+cnt) = coords(3,t);
          ntype(nats+cnt) = t;
          t = list(t);
        enddo
      enddo

      t = head(1);
      do while(t > 0)
        cnt = cnt + 1;
        buffer(1,nats+cnt) = coords(1,t)+Lx;
        buffer(2,nats+cnt) = coords(2,t)+Ly;
        buffer(3,nats+cnt) = coords(3,t)+Lz;
        ntype(nats+cnt) = t;
        t = list(t);
      enddo

      t = head(nx);
      do while(t > 0)
        cnt = cnt + 1;
        buffer(1,nats+cnt) = coords(1,t)-Lx;
        buffer(2,nats+cnt) = coords(2,t)+Ly;
        buffer(3,nats+cnt) = coords(3,t)+Lz;
        ntype(nats+cnt) = t;
        t = list(t);
      enddo

      t = head(nx*(ny-1)+1);
      do while(t > 0)
        cnt = cnt + 1;
        buffer(1,nats+cnt) = coords(1,t)+Lx;
        buffer(2,nats+cnt) = coords(2,t)-Ly;
        buffer(3,nats+cnt) = coords(3,t)+Lz;
        ntype(nats+cnt) = t;
        t = list(t);
      enddo

      t = head(nx*(ny-1)+nx);
      do while(t > 0)
        cnt = cnt + 1;
        buffer(1,nats+cnt) = coords(1,t)-Lx;
        buffer(2,nats+cnt) = coords(2,t)-Ly;
        buffer(3,nats+cnt) = coords(3,t)+Lz;
        ntype(nats+cnt) = t;
        t = list(t);
      enddo

      t = head(nx*ny*(nz-1)+1);
      do while(t > 0)
        cnt = cnt + 1;
        buffer(1,nats+cnt) = coords(1,t)+Lx;
        buffer(2,nats+cnt) = coords(2,t)+Ly;
        buffer(3,nats+cnt) = coords(3,t)-Lz;
        ntype(nats+cnt) = t;
        t = list(t);
      enddo

      t = head(nx*ny*(nz-1)+nx);
      do while(t > 0)
        cnt = cnt + 1;
        buffer(1,nats+cnt) = coords(1,t)-Lx;
        buffer(2,nats+cnt) = coords(2,t)+Ly;
        buffer(3,nats+cnt) = coords(3,t)-Lz;
        ntype(nats+cnt) = t;
        t = list(t);
      enddo

      t = head(nx*ny*(nz-1)+nx*(ny-1)+1);
      do while(t > 0)
        cnt = cnt + 1;
        buffer(1,nats+cnt) = coords(1,t)+Lx;
        buffer(2,nats+cnt) = coords(2,t)-Ly;
        buffer(3,nats+cnt) = coords(3,t)-Lz;
        ntype(nats+cnt) = t;
        t = list(t);
      enddo

      t = head(nx*ny*nz);
      do while(t > 0)
        cnt = cnt + 1;
        buffer(1,nats+cnt) = coords(1,t)-Lx;
        buffer(2,nats+cnt) = coords(2,t)-Ly;
        buffer(3,nats+cnt) = coords(3,t)-Lz;
        ntype(nats+cnt) = t;
        t = list(t);
      enddo

      Nskin = cnt;

      ! And now create a list do everything including the Skin/buffer layer
      dLx = Lx/real(nx,dp); dLy = Ly/real(ny,dp); dLz = Lz/real(nz,dp);
      buffer(1,:) = buffer(1,:)+dLx; 
      buffer(2,:) = buffer(2,:)+dLy; 
      buffer(3,:) = buffer(3,:)+dLz; ! Shift to avoid negative coordinates
      Lx = Lx + 2*dLx; 
      Ly = Ly + 2*dLy; 
      Lz = Lz + 2*dLz;
      nx = nx+2; 
      ny = ny+2; 
      nz = nz+2;

      deallocate(head) 
      deallocate(list)
      allocate(head(nx*ny*nz))
      allocate(list(nats+Nskin))

      do i = 1,nats+Nskin
        cell = 1 + floor(nx*buffer(1,i)/Lx) + floor(ny*buffer(2,i)/Ly)*nx + floor(nz*buffer(3,i)/Lz)*nx*ny;
        list(i) = head(cell);
        head(cell) = i;
      enddo

      do i = 1,nats
        cnt = 0;
        cnt2 = 0;
        !         distvec = 10d5
        do j = -1,1
          do k = -1,1
            do l = -1,1
              Tx = buffer(1,i)+j*dLx;
              Ty = buffer(2,i)+k*dLy;
              Tz = buffer(3,i)+l*dLz;
              cell = 1 + floor(nx*Tx/Lx) + floor(ny*Ty/Ly)*nx + floor(nz*Tz/Lz)*nx*ny;
              t = head(cell);
              do while(t > 0)

                dist = (buffer(1,i)-buffer(1,t))**2 + (buffer(2,i)-buffer(2,t))**2 &
                  + (buffer(3,i)-buffer(3,t))**2

                dist = sqrt(dist)                                 

                if (dist < Rcut)then
                  cnt = cnt + 1;
                  nl%nndist(cnt,i) = dist;
                  nl%nnRx(cnt,i) = buffer(1,t)-dLx;
                  nl%nnRy(cnt,i) = buffer(2,t)-dLy;
                  nl%nnRz(cnt,i) = buffer(3,t)-dLz;
                  nl%nnType(cnt,i) = ntype(t);
                  !                  nl%nnStruct(cnt,i) = ntype(t);                 
                  !                    nl%nnStruct(cnt,i) = ntype(t);
                  !                    nl%nnStructMindist(cnt,i) = dist

                  !                   if (t .le. nats) then
                  cnt2 = cnt2 + 1;
                  nl%nnStruct(cnt2,i) = ntype(t);
                  !                     distvec(ntype(t)) = min(dist,distvec(ntype(t)))                
                  nl%nnStructMindist(cnt2,i) = dist
                  !                   endif
                endif

                t = list(t);
              enddo
            enddo
          enddo

          !           cnt = 0;
          !           do j=1,nats
          !             if(distmat(i,j) .ne.0) then 
          !              cnt = cnt + 1            
          !              nl%nnStructMindist(cnt,i) = distmat(i,j)                                     
          !           enddo
          !           
        enddo

        nl%nrnnlist(i) = cnt;
        nl%nrnnStruct(i) = cnt2;
        !         nl%nrnnStruct(i) = cnt;
      enddo

      deallocate(ntype)
      deallocate(head)
      deallocate(list)    
      deallocate(buffer) !Allocate max buffer atoms

    endif

    deallocate(tmp)

  end subroutine build_nlist

end module neighborlist_latte_mod



