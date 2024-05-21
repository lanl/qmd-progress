!> Hamiltonian module.
!! \ingroup LATTE
!! \brief Routines in this module are used to build the Hamiltonian and Svelap matrix from the system type.
!!
module ham_latte_mod

  use bml
  use tbparams_latte_mod
  use prg_extras_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  public :: get_hindex, get_hindex_coreHalo, get_hsmat, get_hsmat_vect, get_SKBlock, get_SKBlock_vect

contains

  !> Gets the Hamiltonian indices for every atom in the system.
  !! \param spindex Species index for every atom in the system.
  !! \param norbi Number of orbital for every species in the species list.
  !! \param hindex Start and end index for every atom in the system.
  !! \param norb Output parameter corresponding to the number of orbitals of the system.
  !! \note norb is always greater or equal nats.
  !! \note spindex can be gather using the system type:
  !! \verbatim spindex(:) = system%spindex(:) \endverbatim
  !! \note norbi can be gather from the tbparams type as it is a property
  !! that depends strictly on the parametrization;
  !! \verbatim norbi(:) = tbparams%norbi(:) \endverbatim
  !! \todo add verbosity
  !!
  subroutine get_hindex(spindex,norbi,hindex,norb,verbose)
    implicit none
    integer, intent(in) :: spindex(:)
    integer, intent(in) :: norbi(:)
    integer, optional, intent(in) :: verbose
    integer, intent(inout) :: norb
    integer :: nats, cnt, i
    integer, allocatable, intent(inout) :: hindex(:,:)


    nats = size(spindex,dim=1)
    cnt = 1

    if(present(verbose)) then
      if(verbose >= 1) then
        write(*,*)""; write(*,*)"In get hindex ..."
      endif
    endif

    if(.not.allocated(hindex))then
      allocate(hindex(2,nats))
    endif

    cnt = 1
    do i = 1,nats
      hindex(1,i) = cnt
      cnt = cnt + norbi(spindex(i))
      hindex(2,i) = cnt - 1
    enddo

    norb = cnt-1;

    if(present(verbose)) then
      if(verbose >= 1) then
        write(*,*)"Number of orbitals =",norb
      endif
    endif

  end subroutine get_hindex


  !> Gets the Hamiltonian indices for every atom in the system.
  !! \param spindex Species index for every atom in the system.
  !! \param norbi Number of orbital for every species in the species list.
  !! \param hindex Start and end index for every atom in the system.
  !! \param norb Output parameter corresponding to the number of orbitals of the
  !system.
  !! \note norb is always greater or equal nats.
  !! \note spindex can be gather using the system type:
  !! \verbatim spindex(:) = system%spindex(:) \endverbatim
  !! \note norbi can be gather from the tbparams type as it is a property
  !! that depends strictly on the parametrization;
  !! \verbatim norbi(:) = tbparams%norbi(:) \endverbatim
  !! \todo add verbosity
  !!
  subroutine get_hindex_coreHalo(spindex,natsCore,norbi,hindex,norb,norbsCore,verbose)
    implicit none
    integer, intent(in) :: spindex(:)
    integer, intent(in) :: norbi(:), natsCore
    integer, optional, intent(in) :: verbose
    integer, intent(inout) :: norb, norbsCore
    integer :: nats, cnt, i
    integer, allocatable, intent(inout) :: hindex(:,:)


    nats = size(spindex,dim=1)
    cnt = 1

    if(present(verbose))then
      if(verbose >= 1)then
        write(*,*)""; write(*,*)"In get hindex ..."
      endif
    endif

    if(.not.allocated(hindex))then
      allocate(hindex(2,nats))
    endif

    cnt = 1
    do i = 1,nats
      hindex(1,i) = cnt
      cnt = cnt + norbi(spindex(i))
      hindex(2,i) = cnt - 1
      if(i == natsCore)norbsCore = cnt - 1
    enddo

    norb = cnt-1;

    if(present(verbose).and.verbose >= 1)then
      write(*,*)"Number of orbitals in the subsystem =",norb
      write(*,*)"Number of orbitals in the core =",norbsCore
    endif

  end subroutine get_hindex_coreHalo


  !> Constructs Hamiltonian and Overlap Matrix
  !! \brief Construction of the Hamiltonian and Overlap matrix.
  !! \param ham_bml Hamiltonian in bml format.
  !! \param over_bml Overlap in bml format.
  !! \param coordinate Coordinates of the system. This can be getter from the system type:
  !! \verbatim system%coordinate \endverbatim
  !! \param lattice_vector Lattice vectors for the system.
  !! \param spindex Species indices (see system_type).
  !! \param norbi Number of orbital for every species.
  !! \param hindex Contains the Hamiltonian indices for every atom (see get_hindex)
  !! \param onsitesH Onsites energies for every pair of equal type. Allocation:
  !! \verbatim onsitesH(maxints,nsp) \endverbatim
  !! \param onsitesS Same as the Hamiltonian but for the overlap matrix.
  !! In this case the "onsite" elements are equal to 1.0. This was done to maintain
  !! a consistency and have the possibility of generalize the SK block construction.
  !! elements are or the same type.
  !! \param intPairsH,intPairsS See in intPairs_type
  !! \param threshold Threshold value for matrix elements.
  subroutine get_hsmat(ham_bml,over_bml,coordinate,lattice_vector,spindex,&
       norbi,hindex,onsitesH,onsitesS,intPairsH,intPairsS,threshold)
    implicit none
    character(20)                       ::  bml_type, bml_dmode
    integer                             ::  dimi, dimj, i, ii
    integer                             ::  j, jj, nats, norb, mdim
    integer, intent(in)                 ::  hindex(:,:), norbi(:), spindex(:)
    integer                             ::  maxnorbi
    real(dp)                            ::  ra(3), rb(3)
    real(dp), allocatable               ::  block(:,:,:), ham(:,:), over(:,:)
    real(dp), intent(in)                ::  coordinate(:,:), lattice_vector(:,:), onsitesH(:,:), onsitesS(:,:)
    real(dp), intent(in)                ::  threshold
    type(bml_matrix_t), intent(inout)   ::  ham_bml, over_bml
    type(intpairs_type), intent(in)     ::  intPairsH(:,:), intPairsS(:,:)

    nats = size(spindex,dim=1)
    norb = bml_get_N(ham_bml)
    mdim = bml_get_M(ham_bml)
    bml_type = bml_get_type(ham_bml)
    bml_dmode = bml_get_distribution_mode(ham_bml)

    !     if(bml_get_N(ham_bml).LE.0) then  !Carefull we need to clean S and H before rebuilding them!!!
    call bml_deallocate(ham_bml)
    call bml_deallocate(over_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,mdim,ham_bml, &
         bml_dmode)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,mdim,over_bml, &
         bml_dmode)
    !     endif

    maxnorbi = maxval(norbi)

    if(.not.allocated(ham)) then
      allocate(ham(norb,norb))
    endif
    if(.not.allocated(over)) then
      allocate(over(norb,norb))
    endif

    if (.not.allocated(block)) then
      allocate(block(maxnorbi,maxnorbi,nats))
    endif

    !$omp parallel do default(none) &
    !$omp private(ra,rb,dimi,dimj,ii,jj,j) &
    !$omp shared(nats,coordinate,hindex,spindex, intPairsS,intPairsH,threshold,lattice_vector,norbi,onsitesH,onsitesS,ham_bml,over_bml) &
    !$omp shared(block,ham,over)
    do i = 1, nats
      ra(:) = coordinate(:,i)
      dimi = hindex(2,i)-hindex(1,i)+1
      do j = 1, nats
        rb(:) = coordinate(:,j)
        dimj = hindex(2,j)-hindex(1,j)+1
        !Hamiltonian block for a-b atom pair
        call get_SKBlock(spindex(i),spindex(j),coordinate(:,i),&
             coordinate(:,j),lattice_vector,norbi,&
             onsitesH,intPairsH(spindex(i),spindex(j))%intParams,intPairsH(spindex(j),spindex(i))%intParams,block,i)

        do jj=1,dimj
          do ii=1,dimi
            ham(hindex(1,i)-1+ii,hindex(1,j)-1+jj) = block(ii,jj,i)
          enddo
        enddo

        call get_SKBlock(spindex(i),spindex(j),coordinate(:,i),&
             coordinate(:,j),lattice_vector,norbi,&
             onsitesS,intPairsS(spindex(i),spindex(j))%intParams,intPairsS(spindex(j),spindex(i))%intParams,block,i)

        do jj=1,dimj
          do ii=1,dimi
            over(hindex(1,i)-1+ii,hindex(1,j)-1+jj) = block(ii,jj,i)
          enddo
        enddo

        !  write(*,*)spindex(i),spindex(j)
        !  write(*,'(100F10.5)')intPairsS(spindex(i),spindex(j))%intParams
        ! call prg_print_matrix("block",block(:,:,i),1,4,1,4)

      enddo
    enddo

    !$omp end parallel do

    bml_type=bml_get_type(ham_bml) !Get the bml type
    call bml_import_from_dense(bml_type,over,over_bml,threshold,norb) !Dense to dense_bml
    call bml_import_from_dense(bml_type,ham,ham_bml,threshold,norb) !Dense to dense_bml

    if(allocated(ham)) then
      deallocate(ham)
    endif
    if(allocated(over)) then
      deallocate(over)
    endif

    !     call bml_print_matrix("ham_bml",ham_bml,0,6,0,6)
    !     call bml_print_matrix("over_bml",over_bml,0,6,0,6)

  end subroutine get_hsmat

    !> Constructs Hamiltonian and Overlap Matrix
  !! \brief Construction of the Hamiltonian and Overlap matrix.
  !! \param ham_bml Hamiltonian in bml format.
  !! \param over_bml Overlap in bml format.
  !! \param coordinate Coordinates of the system. This can be getter from the system type:
  !! \verbatim system%coordinate \endverbatim
  !! \param lattice_vector Lattice vectors for the system.
  !! \param spindex Species indices (see system_type).
  !! \param norbi Number of orbital for every species.
  !! \param hindex Contains the Hamiltonian indices for every atom (see get_hindex)
  !! \param onsitesH Onsites energies for every pair of equal type. Allocation:
  !! \verbatim onsitesH(maxints,nsp) \endverbatim
  !! \param onsitesS Same as the Hamiltonian but for the overlap matrix.
  !! In this case the "onsite" elements are equal to 1.0. This was done to maintain
  !! a consistency and have the possibility of generalize the SK block construction.
  !! elements are or the same type.
  !! \param intPairsH,intPairsS See in intPairs_type
  !! \param threshold Threshold value for matrix elements.
  subroutine get_hsmat_vect(ham_bml,over_bml,coordinate,lattice_vector,spindex,&
       norbi,hindex,onsitesH,onsitesS,intPairsH,intPairsS,threshold)
    implicit none
    character(20)                       ::  bml_type, bml_dmode
    integer                             ::  dimi, dimj, i, ii
    integer                             ::  j, jj, nats, norb, mdim
    integer, intent(in)                 ::  hindex(:,:), norbi(:), spindex(:)
    integer                             ::  maxnorbi
    real(dp)                            ::  ra(3), rb(3)
    real(dp), allocatable               ::  block(:,:,:), ham(:,:), over(:,:), ham_vect(:,:), over_vect(:,:)
    real(dp), intent(in)                ::  coordinate(:,:), lattice_vector(:,:), onsitesH(:,:), onsitesS(:,:)
    real(dp), intent(in)                ::  threshold
    type(bml_matrix_t), intent(inout)   ::  ham_bml, over_bml
    type(intpairs_type), intent(in)     ::  intPairsH(:,:), intPairsS(:,:)
    real(dp), allocatable               ::  intParams1(:,:,:), intParams2(:,:,:)
    integer, allocatable                ::  norbs_atidx(:)

    nats = size(spindex,dim=1)
    norb = bml_get_N(ham_bml)
    mdim = bml_get_M(ham_bml)
    bml_type = bml_get_type(ham_bml)
    bml_dmode = bml_get_distribution_mode(ham_bml)

    !     if(bml_get_N(ham_bml).LE.0) then  !Carefull we need to clean S and H before rebuilding them!!!
    call bml_deallocate(ham_bml)
    call bml_deallocate(over_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,mdim,ham_bml, &
         bml_dmode)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,mdim,over_bml, &
         bml_dmode)
    !     endif

    maxnorbi = maxval(norbi)

    if(.not.allocated(ham)) then
      allocate(ham(norb,norb))
      allocate(ham_vect(norb,norb))
    endif
    if(.not.allocated(over)) then
      allocate(over(norb,norb))
      allocate(over_vect(norb,norb))
    endif

    if (.not.allocated(block)) then
      allocate(block(maxnorbi,maxnorbi,nats))
    endif

    if (.not.allocated(norbs_atidx)) then
       allocate(norbs_atidx(nats))
    endif
    do i=1,nats
       norbs_atidx(i) = norbi(spindex(i))
    enddo
    ! !$omp parallel do default(none) &
    ! !$omp private(ra,rb,dimi,dimj,ii,jj,j) &
    ! !$omp shared(nats,coordinate,hindex,spindex, intPairsS,intPairsH,threshold,lattice_vector,norbi,onsitesH,onsitesS,ham_bml,over_bml) &
    ! !$omp shared(block,ham,over)
    ! do i = 1, nats
    !   ra(:) = coordinate(:,i)
    !   dimi = hindex(2,i)-hindex(1,i)+1
    !   do j = 1, nats
    !     rb(:) = coordinate(:,j)
    !     dimj = hindex(2,j)-hindex(1,j)+1
    !     !Hamiltonian block for a-b atom pair
    !     call get_SKBlock(spindex(i),spindex(j),coordinate(:,i),&
    !          coordinate(:,j),lattice_vector,norbi,&
    !          onsitesH,intPairsH(spindex(i),spindex(j))%intParams,intPairsH(spindex(j),spindex(i))%intParams,block,i)

    !     do jj=1,dimj
    !       do ii=1,dimi
    !         ham(hindex(1,i)-1+ii,hindex(1,j)-1+jj) = block(ii,jj,i)
    !       enddo
    !     enddo

    !     call get_SKBlock(spindex(i),spindex(j),coordinate(:,i),&
    !          coordinate(:,j),lattice_vector,norbi,&
    !          onsitesS,intPairsS(spindex(i),spindex(j))%intParams,intPairsS(spindex(j),spindex(i))%intParams,block,i)

    !     do jj=1,dimj
    !       do ii=1,dimi
    !         over(hindex(1,i)-1+ii,hindex(1,j)-1+jj) = block(ii,jj,i)
    !       enddo
    !     enddo

    !     !  write(*,*)spindex(i),spindex(j)
    !     !  write(*,'(100F10.5)')intPairsS(spindex(i),spindex(j))%intParams
    !     ! call prg_print_matrix("block",block(:,:,i),1,4,1,4)

    !   enddo
    ! enddo

    !!$omp end parallel do

    allocate(intParams1(nats,16,4))
    allocate(intParams2(nats,16,4))

    !$omp parallel do default(none) &
    !$omp private(ra,rb,dimi,dimj,ii,jj,j,intParams1,intParams2) &
    !$omp shared(nats,coordinate,hindex,spindex, intPairsS,intPairsH,threshold,lattice_vector,norbs_atidx,onsitesH,onsitesS,ham_bml,over_bml) &
    !$omp shared(ham_vect,over_vect)
    do i = 1, nats
       do j = 1,nats
          intParams1(j,:,:) = intPairsH(spindex(i),spindex(j))%intParams(:,1:4)
          intParams2(j,:,:) = intPairsH(spindex(j),spindex(i))%intParams(:,1:4)
       enddo
          
        call get_SKBlock_vect(spindex,coordinate(:,i),coordinate(:,:),lattice_vector,norbs_atidx,&
             onsitesH,intParams1(:,:,:),intParams2(:,:,:),ham_vect(hindex(1,i):hindex(2,i),:),i)

       do j = 1,nats
          intParams1(j,:,:) = intPairsS(spindex(i),spindex(j))%intParams(:,1:4)
          intParams2(j,:,:) = intPairsS(spindex(j),spindex(i))%intParams(:,1:4)
       enddo

       call get_SKBlock_vect(spindex,coordinate(:,i),coordinate(:,:),lattice_vector,norbs_atidx,&
            onsitesS,intParams1(:,:,:),intParams2(:,:,:),over_vect(hindex(1,i):hindex(2,i),:),i)
        
        !  write(*,*)spindex(i),spindex(j)
        !  write(*,'(100F10.5)')intPairsS(spindex(i),spindex(j))%intParams
        ! call prg_print_matrix("block",block(:,:,i),1,4,1,4)

    enddo

    !$omp end parallel do

    ! if(.not.all(abs(ham-ham_vect).lt.1.D-12))then
    !    do i = 1,norb
    !       do j = 1,norb
    !          if(abs(ham(i,j)-ham_vect(i,j)).ge.1.D-12)then
    !             write(*,*)"GET_SKBLOCK_VECT: Vectorized ham differs at i,j = ",i,j,"with difference ",ham(i,j)-ham_vect(i,j)
    !          endif
    !       enddo
    !    enddo
    ! endif
    ! if(.not.all(abs(over-over_vect).lt.1.D-12))then
    !    do i = 1,norb
    !       do j = 1,norb
    !          if(abs(over(i,j)-over_vect(i,j)).ge.1.D-12)then
    !             write(*,*)"GET_SKBLOCK_VECT: Vectorized over differs at i,j = ",i,j,"with difference ",over(i,j)-over_vect(i,j)
    !          endif
    !       enddo
    !    enddo
    ! endif
    ! write(*,*)"ham(1,:) = ",ham(1,:)
    ! write(*,*)"ham_vect(1,:) = ",ham_vect(1,:)
    ! write(*,*)"ham(2,:) = ",ham(2,:)
    ! write(*,*)"ham_vect(2,:) = ",ham_vect(2,:)
    bml_type=bml_get_type(ham_bml) !Get the bml type
    call bml_import_from_dense(bml_type,over_vect,over_bml,threshold,norb) !Dense to dense_bml
    call bml_import_from_dense(bml_type,ham_vect,ham_bml,threshold,norb) !Dense to dense_bml

    if(allocated(ham)) then
       deallocate(ham)
       deallocate(ham_vect)
    endif
    if(allocated(over)) then
      deallocate(over)
      deallocate(over_vect)
    endif
    deallocate(intParams1)
    deallocate(intParams2)

    !     call bml_print_matrix("ham_bml",ham_bml,0,6,0,6)
    !     call bml_print_matrix("over_bml",over_bml,0,6,0,6)

  end subroutine get_hsmat_vect
  !> Function to calculate the bond integral for a given distance
  !! and coefficients set.
  !! \param dr distance between atoms.
  !! \param f parameters (coefficients) for the bond integral.
  real(dp) function bondIntegral(dr,f)
    implicit none
    real(dp) :: rmod
    real(dp) :: polynom
    real(dp) :: rminusr1
    real(dp) :: x
    real(dp), intent(in) :: dr
    real(dp), intent(in) :: f(16)

    if(dr <= f(7))then
      rmod = dr - f(6);
      polynom = rmod*(f(2) + rmod*(f(3) + rmod*(f(4) + f(5)*rmod)));
      x = exp(polynom);
    elseif(dr > f(7).and.dr < f(8))then
      rminusr1 = dr - f(7)
      x = f(9) + rminusr1*(f(10) + rminusr1*(f(11) + rminusr1*(f(12) + rminusr1*(f(13) + rminusr1*f(14)))))
    else
      x = 0
    end if
    bondintegral = f(1)*x

  end function bondintegral

  !> Function to calculate the bond integral for a given distance
  !! and coefficients set.
  !! \param dr distance between atoms.
  !! \param f parameters (coefficients) for the bond integral.
  function bondIntegral_vect(dr,f) result(x)
    implicit none
    real(dp), allocatable :: rmod(:)
    real(dp), intent(in) :: dr(:)
    real(dp), intent(in) :: f(:,:)
    real(dp)             :: x(size(dr))
    logical, allocatable :: low_mask(:),mid_mask(:)
    integer              :: vect_size
    
    vect_size = size(dr)
    
    allocate(low_mask(vect_size))
    allocate(mid_mask(vect_size))
    allocate(rmod(vect_size))    

    mid_mask = dr.gt.f(:,7).and.dr.lt.f(:,8)
    low_mask = dr.le.f(:,7)
    x = 0.0D0
    
    ! Calculate low values
    rmod = dr - f(:,6)
    x = merge(f(:,1)*exp(rmod*(f(:,2) + rmod*(f(:,3) + rmod*(f(:,4) + f(:,5)*rmod)))),0.0D0,low_mask)

    ! Calculate mid values
    rmod = dr - f(:,7)
    x = merge(f(:,1)*(f(:,9) + rmod*(f(:,10) + rmod*(f(:,11) + rmod*(f(:,12) + rmod*(f(:,13) + rmod*f(:,14)))))),x,mid_mask)

    deallocate(low_mask)
    deallocate(mid_mask)
    deallocate(rmod)
    
  end function bondintegral_vect

  !> Standard Slater-Koster sp-parameterization for an atomic block between a pair of atoms
  !! \param sp1 Species index for atom 1. This can be obtained from the
  !! system type as following:
  !! \verbatim sp1 = system%spindex(atom1) \endverbatim
  !! \param sp2 Species index for atom 2.
  !! \param coorda Coordinates for atom 1.
  !! \param coordb Coordinates for atom 2.
  !! \param lattice_vectors Lattice vectors for the system. This can be obtained from
  !! the system type as following:
  !! \verbatim lattice_vectors = system%lattice_vectors \endverbatim
  !! \param norbi Number of orbitals for every species in the system. This can be obtained from
  !! the tbparams type as following:
  !! \verbatim norbi = tbparams%norbi \endverbatim
  !! \param onsites Onsites energies for every pair of equal type. Two different variants
  !! onsitesH and onsitesS will be used as inputs (see get_hsmat routine) Allocation:
  !! \verbatim onsites(maxints,nsp) \endverbatim
  !! \param intParams See intpairs_type.
  !! \param block Output parameter SK block.
  !! \param atnum Input atom number
  subroutine get_SKBlock(sp1,sp2,coorda,coordb,lattice_vectors&
       ,norbi,onsites,intParams,intParamsr,blk,atnum)
    implicit none
    integer                              ::  dimi, dimj, i, nr_shift_X
    integer                              ::  nr_shift_Y, nr_shift_Z
    integer, intent(in)                  ::  norbi(:), sp1, sp2, atnum
    real(dp)                             ::  HPPP, HPPS, HSPS, HSPSR, HSSS
    real(dp)                             ::  L, LBox(3), M, N
    real(dp)                             ::  PPSMPP, PXPX, PXPY, PXPZ
    real(dp)                             ::  PYPX, PYPY, PYPZ, PZPX
    real(dp)                             ::  PZPY, PZPZ, dr, ra(3)
    real(dp)                             ::  rab(3), rb(3), rxb, ryb
    real(dp)                             ::  rzb
    real(dp), allocatable, intent(inout)  ::  blk(:,:,:)
    real(dp), intent(in)                 ::  coorda(:), coordb(:), intParams(:,:), lattice_vectors(:,:)
    real(dp), intent(in)                 ::  onsites(:,:), intParamsr(:,:)

    ra = coorda
    rb = coordb

    dimi= norbi(sp1)
    dimj= norbi(sp2)

    !     write(*,*)atom_type_a, atom_type_b,dimi,dimj

    !!    if(allocated(block))then
    !!      deallocate(block)
    !!    endif

    !!    allocate(block(dimi,dimj))
    blk(:,:,atnum)=0.0_dp

      !     call write_matrix_to_screen("block",block,size(block,dim=1),size(block,dim=2))

      RXb = Rb(1); RYb = Rb(2); RZb = Rb(3)

      ! For cubic lattice
      !     LBox(1) = lattice_vectors(1,1)
      !     LBox(2) = lattice_vectors(2,2)
      !     LBox(3) = lattice_vectors(3,3)

      !Periodic BC shifts in X, Y and Z. Costs a lot extra!
      !do nr_shift_x = -1,1
      !  do nr_shift_y = -1,1
      !    do nr_shift_z = -1,1

            !rb(1) = RXb + nr_shift_x*lattice_vectors(1,1) ! shifts for pbc
            ! rb(1) = rb(1) + nr_shift_y*lattice_vectors(2,1) ! shifts for pbc
            ! rb(1) = rb(1) + nr_shift_z*lattice_vectors(3,1) ! shifts for pbc

            !rb(2) = RYb + nr_shift_y*lattice_vectors(2,2) ! shifts for pbc
            ! rb(2) = rb(2) + nr_shift_x*lattice_vectors(1,2) ! shifts for pbc
            ! rb(2) = rb(2) + nr_shift_z*lattice_vectors(3,2) ! shifts for pbc

            !rb(3) = RZb + nr_shift_z*lattice_vectors(3,3) ! shifts for pbc
            ! rb(3) = rb(3) + nr_shift_y*lattice_vectors(2,3) ! shifts for pbc
            ! rb(3) = rb(3) + nr_shift_x*lattice_vectors(1,3) ! shifts for pbc
      do i = 1,3
         Rab(i) = modulo((Rb(i) - Ra(i)) + 0.5_dp*lattice_vectors(i,i),lattice_vectors(i,i)) - 0.5_dp * lattice_vectors(i,i)
      enddo

            !Rab = Rb-Ra;  ! OBS b - a !!!
            !dR = sqrt(Rab(1)**2+ Rab(2)**2+ Rab(3)**2)
            dR = norm2(Rab)
            
            if(dR.lt.6.5_dp)then
               if(dR .LT.1e-12)then !same position and thus the same type sp1 = sp2
                  do i=1,dimi
                     blk(i,i,atnum) = onsites(i,sp1)
                  enddo
               else
                  
                  L = Rab(1)/dR;  !Direction cosines
                  M = Rab(2)/dR;
                  N = Rab(3)/dR;
                  if(dimi == dimj.and.dimi == 1)then        !s-s  overlap 1 x 1 block
                     HSSS = BondIntegral(dR,intParams(:,1))  !Calculate the s-s bond integral
                     blk(1,1,atnum) = blk(1,1,atnum) + HSSS
                  elseif(dimi < dimj.and.dimi == 1)then    !s-sp overlap 1 x 4 block
                     HSSS = BondIntegral(dR,intParams(:,1))
                     blk(1,1,atnum) = blk(1,1,atnum) + HSSS
                     HSPS = BondIntegral(dR,intParams(:,2))
                     blk(1,2,atnum) = blk(1,2,atnum) + L*HSPS
                     blk(1,3,atnum) = blk(1,3,atnum) + M*HSPS
                     blk(1,4,atnum) = blk(1,4,atnum) + N*HSPS
                  elseif(dimi > dimj.and.dimj == 1)then ! sp-s overlap 4 x 1 block
                     HSSS = BondIntegral(dR,intParams(:,1))
                     blk(1,1,atnum) = blk(1,1,atnum) + HSSS
                     HSPS = BondIntegral(dR,intParams(:,2))
                     blk(2,1,atnum) = blk(2,1,atnum) - L*HSPS
                     blk(3,1,atnum) = blk(3,1,atnum) - M*HSPS
                     blk(4,1,atnum) = blk(4,1,atnum) - N*HSPS
                  elseif(dimi == dimj.and.dimj == 4)then !sp-sp overlap
                     HSSS = BondIntegral(dR,intParams(:,1))
                     HSPS = BondIntegral(dR,intParams(:,2))
                     HSPSR = BondIntegral(dR,intParamsr(:,2))
                     HPPS = BondIntegral(dR,intParams(:,3))
                     HPPP = BondIntegral(dR,intParams(:,4))
                     PPSMPP = HPPS - HPPP
                     PXPX = HPPP + L*L*PPSMPP
                     PXPY = L*M*PPSMPP
                     PXPZ = L*N*PPSMPP
                     PYPX = M*L*PPSMPP
                     PYPY = HPPP + M*M*PPSMPP
                     PYPZ = M*N*PPSMPP
                     PZPX = N*L*PPSMPP
                     PZPY = N*M*PPSMPP
                     PZPZ = HPPP + N*N*PPSMPP
                     blk(1,1,atnum) = blk(1,1,atnum) + HSSS
                     blk(1,2,atnum) = blk(1,2,atnum) + L*HSPS
                     blk(1,3,atnum) = blk(1,3,atnum) + M*HSPS
                     blk(1,4,atnum) = blk(1,4,atnum) + N*HSPS
                     blk(2,1,atnum) = blk(2,1,atnum) - L*HSPSR  !Change spindex
                     blk(2,2,atnum) = blk(2,2,atnum) + PXPX
                     blk(2,3,atnum) = blk(2,3,atnum) + PXPY
                     blk(2,4,atnum) = blk(2,4,atnum) + PXPZ
                     blk(3,1,atnum) = blk(3,1,atnum) - M*HSPSR  !Change spindex
                     blk(3,2,atnum) = blk(3,2,atnum) + PYPX
                     blk(3,3,atnum) = blk(3,3,atnum) + PYPY
                     blk(3,4,atnum) = blk(3,4,atnum) + PYPZ
                     blk(4,1,atnum) = blk(4,1,atnum) - N*HSPSR  !Change spindex
                     blk(4,2,atnum) = blk(4,2,atnum) + PZPX
                     blk(4,3,atnum) = blk(4,3,atnum) + PZPY
                     blk(4,4,atnum) = blk(4,4,atnum) + PZPZ
                  endif
               endif
            endif
         !enddo
      !enddo
   !enddo
   !            write(*,*)"block",dr,block
   !            stop
 end subroutine get_SKBlock

  !> Standard Slater-Koster sp-parameterization for an atomic block between a pair of atoms
  !! \param sp1 Species index for atom 1. This can be obtained from the
  !! system type as following:
  !! \verbatim sp1 = system%spindex(atom1) \endverbatim
  !! \param sp2 Species index for atom 2.
  !! \param coorda Coordinates for atom 1.
  !! \param coordb Coordinates for atom 2.
  !! \param lattice_vectors Lattice vectors for the system. This can be obtained from
  !! the system type as following:
  !! \verbatim lattice_vectors = system%lattice_vectors \endverbatim
  !! \param norbi Number of orbitals for every species in the system. This can be obtained from
  !! the tbparams type as following:
  !! \verbatim norbi = tbparams%norbi \endverbatim
  !! \param onsites Onsites energies for every pair of equal type. Two different variants
  !! onsitesH and onsitesS will be used as inputs (see get_hsmat routine) Allocation:
  !! \verbatim onsites(maxints,nsp) \endverbatim
  !! \param intParams See intpairs_type.
  !! \param block Output parameter SK block.
  !! \param atnum Input atom number
  subroutine get_SKBlock_vect(sp,refcoord,coord,lattice_vectors&
       ,norbs,onsites,intParams1,intParams2,blk,atnum)
    implicit none
    integer                              ::  dimi, dimj, i, nr_shift_X
    integer                              ::  nr_shift_Y, nr_shift_Z, nats, norbsall, mask_size, iorb
    integer, intent(in)                  ::  norbs(:), sp(:), atnum
    real(dp), allocatable                ::  HPPP(:), HPPS(:), HSPS(:), HSSS(:), PPSMPP(:)
    real(dp), allocatable                ::  L(:), M(:), N(:)
    real(dp), allocatable                ::  dr(:), dr_m(:), rab(:,:)
    real(dp), allocatable                ::  dx_m(:), dy_m(:), dz_m(:), blk_m(:,:,:), onsites_m(:)
    real(dp), intent(inout) ::  blk(:,:)
    real(dp), intent(in)                 ::  refcoord(:),coord(:,:), lattice_vectors(:,:)
    real(dp), intent(in)                 ::  onsites(:,:)
    real(dp), intent(in)                 ::  intParams1(:,:,:),intParams2(:,:,:)
    logical, allocatable                 ::  dist_mask(:), onsite_mask(:), calc_mask(:), calcs_mask(:), calcsp_mask(:), param_mask(:,:), calc_mask_for_porbs(:)
    logical, allocatable                 ::  sorb_mask(:), pxorb_mask(:), pyorb_mask(:), pzorb_mask(:)
    logical, allocatable                 ::  sorb_at_mask(:), sporb_at_mask(:)
    integer, allocatable                 ::  atomidx(:), atomidx_m(:), orbidx(:), orbidx_m(:), orbidx_sel(:)
    real(dp), allocatable                ::  intParams(:,:)

    nats = size(coord,dim=2)
    !write(*,*)"GET_SKBLOCK_VECT: nats = ",nats,"when atnum = ",atnum
    norbsall = sum(norbs)
    
    blk(:,:)=0.0_dp

    if(allocated(dr))then
       if(nats.ne.size(dr,dim=1))then
          deallocate(dr)
          deallocate(atomidx)
          deallocate(orbidx)
          deallocate(rab)
          deallocate(dist_mask)
          deallocate(onsite_mask)
          deallocate(calc_mask)
          deallocate(calcs_mask)
          deallocate(calcsp_mask)
          deallocate(param_mask)
          deallocate(sorb_mask)
          deallocate(pxorb_mask)
          deallocate(pyorb_mask)
          deallocate(pzorb_mask)
          deallocate(sorb_at_mask)
          deallocate(sporb_at_mask)
       endif
    endif
    if(.not.allocated(dr))then
       allocate(dr(nats))
       allocate(atomidx(nats))
       atomidx = (/(i,i=1,nats)/)
       allocate(orbidx(norbsall))
       orbidx = (/(i,i=1,norbsall)/)
       allocate(rab(nats,3))
       allocate(dist_mask(nats))
       allocate(onsite_mask(nats))
       allocate(calc_mask(nats))
       allocate(calcs_mask(nats))
       allocate(calcsp_mask(nats))
       allocate(param_mask(nats,16))
       allocate(sorb_at_mask(nats))
       allocate(sporb_at_mask(nats))
       allocate(sorb_mask(norbsall))
       allocate(pxorb_mask(norbsall))
       allocate(pyorb_mask(norbsall))
       allocate(pzorb_mask(norbsall))
       sorb_at_mask = .false.
       sporb_at_mask = .false.
       onsite_mask = .false.
       sorb_mask = .false.
       pxorb_mask = .false.
       pyorb_mask = .false.
       pzorb_mask = .false.
       iorb = 1
       do i = 1,nats
          if(i.eq.atnum)then
             onsite_mask(i) = .true.
             blk(1,iorb) = onsites(1,sp(atnum))
             if(norbs(i).eq.4)then
                blk(2,iorb+1) = onsites(2,sp(atnum))
                blk(3,iorb+2) = onsites(3,sp(atnum))
                blk(4,iorb+3) = onsites(4,sp(atnum))
             endif
          endif
          sorb_mask(iorb) = .true.
          iorb = iorb + 1
          if(norbs(i).eq.1)then
             sorb_at_mask(i) = .true.
          elseif(norbs(i).eq.4)then
             sporb_at_mask(i) = .true.
             pxorb_mask(iorb) = .true.
             pyorb_mask(iorb+1) = .true.
             pzorb_mask(iorb+2) = .true.
             iorb = iorb + 3
          else
             write(*,*)"GET_SKBLOCK_VECT: number of orbitals not 1 or 4 for atom ",i,": ABORT"
             stop
          endif
       enddo
    endif
    
    do i = 1,3
       Rab(:,i) = coord(i,:)
       Rab(:,i) = modulo((Rab(:,i) - refcoord(i) + 0.5_dp*lattice_vectors(i,i)),lattice_vectors(i,i)) - 0.5_dp * lattice_vectors(i,i)
    enddo

    dR(:) = norm2(Rab(:,:),dim=2)

    dist_mask = dr.lt.6.5_dp
    calc_mask = dist_mask.and..not.onsite_mask
    calcsp_mask = calc_mask.and.sporb_at_mask
    calcs_mask = calc_mask.and.sorb_at_mask
    ! First mask using the s orbitals
    ! The s-s term for all atoms is the same in this case
    ! The sorb_mask has the same number of .true. elements as the calc_mask
    orbidx_m = pack(orbidx,mask=sorb_mask)
    orbidx_sel = pack(orbidx_m,mask=calc_mask)
    atomidx_m = pack(atomidx,mask=calc_mask)
    mask_size = size(atomidx_m)
    ! dr_m = pack(dr,calc_mask)
    ! dx_m = pack(rab(:,1),calc_mask)
    ! dy_m = pack(rab(:,2),calc_mask)
    ! dz_m = pack(rab(:,3),calc_mask)
    allocate(intParams(mask_size,16))
    intParams(:,:) = intParams1(atomidx_m(:),:,1)
    blk(1,orbidx_sel) = BondIntegral_vect(dr(atomidx_m(:)),intParams)
    deallocate(intParams)
    deallocate(atomidx_m)
    deallocate(orbidx_sel)
    !deallocate(dr_m)
    !deallocate(dx_m)
    !deallocate(dy_m)
    !deallocate(dz_m)
    ! If the atnum atom is a sp atom, calculate the p*-s elements for
    !   s orbital atoms and sp orbital atoms separately. Do the s orbital
    !   atoms first here.
    if(norbs(atnum).eq.4)then
       ! First the s orbital atoms
       atomidx_m = pack(atomidx,mask=calcs_mask)
       mask_size = size(atomidx_m)
       allocate(intParams(mask_size,16))
       allocate(HSPS(mask_size))
       intParams(:,:) = intParams1(atomidx_m(:),:,2)
       !dr_m = pack(dr,calcs_mask)
       HSPS = BondIntegral_vect(dr(atomidx_m(:)),intParams)
       ! dx_m = pack(rab(:,1),calcs_mask)
       ! dy_m = pack(rab(:,2),calcs_mask)
       ! dz_m = pack(rab(:,3),calcs_mask)
       !L = rab(atomidx_m(:),1)/dr_m
       !M = rab(atomidx_m(:),2)/dr_m
       !N = rab(atomidx_m(:),3)/dr_m
       orbidx_sel = pack(orbidx_m,mask=calcs_mask) ! orbidx_m still is the s orbital mask
       !write(*,*)"GET_SKBLOCK_VECT: S atom S orbital calc mask for atom ",atnum,"= ",orbidx_sel
       if(size(orbidx_sel).ne.size(atomidx_m))then
          write(*,*)"GET_SKBLOCK_VECT: size mismatch of orbital and atom index arrays. Abort."
          stop
       endif
       blk(2,orbidx_sel) = - rab(atomidx_m(:),1)/dr(atomidx_m(:)) * HSPS
       blk(3,orbidx_sel) = - rab(atomidx_m(:),2)/dr(atomidx_m(:)) * HSPS
       blk(4,orbidx_sel) = - rab(atomidx_m(:),3)/dr(atomidx_m(:)) * HSPS
       deallocate(intParams)
       deallocate(HSPS)
       ! deallocate(dr_m)
       ! deallocate(dx_m)
       ! deallocate(dy_m)
       ! deallocate(dz_m)
       deallocate(orbidx_sel)
       deallocate(atomidx_m)
    endif
    ! Now work on the sp orbital atoms
   ! First create some common masked arrays
    atomidx_m = pack(atomidx,mask=calcsp_mask)
    mask_size = size(atomidx_m)
    dr_m = pack(dr,calcsp_mask)
    ! dx_m = pack(rab(:,1),calcsp_mask)
    ! dy_m = pack(rab(:,2),calcsp_mask)
    ! dz_m = pack(rab(:,3),calcsp_mask)
    allocate(L(mask_size))
    allocate(M(mask_size))
    allocate(N(mask_size))
    L = rab(atomidx_m(:),1)/dr_m
    M = rab(atomidx_m(:),2)/dr_m
    N = rab(atomidx_m(:),3)/dr_m
    ! Now work on the p*-s elements first, if relevant
    allocate(intParams(mask_size,16))
    allocate(HSPS(mask_size))
    if(norbs(atnum).eq.4)then
       intParams(:,:) = intParams2(atomidx_m(:),:,2)
       HSPS = BondIntegral_vect(dr_m,intParams)
       orbidx_sel = pack(orbidx_m,mask=calcsp_mask) ! orbidx_m still is the s orbital mask
       if(size(orbidx_sel).ne.size(atomidx_m))then
          write(*,*)"GET_SKBLOCK_VECT: size mismatch of orbital and sp atom index arrays. Abort."
          stop
       endif
       blk(2,orbidx_sel) = - L * HSPS
       blk(3,orbidx_sel) = - M * HSPS
       blk(4,orbidx_sel) = - N * HSPS
       deallocate(orbidx_sel)
    endif
    ! At this point all of the s orbital calculations are done. Move on to p orbitals.
    calc_mask_for_porbs = pack(calc_mask,sporb_at_mask)
    intParams(:,:) = intParams1(atomidx_m(:),:,2)
    HSPS = BondIntegral_vect(dr_m,intParams)
    if(norbs(atnum).eq.4)then ! Calculate extra params if the reference atom is sp type
       allocate(HPPP(mask_size))
       allocate(PPSMPP(mask_size))
       intParams(:,:) = intParams1(atomidx_m(:),:,3)
       PPSMPP = BondIntegral_vect(dr_m,intParams)
       intParams(:,:) = intParams1(atomidx_m(:),:,4)
       HPPP = BondIntegral_vect(dr_m,intParams)
       PPSMPP = PPSMPP - HPPP
    endif
    deallocate(orbidx_m)
    orbidx_m = pack(orbidx,mask=pxorb_mask) ! Select px orbitals
    orbidx_sel = pack(orbidx_m,calc_mask_for_porbs)
    blk(1,orbidx_sel) = + L(:) * HSPS(:)
    if(norbs(atnum).eq.4)then
       blk(2,orbidx_sel) = HPPP + L*L*PPSMPP
       blk(3,orbidx_sel) = M*L*PPSMPP
       blk(4,orbidx_sel) = N*L*PPSMPP
    endif
    ! deallocate(orbidx_m)
    ! deallocate(orbidx_sel)
    ! orbidx_m = pack(orbidx,mask=pyorb_mask) ! Select py orbitals
    ! orbidx_sel = pack(orbidx_m,calc_mask_for_porbs)
    blk(1,orbidx_sel+1) = + M(:) * HSPS(:)
    if(norbs(atnum).eq.4)then
       blk(2,orbidx_sel+1) = L*M*PPSMPP
       blk(3,orbidx_sel+1) = HPPP + M*M*PPSMPP
       blk(4,orbidx_sel+1) = N*M*PPSMPP
    endif
    ! deallocate(orbidx_m)
    ! deallocatE(orbidx_sel)
    ! orbidx_m = pack(orbidx,mask=pzorb_mask) ! Select pz orbitals
    ! orbidx_sel = pack(orbidx_m,calc_mask_for_porbs)
    blk(1,orbidx_sel+2) = + N(:) * HSPS(:)
    if(norbs(atnum).eq.4)then
       blk(2,orbidx_sel+2) = L*N*PPSMPP
       blk(3,orbidx_sel+2) = M*N*PPSMPP
       blk(4,orbidx_sel+2) = HPPP + N*N*PPSMPP
       deallocate(HPPP)
    endif
    deallocate(orbidx_m)
    deallocatE(orbidx_sel)
    deallocate(calc_mask_for_porbs)
    deallocate(HSPS)
    deallocate(intParams)
    deallocate(L)
    deallocate(M)
    deallocate(N)
! if(.false.)then
!  endif   
 end subroutine get_SKBlock_vect

end module ham_latte_mod
