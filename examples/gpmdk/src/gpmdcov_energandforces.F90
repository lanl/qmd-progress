module gpmdcov_EnergAndForces_mod

  use gpmdcov_mod

#ifdef USE_NVTX
    use gpmdcov_nvtx_mod
#endif
    use prg_parallel_mod
  use bml
  use ham_latte_mod
  use tbparams_latte_mod

  private :: get_dH_or_dS_vect_local, get_skblock_vect_local, bondIntegral_vect_local
  private

    integer, parameter :: dp = kind(1.0d0)

    public :: gpmdcov_EnergAndForces
    
  contains
  !> Function to calculate the bond integral for a given distance
  !! and coefficients set.
  !! \param dr distance between atoms.
  !! \param f parameters (coefficients) for the bond integral.
  function bondIntegral_vect_local(dr,f) result(x)
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
    
  end function bondIntegral_vect_local

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
  subroutine get_SKBlock_vect_local(sp,refcoord,coord,lattice_vectors&
       ,norbs,onsites,intParams1,intParams2,blk_out,atnum)
    implicit none
    integer                              ::  dimi, dimj, i, nr_shift_X
    integer                              ::  nr_shift_Y, nr_shift_Z, nats, norbsall, mask_size, iorb
    integer, intent(in)                  ::  norbs(:), sp(:), atnum
    real(dp), allocatable                ::  HPPP(:), HPPS(:), HSPS(:), HSSS(:), PPSMPP(:)
    real(dp), allocatable                ::  L(:), M(:), N(:)
    real(dp), allocatable, save                ::  dr(:), rab(:,:)
    real(dp), allocatable                :: dr_m(:)
    real(dp), allocatable                ::  onsites_m(:)
    real(dp), intent(inout)              ::  blk_out(:,:)
    real(dp), allocatable, save                ::  blk(:,:)
    real(dp), intent(in)                 ::  refcoord(:),coord(:,:), lattice_vectors(:,:)
    real(dp), intent(in)                 ::  onsites(:,:)
    real(dp), intent(in)                 ::  intParams1(:,:,:),intParams2(:,:,:)
    logical, allocatable,save                 ::  dist_mask(:), onsite_mask(:), calc_mask(:), calcs_mask(:)
    logical, allocatable,save                 ::  calcsp_mask(:), param_mask(:,:), calc_mask_for_porbs(:)
    logical, allocatable,save                 ::  sorb_mask(:), pxorb_mask(:), pyorb_mask(:), pzorb_mask(:)
    logical, allocatable,save                 ::  sorb_at_mask(:), sporb_at_mask(:)
    integer, allocatable,save                 ::  atomidx(:), atomidx_m(:), orbidx(:), orbidx_m(:), orbidx_sel(:)
    real(dp), allocatable,save                ::  intParams(:,:)

    nats = size(coord,dim=2)
    norbsall = sum(norbs)

    if(allocated(dr))then
       if(nats.ne.size(dr,dim=1))then
          deallocate(blk)
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
       allocate(blk(4,size(blk_out,dim=2)))
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

    blk(:,:)=0.0_dp

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
    allocate(intParams(mask_size,16))
    intParams(:,:) = intParams1(atomidx_m(:),:,1)
    blk(1,orbidx_sel) = BondIntegral_vect_local(dr(atomidx_m(:)),intParams)
    deallocate(intParams)
    deallocate(atomidx_m)
    deallocate(orbidx_sel)
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
       HSPS = BondIntegral_vect_local(dr(atomidx_m(:)),intParams)
       orbidx_sel = pack(orbidx_m,mask=calcs_mask) ! orbidx_m still is the s orbital mask
       if(size(orbidx_sel).ne.size(atomidx_m))then
          write(*,*)"GET_SKBLOCK_VECT: size mismatch of orbital and atom index arrays. Abort."
          stop
       endif
       blk(2,orbidx_sel) = - rab(atomidx_m(:),1)/dr(atomidx_m(:)) * HSPS
       blk(3,orbidx_sel) = - rab(atomidx_m(:),2)/dr(atomidx_m(:)) * HSPS
       blk(4,orbidx_sel) = - rab(atomidx_m(:),3)/dr(atomidx_m(:)) * HSPS
       deallocate(intParams)
       deallocate(HSPS)
       deallocate(orbidx_sel)
       deallocate(atomidx_m)
    endif
    ! Now work on the sp orbital atoms
   ! First create some common masked arrays
    atomidx_m = pack(atomidx,mask=calcsp_mask)
    mask_size = size(atomidx_m)
    dr_m = pack(dr,calcsp_mask)
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
       HSPS = BondIntegral_vect_local(dr_m,intParams)
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
    HSPS = BondIntegral_vect_local(dr_m,intParams)
    if(norbs(atnum).eq.4)then ! Calculate extra params if the reference atom is sp type
       allocate(HPPP(mask_size))
       allocate(PPSMPP(mask_size))
       intParams(:,:) = intParams1(atomidx_m(:),:,3)
       PPSMPP = BondIntegral_vect_local(dr_m,intParams)
       intParams(:,:) = intParams1(atomidx_m(:),:,4)
       HPPP = BondIntegral_vect_local(dr_m,intParams)
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
    blk(1,orbidx_sel+1) = + M(:) * HSPS(:)
    if(norbs(atnum).eq.4)then
       blk(2,orbidx_sel+1) = L*M*PPSMPP
       blk(3,orbidx_sel+1) = HPPP + M*M*PPSMPP
       blk(4,orbidx_sel+1) = N*M*PPSMPP
    endif
    blk(1,orbidx_sel+2) = + N(:) * HSPS(:)
    if(norbs(atnum).eq.4)then
       blk(2,orbidx_sel+2) = L*N*PPSMPP
       blk(3,orbidx_sel+2) = M*N*PPSMPP
       blk(4,orbidx_sel+2) = HPPP + N*N*PPSMPP
       deallocate(HPPP)
    endif
    blk_out(:,:) = blk(:,:)
    deallocate(blk)
    deallocate(orbidx_m)
    deallocatE(orbidx_sel)
    deallocate(calc_mask_for_porbs)
    deallocate(HSPS)
    deallocate(intParams)
    deallocate(L)
    deallocate(M)
    deallocate(N)
    deallocate(dr_m)
 end subroutine get_SKBlock_vect_local

        !> This routine computes the derivative of H matrix.
  !! \param dx X differential to compute the derivatives
  !! \param coords System coordinates.
  !! \param hindex Contains the Hamiltonian indices for every atom (see get_hindex).
  !! \param spindex Species indices (see system_type).
  !! \param intPairsH See defprg_inition in intPairs_type
  !! \param onsitesH Onsite energies for every orbital of a particular species.
  !! \param symbol System element symbol.
  !! \param lattice_vectors System lattece vectors.
  !! \param norb Number of total orbitals.
  !! \param norbi Number of orbitals for each atomic site.
  !! \param threshold Threshold value for matrix elements.
  !! \param dH0x_bml x derivative of H0.
  !! \param dH0y_bml y derivative of H0.
  !! \param dH0z_bml z derivative of H0.
  !!
  subroutine get_dH_or_dS_vect_local(dx,coords,hindex,spindex,intPairsH,onsitesH,symbol,lattice_vectors, norb, norbi, bml_type, &
       threshold, dH0x_bml,dH0y_bml,dH0z_bml)
    implicit none
    character(2)                       ::  Type_pair(2)
    character(2), intent(in)           ::  symbol(:)
    character(len=*), intent(in)       ::  bml_type
    integer                            ::  IDim, JDim, nats, dimi
    integer                            ::  dimj, i, ii, j
    integer                            ::  jj, l
    integer, intent(in)                ::  hindex(:,:), norb, norbi(:), spindex(:)
    integer                            ::  maxnorbi
    real(dp)                           ::  Rax_m(3), Rax_p(3), Ray_m(3), Ray_p(3)
    real(dp)                           ::  Raz_m(3), Raz_p(3), Rb(3), d, maxblockij
    real(dp), allocatable              ::  Rx(:), Ry(:), Rz(:), blockm(:,:,:)
    real(dp), allocatable              ::  blockp(:,:,:), dh0(:,:), dH0x(:,:), dH0y(:,:), dH0z(:,:)
    real(dp), intent(in)               ::  coords(:,:), dx, lattice_vectors(:,:), onsitesH(:,:)
    real(dp), intent(in)               ::  threshold
    type(bml_matrix_t), intent(inout)  ::  dH0x_bml, dH0y_bml, dH0z_bml
    type(intpairs_type), intent(in)  ::  intPairsH(:,:)
    real(dp), allocatable               ::  intParams1(:,:,:), intParams2(:,:,:)
    real(dp), allocatable               :: dHx_vect(:,:),dHy_vect(:,:),dHz_vect(:,:),ham_vect(:,:)
    integer, allocatable                ::  norbs_atidx(:)
    logical                             ::  test_accuracy
    ! integer, intent(in)                :: nnstruct(:,:)

    test_accuracy = .false.

    nats = size(coords,dim=2)

    if(bml_get_N(dH0x_bml).LT.0)then
      call bml_noinit_matrix(bml_type,bml_element_real,dp,norb,norb,dH0x_bml)
      call bml_noinit_matrix(bml_type,bml_element_real,dp,norb,norb,dH0y_bml)
      call bml_noinit_matrix(bml_type,bml_element_real,dp,norb,norb,dH0z_bml)
    else
      call bml_deallocate(dH0x_bml)
      call bml_deallocate(dH0y_bml)
      call bml_deallocate(dH0z_bml)
      call bml_noinit_matrix(bml_type,bml_element_real,dp,norb,norb,dH0x_bml)
      call bml_noinit_matrix(bml_type,bml_element_real,dp,norb,norb,dH0y_bml)
      call bml_noinit_matrix(bml_type,bml_element_real,dp,norb,norb,dH0z_bml)
    endif

    ! dH0x = zeros(HDIM,HDIM); dH0y = zeros(HDIM,HDIM); dH0z = zeros(HDIM,HDIM);

    if(test_accuracy)then       
       allocate(dH0x(norb,norb))
       allocate(dH0y(norb,norb))
       allocate(dH0z(norb,norb))
       allocate(blockm(maxnorbi,maxnorbi,nats))
       allocate(blockp(maxnorbi,maxnorbi,nats))

       dH0x = 0.0_dp
       dH0y = 0.0_dp
       dH0z = 0.0_dp
    endif
    
    allocate(dHx_vect(norb,norb))
    allocate(dHy_vect(norb,norb))
    allocate(dHz_vect(norb,norb))
    allocate(ham_vect(norb,norb))
    
    allocate(Rx(nats))
    allocate(Ry(nats))
    allocate(Rz(nats))

    Rx = coords(1,:)
    Ry = coords(2,:)
    Rz = coords(3,:)

    maxnorbi = maxval(norbi)

    if(.not.allocated(norbs_atidx))then
       allocate(norbs_atidx(nats))
       do i=1,nats
          norbs_atidx(i) = norbi(spindex(i))
       enddo
    endif

    if(test_accuracy)then
       !$omp parallel do default(none) private(i) &
       !$omp private(Rax_p,Rax_m,Ray_p,Ray_m,Raz_p,Raz_m) &
       !$omp private(dimi,J,Type_pair,dimj,Rb,maxblockij) &
       !$omp shared(nats,RX,RY,RZ,spindex,hindex,lattice_vectors, dx, threshold) &
       !$omp shared(norbi,intPairsH,onsitesH,symbol,dH0x_bml,dH0y_bml,dH0z_bml) &
       !$omp shared(blockm, blockp, dH0x, dH0y, dH0z)
       do I = 1, nats

          Type_pair(1) = symbol(i);
          Rax_p(1) = RX(I)+ dx; Rax_p(2) = RY(I); Rax_p(3) = RZ(I)
          Rax_m(1) = RX(I)- dx; Rax_m(2) = RY(I); Rax_m(3) = RZ(I)
          Ray_p(1) = RX(I); Ray_p(2) = RY(I)+dx; Ray_p(3) = RZ(I)
          Ray_m(1) = RX(I); Ray_m(2) = RY(I)-dx; Ray_m(3) = RZ(I)
          Raz_p(1) = RX(I); Raz_p(2) = RY(I); Raz_p(3) = RZ(I)+dx
          Raz_m(1) = RX(I); Raz_m(2) = RY(I); Raz_m(3) = RZ(I)-dx

          dimi = hindex(2,I)-hindex(1,I)+1;
          do J = 1,nats
             if(J .ne. I)then

                Type_pair(2) = symbol(J);
                Rb(1) = RX(J); Rb(2) = RY(J); Rb(3) = RZ(J)
                dimj = hindex(2,J)-hindex(1,J)+1;

                !! MATLAB code
                !       [fss_sigma,fsp_sigma,fpp_sigma,fpp_pi,Es,Ep,U] = LoadBondIntegralParameters_H(Type_pair); % Used in BondIntegral(dR,fxx_xx)
                !       diagonal(1:2) = [Es,Ep];
                !       dh0 = Slater_Koster_Block(IDim,JDim,Rax_p,Rb,LBox,Type_pair,fss_sigma,fsp_sigma,fpp_sigma,fpp_pi,diagonal);
                !       dH0x(H_INDEX_START(I):H_INDEX_END(I),H_INDEX_START(J):H_INDEX_END(J))  = dH0x(H_INDEX_START(I):H_INDEX_END(I),H_INDEX_START(J):H_INDEX_END(J)) + dh0/(2*dx);
                !       dh0 = Slater_Koster_Block(IDim,JDim,Rax_m,Rb,LBox,Type_pair,fss_sigma,fsp_sigma,fpp_sigma,fpp_pi,diagonal);
                !       dH0x(H_INDEX_START(I):H_INDEX_END(I),H_INDEX_START(J):H_INDEX_END(J))  = dH0x(H_INDEX_START(I):H_INDEX_END(I),H_INDEX_START(J):H_INDEX_END(J)) - dh0/(2*dx);

                call get_SKBlock(spindex(i),spindex(j),Rax_p,&
                     Rb,lattice_vectors,norbi,&
                     onsitesH,intPairsH(spindex(i),spindex(j))%intParams,intPairsH(spindex(j),spindex(i))%intParams,blockp,i)

                maxblockij = 0.0_dp
                do ii=1,dimi
                   do jj=1,dimj
                      maxblockij = max(maxblockij,abs(blockp(ii,jj,i)))
                   enddo
                enddo

                if(maxblockij.gt.0.0_dp)then

                   call get_SKBlock(spindex(i),spindex(j),Rax_m,&
                        Rb,lattice_vectors,norbi,&
                        onsitesH,intPairsH(spindex(i),spindex(j))%intParams,intPairsH(spindex(j),spindex(i))%intParams,blockm,i)

                   blockp(:,:,i) = (blockp(:,:,i) - blockm(:,:,i))/(2.0_dp*dx)

                   do jj=1,dimj
                      do ii=1,dimi
                         dH0x(hindex(1,i)-1+ii,hindex(1,j)-1+jj) = blockp(ii,jj,i)
                      enddo
                   enddo

                   !! MATLAB code
                   !       dh0 = Slater_Koster_Block(IDim,JDim,Ray_p,Rb,LBox,Type_pair,fss_sigma,fsp_sigma,fpp_sigma,fpp_pi,diagonal);
                   !       dH0y(H_INDEX_START(I):H_INDEX_END(I),H_INDEX_START(J):H_INDEX_END(J))  = dH0y(H_INDEX_START(I):H_INDEX_END(I),H_INDEX_START(J):H_INDEX_END(J)) + dh0/(2*dx);
                   !       dh0 = Slater_Koster_Block(IDim,JDim,Ray_m,Rb,LBox,Type_pair,fss_sigma,fsp_sigma,fpp_sigma,fpp_pi,diagonal);
                   !       dH0y(H_INDEX_START(I):H_INDEX_END(I),H_INDEX_START(J):H_INDEX_END(J))  = dH0y(H_INDEX_START(I):H_INDEX_END(I),H_INDEX_START(J):H_INDEX_END(J)) - dh0/(2*dx);

                   call get_SKBlock(spindex(i),spindex(j),Ray_p,&
                        Rb,lattice_vectors,norbi,&
                        onsitesH,intPairsH(spindex(i),spindex(j))%intParams,intPairsH(spindex(j),spindex(i))%intParams,blockp,i)

                   call get_SKBlock(spindex(i),spindex(j),Ray_m,&
                        Rb,lattice_vectors,norbi,&
                        onsitesH,intPairsH(spindex(i),spindex(j))%intParams,intPairsH(spindex(j),spindex(i))%intParams,blockm,i)

                   blockp(:,:,i) = (blockp(:,:,i) - blockm(:,:,i))/(2.0_dp*dx)

                   do jj=1,dimj
                      do ii=1,dimi
                         dH0y(hindex(1,i)-1+ii,hindex(1,j)-1+jj) = blockp(ii,jj,i)
                      enddo
                   enddo

                   !! MATLAB code
                   !       dh0 = Slater_Koster_Block(IDim,JDim,Raz_p,Rb,LBox,Type_pair,fss_sigma,fsp_sigma,fpp_sigma,fpp_pi,diagonal);
                   !       dH0z(H_INDEX_START(I):H_INDEX_END(I),H_INDEX_START(J):H_INDEX_END(J))  = dH0z(H_INDEX_START(I):H_INDEX_END(I),H_INDEX_START(J):H_INDEX_END(J)) + dh0/(2*dx);
                   !       dh0 = Slater_Koster_Block(IDim,JDim,Raz_m,Rb,LBox,Type_pair,fss_sigma,fsp_sigma,fpp_sigma,fpp_pi,diagonal);
                   !       dH0z(H_INDEX_START(I):H_INDEX_END(I),H_INDEX_START(J):H_INDEX_END(J))  = dH0z(H_INDEX_START(I):H_INDEX_END(I),H_INDEX_START(J):H_INDEX_END(J)) - dh0/(2*dx);

                   call get_SKBlock(spindex(i),spindex(j),Raz_p,&
                        Rb,lattice_vectors,norbi,&
                        onsitesH,intPairsH(spindex(i),spindex(j))%intParams,intPairsH(spindex(j),spindex(i))%intParams,blockp,i)

                   call get_SKBlock(spindex(i),spindex(j),Raz_m,&
                        Rb,lattice_vectors,norbi,&
                        onsitesH,intPairsH(spindex(i),spindex(j))%intParams,intPairsH(spindex(j),spindex(i))%intParams,blockm,i)

                   blockp(:,:,i) = (blockp(:,:,i) - blockm(:,:,i))/(2.0_dp*dx)

                   do jj=1,dimj
                      do ii=1,dimi
                         dH0z(hindex(1,i)-1+ii,hindex(1,j)-1+jj) = blockp(ii,jj,i)
                      enddo
                   enddo

                endif
             endif
          enddo
       enddo
       !$omp end parallel do
   endif
    
   allocate(intParams1(nats,16,4))
   allocate(intParams2(nats,16,4))

#ifdef USE_NVTX
      call nvtxStartRange("OMP loop",2)
#endif

   !$omp parallel do default(none) private(i) &
   !$omp private(Rax_p,Rax_m,Ray_p,Ray_m,Raz_p,Raz_m) &
   !$omp private(J) &
   !$omp private(intParams1, intParams2) &
   !$omp shared(nats,RX,RY,RZ,spindex,hindex,lattice_vectors, dx, threshold) &
   !$omp shared(norbs_atidx,intPairsH,onsitesH) &
   !$omp shared(coords,ham_vect) &
   !$omp shared(dHx_vect, dHy_vect, dHz_vect)
   do I = 1, nats
      do j = 1,nats
         intParams1(j,:,:) = intPairsH(spindex(i),spindex(j))%intParams(:,1:4)
         intParams2(j,:,:) = intPairsH(spindex(j),spindex(i))%intParams(:,1:4)
      enddo
      
      Rax_p(1) = RX(I)+ dx; Rax_p(2) = RY(I); Rax_p(3) = RZ(I)
      Rax_m(1) = RX(I)- dx; Rax_m(2) = RY(I); Rax_m(3) = RZ(I)
      Ray_p(1) = RX(I); Ray_p(2) = RY(I)+dx; Ray_p(3) = RZ(I)
      Ray_m(1) = RX(I); Ray_m(2) = RY(I)-dx; Ray_m(3) = RZ(I)
      Raz_p(1) = RX(I); Raz_p(2) = RY(I); Raz_p(3) = RZ(I)+dx
      Raz_m(1) = RX(I); Raz_m(2) = RY(I); Raz_m(3) = RZ(I)-dx
      
      call get_SKBlock_vect_local(spindex,rax_p,coords(:,:),lattice_vectors,norbs_atidx,&
           onsitesH,intParams1(:,:,:),intParams2(:,:,:),dHx_vect(hindex(1,i):hindex(2,i),:),i)

      call get_SKBlock_vect_local(spindex,rax_m,coords(:,:),lattice_vectors,norbs_atidx,&
           onsitesH,intParams1(:,:,:),intParams2(:,:,:),ham_vect(hindex(1,i):hindex(2,i),:),i)

      dHx_vect(hindex(1,i):hindex(2,i),:) = (dHx_vect(hindex(1,i):hindex(2,i),:) &
           - ham_vect(hindex(1,i):hindex(2,i),:))/(2.0_dp*dx)

      call get_SKBlock_vect_local(spindex,ray_p,coords(:,:),lattice_vectors,norbs_atidx,&
           onsitesH,intParams1(:,:,:),intParams2(:,:,:),dHy_vect(hindex(1,i):hindex(2,i),:),i)

      call get_SKBlock_vect_local(spindex,ray_m,coords(:,:),lattice_vectors,norbs_atidx,&
           onsitesH,intParams1(:,:,:),intParams2(:,:,:),ham_vect(hindex(1,i):hindex(2,i),:),i)

      dHy_vect(hindex(1,i):hindex(2,i),:) = (dHy_vect(hindex(1,i):hindex(2,i),:) &
           - ham_vect(hindex(1,i):hindex(2,i),:))/(2.0_dp*dx)

      call get_SKBlock_vect_local(spindex,raz_p,coords(:,:),lattice_vectors,norbs_atidx,&
           onsitesH,intParams1(:,:,:),intParams2(:,:,:),dHz_vect(hindex(1,i):hindex(2,i),:),i)

      call get_SKBlock_vect_local(spindex,raz_m,coords(:,:),lattice_vectors,norbs_atidx,&
           onsitesH,intParams1(:,:,:),intParams2(:,:,:),ham_vect(hindex(1,i):hindex(2,i),:),i)

      dHz_vect(hindex(1,i):hindex(2,i),:) = (dHz_vect(hindex(1,i):hindex(2,i),:) &
           - ham_vect(hindex(1,i):hindex(2,i),:))/(2.0_dp*dx)
   enddo
   
   !$omp end parallel do
   
#ifdef USE_NVTX
      call nvtxEndRange
#endif
      
      call prg_barrierParallel

   if(test_accuracy)then
      if(.not.all(abs(dHx_vect-dH0x).lt.1.D-9))then
         do i = 1,norb
            do j = 1,norb
               if(abs(dHx_vect(i,j)-dH0x(i,j)).ge.1.D-9)then
                  write(*,*)"GET_DH_OR_DS_VECT: Vectorized dHx differs at i,j = ",i,j,"with difference ",dHx_vect(i,j)-dH0x(i,j)
               endif
            enddo
         enddo
      endif
   endif

   call bml_import_from_dense(bml_type,dHx_vect,dH0x_bml,threshold,norb) !Dense to dense_bml
   call bml_import_from_dense(bml_type,dHy_vect,dH0y_bml,threshold,norb) !Dense to dense_bml
   call bml_import_from_dense(bml_type,dHz_vect,dH0z_bml,threshold,norb) !Dense to dense_bml
   
   if (allocated(dH0x)) then
      deallocate(dH0x)
      deallocate(dH0y)
      deallocate(dH0z)
      deallocate(blockm)
      deallocate(blockp)
   endif
    
   deallocate(dHx_vect)
   deallocate(dHy_vect)
   deallocate(dHz_vect)
   deallocate(ham_vect)

    ! stop
  end subroutine get_dH_or_dS_vect_local

  subroutine gpmdcov_EnergAndForces(charges)
    
    use gpmdcov_vars
    use gpmdcov_writeout_mod
    use gpmdcov_kernel_mod
    use gpmdcov_constraints_mod

    Implicit none
    real(dp), intent(in) :: charges(:)
    real(dp), allocatable :: ebandvector(:)
    real(dp), allocatable ::  R1(:), R2(:)
    real(dp) :: dcoords(3), dist
    real(dp) :: smd_total_force(3), smd_total_energy
    real(dp) :: smd_total_energy_allpairs
    real(dp) :: smd_test_force(3), smd_test_energy
    real(dp) :: deltas(3)
    real(dp) :: delta_h, energy_plus, energy_minus, denergy, differ
    type(rankReduceData_t) :: mpimax_in(1), mpimax_out(1)
    integer :: k
    logical :: testsmd
    call gpmdcov_msMem("gpmdcov","Before gpmd_EnergAndForces",lt%verbose,myRank)

    if(.not.allocated(coul_forces)) allocate(coul_forces(3,sy%nats))
    if(.not.allocated(GFPUL))allocate(GFPUL(3,sy%nats))
    if(.not.allocated(GFSCOUL))allocate(GFSCOUL(3,sy%nats))
    if(.not.allocated(SKForce))allocate(SKForce(3,sy%nats))
    if(.not.allocated(collectedforce))allocate(collectedforce(3,sy%nats))
    if(.not.allocated(ebandvector))allocate(ebandvector(gpat%TotalParts))

    GFPUL = 0.0_dp
    GFSCOUL = 0.0_dp
    SKForce = 0.0_dp
    collectedforce = 0.0_dp
    ebandvector = 0.0_dp
    smd_total_force(:) = 0.0_dp
    smd_total_energy = 0.0_dp
    smd_total_energy_allpairs = 0.0_dp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Loop over all the parts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef DO_MPI
    !do ipt= gpat%localPartMin(myRank), gpat%localPartMax(myRank)
    do iptt=1,partsInEachRank(myRank)
      ipt= reshuffle(iptt,myRank)
#else
    do ipt = 1,gpat%TotalParts
#endif
      !Distribute the charges back to the parts.
      do j=1,gpat%sgraph(ipt)%lsize
        jj = gpat%sgraph(ipt)%core_halo_index(j)+1
        syprt(ipt)%net_charge(j) = charges(jj)
      enddo

      norb = syprt(ipt)%estr%norbs

      norb_core = syprt(ipt)%estr%hindex(2,gpat%sgraph(ipt)%llsize)

      if(bml_get_N(aux_bml).gt.0)then
        call bml_deallocate(aux_bml)
        call bml_deallocate(aux1_bml)
        deallocate(row)
      endif

      allocate(row(norb))

      !> Get Electronic energy
      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,nOrb,nOrb,aux1_bml)
      call bml_copy_new(syprt(ipt)%estr%rho,aux_bml)


      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,nOrb,nOrb,rhoat_bml)
      call prg_build_atomic_density(rhoat_bml,tb%numel,syprt(ipt)%estr%hindex,syprt(ipt)%spindex,norb,&
           lt%bml_type)

      call bml_add_deprecated(1.0_dp,aux_bml,-1.0_dp,rhoat_bml,lt%threshold)
      call bml_multiply(aux_bml,syprt(ipt)%estr%ham,aux1_bml,1.0d0, 0.0d0,lt%threshold)
      row=0.0_dp
      call bml_deallocate(rhoat_bml)
      call bml_get_diagonal(aux1_bml,row)

      TRRHOH = 0.0_dp
      do i=1,norb_core
        TRRHOH= TRRHOH+ row(i)
      enddo

      call gpmdcov_message("gpmdcov_EnergAndForces","Energy Band for part =&
      & "//to_string(ipt)//"= "//to_string(TRRHOH),lt%verbose,myRank)

      call bml_deallocate(aux_bml)
      call bml_deallocate(aux1_bml)
      call bml_deallocate(syprt(ipt)%estr%oham)
      deallocate(row)

      syprt(ipt)%estr%eband = TRRHOH
      ebandvector(ipt) = TRRHOH

      dx = 0.0001_dp;

#ifdef USE_NVTX
      call nvtxStartRange("get_dH_and_dS",2)
#endif
      mls_i = mls()
      if(gpmdt%usevectsk)then

         call get_dH_or_dS_vect_local(dx,syprt(ipt)%coordinate,syprt(ipt)%estr%hindex,&
              &syprt(ipt)%spindex,intPairsH,onsitesH,syprt(ipt)%symbol,&
              &syprt(ipt)%lattice_vector, norb, tb%norbi, lt%bml_type, &
              &lt%threshold, dH0x_bml,dH0y_bml,dH0z_bml)
         
         call get_dH_or_dS_vect_local(dx,syprt(ipt)%coordinate,syprt(ipt)%estr%hindex,&
              &syprt(ipt)%spindex,intPairsS,onsitesS,syprt(ipt)%symbol,&
              &syprt(ipt)%lattice_vector, norb, tb%norbi, lt%bml_type, &
              &lt%threshold, dSx_bml,dSy_bml,dSz_bml)
      else
         call get_dH(dx,syprt(ipt)%coordinate,syprt(ipt)%estr%hindex,&
              &syprt(ipt)%spindex,intPairsH,onsitesH,syprt(ipt)%symbol,&
              &syprt(ipt)%lattice_vector, norb, tb%norbi, lt%bml_type, &
              &lt%threshold, dH0x_bml,dH0y_bml,dH0z_bml)

         call get_dS(dx,syprt(ipt)%coordinate,syprt(ipt)%estr%hindex,&
              &syprt(ipt)%spindex,intPairsS,onsitesS,syprt(ipt)%symbol,&
              &syprt(ipt)%lattice_vector, norb, tb%norbi, lt%bml_type, &
              &lt%threshold, dSx_bml,dSy_bml,dSz_bml)
      endif
#ifdef USE_NVTX
      call nvtxEndRange
#endif
      mls_i = mls() - mls_i
      
      mpimax_in(1)%val = mls_i
      mpimax_in(1)%rank = myRank
      
      call maxRankRealParallel(mpimax_in,mpimax_out, 1)
      ! At this point, the answer resides on process root
      
      
      call gpmdcov_msI("gpmdcov_EnergAndForces","Max time for dS and dH is "//to_string(mpimax_out(1)%val)//&
           &" on Rank "//to_string(mpimax_out(1)%rank),lt%verbose,myRank)

    
      if(printRank() == 1 .and. lt%verbose >= 10)then
        call bml_print_matrix("dH0x_bml",dH0x_bml,0,10,0,10)
        call bml_print_matrix("dH0y_bml",dH0y_bml,0,10,0,10)
        call bml_print_matrix("dH0z_bml",dH0z_bml,0,10,0,10)

        call bml_print_matrix("dSx_bml",dSx_bml,0,10,0,10)
        call bml_print_matrix("dSy_bml",dSy_bml,0,10,0,10)
        call bml_print_matrix("dSz_bml",dSz_bml,0,10,0,10)
     endif
#ifdef USE_NVTX
     call nvtxStartRange("get_skforce",3)
#endif
      call get_skforce(syprt(ipt)%nats,syprt(ipt)%estr%rho,dH0x_bml,dH0y_bml,&
           dH0z_bml,syprt(ipt)%estr%hindex,syprt(ipt)%estr%SKForce,lt%threshold)
#ifdef USE_NVTX
      call nvtxEndRange
      call nvtxStartRange("prg_get_pulayforce",4)
#endif
      
      call prg_get_pulayforce(syprt(ipt)%nats,syprt(ipt)%estr%zmat,syprt(ipt)%estr%ham,syprt(ipt)%estr%rho,&
           dSx_bml,dSy_bml,dSz_bml,syprt(ipt)%estr%hindex,syprt(ipt)%estr%FPUL,lt%threshold)
#ifdef USE_NVTX
      call nvtxEndRange
      call nvtxStartRange("get_nonortho_coul_forces",5)
#endif
      !call prg_PulayComponentT(syprt(ipt)%estr%rho,syprt(ipt)%estr%ham,syprt(ipt)%estr%zmat,syprt(ipt)%estr%FPUL,lt%threshold &
      ! &,lt%mdim,lt%bml_type,lt%verbose)


      call get_nonortho_coul_forces(syprt(ipt)%nats, norb, dSx_bml,dSy_bml,dSz_bml,&
           syprt(ipt)%estr%hindex,syprt(ipt)%spindex,syprt(ipt)%estr%rho,syprt(ipt)%net_charge,syprt(ipt)%estr%coul_pot_r,&
           syprt(ipt)%estr%coul_pot_k,tb%hubbardu,syprt(ipt)%estr%FSCOUL,lt%threshold)
#ifdef USE_NVTX
      call nvtxEndRange
#endif
      do i=1,gpat%sgraph(ipt)%llsize
        GFPUL(:,gpat%sgraph(ipt)%core_halo_index(i)+1) = syprt(ipt)%estr%FPUL(:,i)
        GFSCOUL(:,gpat%sgraph(ipt)%core_halo_index(i)+1) = syprt(ipt)%estr%FSCOUL(:,i)
        SKForce(:,gpat%sgraph(ipt)%core_halo_index(i)+1) = syprt(ipt)%estr%SKForce(:,i)
      enddo
      
      call bml_deallocate(dSx_bml)
      call bml_deallocate(dSy_bml)
      call bml_deallocate(dSz_bml)

      call bml_deallocate(dH0x_bml)
      call bml_deallocate(dH0y_bml)
      call bml_deallocate(dH0z_bml)

      !if(.not.kernel%xlbolevel1)then
      call bml_deallocate(syprt(ipt)%estr%rho)
      call bml_deallocate(syprt(ipt)%estr%ham)
      call bml_deallocate(syprt(ipt)%estr%ham0)
      !endif
      !call bml_deallocate(syprt(ipt)%estr%over)
      !call bml_deallocate(syprt(ipt)%estr%zmat)

    enddo

    collectedforce = GFPUL + GFSCOUL + SKForce
    !collectedforce =  GFSCOUL + SKForce

    call gpmdcov_msMem("gpmdcov","Before steered MD (SMD) check",lt%verbose,myRank)
    
    !> Steered MD Force (is using SMD)
    if(gpmdt%usesmd) then

      !> Allocate SMD Arrays
      !! R1 and R2 for xyz coordinates for steered atoms in each pair
      !! directional_smd_force for xyz force components
      if(.not.allocated(R1))allocate(R1(3))
      if(.not.allocated(R2))allocate(R2(3))

      do i = 1,gpmdt%smdnumpairs
        !> Determine coordinates for steered atoms in each pair 
        R1 = sy%coordinate(:,gpmdt%smdatomind1(i))
        R2 = sy%coordinate(:,gpmdt%smdatomind2(i))

        write(*,*) "gpmdcov_EnergAndForces   SMD Pair Number ",i," Atoms ",gpmdt%smdatomind1(i),&
                &" and ",gpmdt%smdatomind2(i)
        !> Call constraints subroutine, harmonic to linear
        !! collectedforce will be updated for steered atoms
        call gpmdcov_constraint_harmonicToLinear(R1, R2, smd_total_force, smd_total_energy, lt%verbose)
     
        !> Update collectedforce to include SMD force for steered atoms
        collectedforce(:,gpmdt%smdatomind1(i)) = collectedforce(:,gpmdt%smdatomind1(i)) + smd_total_force(:)
        collectedforce(:,gpmdt%smdatomind2(i)) = collectedforce(:,gpmdt%smdatomind2(i)) - smd_total_force(:)

        !> Print velocity logging for steered atoms and total SMD energy/force terms
        if(MDstep .gt. 2) then
              call gpmdcov_msI("gpmdcov_EnergAndForces","Velocity of first steered atom " &
                      &//to_string(norm2(sy%velocity(:,gpmdt%smdatomind1(i)))), lt%verbose, myRank)
              call gpmdcov_msI("gpmdcov_EnergAndForces","Velocity of second steered atom " &
                      &//to_string(norm2(sy%velocity(:,gpmdt%smdatomind2(i)))), lt%verbose, myRank)
        endif
        call gpmdcov_msI("gpmdcov_EnergAndForces","SMD Force Magnitude " &
                & // to_string(norm2(smd_total_force)),lt%verbose,myRank)
        call gpmdcov_msI("gpmdcov_EnergAndForces","SMD Total Energy " &
             & // to_string(smd_total_energy),lt%verbose,myRank)
        do k = 1,3
           dcoords(k) = modulo(((R1(k) - R2(k)) + &
                   &0.5_dp*sy%lattice_vector(k,k)),sy%lattice_vector(k,k)) - &
                   &0.5_dp * sy%lattice_vector(k,k)
        enddo
        dist = norm2(dcoords)
        call gpmdcov_msI("gpmdcov_EnergAndForces","SMD distance " &
             & // to_string(dist),lt%verbose,myRank)

        !> Add energy for current steered atom pair to total steering energy
        !! This will then be added to the total energy of the system
        smd_total_energy_allpairs = smd_total_energy_allpairs + smd_total_energy
      enddo

      !> SMD Finite Difference Test
      !! Only test for last SMD atom pair
      !! Set testsmd to true to run finite difference derivative tests
      testsmd = .false.
      if(testsmd) then
              write(*,*) "gpmdcov_EnergAndForces testing SMD"
              call gpmdcov_msI("gpmdcov_EnergAndForces","Testing SMD Derivatives: SMD force (x-direction) " &
                      &//to_string(smd_total_force(1)), lt%verbose, myRank)
              !! Test with three different delta h values
              !! Current testing is only done in the x-direction by updating the x coordinate of R1
              deltas = (/0.001, 0.01, 0.1 /)
              do k=1,3
                  delta_h = deltas(k)
                  R1(1) = R1(1) + delta_h
                  call gpmdcov_constraint_harmonicToLinear(R1, R2, smd_test_force, smd_test_energy, &
                                                           lt%verbose)
                  energy_plus = smd_test_energy
                  R1(1) = R1(1) - 2.0_dp * delta_h
                  call gpmdcov_constraint_harmonicToLinear(R1, R2, smd_test_force, smd_test_energy, &
                                                           lt%verbose)
                  energy_minus = smd_test_energy

                  !! Compute finite difference derivative of energy to compare to output force
                  denergy = (energy_minus - energy_plus)/(2.0_dp*delta_h)
                  call gpmdcov_msI("gpmdcov_EnergAndForces","Testing SMD Derivatives: d energy (x-direction) " &
                          &//to_string(denergy), lt%verbose, myRank)
                  differ = (smd_total_force(1) - denergy)
                  call gpmdcov_msI("gpmdcov_EnergAndForces","Testing SMD Derivatives: difference (force - denergy)  " &
                          &//to_string(differ), lt%verbose, myRank)
                  call gpmdcov_msI("gpmdcov_EnergAndForces","Testing SMD Derivatives: delta h used " &
                          &//to_string(delta_h), lt%verbose, myRank)

                  !> Reset R1(1) for next test
                  R1(1) = R1(1) + delta_h
              enddo
      endif

      !> Deallocate SMD Arrays
      deallocate(R1)
      deallocate(R2)
    endif


    mls_i = mls()

#ifdef DO_MPI
    if (getNRanks() .gt. 1) then
      !       call prg_sumRealReduceN(collectedforce(1,:), sy%nats)
      !       call prg_sumRealReduceN(collectedforce(2,:), sy%nats)
       !       call prg_sumRealReduceN(collectedforce(3,:), sy%nats)
      call prg_barrierParallel
      call prg_sumRealReduceN(collectedforce, sy%nats*3)
      call prg_sumRealReduceN(ebandvector, gpat%TotalParts)
    endif
#endif

    call gpmdcov_msI("gpmdcov_EnergAndForces","MPI rank finished prg_sumRealReduceN &
        &for Forces"//to_string(mls() - mls_i),lt%verbose,myRank)

    coul_forces =  coul_forces_r + coul_forces_k

    !> Get Repulsive energy and forces
    !     call get_PairPot_contrib(sy%coordinate,sy%lattice_vector,sy%spindex,ppot,PairForces,ERep)
    call get_PairPot_contrib_int(sy%coordinate,sy%lattice_vector,nl%nnIx,nl%nnIy,&
         nl%nnIz,nl%nrnnlist,nl%nnType,sy%spindex,ppot,PairForces,ERep)

    !> Get Coulombic energy
    ECoul = 0.0;
    do i = 1,sy%nats
      ECoul = ECoul + charges(i)*(tb%hubbardu(sy%spindex(i))*charges(i) + coul_pot_r(i) + coul_pot_k(i) );
    enddo

    Etot = sum(ebandvector(:)) - 0.5_dp*ECoul  + ERep + smd_total_energy_allpairs

    entropy = 0.0_dp
    if(lt%Entropy)then
        if(lt%MuCalcType == "FromParts" .or. lt%MuCalcType == "Combined")then
          call gpmdcov_getentropy(evalsAll, dvalsAll, beta, Ef, myRank,entropy, lt%verbose)
        else
          write(*,*)"ERROR: Entropy calculation is only valid for MuCalcType= FromParts, or MuCalcType= Combined."
          stop
        endif
    endif
     
    EPOT = Etot + entropy

    if((myRank == 1) .and. (lt%verbose >= 2))then
      write(*,*)"Energy Coulomb = ", ECoul
      write(*,*)"Energy Band =", sum(ebandvector(:))
      write(*,*)"Energy Repulsive = ", ERep
      write(*,*)"Energy Entropic = ", entropy
      write(*,*)"Energy Electronic (Total) =", EPot
    endif

    if(.not.allocated(sy%force))allocate(sy%force(3,sy%nats))

    !TOTAL FORCES
    sy%force =  collectedforce +  PairForces + coul_forces
    !sy%force =  SKForce + GFSCOUL + GFPUL +  PairForces + coul_forces
    !sy%force =  SKForce + GFSCOUL + GFPUL   PairForces + coul_forces
    !sy%force =  coul_forces
    !write(*,*)"FORCESSS",sy%force
    !sy%force =  SKForce + GFSCOUL +  PairForces + coul_forces
    !sy%force =  GFSCOUL 
    !sy%force = SKForce + GFSCOUL + coul_forces + PairForces
    !write(*,*)"FORCES!!!",SKForce,GFSCOUL,coul_forces
    !sy%force =   collectedforce
    !sy%force =   coul_forces

    if(myRank == 1 .and. lt%verbose >= 3)then
      write(*,*)""; write(*,*)"FPUL + FSCOUL + SKForce"
      do i = 1,sy%nats
        write(*,*)"Collected Force",i,collectedforce(1,i),collectedforce(2,i),collectedforce(3,i)
      enddo

      write(*,*)""; write(*,*)"Coulomb forces:"
      do i = 1,sy%nats
        write(*,*)"Coul Force",i,coul_forces(1,i),coul_forces(2,i),coul_forces(3,i)
      enddo

      write(*,*)""; write(*,*)"Repulsive forces:"
      do i = 1,sy%nats
        write(*,*)i,PairForces(1,i),PairForces(2,i),PairForces(3,i)
      enddo

      write(*,*)""; write(*,*) "Total forces"
      do i = 1,sy%nats
        write(*,*)i,sy%force(1,i),sy%force(2,i),sy%force(3,i)
      enddo
    endif

    deallocate(ebandvector)

    call gpmdcov_msMem("gpmdcov","After gpmd_EnergAndForces",lt%verbose,myRank)

  end subroutine gpmdcov_EnergAndForces

  subroutine gpmdcov_getentropy(evals, dvals, beta, mu, rank,entropy, verbose)
    use gpmdcov_writeout_mod
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer                ::  i, norbs
    integer, intent(in)    ::  rank
    real(dp)  ::  fermi
    real(dp), intent(in)   ::  evals(:), dvals(:)
    real(dp), intent(in)   ::  beta
    real(dp), intent(inout)   ::  mu, entropy
    integer, optional, intent(in) :: verbose
    real(dp), allocatable :: fvals(:)

    if (.not.allocated(fvals))then
       allocate(fvals(size(evals)))
    endif
    
    call gpmdcov_msMem("gpmdcov_getentropy","Getting entropic energy contribution ...",verbose,rank)
    
    norbs = size(evals, dim = 1)
    
    entropy = 0.0_dp
    fermi = 0.0_dp
    call gpmdcov_fermifunction(beta,evals,mu,fvals)
    do i = 1,norbs
      fermi = fvals(i)
      if(abs(fermi) > 10d-9 .and. abs(fermi-1.0_dp) > 10d-9)then
        entropy = entropy + (2.0_dp/beta)*dvals(i)*(fermi*log(fermi) + (1.0_dp-fermi)*log(1.0_dp-fermi))
      endif
    enddo

  end subroutine gpmdcov_getentropy

end module gpmdcov_EnergAndForces_mod

