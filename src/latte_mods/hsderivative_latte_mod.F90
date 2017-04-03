!> A module to compute the derivatives of the overlap and 
!! Hamiltonian matrices.
!! @ingroup LATTE 
!! \brief This module will be used to compute the derivatives of H and S.
!!      
module hsderivative_latte_mod

  use prg_openfiles_mod
  use bml    
  use tbparams_latte_mod
  use ham_latte_mod
  use prg_timer_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  public :: get_dH, get_dS

contains 

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
  subroutine get_dH(dx,coords,hindex,spindex,intPairsH,onsitesH,symbol,lattice_vectors, norb, norbi, bml_type, &
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
    real(dp), allocatable              ::  blockp(:,:,:), dh0(:,:)
    real(dp), intent(in)               ::  coords(:,:), dx, lattice_vectors(:,:), onsitesH(:,:)
    real(dp), intent(in)               ::  threshold
    type(bml_matrix_t), intent(inout)  ::  dH0x_bml, dH0y_bml, dH0z_bml
    type(intpairs_type), intent(in)  ::  intPairsH(:,:)
    ! integer, intent(in)                :: nnstruct(:,:)


    write(*,*)"In get_dH ..."

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

    allocate(Rx(nats))
    allocate(Ry(nats))
    allocate(Rz(nats))

    Rx = coords(1,:)
    Ry = coords(2,:)
    Rz = coords(3,:)

    maxnorbi = maxval(norbi)

    if (.not.allocated(blockm)) then
       allocate(blockm(maxnorbi,maxnorbi,nats))
    endif

    if (.not.allocated(blockp)) then
       allocate(blockp(maxnorbi,maxnorbi,nats))
    endif

    call prg_timer_start(dyn_timer,"d calc")

    !$omp parallel do default(none) private(i) &
    !$omp private(Rax_p,Rax_m,Ray_p,Ray_m,Raz_p,Raz_m) &
    !$omp private(dimi,J,Type_pair,dimj,Rb,maxblockij) &    
    !$omp shared(nats,RX,RY,RZ,spindex,hindex,lattice_vectors, dx, threshold) &
    !$omp shared(norbi,intPairsH,onsitesH,symbol,dH0x_bml,dH0y_bml,dH0z_bml) & 
    !$omp shared(blockm, blockp) 
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
              if(abs(blockp(ii,jj,i)).gt.threshold)then 
                call bml_set_element_new(dH0x_bml,hindex(1,i)-1+ii,&
                  hindex(1,j)-1+jj,blockp(ii,jj,i))
              endif
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
              if(abs(blockp(ii,jj,i)).gt.threshold)then 
                call bml_set_element_new(dH0y_bml,hindex(1,i)-1+ii,&
                  hindex(1,j)-1+jj,blockp(ii,jj,i))
              endif
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
              if(abs(blockp(ii,jj,i)).gt.threshold)then 
                call bml_set_element_new(dH0z_bml,hindex(1,i)-1+ii,&
                  hindex(1,j)-1+jj,blockp(ii,jj,i))
              endif
            enddo
          enddo  

         endif 
         endif
      enddo
    enddo
    ! $omp end parallel do
    call prg_timer_stop(dyn_timer,1)

! stop
  end subroutine get_dH

  !> This routine computes the derivative of S matrix.
  !! \param dx X differential to compute the derivatives
  !! \param coords System coordinates.
  !! \param hindex Contains the Hamiltonian indices for every atom (see get_hindex).
  !! \param spindex Species indices (see system_type).
  !! \param intPairsS See defprg_inition in intPairs_type
  !! \param onsitesS Onsite energies for every orbital of a particular species.
  !! \param symbol System element symbol.
  !! \param lattice_vectors System lattece vectors.
  !! \param norb Number of total orbitals.
  !! \param norbi Number of orbitals for each atomic site.
  !! \param threshold Threshold value for matrix elements.
  !! \param dSx_bml x derivative of H0.
  !! \param dSy_bml y derivative of H0.
  !! \param dSz_bml z derivative of H0.
  !!
  subroutine get_dS(dx,coords,hindex,spindex,intPairsS,onsitesS,symbol,lattice_vectors, norb, norbi, bml_type, &
      threshold, dSx_bml,dSy_bml,dSz_bml)
    implicit none 
    character(2)                       ::  Type_pair(2)
    character(2), intent(in)           ::  symbol(:)
    character(len=*), intent(in)       ::  bml_type
    integer                            ::  IDim, JDim, nats, dimi
    integer                            ::  dimj, i, ii, j
    integer                            ::  jj
    integer, intent(in)                ::  hindex(:,:), norb, norbi(:), spindex(:)
    integer                            ::  maxnorbi
    real(dp)                           ::  Rax_m(3), Rax_p(3), Ray_m(3), Ray_p(3)
    real(dp)                           ::  Raz_m(3), Raz_p(3), Rb(3), maxblockij
    real(dp), allocatable              ::  Rx(:), Ry(:), Rz(:), blockm(:,:,:)
    real(dp), allocatable              ::  blockp(:,:,:), dh0(:,:)
    real(dp), intent(in)               ::  coords(:,:), dx, lattice_vectors(:,:), onsitesS(:,:)
    real(dp), intent(in)               ::  threshold
    type(bml_matrix_t), intent(inout)  ::  dSx_bml, dSy_bml, dSz_bml
    type(intpairs_type), intent(in)  ::  intPairsS(:,:)

    write(*,*)"In get_dS ..."

    nats = size(coords,dim=2)

    if(bml_get_N(dSx_bml).LT.0)then 
      call bml_noinit_matrix(bml_type,bml_element_real,dp,norb,norb,dSx_bml)    
      call bml_noinit_matrix(bml_type,bml_element_real,dp,norb,norb,dSy_bml)    
      call bml_noinit_matrix(bml_type,bml_element_real,dp,norb,norb,dSz_bml)          
    else
      call bml_deallocate(dSx_bml)
      call bml_deallocate(dSy_bml)
      call bml_deallocate(dSz_bml)
      call bml_noinit_matrix(bml_type,bml_element_real,dp,norb,norb,dSx_bml)    
      call bml_noinit_matrix(bml_type,bml_element_real,dp,norb,norb,dSy_bml)    
      call bml_noinit_matrix(bml_type,bml_element_real,dp,norb,norb,dSz_bml)          
    endif      

    allocate(Rx(nats))
    allocate(Ry(nats))
    allocate(Rz(nats))

    Rx = coords(1,:)
    Ry = coords(2,:)
    Rz = coords(3,:)

    maxnorbi = maxval(norbi)

    if (.not.allocated(blockm)) then
       allocate(blockm(maxnorbi,maxnorbi,nats))
    endif

    if (.not.allocated(blockp)) then
       allocate(blockp(maxnorbi,maxnorbi,nats))
    endif

    !$omp parallel do default(none) private(i) &
    !$omp private(Rax_p,Rax_m,Ray_p,Ray_m,Raz_p,Raz_m) &
    !$omp private(dimi,J,Type_pair,dimj,Rb,maxblockij) &    
    !$omp shared(nats,RX,RY,RZ,spindex,hindex,lattice_vectors, dx, threshold) &
    !$omp shared(norbi,intPairsS,onsitesS,symbol,dSx_bml,dSy_bml,dSz_bml) &
    !$omp shared(blockm,blockp)
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

          call get_SKBlock(spindex(i),spindex(j),Rax_p,&
            Rb,lattice_vectors,norbi,&
            onsitesS,intPairsS(spindex(i),spindex(j))%intParams,intPairsS(spindex(j),spindex(i))%intParams,blockp,i)

          maxblockij = 0.0_dp
          do ii=1,dimi
            do jj=1,dimj
              maxblockij = max(maxblockij,abs(blockp(ii,jj,i)))
            enddo
          enddo 

          if(maxblockij.gt.0.0_dp)then

          call get_SKBlock(spindex(i),spindex(j),Rax_m,&
            Rb,lattice_vectors,norbi,&
            onsitesS,intPairsS(spindex(i),spindex(j))%intParams,intPairsS(spindex(j),spindex(i))%intParams,blockm,i)

          blockp(:,:,i) = (blockp(:,:,i) - blockm(:,:,i))/(2.0_dp*dx)

          do jj=1,dimj
            do ii=1,dimi
              if(abs(blockp(ii,jj,i)).gt.threshold)then 
                call bml_set_element_new(dSx_bml,hindex(1,i)-1+ii,&
                  hindex(1,j)-1+jj,blockp(ii,jj,i))
              endif
            enddo
          enddo  

          call get_SKBlock(spindex(i),spindex(j),Ray_p,&
            Rb,lattice_vectors,norbi,&
            onsitesS,intPairsS(spindex(i),spindex(j))%intParams,intPairsS(spindex(j),spindex(i))%intParams,blockp,i)

          call get_SKBlock(spindex(i),spindex(j),Ray_m,&
            Rb,lattice_vectors,norbi,&
            onsitesS,intPairsS(spindex(i),spindex(j))%intParams,intPairsS(spindex(j),spindex(i))%intParams,blockm,i)

          blockp(:,:,i) = (blockp(:,:,i) - blockm(:,:,i))/(2.0_dp*dx)

          do jj=1,dimj
            do ii=1,dimi
              if(abs(blockp(ii,jj,i)).gt.threshold)then 
                call bml_set_element_new(dSy_bml,hindex(1,i)-1+ii,&
                  hindex(1,j)-1+jj,blockp(ii,jj,i))
              endif
            enddo
          enddo  

          call get_SKBlock(spindex(i),spindex(j),Raz_p,&
            Rb,lattice_vectors,norbi,&
            onsitesS,intPairsS(spindex(i),spindex(j))%intParams,intPairsS(spindex(j),spindex(i))%intParams,blockp,i)

          call get_SKBlock(spindex(i),spindex(j),Raz_m,&
            Rb,lattice_vectors,norbi,&
            onsitesS,intPairsS(spindex(i),spindex(j))%intParams,intPairsS(spindex(j),spindex(i))%intParams,blockm,i)

          blockp(:,:,i) = (blockp(:,:,i) - blockm(:,:,i))/(2.0_dp*dx)

          do jj=1,dimj
            do ii=1,dimi
              if(abs(blockp(ii,jj,i)).gt.threshold)then 
                call bml_set_element_new(dSz_bml,hindex(1,i)-1+ii,&
                  hindex(1,j)-1+jj,blockp(ii,jj,i))
              endif
            enddo
          enddo  

        endif
        endif
      enddo
    enddo
    !$omp end parallel do

  end subroutine get_dS

end module hsderivative_latte_mod
