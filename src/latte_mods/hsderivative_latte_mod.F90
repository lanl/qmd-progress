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

  public :: get_dH, get_dS, get_dH_or_dS, get_dH_or_dS_vect

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
    real(dp), allocatable              ::  blockp(:,:,:), dh0(:,:), dH0x(:,:), dH0y(:,:), dH0z(:,:)
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

    if (.not.allocated(dH0x)) then
      allocate(dH0x(norb,norb))
    endif
    if (.not.allocated(dH0y)) then
      allocate(dH0y(norb,norb))
    endif
    if (.not.allocated(dH0z)) then
      allocate(dH0z(norb,norb))
    endif

    dH0x = 0.0_dp
    dH0y = 0.0_dp
    dH0z = 0.0_dp

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
    call bml_import_from_dense(bml_type,dH0x,dH0x_bml,threshold,norb) !Dense to dense_bml
    call bml_import_from_dense(bml_type,dH0y,dH0y_bml,threshold,norb) !Dense to dense_bml
    call bml_import_from_dense(bml_type,dH0z,dH0z_bml,threshold,norb) !Dense to dense_bml

    if (allocated(dH0x)) then
      deallocate(dH0x)
    endif
    if (allocated(dH0y)) then
      deallocate(dH0y)
    endif
    if (allocated(dH0z)) then
      deallocate(dH0z)
    endif

    ! stop
  end subroutine get_dH

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
  subroutine get_dH_or_dS(dx,coords,hindex,spindex,intPairsH,onsitesH,symbol,lattice_vectors, norb, norbi, bml_type, &
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
    real(dp), allocatable              ::  dH0x(:,:), dH0y(:,:), dH0z(:,:)
    real(dp), allocatable              ::  H0xm(:,:), H0ym(:,:), H0zm(:,:)
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

    if (.not.allocated(dH0x)) then
      allocate(dH0x(norb,norb))
      allocate(dH0y(norb,norb))
      allocate(dH0z(norb,norb))
      allocate(H0xm(norb,norb))
      allocate(H0ym(norb,norb))
      allocate(H0zm(norb,norb))
    endif

    dH0x = 0.0_dp
    dH0y = 0.0_dp
    dH0z = 0.0_dp
    H0xm = 0.0_dp
    H0ym = 0.0_dp
    H0zm = 0.0_dp

    allocate(Rx(nats))
    allocate(Ry(nats))
    allocate(Rz(nats))

    Rx = coords(1,:)
    Ry = coords(2,:)
    Rz = coords(3,:)

    maxnorbi = maxval(norbi)

    !$omp parallel do default(none) private(i) &
    !$omp private(Rax_p,Rax_m,Ray_p,Ray_m,Raz_p,Raz_m) &
    !$omp private(dimi,J,Type_pair,dimj,Rb,maxblockij) &
    !$omp shared(nats,RX,RY,RZ,spindex,hindex,lattice_vectors, dx, threshold) &
    !$omp shared(norbi,intPairsH,onsitesH,symbol,dH0x_bml,dH0y_bml,dH0z_bml) &
    !$omp shared(dH0x, dH0y, dH0z, H0xm, H0ym, H0zm)
    do I = 1, nats
      do J = 1,nats
         Type_pair(1) = symbol(i);
         Rax_p(1) = RX(I)+ dx; Rax_p(2) = RY(I); Rax_p(3) = RZ(I)
         Rax_m(1) = RX(I)- dx; Rax_m(2) = RY(I); Rax_m(3) = RZ(I)
         Ray_p(1) = RX(I); Ray_p(2) = RY(I)+dx; Ray_p(3) = RZ(I)
         Ray_m(1) = RX(I); Ray_m(2) = RY(I)-dx; Ray_m(3) = RZ(I)
         Raz_p(1) = RX(I); Raz_p(2) = RY(I); Raz_p(3) = RZ(I)+dx
         Raz_m(1) = RX(I); Raz_m(2) = RY(I); Raz_m(3) = RZ(I)-dx
         
         dimi = hindex(2,I)-hindex(1,I)+1;
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
            !       dh0 = Slater_Koster_Block(IDim,JDim,Ray_p,Rb,LBox,Type_pair,fss_sigma,fsp_sigma,fpp_sigma,fpp_pi,diagonal);
            !       dH0y(H_INDEX_START(I):H_INDEX_END(I),H_INDEX_START(J):H_INDEX_END(J))  = dH0y(H_INDEX_START(I):H_INDEX_END(I),H_INDEX_START(J):H_INDEX_END(J)) + dh0/(2*dx);
            !       dh0 = Slater_Koster_Block(IDim,JDim,Ray_m,Rb,LBox,Type_pair,fss_sigma,fsp_sigma,fpp_sigma,fpp_pi,diagonal);
            !       dH0y(H_INDEX_START(I):H_INDEX_END(I),H_INDEX_START(J):H_INDEX_END(J))  = dH0y(H_INDEX_START(I):H_INDEX_END(I),H_INDEX_START(J):H_INDEX_END(J)) - dh0/(2*dx);
            !       dh0 = Slater_Koster_Block(IDim,JDim,Raz_p,Rb,LBox,Type_pair,fss_sigma,fsp_sigma,fpp_sigma,fpp_pi,diagonal);
            !       dH0z(H_INDEX_START(I):H_INDEX_END(I),H_INDEX_START(J):H_INDEX_END(J))  = dH0z(H_INDEX_START(I):H_INDEX_END(I),H_INDEX_START(J):H_INDEX_END(J)) + dh0/(2*dx);
            !       dh0 = Slater_Koster_Block(IDim,JDim,Raz_m,Rb,LBox,Type_pair,fss_sigma,fsp_sigma,fpp_sigma,fpp_pi,diagonal);
            !       dH0z(H_INDEX_START(I):H_INDEX_END(I),H_INDEX_START(J):H_INDEX_END(J))  = dH0z(H_INDEX_START(I):H_INDEX_END(I),H_INDEX_START(J):H_INDEX_END(J)) - dh0/(2*dx);
            
            call get_SKBlock_inplace(spindex(i),spindex(j),Rax_p,&
                 Rb,lattice_vectors,norbi,&
                 onsitesH,intPairsH(spindex(i),spindex(j))%intParams, &
                 intPairsH(spindex(j),spindex(i))%intParams, &
                 dH0x(hindex(1,i):hindex(2,i),hindex(1,j):hindex(2,j)),i)
            
            if(maxval(abs(dH0x(hindex(1,i):hindex(2,i),hindex(1,j):hindex(2,j)))) &
                 .gt.0.0_dp)then
               
               call get_SKBlock_inplace(spindex(i),spindex(j),Rax_m,&
                    Rb,lattice_vectors,norbi,&
                    onsitesH,intPairsH(spindex(i),spindex(j))%intParams, &
                    intPairsH(spindex(j),spindex(i))%intParams, &
                    H0xm(hindex(1,i):hindex(2,i),hindex(1,j):hindex(2,j)),i)
               
               call get_SKBlock_inplace(spindex(i),spindex(j),Ray_p,&
                    Rb,lattice_vectors,norbi,&
                    onsitesH,intPairsH(spindex(i),spindex(j))%intParams, &
                    intPairsH(spindex(j),spindex(i))%intParams, &
                    dH0y(hindex(1,i):hindex(2,i),hindex(1,j):hindex(2,j)),i)
               
               call get_SKBlock_inplace(spindex(i),spindex(j),Ray_m,&
                    Rb,lattice_vectors,norbi,&
                    onsitesH,intPairsH(spindex(i),spindex(j))%intParams, &
                    intPairsH(spindex(j),spindex(i))%intParams, &
                    H0ym(hindex(1,i):hindex(2,i),hindex(1,j):hindex(2,j)),i)
               
               call get_SKBlock_inplace(spindex(i),spindex(j),Raz_p,&
                    Rb,lattice_vectors,norbi,&
                    onsitesH,intPairsH(spindex(i),spindex(j))%intParams, &
                    intPairsH(spindex(j),spindex(i))%intParams, &
                    dH0z(hindex(1,i):hindex(2,i),hindex(1,j):hindex(2,j)),i)
               
               call get_SKBlock_inplace(spindex(i),spindex(j),Raz_m,&
                    Rb,lattice_vectors,norbi,&
                    onsitesH,intPairsH(spindex(i),spindex(j))%intParams, &
                    intPairsH(spindex(j),spindex(i))%intParams, &
                    H0zm(hindex(1,i):hindex(2,i),hindex(1,j):hindex(2,j)),i)
               
            endif
         endif
      enddo
   enddo
   
   !$omp end parallel do

   !call bml_print_matrix("dH0x",dH0x,1,10,1,10)
   !call bml_print_matrix("H0xm",H0xm,1,10,1,10)
   dH0x = (dH0x - H0xm)/(2.0_dp*dx)
   dH0y = (dH0y - H0ym)/(2.0_dp*dx)
   dH0z = (dH0z - H0zm)/(2.0_dp*dx)
            
    call bml_import_from_dense(bml_type,dH0x,dH0x_bml,threshold,norb) !Dense to dense_bml
    call bml_import_from_dense(bml_type,dH0y,dH0y_bml,threshold,norb) !Dense to dense_bml
    call bml_import_from_dense(bml_type,dH0z,dH0z_bml,threshold,norb) !Dense to dense_bml

    if (allocated(dH0x)) then
      deallocate(dH0x)
      deallocate(dH0y)
      deallocate(dH0z)
      deallocate(H0xm)
      deallocate(H0ym)
      deallocate(H0zm)
    endif

    ! stop
  end subroutine get_dH_or_dS

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
    real(dp), allocatable              ::  blockp(:,:,:), dh0(:,:),dSx(:,:),dSy(:,:),dSz(:,:)
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

    if (.not.allocated(dSx)) then
      allocate(dSx(norb,norb))
    endif
    if (.not.allocated(dSy)) then
      allocate(dSy(norb,norb))
    endif
    if (.not.allocated(dSz)) then
      allocate(dSz(norb,norb))
    endif

    allocate(Rx(nats))
    allocate(Ry(nats))
    allocate(Rz(nats))

    dSx = 0.0_dp
    dSy = 0.0_dp
    dSz = 0.0_dp

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
    !$omp shared(blockm,blockp,dSx,dSy,dSz)
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
                dSx(hindex(1,i)-1+ii,hindex(1,j)-1+jj) = blockp(ii,jj,i)
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
                dSy(hindex(1,i)-1+ii,hindex(1,j)-1+jj) = blockp(ii,jj,i)
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
                dSz(hindex(1,i)-1+ii,hindex(1,j)-1+jj) = blockp(ii,jj,i)
              enddo
            enddo

          endif
        endif
      enddo
    enddo
    !$omp end parallel do
    call bml_import_from_dense(bml_type,dSx,dSx_bml,threshold,norb) !Dense to dense_bml
    call bml_import_from_dense(bml_type,dSy,dSy_bml,threshold,norb) !Dense to dense_bml
    call bml_import_from_dense(bml_type,dSz,dSz_bml,threshold,norb) !Dense to dense_bml

    if (allocated(dSx)) then
      deallocate(dSx)
    endif
    if (allocated(dSy)) then
      deallocate(dSy)
    endif
    if (allocated(dSz)) then
      deallocate(dSz)
    endif

  end subroutine get_dS

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
  subroutine get_dH_or_dS_vect(dx,coords,hindex,spindex,intPairsH,onsitesH,symbol,lattice_vectors, norb, norbi, bml_type, &
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
      
      call get_SKBlock_vect(spindex,rax_p,coords(:,:),lattice_vectors,norbs_atidx,&
           onsitesH,intParams1(:,:,:),intParams2(:,:,:),dHx_vect(hindex(1,i):hindex(2,i),:),i)

      call get_SKBlock_vect(spindex,rax_m,coords(:,:),lattice_vectors,norbs_atidx,&
           onsitesH,intParams1(:,:,:),intParams2(:,:,:),ham_vect(hindex(1,i):hindex(2,i),:),i)

      dHx_vect(hindex(1,i):hindex(2,i),:) = (dHx_vect(hindex(1,i):hindex(2,i),:) &
           - ham_vect(hindex(1,i):hindex(2,i),:))/(2.0_dp*dx)

      call get_SKBlock_vect(spindex,ray_p,coords(:,:),lattice_vectors,norbs_atidx,&
           onsitesH,intParams1(:,:,:),intParams2(:,:,:),dHy_vect(hindex(1,i):hindex(2,i),:),i)

      call get_SKBlock_vect(spindex,ray_m,coords(:,:),lattice_vectors,norbs_atidx,&
           onsitesH,intParams1(:,:,:),intParams2(:,:,:),ham_vect(hindex(1,i):hindex(2,i),:),i)

      dHy_vect(hindex(1,i):hindex(2,i),:) = (dHy_vect(hindex(1,i):hindex(2,i),:) &
           - ham_vect(hindex(1,i):hindex(2,i),:))/(2.0_dp*dx)

      call get_SKBlock_vect(spindex,raz_p,coords(:,:),lattice_vectors,norbs_atidx,&
           onsitesH,intParams1(:,:,:),intParams2(:,:,:),dHz_vect(hindex(1,i):hindex(2,i),:),i)

      call get_SKBlock_vect(spindex,raz_m,coords(:,:),lattice_vectors,norbs_atidx,&
           onsitesH,intParams1(:,:,:),intParams2(:,:,:),ham_vect(hindex(1,i):hindex(2,i),:),i)

      dHz_vect(hindex(1,i):hindex(2,i),:) = (dHz_vect(hindex(1,i):hindex(2,i),:) &
           - ham_vect(hindex(1,i):hindex(2,i),:))/(2.0_dp*dx)
   enddo
   
   !$omp end parallel do

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
  end subroutine get_dH_or_dS_vect

end module hsderivative_latte_mod
