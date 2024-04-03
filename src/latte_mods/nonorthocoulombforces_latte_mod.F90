!!> Gets the NonOrtho Coulombic forces.
!! \ingroup PROGRESS
!!
module nonorthocoulombforces_latte_mod

  use bml
  use prg_timer_mod

  implicit none

  private

  integer,parameter :: dp = kind(1.0d0)

  public :: get_nonortho_coul_forces

contains

  !> Nonortho Coulombic Forces
  !!  Coulomb force FSCOUL from nonorthogonality
  !! \param nats Number of atoms.
  !! \param norb Number of orbitals. Its usually the dimension of the Hamiltonian.
  !! \param dSx_bml S derivative in the x direction.
  !! \param dSy_bml S derivative in the y direction.
  !! \param dSz_bml S derivative in the z direction.
  !! \param hindex Start and end index for every atom in the system.
  !! \param spindex Species index. It gives the species index of a particulat atom. See system_type.
  !! \param rho_bml Density matrix in bml format.
  !! \param charges Charges for every atom in the system.
  !! \param Coulomb_Pot_r Coulombic potential (real space contribution).
  !! \param Coulomb_Pot_k Coulombic potential (reciprocal space contribution).
  !! \param hubbardu Hubbard parameter U. This is the onsite e-e potential repulsion. it runs over the species list.
  !! \param FSCOUL Nonortho coulombinc contribution to the force.
  !! \param threshold Threshold value for sparse matrices.
  !!
  subroutine get_nonortho_coul_forces(nats, norb, dSx_bml,dSy_bml,dSz_bml,&
       hindex,spindex,rho_bml,charges,Coulomb_Pot_r,Coulomb_Pot_k,hubbardu,FSCOUL,threshold)
    implicit none
    character(20)                        ::  bml_type
    integer                              ::  I_A, I_B, J_A, J_B
    integer                              ::  count1, i, j, jj
    integer                              ::  norb, norbs
    integer, intent(in)                  ::  nats, hindex(:,:), spindex(:)
    real(dp)                             ::  dQLxdR, dQLydR, dQLzdR, partrace
    real(dp), allocatable                ::  Coulomb_Pot(:), chunk(:,:), chunkx(:,:), chunky(:,:)
    real(dp), allocatable                ::  chunkz(:,:), dDSx(:), dDSy(:), dDSz(:)
    real(dp), allocatable                ::  diagxtmp(:), diagytmp(:), diagztmp(:), row1(:)
    real(dp), allocatable                ::  row2(:), row2x(:), row2y(:), row2z(:)
    real(dp), allocatable, intent(inout)  ::  FSCOUL(:,:)
    real(dp), intent(in)                 ::  Coulomb_Pot_k(:), Coulomb_Pot_r(:), charges(:), hubbardu(:)
    real(dp), intent(in)                 ::  threshold
    type(bml_matrix_t)                   ::  Xtmp_bml, Ytmp_bml, Ztmp_bml, rhot_bml
    type(bml_matrix_t), intent(in)       ::  dSx_bml, dSy_bml, dSz_bml, rho_bml
    real(dp), allocatable                ::  dSx_dense(:,:), dSy_dense(:,:), dSz_dense(:,:)
    real(dp), allocatable                ::  rho_dense(:,:)
    logical, allocatable                 ::  thresh_mask(:,:)

    write(*,*)"In get_nonortho_coul_forces ..."

    if(.not.allocated(FSCOUL))then
      allocate(FSCOUL(3,nats))
    endif

    FSCOUL = 0.0_dp

    allocate(Coulomb_Pot(nats))
    Coulomb_Pot = Coulomb_Pot_r + Coulomb_Pot_k

    norb = bml_get_N(rho_bml)
    bml_type = bml_get_type(dSx_bml)

    allocate(dSx_dense(norb,norb))
    allocate(dSy_dense(norb,norb))
    allocate(dSz_dense(norb,norb))
    allocate(rho_dense(norb,norb))
    allocate(thresh_mask(norb,norb))

    call bml_export_to_dense(dSx_bml,dSx_dense)
    call bml_export_to_dense(dSy_bml,dSy_dense)
    call bml_export_to_dense(dSz_bml,dSz_dense)
    call bml_export_to_dense(rho_bml,rho_dense)

    allocate(dDSX(norb))
    allocate(dDSY(norb))
    allocate(dDSZ(norb))

    dDSX = 0.0_dp
    dDSY = 0.0_dp
    dDSZ = 0.0_dp

    thresh_mask = abs(rho_dense).gt.threshold
    
    !$omp parallel do default(none) private(i) &
    !$omp private(I_A,I_B,j,J_A,J_B,row1,row2,row2x,row2y,row2z,jj) &
    !$omp private(dDSX,dDSY,dDSZ,dQLxdR,dQLydR,dQLzdR,count1) &
    !$omp shared(nats,hindex,norb,rho_dense) &
    !$omp shared(dsx_dense,dsy_dense,dsz_dense,thresh_mask) &
    !$omp shared(FSCOUL,hubbardu,spindex,charges,coulomb_pot)
    do I = 1,nats
      I_A = hindex(1,I);
      I_B = hindex(2,I);

      dDSX = 0.0_dp
      dDSY = 0.0_dp
      dDSZ = 0.0_dp

      do j = I_A,I_B
        row1 = pack(rho_dense(j,:),thresh_mask(j,:))
        row2x = pack(dSx_dense(j,:),thresh_mask(j,:))
        row2y = pack(dSy_dense(j,:),thresh_mask(j,:))
        row2z = pack(dSz_dense(j,:),thresh_mask(j,:))
        dDSX(j) = dDSX(j) + sum(row1(:)*row2x(:))
        dDSY(j) = dDSY(j) + sum(row1(:)*row2y(:))
        dDSZ(j) = dDSZ(j) + sum(row1(:)*row2z(:))
        row1 = pack(rho_dense(:,j),thresh_mask(:,j))
        row2x = pack(dSx_dense(j,:),thresh_mask(:,j))
        row2y = pack(dSy_dense(j,:),thresh_mask(:,j))
        row2z = pack(dSz_dense(j,:),thresh_mask(:,j))
        dDSX(:) = dDSX(:) + unpack(row1(:)*row2x(:),thresh_mask(j,:),0.0_dp)
        dDSY(:) = dDSY(:) + unpack(row1(:)*row2y(:),thresh_mask(j,:),0.0_dp)
        dDSZ(:) = dDSZ(:) + unpack(row1(:)*row2z(:),thresh_mask(j,:),0.0_dp)
      enddo

      do J = 1,nats
        J_A = hindex(1,J);
        J_B = hindex(2,J);
        dQLxdR = 0.0_dp ; dQLydR = 0.0_dp ; dQLzdR = 0.0_dp
        do jj=J_A,J_B
          dQLxdR = dQLxdR + dDSX(jj);
          dQLydR = dQLydR + dDSY(jj);
          dQLzdR = dQLzdR + dDSZ(jj);
        enddo
        FSCOUL(1,I) = FSCOUL(1,I) - &
             dQLxdR*(hubbardu(spindex(J))*charges(J) + Coulomb_Pot(J));
        FSCOUL(2,I) = FSCOUL(2,I) - &
             dQLydR*(hubbardu(spindex(J))*charges(J) + Coulomb_Pot(J));
        FSCOUL(3,I) = FSCOUL(3,I) - &
             dQLzdR*(hubbardu(spindex(J))*charges(J) + Coulomb_Pot(J));
      enddo
    enddo
    !$omp end parallel do

    deallocate(dDSX)
    deallocate(dDSY)
    deallocate(dDSZ)
    deallocate(thresh_mask)
    deallocate(dSx_dense)
    deallocate(dSy_dense)
    deallocate(dSz_dense)
    deallocate(rho_dense)

    call bml_deallocate(rhot_bml)

  end subroutine get_nonortho_coul_forces

end module nonorthocoulombforces_latte_mod
