!!> Gets the NonOrtho Coulombic forces. 
!! \ingroup PROGRESS 
!!
module nonorthocoulombforces_latte_mod

  use bml
  use timer_mod

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

    write(*,*)"In get_nonortho_coul_forces ..."

    if(.not.allocated(FSCOUL))then 
      allocate(FSCOUL(3,nats))
    endif 

    FSCOUL = 0.0_dp

    allocate(Coulomb_Pot(nats))
    Coulomb_Pot = Coulomb_Pot_r + Coulomb_Pot_k

    norb = bml_get_N(rho_bml)
    bml_type = bml_get_type(dSx_bml)

    call bml_noinit_matrix(bml_type,bml_element_real,dp,norb,norb,rhot_bml)          

    call bml_transpose(rho_bml,rhot_bml)

    allocate(row1(norb))
    allocate(row2(norb))
    allocate(row2x(norb))
    allocate(row2y(norb))
    allocate(row2z(norb))

    allocate(chunk(norb,10))
    allocate(chunkx(norb,10))  
    allocate(chunky(norb,10))  
    allocate(chunkz(norb,10))          

    allocate(dDSX(norb))
    allocate(dDSY(norb))
    allocate(dDSZ(norb))

    dDSX = 0.0_dp
    dDSY = 0.0_dp
    dDSZ = 0.0_dp

    !$omp parallel do default(none) private(i) &
    !$omp private(I_A,I_B,j,J_A,J_B,row1,row2,row2x,row2y,row2z,jj) &
    !$omp private(dDSX,dDSY,dDSZ,dQLxdR,dQLydR,dQLzdR,count1) &
    !$omp private(chunk,chunkx,chunky,chunkz) &    
    !$omp shared(nats,hindex,norb,rho_bml,dsx_bml,dsy_bml,dsz_bml,threshold,rhot_bml) &
    !$omp shared(FSCOUL,hubbardu,spindex,charges,coulomb_pot)    
    do I = 1,nats
      I_A = hindex(1,I);
      I_B = hindex(2,I);

      dDSX = 0.0_dp
      dDSY = 0.0_dp
      dDSZ = 0.0_dp

      do j = I_A,I_B
        row1 =0.0_dp; row2x =0.0_dp; row2y =0.0_dp; row2z =0.0_dp        
        call bml_get_row(rho_bml,j,row1)
        call bml_get_row(dSx_bml,j,row2x)        
        call bml_get_row(dSy_bml,j,row2y)        
        call bml_get_row(dSz_bml,j,row2z)        
        do jj=1,norb
          if(abs(row1(jj)).gt.threshold)then
            dDSX(j) = dDSX(j) + row1(jj)*row2x(jj);
            dDSY(j) = dDSY(j) + row1(jj)*row2y(jj);
            dDSZ(j) = dDSZ(j) + row1(jj)*row2z(jj);
          endif
        enddo
      enddo

      count1=0   
      do jj=I_A,I_B
        count1 = count1+1
        chunk(:,count1)=0.0_dp; chunkx(:,count1)=0.0_dp; 
        chunky(:,count1)=0.0_dp; chunkz(:,count1)=0.0_dp
        call bml_get_row(rhot_bml,jj,chunk(:,count1))
        call bml_get_row(dSx_bml,jj,chunkx(:,count1))
        call bml_get_row(dSy_bml,jj,chunky(:,count1))
        call bml_get_row(dSz_bml,jj,chunkz(:,count1))
      enddo

      count1=0
      do jj=I_A,I_B 
        count1 = count1+1
        do j = 1,norb
          if(abs(chunk(j,count1)).gt.threshold)then
            dDSX(j) = dDSX(j) + chunk(j,count1)*chunkx(j,count1)
            dDSY(j) = dDSY(j) + chunk(j,count1)*chunky(j,count1)
            dDSZ(j) = dDSZ(j) + chunk(j,count1)*chunkz(j,count1)
          endif
        enddo          
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
    deallocate(chunk)
    deallocate(chunkx)
    deallocate(chunky)
    deallocate(chunkz)
    deallocate(row1)
    deallocate(row2x)
    deallocate(row2y)
    deallocate(row2z)

    call bml_deallocate(rhot_bml)

  end subroutine get_nonortho_coul_forces

end module nonorthocoulombforces_latte_mod
