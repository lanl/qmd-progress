!> Produces a matrix to get the Pulay Component of the forces. 
!! \ingroup PROGRESS 
!! Please see Niklasson 2008 \cite Niklasson2008
!!
module prg_PulayComponent_mod

  use bml

  implicit none
  
  private

  integer,parameter :: dp = kind(1.0d0)
 
  public :: prg_PulayComponent0, prg_PulayComponentT, prg_get_pulayforce
    
contains

  !> At \f$T=0 K\f$, \f$ P = \rho H \rho \f$
  !!
  !! \param rho_bml Density matrix in bml format. 
  !! \param ham_bml Hamiltonian matrix in bml format.  
  !! \param pcm_bml Pulay matix output in bml format.
  !! \param threshold Threshold for the matrix elements.
  !! \param M Maximum nonzero values per row. 
  !! \param bml_type Bml format type.
  !! \param verbose Verbosity level.
  !! \todo M and bml_type will have to be removed from the input parameter. 
  !!
  subroutine prg_PulayComponent0(rho_bml,ham_bml,pcm_bml,threshold,M,&
    &bml_type,verbose)
  
    implicit none
    
    type(bml_matrix_t), intent(in) :: rho_bml
    type(bml_matrix_t), intent(in) :: ham_bml
    type(bml_matrix_t) :: aux_bml
    type(bml_matrix_t), intent(inout) :: pcm_bml
    integer, intent(in) :: M
    integer :: nOrb, verbose 
    real(dp), intent(in) :: threshold
    character(20), intent(in) :: bml_type

    if(verbose.EQ.1) write(*,*)"In prg_PulayComponent0 ..."
    
    nOrb = bml_get_N(rho_bml)
        
    call bml_zero_matrix(bml_type,bml_element_real,dp,nOrb ,nOrb,aux_bml)
 
    if(bml_get_N(pcm_bml).LE.0)then !If pcm is not allocated  
      call bml_zero_matrix(bml_type,bml_element_real,dp,nOrb,nOrb,pcm_bml)
    else
      call bml_deallocate(pcm_bml) !If pcm is allocated then we set it to 0
      call bml_zero_matrix(bml_type,bml_element_real,dp,nOrb,nOrb,pcm_bml)      
    endif 
 
    call bml_zero_matrix(bml_type,bml_element_real,dp,nOrb ,nOrb,pcm_bml)                 
       
    call bml_multiply(rho_bml, ham_bml, aux_bml, 1.0d0, 1.0d0,threshold) !D*H

    call bml_multiply(aux_bml,rho_bml,pcm_bml, 1.0d0, 1.0d0,threshold)  !(D*H)*D

    call bml_scale(0.5_dp, pcm_bml)

    call bml_deallocate(aux_bml)
    
  end subroutine prg_PulayComponent0
    
  !> At \f$ T > 0K \f$, \f$ P = \rho H S^-1 + S^{-1} H \rho \f$
  !!
  !! \param rho_bml Density matrix in bml format. 
  !! \param ham_bml Hamiltonian matrix in bml format. 
  !! \param Z_bml Congruence transform in bml format.  
  !! \param pcm_bml Pulay matrix output in bml format.
  !! \param threshold Threshold for the matrix elements.
  !! \param M Maximum nonzero values per row. 
  !! \param bml_type Bml format type.
  !! \param verbose Verbosity level.
  !! \todo M and bml_type will have to be removed from the input parameter.   
  !!
  subroutine prg_PulayComponentT(rho_bml,ham_bml,zmat_bml,pcm_bml,threshold &
    &,M,bml_type,verbose)
  
    implicit none
    
    type(bml_matrix_t), intent(in) :: rho_bml
    type(bml_matrix_t), intent(in) :: ham_bml
    type(bml_matrix_t), intent(in) :: zmat_bml        
    type(bml_matrix_t) :: aux_bml
    type(bml_matrix_t) :: aux1_bml    
    type(bml_matrix_t), intent(inout) :: pcm_bml    
    integer, intent(in) :: M
    integer ::  nOrb, verbose 
    real(dp), intent(in) :: threshold
    character(20), intent(in) :: bml_type

    if(verbose.EQ.1) write(*,*)"In prg_PulayComponentT ..."
    
    nOrb = bml_get_N(rho_bml)
    
    if(bml_get_N(pcm_bml).LE.0)then !If pcm is not allocated.
      call bml_zero_matrix(bml_type,bml_element_real,dp,nOrb,nOrb,pcm_bml)
    else
      call bml_deallocate(pcm_bml) !If pcm is allocated then we set it to 0.
      call bml_zero_matrix(bml_type,bml_element_real,dp,nOrb,nOrb,pcm_bml)      
    endif 
        
    call bml_zero_matrix(bml_type,bml_element_real,dp,nOrb ,nOrb,aux_bml)                 
    call bml_zero_matrix(bml_type,bml_element_real,dp,nOrb ,nOrb,aux1_bml)                 
    
    if(bml_get_N(zmat_bml).LT.0)then 
      stop 'ERROR: zmat_bml not allocated'
    endif         
    
    call bml_multiply(zmat_bml,zmat_bml,aux_bml,1.0d0,1.0d0,threshold)  !aux=Z*Z

    call bml_multiply(aux_bml,ham_bml,pcm_bml,1.0d0,1.0d0,threshold)  !Z*Z*H

    call bml_multiply(pcm_bml,rho_bml,aux1_bml,1.0d0,1.0d0,threshold) !aux1=Z*Z*H*D
       
    call bml_scale(0.0_dp, pcm_bml)
    
    call bml_multiply(ham_bml,aux_bml,pcm_bml,1.0d0,1.0d0,threshold)  !H*Z*Z

    call bml_scale(0.0_dp, aux_bml)  !We set aux to 0 to recicle the variable

    call bml_multiply(rho_bml,pcm_bml,aux_bml,1.0d0,1.0d0,threshold) !D*H*Z*Z
    
    call bml_add_deprecated(0.5d0,aux_bml,0.5d0,aux1_bml)

    call bml_copy(aux_bml,pcm_bml)

    call bml_deallocate(aux_bml)
    call bml_deallocate(aux1_bml)    
    
    
  end subroutine prg_PulayComponentT


  !> Pulay Force FPUL from  \f$ 2Tr[ZZ'HD \frac{dS}{dR}] \f$
  !!
  !! \param nats Number of atoms.
  !! \param zmat_bml Congruence transform in bml format.  
  !! \param rho_bml Density matrix.
  !! \param dSx_bml x derivative of S.
  !! \param dSy_bml y derivative of S.
  !! \param dSz_bml z derivative of S.
  !! \param hindex Contains the Hamiltonian indices for every atom (see get_hindex).
  !!
  subroutine prg_get_pulayforce(nats,zmat_bml,ham_bml,rho_bml,&
    dSx_bml,dSy_bml,dSz_bml,hindex,FPUL,threshold)
    implicit none
    real(dp), allocatable, intent(inout) :: FPUL(:,:)
    integer, intent(in) :: nats
    type(bml_matrix_t), intent(in)  ::  dSx_bml, dSy_bml, dSz_bml
    type(bml_matrix_t), intent(in)  ::  rho_bml, ham_bml, zmat_bml
    integer, intent(in)                ::  hindex(:,:)
    integer :: I_A, I_B, i , j, norb
    real(dp), intent(in) :: threshold
    type(bml_matrix_t)  :: Xtmp_bml, Ytmp_bml, Ztmp_bml, SIHD_bml
    type(bml_matrix_t)  :: aux_bml
    real(dp), allocatable :: diagxtmp(:), diagytmp(:), diagztmp(:)
    real(dp) :: partrace

    write(*,*)"In prg_get_pulayforce ..."

    if(.not.allocated(FPUL))then 
      allocate(FPUL(3,nats))
    endif 

    FPUL = 0.0_dp
     
    norb = bml_get_N(rho_bml)

    !SIHD = 2*Z*Z'*H*D;  
    call bml_copy_new(zmat_bml,SIHD_bml) 
    call bml_copy_new(zmat_bml,aux_bml) 
    call bml_transpose(zmat_bml,aux_bml) !Z'
    call bml_multiply(zmat_bml,aux_bml,SIHD_bml,2.0_dp,0.0_dp,threshold) !2*Z*Z'
    call bml_multiply(SIHD_bml,ham_bml,aux_bml,1.0_dp,0.0_dp,threshold) !2*Z*Z'*H
    call bml_multiply(aux_bml,rho_bml,SIHD_bml,1.0_dp,0.0_dp,threshold) !2*Z*Z'*D
    call bml_deallocate(aux_bml)

    call bml_copy_new(rho_bml,Xtmp_bml)
    allocate(diagxtmp(norb))
    call bml_multiply(dSx_bml,SIHD_bml,Xtmp_bml,1.0_dp,0.0_dp,threshold)
    call bml_get_diagonal(Xtmp_bml,diagxtmp)
    call bml_deallocate(Xtmp_bml)

    call bml_copy_new(rho_bml,Ytmp_bml)
    allocate(diagytmp(norb))
    call bml_multiply(dSy_bml,SIHD_bml,Ytmp_bml,1.0_dp,0.0_dp,threshold)
    call bml_get_diagonal(Ytmp_bml,diagytmp)
    call bml_deallocate(Ytmp_bml)

    call bml_copy_new(rho_bml,Ztmp_bml)
    allocate(diagztmp(norb))
    call bml_multiply(dSz_bml,SIHD_bml,Ztmp_bml,1.0_dp,0.0_dp,threshold)
    call bml_get_diagonal(Ztmp_bml,diagztmp)
    call bml_deallocate(Ztmp_bml)
    
    call bml_deallocate(SIHD_bml)

    !$omp parallel do default(none) private(i) &
    !$omp private(I_A,I_B,j,partrace) &
    !$omp shared(hindex,diagxtmp,diagytmp,diagztmp,FPUL,nats)
    do I = 1,nats
      I_A = hindex(1,I);
      I_B = hindex(2,I);

      partrace = 0.0_dp
      do j=I_A,I_B
        partrace = partrace + diagxtmp(j)
      enddo
      FPUL(1,I) = partrace;

      partrace = 0.0_dp
      do j=I_A,I_B
        partrace = partrace + diagytmp(j)
      enddo
      FPUL(2,I) = partrace;

      partrace = 0.0_dp
      do j=I_A,I_B
        partrace = partrace + diagztmp(j)
      enddo
      FPUL(3,I) = partrace;      

    enddo
    !$omp end parallel do
    
    deallocate(diagxtmp)
    deallocate(diagytmp)
    deallocate(diagztmp)

  end subroutine prg_get_pulayforce

end module prg_PulayComponent_mod
