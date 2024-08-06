!> Module to compute the density matrix response and related quantities.
!! \ingroup GPMDK
!! \todo Add the response scf
!! \todo Change name response_SP2 to dm_prt_response
!! \todo Change name response_rs to rs_prt_response
!! \brief More information about the theory can be found at \cite Niklasson2005 and Niklasson2015
module gpmdcov_response_mod

  use bml
  use prg_kernelparser_mod
  use prg_densitymatrix_mod
  use prg_openfiles_mod
  use gpmdcov_nvtx_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.14159265358979323846264338327950_dp

  type, public :: RespData_type
    character(20) :: RespMode
    character(20) :: TypeOfPert
    character(20) :: BmlType
    integer :: Mdim
    real(dp) :: NumThresh
    logical :: ComputeDipole
    logical :: GetResponse
    real(dp) :: FieldIntensity
    real(dp) :: field(3)
  end type RespData_type

  public :: gpmdcov_response_dpdmu

contains

  !>  First-order Canonical Density Matrix Perturbation Theory.
  !! \brief  Assuming known mu0 and representation in H0's eigenbasis Q,
  !! where H0 and P0 become diagonal.
  !! (mu0, eigenvalues e and eigenvectors Q of H0 are assumed known)
  !! Based on PRL 92, 193001 (2004) and PRE 92, 063301 (2015)
  subroutine gpmdcov_response_dpdmu(P1_bml,dPdMu,H1_bml,Norbs,beta,Q_bml,evals,mu0,m,HDIM)

    type(bml_matrix_t), intent(inout)    ::  P1_bml, H1_bml
    type(bml_matrix_t)    ::  dx1_bml
    type(bml_matrix_t) :: Q_bml
    real(dp), parameter   :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0
    integer, intent(in)     :: HDIM, m, Norbs ! = Number of occupied orbitals,m= Number of recursion steps
    real(dp), intent(in)  :: evals(HDIM) !Q and e are eigenvectors and eigenvalues of H0
    real(dp), intent(inout)  :: dPdMu(HDIM) !Q and e are eigenvectors and eigenvalues of H0
    !    real(PREC)              :: P1(HDIM,HDIM) ! Density matrix response
    !    derivative with respect to perturbation H1 to H0
    real(dp), intent(in)  :: beta, mu0 ! Electronic temperature and chemicalpotential
    real(dp), allocatable :: P1(:,:), DX1(:,:)
    real(dp)              :: h_0(HDIM), p_0(HDIM),p_02(HDIM),iD0(HDIM)
    real(dp)              :: cnst, kB, mu1, sumdPdmu
    integer                 :: i, j, k
    character(20) :: bml_type

    kB = 8.61739e-5        ! (eV/K)
    h_0 = evals                ! Diagonal Hamiltonian H0 respresented in the eigenbasis Q
    cnst = beta/(1.D0*2**(m+2)) ! Scaling constant
    p_0 = 0.5D0 + cnst*(h_0-mu0)  ! Initialization for P0 represented in eigenbasis Q
#ifdef USE_NVTX
    call nvtxStartRange("bml_copy_new",1)
#endif
    call bml_copy_new(H1_bml,P1_bml)
    call bml_copy_new(H1_bml,DX1_bml)
#ifdef USE_NVTX
    call nvtxEndRange
#endif
    call bml_scale(-cnst,P1_bml)    !(set mu1 = 0 for simplicity) !Initialization of DM response in Q representation (not diagonal in Q)
    allocate(P1(HDIM,HDIM))
    allocate(DX1(HDIM,HDIM))

    P1 = 0.0_dp
    DX1 = 0.0_dp
#ifdef USE_NVTX
    call nvtxStartRange("bml_export_to_dense",2)
#endif

    call bml_export_to_dense(P1_bml,P1)
    call bml_export_to_dense(DX1_bml,DX1)
#ifdef USE_NVTX
    call nvtxEndRange
#endif
    
!!$acc data copy(P1(1:HDIM,1:HDIM),DX1(1:HDIM,1:HDIM),p_0(1:HDIM))
#ifdef USE_OFFLOAD
    !$omp target enter data map(alloc:P1(1:HDIM,1:HDIM),p_0(1:HDIM),HDIM,j,k)
    !$omp target update to(P1(1:HDIM,1:HDIM),p_0(1:HDIM),HDIM)
#endif
#ifdef USE_NVTX
    call nvtxStartRange("Response Kernel",3)
#endif
    do i = 1,m  ! Loop over m recursion steps
       !p_02 = p_0*p_0
#ifdef USE_OFFLOAD
       !$omp target teams distribute default(none) &
     !!$acc parallel loop deviceptr(P1,DX1,p_0)
     !!$acc parallel loop
     !$omp shared(HDIM,P1,p_0)
       do k = 1,HDIM
          !$omp parallel do
          do j = 1,HDIM
             P1(j,k) = 1.D0/(2.D0*p_0(j)*(p_0(j)-1.D0)+1.D0)*((p_0(j) + p_0(k))*P1(j,k) + 2.D0*(P1(j,k)-(p_0(j) + p_0(k))*P1(j,k))*1.D0/(2.D0*p_0(k)*(p_0(k)-1.0D0)+1.D0)*p_0(k)*p_0(k))
          enddo
          !$omp end parallel do
        !P1(:,k) = 1.D0/(2.D0*(p_0(:)*p_0(:)-p_0(:))+1.D0)*(DX1(:,k) + 2.D0*(P1(:,k)-DX1(:,k))*1.D0/(2.D0*(p_0(k)*p_0(k)-p_0(k))+1.D0)*p_0(k)*p_0(k))
        !P1(k,:) = iD0(k)*(DX1(k,:) + 2.D0*(P1(k,:)-DX1(k,:))*p_0(:))
     enddo
     !!$acc end parallel loop
     !!$acc kernels deviceptr(p_0)
     !!$acc kernels
     !$omp end target teams distribute
     !$omp target
     p_0 = 1.D0/(2.D0*(p_0(:)*p_0(:)-p_0(:))+1.D0)*p_0(:)*p_0(:)
     !$omp end target
#else
       !$omp parallel do default(none) &
       !$omp private(k,j) &
     !$omp shared(HDIM,P1,p_0)
     do k = 1,HDIM
        !$omp simd
        do j = 1,HDIM
           !P1(:,k) = 1.D0/(2.D0*p_0(:)*(p_0(:)-1.D0)+1.D0)*((p_0(:) + p_0(k))*P1(:,k) + 2.D0*(P1(:,k)-(p_0(:) + p_0(k))*P1(:,k))*1.D0/(2.D0*p_0(k)*(p_0(k)-1.0D0)+1.D0)*p_0(k)*p_0(k))
           P1(j,k) = 1.D0/(2.D0*p_0(j)*(p_0(j)-1.D0)+1.D0)*((p_0(j) + p_0(k))*P1(j,k) + 2.D0*(P1(j,k)-(p_0(j) + p_0(k))*P1(j,k))*1.D0/(2.D0*p_0(k)*(p_0(k)-1.0D0)+1.D0)*p_0(k)*p_0(k))
        enddo
        !$omp end simd
        !P1(:,k) = 1.D0/(2.D0*(p_0(:)*p_0(:)-p_0(:))+1.D0)*(DX1(:,k) + 2.D0*(P1(:,k)-DX1(:,k))*1.D0/(2.D0*(p_0(k)*p_0(k)-p_0(k))+1.D0)*p_0(k)*p_0(k))
        !P1(k,:) = iD0(k)*(DX1(k,:) + 2.D0*(P1(k,:)-DX1(k,:))*p_0(:))
     enddo
     !$omp end parallel do
     p_0 = 1.D0/(2.D0*(p_0(:)*p_0(:)-p_0(:))+1.D0)*p_0(:)*p_0(:)
#endif
    enddo
#ifdef USE_NVTX
    call nvtxEndRange
#endif
!!$acc end data
#ifdef USE_OFFLOAD
    !$omp target update from(P1(1:HDIM,1:HDIM))
    !$omp target exit data map(delete:P1(1:HDIM,1:HDIM),DX1(1:HDIM,1:HDIM),p_0(1:HDIM),HDIM,j,k)
#endif
    
    bml_type = bml_get_type(P1_bml)
    call bml_import_from_dense(bml_type,P1,P1_bml,ZERO,HDIM) !Dense to dense_bml

    deallocate(P1)
    deallocate(DX1)

    call bml_deallocate(DX1_bml)

    dPdmu = beta*p_0*(1.D0-p_0)

  end subroutine gpmdcov_response_dpdmu

end module gpmdcov_response_mod
