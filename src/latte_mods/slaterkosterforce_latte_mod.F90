!> Slater-Koster force module.
!! \ingroup LATTE
!! \brief Module to compute the Slater-Koster contribution to the force.
!!
module slaterkosterforce_latte_mod

  use bml

#ifdef USE_NVTX
  use prg_nvtx_mod
#endif
  
  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  public :: get_skforce

contains

  !> Gets the SK contribution to the force.
  !! \note This is computed from: from \f$ Tr[\rho \frac{dH0}{dR}] \f$
  !! \param Nr_atoms Number of atoms.
  !! \param rho_bml Density matrix.
  !! \param dH0x_bml x derivative of H0.
  !! \param dH0y_bml y derivative of H0.
  !! \param dH0z_bml z derivative of H0.
  !! \param hindex Contains the Hamiltonian indices for every atom (see get_hindex).
  !!
  subroutine get_skforce(Nr_atoms,rho_bml,dH0x_bml,dH0y_bml,&
       dH0z_bml,hindex,SKForce,threshold)

    implicit none

    integer, intent(in) :: Nr_atoms
    type(bml_matrix_t), intent(in)  ::  dH0x_bml, dH0y_bml, dH0z_bml
    type(bml_matrix_t), intent(in)  ::  rho_bml
    integer, intent(in)                ::  hindex(:,:)
    real(dp), allocatable, intent(inout) :: SKForce(:,:)
    integer :: I_A, I_B, norbs, i, j, norb
    real(dp), intent(in) :: threshold
    type(bml_matrix_t)  :: Xtmp_bml, Ytmp_bml, Ztmp_bml
    real(dp), allocatable :: diagxtmp(:), diagytmp(:), diagztmp(:)
    real(dp) :: partrace,partracex,partracey,partracez
#ifdef USE_OFFLOAD
    type(c_ptr)                        :: Xtmp_bml_c_ptr, Ytmp_bml_c_ptr, Ztmp_bml_c_ptr
    integer :: ld
    real(c_double), pointer            :: Xtmp_bml_ptr(:,:), Ytmp_bml_ptr(:,:), Ztmp_bml_ptr(:,:)
#endif

    write(*,*)"In get_skforce ..."

    if(.not.allocated(SKForce))then
      allocate(SKForce(3,Nr_atoms))
    endif

    SKForce = 0.0_dp

    norb = bml_get_N(rho_bml)

    ! Slater-Koster Force SKForce from Tr[D*dH0/dR]
    call bml_copy_new(rho_bml,Xtmp_bml)
    call bml_multiply(dH0x_bml,rho_bml,Xtmp_bml,1.0_dp,0.0_dp,threshold)

    call bml_copy_new(rho_bml,Ytmp_bml)
    call bml_multiply(dH0y_bml,rho_bml,Ytmp_bml,1.0_dp,0.0_dp,threshold)

    call bml_copy_new(rho_bml,Ztmp_bml)
    call bml_multiply(dH0z_bml,rho_bml,Ztmp_bml,1.0_dp,0.0_dp,threshold)
    
#ifdef USE_OFFLOAD
    Xtmp_bml_c_ptr = bml_get_data_ptr_dense(Xtmp_bml)
    Ytmp_bml_c_ptr = bml_get_data_ptr_dense(Ytmp_bml)
    Ztmp_bml_c_ptr = bml_get_data_ptr_dense(Ztmp_bml)
    ld = bml_get_ld_dense(Xtmp_bml)
    call c_f_pointer(Xtmp_bml_c_ptr,Xtmp_bml_ptr,shape=[ld,norb])
    call c_f_pointer(Ytmp_bml_c_ptr,Ytmp_bml_ptr,shape=[ld,norb])
    call c_f_pointer(Ztmp_bml_c_ptr,Ztmp_bml_ptr,shape=[ld,norb])
    !$acc enter data copyin(SKForce(1:3,1:Nr_atoms),hindex(1:2,1:Nr_atoms))
    !$acc parallel loop deviceptr(Xtmp_bml_ptr,Ytmp_bml_ptr,Ztmp_bml_ptr) &
    !$acc private(i,j,partracex,partracey,partracez)
    do I = 1,Nr_atoms

      partracex = 0.0_dp
      partracey = 0.0_dp
      partracez = 0.0_dp
      !$acc loop reduction(partracex,partracey,partracez)
      do j=hindex(1,i),hindex(2,i)
        partracex = partracex + Xtmp_bml_ptr(j,j)
        partracey = partracey + Ytmp_bml_ptr(j,j)
        partracez = partracez + Ztmp_bml_ptr(j,j)
     enddo
     !$acc end loop
      SKForce(1,I) = -2.0_dp*partracex;
      SKForce(2,I) = -2.0_dp*partracey;
      SKForce(3,I) = -2.0_dp*partracez;

    enddo
    !$acc end loop
    !$acc exit data copyout(SKForce(1:3,1:Nr_atoms)) delete(hindex(1:2,1:Nr_atoms))
#else
    allocate(diagxtmp(norb))
    call bml_get_diagonal(Xtmp_bml,diagxtmp)
    call bml_deallocate(Xtmp_bml)
    allocate(diagytmp(norb))
    call bml_get_diagonal(Ytmp_bml,diagytmp)
    call bml_deallocate(Ytmp_bml)
    allocate(diagztmp(norb))
    call bml_get_diagonal(Ztmp_bml,diagztmp)
    call bml_deallocate(Ztmp_bml)

    !$omp parallel do default(none) private(i) &
    !$omp private(I_A,I_B,j,partrace) &
    !$omp shared(hindex,diagxtmp,diagytmp,diagztmp,SKForce,Nr_atoms)
    do I = 1,Nr_atoms
      I_A = hindex(1,I);
      I_B = hindex(2,I);

      partrace = 0.0_dp
      do j=I_A,I_B
        partrace = partrace + diagxtmp(j)
      enddo
      SKForce(1,I) = -2.0_dp*partrace;

      partrace = 0.0_dp
      do j=I_A,I_B
        partrace = partrace + diagytmp(j)
      enddo
      SKForce(2,I) = -2.0_dp*partrace;

      partrace = 0.0_dp
      do j=I_A,I_B
        partrace = partrace + diagztmp(j)
      enddo
      SKForce(3,I) = -2.0_dp*partrace;

    enddo
    !$omp end parallel do

    deallocate(diagxtmp)
    deallocate(diagytmp)
    deallocate(diagztmp)
#endif
  end subroutine get_skforce!

end module slaterkosterforce_latte_mod
