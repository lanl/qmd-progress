!> Slater-Koster force module.
!! \ingroup LATTE
!! \brief Module to compute the Slater-Koster contribution to the force.
!!
module slaterkosterforce_latte_mod

  use bml

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
    real(dp) :: partrace

    write(*,*)"In get_skforce ..."

    if(.not.allocated(SKForce))then
      allocate(SKForce(3,Nr_atoms))
    endif

    SKForce = 0.0_dp

    norb = bml_get_N(rho_bml)

    ! Slater-Koster Force SKForce from Tr[D*dH0/dR]
    call bml_copy_new(rho_bml,Xtmp_bml)
    allocate(diagxtmp(norb))
    call bml_multiply(dH0x_bml,rho_bml,Xtmp_bml,1.0_dp,0.0_dp,threshold)
    call bml_get_diagonal(Xtmp_bml,diagxtmp)
    call bml_deallocate(Xtmp_bml)

    call bml_copy_new(rho_bml,Ytmp_bml)
    allocate(diagytmp(norb))
    call bml_multiply(dH0y_bml,rho_bml,Ytmp_bml,1.0_dp,0.0_dp,threshold)
    call bml_get_diagonal(Ytmp_bml,diagytmp)
    call bml_deallocate(Ytmp_bml)

    call bml_copy_new(rho_bml,Ztmp_bml)
    allocate(diagztmp(norb))
    call bml_multiply(dH0z_bml,rho_bml,Ztmp_bml,1.0_dp,0.0_dp,threshold)
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

  end subroutine get_skforce!

end module slaterkosterforce_latte_mod
