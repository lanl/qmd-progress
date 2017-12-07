!> Module to prg_orthogonalize and prg_deorthogonalize any operator.
!! \ingroup PROGRESS
!! Typically the Hamiltonin needs to be prg_orthogonalized:
!! \f$ H_{\mathrm{ortho}} = Z^{\dagger} H Z \f$
!!
!! Also, if the density matrix was obtained from the prg_orthogonalized Hamiltonian,
!! it can be prg_deorthogonalized as:
!! \f$ \rho = Z \rho_{\mathrm{ortho}} Z^{\dagger} \f$
!!
module prg_nonortho_mod

  use bml
  use prg_parallel_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  public :: prg_orthogonalize, prg_deorthogonalize

contains

  !> This routine performs:
  !! \f$ A_{ortho} = Z^{\dagger} A Z \f$
  !!
  !! \param A_bml Matrix to be prg_orthogonalized in bml format.
  !! \param zmat_bml Congruence transform to be used.
  !! \param orthoA_bml Matrix resulting from the orthogonalization.
  !! \param threshold Threshold value to be used in the matrix-matrix operations.
  !! \param bml_type bml format to be used.
  !! \param verbose Verbosity level.
  !!
  subroutine prg_orthogonalize(A_bml,zmat_bml,orthoA_bml,threshold,bml_type,verbose)
    implicit none
    integer :: hdim, mdim
    integer, intent(in) :: verbose
    real(dp), intent(in) :: threshold
    type(bml_matrix_t), intent(inout) :: A_bml
    type(bml_matrix_t), intent(inout) :: zmat_bml
    type(bml_matrix_t) :: AUX_bml
    type(bml_matrix_t), intent(inout) :: OrthoA_bml
    character(len=*), intent(in) :: bml_type

    if(verbose.EQ.1) write(*,*)"In prg_orthogonalize ..."

    hdim= bml_get_N(A_bml)
    mdim= bml_get_M(A_bml)

    !Allocate bml's
    if(bml_get_N(orthoA_bml) .le. 0)then
       call bml_zero_matrix(bml_type,BML_ELEMENT_REAL,dp,HDIM ,MDIM,OrthoA_bml, &
            bml_get_distribution_mode(A_bml))
    endif

    !Do the operations in bml
    call bml_transpose(zmat_bml, aux_bml)

    call bml_multiply(aux_bml, A_bml, OrthoA_bml, 1.0_dp, 0.0_dp,threshold) !Z^t*A

    call bml_multiply(OrthoA_bml, zmat_bml, aux_bml, 1.0_dp, 0.0_dp,threshold) !Z^t*A * Z

    call bml_copy_new(aux_bml, OrthoA_bml)

    call bml_deallocate(aux_bml)

  end subroutine prg_orthogonalize


  !> This routine performs:
  !! \f$ A = Z A_{ortho}  Z^{\dagger} \f$
  !!
  !! \param orthoA_bml Matrix to be prg_deorthogonalized.
  !! \param zmat_bml Congruence transform to be used.
  !! \param A_bml Matrix resulting from the prg_deorthogonalized in bml format.
  !! \param threshold Threshold value to be used in the matrix-matrix operations.
  !! \param bml_type bml format to be used.
  !! \param verbose Verbosity level.
  !!
  subroutine prg_deorthogonalize(orthoA_bml,zmat_bml,a_bml,threshold,bml_type,verbose)
    implicit none
    integer :: HDIM,verbose
    real(dp) :: threshold
    type(bml_matrix_t), intent(inout) :: a_bml
    type(bml_matrix_t), intent(in) :: zmat_bml
    type(bml_matrix_t) :: aux_bml
    type(bml_matrix_t), intent(in) :: orthoA_bml
    character(len=*) :: bml_type

    if(verbose.EQ.1) write(*,*)"In prg_deorthogonalize ..."

    HDIM = bml_get_N(orthoA_bml)

    !Allocate bml's
    if(bml_get_N(a_bml).LE.0) then
      call bml_zero_matrix(bml_type,BML_ELEMENT_REAL,dp,HDIM,HDIM,a_bml, &
        bml_get_distribution_mode(orthoA_bml))
    endif

    call bml_transpose(zmat_bml, aux_bml)

    call bml_multiply(orthoA_bml, aux_bml, a_bml, 1.0_dp, 0.0_dp, threshold) !orthoA*Z^t

! Required when running distributed
#ifdef DO_MPI
     if (getNRanks() > 1 .and. &
         bml_get_distribution_mode(orthoA_bml) == BML_DMODE_DISTRIBUTED) then
         call prg_allGatherParallel(a_bml)
     endif
#endif

    call bml_multiply(zmat_bml, a_bml, aux_bml, 1.0_dp, 0.0_dp, threshold) !Z*orthoA * Z^t

!    call bml_copy(aux_bml, a_bml)
     call bml_copy_new(aux_bml, a_bml)

    call bml_deallocate(aux_bml)

  end subroutine prg_deorthogonalize

end module prg_nonortho_mod
