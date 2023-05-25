module prg_c_interface

  use, intrinsic :: iso_c_binding

  use prg_progress_mod
  use prg_densitymatrix_mod
  use bml_types_m

  implicit none

  public :: prg_version_c
  public :: prg_progress_init_c
  public :: prg_progress_shutdown_c

  ! density matrix module
  public :: prg_build_density_T0_c
  !public :: prg_check_idempotency_c, prg_get_eigenvalues_c
  !public :: prg_get_flevel_c, prg_build_density_T_c, prg_build_atomic_density_c
  !public :: prg_build_density_T_Fermi_c, prg_get_flevel_nt_c, prg_build_density_T_fulldata_c
  !public :: prg_build_density_T_ed_c, prg_toEigenspace_c, prg_toCanonicalspace_c
  !public :: Canon_DM_PRT_c, prg_get_evalsDvalsEvects_c, prg_build_density_fromEvalsAndEvects_c

contains
  subroutine prg_version_c() bind(C, name="prg_version")
     call prg_version()
  end subroutine prg_version_c

  subroutine prg_progress_init_c() bind(C, name="prg_progress_init")
     call prg_progress_init()
  end subroutine prg_progress_init_c

  subroutine prg_progress_shutdown_c() bind(C, name="prg_progress_shutdown")
     call prg_progress_shutdown()
  end subroutine prg_progress_shutdown_c

  ! C-bound wrapper subroutine
  subroutine prg_build_density_T0_c(ham_bml_c, rho_bml_c, threshold, bndfil) bind(C, name="prg_build_density_T0")
    use iso_c_binding, only : c_double, c_ptr
    type(c_ptr), value :: ham_bml_c
    type(c_ptr), value :: rho_bml_c
    real(c_double), intent(in), value :: threshold, bndfil
    
    type(bml_matrix_t) :: ham_bml
    type(bml_matrix_t) :: rho_bml

    ham_bml%ptr = ham_bml_c
    rho_bml%ptr = rho_bml_c

    call prg_build_density_T0(ham_bml, rho_bml, threshold, bndfil)

  end subroutine prg_build_density_T0_c

end module prg_c_interface

