module prg_c_interface

  use iso_c_binding

  use prg_progress_mod
  use prg_densitymatrix_mod
  use bml_types_m

  implicit none

  public :: prg_version_c
  public :: prg_progress_init_c
  public :: prg_progress_shutdown_c

  ! density matrix module
  public :: prg_build_density_T0_c
  public :: prg_check_idempotency_c, prg_get_eigenvalues_c
  public :: prg_get_flevel_c, prg_build_density_T_c, prg_build_atomic_density_c
  public :: prg_build_density_T_Fermi_c, prg_get_flevel_nt_c, prg_build_density_T_fulldata_c
  public :: prg_build_density_T_ed_c, prg_toEigenspace_c, prg_toCanonicalspace_c
  public :: Canon_DM_PRT_c, prg_get_evalsDvalsEvects_c, prg_build_density_fromEvalsAndEvects_c

contains

  ! C wrapper subroutine
  subroutine prg_version_c() bind(C, name="prg_version")
     call prg_version()
  end subroutine prg_version_c

  subroutine prg_progress_init_c() bind(C, name="prg_progress_init")
     call prg_progress_init()
  end subroutine prg_progress_init_c

  subroutine prg_progress_shutdown_c() bind(C, name="prg_progress_shutdown")
     call prg_progress_shutdown()
  end subroutine prg_progress_shutdown_c


!------------------------------------------------
!  Beginning of prg_densitymatrix_mod
!------------------------------------------------

subroutine prg_build_density_T0_c(norbs, ham_bml_c, rho_bml_c, threshold, bndfil, eigenvalues_out) bind(C, name="prg_build_density_T0")
    integer(c_int), value :: norbs
    real(c_double), intent(in), value :: threshold, bndfil
    type(c_ptr), value :: ham_bml_c
    type(c_ptr), value :: rho_bml_c
    type(bml_matrix_t) :: ham_bml
    type(bml_matrix_t) :: rho_bml
    real(c_double) :: eigenvalues_out(norbs)
    real(c_double), allocatable :: eigenvalues(:)

    ham_bml%ptr = ham_bml_c
    rho_bml%ptr = rho_bml_c

    call prg_build_density_T0(ham_bml, rho_bml, threshold, bndfil, eigenvalues)
    eigenvalues_out = eigenvalues
    deallocate(eigenvalues)

end subroutine prg_build_density_T0_c

subroutine prg_build_density_T_c(norbs, ham_bml_c, rho_bml_c, threshold, bndfil, kbt, ef, eigenvalues_out) bind(C, name="prg_build_density_T")
    integer(c_int), value :: norbs
    real(c_double), value :: bndfil
    real(c_double), value :: threshold
    real(c_double), value :: kbt
    real(c_double), value :: ef
    real(c_double) :: eigenvalues_out(norbs)
    real(c_double), allocatable :: eigenvalues(:)
    type(c_ptr), value :: ham_bml_c
    type(bml_matrix_t) :: ham_bml
    type(c_ptr), value :: rho_bml_c
    type(bml_matrix_t) :: rho_bml
    ham_bml%ptr = ham_bml_c
    rho_bml%ptr = rho_bml_c
    call prg_build_density_T(ham_bml, rho_bml, threshold, bndfil, kbt, ef, eigenvalues)
    eigenvalues_out = eigenvalues
    deallocate(eigenvalues)
end subroutine prg_build_density_T_c

subroutine prg_build_density_T_fulldata_c(norbs, ham_bml_c, rho_bml_c, threshold, bndfil, kbt, ef, eigenvalues_out, evects_bml_c, fvals_out) bind(C, name="prg_build_density_T_fulldata")
    integer(c_int), value :: norbs
    real(c_double), value :: bndfil
    real(c_double), value :: threshold
    real(c_double), value :: kbt
    real(c_double), value :: ef
    real(c_double) :: eigenvalues_out(norbs)
    real(c_double) :: fvals_out(norbs)

    real(c_double), allocatable :: fvals(:)
    real(c_double), allocatable :: eigenvalues(:)

    type(c_ptr), value :: ham_bml_c
    type(c_ptr), value :: rho_bml_c
    type(c_ptr), value :: evects_bml_c
    type(bml_matrix_t) :: ham_bml
    type(bml_matrix_t) :: rho_bml
    type(bml_matrix_t) :: evects_bml

    ham_bml%ptr = ham_bml_c
    rho_bml%ptr = rho_bml_c
    evects_bml%ptr = evects_bml_c
    call prg_build_density_T_fulldata(ham_bml, rho_bml, threshold, bndfil, kbt, ef, eigenvalues, evects_bml, fvals)
    fvals_out = fvals
    eigenvalues_out = eigenvalues
    deallocate(eigenvalues)
    deallocate(fvals)

end subroutine prg_build_density_T_fulldata_c

subroutine prg_build_density_T_ed_c(norbs, ham_bml_c, rho_bml_c, evects_bml_c, threshold, bndfil, kbt, ef, evals_out, dvals_out, hindex_out, llsize, verbose) bind(C, name="prg_build_density_T_ed")
    integer(c_int), value :: norbs
    integer(c_int), value :: llsize
    integer(c_int), value :: verbose
    real(c_double), value :: bndfil
    real(c_double), value :: threshold
    real(c_double), value :: kbt
    real(c_double), value :: ef
    real(c_double) :: evals_out(norbs)
    real(c_double) :: dvals_out(norbs)
    integer(c_int) :: hindex_out(norbs, norbs)

    integer(c_int), allocatable  :: hindex(:,:)
    real(c_double), allocatable  :: evals(:)
    real(c_double), allocatable  :: dvals(:)

    type(c_ptr), value :: ham_bml_c
    type(bml_matrix_t) :: ham_bml
    type(c_ptr), value :: rho_bml_c
    type(bml_matrix_t) :: rho_bml
    type(c_ptr), value :: evects_bml_c
    type(bml_matrix_t) :: evects_bml
    ham_bml%ptr = ham_bml_c
    rho_bml%ptr = rho_bml_c
    evects_bml%ptr = evects_bml_c

    call prg_build_density_T_ed(ham_bml, rho_bml, evects_bml, threshold, bndfil, kbt, ef, evals, dvals, hindex, llsize, verbose)

    evals_out = evals
    dvals_out = dvals
    hindex_out = hindex
    deallocate(evals)
    deallocate(dvals)
    deallocate(hindex)

end subroutine prg_build_density_T_ed_c


subroutine prg_get_evalsDvalsEvects_c(norbs, ham_bml_c, threshold, hindex_in, llsize, evals_out, dvals_out, evects_bml_c, verbose) bind(C, name="prg_get_evalsDvalsEvects")
    integer(c_int), value :: norbs, llsize, verbose
    real(c_double), value :: threshold
    integer(c_int)  :: hindex_in(norbs,norbs)
    real(c_double)  :: evals_out(norbs)
    real(c_double)  :: dvals_out(norbs)
    integer(c_int), allocatable  :: hindex(:,:)
    real(c_double), allocatable  :: evals(:)
    real(c_double), allocatable  :: dvals(:)
    type(bml_matrix_t) :: ham_bml
    type(bml_matrix_t) :: evects_bml
    type(c_ptr), value :: ham_bml_c
    type(c_ptr), value :: evects_bml_c

    ham_bml%ptr = ham_bml_c
    evects_bml%ptr = evects_bml_c
    allocate(hindex(2, llsize))
    hindex = hindex_in
    call prg_get_evalsDvalsEvects(ham_bml, threshold, hindex, llsize, evals, dvals, evects_bml, verbose)
    evals_out = evals
    dvals_out = dvals
    deallocate(hindex, evals, dvals)

end subroutine prg_get_evalsDvalsEvects_c


subroutine prg_build_density_fromEvalsAndEvects_c(norbs,evects_bml_c, evals, rho_bml_c, threshold, bndfil, kbt, ef, verbose) bind(C, name="prg_build_density_fromEvalsAndEvects")
    integer(c_int), value :: norbs, verbose
    real(c_double), value :: threshold
    real(c_double), value :: bndfil
    real(c_double), value :: kbt
    real(c_double), value :: ef
    real(c_double) :: evals(norbs)
    type(c_ptr), value :: evects_bml_c
    type(bml_matrix_t) :: evects_bml
    type(c_ptr), value :: rho_bml_c
    type(bml_matrix_t) :: rho_bml

    evects_bml%ptr = evects_bml_c
    rho_bml%ptr = rho_bml_c
    call prg_build_density_fromEvalsAndEvects(evects_bml, evals, rho_bml, threshold, bndfil, kbt, ef, verbose)

end subroutine prg_build_density_fromEvalsAndEvects_c

subroutine prg_build_density_T_Fermi_c(ham_bml_c, rho_bml_c, threshold, kbt, ef, verbose, drho) bind(C, name="prg_build_density_T_Fermi")
    integer(c_int), value  :: verbose
    real(c_double), value :: threshold
    real(c_double), value :: kbt
    real(c_double), value :: ef
    type(c_ptr), value :: ham_bml_c
    type(bml_matrix_t) :: ham_bml
    type(c_ptr), value :: rho_bml_c
    type(bml_matrix_t) :: rho_bml
    real(c_double), value :: drho
    ham_bml%ptr = ham_bml_c
    rho_bml%ptr = rho_bml_c
    call prg_build_density_T_Fermi(ham_bml, rho_bml, threshold, kbt, ef, verbose, drho)
end subroutine prg_build_density_T_Fermi_c


subroutine prg_build_atomic_density_c(rhoat_bml_c, numel, hindex, spindex, norb, bml_type) bind(C, name="prg_build_atomic_density")
    character(c_char), value :: bml_type
    integer(c_int), target :: hindex(norb, norb)
    integer(c_int), value :: norb
    integer(c_int), target :: spindex(norb)
    real(c_double), target :: numel(norb)
    type(c_ptr), value :: rhoat_bml_c
    type(bml_matrix_t) :: rhoat_bml
    rhoat_bml%ptr = rhoat_bml_c
    call prg_build_atomic_density(rhoat_bml, numel, hindex, spindex, norb, bml_type)
end subroutine prg_build_atomic_density_c

subroutine prg_get_flevel_c(norbs, eigenvalues, kbt, bndfil, tol, Ef) bind(C, name="prg_get_flevel")
    integer(c_int) :: norbs
    real(c_double), value :: tol
    real(c_double), value :: bndfil
    real(c_double), target :: eigenvalues(norbs)
    real(c_double), value :: kbt
    real(c_double), value :: Ef
    logical :: err
    call prg_get_flevel(eigenvalues, kbt, bndfil, tol, Ef, err)
end subroutine prg_get_flevel_c

subroutine prg_get_flevel_nt_c(norbs, eigenvalues, kbt, bndfil, tol, ef, verbose) bind(C, name="prg_get_flevel_nt")
    integer(c_int) :: norbs
    real(c_double), value :: bndfil
    real(c_double), value :: kbt
    real(c_double), value :: tol
    real(c_double), target :: eigenvalues(norbs)
    real(c_double), value :: ef
    integer(c_int), optional  :: verbose
    logical :: err
    call prg_get_flevel_nt(eigenvalues, kbt, bndfil, tol, ef, err, verbose)
end subroutine prg_get_flevel_nt_c

subroutine prg_get_eigenvalues_c(norbs, ham_bml_c, eigenvalues_out, verbose) bind(C, name="prg_get_eigenvalues")
    integer(c_int) :: norbs
    integer(c_int), value :: verbose
    real(c_double)     :: eigenvalues_out(norbs)
    real(c_double), allocatable :: eigenvalues(:)
    type(c_ptr), value :: ham_bml_c
    type(bml_matrix_t) :: ham_bml
    ham_bml%ptr = ham_bml_c
    call prg_get_eigenvalues(ham_bml, eigenvalues, verbose)
    eigenvalues_out = eigenvalues
    deallocate(eigenvalues)
end subroutine prg_get_eigenvalues_c

subroutine prg_check_idempotency_c(mat_bml_c, threshold, idempotency) bind(C, name="prg_check_idempotency")
    real(c_double), value :: threshold
    real(c_double), value :: idempotency
    type(c_ptr), value :: mat_bml_c
    type(bml_matrix_t) :: mat_bml
    mat_bml%ptr = mat_bml_c

    call prg_check_idempotency(mat_bml, threshold, idempotency)

end subroutine prg_check_idempotency_c

!
subroutine prg_toEigenspace_c(mat_bml_c, matEig_bml_c, evects_bml_c, threshold, verbose) bind(C, name="prg_toEigenspace")
    integer(c_int), optional  :: verbose
    type(c_ptr), value :: mat_bml_c
    type(bml_matrix_t) :: mat_bml
    type(c_ptr), value :: evects_bml_c
    type(bml_matrix_t) :: evects_bml
    type(c_ptr), value :: matEig_bml_c
    type(bml_matrix_t) :: matEig_bml
    real(c_double), value :: threshold
    mat_bml%ptr = mat_bml_c
    evects_bml%ptr = evects_bml_c
    matEig_bml%ptr = matEig_bml_c
    call prg_toEigenspace(mat_bml, matEig_bml, evects_bml, threshold, verbose)
end subroutine prg_toEigenspace_c

subroutine prg_toCanonicalspace_c(mat_bml_c, matCan_bml_c, evects_bml_c, threshold, verbose) bind(C, name="prg_toCanonicalspace")
    integer(c_int), optional  :: verbose
    type(c_ptr), value :: mat_bml_c
    type(bml_matrix_t) :: mat_bml
    type(c_ptr), value :: evects_bml_c
    type(bml_matrix_t) :: evects_bml
    type(c_ptr), value :: matCan_bml_c
    type(bml_matrix_t) :: matCan_bml
    real(c_double), value :: threshold
    mat_bml%ptr = mat_bml_c
    evects_bml%ptr = evects_bml_c
    matCan_bml%ptr = matCan_bml_c
    call prg_toCanonicalspace(mat_bml, matCan_bml, evects_bml, threshold, verbose)
end subroutine prg_toCanonicalspace_c

subroutine Canon_DM_PRT_c(P1, H1, Nocc, T, Q, e, mu0, m, HDIM) bind(C, name="Canon_DM_PRT")
    integer(c_int), value :: HDIM
    integer(c_int), value :: m
    real(c_double), intent(in) :: H1(HDIM, HDIM), Q(HDIM, HDIM), e(HDIM)
    real(c_double), intent(out) :: P1(HDIM, HDIM)
    real(c_double), value :: T, mu0, Nocc

    call Canon_DM_PRT(P1, H1, Nocc, T, Q, e, mu0, m, HDIM)
end subroutine Canon_DM_PRT_c


!------------------------------------------------
!  end of prg_densitymatrix_mod
!------------------------------------------------


!------------------------------------------------
!  Beginning of prg_charges_mod
!------------------------------------------------



!------------------------------------------------
!  End of prg_charges_mod
!------------------------------------------------



end module prg_c_interface


