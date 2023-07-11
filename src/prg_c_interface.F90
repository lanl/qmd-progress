module prg_c_interface

  use iso_c_binding

  use prg_densitymatrix_mod
  use prg_charges_mod
  use prg_chebyshev_mod
  use prg_ewald_mod
  use prg_dos_mod
  use prg_genz_mod
  use prg_normalize_mod
  use prg_system_mod
  use prg_openfiles_mod
  use prg_progress_mod
  use prg_pulaycomponent_mod
  use prg_sp2_fermi_mod
  use prg_sp2_mod
  use prg_timer_mod
  use bml_types_m

  implicit none

  public :: prg_version_c
  public :: prg_progress_init_c
  public :: prg_progress_shutdown_c

  ! density matrix module
  public :: prg_build_density_T0_c
  public :: prg_check_idempotency_c, prg_get_eigenvalues_c
  public :: prg_get_flevel_c, prg_build_density_T_c, prg_build_atomic_density_c
  public :: prg_build_density_T_fermi_c, prg_get_flevel_nt_c, prg_build_density_T_fulldata_c
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

  subroutine prg_build_density_T_fermi_c(ham_bml_c, rho_bml_c, threshold, kbt, ef, verbose, drho) bind(C, name="prg_build_density_T_fermi")
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
    call prg_build_density_T_fermi(ham_bml, rho_bml, threshold, kbt, ef, verbose, drho)
  end subroutine prg_build_density_T_fermi_c

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
  subroutine prg_get_charges_c(nats, norbs, rho_bml_c, over_bml_c, hindex, charges, numel, spindex, mdimin, threshold) bind(C, name="prg_get_charges")
    integer(c_int), value :: nats, norbs
    integer(c_int), value :: mdimin
    integer(c_int) :: hindex(nats,nats)
    integer(c_int) :: spindex(nats)
    real(c_double) :: charges_out(nats)
    real(c_double) :: numel_in(norbs)
    real(c_double), value :: threshold

    real(c_double), allocatable :: charges(:)
    real(c_double), allocatable :: numel(:)

    type(c_ptr), value :: over_bml_c
    type(bml_matrix_t) :: over_bml
    type(c_ptr), value :: rho_bml_c
    type(bml_matrix_t) :: rho_bml
    over_bml%ptr = over_bml_c
    rho_bml%ptr = rho_bml_c
    allocate(numel(norbs))
    numel = numel_in
    call prg_get_charges(rho_bml, over_bml, hindex, charges, numel, spindex, mdimin, threshold)
    charges_out = charges
    deallocate(charges, numel)

  end subroutine prg_get_charges_c

  subroutine prg_get_hscf_c(nats, ham0_bml_c, over_bml_c, ham_bml_c, spindex, hindex, hubbardu, charges, coulomb_pot_r, coulomb_pot_k, mdimin, threshold) bind(C, name="prg_get_hscf")
    integer(c_int), value :: nats
    integer(c_int), target :: hindex(nats, nats)
    integer(c_int), value :: mdimin
    integer(c_int), target :: spindex(nats)
    real(c_double), target :: charges(nats)
    real(c_double), target :: coulomb_pot_r(nats)
    real(c_double), target :: coulomb_pot_k(nats)
    real(c_double), target :: hubbardu(nats)
    real(c_double), value :: threshold
    type(c_ptr), value :: ham0_bml_c
    type(bml_matrix_t) :: ham0_bml
    type(c_ptr), value :: over_bml_c
    type(bml_matrix_t) :: over_bml
    type(c_ptr), value :: ham_bml_c
    type(bml_matrix_t) :: ham_bml
    ham0_bml%ptr = ham0_bml_c
    over_bml%ptr = over_bml_c
    ham_bml%ptr = ham_bml_c
    call prg_get_hscf(ham0_bml, over_bml, ham_bml, spindex, hindex, hubbardu, charges, coulomb_pot_r, coulomb_pot_k, mdimin, threshold)
  end subroutine prg_get_hscf_c

  subroutine prg_get_hscf_v2_c(nats, ham0_bml_c, over_bml_c, ham_bml_c, spindex, hindex, hubbardu, charges, coulomb_pot_r, coulomb_pot_k, mdimin, threshold) bind(C, name="prg_get_hscf_v2")
    integer(c_int), value :: nats
    integer(c_int) :: hindex(nats, nats)
    integer(c_int), value :: mdimin
    real(c_double), value :: threshold
    integer(c_int) :: spindex(nats)
    real(c_double) :: charges(nats)
    real(c_double) :: coulomb_pot_r(nats)
    real(c_double) :: coulomb_pot_k(nats)
    real(c_double) :: hubbardu(nats)
    type(c_ptr), value :: ham0_bml_c
    type(bml_matrix_t) :: ham0_bml
    type(c_ptr), value :: over_bml_c
    type(bml_matrix_t) :: over_bml
    type(c_ptr), value :: ham_bml_c
    type(bml_matrix_t) :: ham_bml

    ham0_bml%ptr = ham0_bml_c
    over_bml%ptr = over_bml_c
    ham_bml%ptr = ham_bml_c
    call prg_get_hscf_v2(ham0_bml, over_bml, ham_bml, spindex, hindex, hubbardu, charges, coulomb_pot_r, coulomb_pot_k, mdimin, threshold)
  end subroutine prg_get_hscf_v2_c

  !------------------------------------------------
  !  End of prg_charges_mod
  !------------------------------------------------

  !------------------------------------------------
  ! prg_chebyshev_mod
  !------------------------------------------------
  ! chebdata is a datatype defined in fortran, need to convert it into c (TODO)

  !subroutine prg_parse_cheb_c(chebdata, filename) bind(C, name="prg_parse_cheb")
  !  character(c_char), value :: filename
  !  call prg_parse_cheb(chebdata, filename)
  !end subroutine prg_parse_cheb_c

  subroutine prg_build_density_cheb_c(ham_bml_c, rho_bml_c, athr, threshold, ncoeffs, kbt, ef, bndfil, jon, verbose) bind(C, name="prg_build_density_cheb")
    integer(c_int), value :: ncoeffs
    integer(c_int), value :: verbose
    real(c_double), value :: athr
    real(c_double), value :: kbt
    real(c_double), value :: bndfil
    real(c_double), value :: threshold
    real(c_double), value :: ef
    type(c_ptr), value :: ham_bml_c
    type(bml_matrix_t) :: ham_bml
    type(c_ptr), value :: rho_bml_c
    type(bml_matrix_t) :: rho_bml
    integer(c_int), value :: jon
    logical :: jon_l = .False.
    if (jon==1) jon_l = .True.

    ham_bml%ptr = ham_bml_c
    rho_bml%ptr = rho_bml_c
    call prg_build_density_cheb(ham_bml, rho_bml, athr, threshold, ncoeffs, kbt, ef, bndfil, jon_l, verbose)
  end subroutine prg_build_density_cheb_c

  subroutine prg_build_density_cheb_fermi_c(ham_bml_c, rho_bml_c, athr, threshold, ncoeffs, kbt, ef, bndfil, getef, fermitol, jon, npts, trkfunc, verbose) bind(C, name="prg_build_density_cheb_fermi")
    type(c_ptr), value :: ham_bml_c
    type(c_ptr), value :: rho_bml_c
    integer(c_int), value :: npts
    integer(c_int), value :: ncoeffs
    integer(c_int), value :: verbose
    real(c_double), value :: fermitol
    real(c_double), value :: athr
    real(c_double), value :: kbt
    real(c_double), value :: bndfil
    real(c_double), value :: threshold
    real(c_double), value :: ef
    type(bml_matrix_t) :: ham_bml
    type(bml_matrix_t) :: rho_bml
    integer(c_int), value :: getef
    integer(c_int), value :: jon
    integer(c_int), value :: trkfunc
    logical :: getef_l = .False., jon_l = .False., trkfunc_l = .False.

    if (jon == 1) jon_l = .True.
    if (getef == 1) getef_l = .True.
    if (trkfunc == 1) trkfunc_l = .True.

    ham_bml%ptr = ham_bml_c
    rho_bml%ptr = rho_bml_c

    call prg_build_density_cheb_fermi(ham_bml, rho_bml, athr, threshold, ncoeffs, kbt, ef, bndfil, getef_l, fermitol, jon_l, npts, trkfunc_l, verbose)

  end subroutine prg_build_density_cheb_fermi_c

  !------------------------------------------------
  ! prg_dos_mod
  !------------------------------------------------

  subroutine prg_write_tdos_c(nstates, eigenvals, gamma, npts, emin, emax, filename) bind(C, name="prg_write_tdos")
    integer(c_int), value :: nstates, npts
    real(c_double), target :: eigenvals(nstates)
    real(c_double), value :: emax
    real(c_double), value :: emin
    real(c_double), value :: gamma
    character(c_char), value :: filename

    call prg_write_tdos(eigenvals, gamma, npts, emin, emax, filename)
  end subroutine prg_write_tdos_c

  !------------------------------------------------
  ! prg_ewald_mod
  !------------------------------------------------

  subroutine Ewald_Real_Space_Single_latte_c(COULOMBV, I, RXYZ, Box, Nr_elem, DELTAQ, J, U, Element_Pointer, Nr_atoms, COULACC, HDIM, Max_Nr_Neigh) bind(C, name="Ewald_Real_Space_Single_latte")
    integer(c_int),    value :: Nr_atoms
    integer(c_int),    value :: Nr_elem
    integer(c_int),    value :: HDIM
    integer(c_int),    value :: Max_Nr_Neigh
    integer(c_int),    value :: I
    integer(c_int),    value :: J
    integer(c_int),    target :: Element_Pointer(Nr_atoms)
    real(c_double), value :: COULACC
    real(c_double), target :: DELTAQ(Nr_atoms)
    real(c_double), target :: RXYZ(3,Nr_atoms)
    real(c_double), target :: Box(3,3)
    real(c_double), target :: U(Nr_elem)
    real(c_double), intent(out) :: COULOMBV
    call Ewald_Real_Space_Single_latte(COULOMBV, I, RXYZ, Box, Nr_elem, DELTAQ, J, U, Element_Pointer, Nr_atoms, COULACC, HDIM, Max_Nr_Neigh)
  end subroutine Ewald_Real_Space_Single_latte_c

  subroutine Ewald_Real_Space_Single_c(COULOMBV, FCOUL, I, RX, RY, RZ, LBox, DELTAQ, J, U, Element_Type, Nr_atoms, COULACC, TIMERATIO, HDIM, Max_Nr_Neigh) bind(C, name="Ewald_Real_Space_Single")
    integer(c_int),    value :: Nr_atoms
    integer(c_int),    value :: HDIM
    integer(c_int),    value :: Max_Nr_Neigh
    integer(c_int),    value :: I
    integer(c_int),    value :: J
    real(c_double), value :: COULACC
    real(c_double), value :: TIMERATIO
    real(c_double), target :: DELTAQ(Nr_atoms)
    real(c_double), target :: RX(Nr_atoms)
    real(c_double), target :: RY(Nr_atoms)
    real(c_double), target :: RZ(Nr_atoms)
    real(c_double), target :: LBox(3)
    real(c_double), target :: U(Nr_atoms)
    character(c_char), target :: Element_Type(Nr_atoms)
    real(c_double), intent(out) :: COULOMBV
    real(c_double), intent(out) :: FCOUL(3)
    call Ewald_Real_Space_Single(COULOMBV, FCOUL, I, RX, RY, RZ, LBox, DELTAQ, J, U, Element_Type, Nr_atoms, COULACC, TIMERATIO, HDIM, Max_Nr_Neigh)
  end subroutine Ewald_Real_Space_Single_c

  subroutine Ewald_Real_Space_Matrix_latte_c(E, RXYZ, Box, U, Element_Pointer, Nr_atoms, COULACC, nebcoul, totnebcoul, HDIM, Max_Nr_Neigh, Nr_Elem) bind(C, name="Ewald_Real_Space_Matrix_latte")
    integer(c_int),    value :: Nr_atoms
    integer(c_int),    value :: HDIM
    integer(c_int),    value :: Max_Nr_Neigh
    integer(c_int),    value :: Nr_Elem
    real(c_double), value :: COULACC
    real(c_double), intent(out) :: E(Nr_atoms,Nr_atoms)
    real(c_double), target :: RXYZ(3,Nr_atoms)
    real(c_double), target :: Box(3,3)
    real(c_double), target :: U(Nr_elem)
    integer(c_int), target :: Element_Pointer(Nr_atoms)
    integer(c_int),    target :: totnebcoul(Nr_atoms)
    integer(c_int),    target :: nebcoul(4,Max_Nr_Neigh,Nr_atoms)
    call Ewald_Real_Space_Matrix_latte(E, RXYZ, Box, U, Element_Pointer, Nr_atoms, COULACC, nebcoul, totnebcoul, HDIM, Max_Nr_Neigh, Nr_Elem)
  end subroutine Ewald_Real_Space_Matrix_latte_c

  subroutine Ewald_Real_Space_latte_c(COULOMBV, I, RXYZ, Box, DELTAQ, U, Element_Pointer, Nr_atoms, COULACC, nebcoul, totnebcoul, HDIM, Max_Nr_Neigh, Nr_Elem) bind(C, name="Ewald_Real_Space_latte")
    integer(c_int),    value :: Nr_atoms
    integer(c_int),    value :: HDIM
    integer(c_int),    value :: Max_Nr_Neigh
    integer(c_int),    value :: I
    integer(c_int),    value :: Nr_Elem
    real(c_double), value :: COULACC
    real(c_double), target :: RXYZ(3,Nr_atoms)
    real(c_double), target :: Box(3,3)
    real(c_double), target :: DELTAQ(Nr_atoms)
    real(c_double), target :: U(Nr_elem)
    integer(c_int), target :: Element_Pointer(Nr_atoms)
    integer(c_int),    target :: totnebcoul(Nr_atoms)
    integer(c_int),    target :: nebcoul(4,Max_Nr_Neigh,Nr_atoms)
    real(c_double), intent(out) :: COULOMBV
    call Ewald_Real_Space_latte(COULOMBV, I, RXYZ, Box, DELTAQ, U, Element_Pointer, Nr_atoms, COULACC, nebcoul, totnebcoul, HDIM, Max_Nr_Neigh, Nr_Elem)
  end subroutine Ewald_Real_Space_latte_c

  subroutine Ewald_Real_Space_Test_c(COULOMBV, I, RX, RY, RZ, LBox, DELTAQ, U, Element_Type, Nr_atoms, COULACC, nnRx, nnRy, nnRz, nrnnlist, nnType, Max_Nr_Neigh) bind(C, name="Ewald_Real_Space_Test")
    integer(c_int),    value :: Nr_atoms
    integer(c_int),    value :: Max_Nr_Neigh
    integer(c_int),    value :: I
    real(c_double), value :: COULACC
    real(c_double), target :: RX(Nr_atoms)
    real(c_double), target :: RY(Nr_atoms)
    real(c_double), target :: RZ(Nr_atoms)
    real(c_double), target :: LBox(3)
    real(c_double), target :: DELTAQ(Nr_atoms)
    real(c_double), target :: U(Nr_atoms)
    character(c_char), target :: Element_Type(Nr_atoms)
    integer(c_int),    target :: nrnnlist(Nr_atoms)
    integer(c_int),    target :: nnType(Nr_atoms,Max_Nr_Neigh)
    real(c_double), target :: nnRx(Nr_atoms,Max_Nr_Neigh)
    real(c_double), target :: nnRy(Nr_atoms,Max_Nr_Neigh)
    real(c_double), target :: nnRz(Nr_atoms,Max_Nr_Neigh)
    real(c_double), intent(out) :: COULOMBV
    call Ewald_Real_Space_Test(COULOMBV, I, RX, RY, RZ, LBox, DELTAQ, U, Element_Type, Nr_atoms, COULACC, nnRx, nnRy, nnRz, nrnnlist, nnType, Max_Nr_Neigh)
  end subroutine Ewald_Real_Space_Test_c

  subroutine Ewald_Real_Space_c(COULOMBV, FCOUL, I, RX, RY, RZ, LBox, DELTAQ, U, Element_Type, Nr_atoms, COULACC, TIMERATIO, nnRx, nnRy, nnRz, nrnnlist, nnType, HDIM, Max_Nr_Neigh) bind(C, name="Ewald_Real_Space")
    integer(c_int),    value :: Nr_atoms
    integer(c_int),    value :: HDIM
    integer(c_int),    value :: Max_Nr_Neigh
    integer(c_int),    value :: I
    real(c_double), value :: COULACC
    real(c_double), value :: TIMERATIO
    real(c_double), target :: RX(Nr_atoms)
    real(c_double), target :: RY(Nr_atoms)
    real(c_double), target :: RZ(Nr_atoms)
    real(c_double), target :: LBox(3)
    real(c_double), target :: DELTAQ(Nr_atoms)
    real(c_double), target :: U(Nr_atoms)
    character(c_char), target :: Element_Type(Nr_atoms)
    integer(c_int),    target :: nrnnlist(Nr_atoms)
    integer(c_int),    target :: nnType(Nr_atoms,Max_Nr_Neigh)
    real(c_double), target :: nnRx(Nr_atoms,Max_Nr_Neigh)
    real(c_double), target :: nnRy(Nr_atoms,Max_Nr_Neigh)
    real(c_double), target :: nnRz(Nr_atoms,Max_Nr_Neigh)
    real(c_double), intent(out) :: COULOMBV
    real(c_double), intent(out) :: FCOUL(3)
    call Ewald_Real_Space(COULOMBV, FCOUL, I, RX, RY, RZ, LBox, DELTAQ, U, Element_Type, Nr_atoms, COULACC, TIMERATIO, nnRx, nnRy, nnRz, nrnnlist, nnType, HDIM, Max_Nr_Neigh)
  end subroutine Ewald_Real_Space_c

  !------------------------------------------------
  ! prg_genz_mod
  !------------------------------------------------
  !  subroutine prg_parse_ZSP_c(input, filename) bind(C, name="prg_parse_ZSP")
  !    character(c_char), value :: filename
  !    call prg_parse_ZSP(input, filename)
  !  end subroutine prg_parse_ZSP_c
  !
  !  subroutine prg_init_ZSPmat_c(igenz, zk1_bml_c, zk2_bml_c, zk3_bml_c, zk4_bml_c, zk5_bml_c, zk6_bml_c, norb, bml_type, bml_element_type) bind(C, name="prg_init_ZSPmat")
  !    integer(c_int), value :: norb
  !    integer(c_int), value :: igenz
  !    character(c_char), value :: bml_type
  !    character(c_char), optional :: bml_element_type
  !    type(c_ptr), value :: zk1_bml_c
  !    type(bml_matrix_t) :: zk1_bml
  !    type(c_ptr), value :: zk2_bml_c
  !    type(bml_matrix_t) :: zk2_bml
  !    type(c_ptr), value :: zk3_bml_c
  !    type(bml_matrix_t) :: zk3_bml
  !    type(c_ptr), value :: zk4_bml_c
  !    type(bml_matrix_t) :: zk4_bml
  !    type(c_ptr), value :: zk5_bml_c
  !    type(bml_matrix_t) :: zk5_bml
  !    type(c_ptr), value :: zk6_bml_c
  !    type(bml_matrix_t) :: zk6_bml
  !    zk1_bml%ptr = zk1_bml_c
  !    zk2_bml%ptr = zk2_bml_c
  !    zk3_bml%ptr = zk3_bml_c
  !    zk4_bml%ptr = zk4_bml_c
  !    zk5_bml%ptr = zk5_bml_c
  !    zk6_bml%ptr = zk6_bml_c
  !    call prg_init_ZSPmat(igenz, zk1_bml, zk2_bml, zk3_bml, zk4_bml, zk5_bml, zk6_bml, norb, bml_type, bml_element_type)
  !  end subroutine prg_init_ZSPmat_c
  !
  !  subroutine prg_buildZdiag_c(smat_bml_c, zmat_bml_c, threshold, mdimin, bml_type, verbose) bind(C, name="prg_buildZdiag")
  !    character(c_char), value :: bml_type
  !    integer(c_int), value :: mdimin
  !    integer(c_int)  :: verbose
  !    real(c_double), value :: threshold
  !    type(c_ptr), value :: zmat_bml_c
  !    type(bml_matrix_t) :: zmat_bml
  !    type(c_ptr), value :: smat_bml_c
  !    type(bml_matrix_t) :: smat_bml
  !    zmat_bml%ptr = zmat_bml_c
  !    smat_bml%ptr = smat_bml_c
  !    call prg_buildZdiag(smat_bml, zmat_bml, threshold, mdimin, bml_type, verbose)
  !  end subroutine prg_buildZdiag_c
  !
  !  subroutine prg_genz_sp_initialz0_c(smat_bml_c, zmat_bml_c, norb, mdim, bml_type_f, threshold) bind(C, name="prg_genz_sp_initialz0")
  !    character(c_char), value :: bml_type_f
  !    integer(c_int), value :: mdim
  !    integer(c_int), value :: norb
  !    real(c_double), value :: threshold
  !    type(c_ptr), value :: zmat_bml_c
  !    type(bml_matrix_t) :: zmat_bml
  !    type(c_ptr), value :: smat_bml_c
  !    type(bml_matrix_t) :: smat_bml
  !    zmat_bml%ptr = zmat_bml_c
  !    smat_bml%ptr = smat_bml_c
  !    call prg_genz_sp_initialz0(smat_bml, zmat_bml, norb, mdim, bml_type_f, threshold)
  !  end subroutine prg_genz_sp_initialz0_c
  !
  !  subroutine prg_genz_sp_initial_zmat_c(smat_bml_c, zmat_bml_c, norb, mdim, bml_type_f, threshold) bind(C, name="prg_genz_sp_initial_zmat")
  !    character(c_char), value :: bml_type_f
  !    integer(c_int), value :: mdim
  !    integer(c_int), value :: norb
  !    real(c_double), value :: threshold
  !    type(c_ptr), value :: zmat_bml_c
  !    type(bml_matrix_t) :: zmat_bml
  !    type(c_ptr), value :: smat_bml_c
  !    type(bml_matrix_t) :: smat_bml
  !    zmat_bml%ptr = zmat_bml_c
  !    smat_bml%ptr = smat_bml_c
  !    call prg_genz_sp_initial_zmat(smat_bml, zmat_bml, norb, mdim, bml_type_f, threshold)
  !  end subroutine prg_genz_sp_initial_zmat_c
  !
  !  subroutine prg_genz_sp_ref_c(smat_bml_c, zmat_bml, nref, norb, bml_type, threshold) bind(C, name="prg_genz_sp_ref")
  !    integer(c_int), value :: norb
  !    type(c_ptr), value :: smat_bml_c
  !    type(bml_matrix_t) :: smat_bml
  !    smat_bml%ptr = smat_bml_c
  !    call prg_genz_sp_ref(smat_bml, zmat_bml, nref, norb, bml_type, threshold)
  !  end subroutine prg_genz_sp_ref_c

  !------------------------------------------------
  ! prg_grah_mod
  !------------------------------------------------
  ! graph module has subgraph_t data type (TODO)

  !------------------------------------------------
  ! prg_normalize_mod
  !------------------------------------------------

  subroutine prg_normalize_c(h_bml_c) bind(C, name="prg_normalize")
    type(c_ptr),value :: h_bml_c
    type(bml_matrix_t) :: h_bml
    h_bml%ptr = h_bml_c
    call prg_normalize(h_bml)
  end subroutine prg_normalize_c

  subroutine prg_normalize_fermi_c(h_bml_c, h1, hN, mu) bind(C, name="prg_normalize_fermi")
    type(c_ptr),value :: h_bml_c
    type(bml_matrix_t) :: h_bml
    real(c_double), value :: h1
    real(c_double), value :: hN
    real(c_double), value :: mu
    h_bml%ptr = h_bml_c
    call prg_normalize_fermi(h_bml, h1, hN, mu)
  end subroutine prg_normalize_fermi_c

  subroutine prg_normalize_implicit_fermi_c(h_bml_c, cnst, mu) bind(C, name="prg_normalize_implicit_fermi")
    type(c_ptr),value :: h_bml_c
    type(bml_matrix_t) :: h_bml
    real(c_double), value :: cnst
    real(c_double), value :: mu
    h_bml%ptr = h_bml_c
    call prg_normalize_implicit_fermi(h_bml, cnst, mu)
  end subroutine prg_normalize_implicit_fermi_c

  ! gp is a graph_partitioning_t type (todo)
  !subroutine prg_gershgorinReduction_c(gp) bind(C, name="prg_gershgorinReduction")
  !  call prg_gershgorinReduction(gp)
  !end subroutine prg_gershgorinReduction_c

  subroutine prg_normalize_cheb_c(h_bml_c, mu, emin, emax, alpha, scaledmu) bind(C, name="prg_normalize_cheb")
    type(c_ptr),value :: h_bml_c
    type(bml_matrix_t) :: h_bml
    real(c_double), value :: mu
    real(c_double), value :: emin
    real(c_double), value :: emax
    real(c_double), value :: scaledmu
    real(c_double), value :: alpha
    h_bml%ptr = h_bml_c
    call prg_normalize_cheb(h_bml, mu, emin, emax, alpha, scaledmu)
  end subroutine prg_normalize_cheb_c

  !------------------------------------------------
  ! prg_sys_mod
  !------------------------------------------------

  subroutine prg_get_nameandext_c(fullfilename, filename, ext) bind(C, name="prg_get_nameandext")
    character(c_char), value :: fullfilename
    character(c_char), value :: filename
    character(c_char), value :: ext
    call prg_get_nameandext(fullfilename, filename, ext)
  end subroutine prg_get_nameandext_c

  !  subroutine prg_parse_system_c(system, filename, extin) bind(C, name="prg_parse_system")
  !
  !    type(system_type), intent(out)  ::  system
  !    character(c_char) :: extin
  !    character(c_char), value :: filename
  !    call prg_parse_system(system, filename, extin)
  !  end subroutine prg_parse_system_c
  !
  !  subroutine prg_destroy_system_c(sy) bind(C, name="prg_destroy_system")
  !
  !    type(system_type) ::  sy
  !    call prg_destroy_system(sy)
  !  end subroutine prg_destroy_system_c
  !
  !  subroutine prg_destroy_estr_c(estr) bind(C, name="prg_destroy_estr")
  !
  !    type(estruct_type) ::  estr
  !    call prg_destroy_estr(estr)
  !  end subroutine prg_destroy_estr_c
  !
  !  subroutine prg_write_system_c(system, filename, extin) bind(C, name="prg_write_system")
  !    type(system_type) ::  system
  !    character(c_char) :: extin
  !    character(c_char), value :: filename
  !    call prg_write_system(system, filename, extin)
  !  end subroutine prg_write_system_c

  !  subroutine prg_write_trajectory_c(system, iter, each, prg_deltat, filename, extension) bind(C, name="prg_write_trajectory")
  !    type(system_type), intent(out)  ::  system
  !    character(c_char), value :: filename
  !    character(c_char), value :: extension
  !    integer(c_int), value :: iter
  !    integer(c_int), value :: each
  !    real(c_double), value :: prg_deltat
  !    call prg_write_trajectory(system, iter, each, prg_deltat, filename, extension)
  !  end subroutine prg_write_trajectory_c
  !
  !  subroutine prg_write_trajectoryandproperty_c(system, iter, each, prg_deltat, scalarprop, filename, extension) bind(C, name="prg_write_trajectoryandproperty")
  !    character(c_char), value :: filename
  !    character(c_char), value :: extension
  !    integer(c_int), value :: iter
  !    integer(c_int), value :: each
  !    real(c_double), value :: prg_deltat
  !    real(c_double), target :: scalarprop(:)
  !    call prg_write_trajectoryandproperty(system, iter, each, prg_deltat, scalarprop, filename, extension)
  !  end subroutine prg_write_trajectoryandproperty_c
  !
  !  subroutine prg_make_random_system_c(system, nats, seed, lx, ly, lz) bind(C, name="prg_make_random_system")
  !    integer(c_int), value :: nats
  !    integer(c_int), value :: seed
  !    real(c_double), value :: lx
  !    real(c_double), value :: ly
  !    real(c_double), value :: lz
  !    call prg_make_random_system(system, nats, seed, lx, ly, lz)
  !  end subroutine prg_make_random_system_c
  !
  !  subroutine prg_parameters_to_vectors_c(abc_angles, lattice_vector) bind(C, name="prg_parameters_to_vectors")
  !    real(c_double), target :: abc_angles(2,3)
  !    real(c_double), intent(out) :: lattice_vector(3,3)
  !    call prg_parameters_to_vectors(abc_angles, lattice_vector)
  !  end subroutine prg_parameters_to_vectors_c
  !
  !  subroutine prg_vectors_to_parameters_c(lattice_vector, abc_angles) bind(C, name="prg_vectors_to_parameters")
  !    real(c_double), target :: lattice_vector(3,3)
  !    real(c_double), intent(out) :: abc_angles(2,3)
  !    call prg_vectors_to_parameters(lattice_vector, abc_angles)
  !  end subroutine prg_vectors_to_parameters_c
  !
  !  subroutine prg_get_origin_c(coords, origin) bind(C, name="prg_get_origin")
  !    real(c_double), allocatable,  :: origin(:)
  !    real(c_double), target :: coords(:,:)
  !    call prg_get_origin(coords, origin)
  !  end subroutine prg_get_origin_c
  !
  !  subroutine prg_get_distancematrix_c(coords, dmat) bind(C, name="prg_get_distancematrix")
  !    real(c_double), target :: coords(:,:)
  !    real(c_double), intent(out), allocatable :: dmat(:,:)
  !    call prg_get_distancematrix(coords, dmat)
  !  end subroutine prg_get_distancematrix_c
  !
  !  subroutine prg_translateandfoldtobox_c(coords, lattice_vectors, origin, verbose) bind(C, name="prg_translateandfoldtobox")
  !    integer(c_int) :: verbose
  !    real(c_double), allocatable,  :: origin(:)
  !    real(c_double), allocatable,  :: coords(:,:)
  !    real(c_double), target :: lattice_vectors(:,:)
  !    call prg_translateandfoldtobox(coords, lattice_vectors, origin, verbose)
  !  end subroutine prg_translateandfoldtobox_c
  !
  !  subroutine prg_centeratbox_c(coords, lattice_vectors, verbose) bind(C, name="prg_centeratbox")
  !    integer(c_int) :: verbose
  !    real(c_double), allocatable,  :: coords(:,:)
  !    real(c_double), target :: lattice_vectors(:,:)
  !    call prg_centeratbox(coords, lattice_vectors, verbose)
  !  end subroutine prg_centeratbox_c
  !
  !  subroutine prg_wraparound_c(coords, lattice_vectors, index, verbose) bind(C, name="prg_wraparound")
  !    integer(c_int), value :: index
  !    integer(c_int) :: verbose
  !    real(c_double), allocatable,  :: coords(:,:)
  !    real(c_double), target :: lattice_vectors(:,:)
  !    call prg_wraparound(coords, lattice_vectors, index, verbose)
  !  end subroutine prg_wraparound_c
  !
  !  subroutine prg_translatetogeomcandfoldtobox_c(coords, lattice_vectors, origin) bind(C, name="prg_translatetogeomcandfoldtobox")
  !    real(c_double), allocatable,  :: origin(:)
  !    real(c_double), allocatable,  :: coords(:,:)
  !    real(c_double), target :: lattice_vectors(:,:)
  !    call prg_translatetogeomcandfoldtobox(coords, lattice_vectors, origin)
  !  end subroutine prg_translatetogeomcandfoldtobox_c
  !
  !  subroutine prg_replicate_c(coords, symbols, lattice_vectors, nx, ny, nz) bind(C, name="prg_replicate")
  !    integer(c_int), value :: nx
  !    integer(c_int), value :: ny
  !    integer(c_int), value :: nz
  !    character(c_char), allocatable,  :: symbols(:)
  !    real(c_double), allocatable,  :: coords(:,:)
  !    real(c_double), target :: lattice_vectors(:,:)
  !    call prg_replicate(coords, symbols, lattice_vectors, nx, ny, nz)
  !  end subroutine prg_replicate_c
  !
  !  subroutine prg_replicate_system_c(sy, syf, nx, ny, nz) bind(C, name="prg_replicate_system")
  !    integer(c_int), value :: nx
  !    integer(c_int), value :: ny
  !    integer(c_int), value :: nz
  !    call prg_replicate_system(sy, syf, nx, ny, nz)
  !  end subroutine prg_replicate_system_c
  !
  !  subroutine prg_cleanuprepeatedatoms_c(nats, coords, symbols, verbose) bind(C, name="prg_cleanuprepeatedatoms")
  !    integer(c_int), value :: nats
  !    integer(c_int) :: verbose
  !    real(c_double), allocatable,  :: coords(:,:)
  !    character(c_char), allocatable,  :: symbols(:)
  !    call prg_cleanuprepeatedatoms(nats, coords, symbols, verbose)
  !  end subroutine prg_cleanuprepeatedatoms_c
  !
  !  subroutine prg_get_recip_vects_c(lattice_vectors, recip_vectors, volr, volk) bind(C, name="prg_get_recip_vects")
  !    real(c_double), allocatable,  :: recip_vectors(:,:)
  !    real(c_double), target :: lattice_vectors(:,:)
  !    real(c_double), value :: volk
  !    real(c_double), value :: volr
  !    call prg_get_recip_vects(lattice_vectors, recip_vectors, volr, volk)
  !  end subroutine prg_get_recip_vects_c
  !
  !  subroutine prg_get_dihedral_c(coords, id1, id2, id3, id4, dihedral) bind(C, name="prg_get_dihedral")
  !    real(c_double), target :: coords(:,:)
  !    real(c_double), intent(out) :: dihedral
  !    integer(c_int), value :: id1
  !    integer(c_int), value :: id2
  !    integer(c_int), value :: id3
  !    integer(c_int), value :: id4
  !    call prg_get_dihedral(coords, id1, id2, id3, id4, dihedral)
  !  end subroutine prg_get_dihedral_c
  !
  !  subroutine prg_get_covgraph_c(sy, nnStruct, nrnnstruct, bml_type, factor, gcov_bml_c, mdimin, verbose) bind(C, name="prg_get_covgraph")
  !    character(c_char), value :: bml_type
  !    integer(c_int), value :: mdimin
  !    integer(c_int), target :: nnStruct(:,:)
  !    integer(c_int), target :: nrnnstruct(:)
  !    integer(c_int)  :: verbose
  !    real(c_double), value :: factor
  !    type(c_ptr), value :: gcov_bml_c
  !    type(bml_matrix_t) :: gcov_bml
  !    gcov_bml%ptr = gcov_bml_c
  !    call prg_get_covgraph(sy, nnStruct, nrnnstruct, bml_type, factor, gcov_bml, mdimin, verbose)
  !  end subroutine prg_get_covgraph_c
  !
  !  subroutine prg_get_covgraph_h_c(sy, nnStruct, nrnnstruct, rcut, graph_h, mdimin, verbose) bind(C, name="prg_get_covgraph_h")
  !    integer(c_int), target :: nnStruct(:,:)
  !    integer(c_int), target :: nrnnstruct(:)
  !    integer(c_int), value :: mdimin
  !    integer(c_int)  :: verbose
  !    real(c_double), value :: rcut
  !    integer(c_int), allocatable,  :: graph_h(:,:)
  !    call prg_get_covgraph_h(sy, nnStruct, nrnnstruct, rcut, graph_h, mdimin, verbose)
  !  end subroutine prg_get_covgraph_h_c
  !
  !  subroutine prg_get_subsystem_c(sy, lsize, indices, sbsy, verbose) bind(C, name="prg_get_subsystem")
  !    integer(c_int), target :: indices(:)
  !    integer(c_int), value :: lsize
  !    integer(c_int)  :: verbose
  !    call prg_get_subsystem(sy, lsize, indices, sbsy, verbose)
  !  end subroutine prg_get_subsystem_c
  !
  !  subroutine prg_destroy_subsystems_c(sbsy, verbose) bind(C, name="prg_destroy_subsystems")
  !    integer(c_int)  :: verbose
  !    call prg_destroy_subsystems(sbsy, verbose)
  !  end subroutine prg_destroy_subsystems_c
  !
  !  subroutine prg_molpartition_c(sy, npart, nnStructMindist, nnStruct, nrnnstruct, hetatm, gp, verbose) bind(C, name="prg_molpartition")
  !    character(c_char), value :: hetatm
  !    integer(c_int), target :: nnStruct(:,:)
  !    integer(c_int), target :: nrnnstruct(:)
  !    integer(c_int), value :: npart
  !    integer(c_int)  :: verbose
  !    real(c_double), target :: nnStructMindist(:,:)
  !    call prg_molpartition(sy, npart, nnStructMindist, nnStruct, nrnnstruct, hetatm, gp, verbose)
  !  end subroutine prg_molpartition_c
  !
  !  subroutine prg_get_partial_atomgraph_c(rho_bml_c, hindex, gch_bml_c, threshold, verbose) bind(C, name="prg_get_partial_atomgraph")
  !    integer(c_int), target :: hindex(:,:)
  !    integer(c_int)  :: verbose
  !    real(c_double), value :: threshold
  !    type(c_ptr), value :: rho_bml_c
  !    type(bml_matrix_t) :: rho_bml
  !    type(c_ptr), value :: gch_bml_c
  !    type(bml_matrix_t) :: gch_bml
  !    rho_bml%ptr = rho_bml_c
  !    gch_bml%ptr = gch_bml_c
  !    call prg_get_partial_atomgraph(rho_bml, hindex, gch_bml, threshold, verbose)
  !  end subroutine prg_get_partial_atomgraph_c
  !
  !  subroutine prg_collect_graph_p_c(rho_bml_c, nc, nats, hindex, chindex, graph_p, threshold, mdimin, verbose) bind(C, name="prg_collect_graph_p")
  !    integer(c_int), allocatable,  :: graph_p(:,:)
  !    integer(c_int), target :: chindex(:)
  !    integer(c_int), target :: hindex(:,:)
  !    integer(c_int), value :: nats
  !    integer(c_int), value :: nc
  !    integer(c_int), value :: mdimin
  !    integer(c_int)  :: verbose
  !    real(c_double), value :: threshold
  !    type(c_ptr), value :: rho_bml_c
  !    type(bml_matrix_t) :: rho_bml
  !    rho_bml%ptr = rho_bml_c
  !    call prg_collect_graph_p(rho_bml, nc, nats, hindex, chindex, graph_p, threshold, mdimin, verbose)
  !  end subroutine prg_collect_graph_p_c
  !
  !  subroutine prg_merge_graph_c(graph_p, graph_h) bind(C, name="prg_merge_graph")
  !    integer(c_int),              target :: graph_h(:,:)
  !    integer(c_int), target :: graph_p(:,:)
  !    call prg_merge_graph(graph_p, graph_h)
  !  end subroutine prg_merge_graph_c
  !
  !  subroutine prg_merge_graph_adj_c(graph_p, graph_h, xadj, adjncy) bind(C, name="prg_merge_graph_adj")
  !    integer(c_int), allocatable,  :: adjncy(:)
  !    integer(c_int), allocatable,  :: graph_h(:,:)
  !    integer(c_int), allocatable,  :: graph_p(:,:)
  !    integer(c_int), allocatable,  :: xadj(:)
  !    call prg_merge_graph_adj(graph_p, graph_h, xadj, adjncy)
  !  end subroutine prg_merge_graph_adj_c
  !
  !  subroutine prg_adj2bml_c(xadj, adjncy, bml_type, g_bml_c) bind(C, name="prg_adj2bml")
  !    character(c_char), value :: bml_type
  !    integer(c_int), target :: adjncy(:)
  !    integer(c_int), target :: xadj(:)
  !    type(c_ptr), value :: g_bml_c
  !    type(bml_matrix_t) :: g_bml
  !    g_bml%ptr = g_bml_c
  !    call prg_adj2bml(xadj, adjncy, bml_type, g_bml)
  !  end subroutine prg_adj2bml_c
  !
  !  subroutine prg_graph2bml_c(graph, bml_type, g_bml_c) bind(C, name="prg_graph2bml")
  !    character(c_char), value :: bml_type
  !    integer(c_int), allocatable,  :: graph(:,:)
  !    type(c_ptr), value :: g_bml_c
  !    type(bml_matrix_t) :: g_bml
  !    g_bml%ptr = g_bml_c
  !    call prg_graph2bml(graph, bml_type, g_bml)
  !  end subroutine prg_graph2bml_c
  !
  !  subroutine prg_graph2vector_c(graph, vector, maxnz) bind(C, name="prg_graph2vector")
  !    integer(c_int), target :: graph(:,:)
  !    integer(c_int), allocatable :: vector(:)
  !    integer(c_int), value :: maxnz
  !    call prg_graph2vector(graph, vector, maxnz)
  !  end subroutine prg_graph2vector_c
  !
  !  subroutine prg_vector2graph_c(vector, graph, maxnz) bind(C, name="prg_vector2graph")
  !    integer(c_int), target :: graph(:,:)
  !    integer(c_int), allocatable,  :: vector(:)
  !    integer(c_int), value :: maxnz
  !    call prg_vector2graph(vector, graph, maxnz)
  !  end subroutine prg_vector2graph_c
  !
  !  subroutine prg_sortadj_c(xadj, adjncy) bind(C, name="prg_sortadj")
  !    integer(c_int), target :: xadj(:)
  !    integer(c_int), allocatable,  :: adjncy(:)
  !    call prg_sortadj(xadj, adjncy)
  !  end subroutine prg_sortadj_c


  !------------------------------------------------
  ! prg_sp2_mod
  !------------------------------------------------

  subroutine prg_sp2_fermi_init_c(h_bml_c, nsteps, nocc, tscale, threshold, occErrLimit, traceLimit, x_bml_c, mu, beta, h1, hN, sgnlist) bind(C, name="prg_sp2_fermi_init")
    type(c_ptr), value :: h_bml_c
    type(bml_matrix_t) :: h_bml
    type(c_ptr), value :: x_bml_c
    type(bml_matrix_t) :: x_bml
    integer(c_int), value :: nsteps
    integer(c_int), target :: sgnlist(nsteps)
    real(c_double), value :: nocc
    real(c_double), value :: tscale
    real(c_double), value :: threshold
    real(c_double), value :: occErrLimit
    real(c_double), value :: traceLimit
    real(c_double), value :: mu
    real(c_double), value :: beta
    real(c_double), value :: h1
    real(c_double), value :: hN
    h_bml%ptr = h_bml_c
    x_bml%ptr = x_bml_c
    call prg_sp2_fermi_init(h_bml, nsteps, nocc, tscale, threshold, occErrLimit, traceLimit, x_bml, mu, beta, h1, hN, sgnlist)
  end subroutine prg_sp2_fermi_init_c

  subroutine prg_sp2_fermi_init_norecs_c(h_bml_c, nsteps, nocc, tscale, threshold, occErrLimit, traceLimit, x_bml_c, mu, beta, h1, hN, sgnlist, verbose) bind(C, name="prg_sp2_fermi_init_norecs")
    type(c_ptr), value :: h_bml_c
    type(bml_matrix_t) :: h_bml
    type(c_ptr), value :: x_bml_c
    type(bml_matrix_t) :: x_bml
    integer(c_int), value :: nsteps
    integer(c_int), target :: sgnlist(nsteps)
    real(c_double), value :: nocc
    real(c_double), value :: tscale
    real(c_double), value :: threshold
    real(c_double), value :: occErrLimit
    real(c_double), value :: traceLimit
    real(c_double), value :: mu
    real(c_double), value :: beta
    real(c_double), value :: h1
    real(c_double), value :: hN
    integer(c_int), optional :: verbose
    h_bml%ptr = h_bml_c
    x_bml%ptr = x_bml_c
    call prg_sp2_fermi_init_norecs(h_bml, nsteps, nocc, tscale, threshold, occErrLimit, traceLimit, x_bml, mu, beta, h1, hN, sgnlist, verbose)
  end subroutine prg_sp2_fermi_init_norecs_c

  subroutine prg_sp2_fermi_c(h_bml_c, osteps, nsteps, nocc, mu, beta, h1, hN, sgnlist, threshold, eps, traceLimit, x_bml_c) bind(C, name="prg_sp2_fermi")
    type(c_ptr), value :: h_bml_c
    type(bml_matrix_t) :: h_bml
    type(c_ptr), value :: x_bml_c
    type(bml_matrix_t) :: x_bml
    integer(c_int), value :: osteps
    integer(c_int), value :: nsteps
    integer(c_int), target :: sgnlist(nsteps)
    real(c_double), value :: nocc
    real(c_double), value :: threshold
    real(c_double), value :: eps
    real(c_double), value :: traceLimit
    real(c_double), value :: beta
    real(c_double), value :: h1
    real(c_double), value :: hN
    real(c_double), value :: mu
    h_bml%ptr = h_bml_c
    x_bml%ptr = x_bml_c
    call prg_sp2_fermi(h_bml, osteps, nsteps, nocc, mu, beta, h1, hN, sgnlist, threshold, eps, traceLimit, x_bml)
  end subroutine prg_sp2_fermi_c

  subroutine prg_sp2_entropy_function_c(mu, h1, hN, nsteps, sgnlist, GG, ee) bind(C, name="prg_sp2_entropy_function")
    real(c_double), value :: mu
    real(c_double), value :: h1
    real(c_double), value :: hN
    integer(c_int), value :: nsteps
    integer(c_int), target :: sgnlist(nsteps)
    real(c_double) :: GG(1001)
    real(c_double) :: ee(1001)
    real(c_double), allocatable :: GG_tmp(:), ee_tmp(:)
    call prg_sp2_entropy_function(mu, h1, hN, nsteps, sgnlist, GG_tmp, ee_tmp)
    GG = GG_tmp
    ee = ee_tmp
    deallocate(GG_tmp, ee_tmp)
  end subroutine prg_sp2_entropy_function_c

  !------------------------------------------------
  ! prg_pulaycomponent_mod
  !------------------------------------------------
  subroutine prg_PulayComponent0_c(rho_bml_c, ham_bml_c, pcm_bml_c, threshold, M, bml_type, verbose) bind(C, name="prg_PulayComponent0")
    type(c_ptr), value :: rho_bml_c
    type(bml_matrix_t) :: rho_bml
    type(c_ptr), value :: ham_bml_c
    type(bml_matrix_t) :: ham_bml
    type(c_ptr), value :: pcm_bml_c
    type(bml_matrix_t) :: pcm_bml
    integer(c_int), value :: M
    integer(c_int), value :: verbose
    real(c_double), value :: threshold
    character(c_char), value :: bml_type
    rho_bml%ptr = rho_bml_c
    ham_bml%ptr = ham_bml_c
    pcm_bml%ptr = pcm_bml_c
    call prg_PulayComponent0(rho_bml, ham_bml, pcm_bml, threshold, M, bml_type, verbose)
  end subroutine prg_PulayComponent0_c

  subroutine prg_PulayComponentT_c(rho_bml_c, ham_bml_c, zmat_bml_c, pcm_bml_c, threshold, M, bml_type, verbose) bind(C, name="prg_PulayComponentT")
    type(c_ptr), value :: rho_bml_c
    type(bml_matrix_t) :: rho_bml
    type(c_ptr), value :: ham_bml_c
    type(bml_matrix_t) :: ham_bml
    type(c_ptr), value :: zmat_bml_c
    type(bml_matrix_t) :: zmat_bml
    type(c_ptr), value :: pcm_bml_c
    type(bml_matrix_t) :: pcm_bml
    integer(c_int), value :: M
    integer(c_int), value :: verbose
    real(c_double), value :: threshold
    character(c_char), value :: bml_type
    rho_bml%ptr = rho_bml_c
    ham_bml%ptr = ham_bml_c
    zmat_bml%ptr = zmat_bml_c
    pcm_bml%ptr = pcm_bml_c
    call prg_PulayComponentT(rho_bml, ham_bml, zmat_bml, pcm_bml, threshold, M, bml_type, verbose)
  end subroutine prg_PulayComponentT_c

  !  subroutine prg_get_pulayforce_c(nats, zmat_bml_c, ham_bml_c, rho_bml_c, dSx_bml_c, dSy_bml_c, dSz_bml_c, hindex, FPUL, threshold) bind(C, name="prg_get_pulayforce")
  !    real(c_double), allocatable,  :: FPUL(:,:)
  !    integer(c_int), value :: nats
  !    type(c_ptr), value :: dSx_bml_c
  !    type(bml_matrix_t) :: dSx_bml
  !    type(c_ptr), value :: dSy_bml_c
  !    type(bml_matrix_t) :: dSy_bml
  !    type(c_ptr), value :: dSz_bml_c
  !    type(bml_matrix_t) :: dSz_bml
  !    type(c_ptr), value :: rho_bml_c
  !    type(bml_matrix_t) :: rho_bml
  !    type(c_ptr), value :: ham_bml_c
  !    type(bml_matrix_t) :: ham_bml
  !    type(c_ptr), value :: zmat_bml_c
  !    type(bml_matrix_t) :: zmat_bml
  !    integer(c_int), target :: hindex(:,:)
  !    real(c_double), value :: threshold
  !    dSx_bml%ptr = dSx_bml_c
  !    dSy_bml%ptr = dSy_bml_c
  !    dSz_bml%ptr = dSz_bml_c
  !    rho_bml%ptr = rho_bml_c
  !    ham_bml%ptr = ham_bml_c
  !    zmat_bml%ptr = zmat_bml_c
  !    call prg_get_pulayforce(nats, zmat_bml, ham_bml, rho_bml, dSx_bml, dSy_bml, dSz_bml, hindex, FPUL, threshold)
  !  end subroutine prg_get_pulayforce_c

  !------------------------------------------------
  ! prg_sp2_mod
  !------------------------------------------------

  subroutine prg_sp2_basic_c(h_bml_c, rho_bml_c, threshold, bndfil, minsp2iter, maxsp2iter, sp2conv, idemtol, verbose) bind(C, name="prg_sp2_basic")
    integer(c_int), value :: minsp2iter
    integer(c_int), value :: maxsp2iter
    integer(c_int), value :: verbose
    real(c_double), value :: bndfil
    real(c_double), value :: threshold
    real(c_double), value :: idemtol
    character(c_char), value :: sp2conv
    type(c_ptr), value :: rho_bml_c
    type(bml_matrix_t) :: rho_bml
    type(c_ptr), value :: h_bml_c
    type(bml_matrix_t) :: h_bml
    rho_bml%ptr = rho_bml_c
    h_bml%ptr = h_bml_c
    call prg_sp2_basic(h_bml, rho_bml, threshold, bndfil, minsp2iter, maxsp2iter, sp2conv, idemtol, verbose)
  end subroutine prg_sp2_basic_c

  subroutine prg_sp2_basic_tcore_c(h_bml_c, rho_bml_c, rhofull_bml_c, threshold, bndfil, minsp2iter, maxsp2iter, sp2conv, idemtol, verbose) bind(C, name="prg_sp2_basic_tcore")
    integer(c_int), value :: minsp2iter
    integer(c_int), value :: maxsp2iter
    integer(c_int), value :: verbose
    real(c_double), value :: bndfil
    real(c_double), value :: threshold
    real(c_double), value :: idemtol
    character(c_char), value :: sp2conv
    type(c_ptr), value :: rho_bml_c
    type(bml_matrix_t) :: rho_bml
    type(c_ptr), value :: rhofull_bml_c
    type(bml_matrix_t) :: rhofull_bml
    type(c_ptr), value :: h_bml_c
    type(bml_matrix_t) :: h_bml
    rho_bml%ptr = rho_bml_c
    rhofull_bml%ptr = rhofull_bml_c
    h_bml%ptr = h_bml_c
    call prg_sp2_basic_tcore(h_bml, rho_bml, rhofull_bml, threshold, bndfil, minsp2iter, maxsp2iter, sp2conv, idemtol, verbose)
  end subroutine prg_sp2_basic_tcore_c

  subroutine prg_sp2_alg2_c(h_bml_c, rho_bml_c, threshold, bndfil, minsp2iter, maxsp2iter, sp2conv, idemtol, verbose) bind(C, name="prg_sp2_alg2")
    integer(c_int), value :: minsp2iter
    integer(c_int), value :: maxsp2iter
    real(c_double), value :: threshold
    real(c_double), value :: bndfil
    real(c_double), value :: idemtol
    character(c_char), value :: sp2conv
    type(c_ptr),value :: h_bml_c
    type(bml_matrix_t) :: h_bml
    type(c_ptr),value :: rho_bml_c
    type(bml_matrix_t) :: rho_bml
    integer(c_int)  :: verbose
    h_bml%ptr = h_bml_c
    rho_bml%ptr = rho_bml_c
    call prg_sp2_alg2(h_bml, rho_bml, threshold, bndfil, minsp2iter, maxsp2iter, sp2conv, idemtol, verbose)
  end subroutine prg_sp2_alg2_c

  subroutine prg_sp2_alg2_genseq_c(h_bml_c, rho_bml_c, threshold, bndfil, minsp2iter, maxsp2iter, sp2conv, idemtol, pp, icount, vv, verbose) bind(C, name="prg_sp2_alg2_genseq")
    integer(c_int), value :: minsp2iter
    integer(c_int), value :: maxsp2iter
    integer(c_int), value :: icount
    integer(c_int), target :: pp(:)
    real(c_double), value :: threshold
    real(c_double), value :: bndfil
    real(c_double), value :: idemtol
    real(c_double), target :: vv(:)
    character(c_char), value :: sp2conv
    type(c_ptr),value :: h_bml_c
    type(bml_matrix_t) :: h_bml
    type(c_ptr),value :: rho_bml_c
    type(bml_matrix_t) :: rho_bml
    integer(c_int)  :: verbose
    h_bml%ptr = h_bml_c
    rho_bml%ptr = rho_bml_c
    call prg_sp2_alg2_genseq(h_bml, rho_bml, threshold, bndfil, minsp2iter, maxsp2iter, sp2conv, idemtol, pp, icount, vv, verbose)
  end subroutine prg_sp2_alg2_genseq_c

  subroutine prg_sp2_alg2_seq_c(h_bml_c, rho_bml_c, threshold, pp, icount, vv, verbose) bind(C, name="prg_sp2_alg2_seq")
    integer(c_int), value :: icount
    integer(c_int), target :: pp(:)
    real(c_double), value :: threshold
    real(c_double), target :: vv(:)
    type(c_ptr),value :: h_bml_c
    type(bml_matrix_t) :: h_bml
    type(c_ptr),value :: rho_bml_c
    type(bml_matrix_t) :: rho_bml
    integer(c_int)  :: verbose
    h_bml%ptr = h_bml_c
    rho_bml%ptr = rho_bml_c
    call prg_sp2_alg2_seq(h_bml, rho_bml, threshold, pp, icount, vv, verbose)
  end subroutine prg_sp2_alg2_seq_c

  subroutine prg_prg_sp2_alg2_seq_inplace_c(rho_bml_c, threshold, pp, icount, vv, mineval, maxeval, verbose) bind(C, name="prg_prg_sp2_alg2_seq_inplace")
    integer(c_int), value :: icount
    integer(c_int), target :: pp(:)
    real(c_double), value :: threshold
    real(c_double), target :: vv(:)
    real(c_double)  :: mineval
    real(c_double)  :: maxeval
    type(c_ptr),value :: rho_bml_c
    type(bml_matrix_t) :: rho_bml
    integer(c_int)  :: verbose
    rho_bml%ptr = rho_bml_c
    call prg_prg_sp2_alg2_seq_inplace(rho_bml, threshold, pp, icount, vv, mineval, maxeval, verbose)
  end subroutine prg_prg_sp2_alg2_seq_inplace_c

  subroutine prg_sp2_alg1_c(h_bml_c, rho_bml_c, threshold, bndfil, minsp2iter, maxsp2iter, sp2conv, idemtol, verbose) bind(C, name="prg_sp2_alg1")
    integer(c_int), value :: minsp2iter
    integer(c_int), value :: maxsp2iter
    integer(c_int)  :: verbose
    real(c_double), value :: threshold
    real(c_double), value :: bndfil
    real(c_double), value :: idemtol
    character(c_char), value :: sp2conv
    type(c_ptr),value :: h_bml_c
    type(bml_matrix_t) :: h_bml
    type(c_ptr),value :: rho_bml_c
    type(bml_matrix_t) :: rho_bml
    h_bml%ptr = h_bml_c
    rho_bml%ptr = rho_bml_c
    call prg_sp2_alg1(h_bml, rho_bml, threshold, bndfil, minsp2iter, maxsp2iter, sp2conv, idemtol, verbose)
  end subroutine prg_sp2_alg1_c

  subroutine prg_sp2_alg1_seq_c(h_bml_c, rho_bml_c, threshold, pp, icount, vv) bind(C, name="prg_sp2_alg1_seq")
    integer(c_int), value :: icount
    integer(c_int), target :: pp(:)
    real(c_double), value :: threshold
    real(c_double), target :: vv(:)
    type(c_ptr),value :: h_bml_c
    type(bml_matrix_t) :: h_bml
    type(c_ptr),value :: rho_bml_c
    type(bml_matrix_t) :: rho_bml
    h_bml%ptr = h_bml_c
    rho_bml%ptr = rho_bml_c
    call prg_sp2_alg1_seq(h_bml, rho_bml, threshold, pp, icount, vv)
  end subroutine prg_sp2_alg1_seq_c

  subroutine prg_prg_sp2_alg1_seq_inplace_c(rho_bml_c, threshold, pp, icount, vv, mineval, maxeval) bind(C, name="prg_prg_sp2_alg1_seq_inplace")
    integer(c_int), value :: icount
    integer(c_int), target :: pp(:)
    real(c_double), value :: threshold
    real(c_double), target :: vv(:)
    real(c_double), value :: mineval
    real(c_double), value :: maxeval
    type(c_ptr),value :: rho_bml_c
    type(bml_matrix_t) :: rho_bml
    rho_bml%ptr = rho_bml_c
    call prg_prg_sp2_alg1_seq_inplace(rho_bml, threshold, pp, icount, vv, mineval, maxeval)
  end subroutine prg_prg_sp2_alg1_seq_inplace_c

  subroutine prg_sp2_submatrix_c(ham_bml_c, rho_bml_c, threshold, pp, icount, vv, mineval, maxeval, core_size) bind(C, name="prg_sp2_submatrix")
    integer(c_int), value :: icount
    integer(c_int), target :: pp(:)
    integer(c_int), value :: core_size
    real(c_double), value :: threshold
    real(c_double), target :: vv(:)
    real(c_double), value :: mineval
    real(c_double), value :: maxeval
    type(c_ptr),value :: ham_bml_c
    type(bml_matrix_t) :: ham_bml
    type(c_ptr),value :: rho_bml_c
    type(bml_matrix_t) :: rho_bml
    ham_bml%ptr = ham_bml_c
    rho_bml%ptr = rho_bml_c
    call prg_sp2_submatrix(ham_bml, rho_bml, threshold, pp, icount, vv, mineval, maxeval, core_size)
  end subroutine prg_sp2_submatrix_c

  subroutine prg_sp2_submatrix_inplace_c(rho_bml_c, threshold, pp, icount, vv, mineval, maxeval, core_size) bind(C, name="prg_sp2_submatrix_inplace")
    integer(c_int), value :: icount
    integer(c_int), target :: pp(:)
    integer(c_int), value :: core_size
    real(c_double), value :: threshold
    real(c_double), target :: vv(:)
    real(c_double), value :: mineval
    real(c_double), value :: maxeval
    type(c_ptr),value :: rho_bml_c
    type(bml_matrix_t) :: rho_bml
    rho_bml%ptr = rho_bml_c
    call prg_sp2_submatrix_inplace(rho_bml, threshold, pp, icount, vv, mineval, maxeval, core_size)
  end subroutine prg_sp2_submatrix_inplace_c

  !------------------------------------------------
  ! prg_timer_mod
  !------------------------------------------------

  subroutine timer_prg_init_c() bind(C, name="timer_prg_init")
    call timer_prg_init()
  end subroutine timer_prg_init_c

  subroutine prg_timer_shutdown_c() bind(C, name="prg_timer_shutdown")
    call prg_timer_shutdown()
  end subroutine prg_timer_shutdown_c

  subroutine prg_timer_start_c(itimer, tag) bind(C, name="prg_timer_start")
    integer(c_int) :: itimer
    character(c_char), value :: tag
    call prg_timer_start(itimer, tag)
  end subroutine prg_timer_start_c

  subroutine prg_timer_stop_c(itimer, verbose) bind(C, name="prg_timer_stop")
    integer(c_int) :: itimer, verbose
    call prg_timer_stop(itimer, verbose)
  end subroutine prg_timer_stop_c

  subroutine prg_timer_collect_c() bind(C, name="prg_timer_collect")
    call prg_timer_collect()
  end subroutine prg_timer_collect_c

  subroutine prg_timer_results_c() bind(C, name="prg_timer_results")
    call prg_timer_results()
  end subroutine prg_timer_results_c

  subroutine prg_print_date_and_time_c(tag) bind(C, name="prg_print_date_and_time")
    character(c_char), value :: tag
    call prg_print_date_and_time(tag)
  end subroutine prg_print_date_and_time_c



  !------------------------------------------------
  ! prg_openfiles_mod
  !------------------------------------------------

  subroutine prg_open_file_c(io, name) bind(C, name="prg_open_file")
    character(c_char), value :: name
    integer(c_int), value :: io
    call prg_open_file(io, name)
  end subroutine prg_open_file_c

  subroutine prg_open_file_to_read_c(io, name) bind(C, name="prg_open_file_to_read")
    character(c_char), value :: name
    integer(c_int), value :: io
    call prg_open_file_to_read(io, name)
  end subroutine prg_open_file_to_read_c

  !------------------------------------------------
  ! prg_xx_mod
  !------------------------------------------------



end module prg_c_interface
