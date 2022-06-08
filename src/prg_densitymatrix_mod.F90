!> Module to obtain the density matrix by diagonalizing an orthogonalized Hamiltonian.
!!
!! \ingroup PROGRESS
!!
module prg_densitymatrix_mod

  use bml
  use prg_parallel_mod

  implicit none

  private  !Everything is private by default

  integer, parameter :: dp = kind(1.0d0)

  public :: prg_build_density_T0, prg_check_idempotency, prg_get_eigenvalues
  public :: prg_get_flevel, prg_build_density_T, prg_build_atomic_density
  public :: prg_build_density_T_Fermi, prg_get_flevel_nt, prg_build_density_T_fulldata
  public :: prg_build_density_T_ed, prg_toEigenspace, prg_toCanonicalspace
  public :: Canon_DM_PRT, prg_get_evalsDvalsEvects, prg_build_density_fromEvalsAndEvects

contains

  !> Builds the density matrix from \f$ H_0 \f$ for zero electronic temperature.
  !! \f$ \rho = C \Theta(\mu I - \epsilon) C^{\dagger} \f$
  !! Where,\f$ C \f$ is the matrix eigenvector and \f$ \epsilon \f$ is the matrix eigenvalue.
  !! \f$ \Theta() \f$  is the Heaviside function.
  !! \param ham_bml Input Orthogonalized Hamiltonian matrix.
  !! \param rho_bml Output density matrix.
  !! \param threshold Threshold for sparse matrix algebra.
  !! \param bndfil Filing factor.
  !! \param eigenvalues_out Output the eigenvalues.
  !! \warning This does not solve the generalized eigenvalue problem.
  !! The Hamiltonian that comes in has to be preprg_orthogonalized.
  !!
  subroutine prg_build_density_T0(ham_bml, rho_bml, threshold, bndfil, eigenvalues_out)

    character(20)                      ::  bml_type
    integer                            ::  i, norb
    real(8), intent(in)                ::  bndfil, threshold
    real(dp)                           ::  nocc
    real(dp), allocatable              ::  eigenvalues(:)
    real(dp), allocatable, optional, intent(out)  ::  eigenvalues_out(:)
    type(bml_matrix_t)                 ::  aux1_bml, aux_bml, eigenvectors_bml, occupation_bml
    type(bml_matrix_t), intent(in)     ::  ham_bml
    type(bml_matrix_t), intent(inout)  ::  rho_bml
    character(20) :: bml_dmode

    if (getNRanks().gt.1)then
      if (printRank() .eq. 1) print*,'prg_build_density_T0: BML_DMODE_DISTRIBUTED'
      bml_dmode = BML_DMODE_DISTRIBUTED
    else
      print*,'prg_build_density_T0: BML_DMODE_SEQUENTIAL'
      bml_dmode = BML_DMODE_SEQUENTIAL
    endif

    norb = bml_get_n(ham_bml)
    bml_type = bml_get_deep_type(ham_bml)

    allocate(eigenvalues(nOrb))

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,eigenvectors_bml,bml_dmode)

    call bml_diagonalize(ham_bml,eigenvalues,eigenvectors_bml)
    if(present(eigenvalues_out))then
      if(allocated(eigenvalues_out))deallocate(eigenvalues_out)
      allocate(eigenvalues_out(nOrb))
      eigenvalues_out = eigenvalues
    endif

    nocc = norb*bndfil

    do i=1,norb    !Reusing eigenvalues to apply the theta function.
      if(real(i)-nocc < 0.0001_dp) then
        eigenvalues(i) = 2.0_dp
      elseif(abs(real(i)-real(nocc)) < 0.0001_dp) then
        eigenvalues(i) = 2.0_dp
      else
        eigenvalues(i) = 0.0_dp
      endif
    enddo
    if(abs(nocc - int(nocc)) > 0.01_dp)then
      eigenvalues(int(nocc)+1) = 1.0_dp
    endif

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,occupation_bml,bml_dmode)
    call bml_set_diagonal(occupation_bml, eigenvalues) !eps(i,i) = eps(i)

    deallocate(eigenvalues)

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,aux_bml,bml_dmode)
    call bml_multiply(eigenvectors_bml, occupation_bml, aux_bml, 1.0_dp, 0.0_dp, threshold)
    call bml_deallocate(occupation_bml)

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,aux1_bml,bml_dmode)
    call bml_transpose(eigenvectors_bml, aux1_bml)
    call bml_deallocate(eigenvectors_bml)

    call bml_multiply(aux_bml, aux1_bml, rho_bml, 1.0_dp, 0.0_dp, threshold)
    call bml_deallocate(aux_bml)
    call bml_deallocate(aux1_bml)

  end subroutine prg_build_density_T0

  !> Builds the density matrix from \f$ H_0 \f$ for electronic temperature T.
  !! \f$ \rho = C f(\mu I - \epsilon) C^{\dagger} \f$
  !! Where,\f$ C \f$ is the matrix eigenvector and \f$ \epsilon \f$ is the matrix eigenvalue.
  !! \f$ f \f$  is the Fermi function.
  !! \param ham_bml Input Orthogonalized Hamiltonian matrix.
  !! \param rho_bml Output density matrix,
  !! \param threshold Threshold for sparse matrix algebra.
  !! \param bndfil Filing factor.
  !! \param kbt Electronic temperature.
  !! \param ef Fermi level.
  !! \param eigenvalues_out Output the eigenvalues.
  !! \warning This does not solve the generalized eigenvalue problem.
  !! The Hamiltonian that comes in has to be preorthogonalized.
  !!
  subroutine prg_build_density_T(ham_bml, rho_bml, threshold, bndfil, kbt, ef, eigenvalues_out)

    character(20)                      ::  bml_type
    integer                            ::  i, norb
    real(dp), intent(in)               ::  bndfil, threshold, kbt
    real(dp), intent(inout)            ::  ef
    real(dp)                           ::  nocc, fleveltol, mlsi, efOld
    real(dp), allocatable              ::  eigenvalues(:)
    real(dp), allocatable, optional, intent(out)  ::  eigenvalues_out(:)
    type(bml_matrix_t)                 ::  aux1_bml, aux_bml, eigenvectors_bml, occupation_bml
    type(bml_matrix_t), intent(in)     ::  ham_bml
    type(bml_matrix_t), intent(inout)  ::  rho_bml
    logical :: err

    if (printRank() .eq. 1) then
      write(*,*)"In get_density_t ..."
    endif

    norb = bml_get_n(ham_bml)
    bml_type = bml_get_deep_type(ham_bml)

    allocate(eigenvalues(nOrb))
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,eigenvectors_bml)

    call bml_diagonalize(ham_bml,eigenvalues,eigenvectors_bml)

    if(present(eigenvalues_out))then
      if(allocated(eigenvalues_out))deallocate(eigenvalues_out)
      allocate(eigenvalues_out(nOrb))
      eigenvalues_out = eigenvalues
    endif

    fleveltol = 1.0e-11

    efOld = ef
    call prg_get_flevel_nt(eigenvalues,kbt,bndfil,fleveltol,ef,err)
    if(err)call prg_get_flevel(eigenvalues,kbt,bndfil,fleveltol,ef,err)
    if(err)then
      write(*,*)"WARNING: Ef/Chemical potential search failed. We'll use the previous one to proceed"
      ef = efOld
    endif

    nocc = norb*bndfil

    do i=1,norb   !Reusing eigenvalues to apply the theta function.
      eigenvalues(i) = 2.0_dp*fermi(eigenvalues(i),ef,kbt)
    enddo

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,occupation_bml)
    call bml_set_diagonal(occupation_bml, eigenvalues) !eps(i,i) = eps(i)

    deallocate(eigenvalues)

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,aux_bml)
    call bml_multiply(eigenvectors_bml, occupation_bml, aux_bml, 1.0_dp, 0.0_dp, threshold)
    call bml_deallocate(occupation_bml)

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,aux1_bml)
    call bml_transpose(eigenvectors_bml, aux1_bml)
    call bml_deallocate(eigenvectors_bml)

    call bml_multiply(aux_bml, aux1_bml, rho_bml, 1.0_dp, 0.0_dp, threshold)
    call bml_deallocate(aux_bml)
    call bml_deallocate(aux1_bml)

  end subroutine prg_build_density_T


  !> Builds the density matrix from \f$ H_0 \f$ for electronic temperature T.
  !! \f$ \rho = C f(\mu I - \epsilon) C^{\dagger} \f$
  !! Where,\f$ C \f$ is the matrix eigenvector and \f$ \epsilon \f$ is the matrix eigenvalue.
  !! \f$ f \f$  is the Fermi function.
  !! \param ham_bml Input Orthogonalized Hamiltonian matrix.
  !! \param rho_bml Output density matrix,
  !! \param threshold Threshold for sparse matrix algebra.
  !! \param bndfil Filing factor.
  !! \param kbt Electronic temperature.
  !! \param ef Fermi level.
  !! \param eigenvalues_out Output the eigenvalues.
  !! \param evects_bml Output the eigenvectors.
  !! \param fvals Output the occupancies.
  !! \warning This does not solve the generalized eigenvalue problem.
  !! The Hamiltonian that comes in has to be preorthogonalized.
  !!
  subroutine prg_build_density_T_fulldata(ham_bml, rho_bml, threshold, bndfil, kbt, ef, eigenvalues_out&
       &, evects_bml, fvals)

    character(20)                      ::  bml_type
    integer                            ::  i, norb
    real(dp), intent(in)               ::  bndfil, threshold, kbt
    real(dp), intent(inout)            ::  ef
    real(dp)                           ::  nocc, fleveltol, efOld
    real(dp), allocatable              ::  eigenvalues(:)
    real(dp), allocatable, intent(inout)  ::  eigenvalues_out(:), fvals(:)
    type(bml_matrix_t)                 ::  aux1_bml, aux_bml, occupation_bml
    type(bml_matrix_t), intent(in)     ::  ham_bml
    type(bml_matrix_t), intent(inout)  ::  rho_bml, evects_bml
    logical :: err

    if (printRank() .eq. 1) then
      write(*,*)"In get_density_t_fulldata ..."
    endif

    norb = bml_get_n(ham_bml)
    bml_type = bml_get_deep_type(ham_bml)

    allocate(eigenvalues(nOrb))
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,evects_bml)

    call bml_diagonalize(ham_bml,eigenvalues,evects_bml)

    if(.not.allocated(eigenvalues_out))then
      allocate(eigenvalues_out(nOrb))
      allocate(fvals(nOrb))
      call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,evects_bml)
    endif

    eigenvalues_out = eigenvalues

    fleveltol = 1.0e-11
    fvals = 0.0_dp

    call prg_get_flevel_nt(eigenvalues,kbt,bndfil,fleveltol,ef,err)
    if(err)call prg_get_flevel(eigenvalues,kbt,bndfil,fleveltol,ef,err)
    if(err)then
      write(*,*)"WARNING: Ef/Chemical potential search failed. We'll use the previous one to proceed"
      ef = efOld
    endif

    nocc = norb*bndfil

    do i=1,norb   !Aapply Fermi function.
      fvals(i) = 2.0_dp*fermi(eigenvalues(i),ef,kbt)
    enddo

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,occupation_bml)
    call bml_set_diagonal(occupation_bml, fvals) !eps(i,i) = eps(i)

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,aux_bml)
    call bml_multiply(evects_bml, occupation_bml, aux_bml, 1.0_dp, 0.0_dp, threshold)
    call bml_deallocate(occupation_bml)

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,aux1_bml)
    call bml_transpose(evects_bml, aux1_bml)

    call bml_multiply(aux_bml, aux1_bml, rho_bml, 1.0_dp, 0.0_dp, threshold)
    call bml_deallocate(aux_bml)
    call bml_deallocate(aux1_bml)

  end subroutine prg_build_density_T_fulldata


  !> Builds the density matrix from \f$ H_0 \f$ for electronic temperature T.
  !! \f$ \rho = C f(\mu I - \epsilon) C^{\dagger} \f$
  !! Where,\f$ C \f$ is the matrix eigenvector and \f$ \epsilon \f$ is the matrix eigenvalue.
  !! \f$ f \f$  is the Fermi function.
  !! \param ham_bml Input Orthogonalized Hamiltonian matrix.
  !! \param rho_bml Output density matrix,
  !! \param threshold Threshold for sparse matrix algebra.
  !! \param bndfil Filing factor.
  !! \param kbt Electronic temperature.
  !! \param ef Fermi level.
  !! \param evals Output the eigenvalues.
  !! \param fvals Output the occupancies.
  !! \param dvals Contribution to population from every evect.
  !! \param core_indices Indices in the core.
  !! \warning This does not solve the generalized eigenvalue problem.
  !! The Hamiltonian that comes in has to be orthogonalized.
  !!
  subroutine prg_build_density_T_ed(ham_bml, rho_bml, evects_bml, threshold, bndfil, kbt, ef, evals&
       &, dvals, hindex, llsize, verbose)

    character(20)                      ::  bml_type
    integer, allocatable, intent(in)   ::  hindex(:,:)
    integer                            ::  i, j, jj, k, kk, l, norb, norbCore
    integer, intent(in)                ::  llsize, verbose
    real(dp), intent(in)               ::  bndfil, threshold, kbt
    real(dp), intent(inout)            ::  ef
    real(dp)                           ::  nocc
    real(dp), allocatable              ::  eigenvalues(:), row(:),fvals(:),ham(:,:),evects(:,:)
    real(dp), allocatable, intent(inout)  ::  evals(:), dvals(:)
    type(bml_matrix_t)                 ::  aux1_bml, aux_bml, occupation_bml
    type(bml_matrix_t), intent(in)     ::  ham_bml
    type(bml_matrix_t), intent(inout)  ::  rho_bml, evects_bml

    if (printRank() .eq. 1 .and. verbose >= 1) then
      write(*,*)"In get_density_t_efd ..."
    endif

    norb = bml_get_n(ham_bml)
    bml_type = bml_get_type(ham_bml)

    if(allocated(evals))then
      deallocate(evals)
      deallocate(dvals)
    endif
    allocate(evals(norb))
    allocate(dvals(norb))
    allocate(fvals(norb))
    allocate(ham(norb,norb))
    allocate(evects(norb,norb))
    call bml_export_to_dense(ham_bml,ham)
    if(bml_get_n(evects_bml) < 0)then
      call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,evects_bml)
    endif
    call Eig(ham,evects,evals,'V',norb)
    call bml_import_from_dense(bml_type,evects,evects_bml,0.0_dp,norb)
    !call bml_diagonalize(ham_bml,evals,evects_bml)
    write(*,*)"ham",ham(1,1),ham(1,2),ham(1,3)
    deallocate(evects)
    deallocate(ham)
    call bml_print_matrix("evects",evects_bml,0,4,0,4)
    write(*,*)"Evals",evals(1),evals(2),evals(3)
    nocc = norb*bndfil

    do i=1,norb   !Apply Fermi function.
      fvals(i) = 2.0_dp*fermi(evals(i),ef,kbt)
    enddo

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,occupation_bml)
    call bml_set_diagonal(occupation_bml, fvals) !eps(i,i) = eps(i)

    deallocate(fvals)

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,aux_bml)
    call bml_multiply(evects_bml, occupation_bml, aux_bml, 1.0_dp, 0.0_dp, threshold)
    call bml_deallocate(occupation_bml)

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,aux1_bml)
    call bml_transpose(evects_bml, aux1_bml)

    allocate(row(norb))
    dvals = 0.0_dp
    do i = 1,norb
      call bml_get_row(aux1_bml,i,row)
      do k = 1,llsize
        do l = hindex(1,k),hindex(2,k)
          dvals(i) = dvals(i) + row(l)**2
        enddo
      enddo
    enddo
    deallocate(row)

    call bml_multiply(aux_bml, aux1_bml, rho_bml, 1.0_dp, 0.0_dp, threshold)

    call bml_deallocate(aux_bml)
    call bml_deallocate(aux1_bml)

  end subroutine prg_build_density_T_ed


  !> Gets the eigenvalues and eigenvectors and the core contribution to each
  !! eigenvalue (dvals) .
  !! \param ham_bml Input Orthogonalized Hamiltonian matrix.
  !! \param threshold Threshold for sparse matrix algebra.
  !! \param hindex Start and end index for every atom in the system.
  !! \param llsize Number of atoms in the core of the graph part.
  !! \param evals Eigenvalues.
  !! \param dvals Contribution to population from every evect.
  !! \param evects_bml Output the eigenvectors.
  !! \warning This does not solve the generalized eigenvalue problem.
  !! The Hamiltonian that comes in has to be orthogonalized.
  !!
  subroutine prg_get_evalsDvalsEvects(ham_bml, threshold, hindex,&
       & llsize, evals, dvals, evects_bml, verbose)

    character(20)                      ::  bml_type
    integer, allocatable, intent(in)   ::  hindex(:,:)
    integer                            ::  i, k, l, norb
    integer, intent(in)                ::  llsize, verbose
    real(dp), intent(in)               ::  threshold
    real(dp), allocatable              ::  row(:),ham(:,:),evects(:,:)
    real(dp), allocatable, intent(inout)  ::  evals(:), dvals(:)
    type(bml_matrix_t), intent(in)     ::  ham_bml
    type(bml_matrix_t), intent(inout)  ::  evects_bml
    type(bml_matrix_t)                 ::  aux1_bml

    if (printRank() .eq. 1 .and. verbose >= 1) then
      write(*,*)"In get_evalsDvalsEvects ..."
    endif

    norb = bml_get_n(ham_bml)
    bml_type = bml_get_type(ham_bml)

    if(allocated(evals))then
      deallocate(evals)
      deallocate(dvals)
    endif
    allocate(evals(norb))
    allocate(dvals(norb))
    allocate(ham(norb,norb))
    allocate(evects(norb,norb))
    call bml_export_to_dense(ham_bml,ham)
    if(bml_get_n(evects_bml) < 0)then
      call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,evects_bml)
    endif
    call Eig(ham,evects,evals,'V',norb)
    call bml_import_from_dense(bml_type,evects,evects_bml,0.0_dp,norb)
    !call bml_diagonalize(ham_bml,evals,evects_bml)
    deallocate(evects)
    deallocate(ham)

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,aux1_bml)
    call bml_transpose(evects_bml, aux1_bml)

    allocate(row(norb))
    dvals = 0.0_dp
    do i = 1,norb
      call bml_get_row(aux1_bml,i,row)
      do k = 1,llsize
        do l = hindex(1,k),hindex(2,k)
          dvals(i) = dvals(i) + row(l)**2
        enddo
      enddo
    enddo
    deallocate(row)
  end subroutine prg_get_evalsDvalsEvects


  !> Builds the density matrix from the evects and evals for electronic temperature T.
  !! \param evects_bml Eigenvectors.
  !! \param evals Eigenvalues.
  !! \param threshold Threshold for sparse matrix algebra.
  !! \param bndfil Filing factor.
  !! \param kbt Electronic temperature.
  !! \param ef Fermi level.
  !! \param verbose Verbosity level.
  !!
  subroutine prg_build_density_fromEvalsAndEvects(evects_bml, evals, rho_bml, threshold,&
       & bndfil, kbt, ef, verbose)

    character(20)                      ::  bml_type
    integer                            ::  i, norb
    integer, intent(in)                ::  verbose
    real(dp), intent(in)               ::  bndfil, threshold, kbt
    real(dp), intent(inout)            ::  ef
    real(dp)                           ::  nocc
    real(dp), intent(in)               ::  evals(*)
    real(dp), allocatable              ::  fvals(:)
    type(bml_matrix_t)                 ::  aux1_bml, aux_bml, occupation_bml
    type(bml_matrix_t), intent(in)     ::  evects_bml
    type(bml_matrix_t), intent(inout)  ::  rho_bml

    if (printRank() .eq. 1 .and. verbose >= 1) then
      write(*,*)"In get_density_fromEvalsAndEvects ..."
    endif

    norb = bml_get_N(evects_bml)
    bml_type = bml_get_type(evects_bml)

    allocate(fvals(norb))

    nocc = norb*bndfil

    do i=1,norb   !Apply Fermi function.
      fvals(i) = 2.0_dp*fermi(evals(i),ef,kbt)
    enddo

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,occupation_bml)
    call bml_set_diagonal(occupation_bml, fvals) !eps(i,i) = eps(i)

    deallocate(fvals)

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,aux_bml)
    call bml_multiply(evects_bml, occupation_bml, aux_bml, 1.0_dp, 0.0_dp,threshold)
    call bml_deallocate(occupation_bml)

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,aux1_bml)
    call bml_transpose(evects_bml, aux1_bml)

    call bml_multiply(aux_bml, aux1_bml, rho_bml, 1.0_dp, 0.0_dp, threshold)

    call bml_deallocate(aux_bml)
    call bml_deallocate(aux1_bml)

  end subroutine prg_build_density_fromEvalsAndEvects



  !> Builds the density matrix from \f$ H_0 \f$ for electronic temperature T.
  !! \f$ \rho = C f(\mu I - \epsilon) C^{\dagger} \f$
  !! Where,\f$ C \f$ is the matrix eigenvector and \f$ \epsilon \f$ is the matrix eigenvalue.
  !! \f$ f \f$  is the Fermi function. In this routine the Fermi level is passed as an argument.
  !! \param ham_bml Input Orthogonalized Hamiltonian matrix.
  !! \param rho_bml Output density matrix,
  !! \param threshold Threshold for sparse matrix algebra.
  !! \warning This does not solve the generalized eigenvalue problem.
  !! The Hamiltonian that comes in has to be preorthogonalized.
  !!
  subroutine prg_build_density_T_Fermi(ham_bml, rho_bml, threshold, kbt, ef, verbose, drho)

    character(20)                      ::  bml_type
    integer                            ::  i, norb, mdim
    integer, optional, intent(in)      ::  verbose
    real(dp), intent(in)               ::  threshold, kbt
    real(dp), intent(in)               ::  ef
    real(dp)                           ::  fleveltol
    real(dp), allocatable              ::  eigenvalues(:)
    real(dp), allocatable              ::  nocc(:)
    type(bml_matrix_t)                 ::  aux1_bml, aux_bml, eigenvectors_bml, occupation_bml
    type(bml_matrix_t), intent(in)     ::  ham_bml
    type(bml_matrix_t), intent(inout)  ::  rho_bml
    real(dp), optional, intent(inout)  ::  drho ! gradient of density w.r.t. ef

    if (printRank() .eq. 1) then
      if(present(verbose))then
        if(verbose >= 1)then
          write(*,*)"In get_density_t_Fermi ..."
          write(*,*)"Ef = ",ef
        endif
      endif
    endif

    norb = bml_get_n(ham_bml)
    mdim = bml_get_M(ham_bml)
    bml_type = bml_get_deep_type(ham_bml)

    allocate(eigenvalues(nOrb),nocc(norb))
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,mdim,eigenvectors_bml)

    call bml_diagonalize(ham_bml,eigenvalues,eigenvectors_bml)

    do i=1,norb   !Reusing eigenvalues to apply the theta function.
      nocc(i) = 2.0_dp*fermi(eigenvalues(i),ef,kbt)
    enddo

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,mdim,occupation_bml)
    call bml_set_diagonal(occupation_bml, nocc) !eps(i,i) = eps(i)

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,mdim,aux_bml)
    call bml_multiply(eigenvectors_bml, occupation_bml, aux_bml, 1.0_dp, 0.0_dp, threshold)
    call bml_deallocate(occupation_bml)

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,mdim,aux1_bml)
    call bml_transpose(eigenvectors_bml, aux1_bml)

    call bml_multiply(aux_bml, aux1_bml, rho_bml, 1.0_dp, 0.0_dp, threshold)

    if (present(drho)) then
      do i=1,norb   !Reusing eigenvalues to apply the theta function.
        nocc(i) = 2.0_dp * dfermi(eigenvalues(i),ef,kbt)
      enddo

      call bml_zero_matrix(bml_type,bml_element_real,dp,norb,mdim,occupation_bml)
      call bml_set_diagonal(occupation_bml, nocc) !eps(i,i) = eps(i)

      call bml_multiply(eigenvectors_bml, occupation_bml, aux_bml, 1.0_dp, 0.0_dp, threshold)

      call bml_zero_matrix(bml_type,bml_element_real,dp,norb,mdim,aux1_bml)
      call bml_transpose(eigenvectors_bml, aux1_bml)

      call bml_multiply(aux_bml, aux1_bml, occupation_bml, 1.0_dp, 0.0_dp, threshold)
      drho = bml_trace(occupation_bml)
      call bml_deallocate(occupation_bml)
    endif

    deallocate(nocc)
    deallocate(eigenvalues)
    call bml_deallocate(eigenvectors_bml)
    call bml_deallocate(aux_bml)
    call bml_deallocate(aux1_bml)

  end subroutine prg_build_density_T_Fermi

  !> Builds the atomic density matrix.
  !! \f$ \rho_{ii} = mathcal{Z}_{ii} \f$
  !! Where,\f$ mathcal{Z}_{ii} \f$ is the number of electrons for orbital i.
  !! \param rhoat Output atomic diagonal density matrix,
  !! \param hindex Start and end index for every atom in the system.
  !! \param numel Number of electrons per specie. It runs over the specie index.
  !! \param spindex Specie index.
  !! \param norbs Number of orbitals.
  !!
  subroutine prg_build_atomic_density(rhoat_bml,numel,hindex,spindex,norb,bml_type)

    character(len=*), intent(in)          ::  bml_type
    integer                            ::  i, index, n_orb, nats
    integer, intent(in)                ::  hindex(:,:), norb, spindex(:)
    real(dp)                           ::  occ
    real(dp), allocatable              ::  d_atomic(:), rhoat(:)
    real(dp), intent(in)               ::  numel(:)
    type(bml_matrix_t), intent(inout)  ::  rhoat_bml

    nats = size(hindex,dim=2)

    allocate(rhoat(norb))

    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,rhoat_bml)

    index = 0;

    do i = 1,nats
      n_orb = hindex(2,i)-hindex(1,i) + 1;
      if(n_orb == 1)then
        index = index + 1;
        rhoat(index) = numel(spindex(i));
      else
        if(numel(spindex(i)) <= 2)then
          index = index + 1;
          rhoat(index) = numel(spindex(i));
          index = index + 1;
          rhoat(index) = 0.0_dp;
          index = index + 1;
          rhoat(index) = 0.0_dp;
          index = index + 1;
          rhoat(index) = 0.0_dp;
        else
          index = index + 1;
          rhoat(index) = 2.0_dp;

          index = index + 1;
          occ = (numel(spindex(i))-2.0_dp)/3.0_dp;
          rhoat(index) = occ;
          index = index + 1;
          rhoat(index) = occ;
          index = index + 1;
          rhoat(index) = occ;
        endif
      endif
    enddo

    call bml_set_diagonal(rhoat_bml,rhoat,0.0_dp)

    deallocate(rhoat)

  end subroutine prg_build_atomic_density

  !> Routine to compute the Fermi level given a set of eigenvalues and a temperature.
  !! It applies the Bisection method over the function:
  !! \f$ g(\mu) = \sum_k 2 f(\epsilon_k - \mu) - N = 0 \f$
  !! Where \f$ f(\epsilon_k - \mu) = \frac{1}{1+\exp{(\epsilon_k - \mu)/(k_bT)}}\f$.
  !! \param eigenvalues Eigenvalues of the system (\f$ \{ \epsilon_k \} \f$).
  !! \param kbt Temperature times the Boltzmann's constant  (\f$ k_bT  \f$).
  !! \param bndfil Filing factor (\f$ N_{el}/(2*N_{orbs})\f$).
  !! \param tol Tolerance for the bisection method.
  !! \param Ef Fermi level (\f$ \mu \f$).
  !! \param err Error logical variable
  subroutine prg_get_flevel(eigenvalues,kbt,bndfil,tol,Ef,err)

    integer                  ::  i, j, k, m
    integer                  ::  norb
    real(dp)                 ::  Ft1, Ft2, Kb, Prod
    real(dp)                 ::  T, step, tol, nel
    real(dp), intent(in)     ::  bndfil, eigenvalues(:), kbt
    real(dp), intent(inout)  ::  Ef
    logical, intent(inout)   ::  err

    norb = size(eigenvalues,dim=1)
    nel = bndfil*2.0_dp*norb
    Ef=eigenvalues(1)
    step=abs((eigenvalues(norb)-eigenvalues(1)))
    Ft1=0.0_dp
    Ft2=0.0_dp
    prod=0.0_dp

    !Sum of the occupations
    do i=1,norb
      ft1 = ft1 + 2.0_dp*fermi(eigenvalues(i),ef,kbt)
    enddo
    ft1=ft1-nel

    err = .false.
    do m=1,1000001

      if(m.gt.1000000)then
        err = .true.
        write(*,*) "WARNING: Bisection method in prg_get_flevel not converging ..."
      endif

      if(abs(ft1).lt.tol)then !tolerance control
        return
      endif

      ef = ef + step

      ft2=0.0_dp

      !New sum of the occupations
      do i=1,norb
        ft2 = ft2 + 2.0_dp*fermi(eigenvalues(i),ef,kbt)
      enddo

      ft2=ft2-nel

      !Product to see the change in sign.
      prod = ft2*ft1

      if(prod.lt.0)then
        ef=ef-step
        step=step/2.0_dp !If the root is inside we shorten the step.
      else
        ft1=ft2  !If not, Ef moves forward.
      endif

    enddo

  end subroutine prg_get_flevel

  !> Routine to compute the Fermi level given a set of eigenvalues and a temperature.
  !! It applies the Newton-Raphson method over the function:
  !! \f$ g(\mu) = \sum_k 2 f(\epsilon_k - \mu) - N = 0 \f$
  !! Where \f$ f(\epsilon_k - \mu) = \frac{1}{1+\exp{(\epsilon_k - \mu)/(k_bT)}}\f$.
  !! \param eigenvalues Eigenvalues of the system (\f$ \{ \epsilon_k \} \f$).
  !! \param kbt Temperature times the Boltzmann's constant  (\f$ k_bT  \f$).
  !! \param bndfil Filing factor (\f$ N_{el}/(2*N_{orbs})\f$).
  !! \param tol Tolerance for the bisection method.
  !! \param Ef Fermi level (\f$ \mu \f$).
  !! \param err Error logical variable
  !!
  subroutine prg_get_flevel_nt(eigenvalues,kbt,bndfil,tol,ef,err,verbose)

    integer                  ::  i, m
    integer                  ::  norb
    real(dp)                 ::  ef0, f1, f2, nel
    real(dp)                 ::  step
    real(dp), intent(in)     ::  bndfil, kbt
    real(dp), intent(in)     ::  tol, eigenvalues(:)
    real(dp), intent(inout)  ::  ef
    integer, optional, intent(in) :: verbose
    logical, intent(inout) :: err

    norb = size(eigenvalues, dim=1)
    nel = bndfil*2.0_dp*real(norb,dp)
    ef0 = ef
    f1 = 0.0_dp
    f2 = 0.0_dp
    step = 0.1_dp
    err = .false.

    !$omp parallel do default(none) private(i) &
    !$omp shared(eigenvalues,kbt,ef,norb) &
    !$omp reduction(+:f1)
    do i=1,norb
      f1 = f1 + 2.0_dp*fermi(eigenvalues(i),ef,kbt)
    enddo
    !$omp end parallel do

    f1=f1-nel
    ef = ef0 + step

    f2 = 0.0_dp

    !$omp parallel do default(none) private(i) &
    !$omp shared(eigenvalues,kbt,ef,norb) &
    !$omp reduction(+:f2)
    do i=1,norb
      f2 = f2 + 2.0_dp*fermi(eigenvalues(i),ef,kbt)
    enddo
    !$omp end parallel do

    f2=f2-nel
    ef0 = ef
    ef = -f2*step/(f2-f1) + ef0
    f1 = f2
    step = ef - ef0

    do m = 1,1000001
      if(m.gt.1000000)then
        write(*,*) "WARNING: Newton method in prg_get_flevel_nt is not converging ..."
        err = .true.
        exit
      endif

      !New sum of the occupations
      f2 = 0.0_dp
      !$omp parallel do default(none) private(i) &
      !$omp shared(eigenvalues,ef,kbt,norb) &
      !$omp reduction(+:f2)
      do i=1,norb
        f2 = f2 +  2.0_dp*fermi(eigenvalues(i),ef,kbt)
      enddo
      !$omp end parallel do

      f2=f2-nel
      ef0 = ef
      ef = -f2*step/(f2-f1) + ef0
      f1 = f2
      step = ef - ef0
      if(abs(f1).lt.tol)then !tolerance control
        return
      endif
    enddo

  end subroutine prg_get_flevel_nt


  !> Gets the eigenvalues of the Orthogonalized Hamiltonian.
  !! \param ham_bml Input Orthogonalized Hamiltonian matrix.
  !! \param eigenvalues Output eigenvalues of the system.
  !! \param verbose Verbosity level.
  !!
  subroutine prg_get_eigenvalues(ham_bml,eigenvalues,verbose)

    character(20)                        ::  bml_type
    integer                              ::  i, norb
    integer, intent(in)                  ::  verbose
    real(dp)                             ::  nocc
    real(dp), allocatable                ::  aux(:,:)
    real(dp), allocatable, intent(inout)  ::  eigenvalues(:)
    type(bml_matrix_t)                   ::  aux_bml, eigenvectors_bml
    type(bml_matrix_t), intent(in)       ::  ham_bml

    norb = bml_get_n(ham_bml)
    bml_type = bml_get_deep_type(ham_bml)

    if(verbose.ge.1)write(*,*)"In prg_get_eigenvalues ..."
    if(verbose.ge.1)write(*,*)"Number of states =",norb

    if(.not.allocated(eigenvalues))allocate(eigenvalues(norb))

    !Ensure dense type to diagonalize
    if(trim(bml_type).ne."dense")then
      allocate(aux(norb,norb))
      call bml_export_to_dense(ham_bml,aux)
      call bml_zero_matrix(bml_matrix_dense,bml_element_real,dp,norb,norb,aux_bml)
      call bml_import_from_dense(bml_matrix_dense,aux,aux_bml,0.0_dp,norb)
      deallocate(aux)
    else
      call bml_copy_new(ham_bml,aux_bml)
    endif

    call bml_zero_matrix(bml_matrix_dense,bml_element_real,dp,norb,norb,eigenvectors_bml)

    call bml_diagonalize(aux_bml,eigenvalues,eigenvectors_bml)

    call bml_deallocate(eigenvectors_bml)
    call bml_deallocate(aux_bml)

  end subroutine prg_get_eigenvalues


  !> To check the idempotency error of a matrix.
  !! This is calculated as the Frobenius norm of \f$ (A - A^2) \f$
  !! \param mat_bml Some bml matrix
  !! \param idempotency (Output value of the idempotency error)
  !!
  subroutine prg_check_idempotency(mat_bml, threshold, idempotency)

    character(20)                   ::  bml_type
    integer                         ::  N, i, j
    real(dp), intent(in)            ::  threshold
    real(dp), intent(out)           ::  idempotency
    type(bml_matrix_t)              ::  aux_bml
    type(bml_matrix_t), intent(in)  ::  mat_bml

    call bml_copy_new(mat_bml, aux_bml)

    call bml_multiply(mat_bml, mat_bml, aux_bml, 1.0_dp, 0.0_dp, threshold) !A^2
    call bml_add_deprecated(1.0_dp,aux_bml,-1.0_dp,mat_bml, threshold) !A^2 - A

    idempotency = bml_fnorm(aux_bml)

    call bml_deallocate(aux_bml)

  end subroutine prg_check_idempotency

  !> Gives the Fermi distribution value for energy e.
  !! \param e Energy.
  !! \param ef Fermi energy.
  !!
  real(dp) function fermi(e,ef,kbt)

    real(dp), intent(in) :: e, ef, kbt

    if ((e-ef)/kbt > 100.0_dp) then
      fermi = 0.0_dp
    else
      fermi = 1.0_dp/(1.0_dp+exp((e-ef)/(kbt)))
    endif

  end function fermi

  !> Gives the gradient of Fermi distribution w.r.t. ef.
  !! \param e Energy.
  !! \param ef Fermi energy.
  !!
  real(dp) function dfermi(e,ef,kbt)

    real(dp), intent(in) :: e, ef, kbt

    dfermi = 0.0_dp
    if ((e-ef)/kbt > 100.0_dp) then
      if (abs(e-ef)<0.001_dp) then
        dfermi = - 1000.0_dp
      endif
    else
      dfermi = 1.0_dp/(1.0_dp+exp((e-ef)/(kbt)))
      dfermi = dfermi * dfermi * exp((e-ef)/(kbt)) * 1.0_dp/kbt
    endif

  end function dfermi

  !> Change an operator into the eigenspace representation
  !! \brief This routine performs
  !! \f$ O_{eig} = U^{T} O U \f$, where operator U is the unitary transformations
  !! constructed from the eigenvectors of a Hamiltonian
  !! \param mat_bml Operator to be transformed
  !! \param matEig_bml Output operator after transformation
  !! \param evects_bml Eigenvectors matrix (U)
  !! \param threshold Threshold parameter
  !! \verbose Verbosity level
  !!
  subroutine prg_toEigenspace(mat_bml,matEig_bml,evects_bml,threshold,verbose)

    integer, optional, intent(in)   :: verbose
    type(bml_matrix_t), intent(in)  :: mat_bml, evects_bml
    type(bml_matrix_t), intent(inout) :: matEig_bml
    type(bml_matrix_t)              ::  aux_bml
    real(dp), intent(in)            :: threshold
    integer                         :: hdim, mdim
    character(20)                   :: bml_type

    if(verbose.eq.1) write(*,*)"In prg_toEigenspace ..."

    hdim= bml_get_N(mat_bml)
    mdim= bml_get_M(mat_bml)
    bml_type = bml_get_type(mat_bml)
    !Allocate bml's
    if(bml_get_N(matEig_bml) < 0)then
      call bml_zero_matrix(bml_type,bml_element_real,dp,hdim ,mdim,matEig_bml, &
           & bml_get_distribution_mode(mat_bml))
    endif
    !Do the operations in bml
    call bml_transpose(evects_bml, aux_bml)

    call bml_multiply(aux_bml, mat_bml, matEig_bml, 1.0_dp, 0.0_dp,threshold) !U^t*O

    call bml_multiply(matEig_bml, evects_bml, aux_bml, 1.0_dp, 0.0_dp,threshold) !U^t*O * U

    call bml_copy(aux_bml, matEig_bml)

    call bml_deallocate(aux_bml)

  end subroutine prg_toEigenspace

  !> Change an operator into the eigenspace representation
  !! \brief This routine performs
  !! \f$ O_{eig} = U O U^{t} \f$, where operator U is the unitary
  !transformations
  !! constructed from the eigenvectors of a Hamiltonian
  !! \param mat_bml Operator to be transformed
  !! \param matEig_bml Output operator after transformation
  !! \param evects_bml Eigenvectors matrix (U)
  !! \param threshold Threshold parameter
  !! \verbose Verbosity level
  !!
  subroutine prg_toCanonicalspace(mat_bml,matCan_bml,evects_bml,threshold,verbose)

    integer, optional, intent(in)   :: verbose
    type(bml_matrix_t), intent(in)  :: mat_bml, evects_bml
    type(bml_matrix_t), intent(inout) :: matCan_bml
    type(bml_matrix_t)              ::  aux_bml
    real(dp), intent(in)            :: threshold
    integer                         :: hdim, mdim
    character(20)                   :: bml_type

    if(verbose.eq.1) write(*,*)"In prg_toCacninicalspace ..."

    hdim= bml_get_N(mat_bml)
    mdim= bml_get_M(mat_bml)
    bml_type = bml_get_type(mat_bml)

    !Allocate bml's
    if(bml_get_N(matCan_bml) < 0)then
      call bml_zero_matrix(bml_type,bml_element_real,dp,hdim ,mdim,matCan_bml,&
           & bml_get_distribution_mode(mat_bml))
    endif

    !Do the operations in bml

    call bml_transpose(evects_bml, aux_bml)
    call bml_multiply(mat_bml, aux_bml, matCan_bml, 1.0_dp, 0.0_dp,threshold) !O*U^t

    call bml_multiply(evects_bml,matCan_bml, aux_bml, 1.0_dp, 0.0_dp,threshold) !U*O * U^t

    call bml_copy(aux_bml, matCan_bml)

    call bml_deallocate(aux_bml)

  end subroutine prg_toCanonicalspace


  subroutine  Canon_DM_PRT(P1,H1,Nocc,T,Q,e,mu0,m,HDIM)

    implicit none
    integer, parameter      :: PREC = 8
    real(PREC), parameter   :: ONE = 1.D0, TWO = 2.D0, ZERO = 0.D0
    integer, intent(in)     :: HDIM, m ! Nocc = Number of occupied orbitals, m
    real(PREC), intent(in)  :: H1(HDIM,HDIM), Q(HDIM,HDIM), e(HDIM), nocc !
    real(PREC), intent(out) :: P1(HDIM,HDIM) ! Density matrix response derivative
    real(PREC), intent(in)  :: T, mu0 ! Electronic temperature and chemicalpotential
    real(PREC)              :: X(HDIM,HDIM), DX1(HDIM,HDIM), Y(HDIM,HDIM) !Temporary matrices
    real(PREC)              :: h_0(HDIM), p_0(HDIM), dPdmu(HDIM), p_02(HDIM),iD0(HDIM)
    real(PREC)              :: beta, cnst, kB, mu1, sumdPdmu
    integer                 :: i, j, k, norbsCore

    kB = 8.61739e-5        ! (eV/K)
    beta = 1.D0/(kB*T)     ! Temp in Kelvin
    h_0 = e                ! Diagonal Hamiltonian H0 respresented in the
    cnst = beta/(1.D0*2**(m+2)) ! Scaling constant
    p_0 = 0.5D0 + cnst*(h_0-mu0)  ! Initialization for P0 represented in

    !  call MMult(ONE,Q,H1,ZERO,X,'T','N',HDIM)      ! Main cost are the
    !  transformation
    !  call MMult(ONE,X,Q,ZERO,Y,'N','N',HDIM)       ! to the eigenbasis of H1
    !  P1 = -cnst*Y    !(set mu1 = 0 for simplicity) ! Initialization of DM response
    !  (not diagonal in Q)
    P1 = -cnst*H1    !(set mu1 = 0 for simplicity) ! Initialization of DM response


    do i = 1,m  ! Loop over m recursion steps
      p_02 = p_0*p_0
      do j = 1,HDIM
        do k = 1,HDIM
          DX1(k,j) = p_0(k)*P1(k,j) + P1(k,j)*p_0(j)
        enddo
      enddo
      iD0 = 1.D0/(2.D0*(p_02-p_0)+1.D0)
      p_0 = iD0*p_02
      do j = 1,HDIM
        do k = 1,HDIM
          P1(k,j) = iD0(k)*(DX1(k,j) + 2.D0*(P1(k,j)-DX1(k,j))*p_0(j))
        enddo
      enddo
    enddo
    dPdmu = beta*p_0*(1.D0-p_0)
    mu1 = 0.D0
    sumdPdmu = 0.D0
    !do i = 1,HDIM
    do i = 1,norbsCore
      mu1 = mu1 + P1(i,i)
      sumdPdmu = sumdPdmu + dPdmu(i)
    enddo
    write(*,*)"mu1",mu1

    !mu1 = -mu1/SUM(dPdmu)
    mu1 = -mu1/sumdPdmu
    do i = 1,HDIM
      P1(i,i) = P1(i,i) + mu1*dPdmu(i)  ! Trace correction by adding (dP/dmu)*(dmu/dH1) to dP/dH1
    enddo

    !  call MMult(ONE,Q,P1,ZERO,X,'N','N',HDIM) ! and back transformation of P1,
    !  with a total of
    !  call MMult(ONE,X,Q,ZERO,P1,'N','T',HDIM) ! 4 matrix-matrix multiplication
    !  dominates the cost

  end subroutine Canon_DM_PRT

  subroutine Eig(A,Q,ee,TYPE,HDIM)

    implicit none
    integer,      parameter       :: PREC = 8
    integer,      intent(in)      :: HDIM
    integer    :: INFO
    integer    :: DIAG_LWORK
    real(PREC) :: DIAG_WORK(1+6*HDIM+2*HDIM*HDIM)
    integer    :: DIAG_LIWORK
    integer    :: DIAG_IWORK(3+5*HDIM)
    real(PREC), intent(in)      :: A(HDIM,HDIM)
    real(PREC), intent(out)     :: Q(HDIM,HDIM), ee(HDIM)
    real(PREC) :: EVECS(HDIM,HDIM), EVALS(HDIM)
    character(LEN=1), parameter :: JOBZ = "V",  UPLO = "U"
    character(1), intent(in)    :: TYPE  ! 'N'(for eigenvaules only) or 'V'(otherwise)

    DIAG_LWORK = 1+6*HDIM+2*HDIM*HDIM
    DIAG_LIWORK = 3 + 5*HDIM

    EVECS = A
    CALL DSYEVD(JOBZ, UPLO, HDIM, EVECS, HDIM, EVALS, &
         DIAG_WORK, DIAG_LWORK, DIAG_IWORK, DIAG_LIWORK, INFO)

    Q = EVECS
    ee = EVALS

  end subroutine Eig


end module prg_densitymatrix_mod
