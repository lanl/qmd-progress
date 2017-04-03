!> Module to obtain the density matrix by diagonalizing an prg_orthogonalized Hamiltonian.
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
  public :: prg_get_flevel, prg_build_density_T, prg_build_atomic_density, prg_build_density_T_Fermi

contains

  !> Builds the density matrix from \f$ H_0 \f$ for zero electronic temperature.
  !! \f$ \rho = C \Theta(\mu I - \epsilon) C^{\dagger} \f$
  !! Where,\f$ C \f$ is the matrix eigenvector and \f$ \epsilon \f$ is the matrix eigenvalue.
  !! \f$ \Theta() \f$  is the Heaviside function.
  !! \param ham_bml Input Orthogonalized Hamiltonian matrix.
  !! \param rho_bml Output density matrix,
  !! \param threshold Threshold for sparse matrix algebra.
  !! \param bndfil Filing factor.
  !! \warning This does not solve the generalized eigenvalue problem.
  !! The Hamiltonian that comes in has to be preprg_orthogonalized.
  !!
  subroutine prg_build_density_T0(ham_bml, rho_bml, threshold, bndfil)
    implicit none
    character(20)                      ::  bml_type
    integer                            ::  i, norb
    real(8), intent(in)                ::  bndfil, threshold
    real(dp)                           ::  nocc
    real(dp), allocatable              ::  eigenvalues(:)
    type(bml_matrix_t)                 ::  aux1_bml, aux_bml, eigenvectors_bml, occupation_bml
    type(bml_matrix_t), intent(in)     ::  ham_bml
    type(bml_matrix_t), intent(inout)  ::  rho_bml

    norb = bml_get_n(ham_bml)
    bml_type = bml_get_type(ham_bml)

    allocate(eigenvalues(nOrb))
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,eigenvectors_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,occupation_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,aux_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,aux1_bml)

    call bml_diagonalize(ham_bml,eigenvalues,eigenvectors_bml)

    nocc = norb*bndfil

    do i=1,norb    !Reusing eigenvalues to apply the theta function.
      if(i.LE.nocc) then
        eigenvalues(i) = 2.0_dp
      else
        eigenvalues(i) = 0.0_dp
      endif
    enddo
    if(abs(nocc - int(nocc)).gt.0)then
      eigenvalues(int(nocc)+1) = 1.0_dp
    endif

    call bml_set_diagonal(occupation_bml, eigenvalues) !eps(i,i) = eps(i)


    call bml_multiply(eigenvectors_bml, occupation_bml, aux_bml, 1.0_dp, 0.0_dp, threshold)
    call bml_transpose(eigenvectors_bml, aux1_bml)
    call bml_multiply(aux_bml, aux1_bml, rho_bml, 1.0_dp, 0.0_dp, threshold)

    deallocate(eigenvalues)
    call bml_deallocate(eigenvectors_bml)
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
  !! \warning This does not solve the generalized eigenvalue problem.
  !! The Hamiltonian that comes in has to be preprg_orthogonalized.
  !!
  subroutine prg_build_density_T(ham_bml, rho_bml, threshold, bndfil, kbt, ef)
    implicit none
    character(20)                      ::  bml_type
    integer                            ::  i, norb
    real(8), intent(in)                ::  bndfil, threshold, kbt
    real(8), intent(inout)             ::  ef
    real(dp)                           ::  nocc, fleveltol
    real(dp), allocatable              ::  eigenvalues(:)
    type(bml_matrix_t)                 ::  aux1_bml, aux_bml, eigenvectors_bml, occupation_bml
    type(bml_matrix_t), intent(in)     ::  ham_bml
    type(bml_matrix_t), intent(inout)  ::  rho_bml

    if (printRank() .eq. 1) then
      write(*,*)"In get_density_t ..."
    endif

    norb = bml_get_n(ham_bml)
    bml_type = bml_get_type(ham_bml)

    allocate(eigenvalues(nOrb))
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,eigenvectors_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,occupation_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,aux_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,aux1_bml)

    call bml_diagonalize(ham_bml,eigenvalues,eigenvectors_bml)

    fleveltol = 1.0e-5

    call prg_get_flevel(eigenvalues,kbt,bndfil,fleveltol,ef)

    nocc = norb*bndfil

    do i=1,norb   !Reusing eigenvalues to apply the theta function.
      eigenvalues(i) = 2.0_dp*fermi(eigenvalues(i),ef,kbt)
    enddo

    call bml_set_diagonal(occupation_bml, eigenvalues) !eps(i,i) = eps(i)

    call bml_multiply(eigenvectors_bml, occupation_bml, aux_bml, 1.0_dp, 0.0_dp, threshold)
    call bml_transpose(eigenvectors_bml, aux1_bml)
    call bml_multiply(aux_bml, aux1_bml, rho_bml, 1.0_dp, 0.0_dp, threshold)

    deallocate(eigenvalues)
    call bml_deallocate(eigenvectors_bml)
    call bml_deallocate(aux_bml)
    call bml_deallocate(aux1_bml)
    call bml_deallocate(occupation_bml)

  end subroutine prg_build_density_T


  !> Builds the density matrix from \f$ H_0 \f$ for electronic temperature T.
  !! \f$ \rho = C f(\mu I - \epsilon) C^{\dagger} \f$
  !! Where,\f$ C \f$ is the matrix eigenvector and \f$ \epsilon \f$ is the matrix eigenvalue.
  !! \f$ f \f$  is the Fermi function. In this routine the Fermi level is passed as an argument.
  !! \param ham_bml Input Orthogonalized Hamiltonian matrix.
  !! \param rho_bml Output density matrix,
  !! \param threshold Threshold for sparse matrix algebra.
  !! \warning This does not solve the generalized eigenvalue problem.
  !! The Hamiltonian that comes in has to be preprg_orthogonalized.
  !!
  subroutine prg_build_density_T_Fermi(ham_bml, rho_bml, threshold, kbt, ef,verbose)
    implicit none
    character(20)                      ::  bml_type
    integer                            ::  i, norb
    integer, optional, intent(in)      ::  verbose
    real(8), intent(in)                ::  threshold, kbt
    real(8), intent(in)                ::  ef
    real(dp)                           ::  nocc, fleveltol
    real(dp), allocatable              ::  eigenvalues(:)
    type(bml_matrix_t)                 ::  aux1_bml, aux_bml, eigenvectors_bml, occupation_bml
    type(bml_matrix_t), intent(in)     ::  ham_bml
    type(bml_matrix_t), intent(inout)  ::  rho_bml

    if (printRank() .eq. 1) then
      if(present(verbose).and.verbose >= 1)then
        write(*,*)"In get_density_t_Fermi ..."
        write(*,*)"Ef = ",ef
      endif
    endif

    norb = bml_get_n(ham_bml)
    bml_type = bml_get_type(ham_bml)

    allocate(eigenvalues(nOrb))
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,eigenvectors_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,occupation_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,aux_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,aux1_bml)

    call bml_diagonalize(ham_bml,eigenvalues,eigenvectors_bml)

    do i=1,norb   !Reusing eigenvalues to apply the theta function.
      eigenvalues(i) = 2.0_dp*fermi(eigenvalues(i),ef,kbt)
    enddo

    call bml_set_diagonal(occupation_bml, eigenvalues) !eps(i,i) = eps(i)

    call bml_multiply(eigenvectors_bml, occupation_bml, aux_bml, 1.0_dp, 0.0_dp, threshold)
    call bml_transpose(eigenvectors_bml, aux1_bml)
    call bml_multiply(aux_bml, aux1_bml, rho_bml, 1.0_dp, 0.0_dp, threshold)

    deallocate(eigenvalues)
    call bml_deallocate(eigenvectors_bml)
    call bml_deallocate(occupation_bml)
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
    implicit none
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
  !! \param kbt Temperature times the Boltzmans's constant  (\f$ k_bT  \f$).
  !! \param bndfil Filing factor (\f$ N_{el}/(2*N_{orbs})\f$).
  !! \param tol Tolerance for the bisection method.
  !! \param Ef Fermi level (\f$ \mu \f$).
  !!
  subroutine prg_get_flevel(eigenvalues,kbt,bndfil,tol,Ef)
    implicit none
    integer                  ::  i, j, k, m
    integer                  ::  norb
    real(dp)                 ::  Ft1, Ft2, Kb, Prod
    real(dp)                 ::  T, step, tol, nel
    real(dp), intent(in)     ::  bndfil, eigenvalues(:), kbt
    real(dp), intent(inout)  ::  Ef

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

    do m=1,1000001

      if(m.GT.1000000)then
        stop "Bisection method in prg_get_flevel not converging ..."
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


  !> Gets the eigenvalues of the Orthogonalized Hamiltonian.
  !! \param ham_bml Input Orthogonalized Hamiltonian matrix.
  !! \param eigenvalues Output eigenvalues of the system.
  !! \param verbose Verbosity level.
  !!
  subroutine prg_get_eigenvalues(ham_bml,eigenvalues,verbose)
    implicit none
    character(20)                        ::  bml_type
    integer                              ::  i, norb
    integer, intent(in)                  ::  verbose
    real(dp)                             ::  nocc
    real(dp), allocatable                ::  aux(:,:)
    real(dp), allocatable, intent(inout)  ::  eigenvalues(:)
    type(bml_matrix_t)                   ::  aux_bml, eigenvectors_bml
    type(bml_matrix_t), intent(in)       ::  ham_bml

    norb = bml_get_n(ham_bml)
    bml_type = bml_get_type(ham_bml)

    if(verbose.GE.1)write(*,*)"In prg_get_eigenvalues ..."
    if(verbose.GE.1)write(*,*)"Number of states =",norb

    if(.not.allocated(eigenvalues))allocate(eigenvalues(norb))

    !Ensure dense type to diagonalize
    if(trim(bml_type).ne."dense")then
      allocate(aux(norb,norb))
      call bml_convert_to_dense(ham_bml,aux)
      call bml_zero_matrix(bml_matrix_dense,bml_element_real,dp,norb,norb,aux_bml)
      call bml_convert_from_dense(bml_matrix_dense,aux,aux_bml,0.0_dp,norb)
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
    implicit none
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
    implicit none
    real(dp), intent(in) :: e, ef, kbt

    fermi = 1.0_dp/(1.0_dp+exp((e-ef)/(kbt)))

  end function fermi

end module prg_densitymatrix_mod
