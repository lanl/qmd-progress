!> The prg_normalize module.
    !
    !

module prg_normalize_mod

  use bml
  use prg_graph_mod
  use prg_parallel_mod

  implicit none

  private  !Everything is private by default

  integer, parameter :: dp = kind(1.0d0)

  public :: prg_normalize
  public :: prg_normalize_fermi
  public :: prg_normalize_cheb
  public :: prg_gershgorinReduction

contains

  !> Normalize a Hamiltonian matrix prior to running the SP2 algorithm.
  !!
  !! X0 = (e_max * I - H) / (e_max - e_min)
  !!
  !! where e_max and e_min are obtained sing the Gershgorin circle theorem.
  !!
  !! \param h_bml Input/Output Hamiltonian matrix
  subroutine prg_normalize(h_bml)

    type(bml_matrix_t),intent(inout) :: h_bml

    real(dp) :: alpha, beta, maxMinusMin
    real(dp), allocatable :: gbnd(:)

    allocate(gbnd(2))

    call bml_gershgorin(h_bml, gbnd)
    maxMinusMin = gbnd(2) - gbnd(1)
    alpha = -1.00_dp / maxMinusMin
    beta = gbnd(2) /maxMinusMin
    call bml_scale_add_identity(h_bml, alpha, beta, 0.00_dp)

    deallocate(gbnd)

  end subroutine prg_normalize

  !> Normalize a Hamiltonian matrix prior to running the truncated SP2 algorithm.
  !!
  !! X0 = ((hN-mu) * I - H) / (hN - h1)
! X0 = (hN*I-H0-mu*I)/(hN-h1)
  !!
  !! where h1 and hN are scaled Gershgorin bounds.
  !!
  !! \param H_bml Hamiltonian matrix
  !! \param h1    Scaled minimum Gershgorin bound.
  !! \param hN    Scaled maximum Gershgorin bound.
  !! \param mu    Chemical potential
  subroutine prg_normalize_fermi(h_bml, h1, hN, mu)

    type(bml_matrix_t),intent(inout) :: h_bml
    real(dp), intent(in) :: h1, hN, mu

    real(dp) :: alpha, beta, maxMinusMinInverse

    maxMinusMinInverse = 1.0_dp/(hN - h1)
    alpha = -1.0_dp
    beta = hN - mu
    call bml_scale_add_identity(h_bml, alpha, beta, 0.00_dp)
    call bml_scale(maxMinusMinInverse, h_bml)

  end subroutine prg_normalize_fermi

  !> Determine gershgorin bounds across all parts, local and distributed.
  subroutine prg_gershgorinReduction(gp)

    implicit none

    type(graph_partitioning_t), intent(inout) :: gp

    integer :: i, nRanks
    real(dp) :: localVal

    nRanks = getNRanks()
    if (nRanks .eq. 1) return

    localVal = gp%mineval
    call prg_minRealReduce(localVal)
    gp%mineval = localVal

    localVal = gp%maxeval
    call prg_maxRealReduce(localVal)
    gp%maxeval = localVal

  end subroutine prg_gershgorinReduction

  !> Normalize a Hamiltonian matrix prior to running the Chebyshev algorithm.
  !!
  !! X0 = 2*(H - e_min*I) / (e_max - e_min) - I
  !!
  !! where e_max and e_min are obtained sing the Gershgorin circle theorem.
  !!
  !! \param h_bml Input/Output Hamiltonian matrix
  subroutine prg_normalize_cheb(h_bml,mu,alpha,scaledmu)

    type(bml_matrix_t),intent(inout) :: h_bml

    real(dp) :: beta, maxMinusMin
    real(dp), intent(in) :: mu
    real(dp), intent(inout) :: scaledmu, alpha
    real(dp), allocatable :: gbnd(:)

    allocate(gbnd(2))

    call bml_gershgorin(h_bml, gbnd)
    maxMinusMin = gbnd(2) - gbnd(1)
    alpha = 2.00_dp / maxMinusMin
    beta = -1.00_dp*(2.00_dp*gbnd(1)/maxMinusMin + 1.00_dp)
    call bml_scale_add_identity(h_bml, alpha, beta, 0.00_dp)
    scaledmu = alpha*mu + beta

    deallocate(gbnd)

  end subroutine prg_normalize_cheb

end module prg_normalize_mod
