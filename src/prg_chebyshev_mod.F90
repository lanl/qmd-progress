!> Module to obtain the density matrix by applying a Chebyshev
!! polynomial expansion.
!!
!! \ingroup PROGRESS
!!
!!  See Amparo Gil 2007 \cite Amparo2007 ,
!!  See Silver et al \cite Silver1996 ,
!!  See Weisse et al \cite Weisse2006
!!
!! \todo Add the power method in BML to get a better estimate of the spectral boundaries.
!!
module prg_Chebyshev_mod

  use bml
  use prg_normalize_mod
  use prg_densitymatrix_mod
  use prg_openfiles_mod

  implicit none

  private  !Everything is private by default

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.14159265358979323846264338327950_dp

  public :: prg_build_density_cheb

contains

  !> Builds the density matrix from \f$ H_0 \f$ for a Fermi function
  !! approximated with a Chebyshev polynomial expansion.
  !!
  !! \f$ \rho_{n+1} = b_{n+1}T_{n+1} + \rho_{n} \f$
  !! Where,\f$ T_n \f$ is the nth Chebyshev polynomial and
  !! \f$ b_{n} \f$  is the nth coefficient of the expansion for the Fermi function.
  !! \param ham_bml Input Orthogonalized Hamiltonian matrix.
  !! \param rho_bml Output density matrix,
  !! \param threshold Threshold for sparse matrix algebra.
  !! \param ncoeffs Number of Chebyshev coefficients.
  !! \param kbt Electronic temperature in the energy units of the Hamiltonian.
  !! \param ef Fermi level in the energy units of the Hamiltonian.
  !! \param verbose Verbosity level.
  !!
  subroutine prg_build_density_cheb(ham_bml, rho_bml, threshold, ncoeffs, kbt, ef, verbose)
  
    character(20)                      ::  bml_type
    integer                            ::  N, i, j, norb, io
    integer, intent(in)                ::  ncoeffs, verbose
    real(dp)                           ::  alpha, de, derf, error
    real(dp)                           ::  maxder, maxderf, maxdiff0, maxdiff1
    real(dp)                           ::  mycoeff, nocc, pi, scaledef
    real(dp)                           ::  scaledkbt, occ
    real(dp), allocatable              ::  coeffs(:), coeffs1(:), domain(:), domain0(:)
    real(dp), allocatable              ::  domain2(:), gbnd(:), tn(:), tnm1(:)
    real(dp), allocatable              ::  tnp1(:), trace(:)
    real(dp), intent(in)               ::  ef, kbt, threshold
    type(bml_matrix_t)                 ::  aux1_bml, aux_bml, eigenvectors_bml, occupation_bml
    type(bml_matrix_t)                 ::  tn_bml, tnm1_bml, tnp1_bml, x_bml
    type(bml_matrix_t), intent(in)     ::  ham_bml
    type(bml_matrix_t), intent(inout)  ::  rho_bml

    write(*,*)"Building rho via Chebyshev expansion of the Fermi function ..."

    norb = bml_get_n(ham_bml)
    bml_type = bml_get_type(ham_bml)

    call bml_copy_new(ham_bml,x_bml)

    call bml_gershgorin(x_bml, gbnd)

    call prg_normalize_cheb(x_bml,ef,alpha,scaledef)

    de = 0.01_dp !This energy step can be set smaller if needed
    N = floor((gbnd(2)-gbnd(1))/de)

    ! Function with Real domain to keep track of the expansion
    allocate(domain0(N))
    allocate(domain(N))
    allocate(domain2(N))

    ! Chebyshev polynomial for recursion applied to the tracking function
    allocate(tnp1(N))
    allocate(tnm1(N))
    allocate(tn(N))

    ! Chebyshev coefficients for the expansion
    allocate(coeffs(ncoeffs))
    allocate(coeffs1(ncoeffs))

    scaledef = ef  !For this version there is no scaled ef
    scaledkbt = kbt !For this version there is no scaled kbT

    call prg_get_chebcoeffs(scaledkbt,scaledef,ncoeffs,coeffs,gbnd(1),gbnd(2))

    ! Prepare bml matrices for recursion.
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,tnp1_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,tn_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,tnm1_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,aux1_bml)

    ! Set the domain for the tracking function
    do i=1,N
       domain0(i) = gbnd(1) + i*de
       domain0(i) = 2.0_dp*(domain0(i) - gbnd(1))/(gbnd(2)-gbnd(1)) - 1.0_dp
       tnp1(i) = 0.0_dp
       tnm1(i) = 1.0_dp
       tn(i) = domain0(i)
       domain(i) = 0.0_dp
    enddo

    !First step of recursion ...
    mycoeff = jackson(ncoeffs,1)*coeffs(1) !Application of the Jackson kernel
    call bml_add_identity(tnm1_bml, 1.0_dp, threshold) !T0
    call bml_scale(mycoeff,tnm1_bml,aux1_bml) !Rho(0) = coeffs(1)*T0
    domain = 0.0_dp + mycoeff*tnm1

    !Second step of recursion ...
    mycoeff = jackson(ncoeffs,2)*coeffs(2)
    call bml_copy(x_bml,tn_bml) !T1
    call bml_add_deprecated(1.0_dp,aux1_bml,mycoeff,tn_bml,threshold) !Rho(1) = Rho(0) + coeffs(2)*T(1)
    domain = domain + mycoeff*tn

    !Third to n-1 step of recursion ...
    if(verbose >= 1) write(*,*)"Chebyshev recursion ..."
    do i=2,ncoeffs-1

       mycoeff = coeffs(i+1)

       call bml_copy(tnm1_bml,tnp1_bml)
       tnp1 = tnm1
       call bml_multiply(x_bml,tn_bml,tnp1_bml,2.0_dp,-1.0_dp) !T(n+1) = 2xT(n) - T(n-1)
       tnp1 = 2.0_dp*domain0*tn - tnp1
       mycoeff = mycoeff*jackson(ncoeffs,i)

       call bml_add_deprecated(1.0_dp,aux1_bml,mycoeff,tnp1_bml,threshold) !Rho(n+1) = Rho(n) + b(n+1)*T(n+1)
       domain = domain + mycoeff*tnp1

       call bml_copy(tn_bml,tnm1_bml)
       call bml_copy(tnp1_bml,tn_bml)
       tnm1 = tn
       tn = tnp1

       occ = 2.0_dp*bml_trace(aux1_bml)
       if(verbose >= 1) write(*,*)"step, occ",i,occ

    enddo

    if(verbose >= 2) then
       call prg_open_file(io,"fermi_approx.out")
       write(io,*)"# Energy, FApprox, Fermi"
       do i=1,N
          write(io,*)gbnd(1) + i*de,domain(i),fermi(gbnd(1) + i*de, scaledef, scaledkbt)
       enddo
    endif

    if(verbose >= 2) then
      maxder = absmaxderivative(domain,de)
      write(*,*)"TargetKbt =",scaledkbt
      write(*,*)"AbsMaxDerivative =",maxder
      write(*,*)"kbT = 1/(4*AbsMaxDerivative) =",1.0_dp/(4.0_dp*maxder)
    endif

    call bml_copy(aux1_bml,rho_bml)
    call bml_scale(2.0d0, rho_bml)

    write(*,*)"TotalOccupation =", bml_trace(rho_bml)

  end subroutine prg_build_density_cheb

  !> Evaluates the Jackson Kernel Coefficient.
  !! \param ncoeffs Number of Chebyshev polynomial.
  !! \param i Coefficient number i.
  !!
  real(dp) function jackson(ncoeffs,i)

    integer, intent(in) :: ncoeffs,i

    if(i == 1)then
       jackson = 1.0_dp
    endif

    if(i == 2)then
       jackson = real(ncoeffs)*cos(pi/(real(ncoeffs+1))) &
            +   sin(pi/(real(ncoeffs+1)))*1.0_dp/tan(pi/(real(ncoeffs+1)))
       jackson = jackson/real(ncoeffs+1)
    endif

    if(i > 2)then
       jackson = real(ncoeffs - i + 1)*cos(pi*i/(real(ncoeffs+1))) &
            +   sin(pi*i/(real(ncoeffs+1)))*1.0_dp/tan(pi/(real(ncoeffs+1)))
       jackson = jackson/real(ncoeffs + 1)
    endif

    !       jackson = 1.0_dp
  end function jackson

  !> Gets the coefficients of the Chebyshev expansion.
  !! \param kbt Electronic temperature
  !! \param ef Fermi level
  !! \param ncoeffs Number of Chebyshev coefficients
  !! \param coeffs Output vector for the Chebyshev coefficients
  !! \param emin lowest boundary for the eigenvalues of H
  !! \param emax highest boundary for the eigenvalues of H
  !!
  !! \note coeffs(1) = C0 !
  !!
  subroutine prg_get_chebcoeffs(kbt,ef,ncoeffs,coeffs,emin,emax)

    integer                  ::  i, j, npts, r
    integer, intent(in)      ::  ncoeffs
    real(dp)                 ::  Int, Kr, Kr0, fapp
    real(dp)                 ::  x, xj
    real(dp), intent(in)     ::  emax, emin, kbt, ef
    real(dp), intent(inout)  ::  coeffs(:)

    npts= 500 !Length of the discretization

    !Get coefficient of nth cheb expansion
    Kr = 0.5_dp*real(npts+1.0d0)
    Kr0 = real(npts+1.0d0)

    coeffs = 0.0d0

    do r = 0,ncoeffs-1

       Int = 0.0d0
       do j=0,npts
          xj = cos((j+0.5_dp)*pi/(npts + 1))
          x = (emax-emin)*(xj + 1.0d0)/2.0d0 + emin
          Int = Int + Tr(r,xj)*fermi(x,ef,kbt)
       enddo

       if(r == 0) then
          coeffs(r+1) = Int/Kr0
       else
          coeffs(r+1) = Int/Kr
       endif

    enddo

  end subroutine prg_get_chebcoeffs

  !> Chebyshev polynomial obtained by recursion.
  !! \param r rth polynomial
  !! \param x argument the evaluate the polynomial.
  !!
  real(dp) function Tr(r,x)

    real(dp), intent(in) :: x
    integer, intent(in)  :: r

    Tr = cos(r*acos(x))

  end function Tr

  !> Gives the Fermi distribution value for energy e.
  !! \param e Energy.
  !! \param ef Fermi energy.
  !!
  real(dp) function fermi(e,ef,kbt)

    real(dp), intent(in) :: e, ef, kbt

    fermi = 1.0_dp/(1.0_dp+exp((e-ef)/(kbt)))

  end function fermi

  !> Gets the absolute maximum of the derivative of a function
  !! \param func
  !! \param de Energy step
  !!
  real(dp) function absmaxderivative(func,de)

    real(dp), intent(in) :: func(:), de
    integer :: j

    absmaxderivative = -10000.0d0

    do j=1,size(func, dim=1)-1
       if(abs(func(j+1) - func(j))/de > absmaxderivative) &
       absmaxderivative = abs(func(j+1) - func(j))/de
    enddo

  end function absmaxderivative

end module prg_Chebyshev_mod
