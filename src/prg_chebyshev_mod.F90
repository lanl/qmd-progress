!> Module to obtain the density matrix by applying a Chebyshev
!! polynomial expansion.
!!
!! \ingroup PROGRESS
!!
!!  See Amparo Gil 2007 \cite Amparo2007 ,
!!  See Silver et al \cite Silver1996 ,
!!  See Weisse et al \cite Weisse2006
!!
module prg_Chebyshev_mod

  use bml
  use prg_normalize_mod
  use prg_densitymatrix_mod
  use prg_openfiles_mod
  use prg_extras_mod
  use prg_kernelparser_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 3.14159265358979323846264338327950_dp

  !> General Cheb solver type
  !!
  type, public :: Chebdata_type
     character(100)   ::  flavor
     character(100)   ::  bml_type, jobname
     integer          ::  mdim, ncoeffs, ndim, verbose
     integer          ::  npts
     real(dp)         ::  atr, bndfil, ef, estep
     real(dp)         ::  fermitol, kbt, threshold
     logical          ::  getef, jon, trkfunc
  end type Chebdata_type

  public :: prg_build_density_cheb, prg_build_density_cheb_fermi
  public :: prg_parse_cheb

contains

  !> Chebyshev parser.
  !! This module is used to parse all the input variables for the cheb
  !! electronic structure solver.
  !! Adding a new input keyword to the parser:
  !! - If the variable is real, we have to increase nkey_re.
  !! - Add the keyword (character type) in the keyvector_re vector.
  !! - Add a default value (real type) in the valvector_re.
  !! - Define a new variable and pass the value through valvector_re(num)
  !! where num is the position of the new keyword in the vector.
  !!
  subroutine prg_parse_cheb(chebdata,filename)

    implicit none
    type(chebdata_type), intent(inout) :: chebdata
    integer, parameter :: nkey_char = 3, nkey_int = 5, nkey_re = 7, nkey_log = 3
    character(len=*) :: filename

    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=50) :: &
         'JobName=', 'BMLType=', 'Flavor=' ]
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
         'MyJob'   , 'Dense'   , 'Alg1' ]

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
         'MDim=', 'NDim=', 'NCoeffs=','Verbose=','NPoints=']
    integer :: valvector_int(nkey_int) = (/ &
         -1   ,    0    ,     50      ,  0, 500 /)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
         'NumThresh=','FermiTol=','BndFil=','EStep=',"ATr=","Kbt=","Ef=" ]
    real(dp) :: valvector_re(nkey_re) = (/&
         0.0    ,   0.00000001   ,0.0, 0.01, 0.0, 0.0, 0.0 /)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=50) :: &
         'GetEf=', 'Jackson=','TRKFunction=']
    logical :: valvector_log(nkey_log) = (/&
         .false., .false., .false./)

    !Start and stop characters
    character(len=50), parameter :: startstop(2) = [character(len=50) :: &
         'CHEB{', '}']

    call prg_parsing_kernel(keyvector_char,valvector_char&
         ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
         keyvector_log,valvector_log,trim(filename),startstop)

    !Characters
    chebdata%JobName = valvector_char(1)

    if(valvector_char(2) == "Dense")then
       chebdata%bml_type = BML_MATRIX_DENSE
    elseif(valvector_char(2) == "Ellpack")then
       chebdata%bml_type = BML_MATRIX_ELLPACK
    endif
    chebdata%flavor = valvector_char(3)

    !Reals
    chebdata%threshold = valvector_re(1)
    chebdata%fermitol = valvector_re(2)
    chebdata%bndfil = valvector_re(3)
    chebdata%estep = valvector_re(4)
    chebdata%atr = valvector_re(5)
    chebdata%kbt = valvector_re(6)
    chebdata%ef = valvector_re(7)

    !Logicals
    chebdata%getef = valvector_log(1)
    chebdata%jon = valvector_log(2)
    chebdata%trkfunc = valvector_log(3)

    !Integers
    chebdata%mdim = valvector_int(1)
    chebdata%ndim = valvector_int(2)
    chebdata%ncoeffs = valvector_int(3)
    chebdata%verbose = valvector_int(4)
    chebdata%npts = valvector_int(5)

  end subroutine prg_parse_cheb

  !> Builds the density matrix from \f$ H_0 \f$ for a Fermi function
  !! approximated with a Chebyshev polynomial expansion.
  !!
  !! \f$ \rho_{n+1} = b_{n+1}T_{n+1} + \rho_{n} \f$
  !! Where,\f$ T_n \f$ is the nth Chebyshev polynomial and
  !! \f$ b_{n} \f$  is the nth coefficient of the expansion for the Fermi function.
  !! In the sparse version (when ellpack is used) the threshold can be varied
  !! linearly with the polynomial degree. The function is the following:
  !! \f$ Thresh(n) = Thresh_0 [a_{thr} (n-1) + (1-a_{thr})]  \f$
  !! \param ham_bml Input Orthogonalized Hamiltonian matrix.
  !! \param rho_bml Output density matrix.
  !! \param athr Threshold linear increasing constant.
  !! \param threshold Threshold for sparse matrix algebra.
  !! \param ncoeffs Number of Chebyshev coefficients.
  !! \param kbt Electronic temperature in the energy units of the Hamiltonian.
  !! \param ef Fermi level in the energy units of the Hamiltonian.
  !! \param bndfil Band filing factor.
  !! \param verbose Verbosity level.
  !!
  subroutine prg_build_density_cheb(ham_bml, rho_bml, athr, threshold, ncoeffs, &
       kbt, ef, bndfil, jon, verbose)
    character(20)                      ::  bml_type
    integer                            ::  npts, enpts, i, io
    integer                            ::  norb, mdim
    integer, intent(in)                ::  ncoeffs, verbose
    real(dp)                           ::  alpha, de, maxder, mls_I
    real(dp)                           ::  mycoeff, occ, scaledef, scaledkbt
    real(dp)                           ::  threshold1
    real(dp), allocatable              ::  coeffs(:), coeffs1(:), domain(:), domain0(:)
    real(dp), allocatable              ::  domain2(:), gbnd(:), tn(:), tnm1(:)
    real(dp), allocatable              ::  tnp1(:), tnp1_dense(:,:), tracesT(:)
    real(dp), intent(in)               ::  athr, bndfil, kbt, threshold
    real(dp), intent(in)               ::  ef
    type(bml_matrix_t)                 ::  aux_bml, tn_bml, tnm1_bml, tnp1_bml
    type(bml_matrix_t)                 ::  x_bml
    type(bml_matrix_t), intent(in)     ::  ham_bml
    type(bml_matrix_t), intent(inout)  ::  rho_bml
    logical, intent(in)                ::  jon

    if(verbose >= 1)write(*,*)"Building rho via Chebyshev expansion of the Fermi function ..."

    norb = bml_get_n(ham_bml) !Get the number of orbitals from H
    bml_type = bml_get_type(ham_bml) !Get the bml type

    if(verbose >= 1)mls_I = mls()
    call bml_copy_new(ham_bml,x_bml)
    call bml_gershgorin(x_bml, gbnd)
    call prg_normalize_cheb(x_bml,ef,gbnd(1),gbnd(2),alpha,scaledef)
    if(verbose >= 1)write(*,*)"Time for gershgorin and normalize",mls()-mls_I
    de = 0.01_dp !This energy step can be set smaller if needed
    enpts = floor((gbnd(2)-gbnd(1))/de)

    ! Defining a function with "Real domain" to keep track of the expansion
    allocate(domain0(enpts))
    allocate(domain(enpts))
    allocate(domain2(enpts))

    ! Chebyshev polynomial for recursion applied to the tracking function
    allocate(tnp1(enpts))
    allocate(tnm1(enpts))
    allocate(tn(enpts))
    allocate(tnp1_dense(norb,norb))

    ! Chebyshev coefficients for the expansion
    allocate(coeffs(ncoeffs))
    allocate(coeffs1(ncoeffs))
    allocate(tracesT(ncoeffs))
    tracesT = 0.0_dp

    ! First computation of the Chebyshev coefficients (Non-Ef)
    if(verbose >= 1)mls_I = mls()
    call prg_get_chebcoeffs(npts,kbt,ef,ncoeffs,coeffs,gbnd(1),gbnd(2))
    if(verbose >= 1)write(*,*)"Time for prg_get_chebcoeffs",mls()-mls_I

    ! Prepare bml matrices for recursion.
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,tnp1_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,tn_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,tnm1_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,norb,aux_bml)

    ! Set the domain for the tracking function
    if(verbose >= 1)mls_I = mls()
    do i=1,enpts
       domain0(i) = gbnd(1) + i*de
       domain0(i) = 2.0_dp*(domain0(i) - gbnd(1))/(gbnd(2)-gbnd(1)) - 1.0_dp
       tnp1(i) = 0.0_dp
       tnm1(i) = 1.0_dp
       tn(i) = domain0(i)
       domain(i) = 0.0_dp
    enddo
    if(verbose >= 1)write(*,*)"Time for setting the tracking function",mls()-mls_I

    !First step of recursion ...
    if(verbose >= 1)mls_I = mls()
    mycoeff = jackson(ncoeffs,1,jon)*coeffs(1) !Application of the Jackson kernel
    call bml_add_identity(tnm1_bml, 1.0_dp, threshold) !T0
    call bml_scale(mycoeff,tnm1_bml,aux_bml) !Rho(0) = coeffs(1)*T0
    !Second step of recursion ...
    call bml_copy(x_bml,tn_bml) !T1
    call bml_add_deprecated(1.0_dp,aux_bml,mycoeff,tn_bml,threshold) !Rho(1) = Rho(0) + coeffs(2)*T(1)

    tracesT(1) = bml_trace(tnm1_bml)
    domain = 0.0_dp + mycoeff*tnm1
    mycoeff = jackson(ncoeffs,2,jon)*coeffs(2)
    tracesT(2) = bml_trace(tn_bml)
    domain = domain + mycoeff*tn

    !Third to n-1 step of recursion ...
    if(verbose >= 1) write(*,*)"Chebyshev recursion ..."
    do i=2,ncoeffs-1

       mycoeff = coeffs(i+1)*jackson(ncoeffs,i+1,jon)

       call bml_copy(tnm1_bml,tnp1_bml)
       tnp1 = tnm1

       threshold1 =  threshold*(athr*real(i-1) + (1.0_dp-athr))

       call bml_multiply(x_bml,tn_bml,tnp1_bml,2.0_dp,-1.0_dp,threshold1) !T(n+1) = 2xT(n) - T(n-1)
       tracesT(i+1) = bml_trace(tnp1_bml)

       if(verbose >= 3)then
          write(*,*)"Time for mult",mls()-mls_I
          write(*,*)"Coeff",abs(mycoeff)
          write(*,*)"Bandwidth of Tn, Threshold",bml_get_bandwidth(tnp1_bml),threshold1
          write(*,*)"Sparsity of Tn",bml_get_sparsity(tnp1_bml,threshold)
       endif

       tnp1 = 2.0_dp*domain0*tn - tnp1

       call bml_add_deprecated(1.0_dp,aux_bml,mycoeff,tnp1_bml,threshold) !Rho(n+1) = Rho(n) + b(n+1)*T(n+1)
       domain = domain + mycoeff*tnp1

       call bml_copy(tn_bml,tnm1_bml)
       call bml_copy(tnp1_bml,tn_bml)
       tnm1 = tn
       tn = tnp1

    enddo
    if(verbose >= 1)write(*,*)"Time for recursion",mls()-mls_I

    if(verbose >= 2) then
       call prg_open_file(io,"fermi_approx.out")
       write(io,*)"# Energy, FApprox, Fermi"
       do i=1,enpts
          write(io,*)gbnd(1) + i*de,domain(i),fermi(gbnd(1) + i*de, ef, kbt)
       enddo
    endif

    if(verbose >= 2) then
       maxder = absmaxderivative(domain,de)
       write(*,*)"TargetKbt =",scaledkbt
       write(*,*)"AbsMaxDerivative =",maxder
       write(*,*)"kbT = 1/(4*AbsMaxDerivative) =",1.0_dp/(4.0_dp*maxder)
    endif

    call bml_copy(aux_bml,rho_bml)
    call bml_scale(2.0d0, rho_bml)

    if(verbose >= 1)write(*,*)"TotalOccupation =", bml_trace(rho_bml)

  end subroutine prg_build_density_cheb

  !> Builds the density matrix from \f$ H_0 \f$ for a Fermi function
  !! approximated with a Chebyshev polynomial expansion.
  !! In this case the self-consistent recursion is applied to converge
  !! to the correct number of electrons and obtain the Fermi level.
  !!
  !! \f$ \rho_{n+1} = b_{n+1}T_{n+1} + \rho_{n} \f$
  !! Where,\f$ T_n \f$ is the nth Chebyshev polynomial and
  !! \f$ b_{n} \f$  is the nth coefficient of the expansion for the Fermi function.
  !! In the sparse version (when ellpack is used) the threshold can be varied
  !! linearly with the polynomial degree. The function is the following:
  !! \f$ Thresh(n) = Thresh_0 [a_{thr} (n-1) + (1-a_{thr})]  \f$
  !! \param ham_bml Input Orthogonalized Hamiltonian matrix.
  !! \param rho_bml Output density matrix.
  !! \param athr Threshold linear increasing constant.
  !! \param threshold Threshold for sparse matrix algebra.
  !! \param ncoeffs Number of Chebyshev coefficients.
  !! \param kbt Electronic temperature in the energy units of the Hamiltonian.
  !! \param ef Fermi level in the energy units of the Hamiltonian.
  !! \param bndfil Band filing factor.
  !! \param npts Number of energy points to compute the coefficients
  !! \param verbose Verbosity level.
  !!
  subroutine prg_build_density_cheb_fermi(ham_bml, rho_bml, athr, threshold, ncoeffs, &
       kbt, ef, bndfil, getef, fermitol, jon, npts, trkfunc, verbose)

    character(20)                      ::  bml_type
    integer                            ::  npts, enpts, i, io
    integer                            ::  norb, mdim
    integer, intent(in)                ::  ncoeffs, verbose
    real(dp)                           ::  alpha, de, maxder, mls_I, mls_R
    real(dp)                           ::  mycoeff, occ, scaledef
    real(dp)                           ::  threshold1, fermitol
    real(dp), allocatable              ::  coeffs(:), coeffs1(:), domain(:), domain0(:)
    real(dp), allocatable              ::  domain2(:), gbnd(:), tn(:), tnm1(:)
    real(dp), allocatable              ::  tnp1(:), tnp1_dense(:,:), tracesT(:)
    real(dp), intent(in)               ::  athr, bndfil, kbt, threshold
    real(dp), intent(out)              ::  ef
    type(bml_matrix_t)                 ::  aux_bml, tn_bml, tnm1_bml, tnp1_bml
    type(bml_matrix_t)                 ::  x_bml
    type(bml_matrix_t), intent(in)     ::  ham_bml
    type(bml_matrix_t), intent(inout)  ::  rho_bml
    logical, intent(in)                ::  getef, jon, trkfunc

    if(verbose >= 1)write(*,*)"Building rho via Chebyshev expansion of the Fermi function ..."
    ef = 0.0_dp
    write(*,*)"ef",ef
    norb = bml_get_n(ham_bml) !Get the number of orbitals from H
    mdim = bml_get_m(ham_bml)
    bml_type = bml_get_type(ham_bml) !Get the bml type
    mdim = bml_get_m(ham_bml)
    if(verbose >= 1)mls_I = mls()
    call bml_copy_new(ham_bml,x_bml)
    allocate(gbnd(2))
    call bml_gershgorin(x_bml, gbnd)
    if(verbose >= 1)write(*,*)"Estimation for emin, emax", gbnd(1), gbnd(2)
    call prg_normalize_cheb(x_bml,ef,gbnd(1),gbnd(2),alpha,scaledef)
    if(verbose >= 1)write(*,*)"Time for gershgorin and normalize",mls()-mls_I

    if(trkfunc)then
       de = 0.001_dp !This energy step can be set smaller if needed
       enpts = floor((gbnd(2)-gbnd(1))/de)

       ! Defining a function with "Real domain" to keep track of the expansion
       allocate(domain0(enpts))
       allocate(domain(enpts))
       allocate(domain2(enpts))

       ! Chebyshev polynomial for recursion applied to the tracking function
       allocate(tnp1(enpts))
       allocate(tnm1(enpts))
       allocate(tn(enpts))
       allocate(tnp1_dense(norb,norb))
    endif

    ! Chebyshev coefficients for the expansion
    allocate(coeffs(ncoeffs))
    allocate(coeffs1(ncoeffs))
    allocate(tracesT(ncoeffs))
    tracesT = 0.0_dp

    ! First computation of the Chebyshev coefficients (Non-Ef)
    if(verbose >= 1)mls_I = mls()
    call prg_get_chebcoeffs(npts,kbt,ef,ncoeffs,coeffs,gbnd(1),gbnd(2))
    if(verbose >= 1)write(*,*)"Time for prg_get_chebcoeffs",mls()-mls_I

    ! Prepare bml matrices for recursion.
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,mdim,tnp1_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,mdim,tn_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,mdim,tnm1_bml)
    call bml_zero_matrix(bml_type,bml_element_real,dp,norb,mdim,aux_bml)

    ! Set the domain for the tracking function
    if(trkfunc)then
       mls_I = mls()
       do i=1,enpts
          domain0(i) = gbnd(1) + i*de
          domain0(i) = 2.0_dp*(domain0(i) - gbnd(1))/(gbnd(2)-gbnd(1)) - 1.0_dp
          tnp1(i) = 0.0_dp
          tnm1(i) = 1.0_dp
          tn(i) = domain0(i)
          domain(i) = 0.0_dp
       enddo
       write(*,*)"Time for setting the tracking function",mls()-mls_I
    endif

    if(getef)then
       if(verbose >= 1)write(*,*)"Computing Ef ..."
       !Compute Ts to get the traces
       if(verbose >= 1)mls_I = mls()
       call bml_add_identity(tnm1_bml, 1.0_dp, threshold) !T0
       tracesT(1) = norb
       call bml_copy(x_bml,tn_bml) !T1
       tracesT(2) = bml_trace(tn_bml)
       do i=2,ncoeffs-1
          call bml_copy(tnm1_bml,tnp1_bml)
          threshold1 =  threshold*(athr*real(i-1) + (1.0_dp-athr))
          call bml_multiply(x_bml,tn_bml,tnp1_bml,2.0_dp,-1.0_dp,threshold1) !T(n+1) = 2xT(n) - T(n-1)
          tracesT(i+1) = bml_trace(tnp1_bml)
          call bml_copy(tn_bml,tnm1_bml)
          call bml_copy(tnp1_bml,tn_bml)
       enddo
       if(verbose >= 1)write(*,*)"Time for computing traces",mls()-mls_I

       ! Clear bml matrices
       call bml_deallocate(tnp1_bml)
       call bml_deallocate(tn_bml)
       call bml_deallocate(tnm1_bml)
       call bml_zero_matrix(bml_type,bml_element_real,dp,norb,mdim,tnp1_bml)
       call bml_zero_matrix(bml_type,bml_element_real,dp,norb,mdim,tn_bml)
       call bml_zero_matrix(bml_type,bml_element_real,dp,norb,mdim,tnm1_bml)

       mls_I = mls()
       call prg_get_chebcoeffs_fermi_nt(npts,kbt,ef,tracesT,ncoeffs,coeffs,gbnd(1),&
            gbnd(2),bndfil,norb,fermitol,jon,verbose)
       if(verbose >= 1)write(*,*)"Time for prg_get_chebcoeffs_fermi_nt",mls()-mls_I
       if(verbose >= 1)write(*,*)"Converged Ef",ef
       call prg_get_chebcoeffs(npts,kbt,ef,ncoeffs,coeffs,gbnd(1),gbnd(2))
    endif

    !First step of recursion ...
    if(verbose >= 1)mls_R = mls()
    mycoeff = jackson(ncoeffs,1,jon)*coeffs(1) !Application of the Jackson kernel
    call bml_add_identity(tnm1_bml, 1.0_dp, threshold) !T0
    call bml_scale(mycoeff,tnm1_bml,aux_bml) !Rho(0) = coeffs(1)*T0
    if(trkfunc)then
       tracesT(1) = bml_trace(tnm1_bml)
       domain = 0.0_dp + mycoeff*tnm1
    endif

    !Second step of recursion ...
    mycoeff = jackson(ncoeffs,2,jon)*coeffs(2)
    call bml_copy(x_bml,tn_bml) !T1
    call bml_add_deprecated(1.0_dp,aux_bml,mycoeff,tn_bml,threshold) !Rho(1) = Rho(0) + coeffs(2)*T(1)
    if(trkfunc)then
       tracesT(2) = bml_trace(tn_bml)
       domain = domain + mycoeff*tn
    endif

    !Third to n-1 step of recursion ...
    if(verbose >= 1) write(*,*)"Chebyshev recursion ..."
    do i=2,ncoeffs-1

       mycoeff = coeffs(i+1)*jackson(ncoeffs,i+1,jon)

       call bml_copy(tnm1_bml,tnp1_bml)
       if(trkfunc)tnp1 = tnm1

       if(verbose >= 4)mls_I = mls()
       threshold1 =  threshold*(athr*real(i-1) + (1.0_dp-athr))

       call bml_multiply(x_bml,tn_bml,tnp1_bml,2.0_dp,-1.0_dp,threshold1) !T(n+1) = 2xT(n) - T(n-1)
       tracesT(i+1) = bml_trace(tnp1_bml)

       if(verbose >= 4)then
          write(*,*)"Time for mult",mls()-mls_I
          write(*,*)"Coeff",i,abs(mycoeff)
          write(*,*)"Bandwidth of Tn, Threshold",bml_get_bandwidth(tnp1_bml),threshold1
          write(*,*)"Sparsity of Tn",bml_get_sparsity(tnp1_bml,threshold)
       endif

       if(trkfunc)tnp1 = 2.0_dp*domain0*tn - tnp1

       call bml_add_deprecated(1.0_dp,aux_bml,mycoeff,tnp1_bml,threshold) !Rho(n+1) = Rho(n) + b(n+1)*T(n+1)
       if(trkfunc)domain = domain + mycoeff*tnp1

       call bml_copy(tn_bml,tnm1_bml)
       call bml_copy(tnp1_bml,tn_bml)
       if(trkfunc)then
          tnm1 = tn
          tn = tnp1
       endif

       if(verbose >= 4)then
          occ = 2.0_dp*bml_trace(aux_bml)
          if(verbose >= 1) write(*,*)"Step, Occupation",i,occ
          write(*,*)"Time for", i,"recursion",mls()-mls_I
       endif
    enddo
    if(verbose >= 1)write(*,*)"Time for recursion",mls()-mls_R

    if(trkfunc) then
       call prg_open_file(io,"fermi_approx.out")
       write(io,*)"# Energy, FApprox, Fermi"
       do i=1,enpts
          write(io,*)gbnd(1) + i*de,domain(i),fermi(gbnd(1) + i*de, ef, kbt)
       enddo
       close(io)
    endif

    if(verbose >= 2) then
       maxder = absmaxderivative(domain,de)
       write(*,*)"TargetKbt =",kbt
       write(*,*)"AbsMaxDerivative =",maxder
       write(*,*)"Actual Kbt = 1/(4*AbsMaxDerivative) =",1.0_dp/(4.0_dp*maxder)
    endif

    call bml_copy(aux_bml,rho_bml)
    call bml_scale(2.0d0, rho_bml)

    if(verbose >= 1)write(*,*)"TotalOccupation =", bml_trace(rho_bml)

    if(trkfunc)then
       deallocate(domain0)
       deallocate(domain)
       deallocate(domain2)
       deallocate(tnp1)
       deallocate(tnm1)
       deallocate(tn)
       deallocate(tnp1_dense)
    endif

    deallocate(coeffs)
    deallocate(coeffs1)
    deallocate(tracesT)
    call bml_deallocate(tnp1_bml)
    call bml_deallocate(tn_bml)
    call bml_deallocate(tnm1_bml)
    call bml_deallocate(aux_bml)
    call bml_deallocate(x_bml)

  end subroutine prg_build_density_cheb_fermi

  !> Evaluates the Jackson Kernel Coefficients.
  !! \param ncoeffs Number of Chebyshev polynomial.
  !! \param i Coefficient number i.
  !!
  real(dp) function jackson(ncoeffs,i,jon)

    integer, intent(in) :: ncoeffs,i
    logical, intent(in) :: jon

    if(jon)then
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
    else
       jackson = 1.0_dp
    endif

  end function jackson

  !> Gets the coefficients of the Chebyshev expansion.
  !! \param npts Number of points for discretization.
  !! \param kbt Electronic temperature.
  !! \param ef Fermi level.
  !! \param ncoeffs Number of Chebyshev coefficients.
  !! \param coeffs Output vector for the Chebyshev coefficients.
  !! \param emin lowest boundary for the eigenvalues of H.
  !! \param emax highest boundary for the eigenvalues of H.
  !!
  subroutine prg_get_chebcoeffs(npts,kbt,ef,ncoeffs,coeffs,emin,emax)

    integer                  ::  i,j
    integer, intent(in)      ::  ncoeffs, npts
    real(dp)                 ::  Int, Kr, Kr0, x
    real(dp)                 ::  xj
    real(dp), intent(in)     ::  ef, emax, emin, kbt
    real(dp), intent(inout)  ::  coeffs(:)

    !Get coefficient of nth cheb expansion.
    Kr = 0.5_dp*real(npts+1.0d0)
    Kr0 = real(npts+1.0d0)

    coeffs = 0.0d0

    !$omp parallel do default(none) private(i) &
    !$omp private(j,x,xj,Int) &
    !$omp shared(emin,emax,npts,ef,kbt,Kr,Kr0,coeffs,ncoeffs)
    do i = 0,ncoeffs-1

       Int = 0.0d0
       do j=0,npts
          xj = cos((j+0.5_dp)*pi/(npts + 1))
          x = (emax-emin)*(xj + 1.0d0)/2.0d0 + emin
          Int = Int + Tr(i,xj)*fermi(x,ef,kbt)
       enddo

       if(i == 0) then
          coeffs(i+1) = Int/Kr0
       else
          coeffs(i+1) = Int/Kr
       endif

    enddo
    ! $omp end parallel do

  end subroutine prg_get_chebcoeffs

  !> Gets the coefficients of the Chebyshev expansion with Ef computation.
  !! \brief In this case we are applying the bisection method to find the root.
  !! \param npts Number of points for the discretization.
  !! \param kbt Electronic temperature.
  !! \param ef Fermi level.
  !! \param tracesT Input traces for matrix polynomials.
  !! \param ncoeffs Number of Chebyshev coefficients.
  !! \param coeffs Output vector for the Chebyshev coefficients.
  !! \param emin lowest boundary for the eigenvalues of H.
  !! \param emax highest boundary for the eigenvalues of H.
  !! \param tol Tolerance for the bisection method.
  !! \param verbose Verbosity level.
  !!
  subroutine prg_get_chebcoeffs_fermi_bs(npts,kbt,ef,tracesT,ncoeffs,coeffs,emin,&
       emax,bndfil,norb,tol,jon,verbose)

    integer                  ::  i, m
    integer, intent(in)      ::  ncoeffs, norb, verbose, npts
    real(dp)                 ::  f1, f2, nel
    real(dp)                 ::  prod, step
    real(dp), intent(in)     ::  bndfil, emax, emin, kbt
    real(dp), intent(in)     ::  tol, tracesT(:)
    real(dp), intent(inout)  ::  coeffs(:), ef
    logical, intent(in)      ::  jon

    nel = bndfil*2.0_dp*real(norb,dp)
    ef=emin
    step=abs(emax-emin)
    f1=0.0_dp
    f2=0.0_dp
    prod=0.0_dp

    if(verbose >= 1)write(*,*)"In prg_get_chebcoeffs_fermi ..."

    !Sum of the occupations
    do i=1,ncoeffs
       f1 = f1 + 2.0_dp*jackson(ncoeffs,i,jon)*coeffs(i)*tracesT(i)
    enddo
    f1=f1-nel

    do m=1,1000001

       if(m.gt.1000000)then
          stop "Bisection method in prg_get_chebcoeffs_fermi is not converging ..."
       endif

       if(verbose >= 2)write(*,*)"ef,f(ef)",ef,f1
       if(abs(f1).lt.tol)then !tolerance control
          return
       endif

       ef = ef + step

       f2=0.0_dp
       call prg_get_chebcoeffs(npts,kbt,ef,ncoeffs,coeffs,emin,emax)
       !New sum of the occupations
       do i=1,ncoeffs
          f2 = f2 + 2.0_dp*jackson(ncoeffs,i,jon)*coeffs(i)*tracesT(i)
       enddo
       f2=f2-nel

       !Product to see the change in sign.
       prod = f2*f1

       if(prod.lt.0)then
          ef=ef-step
          step=step/2.0_dp !If the root is inside we shorten the step.
       else
          f1=f2  !If not, Ef moves forward.
       endif

    enddo

  end subroutine prg_get_chebcoeffs_fermi_bs

  !> Gets the coefficients of the Chebyshev expansion with Ef computation.
  !! \brief In this case the Newton-Raphson method is applied to find the root.
  !! \param npst Number of points for the discretization.
  !! \param kbt Electronic temperature.
  !! \param ef Fermi level.
  !! \param tracesT Input traces for matrix polynomials.
  !! \param ncoeffs Number of Chebyshev coefficients.
  !! \param coeffs Output vector for the Chebyshev coefficients.
  !! \param emin lowest boundary for the eigenvalues of H.
  !! \param emax highest boundary for the eigenvalues of H.
  !! \param bndfil Band filing factor.
  !! \param norb Number of orbitals.
  !! \param tol Tolerance for NR method.
  !! \param verbose Verbosity level.
  !!
  subroutine prg_get_chebcoeffs_fermi_nt(npts,kbt,ef,tracesT,ncoeffs,coeffs,emin,emax &
       ,bndfil,norb,tol,jon,verbose)

    integer                  ::  i, m
    integer, intent(in)      ::  ncoeffs, norb, verbose, npts
    real(dp)                 ::  ef0, f1, f2, nel
    real(dp)                 ::  step
    real(dp), intent(in)     ::  bndfil, emax, emin, kbt
    real(dp), intent(in)     ::  tol, tracesT(:)
    real(dp), intent(inout)  ::  coeffs(:), ef
    logical, intent(in)      ::  jon

    nel = bndfil*2.0_dp*real(norb,dp)
    ef0 = ef
    f1 = 0.0_dp
    f2 = 0.0_dp
    step = 0.1_dp

    if(verbose >= 1) write(*,*)"In prg_get_chebcoeffs_fermi_nt ..."

    !Sum of the occupations
    !$omp parallel do default(none) private(i) &
    !$omp shared(jon,coeffs,ncoeffs,tracesT) &
    !$omp reduction(+:f1)
    do i=1,ncoeffs
       f1 = f1 + 2.0_dp*jackson(ncoeffs,i,jon)*coeffs(i)*tracesT(i)
    enddo
    !$omp end parallel do

    f1=f1-nel
    ef = ef0 + step
    call prg_get_chebcoeffs(npts,kbt,ef,ncoeffs,coeffs,emin,emax)
    f2 = 0.0_dp

    !$omp parallel do default(none) private(i) &
    !$omp shared(jon,coeffs,ncoeffs,tracesT) &
    !$omp reduction(+:f2)
    do i=1,ncoeffs
       f2 = f2 + 2.0_dp*jackson(ncoeffs,i,jon)*coeffs(i)*tracesT(i)
    enddo
    !$omp end parallel do

    f2=f2-nel
    ef0 = ef
    ef = -f2*step/(f2-f1) + ef0
    f1 = f2
    step = ef - ef0

    do m = 1,1000001
       if(m.gt.1000000)then
          stop "Newton method in prg_get_chebcoeffs_fermi_nt is not converging ..."
       endif
       !New sum of the occupations
       f2 = 0.0_dp
       call prg_get_chebcoeffs(npts,kbt,ef,ncoeffs,coeffs,emin,emax)

       !$omp parallel do default(none) private(i) &
       !$omp shared(jon,coeffs,ncoeffs,tracesT) &
       !$omp reduction(+:f2)
       do i=1,ncoeffs
          f2 = f2 + 2.0_dp*jackson(ncoeffs,i,jon)*coeffs(i)*tracesT(i)
       enddo
       !$omp end parallel do

       f2=f2-nel
       if(verbose >= 2) write(*,*)"ef,f(ef)",ef,f2
       ef0 = ef
       ef = -f2*step/(f2-f1) + ef0
       f1 = f2
       step = ef - ef0
       if(abs(f1).lt.tol)then !tolerance control
          return
       endif
    enddo

  end subroutine prg_get_chebcoeffs_fermi_nt

  !> Chebyshev polynomial obtained by recursion.
  !! \param r rth polynomial.
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

  !> Gets the absolute maximum of the derivative of a function.
  !! \param func.
  !! \param de Energy step.
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
