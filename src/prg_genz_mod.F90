!> To produce a matrix \f$Z\f$ which is needed to orthogonalize \f$H\f$.
!! \ingroup PROGRESS
!! \brief \f$ H_{orth} = Z^{\dagger}HZ \f$
!! See Negre 2016 \cite Negre2016
!!
module prg_genz_mod

  use prg_openfiles_mod
  use bml
  use prg_kernelparser_mod
  use prg_parallel_mod
  use prg_extras_mod

  implicit none

  private  !Everything is private by default

  integer, parameter :: dp = kind(1.0d0)

  public :: prg_buildZdiag, prg_buildZsparse, prg_parse_ZSP, prg_init_ZSPmat
  public :: prg_genz_sp_initial_zmat, prg_genz_sp_ref, prg_genz_sp_initialz0

  !> Input for the genz driver.
  !! This type controlls all the variables that are needed by genz
  !!
  type, public :: genZSPinp  !< The ZSpinp data type

     !> To have different levels of verbose
     integer :: verbose

     !> !Lentgth of the "firsts iteration period".
     integer :: nfirst

     !> !Initial number of recursive refinements.
     integer :: nrefi

     !> !Initial number of recursive refinements.
     integer :: nreff

     !> Initial threshold value.
     real(dp) :: numthresi

     !> Final threshold value.
     real(dp) :: numthresf

     !> If we want to do XL integration scheme for Z.
     logical :: integration

     !> To keep track of the genz iterations
     integer :: igenz

     !> Logical variable to compute in sparse or dense mode
     logical :: ZSP

     !> Max nonzero elements per row for every row see \cite Mniszewski2015 .
     integer :: mdim

     !> Matrix format (Dense or Ellpack).
     character(20) :: bml_type

  end type genZSPinp

contains

  !> The parser for genz solver.
  !!
  subroutine prg_parse_ZSP(input,filename)
    type(genZSPinp), intent(inout) :: input
    integer, parameter :: nkey_char = 1, nkey_int = 5, nkey_re = 2, nkey_log = 2
    character(len=*) :: filename

    !Library of keywords with the respective defaults.
    character(len=50), parameter :: keyvector_char(nkey_char) = [character(len=100) :: &
         'BMLType=']
    character(len=100) :: valvector_char(nkey_char) = [character(len=100) :: &
         'Dense']

    character(len=50), parameter :: keyvector_int(nkey_int) = [character(len=50) :: &
         'Verbose=','NFirst=','NRefI=','NRefF=','MDim=']
    integer :: valvector_int(nkey_int) = (/ &
         0,10,3,1,-1 /)

    character(len=50), parameter :: keyvector_re(nkey_re) = [character(len=50) :: &
         'NumthreshI=','NumthreshF=' ]
    real(dp) :: valvector_re(nkey_re) = (/&
         1.0e-8, 1.0e-5  /)

    character(len=50), parameter :: keyvector_log(nkey_log) = [character(len=100) :: &
         'ZSP=','Int=']
    logical :: valvector_log(nkey_log) = (/&
         .false., .true. /)

    !Start and stop characters
    character(len=50), parameter :: startstop(2) = [character(len=50) :: &
         'ZSP{', '}']

    call prg_parsing_kernel(keyvector_char,valvector_char&
         ,keyvector_int,valvector_int,keyvector_re,valvector_re,&
         keyvector_log,valvector_log,trim(filename),startstop)

    !Characters
    if(valvector_char(1) == "Dense")then
       input%bml_type = bml_matrix_dense
    elseif(valvector_char(1) == "Ellpack")then
       input%bml_type = bml_matrix_ellpack
    endif

    !Integers
    input%verbose = valvector_int(1)
    input%nfirst = valvector_int(2)
    input%nrefi =  valvector_int(3)
    input%nreff = valvector_int(4)
    input%mdim = valvector_int(5)

    !Logicals
    input%zsp = valvector_log(1)
    input%integration = valvector_log(2)

    !Reals
    input%numthresi = valvector_re(1)
    input%numthresf = valvector_re(2)

  end subroutine prg_parse_ZSP

  !> Initiates the matrices for the Xl integration of Z
  !! \param self input zsp variables
  !! \param zk1_bml-zk6_bml history record of the previous Z matrices.
  !! \param norb number of orbitals.
  !! \param bml_type the bml format we are passing.
  !!
  subroutine prg_init_ZSPmat(igenz,zk1_bml,zk2_bml,zk3_bml&
       ,zk4_bml,zk5_bml,zk6_bml,norb,bml_type,bml_element_type)
    integer :: norb, igenz
    character(20) :: bml_type, my_bml_element_type
    character(20), optional :: bml_element_type
    type(bml_matrix_t) :: zk1_bml
    type(bml_matrix_t) :: zk2_bml
    type(bml_matrix_t) :: zk3_bml
    type(bml_matrix_t) :: zk4_bml
    type(bml_matrix_t) :: zk5_bml
    type(bml_matrix_t) :: zk6_bml

    if(present(bml_element_type))then
      my_bml_element_type = bml_element_type
    else
      my_bml_element_type = bml_element_real
    endif

    igenz = 0

    if(bml_get_n(zk1_bml).le.0)then
       call bml_zero_matrix(bml_type,my_bml_element_type,dp,norb,norb,zk1_bml)
       call bml_zero_matrix(bml_type,my_bml_element_type,dp,norb,norb,zk2_bml)
       call bml_zero_matrix(bml_type,my_bml_element_type,dp,norb,norb,zk3_bml)
       call bml_zero_matrix(bml_type,my_bml_element_type,dp,norb,norb,zk4_bml)
       call bml_zero_matrix(bml_type,my_bml_element_type,dp,norb,norb,zk5_bml)
       call bml_zero_matrix(bml_type,my_bml_element_type,dp,norb,norb,zk6_bml)
    endif

  end subroutine prg_init_ZSPmat

  !> Usual subroutine involving diagonalization.
  !! \f$ Z=U\sqrt{s}U^{\dagger} \f$, where \f$ U \f$ = eigenvectors and \f$ s \f$ = eigenvalues.
  !! The purpose of this subroutine is to have an exact way of computing
  !! z for comparing with the sparse approach.
  !! \param smat_bml Overlap matrix in bml format.
  !! \param zmat_bml Congruence transform in bml format.
  !! \param threshold Threshold value to use, in this case, only in the backtransformation to ellpack format.
  !! \param mdim Maximun nonzero to use, in this case, only in the backtransformation to ellpack format.
  !! \param bml_type the bml type we are passing.
  !!
  subroutine prg_buildZdiag(smat_bml,zmat_bml,threshold,mdimin,bml_type,verbose)
    !     use extras
    real(dp)                           ::  err_check
    character(len=*), intent(in)       ::  bml_type
    character(20)                      ::  bml_element_type
    integer                            ::  i, j, mdim, norb
    integer, intent(in)                ::  mdimin
    integer, optional, intent(in)      ::  verbose
    real(dp)                           ::  mls_i
    real(dp)                           ::  invsqrt, threshold
    real(dp), allocatable              ::  nono_evals(:), nonotmp(:,:), smat(:,:), umat(:,:)
    real(dp), allocatable              ::  zmat(:,:)
    type(bml_matrix_t)                 ::  nonotmp_bml, saux_bml, umat_bml, umat_t_bml
    type(bml_matrix_t)                 ::  zmat_bml
    type(bml_matrix_t), intent(inout)  ::  smat_bml

    if(present(verbose))then 
       if(verbose >= 1)then
          write(*,*)""; write(*,*)"In buildzdiag ..."
       endif
    endif

    if(bml_get_precision(smat_bml) == 1 .or.&
      &bml_get_precision(smat_bml) == 2)then
      bml_element_type = "real"
    elseif(bml_get_precision(smat_bml) == 3 .or.&
      &bml_get_precision(smat_bml) == 4)then
      bml_element_type = "complex"
    endif

    norb = bml_get_n(smat_bml)
    mdim = bml_get_m(smat_bml)

    !Allocate temporary matrices.
    allocate(nono_evals(norb))
    allocate(umat(norb,norb))
    allocate(nonotmp(norb,norb))
    allocate(zmat(norb,norb))
    allocate(smat(norb,norb))
    !To bml dense. this is done because the diagonalization
    !it is only implemented for bml_dense. In future versions of bml
    !the api should do this automatically.
    call bml_export_to_dense(smat_bml,smat) !my_bml_type to dense

    call bml_zero_matrix(bml_matrix_dense,bml_element_type,dp,norb,norb,saux_bml) !Allocate bml dense

    call bml_import_from_dense(bml_matrix_dense,smat,saux_bml,threshold,mdim) !Dense to dense_bml

    !call bml_print_matrix("Smat_bml",smat_bml,0,6,0,6)
    ! call bml_print_matrix("Smat",saux_bml,0,6,0,6)

    !Reseting zmat to make it bml dense. Same reason as before.
    call bml_deallocate(zmat_bml)
    call bml_zero_matrix(bml_matrix_dense,bml_element_type,dp,norb,norb,zmat_bml)

    !Auxiliary matrices.
    call bml_zero_matrix(bml_matrix_dense,bml_element_type,dp,norb,norb,umat_bml)
    call bml_zero_matrix(bml_matrix_dense,bml_element_type,dp,norb,norb,umat_t_bml)
    call bml_zero_matrix(bml_matrix_dense,bml_element_type,dp,norb,norb,nonotmp_bml)

    !Eigenvectors and eigenalues of the overlap s.
    mls_i = mls()
    call bml_diagonalize(saux_bml,nono_evals,umat_bml)
    if(present(verbose))then
      if(verbose >= 1)then
         write(*,*)"Time for S diag = "//to_string(mls() - mls_i)//" ms"
      endif
    endif

    call bml_export_to_dense(umat_bml, umat)

    nonotmp=0.0_dp

    !Doing u s^-1/2
    do i = 1, norb

       if(nono_evals(i).lt.0.0_dp) stop 'matrix s has a 0 eigenvalue'

       invsqrt = 1.0_dp/sqrt(nono_evals(i))

       do j = 1, norb
          nonotmp(j,i) = umat(j,i) * invsqrt
       end do

    end do

    !Computing u^dag
    call bml_transpose(umat_bml,umat_t_bml)

    ! #ifdef DO_MPI
    !     if (getNRanks() .gt. 1 .and. &
    !         bml_get_distribution_mode(umat_t_bml) .eq. BML_DMODE_DISTRIBUTED) then
    !         call prg_allGatherParallel(umat_t_bml)
    !     endif
    ! #endif

    call bml_import_from_dense(BML_MATRIX_DENSE, nonotmp, nonotmp_bml)

    !Doing u s^-1/2 u^t
    !!call bml_multiply(nonotmp_bml,umat_t_bml, zmat_bml, 1.0_dp, 1.0_dp,threshold)
    call bml_multiply(nonotmp_bml,umat_t_bml, zmat_bml, 1.0_dp, 0.0_dp,threshold)

    !If the original type was ellpack then we convert back from
    !dense to ellpack. This is done just to be able to test ellpack with sp2 and buildzdiag.
    if(bml_type.eq."ellpack")then
       call bml_export_to_dense(zmat_bml, zmat)!Dense_bml to dense.
       call bml_deallocate(zmat_bml)
       call bml_import_from_dense(bml_matrix_ellpack,zmat,zmat_bml,threshold,mdim, bml_get_distribution_mode(smat_bml)) !Dense to ellpack_bml.
    endif

    !!To check for the accuracy of the approximation (prg_delta). this is done using matmul
    !!so its very inefficient. Only uncomment for debugging purpose.
    !call bml_export_to_dense(zmat_bml, zmat)
    !call prg_delta(zmat,smat,norb,err_check)
    !write(*,*)"err", err_check, norb
    !stop

    deallocate(nonotmp)
    deallocate(nono_evals)
    deallocate(umat)
    deallocate(smat)
    deallocate(zmat)
    call bml_deallocate(umat_bml)
    call bml_deallocate(umat_t_bml)
    call bml_deallocate(saux_bml)
    call bml_deallocate(nonotmp_bml)

  end subroutine prg_buildZdiag

  !> Inverse factorization using Niklasson's algorithm.
  !! \param smat_bml overlap matrix
  !! \param zmat_bml congruence transform to be updated or computed. (bml format)
  !! \param igenz counter to keep track of the calls to this subroutine.
  !! \param mdim dimension of the maxnonzero per row.
  !! \param zk1_bml-zk6_bml: history of the past congruence transforms.
  !! \param nfirst first pre septs with nrefi and thresholdi.
  !! \param nrefi number of refinement iterations for the firsts "nfirst" steps.
  !! \param nreff number of refinement iterations for the rest of the steps.
  !! \param integration if we want to apply xl integration scheme for z (default is always .true.)
  !! \param verbose to print extra information.
  !!
  subroutine prg_buildZsparse(smat_bml,zmat_bml,igenz,mdim,&
       bml_type,zk1_bml,zk2_bml,zk3_bml&
       ,zk4_bml,zk5_bml,zk6_bml,nfirst,nrefi,nreff&
       ,thresholdi,thresholdf,integration,verbose)
    !use extras
    !real(dp)              ::  err_check
    !real(dp), allocatable  ::  smat(:,:), zmat(:,:)
    character(20)            ::  bml_type
    integer                  ::  KK, i, igenz, mdim
    integer                  ::  nfirst, norb, nref, nreff
    integer                  ::  nrefi, verbose
    logical                  ::  integration
    real(dp)                  ::  sec_i, sec_ii
    real(dp)                 ::  threshold, thresholdf, thresholdi
    type(bml_matrix_t)       ::  smat_bml, zk1_bml, zk2_bml, zk3_bml
    type(bml_matrix_t)       ::  zk4_bml, zk5_bml, zk6_bml, zmat_bml

    norb = bml_get_n(smat_bml)
    if(mdim<0) mdim = norb

    KK=6 !Total number of stored z matrices for the xl integration scheme.

    if(verbose.eq.1)sec_ii=mls() !Gets the actual time in ms.

    ! Number of refinements and threshold value for the first "nfirst" iterations.
    if(igenz.lt.nfirst)then
       nref=nrefi; threshold = thresholdi !For the firsts iterations.
    else
       nref=nreff; threshold = thresholdf !For the following iterations.
    end if

    if(verbose.eq.1)sec_i=mls() !Firs calculation of z using the graph approach.
    if(igenz.eq.1) call prg_genz_sp_initialz0(smat_bml,zmat_bml,norb,mdim,bml_type,threshold)
    if(verbose.eq.1)write(*,*)"Time for prg_initial estimate "//to_string(mls()-sec_i)//" ms"

    if(verbose.eq.1)sec_i=mls()! integration scheme.
    if(integration)then
       call prg_genz_sp_int(zmat_bml,zk1_bml,zk2_bml,zk3_bml&
            &,zk4_bml,zk5_bml,zk6_bml,igenz,norb,bml_type,threshold)
    end if
    if(verbose.eq.1)write(*,*)"Time for xl scheme "//to_string(mls()-sec_i)//" ms"

    if(verbose.eq.1)sec_i=mls()! Refinement.
    call prg_genz_sp_ref(smat_bml,zmat_bml,nref,mdim,bml_type,threshold)
    if(verbose.eq.1)write(*,*)"Time for prg_genz_sp_ref "//to_string(mls()-sec_i)//" ms"

    !     call bml_export_to_dense()
    !     call prg_delta(xmat,smat,norb,err_check)  !to check for the accuracy of the approximation (prg_delta)
    !     write(*,*)"err", err_check, norb
    !     stop

    if(verbose.eq.1)write(*,*)"Time for prg_buildZsparse "//to_string(igenz)//" "//to_string(mls()-sec_ii)//" ms"

  end subroutine prg_buildZsparse

  !> Initial estimation of Z.
  !! \note Most of the operations are done in pure dense format.
  !! The purpose of this subroutine is to have an exact way of computing
  !! z for comparing with the sparse approach.
  !! \param smat_bml Overlap matrix in bml format.
  !! \param zmat_bml Congruence transform in bml format.
  !! \param norb Congruence transform in bml format.
  !! \param mdim Congruence transform in bml format.
  !! \param bml_type_f The bml final type of zmat_bml.
  !! \param threshold Threshold value to use, in this case, only in the backtransformation to ellpack format.
  !!
  subroutine prg_genz_sp_initialz0(smat_bml,zmat_bml,norb,mdim,bml_type_f,threshold)
    character(20)                   ::  bml_type, bml_type_f, bml_element_type
    integer                         ::  i, ii, j, jj
    integer                         ::  k, l, mdim, norb
    integer, allocatable            ::  kk(:)
    real(dp)                        ::  err_check, invsqrt, threshold, thresholdz0
    real(dp), allocatable           ::  sitmp(:,:), smat(:,:), sthres(:,:), stmp(:,:)
    real(dp), allocatable           ::  stmp_evals(:), utmp(:,:), zmat(:,:), ztmp(:,:)
    type(bml_matrix_t)              ::  sitmp_bml, sthres_bml, stmp_bml, utmp_bml
    type(bml_matrix_t)              ::  xtmp_bml, zmat_bml
    type(bml_matrix_t), intent(in)  ::  smat_bml

    !When this subroutine is parallelized, the following threshold can be set to larger (coarser) values.
    !The reason for this is to reduce the memory per thread.
    thresholdz0= 1.0d-10 !This could be passed as an input.

    bml_type= bml_matrix_dense !All the operations are performed in bml_dense.
    ! bml_element_type = bml_get_element_type(smat_bml)
    if(bml_get_precision(smat_bml) == 1 .or.&
      &bml_get_precision(smat_bml) == 2)then
      bml_element_type = "real"
    elseif(bml_get_precision(smat_bml) == 3 .or.&
      &bml_get_precision(smat_bml) == 4)then
      bml_element_type = "complex"
    endif

    !Part of the operations are still done in pure dense format.
    allocate(zmat(norb,norb))
    zmat = 0.0_dp

    allocate(sthres(norb,norb))

    allocate(smat(norb,norb))
    call bml_export_to_dense(smat_bml,smat)

    call bml_zero_matrix(bml_type,bml_element_type,dp,norb,norb,sthres_bml) !a thresholded version of s

    call bml_import_from_dense(bml_type,smat,sthres_bml,thresholdZ0,mdim)

    call bml_export_to_dense(sthres_bml,sthres) !We will use the dense to extract the small blocks.

    allocate(kk(norb)) !The nonzero positions in the row

    !NOTE: The following loop must be OMP parallelized as parallelization is trivial, however, as BML
    !operations are taking all the threads this is not possible. BML needs to recognize automatically
    !when "threads" are requested by the hosting code.
    do i = 1, norb !Z0 is the prg_initial guess for the iterative refinement.

       jj=0
       kk=0
       do j=1,norb
          if(abs(sthres(i,j)).gt.thresholdZ0)then
             jj=jj+1
             kk(jj)=j
          endif
       enddo
       k=jj

       !The followings will be dense matrices since they are the small blocks.
       call bml_zero_matrix(bml_type,bml_element_type,dp, k,k,stmp_bml)
       call bml_zero_matrix(bml_type,bml_element_type,dp, k,k,sitmp_bml)
       call bml_zero_matrix(bml_type,bml_element_type,dp, k,k,utmp_bml)
       call bml_zero_matrix(bml_type,bml_element_type,dp, k,k,xtmp_bml)

       allocate(stmp(k,k),sitmp(k,k))
       allocate(utmp(k,k),stmp_evals(k),ztmp(k,k))

       stmp = 0.0_dp !Small dense block to be extracted from s.
       stmp_evals = 0.0_dp !eigenvalues for the small dense s matrices.
       sitmp = 0.0_dp
       utmp = 0.0_dp !Matrix containing the eigenvectors of small s
       ztmp = 0.0_dp

       do j = 1, k !Extracting the small s matrices
          do l = 1, k
             stmp(l,j) = smat(kk(l),kk(j))
          end do
       end do

       call bml_import_from_dense(bml_type, stmp, stmp_bml)
       call bml_diagonalize(stmp_bml,stmp_evals,utmp_bml)
       call bml_export_to_dense(utmp_bml,utmp)

       do j = 1, k  !Applying the inverse sqrt function to the eigenvalues.
          invsqrt = 1.0_dp/sqrt(stmp_evals(j))
          do l = 1, k
             sitmp(l,j) = utmp(l,j) * invsqrt
          end do
       end do

       utmp=transpose(utmp)

       call bml_import_from_dense(bml_type , utmp , utmp_bml)
       call bml_import_from_dense(bml_type, sitmp, sitmp_bml)
       call bml_multiply(sitmp_bml,utmp_bml,xtmp_bml,1.0_dp,0.0_dp)
       call bml_export_to_dense(xtmp_bml,ztmp)

       do l = 1, k  !Reconstructing the large Z0 based in the small Z.
          do j = 1, k
             if(kk(j).eq.i)then !For the row corresponding to the extracted dense block.
                zmat(i,kk(l)) = ztmp(j,l)
             end if
          end do
       end do

       deallocate(stmp,sitmp)  !Deallocate temporary matrices
       deallocate(utmp)
       deallocate(stmp_evals)
       deallocate(ztmp)
       call bml_deallocate(sitmp_bml)
       call bml_deallocate(stmp_bml)
       call bml_deallocate(utmp_bml)

    end do

    call bml_zero_matrix(bml_type_f,bml_element_type,dp,norb,norb,zmat_bml)
    call bml_import_from_dense(bml_type_f,zmat, zmat_bml,threshold,mdim) !Converting to bml format

    ! call prg_delta(zmat,smat,norb,err_check)  !to check for the accuracy of the approximation (prg_delta)
    ! call sparsity(smat,norb,spa)
    ! write(*,*)"err", err_check, norb
    ! stop

    deallocate(sthres)
    deallocate(kk)
    call bml_deallocate(sthres_bml)

  end subroutine prg_genz_sp_initialz0

  !> Initial estimation of Z.
  !! \note Most of the operations are done in pure dense format.
  !! The purpose of this subroutine is to have an exact way of computing
  !! z for comparing with the sparse approach.
  !! \param smat_bml Overlap matrix in bml format.
  !! \param zmat_bml Congruence transform in bml format.
  !! \param norb Congruence transform in bml format.
  !! \param mdim Congruence transform in bml format.
  !! \param bml_type_f The bml final type of zmat_bml.
  !! \param threshold Threshold value to use, in this case, only in the backtransformation to ellpack format.
  !!
  subroutine prg_genz_sp_initial_zmat(smat_bml,zmat_bml,norb,mdim,bml_type_f,threshold)
    !     use extras
    implicit none
    character(20)                   ::  bml_type, bml_type_f, bml_element_type
    integer                         ::  i, ii, j, jj
    integer                         ::  k, l, mdim, norb
    integer, allocatable            ::  kk(:)
    real(dp), intent(in)            ::  threshold
    real(dp)                        ::  err_check, invsqrt
    real(dp)                        ::  thresholdz0
    real(dp), allocatable           ::  sitmp(:,:), smat(:,:), sthres(:,:)
    real(dp), allocatable           ::  stmp(:,:)
    real(dp), allocatable           ::  stmp_evals(:), utmp(:,:)
    real(dp), allocatable           ::  zmat(:,:), ztmp(:,:)
    type(bml_matrix_t)              ::  sitmp_bml, sthres_bml, stmp_bml
    type(bml_matrix_t)              ::  xtmp_bml, zmat_bml, utmp_bml
    type(bml_matrix_t), intent(in)  ::  smat_bml

    !When this subroutine is parallelized, the following threshold can be set to
    !larger (coarser) values.
    !The reason for this is to reduce the memory per thread.
    !thresholdz0= 1.0d-10 !This could be passed as an input.

    bml_type= bml_matrix_dense !All the operations are performed in bml_dense.
    ! bml_element_type = bml_get_element_type(smat_bml)
    if(bml_get_precision(smat_bml) == 1 .or.&
      &bml_get_precision(smat_bml) == 2)then
      bml_element_type = "real"
    elseif(bml_get_precision(smat_bml) == 3 .or.&
      &bml_get_precision(smat_bml) == 4)then
      bml_element_type = "complex"
    endif

    !Part of the operations are still done in pure dense format.
    allocate(zmat(norb,norb))
    allocate(sthres(norb,norb))

    allocate(smat(norb,norb))
    call bml_export_to_dense(smat_bml,smat)

    ! a thresholded version of s
    call bml_zero_matrix(bml_type,bml_element_type,dp,norb,norb,sthres_bml)

    !call bml_import_from_dense(bml_type,smat,sthres_bml,thresholdZ0,mdim)
    call bml_import_from_dense(bml_type,smat,sthres_bml,threshold,mdim)

    ! We will use the dense to extract the small blocks.
    call bml_export_to_dense(sthres_bml,sthres)

    allocate(kk(norb)) !The nonzero positions in the row

    !NOTE: The following loop must be OMP parallelized as parallelization is
    !trivial, however, as BML
    !operations are taking all the threads this is not possible. BML needs to
    !recognize automatically
    !when "threads" are requested by the hosting code.
    do i = 1, norb !Z0 is the prg_initial guess for the iterative refinement.
       !do i = bml_getLocalRowMin(smat_bml, getMyRank()), &
       !       bml_getLocalRowMax(smat_bml, getMyRank())

       jj=0
       kk=0
       do j=1,norb
          !if(abs(sthres(i,j)).gt.thresholdZ0)then
          if(abs(sthres(i,j)).gt.threshold)then
             jj=jj+1
             kk(jj)=j
          endif
       enddo
       k=jj

       !The followings will be dense matrices since they are the small blocks.
       call bml_zero_matrix(bml_type,bml_element_type,dp, k,k,stmp_bml)
       call bml_zero_matrix(bml_type,bml_element_type,dp, k,k,sitmp_bml)
       call bml_zero_matrix(bml_type,bml_element_type,dp, k,k,utmp_bml)
       call bml_zero_matrix(bml_type,bml_element_type,dp, k,k,xtmp_bml)

       allocate(stmp(k,k),sitmp(k,k))
       allocate(utmp(k,k),stmp_evals(k),ztmp(k,k))

       stmp = 0.0_dp !Small dense block to be extracted from s.
       stmp_evals = 0.0_dp !eigenvalues for the small dense s matrices.
       sitmp = 0.0_dp
       utmp = 0.0_dp !Matrix containing the eigenvectors of small s
       ztmp = 0.0_dp

       do j = 1, k !Extracting the small s matrices
          do l = 1, k
             stmp(l,j) = smat(kk(l),kk(j))
          end do
       end do

       call bml_import_from_dense(BML_MATRIX_DENSE, stmp, stmp_bml)
       call bml_diagonalize(stmp_bml,stmp_evals,utmp_bml)
       call bml_export_to_dense(utmp_bml,utmp)

       do j = 1, k  !Applying the inverse sqrt function to the eigenvalues.
          invsqrt = 1.0_dp/sqrt(stmp_evals(j))
          do l = 1, k
             sitmp(l,j) = utmp(l,j) * invsqrt
          end do
       end do

       utmp=transpose(utmp)

       call bml_import_from_dense(bml_type, utmp, utmp_bml)
       call bml_import_from_dense(bml_type, sitmp, sitmp_bml)
       call bml_multiply(sitmp_bml,utmp_bml,xtmp_bml,1.0_dp,0.0_dp)
       call bml_export_to_dense(xtmp_bml,ztmp)

       do l = 1, k  !Reconstructing the large Z0 based in the small Z.
          do j = 1, k
             if(kk(j).eq.i)then !For the row corresponding to the extracted dense
                !block.
                zmat(i,kk(l)) = ztmp(j,l)
             end if
          end do
       end do

       deallocate(stmp,sitmp)  !Deallocate temporary matrices
       deallocate(utmp)
       deallocate(stmp_evals)
       deallocate(ztmp)
       call bml_deallocate(sitmp_bml)
       call bml_deallocate(stmp_bml)
       call bml_deallocate(utmp_bml)

    end do

    !call bml_zero_matrix(bml_type_f,bml_element_type,bml_element_precision,norb,norb,zmat_bml)
    call bml_import_from_dense(bml_type_f,zmat,zmat_bml,threshold,mdim, &
         bml_get_distribution_mode(smat_bml))
    !Converting to bml format

    ! call prg_delta(zmat,smat,norb,err_check)  !to check for the accuracy of
    ! the approximation (prg_delta)
    ! call sparsity(smat,norb,spa)
    ! write(*,*)"err", err_check, norb
    ! stop
    deallocate(sthres)
    deallocate(kk)
    call bml_deallocate(sthres_bml)

  end subroutine prg_genz_sp_initial_zmat

  !> Inverse factorization using Niklasson's algorithm.
  !! \param smat_bml overlap matrix
  !! \param zmat_bml congruence transform to be updated or computed. (bml format)
  !! \param mdim dimension of the maxnonzero per row.
  !! \param zk1_bml-zk6_bml: history of the past congruence transforms.
  !! \param igenz counter to keep track of the calls to this subroutine.
  !! \param norb Congruence transform in bml format.
  !! \param bml_type_f The bml final type of zmat_bml.
  !! \param threshold Threshold value to use.
  !!
  subroutine prg_genz_sp_int(zmat_bml,zk1_bml,zk2_bml,zk3_bml&
       &,zk4_bml,zk5_bml,zk6_bml,igenz,norb,bml_type&
       &,threshold)
    integer :: igenz,norb,KK
    character(20) :: bml_element_type
    real(dp) :: alpha, kappa, c0, c1, c2, c3, c4, c5
    real(dp) :: threshold
    type(bml_matrix_t) :: zmat_bml
    type(bml_matrix_t) :: zk1_bml,zk2_bml,zk3_bml,zk4_bml,zk5_bml,zk6_bml
    character(20) :: bml_type

    KK=6
    alpha=0.0180_dp
    kappa=1.82_dp

    ! bml_element_type = bml_get_element_type(zmat_bml)
    if(bml_get_precision(zmat_bml) == 1 .or.&
      &bml_get_precision(zmat_bml) == 2)then
      bml_element_type = "real"
    elseif(bml_get_precision(zmat_bml) == 3 .or.&
      &bml_get_precision(zmat_bml) == 4)then
      bml_element_type = "complex"
    endif

    !The following constants are the original constants premultiplied by alpha.

    C0=-0.1080_dp
    C1=0.2520_dp
    C2=-0.1440_dp
    C3=-0.0540_dp
    C4=0.0720_dp
    C5=-0.0180_dp


    if(igenz.eq.kk)then

       call bml_zero_matrix(bml_type,bml_element_type,dp,norb ,norb,zk1_bml)
       call bml_zero_matrix(bml_type,bml_element_type,dp,norb ,norb,zk2_bml)
       call bml_zero_matrix(bml_type,bml_element_type,dp,norb ,norb,zk3_bml)
       call bml_zero_matrix(bml_type,bml_element_type,dp,norb ,norb,zk4_bml)
       call bml_zero_matrix(bml_type,bml_element_type,dp,norb ,norb,zk5_bml)
       call bml_zero_matrix(bml_type,bml_element_type,dp,norb ,norb,zk6_bml)

       call bml_copy(zmat_bml,zk1_bml);call bml_copy(zmat_bml,zk2_bml);call bml_copy(zmat_bml,zk3_bml)
       call bml_copy(zmat_bml,zk4_bml);call bml_copy(zmat_bml,zk5_bml);call bml_copy(zmat_bml,zk6_bml)

    end if

    if(igenz.ge.kk)then !Here we change z by applying the t-r integration scheme

       call bml_add_deprecated(1.0_dp,zmat_bml,-1.0_dp,zk6_bml,threshold)  !Z(t)-Z~(t)
       call bml_scale(kappa,zmat_bml)
       call bml_add_deprecated(1.0_dp,zmat_bml,2.0_dp,zk6_bml,threshold)   !Z(t)+2Z~(t)
       call bml_add_deprecated(1.0_dp,zmat_bml,-1.0_dp,zk5_bml,threshold)  !Z(t)-Z~(t-dt)

       !Dissipation force term:
       call bml_add_deprecated(1.0_dp,zmat_bml,c0,zk6_bml,threshold) !Z(t)+c0*Z~(t)
       call bml_add_deprecated(1.0_dp,zmat_bml,c1,zk5_bml,threshold)
       call bml_add_deprecated(1.0_dp,zmat_bml,c2,zk4_bml,threshold)
       call bml_add_deprecated(1.0_dp,zmat_bml,c3,zk3_bml,threshold)
       call bml_add_deprecated(1.0_dp,zmat_bml,c4,zk2_bml,threshold)
       call bml_add_deprecated(1.0_dp,zmat_bml,c5,zk1_bml,threshold) !Z(t)+c0*Z~(t-5*dt)

    end if

    if(igenz.ge.kk)then !Here we are shifting the z matrices.

       call bml_copy(zk2_bml,zk1_bml) !Z~(t-5*dt)=Z~(t-4*dt)
       call bml_copy(zk3_bml,zk2_bml) !Z~(t-4*dt)=Z~(t-3*dt)
       call bml_copy(zk4_bml,zk3_bml)
       call bml_copy(zk5_bml,zk4_bml)
       call bml_copy(zk6_bml,zk5_bml) !Z~(t-dt)=Z~(t)
       call bml_copy(zmat_bml,zk6_bml) !Z~(t)=Z~(t+dt)

    end if

  end subroutine prg_genz_sp_int

  !> Iterative refinement.
  !! \param smat_bml overlap matrix
  !! \param zmat_bml congruence transform to be updated or computed. (bml format)
  !! \param nref Number of refinement iterations.
  !! \param bml_type_f The bml final type of zmat_bml.
  !! \param threshold Threshold value to use.
  !! \param verbose to print extra information.
  !!
  subroutine prg_genz_sp_ref(smat_bml,zmat_bml,nref,norb,bml_type,threshold)
    integer :: k
    integer, intent(inout) :: norb
    integer, intent(in) :: NREF
    real(dp) :: err_check, sec_i
    real(dp), allocatable :: smat(:,:), zmat(:,:)
    real(dp),intent(in) :: threshold
    type(bml_matrix_t) :: temp_bml
    type(bml_matrix_t) :: temp1_bml
    type(bml_matrix_t) :: temp2_bml
    type(bml_matrix_t) :: idscaled_bml
    type(bml_matrix_t) :: xmat_t_bml
    type(bml_matrix_t),intent(inout) :: zmat_bml
    type(bml_matrix_t) :: aux_bml
    type(bml_matrix_t), intent(in) :: smat_bml
    character(20),intent(in) :: bml_type
    character(20) :: bml_element_type

    norb = bml_get_n(smat_bml)

    ! bml_element_type = bml_get_element_type(smat_bml)
    if(bml_get_precision(smat_bml) == 1 .or.&
      &bml_get_precision(smat_bml) == 2)then
      bml_element_type = "real"
    elseif(bml_get_precision(smat_bml) == 3 .or.&
      &bml_get_precision(smat_bml) == 4)then
      bml_element_type = "complex"
    endif

    call bml_zero_matrix(bml_type,bml_element_type,dp,norb,norb,idscaled_bml)

    call bml_add_identity(idscaled_bml, 1.0_dp, threshold)  !1.0 [0] + 1.0 I
    call bml_scale(1.8750_dp,idscaled_bml) ! 1.875*I

    call bml_noinit_matrix(bml_type,bml_element_type,dp,norb ,norb,temp_bml)
    call bml_noinit_matrix(bml_type,bml_element_type,dp,norb ,norb,temp2_bml)

    sec_i=mls() !Firs calculation of z using the graph approach.
    do k = 1, NREF !Iterative refinement

       !Enforcing symmetry (in bml).
       call bml_transpose(zmat_bml, xmat_t_bml) !Z^t
       call bml_add_deprecated(0.50_dp,zmat_bml, 0.50_dp, xmat_t_bml,threshold) !(Z^t+Z)/2

       call bml_multiply(smat_bml,zmat_bml,temp_bml, 1.0_dp, 0.0_dp,threshold)  !S*Z

       call bml_multiply(zmat_bml,temp_bml, temp2_bml, 1.0_dp, 0.0_dp,threshold)  !X = Z^t*S*Z

       call bml_multiply(temp2_bml, temp2_bml, temp_bml, 1.0_dp, 0.0_dp,threshold) !X*X

       call bml_scale(0.3750_dp, temp_bml)
       call bml_scale(-1.250_dp, temp2_bml)

       !Temp = 1.875*I - 1.25*X + 0.375*X^2
       call bml_add_deprecated(1.0_dp,temp2_bml,1.0_dp, idscaled_bml,threshold)
       call bml_add_deprecated(1.0_dp,temp_bml,1.0_dp, temp2_bml,threshold)

       call bml_multiply(zmat_bml,temp_bml,temp2_bml, 1.0_dp, 0.0_dp,threshold) !Z*Temp

       call bml_copy(temp2_bml,zmat_bml)

    end do
    call bml_deallocate(temp_bml)
    call bml_deallocate(temp1_bml)
    call bml_deallocate(temp2_bml)

    write(*,*)"Time for ref loop "//to_string(mls()-sec_i)//" ms"


    call bml_deallocate(idscaled_bml)
    call bml_deallocate(xmat_t_bml)

  end subroutine prg_genz_sp_ref

end module prg_genz_mod
