!> The SP2 module.
!! \ingroup PROGRESS
!!
!! \author S. Mniszewski
!! (smn@lanl.gov)
!!
!! \brief This subroutine implements Niklasson's SP2 density matrix purification
!! algorithm.
!!
module sp2_mod

  use bml   
  use normalize_mod
  use timer_mod
  use parallel_mod
  
  implicit none
 
  private  !Everything is private by default
  
  integer, parameter :: dp = kind(1.0d0)

  public :: sp2_basic
  public :: sp2_alg1
  public :: sp2_alg1_genseq
  public :: sp2_alg1_seq
  public :: sp2_alg1_seq_inplace
  public :: sp2_alg2
  public :: sp2_alg2_genseq
  public :: sp2_alg2_seq
  public :: sp2_alg2_seq_inplace
  public :: sp2_submatrix
  public :: sp2_submatrix_inplace
  
contains

  !> Calculates the density matrix from a Hamiltonian matrix by
  !! purification. The method implemented here is the very first verion of the
  !! SP2 method.
  !!
  !! \param h_bml Input Hamiltonian matrix
  !! \param rho_bml Output density matrix
  !! \param threshold Threshold for sparse matrix algebra
  !! \param bndfil Bond 
  !! \param minsp2iter Minimum sp2 iterations
  !! \param maxsp2iter Maximum SP2 iterations
  !! \param sp2conv Convergence type
  !! \param idemtol Idempotency tolerance
  !! \param verbose A verbosity level
  subroutine sp2_basic(h_bml,rho_bml,threshold,bndfil,minsp2iter,maxsp2iter &
      ,sp2conv,idemtol,verbose)

    implicit none
    integer :: hdim,iter,minsp2iter,maxsp2iter
    integer :: verbose
    real(dp) :: trx, occ, trx2, limdiff
    real(dp) :: gershfact, bndfil, maxeval, maxminusmin
    real(dp) :: threshold, mls_i, mls_ii
    real(dp) :: idemtol
    real(dp), allocatable :: trace(:)
    character(len=*), intent(in) :: sp2conv
    type(bml_matrix_t), intent(inout) :: rho_bml
    type(bml_matrix_t), intent(in) :: h_bml
    type(bml_matrix_t) :: x2_bml


    hdim = bml_get_n(h_bml)  !to get it from the bml matrix

    ! we're also using niklasson's scheme to determine convergence    
    occ = bndfil*real(hdim,dp)

    ! Normalize
    call bml_copy(h_bml, rho_bml)
    call normalize(rho_bml)

    allocate(trace(2))

    ! X2 <- X
    call bml_copy_new(rho_bml, x2_bml)

    do iter = 0, maxsp2iter

      trx = bml_trace(rho_bml)

      ! X2 <- X * X
      call bml_multiply_x2(rho_bml, x2_bml, threshold, trace)

      trx2 = trace(2)

#ifdef DO_MPI_BLOCK
      !< Trace reduction
      if (getNRanks() > 1) then
          call sumRealReduce2(trx, trx2)
      endif
#endif

      if(verbose.gt.1 .and. printRank() .eq. 1) write(*,*)"sp2iter", iter, occ, trx, abs(occ-trx)

      if(trx - occ .le. 0.0_dp) then

        ! X <- 2 * X - X2
        call bml_add_deprecated(2.00_dp, rho_bml, -1.00_dp, x2_bml, threshold)

        trx = 2.0_dp * trx - trx2

      else

        ! X <- X2
        call bml_copy(x2_bml,rho_bml) !x2 = d          

        trx = trx2

      end if

      if(abs(occ-trx) .lt. idemtol .and. iter .gt. minsp2iter) exit

      if(iter .eq. maxsp2iter) then
        write(*,*) "sp2 purification is not converging: stop!"
        error stop
      end if

#ifdef DO_MPI_BLOCK
        !< Send new matrix to all ranks
        if (getNRanks() > 1) then
          !call allGatherParallel(rho_bml)
        endif
#endif

  end do

  call bml_scale(2.0d0, rho_bml) !d = 2*d

  call bml_deallocate(x2_bml)
  deallocate(trace)

end subroutine sp2_basic

  !! Calculate the density matrix from a Hamiltonian matrix by
  !! purification. The method implemented here is the SP2 method
  !! [A. Niklasson, Phys. Rev. B, 66, 155115 (2002)].
  !!
  !! \param h_bml Input Hamiltonian matrix
  !! \param rho_bml Output density matrix
  !! \param threshold Threshold for sparse matrix algebra
  !! \param bndfil Bond 
  !! \param minsp2iter Minimum sp2 iterations
  !! \param maxsp2iter Maximum SP2 iterations
  !! \param sp2conv Convergence type
  !! \param idemtol Idempotency tolerance
  !! \param verbose Verbosity level
subroutine sp2_alg2(h_bml, rho_bml, threshold, bndfil, &
                      minsp2iter, maxsp2iter, sp2conv, idemtol,verbose)

    integer, intent(in) :: minsp2iter, maxsp2iter
    real(dp), intent(in) :: threshold, bndfil, idemtol
    character(len=*), intent(in) :: sp2conv
    type(bml_matrix_t),intent(in) :: h_bml
    type(bml_matrix_t),intent(inout) :: rho_bml

    integer, optional, intent(in) :: verbose
    integer :: HDIM, iter, breakloop
    real(dp) :: trx, trx2, trxOld, tr2xx2
    real(dp) :: occ, limdiff
    real(dp) :: idemperr, idemperr1, idemperr2
    real(dp), allocatable :: trace(:)

    type(bml_matrix_t) :: x2_bml

    idemperr = 0.0_dp
    idemperr1 = 0.0_dp
    idemperr2 = 0.0_dp

    HDIM = bml_get_N(h_bml)
    occ = bndfil*FLOAT(HDIM)
    if (printRank() .eq. 1) write(*,*) 'OCC = ', occ

    !! Normalize
    call bml_copy(h_bml, rho_bml)
    call normalize(rho_bml)

    allocate(trace(2))

    trx = bml_trace(rho_bml)

    iter = 0
    breakloop = 0

    ! X2 <- X
    call bml_copy_new(rho_bml, x2_bml)

    do while (breakloop .eq. 0 .and. iter .lt. maxsp2iter)

      iter = iter + 1

      ! X2 <- X * X
      call bml_multiply_x2(rho_bml, x2_bml, threshold, trace)

      trx2 = trace(2)

#ifdef DO_MPI_BLOCK
      !< Trace reduction
      if (getNRanks() > 1) then
          call sumRealReduce2(trx, trx2)
      endif
#endif
      
      if(present(verbose))then
        if(verbose.GE.1) then 
          if (printRank() .eq. 1) write(*,*) 'iter = ', iter, 'trx = ', trx, ' trx2 = ', trx2
        endif
      endif  

      tr2xx2 = 2.0_dp*trx - trx2
      trXOld = trx
      limDiff = abs(trx2 - occ) - abs(tr2xx2 - occ)

      if (limdiff .ge. idemtol) then

        ! X <- 2 * X - X2
        call bml_add_deprecated(2.0_dp, rho_bml, -1.0_dp, x2_bml, threshold)

        trx = 2.0_dp * trx - trx2

      elseif(limdiff .lt. -idemtol) then

        ! X <- X2
        call bml_copy(x2_bml, rho_bml)

        trx = trx2

      else

        iter = iter - 1
        trx = trxOld
        breakloop = 1

      end if

      idemperr2 = idemperr1
      idemperr1 = idemperr
      idemperr = abs(trx - trxOld)

      if (sp2conv .eq. "Rel" .and. iter .ge. minsp2iter .and. &
         (idemperr .ge. idemperr2 .or. idemperr .lt. idemtol)) then
         breakloop = 1
      end if

      if (iter .eq. maxsp2iter) then
         write(*,*) "SP2 purification is not converging: STOP!"
         stop
      end if

#ifdef DO_MPI_BLOCK
        !< Send new matrix to all ranks
        if (getNRanks() > 1) then
          call allGatherParallel(rho_bml)
        endif
#endif

     end do

     ! X <- 2 * X
     call bml_scale(2.0_dp, rho_bml)       !D = 2.0_dp*D

     call bml_deallocate(x2_bml)
     deallocate(trace)

  end subroutine sp2_alg2

  !! Perform SP2 algorithm, generate sequence, and calculate norm.
  !!
  !! \param h_bml Input Hamiltonian matrix
  !! \param rho_bml Output density matrix
  !! \param threshold Threshold for sparse matrix algebra
  !! \param bndfil Bond 
  !! \param minsp2iter Minimum sp2 iterations
  !! \param maxsp2iter Maximum SP2 iterations
  !! \param sp2conv Convergence type
  !! \param idemtol Idempotency tolerance
  !! \param pp Vector containing sequence of 0s and 1s
  !! \param icount Sequence count
  !! \param vv Vector for sum of squares per iteration
  subroutine sp2_alg2_genseq(h_bml, rho_bml, threshold, bndfil, minsp2iter, &
                             maxsp2iter, sp2conv, idemtol, pp, icount, vv, verbose)
    implicit none
    integer, intent(in) :: minsp2iter, maxsp2iter
    integer, intent(inout) :: icount
    integer, intent(inout) :: pp(:)
    real(dp), intent(in) :: threshold, bndfil, idemtol
    real(dp), intent(inout) :: vv(:)
    character(len=*), intent(in) :: sp2conv
    type(bml_matrix_t),intent(in) :: h_bml
    type(bml_matrix_t),intent(inout) :: rho_bml

    integer, optional, intent(in) :: verbose
    integer :: HDIM, iter, breakloop, myRank
    real(dp) :: trx, trx2, trxOld, tr2xx2
    real(dp) :: occ, limdiff, ssum
    real(dp) :: idemperr, idemperr1, idemperr2
    real(dp), allocatable :: trace(:)

    type(bml_matrix_t) :: x2_bml

    idemperr = 0.0_dp
    idemperr1 = 0.0_dp
    idemperr2 = 0.0_dp 

    myRank = getMyRank()

    HDIM = bml_get_N(h_bml)
    occ = bndfil*FLOAT(HDIM)
    if (printRank() .eq. 1) write(*,*) 'OCC = ', occ

    !! Normalize
    call bml_copy(h_bml, rho_bml)
    call normalize(rho_bml)

#ifdef DO_MPI
        !< Send new matrix to all ranks
        if (getNRanks() .gt. 1 .and. &
            bml_get_distribution_mode(rho_bml) .eq.  BML_DMODE_DISTRIBUTED) then
          call allGatherParallel(rho_bml)
        endif
#endif

    allocate(trace(2))

    iter = 0
    breakloop = 0

    call bml_copy_new(rho_bml, x2_bml)

    do while (breakloop .eq. 0 .and. iter .lt. maxsp2iter)

      iter = iter + 1

      ! X2 <- X * X
      call bml_multiply_x2(rho_bml, x2_bml, threshold, trace)
      ssum = bml_sum_squares2(rho_bml, x2_bml, 1.0_dp, -1.0_dp, threshold)
      trx = trace(1)
      trx2 = trace(2)

#ifdef DO_MPI
      !< Trace and trace norm reduction
      if (getNRanks() .gt. 1 .and. & 
          bml_get_distribution_mode(rho_bml) .eq.  BML_DMODE_DISTRIBUTED) then
          call sumRealReduce3(trx, trx2, ssum)
      endif
#endif

      vv(iter) = sqrt(ssum)

      if(present(verbose))then
        if(verbose.GE.1) then 
          if (printRank() .eq. 1) write(*,*) 'iter = ', iter, 'trx = ', trx, ' trx2 = ', trx2
        endif
      endif
        
      tr2xx2 = 2.0_dp * trx - trx2
      trXOld = trx
      limDiff = abs(trx2 - occ) - abs(tr2xx2 - occ)

      if (limdiff .ge. idemtol) then

        ! X <- 2 * X - X2
        call bml_add_deprecated(2.0_dp, rho_bml, -1.0_dp, x2_bml, threshold)

        trx = 2.0_dp * trx - trx2

        pp(iter) = 0

      elseif(limdiff .lt. -idemtol) then

        ! X <- X2
        call bml_copy(x2_bml, rho_bml)

        trx = trx2

        pp(iter) = 1

      else

        iter = iter - 1
        trx = trxOld
        breakloop = 1

      end if

      idemperr2 = idemperr1
      idemperr1 = idemperr
      idemperr = abs(trx - trxOld)

      if (sp2conv .eq. "Rel" .and. iter .ge. minsp2iter .and. &
         (idemperr .ge. idemperr2 .or. idemperr .lt. idemtol)) then
         breakloop = 1
      end if

      if (iter .eq. maxsp2iter) then
         write(*,*) "SP2 purification is not converging: STOP!"
         stop
      end if

#ifdef DO_MPI
        !< Send new matrix to all ranks
        if (getNRanks() .gt. 1 .and. &
            bml_get_distribution_mode(rho_bml) .eq.  BML_DMODE_DISTRIBUTED) then
          call allGatherParallel(rho_bml)
        endif
#endif

     end do

     icount = iter

     ! X = 2 * X
     call bml_scale(2.0_dp, rho_bml)       !D = 2.0_dp*D

     call bml_deallocate(x2_bml)
     deallocate(trace)

  end subroutine sp2_alg2_genseq

  !! Perform SP2 algorithm given a sequence and calculate norm.
  !!
  !! \param h_bml Input Hamiltonian matrix
  !! \param rho_bml Output density matrix
  !! \param threshold Threshold for sparse matrix algebra
  !! \param pp Vector containing sequence of 0s and 1s
  !! \param icount Sequence count
  !! \param vv Vector of sum of squares per iteration
  subroutine sp2_alg2_seq(h_bml, rho_bml, threshold, pp, icount, vv, verbose)

    integer, intent(inout) :: icount
    integer, intent(inout) :: pp(:)
    real(dp), intent(in) :: threshold
    real(dp), intent(inout) :: vv(:)
    type(bml_matrix_t),intent(in) :: h_bml
    type(bml_matrix_t),intent(inout) :: rho_bml

    integer, optional, intent(in) :: verbose
    integer :: iter
    real(dp) :: occ, trx, trx2, ssum
    real(dp), allocatable :: trace(:)

    type(bml_matrix_t) :: x2_bml

    !! Normalize
    call bml_copy(h_bml, rho_bml)
    call normalize(rho_bml)

    allocate(trace(2))

    trx = bml_trace(rho_bml)

    iter = 0

    ! X2 <- X
    call bml_copy_new(rho_bml, x2_bml)

    do while (iter .lt. icount)

      iter = iter + 1

      ! X2 <- X * X
      call bml_multiply_x2(rho_bml, x2_bml, threshold, trace)
      ssum = bml_sum_squares2(rho_bml, x2_bml, 1.0_dp, -1.0_dp, threshold)
      trx2 = trace(2)

#ifdef DO_MPI_BLOCK
      !< Trace and trace norm reduction
      if (getNRanks() > 1) then
          call sumRealReduce3(trx, trx2, ssum)
      endif
#endif

      vv(iter) = sqrt(ssum)

      if(present(verbose))then
        if(verbose.GE.1) then 
          if (printRank() .eq. 1) write(*,*) 'iter = ', iter, 'trx = ', trx, ' trx2 = ', trx2
        endif
      endif
       
      if (pp(iter) .eq. 0) then

        ! X <- 2 * X - X2
        call bml_add_deprecated(2.0_dp, rho_bml, -1.0_dp, x2_bml, threshold)

        trx = 2.0_dp * trx - trx2

      else

        ! X <- X2
        call bml_copy(x2_bml, rho_bml)

        trx = trx2

      end if

#ifdef DO_MPI_BLOCK
        !< Send new matrix to all ranks
        if (getNRanks() > 1) then
          call allGatherParallel(rho_bml)
        endif
#endif

     end do

     ! X <- 2 * X
     call bml_scale(2.0_dp, rho_bml)       !D = 2.0_dp*D

     call bml_deallocate(x2_bml)
     deallocate(trace)
  
  end subroutine sp2_alg2_seq
  
  !! Perform SP2 algorithm given a sequence and calculate norm.
  !!
  !! \param rho_bml Input Hamiltonian/Output density matrix
  !! \param threshold Threshold for sparse matrix algebra
  !! \param pp Vector containing sequence of 0s and 1s
  !! \param icount Sequence count
  !! \param vv Vector of sum of squares per iteration
  !! \param mineval Min value used for normalization (optional)
  !! \param maxeval Max value used for normalization (optional)
  subroutine sp2_alg2_seq_inplace(rho_bml, threshold, pp, icount, vv, &
       mineval, maxeval, verbose)

    integer, intent(inout) :: icount
    integer, intent(inout) :: pp(:)
    real(dp), intent(in) :: threshold
    real(dp), intent(inout) :: vv(:)
    real(dp), optional, intent(in) :: mineval, maxeval
    type(bml_matrix_t),intent(inout) :: rho_bml

    integer, optional, intent(in) :: verbose
    integer :: iter
    real(dp) :: occ, trx, trx2, ssum
    real(dp), allocatable :: trace(:)

    type(bml_matrix_t) :: x2_bml

    !! Normalize
    call bml_normalize(rho_bml, mineval, maxeval)

    allocate(trace(2))

    trx = bml_trace(rho_bml)

    iter = 0

    ! X2 <- X
    call bml_copy_new(rho_bml, x2_bml)

    do while (iter .lt. icount)

      iter = iter + 1

      ! X2 <- X * X
      call bml_multiply_x2(rho_bml, x2_bml, threshold, trace)
      ssum = bml_sum_squares2(rho_bml, x2_bml, 1.0_dp, -1.0_dp, threshold)
      trx2 = trace(2)

#ifdef DO_MPI_BLOCK
      !< Trace and trace norm reduction
      if (getNRanks() > 1) then
          call sumRealReduce3(trx, trx2, ssum)
      endif
#endif

      vv(iter) = sqrt(ssum)

      if(present(verbose))then
        if(verbose.GE.1) then 
          if (printRank() .eq. 1) write(*,*) 'iter = ', iter, 'trx = ', trx, ' trx2 = ', trx2
        endif
      endif
        
      if (pp(iter) .eq. 0) then

        ! X <- 2 * X - X2
        call bml_add_deprecated(2.0_dp, rho_bml, -1.0_dp, x2_bml, threshold)

        trx = 2.0_dp * trx - trx2

      else

        ! X <- X2
        call bml_copy(x2_bml, rho_bml)

        trx = trx2

      end if

#ifdef DO_MPI_BLOCK
        !< Send new matrix to all ranks
        if (getNRanks() > 1) then
          call allGatherParallel(rho_bml)
        endif
#endif

     end do

     ! X <- 2 * X
     call bml_scale(2.0D0, rho_bml)       !D = 2.0_dp*D

     call bml_deallocate(x2_bml)
     deallocate(trace)

  end subroutine sp2_alg2_seq_inplace

  !! Perform SP2 algorithm.
  !!
  !! \param h_bml Input Hamiltonian matrix
  !! \param rho_bml Output density matrix
  !! \param threshold Threshold for sparse matrix algebra
  !! \param bndfil Bond 
  !! \param minsp2iter Minimum sp2 iterations
  !! \param maxsp2iter Maximum SP2 iterations
  !! \param sp2conv Convergence type
  !! \param idemtol Idempotency tolerance
  subroutine sp2_alg1(h_bml, rho_bml, threshold, bndfil, minsp2iter, maxsp2iter, &
                          sp2conv, idemtol, verbose)

    integer, intent(in) :: minsp2iter, maxsp2iter
    integer, optional, intent(in) :: verbose
    real(dp), intent(in) :: threshold, bndfil, idemtol
    character(len=*), intent(in) :: sp2conv
    type(bml_matrix_t),intent(in) :: h_bml
    type(bml_matrix_t),intent(inout) :: rho_bml

    integer :: HDIM, iter, breakloop
    real(dp) :: trx, trx2, trxOld, tr2xx2
    real(dp) :: occ, limdiff
    real(dp) :: idemperr, idemperr1, idemperr2

    type(bml_matrix_t) :: x2_bml

    idemperr = 0.0_dp
    idemperr1 = 0.0_dp
    idemperr2 = 0.0_dp

    ! We're also using Niklasson's scheme to determine convergence

    HDIM = bml_get_N(h_bml)
    occ = bndfil*FLOAT(HDIM)
    if (printRank() .eq. 1) write(*,*) 'OCC = ', occ

    !! Normalize
    call bml_copy(h_bml, rho_bml)
    call normalize(rho_bml)

    trx = bml_trace(rho_bml)

    iter = 0
    breakloop = 0

    call bml_copy_new(rho_bml, x2_bml)

    do while (breakloop .eq. 0 .and. iter .lt. maxsp2iter)

      iter = iter + 1

      ! X2 <- X
      ! X2 <- X - X * X
      call bml_copy(rho_bml, x2_bml)
      call bml_multiply(rho_bml, rho_bml, x2_bml, -1.0_dp, 1.0_dp, threshold)

      trx2 = bml_trace(x2_bml)

#ifdef DO_MPI_BLOCK
      !< Trace and trace norm reduction
      if (getNRanks() > 1) then
          call sumRealReduce2(trx, trx2)
      endif
#endif

      if(present(verbose).and.verbose.ge.10)then 
        if (printRank() .eq. 1) write(*,*) 'iter = ', iter, 'trx = ', trx, ' trx2 = ', trx2
      endif  

      limdiff = abs(trx - trx2 - occ) - abs(trx + trx2 - occ)

      if (limdiff .ge. idemtol) then

        ! X <- X + (X - X * X) <- 2 * X - X * X
        call bml_add_deprecated(1.0_dp, rho_bml, 1.0_dp, x2_bml, threshold)

        trx = trx + trx2

      elseif(limdiff .lt. -idemtol) then

        ! X <- X - (X - X * X) <- X * X
        call bml_add_deprecated(1.0_dp, rho_bml, -1.0_dp, x2_bml, threshold)

        trx = trx - trx2

      else

        iter = iter - 1
        breakloop = 1

      end if

      idemperr2 = idemperr1
      idemperr1 = idemperr
      idemperr = abs(trx2)

      if (sp2conv .eq. "Rel" .and. iter .ge. minsp2iter .and. &
         (idemperr2 .le. idemperr .or. idemperr .lt. idemtol)) then
         breakloop = 1
      end if

      if (iter .eq. maxsp2iter) then
         write(*,*) "SP2 purification is not converging: STOP!"
         stop
      end if

#ifdef DO_MPI_BLOCK
        !< Send new matrix to all ranks
        if (getNRanks() > 1) then
          call allGatherParallel(rho_bml)
        endif
#endif

     end do

     ! X <- 2 * X
     call bml_scale(2.0_dp, rho_bml)       !D = 2.0_dp*D

     call bml_deallocate(x2_bml)

  end subroutine sp2_alg1

  !! Perform SP2 algorithm, generate sequence, and calculate norm.
  !!
  !! \param h_bml Input Hamiltonian matrix
  !! \param rho_bml Output density matrix
  !! \param threshold Threshold for sparse matrix algebra
  !! \param bndfil Bond 
  !! \param minsp2iter Minimum sp2 iterations
  !! \param maxsp2iter Maximum SP2 iterations
  !! \param sp2conv Convergence type
  !! \param idemtol Idempotency tolerance
  !! \param pp Vector containing sequence of 0s and 1s
  !! \param icount Sequence count
  !! \param vv Vector of sqrt of TrNorm per iteration
  subroutine sp2_alg1_genseq(h_bml, rho_bml, threshold, bndfil, &
                          minsp2iter, maxsp2iter, sp2conv, idemtol, &
                          pp, icount, vv)

    integer, intent(in) :: minsp2iter, maxsp2iter
    integer, intent(inout) :: icount
    integer, intent(inout) :: pp(:)
    real(dp), intent(in) :: threshold, bndfil, idemtol
    real(dp), intent(inout) :: vv(:)
    character(len=*), intent(in) :: sp2conv
    type(bml_matrix_t),intent(in) :: h_bml
    type(bml_matrix_t),intent(inout) :: rho_bml

    integer :: HDIM, iter, breakloop
    real(dp) :: trx, trx2, trxOld, tr2xx2
    real(dp) :: occ, limdiff, ssum
    real(dp) :: idemperr, idemperr1, idemperr2
    type(bml_matrix_t) :: x2_bml

    idemperr = 0.0_dp
    idemperr1 = 0.0_dp
    idemperr2 = 0.0_dp

    ! We're also using Niklasson's scheme to determine convergence

    HDIM = bml_get_N(h_bml)
    occ = bndfil*FLOAT(HDIM)
    if (printRank() .eq. 1) write(*,*) 'OCC = ', occ

    !! Normalize
    call bml_copy(h_bml, rho_bml)
    call normalize(rho_bml)

    trx = bml_trace(rho_bml)

    iter = 0
    breakloop = 0

    call bml_copy_new(rho_bml, x2_bml)

    do while (breakloop .eq. 0 .and. iter .lt. maxsp2iter)

      iter = iter + 1

      ! X2 <- X
      ! X2 <- X - X * X
      call bml_copy(rho_bml, x2_bml)
      call bml_multiply(rho_bml, rho_bml, x2_bml, -1.0_dp, 1.0_dp, threshold)
      ssum = bml_sum_squares(x2_bml)
      trx2 = bml_trace(x2_bml)

#ifdef DO_MPI_BLOCK
      !< Trace and trace norm reduction
      if (getNRanks() > 1) then
          call sumRealReduce3(trx, trx2, ssum)
      endif
#endif

      vv(iter) = sqrt(ssum)

      if (printRank() .eq. 1) write(*,*) 'iter = ', iter, 'trx = ', trx, ' trx2 = ', trx2

      limdiff = abs(trx - trx2 - occ) - abs(trx + trx2 - occ)

      if (limdiff .ge. idemtol) then

        ! X <- X + (X - X * X) <- 2 * X - X * X
        call bml_add_deprecated(1.0_dp, rho_bml, 1.0_dp, x2_bml, threshold)

        trx = trx + trx2

        pp(iter) = 0

      elseif(limdiff .lt. -idemtol) then

        ! X <- X - (X - X * X) <- X * X
        call bml_add_deprecated(1.0_dp, rho_bml, -1.0_dp, x2_bml, threshold)

        trx = trx - trx2

        pp(iter) = 1

      else

        iter = iter - 1
        breakloop = 1

      end if

      idemperr2 = idemperr1
      idemperr1 = idemperr
      idemperr = abs(trx2)

      if (sp2conv .eq. "Rel" .and. iter .ge. minsp2iter .and. &
         (idemperr2 .le. idemperr .or. idemperr .lt. idemtol)) then
         breakloop = 1
      end if

      if (sp2conv .eq. "Abs" .and. abs(limdiff) .lt. idemtol) exit

      if (iter .eq. maxsp2iter) then
         write(*,*) "SP2 purification is not converging: STOP!"
         stop
      end if

#ifdef DO_MPI_BLOCK
        !< Send new matrix to all ranks
        if (getNRanks() > 1) then
          call allGatherParallel(rho_bml)
        endif
#endif

     end do

     icount = iter

     ! X <- 2 * X
     call bml_scale(2.0_dp, rho_bml)       !D = 2.0_dp*D

     call bml_deallocate(x2_bml)

  end subroutine sp2_alg1_genseq

  !! Perform SP2 algorithm using sequence and calculate norm.
  !!
  !! \param h_bml Input Hamiltonian matrix
  !! \param rho_bml Output density matrix
  !! \param threshold Threshold for sparse matrix algebra
  !! \param pp Vector containing sequence of 0s and 1s
  !! \param icount Sequence count
  !! \param Vector of sum of squares per iteration
  subroutine sp2_alg1_seq(h_bml, rho_bml, threshold, pp, icount, vv)

    integer, intent(inout) :: icount
    integer, intent(inout) :: pp(:)
    real(dp), intent(in) :: threshold
    real(dp), intent(inout) :: vv(:)
    type(bml_matrix_t),intent(in) :: h_bml
    type(bml_matrix_t),intent(inout) :: rho_bml

    integer :: iter, breakloop
    real(dp) :: trx, trx2, ssum
    type(bml_matrix_t) :: x2_bml

    ! We're also using Niklasson's scheme to determine convergence

    !! Normalize
    call bml_copy(h_bml, rho_bml)
    call normalize(rho_bml)

    trx = bml_trace(rho_bml)

    iter = 0
    breakloop = 0

    call bml_copy_new(rho_bml, x2_bml)

    do while (iter .lt. icount)

      iter = iter + 1

      ! X2 <- X
      ! X2 <- X - X * X
      call bml_copy(rho_bml, x2_bml)
      call bml_multiply(rho_bml, rho_bml, x2_bml, -1.0_dp, 1.0_dp, threshold)
      ssum = bml_sum_squares(x2_bml)

      trx2 = bml_trace(x2_bml)

#ifdef DO_MPI_BLOCK
      !< Trace and trace norm reduction
      if (getNRanks() > 1) then
          call sumRealReduce3(trx, trx2, ssum)
      endif
#endif

      vv(iter) = sqrt(ssum)

      if (printRank() .eq. 1) write(*,*) 'iter = ', iter, 'trx = ', trx, ' trx2 = ', trx2

      if (pp(iter) .eq. 0) then

        ! X <- X + (X - X * X) <- 2 * X - X * X
        call bml_add_deprecated(1.0_dp, rho_bml, 1.0_dp, x2_bml, threshold)

        trx = trx + trx2

      else

        ! X <- X - (X - X * X) <- X * X
        call bml_add_deprecated(1.0_dp, rho_bml, -1.0_dp, x2_bml, threshold)

        trx = trx - trx2

      end if

#ifdef DO_MPI_BLOCK
        !< Send new matrix to all ranks
        if (getNRanks() > 1) then
          call allGatherParallel(rho_bml)
        endif
#endif

     end do

     ! X <- 2 * X
     call bml_scale(2.0_dp, rho_bml)       !D = 2.0_dp*D

     call bml_deallocate(x2_bml)

  end subroutine sp2_alg1_seq

  !! Perform SP2 algorithm using sequence and calculate norm.
  !!
  !! \param rho_bml Input Hamiltonian/Output density matrix
  !! \param threshold Threshold for sparse matrix algebra
  !! \param pp Vector containing sequence of 0s and 1s
  !! \param icount Sequence count
  !! \param Vector of sum of squares per iteration
  !! \param mineval Min value used for normalization (optional)
  !! \param maxeval Max value used for normalization (optional)
  subroutine sp2_alg1_seq_inplace(rho_bml, threshold, pp, icount, vv, &
      mineval, maxeval)

    integer, intent(inout) :: icount
    integer, intent(inout) :: pp(:)
    real(dp), intent(in) :: threshold
    real(dp), intent(inout) :: vv(:)
    real(dp), intent(in) :: mineval, maxeval
    type(bml_matrix_t),intent(inout) :: rho_bml

    integer :: iter, breakloop
    real(dp) :: trx, trx2, ssum
    type(bml_matrix_t) :: x2_bml

    ! We're also using Niklasson's scheme to determine convergence

    !! Normalize
    call bml_normalize(rho_bml, mineval, maxeval)

    trx = bml_trace(rho_bml)

    iter = 0
    breakloop = 0

    call bml_copy_new(rho_bml, x2_bml)

    do while (iter .lt. icount)

      iter = iter + 1

      ! X2 <- X
      ! X2 <- X - X * X
      call bml_copy(rho_bml, x2_bml)
      call bml_multiply(rho_bml, rho_bml, x2_bml, -1.0_dp, 1.0_dp, threshold)
      ssum = bml_sum_squares(x2_bml)

      trx2 = bml_trace(x2_bml)

#ifdef DO_MPI_BLOCK
      !< Trace and trace norm reduction
      if (getNRanks() > 1) then
          call sumRealReduce3(trx, trx2, ssum)
      endif
#endif

      vv(iter) = sqrt(ssum)

      if (pp(iter) .eq. 0) then

        ! X <- X + (X - X * X) <- 2 * X - X * X
        call bml_add_deprecated(1.0_dp, rho_bml, 1.0_dp, x2_bml, threshold)

        trx = trx + trx2

      else

        ! X <- X - (X - X * X) <- X * X
        call bml_add_deprecated(1.0_dp, rho_bml, -1.0_dp, x2_bml, threshold)

        trx = trx - trx2

      end if

#ifdef DO_MPI_BLOCK
        !< Send new matrix to all ranks
        if (getNRanks() > 1) then
          call allGatherParallel(rho_bml)
        endif
#endif

     end do

     ! X <- 2 * X
     call bml_scale(2.0_dp, rho_bml)       !D = 2.0_dp*D

     call bml_deallocate(x2_bml)

  end subroutine sp2_alg1_seq_inplace

  !> Perform SP2 algorithm using sequence and calculate norm
  !! for a submatrix.
  !!
  !! \param rho_bml Input Hamiltonian/Output density matrix
  !! \param threshold Threshold for sparse matrix algebra
  !! \param pp Vector containing sequence of 0s and 1s
  !! \param icount Sequence count
  !! \param vv Vector of sum of squares per iteration
  !! \param mineval Min value used for normalization (optional)
  !! \param maxeval Max value used for normalization (optional)
  !! \param core_size Number of core rows
  subroutine sp2_submatrix(ham_bml, rho_bml, threshold, pp, icount, vv, &
      mineval, maxeval, core_size)

    integer, intent(in) :: icount
    integer, intent(in) :: pp(:)
    integer, intent(in) :: core_size
    real(dp), intent(in) :: threshold
    real(dp), intent(inout) :: vv(:)
    real(dp), intent(in) :: mineval, maxeval
    type(bml_matrix_t),intent(in) :: ham_bml
    type(bml_matrix_t),intent(inout) :: rho_bml

    integer :: iter
    real(dp) :: trx, trx2, factor, ssum
    type(bml_matrix_t) :: x2_bml

    ! We're also using Niklasson's scheme to determine convergence

    !! Normalize
    call bml_copy(ham_bml, rho_bml)
    call bml_normalize(rho_bml, mineval, maxeval)

    trx = bml_trace(rho_bml)

    call bml_copy_new(rho_bml, x2_bml)

    do iter = 1, icount

      ! X2 <- X
      ! X2 <- X - X * X
      call bml_copy(rho_bml, x2_bml)

      call bml_multiply(rho_bml, rho_bml, x2_bml, -1.0_dp, 1.0_dp, threshold)

      vv(iter) = vv(iter) + bml_sum_squares_submatrix(x2_bml, core_size)

      trx2 = bml_trace(x2_bml)

!      if (printRank() .eq. 1) write(*,*) 'iter = ', iter, &
!        'trx = ', trx, ' trx2 = ', trx2, ' vv = ', vv(iter)

      factor = 1.0_dp - 2.0_dp * pp(iter)

      ! X <- X + (X - X * X) <- 2 * X - X * X  when pp(iter) == 0
      ! X <- X - (X - X * X) <- X * X          when pp(iter) == 1
      call bml_add_deprecated(1.0_dp, rho_bml, factor, x2_bml, threshold)

      trx = trx + factor * trx2

     end do

     ! X <- 2 * X
     call bml_scale(2.0_dp, rho_bml)       !D = 2.0_dp*D

     call bml_deallocate(x2_bml)

  end subroutine sp2_submatrix

  !! Perform SP2 algorithm using sequence and calculate norm
  !! for a submatrix.
  !!
  !! \param rho_bml Input Hamiltonian/Output density matrix
  !! \param threshold Threshold for sparse matrix algebra
  !! \param pp Vector containing sequence of 0s and 1s
  !! \param icount Sequence count
  !! \param vv Vector of sum of squares per iteration
  !! \param mineval Min value used for normalization (optional)
  !! \param maxeval Max value used for normalization (optional)
  !! \param core_size Number of core rows
  subroutine sp2_submatrix_inplace(rho_bml, threshold, pp, icount, vv, &
      mineval, maxeval, core_size)

    integer, intent(inout) :: icount
    integer, intent(inout) :: pp(:)
    integer, intent(in) :: core_size
    real(dp), intent(in) :: threshold
    real(dp), intent(inout) :: vv(:)
    real(dp), intent(in) :: mineval, maxeval
    type(bml_matrix_t),intent(inout) :: rho_bml

    integer :: iter
    real(dp) :: trx, trx2, factor, ssum
    type(bml_matrix_t) :: x2_bml

    ! We're also using Niklasson's scheme to determine convergence

    !! Normalize
    call bml_normalize(rho_bml, mineval, maxeval)

    trx = bml_trace(rho_bml)

    call bml_copy_new(rho_bml, x2_bml)

    do iter = 1, icount

      ! X2 <- X
      ! X2 <- X - X * X
      call bml_copy(rho_bml, x2_bml)

      call bml_multiply(rho_bml, rho_bml, x2_bml, -1.0_dp, 1.0_dp, threshold)

      vv(iter) = vv(iter) + bml_sum_squares_submatrix(x2_bml, core_size)

      trx2 = bml_trace(x2_bml)

      factor = 1.0_dp - 2.0_dp * pp(iter)

      ! X <- X + (X - X * X) <- 2 * X - X * X  when pp(iter) == 0
      ! X <- X - (X - X * X) <- X * X          when pp(iter) == 1
      call bml_add_deprecated(1.0_dp, rho_bml, factor, x2_bml, threshold)

      trx = trx + factor * trx2

     end do

     ! X <- 2 * X
     call bml_scale(2.0_dp, rho_bml)       !D = 2.0_dp*D

     call bml_deallocate(x2_bml)

  end subroutine sp2_submatrix_inplace

end module sp2_mod
