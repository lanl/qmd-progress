!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 2010.  Los Alamos National Security, LLC. This material was    !
! produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos !
! National Laboratory (LANL), which is operated by Los Alamos National     !
! Security, LLC for the U.S. Department of Energy. The U.S. Government has !
! rights to use, reproduce, and distribute this software.  NEITHER THE     !
! GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY,     !
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS         !
! SOFTWARE.  If software is modified to produce derivative works, such     !
! modified software should be clearly marked, so as not to confuse it      !
! with the version available from LANL.                                    !
!                                                                          !
! Additionally, this program is free software; you can redistribute it     !
! and/or modify it under the terms of the GNU General Public License as    !
! published by the Free Software Foundation; version 2.0 of the License.   !
! Accordingly, this program is distributed in the hope that it will be     !
! useful, but WITHOUT ANY WARRANTY; without even the implied warranty of   !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General !
! Public License for more details.                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> The parallel module.
!! \ingroup PROGRESS
    !
    !

module parallel_mod

  use bml

#ifdef DO_MPI
  use MPI
#endif

#ifdef DO_MPI
#ifdef SINGLE
#define REAL_MPI_TYPE MPI_FLOAT
#else
#define REAL_MPI_TYPE MPI_DOUBLE
#endif
#endif

  implicit none

  private  !Everything is private by default

  integer, parameter :: dp = kind(1.0d0)

  integer :: myRank, nRanks
  integer :: ierr, reqCount
  integer, allocatable :: requestList(:), rUsed(:)

  public :: rankReduceData_t
  public :: getNRanks
  public :: getMyRank
  public :: printRank
  public :: initParallel 
  public :: shutdownParallel
  public :: barrierParallel
  public :: sendReceiveParallel
  public :: isendParallel
  public :: sendParallel
  public :: irecvParallel
  public :: recvParallel
  public :: sumIntParallel
  public :: sumRealParallel
  public :: maxIntParallel
  public :: maxRealParallel
  public :: minIntParallel
  public :: minRealParallel
  public :: minRealReduce
  public :: maxRealReduce
  public :: maxIntReduce2
  public :: sumIntReduce2
  public :: sumIntReduceN
  public :: sumRealReduce
  public :: sumRealReduce2
  public :: sumRealReduce3
  public :: sumRealReduceN
  public :: minRankRealParallel
  public :: maxRankRealParallel
  public :: bcastParallel
  public :: allGatherRealParallel
  public :: allGatherIntParallel
  public :: allGatherVRealParallel
  public :: allGatherVIntParallel
  public :: allSumRealReduceParallel
  public :: allSumIntReduceParallel
  public :: allGatherParallel
  public :: wait
  
  !> Data structure for rection over MPI ranks
  type rankReduceData_t

    !> Data value
    real(dp) :: val

    !> MPI rank
    integer :: rank

  end type rankReduceData_t

  contains

  !
  ! Return total number of ranks
  !
  function getNRanks() result(nR)
    integer :: nR

    nR = nRanks
    return

  end function getNRanks

  !
  ! Return the local rank
  !
  function getMyRank() result(mR)
    integer :: mR

    mR = myRank

    return

  end function getMyRank

  !
  ! Check if this is rank 0 for printing
  !
  function printRank() result(pR)
    integer :: pR

    pR = 0

    if (myRank .eq. 0) then
      pR = 1
    endif

    return

  end function printRank

  !
  ! Initialize MPI
  !
  subroutine initParallel()

#ifdef DO_MPI
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, myRank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nRanks, ierr)

    call bml_initF(MPI_COMM_WORLD)

!    if (printRank() .eq. 1) then
        write(*,*) "MPI started in progress, rank ", myRank, " out of ", nRanks, " ranks"
!    endif

#else
    myRank = 0
    nRanks = 1
#endif

    allocate(requestList(nRanks))
    allocate(rUsed(nRanks))
    rUsed = 0

  end subroutine initParallel

  !
  ! Shutdown MPI
  !
  subroutine shutdownParallel()

    deallocate(requestList)
    deallocate(rUsed)
  
#ifdef DO_MPI
    call bml_shutdownF()

    call MPI_Finalize(ierr)
#endif

  end subroutine shutdownParallel

  !
  ! Save request
  !
  function saveRequest(irequest)

    integer :: saveRequest
    integer,intent(in) :: irequest
    integer :: i

#ifdef DO_MPI
    do i = 1, nRanks
      if (rUsed(i) == 0) then
        requestList(i) = irequest
        rUsed(i) = 1

        saveRequest = i
        return
      endif
    enddo
#endif

    saveRequest = -1
    return

  end function saveRequest

  !
  ! Barrier - all ranks synchronized
  !
  subroutine barrierParallel()

#ifdef DO_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

  end subroutine barrierParallel

  !
  ! Coordinated send and receive
  !
  subroutine sendReceiveParallel(sendBuf, sendLen, dest, recvBuf, recvLen, source, nreceived)

    real(dp), intent(in) :: sendBuf(*)   
    real(dp), intent(out) :: recvBuf(*)    
    integer, intent(in) :: sendLen, recvLen, dest, source
    integer, intent(out) :: nReceived
#ifdef DO_MPI
    integer :: mstatus(MPI_STATUS_SIZE)

    call MPI_Sendrecv(sendBuf, sendLen, REAL_MPI_TYPE, dest, 0, &
                      recvBuf, recvLen, REAL_MPI_TYPE, source, 0, &
                      MPI_COMM_WORLD, mstatus, ierr)
    call MPI_Get_count(mstatus, REAL_MPI_TYPE, nreceived, ierr)
#else
    recvBuf(1:sendLen) = sendBuf(1:sendLen)
    nreceived = sendLen
#endif

  end subroutine sendReceiveParallel

  !
  ! Non-blocking send
  !
  subroutine isendParallel(sendBuf, sendLen, dest)

    real(dp),intent(in) :: sendBuf(*)
    integer,intent(in) :: sendLen, dest
    integer :: request

#ifdef DO_MPI
    call MPI_ISend(sendBuf, sendLen, REAL_MPI_TYPE, dest, &
                         0, MPI_COMM_WORLD, request, ierr)
#endif

  end subroutine isendParallel

  !
  ! Blocking send
  !
  subroutine sendParallel(sendBuf, sendLen, dest)

    real(dp),intent(in) :: sendBuf(*)
    integer,intent(in) :: sendLen, dest

#ifdef DO_MPI
    call MPI_Send(sendBuf, sendLen, REAL_MPI_TYPE, dest, &
                        0, MPI_COMM_WORLD, ierr) 
#endif

  end subroutine sendParallel

  !
  ! Non-blocking receive from any other rank
  !
  subroutine irecvParallel(recvBuf, recvLen, rind)

    real(dp) :: recvBuf(*)
    integer,intent(in) :: recvLen
    integer :: rind
#ifdef DO_MPI
    integer :: request 

    call MPI_IRecv(recvBuf, recvLen, REAL_MPI_TYPE, MPI_ANY_SOURCE, &
                         0, MPI_COMM_WORLD, request, ierr)
    rind = saveRequest(request)
#endif

  end subroutine irecvParallel

  !
  ! Blocking receive from any other rank
  !
  subroutine recvParallel(recvBuf, recvLen)

    real(dp) :: recvBuf(*)
    integer,intent(in) :: recvLen
#ifdef DO_MPI
    integer :: mstatus(MPI_STATUS_SIZE)

    call MPI_Recv(recvBuf, recvLen, REAL_MPI_TYPE, MPI_ANY_SOURCE, &
                         0, MPI_COMM_WORLD, mstatus, ierr)
#endif

  end subroutine recvParallel

  !
  ! Integer sum reduce
  !
  subroutine sumIntParallel(sendBuf, recvBuf, icount)

    integer,intent(in) :: sendBuf(*)
    integer :: recvBuf(*)
    integer,intent(in) :: icount
    integer :: i

#ifdef DO_MPI
    call MPI_AllReduce(sendBuf, recvBuf, icount, MPI_INT, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
#else
    do i = 1, icount
      recvBuf(i) = sendBuf(i)
    enddo
#endif

  end subroutine sumIntParallel

  !
  ! Real sum reduce
  !
  subroutine sumRealParallel(sendBuf, recvBuf, icount)

    real(dp),intent(in) :: sendBuf(*)
    real(dp),intent(out) :: recvBuf(*)
    integer,intent(in) :: icount
    integer :: i

#ifdef DO_MPI
    call MPI_AllReduce(sendBuf, recvBuf, icount, REAL_MPI_TYPE, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)                    
#else
    do i = 1, icount
      recvBuf(i) = sendBuf(i)
    enddo
#endif

  end subroutine sumRealParallel

  !
  ! Integer max reduce
  !
  subroutine maxIntParallel(sendBuf, recvBuf, icount)

    integer,intent(in) :: sendBuf(*)
    integer,intent(out) :: recvBuf(*)
    integer,intent(in) :: icount
    integer :: i

#ifdef DO_MPI
    call MPI_AllReduce(sendBuf, recvBuf, icount, MPI_INT, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)
#else
    do i = 1, icount
      recvBuf(i) = sendBuf(i)
    enddo
#endif

  end subroutine maxIntParallel

  !
  ! Real max reduce
  !
  subroutine maxRealParallel(sendBuf, recvBuf, icount)

    real(dp),intent(in) :: sendBuf(*)
    real(dp),intent(out) :: recvBuf(*)
    integer,intent(in) :: icount
    integer :: i

#ifdef DO_MPI
    call MPI_AllReduce(sendBuf, recvBuf, icount, REAL_MPI_TYPE, &
                       MPI_MAX, MPI_COMM_WORLD, ierr)
#else
    do i = 1, icount
      recvBuf(i) = sendBuf(i)
    enddo
#endif

  end subroutine maxRealParallel

  !
  ! Integer min reduce
  !
  subroutine minIntParallel(sendBuf, recvBuf, icount)

    integer,intent(in) :: sendBuf(*)
    integer,intent(out) :: recvBuf(*)
    integer,intent(in) :: icount
    integer :: i

#ifdef DO_MPI
    call MPI_AllReduce(sendBuf, recvBuf, icount, MPI_INT, &
                       MPI_MIN, MPI_COMM_WORLD, ierr)
#else
    do i = 1, icount
      recvBuf(i) = sendBuf(i)
    enddo
#endif

  end subroutine minIntParallel
                                    
  !
  ! Real min reduce
  !
  subroutine minRealParallel(sendBuf, recvBuf, icount)

    real(dp),intent(in) :: sendBuf(*)
    real(dp),intent(out) :: recvBuf(*)
    integer,intent(in) :: icount
    integer :: i

#ifdef DO_MPI
    call MPI_AllReduce(sendBuf, recvBuf, icount, REAL_MPI_TYPE, &
                       MPI_MIN, MPI_COMM_WORLD, ierr)
#else
    do i = 1, icount
      recvBuf(i) = sendBuf(i)
    enddo
#endif

  end subroutine minRealParallel

  !
  ! Real min reduce for 1 value
  !
  subroutine minRealReduce(rvalue)

    real(dp),intent(inout) :: rvalue
    real(dp) :: sLocal(1),sGlobal(1)

    sLocal(1) = rvalue

    call minRealParallel(sLocal, sGlobal, 1);

    rvalue = sGlobal(1)

  end subroutine minRealReduce

  !
  ! Real max reduce for 1 value
  !
  subroutine maxRealReduce(rvalue)

    real(dp),intent(inout) :: rvalue
    real(dp) :: sLocal(1), sGlobal(1)

    sLocal(1) = rvalue

    call maxRealParallel(sLocal, sGlobal, 1);

    rvalue = sGlobal(1)

  end subroutine maxRealReduce

  !
  ! Integer max reduce for 2 values
  !
  subroutine maxIntReduce2(value1, value2)

    integer,intent(inout) :: value1, value2
    integer :: sLocal(2), sGlobal(2)

    sLocal(1) = value1
    sLocal(2) = value2

    call maxIntParallel(sLocal, sGlobal, 2);

    value1 = sGlobal(1)
    value2 = sGlobal(2)

  end subroutine maxIntReduce2

  !
  ! Integer sum reduce for 2 values
  !
  subroutine sumIntReduce2(value1, value2)

    integer,intent(inout) :: value1, value2
    integer :: sLocal(2), sGlobal(2)

    sLocal(1) = value1
    sLocal(2) = value2

    call sumIntParallel(sLocal, sGlobal, 2);

    value1 = sGlobal(1)
    value2 = sGlobal(2)

  end subroutine sumIntReduce2

  !
  ! Real sum reduce for 1 value
  !
  subroutine sumRealReduce(value1)

    real(dp),intent(inout) :: value1
    real(dp):: sLocal(1), sGlobal(1)

    sLocal(1) = value1

    call sumRealParallel(sLocal, sGlobal, 1);

    value1 = sGlobal(1)

  end subroutine sumRealReduce

  !
  ! Real sum reduce for 2 values
  !
  subroutine sumRealReduce2(value1, value2)

    real(dp),intent(inout) :: value1, value2
    real(dp):: sLocal(2), sGlobal(2)

    sLocal(1) = value1
    sLocal(2) = value2

    call sumRealParallel(sLocal, sGlobal, 2);

    value1 = sGlobal(1)
    value2 = sGlobal(2)

  end subroutine sumRealReduce2

  !
  ! Real sum reduce for 3 values
  !
  subroutine sumRealReduce3(value1, value2, value3)

    real(dp),intent(inout) :: value1, value2, value3
    real(dp):: sLocal(3), sGlobal(3)

    sLocal(1) = value1
    sLocal(2) = value2
    sLocal(3) = value3

    call sumRealParallel(sLocal, sGlobal, 3);

    value1 = sGlobal(1)
    value2 = sGlobal(2)
    value3 = sGlobal(3)

  end subroutine sumRealReduce3

  !
  ! Real sum reduce for N values
  !
  subroutine sumRealReduceN(valueVec, N)

    real(dp),intent(inout) :: valueVec(N)
    integer, intent(in) :: N

    real(dp), allocatable :: sGlobal(:)

    allocate(sGlobal(N))
    
    call sumRealParallel(valueVec, sGlobal, N);

    valueVec = sGlobal
   
    deallocate(sGlobal)

  end subroutine sumRealReduceN

  
  !
  ! Real sum reduce for Int values
  !
  subroutine sumIntReduceN(valueVec, N)

    integer, intent(inout) :: valueVec(N)
    integer, intent(in) :: N

    integer, allocatable :: sGlobal(:)

    allocate(sGlobal(N))

    call sumIntParallel(valueVec, sGlobal, N);

    valueVec = sGlobal

    deallocate(sGlobal)

  end subroutine sumIntReduceN
  
  !
  ! Wrapper for MPI_Allreduce real min with rank.
  !
  subroutine minRankRealParallel(sendBuf, recvBuf, icount)

    type(rankReduceData_t), intent(in) :: sendBuf(*)
    type(rankReduceData_t), intent(out) :: recvBuf(*)
    integer, intent(in) :: icount

#ifdef DO_MPI
   call MPI_Allreduce(sendBuf, recvBuf, icount, MPI_DOUBLE_INT, & 
                      MPI_MINLOC, MPI_COMM_WORLD, ierr)
#else
   integer :: i

   do i = 1, icount
      recvBuf(i)%val = sendBuf(i)%val
      recvBuf(i)%rank = sendBuf(i)%rank
   enddo
#endif

  end subroutine minRankRealParallel

  !
  ! Wrapper for MPI_Allreduce real max with rank.
  !
  subroutine maxRankRealParallel(sendBuf, recvBuf, icount)

    type(rankReduceData_t), intent(in) :: sendBuf(*)
    type(rankReduceData_t), intent(out) :: recvBuf(*)
    integer, intent(in) :: icount

#ifdef DO_MPI
   call MPI_Allreduce(sendBuf, recvBuf, icount, MPI_DOUBLE_INT, &
                      MPI_MAXLOC, MPI_COMM_WORLD, ierr)
#else
   integer :: i

   do i = 1, icount
      recvBuf(i)%val = sendBuf(i)%val
      recvBuf(i)%rank = sendBuf(i)%rank
   enddo
#endif

  end subroutine maxRankRealParallel

  !
  ! Wrapper for MPI broadcast
  !
  subroutine bcastParallel(buf, blen, root)

    character, intent(in) :: buf(*)
    integer, intent(in) :: blen, root

#ifdef DO_MPI
   call MPI_Bcast(buf, blen, MPI_BYTE, root, MPI_COMM_WORLD, ierr)
#endif

  end subroutine bcastParallel

  !
  ! Wrapper for real MPI_AllGather
  !
  subroutine allGatherRealParallel(sendBuf, sendLen, recvBuf, recvLen)

    real(dp), intent(in) :: sendBuf(*)
    real(dp), intent(out) :: recvBuf(*)
    integer, intent(in) :: sendLen, recvLen
    
#ifdef DO_MPI
    call MPI_Allgather(sendBuf, sendLen, REAL_MPI_TYPE, recvBuf, recvLen, &
                       REAL_MPI_TYPE, MPI_COMM_WORLD, ierr)
#endif

  end subroutine allGatherRealParallel

  !
  ! Wrapper for integer MPI_AllGather
  !
  subroutine allGatherIntParallel(sendBuf, sendLen, recvBuf, recvLen)

    integer, intent(in) :: sendBuf(*)
    integer, intent(out) :: recvBuf(*)
    integer, intent(in) :: sendLen, recvLen

#ifdef DO_MPI
  call MPI_Allgather(sendBuf, sendLen, MPI_INT, recvBuf, recvLen, &
                     MPI_INT, MPI_COMM_WORLD, ierr)
#endif

  end subroutine allGatherIntParallel

  !
  ! Wrapper for real MPI_AllGatherV
  !
  subroutine allGatherVRealParallel(sendBuf, sendLen, recvBuf, recvLen, recvDispl)

    real(dp), intent(in) :: sendBuf(*)
    real(dp), intent(out) :: recvBuf(*)
    integer, intent(in) :: sendLen
    integer, intent(in) :: recvLen(*)
    integer, intent(in) :: recvDispl(*)

#ifdef DO_MPI
    call MPI_Allgatherv(sendBuf, sendLen, REAL_MPI_TYPE, recvBuf, &
                        recvLen, recvDispl, REAL_MPI_TYPE, &
                        MPI_COMM_WORLD, ierr)
#endif


  end subroutine allGatherVRealParallel

  !
  ! Wrapper for integer MPI_AllGatherV
  !
  subroutine allGatherVIntParallel(sendBuf, sendLen, recvBuf, recvLen, recvDispl)

    integer, intent(in) :: sendBuf(*)
    integer, intent(out) :: recvBuf(*)
    integer, intent(in) :: sendLen 
    integer, intent(in) :: recvLen(*) 
    integer, intent(in) :: recvDispl(*)

#ifdef DO_MPI
    call MPI_Allgatherv(sendBuf, sendLen, MPI_INT, recvBuf, recvLen, &
                        recvDispl, MPI_INT, MPI_COMM_WORLD, ierr)
#endif

  end subroutine allGatherVIntParallel

  !
  ! Wrapper for real MPI_AllReduce
  !
  subroutine allSumRealReduceParallel(buf, buflen)

    real(dp), intent(inout) :: buf(*)
    integer, intent(in) :: buflen

#ifdef DO_MPI
    call MPI_Allreduce(MPI_IN_PLACE, buf, buflen, REAL_MPI_TYPE, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

  end subroutine allSumRealReduceParallel

  !
  ! Wrapper for integer MPI_AllReduce
  !
  subroutine allSumIntReduceParallel(buf, buflen)

    integer, intent(inout) :: buf(*)
    integer, intent(in) :: buflen

#ifdef DO_MPI
    call MPI_Allreduce(MPI_IN_PLACE, buf, buflen, MPI_INT, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

  end subroutine allSumIntReduceParallel

  !
  ! Wrapper for bml_allGatherVParallel
  !
  subroutine allGatherParallel(a)

    type (bml_matrix_t), intent(inout) :: a

#ifdef DO_MPI
    call bml_allGatherVParallel(a)
#endif

  end subroutine allGatherParallel


  !
  ! Wrapper for MPI_wait
  !
  subroutine wait()

    integer :: status(3)
    integer :: request, ierr

#ifdef DO_MPI
!     call MPI_WAIT(request, status, ierr)
#endif

  end subroutine wait
    

  
end module parallel_mod
