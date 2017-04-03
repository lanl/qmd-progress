!> The timer module.
!! \ingroup PROGRESS
!!
!! Example use of dynamic timing: 
!!
!!    call timer_prg_init()
!!     
!!     call prg_timer_start(dyn_timer,"timer_tag")
!!
!!     .... code lines ... 
!!
!!     call prg_timer_stop(dyn_timer,1)  
!!
!!
!! This will write the time it takes to execute "code lines" and it will name it "timer_tag"
!!
    !
    ! Timer routines.
    !

module prg_timer_mod

  use prg_parallel_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  public :: timer_status_t
  public :: timer_prg_init
  public :: prg_timer_shutdown
  public :: prg_timer_start
  public :: prg_timer_stop
  public :: prg_timer_collect
  public :: prg_timer_results
  public :: time2milliseconds
  public :: prg_print_date_and_time

  integer, public :: loop_timer, sp2_timer, genx_timer
  integer, public :: part_timer, subgraph_timer, deortho_timer
  integer, public :: ortho_timer, zdiag_timer, graphsp2_timer
  integer, public :: subind_timer, subext_timer, subsp2_timer
  integer, public :: suball_timer, bmult_timer, badd_timer
  integer, public :: dyn_timer, mdloop_timer, buildz_timer
  integer, public :: realcoul_timer, recipcoul_timer, pairpot_timer
  integer, public :: halfverlet_timer, pos_timer, nlist_timer

  !> Timer status type
  type timer_status_t
  
    !> Timer name
    character(LEN=20) :: tname

    !> Start time
    integer :: tstart

    !> Current total time
    integer :: ttotal

    !> Current call count
    integer :: tcount

    !> Rank with min value
    integer :: minRank

    !> Rank with max value
    integer :: maxRank

    !> Sum time - total time in secs
    real(dp):: tsum

    !> Minimum value over all ranks
    real(dp) :: minValue

    !> Maximum value over all ranks
    real(dp) :: maxValue

    !> Average value over all ranks
    real(dp) :: tavg

    !> Stdev across all ranks
    real(dp) :: tstdev

    !> Percent of time across all timers
    real(dp) :: tpercent

  end type timer_status_t

  !
  ! Adding a new timer requires the following.
  !
  ! integer :: new_timer
  !
  ! In prg_init_timer, increment the timer count, specify a number
  ! and name for the new timer.
  !
  ! ...
  !
  ! ! Increment when adding a new timer
  ! num_timers = 4
  !
  !  ...
  !
  ! ! Timer handles, names, and counters
  !  loop_timer = 1
  !  sp2_timer  = 2
  !  genx_timer = 3
  !  new_timer  = 4
  !
  !  ptimer(loop_timer)%tname = "Loop"
  !  ptimer(sp2_timer)%tname  = "  SP2"
  !  ptimer(genx_timer)%tname = "  GenX" 
  !  ptimer(new_timer)%tname  = "  New"
  !

  integer :: tstart_clock, tstop_clock, tclock_rate, tclock_max
  integer :: num_timers

  type (timer_status_t), allocatable :: ptimer(:)

  private :: int2char

contains

  !> Initialize timers
  subroutine timer_prg_init()

    integer :: i

    ! Increment when adding a new timer
    num_timers = 24 

    allocate(ptimer(num_timers))

    ! Timer handles, names, and counters
    loop_timer = 1
    subgraph_timer = 2
    sp2_timer  = 3
    genx_timer = 4
    part_timer = 5
    deortho_timer  = 6
    ortho_timer  = 7
    zdiag_timer  = 8
    graphsp2_timer  = 9
    subind_timer = 10
    subext_timer = 11
    subsp2_timer = 12
    suball_timer = 13
    bmult_timer = 14
    badd_timer = 15
    dyn_timer = 16
    mdloop_timer = 17
    buildz_timer = 18
    realcoul_timer = 19
    recipcoul_timer = 20
    pairpot_timer = 21
    halfverlet_timer = 22
    pos_timer = 23
    nlist_timer = 24
    
    ptimer(loop_timer)%tname       = "Loop"
    ptimer(subgraph_timer)%tname   = " Subgraph"
    ptimer(sp2_timer)%tname        = "  SP2"
    ptimer(genx_timer)%tname       = "  GenX"
    ptimer(part_timer)%tname       = "  Part"
    ptimer(deortho_timer)%tname    = "  Deortho"
    ptimer(ortho_timer)%tname      = "  Ortho"
    ptimer(zdiag_timer)%tname      = "  Zdiag"
    ptimer(graphsp2_timer)%tname   = "  GraphSP2"
    ptimer(subind_timer)%tname     = "   SubInd"
    ptimer(subext_timer)%tname     = "   SubExt"
    ptimer(subsp2_timer)%tname     = "   SubSP2"
    ptimer(suball_timer)%tname     = "   SubAll"
    ptimer(bmult_timer)%tname      = "    BMult"
    ptimer(badd_timer)%tname       = "    BAdd"
    ptimer(dyn_timer)%tname        = "        " !Reserved for dynamic timing  
    ptimer(mdloop_timer)%tname     = "MDLoop"   
    ptimer(buildz_timer)%tname     = "  BuildZ"   
    ptimer(realcoul_timer)%tname   = "  RealCoul"   
    ptimer(recipcoul_timer)%tname  = "  RecipCoul"   
    ptimer(pairpot_timer)%tname    = "  PairPot"   
    ptimer(halfverlet_timer)%tname = "  HalfVerlet"   
    ptimer(pos_timer)%tname        = "  Pos"   
    ptimer(nlist_timer)%tname      = "  NList"   
    
    do i = 1, num_timers
      ptimer(i)%ttotal = 0
      ptimer(i)%tcount = 0
    end do

  end subroutine timer_prg_init

  !> Get timer id
  subroutine prg_timer_getid()
  
  end subroutine prg_timer_getid
  
  !> Done with timers
  subroutine prg_timer_shutdown()

    deallocate(ptimer)

  end subroutine prg_timer_shutdown

  !> Start Timing
  !!
  !! \param itimer The index of the timer to start.
  !! \param tag Optional parameter to retag the timer on the fly.
  subroutine prg_timer_start(itimer,tag)

    integer, intent(in) :: itimer
    character(len=*), intent(in), optional :: tag 
    
    if(present(tag))then 
      ptimer(itimer)%tname = tag
    endif
    
    call system_clock(tstart_clock, tclock_rate, tclock_max)
    ptimer(itimer)%tstart = tstart_clock

  end subroutine prg_timer_start

  !> Stop timing
  !!
  !! \param itimer The index of the timer to stop.
  !! \param verbose Optional parameters to print partial times.
  subroutine prg_timer_stop(itimer,verbose)

    integer, intent(in) :: itimer
    integer :: tprg_delta
    integer, intent(in), optional :: verbose

    call system_clock(tstop_clock, tclock_rate, tclock_max)
    tprg_delta = tstop_clock - ptimer(itimer)%tstart
    if(present(verbose))then
      if(verbose.GT.0)then 
        write(*,*)"Time for ",ptimer(itimer)%tname,"=",tprg_delta
      endif
    endif
    ptimer(itimer)%ttotal = ptimer(itimer)%ttotal + tprg_delta
    ptimer(itimer)%tcount = ptimer(itimer)%tcount + 1

  end subroutine prg_timer_stop  

  ! Collect timer results
  !
  subroutine prg_timer_collect()

  integer :: i
  real(dp) :: temp
  real(dp), allocatable :: sendBuf(:), recvBuf(:)
  type(rankReduceData_t), allocatable :: reduceSendBuf(:)
  type(rankReduceData_t), allocatable :: reduceRecvBuf(:)

  real(dp) :: rranks

  allocate(sendBuf(num_timers))
  allocate(recvBuf(num_timers))

  rranks = float(getNRanks())

  !! Determine average of each timer across ranks
  do i = 1, num_timers
    sendBuf(i) = float(ptimer(i)%ttotal)/float(tclock_rate)
  enddo
  call sumRealParallel(sendBuf, recvBuf, num_timers);

  do i = 1, num_timers
    ptimer(i)%tavg = recvBuf(i) / rranks
  enddo

  !! Determine min and max across ranks and which rank
  allocate(reduceSendBuf(num_timers))
  allocate(reduceRecvBuf(num_timers))

  do i = 1, num_timers
    reduceSendBuf(i)%val = float(ptimer(i)%ttotal)/float(tclock_rate)
    reduceSendBuf(i)%rank = getMyRank()
 enddo 
  call minRankRealParallel(reduceSendBuf, reduceRecvBuf, num_timers);
  do i = 1, num_timers
    ptimer(i)%minValue = reduceRecvBuf(i)%val
    ptimer(i)%minRank = reduceRecvBuf(i)%rank
  enddo 
  call maxRankRealParallel(reduceSendBuf, reduceRecvBuf, num_timers);
  do i = 1, num_timers
    ptimer(i)%maxValue = reduceRecvBuf(i)%val
    ptimer(i)%maxRank = reduceRecvBuf(i)%rank
  enddo

  deallocate(reduceSendBuf)
  deallocate(reduceRecvBuf)

  !! Determine standard deviation
  do i = 1, num_timers
    temp = float(ptimer(i)%ttotal)/float(tclock_rate) - ptimer(i)%tavg
    sendBuf(i) = temp * temp;
  enddo
  call sumRealParallel(sendBuf, recvBuf, num_timers);
  do i = 1, num_timers
    ptimer(i)%tstdev = sqrt(recvBuf(i) / rranks)
  enddo

  deallocate(sendBuf)
  deallocate(recvBuf)

  end subroutine prg_timer_collect

  ! Print performance results
  !
  subroutine prg_timer_results()

    integer :: i

    ! Collect results across all ranks
    call prg_timer_collect()

    ! Print timer results
    if (printRank() .eq. 1) then

      write(*,*) ""
      write(*,*) "Timings for Rank ", getMyRank()
      write(*,*) "Timer                 # Calls  Avg/Call (s)     Total (s)       % Time"
      write(*,*) ""

      do i = 1, num_timers
        if (ptimer(i)%tcount .gt. 0) then
    !!      ptimer(i)%tavg = (float(ptimer(i)%ttotal)/float(tclock_rate))/float(ptimer(i)%tcount)
          ptimer(i)%tsum = float(ptimer(i)%ttotal)/float(tclock_rate)
          ptimer(i)%tpercent = (ptimer(i)%tsum / ptimer(1)%tsum) * 100.0
          write(*,10) ptimer(i)%tname, ptimer(i)%tcount, ptimer(i)%tsum/float(ptimer(i)%tcount), ptimer(i)%tsum, ptimer(i)%tpercent
10        format(A23, I6, 3G16.6)
        end if
      end do

      write(*,*) ""
       write(*,*) "Timing Statistics Across ", getNRanks(), " Ranks:"
      write(*,*) "Timer                      Rank: Min(s)        Rank: Max(s)            Avg(s)        Stdev(s)"
      write(*,*)

      do i = 1, num_timers
        if (ptimer(i)%tcount > 0) then
          write(*, 20) ptimer(i)%tname, &
                       ptimer(i)%minRank, ptimer(i)%minValue, &
                       ptimer(i)%maxRank, ptimer(i)%maxValue, &
                       ptimer(i)%tavg, ptimer(i)%tstdev 
20        format(A23,2X,I4,G16.6,I4,3G16.6) 
        endif
      enddo
    endif

  end subroutine prg_timer_results
 
  function time2milliseconds() result(mls)

    real(8) :: mls
    integer :: timevector(8)

    call date_and_time(values=timevector)
    mls = timevector(5)*60*60*1000 + timevector(6)*60*1000 + &
        timevector(7)*1000 + timevector(8)

  end function time2milliseconds

  subroutine prg_print_date_and_time(tag)

    implicit none
    
    character(len=*), intent(in) :: tag
    character(2) :: monthchar, daychar,hourchar,minchar,secchar    
    integer :: sec, mins, hour, day, month, year
    integer :: timevector(8)

    call date_and_time(values=timevector)
    
    year = timevector(1); month = timevector(2); day = timevector(3)
    hour = timevector(5); mins = timevector(6); sec = timevector(7)
    
    monthchar = int2char(month); daychar = int2char(day)
    hourchar = int2char(hour); minchar = int2char(mins); secchar = int2char(sec)
      
    write(*,'(a2,a,x,A2,a1,A2,a1,i4,x,a2,x,A2,a1,A2,a1,A2)')"# ", &
      trim(tag),monthchar,"/" &
      ,daychar,"/",year, "at", hourchar,":",minchar,":",secchar

  end subroutine prg_print_date_and_time
  
  function int2char(ival)
  
    implicit none

    integer, intent(in) :: ival 
    character(2) :: int2char, myintchar
    
    if ((ival/10) .lt. 1) then
      write(myintchar,'(I2)') ival
      myintchar="0"//trim(adjustl(myintchar))
    else
      write(myintchar,'(I2)') ival
    endif  
    
    int2char = myintchar
    
  end function int2char

end module prg_timer_mod
