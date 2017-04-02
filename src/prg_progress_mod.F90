!> The progress module.
!! \ingroup PROGRESS 
    !
    !

module prg_progress_mod

  use bml
  use prg_parallel_mod
  use prg_timer_mod

  implicit none

  private  !Everything is private by default

  integer, parameter :: dp = kind(1.0d0)

  public :: progress_init
  public :: progress_shutdown

  contains

  !> Initialize progress.
  subroutine progress_init()

    ! Initialize MPI
    call initParallel()

    ! Initialize timers
    call timer_init()
    call timer_start(loop_timer)

  end subroutine progress_init

  !> Shutdown progress.
  subroutine progress_shutdown()

    ! Timer report and finalize
    call timer_stop(loop_timer)
    call timer_results()
    call timer_shutdown()

    ! Finalize MPI
    call shutdownParallel()

  end subroutine progress_shutdown

end module prg_progress_mod
