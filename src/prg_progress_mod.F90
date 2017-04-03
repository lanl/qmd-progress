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

  public :: prg_progress_init
  public :: prg_progress_shutdown

  contains

  !> Initialize progress.
  subroutine prg_progress_init()

    ! Initialize MPI
    call prg_initParallel()

    ! Initialize timers
    call timer_prg_init()
    call prg_timer_start(loop_timer)

  end subroutine prg_progress_init

  !> Shutdown progress.
  subroutine prg_progress_shutdown()

    ! Timer report and finalize
    call prg_timer_stop(loop_timer)
    call prg_timer_results()
    call prg_timer_shutdown()

    ! Finalize MPI
    call prg_shutdownParallel()

  end subroutine prg_progress_shutdown

end module prg_progress_mod
