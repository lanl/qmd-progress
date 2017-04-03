!> The homolumo module.
!! \ingroup PROGRESS
    !
    !

module prg_homolumo_mod

  use bml

  implicit none

  private  !Everything is private by default

  integer, parameter :: dp = kind(1.0d0)

  public :: prg_homolumogap
  public :: prg_sp2sequence

contains

  !! Calculate the homo-lumo gap
  subroutine prg_homolumogap(vv, imax, pp, mineval, maxeval, ehomo, elumo, &
      egap, verbose)

    integer, intent(in) :: imax
    integer, intent(in), optional :: verbose
    integer, intent(in) :: pp(:)
    real(dp), intent(in) :: vv(:)
    real(dp), intent(in) :: mineval, maxeval
    real(dp), intent(inout) :: ehomo, elumo, egap

    integer :: i, j
    real(dp) :: precomp, x_a, x_b, y_a, y_b, hgamma

    x_a = 0.0_dp
    x_b = 1.0_dp
    y_a = 0.0_dp
    y_b = 0.0_dp

    i = imax
    do while (vv(i) .le. 0.0)
      i=i-1      
    enddo

    hgamma = 6.0_dp - 4.0_dp * sqrt(2.0_dp)
    hgamma = hgamma * (1.0_dp - hgamma)
    
    if(present(verbose).and.verbose.GE.1) write(*,*)"In prg_homolumogap ..."

    do while (vv(i) .lt. hgamma)

        precomp = sqrt(1.0_dp - 4.0_dp * vv(i))

        y_a = 0.5_dp * (1.0_dp + precomp)
        y_b = 0.5_dp * (1.0_dp - precomp)

        ! write(*,*)"ya,yb",y_a,y_b
        ! write(*,*)"y_a = ",y_a,"y_a = ",y_b,"vv(i) = ",vv(i)
        do j = i-1, 1, -1
            if (pp(j) .gt. 0) then
                y_a = sqrt(y_a)
                y_b = sqrt(y_b)
            else
                y_a = 1.0_dp - sqrt(1.0_dp - y_a)
                y_b = 1.0_dp - sqrt(1.0_dp - y_b)
            endif
       enddo

       x_a = max(x_a, y_a)
       x_b = min(x_b, y_b)

       if(present(verbose).and.verbose.GE.2) write(*,*)"x_a = ",x_a,"x_b = ",x_b

       i = i - 1
       if (i .lt. 1) then
           write(*,*) "prg_homolumogap error: i < 1, i = ", i
       endif
    enddo

    ehomo = maxeval - x_a * (maxeval - mineval)
    elumo = maxeval - x_b * (maxeval - mineval)

    egap = elumo - ehomo

  end subroutine prg_homolumogap

  !! Calculate the SP2 sequence given thw min/max evals and homo amd lumo
  !! values. Determine sequence of X^2 and 2X-X^2 for SP2.
  !!
  !! pp(i) = 0 -> 2X-X^2
  !! pp(i) = 1 -> x^2
  !!
  subroutine prg_sp2sequence(pp, imax, mineval, maxeval, ehomo, elumo, errlimit, verbose)

    integer, intent(inout) :: imax
    integer, intent(inout) :: pp(:) 
    real(dp), intent(in) :: mineval, maxeval, ehomo, elumo 
    real(dp), intent(in) :: errlimit
    integer, intent(in), optional :: verbose

    integer :: it
    real(dp) :: eh, el, error, sgm
   
    eh = (maxeval - ehomo) / (maxeval - mineval)
    el = (maxeval - elumo) / (maxeval - mineval)
    error = 1.0_dp

    it = 0
    
    if(present(verbose).and.verbose.GE.1) write(*,*)"In prg_sp2sequence ..."

    do while (error .gt. errlimit)
        it = it + 1
        
        if(present(verbose).and.verbose.GE.2) write(*,*)"error = ",error
        
        if ((abs(1.0_dp - eh * eh) + abs(el * el)) .lt. &
            (abs(1.0_dp - (2.0_dp * eh - eh * eh) + &
             abs(2.0_dp * el - el * el)))) then
            pp(it) = 1
        else
            pp(it) = 0
        endif  

        sgm = 1.0_dp - 2.0_dp * pp(it)
        eh = eh + sgm * (eh - eh * eh)
        el = el + sgm * (el - el * el)
        error = abs(1.0_dp - eh) + abs(el)

        if (it .ge. 100) then
            error = 0.0_dp
            write(*,*) "prg_sp2sequence error: SP2 not converging"
        endif

    enddo

    imax = it

  end subroutine prg_sp2sequence

end module prg_homolumo_mod
