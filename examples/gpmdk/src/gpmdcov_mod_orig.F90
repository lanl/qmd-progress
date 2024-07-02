!> Auxiliary modules.
!! \brief This module will be used to have auxiliary routines.
!!
!! @ingroup GPMDCOV
!!
!!
module gpmdcov_mod_orig

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  public :: gpmdcov_musearch

contains


  !> Perform musearch.
  !! \param evals Eigenvalues of the system.
  !! \para fvals Occupations of the states.
  !! \param dvals Partial traces.
  !! \prarm beta Electronic temperature.
  !! \param maxMuIter Maximum numner of iterations.
  !! \param occTol Occupation tolerance.
  !! \param mu Chemical potential.
  !!
  subroutine gpmdcov_musearch(evals, fvals, dvals, beta, nocc, maxMuIter, occTol, mu, verbose)
    implicit none
    integer                ::  i, j, norbs
    integer, intent(in)    ::  maxMuIter
    real(dp)               ::  den, occErr, mu0, occ
    real(dp), intent(in)   ::  evals(:), dvals(:)
    real(dp), intent(inout)   ::  fvals(:)
    real(dp), intent(in)   ::  occTol, beta, nocc
    real(dp), intent(inout)   ::  mu
    integer, optional, intent(in) :: verbose


    norbs = size(evals, dim = 1)

    mu0 = mu
    write(*,*)"NOCC",nocc, norbs, mu0

    do j = 1, norbs
      fvals(j) = 1.0_dp/(exp(beta*(evals(j)-mu0))+1.0_dp)
    end do

    occ = 0.0_dp
    do j = 1, norbs
      !                write(*,*)"fvals, dvals",j,fvals(j),dvals(j), evals(j)
      occ = occ + fvals(j)*dvals(j)
    end do

    !stop
    occErr = abs(occ-nocc)

    if (present(verbose)) then
      if (verbose >= 1) write(*,*)"In gpmdcov_musearch ..."
    endif

    do i = 1, maxMuIter
      do j = 1, norbs
        fvals(j) = 1.0_dp/(exp(beta*(evals(j)-mu0))+1.0_dp)
      end do

      occ = 0.0_dp
      den = 0.0_dp
      do j = 1, norbs
        occ = occ + fvals(j)*dvals(j)
        den = den + beta*fvals(j)*(1.0_dp-fvals(j))*dvals(j)
      end do

      occErr = abs(occ - nocc)
      mu0 = mu0 + (nocc - occ)/den
      mu = mu0

      write(*,*)"occ, occErr, nocc, den, mu", occ, occErr, nocc, den, mu
      if (i .eq. maxMuIter) then
        write(*,*) "Max mu iters reached"
        stop
      else
        if(occErr .lt. occTol) exit
      end if

    end do

  end subroutine gpmdcov_musearch

end module gpmdcov_mod_orig
