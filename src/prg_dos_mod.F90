!> A module to compute the Density of state (DOS) and lDOS.
!! \brief This module will be used to compute DOS and lDOS.
!!
!! @todo Add LDOS.
!! @ingroup PPROGRESS
!!
!!
module prg_dos_mod

  use prg_openfiles_mod
  use prg_ptable_mod

  implicit none

  private

  integer, parameter :: dp = kind(1.0d0)

  public :: prg_write_tdos

contains

  !> Writes the total DOS into a file.
  !! \f$ DOS(\epsilon) = \sum_k L(\epsilon - \epsilon_k) \f$
  !! Where \f$ \int_{-\infty}^{\infty} DOS(\epsilon) = Nstates \f$
  !! \note DOS is NOT shifted respect to Ef.
  !! \param eigenvals Eigenvalues of the system.
  !! \para gamma Lorentzian width.
  !! \param npts Number of energy points.
  !! \param emin Minimum energy value.
  !! \param emax Maximum energy value.
  !! \param filename Filename to write the DOS.
  !!
  subroutine prg_write_tdos(eigenvals, gamma, npts, emin, emax, filename)
    implicit none
    character(len=*),intent(in)       ::  filename
    integer                ::  i, io
    integer, intent(in)    ::  npts
    real(dp)               ::  de
    real(dp), allocatable  ::  loads(:)
    real(dp), intent(in)   ::  eigenvals(:), emax, emin, gamma

    call prg_open_file(io,filename)

    de = (emax-emin)/real(npts)

    allocate(loads(size(eigenvals)))

    loads = 1.0_dp
    write(io,*)"#  Energy    DOS"
    do i = 1, npts
      write(io,*) emin + de*i,lorentz(emin + de*i, eigenvals, loads, gamma)
    end do

    close(io)

    call prg_open_file(io, "eigenvals")

    write(io,*)"#  i   Eval"
    do i = 1, size(eigenvals,dim=1)
      write(io,*) i, eigenvals(i)
    end do

    close(io)

  end subroutine prg_write_tdos


  !> Lorentzian Function
  !! \brief Computes:
  !! \f$ L(\epsilon) = \sum_{k} \frac{\omega(k)\Gamma}{2 \pi}\frac{1}{(\epsilon - \epsilon_k)^2 + (\Gamma/2)^2} \f$
  !! \param energy Energy point.
  !! \param eigenvals Eigenvalues of the system.
  !! \param Gamma Lorentz function broadening.
  !!
  real(dp) function lorentz(energy, eigenvals, loads, Gamma)
    implicit none
    integer               ::  Nstates, k
    real(dp)              ::  auxfactor, auxterm, pi
    real(dp), intent(in)  ::  Gamma, eigenvals(:), energy, loads(:)

    Nstates = size(eigenvals,dim=1)
    pi = 3.14159265358979323846264338327950_dp

    !Lorentz parameters
    auxfactor = Gamma/(2.0_dp*pi)
    auxterm = (Gamma/2.0_dp)**2
    lorentz = 0.0_dp

    do k = 1, Nstates
      lorentz = lorentz + loads(k)/((energy-eigenvals(k))**2 + auxterm)
    end do

    lorentz = auxfactor*lorentz

  end function lorentz


end module prg_dos_mod
