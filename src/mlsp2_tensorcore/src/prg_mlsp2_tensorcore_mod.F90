!> The SP2 Tensor Core module.
!! \ingroup PROGRESS
!!
!! \brief This subroutine implements SP2 density matrix
!!  purification algorithm using Tenso Core devices.
!!
module prg_mlsp2_tensorcore_mod

  use, intrinsic :: iso_c_binding !C interface

  implicit none
  public :: prg_mlsp2_tensorcore_f , prg_mlsp2_tensorcore_C

  interface

    subroutine prg_mlsp2_tensorcore_C(N,H,D,beta,mu,verbose) bind(C, name="prg_mlsp2_tensorcore")
      import :: C_PTR, C_INT, C_FLOAT, C_CHAR, C_DOUBLE
      type(C_PTR), value :: D
      type(C_PTR), value :: H
      real(C_DOUBLE), value, intent(in) :: beta, mu
      integer(C_INT), value, intent(in) :: N, verbose
    end subroutine prg_mlsp2_tensorcore_C

  end interface

contains

  !> Calculates the density matrix from a Hamiltonian matrix by
  !! purification.
  !!
  !! \param N Number of orbitals (Size of the Hamiltonian).
  !! \param H Input Hamiltonian matrix.
  !! \param D Output density matrix.
  !! \param threshold Threshold for sparse matrix algebra.
  !! \param bndfil Band filling fraction.
  !! \param minsp2iter Minimum sp2 iterations.
  !! \param maxsp2iter Maximum SP2 iterations.
  !! \param sp2conv Convergence type.
  !! \param idemtol Idempotency tolerance.
  !! \param verbose A verbosity level.
  subroutine prg_mlsp2_tensorcore_f(N,H,D,beta,mu,verbose)
    integer(C_INT), intent(in) :: N, verbose
    real(C_DOUBLE), target, intent(inout) :: D(*)
    real(C_DOUBLE), target, intent(in) :: H(*)
    real(C_DOUBLE), intent(in) :: beta, mu

    !Call the interface
    call prg_mlsp2_tensorcore_C(N,c_loc(H),c_loc(D),beta,mu,verbose)

  end subroutine prg_mlsp2_tensorcore_f

end module prg_mlsp2_tensorcore_mod
