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

!> The normalize module.
    !
    !

module normalize_mod

  use bml
  use graph_mod
  use parallel_mod

  implicit none

  private  !Everything is private by default

  integer, parameter :: dp = kind(1.0d0)

  public :: normalize
  public :: gershgorinReduction

contains

  !> Normalize a Hamiltonian matrix prior to running the SP2 algorithm.
  !!
  !! X0 = (e_max * I - H) / (e_max - e_min)
  !!
  !! where e_max and e_min are obtained sing the Gershgorin circle theorem.
  !!
  !! \param H_bml Input/Output Hamiltonian matrix
  subroutine normalize(H_bml)

    type(bml_matrix_t),intent(inout) :: H_bml

    real(dp) :: alpha, beta, maxMinusMin
    real(dp), allocatable :: gbnd(:)

    allocate(gbnd(2))

    call bml_gershgorin(H_bml, gbnd)
    maxMinusMin = gbnd(2) - gbnd(1)
    alpha = -1.00_dp / maxMinusMin
    beta = gbnd(2) /maxMinusMin 
    call bml_scale_add_identity(H_bml, alpha, beta, 0.00_dp)

    deallocate(gbnd)

  end subroutine normalize

  !> Determine gershgorin bounds across all parts, local and distributed.
  subroutine gershgorinReduction(gp)

    implicit none

    type(graph_partitioning_t), intent(inout) :: gp

    integer :: i, nRanks
    real(dp) :: localVal

    nRanks = getNRanks()
    if (nRanks .eq. 1) return

    localVal = gp%mineval
    call minRealReduce(localVal)
    gp%mineval = localVal

    localVal = gp%maxeval
    call maxRealReduce(localVal)
    gp%maxeval = localVal
 
  end subroutine gershgorinReduction

end module normalize_mod
