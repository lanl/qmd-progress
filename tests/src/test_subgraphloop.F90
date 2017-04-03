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

!> The test_subgraphloop module.
    !
    !
module test_prg_subgraphloop_mod

  use prg_graph_mod
  use prg_sp2_mod
  use prg_subgraphloop_mod
  use prg_homolumo_mod
  use prg_timer_mod
  use bml
  use omp_lib

  implicit none

  private  !Everything is private by default

  integer, parameter :: dp = kind(1.0d0)

  public :: test_subgraphloop

  contains

  !! Test subgraph SP2 loop.
  !!
  !! \param mtxFile Input Hamiltonian file
  !! \param N Number of rows/orbitals
  !! \param M Max number of non-zeroes per row
  !! \param threshold Threshold for sparse matrix algebra
  !! \param bndfil Bond 
  !! \param minsp2iter Minimum sp2 iterations
  !! \param maxsp2iter Maximum SP2 iterations
  !! \param sp2conv Convergence type
  !! \param idemtol Idempotency tolerance
  !! \param gthreshold Graph threshold
  !! \param errlimit Error limit for sequence calculation
  !! \param nodesPerPart Number of nodes per partition
  subroutine test_subgraphloop(h_bml, rho_bml, sthreshold, bndfil, &
    minsp2iter, maxsp2iter, sp2conv, idemtol, gthreshold, errlimit, &
    nodesPerPart)

    integer, intent(in) :: nodesPerPart
    integer, intent(in) :: minsp2iter, maxsp2iter
    real(dp), intent(in) :: sthreshold, gthreshold, bndfil, idemtol
    real(dp), intent(in) :: errlimit
    character(len=*), intent(in) :: sp2conv
    type(bml_matrix_t), intent(in) :: h_bml
    type(bml_matrix_t), intent(inout) :: rho_bml

    type(bml_matrix_t) :: g_bml
    type(graph_partitioning_t) :: gp

    integer :: N, M, i, icount
    integer :: pp(100)
    real(dp) :: vv(100)
    real(dp) :: mineval, maxeval, traceMult
    real(dp) :: ehomo, elumo, egap
    real(dp), allocatable :: gbnd(:)

    N = bml_get_N(h_bml)
    M = bml_get_M(h_bml)

    !! Create the same size matrix for the connectivity graph
    call bml_zero_matrix(BML_MATRIX_ELLPACK, BML_ELEMENT_REAL, dp, &
          N, M, g_bml);

    !! Run regular SP2 
    call prg_timer_start(sp2_timer)
    call prg_sp2_alg2_genseq(h_bml, g_bml, sthreshold, bndfil, minsp2iter, &
           maxsp2iter, sp2conv, idemtol, pp, icount, vv) 
    call prg_timer_stop(sp2_timer)

    !! List fnorm per iteration
    do i = 1, icount
     write(*,*) "iter = ", i, " vv = ", vv(i) 
    enddo

    !! Calculate Trace[HP]
    !!traceMult = bml_traceMult(h_bml, g_bml)
    !write(*,*) "Trace[HP] for SP2 = ", traceMult

    !! Calculate Homo-Lumo gap
    allocate(gbnd(2))
    call bml_gershgorin(h_bml, gbnd)
    mineval = gbnd(1)
    maxeval = gbnd(2)
    deallocate(gbnd)
    write(*,*) "Gershgorin: mineval = ", mineval, " maxeval = ", maxeval
    call prg_homolumogap(vv, icount, pp, mineval, maxeval, ehomo, elumo, egap)
    write(*,*) "Homo-lumo: ehomo = ", ehomo, " elumo = ", elumo, &
           " egap = ", egap
           
 
    !! Start graph partitioning SP2 

    !! Threshold the graph
    call bml_threshold(g_bml, gthreshold)

    !! Create graph partitioning of equal parts
    call prg_equalPartition(gp, nodesPerPart, N)
    gp%mineval = mineval
    gp%maxeval = maxeval

    !! Calculate SP2 sequence
    call prg_sp2sequence(gp%pp, gp%maxIter, mineval, maxeval, ehomo, elumo, &
          errlimit)
    write(*,*) "SP2Sequence: Max iterations = ", gp%maxIter

    !! Run subgraph SP2 loop
    call prg_timer_start(graphsp2_timer)
    call prg_subgraphSP2Loop(h_bml, g_bml, rho_bml, gp, sthreshold) 
    call prg_timer_stop(graphsp2_timer)

    !! Calculate fnorm across subgraphs per iteration
    call prg_fnormGraph(gp)
    
    !! Calculate homo-lumo gap
    call prg_homolumogap(gp%vv, gp%maxIter, gp%pp, gp%mineval, gp%maxeval, ehomo, elumo, &
         egap)
    write(*,*) "Homo-lumo: ehomo = ", ehomo, " elumo = ", elumo, &
           " egap = ", egap

    !! Calculate trace[HP]
    !traceMult = bml_traceMult(h_bml, rho_bml)
    !write(*,*) "Trace[HP] for subgraph SP2 = ", traceMult

    call bml_deallocate(g_bml)

    call prg_destroyGraphPartitioning(gp)

  end subroutine test_subgraphloop

end module test_prg_subgraphloop_mod
