!> Solver for computing the density matrix
subroutine gpmd_RhoSolver(orthoh_bml,orthop_bml)

  type(bml_matrix_t), intent(in) :: orthoh_bml
  type(bml_matrix_t), intent(inout) :: orthop_bml

  if(lt%verbose >= 1 .and. myRank == 1) write(*,*)"In solver ..."
  if(lt%verbose >= 1 .and. myRank == 1) call prg_timer_start(dyn_timer,"Solver")

  if(lt%method.EQ."GSP2")then
    call prg_timer_start(graphsp2_timer)
    call prg_subgraphSP2Loop(orthoh_bml, g_bml, orthop_bml, gp, lt%threshold)
    call prg_timer_stop(graphsp2_timer)
    ! call prg_sp2_alg1_seq(orthoh_bml,orthop_bml,lt%threshold, gp%pp, gp%maxIter, gp%vv)
  elseif(lt%method.EQ."SP2")then
    call prg_sp2_alg2(orthoh_bml,orthop_bml,lt%threshold, bndfil, sp2%minsp2iter, sp2%maxsp2iter &
      ,sp2%sp2conv,sp2%sp2tol,lt%verbose)
  elseif(lt%method.EQ."Diag")then
    ! call build_density_t0(orthoh_bml,orthop_bml,lt%threshold,bndfil)
    ! call prg_build_density_T(orthoh_bml,orthop_bml,lt%threshold,bndfil, 0.1_dp, Ef)
    call prg_build_density_T_Fermi(orthoh_bml,orthop_bml,lt%threshold, 0.1_dp, Ef)
    if(lt%verbose >= 1 .and. myRank == 1) write(*,*)"ipt =",ipt,"Ef =",Ef
  else
    stop"No valid Method in LATTE parameters"
  endif

  if(lt%verbose >= 1 .and. myRank == 1) call prg_timer_stop(dyn_timer,1)

  ! #ifdef DO_MPI_BLOCK
  !     call prg_allGatherParallel(orthop_bml)
  ! #endif

  if(lt%verbose >= 2 .and. myRank == 1)then
    call bml_print_matrix("orthop_bml",orthop_bml,0,6,0,6)
  endif

end subroutine gpmd_RhoSolver
