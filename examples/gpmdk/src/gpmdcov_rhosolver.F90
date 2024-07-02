module gpmdcov_RhoSolver_mod

  contains 
 !> Solver for computing the density matrix
  subroutine gpmdcov_RhoSolver(orthoh_bml,orthop_bml,evects_bml)
    use gpmdcov_vars
    use gpmdcov_mod
    type(bml_matrix_t), intent(in) :: orthoh_bml
    type(bml_matrix_t), intent(inout) :: orthop_bml
    type(bml_matrix_t), intent(inout) :: evects_bml

    if(lt%verbose >= 1 .and. myRank == 1) write(*,*)"starting solver ..."

    if(lt%method.EQ."GSP2")then
      call prg_subgraphSP2Loop(orthoh_bml, g_bml, orthop_bml, gp, lt%threshold)
      ! call prg_sp2_alg1_seq(orthoh_bml,orthop_bml,lt%threshold, gp%pp, gp%maxIter, gp%vv)
    elseif(lt%method.EQ."SP2")then
      call prg_sp2_alg2(orthoh_bml,orthop_bml,lt%threshold, bndfil, sp2%minsp2iter, sp2%maxsp2iter &
           ,sp2%sp2conv,sp2%sp2tol,lt%verbose)
    elseif(lt%method.EQ."Diag")then
      ! call build_density_t0(orthoh_bml,orthop_bml,lt%threshold,bndfil)
      ! call prg_build_density_T(orthoh_bml,orthop_bml,lt%threshold,bndfil, 0.1_dp, Ef)
      !call prg_build_density_T_Fermi(orthoh_bml,orthop_bml,lt%threshold, 0.1_dp, Ef)
      call prg_build_density_T_Fermi(orthoh_bml,orthop_bml,lt%threshold, 0.1_dp, Ef)
      if(lt%verbose >= 1 .and. myRank == 1) write(*,*)"ipt =",ipt,"Ef =",Ef
    elseif(lt%method.EQ."DiagEf")then
      if(bml_get_n(evects_bml) < 0) call bml_deallocate(evects_bml)
      call bml_zero_matrix(lt%bml_type,bml_element_real,dp,norb,norb,evects_bml)      
      call prg_build_density_T_ed(orthoh_bml,orthop_bml,evects_bml,lt%threshold,bndfil, lt%kbt, Ef, &
           & syprt(ipt)%estr%evals, syprt(ipt)%estr%dvals, syprt(ipt)%estr%hindex, gpat%sgraph(ipt)%llsize,lt%verbose)
    else
      stop "No valid Method in LATTE parameters"
    endif



    if(lt%verbose >= 2 .and. myRank == 1)then
      call bml_print_matrix("orthop_bml",orthop_bml,0,6,0,6)
    endif
    if(lt%verbose >= 1 .and. myRank == 1) write(*,*)"leaving solver ..."

  end subroutine gpmdcov_RhoSolver

end module gpmdcov_RhoSolver_mod
