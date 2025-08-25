subroutine gpmdcov_buildz(overin_bml,zmatout_bml)
  use gpmdcov_vars
#ifdef USE_NVTX
    use prg_nvtx_mod
#endif

  type(bml_matrix_t), intent(inout) :: zmatout_bml
  type(bml_matrix_t), intent(inout) :: overin_bml

  igenz = igenz + 1
  

  if(lt%zmat == "ZSP")then !Congruence transformation.

#ifdef USE_NVTX
    call nvtxStartRange("prg_buildzsparse",5)
#endif
    call prg_buildzsparse(overin_bml,zmatout_bml,igenz,lt%mdim,&
         lt%bml_type, zk1_bml,zk2_bml,zk3_bml&
         ,zk4_bml,zk5_bml,zk6_bml,zsp%nfirst,zsp%nrefi,zsp%nreff,&
         zsp%numthresi,zsp%numthresf,zsp%integration,zsp%verbose)

  else

#ifdef USE_NVTX
    call nvtxStartRange("prg_buildzdiag",5)
#endif
    !Build Z matrix using diagonalization (usual method).
    call prg_buildzdiag(overin_bml,zmatout_bml,lt%threshold,lt%mdim,lt%bml_type,0,err_status)
    if(err_status)then
      if(.not. lib_mode)then
        stop "ERROR: Possible non positive definite overlap"
        return
      endif
    endif


 endif
 
#ifdef USE_NVTX
    call nvtxEndRange
#endif

end subroutine gpmdcov_buildz
