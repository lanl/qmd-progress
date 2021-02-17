  !> Finalize gpmd
  subroutine gpmdcov_Finalize()

    use gpmdcov_vars

    deallocate(gbnd)

    call bml_deallocate(ham_bml)
    call bml_deallocate(ham0_bml)
    !     call bml_deallocate(orthoh_bml)
    call bml_deallocate(g_bml)
    call bml_deallocate(rho_bml)
    call bml_deallocate(over_bml)
    call bml_deallocate(zmat_bml)

    ! Progress is done
    call prg_progress_shutdown()

  end subroutine gpmdcov_Finalize

