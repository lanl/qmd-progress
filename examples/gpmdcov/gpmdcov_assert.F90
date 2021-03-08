module gpmdcov_assert_mod

 use prg_extras_mod

contains

  subroutine gpmdcov_assert_input(rank) 
    use gpmdcov_vars
    integer, intent(in) :: rank

    if(rank == 1)then 
    if(lt%method == "Diag" .and. lt%mucalctype == "FromParts")then 
      write(*,*)""
      write(*,*)"ERROR: Incompatible parameters set"
      write(*,*)"Keyword Method= Diag is incompatible with MuCalcType= FromParts"
      write(*,*)"Set Method= DiagEf instead"
      write(*,*)""
      stop 
    endif 

    if(lt%method == "Diag" .and. lt%mucalctype == "Combined")then
      write(*,*)""
      write(*,*)"ERROR: Incompatible parameters set"
      write(*,*)"Keyword Method= Diag is incompatible with MuCalcType= Combined"
      write(*,*)"Set Method= DiagEf instead"
      write(*,*)""
      stop
    endif
 
    if(lt%method == "Diag" .and. lt%entropy .eqv. .true.)then
      write(*,*)""
      write(*,*)"ERROR: Incompatible parameters set"
      write(*,*)"Keyword Method= Diag is incompatible with Entropy= T"
      write(*,*)"Set Method= DiagEf instead"
      write(*,*)""
      stop
    endif
    endif

  end subroutine gpmdcov_assert_input

end module gpmdcov_assert_mod
