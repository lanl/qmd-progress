module prg_sp2_tensorcore_mod 

  use, intrinsic :: iso_c_binding !C interface

  implicit none
  public :: prg_sp2_tensorcore_f , prg_sp2_tensorcore_C
 
  interface 

    subroutine prg_sp2_tensorcore_C(N,H,D,eps,bndfil,minsp2iter,maxsp2iter,&
              & sp2conv,idemtol,verbose) bind(C, name="prg_sp2_tensorcore")
      import :: C_PTR, C_INT, C_FLOAT, C_CHAR 
      type(C_PTR), value :: D
      type(C_PTR), value :: H
      real(C_FLOAT), value, intent(in) :: eps, idemtol, bndfil
      integer(C_INT), value, intent(in) :: N, minsp2iter, maxsp2iter, verbose
      character(C_CHAR), intent(in)  :: sp2conv(*)
    end subroutine prg_sp2_tensorcore_C
     
  end interface 

  contains 
  
  subroutine prg_sp2_tensorcore_f(N,H,D,eps,bndfil,minsp2iter,maxsp2iter,&
            & sp2conv,idemtol,verbose)
    integer(C_INT), intent(in) :: N, minsp2iter, maxsp2iter, verbose
    real(C_DOUBLE), target :: D(*)
    real(C_FLOAT), target, intent(in) :: H(*)
    real(C_FLOAT), intent(in) :: eps, bndfil, idemtol
    character(C_CHAR), intent(in) :: sp2conv
    
    !Call the interface
    write(*,*)"Entering prg_sp2_tensorcore_C ..."
    call prg_sp2_tensorcore_C(N,c_loc(H),c_loc(D),eps,bndfil,minsp2iter,maxsp2iter,sp2conv,idemtol,verbose)
 
  end subroutine prg_sp2_tensorcore_f
 
end module prg_sp2_tensorcore_mod

