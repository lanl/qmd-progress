module prg_conjgrad_mod.F90

use bml 

implicit none 

integer, parameter :: dp = kind(1.0d0)

public :: prg_conjgrad

contains 

subroutine prg_conjgrad(x_bml, A_bml, b_bml, cg_tol, threshold)

implicit none 

type(bml_matrix_t), intent(in) :: A_bml, b_bml 
type(bml_matrix_t), intent(inout) :: x_bml
real(dp), intent(in) :: cg_tol, threshold 

type(bml_matrix_t) :: r_bml, p_bml, w_bml 

integer :: k 

end subroutine prg_conjgrad

end module prg_conjgrad_mod.F90 
