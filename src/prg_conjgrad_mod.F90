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
real(dp), allocatable :: alpha(:), beta(:)
integer :: N, M, k
real :: r_norm

bml_type = bml_get_type(h_bml)
N = bml_get_N(h_bml)
M = bml_get_M(h_bml)

allocate(alpha(N))
allocate(beta(N))

call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, r_bml)
call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, p_bml)
call bml_zero_matrix(bml_type, bml_element_real, dp, N, M, w_bml)

! r = b - Ax
call bml_multiply_AB(A_bml,x_bml,r_bml,threshold)
call bml_add(r_bml,b_bml,-1.0_dp,1.0_dp,threshold)
r = bml_fnorm(r_bml)
k = 0

do while (r .gt. cg_tol) 
        k = k+1
        if (k .eq. 1) then 
                call bml_copy(r_bml,p_bml)
        else
        endif 
enddo 

deallocate(alpha)
deallocate(beta)

call bml_deallocate(r_bml)
call bml_deallocate(p_bml)
call bml_deallocate(w_bml)

end subroutine prg_conjgrad

end module prg_conjgrad_mod.F90 
