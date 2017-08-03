!!comments to describe what it does


module prg_evolvedens_mod
  use bml
  implicit none

  


  private  !Everything is private by default 
  integer, parameter :: dp = kind(1.0d0)


  public::prg_get_sparsity_cplxmat, prg_get_sparsity_realmat

contains 



 subroutine prg_get_sparsity_cplxmat(matrix_type,element_type,sparsity_threshold,a_dense)
   implicit none
    character(len=*), intent(in) :: matrix_type ,element_type
    complex(8),intent(in)::a_dense(:,:)
    type(bml_matrix_t):: a
    real(8),intent(in)::sparsity_threshold
    real(8)::convert_threshold
    integer asize
    character(len=20):: test_type

    asize=SIZE(a_dense,1)
    convert_threshold=1.0d0
    !sparsity_threshold = 1.0d-5
    call bml_zero_matrix(matrix_type, element_type,1,asize,asize, a)
    call bml_convert_from_dense(matrix_type, a_dense, a, convert_threshold)
    write(*,*)"thr,sparsity ",sparsity_threshold, bml_get_sparsity(a, sparsity_threshold)

    test_type=bml_get_type(a)
    write(*,*)"matrix type  ", test_type


 end subroutine prg_get_sparsity_cplxmat


 subroutine prg_get_sparsity_realmat(matrix_type,element_type,sparsity_threshold,a_dense)
   implicit none
    character(len=*), intent(in) :: matrix_type ,element_type
    real(8),intent(in)::a_dense(:,:)
    type(bml_matrix_t):: a
    real(8),intent(in)::sparsity_threshold
    real(8)::convert_threshold
    integer asize


    asize=SIZE(a_dense,1)
    convert_threshold=1.0d0
    !sparsity_threshold = 1.0d-5
    call bml_zero_matrix(matrix_type, element_type,1,asize,asize, a)
    call bml_convert_from_dense(matrix_type, a_dense, a, convert_threshold)
    write(*,*)"thr,sparsity ", sparsity_threshold,bml_get_sparsity(a, sparsity_threshold)


 end subroutine prg_get_sparsity_realmat

end module prg_evolvedens_mod
