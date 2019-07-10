
!!
program main

  use bml

  !progress lib modes
  use prg_implicit_fermi_mod
  use prg_conjgrad_mod 
  !use hamiltonian_mod

  implicit none
   
  type(bml_matrix_t) :: A_bml, x_bml, b_bml
  real :: cgtol,threshold
  character(20) :: bml_type

  bml_type = "ellpack"
  cgtol = 10^-6
  threshold = 10^-9

  call bml_random_matrix(bml_type,bml_element_real,dp,10,3,A_bml)
  call bml_random_matrix(bml_type,bml_element_real,dp,10,3,b_bml)
  call bml_random_matrix(bml_type,bml_element_real,dp,10,3,x_bml)
  
  call prg_conjgrad_mod(A_bml,x_bml,b_bml,cgtol,threshold) 

  call bml_deallocate(A_bml)
  call bml_deallocate(b_bml)
  call bml_deallocate(x_bml)
end program main
