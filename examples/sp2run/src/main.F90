!> High-level program to compute the density matrix with the SP2 method.
!!
program sp2run

  !BML lib.
  use bml

  !PROGRESS lib modes
  use prg_sp2_mod
  use prg_sp2parser_mod

  implicit none
  integer, parameter                ::  dp = kind(1.0d0)
  integer                           ::  norb
  type(sp2data_type)                ::  sp2
  type(bml_matrix_t)                ::  ham_bml,rho_bml

  write(*,*)"Reading input sp2 parser ..."

  !> parsing input file.
  call prg_parse_sp2(sp2,"input.in") !Reads the input sp2 parameters.

  !> Get Hamiltonian.
  norb=sp2%ndim
  call bml_zero_matrix(sp2%bml_type,bml_element_real,dp,norb,norb,ham_bml)
  call bml_zero_matrix(sp2%bml_type,bml_element_real,dp,norb,norb,rho_bml)
  call bml_read_matrix(ham_bml,"ham.mtx")

  !> Do sp2.
  write(*,*)sp2%flavor
  if(sp2%flavor .eq. "Basic") then
    call prg_sp2_basic(ham_bml,rho_bml,sp2%threshold,sp2%bndfil,sp2%minsp2iter,sp2%maxsp2iter &
         ,sp2%sp2conv,sp2%sp2tol,sp2%verbose)
  elseif(sp2%flavor .eq. "Alg1") then
    call prg_sp2_alg1(ham_bml,rho_bml,sp2%threshold,sp2%bndfil,sp2%minsp2iter,sp2%maxsp2iter &
         ,sp2%sp2conv,sp2%sp2tol,sp2%verbose)
  elseif(sp2%flavor .eq. "Alg2") then
    call prg_sp2_alg2(ham_bml,rho_bml,sp2%threshold,sp2%bndfil,sp2%minsp2iter,sp2%maxsp2iter &
         ,sp2%sp2conv,sp2%sp2tol,sp2%verbose)
  else
    write(*,*)"SP2 algorithm variant/flavor not recognized"
  endif

  call bml_print_matrix("rho_bml",rho_bml,0,6,0,6)

end program sp2run
