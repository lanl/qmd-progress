!> High-level program to construct a model Hamiltonian
!!
program hmodel

  !BML lib.
  use bml

  !PROGRESS lib modes
  use prg_modelham_mod
  use prg_system_mod
  use prg_densitymatrix_mod
  use prg_dos_mod
  use prg_sp2_mod
  use prg_timer_mod
  use prg_extras_mod

  implicit none
  integer, parameter ::  dp = kind(1.0d0)
  integer ::  norbs,seed
  integer ::  verbose
  type(bml_matrix_t) ::  ham_bml,rho_bml,rhos_bml,evects_bml,aux_bml
  type(mham_type) ::  mham
  type(system_type) ::  sys
  real(dp) ::  threshold, bndfil
  real(dp), allocatable :: trace(:)
  real(dp), allocatable :: eigenvalues(:)
  real(dp) :: ef,sparsity,dec,mlsi

  !Parsing input file.
  call prg_parse_mham(mham,"input.in") !Reads the input for modelham

  norbs=mham%norbs
  call bml_zero_matrix(mham%bml_type,bml_element_real,dp,norbs,norbs,ham_bml)
  call bml_zero_matrix(mham%bml_type,bml_element_real,dp,norbs,norbs,rho_bml)
  call bml_zero_matrix(mham%bml_type,bml_element_real,dp,norbs,norbs,evects_bml)
  call bml_zero_matrix(mham%bml_type,bml_element_real,dp,norbs,norbs,rhos_bml)
  call bml_zero_matrix(mham%bml_type,bml_element_real,dp,norbs,norbs,aux_bml)

  seed = 1000
  verbose = 1
  threshold = 1.0d-5
  bndfil = 0.5d0
  
  allocate(eigenvalues(norbs))
  
  !Constructng the Hamiltonian
  call prg_twolevel_model(mham%ea, mham%eb, mham%dab, mham%daiaj, mham%dbibj, &
  &mham%dec, mham%rcoeff, ham_bml, verbose)

  call bml_print_matrix("ham_bml",ham_bml,0,10,0,10)
  
  sparsity = bml_get_sparsity(ham_bml, 1.0D-5)  
  write(*,*)"Sparsity Ham=",sparsity

  !Computing the density matrix with diagonalization
  mlsi = mls()
  call prg_build_density_T0(ham_bml, rho_bml, threshold, bndfil, eigenvalues)
  write(*,*)"Time for prg_build_density_T0",mls()-mlsi

  sparsity = bml_get_sparsity(rho_bml, 1.0D-5)  
  write(*,*)"Sparsity Rho=",sparsity
  
  !Getting the fermi level
  ef = (eigenvalues(int(norbs/2)+1) + eigenvalues(int(norbs/2)))/2 
  eigenvalues = eigenvalues - ef 
  
  !Writting the total DOS
  call prg_write_tdos(eigenvalues, 0.05d0, 10000, -20.0d0, 20.0d0, "tdos.dat")

  !Solving for Rho using SP2
  mlsi = mls()
  call prg_sp2_basic(ham_bml,rhos_bml,threshold,bndfil,15,100 &
       ,"Rel",1.0D-5,20)
  write(*,*)"Time for prg_sp2_alg1",mls()-mlsi
  call bml_print_matrix("rho_bml",rho_bml,0,10,0,10)
  call bml_print_matrix("rhos_bml",rhos_bml,0,10,0,10)

  call bml_copy(rhos_bml,aux_bml)
  call bml_add(aux_bml,rho_bml,1.0d0,-1.0d0,threshold)
  write(*,*)"|DM_sp2-DM_diag|",bml_fnorm(aux_bml)

  call bml_multiply(rhos_bml, rhos_bml, aux_bml, 0.5_dp, 0.0_dp, threshold)
  call bml_print_matrix("rhos_bml^2",aux_bml,0,10,0,10)
  call bml_add(aux_bml,rhos_bml,1.0d0,-1.0d0,threshold)
  write(*,*)"|DM_sp2-DM_sp2^2|",bml_fnorm(aux_bml)
  
  call bml_multiply(ham_bml,rhos_bml,aux_bml,1.0_dp,0.0_dp,threshold)
  call bml_multiply(rhos_bml,ham_bml,aux_bml,1.0_dp,-1.0_dp,threshold)
  write(*,*)"|DM_sp2*H-H*DM_sp2|",bml_fnorm(aux_bml)

  call bml_write_matrix(ham_bml, "hamiltonian.mtx")
  
end
