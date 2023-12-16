!> Graph-based aproach driver.
!!
program gpsolve

  use bml

  !PROGRESS lib modes
  use prg_modelham_mod
  use prg_system_mod
  use prg_timer_mod
  use prg_extras_mod
  use prg_partition_mod
  use prg_graph_mod
  use prg_parallel_mod
  use prg_progress_mod
  use prg_densitymatrix_mod
  use prg_subgraphloop_mod
  use prg_graphsolver_mod

  implicit none
  integer, parameter ::  dp = kind(1.0d0)
  integer ::  norbs,seed,nnodes,ipt,inorbs
  integer ::  verbose, i, tnnz,nparts,vsize(2)
  integer, allocatable :: xadj(:), adjncy(:), vector(:)
  type(bml_matrix_t) ::  ham_bml,rho_bml,rhos_bml,evects_bml,aux_bml,g_bml
  type(mham_type) ::  mham
  type(system_type) ::  sys
  type(system_type), allocatable    ::  syprt(:)
  real(dp) ::  threshold, bndfil, maxCH, pnorm=6, threshold_g, sparsityRho
  real(dp), allocatable :: trace(:)
  real(dp), allocatable :: eigenvalues(:)
  real(dp) :: ef,sparsity,dec,mlsi,mlsii,smooth_maxCH, sumCubes, timeReg, timeGP
  type(graph_partitioning_t) ::  gpat
  integer, allocatable  ::  part(:), core_count(:), Halo_count(:,:),CH_count(:)
  character(20) :: bml_type

  call prg_progress_init()

  if(printRank() == 1)call prg_version()

  norbs=  2162
  !norbs= 8648
  !norbs= 12972
  norbs= 19458
  !norbs = 29187
  bml_type = "ellpack"
  verbose = 1
  threshold = 1.0d-6
  threshold_g = 1.0d-2
  bndfil = 0.649398_dp
  Ef = 0.0_dp
  nparts = 4

  if(printRank() == 1)write(*,*)"threshold_g=",threshold_g

  ! Allocate bml matrices
  ! Reading the Hamiltonian
  call bml_zero_matrix(bml_type,bml_element_real,dp,norbs,norbs,ham_bml)
  call bml_read_matrix(ham_bml,"oham.mtx")
  write(*,*)"done reading matrix..."

  call bml_zero_matrix(bml_type,bml_element_real,dp,norbs,norbs,rho_bml)
  call bml_zero_matrix(bml_type,bml_element_real,dp,norbs,norbs,aux_bml)
  call bml_zero_matrix(bml_type,bml_element_real,dp,norbs,norbs,g_bml)

  call bml_threshold(ham_bml, threshold)
  if(printRank() == 1)call bml_print_matrix("ham_bml",ham_bml,0,10,0,10)
  sparsity = bml_get_sparsity(ham_bml, 1.0D-2)
  if(printRank() == 1)write(*,*)"Sparsity Ham=",sparsity

  ! Construct the density matrix from diagonalization of full matrix to compare with
  if(printRank() == 1)mlsi = mls()
  call prg_build_density_T_Fermi(ham_bml,aux_bml,threshold, 0.1_dp, Ef)
  timeReg = mls() - mlsi
  if(printRank() == 1)write(*,*)"Total time full diag =",mls()-mlsi
  if(printRank() == 1)call bml_print_matrix("aux_bml",aux_bml,0,10,0,10)

  call bml_copy(aux_bml,g_bml)
  call bml_threshold(g_bml, threshold_g)

  ! Call API
  mlsi = mls()
  call prg_build_densityGP_T0(ham_bml, g_bml, rho_bml, threshold, bndfil, Ef, nparts, 10)
  timeGP = mls() - mlsi
  sparsityRho = bml_get_sparsity(rho_bml, threshold)
  if(printRank() == 1)write(*,*)"Sparsity Rho=",sparsityRho

  if(printRank() == 1)call bml_print_matrix("rhoGP",rho_bml,0,10,0,10)
  if(printRank() == 1)call bml_print_matrix("rho",aux_bml,0,10,0,10)
  if(printRank() == 1)call bml_add(aux_bml,rho_bml,1.0d0,-1.0d0,threshold)
  if(printRank() == 1)write(*,*)"|rhoGP-rho|",bml_fnorm(aux_bml)
  if(printRank() == 1)write(*,*)"|rhoGP-rho|/((1-sparsityRho)*N2)",bml_fnorm(aux_bml)/((1.0_dp - sparsityRho)*norbs*norbs)
  if(printRank() == 1)write(*,*)"Time for solving rho GB", timeGP

  call  prg_shutdownParallel()

end program gpsolve
