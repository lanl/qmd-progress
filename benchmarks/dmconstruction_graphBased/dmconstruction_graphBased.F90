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
  integer ::  verbose, i, tnnz,nparts,vsize(2), myRank
  integer, allocatable :: xadj(:), adjncy(:), vector(:)
  type(bml_matrix_t) ::  ham_bml,rho_bml,rhos_bml,evects_bml,aux_bml,g_bml
  type(mham_type) ::  mham
  type(system_type) ::  sys
  type(system_type), allocatable    ::  syprt(:)
  real(dp) ::  threshold, bndfil, maxCH, pnorm=6, threshold_g
  real(dp), allocatable :: trace(:)
  real(dp), allocatable :: eigenvalues(:)
  real(dp) :: ef,sparsity,dec,mlsi,mlsii,smooth_maxCH, sumCubes
  type(graph_partitioning_t) ::  gpat
  integer, allocatable  ::  part(:), core_count(:), Halo_count(:,:),CH_count(:)

  ! Parsing input file
  call prg_parse_mham(mham,"input.in") !Reads the input for modelham

  ! Allocate bml matrices
  norbs=mham%norbs
  call bml_zero_matrix(mham%bml_type,bml_element_real,dp,norbs,norbs,ham_bml)
  call bml_zero_matrix(mham%bml_type,bml_element_real,dp,norbs,norbs,rho_bml)
  call bml_zero_matrix(mham%bml_type,bml_element_real,dp,norbs,norbs,aux_bml)
  call bml_zero_matrix(mham%bml_type,bml_element_real,dp,norbs,norbs,g_bml)

  ! Constructing the Hamiltonian
  seed = 1000
  verbose = 1
  threshold = 1.0d-20
  bndfil = 0.5d0
  Ef = 0.0
  nparts = 10

  call prg_twolevel_model(mham%ea, mham%eb, mham%dab, mham%daiaj, mham%dbibj, &
       &mham%dec, mham%rcoeff, mham%reshuffle, mham%seed, ham_bml, verbose)
  call bml_threshold(ham_bml, threshold)
  if(myRank == 1)call bml_print_matrix("ham_bml",ham_bml,0,10,0,10)
  sparsity = bml_get_sparsity(ham_bml, 1.0D-5)
  if(myRank == 1)write(*,*)"Sparsity Ham=",sparsity

  ! Construct the graph out ot H^2 and apply threshold
  threshold_g = 1.0d-20
  call bml_multiply_x2(ham_bml,g_bml,threshold_g,trace)
  call bml_threshold(g_bml, threshold_g)

  ! Call API 
  mlsi = mls()
  call prg_build_densityGP_T0(ham_bml, g_bml, rho_bml, threshold, bndfil, Ef, nparts)
  if(myRank == 1)write(*,*)"Total time graph =",mls()-mlsi
  
  ! Construct the density matrix from diagonalization of full matrix to compare with
  if(myRank == 1)mlsi = mls()
  call prg_build_density_T_Fermi(ham_bml,aux_bml,threshold, 0.1_dp, Ef)
  if(myRank == 1)write(*,*)"Total time full diag =",mls()-mlsi

  call bml_print_matrix("rhoGP",rho_bml,0,10,0,10)
  call bml_print_matrix("rho",aux_bml,0,10,0,10)
  call bml_add(aux_bml,rho_bml,1.0d0,-1.0d0,threshold)
  write(*,*)"|rhoGP-rho|",bml_fnorm(aux_bml)


end program gpsolve
