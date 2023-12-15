!> Graph-based aproach driver.
!!
program gpsolve

  use bml
  use aux_mod
  use MPI

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
  use ham_latte_mod
  use tbparams_latte_mod
  use ppot_latte_mod
  use prg_ptable_mod
  use prg_genz_mod
  use prg_nonortho_mod

  integer, parameter ::  dp = kind(1.0d0)
  character(2), allocatable         ::  TypeA(:,:), TypeB(:,:)
  character(3), allocatable         ::  intKind(:)
  integer                           ::  myRank, nel, norbs, ierr
  integer, allocatable              ::  hindex(:,:)
  real(dp)                          ::  bndfil, ef, mlsgraph, mlsi
  real(dp)                          ::  mlsreg, sparsity, threshold
  real(dp), allocatable             ::  onsitesH(:,:), onsitesS(:,:), origin(:), trace(:)
  type(bioham_type)                 ::  bioham
  type(bml_matrix_t)                ::  aux_bml, g_bml, ham_bml, oham_bml
  type(bml_matrix_t)                ::  over_bml, rho_bml, zmat_bml
  type(gppar_type)                  ::  gppar
  type(intpairs_type), allocatable  ::  intPairsH(:,:), intPairsS(:,:)
  type(ppot_type), allocatable      ::  ppot(:,:)
  type(system_type)                 ::  sy, syf
  type(tbparams_type)               ::  tb
  real(dp), allocatable :: umat(:,:)


  if(myRank == 1)call prg_version()

  ! Parsing input file
  call prg_parse_bioham(bioham,"input_bio.in") !Reads the input for modelham

  ! Parsing graph related parameters
  call prg_parse_gppar(gppar,"input_bio.in") !Reads the input for modelham

  ! Reading the system
  call prg_parse_system(sy,"prot.pdb")

  call prg_replicate_system(sy,syf,bioham%replicatex,bioham%replicatey,bioham%replicatez)

  call prg_destroy_system(sy)

  ! Center sytem inside the box and fold it by the lattice_vectors.
  allocate(origin(3)); origin = 0.0_dp
  call prg_translateandfoldtobox(syf%coordinate,syf%lattice_vector,origin)
  !call prg_barrierParallel()

  call prg_write_system(syf,"protf.pdb")

  ! Constructing the Hamiltonian
  call load_latteTBparams(tb,syf%splist,bioham%parampath)

  ! Get the mapping of the Hamiltonian index with the atom index
  allocate(hindex(2,syf%nats))

  ! Bond integrals parameters for LATTE Hamiltonian.
  call load_bintTBparamsH(syf%splist,tb%onsite_energ,&
       typeA,typeB,intKind,onsitesH,onsitesS,intPairsH,intPairsS,bioham%parampath)

  if(myRank == 1)call write_bintTBparamsH(typeA,typeB,&
       intKind,intPairsH,intPairsS,adjustl(trim(bioham%jobname))//"_mybondints.nonorth")

  ! Load Pair potentials for LATTE TB.
  call load_PairPotTBparams(bioham%parampath,syf%splist,ppot)

  call get_hindex(syf%spindex,tb%norbi,hindex,norbs)

  call bml_zero_matrix(BML_MATRIX_DENSE,bml_element_real,dp,norbs,norbs,ham_bml)
  call bml_zero_matrix(BML_MATRIX_DENSE,bml_element_real,dp,norbs,norbs,over_bml)

  call get_hsmat(ham_bml,over_bml,syf%coordinate,&
       syf%lattice_vector,syf%spindex,&
       tb%norbi,hindex,onsitesH,onsitesS,intPairsH,intPairsS,bioham%threshold)

  if(bioham%mdim == 0) bioham%mdim = norbs
  ! Get occupation based on last shell population.
  nel = sum(element_numel(syf%atomic_number))
  bndfil = nel/(2.0_dp*real(norbs,dp))

  call bml_threshold(ham_bml,bioham%threshold)

  if(myRank == 1)call bml_print_matrix("ham_bml",ham_bml,0,10,0,10)

  sparsity = bml_get_sparsity(ham_bml,bioham%threshold)
  write(*,*)"Sparsity Ham=",sparsity
  
  write(*,*)"Orthogonalize Hamiltonian..."
  call bml_zero_matrix(BML_MATRIX_DENSE,bml_element_real,dp,norbs,norbs,zmat_bml)
  call prg_buildzdiag(over_bml,zmat_bml,bioham%threshold,norbs,BML_MATRIX_DENSE,10)
  
  write(*,*)"Orthogonalize Hamiltonian..."
  call bml_zero_matrix(BML_MATRIX_DENSE,bml_element_real,dp,norbs,norbs,oham_bml)
  call prg_orthogonalize(ham_bml,zmat_bml,oham_bml,&
       bioham%threshold,bioham%bml_type,bioham%verbose)

      ! Construct the density matrix from diagonalization of full matrix to compare with
  call bml_zero_matrix(BML_MATRIX_DENSE,bml_element_real,dp,norbs,norbs,rho_bml)
  if(printRank() == 1)mlsi = mls()
  call prg_build_density_T_Fermi(ham_bml,rho_bml,threshold, 0.1_dp, Ef)
  timeReg = mls() - mlsi
  if(printRank() == 1)write(*,*)"Total time full diag =",mls()-mlsi


  call bml_write_matrix(oham_bml,"oham.mtx")
  write(*,*)"Norbs",Norbs
  write(*,*)"Nel",nel
  write(*,*)"bndfil",bndfil
  write(*,*)"Time for solving rho",timeReg

end program gpsolve
