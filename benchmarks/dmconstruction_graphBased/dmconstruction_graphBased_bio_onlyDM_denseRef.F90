!> Graph-based aproach driver.
!!
program gpsolve

  use bml
  use aux_mod
#ifdef DO_MPI
  use MPI
#endif

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

  call prg_progress_init()

  myRank = getMyRank() + 1

  if(myRank == 1)call prg_version()

  call prg_barrierParallel()

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
  call prg_barrierParallel()

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
  call prg_barrierParallel()

  sparsity = bml_get_sparsity(ham_bml,bioham%threshold)
  if(myRank == 1)write(*,*)"Sparsity Ham=",sparsity
  call prg_barrierParallel()

  ! Construct the graph out ot S^2 and apply threshold
  if(myRank == 1)write(*,*)"Build graph..."
  call bml_zero_matrix(BML_MATRIX_DENSE,bml_element_real,dp,norbs,norbs,g_bml)
  call bml_multiply_x2(over_bml,g_bml,bioham%threshold,trace)
  call bml_threshold(g_bml,gppar%threshold)
  call prg_barrierParallel()

  ! Orthogonalize Hamiltonian using overlap matrix
  if(myRank == 1)write(*,*)"Orthogonalize Hamiltonian..."
  call bml_zero_matrix(BML_MATRIX_DENSE,bml_element_real,dp,norbs,norbs,zmat_bml)
  call prg_buildzdiag(over_bml,zmat_bml,bioham%threshold,bioham%mdim,BML_MATRIX_DENSE,10)
  call prg_barrierParallel()
  if(myRank == 1)write(*,*)"prg_buildzdiag done..."

  call bml_zero_matrix(BML_MATRIX_DENSE,bml_element_real,dp,norbs,norbs,oham_bml)
  call prg_orthogonalize(ham_bml,zmat_bml,oham_bml,&
       bioham%threshold,bioham%bml_type,bioham%verbose)
  call prg_barrierParallel()
  if(myRank == 1)write(*,*)"Hamiltonian ready!"

  ! Construct the density matrix from diagonalization of full matrix to compare with
  call bml_zero_matrix(BML_MATRIX_DENSE,bml_element_real,dp,norbs,norbs,aux_bml)
  if(myRank == 1)write(*,*)"Diagonalize replicated H for reference..."
  mlsi = mls()
  call prg_build_density_T_Fermi(oham_bml,aux_bml,bioham%threshold, 0.1_dp, Ef)
  mlsreg = mls()-mlsi

  allocate(umat(norbs,norbs))

  ! Call graph-based solver
  if(myRank == 1)write(*,*)"Distributed solver..."
  Ef = 0.0_dp

  call bml_export_to_dense(g_bml, umat)
  call bml_deallocate(g_bml)
  call bml_zero_matrix(BML_MATRIX_ELLPACK,bml_element_real,dp,norbs,norbs,g_bml)
  call bml_import_from_dense(BML_MATRIX_ELLPACK, umat, g_bml, bioham%threshold, norbs)

  call bml_export_to_dense(oham_bml, umat)
  call bml_deallocate(oham_bml)
  call bml_zero_matrix(BML_MATRIX_ELLPACK,bml_element_real,dp,norbs,norbs,oham_bml)
  call bml_import_from_dense(BML_MATRIX_ELLPACK, umat, oham_bml, bioham%threshold, norbs)

  call bml_zero_matrix(BML_MATRIX_ELLPACK,bml_element_real,dp,norbs,norbs,rho_bml)

  if(myRank == 1)write(*,*)"threshold..."
  call bml_threshold(g_bml,gppar%threshold)
  mlsi = mls()
  if(myRank == 1)write(*,*)"prg_build_densityGP_T0..."
  call prg_build_densityGP_T0(oham_bml, g_bml, rho_bml, gppar%threshold, bndfil, Ef, gppar%numparts, 10)
  mlsgraph = mls()-mlsi

  if(myRank == 1)call bml_print_matrix("rhoGP",rho_bml,0,10,0,10)
  if(myRank == 1)call bml_print_matrix("rho",aux_bml,0,10,0,10)

  call bml_export_to_dense(rho_bml, umat)
  call bml_deallocate(rho_bml)
  call bml_zero_matrix(BML_MATRIX_DENSE,bml_element_real,dp,norbs,norbs,rho_bml)
  call bml_import_from_dense(BML_MATRIX_DENSE, umat, rho_bml, bioham%threshold, norbs)
  call bml_add(aux_bml,rho_bml,1.0d0,-1.0d0,threshold)

  if(myRank == 1)then
    write(*,*)"Number of replicas =",bioham%replicatex*bioham%replicatey*bioham%replicatez
    write(*,*)"Number of atoms =",syf%nats
    write(*,*)"Number of orbitals =",norbs
    write(*,*)"Time for graphSolver =",mlsgraph
    write(*,*)"Time for regularSolver =",mlsreg
    write(*,*)"SpeedupRho=",mlsreg/mlsgraph
    write(*,*)"ErrorRho=",bml_fnorm(aux_bml)/real(norbs*norbs,dp)
  endif


end program gpsolve
