!> Graph-based aproach driver.
!!
program biosolve

  use bml
  use aux_mod
  use prg_sp2_mod

  !PROGRESS lib modes
  use prg_modelham_mod
  use prg_system_mod
  use prg_timer_mod
  use prg_extras_mod
  use prg_parallel_mod
  use prg_progress_mod
  use prg_densitymatrix_mod
  use ham_latte_mod
  use tbparams_latte_mod
  use ppot_latte_mod
  use prg_ptable_mod
  use prg_genz_mod
  use prg_nonortho_mod

  integer, parameter ::  dp = kind(1.0d0)
  character(2), allocatable         ::  TypeA(:,:), TypeB(:,:)
  character(3), allocatable         ::  intKind(:)
  character(20)                     ::  arg
  integer                           ::  myRank, nel, norbs, nnz, nreps
  integer, allocatable              ::  hindex(:,:)
  real(dp)                          ::  bndfil, ef, mlssp2, mlsi
  real(dp)                          ::  mlsdiag, sparsity, threshold
  real(dp), allocatable             ::  onsitesH(:,:), onsitesS(:,:), origin(:), trace(:)
  type(bioham_type)                 ::  bioham
  type(bml_matrix_t)                ::  aux_bml, ham_bml, oham_bml
  type(bml_matrix_t)                ::  over_bml, rho_bml
  type(intpairs_type), allocatable  ::  intPairsH(:,:), intPairsS(:,:)
  type(ppot_type), allocatable      ::  ppot(:,:)
  type(system_type)                 ::  sy, syf
  type(tbparams_type)               ::  tb
  real(dp) ::  tol
  real(dp), allocatable :: eigenvalues(:)

#ifdef CRAY_SDK
  integer iargc
#endif

  call prg_progress_init()
  myRank = getMyRank() + 1

  ! Parsing input file
  call prg_parse_bioham(bioham,"input.in") !Reads the input for modelham
  nreps=1
  if (iargc().gt.1)then
    call getarg(2,arg)
    read(arg,*)nreps
  endif

  ! Reading the system
  call prg_parse_system(sy,"prot.pdb")

  call prg_replicate_system(sy,syf,bioham%replicatex,bioham%replicatey,bioham%replicatez)
  call prg_destroy_system(sy)

  ! Center sytem inside the box and fold it by the lattice_vectors.
  allocate(origin(3)); origin = 0.0_dp
  call prg_translateandfoldtobox(syf%coordinate,syf%lattice_vector,origin)

  call prg_write_system(syf,"protf.pdb")

  ! Constructing the Hamiltonian
  call load_latteTBparams(tb,syf%splist,bioham%parampath)

  ! Get the mapping of the Hamiltonian index with the atom index
  allocate(hindex(2,syf%nats))

  ! Bond integrals parameters for LATTE Hamiltonian.
  call load_bintTBparamsH(syf%splist,tb%onsite_energ,&
       typeA,typeB,intKind,onsitesH,onsitesS,intPairsH,intPairsS,bioham%parampath)

  call write_bintTBparamsH(typeA,typeB,&
       intKind,intPairsH,intPairsS,adjustl(trim(bioham%jobname))//"_mybondints.nonorth")

  ! Load Pair potentials for LATTE TB.
  call load_PairPotTBparams(bioham%parampath,syf%splist,ppot)

  call get_hindex(syf%spindex,tb%norbi,hindex,norbs)

  call bml_zero_matrix(bioham%bml_type,bml_element_real,dp,norbs,norbs,ham_bml)
  call bml_zero_matrix(bioham%bml_type,bml_element_real,dp,norbs,norbs,over_bml)

  call get_hsmat(ham_bml,over_bml,syf%coordinate,&
       syf%lattice_vector,syf%spindex,&
       tb%norbi,hindex,onsitesH,onsitesS,intPairsH,intPairsS,bioham%threshold)

  if(bioham%mdim == 0) bioham%mdim = norbs
  ! Get occupation based on last shell population.
  nel = sum(element_numel(syf%atomic_number))
  if(myRank == 1)write(*,*)"N electrons = ",nel
  bndfil = nel/(2.0_dp*real(norbs,dp))
  if(myRank == 1)write(*,*)"bndfil = ",bndfil

  call bml_threshold(ham_bml,bioham%threshold)

  if(myRank == 1)call bml_print_matrix("ham_bml",ham_bml,0,10,0,10)
  sparsity = bml_get_sparsity(ham_bml,bioham%threshold)
  if(myRank == 1)write(*,*)"Sparsity Ham=",sparsity

  call bml_zero_matrix(bioham%bml_type,bml_element_real,dp,norbs,norbs,aux_bml)

  call prg_buildzdiag(over_bml,aux_bml,bioham%threshold,bioham%mdim,bioham%bml_type)
  if(myRank == 1)call bml_print_matrix("zmat",aux_bml,0,10,0,10)
  call bml_deallocate(over_bml)

  call bml_zero_matrix(bioham%bml_type,bml_element_real,dp,norbs,norbs,oham_bml)
  call prg_orthogonalize(ham_bml,aux_bml,oham_bml,&
       bioham%threshold,bioham%bml_type,bioham%verbose)
  call bml_deallocate(ham_bml)

  call bml_zero_matrix(bioham%bml_type,bml_element_real,dp,norbs,norbs,rho_bml)
  ! Call SP2
  tol = 2.0D-5*norbs*bndfil
  do i =1, nreps
    mlsi = mls()
    call prg_sp2_alg1(oham_bml, rho_bml, bioham%threshold, bndfil, 15,100, "Rel", tol, 20)
    mlssp2 = mls()-mlsi
    if (printRank() .eq. 1)write(*,*)"Time_for_prg_sp2_alg1",mlssp2
  enddo

  call bml_deallocate(aux_bml)

  if(myRank == 1)write(*,*)"Solve dense eigenvalue problem..."
  call bml_zero_matrix(bioham%bml_type,bml_element_real,dp,norbs,norbs,aux_bml)
  allocate(eigenvalues(norbs))

  ! Construct the density matrix from diagonalization of full matrix to compare with
  do i =1, nreps
    mlsi = mls()
    call prg_build_density_T0(oham_bml,aux_bml,bioham%threshold, bndfil, eigenvalues)
    mlsdiag = mls()-mlsi
    if (printRank() .eq. 1)write(*,*)"Time_for_prg_build_density_T0",mlsdiag
  enddo

  sparsity = bml_get_sparsity(rho_bml,bioham%threshold)

  if(myRank == 1)call bml_print_matrix("rhoSP2",rho_bml,0,10,0,10)
  if(myRank == 1)call bml_print_matrix("rhoDIAG",aux_bml,0,10,0,10)
  threshold = 0.D0
  call bml_add(aux_bml,rho_bml,1.0d0,-1.0d0,threshold)

  if(myRank == 1)then
    write(*,*)"DM sparsity              = ",sparsity
    write(*,*)"Threshold                = ",bioham%threshold
    write(*,*)"Number of replicas       = ",bioham%replicatex*bioham%replicatey*bioham%replicatez
    write(*,*)"Number of atoms          = ",syf%nats
    write(*,*)"Number of orbitals       = ",norbs
    write(*,*)"Time for SP2             = ",mlssp2
    write(*,*)"Time for Diagonalization = ",mlsdiag
    write(*,*)"Speedup                  = ",mlsdiag/mlssp2
    write(*,*)"Error DM                 = ",bml_fnorm(aux_bml)/real(norbs*norbs,dp)
  endif

  call bml_deallocate(aux_bml)
  call bml_deallocate(rho_bml)

end program biosolve
