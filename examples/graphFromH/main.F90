!> Constructs a Hamiltonian and get the adjacency
!! matrix in bml format.
!!
program gpsolve

  !BML lib.
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

  ! Initialize progress MPI and get ranks
  call prg_progress_init()
   myRank = getMyRank() + 1

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
  call prg_twolevel_model(mham%ea, mham%eb, mham%dab, mham%daiaj, mham%dbibj, &
       &mham%dec, mham%rcoeff, mham%reshuffle, mham%seed, ham_bml, verbose)
  call bml_threshold(ham_bml, threshold)
  if(myRank == 1)call bml_print_matrix("ham_bml",ham_bml,0,10,0,10)
  sparsity = bml_get_sparsity(ham_bml, 1.0D-5)
  if(myRank == 1)write(*,*)"Sparsity Ham=",sparsity
  
  ! if(myRank == 1)mlsi = mls()
  !call prg_build_density_T_Fermi(ham_bml,aux_bml,threshold, 0.1_dp, Ef)
  !if(myRank == 1)write(*,*)"Total time full diag =",mls()-mlsi

  ! Construct the graph out ot H^2 and apply threshold
  threshold_g = 1.0d-8
  call bml_multiply_x2(ham_bml,g_bml,threshold_g,trace)
  call bml_threshold(g_bml, threshold_g)

  ! Doing graph partitioning using Metis
  allocate(xadj(norbs+1)) ! Adjacency in csr format

  tnnz = 0 
  do i=1,norbs
  !tnnz = tnnz + bml_get_row_bandwidth(ham_bml,i)
    tnnz = tnnz + bml_get_row_bandwidth(g_bml,i)
  enddo

  allocate(adjncy(tnnz+1))

  call bml_adjacency(g_bml, xadj, adjncy, 1)

  nparts = 32
  nnodes = norbs

  allocate(part(nnodes))
  allocate(core_count(nparts))
  allocate(CH_count(nparts))
  allocate(Halo_count(nparts, nnodes))

  call prg_metisPartition(gpat, nnodes, nnodes, xadj, adjncy, nparts, part, core_count,&
       CH_count, Halo_count, sumCubes, maxCH, smooth_maxCH, pnorm)

  write(*,*)"gpat = ",gpat%pname
  write(*,*)"totalParts = ",gpat%totalParts
  write(*,*)"totalNodes = ",gpat%totalNodes
  write(*,*)"partMin = ",gpat%localPartMin(1)
  write(*,*)"partMax = ",gpat%localPartMax(1)

  allocate(syprt(gpat%TotalParts))

  call prg_wait

  if(myRank == 1)mlsi = mls()

  ! Extract core halo indices from partitions
  do i=1,gpat%TotalParts
    call bml_matrix2submatrix_index(g_bml,&
         gpat%sgraph(i)%nodeInPart,gpat%nnodesInPart(i),&
         gpat%sgraph(i)%core_halo_index, &
         vsize,.true.)
    gpat%sgraph(i)%lsize = vsize(1)
    gpat%sgraph(i)%llsize = vsize(2)
    if(myRank == 1)write(*,*)"Core and halo",vsize(1),vsize(2)
  enddo

  ! Construct DM for every part
  do i = gpat%localPartMin(myRank), gpat%localPartMax(myRank)
  !  do i=1,gpat%TotalParts ! If no MPI
    call bml_matrix2submatrix_index(g_bml,&
         gpat%sgraph(i)%nodeInPart,gpat%nnodesInPart(i),&
         gpat%sgraph(i)%core_halo_index, &
         vsize,.true.)
    gpat%sgraph(i)%lsize = vsize(1)
    gpat%sgraph(i)%llsize = vsize(2)
    inorbs = vsize(1)
    if(myRank == 1)write(*,*)"Core and CH sizes:",gpat%sgraph(i)%llsize,gpat%sgraph(i)%lsize
    call bml_zero_matrix("dense",bml_element_real,dp,inorbs,inorbs,syprt(i)%estr%ham)
    if(allocated(vector))deallocate(vector)
    allocate(vector(gpat%sgraph(i)%lsize))
    vector(:) = gpat%sgraph(i)%core_halo_index(1:gpat%sgraph(i)%lsize)
    call bml_matrix2submatrix(ham_bml, syprt(i)%estr%ham, vector, &
         & inorbs)

    !Computing the density matrix with diagonalization
    mlsii = mls()
    call bml_zero_matrix("dense",bml_element_real,dp,inorbs,inorbs,syprt(i)%estr%rho)
    call prg_build_density_T_Fermi(syprt(i)%estr%ham,syprt(i)%estr%rho,0.0_dp, 0.1_dp, Ef)
    write(*,*)"Time for prg_build_density_T0",mls()-mlsii
    !call bml_print_matrix("syprt(i)%estr%rho",syprt(i)%estr%rho,0,5,0,5)
    call bml_submatrix2matrix(syprt(i)%estr%rho,rho_bml,vector,gpat%sgraph(i)%lsize,gpat%sgraph(i)%llsize,threshold)
  enddo
 
  ! Collect the small DM matrices into its full system corresponding 
  call prg_partOrdering(gpat)
  call prg_collectMatrixFromParts(gpat, rho_bml)
  call prg_wait

  if(myRank == 1)write(*,*)"Total time graph =",mls()-mlsi
  
  ! Construct the density matrix from diagonalization of full matrix to compare with
  if(myRank == 1)mlsi = mls()
  call prg_build_density_T_Fermi(ham_bml,aux_bml,threshold, 0.1_dp, Ef)
  if(myRank == 1)write(*,*)"Total time full diag =",mls()-mlsi

  call bml_print_matrix("rhoGP",rho_bml,0,10,0,10)
  call bml_print_matrix("rho",aux_bml,0,10,0,10)
  call bml_add(aux_bml,rho_bml,1.0d0,-1.0d0,threshold)
  write(*,*)"|rhoGP-rho|",bml_fnorm(aux_bml)

  call prg_progress_shutdown

end program gpsolve
